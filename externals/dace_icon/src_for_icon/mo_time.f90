!
!+ Define derived type to hold time+date information and operations thereon
!
MODULE mo_time
!
! Description:
!   Define derived type t_time to hold time+date information.
!   (seconds to centuries)
!   Define operations on derived type t_time:
!     - initialisation
!     - comparison  (<,>,==,..)
!     - arithmetics (+,-,*,max,min,...)
!     - derive various integer   representations (yyyymmddhh, ...)
!     - derive various character representations (yyyymmddhh, ...)
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
! V1_4         2009/03/26 Harald Anlauf
!  Add SEQUENCE attribute to t_time
! V1_5         2009/05/25 Andreas Rhodin
!  new function ihhhmm (integer representation: arbitrary number of h+min)
! V1_8         2009/12/09 Harald Anlauf
!  minval_time: add optional argument mask
! V1_9         2010/04/20 Harald Anlauf
!  hours: make conversion numerically symmetric w.r.t. time=0
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Andreas Rhodin
!  new subroutine solar (azimuth and zenith angle)
! V1_15        2011/12/06 Harald Anlauf
!  mo_time: declare more subroutines as pure / elemental
! V1_19        2012-04-16 Andreas Rhodin
!  cosmetic changes
! V1_22        2013-02-13 Andreas Rhodin
!  new conversion routines t_time <-> 'modified Julan date' (STD operator)
! V1_23        2013-03-26 Andreas Rhodin
!  new functions: iddhhmmss, ddhhmmss, chhhmm
! V1_42        2015-06-08 Andreas Rhodin
!  derive t_time from GRADS date/time string
! V1_45        2015-12-15 Andreas Rhodin
!  fix function chhhmm
! V1_47        2016-06-06 Andreas Rhodin
!  new specific function integer_times_time
! V1_48        2016-10-06 Harald Anlauf
!  new function localdate: Derive local (mean solar) date/time for longitude
!
! Code Description:
! Language: Fortran 95.
! Software Standards:
!
! Authors:
! Andreas Rhodin  DWD  2001-2008
! Harald Anlauf   DWD  2007-2008
!------------------------------------------------------------------------------
  !=============
  ! Modules used
  !=============
  use mo_kind,      only: wp, dp, i8
  use mo_exception, only: finish
  use mo_mpi_dace,  only: p_bcast
  implicit none
  private
  !-----------------------------------------------
  ! Data type 't_time' to hold time yyyymmddhhmmss
  !-----------------------------------------------
  public :: t_time
  public :: zero_time
  public :: invalid_time
  !--------------------------------------------
  ! derive integer representation from 't_time'
  !--------------------------------------------
  public :: iyyyy      ! years
  public :: imm        ! months
  public :: idd        ! days
  public :: ihh        ! hours
  public :: imi        ! minutes
  public :: iss        ! seconds
  public :: i_time     ! yyyy, mm, dd, hh, mi, ss
  public :: iddmmyyhh
  public :: ihhmmss
  public :: ihhmm
  public :: ihhhmm
  public :: ihhhmmss
  public :: iyyyymmdd
  public :: iyyddd
  public :: days_year  ! days in current year
  public :: frac_year  ! elapsed fraction of current year
  !---------------------------------------------
  ! derive real(wp) representation from 't_time'
  !---------------------------------------------
  public :: days
  public :: hours
  public :: minutes
  public :: seconds
  public :: mjd     ! modified Julian date (GNSSGB operator)
  !----------------------------------------------
  ! derive character representation from 't_time'
  !----------------------------------------------
  public :: cddmmyyhh
  public :: chhmmss
  public :: chhmm
  public :: chhhmm
  public :: chhhmmss
  public :: cddhhmmss
  public :: cdddhhmmss
  public :: cyyyymm
  public :: cyyyymmdd
  public :: cyymmdd
  public :: cyyyymmddhh
  public :: cyyyymmddhhmm
  public :: cyyyymmddhhmmss
  public :: cyyddd
  public :: chhzddmmmyyyy  ! GRADS initial time
  public :: cvvkk          ! GRADS time increment
  public :: cdate          ! 0000-00-00
  public :: ctime          ! 00:00:00
  public :: date_tmpl
  !--------------------------------------------
  ! derive 't_time' from integer representation
  !--------------------------------------------
  public :: time_ddmmyyhh
  public :: time_yyyymmddhh
  public :: time_yyyymmdd_hhmmss
  public :: init_time
  public :: now                  ! returns actual date
  !----------------------------------------------
  ! derive 't_time' from character representation
  !----------------------------------------------
  public :: time_cyyyymmddhhmmss
  public :: time_cyyyymmddhhmm
  public :: time_cyyyymmddhh
  public :: time_c
  public :: time_grads        ! from GRADS date/time string
  !------------------------------------------
  ! derive 't_time' from modified Julian date
  !------------------------------------------
  public :: time_mjd
  !-----------------------
  ! operations on 't_time'
  !-----------------------
  public :: print
  public :: max
  public :: min
  public :: minval
  public :: maxval
  public :: operator (+)
  public :: operator (-)
  public :: operator (*)
  public :: operator (/)
  public :: operator (>)
  public :: operator (<)
  public :: operator (>=)
  public :: operator (<=)
  public :: operator (==)
  public :: operator (/=)
  public :: sign
  public :: solar         ! calculate solar parameters: azimuth, zenith angle
  public :: p_bcast       ! MPI broadcast
  public :: localdate     ! Derive local (mean solar) date/time for longitude
!==============================================================================
  !-----------------------------------------------
  ! Data type 't_time' to hold time yyyymmddhhmmss
  !-----------------------------------------------
  type t_time
    SEQUENCE
    integer :: days = 0 ! yyyymmdd
    integer :: secs = 0 ! hhmmss
  end type t_time
  !------------------------
  ! Zero time, invalid time
  !------------------------
  type (t_time) ,parameter :: zero_time    = t_time (      0 ,      0 )
  type (t_time) ,parameter :: invalid_time = t_time (-huge(0),-huge(0))
!==============================================================================
  interface operator (+)
    module procedure add_times
  end interface
!------------------------------------------------------------------------------
  interface operator (-)
    module procedure subtr_times
    module procedure minus_time
  end interface
!------------------------------------------------------------------------------
  interface operator (*)
    module procedure real_times_time
    module procedure integer_times_time
  end interface
!------------------------------------------------------------------------------
  interface operator (/)
    module procedure time_over_real
  end interface
!------------------------------------------------------------------------------
  interface operator (>)
    module procedure time_gt_time
  end interface
!------------------------------------------------------------------------------
  interface operator (<)
    module procedure time_lt_time
  end interface
!------------------------------------------------------------------------------
  interface operator (>=)
    module procedure time_ge_time
  end interface
!------------------------------------------------------------------------------
  interface operator (<=)
    module procedure time_le_time
  end interface
!------------------------------------------------------------------------------
  interface operator (==)
    module procedure time_eq_time
  end interface
!------------------------------------------------------------------------------
  interface operator (/=)
    module procedure time_ne_time
  end interface
!------------------------------------------------------------------------------
  interface sign
    module procedure sign_int_time
  end interface
!------------------------------------------------------------------------------
  interface print
    module procedure print_time
  end interface
!------------------------------------------------------------------------------
  interface max
    module procedure max_time
  end interface
!------------------------------------------------------------------------------
  interface minval
    module procedure minval_time_masked
  end interface
!------------------------------------------------------------------------------
  interface maxval
    module procedure maxval_time
  end interface
!------------------------------------------------------------------------------
  interface min
    module procedure min_time
  end interface
!------------------------------------------------------------------------------
  interface p_bcast
    module procedure bcast_time
  end interface
!------------------------------------------------------------------------------
  interface init_time
    module procedure init_time
    module procedure init_times
    module procedure init_time_from_date
  end interface init_time
!==============================================================================
  integer          ,parameter :: mjd_offset = 2400001
  character(len=3) ,parameter :: mmm(12)    = &
  (/'jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec'/)
!==============================================================================
  integer,          parameter :: n_tmpl_forms = 4
  character(len=20)           :: date_tmpl_forms(n_tmpl_forms)
  data date_tmpl_forms(1) /'_YYMMDDHH_'/
  data date_tmpl_forms(2) /'_YYYYMMDDHH_'/
  data date_tmpl_forms(3) /'_YYYYMMDDHHMM_'/
  data date_tmpl_forms(4) /'_YYYYMMDDHHMMSS_'/
!==============================================================================
contains
!==============================================================================
  subroutine init_time (time, yyyy, mo, dd, hh, mi, ss, yyyymmdd, hhmmss)
  type (t_time) ,intent(out)          :: time
  integer       ,intent(in) ,optional :: yyyy
  integer       ,intent(in) ,optional :: mo
  integer       ,intent(in) ,optional :: dd
  integer       ,intent(in) ,optional :: hh
  integer       ,intent(in) ,optional :: mi
  integer       ,intent(in) ,optional :: ss
  integer       ,intent(in) ,optional :: yyyymmdd
  integer       ,intent(in) ,optional :: hhmmss
    integer :: y,d,m
    time% days = 0
    if (present (yyyy) .or. present (mo)) then
      if (.not. present (yyyy)  .or. &
          .not. present (mo) .or. &
          .not. present (dd)) then
!         write(0,*) 'init_time: year, month and day must all be present.'
!         stop
        call finish ('init_time','year, month and day must all be present.')
      endif
      time% days = julday (mo, dd, yyyy)
    else if (present( yyyymmdd )) then
      y =      yyyymmdd / 10000
      m = mod (yyyymmdd /   100 ,100)
      d = mod (yyyymmdd         ,100)
      time% days = julday (m, d, y)
    else
      if (present (dd)) time% days = dd
    endif
    time% secs = 0
    if(present( hh     )) time% secs = time% secs + hh * 3600
    if(present( mi     )) time% secs = time% secs + mi *   60
    if(present( ss     )) time% secs = time% secs + ss
    time% days = time% days + time% secs / 86400
    time% secs =        mod ( time% secs , 86400 )
    if(present( hhmmss )) then
      time% secs = mod (hhmmss / 10000 , 24) * 3600 &
                 + mod (hhmmss /   100 ,100) *   60 &
                 + mod (hhmmss        , 100)
      time% days = time% days + hhmmss / 240000
    endif
  end subroutine init_time
!------------------------------------------------------------------------------
  subroutine init_times (time, yyyy, mo, dd, hh, mi, ss, yyyymmdd, hhmmss)
  type (t_time) ,intent(out)          :: time (:)
  integer       ,intent(in) ,optional :: yyyy
  integer       ,intent(in) ,optional :: mo
  integer       ,intent(in) ,optional :: dd
  integer       ,intent(in) ,optional :: hh
  integer       ,intent(in) ,optional :: mi
  integer       ,intent(in) ,optional :: ss
  integer       ,intent(in) ,optional :: yyyymmdd
  integer       ,intent(in) ,optional :: hhmmss
    if (size(time) == 0) return
    call init_time (time(1), yyyy, mo, dd, hh, mi, ss, yyyymmdd, hhmmss)
    time(2:) = time(1)
  end subroutine init_times
!------------------------------------------------------------------------------
  subroutine init_time_from_date (time, date, status)
    type(t_time) ,intent(out)           :: time
    integer(i8)  ,intent(in)            :: date       ! yyyymmddhhmmss
    integer      ,intent(out), optional :: status     ! Error status
    !------------------------------------------------------
    ! Derive t_time from date represented as yyyymmddhhmmss
    !------------------------------------------------------
    integer :: yyyy, mo, dd, hh, mi, ss, stat
    if (date > 99991231235959_i8 .or. date < 10000000000000_i8) then
       stat = -1
    else
       stat = 0
    end if
    if (present (status)) then
       status = stat
    else if (stat /= 0) then
       call finish ('init_time','yyyymmddhhmmss out of range')
    end if
    yyyy =      date / 10000000000_i8
    mo   = mod (date / 100000000_i8   ,100_i8)
    dd   = mod (date / 1000000_i8     ,100_i8)
    hh   = mod (date / 10000_i8       ,100_i8)
    mi   = mod (date / 100_i8         ,100_i8)
    ss   = mod (date                  ,100_i8)
    call init_time (time, yyyy, mo, dd, hh, mi, ss)
  end subroutine init_time_from_date
!------------------------------------------------------------------------------
  function now (local)
  type (t_time)     :: now   ! actual date
  logical, optional :: local ! set .true. for local time instead of UTC
  !--------------------
  ! returns actual date
  !--------------------
    integer       :: i(8)
    type (t_time) :: diff
    logical       :: utc
    call date_and_time (values=i)
    call init_time (now, i(1),i(2),i(3),i(5),i(6),i(7))
    utc = .true.
    if (present(local)) utc = .not. local
    if (utc .and. i(4)/=0) then ! difference to UTC
      call init_time (diff, mi=i(4))
      now = now - diff
    end if
  end function now
!------------------------------------------------------------------------------
  elemental function time_mjd (mjd) result (time)
  type (t_time)             :: time
  real(dp)      ,intent(in) :: mjd
    time% days =        mjd_offset + int(mjd)
    time% secs = nint( (mjd        - int(mjd)) * 86400)
  end function time_mjd
!------------------------------------------------------------------------------
  function time_ddmmyyhh (ddmmyyhh) result (time)
  type (t_time)              :: time
  integer       ,intent(in)  :: ddmmyyhh
    integer :: dd, mm, yy, hh
    dd =      ddmmyyhh / 1000000
    mm = mod (ddmmyyhh / 10000  ,100)
    yy = mod (ddmmyyhh / 100    ,100)
    hh = mod (ddmmyyhh          ,100)
    if((mm/=0).and.(dd/=0)) then
      if(yy > 50) then
        yy=yy+1900
      else
        yy=yy+2000
      endif
      time% days = julday (mm, dd, yy)
    else if ((mm==0).and.(yy==0)) then
      time% days = dd
    else
      call finish ('time_ddmmyyhh','invalid input')
    endif
    time% secs = hh * 3600
  end function time_ddmmyyhh
!------------------------------------------------------------------------------
  elemental  function time_yyyymmddhh (yyyymmddhh) result (time)
  type (t_time)           :: time
  integer      ,intent(in):: yyyymmddhh
    integer :: yyyy, mm, dd, hh
    hh   = mod (yyyymmddhh        ,100)
    dd   = mod (yyyymmddhh/100    ,100)
    mm   = mod (yyyymmddhh/10000  ,100)
    yyyy =      yyyymmddhh/1000000
    time% days = julday (mm, dd, yyyy)
    time% secs = 3600*hh
 end function time_yyyymmddhh
!------------------------------------------------------------------------------
  elemental  function time_yyyymmdd_hhmmss (yyyymmdd, hhmmss) result (time)
  type (t_time)                        :: time
  integer       ,intent(in)            :: yyyymmdd
  integer       ,intent(in)  ,optional :: hhmmss
    integer :: yyyy, mm, dd
    integer :: hh, mi, ss
    dd   = mod (yyyymmdd      ,100)
    mm   = mod (yyyymmdd/100  ,100)
    yyyy =      yyyymmdd/10000
    time% days = julday (mm, dd, yyyy)
    if(present(hhmmss)) then
      ss = mod (hhmmss      ,100)
      mi = mod (hhmmss/100  ,100)
      hh =      hhmmss/10000
      time% secs = 3600*hh+60*mi+ss
    endif
 end function time_yyyymmdd_hhmmss
!------------------------------------------------------------------------------
  elemental function time_cyyyymmddhhmmss (cyyyymmddhhmmss) result (time)
  type (t_time)                :: time
  character(len=14),intent(in) :: cyyyymmddhhmmss
    integer :: yyyymmdd
    integer :: hhmmss
    read(cyyyymmddhhmmss,'(i8,i6)') yyyymmdd,hhmmss
    time   = time_yyyymmdd_hhmmss (yyyymmdd, hhmmss)
  end function time_cyyyymmddhhmmss
!------------------------------------------------------------------------------
  elemental function time_cyyyymmddhhmm (cyyyymmddhhmm) result (time)
  type (t_time)                :: time
  character(len=12),intent(in) :: cyyyymmddhhmm
    integer :: yyyymmdd
    integer :: hhmm
    integer :: hhmmss
    read(cyyyymmddhhmm,'(i8,i4)') yyyymmdd,hhmm
    hhmmss = hhmm * 100
    time   = time_yyyymmdd_hhmmss (yyyymmdd, hhmmss)
  end function time_cyyyymmddhhmm
!------------------------------------------------------------------------------
  elemental function time_cyyyymmddhh (cyyyymmddhh) result (time)
  type (t_time)                :: time
  character(len=10),intent(in) :: cyyyymmddhh
    integer :: yyyymmdd
    integer :: hhmmss
    integer :: hh
    read(cyyyymmddhh,'(i8,i2)') yyyymmdd,hh
    hhmmss = hh * 10000
    time   = time_yyyymmdd_hhmmss (yyyymmdd, hhmmss)
  end function time_cyyyymmddhh
!------------------------------------------------------------------------------
  elemental function time_c (str) result (t)
  type (t_time)               :: t
  character(len=*),intent(in) :: str
    if (verify(trim(str),'0123456789') > 0) then
      t = invalid_time
      return
    end if
    select case (len_trim(str))
    case(4)
      t = time_cyyyymmddhh(trim(str)//'010100')
    case(6)
      t = time_cyyyymmddhh(trim(str)//'0100')
    case(8)
      t = time_cyyyymmddhh(trim(str)//'00')
    case(10)
      t = time_cyyyymmddhh(trim(str))
    case(12)
      t = time_cyyyymmddhhmm(trim(str))
    case(14)
      t = time_cyyyymmddhhmmss(trim(str))
    case default
      t = invalid_time
    end select
  end function time_c
!==============================================================================
  elemental function iddmmyyhh (time)
  integer                   :: iddmmyyhh
  type (t_time) ,intent(in) :: time
    integer :: hh, yy, mm, dd
    hh =         time% secs/3600
    call caldat (time% days, mm, dd, yy)
    yy = mod(yy,100)
    iddmmyyhh = dd*1000000+mm*10000+yy*100+hh
  end function iddmmyyhh
!------------------------------------------------------------------------------
  elemental function cddmmyyhh (time)
  character(len=8)             :: cddmmyyhh
  type (t_time)    ,intent(in) :: time
    call i_to_c(cddmmyyhh,iddmmyyhh (time))
  end function cddmmyyhh
!------------------------------------------------------------------------------
  subroutine i_time (time, yyyy, mm, dd, hh, mi, ss)
  type (t_time) ,intent(in)  :: time
  integer       ,intent(out) :: yyyy, mm, dd, hh, mi,ss
    ss = mod (time% secs     ,60)
    mi = mod (time% secs/60  ,60)
    hh =      time% secs/3600
    call caldat (time% days, mm, dd, yyyy)
  end subroutine i_time
!------------------------------------------------------------------------------
  elemental function ihhmmss (time)
  integer                   :: ihhmmss
  type (t_time) ,intent(in) :: time
    integer :: hh, mm, ss
    ss = mod (time% secs     ,60)
    mm = mod (time% secs/60  ,60)
    hh =      time% secs/3600
    ihhmmss = 10000*hh+100*mm+ss
  end function ihhmmss
!------------------------------------------------------------------------------
  elemental function iddhhmmss (time)
  integer                   :: iddhhmmss
  type (t_time) ,intent(in) :: time
    integer :: hh, mm, ss, dd
    ss = mod (time% secs     ,60)
    mm = mod (time% secs/60  ,60)
    hh =      time% secs/3600
    dd =      time% days
    iddhhmmss = 1000000*dd+10000*hh+100*mm+ss
  end function iddhhmmss
!------------------------------------------------------------------------------
  elemental function ihhmm (time)
  integer                   :: ihhmm
  type (t_time) ,intent(in) :: time
    integer :: hh, mm
    mm = mod (time% secs/60  ,60)
    hh =      time% secs/3600
    ihhmm = 100*hh+mm
  end function ihhmm
!------------------------------------------------------------------------------
  elemental function ihhhmm (time)
  integer                   :: ihhhmm
  type (t_time) ,intent(in) :: time
    integer :: hh, mm
    mm = mod (time% secs/60  ,60)
    hh =      time% secs/3600
    ihhhmm = 100*hh+mm + 2400*time%days
  end function ihhhmm
!------------------------------------------------------------------------------
  elemental function ihhhmmss (time)
    integer                  :: ihhhmmss
    type(t_time) ,intent(in) :: time
    integer :: hh, mm, ss
    ss = mod (time% secs     ,60)
    mm = mod (time% secs/60  ,60)
    hh =      time% secs/3600
    ihhhmmss = 10000*hh + 100*mm + ss + 240000*time%days
  end function ihhhmmss
!------------------------------------------------------------------------------
  elemental function ihh (time)
  integer                   :: ihh
  type (t_time) ,intent(in) :: time
    ihh =      time% secs/3600
  end function ihh
!------------------------------------------------------------------------------
  elemental function imi (time)
  integer                   :: imi
  type (t_time) ,intent(in) :: time
    imi = mod (time% secs/60, 60)
  end function imi
!------------------------------------------------------------------------------
  elemental function iss (time)
  integer                   :: iss
  type (t_time) ,intent(in) :: time
    iss = mod (time% secs, 60)
  end function iss
!------------------------------------------------------------------------------
  elemental function days (time)
  real(wp)                  :: days
  type (t_time) ,intent(in) :: time
    days = time% days + time% secs / 86400._wp
  end function days
!------------------------------------------------------------------------------
  elemental function mjd (time)
  real(dp)                  :: mjd
  type (t_time) ,intent(in) :: time
    mjd = (time% days - mjd_offset) + real(time% secs,dp) / 86400
  end function mjd
!------------------------------------------------------------------------------
  elemental function hours (time)
  real(wp)                  :: hours
  type (t_time) ,intent(in) :: time
    type(t_time) :: abs_t   ! "Absolute value" of argument
    if (time% days >= 0) then
       hours =     time% days * 24._wp +  time% secs / 3600._wp
    else
       abs_t = - time
       hours = - (abs_t% days * 24._wp + abs_t% secs / 3600._wp)
    end if
  end function hours
!------------------------------------------------------------------------------
  elemental function minutes (time)
  real(wp)                  :: minutes
  type (t_time) ,intent(in) :: time
    minutes = time% days * 1440._wp + time% secs / 60._wp
  end function minutes
!------------------------------------------------------------------------------
  elemental function seconds (time)
  real(wp)                  :: seconds
  type (t_time) ,intent(in) :: time
    seconds = time% days * 86400._wp + time% secs
  end function seconds
!------------------------------------------------------------------------------
  elemental function iyyyymmdd (time)
  integer                   :: iyyyymmdd
  type (t_time) ,intent(in) :: time
    integer :: yy, mm, dd
    call caldat (time% days, mm, dd, yy)
    iyyyymmdd = yy*10000+mm*100+dd
  end function iyyyymmdd
!------------------------------------------------------------------------------
  elemental function iyyyy (time)
  integer                   :: iyyyy
  type (t_time) ,intent(in) :: time
    integer :: mm, dd
    call caldat (time% days, mm, dd, iyyyy)
  end function iyyyy
!------------------------------------------------------------------------------
  elemental function imm (time)
  integer                   :: imm
  type (t_time) ,intent(in) :: time
    integer :: dd, yy
    call caldat (time% days, imm, dd, yy)
  end function imm
!------------------------------------------------------------------------------
  elemental function idd (time)
  integer                   :: idd
  type (t_time) ,intent(in) :: time
    integer :: mm, yy
    call caldat (time% days, mm, idd, yy)
  end function idd
!------------------------------------------------------------------------------
  elemental function iyymmdd (time)
  integer                   :: iyymmdd
  type (t_time) ,intent(in) :: time
    iyymmdd = mod (iyyyymmdd (time),1000000)
  end function iyymmdd
!------------------------------------------------------------------------------
  elemental function chhmmss (time)
  character(len=6)             :: chhmmss
  type (t_time)    ,intent(in) :: time
    call i_to_c (chhmmss, ihhmmss (time))
  end function chhmmss
!------------------------------------------------------------------------------
  elemental function cddhhmmss (time)
  character(len=8)             :: cddhhmmss
  type (t_time)    ,intent(in) :: time
    call i_to_c (cddhhmmss, iddhhmmss (time))
  end function cddhhmmss
!------------------------------------------------------------------------------
  elemental function cdddhhmmss (time)
    character(len=9)             :: cdddhhmmss
    type(t_time)     ,intent(in) :: time
    call i_to_c (cdddhhmmss, iddhhmmss (time))
  end function cdddhhmmss
!------------------------------------------------------------------------------
  elemental function chhmm (time)
  character(len=4)             :: chhmm
  type (t_time)    ,intent(in) :: time
    call i_to_c (chhmm, ihhmm (time))
  end function chhmm
!------------------------------------------------------------------------------
  elemental function chhhmm (time)
  character(len=5)             :: chhhmm
  type (t_time)    ,intent(in) :: time
    call i_to_c (chhhmm, ihhhmm (time))
  end function chhhmm
!------------------------------------------------------------------------------
  elemental function chhhmmss (time)
    character(len=7)            :: chhhmmss
    type(t_time)    ,intent(in) :: time
    call i_to_c (chhhmmss, ihhhmmss (time))
  end function chhhmmss
!------------------------------------------------------------------------------
  elemental function cyyyymm (time)
  character(len=6)             :: cyyyymm
  type (t_time)    ,intent(in) :: time
    if (time% days /= 0) then
      call i_to_c (cyyyymm, iyyyymmdd (time) / 100)
    else
      call i_to_c (cyyyymm, 0)
    endif
  end function cyyyymm
!------------------------------------------------------------------------------
  elemental function cyyyymmdd (time)
  character(len=8)             :: cyyyymmdd
  type (t_time)    ,intent(in) :: time
    if (time% days /= 0) then
      call i_to_c (cyyyymmdd, iyyyymmdd (time))
    else
      call i_to_c (cyyyymmdd, 0)
    endif
  end function cyyyymmdd
!------------------------------------------------------------------------------
  elemental function cyyyymmddhh (time)
  character(len=10)            :: cyyyymmddhh
  type (t_time)    ,intent(in) :: time
    cyyyymmddhh(1: 8) = cyyyymmdd (time)
    call i_to_c (cyyyymmddhh(9:10), time% secs/3600)
  end function cyyyymmddhh
!------------------------------------------------------------------------------
  elemental function cyyyymmddhhmm (time)
  character(len=12)            :: cyyyymmddhhmm
  type (t_time)    ,intent(in) :: time
    cyyyymmddhhmm(1: 8) = cyyyymmdd (time)
    call i_to_c (cyyyymmddhhmm( 9:10),     time% secs/3600)
    call i_to_c (cyyyymmddhhmm(11:12), mod(time% secs/  60,60))
  end function cyyyymmddhhmm
!------------------------------------------------------------------------------
  elemental function cyyyymmddhhmmss (time)
  character(len=14)            :: cyyyymmddhhmmss
  type (t_time)    ,intent(in) :: time
    cyyyymmddhhmmss(1: 8) = cyyyymmdd (time)
    call i_to_c (cyyyymmddhhmmss( 9:10),     time% secs/3600)
    call i_to_c (cyyyymmddhhmmss(11:12), mod(time% secs/  60,60))
    call i_to_c (cyyyymmddhhmmss(13:14), mod(time% secs     ,60))
  end function cyyyymmddhhmmss
!------------------------------------------------------------------------------
  elemental function chhzddmmmyyyy (time) result (c)
  character(len=12)            :: c                    ! GRADS initial time
  type (t_time)    ,intent(in) :: time
    character(len=10) :: c10
    integer           :: i
    c10 = cyyyymmddhh (time)
    i   = imm         (time)
    c (1:2)  = c10 (9:10)
    c (3:3)  = 'z'
    c (4:5)  = c10 (7:8)
    c (6:8)  = mmm (i)
    c (9:12) = c10 (1:4)
  end function chhzddmmmyyyy
!------------------------------------------------------------------------------
  function time_grads (grads) result (time)
  character(len=*) ,intent(in) :: grads
  type(t_time)                 :: time
  !--------------------------------------------------------------------
  ! convert from GRADS date/time string: hh[:mm]Zddmmmyyyy
  !
  ! where:
  ! hh   = hour (two digit integer)
  ! mm   = minute (two digit integer)
  ! dd   = day (one or two digit integer)
  ! mmm  = 3-character month
  ! yyyy = year (may be a two or four digit integer;
  !              2 digits implies a year between 1950 and 2049)
  !
  ! If not specified, hh defaults to 00, mm defaults to 00, and dd
  ! defaults to 1. The month and year must be specified. No intervening
  ! blanks are allowed in the GrADS absolute date/time format.
  !--------------------------------------------------------------------

    integer :: yyyy,mo,dd,hh,mi
    integer :: im, id, iz, ib, j

    id = index (grads,':')
    iz = index (grads,'z'); if (iz==0) iz = index (grads,'Z')
    ib = index (grads,' '); if (ib==0) ib = len   (grads) + 1

    im = 0
    do j = 1, 12
      im = index (grads,mmm(j))
      if (im/=0) exit
    end do

    hh = 0
    if (id==0) id = iz
    select case (id)
    case (0)
    case (2)
      read (grads(1:id-1),'(i1)') hh
    case (3)
      read (grads(1:id-1),'(i2)') hh
    case default
      call finish ('time_grads',"invalid position of ':' or 'z'")
    end select

    mi = 0
    select case (iz-id)
    case (0,1)
    case (2)
      read (grads(id+1:iz-1),'(i1)') mi
    case (3)
      read (grads(id+1:iz-1),'(i2)') mi
    end select

    dd = 1
    select case (im-iz)
    case (0,1)
    case (2)
      read (grads(iz+1:im-1),'(i1)') dd
    case (3)
      read (grads(iz+1:im-1),'(i2)') dd
    end select

    if (im==0) call finish ('time_grads',"invalid string for month")
    mo = j

    select case (ib-im)
    case (5)
      read (grads(im+3:ib-1),'(i2)') yyyy; yyyy=yyyy + 1900
    case (7)
      read (grads(im+3:ib-1),'(i4)') yyyy
    case default
      call finish ('time_grads',"invalid string for 'yy' or 'yyyy'")
    end select

    call init_time (time, yyyy, mo, dd, hh, mi)

  end function time_grads
!------------------------------------------------------------------------------
  elemental function cvvkk  (time) result (c)
  character(len=4)             :: c                    ! GRADS time increment
  type (t_time)    ,intent(in) :: time
    if (time% days >= 10000) then
      call i_to_c (c(1:2),time% days/10000)
      c(3:4) = 'yr'
    else if (time% days >= 100) then
      call i_to_c (c(1:2),time% days/100)
      c(3:4) = 'mo'
    else if (time% days >= 1) then
      call i_to_c (c(1:2),time% days)
      c(3:4) = 'dy'
    else if (time% secs >= 3600) then
      call i_to_c (c(1:2),time% secs/3600)
      c(3:4) = 'hr'
    else
      call i_to_c (c(1:2),time% secs/60)
      c(3:4) = 'mn'
    endif
  end function cvvkk
!------------------------------------------------------------------------------
  elemental function cyymmdd (time)
  character(len=6)             :: cyymmdd
  type (t_time)    ,intent(in) :: time
    call i_to_c (cyymmdd, iyymmdd (time))
  end function cyymmdd
!------------------------------------------------------------------------------
  elemental function iyyddd (time)
  integer                   :: iyyddd
  type (t_time) ,intent(in) :: time
    integer i, yy, mm, dd
    call caldat (time% days, mm, dd, yy)
    i = time% days - julday  (1, 1, yy) + 1
    iyyddd = yy * 1000 + i
  end function iyyddd
!------------------------------------------------------------------------------
  elemental function cyyddd (time)
  character(len=5)             :: cyyddd
  type (t_time)    ,intent(in) :: time
    call i_to_c (cyyddd, iyyddd (time))
  end function cyyddd
!------------------------------------------------------------------------------
  elemental function cdate (time)
  character(len=10)            :: cdate
  type (t_time)    ,intent(in) :: time
    cdate(5:5) = '-'
    cdate(8:8) = '-'
    call i_to_c (cdate(1: 4), iyyyy (time))
    call i_to_c (cdate(6: 7), imm   (time))
    call i_to_c (cdate(9:10), idd   (time))
  end function cdate
!------------------------------------------------------------------------------
  elemental function ctime (time)
  character(len=8)             :: ctime
  type (t_time)    ,intent(in) :: time
    integer hhmmss
    hhmmss = ihhmmss (time)
    ctime(3:3) = ':'
    ctime(6:6) = ':'
    call i_to_c (ctime(1:2),      hhmmss / 10000    )
    call i_to_c (ctime(4:5), mod (hhmmss / 100 ,100))
    call i_to_c (ctime(7:8), mod (hhmmss       ,100))
  end function ctime
!==============================================================================
  elemental function add_times (time1, time2) result (time)
  type (t_time)              :: time
  type (t_time) ,intent(in)  :: time1
  type (t_time) ,intent(in)  :: time2
    integer :: secd, seconds, days
    secd = 60*60*24
    seconds = time1% secs + time2% secs
    days    = time1% days + time2% days
    if (seconds >= 0) then
      days = days + seconds/secd
      seconds = mod (seconds,secd)
    else
      days = days + seconds/secd - 1
      seconds = secd - mod (-seconds,secd)
    endif
    time% days = days
    time% secs = seconds
  end function add_times
!------------------------------------------------------------------------------
  elemental function subtr_times (time1, time2) result (time)
  type (t_time)              :: time
  type (t_time) ,intent(in)  :: time1
  type (t_time) ,intent(in)  :: time2
  integer :: secd, seconds, days
    secd = 60*60*24
    seconds = time1% secs - time2% secs
    days    = time1% days - time2% days
    if (seconds <  0) then
      days = days + seconds/secd - 1
      seconds = secd - mod (-seconds,secd)
    endif
    if (seconds >= 0) then
      days = days + seconds/secd
      seconds = mod (seconds,secd)
    endif
    time% days = days
    time% secs = seconds
  end function subtr_times
!------------------------------------------------------------------------------
  elemental function minus_time (time1) result (time)
  type (t_time)              :: time
  type (t_time) ,intent(in)  :: time1
  integer :: secd, seconds, days
    secd = 60*60*24
    seconds = - time1% secs
    days    = - time1% days
    if (seconds >= 0) then
      days = days + seconds/secd
      seconds = mod (seconds,secd)
    else
      days = days + seconds/secd - 1
      seconds = secd - mod (-seconds,secd)
    endif
    time% days = days
    time% secs = seconds
  end function minus_time
!==============================================================================
  elemental function time_gt_time (time1, time2) result (y)
  logical                    :: y
  type (t_time) ,intent(in)  :: time1
  type (t_time) ,intent(in)  :: time2
  if (time1% days > time2% days) then
    y = .true.
  else if (time1% days < time2% days) then
    y = .false.
  else if (time1% secs > time2% secs) then
    y = .true.
  else
    y = .false.
  endif
  end function time_gt_time
!------------------------------------------------------------------------------
  elemental function time_le_time (time1, time2) result (y)
  logical                    :: y
  type (t_time) ,intent(in)  :: time1
  type (t_time) ,intent(in)  :: time2
    y = .not. time1 > time2
  end function time_le_time
!------------------------------------------------------------------------------
  elemental function time_lt_time (time1, time2) result (y)
  logical                    :: y
  type (t_time) ,intent(in)  :: time1
  type (t_time) ,intent(in)  :: time2
    y = time2 > time1
  end function time_lt_time
!------------------------------------------------------------------------------
  elemental function time_ge_time (time1, time2) result (y)
  logical                    :: y
  type (t_time) ,intent(in)  :: time1
  type (t_time) ,intent(in)  :: time2
    y = .not. time2 > time1
  end function time_ge_time
!------------------------------------------------------------------------------
  elemental function time_eq_time (time1, time2) result (y)
  logical                    :: y
  type (t_time) ,intent(in)  :: time1
  type (t_time) ,intent(in)  :: time2
    y = time1% days == time2% days .and. time1% secs == time2% secs
  end function time_eq_time
!------------------------------------------------------------------------------
  elemental function time_ne_time (time1, time2) result (y)
  logical                    :: y
  type (t_time) ,intent(in)  :: time1
  type (t_time) ,intent(in)  :: time2
    y = time1% days /= time2% days .or. time1% secs /= time2% secs
  end function time_ne_time
!------------------------------------------------------------------------------
  elemental function sign_int_time (int, time) result (sign)
  integer                    :: sign
  type (t_time) ,intent(in)  :: time
  integer       ,intent(in)  :: int
    sign = abs (int)
    if (time <  zero_time) sign = - sign
    if (time == zero_time) sign =   0
  end function sign_int_time
!------------------------------------------------------------------------------
  elemental function max_time (time1, time2) result (y)
  type (t_time)              :: y
  type (t_time) ,intent(in)  :: time1
  type (t_time) ,intent(in)  :: time2
    if (time1 > time2) then
      y = time1
    else
      y = time2
    endif
  end function max_time
!------------------------------------------------------------------------------
  elemental function min_time (time1, time2) result (y)
  type (t_time)              :: y
  type (t_time) ,intent(in)  :: time1
  type (t_time) ,intent(in)  :: time2
    if (time1 > time2) then
      y = time2
    else
      y = time1
    endif
  end function min_time
!------------------------------------------------------------------------------
  pure function minval_time (x) result (y)
  type (t_time) ,intent(in)  :: x (:)
  type (t_time)              :: y
    integer       :: i
    if (size(x)==0) return
    y = x(1)
    do i=2,size(x)
      if (x(i) < y) y = x(i)
    end do
  end function minval_time
!------------------------------------------------------------------------------
  pure function minval_time_masked (x, mask) result (y)
  type (t_time)     ,intent(in) :: x   (:)
  logical, optional ,intent(in) :: mask(:)
  type (t_time)                 :: y
    if (present (mask)) then
       y = minval_time (pack (x, mask))
    else
       y = minval_time (x)
    end if
  end function minval_time_masked
!------------------------------------------------------------------------------
  pure function maxval_time (x) result (y)
  type (t_time) ,intent(in)  :: x (:)
  type (t_time)              :: y
    integer       :: i
    if (size(x)==0) return
    y = x(1)
    do i=2,size(x)
      if (x(i) > y) y = x(i)
    end do
  end function maxval_time
!------------------------------------------------------------------------------
  subroutine bcast_time (time, source, comm)
  type (t_time) ,intent(inout)        :: time
  integer       ,intent(in)           :: source
  integer       ,intent(in) ,optional :: comm
    call p_bcast (time% days, source, comm)
    call p_bcast (time% secs, source, comm)
  end subroutine bcast_time
!==============================================================================
  elemental function real_times_time (r, x) result (y)
  type (t_time)              :: y
  real (wp)     ,intent(in)  :: r
  type (t_time) ,intent(in)  :: x
    integer, parameter :: secd = 60*60*24
    real(wp)           :: days
    days    = r * x% days
    y% days = int (days)
    y% secs = nint (r * x% secs + secd * (days - y% days))
    if (y% secs >= 0) then
      y% days = y% days + y% secs/secd
      y% secs = mod (y% secs,secd)
    else
      y% days = y% days + y% secs/secd - 1
      y% secs = secd - mod (-y% secs,secd)
    endif
  end function real_times_time
!------------------------------------------------------------------------------
  elemental function integer_times_time (i, x) result (y)
  type (t_time)              :: y
  integer       ,intent(in)  :: i
  type (t_time) ,intent(in)  :: x
    integer, parameter :: secd = 60*60*24
    y% days = i * x% days
    y% secs = i * x% secs
    if (y% secs >= 0) then
      y% days = y% days + y% secs/secd
      y% secs = mod (y% secs,secd)
    else
      y% days = y% days + y% secs/secd - 1
      y% secs = secd - mod (-y% secs,secd)
    endif
  end function integer_times_time
!------------------------------------------------------------------------------
  elemental function time_over_real (x, r) result (y)
  type (t_time)              :: y
  type (t_time) ,intent(in)  :: x
  real (wp)     ,intent(in)  :: r
    integer, parameter :: secd = 60*60*24
    real(wp)           :: days
    days    = x% days / r
    y% days = int (days)
    y% secs = nint (x% secs / r + secd * (days - y% days))
    if (y% secs >= 0) then
      y% days = y% days + y% secs/secd
      y% secs = mod (y% secs,secd)
    else
      y% days = y% days + y% secs/secd - 1
      y% secs = secd - mod (-y% secs,secd)
    endif
  end function time_over_real
!------------------------------------------------------------------------------
  subroutine print_time (time, iunit, intend)
  type (t_time)    ,intent(in)           :: time   ! variable to print
  integer          ,intent(in) ,optional :: iunit  ! unit   (default=6)
  character(len=*) ,intent(in) ,optional :: intend ! intend string ('')
  !-------------------------
  ! print t_time (formatted)
  ! (for debugging)
  !-------------------------
    integer           :: iu, n
    character(len=32) :: c
    !--------------------
    ! optional parameters
    !--------------------
    n=0
    c=''
    if (present(intend)) then
      c = intend
      n = len_trim(c)
      if(n>0) then
        if(c(n:n)=='.') c(n:n)=' '
      endif
    endif
    iu = 6; if (present(iunit)) iu = iunit
    !------
    ! print
    !------
    write(iu,"(a,' (t_time)')") c(:n)
    write(iu,"(a,' days          :',i8)")  c(:n), time% days
    write(iu,"(a,' secs          :',i14)") c(:n), time% secs
    if (time% days > 0 .and. time% secs >= 0 .and. time% secs < 86400) then
      write(iu,"(a,' yyyymmddhhmmss:',a8,a6)") c(:n), cyyyymmdd(time),&
                                                       chhmmss (time)
    endif
  end subroutine print_time
!==============================================================================
  elemental function days_year (time) result (days)
  integer                   :: days
  type (t_time) ,intent(in) :: time
    integer :: yy, mm, dd
    call caldat (time% days, mm, dd, yy)
    days = julday (1,1,yy+1) - julday (1,1,yy)
  end function days_year
!==============================================================================
  elemental function frac_year (time) result (frac)
  real(kind=wp)             :: frac
  type (t_time) ,intent(in) :: time
    integer :: iy, jd0
    iy = iyyyy(time)
    jd0 = julday(1, 1, iy)
    frac = (time%days - jd0) + time%secs / 86400._wp
    frac = frac / days_year(time)
  end function frac_year
!------------------------------------------------------------------------------
! elemental function days_month (time, mm) result (days)
! integer                             :: days
! type (t_time) ,intent(in)           :: time
! integer       ,intent(in) ,optional :: mm
!   integer :: yy, mo, dd
!   call caldat (time% days, mo, dd, yy)
!   if(present(mm)) mo = mm
!   days = julday (1,mo,yy+1) - julday (1,mo,yy)    ! This is wrong!?
! end function days_month
!==============================================================================
  elemental function julday (mm, id, iyyy)
  integer             :: julday
  integer, intent(in) :: mm, id, iyyy
  !-----------------------------------------------------------------------
  ! In this routine JULDAY returns the Julian day number that begins at
  ! noon of the calender date specified by month MM, day ID and year IYYY,
  ! all integer variables. Positive year signifies A.D.; negative,
  ! B.C. Remember that the year after 1 B.C. was 1 A.D.
  !-----------------------------------------------------------------------
    integer, parameter  :: igreg = 15 + 31 * (10 + 12 * 1582)
    integer             :: ja, jm, jy
    jy=iyyy
!   if (jy == 0) call finish ('julday','there is no year zero')
    if (jy <= 0) jy = jy+1
    if (mm >  2) then
      jm = mm + 1
    else
      jy = jy -  1
      jm = mm + 13
    endif
    julday = int( 365.25 * jy) + int( 30.6001 *jm) + id + 1720995
    if (id + 31 * (mm + 12 * iyyy) >= IGREG) then
      ja     = int (0.01*jy)
      julday = julday + 2 - ja + int (0.25 * ja)
    endif
!print *,'julday: julday, mm, id, iyyy =',julday, mm, id, iyyy
  end function julday
! adapted from (C) Copr. 1986-92 Numerical Recipes Software
!------------------------------------------------------------------------------
  elemental subroutine caldat (julian, mm, id, iyyy)
  integer, intent(in)  :: julian
  integer, intent(out) :: mm, id, iyyy
  !-----------------------------------------------------------------------
  ! Inverse of the function JULDAY. Here JULIAN is input as a Julian day
  ! number, and the routine outputs MM, ID, and IYYY as the month, day and
  ! year on which the specified Julian day started at noon.
  !-----------------------------------------------------------------------
    integer, parameter   :: igreg = 2299161
    integer              :: ja, jalpha, jb, jc, jd, je
    if(julian >= igreg) then
      jalpha = int((( julian - 1867216) - 0.25) / 36524.25)
      ja     = julian + 1 + jalpha - int (0.25 * jalpha)
    else
      ja = julian
    endif
    jb = ja + 1524
    jc = int (6680. + ((jb - 2439870) - 122.1) / 365.25)
    jd = 365 * jc + int(0.25 * jc)
    je = int ((jb - jd) / 30.6001)
    id = jb - jd - int (30.6001 * je)
    mm = je - 1
    if (mm > 12) mm = mm - 12
    iyyy = jc - 4715
    if (mm   >  2) iyyy = iyyy - 1
    if (iyyy <= 0) iyyy = iyyy - 1
!print *,'caldat: julian, mm, id, iyyy=',julian, mm, id, iyyy
  end subroutine caldat
! adapted from (C) Copr. 1986-92 Numerical Recipes Software
!------------------------------------------------------------------------------
  elemental subroutine i_to_c (c,i)
  character(len=*) ,intent(out) :: c
  integer          ,intent(in)  :: i
  !------------------------------------------------------
  ! convert an non-negative integer to a character string
  ! right justified, padded with leading zeros
  !------------------------------------------------------
    integer, parameter :: i0 = ichar('0')
    integer :: j,n
    n  = i
    if (n >= 0) then
      do j=len(c),1,-1
        c(j:j)=char(mod(n,10)+i0)
        n=n/10
      end do
    end if
    if (n/=0) c = repeat ('*',len(c))
  end subroutine i_to_c
!==============================================================================
  elemental subroutine solar (time, phi, lam, azimuth, zenith)
  !---------------------------------
  ! calculate solar parameters:
  !   azimuth  azimuth of sun
  !   h        solar zenith angle
  !---------------------------------
  type(t_time) ,intent(in)  :: time    ! date & time
  real(wp)     ,intent(in)  :: phi     ! latitude            (degree)
  real(wp)     ,intent(in)  :: lam     ! longitude           (degree)
  real(wp)     ,intent(out) :: azimuth ! solar azimuth       (degree)
  real(wp)     ,intent(out) :: zenith  ! solar zenith  angle (degree)

    real(wp)        :: jd
    real(wp)        :: L, g, gamma, alpha,dec, t0, n
    real(wp)        :: tstar, tfrueh, taus, ecobli, h

    real(wp) ,parameter :: degrad = 57.295779513_wp
    real(wp) ,parameter :: raddeg = 0.0174532925_wp

    !------------
    ! Julian date
    !------------
    jd = time% days
    n = jd - 2451545._wp

    !-------------------------------
    ! mittlere Ekliptikale der Sonne
    !-------------------------------
    L = modulo (280.460_wp + 0.9856474_wp * n, 360.0_wp)

    !-------------------
    ! mittlere Annomalie
    !-------------------
    g = modulo (357.528_wp + 0.9856003_wp * n, 360.0_wp)

    !-----------------------------------------
    ! Mittelpunktsgleichung ekliptische Laenge
    !-----------------------------------------
    gamma = modulo (L + 1.915_wp  * sin(g*raddeg)       &
                      + 0.020e0_wp* sin(2._wp*g*raddeg) &
                   , 360.0_wp)
    ecobli = 23.439_wp - 4.0e-07_wp * n

    !--------------
    ! Rektaszension
    !--------------
    alpha = degrad * atan ((cos (ecobli*raddeg) * sin (gamma*raddeg)) &
          / cos (gamma*raddeg))
    if (cos(gamma*raddeg) < 0.0_wp) alpha = alpha + 180.0_wp

    !------------
    ! Deklination
    !------------
    dec = degrad * asin(sin(ecobli*raddeg) * sin(gamma*raddeg))

    !-------------------------
    ! Julianisches Jahrhundert
    !-------------------------
    t0    = (jd - 2451545._wp) / 36525.0_wp

    !----------
    ! Sternzeit
    !----------
    tstar = modulo (                                                         &
            6.697376_wp + 2400.05134_wp*t0 + 1.002738_wp*time% secs/3600._wp &
            , 24.0_wp)

    !------------------------------------
    ! Stundenwinkel des Fruehlingspunktes
    !------------------------------------
    tfrueh = modulo (tstar * 15.0_wp + lam, 360.0_wp)

    !------------------------
    ! Stundenwinkel der Sonne
    !------------------------
    taus = tfrueh - alpha
    if (taus < -180.0_wp) taus = taus + 360.0_wp
    if (taus >  180.0_wp) taus = taus - 360.0_wp

    !--------
    ! Azimuth
    !--------
    azimuth = degrad * atan(sin(taus*raddeg)       &
              /(cos(taus*raddeg)*sin(phi*raddeg)   &
                - tan(dec*raddeg)* cos(phi*raddeg)))
    if ((cos(taus*raddeg)*sin(phi*raddeg) -       &
         tan(dec *raddeg)*cos(phi*raddeg)) < 0._wp) azimuth = azimuth + 180.0_wp
    if (azimuth < -180.0_wp)                        azimuth = azimuth + 360.0_wp
    if (azimuth >  180.0_wp)                        azimuth = azimuth - 360.0_wp

    !-------------
    ! Hoehenwinkel
    !-------------
    h = degrad * asin(cos(dec*raddeg)*cos(taus*raddeg)*cos(phi*raddeg) &
               + sin(dec*raddeg)*sin(phi*raddeg))
    !-------------------------------------------------------
    ! return solar zenith angle, same convention as in RTTOV
    !-------------------------------------------------------
    zenith = 90._wp - h

  end subroutine solar
!==============================================================================
  elemental function localdate (time, lon)
    type(t_time), intent(in)  :: time           ! Reference time (UT1)
    real(wp),     intent(in)  :: lon            ! Longitude [deg]
    type(t_time)              :: localdate      ! Local (mean solar) date/time
    !--------------------------------------------------
    ! Derive local (mean solar) date/time for longitude
    !--------------------------------------------------
    real(wp) :: dlon      ! Longitude reduced to range [-180°,180°]
    real(wp) :: dt        ! Time shift [s]

    dlon = modulo (lon, 360._wp)
    if (dlon > 180._wp) dlon = dlon - 360._wp
    dt   = dlon * (86400._wp / 360._wp)
    localdate = time + t_time (days=0, secs=nint (dt))
  end function localdate
!==============================================================================
  subroutine date_tmpl(str, idfm, ind0, ind1, t)
    character(len=*), intent(inout)          :: str
    integer,          intent(out),  optional :: idfm
    integer,          intent(out),  optional :: ind0
    integer,          intent(out),  optional :: ind1
    type(t_time),     intent(in),   optional :: t

    character(len=14) :: dstr
    character(len=len(str)) :: aux
    integer           :: j, ii, i0, i1

    ii = 0
    do j = 1, n_tmpl_forms
      i0 = index(str, trim(date_tmpl_forms(j)))
      if (i0 > 0) then
        ii = j
        i1 = i0-1 + len_trim(date_tmpl_forms(j))
        exit
      end if
    end do
    if (present(idfm)) idfm = ii
    if (present(ind0)) ind0 = i0
    if (present(ind1)) ind1 = i1

    if (present(t) .and. ii > 0) then
      select case(ii)
      case(1)
        dstr = cyyyymmddhh(t)
        dstr = dstr(3:)
      case(2)
        dstr = cyyyymmddhh(t)
      case(3)
        dstr = cyyyymmddhhmm(t)
      case(4)
        dstr = cyyyymmddhhmmss(t)
      end select
      if (len_trim(str) > i1) then
        aux = str(i1+1:)
      else
        aux = ''
      end if
      str = str(1:i0-1)//trim(dstr)//trim(aux)
    end if

  end subroutine date_tmpl

end module mo_time
