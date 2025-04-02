!
!+ read and write GRADS data sets
!
MODULE mo_grads
!
! Description:
!   Read and write GRADS data sets:
!   IEEE binary files (little/big endian transparent)
!   and corresponding .ctl GRADS meta data files
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
! V1_8         2009/12/09 Harald Anlauf
!  write_var3, read_var2, read_var3: print name of variable before crashing
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Andreas Rhodin
!  changed 3dvar revision numbers (CVS->SVN)
! V1_26        2013/06/27 Andreas Rhodin
!  new routines set_path set_time, record
! V1_31        2014-08-21 Andreas Rhodin
!  enhancements for station data files
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Andreas Rhodin  MPIfM/DWD  2000-2007  original source
! Harald Anlauf   DWD        2007-2008  fixes for various compiler bugs
!==============================================================================
#include "tr15581.incf"
  !-------------
  ! used modules
  !-------------
  use mo_kind,          only: sp, wp, i4         ! single,working precision
  use mo_endian,        only: little,           &! returns .true. on LEmachine
                              flip               ! flip single precision field
  use mo_fortran_units, only: get_unit_number, & ! obtain Fortran unit number
                              return_unit_number ! release Fortran unit number
  use mo_exception,     only: finish             ! abort on errror condition
  use mo_transform,     only: gauaw              ! calculate Gaussian latitudes
  use mo_dace_string,   only: tolower,          &! Lowercase of string
                              split              ! Split text string into array
  use mo_time,          only: t_time,           &! time derived type
                              chhzddmmmyyyy,    &! GRADS initial time
                              cvvkk              ! GRADS time increment
  implicit none
!==============================================================================
  !----------------
  ! public entities
  !----------------
  private
  public :: t_ctl, &   ! data type holding .ctl file information
            t_var, &   ! component of t_ctl holding information on datasets
            t_stat     ! component of t_ctl (station data header)
  public :: read_ctl   ! read .ctl file information
  public :: write_ctl  ! write .ctl file
  public :: init_ctl   ! preset .ctl file
  public :: set_path   ! set path and dataset entries in .ctl file
  public :: set_time   ! set time and increment
  public :: add_var    ! add variable in .ctl file
  public :: write_stat ! write station data to grads file
  public :: write_var  ! write variable to grads file
  public :: record     ! determine record number in GRADS file
  public :: read_var   ! read variable from grads file
  public :: print_stat ! print content of a station data file
  public :: c_month    ! character string encoding for months
  public :: mmm2mm     ! month character string encoding to integer
  public :: destruct   ! destruct t_ctl data type
!==============================================================================
  !-----------
  ! data types
  !-----------
  type t_var
    character(len=32)         :: name          = ''
    integer                   :: levels        = 1
    character(len=16)         :: units         = '99'
    character(len=40)         :: comment       = 'no comment'
    integer                   :: records       = 0
  end type t_var

  type t_stat
    character(len=8) :: id   = ''
    real(sp)         :: lat  = -huge(1._sp)
    real(sp)         :: lon  = -huge(1._sp)
    real(sp)         :: t    = 0._sp
    integer(i4)      :: nlev = 0
    integer(i4)      :: flag = 0
  end type t_stat

  type t_ctl
    character(len=128)        :: path          = '' ! actual path+filename
    character(len=128)        :: dset          = '' ! filename written to .ctl
    character(len=128)        :: stnmap        = ''
    character(len=128)        :: title         = ''
    character(len=8)          :: dtype         = ''
    type(t_stat)              :: stat
    integer                   :: irecstat      = 0
    integer                   :: iunit         = -1 ! Fortran unit number
    logical                   :: yrev          = .false.
    logical                   :: zrev          = .false.
    logical                   :: sequential    = .false.
    logical                   :: big_endian    = .false.
    logical                   :: little_endian = .false.
    logical                   :: byteswapped   = .false.
    integer                   :: fileheader    = 0
    integer                   :: theader       = 0
    integer                   :: xyheader      = 0
    real(sp)                  :: undef         = -huge(1._sp)
    integer                   :: xdefn         = 1
    integer                   :: ydefn         = 1
    integer                   :: zdefn         = 1
    integer                   :: tdefn         = 1
    character(len=8)          :: xdef          = 'linear'
    character(len=8)          :: ydef          = 'linear'
    character(len=8)          :: tdef          = 'linear'
    character(len=8)          :: zdef          = 'linear'
    real(sp)                  :: xdefi(2)      = (/0.,1./)
    real(sp)                  :: ydefi(2)      = (/0.,1./)
    real(sp)                  :: zdefi(2)      = (/0.,1./)
#if defined(__GFORTRAN__)
    ! Work around gfortran bug #31487
    character(len=16)         :: tdefi(2)      = (/'0z1jan0000      ',&
                                                   '1hr             '/)
#else
    character(len=16)         :: tdefi(2)      = (/'0z1jan0000','1hr       '/)
#endif
    integer                   :: vars          = 0
    integer                   :: svars         = 0
    integer                   :: lvars         = 0
    integer                   :: records       = 0
    type (t_var)     ,pointer :: var(:)        => NULL()
    real(wp)         ,pointer :: xlevels(:)    => NULL()
    real(wp)         ,pointer :: ylevels(:)    => NULL()
    real(wp)         ,pointer :: zlevels(:)    => NULL()
    character(len=78),pointer :: comment(:)    => NULL()
  end type t_ctl
  !----------
  ! constants
  !----------
  character(len=3) ,parameter :: c_month (12) = &
                                 (/'jan','feb','mar','apr','may','jun',&
                                   'jul','aug','sep','oct','nov','dec'/)
  character(len=3) ,parameter :: u_month (12) = &
                                 (/'JAN','FEB','MAR','APR','MAY','JUN',&
                                   'JUL','AUG','SEP','OCT','NOV','DEC'/)
  !-----------
  ! interfaces
  !-----------
  interface write_var
    module procedure write_var2
    module procedure write_var3
  end interface write_var

  interface read_var
    module procedure read_var2
    module procedure read_var3
  end interface read_var

  interface init_ctl
    module procedure init_ctl_
  end interface init_ctl

  interface destruct
    module procedure destruct_ctl
  end interface destruct
!------------------------------------------------------------------------------
contains
!------------------------------------------------------------------------------
  function mmm2mm (c) result (i)
  character (len=3) ,intent(in) :: c
  integer                       :: i
  !-------------------------------------------
  ! month character string encoding to integer
  !-------------------------------------------
    integer :: k
    i = 0
    do k=1,12
      if(c==c_month(k) .or. c==u_month(k)) then
        i = k
        exit
      endif
    end do
  end function mmm2mm
!------------------------------------------------------------------------------
  subroutine set_path (ctl, file)
  type (t_ctl)         ,intent(inout) :: ctl
  character(len=*)     ,intent(in)    :: file
  !-------------------------------------------------
  ! set dset      (file name in .ctl, may contain ^)
  ! pathname used (without ^)
  !-------------------------------------------------
    integer :: idx

    idx = index (file, '^')
    if (idx > 0) then
      ctl% dset = file(idx:)
      ctl% path = file(:idx-1)//file(idx+1:)
    else
      ctl% dset = file
      ctl% path = file
    endif
  end subroutine set_path
!------------------------------------------------------------------------------
  subroutine set_time (ctl, tdefn, tdefi, tdefd, t_defi, t_defd)
  type (t_ctl)         ,intent(inout)        :: ctl
  integer              ,intent(in) ,optional :: tdefn     !no.time slices
  character(len=*)     ,intent(in) ,optional :: tdefi     !initial time
  character(len=*)     ,intent(in) ,optional :: tdefd     !time increment
  type(t_time)         ,intent(in) ,optional :: t_defi    !initial time
  type(t_time)         ,intent(in) ,optional :: t_defd    !time increment
  !---------------------------
  ! set time info in .ctl file
  !---------------------------
    if (present(tdefn )) ctl% tdefn    = tdefn
    if (present(tdefi )) ctl% tdefi(1) = tdefi
    if (present(tdefd )) ctl% tdefi(2) = tdefd
    if (present(t_defi)) ctl% tdefi(1) = chhzddmmmyyyy (t_defi)
    if (present(t_defd)) ctl% tdefi(2) = cvvkk         (t_defd)
  end subroutine set_time
!------------------------------------------------------------------------------
  subroutine init_ctl_(ctl, file, dtype, title,                       &
                       nx, ny, ngl, ke, di, dj, lo1, la1, xlev, zlev, &
                       tdefn, tdefi, tdefd, undef, yrev, zrev,        &
                       surdat, levdat, surcom, levcom, t_defi, t_defd,&
                       tmpl, comment)
  type (t_ctl)         ,intent(inout)        :: ctl
  character(len=*)     ,intent(in) ,optional :: file
  character(len=*)     ,intent(in) ,optional :: dtype
  character(len=*)     ,intent(in) ,optional :: title
  integer              ,intent(in) ,optional :: nx, ny, ngl, ke
  real(wp)             ,intent(in) ,optional :: di, dj, lo1, la1
  real(wp)             ,intent(in) ,optional :: xlev (:)  !Disguised longitudes
  real(wp)             ,intent(in) ,optional :: zlev (:)  !Pressure levels
  integer              ,intent(in) ,optional :: tdefn     !no.time slices
  character(len=*)     ,intent(in) ,optional :: tdefi     !initial time
  character(len=*)     ,intent(in) ,optional :: tdefd     !time increment
  real(wp)             ,intent(in) ,optional :: undef     !undefined value
  logical              ,intent(in) ,optional :: yrev
  logical              ,intent(in) ,optional :: zrev
  character(len=*)     ,intent(in) ,optional :: surdat(:) !station surface data
  character(len=*)     ,intent(in) ,optional :: levdat(:) !station level   data
  character(len=*)     ,intent(in) ,optional :: surcom(:) !station s. comments
  character(len=*)     ,intent(in) ,optional :: levcom(:) !station l. comments
  type(t_time)         ,intent(in) ,optional :: t_defi    !initial time
  type(t_time)         ,intent(in) ,optional :: t_defd    !time increment
  type(t_ctl)          ,intent(in) ,optional :: tmpl      !template
  character(len=*)     ,intent(in) ,optional :: comment(:)!comment lines
    integer :: i, ic
    call destruct (ctl)
    ctl% big_endian = .true.
    ctl% records    = 0
    !--------------------------------
    ! set GRADS binary data file name
    !--------------------------------
    if (present(file   )) call set_path (ctl, file)
    if (present(comment)) then
      allocate (ctl% comment (size(comment)))
      ctl% comment = comment
    endif
    !------------------------------
    ! take parameters from template
    !------------------------------
    if (present(tmpl )) then
      ctl% dtype    = tmpl% dtype
      ctl% title    = tmpl% title
      ctl% undef    = tmpl% undef
      ctl% yrev     = tmpl% yrev
      ctl% zrev     = tmpl% zrev

      ctl% tdefn    = tmpl% tdefn
      ctl% tdefi    = tmpl% tdefi
      ctl% xdefn    = tmpl% xdefn
      ctl% xdefi    = tmpl% xdefi
      ctl% ydefn    = tmpl% ydefn
      ctl% ydefi    = tmpl% ydefi
      ctl% zdefn    = tmpl% zdefn
      ctl% zdefi    = tmpl% zdefi

      if (associated (tmpl% xlevels)) then
        allocate (ctl% xlevels (size (tmpl% xlevels)))
        ctl% xlevels = tmpl% xlevels
        ctl% xdef    = 'levels'
      endif

      if (associated(tmpl% ylevels)) then
        allocate (ctl% ylevels (size(tmpl% ylevels)))
        ctl% ylevels = tmpl% ylevels
        ctl% ydef    = 'levels'
      endif

      if (associated(tmpl% zlevels)) then
        allocate (ctl% zlevels (size(tmpl% zlevels)))
        ctl% zlevels = tmpl% zlevels
        ctl% zdef    = 'levels'
      endif
    endif
    !---------------
    ! set some flags
    !---------------
    if (present(dtype)) ctl% dtype    = dtype
    if (present(title)) ctl% title    = title
    if (present(undef)) ctl% undef    = undef
    if (present(yrev )) ctl% yrev     = yrev
    if (present(zrev )) ctl% zrev     = zrev
    !------------------
    ! set gaussian grid
    !------------------
    if (present(ngl  )) then
      if(ngl>0) then
        ctl% ydef = 'levels'
        if (associated(ctl% ylevels)) deallocate (ctl% ylevels)
        allocate (ctl% ylevels(ngl))
        ctl% ydefn    = ngl
        ctl% xdefn    = ngl * 2
        ctl% xdefi(1) = 0._wp
        ctl% xdefi(2) = 360._wp/ctl% xdefn
        call gauaw (ga=ctl% ylevels)
        ctl% ylevels  = - 90._wp/asin(1._wp) * asin (ctl% ylevels)
        if (ctl% yrev) ctl% ylevels  = - ctl% ylevels
        ctl% ydefi(1) = ctl% ylevels(1)
      endif
    endif
    !-------------------------------------------
    ! set regular grid, explicit dimensions, etc
    !-------------------------------------------
    call set_time (ctl, tdefn, tdefi, tdefd, t_defi, t_defd)
    if (present(nx    )) ctl% xdefn    = nx
    if (present(lo1   )) ctl% xdefi(1) = lo1
    if (present(di    )) ctl% xdefi(2) = di
    if (present(ny    )) ctl% ydefn    = ny
    if (present(la1   )) ctl% ydefi(1) = la1
    if (present(dj    )) ctl% ydefi(2) = dj
    if (present(ke    )) ctl% zdefn    = ke
    !----------------------------------------------------
    ! Handle nonequidistant x-axis (disguised longitudes)
    !----------------------------------------------------
    if (present (xlev)) then
       ctl% xdefn   = size (xlev); if (present (nx)) ctl% xdefn = nx
       ctl% xdef    = 'levels'
       if (associated (ctl% xlevels)) deallocate (ctl% xlevels)
       allocate (ctl% xlevels (ctl% xdefn))
       ctl% xlevels = xlev(:ctl% xdefn)
    endif
    !-----------------------------------
    ! set nonequidistant vertical levels
    !-----------------------------------
    if (present(zlev)) then
      ctl% zdefn   = size (zlev); if (present(ke)) ctl% zdefn = ke
      ctl% zdef    = 'levels'
      if (associated(ctl% zlevels)) deallocate (ctl% zlevels)
      allocate (ctl% zlevels (ctl% zdefn))
      ctl% zlevels = zlev(:ctl% zdefn)
!     if (zlev(1) < zlev(ctl% zdefn)) then   ! enforce bottom to top
!       ctl% zlevels = zlev(ctl% zdefn:1:-1)
!     else
!       ctl% zlevels = zlev(:ctl% zdefn)
!     endif
    endif
    !------------------------------------
    ! define surface variables if present
    !------------------------------------
    if (present(surdat)) then
      ctl% dtype = 'station'
      ic = 0; if (present(surcom)) ic = size (surcom)
      do i=1,size(surdat)
        if(i<=ic) then
          call add_var (ctl, surdat(i), levels=0, comment=surcom(i))
        else
          call add_var (ctl, surdat(i), levels=0)
        endif
      end do
    endif
    !--------------------------------
    ! define surface level if present
    !--------------------------------
    if (present(levdat)) then
      ctl% dtype = 'station'
      ic = 0; if (present(levcom)) ic = size (levcom)
      do i=1,size(levdat)
        if(i<=ic) then
          call add_var (ctl, levdat(i), comment=levcom(i))
        else
          call add_var (ctl, levdat(i))
        endif
      end do
    endif
    !----------------------------------
    ! special settings for station data
    !----------------------------------
    select case (ctl% dtype)
    case ('station','STATION')
      i = index(ctl% dset, '.dat ')
      if (i/=0) then
        ctl% stnmap = ctl% dset(1:i)//'map'
      else
        ctl% stnmap = trim(ctl% dset)//'.map'
      endif
      ctl% little_endian =       little()
      ctl% big_endian    = .not. little()
      ctl% zdefn         = 1
      ctl% records       = 0
    end select
    !--------------------------------------------------------
    ! Patch:
    ! write grads files with the native endian
    ! (conversion currently doesnt work with the NAG compiler
    !--------------------------------------------------------
    ctl% little_endian =       little()
    ctl% big_endian    = .not. little()
  end subroutine init_ctl_
!------------------------------------------------------------------------------
  subroutine destruct_ctl (ctl)
  type (t_ctl) ,intent(inout) :: ctl
    integer :: ierr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! something goes wrong with the NAG compiler
!! obviously there is no default-initialization for t_ctl declared
!! in other modules. Therefore the stat parameter is used
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    type (t_ctl) :: empty_ctl
    ierr = 0
    if (associated(ctl% var)) then
      deallocate (ctl% var     ,stat=ierr)
      if (ierr>0) write (0,*)'destruct_ctl: error: deallocate (var)'
    endif
    if (associated(ctl% xlevels)) then
       deallocate (ctl% xlevels ,stat=ierr)
       if (ierr>0) write (0,*)'destruct_ctl: error: deallocate (xlevels)'
    endif
    if (associated(ctl% ylevels)) then
      deallocate (ctl% ylevels ,stat=ierr)
      if (ierr>0) write (0,*)'destruct_ctl: error: deallocate (ylevels)'
    endif
    if (associated(ctl% zlevels)) then
      deallocate (ctl% zlevels ,stat=ierr)
      if (ierr>0) write (0,*)'destruct_ctl: error: deallocate (zlevels)'
    endif
    if (associated(ctl% comment)) then
      deallocate (ctl% comment ,stat=ierr)
      if (ierr>0) write (0,*)'destruct_ctl: error: deallocate (comment)'
    endif
    ctl = empty_ctl
  end subroutine destruct_ctl
!==============================================================================
  subroutine add_var (ctl, name, levels, comment)
  !---------------------------
  ! add a variable to the list
  !---------------------------
  type (t_ctl)     ,intent(inout)        :: ctl
  character(len=*) ,intent(in)           :: name
  integer          ,intent(in) ,optional :: levels
  character(len=*) ,intent(in) ,optional :: comment
    type (t_var) ,allocatable  :: var(:)
    integer                    :: n, i
    !---------------------------------------------
    ! check if variable name is already registered
    !---------------------------------------------
    do i=1,ctl%vars
      if (ctl% var(i)% name == name) then
        if (present(levels)) then
          if (ctl% var(i)% levels /= levels) call finish('add_var',&
            'inconsistent number of levels: '//ctl% var(i)% name)
        endif
        return
      endif
    end do
    !--------------
    ! increase list
    !--------------
    n = ctl%vars + 1
    if (associated(ctl% var)) then
      allocate (var(n-1));    var = ctl% var;       deallocate (ctl% var)
      allocate (ctl% var(n)); ctl% var(:n-1) = var; deallocate (var)
    else
      allocate (ctl% var(n))
    endif
    ctl% vars = n
    !----------------------
    ! insert new list entry
    !----------------------
    ctl% var(n)% name    = name
    ctl% var(n)% levels  = ctl% zdefn
    if (present(levels))  ctl% var(n)% levels  = levels
    if (ctl% var(n)% levels == 0) then
      ctl% svars = ctl% svars + 1
    else
      ctl% lvars = ctl% lvars + 1
    endif
    if (present(comment)) ctl% var(n)% comment = comment
    ctl% var(n)% records = ctl% records + 1
    ctl% records = ctl% records + ctl% var(n)% levels
  end subroutine add_var
!------------------------------------------------------------------------------
  subroutine write_stat (ctl, id, lat, lon, z, surf, lev, new)
  !-----------------------------------
  ! write surface data to a grads file
  !-----------------------------------
  type(t_ctl)        ,intent(inout) :: ctl     ! .ctl file information
  character(len=*)   ,intent(in)    :: id      ! station id
  real(wp)           ,intent(in)    :: lat     ! latitude  [degree]
  real(wp)           ,intent(in)    :: lon     ! longitude [degree]
!+real(wp) ,optional ,intent(in)    :: t       ! time slice
  real(wp) ,optional ,intent(in)    :: z       ! level
  real(wp) ,optional ,intent(in)    :: surf(:) ! surface data
  real(wp) ,optional ,intent(in)    :: lev (:) ! level data
  logical  ,optional ,intent(in)    :: new     ! new station indicator
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! time handling was not properly implemented. Temporarily disabled.
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    type(t_stat)     :: empty_header
    integer          :: irecl
    logical          :: new_station
    character(len=7) :: status
    integer          :: iveri
    logical          :: exist
    !----------
    ! open file
    !----------
    if (ctl% iunit < 0) then
      inquire (iolength= irecl) 1._sp
      inquire (file=ctl% path, exist=exist)
      ctl% iunit = get_unit_number()
      status =                        'old'
      if (ctl% records == 0) status = 'replace'
      if (.not.exist)        status = 'new'
      open (ctl% iunit ,file=ctl% path ,              &
            access='direct' ,recl=irecl, status=status)
    endif
    !------------------------------------------------
    ! write empty header if new time slice is reached
    !------------------------------------------------
!+    if(present(t)) then
!+      if (ctl% stat% t /= t .and. ctl% records/=0) then
!+        ctl% stat    = empty_header
!+        ctl% stat% t = t
!+        call write_header (ctl% stat, ctl% records)
!+        ctl% records  = ctl% records + 7
!+        ctl% irecstat = ctl% records
!+      endif
!+    endif
    !-------------------------------------------------
    ! new id or new surface data indicates new station
    ! print warning for double nonumeric station id
    !-------------------------------------------------
    new_station = (id /= ctl% stat% id)
    if (.not. new_station .and. present(surf)) then
      iveri = verify (id, ' 0123456789')
      if (iveri==0) write (0,*) 'write_stat: WARNING: double surface data:',id
      new_station = .true.
    endif
    if (present(new)) new_station = new
    if (new_station) then
      !--------------------------------------
      ! new station: store header information
      !--------------------------------------
      ctl% stat% id   = id
      ctl% stat% lat  = lat
      ctl% stat% lon  = lon
      ctl% stat% nlev = 0
      ctl% stat% flag = 0
!+    if(present(t)) ctl% stat% t = t
      !------------------------------
      ! keep space for station header
      !------------------------------
      ctl% irecstat = ctl% records
      ctl% records  = ctl% records + 7
    endif
    !-------------------
    ! write surface data
    !-------------------
    if (present(surf)) then
!     if (ctl% stat% flag == 1) call finish('write_stat','double surface data')
!     if (ctl% stat% nlev >  0) call finish('write_stat','surf data after lev')
      ctl% stat% flag = 1
      ctl% stat% nlev = 1
      call write_data (surf)
    endif
    !-----------------
    ! write level data
    !-----------------
    if (present(lev)) then
      if(.not.present(z))call finish('write_stat','z not present for lev data')
      ctl% stat% nlev = ctl% stat% nlev + 1
      call write_data ((/z/))
      call write_data (lev)
    endif
    !-------------
    ! write header
    !-------------
    call write_header (ctl% stat, ctl% irecstat)
    !------------------------------------------
    ! write empty header in advance, close file
    !------------------------------------------
    call write_header (empty_header, ctl% records)
  contains
    !-------------
    ! write header
    !-------------
    subroutine write_header (header, irec)
    type (t_stat) :: header
    integer       :: irec
      write (ctl% iunit, rec = irec + 1) header% id(1:4)
      write (ctl% iunit, rec = irec + 2) header% id(5:8)
      write (ctl% iunit, rec = irec + 3) header% lat
      write (ctl% iunit, rec = irec + 4) header% lon
      write (ctl% iunit, rec = irec + 5) header% t
      write (ctl% iunit, rec = irec + 6) header% nlev
      write (ctl% iunit, rec = irec + 7) header% flag
    end subroutine write_header
    !-----------
    ! write data
    !-----------
    subroutine write_data (data)
    real(wp) :: data (:)
      integer  :: j
      do j=1,size(data)
        ctl% records = ctl% records + 1
        write (ctl% iunit, rec = ctl% records) real(data(j),sp)
      end do
    end subroutine write_data
  end subroutine write_stat
!------------------------------------------------------------------------------
  subroutine print_stat (ctlfile)
  character(len=*) ,intent(in) :: ctlfile
  !-------------------------------------
  ! print content of a station data file
  !-------------------------------------
    integer               :: irecl
    integer               :: iunit
    integer               :: ios
    type (t_ctl)          :: ctl
    real(wp)              :: lev  (1)
    real(wp) ,allocatable :: svar (:)
    real(wp) ,allocatable :: lvar (:)
    integer               :: flag
    integer               :: nlev
    integer               :: ilev
    !--------------
    ! read ctl-file
    !--------------
    call read_ctl (ctl, ctlfile, ios)
    if (ios /= 0) &
      call finish ('print_stat','cannot open ctl-file: '//trim(ctlfile))
    !---------------
    ! write ctl-file
    !---------------
    write(6,'()')
    call write_ctl(ctl, unit=6)
    write(6,'(a,i4)') 'svars' ,     ctl% svars
    write(6,'(a,i4)') 'lvars' ,     ctl% lvars
    write(6,'(a,a)')  'path  ',trim(ctl% path)
    write(6,'()')
    !-------
    ! checks
    !-------
    select case (ctl% dtype)
    case ('station','STATION')
    case default
      call finish ('print_stat','no station data file: '//trim(ctlfile))
    end select
    if (ctl% svars + ctl% lvars /= ctl% vars) then
      write(0,*) 'print_stat: ctl% svars + ctl% lvars /= ctl% vars:',&
                              ctl% svars,  ctl% lvars,   ctl% vars
      call finish ('print_stat','ctl% svars + ctl% lvars /= ctl% vars')
    endif
    allocate (svar (ctl% svars))
    allocate (lvar (ctl% lvars))
    !----------
    ! open file
    !----------
    inquire (iolength= irecl) 1._sp
    iunit = get_unit_number()
    write(6,*) 'opening file: >'//trim(ctl% path)//'<'
    open  (iunit ,file= ctl% path ,access= 'direct' ,recl= irecl, status='old')
    ctl% records = 0
    do
      call read_header  (ctl% stat, ctl% records, ios)
      if (ios/=0) then
        write(6,*) 'print_stat: iostat =',ios
        exit
      endif
      call print_header (ctl% stat)
      flag = ctl% stat% flag
      select case (flag)
      case (1,0)
      case default
        write(0,*) 'print_stat: station data flag =',ctl% stat% flag
        call finish ('print_stat','station data flag /= 0 or 1')
      end select
      nlev = ctl% stat% nlev
      do ilev = 1,nlev
        if (flag==1) then
          write(6,*) 'reading surface data',ilev
          call read_data (svar)
          write(6,'(a,(20f8.2))') 'surface data',svar
          flag = 0
        else
          write(6,*) 'reading level data',ilev
          call read_data (lev)
          write(6,'(a,f12.0)') 'level ',lev
          call read_data (lvar)
          write(6,'(a,(20f8.2))') 'level data',lvar
        endif
      end do
    end do
    close (iunit)
    call return_unit_number(iunit)
    deallocate (svar, lvar)
  contains

    subroutine read_header (header, irec, iostat)
    type (t_stat) ,intent(out)            :: header
    integer       ,intent(inout)          :: irec
    integer       ,intent(out)  ,optional :: iostat
      integer :: ios
      read (iunit, rec = irec + 1, iostat = ios) header% id(1:4)
      if (ios/=0) then
        write(6,*)'read_header: ERROR reading station data header'
        if (present (iostat)) then
          iostat = ios
          return
        else
           call finish ('read_header','ERROR reading station data header')
        endif
      endif
      read (iunit, rec = irec + 2) header% id(5:8)
      read (iunit, rec = irec + 3) header% lat
      read (iunit, rec = irec + 4) header% lon
      read (iunit, rec = irec + 5) header% t
      read (iunit, rec = irec + 6) header% nlev
      read (iunit, rec = irec + 7) header% flag
      irec = irec + 7
      if (present(iostat)) iostat = 0
    end subroutine read_header

    subroutine print_header (header)
    type (t_stat) ,intent(in) :: header
      write (6,'(a,a)')     'id  : ',header% id
      write (6,'(a,f10.3)') 'lat : ',header% lat
      write (6,'(a,f10.3)') 'lon : ',header% lon
      write (6,'(a,f10.3)') 't   : ',header% t
      write (6,'(a,i6)')    'nlev: ',header% nlev
      write (6,'(a,i6)')    'flag: ',header% flag
    end subroutine print_header

    subroutine read_data (data)
    real(wp) ,intent(out) :: data (:)
      real(sp) :: xsp
      integer  :: j
      do j=1,size(data)
        ctl% records = ctl% records + 1
        read (iunit, rec = ctl% records) xsp
        data(j) = xsp
      end do
    end subroutine read_data

  end subroutine print_stat
!------------------------------------------------------------------------------
  subroutine write_var2 (ctl, x, name, t, z, comment, yrev, zrev, iostat)
  !---------------------------------
  ! write a data set to a grads file
  !---------------------------------
  type (t_ctl)     ,intent(inout)        :: ctl     ! grads ctl structure
  real(wp)         ,intent(in)           :: x (:,:) ! field to write
  character(len=*) ,intent(in)           :: name    ! name of data set
  integer          ,intent(in) ,optional :: t       ! time slice
  integer          ,intent(in) ,optional :: z       ! level (GRADS order: b->t)
  character(len=*) ,intent(in) ,optional :: comment ! data set description
  logical          ,intent(in) ,optional :: yrev    ! flip n-s direction
  logical          ,intent(in) ,optional :: zrev    ! passed to write_var3
  integer          ,intent(out),optional :: iostat  ! error return argument

    integer               :: irecl, it, irec, i, j, nx, ny, z0, ios
    logical               :: yf
    real(sp)              :: buf (size(x,1),size(x,2))
    integer(i4)           :: ibuf(size(x,1),size(x,2)) ! aux. integer buffer
    real(wp) ,allocatable :: xx (:,:,:)
    character(len=7)      :: status
    logical               :: exist
    !---------------------------------------------
    ! pass to write_var3 in case of vertical slice
    !---------------------------------------------
    nx = size(x,1)
    ny = size(x,2)
    if ((ctl% xdefn /= nx  .or.  &
         ctl% ydefn /= ny) .and. &
         ctl% xdefn == 1   .and. &
         ctl% ydefn == nx  .and. &
         ctl% zdefn == ny        )  then
      allocate (xx (1,nx,ny))
      xx(1,:,:) = x
      call  write_var3 (ctl, xx, name, t, z, comment, yrev, zrev, iostat)
      deallocate (xx)
      return
    endif
    !------------------------------
    ! process (optional) parameters
    !------------------------------
    it = 1         ;if(present(t))     it = t
    yf = ctl% yrev ;if(present(yrev))  yf = yrev .neqv. ctl% yrev
    z0 = 0         ;if(present(z))     z0 = z - 1
    if(present(iostat)) iostat = 0
    if(present(zrev)) call finish('write_var2','zrev not processed')
    !---------------------
    ! check shape of array
    !---------------------
    if (any (shape(x) /= (/ctl% xdefn, ctl% ydefn/))) then
      if(present(iostat)) then
        iostat = 1
        return
      endif
      print *,'write_var2: finish, shape(x) /= (/xdefn,ydefn/):',name
      print *,'write_var2: finish,',shape(x),'/=',ctl% xdefn, ctl% ydefn
      print *,'write_var2: ctl% zdefn',ctl% zdefn
      call finish('write_var2','shape(x) /= (/xdefn,ydefn/): '//name)
    endif
    !---------------------------
    ! add entry to ctl structure
    !---------------------------
    inquire (file=ctl% path, exist=exist)
    status =                        'old'
    if (ctl% records == 0) status = 'replace'
    if (.not.exist)        status = 'new'
    if (it==1 .and. .not. present(z)) &
      call add_var (ctl, name, levels=1, comment=comment)
    ctl% tdefn = max (ctl% tdefn, it)
    !------------------------
    ! determine record number
    !------------------------
    irec = 0
    do i=1,ctl% vars
      if(ctl%var(i)% name == name) then
        irec = ctl%var(i)% records
        exit
      endif
    end do
    if (irec == 0) then
      print *,'write_var2: finish, variable ',trim(name),' is not in ctl list.'
      call write_ctl(ctl, unit=6)
      call finish('write_var2','variable is not in ctl list: '//name)
    endif
    irec = irec + z0 + (it-1) * ctl% records
    !----------------------------------------
    ! convert to storage precision, flip axis
    !----------------------------------------
    if(yf) then
      do i=1,ny
        buf(:,i) = x(:,ny-i+1)
      end do
    else
      buf = x
    endif
    !----------
    ! open file
    !----------
    if (ctl% iunit < 0) then
      inquire (iolength = irecl) ibuf
      ctl% iunit = get_unit_number()
      open (ctl% iunit ,file=ctl% path ,                          &
            access='direct' ,recl=irecl, status=status, iostat=ios)
      if(ios/=0) call finish('write_var2','cannot open '//trim(ctl%path))
    endif
    !-------------------------------------------------------
    ! Work around "feature" of the Intel compiler with -O0
    ! when converting reals from little endian to big endian
    !-------------------------------------------------------
    do j = 1, size (x,2)
       ibuf(:,j) = transfer (buf(:,j), ibuf)
    end do
    !-------------------
    ! flip storage order
    !-------------------
    if (byteswapped(ctl)) call flip (ibuf)
    !---------------
    ! write data set
    !---------------
    write (ctl% iunit,rec=irec) ibuf
  end subroutine write_var2
!------------------------------------------------------------------------------
  subroutine write_var3 (ctl, x, name, t, z, comment, yrev, zrev, iostat)
  !---------------------------------
  ! write a data set to a grads file
  !---------------------------------
  type (t_ctl)     ,intent(inout)        :: ctl       ! grads ctl structure
  real(wp)         ,intent(in)           :: x (:,:,:) ! field to write
  character(len=*) ,intent(in)           :: name      ! name of data set
  integer          ,intent(in) ,optional :: t         ! time slice
  integer          ,intent(in) ,optional :: z         ! lowestlevel(GRADSorder)
  character(len=*) ,intent(in) ,optional :: comment   ! data set description
  logical          ,intent(in) ,optional :: yrev      ! x passed N   to S
  logical          ,intent(in) ,optional :: zrev      ! x passed top to bot.
  integer          ,intent(out),optional :: iostat

    integer          :: irecl, it, irec, i, j, ny, nz, z0, ios
    logical          :: yf, zf
    real(sp)         :: buf (size(x,1),size(x,2),size(x,3))
    integer(i4)      :: ibuf(size(x,1),size(x,2))        ! aux. integer buffer
    character(len=7) :: status
    logical          :: exist

    it = 1         ;if(present(t))     it = t
    yf = ctl% yrev ;if(present(yrev))  yf = yrev .neqv. ctl% yrev
    zf = ctl% zrev ;if(present(zrev))  zf = zrev .neqv. ctl% zrev
    z0 = 0         ;if(present(z))     z0 = z - 1
    ny = size(x,2)
    nz = size(x,3)
    if(present(iostat)) iostat = 0
    !---------------------
    ! check shape of array
    !---------------------
    if (any (shape(x(:,:,1)) /= (/ctl% xdefn, ctl% ydefn/))) then
      if(present(iostat)) then
        iostat = 1
        return
      endif
      write(6,*) 'write_var3: shape(x) /= (/xdefn,ydefn/)'
      write(6,*) '            shape(x)        =',shape(x)
      write(6,*) '            (/xdefn,ydefn/) =',ctl% xdefn, ctl% ydefn
      call finish('write_var3','shape(x) /= (/xdefn,ydefn/): '//name)
    endif
    !---------------------------
    ! add entry to ctl structure
    !---------------------------
    inquire (file=ctl% path, exist=exist)
    status =                        'old'
    if (ctl% records == 0) status = 'replace'
    if (.not.exist)        status = 'new'
    if (it==1 .and. .not. present(z)) &
      call add_var (ctl, name, levels=size(x,3), comment=comment)
    ctl% tdefn = max (ctl% tdefn, it)
    !------------------------
    ! determine record number
    !------------------------
    irec = 0
    do i=1,ctl% vars
      if(ctl%var(i)% name == name) then
        irec = ctl%var(i)% records
        exit
      endif
    end do
    if (irec == 0) call finish('write_var3',                        &
                               'variable is not in ctl list: '//name)
    irec = irec + z0 + (it-1) * ctl% records
    !----------------------------------------
    ! convert to storage precision, flip axes
    !----------------------------------------
    if(yf) then
      do i=1,ny
        buf(:,i,:) = x(:,ny-i+1,:)
      end do
    else
      buf = x
    endif
    if (zf) then
      buf(:,:,:) = buf(:,:,nz:1:-1)
    endif
    !----------
    ! open file
    !----------
    if (ctl% iunit < 0) then
      inquire (iolength = irecl) ibuf(:,:)
      ctl% iunit = get_unit_number()
      open (ctl% iunit ,file=ctl% path , &
            access='direct', recl=irecl, &
#ifdef HAVE_ASYNC_IO
            asynchronous='yes',          &
#endif
            status=status, iostat=ios    )
      if(ios/=0) call finish('write_var3','cannot open '//trim(ctl%path))
    endif
    !---------------
    ! write data set
    !---------------
    do i=1,size(x,3)
      !-------------------------------------------------------
      ! Work around "feature" of the Intel compiler with -O0
      ! when converting reals from little endian to big endian
      !-------------------------------------------------------
      do j = 1, size (x,2)
         ibuf(:,j) = transfer (buf(:,j,i), ibuf)
      end do
      !-------------------
      ! flip storage order
      !-------------------
      if (byteswapped(ctl)) call flip (ibuf)
      write (ctl% iunit,rec=irec) ibuf(:,:)
      irec=irec+1
    end do
  end subroutine write_var3
!------------------------------------------------------------------------------
  function record (ctl, name, t, z, index)
  integer :: record
  type (t_ctl)     ,intent(in)            :: ctl     ! grads ctl structure
  character(len=*) ,intent(in)            :: name    ! name of data set
  integer          ,intent(in)  ,optional :: t       ! time slice
  integer          ,intent(in)  ,optional :: z       ! level (GRADS order: b->t)
  integer          ,intent(out) ,optional :: index   ! var index in .ctl
  !------------------------
  ! determine record number
  !------------------------
    integer :: i, irec, it, z0

    it   = 1; if(present(t)) it = t
    z0   = 0; if(present(z)) z0 = z - 1
    irec = 0
    i    = 0
    do i=1, ctl% vars
      if (ctl%var(i)% name == name) then
        irec = ctl%var(i)% records
        exit
      endif
    end do
    if (irec /= 0) irec = irec + z0 + (it-1) * ctl% records
    record = irec
    if (present(index)) index = i
  end function record
!------------------------------------------------------------------------------
  subroutine read_var2 (ctl, x, name, t, z, yrev, iostat)
  !------------------------------
  ! Read 2d field from GrADS file
  !------------------------------
  type (t_ctl)     ,intent(in)           :: ctl     ! grads ctl structure
  real(wp)         ,intent(out)          :: x (:,:) ! field to read
  character(len=*) ,intent(in)           :: name    ! name of data set
  integer          ,intent(in) ,optional :: t       ! time slice
  integer          ,intent(in) ,optional :: z       ! level (GRADS order: b->t)
  logical          ,intent(in) ,optional :: yrev    ! flip n-s direction
  integer          ,intent(out),optional :: iostat  ! io stat return variable

    integer      :: iunit, irecl, irec, i, j, ny
    real(sp)     :: buf (size(x,1),size(x,2))
    integer(i4)  :: ibuf(size(x,1),size(x,2))       ! aux. integer buffer
    logical      :: yf
    yf = ctl% yrev ;if(present(yrev))  yf = yrev .neqv. ctl% yrev
    if(present(iostat)) iostat = 0
    !---------------------
    ! check shape of array
    !---------------------
    if (any (shape(x) /= (/ctl% xdefn, ctl% ydefn/))) then
      write(6,*) 'read_var2: shape(x) /= (/xdefn,ydefn/)'
      write(6,*) '           shape(x)        =',shape(x)
      write(6,*) '           (/xdefn,ydefn/) =',ctl% xdefn, ctl% ydefn
      call finish('read_var2','shape(x) /= (/xdefn,ydefn/): '//name)
    endif
    !------------------------
    ! determine record number
    !------------------------
    irec = record (ctl, name, t, z, index=i)
    if (irec == 0) then
      !---------------------
      ! variable not present
      !---------------------
      if(present(iostat)) then
        iostat = 1
        return
      else
        write(6,*) 'read_var2: variable is not in ctl list.'
        write(6,*) '           variable =',name
        do i=1,ctl% vars
          write(6,*) '          listentry =',ctl%var(i)% name,&
                                             ctl%var(i)% name == name
        end do
      endif
      call finish('read_var2','variable is not in ctl list.')
    endif
    !--------------
    ! read data set
    !--------------
    inquire (iolength = irecl) ibuf
    iunit = get_unit_number()
    open (iunit ,file=ctl%path ,access='direct' ,recl=irecl ,status='old')
    !------------------------------------------------------
    ! Work around "feature" of the Intel compiler with -O0:
    ! Read data into temporary array of integers
    ! (instead of directly to array of reals),
    ! flip byte order of integer *before* transfer to real.
    !------------------------------------------------------
    read(iunit,rec=irec) ibuf
    close(iunit)
    call return_unit_number(iunit)
    if (byteswapped(ctl)) call flip (ibuf)
    do j = 1, size (x,2)
       buf(:,j) = transfer (ibuf(:,j), buf)
    end do
    !-------------------
    ! flip storage order
    !-------------------
    if(yf) then
      ny = size(x,2)
      do i=1,ny
        x(:,ny-i+1) = buf(:,i)
      end do
    else
      x = buf
    endif
  end subroutine read_var2
!------------------------------------------------------------------------------
  subroutine read_var3 (ctl, x, name, t, z, yrev, zrev, iostat)
  !------------------------------
  ! Read 3d field from GrADS file
  !------------------------------
  type (t_ctl)     ,intent(in)           :: ctl       ! grads ctl structure
  real(wp)         ,intent(out)          :: x (:,:,:) ! field to read
  character(len=*) ,intent(in)           :: name      ! name of data set
  integer          ,intent(in) ,optional :: t         ! time slice
  integer          ,intent(in) ,optional :: z         ! lowestlevel(GRADSorder)
  logical          ,intent(in) ,optional :: yrev      ! flip n-s   direction
  logical          ,intent(in) ,optional :: zrev      ! flip vert. direction
  integer          ,intent(out),optional :: iostat    ! io stat return variable

    integer      :: iunit, irecl, irec, i, j, ny, nz
    real(sp)     :: buf (size(x,1),size(x,2),size(x,3))
    integer(i4)  :: ibuf(size(x,1),size(x,2))         ! aux. integer buffer
    logical      :: yf, zf
    yf = ctl% yrev ;if(present(yrev))  yf = yrev .neqv. ctl% yrev
    zf = ctl% zrev ;if(present(zrev))  zf = zrev .neqv. ctl% zrev
    if(present(iostat)) iostat = 0
    !------------------------
    ! determine record number
    !------------------------
    irec = record (ctl, name, t, z, index=i)
    if (irec == 0) then
      !---------------------
      ! variable not present
      !---------------------
      if(present(iostat)) then
        iostat = 1
        return
      else
        write(6,*) 'read_var3: variable is not in ctl list: ',trim(name)
        do i=1,ctl% vars
          write(6,*) '          listentry =',ctl%var(i)% name,&
                                             ctl%var(i)% name == name
        end do
      endif
      call finish('read_var3','variable is not in ctl list: '//trim(name))
    endif
    !---------------------
    ! check shape of array
    !---------------------
    if (any (shape(x) /= (/ctl% xdefn, ctl% ydefn, ctl%var(i)% levels/))) then
      write(6,*) 'read_var3: shape(x) /= (/xdefn,ydefn,var%levels/)'
      write(6,*) '           shape(x) = ',shape(x)
      write(6,*) '           ctl      : ',&
                 ctl% xdefn, ctl% ydefn, ctl%var(i)%levels
      call finish('read_var3','shape(x) /= (/xdefn,ydefn,var%levels/): '//name)
    endif
    !--------------
    ! read data set
    !--------------
    inquire (iolength = irecl) ibuf(:,:)
    iunit = get_unit_number()
    open (iunit ,file=ctl%path ,access='direct' ,recl=irecl ,status='old')
    do i=1,size(x,3)
      !------------------------------------------------------
      ! Work around "feature" of the Intel compiler with -O0:
      ! Read data into temporary array of integers
      ! (instead of directly to array of reals),
      ! flip byte order of integer *before* transfer to real.
      !------------------------------------------------------
      read(iunit,rec=irec) ibuf(:,:)
      irec=irec+1
      if (byteswapped(ctl)) call flip (ibuf)
      do j = 1, size (x,2)
         buf(:,j,i) = transfer (ibuf(:,j), buf)
      end do
    end do
    close(iunit)
    call return_unit_number(iunit)
    !-------------------
    ! flip storage order
    !-------------------
    if (zf) then
      nz = size(x,3)
      buf(:,:,:) = buf(:,:,nz:1:-1)
    endif
    if(yf) then
      ny = size(x,2)
      do i=1,ny
        x(:,ny-i+1,:) = buf(:,i,:)
      end do
    else
      x = buf
    endif
  end subroutine read_var3
!------------------------------------------------------------------------------
  subroutine write_ctl (ctl, unit)
  type (t_ctl) ,intent(inout)        :: ctl
  integer      ,intent(in) ,optional :: unit
    integer            :: i, ios, iu
    character(len=128) :: file
    !----------------------
    ! open file if required
    !----------------------
    if (present(unit)) then
      iu = unit
    else
      iu = get_unit_number()
      i = index(ctl% path, '.dat ')
      if (i/=0) then
        file = ctl%path(1:i)//'ctl'
      else
        file = trim(ctl%path)//'.ctl'
      endif
      open (iu,file=file,action='write',iostat=ios)
      if(ios/=0) call finish('write_ctl','cannot open '//trim(file))
    endif
    !--------------
    ! write ctl ...
    !--------------
    if(iu==6) then
      write(iu,'(a,1x,a)')                    'path ' ,trim(ctl% path)
    endif
    write(iu,'(a,1x,a)')                      'dset ' ,trim(ctl% dset)
    if (associated (ctl% comment)) then
      do i = 1, size (ctl% comment)
        if(ctl% comment(i)/=' ') &
                         write(iu,'(a,1x,a)') '*'     ,trim(ctl% comment(i))
      end do
    end if
    if(ctl% dtype  /='') write(iu,'(a,1x,a)') 'dtype' ,trim(ctl% dtype)
    if(ctl% stnmap /='') write(iu,'(a,1x,a)') 'stnmap',trim(ctl% stnmap)
    if(ctl% title  /='') write(iu,'(a,1x,a)') 'title' ,trim(ctl% title)

    if(ctl% sequential .or. ctl% yrev .or. ctl% little_endian .or. &
       ctl% big_endian .or. ctl% byteswapped) then
      write                       (iu,'(a)',advance='no') 'options'
      if(ctl% sequential)    write(iu,'(a)',advance='no') ' sequential'
      if(ctl% yrev)          write(iu,'(a)',advance='no') ' yrev'
      if(ctl% zrev)          write(iu,'(a)',advance='no') ' zrev'
      if(ctl% little_endian) write(iu,'(a)',advance='no') ' little_endian'
      if(ctl% big_endian)    write(iu,'(a)',advance='no') ' big_endian'
      if(ctl% byteswapped)   write(iu,'(a)',advance='no') ' byteswapped'
      write                    (iu,'()')
    endif

    if(ctl% fileheader/=0) write(iu,'(a,1x,i8)') 'fileheader', ctl% fileheader
    if(ctl% theader   /=0) write(iu,'(a,1x,i8)') 'theader   ', ctl% theader
    if(ctl% xyheader  /=0) write(iu,'(a,1x,i8)') 'xyheader  ', ctl% xyheader

    write(iu,'(a,1x,g12.6)')          'undef',ctl% undef
    select case (ctl% dtype)
    !---------------------------------------------
    ! write xdef, ydef, zdef for gridded data only
    !---------------------------------------------
    case ('')
      select case (ctl%xdef)
      case ('linear','LINEAR')
         write(iu,'(a,1x,i6,1x,a,2g13.6)') 'xdef' ,ctl%xdefn,ctl%xdef,&
                                                 ctl%xdefi(1),ctl%xdefi(2)
      case ('levels','LEVELS')
         write(iu,'(a,1x,i6,1x,a,2g13.6)') 'xdef' ,ctl%xdefn,ctl%xdef
         write(iu,'(15x,6f8.3)') real (ctl%xlevels,sp)
      end select
      select case (ctl%ydef)
      case ('linear','LINEAR')
        write(iu,'(a,1x,i6,1x,a,2g13.6)') 'ydef' ,ctl%ydefn,ctl%ydef,&
                                                 ctl%ydefi(1),ctl%ydefi(2)
      case ('levels','LEVELS')
        write(iu,'(a,1x,i6,1x,a,2g13.6)') 'ydef' ,ctl%ydefn,ctl%ydef
        write(iu,'(15x,6f8.3)') real(ctl%ylevels,sp)
      end select
      select case (ctl%zdef)
      case ('linear','LINEAR')
        write(iu,'(a,1x,i6,1x,a,2g13.6)') 'zdef' ,ctl%zdefn,ctl%zdef,&
                                                 ctl%zdefi(1),ctl%zdefi(2)
      case ('levels','LEVELS')
        write(iu,'(a,1x,i6,1x,a,2g13.6)') 'zdef' ,ctl%zdefn,ctl%zdef
        write(iu,'(15x,6f9.3)') real(ctl%zlevels,sp)
      end select
    end select

    write(iu,'(a,1x,i6,3(1x,a))')     'tdef' ,ctl%tdefn,ctl%tdef,&
                                             ctl%tdefi(1),ctl%tdefi(2)
    write(iu,'(a,1x,i4)')             'vars' ,ctl%vars
    do i=1,ctl%vars
      write(iu,"(a,1x,i4,1x,a,1x,a)") &
        ctl%var(i)%name, ctl%var(i)%levels,         &
        trim(ctl%var(i)%units),trim(ctl%var(i)%comment)
    end do
    write(iu,'(a)')                   'endvars'
    !------------------------
    ! close files if required
    !------------------------
    if (.not.present(unit)) then
      close (iu)
      call return_unit_number(iu)
    endif
    if (ctl% iunit >= 0) then
      close (ctl% iunit)
      call return_unit_number (ctl% iunit)
      ctl% iunit = -1
    endif
  end subroutine write_ctl
!------------------------------------------------------------------------------
  subroutine read_ctl (ctl, file, iostat, debug)
  type (t_ctl)     ,intent(inout)         :: ctl
  character(len=*) ,intent(in)            :: file
  integer          ,intent(out) ,optional :: iostat
  logical          ,intent(in)  ,optional :: debug
    type(t_ctl)                 :: empty
    integer                     :: iu, ios, idx, idx2, i, n
    character(len=128)          :: line
    character(len=16)           :: code
    character(len=128)          :: path
    character(len=16)           :: words(32)
    character(len=78)  ,pointer :: comment(:)
    logical                     :: dbg
    dbg = .false.; if (present (debug)) dbg = debug
    iu = get_unit_number()
    idx = index(file, '^')
    if (idx>0) then
      path = file (:idx-1)//file(idx+1:)
    else
      path = file
    endif
    open (iu,file=trim(path),action='read',status='old',iostat=ios)
    if(present(iostat)) iostat = ios
    if(ios/=0 .and. .not.present(iostat)) &
      call finish('read_ctl','cannot open '//trim(path))
    if(associated(ctl% var    )) deallocate (ctl% var)
    if(associated(ctl% ylevels)) deallocate (ctl% ylevels)
    ctl = empty
    do
      read(iu,'(a)',iostat=ios) line
      if(ios/=0)         exit
      if (dbg) write (0,*) 'read_ctl: line=', trim (line)
      if(line(1:1)=='*') then
        if (associated(ctl% comment)) then
          comment => ctl% comment
          allocate (ctl% comment (size(comment)+1))
          ctl% comment(1:size(comment)  ) = comment
          ctl% comment(  size(comment)+1) = line(3:)
          deallocate (comment)
        else
          allocate (ctl% comment(1))
          ctl% comment(1) = line(3:)
        endif
        cycle
      endif
      idx=index(line,' ')
      if(idx>0) then
        code = tolower (line(:idx-1))
        select case (code)
        case ('dset','DSET')
          idx2 = verify (line(idx:),' ')
          ctl% dset = line(idx+idx2-1:)
          if (ctl% dset(1:1)=='^') then
            idx2 = index (path,'/',.true.)
            if (idx2 > 0) then
              ctl% path = path(:idx2)//ctl% dset(2:)
            else
              ctl% path = ctl% dset(2:)
            endif
          else
            ctl% path = ctl% dset
          endif
        case ('undef','UNDEF')
          read(line(idx+1:),*) ctl% undef
        case ('dtype','DTYPE')
          read(line(idx+1:),*) ctl% dtype
        case ('stnmap','STNMAP')
          read(line(idx+1:),*) ctl% stnmap
        case ('title','TITLE')
          read(line(idx+1:),*) ctl% title
        case ('options','OPTIONS')
          ctl% sequential    = index (line,'sequential')    /= 0
          ctl% yrev          = index (line,'yrev')          /= 0
          ctl% zrev          = index (line,'zrev')          /= 0
          ctl% little_endian = index (line,'little_endian') /= 0
          ctl% big_endian    = index (line,'big_endian')    /= 0
          ctl% byteswapped   = index (line,'byteswapped')   /= 0
        case ('xdef','XDEF')
          read(line(idx+1:),*) ctl%xdefn,ctl%xdef,ctl%xdefi(1),ctl%xdefi(2)
        case ('ydef','YDEF')
          read(line(idx+1:),*) ctl%ydefn,ctl%ydef
          select case (ctl%ydef)
          case ('linear','LINEAR')
            read(line(idx+1:),*) ctl%ydefn,ctl%ydef,ctl%ydefi(1),ctl%ydefi(2)
          case ('levels','LEVELS')
            allocate (ctl% ylevels (ctl%ydefn))
            ! Handle level values on the line after declaration
            call split (words, line(idx+1:), n)
            if (dbg) write (0,*) "read_ctl: n =", n
            if (n >= 3) then
               ! Get the whole bunch...
               do i = 3, n
                  if (dbg) write (0,*) "read_ctl: next ylevel:", words(i)
                  read (words(i),*) ctl% ylevels(i-2)
               end do
               read(iu,*) ctl% ylevels(n-2+1:)
            else
               read(iu,*) ctl% ylevels(:)
            end if
          end select
        case ('zdef','ZDEF')
          read(line(idx+1:),*) ctl%zdefn,ctl%zdef
          select case (ctl%zdef)
          case ('linear','LINEAR')
            read(line(idx+1:),*) ctl%zdefn,ctl%zdef,ctl%zdefi(1),ctl%zdefi(2)
          case ('levels','LEVELS')
            allocate (ctl% zlevels (ctl%zdefn))
            read(iu,*) ctl% zlevels
          end select
        case ('tdef','TDEF')
          read(line(idx+1:),*) ctl%tdefn,ctl%tdef,ctl%tdefi(1),ctl%tdefi(2)
        case ('vars','VARS')
          read(line(idx+1:),*) ctl%vars
          allocate (ctl%var(ctl%vars))
          ctl% records = 0
          ctl% svars   = 0
          ctl% lvars   = 0
          do i=1,ctl%vars
            read(iu,'(a)',iostat=ios) line
            if(ios/=0) exit
            ctl%var(i)%units     = '99'
            ctl% var(i)% comment = ''
            read(line,*,iostat=ios) ctl%var(i)%name, ctl%var(i)%levels
            if (ctl%var(i)%levels == 0) then
              ctl%svars = ctl%svars + 1
            else
              ctl%lvars = ctl%lvars + 1
            endif
            ctl%var(i)%levels = max(1,ctl%var(i)%levels)
            idx =       verify(line         ,' ')  ! start of name
            idx = idx + index (line(idx +1:),' ')  ! end of name
            idx = idx + verify(line(idx +1:),' ')  ! start of levels
            idx = idx + index (line(idx +1:),' ')  ! end of levels
            idx = idx + verify(line(idx +1:),' ')  ! start of units
            idx2= idx + index (line(idx +1:),' ')  ! end of units
            if(idx2>idx) ctl%var(i)%units = line(idx:idx2-1)
            idx = idx2+ verify(line(idx2+1:),' ')  ! start of comment
            ctl% var(i)% comment = line(idx:)
            ctl% var(i)% records = ctl% records + 1
            ctl% records = ctl% records + ctl% var(i)% levels
          end do
        case ('endvars','ENDVARS')
          exit
        case ('index')
           print *, 'read_ctl: *** INDEX dataset declaration not supported ***'
           print *, 'read_ctl: ',trim(line)
           cycle
        case default
          print *,'read_ctl: UNKNOWN LINE:',trim(line)
        end select
      endif
    end do
    close (iu)
    call return_unit_number(iu)
  end subroutine read_ctl
!------------------------------------------------------------------------------
  function byteswapped (ctl)
  type (t_ctl) ,intent(in) :: ctl
  logical                  :: byteswapped
    byteswapped = ctl% byteswapped
    if (ctl% big_endian)    byteswapped = little()
    if (ctl% little_endian) byteswapped = .not. little()
  end function byteswapped
!------------------------------------------------------------------------------
end module mo_grads
