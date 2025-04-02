!
!+ Read blacklist to reject suspicious observations
!
MODULE mo_blacklist
!
! Description:
!   This module holds routines 'read_blacklists' to read blacklisting
!   data and 'black_entry' to access the data concerning a specific
!   station. The Names and Paths of the blacklist files to be read are
!   specified by namelist '/BLACKLIST/'. The information is kept in a
!   private module variable of type 't_black'.
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
! V1_8         2009/12/09 Andreas Rhodin
!  precautions to read very old blacklist format (geop. and wind only)
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_11        2010/06/16 Harald Anlauf
!  /BLACKLIST/: add "verbose" for verbosity of blacklist processing
! V1_13        2011/11/01 Andreas Rhodin
!  allow multiple blacklist entries per station id; implement whitelist
! V1_15        2011/12/06 Andreas Rhodin
!  fix bugs in whitelist processing
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Andreas Rhodin DWD   2002-2008  original code
!------------------------------------------------------------------------------

  !=============
  ! Modules used
  !=============
  use mo_kind,          only: i8,wp               ! integer*8 kind parameter
  use mo_endian,        only: little,            &! .true. for little endian
                              flip                ! flip integer field
  use mo_namelist,      only: position_nml,      &! position namelist
                              POSITIONED,        &! return value
                              nnml                ! namelist file unit number
  use mo_fortran_units, only: get_unit_number,   &! obtain a free unit number
                              return_unit_number
  use mo_exception,     only: finish, message     ! abort, warning routine
  use mo_mpi_dace,      only: dace,              &! MPI group info
                              p_bcast             ! generic MPI broadcast
  use mo_run_params,    only: input,             &! path for input
                              path_file           ! add a path to a file name
  use mo_obs_tables,    only: obstyp              ! relate names to obstypes
  use mo_usstd,         only: h_p_usstd           ! h(gpm) from p US std.atm
  use mo_wigos,         only: t_wsi,             &! WIGOS id data type
                              wsi_mode,          &! WIGOS station id mode
                              wsi_prefix,        &! Prefix of "fake station id"
                              wsi_decode,        &! Derive WSI from given string
                              wsi_encode          ! statid,stathash from t_wsi
  implicit none
  !================
  ! Public entities
  !================
  private
!EOX
!
! !PUBLIC TYPES:
!
  public :: t_black         ! blacklist entry data type
  public :: t_white         ! whitelist entry data type
  public :: p_white         ! pinter to t_white array
  public :: t_blacklists    ! container for lists of type t_black
!
! !PUBLIC MEMBER FUNCTIONS:
!
  public :: read_blacklists ! read blacklists
  public :: black_entry     ! find entry in blacklist
  public :: clean_blacklist ! deallocate blacklists
  public :: empty_blacklist ! create empty blacklists
  public :: set_mask        ! set mask in blacklist entry
!EOP
  public :: lists           ! blacklists
  public :: verbose         ! verbosity level of blacklisting
  !----------------------
  ! data type definitions
  !----------------------
!BOP
!
! !DATATYPE: t_black
!
! !DESCRIPTION:
!
! The Data read from the blacklisting files is stored in an array of data type
! {\bf t\_blacklists}.
!
! Each element of the array holds information on a specific observation
! type. Observation type number is stored in the component {\bf type} of
! type {\bf t\_black}. The blacklisting information is stored in the
! array component {\bf list} whose size is given by {\bf n}.
!
! The components of {\bf t\_black} reflect the information given in the
! blacklisting files: {\bf statid} holds the name of the station (with
! {\tt "."} used as wildcards). {\bf id} and {\bf mask} are integer
! representations of the content of {\bf statid} used internally for
! efficient wildcard processing. Components {\bf pg}, {\bf pw}, {\bf
! pt}, and {\bf pd} hold pressure limits for height, wind, temperature,
! and dewpoint, respectively.

  !------------------------
  ! Derived type definition
  !------------------------
  type t_black
    character(len=8)  :: statid = ''    ! station id (wildcard)
    integer(i8)       :: mask   = 0     ! mask (for wildcards)
    integer(i8)       :: id     = 0     ! station id (masked)
    integer           :: len    = 0     ! len_trim (station id)
    integer           :: pg(2)  = 0     ! pressure limits for height      (Pa)
    integer           :: pw(2)  = 0     ! pressure limits for wind        (Pa)
    integer           :: pt(2)  = 0     ! pressure limits for temperature (Pa)
    integer           :: pd(2)  = 0     ! pressure limits for dewpoint    (Pa)
    real(wp)          :: zg(2)  = 0._wp ! height   limits for height      (gpm)
    real(wp)          :: zw(2)  = 0._wp ! height   limits for wind        (gpm)
    real(wp)          :: zt(2)  = 0._wp ! height   limits for temperature (gpm)
    real(wp)          :: zd(2)  = 0._wp ! height   limits for dewpoint    (gpm)
    logical           :: lg     =.false.! geopotental blacklisted
    logical           :: lw     =.false.! wind        blacklisted
    logical           :: lt     =.false.! temperature blacklisted
    logical           :: ld     =.false.! dewpoint    blacklisted
    logical           :: mwl    =.false.! missing whitelist entry
    real              :: bcg    = 0.    ! geopotental bias correction
  end type t_black

  type t_white
    character(len=8)  :: statid   = ''    ! station id
    integer(i8)       :: id       = 0     ! station id
  end type t_white

  type p_white
    integer                :: codetype = 0       ! code type
    integer                :: n        =  0      ! number of items in whitelist
    type(t_white) ,pointer :: p  (:)   => NULL() ! whitelist (codetype)
  end type p_white

  type t_blacklists
    integer                :: type     =  0      ! observation type number
    integer                :: n        =  0      ! number of items in blacklist
    integer                :: nw       =  0      ! number of items in whitelist
    type(t_black) ,pointer :: list (:) => NULL() ! blacklist
    type(p_white) ,pointer :: wlist(:) => NULL() ! whitelist (codetype)
  end type t_blacklists

!EOP
  !=================
  ! Module variables
  !=================
  !----------
  ! constants
  !----------
  integer ,parameter :: mblen = 10000 ! estimated max length of blacklists
  integer ,parameter :: mty   =    20 ! estimated max number of observ. types
  logical ,save      :: bigend= .true.! true for big endian
  !----------
  ! variables
  !----------
  type(t_blacklists) ,pointer :: lists (:) => NULL()
  integer                     :: nlists    = 0
  !---------
  ! namelist
  !---------
!BOP
!
! !NAMELIST: /BLACKLIST/
!
! !DESCRIPTION:
!
! By namelist group {\bf /BLACKLIST/} the path and name of the blacklist
! files are specified by variables {\bf path} and {\bf name},
! respectively. The namelist group may be called repeatedly for
! different blacklisting files. In this case the path defined in
! previous calls remains valid. The logical variable {\bf default}
! may be set to {\bf .false.} to supress reading of the default
! blacklist (\/io\/input\/blacklist\_yyyymmddhh).
!
! !DEFINITION:
!
  character(len=128) :: file
  character(len=128) :: path
  logical            :: default
  integer            :: verbose

  namelist /BLACKLIST/ path, file, default, verbose
!
!EOP

!!        T P.G.U P.G.O P.W.U P.W.O P.T.U P.T.O P.D.U P.D.O

!0010     1  1100     0     0     0     0     0     0     0
!04045    1  1100     0     0     0     0     0     0     0
!0XRA6    1  1100     0     0     0     0     0     0     0

!......   2     0     0    59    10    59    10    59    10
!034BATBA 2     0     0     0     0  1100     0     0     0
!0EVEIEBA 2     0     0     0     0  1100     0     0     0
!ABX...   2     0     0  1100     0     0     0     0     0
!ACHUIEBA 2     0     0     0     0  1100     0     0     0


!14904    4  1100     0     0     0     0     0     0     0
!14914    4  1100     0     0     0     0     0     0     0

!03023    5  1100     0     0     0  1100     0     0     0
!ZDLP     5  1100     0     0     0  1100     0     0     0

!13586    6     0     0  1100     0     0     0     0     0
!15480    6     0     0  1100     0     0     0     0     0
!17240    6     0     0  1100     0     0     0     0     0
!31960    6     0     0  1100     0     0     0     0     0

! 1 SYNOP
! 2 AIREP
! 3 SATOB
! 4 BOJE
! 5 TEMP
! 6 PILOT
! 7 SATEM


contains
!-------------------------------------------------------------------

!BOP
! !IROUTINE: read_blacklists
!
! !DESCRIPTION:
!
! Blacklisting files are read from the files specified by namelist {\bf
! /BLACKLIST/}. In a parallel environment the namelists and blacklisting
! files are read by the I/O processor element only and then broadcasted
! to the other processor elements. The information is stored in a module
! variable and may be accessed by calling subroutine {\bf black\_entry}
! later.
!
! !INTERFACE:
!
  subroutine read_blacklists

!EOP
    !
    ! local variables
    !
    integer           :: i,j,k   ! indices
    logical           :: first   ! flag for first namelist group
    integer           :: ierr    ! namelist position status variable
    integer           :: iu, iuu ! Fortran unit numbers
    integer           :: typ     ! observation type number
    integer           :: code    ! observation code number
    type(t_black)     :: bl      ! one line of the blacklist
    integer           :: nl      ! number of lines read (blacklist)
    type(t_white)     :: wl      ! one line of the whitelist
    integer           :: nlw     ! number of lines read (whitelist)
    integer           :: ios     ! I/O status variable
    character(len=96) :: line    ! line to read
    integer           :: iname   ! 5 or 8 character for station name
    character(len=32) :: wsistr  ! WIGOS station ID as text
    type(t_wsi)       :: wsi     ! WIGOS station ID (decoded)
    integer(i8)       :: wsihash ! WIGOS station ID (internal hash)
    integer           :: status  ! status value
    logical           :: ok      ! status flag
    character(len=10) :: statid  ! "Fake" station ID
    !----------------------
    ! delete old blacklists
    !----------------------
    call clean_blacklist
    bigend = .not.little()
    !---------------------------------
    ! Read blacklists on I/O processor
    !---------------------------------
    if (dace% lpio) then
      allocate (lists(mty))
      nlists = 0
      !---------
      ! printout
      !---------
      write (6,'(a)') repeat ('=',79)
      write (6,'(a)')
      write (6,'(a)') '  Reading Blacklists'
      write (6,'(a)')
      !---------------------------------
      ! read namelist groups /BLACKLIST/
      !---------------------------------
      path    = input
      verbose = 0
      default = .true.
      first   = .true.
      iu = get_unit_number()
      do  ! read namelist group
        call position_nml ('BLACKLIST', lrewind=first, status=ierr)
        first = .false.
        if (ierr == POSITIONED) then
          !------------------------------------------------------------
          ! namelist group positioned: read path/filename from namelist
          !------------------------------------------------------------
          file = ''
#if defined(__ibm__)
          read (nnml ,nml=BLACKLIST, iostat=ierr)
          if (ierr/=0) call finish ('nml_run_flags',               &
                                    'ERROR in namelist /BLACKLIST/')
#else
          read (nnml ,nml=BLACKLIST)
#endif
          if (file == '') cycle
        else
          !----------------------------------------------------------------
          ! no further namelist group positioned: set default path/filename
          !----------------------------------------------------------------
          if (.not. default) exit
          default = .false.
          path = input
          file = 'blklsttmp'
          write (6,*) 'reading default blacklist: ',trim(file)
        endif
        !----------
        ! open file
        !----------
        if (file=='from_namelist') then
          iuu = nnml
          write(6,'(a,a,a,i3)') ' reading: ',trim(file),' ,unit=',iuu
        else
          file = path_file (path, file)
          write(6,'(a,a)') ' reading: ',trim(file)
          iuu = iu
          open (iuu,file=file, status='old', iostat=ios)
          if (ios /= 0) call finish ('read_blacklists', &
                        'blacklist not present: '//trim(file))
        endif
        !-----------
        ! read lines
        !-----------
        nl    = 0
        nlw   = 0
        iname = 8
        do
          read (iuu, '(a)' ,end=999) line
          read (line(1:iname),'(a)') bl%statid
          if (bl%statid(1:6) == 'IIIII '  ) iname = 5
          if (bl%statid(1:8) == 'IIIIIIII') iname = 8
          if (bl%statid(1:5) == 'IIIII'   ) cycle
          if (bl%statid      == ''        ) cycle
          if (bl%statid      == '-EOF-'   ) exit
          if (bl%statid(1:8) == 'WHITELIS') exit  !+ currently ignore whitelist
          if (bl%statid(1:1) == wsi_prefix) then
            if (wsi_mode     == 0         ) cycle ! ignore WIGOS ID entries
            if (iname        /= 8         ) cycle
            wsistr = ""
            if (line(61:) == "") then
              ! Short line w/o WIGOS identifier
              read(line(9:),*,iostat=ios) typ, bl%pg, bl%pw, bl%pt, bl%pd
            else
              read(line(9:),*,iostat=ios) typ, bl%pg, bl%pw, bl%pt, bl%pd, wsistr
            end if
            ok = (ios == 0)
            if (ok) then
              select case (wsistr(1:2))
              case ("  ")
                !------------------------------------------------------------
                ! Special case for blacklisting based on issuer of identifier
                ! or an original fake station identifier by ECMWF.
                ! Handle complete blacklisting of all stations with WSIs.
                !------------------------------------------------------------
                if (wsi_mode        == 3        ) ok = .false.
                if (wsi_mode        /= 1  .and. &
                    bl% statid(5:8) /= "...."   ) ok = .false.
                if (bl% statid(2:8) == ".......") ok = .true.
              case ("0-")
                call wsi_decode (wsistr, wsi, status)
                ok = (status == 0) .and. wsi% wigii < 20000
                if (ok) then
                  call wsi_encode (wsi, statid, wsihash)
                  if (bl% statid(2:8) == ".......") then
                    bl% statid = statid
!                 if (verbose > 0)                                      &
                    write(*,*) "mapping: ", trim (wsistr), " -> ", statid
                  else if (statid /= bl% statid) then
                    write(*,*) "read_blacklists: ", statid, " /= ", bl% statid
                  end if
                  if (wsi_mode /= 1) bl% statid = statid
                end if
              case default
                ok = .false.
              end select
            end if
            if (.not. ok) then
              write(0,'(a)') "read_blacklists:WARNING: invalid line skipped!"
              write(0,'(a)') trim (line)
              cycle
            end if
          else
            read (line(iname+1:60),*,iostat=ios) typ, bl%pg, bl%pw, bl%pt, bl%pd
            if (ios/=0) then    ! read very old format (geop. and wind only)
              read (line(iname+1:60),*) typ, bl%pg, bl%pw
              bl%pt = 0
              bl%pd = 0
            endif
          end if
          bl%len = len_trim(bl%statid)
          bl%pg  = bl%pg * 100
          bl%pw  = bl%pw * 100
          bl%pt  = bl%pt * 100
          bl%pd  = bl%pd * 100
          nl = nl + 1
          !---------------------
          ! derive logical flags
          !---------------------
          bl%lg = (bl%pg(1) /= 0 .and. bl%pg(1) >= bl%pg(2))
          bl%lw = (bl%pw(1) /= 0 .and. bl%pw(1) >= bl%pw(2))
          bl%lt = (bl%pt(1) /= 0 .and. bl%pt(1) >= bl%pt(2))
          bl%ld = (bl%pd(1) /= 0 .and. bl%pd(1) >= bl%pd(2))
          !---------------------
          !compute height levels
          !---------------------
          bl% zg = [-999._wp,99999._wp]
          bl% zw = [-999._wp,99999._wp]
          bl% zt = [-999._wp,99999._wp]
          bl% zd = [-999._wp,99999._wp]
          if (bl% lg) then
            if (bl% pg(1) < 1100) bl% zg(1) = h_p_usstd(real(bl% pg(1),wp))
            if (bl% pg(2) >    0) bl% zg(2) = h_p_usstd(real(bl% pg(2),wp))
          end if
          if (bl% lw) then
            if (bl% pw(1) < 1100) bl% zw(1) = h_p_usstd(real(bl% pw(1),wp))
            if (bl% pw(2) >    0) bl% zw(2) = h_p_usstd(real(bl% pw(2),wp))
          end if
          if (bl% lt) then
            if (bl% pt(1) < 1100) bl% zt(1) = h_p_usstd(real(bl% pt(1),wp))
            if (bl% pt(2) >    0) bl% zt(2) = h_p_usstd(real(bl% pt(2),wp))
          end if
          if (bl% ld) then
            if (bl% pd(1) < 1100) bl% zd(1) = h_p_usstd(real(bl% pd(1),wp))
            if (bl% pd(2) >    0) bl% zd(2) = h_p_usstd(real(bl% pd(2),wp))
          end if
          !----------
          ! find type
          !----------
          j = 0
          do i=1,nlists
            if (lists(i)%type == typ) then; j=i; exit; endif
          end do
          !--------------------------------------
          ! allocate new list if type not present
          !--------------------------------------
          if (j==0) then
            nlists = nlists + 1
            if (nlists > mty) call finish ('read_blacklists','nlists > mty')
            j = nlists
            lists(j)% type = typ
            lists(j)% n    = 0
            lists(j)% nw   = 0
            allocate (lists(j)%  list (mblen))
            allocate (lists(j)% wlist (mty)  )
          endif
          !------------------------------------------
          ! derive integer representation of wildcard
          !------------------------------------------
          call set_mask (bl)
          !-----------
          ! store line
          !-----------
          lists (j)% n = lists (j)% n+ 1
          if (lists (j)% n > mblen) call finish ('read_blacklists','n > mblen')
          lists (j)% list (lists (j)% n) = bl
        end do
        !---------------
        ! read whitelist
        !---------------
        nlw   = 0
        if (bl%statid(1:8) == 'WHITELIS') then
          iname = 8
          do
            read (iuu, '(a)' ,end=999) line
            read (line(1:iname),'(a)') wl% statid
            if (wl%statid == ''      ) cycle
            if (wl%statid == '-EOF-' ) exit
            if (wl%statid(1:1) == wsi_prefix) then
              if (wsi_mode     == 0         ) cycle ! ignore WIGOS ID entries
              wsistr = ""
              if (line(20:) == "") then
                 ! Short line w/o WIGOS identifier
                 read(line(9:),*,iostat=ios) typ, code
              else
                 read(line(9:),*,iostat=ios) typ, code, wsistr
              end if
              ok = (ios == 0)
              if (ok) then
                select case (wsistr(1:2))
                case ("  ")
                  !------------------------------------------------------
                  ! Special case of whitelisting using given fake station
                  ! identifier by ECMWF.  No wildcarding allowed here.
                  !------------------------------------------------------
                  if (wsi_mode        /= 1            ) ok = .false.
                  if (wsi_mode        == 1   .and.    &
                      index (wl% statid(2:8), ".") > 0) ok = .false.
                case ("0-")
                  call wsi_decode (wsistr, wsi, status)
                  ok = (status == 0) .and. wsi% wigii < 20000
                  if (ok) then
                    call wsi_encode (wsi, statid, wsihash)
                    if (wl% statid(2:8) == ".......") then
                      wl% statid = statid
!                     if (verbose > 0)                                      &
                        write(*,*) "mapping: ", trim (wsistr), " -> ", statid
                    else if (statid /= wl% statid) then
                      write(*,*) "read whitelist: ", statid, " /= ", wl% statid
                    end if
                    if (wsi_mode /= 1) wl% statid = statid
                  end if
                case default
                  ok = .false.
                end select
              end if
              if (.not. ok) then
                write(0,'(a)') "read whitelist:WARNING: invalid line skipped!"
                write(0,'(a)') trim (line)
                cycle
              end if
            else
              read (line(iname+1:),*) typ, code
            end if
            nlw = nlw + 1
            !----------
            ! find type
            !----------
            j   = 0
            do i=1,nlists
              if (lists(i)%type == typ) then; j=i; exit; endif
            end do
            !--------------------------------------
            ! allocate new list if type not present
            !--------------------------------------
            if (j==0) then
              nlists = nlists + 1
              if (nlists > mty) call finish ('read_blacklists','nlists > mty')
              j = nlists
              lists(j)% type = typ
              lists(j)% n    = 0
              lists(j)% nw   = 0
              allocate (lists(j)%  list (mblen))
              allocate (lists(j)% wlist (mty)  )
            endif
            !--------------
            ! find codetype
            !--------------
            k = 0
            do i=1,mty
              if (lists(j)% wlist(i)% codetype == code .or. &
                  lists(j)% wlist(i)% codetype == 0         ) then
                k=i
                exit
              endif
            end do
            if (k==0) call finish('read_blacklists','more than mty whitelists')
            !--------------------------------------------
            ! allocate whitelist for codetype if required
            !--------------------------------------------
            if (lists(j)% wlist(k)% codetype == 0) then
              lists(j)% wlist(k)% codetype = code
              allocate (lists(j)% wlist(k)% p (mblen))
              lists(j)% nw = k
            endif
            !--------------------------------------------
            ! derive integer representation of station id
            !--------------------------------------------
            wl% id = transfer  (wl% statid, wl% id)
            if (bigend) call flip (wl% id)
            !-----------
            ! store line
            !-----------
            lists(j)% wlist(k)% n = lists(j)% wlist(k)% n + 1
            i = lists(j)% wlist(k)% n
            if (i > mblen) &
              call finish ('read_blacklists','more than mblen whitelists')
            lists(j)% wlist(k)% p(i) = wl
          end do
        endif
999     continue ! End Of File
        write (6,*) nl , ' blacklist lines read.'
        write (6,*) nlw, ' whitelist lines read.'
        if (iuu /= nnml) then
          close (iuu)
        else
          backspace (iuu)
        endif
      end do  ! read namelist group
      !---------
      ! printout
      !---------
      write (6,'(a)')
      if (nlists==0) then
        write (6,'(a)') '  no blacklists read'
      else
        write (6,'(a,i3)') '  type blacklist-entries'
        do i=1,nlists
          write (6,"(2x,a8,' (',i2,')  ',i5)") obstyp(lists(i)%type)% name,  &
                                                      lists(i)%type,         &
                                                      lists(i)%n
        end do
      endif
      write (6,'(a)')
      !--------
      ! cleanup
      !--------
      call return_unit_number (iu)
    endif
    !---------------------
    ! broadcast blacklists
    !---------------------
    call bcast_blacklist (lists, nlists, dace% pio)
    call p_bcast         (verbose,       dace% pio)
  end subroutine read_blacklists
!------------------------------------------------------------------------------
  subroutine set_mask (bl)
  type (t_black), intent (inout) :: bl
  !------------------------------------------
  ! derive integer representation of wildcard
  !------------------------------------------
    integer(i8) :: m     ! mask variable
    integer     :: i     ! index variable
    bl% id   = transfer  (bl% statid, bl% id)
    bl%len   = len_trim  (bl% statid)
    if(bigend) call flip (bl% id)
    bl% mask = not(  0_i8)
    m        = not(255_i8)
    do i=1,8
      if(bl% statid(i:i)=='.') then
        bl% id   = iand (bl% id,   m)
        bl% mask = iand (bl% mask, m)
      endif
      m = ior( ishft(m, 8), 255_i8)
    end do
  end subroutine set_mask
!------------------------------------------------------------------------------
  subroutine clean_blacklist
  !----------------------
  ! deallocate blacklists
  !----------------------
    integer :: i, j
    if (associated (lists)) then
      do i=1,nlists
        if (associated (lists(i)%  list)) deallocate (lists(i)% list)
        if (associated (lists(i)% wlist)) then
          do j = 1, lists(i)% nw
            if (associated (lists(i)% wlist(j)% p)) &
                deallocate (lists(i)% wlist(j)% p)
          end do
          deallocate       (lists(i)% wlist)
        endif
      end do
      deallocate (lists)
    endif
    nlists = 0
  end subroutine clean_blacklist
!------------------------------------------------------------------------------
  subroutine empty_blacklist
  !-----------------------
  ! create empty blacklist
  !-----------------------
    integer :: i
    call clean_blacklist
    allocate (lists(mty))
    nlists = mty
    bigend = .not.little()
    do i = 1, mty
      allocate (lists(i)%  list (0))
      allocate (lists(i)% wlist (0))
      lists(i)% type = i
      lists(i)% n    = 0
      lists(i)% nw   = 0
    end do
  end subroutine empty_blacklist
!BOP---------------------------------------------------------------------------
! !IROUTINE: black_entry
!
! !DESCRIPTION:
!
! Access to the blacklisting information is provided by function {\bf
! black\_entry}.  The argument {\bf found} returns the number of
! matching blacklist entries.  If found>0, the routine returns the
! blacklist entries {\bf entry} matching the given station
! identifier {\bf statid} and observation type {\bf typ}.
! If found<0, the size of array entry was insufficient.
!
! !INTERFACE:
!
  function black_entry (statid, typ, code, entry) result (found)
  integer                    :: found     ! number of matching entries
  character(len=8)           :: statid    ! station identifier
  integer       ,intent(in)  :: typ       ! observation type
  integer       ,intent(in)  :: code      ! observation code type
  type(t_black) ,intent(out) :: entry (:) ! blacklist entry

!EOP
    !----------------
    ! local variables
    !----------------
    integer                 :: i,j,len,mb,n
    integer(i8)             :: id
    type (t_black) ,pointer :: e
#if defined (__SX__)
    logical,    allocatable :: mask(:)     ! Auxiliary mask
    integer,    allocatable :: idx(:)      ! Auxiliary index vector
    integer                 :: n_idx, k
#endif
    !----------------------
    ! executable statements
    !----------------------
    mb = size (entry)
    found = 0
    len = len_trim (statid)
    id  = transfer (statid, id)
    if(bigend) call flip (id)
    entry% pg (1) = 0
    entry% pw (1) = 0
    entry% pt (1) = 0
    entry% pd (1) = 0
    entry% pg (2) = 999999
    entry% pw (2) = 999999
    entry% pt (2) = 999999
    entry% pd (2) = 999999
    entry% zg (1) = 0
    entry% zw (1) = 0
    entry% zt (1) = 0
    entry% zd (1) = 0
    entry% zg (2) = 999999
    entry% zw (2) = 999999
    entry% zt (2) = 999999
    entry% zd (2) = 999999
    !++++++++++++++++++++++++++++++++++++++++++++
    ! default initialisation does not work on SX6
    !++++++++++++++++++++++++++++++++++++++++++++
!   if (entry%lg .or. entry%lw .or. entry%lt .or. entry%ld) &
!     call finish('black_entry','default initialisation error')
    entry% lg  = .false.
    entry% lw  = .false.
    entry% lt  = .false.
    entry% ld  = .false.
    entry% mwl = .false.
    do i=1,nlists
      if (lists(i)% type == typ) then

        !-----------------
        ! handle blacklist
        !-----------------
        n = lists(i)% n
        if (n > 0) then
#if defined (__SX__)
          !----------------------------------------------
          ! Vectorized scan using an auxiliary index list
          !----------------------------------------------
          allocate (mask(n))
          mask = (lists(i)% list(:n)% len == len                               &
            .and. lists(i)% list(:n)% id  == iand (lists(i)% list(:n)% mask, id))
          n_idx = count (mask)
          if (n_idx > 0) then
            allocate (idx(n_idx))
            idx = pack ( (/ (k, k=1,n) /), mask)
            !------------------------------------------
            ! Now we walk only through the list of hits
            !------------------------------------------
            do k = 1, n_idx
              j = idx(k)
#else
            !---------------------------
            ! Slow scalar blacklist scan
            !---------------------------
          if (.true.) then
            do j=1,lists(i)% n
#endif
              e => lists(i)% list(j)
              if     (e% len  == len               &
                .and. e% id   == iand (e% mask, id)) then
                found = found + 1
                if (found <= mb) entry (found) = e
              endif
            end do ! list of hits
          endif    ! n_idx > 0
        endif      ! n     > 0
        !-----------------
        ! handle whitelist
        !-----------------
        n = lists(i)% nw
        if (n > 0) then ! obstype   whitelist present
          do j = 1, n   ! codetypes whitelists
            if   (lists(i)% wlist(j)% codetype == code) then
              if (lists(i)% wlist(j)% n        >  0   ) then
                if (all (lists(i)% wlist(j)% p% id /= id)) then
                  found = max (found, 1)
                  entry% mwl = .true.
                endif   ! whitelist entry missing
              endif     ! whitelist present
              exit
            endif       ! codetype  matches
          end do        ! codetypes whitelists
        endif           ! obstype   whitelist present
        !----------------------------------------------
        ! check for sufficient size of output parameter
        !----------------------------------------------
        if (found > mb) then
          call message ('black_entry','increase number of entries mb !')
          found = - found
        endif
        exit
      endif  ! obstype matches
    end do   ! nlists (obstypes)
  end function black_entry
!==============================================================================
  subroutine bcast_blacklist (blacklist, n, source)
  type (t_blacklists) ,pointer       :: blacklist(:) ! blacklists to broadcast
  integer             ,intent(inout) :: n            ! number of blacklists
  integer             ,intent(in)    :: source       ! source processor element
    !----------------
    ! local variables
    !----------------
    integer                 :: s_black      ! size of data type t_black
    integer                 :: s_blacklists ! size of data type t_blacklists
    integer                 :: s_white      ! size of data type t_white
    type (t_black)          :: d_black      ! dummy data type
    type (t_blacklists)     :: d_blacklists ! dummy data type
    type (t_white)          :: d_white      ! dummy data type
    integer                 :: i, j         ! blacklist index variable
    type (t_white) ,pointer :: p (:)

#if (defined (__GFORTRAN__) && (__GNUC__ >= 10)) || defined (NAGFOR)
    !----------------------------------------------
    ! include interfaces for external bcast routine
    !----------------------------------------------
    interface
       subroutine p_bcast_derivedtype (buffer, count, source, comm)
         type(*) ,intent(inout)     :: buffer(*)    ! variable to bcast
         integer ,intent(in)        :: count        ! len(byte) of variable
         integer ,intent(in)        :: source       ! source processor index
         integer ,intent(in)        :: comm         ! communicator
       end subroutine p_bcast_derivedtype
    end interface
#endif

    !--------------------------
    ! determine data type sizes
    !--------------------------
    s_blacklists = size (transfer(d_blacklists ,(/' '/)))
    s_black      = size (transfer(d_black      ,(/' '/)))
    s_white      = size (transfer(d_white      ,(/' '/)))
    !--------------------
    ! broadcast container
    !--------------------
    call p_bcast (n, source)
    if (dace% pe /= source) allocate (blacklist (n))
    call p_bcast_derivedtype (blacklist, n*s_blacklists, &
                              source, (dace% comm)       )
    !----------------
    ! send components
    !----------------
    do i = 1, n
      if (dace% pe /= source) then
        allocate (blacklist(i)%  list (blacklist(i)% n ))
        allocate (blacklist(i)% wlist (blacklist(i)% nw))
      endif
      call p_bcast_derivedtype (blacklist(i)% list,        &
                                blacklist(i)% n * s_black, &
                                source, (dace% comm)       )
      call p_bcast (blacklist(i)% wlist(1:blacklist(i)% nw)% n, &
                                  source, dace% comm            )
      call p_bcast (blacklist(i)% wlist(1:blacklist(i)% nw)% codetype, &
                                  source, dace% comm                   )
      do j = 1, blacklist(i)% nw
        if (blacklist(i)% wlist(j)% n > 0) then
          if (dace% pe /= source) then
            allocate (blacklist(i)% wlist(j)% p (blacklist(i)% wlist(j)% n))
          else
            p => blacklist(i)% wlist(j)% p
            allocate (blacklist(i)% wlist(j)% p (blacklist(i)% wlist(j)% n))
            blacklist(i)% wlist(j)% p = p (:blacklist(i)% wlist(j)% n)
            deallocate (p)
            write (6,*) blacklist(i)% wlist(j)% n,           &
                        ' whitelist entries for type/code:', &
                        blacklist(i)% type,                  &
                        blacklist(i)% wlist(j)% codetype
          endif
          call p_bcast_derivedtype (blacklist(i)% wlist(j)% p,          &
                                    blacklist(i)% wlist(j)% n * s_white,&
                                    source, (dace% comm)                )
        endif
      end do
    end do
  end subroutine bcast_blacklist
!==============================================================================
end module mo_blacklist
