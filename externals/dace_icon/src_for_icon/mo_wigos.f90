!
!+ Handling of WIGOS station identifiers
!
MODULE mo_wigos
!
! Description:
!   This module provides support for WIGOS Station Identifiers (WSIs):
!   data type 't_wigos',
!   conversion routines to and from text representations,
!   generation of "fake station IDs" (ECMWF-type or other variants),
!   and an interface to message digest algorithm MD5
!   for the derivation of station ID hashes.
!
! Current Maintainer: DWD, Harald Anlauf
!    phone: +49 69 8062 4941
!    fax:   +49 69 8062 3721
!    email: harald.anlauf@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V2_11        2022-04-14 Harald Anlauf
!
! Code Description:
! Language: Fortran.
! Software Standards:
!
! Authors:
! Harald Anlauf  DWD   2022       original code
!------------------------------------------------------------------------------
  !-------------
  ! Modules used
  !-------------
  use mo_kind,       only: i8              ! 8-byte integer kind
  use mo_exception,  only: finish          ! exit on error
  use mo_namelist,   only: position_nml,  &! position namelist
                           nnml,          &! namelist Fortran unit number
                           POSITIONED      ! ok    code from position_nml
  use mo_mpi_dace,   only: dace,          &! MPI group info
                           p_bcast         ! broadcast routine
  implicit none
  !----------------
  ! Public entities
  !----------------
  private
  public :: t_wsi           ! Derived type holding WIGOS station identifier
  public :: wsi_encode      ! Derive t_spot%statid, t_spot%stathash from t_wsi
  public :: wsi_decode      ! Derive t_wsi from given string
  public :: wsi_to_text     ! Convert t_wsi to text format
  public :: wsi_mode        ! WSI handling mode
  public :: wsi_verbose     ! Verbosity level
  public :: wsi_from_statid ! Derive WSI from short station name as fallback?
  public :: wsi_prefix      ! WSI prefix character
  public :: wsi_digesttype  ! WSI message digest type
  public :: read_wigos_nml  ! Read WIGOS namelist
  public :: operator (==)   ! Same WSI?
  public :: dace_md5
  !============================================================================
  !------------------------------------
  ! Derived types and module variables:
  !------------------------------------
  type t_wsi
    logical       :: valid = .false. ! WSI is valid
    integer       :: wigis = -1      ! WIGOS identifier series
    integer       :: wigii = -1      ! WIGOS issuer of identifier
    integer       :: wigin = -1      ! WIGOS issue number
    character(16) :: wigli = ""      ! WIGOS local identifier
  end type t_wsi
  !============================================================================
  ! Namelist variables
  !-------------------
  integer      :: wsi_mode        = 0       ! 0:off, 1:ECMWF, >=2:DACE scheme
  integer      :: wsi_verbose     = 0       ! verbosity level (0,1,2)
  integer      :: wsi_digestlen   = -1      ! automatic (scheme-dep.) or value
  integer      :: wsi_encbase     = 16      ! encoding base: 16/32/64
  logical      :: wsi_from_statid = .false. ! derive from short station name?
  protected    :: wsi_mode, wsi_verbose
  !-------------------
  ! Reserved variables
  !-------------------
  character(1) :: wsi_prefix      = "_"     ! (fixed)
  character(3) :: wsi_digesttype  = "md5"   ! (fixed)
  protected    :: wsi_prefix
  protected    :: wsi_digesttype
  !----------------------------------------------------------------------------
  namelist /wigos/ wsi_mode,      wsi_verbose, &
                   wsi_digestlen, wsi_encbase, &
                   wsi_from_statid
  !----------------------------------------------------------------------------
  interface operator (==)
     module procedure equal_wsi
  end interface operator (==)
  !----------------------------------------------------------------------------
  interface
     subroutine dace_md5 (text, len, base,      &
                          digest, string, status) bind(c,name="dace_md5")
       use, intrinsic :: iso_c_binding, only : c_char, c_int
       character(kind=c_char), intent(in)  :: text(*)
       integer(c_int), value,  intent(in)  :: len
       integer(c_int), value,  intent(in)  :: base
       character(kind=c_char), intent(out) :: digest(16)
       character(kind=c_char), intent(out) :: string(32)
       integer(c_int),         intent(out) :: status
     end subroutine dace_md5
  end interface
  !============================================================================
contains
  !============================================================================
  subroutine wsi_to_text (wsi, s)
    !-----------------------------
    ! Convert t_wsi to text format
    !-----------------------------
    type(t_wsi),  intent(in)  :: wsi
    character(*), intent(out) :: s
    character(32) :: t
    if (wsi% valid) then
       write (t,'(i0,"-",i0,"-",i0,"-",a)') &
            wsi% wigis, wsi% wigii, wsi% wigin, wsi% wigli
       s = t
    else
       s = repeat ("*", len (s))
    end if
  end subroutine wsi_to_text
  !----------------------------------------------------------------------------
  subroutine wsi_encode (wsi, statid, stathash)
    !-------------------------------------------------
    ! Derive t_spot%statid, t_spot%stathash from t_wsi
    !-------------------------------------------------
    type(t_wsi),  intent(in)  :: wsi
    character(*), intent(out) :: statid
    integer(i8),  intent(out) :: stathash
    character(32) :: t
    character(16) :: fsi
    character(20) :: ecid
    character(4)  :: echash
    character(16) :: digest
    character(32) :: string
    integer       :: slen, status, encbase

    if (.not. wsi% valid) goto 999
    select case (wsi_mode)
    case (1)   ! ECMWF
       write (ecid,'(i0,a)') wsi% wigin, wsi% wigli
       slen = len_trim (ecid)
       call dace_md5 (ecid, slen, 16, digest, string, status)
       if (status /= 0) goto 999
       echash  = string(29:32)
       encbase = 16
    case (2,3) ! DACE
       encbase = wsi_encbase
    end select

    call wsi_to_text (wsi, t)
    slen = len_trim (t)
    call dace_md5 (t, slen, encbase, digest, string, status)
    if (status /= 0) goto 999

    select case (wsi_mode)
    case (1)
       write (fsi,'(a,i3.3,a)') wsi_prefix, wsi% wigii, echash
    case (2)
       write (fsi,'(a,i3.3,a)') wsi_prefix, wsi% wigii, string(1:wsi_digestlen)
    case (3)
       fsi = wsi_prefix // string(1:wsi_digestlen)
    end select
    statid   = fsi
    stathash = ieor (transfer (digest(1: 8), stathash), &
                     transfer (digest(9:16), stathash)  )
    return

999 continue
    statid   = repeat ('*', len (statid))
    stathash = 0
    stop 999
  end subroutine wsi_encode
  !----------------------------------------------------------------------------
  subroutine wsi_decode (str, wsi, status)
    !-------------------------------
    ! Derive t_wsi from given string
    !-------------------------------
    character(*), intent(in)  :: str
    type(t_wsi),  intent(out) :: wsi
    integer,      intent(out) :: status
    integer :: slen, i2, i3, lidlen
    status = 1
    slen = len_trim (str)
    ! min. length id: 0-1-0-X
    if (slen     <  7   ) return
    ! Verify presence of identifier series
    if (str(1:2) /= "0-") return
    wsi% wigis = 0

    ! Check and decode issuer of identifier (length: i2-1)
    i2 = index (str(3   :slen), "-")
    if (i2 <= 1 .or. i2 > 6) return
    wsi% wigii = atoi (str(3:1+i2))
    if (wsi% wigii < 0 .or. wsi% wigii > 65535) return

    ! Check and decode issue number (length: i3-1)
    i3 = index (str(3+i2:slen), "-")
    if (i3 <= 1 .or. i3 > 6) return
    wsi% wigin = atoi (str(3+i2:1+i2+i3))
    if (wsi% wigin < 0 .or. wsi% wigin > 65535) return

    ! Valid length of local identifier?
    lidlen = slen - (2+i2+i3)
    if (lidlen <= 0 .or. lidlen > 16) return
    wsi% wigli = str(3+i2+i3:2+i2+i3+lidlen)

    ! Gross check: local identifier must not start with blank or "-"
    select case (wsi% wigli(1:1))
    case (' ', '-')
       return
    end select

    wsi% valid = .true.
    status = 0
!   print *, "wsi:", wsi% wigis, wsi% wigii, wsi% wigin, "|", wsi% wigli, "|"

  contains
    !----------------------------
    ! ASCII to integer conversion
    ! Returns -1 on error.
    !----------------------------
    integer function atoi (s) result (res)
      character(*), intent(in) :: s
      integer :: i, k
      res = 0
      do i = 1, len (s)
         k = ichar (s(i:i)) - ichar ('0')
         if (k < 0 .or. k > 9) then
            res = -1
            return
         end if
         res = 10*res + k
      end do
    end function atoi
  end subroutine wsi_decode
  !----------------------------------------------------------------------------
  subroutine read_wigos_nml ()
    !----------------------
    ! read namelist /WIGOS/
    !----------------------
    integer :: ierr
    logical :: nml_read = .false.   ! namelist not yet read

    !----------------------------
    ! call this routine only once
    !----------------------------
    if (nml_read) return
    nml_read = .true.

    !-------------
    ! set defaults
    !-------------
    wsi_mode        = 0       ! 0:off, 1:ECMWF, >=2:DACE scheme
    wsi_verbose     = 0       ! verbosity level (0,1,2)
    wsi_digestlen   = -1      ! automatic (scheme-dep.) or value
    wsi_encbase     = 16      ! encoding base: 16/32/64
    wsi_from_statid = .false. ! derive from short station name?

    !--------------
    ! read namelist
    !--------------
    if (dace% lpio) then
      write(6,'()')
      write(6,'(a)') repeat('-',79)
      write(6,'()')
      call position_nml ('WIGOS', status=ierr)
      select case (ierr)
      case (POSITIONED)
        write(6,*) 'Namelist /WIGOS/:'
#if defined(__ibm__)
        read (nnml ,nml=WIGOS, iostat=ierr)
        if (ierr/=0) call finish ('read_wigos_nml',          &
                                  'ERROR in namelist /WIGOS/')
#else
        read (nnml ,nml=WIGOS)
#endif
      case default
        write(6,*) 'Namelist /WIGOS/ not present, defaults used'
      end select
      !---------------------------
      ! Use mode-specific defaults
      !---------------------------
      select case (wsi_mode)
      case (1)
         wsi_digestlen = 4
         wsi_encbase   = 16
      case (2)
         if (wsi_digestlen < 0) wsi_digestlen = 4
      case (3)
         if (wsi_digestlen < 0) wsi_digestlen = 7
      case default
         wsi_from_statid = .false.
      end select
      !---------
      ! printout
      !---------
      write(6,'()')
      write(6,'(a,i6)') ' wsi_mode        =',wsi_mode
      write(6,'(a,i6)') ' wsi_verbose     =',wsi_verbose
      write(6,'(a,i6)') ' wsi_digestlen   =',wsi_digestlen
      write(6,'(a,i6)') ' wsi_encbase     =',wsi_encbase
      write(6,'(a,l6)') ' wsi_from_statid =',wsi_from_statid
      !-------
      ! Checks
      !-------
      if (wsi_mode    <  0  .or.   wsi_mode   >  3)             &
           call finish ("read_wigos_nml","wsi_mode must be 0..3")
      if (wsi_encbase /= 16 .and. wsi_encbase /= 32 .and. wsi_encbase /= 64) &
           call finish ("read_wigos_nml","wsi_encbase must be 16 or 32 or 64")
      select case (wsi_mode)
      case (2)
         if (wsi_digestlen < 4 .or. wsi_digestlen > 6) &
              call finish ("read_wigos_nml","wsi_digestlen out of range")
      case (3)
         if (wsi_digestlen < 6 .or. wsi_digestlen > 9) &
              call finish ("read_wigos_nml","wsi_digestlen out of range")
      end select
    end if
    call p_bcast (wsi_mode       , dace% pio)
    call p_bcast (wsi_verbose    , dace% pio)
    call p_bcast (wsi_digestlen  , dace% pio)
    call p_bcast (wsi_encbase    , dace% pio)
    call p_bcast (wsi_from_statid, dace% pio)

  end subroutine read_wigos_nml
  !----------------------------------------------------------------------------
  elemental function equal_wsi (a, b) result (res)
    type(t_wsi), intent(in) :: a, b
    logical                 :: res
    res = a% valid             .and. &
          b% valid             .and. &
          a% wigis == b% wigis .and. &
          a% wigii == b% wigii .and. &
          a% wigin == b% wigin .and. &
          a% wigli == b% wigli
  end function equal_wsi
  !----------------------------------------------------------------------------
end module mo_wigos
