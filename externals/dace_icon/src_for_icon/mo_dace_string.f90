!
!+ Character string conversion utilities
!
MODULE mo_dace_string
!
! Description:
!   This module holds character string conversion utilities.
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
! V1_4         2009/03/26 Andreas Rhodin
!  New function char5: convert integer to character(len=5)
! V1_5         2009/05/25 Harald Anlauf
!  tolower, toupper, underscore: optimize for NEC SX where LEN_TRIM is slow
! V1_8         2009/12/09 Andreas Rhodin
!  new function char4  ! Conversion: integer -> char(len=4)
! V1_9         2010/04/20 Harald Anlauf
!  unquote_string: remove quotes around string
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Andreas Rhodin
!  new function char1 (INTEGER -> CHARACTER(LEN=1))
! V1_19        2012-04-16 Andreas Rhodin
!  new routine eval_string: add or remove words from a string
! V1_23        2013-03-26 Andreas Rhodin
!  new function 'replace' (replace character in string)
! V1_26        2013/06/27 Harald Anlauf
!  Changes for GRIB2/GRIB_API/ICON
! V1_42        2015-06-08 Harald Anlauf
!  byte2hex: extend to array of char
! V1_47        2016-06-06 Harald Anlauf
!  new conversion 'charn'
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Andreas Rhodin  DWD  2001-2008
! Harald Anlauf   DWD  2008
!------------------------------------------------------------------------------
  USE mo_exception, ONLY: finish
  IMPLICIT NONE
  PRIVATE
  !----------------
  ! Public entities
  !----------------
  PUBLIC :: tolower        ! Conversion   : 'ABCXYZ' -> 'abcxyz'
  PUBLIC :: toupper        ! Conversion   : 'abcxyz' -> 'ABCXYZ'
  PUBLIC :: underscore     ! Conversion   : 'a b c ' -> 'a_b_c'
  PUBLIC :: replace        ! Conversion   : replace characters by others
  PUBLIC :: strip          ! Conversion   : strip a string, i.e. remove characters ar start/end
  PUBLIC :: char1          ! Conversion   : INTEGER  -> CHAR (LEN=1)
  PUBLIC :: char2          ! Conversion   : INTEGER  -> CHAR (LEN=2)
  PUBLIC :: char3          ! Conversion   : INTEGER  -> CHAR (LEN=3)
  PUBLIC :: char4          ! Conversion   : INTEGER  -> CHAR (LEN=4)
  PUBLIC :: char5          ! Conversion   : INTEGER  -> CHAR (LEN=5)
  PUBLIC :: charn          ! Conversion   : INTEGER  -> CHAR (LEN=n)
  PUBLIC :: byte2hex       ! Convert CHAR to 'hexadecimal' CHAR(LEN=2)
  PUBLIC :: separator      ! Format string: (/"-----...-----"/)
  PUBLIC :: split          ! Put words of a string into a character-array
  PUBLIC :: split2         ! Put words of a string into a character-array
  PUBLIC :: intstr2array   ! Convert a string with an integer list to an integer/logical array
  PUBLIC :: concat         ! concatenate array of chars to string
  PUBLIC :: unquote_string ! Remove quotes around string
  PUBLIC :: eval_string    ! add/remove words from string
  PUBLIC :: decode_uuid    ! Decode uuid (hexadecimal string -> char array)
  !-----------------
  ! module variables
  !-----------------
  character(len=*) ,parameter :: separator = '("'//repeat('-',79)//'")'
  !-----------------------------------------------------------------
  ! number of words and characters handled by subroutine eval_string
  !-----------------------------------------------------------------
  integer ,parameter :: nw = 128   ! number of words
  integer ,parameter :: nc =  16   ! number of characters
  !-----------
  ! interfaces
  !-----------
  interface byte2hex
     module procedure byte2hex      ! Convert single char and
     module procedure byte2hex_1d   ! .. 1d-array of char to hex strings
  end interface byte2hex
!==============================================================================
CONTAINS
!==============================================================================
  pure FUNCTION tolower (upper)
    !-----------------------------------
    ! Conversion: Uppercase -> Lowercase
    !-----------------------------------
    CHARACTER(LEN=*)              ,INTENT(in) :: upper
    CHARACTER(LEN=LEN_TRIM(upper))            :: tolower

    INTEGER            :: i
    INTEGER ,PARAMETER :: idel = ICHAR('a')-ICHAR('A')

    DO i=1,LEN (tolower)
      IF (ICHAR(upper(i:i)) >= ICHAR('A') .AND. &
          ICHAR(upper(i:i)) <= ICHAR('Z')) THEN
        tolower(i:i) = CHAR( ICHAR(upper(i:i)) + idel )
      ELSE
        tolower(i:i) = upper(i:i)
      END IF
    END DO

  END FUNCTION tolower
!------------------------------------------------------------------------------
  pure FUNCTION toupper (lower)
    !-----------------------------------
    ! Conversion: Lowercase -> Uppercase
    !-----------------------------------
    CHARACTER(LEN=*)              ,INTENT(in) :: lower
    CHARACTER(LEN=LEN_TRIM(lower))            :: toupper

    INTEGER            :: i
    INTEGER ,PARAMETER :: idel = ICHAR('a')-ICHAR('A')

    DO i=1,LEN (toupper)
      IF (ICHAR(lower(i:i)) >= ICHAR('a') .AND. &
          ICHAR(lower(i:i)) <= ICHAR('z')) THEN
        toupper(i:i) = CHAR( ICHAR(lower(i:i)) - idel )
      ELSE
        toupper(i:i) = lower(i:i)
      END IF
    END DO

  END FUNCTION toupper
!------------------------------------------------------------------------------
  pure FUNCTION replace (string, c, r)
    !---------------------------------------
    ! Conversion: replace character c with r
    !---------------------------------------
    CHARACTER(LEN=*)    ,INTENT(in) :: string
    CHARACTER           ,INTENT(in) :: c
    CHARACTER           ,INTENT(in) :: r
    CHARACTER(LEN=LEN_TRIM(string)) :: replace

    INTEGER            :: i

    DO i=1,LEN (replace)
      IF (string(i:i) == c) then
        replace(i:i) = r
      ELSE
        replace(i:i) = string(i:i)
      END IF
    END DO

  END FUNCTION replace
!------------------------------------------------------------------------------
  pure function strip(string, set, se)
    character(len=*), intent(in)           :: string  ! string to be modified
    character(len=*), intent(in)           :: set     ! characters to stripped from string
    integer,          intent(in), optional :: se      ! flags:
                                                      ! 1: strip start of string
                                                      ! 2: strip end   of string
                                                      ! 3: strip start and end of string
    character(len=len_trim(string))        :: strip

    integer :: i, se_loc

    strip = string
    if (set == '') return
    
    if (present(se)) then
      se_loc = se
    else
      se_loc = 3
    end if

    if (iand(se_loc, 1) > 0) then
      do i = 1, len_trim(strip)
        if (scan(strip, set) == i) then
          strip(i:i) = ' '
        else
          exit
        end if
      end do
      strip = adjustl(strip)
    end if
    if (iand(se_loc, 2) > 0) then
      do i = len_trim(strip), 1, -1
        if (scan(strip, set, back=.true.) == i) then
          strip(i:i) = ' '
        else
          exit
        end if
      end do
      strip = trim(strip)
    end if

  end function strip
!------------------------------------------------------------------------------
  pure FUNCTION underscore (blanks)
    !----------------------------------------------------------
    ! Conversion: string with blanks -> string with underscores
    !----------------------------------------------------------
    CHARACTER(LEN=*)    ,INTENT(in) :: blanks
    CHARACTER(LEN=LEN_TRIM(blanks)) :: underscore

    INTEGER            :: i

    DO i=1,LEN (underscore)
      IF (blanks(i:i) == ' ') then
        underscore(i:i) = '_'
      ELSE
        underscore(i:i) = blanks(i:i)
      END IF
    END DO

  END FUNCTION underscore
!------------------------------------------------------------------------------
  pure FUNCTION char2 (i, zero)
    !----------------------------------------
    ! Conversion: INTEGER -> CHARACTER(LEN=2)
    !----------------------------------------
    CHARACTER(LEN=2)                       :: char2 ! result
    INTEGER          ,INTENT(in)           :: i     ! argument
    CHARACTER        ,INTENT(in) ,OPTIONAL :: zero  ! padding instead of '0'

    INTEGER ,PARAMETER :: i0 = ICHAR ('0')

    IF (i>99 .OR. i<0) THEN
      char2 = '**'
    ELSE
      char2(1:1) = CHAR(    i/10  + i0)
      char2(2:2) = CHAR(MOD(i,10) + i0)
    ENDIF

    IF(PRESENT(zero)) THEN
      IF(char2(1:1) /= '0') return; char2(1:1) = zero
      IF(char2(2:2) /= '0') return; char2(2:2) = zero
    ENDIF
  END FUNCTION char2
!------------------------------------------------------------------------------
  pure FUNCTION char1 (i)
    !----------------------------------------
    ! Conversion: INTEGER -> CHARACTER(LEN=1)
    !----------------------------------------
    CHARACTER(LEN=1)                       :: char1 ! result
    INTEGER          ,INTENT(in)           :: i     ! argument

    INTEGER ,PARAMETER :: i0 = ICHAR ('0')

    IF (i>9 .OR. i<0) THEN
      char1 = '*'
    ELSE
      char1 = CHAR(i + i0)
    ENDIF

  END FUNCTION char1
!------------------------------------------------------------------------------
  pure FUNCTION char3 (i, zero)
    !----------------------------------------
    ! Conversion: INTEGER -> CHARACTER(LEN=3)
    !----------------------------------------
    CHARACTER(LEN=3)                       :: char3 ! result
    INTEGER          ,INTENT(in)           :: i     ! argument
    CHARACTER        ,INTENT(in) ,OPTIONAL :: zero  ! padding instead of '0'

    INTEGER ,PARAMETER :: i0 = ICHAR ('0')

    IF (i>999 .OR. i<0) THEN
      char3 = '***'
    ELSE
      char3(1:1) = CHAR(    i/100     + i0)
      char3(2:2) = CHAR(MOD(i/10 ,10) + i0)
      char3(3:3) = CHAR(MOD(i    ,10) + i0)
    ENDIF

    IF(PRESENT(zero)) THEN
      IF(char3(1:1) /= '0') return; char3(1:1) = zero
      IF(char3(2:2) /= '0') return; char3(2:2) = zero
      IF(char3(3:3) /= '0') return; char3(3:3) = zero
    ENDIF
  END FUNCTION char3
!------------------------------------------------------------------------------
  pure FUNCTION char4 (i, zero)
    !----------------------------------------
    ! Conversion: INTEGER -> CHARACTER(LEN=5)
    !----------------------------------------
    CHARACTER(LEN=4)                       :: char4 ! result
    INTEGER          ,INTENT(in)           :: i     ! argument
    CHARACTER        ,INTENT(in) ,OPTIONAL :: zero  ! padding instead of '0'

    INTEGER ,PARAMETER :: i0 = ICHAR ('0')

    IF (i>9999 .OR. i<0) THEN
      char4 = '****'
    ELSE
      char4(1:1) = CHAR(    i/1000     + i0)
      char4(2:2) = CHAR(MOD(i/100 ,10) + i0)
      char4(3:3) = CHAR(MOD(i/10  ,10) + i0)
      char4(4:4) = CHAR(MOD(i     ,10) + i0)
    ENDIF

    IF(PRESENT(zero)) THEN
      IF(char4(1:1) /= '0') return; char4(1:1) = zero
      IF(char4(2:2) /= '0') return; char4(2:2) = zero
      IF(char4(3:3) /= '0') return; char4(3:3) = zero
      IF(char4(4:4) /= '0') return; char4(4:4) = zero
    ENDIF
  END FUNCTION char4
!------------------------------------------------------------------------------
  pure FUNCTION char5 (i, zero)
    !----------------------------------------
    ! Conversion: INTEGER -> CHARACTER(LEN=5)
    !----------------------------------------
    CHARACTER(LEN=5)                       :: char5 ! result
    INTEGER          ,INTENT(in)           :: i     ! argument
    CHARACTER        ,INTENT(in) ,OPTIONAL :: zero  ! padding instead of '0'

    INTEGER ,PARAMETER :: i0 = ICHAR ('0')

    IF (i>99999 .OR. i<0) THEN
      char5 = '*****'
    ELSE
      char5(1:1) = CHAR(    i/10000     + i0)
      char5(2:2) = CHAR(MOD(i/1000 ,10) + i0)
      char5(3:3) = CHAR(MOD(i/100  ,10) + i0)
      char5(4:4) = CHAR(MOD(i/10   ,10) + i0)
      char5(5:5) = CHAR(MOD(i      ,10) + i0)
    ENDIF

    IF(PRESENT(zero)) THEN
      IF(char5(1:1) /= '0') return; char5(1:1) = zero
      IF(char5(2:2) /= '0') return; char5(2:2) = zero
      IF(char5(3:3) /= '0') return; char5(3:3) = zero
      IF(char5(4:4) /= '0') return; char5(4:4) = zero
      IF(char5(5:5) /= '0') return; char5(5:5) = zero
    ENDIF
  END FUNCTION char5
!------------------------------------------------------------------------------
  pure FUNCTION charn (i, n, zero)
    !----------------------------------------
    ! Conversion: INTEGER -> CHARACTER(LEN=n)
    !----------------------------------------
    INTEGER          ,INTENT(in)           :: i     ! argument
    INTEGER          ,INTENT(in)           :: n     ! string length
    CHARACTER        ,INTENT(in) ,OPTIONAL :: zero  ! padding instead of '0'
    CHARACTER(LEN=n)                       :: charn ! result

    INTEGER ,PARAMETER :: i0 = ICHAR ('0')
    integer            :: k, m

    IF (i >= 0 .and. i < 10**n) THEN
       m = i
       do k = n, 1, -1
          charn(k:k) = CHAR (MOD (m, 10) + i0)
          m          = m / 10
       end do
    else
       charn = repeat ("*", n)
    ENDIF

    IF (PRESENT (zero)) THEN
       do k = 1, n-1
          if (charn(k:k) /= '0') return
          charn(k:k) = zero
       end do
    ENDIF
  END FUNCTION charn
!------------------------------------------------------------------------------
  pure function byte2hex (c) result (hex)
    !--------------------------------------------
    ! Convert single byte to 'hexadecimal' string
    !--------------------------------------------
    character(len=1), intent(in) :: c
    character(len=2)             :: hex
    integer :: x
    x = iachar (c)
    hex(1:1) = nibble (      x / 16)
    hex(2:2) = nibble (iand (x,  15))
  contains
    pure function nibble (x)
      !----------------------------------------------
      ! Convert half-byte ('nibble') to 'hexadecimal'
      !----------------------------------------------
      integer, intent(in) :: x
      character           :: nibble
      select case (x)
      case (0:9)
         nibble = achar (iachar ('0') + x)
      case default
         nibble = achar (iachar ('a') - 10 + x)
      end select
    end function nibble
  end function byte2hex
!------------------------------------------------------------------------------
  pure function byte2hex_1d (c) result (hex)
    !-------------------------------------------
    ! Convert byte array to 'hexadecimal' string
    !-------------------------------------------
    character(len=1)          ,intent(in) :: c(:)
    character(len=2*size (c))             :: hex
    integer :: k
    do k = 1, size (c)
       hex(2*k-1:2*k) = byte2hex(c(k))
    end do
  end function byte2hex_1d
!------------------------------------------------------------------------------
  pure subroutine split (array, string, len)
  character (len=*) ,intent(out)          :: array (:)
  character (len=*) ,intent(in)           :: string
  integer           ,intent(out)          :: len
  !------------------------------------------------------------------------
  ! Put each word (seperated by blanks) of STRING into an element of ARRAY.
  ! In LEN return the number of words returned or -1 in case of error
  ! (dimension of ARRAY too small).
  !------------------------------------------------------------------------

    integer            :: l, i, n, i1, i2, i3

    l = len_trim (string)
    n = size     (array)

    array = ''

    i1  = 0                             ! start position in STRING
    i   = 0                             ! counter for number of words
    len = -1                            ! (error) return argument
    do
      i2 = index(string(i1+1:),' ')     ! position of blank
      if (i2 /= 1)  then                ! skip leading blank
        i3 = l                          ! i3 = end of string
        if (i2 > 1) i3 = i1+i2-1        !      or last char before blank
        if (i3-i1>0) then               ! word bounded by i1..i3
          i = i + 1                     ! next word
          if (i > n) return             ! error return, array too small
          array (i) = string (i1+1:i3)  ! store word
        endif
        if (i2 == 0) exit               ! exit loop, no further blank
      endif
      i1 = i1 + i2                      ! skip word
    end do
    len = i                             ! return number of words

  end subroutine split
!------------------------------------------------------------------------------
  subroutine split2 (string, sep, array, n, status)
    character (len=*) ,intent(in)            :: string
    character (len=*) ,intent(in),  optional :: sep
    character (len=*) ,intent(out), optional :: array (:)
    integer           ,intent(out), optional :: n
    integer           ,intent(out), optional :: status
    !------------------------------------------------------------------------
    ! Put each word of STRING into an element of ARRAY.
    ! In N return the number of words
    !------------------------------------------------------------------------

    integer            :: l, i, n_arr, nsep
    integer            :: i1, i2, i3
    character(len=300) :: sep_ = ' '

    if (present(sep)) then
      nsep = min(len(sep),len(sep_))
      sep_ = sep(1:nsep)
    else
      nsep = 1
      sep_ = ' '
    end if

    if (present(array)) then
      n_arr = size(array)
      array = ''
    else
      n_arr = 0
    end if

    if (present(status)) status = 0

    l  = len_trim (string)
    i1 = 0                             ! start position in STRING
    i  = 0                             ! counter for number of words
    do
      i2 = scan(string(i1+1:l),sep_(1:nsep)) ! position of seperator
      if (i2 /= 1)  then                     ! skip leading separators
        if (i2 > 1) then
          i3 = i1+i2-1                       ! last char before seperator ...
        else
          i3 = l                             !    ... or i3 = end of string
        end if
        if (i3-i1>0) then                    ! word bounded by i1..i3
          i = i + 1                          ! next word
          if (i <= n_arr) then
            if (i3-i1 <= len(array(i))) then
              array(i) = string (i1+1:i3)
            else
              if (present(status)) then
                status = 2
                exit
              else
                stop 2
              end if
            end if
          else if (present(array)) then
            if (present(status)) then
              status = 1
              exit
            else
              stop 1
            end if
          end if
        endif
        if (i2 == 0) exit               ! exit loop, no further blank
      endif
      i1 = i1 + i2                      ! skip word
    end do
    if (present(n)) n = i

  end subroutine split2

!------------------------------------------------------------------------------
  subroutine intstr2array(str, sep, n, imin, imax, iarr, p_iarr, larr, p_larr, status, fstr, ldeb, pe)
    character(len=*), intent(in)                     :: str
    character(len=*), intent(in),           optional :: sep
    integer,          intent(out),          optional :: n           ! number of integer values (min. size of iarr)
    integer,          intent(out),          optional :: imin        ! minimum integer number (lower bound of larr)
    integer,          intent(out),          optional :: imax        ! maximum integer number (upper bound of larr)
    integer,          intent(out),          optional :: iarr(:)     ! array of integer values contained in string
    integer,          intent(out), pointer, optional :: p_iarr(:)   ! array of integer values contained in string
    logical,          intent(out),          optional :: larr(:)     ! logical array (imin:imax), .true. if index in contained in string
    logical,          intent(out), pointer, optional :: p_larr(:)   ! logical array (imin:imax), .true. if index in contained in string
    integer,          intent(out),          optional :: status      ! return code
    character(len=*), intent(out),          optional :: fstr        ! substring, that could not be interpreted
    logical,          intent(in),           optional :: ldeb
    integer,          intent(in),           optional :: pe
    !-----------------------------------------------------------------------------------
    ! Converts strings, that contain a list of integers into an integer or logical array
    ! The values might be separated by blanks, commas or semicolons. ":" and "-" are
    ! interpreted as value ranges. For example: "1, 2, 5-11;20 21, 25:27" results in
    ! iarr = (/1, 2, 5, 6, 7, 8, 9, 10, 11, 20, 21, 25, 26, 27/)
    !-----------------------------------------------------------------------------------
    integer,           parameter   :: n_alloc = 100
    character(len=100)             :: msg     = ''
    character(len=20), allocatable :: sub(:)
    integer                        :: n_sub ! number of substrings
    integer,           pointer     :: ip(:) => null()
    integer                        :: n_ip  ! size of ip
    integer                        :: n_i   ! number of integer values
    integer                        :: i1, i2, i, j, stat, ib, i_range
    integer                        :: imax_, imin_
    logical,           pointer     :: lp(:) => null()
    logical                        :: ldeb_
    integer                        :: pe_

    ! set default values
    if (present(n     )) n      =  0
    if (present(imax  )) imax   =  0
    if (present(iarr  )) iarr   =  0
    if (present(p_iarr)) p_iarr => null()
    if (present(larr  )) larr   =  .false.
    if (present(p_larr)) p_larr => null()
    if (present(status)) status =  0
    if (present(fstr  )) fstr   =  ''
    if (present(ldeb)) then
      ldeb_ = ldeb
    else
      ldeb_ = .false.
    end if
    if (present(pe)) then
      pe_ = pe
    else
      pe_ = -1
    end if

    if (ldeb_) then
      write(1000+pe_,*) 'intstr2array str',len_trim(str),trim(str)
      flush 1000+pe_
    end if

    ! split string
    call split2(str, sep=sep, n=n_sub, status=stat)
    if (stat /= 0) then
      call error(5,'split2 failed (1): "'//trim(str)//'"')
      return
    end if
    if (n_sub <= 0) RETURN
    allocate(sub(n_sub))
    call split2(str, sep=sep, array=sub, n=n_sub, status=stat)
    if (stat /= 0) then
      call error(6,'split2 failed (2): "'//trim(str)//'"')
      return
    end if
    if (ldeb_) then
      write(1000+pe_,*) 'intstr2array n_sub',n_sub
      do i = 1, n_sub
        write(1000+pe_,*) 'intstr2array sub',i, trim(sub(i))
        flush 1000+pe_
      end do
    end if

    ! interpret substrings
    n_ip = n_sub
    allocate(ip(n_ip))
    n_i = 0
    do i = 1, n_sub
      i_range = scan(sub(i), '-:')
      if (i_range > 0) then
        if (trim(sub(i)(1:i_range-1)) /= '') then
          read(sub(i)(1:i_range-1),*,iostat=stat) i1
          if (stat /= 0) then
            call error(stat)
            i1 = huge(i1)
          end if
        else
          i1 = 1
        end if
        if (trim(sub(i)(i_range+1:)) /= '') then
          read(sub(i)(i_range+1:),*,iostat=stat) i2
          if (stat /= 0) then
            call error(stat)
            i2 = -huge(i2)
          end if
        else
          call error(1)
          i2 = -huge(i2)
        end if
      else
        read(sub(i),*,iostat=stat) i1
        if (stat /= 0) then
          call error(stat)
          i1 = huge(i1)
          i2 = -huge(i2)
        else
          i2 = i1
        end if
      end if
      do j = i1, i2
        if (n_i + 1 > n_ip) call realloc
        n_i = n_i + 1
        ip(n_i) = j
      end do
    end do
    imax_ = maxval(ip(1:n_i))
    imin_ = minval(ip(1:n_i))

    if (ldeb_) then
      write(1000+pe_,*) 'intstr2array n_i',n_i
      flush 1000+pe_
    end if

    ! fill optional output variables
    if (present(n   )) n    = n_i
    if (present(imax)) imax = imax_
    if (present(imin)) imin = imin_
    if (present(iarr)) then
      j = min(size(iarr),n_i)
      iarr(1:j) = ip(1:j)
      if (size(iarr) < n_i) then
        write(msg, '("iarr too small: size(iarr)=",I8," n_i=",I8)') size(iarr),n_i
        call error(2, trim(msg))
      end if
    end if
    if (present(larr) .or. present(p_larr)) then
      allocate(lp(imin_:imax_))
      lp(:)         = .false.
      lp(ip(1:n_i)) = .true.
    end if
    if (present(larr)) then
      ib =lbound(larr,1)
      i1 = imin_
      if (ib < imin_) then
        larr(:i1-1) = .false.
      else if (ib > imin_) then
        i1 = ib
        if (present(status)) call error(3,'lower bound of larr too big')
      end if
      ib =ubound(larr,1)
      i2 = imax_
      if (ib > imax_) then
        larr(i2+1:) = .false.
      else if (ib < imax_) then
        i2 = ib
        if (present(status)) call error(4,'upper bound of larr too small')
      end if
      larr(i1:i2) = lp(i1:i2)
    end if
    if (present(p_larr)) then
      p_larr => lp
    else if (associated(lp)) then
      deallocate(lp)
    end if
    if (present(p_iarr)) then
      p_iarr => ip
    else
      deallocate(ip)
    end if

  contains

    subroutine realloc
      integer, pointer :: ip_(:) => null()
      allocate(ip_(n_ip+n_alloc))
      ip_(1:n_ip) = ip(1:n_ip)
      deallocate (ip)
      ip => ip_
      n_ip = size(ip)
    end subroutine realloc

    subroutine error(stat,msg)
      integer,          intent(in)           :: stat
      character(len=*), intent(in), optional :: msg
      if (present(status)) status = stat
      if (present(fstr  )) then
        if (present(msg)) then
          if (fstr == '') then
            fstr = trim(msg)
          else
            fstr = trim(fstr)//' and '//trim(msg)
          end if
        else
          if (fstr == '') then
            fstr = 'unable to interpret:'
          else
            fstr = trim(fstr)//' and'
          end if
          fstr = trim(fstr)//' "'//trim(sub(i))//'"'
        end if
      end if
    end subroutine error

  end subroutine intstr2array
!------------------------------------------------------------------------------
  subroutine concat (string, array)
  character(len=*) ,intent(out) :: string
  character(len=*) ,intent(in)  :: array (:)
  !-----------------------------------------------------------
  ! concatenate an array of words (array of character strings)
  ! words are separated by single blanks
  !-----------------------------------------------------------
  integer :: i, j, l
  string = ''
  j = 0
  do i = 1, size (array)
    l = len_trim (array(i))
    if (j+l > len(string)) call finish('concat','string is too small')
    string (j+1:j+l) = array (i) (1:l)
    j = j + l + 1
  end do
  end subroutine concat
!------------------------------------------------------------------------------
  subroutine unquote_string (s)
    !----------------------------
    ! Remove quotes around string
    !----------------------------
    character(len=*), intent(inout) :: s

    integer   :: k
    character :: c

    if (len (s) == 0) return
    s = adjustl (s)
    c = s(1:1)
    if (c /= "'" .and. c /= '"') return
    k = index (s(2:), c)
    if (k > 0) then
       s(k+1:k+1) = " "
    end if
    s = s(2:)
  end subroutine unquote_string
!------------------------------------------------------------------------------
  subroutine eval_string (dest, source, mode)
  character(len=*) ,intent(inout)        :: dest   ! destination string
  character(len=*) ,intent(in)           :: source ! source string
  integer          ,intent(in) ,optional :: mode   ! -1,0,1: start with: - = +
  !---------------------------------------------
  ! add or remove words from a string.
  ! special characters accounted for in source:
  !   = : dest is set to subsequent words in string
  !   + : subsequent words are added to dest
  !   - : subsequent words are removed from dest
  !   % : empty word; '= %' will delete dest
  !---------------------------------------------
    integer           :: m
    integer           :: nd, ns
    integer           :: id, is, i
    character(len=nc) :: d (nw)
    character(len=nw) :: s (nw)

    !--------------------------------------------------------------
    ! start with mode '=' if not given by optional parameter 'mode'
    !--------------------------------------------------------------
    m = 0; if (present(mode)) m = mode
    !-------------------------
    ! pack strings into arrays
    !-------------------------
    call split (d, dest,   nd)
    call split (s, source, ns)
    if (nd<0.or.ns<0) &
      call finish ('eval_string', 'array too small, increase nw!')
    !----------------------------
    ! loop over words in 'source'
    !----------------------------
    do is = 1, ns
      select case (s(is))
      !--------------
      ! change mode ?
      !--------------
      case ('=')
        m =  0
      case ('+')
        m =  1
      case ('-')
        m = -1
      case default
        !--------------------------------------------------------
        ! for mode '=' remove content of 'dest' and switch to '+'
        !--------------------------------------------------------
        if (m == 0) then
          m  = 1  ! switch to +
          nd = 0  ! remove 'dest'
        endif
        !-----------------------------------------------------------
        ! special character: % (empty word, = % deletes destination)
        !-----------------------------------------------------------
        if (s(is) == '%') cycle
        !-------------------------------------------------
        ! search for matching words in 'source' and 'dest'
        !-------------------------------------------------
        id = 0
        do i = 1, nd
          if (d(i) == s(is)) then
            id = i
            exit
          endif
        end do
        select case (m)
        case (1)
          !----------------------------------
          ! add word if not present in 'dest'
          !----------------------------------
          if (id==0) then
            if (nd==nw) &
              call finish ('eval_string', 'array too small, increase nw!')
            nd = nd + 1
            d(nd) = s(is)
          endif
        case (-1)
          !---------------------------------
          ! remove word if present in 'dest'
          !---------------------------------
          if (id>0) then
            d(id:nd-1) = d(id+1:nd)
            nd = nd - 1
          endif
        case default
          call finish ('eval_string','invalid value for mode')
        end select
      end select
    end do
    !-----------------------
    ! pack array into string
    !-----------------------
    call concat (dest, d(1:nd))

  end subroutine eval_string
!------------------------------------------------------------------------------
  subroutine decode_uuid (uuid_str, uuid)
    character(len=*), intent(in)  :: uuid_str   ! uuid encoded as string
    character(len=1), intent(out) :: uuid(:)    ! decoded uuid

    integer          :: i, j, l, n, b
    character(len=2) :: buf

    uuid(:) = achar (0)
    l = verify (uuid_str, "0123456789ABCDEFabcdef-")
    if (l > 0) then
       write (0,*) "Warning: invalid character in uuid: '", uuid_str(l:l),"'"
       return
    end if
    n = len  (uuid_str)
    i = 1
    j = 0
    do while (i < n)
       buf = uuid_str(i:i+1)
       if (buf(1:1) == "-") then
          i = i + 1                     ! Skip over dashes
          cycle
       end if
       i = i + 2
       read (buf,'(Z2)') b
       j = j + 1
       if (j > size (uuid)) call finish ("decode_uuid", "uuid input too long!")
       uuid(j) = achar (b)
    end do
    if (i == n) call finish ("decode_uuid", "uuid bad length")
  end subroutine decode_uuid
!------------------------------------------------------------------------------
END MODULE mo_dace_string
