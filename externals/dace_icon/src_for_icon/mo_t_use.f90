!
!+ Handle the 'flags' set for a report or datum.
!
MODULE mo_t_use
!
! Description:
!   This module handles the 'flags' set for a Report or Datum.
!
!   The actual settings are kept in variables of type 't_use'. This
!   datatype holds i) the 'state' (active, passive,..) of the report, ii)
!   the 'check' (first guess check, etc.) which caused the current state
!   of the report, and iii) a bit field 'flags' which indicates which of
!   the checks had a positive result (generally causing a degradation of
!   the state flag.
!
!   Constants are provided for valid values of the components state
!   (STAT_...) and check (CHK_...) .
!
!   Routines to decrease (decr_use) or increase (incr_use) the status
!   variable are provided.
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
!  ifdef USE_SEQUENCE
! V1_8         2009/12/09 Andreas Rhodin
!  new flag: CHK_BIASCOR = 3  (no bias correction possible)
! V1_9         2010/04/20 Andreas Rhodin
!  define status flag STAT_DEFAULT
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Detlef Pingel
!  changes for IASI, bias correction, cloud check
! V1_15        2011/12/06 Andreas Rhodin
!  define new flag "OPERATOR" (observation operator not applicable)
!  replace unused status flag value STAT_USED by STAT_OBS_ONLY
! V1_20        2012-06-18 Andreas Rhodin
!  reverse_bits: bugfix for btest(bit#32)
! V1_26        2013/06/27 Andreas Rhodin
!  replace STAT_ACTIVE_0 by STAT_NOTACTIVE
! V1_35        2014-11-07 Andreas Rhodin
!  remove CHK_MAXPR, introduce CHK_FG_LBC (fg check vs lateral bc. in COSMO)
! V1_47        2016-06-06 Andreas Rhodin
!  reverse_code: generic routine (vector or scalar);
!  account for status ST_DISMISS in DACE<->feedback translation table
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Andreas Rhodin  DWD  2004-2008  original source
!------------------------------------------------------------------------------

  !-------------
  ! Modules used
  !-------------
  use mo_kind,        only: i1, i4           ! kind parameters
  use mo_dace_string, only: toupper          ! convert string
  use mo_fdbk_tables, only: FL_OBSTYPE     ,&! Flag values: passive report type
                            FL_BLACKLIST   ,&! blacklist (or not whitelist)
                          ! FL_SUSP_LOCT   ,&! suspicious location or date/time
                            FL_TIME        ,&! time     not in valid range
                            FL_AREA        ,&! location not in valid area
                            FL_HEIGHT      ,&! '' not in valid height range
                            FL_SURF        ,&! incorrect surface (land,ice,etc)
                            FL_CLOUD       ,&! cloud check
                            FL_PRACTICE    ,&! bad reporting practice/insf.data
                            FL_DATASET     ,&! dataset quality flags
                            FL_REDUNDANT   ,&! redundant report
                          ! FL_FLIGHTTRACK ,&! flight track error flag
                            FL_MERGE       ,&! merged reports (e.g. TEMP ABCD)
                            FL_THIN        ,&! thinning
                            FL_RULE        ,&! complex rule
                            FL_OBS_ERR     ,&! observation error too large
                            FL_GROSS       ,&! gross error flag
                            FL_NO_BIASCOR  ,&! no bias correction possible
                            FL_FG          ,&! obs-fg check
                            FL_FG_LBC      ,&! obs-fg check vs. lateral bc.
                            FL_NO_OBS      ,&! no observations in report
                            FL_NONE        ,&! no flag set
                            FL_OPERATOR    ,&! observ. operator not applicable
                            ST_OBS_ONLY    ,&! no model equivalent available
                            ST_ACCEPTED    ,&! active and VQC accepted
                            ST_ACTIVE      ,&! used in the assimilation
                            ST_MERGED      ,&! not used, merged into multilevel
                            ST_PASSIVE     ,&! not used, only monitored
                            ST_REJECTED    ,&! not used, suspicious quality
                            ST_PAS_REJ     ,&! passive and rejected
                            ST_DISMISS       ! dismiss observation
  implicit none

  !----------------
  ! Public entities
  !----------------
  private
  public :: t_use        ! data type to hold the state of a report
  public :: use_0        ! default values of type use_0
  public :: decr_use     ! routine to decrease the state of a report or datum
  public :: incr_use     ! routine to increase the state of a report or datum
  public :: n_chk        ! number of states defined
  public :: n_stat       ! number of use-flag values defined
  public :: chk          ! mnemonics, default states, codes for 'checks'
  public :: t_chk        ! data type to hold chk
  public :: stats        ! table of      'states'
  public :: stat_mnem    ! mnemonics for 'states'
  public :: stat_key     ! key       for 'states'
  public :: chk_key      ! key       for 'checks'
  public :: change_code  ! convert code     (3dvar to feedback file convention)
  public :: change_bits  ! convert bit-field(3dvar to feedback file convention)
  public :: reverse_code ! convert code back     (feedback to 3dvar convention)
  public :: reverse_bits ! convert bit-field back(feedback to 3dvar convention)
  public ::             &! constants for valid state values.
            STAT_FORGET, STAT_DISMISS, STAT_OBS_ONLY,                   &
            STAT_PAS_REJ, STAT_PASSIVE, STAT_REJECTED,                  &
            STAT_NOTACTIVE, STAT_ACTIVE_0I, STAT_ACTIVE, STAT_ACTIVE_1I, &
            STAT_ACTIVE_1,                                              &
            STAT_ACCEPTED,                                              &
            STAT_INVALID, STAT_ABORT, STAT_MERGED,                      &
            STAT_DEFAULT
  public ::             &! constants for valid checks:
            CHK_NONE, CHK_BIASCOR, CHK_SUBTYP, CHK_TIME, CHK_AREA, CHK_HEIGHT, &
            CHK_CORR, CHK_CORRERR, CHK_REDUNDANT, CHK_DBLERR, CHK_INSDAT,      &
            CHK_OPERATOR, CHK_THIN, CHK_MERGE, CHK_FG,                         &
            CHK_RULE, CHK_FG_LBC, CHK_DOMAIN,                                  &
            CHK_SURF, CHK_DATASET, CHK_CLOUD, CHK_QI, CHK_NOTUSED,             &
            CHK_BLACKLIST, CHK_WHITELIST, CHK_FINAL, CHK_NOIMPL, CHK_OBS_ERR,  &
            CHK_CONSIST,CHK_GROSS,CHK_NO_OBS

  !===========
  ! Interfaces
  !===========
  interface reverse_code
    module procedure reverse_code_0  ! scalar version
    module procedure reverse_code_1  ! vector version
  end interface

  !=====================
  ! Constant definitions
  !=====================
  !-----------------------------------------------
  ! status of a report as result of various checks
  !-----------------------------------------------
  integer, parameter :: STAT_DEFAULT   =  0 ! to be replaced by other value
  integer, parameter :: STAT_ABORT     =  1 ! severe error, abort the program
  integer, parameter :: STAT_FORGET    =  2 ! don't use the report
  integer, parameter :: STAT_MERGED    =  3 !
  integer, parameter :: STAT_DISMISS   =  4 ! don't use, write to report file
  integer, parameter :: STAT_OBS_ONLY  =  5 !
  integer, parameter :: STAT_PAS_REJ   =  6 ! report is monitored + rejected
  integer, parameter :: STAT_PASSIVE   =  7 ! report is monitored
  integer, parameter :: STAT_REJECTED  =  8 ! rejected
  integer, parameter :: STAT_NOTACTIVE =  9 ! active, low weight in VQC
  integer, parameter :: STAT_ACTIVE_0I = 10 ! .. (deprecated, do not use!)
  integer, parameter :: STAT_ACTIVE    = 11 ! report is active
  integer, parameter :: STAT_ACTIVE_1I = 12 ! .. (deprecated, do not use!)
  integer, parameter :: STAT_ACTIVE_1  = 13 ! active, high weight in VQC
  integer, parameter :: STAT_ACCEPTED  = 14 ! accepted observation
  integer, parameter :: STAT_INVALID   = 15 ! invalid value
  integer, parameter :: n_stat         = 15 ! number of distinct states

  type t_chk
    character(len=12) :: mnem          ! mnemonic
    integer           :: default_state ! default state
    integer           :: code          ! code in feedback file
    logical           :: back          ! choice for backwards translation
  end type t_chk

  type (t_chk) ,parameter :: stats (n_stat) = &
    (/  t_chk('ABORT       ',STAT_ABORT    , -1          , .false. ),&
        t_chk('FORGET      ',STAT_FORGET   , -1          , .false. ),&
        t_chk('MERGED      ',STAT_MERGED   , ST_MERGED   , .true.  ),&
        t_chk('DISMISS     ',STAT_DISMISS  , ST_DISMISS  , .true.  ),&
        t_chk('OBS_ONLY    ',STAT_OBS_ONLY , ST_OBS_ONLY , .true.  ),&
        t_chk('PAS_REJ     ',STAT_PAS_REJ  , ST_PAS_REJ  , .true.  ),&
        t_chk('PASSIVE     ',STAT_PASSIVE  , ST_PASSIVE  , .true.  ),&
        t_chk('REJECTED    ',STAT_REJECTED , ST_REJECTED , .true.  ),&
        t_chk('NOTACTIVE   ',STAT_NOTACTIVE, ST_PASSIVE  , .false. ),&
        t_chk('ACTIVE_0I   ',STAT_ACTIVE_0I, ST_ACTIVE   , .false. ),&
        t_chk('ACTIVE      ',STAT_ACTIVE   , ST_ACTIVE   , .true.  ),&
        t_chk('ACTIVE_1I   ',STAT_ACTIVE_1I, ST_ACTIVE   , .false. ),&
        t_chk('ACTIVE_1    ',STAT_ACTIVE_1 , ST_ACTIVE   , .false. ),&
        t_chk('ACCEPTED    ',STAT_ACCEPTED , ST_ACCEPTED , .true.  ),&
        t_chk('INVALID     ',STAT_INVALID  , -1          , .true.  )/)

  character(len=12) ,parameter :: stat_mnem (0:n_stat) = &
    (/'DEFAULT     ',&
      'ABORT       ',&
      'FORGET      ',&
      'MERGED      ',&
      'DISMISS     ',&
      'OBS_ONLY    ',&
      'PAS_REJ     ',&
      'PASSIVE     ',&
      'REJECTED    ',&
      'NOTACTIVE   ',&
      'ACTIVE_0I   ',&
      'ACTIVE      ',&
      'ACTIVE_1I   ',&
      'ACTIVE_1    ',&
      'ACCEPTED    ',&
      'INVALID     '/)

  !----------------------------------------------------------------
  ! list of events (checks) which may change the status of a report
  !----------------------------------------------------------------
  integer, parameter :: CHK_NONE       =  1 !
  integer, parameter :: CHK_NOIMPL     =  2 ! not implemented
  integer, parameter :: CHK_BIASCOR    =  3 ! no bias correction possible
  integer, parameter :: CHK_SUBTYP     =  4 ! valid report subtype
  integer, parameter :: CHK_CORR       =  5 ! correction
  integer, parameter :: CHK_CORRERR    =  6 ! erroneous correction
  integer, parameter :: CHK_MERGE      =  7 ! merge reports (e.g. TEMP ABCD)
  integer, parameter :: CHK_REDUNDANT  =  8 ! double occurence of station
  integer, parameter :: CHK_DBLERR     =  9 ! erroneous double occurence
  integer, parameter :: CHK_OPERATOR   = 10 ! observation operator invalid
  integer, parameter :: CHK_INSDAT     = 11 ! insufficient data
  integer, parameter :: CHK_TIME       = 12 ! time in valid range
  integer, parameter :: CHK_AREA       = 13 ! location in valid area
  integer, parameter :: CHK_HEIGHT     = 14 ! location in valid height range
  integer, parameter :: CHK_THIN       = 15 ! thinning
  integer, parameter :: CHK_RULE       = 16 ! complex rule
  integer, parameter :: CHK_FG_LBC     = 17 ! max. number of processed reports
  integer, parameter :: CHK_GROSS      = 18 ! gross error
  integer, parameter :: CHK_SURF       = 19 ! correct surface ?
  integer, parameter :: CHK_DATASET    = 20 ! dataset flags (e.g.1D-Var flag)
  integer, parameter :: CHK_CLOUD      = 21 ! rejected due to cloud flag
  integer, parameter :: CHK_FG         = 22 ! obs-fg check
  integer, parameter :: CHK_OBS_ERR    = 23 ! obs/fg error check
  integer, parameter :: CHK_BLACKLIST  = 24 ! blacklist
  integer, parameter :: CHK_WHITELIST  = 25 ! whitelist
  integer, parameter :: CHK_NO_OBS     = 26 ! no observation in report
  integer, parameter :: CHK_QI         = 27 ! quality index
  integer, parameter :: CHK_NOTUSED    = 28 ! not used
  integer, parameter :: CHK_FINAL      = 29 ! final acceptance of datum
  integer, parameter :: CHK_DOMAIN     = 30 ! out of domain
  integer, parameter :: CHK_CONSIST    = 31 ! consistency in obs. profile
  integer, parameter :: n_chk          = 31 ! max value for CHK_..

  type (t_chk) ,parameter :: chk (n_chk) = &
    (/  t_chk('NONE        ',STAT_FORGET   , FL_NONE      , .true.  ),&
        t_chk('NOIMPL      ',STAT_ABORT    , FL_PRACTICE  , .false. ),&
        t_chk('BIASCOR     ',STAT_PASSIVE  , FL_NO_BIASCOR, .true.  ),&
        t_chk('SUBTYP      ',STAT_DISMISS  , FL_OBSTYPE   , .false. ),&
        t_chk('CORR        ',STAT_DISMISS  , FL_PRACTICE  , .false. ),&
        t_chk('CORRERR     ',STAT_DISMISS  , FL_PRACTICE  , .false. ),&
        t_chk('MERGE       ',STAT_MERGED   , FL_MERGE     , .true.  ),&
        t_chk('REDUNDANT   ',STAT_DISMISS  , FL_REDUNDANT , .true.  ),&
        t_chk('DBLERR      ',STAT_DISMISS  , FL_REDUNDANT , .false. ),&
        t_chk('OPERATOR    ',STAT_OBS_ONLY , FL_OPERATOR  , .true.  ),&
        t_chk('INSDAT      ',STAT_DISMISS  , FL_PRACTICE  , .false. ),&
        t_chk('TIME        ',STAT_FORGET   , FL_TIME      , .true.  ),&
        t_chk('AREA        ',STAT_DISMISS  , FL_AREA      , .true.  ),&
        t_chk('HEIGHT      ',STAT_DISMISS  , FL_HEIGHT    , .true.  ),&
        t_chk('THIN        ',STAT_FORGET   , FL_THIN      , .true.  ),&
        t_chk('RULE        ',STAT_DISMISS  , FL_RULE      , .true.  ),&
        t_chk('FG_LBC      ',STAT_REJECTED , FL_FG_LBC    , .true.  ),&
        t_chk('GROSS       ',STAT_REJECTED , FL_GROSS     , .true.  ),&
        t_chk('SURF        ',STAT_DISMISS  , FL_SURF      , .true.  ),&
        t_chk('DATASET     ',STAT_DISMISS  , FL_DATASET   , .true.  ),&
        t_chk('CLOUD       ',STAT_REJECTED , FL_CLOUD     , .true.  ),&
        t_chk('FG          ',STAT_REJECTED , FL_FG        , .true.  ),&
        t_chk('OBS_ERR     ',STAT_PASSIVE  , FL_OBS_ERR   , .true.  ),&
        t_chk('BLACKLIST   ',STAT_REJECTED , FL_BLACKLIST , .true.  ),&
        t_chk('WHITELIST   ',STAT_REJECTED , FL_BLACKLIST , .false. ),&
        t_chk('NO_OBS      ',STAT_REJECTED , FL_NO_OBS    , .true.  ),&
        t_chk('QI          ',STAT_DISMISS  , FL_DATASET   , .false. ),&
        t_chk('NOTUSED     ',STAT_DISMISS  , FL_OBSTYPE   , .true.  ),&
        t_chk('FINAL       ',STAT_ACCEPTED , FL_NONE      , .false. ),&
        t_chk('DOMAIN      ',STAT_FORGET   , FL_AREA      , .false. ),&
        t_chk('CONSIST     ',STAT_DISMISS  , FL_PRACTICE  , .true.  )/)
  !---------------------
  ! Data type definition
  !---------------------
  type t_use
#ifdef USE_SEQUENCE
    SEQUENCE
#endif
    integer(i1) :: state = STAT_INVALID ! status of the report
    integer(i1) :: check = 0            ! event which lead to the actual status
    integer(i4) :: flags = 0            ! bit field of events
  end type t_use

  !----------------
  ! module variable
  !----------------
  type(t_use), parameter :: use_0 = t_use(STAT_INVALID,0,0)! default settings

!==============================================================================
contains
!==============================================================================
  function stat_key (mnem) result (key)
  character(len=*) ,intent(in) :: mnem
  integer                      :: key
    character(len=12) :: upper
    integer           :: i
    upper = toupper (mnem)
    do i = 0, n_stat
      if (upper == stat_mnem(i)) then
        key = i
        return
      endif
    end do
    key = -1
  end function stat_key
!------------------------------------------------------------------------------
  function chk_key (mnem) result (key)
  character(len=*) ,intent(in) :: mnem
  integer                      :: key
    character(len=12) :: upper
    integer           :: i
    upper = toupper (mnem)
    do i = 1, size (chk)
      if (upper == chk(i)% mnem) then
        key = i
        return
      endif
    end do
    key = -1
  end function chk_key
!------------------------------------------------------------------------------
  subroutine change_bits (bits, table)
  integer ,intent(inout) :: bits  (:)
  integer ,intent(in)    :: table (:)
  !------------------------------------------------------
  ! convert bit-field (3dvar to feedback file convention)
  !------------------------------------------------------
    integer :: i,n,tmp
    do n=1,size(bits)
      tmp     = bits(n)
      bits(n) = 0
      do i=1,size(table)
        if (btest(tmp,i)) then
          if (table(i)>=0) then
            bits(n) = ibset(bits(n),table(i))
          endif
        endif
      end do
    end do
  end subroutine change_bits
!------------------------------------------------------------------------------
  subroutine reverse_bits (bits, table, unknown)
  integer     ,intent(inout)        :: bits  (:) ! bit-field to convert
  type(t_chk) ,intent(in)           :: table (:) ! translation table
  integer     ,intent(in) ,optional :: unknown   ! to be used for unknown code
  !-----------------------------------------------------------
  ! convert bit-field back (feedback file to 3dvar convention)
  !-----------------------------------------------------------
    integer :: n   ! bit-field array loop index
    integer :: i   ! bit number in 3dvar
    integer :: tmp ! code in feedback file
    integer :: j   ! bit number in feedback file
    do n=1,size(bits)
      tmp     = bits(n)
      bits(n) = 0
      do i=1,size(table)
        j = table(i)% code
        if (j < 0 .or. j >= bit_size(tmp) .or. .not. table(i)% back) cycle
        if (btest(tmp,j)) then
          bits(n) = ibset (bits(n),i)
          tmp     = ibclr (tmp    ,j)
        endif
      end do
      if (tmp /= 0 .and. present (unknown)) bits(n) = ibset (bits(n),unknown)
    end do
  end subroutine reverse_bits
!-----------------------------------------------------------------------------
  subroutine change_code (code, table)
  integer ,intent(inout) :: code  (:)
  integer ,intent(in)    :: table (:)
  !-------------------------------------------------
  ! convert code (3dvar to feedback file convention)
  !-------------------------------------------------
    integer :: n
    do n=1,size(code)
      code(n) = table(code(n))
    end do
  end subroutine change_code
!-----------------------------------------------------------------------------
  subroutine reverse_code_0 (code, table, unknown)
  integer     ,intent(inout)        :: code      ! code to convert
  type(t_chk) ,intent(in)           :: table (:) ! translation table
  integer     ,intent(in) ,optional :: unknown   ! to be used for unknown code
  !-----------------------------------------------------
  ! convert code back (feedback file to DACE convention)
  !-----------------------------------------------------
    integer :: c(1)
    c = code
    call reverse_code_1 (c, table, unknown)
    code = c(1)
  end subroutine reverse_code_0
!-----------------------------------------------------------------------------
  subroutine reverse_code_1 (code, table, unknown)
  integer     ,intent(inout)        :: code  (:) ! code to convert
  type(t_chk) ,intent(in)           :: table (:) ! translation table
  integer     ,intent(in) ,optional :: unknown   ! to be used for unknown code
  !-----------------------------------------------------
  ! convert code back (feedback file to DACE convention)
  !-----------------------------------------------------
    integer :: n   ! bit-field array loop index
    integer :: i   ! bit number in 3dvar
    integer :: c   ! code in feedback file (source)
    if (size (table) == 0) return
    do n=1,size(code)
      c = code(n)
      code(n) = -1
      do i = 1, size (table)
        if (table(i)% back .and. table(i)% code == c) then
          code(n) = i
          exit
        endif
      end do
      if (code(n) == -1 .and. present (unknown)) code(n) = unknown
      if (code(n) == -1) code(n) = -c
    end do
  end subroutine reverse_code_1
!==============================================================================
  elemental subroutine decr_use (u, state, check, lflag)
  type (t_use) ,intent(inout)        :: u     ! use state variable
  integer      ,intent(in) ,optional :: state ! state (active, passive,..) to set
  integer      ,intent(in)           :: check ! check which caused the state
  logical      ,intent(in) ,optional :: lflag ! set flag ?
  !-----------------------------------------
  ! decrease the state of a report or datum.
  !-----------------------------------------
    logical :: lf
    integer :: st

    st = STAT_DEFAULT
    if (present(state))     st = state
    if (st == STAT_DEFAULT) st = chk(check)% default_state

    if (u% state == STAT_PASSIVE .and. &
        st       == STAT_REJECTED      ) st = STAT_PAS_REJ
    if (st       == STAT_PASSIVE .and. &
        u% state == STAT_REJECTED      ) st = STAT_REJECTED
    if (u% state > st) then
      u% state = st
      u% check = check
    endif
    lf = st <= STAT_NOTACTIVE; if (present(lflag)) lf = lflag
    if (lf) u% flags = ibset (u% flags, check)
  end subroutine decr_use
!------------------------------------------------------------------------------
  subroutine incr_use (u, state, check)
  type (t_use) ,intent(inout) :: u     ! use state variable
  integer      ,intent(in)    :: state ! state (active, passive,..) to set
  integer      ,intent(in)    :: check ! check which caused the state
  !-----------------------------------------
  ! increase the state of a report or datum.
  !-----------------------------------------
    if (u% state < state .or. u% state == STAT_INVALID) then
      u% state  = state
      u% check = check
    endif
    u% flags = ibset (u% flags, check)
  end subroutine incr_use
!==============================================================================
end module mo_t_use
