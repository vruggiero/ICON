!
!+ 3DVAR/COSMO tables relating names(mnemonics) + numerical values or bit flags
!
! $Id$
!-------------------------------------------------------------------------------
!
MODULE mo_t_table
!
!------------------------------------------------------------------------------
! Description:
!   Derived data types and operators to maintain tables relating names
!   (mnemonics) and numerical values or bit flags. These tables are used
!   in the interface routines for processing the feedback file.
!   This module is used commonly by the COSMO model and 3DVAR program packages !
!
! Current Code Owners:
!    For DWD 3DVAR:                        For COSMO:
!    DWD, Andreas Rhodin                   DWD, Christoph Schraff
!    phone: +49 69 8062 2722               phone: +49 69 8062 2725
!    fax:   +49 69 8062 3721               fax:   +49 69 8062 3721
!    email: andreas.rhodin@dwd.de          email: christoph.schraff@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_1         2008/11/05 Andreas Rhodin
!  First operational 3D-Var release
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Andreas Rhodin
!  fix spelling errors
! V1_18        2012/01/06 Andreas Rhodin
!  Extended comments and documentation
! V1_22        2013-02-13 Christoph Schraff
!  function bit1_pos_nr included in mo_t_table
! V1_26        2013/06/27 Andreas Rhodin
!   new function text_value to derive comment line from t_table
! V1_27        2013-11-08 Andreas Rhodin
!  cleanup
! V1_29        2014/04/02 Andreas Rhodin
!  use longtable instead of tabular for tables in feedback file doku
! V1_37        2014-12-23 Andreas Rhodin
!  homogenize interface with COSMO V5.1
! V1_45        2015-12-15 Andreas Rhodin
!  cleanup
!
! CAUTION: This module is used commonly by the 3DVAR and COSMO main programs.!!!
!!!        Therefore, anybody wanting to introduce a modification to this    !!!
!!!        module in the context of either of these programs must consult    !!!
!!!        the 'current code owner' of this module for the other program,    !!!
!!!        in order to allow for checking that the modification will comply  !!!
!!!        with both program packages. This must be done before the          !!!
!!!        modification is put into the Version Control System (VCS).        !!!
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Andreas Rhodin  DWD  2007  original source to be used in 3DVAR/LETKF/COSMO
!------------------------------------------------------------------------------
implicit none

!----------------
! Public entities
!----------------
private
public :: t_table          ! table       derived type
public :: t_entry          ! table entry derived type
public :: position         ! find position of table entry (name or value given)
public :: name_value       ! find name        of table entry (value given)
public :: text_value       ! find description of table entry (value given)
public :: value_name       ! find value       of table entry (name given)
public :: INVALID_VALUE    ! returned by 'value_name' if entry is not present
public :: init_table       ! initialise the table table
public :: bit1_pos_nr      ! first bit = 1 in a word searched in a given order
public :: NLEN, ULEN, DLEN ! length of character strings used in t_table

  !----------------------------------------------------------------------
  ! Constants:
  ! Definition of charachter string lengths for table entries
  !            of invalid value returned by routine 'value_name'
  !             for the case that the requested table entry is not present
  !----------------------------------------------------------------------
  integer,   parameter :: NLEN          =   16 ! length of entry name
  integer,   parameter :: ULEN          =   12 ! length of units field
  integer,   parameter :: DLEN          =   64 ! length of entry description
  integer,   parameter :: INVALID_VALUE = -999 ! returned by 'value_name'
  character, parameter :: BS            = '\\'(1:1) ! portable backslash

  !------------------------------------------
  ! Derived type definition for table entries
  !------------------------------------------
  type t_entry
    integer             :: value        ! numerical value or bit number
    character(len=NLEN) :: name         ! associated name (mnemonic)
    character(len=ULEN) :: units        ! units of the numerical value
    character(len=DLEN) :: description  ! some text
  end type t_entry

  !--------------------------------------
  ! Derived type definition for the table
  !--------------------------------------
  type t_table
    type (t_entry) ,pointer :: e(:)    ! list of table entries
    character(len=16)       :: name    ! name of table
    character(len=128)      :: caption ! caption (for LaTEX doc)
    integer                 :: n       ! number of entries
    integer                 :: first   ! smallest value or bit number in table
    integer                 :: last    ! largest  value or bit number in table
    logical                 :: bit     ! used as bit-table
  end type t_table

  !-------------------------------------------------------
  ! interface for subroutine position:
  !   determines the position of the entry in the table
  !   either from the name or from the value or bit number
  !-------------------------------------------------------
  interface position
    module procedure position_name
    module procedure position_value
  end interface position

!==============================================================================
contains
!------------------------------------------------------------------------------
  elemental function position_name (table, name)
  type (t_table)   ,intent(in) :: table
  character(len=*) ,intent(in) :: name
  integer                      :: position_name
  !----------------------------------------------------------------
  ! determines the position of the entry in the table from its name
  !----------------------------------------------------------------
    integer :: i
    position_name = 0
    do i = 1, size(table% e)
      if (table% e(i)% name == name) then
        position_name = i
        exit
      endif
    end do
  end function position_name
!------------------------------------------------------------------------------
  elemental function position_value (table, value)
  type (t_table) ,intent(in) :: table
  integer        ,intent(in) :: value
  integer                    :: position_value
  !--------------------------------------------------
  ! determines the position of the entry in the table
  ! from its value or bit number
  !--------------------------------------------------
    integer :: i
    position_value = 0
    do i = 1, size(table% e)
      if (table% e(i)% value == value) then
        position_value = i
        exit
      endif
    end do
  end function position_value
!------------------------------------------------------------------------------
  elemental function value_name (table, name)
  type (t_table)   ,intent(in) :: table
  character(len=*) ,intent(in) :: name
  integer                      :: value_name
  !------------------------------------------------------------------
  ! determines the value or bit number of a table entry from its name
  !------------------------------------------------------------------
    integer :: i
    value_name = INVALID_VALUE
    do i = 1, size(table% e)
      if (table% e(i)% name == name) then
        value_name = table% e(i)% value
        exit
      endif
    end do
  end function value_name
!------------------------------------------------------------------------------
  elemental function name_value (table, value)
  type (t_table)      ,intent(in) :: table
  integer             ,intent(in) :: value
  character(len=NLEN)             :: name_value
  !------------------------------------------------------------------
  ! determines the name of a table entry from its value or bit number
  !------------------------------------------------------------------
    integer :: i
    name_value = ''
    do i = 1, size(table% e)
      if (table% e(i)% value == value) then
        name_value = table% e(i)% name
        exit
      endif
    end do
  end function name_value
!------------------------------------------------------------------------------
  elemental function text_value (table, value)
  type (t_table)      ,intent(in) :: table
  integer             ,intent(in) :: value
  character(len=DLEN)             :: text_value
  !------------------------------------------------------------------
  ! determines the name of a table entry from its value or bit number
  !------------------------------------------------------------------
    integer :: i
    text_value = ''
    do i = 1, size(table% e)
      if (table% e(i)% value == value) then
        text_value = table% e(i)% description
        exit
      endif
    end do
  end function text_value
!==============================================================================
  subroutine init_table (table, entries, name, caption, bit, latex, csv)
  type (t_table)               ,pointer  :: table      ! table
  type (t_entry)   ,intent(in) ,target   :: entries(:) ! table entries
  character(len=*) ,intent(in)           :: name       ! name of table
  character(len=*) ,intent(in)           :: caption    ! caption (for TeX)
  logical          ,intent(in)           :: bit        ! used as bit flag
  logical          ,intent(in) ,optional :: latex      ! write LaTeX file
  logical          ,intent(in) ,optional :: csv        ! write CSV file
  !---------------------------------------------------------------
  ! Initialise the table:
  !
  !   1) allocate pointer 'table'
  !   2) link table with table entries
  !   3) set NAME, CAPTION and BIT components
  !   4) optionally write a latex file tab.NAME.tex
  !   5) optionally write a CSV file tab.NAME.csv
  !      (CSV - comma separated values,
  !             ASCII file with several values per row)
  !
  !   Bit should be set to .true. if the table entries are to be
  !   interpreted as bit positions.
  !
  !   In the LaTeX output the following characters are replaced:
  !     '_'  by  '\_'
  !     '%'  by  '\%'
  !     '#'  by  '\#'
  !     '~' by   '\hfill '
  !-------------------------------------------------------------

    !-------------------------------------------------
    ! local variables (used for LaTeX processing only)
    !-------------------------------------------------
    integer           :: i       ! table entry index
    logical           :: lunits  ! true if 'units' is given for an entry
    character(len=49) :: legend  ! legend of table (1st row)
    character(len=5)  :: cunits  ! 'units' as character string
    character(len=5)  :: cvalue  ! value or bit number as character string
    logical           :: tex     ! local copy of argument 'latex'
    logical           :: csval   ! local copy of argument 'csv'

    !-----------------------------------
    !   1) allocate pointer 'table'
    !   2) link table with table entries
    !   3) set remaining components
    !-----------------------------------
    allocate (table)
    table% e       => entries
    table% name    =  name
    table% caption =  caption
    table% n       =  size   (entries)
    table% first   =  minval (entries% value)
    table% last    =  maxval (entries% value)
    table% bit     =  bit

    !----------------------------------------------
    ! 4) optionally write a latex file tab.NAME.tex
    !----------------------------------------------
    tex  = .false.; if (present(latex)) tex = latex
    if (tex) then

      !---------------
      ! prepare header
      !---------------
      lunits = any (entries% units /= '')
      cunits = ''      ;if(lunits) cunits = 'units'
      cvalue = 'value' ;if(bit)    cvalue = 'bit'//BS//'#'
      legend = cvalue//' & name & '//cunits//'  & description '//BS//BS

      !------------------------
      ! open file, write header
      !------------------------
      open  (1,file = 'tab.'//trim(table% name)//'.tex')
      write (1,'(a)') BS//'begin{center}'
      write (1,'(a)') BS//'begin{longtable}{rlcl}'
      write (1,'(a)') trim(legend)
      write (1,'(a)') ' & & & '//BS//BS
      write (1,'(a)') BS//'hline'
      write (1,'(a)') ' & & & '//BS//BS
      !---------------------------
      ! write rows (table entries)
      !---------------------------
      do i = 1, size (table% e)
        write (1,'(i5," '//BS//'quad & ",a," & ",a," & ",a," '//BS//BS//'")') &
                    table% e(i)% value,                                       &
          trim(fix_tex(table% e(i)% name)),                                   &
          trim(fix_tex(table% e(i)% units)),                                  &
          trim(fix_tex(table% e(i)% description))
      end do
      !--------------
      ! write trailer
      !--------------
      write (1,'(a)') ' & & &'//BS//BS
      write (1,'(a)') BS//'hline'
      write (1,'(a)') BS//'caption{'
      write (1,'(a)') trim(table% caption)//'}'
      write (1,'(a)') BS//'label{tab.'//trim(table% name)//'}'
      write (1,'(a)') BS//'end{longtable}'
      write (1,'(a)') BS//'end{center}'
    endif
    !----------------------------------------------
    ! 5) optionally write a CSV file tab.NAME.csv
    !----------------------------------------------
    csval  = .false.; if (present(csv)) csval = csv
    if (csval) then
      open  (1,file = 'tab.'//trim(table% name)//'.csv')

      !---------------------------------------
      ! write first row with table description
      !---------------------------------------
      write (1,'(i5,",""",a,""",""",a,""",""",a,"""")') &
           -999,                                        &
           trim(table% name),                           &
           '',                                          &
           trim(table% caption)

      !---------------------------
      ! write rows (table entries)
      !---------------------------
      do i = 1, size (table% e)
        write (1,'(i5,",""",a,""",""",a,""",""",a,"""")')   &
               table% e(i)% value,                          &
          trim(table% e(i)% name),                          &
          trim(table% e(i)% units),                         &
          trim(table% e(i)% description)
      end do

      close(1)
    end if
  !----------------------------------------------------------------------------
  contains
    function fix_tex(string)
    character(len=*)              :: string
    character(len=len(string)+20) :: fix_tex
    !-------------------------------------------------------------
    !   In the LaTeX output the following characters are replaced:
    !     '_'  by  '\_'
    !     '%'  by  '\%'
    !     '#'  by  '\#'
    !     '~' by   '\hfill '
    !-------------------------------------------------------------
      integer :: lout
      integer :: i,j
      lout = len (fix_tex)
      j = 0
      fix_tex = ''
      do i = 1, len(string)
        j=j+1
        if(j>LOUT) exit
        if (string(i:i) == '~') then
          if (j<LOUT-6) then
            fix_tex(j:j+6) = BS//'hfill '
            j=j+6
            cycle
          endif
        endif
        select case (string(i:i))
        case ('_','#','%')
          if (j<LOUT) then
            fix_tex(j:j) = BS
            j=j+1
          else
            exit
          endif
        end select
        fix_tex(j:j) = string(i:i)
      end do
    end function fix_tex
  end subroutine init_table
!==============================================================================
  integer function bit1_pos_nr ( iword, norder, iorder )

  integer ,intent(in) :: iword          ! input word
  integer ,intent(in) :: norder         ! number of bits to be searched through
  integer ,intent(in) :: iorder(norder) ! order  of bits to be searched through
  !--------------------------------------------------------
  ! searches through the bits in 'iword' in the order given
  ! by the vector 'iorder(norder)' and returns the position
  ! of the first bit which is found to be equal to 1;
  ! if all bits are zero, return value is the last bit
  ! according to the order given by 'iorder'
  ! (this function is used to get 'check' from 'flags')
  !--------------------------------------------------------
  integer  :: jj     ! loop index
  integer  :: ibps   ! bit position
  !---------------------------------- End of header -------

  bit1_pos_nr = iorder(norder)
  do jj = 1 , norder
    ibps = iorder(jj)
!   if (iand( ishft( iword,-ibps ), 1 ) == 1) then
    if (btest( iword, ibps )) then
      bit1_pos_nr = ibps
                                                                           exit
    endif
  enddo
  end function bit1_pos_nr

!------------------------------------------------------------------------------
end module mo_t_table
