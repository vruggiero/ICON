!
!+ Interface to the DWD BUFRX routines
!
MODULE mo_bufr_dwd
!
! Description:
!   Defines data type 't_bufr' to hold a BUFR record in encoded
!   or decoded representation. Defines operations on this data type
!   (interfaces to the DWD BUFRX routines). Automatically
!   allocates pointer array with appropriate size.
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
! V1_4         2009/03/26 Harald Anlauf
!  Preprocessor macro NOBUFR to disable linking to BUFR library
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Andreas Rhodin
!  provide double precision scale factor; increase accuracy of printout
! V1_22        2013-02-13 Andreas Rhodin
!  consistently use NF_FILL_FLOAT as real invalid value indicator
! V1_50        2017-01-09 Harald Anlauf
!  implement BUFRx encoder
! V1_51        2017-02-24 Michael Bender
!  bufr_encode: work around possible gfortran problem
!
! Code Description:
! Language: Fortran 95.
! Software Standards:
!
! Authors:
! Andreas Rhodin  DWD 2002  original code, ECMWF EMOS  library
! Andreas Rhodin  DWD 2003  change to      DWD   BUFR3 library
! Harald Anlauf   DWD 2008  change to      DWD   BUFRX library
!-------------------------------------------------------------
#if defined (__ICON__) && !defined (NOBUFR)
#define NOBUFR                   /* disable support BUFR outside DACE */
#endif
!=============
! Modules used
!=============
use mo_kind,      only: sp, dp   ! single, double precision kind parameter
use mo_exception, only: finish   ! abort on error
implicit none
!==============================================================================
!================
! Public entities
!================
  private
  !---------------------------------------------------------
  ! parameters for compatibility with old DWD raw I/O format
  !---------------------------------------------------------
  public :: raw                      ! write (cuegin compatible) DWD format
  public :: max_chars                ! max. no. of bytes to skip (>=8)
  !----------------------
  ! Data type definitions
  !----------------------
  public :: t_sec0                   ! Section 0
  public :: t_sec1                   ! Section 1
  public :: t_sec2                   ! Section 2
  public :: t_sec3                   ! Section 3
  public :: t_sec4                   ! Section 4
  public :: t_bufr                   ! Container
  !------------------------------------
  ! Routines acting on data type t_bufr
  !------------------------------------
  public :: bufr_open_file           ! open the BUFR file
  public :: bufr_close_file          ! close the BUFR file
  public :: bufr_read_bufr           ! read  a BUFR record
  public :: bufr_write_bufr          ! write a BUFR record
  public :: bufr_destroy             ! release memory of a BUFR record
  public :: bufr_get_sections        ! decode sections 1 to 4
  public :: bufr_get_data            ! decode, get data and meta-data
  public :: bufr_print_sections      ! print sections 0-4
  public :: bufr_print_subset        ! print body of a BUFR record
  public :: bufr_get_character       ! get character value
  public :: bufr_allocate            ! (re)allocate array components of t_bufr
  public :: bufr_deallocate          ! deallocate array components of t_bufr
  public :: bufr_encode              ! encode BUFR sections 0-4
  !--------------------------------------------------
  ! Interfaces to the low level DWD BUFR (C-)routines
  !--------------------------------------------------
  public :: bufrx_open_file          ! open the BUFR file
  public :: bufrx_close_file         ! close the BUFR file
  public :: bufrx_read_bufr          ! read  a BUFR record
  public :: bufrx_write_bufr         ! write a BUFR record
  public :: bufrx_destroy            ! release memory of a BUFR record
! public :: bufrx_get_sections       ! decode sections 1 to 4
  public :: bufrx_decode             ! decode the data
  public :: bufrx_get_data           ! actually get the data
  public :: bufrx_get_entry_descs    ! get meta-data
  public :: bufrx_get_character_value! get character datum
  public :: bufrx_set_section2_mode  ! set section 2 mode
  public :: bufrx_encode             ! encode the data
  public :: bufr_get_entry_texts     ! get text description
  public :: bufr_get_entry_units     ! get units description
  public :: bufrx_get_error_message  ! get error message
  !---------------------------------
  ! constant definitions for 'itype'
  !---------------------------------
  public :: ty_data       ! Normale Daten
  public :: ty_assoc      ! Associated field ("Quality bits")
  public :: ty_unknown    ! Unbekanntes Datenelement, 206YYY-Deskriptor
  public :: ty_char       ! Character-Daten, 206YYY-Deskriptor
  public :: ty_ref_ch     ! Geaenderter Referenzwert 204YYY-Deskriptor
  public :: ty_qual       ! Quality Assessment information
  public :: ty_loop_start ! Indikator fuer Beginn eines Loop
  public :: ty_loop_end   ! Indikator fuer Ende eines Loops
  public :: ty_loop_cnt   ! Indikator fuer Stand des Loop-Counters
  public :: ty_2XXYYY     ! Indikator fuer Deskriptor 2XXYYY (mit XX != 5)
  public :: ty_seq_start  ! Indikator fuer Start einer Sequenz
  public :: ty_seq_end    ! Indikator fuer Ende einer Sequenz
  !-------------------
  ! invalid data value
  !-------------------
  public :: inv_bufr
  !----------------------------------------------------------------------
  ! missing value indicators taken from previoud (EMOS library) interface
  !----------------------------------------------------------------------
  public :: nvind ! missing value indicator (integer)
!==============================================================================
!====================================
! Interfaces to the DWD BUFR routines
!====================================
  interface

    !-------------------
    ! open the BUFR file
    !-------------------
    subroutine bufrx_open_file (filename, output, raw, iunit, ier)
    implicit none
    character(*) ,intent(in)  :: filename  ! name of file to open
    integer      ,intent(in)  :: output    ! 0:input, 1:output
    integer      ,intent(in)  :: raw       ! 1:raw, 0:DWD format(cuegin comp.)
    integer      ,intent(out) :: iunit     ! file handle
    integer      ,intent(out) :: ier       ! error code: 0=OK
    end subroutine bufrx_open_file

    !--------------------
    ! close the BUFR file
    !--------------------
    subroutine bufrx_close_file (iunit, status, ier)
    implicit none
    integer          ,intent(in)  :: iunit     ! file handle
    character(len=*) ,intent(in)  :: status    ! 'del..' or 'DEL..' for delete
    integer          ,intent(out) :: ier       ! error flag: 0=OK
    end subroutine bufrx_close_file

    !-------------------
    ! read a BUFR record
    !-------------------
    subroutine bufrx_read_bufr (iunit, max_chars, ihandle, ieof, ier)
    implicit none
    integer       ,intent(in)  :: iunit     ! file handle
    integer       ,intent(in)  :: max_chars ! max. no. of bytes to skip (>=8)
    integer       ,intent(out) :: ihandle   ! BUFR record handle
    integer       ,intent(out) :: ieof      ! 0:OK, 1:EOF
    integer       ,intent(out) :: ier       ! error code: 0=OK
    end subroutine bufrx_read_bufr

    !--------------------
    ! write a BUFR record
    !--------------------
    subroutine bufrx_write_bufr (iunit, ihandle, ier)
    implicit none
    integer       ,intent(in)  :: iunit     ! file handle
    integer       ,intent(in)  :: ihandle   ! BUFR record handle
    integer       ,intent(out) :: ier       ! error code: 0=OK
    end subroutine bufrx_write_bufr

    !--------------------------------
    ! release memory of a BUFR record
    !--------------------------------
    subroutine bufrx_destroy (ihandle, ier)
    implicit none
    integer       ,intent(in)  :: ihandle   ! BUFR record handle
    integer       ,intent(out) :: ier       ! error flag: 0=OK
    end subroutine bufrx_destroy

    !-----------------------
    ! decode sections 1 to 4
    !-----------------------
!   !+++++++++++++++++++++++++++++++++++++++++++++++++++++
!   ! not activated because the Fortran Standart cuurently
!   ! prohibits to use the actual module 'use mo_bufr_dwd'
!   ! for a clean implementation wait for Fortran 200X
!   !+++++++++++++++++++++++++++++++++++++++++++++++++++++
!   subroutine bufrx_get_sections (ihandle, ibufsec0, ibufsec1, ibufsec2,&
!                                           ibufsec3, ibufsec4, ier)
!   use mo_bufr_dwd
!   implicit none
!   type(t_sec0) ,intent(in)  :: ihandle
!   type(t_sec1) ,intent(out) :: ibufsec0 ( 2)
!   type(t_sec2) ,intent(out) :: ibufsec1 (15)
!   type(t_sec3) ,intent(out) :: ibufsec2 (12)
!   type(t_sec4) ,intent(out) :: ibufsec3 ( 5)
!   type(t_sec5) ,intent(out) :: ibufsec4 ( 1)
!   integer      ,intent(out) :: ier
!   end subroutine bufrx_get_sections

    !----------------
    ! decode the data
    !----------------
    subroutine bufrx_decode (ihandle, mnem_list, mnem_list_len,             &
                                      idesc_list, idesc_list_len, ioptions, &
                                      max_subset_len, num_entry_descs, ier  )
    implicit none
    integer          ,intent(in)  :: ihandle         ! BUFR record handle
    character(len=*) ,intent(in)  :: mnem_list       ! mnemonic list
    integer          ,intent(in)  :: mnem_list_len   ! number of mnemonics
    integer          ,intent(in)  :: idesc_list(*)   ! descriptor list
    integer          ,intent(in)  :: idesc_list_len  ! size of idesc_list
    integer          ,intent(in)  :: ioptions        ! 1: no loop information
    integer          ,intent(out) :: max_subset_len  ! max. length of a subset
    integer          ,intent(out) :: num_entry_descs ! no. of entry descriptors
    integer          ,intent(out) :: ier             ! error flag: 0=OK
    end subroutine bufrx_decode

    !-------------
    ! get the data
    !-------------
    subroutine bufrx_get_data (ihandle, ibufdat, idescidx, nbufdat, &
                                        dim1, dim2, ier)
    implicit none
    integer ,intent(in)  :: ihandle              ! BUFR record handle
    integer ,intent(in)  :: dim1                 ! at least num_subsets
    integer ,intent(in)  :: dim2                 ! at least max_subset_len
    integer ,intent(out) :: ibufdat  (dim1,dim2) ! data
    integer ,intent(out) :: idescidx (dim1,dim2) ! index data
    integer ,intent(out) :: nbufdat  (dim1)      ! number of entries in subset
    integer ,intent(out) :: ier                  ! error flag: 0=OK
    end subroutine bufrx_get_data

    !------------------
    ! get the meta data
    !------------------
    subroutine bufrx_get_entry_descs (ihandle, itype, ifxy, is_char, iscale, &
                                      scale, nbits, irefval, ymnem, ier)
    use mo_kind, only: sp ! single precision kind parameter
    implicit none
    integer          ,intent(in)  :: ihandle     ! BUFR record handle
    integer          ,intent(out) :: itype   (*) ! type of entry-descriptor
    integer          ,intent(out) :: ifxy    (*) ! fxxyyy as decimal number
    integer          ,intent(out) :: is_char (*) ! flag for character
    integer          ,intent(out) :: iscale  (*) ! -exponent of scalingfactor
    real(sp)         ,intent(out) :: scale   (*) ! scaling factor
    integer          ,intent(out) :: nbits   (*) ! number of bits in BUFR
    integer          ,intent(out) :: irefval (*) ! reference value
    character(len=*) ,intent(out) :: ymnem   (*) ! Mnemonic
    integer          ,intent(out) :: ier         ! error flag: 0=OK
    end subroutine bufrx_get_entry_descs

    !--------------------
    ! get character datum
    !--------------------
    subroutine bufrx_get_character_value (ihandle, ibufdat, ychar, ier)
    implicit none
    integer          ,intent(in)  :: ihandle       ! BUFR record handle
    integer          ,intent(in)  :: ibufdat       ! integer BUFR value
    character(len=*) ,intent(out) :: ychar         ! character value
    integer          ,intent(out) :: ier           ! error flag: 0=OK
    end subroutine bufrx_get_character_value

    !--------------------------
    ! set section 2 mode
    ! mode=0: size(ibufsec2)=12
    ! mode=1: size(ibufsec2)=20
    !--------------------------
    subroutine bufrx_set_section2_mode (mode)
    implicit none
    integer          ,intent(in)  :: mode
    end subroutine bufrx_set_section2_mode

    !----------------
    ! encode the data
    !----------------
    subroutine bufrx_encode (ibufsec0, ibufsec1, ibufsec2, ibufsec3,      &
                             ibufsec4, idescr, nbufdat, ibufdat, ybufdat, &
                             setmne, ybufmne, ihandle, ier  )
    implicit none
    integer      ,intent(inout) :: ibufsec0(*)              ! Section 0 data (1:2)
    integer      ,intent(inout) :: ibufsec1(*)              ! Section 1 data (1:17)
    integer      ,intent(inout) :: ibufsec2(*)              ! Section 2 data (1:20)
    integer      ,intent(inout) :: ibufsec3(*)              ! Section 3 data (1:5)
    integer      ,intent(out)   :: ibufsec4(*)              ! Section 4 size (1:1)
    integer      ,intent(in)    :: idescr  (*)              ! descriptor list
    integer      ,intent(in)    :: nbufdat (*)              ! descriptor list
    integer      ,intent(in)    :: ibufdat (*)              ! input data to encode
    character(*) ,intent(in)    :: ybufdat (ibufsec3(3),*)  ! character data
    integer      ,intent(in)    :: setmne                   ! flag to set mnemonics
    character(*) ,intent(out)   :: ybufmne (ibufsec3(3),*)  ! mnemonic list
    integer      ,intent(out)   :: ihandle                  ! BUFR record handle
    integer      ,intent(out)   :: ier                      ! error flag: 0=OK
    end subroutine bufrx_encode

    !---------------------------------
    ! get description of error message
    !---------------------------------
    subroutine bufrx_get_error_message (text)
    implicit none
    character(len=*) ,intent(out) :: text         ! description of error
    end subroutine bufrx_get_error_message

  end interface
!==============================================================================
!==========================================
! data type definitions for sections 1 to 4
!==========================================

  !----------
  ! Section 0
  !----------
  type t_sec0
    integer ::        &
      num_bytes      ,&! ( 1) (*) BUFR-Laenge in Bytes
      edition          ! ( 2)     BUFR-Editions-Nummer
  end type t_sec0

  !----------
  ! Section 1 (layout according to BUFR edition 4)
  !----------
  type t_sec1
    integer ::        &
      len            ,&! ( 1) (*)  Laenge der Section 1 in Bytes
      table          ,&! ( 2)      BUFR-Master-Tabelle
      center         ,&! ( 3)      Erzeugendes Zentrum
      sub_center     ,&! ( 4)      Erzeugendes Sub-Zentrum
      update         ,&! ( 5)      Update-Sequenz-Nummer
      s2_present     ,&! ( 6)      Section 2-FLAG,
                       !           vorhanden                           = 128
                       !           nicht vorhanden                     =   0
      message_type   ,&! ( 7)      Daten-Kategorie nach BUFR-Tabelle A
      sub_category   ,&! ( 8)      Data sub-category  (common code table 13)
      message_subtype,&! ( 9)      Daten-Sub-Kategorie (lokale Tabelle)
      version        ,&! (10) (**) Versionsnummer der BUFR-Master-Tabellen
      local_version  ,&! (11) (**) Versionsnummer der lokalen Tabellen
      year           ,&! (12)      BUFR3: Jahr im Jahrhundert
                       !           100 ist Jahr 2000 im 20.Jahrhundert
                       !           converted to ccyy as defined in BUFR4
      month          ,&! (13)      Monat            (most typical date)
      day            ,&! (14)      Tag              (most typical date)
      hour           ,&! (15)      Stunde           (most typical date)
      minute         ,&! (16)      Minute           (most typical date)
      second           ! (17)      Sekunde
  end type t_sec1

  !----------
  ! Section 2
  !----------
  type t_sec2
    integer ::        &
      len            ,&! ( 1) (*)  Laenge der Section 2 in Bytes
      reserved       ,&! ( 2)      reserviert
      day_back       ,&! ( 3)      Tag des Backup-Files
      type_back      ,&! ( 4)      Art des Backup-Files
                       !           fuer ASCII                           =  0
                       !           fuer Binaer                          =  1
      record         ,&! ( 5)      Record-Nummer auf Backup-file
      machine        ,&! ( 6)      Markierung fuer Rechner,
                       !           fuer RUS1                            =  1
                       !           fuer RUS2                            =  0
      year_dec       ,&! ( 7)      Entschluesselungsjahr im Jahrhundert
      month_dec      ,&! ( 8)      Entschluesselungsmonat
      day_dec        ,&! ( 9)      Entschluesselungstag
      hour_dec       ,&! (10)      Entschluesselungsstunde
      min_dec        ,&! (11)      Entschluesselungsminute
      bank           ,&! (12)      Datenbank-Kennziffer
      center         ,&! (13)      generating center
      sequence       ,&! (14)      conversion sequence number
      year           ,&! (15)      year   of conversion
      month          ,&! (16)      month  of conversion
      day            ,&! (17)      year   of conversion
      hour           ,&! (18)      hour   of conversion
      minute         ,&! (19)      minute of conversion
      second           ! (20)      second of conversion
                       !           =-1E9, -> IKZ ist undefiniert
  end type t_sec2

  !----------
  ! Section 3
  !----------
  type t_sec3
    integer ::        &
      len            ,&! ( 1) (*)  Laenge der Section 3 in Bytes
      reserved       ,&! ( 2)      reserviert
      num_subsets    ,&! ( 3)      Anzahl einzelner Meldungen
      flag           ,&! ( 4)      Markierung (128 oder 0),
                       !           BIT1=1, -> Beobachtungen
                       !               =0, -> anderes
                       !           BIT2=1, -> komprimiert
                       !               =0, -> nicht komprimiert
      num_desc         ! ( 5)      Anzahl(k) der Deskriptoren
  end type t_sec3

  !----------
  ! Section 4
  !----------
  type t_sec4
    integer ::        &
    len                ! ( 1) (*)  Laenge der Section 4 in Bytes
  end type t_sec4

  !===================================
  ! data type holding temporary fields
  ! for BUFR decoding/encoding
  !===================================
  type t_bufr
    !------------------------------------------
    ! handles required by the DWD-BUFR routines
    !------------------------------------------
    integer          :: iunit   = -1             ! input  file handle
    integer          :: ounit   = -1             ! output file handle
    integer          :: ihandle = -1             ! BUFR record handle
    !----------------------------------------
    ! sections returned by bufrx_get_sections
    !----------------------------------------
    type (t_sec0)    :: sec0
    type (t_sec1)    :: sec1
    type (t_sec2)    :: sec2
    type (t_sec3)    :: sec3
    type (t_sec4)    :: sec4
    !-------------------------------------
    ! array sizes returned by bufrx_decode
    !-------------------------------------
    integer          :: max_subset_len           ! max. length of a subset
    integer          :: num_entry_descs          ! no. of entry descriptors
    !--------------------------------
    ! buffer sizes actually allocated
    !--------------------------------
    integer          :: dim1 = 0                 ! at least sec3% num_subsets
    integer          :: dim2 = 0                 ! at least max_subset_len
    integer          :: dim3 = 0                 ! at least num_entry_descs
    !-----------------------------
    ! data set by bufrx_get_data
    !-----------------------------
    integer ,pointer :: ibufdat  (:,:)=>NULL()   ! (dim1,dim2) data
    integer ,pointer :: idescidx (:,:)=>NULL()   ! (dim1,dim2) index data
    integer ,pointer :: nbufdat  (:)  =>NULL()   ! (dim1) number of entries
                                                 !        in subset
    !----------------------------------
    ! input data needed by bufrx_encode
    !----------------------------------
    character(len=64),pointer:: ybufdat(:,:)=>NULL() !(dim1,dim2) Character data
    !------------------------------------
    ! diagnostic data set by bufrx_encode
    !------------------------------------
    character(len=8), pointer:: ybufmne(:,:)=>NULL() !(dim1,dim2) Mnemonics of
                                                     !            encoded data
    !---------------------------------------
    ! meta-data set by bufrx_get_entry_descs
    !---------------------------------------
    integer          ,pointer:: itype  (:)=>NULL()!(dim3) type of entry-descr.
    integer          ,pointer:: ifxy   (:)=>NULL()!(dim3) fxxyyy as decimal no.
    integer          ,pointer:: is_char(:)=>NULL()!(dim3) flag for character
    integer          ,pointer:: iscale (:)=>NULL()!(dim3) -exponent of scaling.
    real(sp)         ,pointer:: scale  (:)=>NULL()!(dim3) scaling factor
    real(dp)         ,pointer:: dscale (:)=>NULL()!(dim3) " for higher accuracy
    integer          ,pointer:: nbits  (:)=>NULL()!(dim3) number of bits inBUFR
    integer          ,pointer:: irefval(:)=>NULL()!(dim3) reference value
    character(len=8) ,pointer:: ymnem  (:)=>NULL()!(dim3) Mnemonic
    character(len=64),pointer:: ytext  (:)=>NULL()!(dim3) Description
    character(len=16),pointer:: yunit  (:)=>NULL()!(dim3) Units

  end type t_bufr

!==============================================================================
!=================
! module variables
!=================
  integer            :: raw       = 0           ! write DWD format
  integer            :: max_chars = 256         ! max. no. of bytes to skip
  integer ,parameter :: inv_bufr  = -1000000000 ! invalid value in bufr message
  !---------------------------------
  ! constant definitions for 'itype'
  !---------------------------------
  integer ,parameter :: ty_data       =  0 ! Normale Daten
  integer ,parameter :: ty_assoc      =  1 ! Associated field ("Quality bits")
  integer ,parameter :: ty_unknown    =  2 ! Unbekanntes Datenelement
  integer ,parameter :: ty_char       =  3 ! Character-Daten, 206YYY-Deskriptor
  integer ,parameter :: ty_ref_ch     =  4 ! Geaenderter Referenzwert
  integer ,parameter :: ty_qual       =  5 ! Quality Assessment information
  integer ,parameter :: ty_loop_start =  6 ! Indikator fuer Beginn eines Loop
  integer ,parameter :: ty_loop_end   =  7 ! Indikator fuer Ende eines Loops
  integer ,parameter :: ty_loop_cnt   =  8 ! Indikator fuer Loop-Counter
  integer ,parameter :: ty_2XXYYY     =  9 ! Deskriptor 2XXYYY (mit XX != 5)
  integer ,parameter :: ty_seq_start  = 10 ! Indikator fuer Start einer Sequenz
  integer ,parameter :: ty_seq_end    = 11 ! Indikator fuer Ende einer Sequenz
  !-----------------------------------------------------------------------
  ! Missing value indicators taken from previoud (EMOS library) interface.
  ! Subject to change ?
  !-----------------------------------------------------------------------
  integer  ,parameter :: nvind = 2147483647 ! missing value indicator
!==============================================================================
contains
!====================================
! Routines acting on data type t_bufr
!====================================
!------------------------------------------------------------------------------
#ifndef NOBUFR
  !-------------------------------------------------
  ! allocate pointer components with sufficient size
  !-------------------------------------------------
  subroutine reallocate (bufr, dim1, dim2, dim3)

  type(t_bufr) ,intent(inout) :: bufr
  integer      ,intent(in)    :: dim1
  integer      ,intent(in)    :: dim2
  integer      ,intent(in)    :: dim3

    if (bufr% dim1 < dim1) then
      if (associated (bufr% nbufdat) ) deallocate (bufr% nbufdat)
      allocate (bufr% nbufdat (dim1))
    endif

    if (bufr% dim1 /= dim1 .or. bufr% dim2 < dim2) then
      if (associated (bufr% ibufdat) ) deallocate (bufr% ibufdat)
      if (associated (bufr% idescidx)) deallocate (bufr% idescidx)
      if (associated (bufr% ybufdat) ) deallocate (bufr% ybufdat)
      if (associated (bufr% ybufmne) ) deallocate (bufr% ybufmne)
      allocate (bufr% ibufdat  (dim1,dim2))
      allocate (bufr% idescidx (dim1,dim2))
      allocate (bufr% ybufdat  (dim1,dim2))
      allocate (bufr% ybufmne  (dim1,dim2))
    endif

    if (bufr% dim3 < dim3) then
      if (associated (bufr% itype)  ) deallocate (bufr% itype)
      if (associated (bufr% ifxy)   ) deallocate (bufr% ifxy)
      if (associated (bufr% is_char)) deallocate (bufr% is_char)
      if (associated (bufr% iscale) ) deallocate (bufr% iscale)
      if (associated (bufr% scale)  ) deallocate (bufr% scale)
      if (associated (bufr% dscale) ) deallocate (bufr% dscale)
      if (associated (bufr% nbits)  ) deallocate (bufr% nbits)
      if (associated (bufr% irefval)) deallocate (bufr% irefval)
      if (associated (bufr% ymnem)  ) deallocate (bufr% ymnem)
      if (associated (bufr% ytext)  ) deallocate (bufr% ytext)
      if (associated (bufr% yunit)  ) deallocate (bufr% yunit)
      allocate (bufr% itype   (dim3))
      allocate (bufr% ifxy    (dim3))
      allocate (bufr% is_char (dim3))
      allocate (bufr% iscale  (dim3))
      allocate (bufr% scale   (dim3))
      allocate (bufr% dscale  (dim3))
      allocate (bufr% nbits   (dim3))
      allocate (bufr% irefval (dim3))
      allocate (bufr% ymnem   (dim3))
    endif

    bufr% dim1 = dim1
    bufr% dim2 = dim2
    bufr% dim3 = dim3

  end subroutine reallocate

#endif
!------------------------------------------------------------------------------
#ifndef NOBUFR
  !------------------------------
  ! deallocate pointer components
  !------------------------------
  subroutine deallocate (bufr)
  type(t_bufr) ,intent(inout) :: bufr

    if (associated (bufr% nbufdat) ) deallocate (bufr% nbufdat)
    if (associated (bufr% ibufdat) ) deallocate (bufr% ibufdat)
    if (associated (bufr% idescidx)) deallocate (bufr% idescidx)
    if (associated (bufr% ybufdat) ) deallocate (bufr% ybufdat)
    if (associated (bufr% ybufmne) ) deallocate (bufr% ybufmne)
    if (associated (bufr% itype)   ) deallocate (bufr% itype)
    if (associated (bufr% ifxy)    ) deallocate (bufr% ifxy)
    if (associated (bufr% is_char) ) deallocate (bufr% is_char)
    if (associated (bufr% iscale)  ) deallocate (bufr% iscale)
    if (associated (bufr% scale)   ) deallocate (bufr% scale)
    if (associated (bufr% dscale)  ) deallocate (bufr% dscale)
    if (associated (bufr% nbits)   ) deallocate (bufr% nbits)
    if (associated (bufr% irefval) ) deallocate (bufr% irefval)
    if (associated (bufr% ymnem)   ) deallocate (bufr% ymnem)
    if (associated (bufr% ytext)   ) deallocate (bufr% ytext)
    if (associated (bufr% yunit)   ) deallocate (bufr% yunit)

    bufr% dim1 = 0
    bufr% dim2 = 0
    bufr% dim3 = 0

  end subroutine deallocate

#endif
!------------------------------------------------------------------------------
  subroutine bufr_print_sections (bufr)
  type(t_bufr) ,intent(in) :: bufr

  write(6,'()')
  write(6,'(a)') ' Section 0'
  write(6,'()')
  write(6,'(i8,a)')bufr%sec0% num_bytes,      ' BUFR-Laenge in Bytes'
  write(6,'(i8,a)')bufr%sec0% edition,        ' BUFR-Editions-Nummer'
  write(6,'()')
  write(6,'(a)') ' Section 1'
  write(6,'()')
  write(6,'(i8,a)')bufr%sec1% len,            ' Laenge der Section 1 in Bytes'
  write(6,'(i8,a)')bufr%sec1% table,          ' BUFR-Master-Tabelle'
  write(6,'(i8,a)')bufr%sec1% sub_center,     ' Erzeugendes Sub-Zentrum'
  write(6,'(i8,a)')bufr%sec1% center,         ' Erzeugendes Zentrum'
  write(6,'(i8,a)')bufr%sec1% update,         ' Update-Sequenz-Nummer'
  write(6,'(i8,a)')bufr%sec1% s2_present,     ' Section 2-FLAG'
  write(6,'(i8,a)')bufr%sec1% message_type,   ' Daten-Kategorie'
  write(6,'(i8,a)')bufr%sec1% message_subtype,' Daten-Sub-Kategorie'
  write(6,'(i8,a)')bufr%sec1% version,        ' Versionsnummer der Tabellen'
  write(6,'(i8,a)')bufr%sec1% local_version,  ' Versionsnummer lokale Tab.'
  write(6,'(i8,a)')bufr%sec1% year,           ' Jahr im Jahrhundert'
  write(6,'(i8,a)')bufr%sec1% month,          ' Monat'
  write(6,'(i8,a)')bufr%sec1% day,            ' Tag'
  write(6,'(i8,a)')bufr%sec1% hour,           ' Stunde'
  write(6,'(i8,a)')bufr%sec1% minute,         ' Minute'
  write(6,'()')
  write(6,'(a)') ' Section 2'
  write(6,'()')
  write(6,'(i8,a)')bufr%sec2% len,      ' Laenge der Section 2 in Bytes'
  write(6,'(i8,a)')bufr%sec2% reserved, ' reserviert'
  write(6,'(i8,a)')bufr%sec2% day_back, ' Tag des Backup-Files'
  write(6,'(i8,a)')bufr%sec2% type_back,' Art des Backup-Files'
  write(6,'(i8,a)')bufr%sec2% record,   ' Record-Nummer auf Backup-file'
  write(6,'(i8,a)')bufr%sec2% machine,  ' Markierung fuer Rechner'
  write(6,'(i8,a)')bufr%sec2% year_dec, ' Entschluesselungsjahr'
  write(6,'(i8,a)')bufr%sec2% month_dec,' Entschluesselungsmonat'
  write(6,'(i8,a)')bufr%sec2% day_dec,  ' Entschluesselungstag'
  write(6,'(i8,a)')bufr%sec2% hour_dec, ' Entschluesselungsstunde'
  write(6,'(i8,a)')bufr%sec2% min_dec,  ' Entschluesselungsminute'
  write(6,'(i8,a)')bufr%sec2% bank,     ' Datenbank-Kennziffer'
  write(6,'()')
  write(6,'(a)') ' Section 3'
  write(6,'()')
  write(6,'(i8,a)')bufr%sec3% len,        ' Laenge der Section 3 in Bytes'
  write(6,'(i8,a)')bufr%sec3% reserved,   ' reserviert'
  write(6,'(i8,a)')bufr%sec3% num_subsets,' Anzahl einzelner Meldungen'
  write(6,'(i8,a)')bufr%sec3% flag,       ' BIT1:Beobachtungen, 2:komprimiert'
  write(6,'(i8,a)')bufr%sec3% num_desc,   ' Anzahl(k) der Deskriptoren'
  write(6,'()')
  write(6,'(a)') ' Section 4'
  write(6,'()')
  write(6,'(i8,a)')bufr%sec4% len,        ' Laenge der Section 4 in Bytes'
  write(6,'()')

  end subroutine bufr_print_sections
!------------------------------------------------------------------------------
  subroutine bufr_print_subset (bufr, is)
  type (t_bufr) ,intent(in) :: bufr ! BUFR record
  integer       ,intent(in) :: is   ! subset index
  !-------------------------------------
  ! print body (subset) of a BUFR record
  !-------------------------------------
    integer           :: ie, id  ! indices
    integer           :: iv      ! decoded integer   value
    real(dp)          :: rv      ! decoded real      value
    character(len=40) :: cv      ! decoded character value

    !------------------
    ! loop over entries
    !------------------
    do ie = 1, bufr% nbufdat (is)
      id = bufr% idescidx (is, ie)
      iv = bufr% ibufdat  (is, ie)

      if (bufr% is_char(id) /= 0) then
        !----------------------
        ! print character value
        !----------------------
        call bufr_get_character (bufr, iv, cv)
        write(6,'(i6,1x,i2,1x,a8,1x,a30,1x,a,a40)') &
          bufr% ifxy (id), bufr% itype (id), bufr% ymnem (id), cv, &
          '['//bufr% yunit (id)//'] ', bufr% ytext (id)
      elseif (iv == inv_bufr) then
        !--------------------
        ! print invalid entry
        !--------------------
        write(6,'(i6,1x,i2,1x,a8,1x,a10,21x,a,a40)') &
          bufr% ifxy (id), bufr% itype (id), bufr% ymnem (id), &
          '   invalid', '['//bufr% yunit (id)//'] ', bufr% ytext (id)
      else
        !--------------------
        ! print numeric value
        !--------------------
        rv = bufr% dscale (id) * iv
        write(6,'(i6,1x,i2,1x,a8,1x,i10,f20.8,1x,a,a40)') &
          bufr% ifxy (id), bufr% itype (id), bufr% ymnem (id), iv, rv, &
          '['//bufr% yunit (id)//'] ', bufr% ytext (id)
      endif
    end do

  end subroutine bufr_print_subset
!------------------------------------------------------------------------------
  !-------------------
  ! open the BUFR file
  !-------------------
  subroutine bufr_open_file (bufr, filename, output, iunit, ier)
  type(t_bufr)     ,intent(inout) :: bufr
  character(len=*) ,intent(in)    :: filename ! name of file to open
  integer ,optional,intent(in)    :: output   ! 0:input, 1:output
  integer ,optional,intent(out)   :: iunit    ! file handle
  integer ,optional,intent(out)   :: ier      ! error flag: 0=OK

#ifdef NOBUFR
call finish ('bufr_open_file','BUFRx library not linked')
#else

    integer :: out, ie, iu
    out = 0; if (present(output)) out = output

    call bufrx_open_file (filename, out, raw, iu, ie)

    select case (out)
    case (0)
      bufr% iunit = iu
    case (1)
      bufr% ounit = iu
    case default
      call finish ('bufr_open_file','output is not 0 or 1')
    end select

    if (present(iunit)) iunit = iu
    if (present(ier)) then
      ier = ie
    else if (ie /= 0) then
      call error ('bufr_open_file (bufrx_open_file)')
    endif

#endif

  end subroutine bufr_open_file
!------------------------------------------------------------------------------
  !--------------------
  ! close the BUFR file
  !--------------------
  subroutine bufr_close_file (bufr, status, output, ier)
  type(t_bufr)     ,intent(inout)         :: bufr    ! BUFR variable
  character(len=*) ,intent(in)  ,optional :: status  ! 'del' to delete file
  integer          ,intent(in)  ,optional :: output  ! 0,1 for in/output only
  integer          ,intent(out) ,optional :: ier     ! 0: OK

#ifdef NOBUFR
call finish ('bufr_close_file','BUFRx library not linked')
#else

    character(len=3) :: stat
    integer          :: ie
    logical          :: i,o,d

    i    = bufr% iunit /= -1  ! close the input file  ?
    o    = bufr% ounit /= -1  ! close the output file ?
    d    = .true.             ! delete the BUFR variable

    stat = ''; if(present(status)) stat = status

    if (present (output)) then
      d = .false.
      select case (output)
      case (1)
        i = .false.
      case (0)
        o = .false.
      case default
        call finish ('bufr_close_file','output is not 0 or 1')
      end select
    endif

    if (bufr% iunit /= -1) then
      call bufrx_close_file (bufr% iunit, stat, ie)
      bufr% iunit = -1
      if (present(ier)) then
        ier = ie
      else if (ie /= 0) then
        call error ('bufr_close_file (input)')
      endif
    endif

    if (bufr% ounit /= -1) then
      call bufrx_close_file (bufr% ounit, stat, ie)
      bufr% ounit = -1
      if (present(ier)) then
        ier = ie
      else if (ie /= 0) then
        call error ('bufr_close_file (output)')
      endif
    endif

    call deallocate (bufr)

#endif

  end subroutine bufr_close_file
!------------------------------------------------------------------------------
  !-------------------
  ! read a BUFR record
  !-------------------
  subroutine bufr_read_bufr (bufr, ieof, ier)
  type(t_bufr)      ,intent(inout) :: bufr
  integer ,optional ,intent(out)   :: ieof      ! 0:OK, 1:EOF
  integer ,optional ,intent(out)   :: ier       ! error flag: 0=OK

#ifdef NOBUFR
call finish ('bufr_read_bufr','BUFRx library not linked')
#else

    integer :: i, ie

    call bufrx_read_bufr (bufr% iunit, max_chars, bufr% ihandle, i, ie)

    if (present(ieof)) then
      ieof = i
    else if (i/=0) then
      call finish ('bufr_read_bufr','EOF reached')
    endif

    if (present(ier)) then
      ier = ie
    else if (ie/=0) then
      call error ('bufr_read_bufr (bufrx_read_bufr)')
    endif

#endif

  end subroutine bufr_read_bufr
!------------------------------------------------------------------------------
  !--------------------
  ! write a BUFR record
  !--------------------
  subroutine bufr_write_bufr (bufr, ier)
  type(t_bufr)      ,intent(in)  :: bufr
  integer ,optional ,intent(out) :: ier       ! error flag: 0=OK

#ifdef NOBUFR
call finish ('bufr_write_bufr','BUFRx library not linked')
#else

    integer :: ie

    call bufrx_write_bufr (bufr% ounit, bufr% ihandle, ie)

    if (present(ier)) then
      ier = ie
    else if (ie/=0) then
      call error ('bufr_write_bufr (bufrx_write_bufr)')
    endif

#endif

  end subroutine bufr_write_bufr
!------------------------------------------------------------------------------
  !--------------------------------
  ! release memory of a BUFR record
  !--------------------------------
  subroutine bufr_destroy (bufr, ier)
  type(t_bufr)  ,intent(inout)           :: bufr
  integer       ,intent(out)   ,optional :: ier

#ifdef NOBUFR
call finish ('bufr_destroy','BUFRx library not linked')
#else

    integer :: ie

    call bufrx_destroy (bufr% ihandle, ie)
    bufr% ihandle = -1

    if (present(ier)) then
      ier = ie
    else if (ie/=0) then
      call error ('bufr_destroy (bufrx_destroy)')
    endif

#endif

  end subroutine bufr_destroy
!------------------------------------------------------------------------------
  !-----------------------
  ! decode sections 1 to 4
  !-----------------------
  subroutine bufr_get_sections (bufr, ier)
  type(t_bufr) ,intent(inout)         :: bufr
  integer      ,intent(out) ,optional :: ier

#ifdef NOBUFR
call finish ('bufr_get_sections','BUFRx library not linked')
#else

    integer :: ie, tmp
    integer :: sec1(17)

    call bufrx_get_sections (bufr% ihandle, bufr% sec0,       sec1,   &
                                bufr% sec2, bufr% sec3, bufr% sec4, ie)

    if (present (ier)) then
      ier = ie
    else if (ie /=0) then
      call error ('bufr_get_sections (bufrx_get_sections)')
    endif

    !------------------------------------------
    ! account for different layout of section 1
    ! in BUFR edition 3 and 4
    !------------------------------------------
    select case (bufr% sec0% edition)
    case (3)
      !-----------------------------------
      ! new in BUFR4: second, sub-category
      ! set to defaults if BUFR3 was read
      !-----------------------------------
      sec1 (  17) =   0         ! second
      sec1 (9:16) = sec1 (8:15)
      sec1 (8   ) = 255         ! data sub-category = undefined
      bufr% sec1 = transfer (sec1,bufr% sec1)
      !---------------------------------------------------
      ! change year in the sentury (BUFR3) to year (BUFR4)
      !---------------------------------------------------
      select case (bufr% sec1% year)
      case (  1:20)
        bufr% sec1% year = bufr% sec1% year + 2000
      case ( 21:99)
        bufr% sec1% year = bufr% sec1% year + 1900
      case (100   )
        bufr% sec1% year =                    2000
      case default
        call finish ('bufr_get_sections','invalid year')
      end select
      !--------------------------------------------------
      ! changed sequence of center,subcenter in edition 4
      !--------------------------------------------------
      tmp = bufr% sec1% sub_center
      bufr% sec1% sub_center = bufr% sec1% center
      bufr% sec1% center     = tmp
    case (4)
      bufr% sec1 = transfer (sec1,bufr% sec1)
    case default
      call finish ('bufr_get_sections','invalid section number (not 3,4):')
    end select

#endif

  end subroutine bufr_get_sections
!------------------------------------------------------------------------------
  !----------------
  ! decode the data
  !----------------
  subroutine bufr_get_data (bufr, mnem_list, ier)
  type(t_bufr)     ,intent(inout)           :: bufr
  character(len=8) ,intent(in)    ,optional :: mnem_list (:)
  integer          ,intent(out)   ,optional :: ier

#ifdef NOBUFR
call finish ('bufr_get_data','BUFRx library not linked')
#else

    integer                       :: i
    !----------------------------------
    ! flags for subroutine bufrx_decode
    !----------------------------------
    character(len=256) :: ml            ! mnemonic list
    integer            :: mll           ! number of mnemonics
    integer            :: ioptions = 0  ! 1: no loop information
    integer            :: desc_list(1)
    integer            :: desc_list_len
    integer            :: ie

    !---------------------------
    ! process optional arguments
    !---------------------------
    mll = 0
    if (present (mnem_list)) then
      do i=1,size(mnem_list)
        if (mnem_list(i)/=' ') then
          mll = mll + 1
          if (mll > len(ml)/8) &
            call finish('bufr_get_data','mnem_list is too long')
          ml((mll-1)*8+1:(mll-1)*8+8) = mnem_list (i)
        endif
      end do
    endif
    desc_list     = 0
    desc_list_len = 0

    !----------------
    ! decode the data
    !----------------
    call bufrx_decode (bufr% ihandle, ml(1:8), mll,                   &
                       desc_list, desc_list_len, ioptions,            &
                       bufr% max_subset_len, bufr% num_entry_descs, ie)

    if (ie /=0) then
      if (present (ier)) then
        ier = ie
        return
      else
        call error ('bufr_get_data (bufrx_decode)')
      endif
    endif

    !-----------------
    ! allocate buffers
    !-----------------
    call reallocate (bufr, bufr% sec3% num_subsets,   &
                           bufr%       max_subset_len,&
                           bufr%       num_entry_descs)

    !-------------
    ! get the data
    !-------------
    call bufrx_get_data (bufr% ihandle, &
                         bufr% ibufdat, bufr% idescidx, bufr% nbufdat, &
                         bufr% dim1, bufr% dim2, ie)

    if (ie /=0) then
      if (present (ier)) then
        ier = ie
        return
      else
        call error ('bufr_get_data (bufrx_get_data)')
      endif
    endif

    !------------------
    ! get the meta data
    !------------------
    call bufrx_get_entry_descs (bufr% ihandle, bufr% itype,   bufr% ifxy, &
                                bufr% is_char, bufr% iscale,  bufr% scale,&
                                bufr% nbits,   bufr% irefval, bufr% ymnem,&
                                ie)
    if (present (ier)) then
      ier = ie
    else if (ie /=0) then
      call error ('bufr_get_data (bufrx_get_entry_descs)')
    endif
    bufr% dscale = 10._dp**(-bufr% iscale)

#endif

  end subroutine bufr_get_data
!------------------------------------------------------------------------------
  !----------------
  ! decode the data
  !----------------
  subroutine bufr_get_character (bufr, ibufdat, ychar, ier)
  type(t_bufr)     ,intent(in)            :: bufr
  integer          ,intent(in)            :: ibufdat
  character(len=*) ,intent(out)           :: ychar
  integer          ,intent(out) ,optional :: ier

#ifdef NOBUFR
call finish ('bufr_get_character','BUFRx library not linked')
#else

    integer :: ie

    call bufrx_get_character_value (bufr% ihandle, ibufdat, ychar, ie)

    if (present(ier)) then
      ier = ie
    else if(ie/=0) then
      call error ('bufr_get_character (bufrx_get_character_value)')
    endif

#endif

  end subroutine bufr_get_character
!------------------------------------------------------------------------------
  subroutine bufr_get_entry_texts (bufr, ier)
  type(t_bufr) ,intent(inout)          :: bufr
  integer      ,intent(out)  ,optional :: ier

#ifdef NOBUFR
call finish ('bufr_get_entry_texts','BUFRx library not linked')
#else

    integer :: ie

    if (.not.associated (bufr% ytext)) &
      allocate (bufr% ytext (bufr% num_entry_descs))
    if (size(bufr% ytext) < bufr% num_entry_descs) then
      if (associated (bufr% ytext)) deallocate (bufr% ytext)
      allocate (bufr% ytext (bufr% num_entry_descs))
    endif

    call bufrx_get_entry_texts (bufr% ihandle, bufr% ytext, ie)

    if (present (ier)) then
      ier = ie
    else if (ie/=0) then
      call error ('bufr_get_entry_texts (bufrx_get_entry_texts)')
    endif

#endif

  end subroutine bufr_get_entry_texts
!------------------------------------------------------------------------------
  subroutine bufr_get_entry_units (bufr, ier)
  type(t_bufr) ,intent(inout)          :: bufr
  integer      ,intent(out)  ,optional :: ier

#ifdef NOBUFR
call finish ('bufr_get_entry_units','BUFRx library not linked')
#else

    integer :: ie

    if (.not.associated (bufr% yunit)) &
      allocate (bufr% yunit (bufr% num_entry_descs))
    if (size(bufr% yunit) < bufr% num_entry_descs) then
      if (associated (bufr% yunit)) deallocate (bufr% yunit)
      allocate (bufr% yunit (bufr% num_entry_descs))
    endif

    call bufrx_get_entry_units (bufr% ihandle, bufr% yunit, ie)

    if (present (ier)) then
      ier = ie
    else if (ie/=0) then
      call error ('bufr_get_entry_units (bufrx_get_entry_units)')
    endif

#endif

  end subroutine bufr_get_entry_units
!------------------------------------------------------------------------------
  subroutine bufr_allocate (bufr, num_subsets, subset_len, num_entry_descs)
  type(t_bufr) ,intent(inout) :: bufr
  integer      ,intent(in)    :: num_subsets      ! no. subsets
  integer      ,intent(in)    :: subset_len       ! >= max. subset length
  integer      ,intent(in)    :: num_entry_descs  ! no. entry descriptors

#ifdef NOBUFR
call finish ('bufr_allocate','BUFRx library not linked')
#else
    call reallocate (bufr, dim1=num_subsets,   &
                           dim2=subset_len,    &
                           dim3=num_entry_descs)
#endif

  end subroutine bufr_allocate
!------------------------------------------------------------------------------
  subroutine bufr_deallocate (bufr)
  type(t_bufr) ,intent(inout) :: bufr

#ifdef NOBUFR
call finish ('bufr_deallocate','BUFRx library not linked')
#else
    call deallocate (bufr)
#endif

  end subroutine bufr_deallocate
!------------------------------------------------------------------------------
  subroutine bufr_encode (bufr, idescr, ier)
    type(t_bufr) ,intent(inout)          :: bufr
    integer      ,intent(in)             :: idescr(:) !(dim3) descriptors
    integer      ,intent(out)  ,optional :: ier

#ifdef NOBUFR
call finish ('bufr_encode','BUFRx library not linked')
#else

    integer              :: ie
    integer              :: ms0, ms1, ms2, ms3,ms4
    integer, allocatable :: isec0(:), isec1(:), isec2(:), isec3(:), isec4(:)

    ms0 = size (transfer (bufr% sec0,[1]))
    ms1 = size (transfer (bufr% sec1,[1]))
    ms2 = size (transfer (bufr% sec2,[1]))
    ms3 = size (transfer (bufr% sec3,[1]))
    ms4 = size (transfer (bufr% sec4,[1]))
    allocate (isec0(max(ms0,2)), isec1(max(ms1,17)), isec2(max(ms2,20)) &
                               , isec3(max(ms3, 5)), isec4(max(ms4, 1)) )

    isec0         = inv_bufr
    isec1         = inv_bufr
    isec2         = inv_bufr
    isec3         = inv_bufr
    isec4         = inv_bufr
    isec0(1:ms0)  = transfer (bufr% sec0, isec0)
    isec1(1:ms1)  = transfer (bufr% sec1, isec1)
    isec2(1:ms2)  = transfer (bufr% sec2, isec2)
    isec3(1:ms3)  = transfer (bufr% sec3, isec3)
    isec4(1:ms4)  = transfer (bufr% sec4, isec4)
    bufr% ybufmne = ""

    call bufrx_set_section2_mode (1)

#ifdef __GFORTRAN__     /* work around possible gfortran problem */
    if (isec3(1) > 0) then
       print *, "isec3=", isec3
       write(*,*) size(bufr%ibufdat)
       write(*,*) size(bufr%idescidx)
       write(*,*) size(bufr%nbufdat)
    end if
#endif
    call bufrx_encode (isec0, isec1, isec2, isec3, isec4, idescr,   &
                       bufr% nbufdat, bufr% ibufdat, bufr% ybufdat, &
                       1, bufr% ybufmne, bufr% ihandle, ie          )

    if (present (ier)) then
      ier = ie
    else if (ie/=0) then
      call error ('bufr_encode (bufrx_encode)')
    endif

    bufr% sec0 = transfer (isec0(1:ms0), bufr% sec0)
    bufr% sec1 = transfer (isec1(1:ms1), bufr% sec1)
    bufr% sec2 = transfer (isec2(1:ms2), bufr% sec2)
    bufr% sec3 = transfer (isec3(1:ms3), bufr% sec3)
    bufr% sec4 = transfer (isec4(1:ms4), bufr% sec4)

#endif

  end subroutine bufr_encode
!------------------------------------------------------------------------------
#ifndef NOBUFR

  subroutine error (name)
  character(len=*) ,intent(in) :: name
  !-----------------------------------------
  ! abort in case of an error:
  ! call finish with BUFRX error description
  !-----------------------------------------
    character(len=256) :: err_text
    err_text = ''
    call bufrx_get_error_message (err_text)
    call finish (name, trim(err_text))
  end subroutine error

#endif
!------------------------------------------------------------------------------
end module mo_bufr_dwd
