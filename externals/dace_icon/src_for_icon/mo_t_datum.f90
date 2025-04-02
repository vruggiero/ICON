!
!+ Derived type definition and operations for a single observed entity
!
MODULE mo_t_datum
!
! Description:
!   Derived type definition and operations for a single observed entity.
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
!  #ifdef USE_SEQUENCE
! V1_5         2009/05/25 Andreas Rhodin
!  changes for TEMP/PILOT level selection
! V1_7         2009/08/24 Andreas Rhodin
!  define constant PRC_ADD_1  (for artificial data generation)
! V1_8         2009/12/09 Andreas Rhodin
!  new component t_datum% mon_pos, mon_rec (position in monitoring file)
! V1_9         2010/04/20 Andreas Rhodin
!  cleanup: dont pass "stat_mnem" to other modules
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Harald Anlauf
!  Fix QC flag handling; properly merge partial TEMP reports; technical changes
! V1_15        2011/12/06 Andreas Rhodin
!  implement data structures to trace invalid observation operators.
!  replace unused status flag value STAT_USED by STAT_OBS_ONLY.
! V1_19        2012-04-16 Andreas Rhodin
!  read spec_index from fof-file, store on body (for RADAR)
!  revised IASI diagnostics in feedback file
! V1_22        2013-02-13 Andreas Rhodin
!  changes for RADAR and STD operator
! V1_26        2013/06/27 Andreas Rhodin
!  replace STAT_ACTIVE_0 by STAT_NOTACTIVE
! V1_45        2015-12-15 Andreas Rhodin
!  replace pl_width (integer) with plev_width (real)
! V1_46        2016-02-05 Andreas Rhodin
!  temporary components in t_datum, t_obs for COMET (italian met. center)
! V1_50        2017-01-09 Harald Anlauf
!  Preparations for high-resolution BUFR TEMPs
!
! Code Description:
! Language: Fortran 95.
! Software Standards:
!
! Author:
! A.Rhodin  DWD  2003-2008
!==============================================================================

  use mo_kind,        only: sp,i1,i2 ,wp  ! precision kind parameter
  use mo_netcdf_param,only: NF_FILL_FLOAT ! missing value indicator
  use mo_t_use,       only: t_use,       &! status data type
                            use_0         ! default values of type t_use
  use mo_obs_rules,   only: t_set,       &! data type
                            empty_set     ! default values of type t_set
  implicit none

  !================
  ! public entities
  !================
  private
  public :: t_datum      ! data type for a single entity
  public :: inv_datum    ! invalid datum
  public :: set_datum    ! set observed datum from BUFR-message
  public :: set_qbits    ! set qc information from BUFR-message
  public :: cmp_datum    ! compare data
  public :: merge_datum  ! merge data
  public :: count_datum  ! count relative data content
  public :: qbit_conv    ! convention: -1=none 0=old(buggy)1=old(DWD)2=new(WMO)
  public :: print        ! generic print routine
  public :: rvind        ! real missing value indicator
  public :: QC_OK,QC_MISS,QC_NOUSE,QC_BLACK,QC_INCON,QC_CLIM,QC_DOUBLE,&
            QC_DENSE,QC_INTPOL
  public :: SRC_OBS, SRC_CB1, SRC_CB2, SRC_CB3, SRC_CDB, &
            SRC_DOK, SRC_DER, SRC_NON

  !==========
  ! constants
  !==========
  !-------------------------------------------------------------------
  ! invalid value used in observation processing
  ! same as NetCDF NF_FILL_FLOAT
  ! +0._sp: saveguard for single-precision constant compared to double
  !-------------------------------------------------------------------
  real(sp) ,parameter :: rvind = NF_FILL_FLOAT+0._sp ! missing value indicator
  !----------------------------------------------
  ! quality control flags (component t_datum% qc)
  !----------------------------------------------
  integer(i2) ,parameter :: QC_OK     =   0 ! OK
  integer(i2) ,parameter :: QC_MISS   =   1 ! missing value
  integer(i2) ,parameter :: QC_DB     =   2 ! flagged bad by databank q_bits
  integer(i2) ,parameter :: QC_TIME   =   4 ! out of time frame
  integer(i2) ,parameter :: QC_NOUSE  =   8 ! do not use flag set
  integer(i2) ,parameter :: QC_BLACK  =  16 ! blacklisted observation
  integer(i2) ,parameter :: QC_INCON  =  32 ! inconsistent observation
  integer(i2) ,parameter :: QC_CLIM   =  64 ! climatological range error
  integer(i2) ,parameter :: QC_DOUBLE = 128 ! double occurence of observation
  integer(i2) ,parameter :: QC_DENSE  = 256 ! data density to high
  integer(i2) ,parameter :: QC_INTPOL = 512 ! interpolation error
  integer     ,save      :: qbit_conv =   0 ! qbit indicator convention (BUFR)
!  !--------------------------------
!  ! use flag (component t_datum% u)
!  !--------------------------------
!  integer(i1) ,parameter :: XXX_NOT = 0 ! dont use this datum
!  integer(i1) ,parameter :: XXX_USE = 1 ! use datum (other than ASS or MON)
!  integer(i1) ,parameter :: XXX_ASS = 2 ! assimilate datum
!  integer(i1) ,parameter :: XXX_MON = 4 ! monitor datum
  !-------------------------------------
  ! source flag (component t_datum% src)
  !-------------------------------------
  integer(i1) ,parameter :: SRC_OBS =  0 ! observed     (BUFR)
  integer(i1) ,parameter :: SRC_CB1 =  1 ! corrected    (BUFR)
  integer(i1) ,parameter :: SRC_CB2 =  2 ! corrected    (BUFR)
  integer(i1) ,parameter :: SRC_CB3 =  3 ! corrected    (BUFR)
  integer(i1) ,parameter :: SRC_CDB =  4 ! corrected    (data base)
  integer(i1) ,parameter :: SRC_DOK =  8 ! derived (OK) (3D-Var)
  integer(i1) ,parameter :: SRC_DER = 16 ! derived (?)  (3D-Var)
  integer(i1) ,parameter :: SRC_NON = 32 ! not set
  !=================
  ! type definitions
  !=================
  !--------------------------------------
  ! Data type for a single observed datum
  !--------------------------------------
  type t_datum
#ifdef USE_SEQUENCE
    SEQUENCE
#endif
    character(len=8) :: mn         =''       ! mnemonic
    real(sp)         :: o          = rvind   ! observed   value
    real(sp)         :: bg         = rvind   ! background value
    real(sp)         :: ac         = -1._sp  ! accuracy (same dimension as 'o')
    real(sp)         :: bc         =  0._sp  ! bias correction
    real(sp)         :: wqc        =  0._sp  ! VQC weight
    real(sp)         :: eo         =  0._sp  ! observational error
    real(sp)         :: eb         =  0._sp  ! background    error
    real(sp)         :: plev       = -1._sp  ! observation height (Pa)
    real(sp)         :: lat        = rvind   ! latitude  (degree)
    real(sp)         :: lon        = rvind   ! longitude (degree)
    real(sp)         :: obs_par(2) = rvind   ! observation specific parameters
    real(sp)         :: plev_width = -1._sp  ! width of sensitivity (log p)
    integer(i2)      :: qc         = QC_MISS ! quality control flag
    integer(i2)      :: mon_pos    = -1      ! position in monitoring file
    integer(i2)      :: lev_sig    =  0      ! level significance
    integer(i2)      :: lev_typ    = -1      ! level type
    integer(i2)      :: spec_index = -1      ! observation specific index
    integer(i2)      :: secs       =  0      ! time offset (seconds)
    type(t_use)      :: use                  ! status flags
    type(t_set)      :: set                  ! namelist settings
    integer(i1)      :: src        = SRC_NON ! source (observed, derived)
    integer(i1)      :: pcc        = -1      ! per cent confidence
    integer(i1)      :: ch_bl      =  0      ! channel blacklist
    integer(i1)      :: op_na      =  0      ! obs.op. not applicable
    !+++++++++++++++++++++++++++++++++
    ! preliminary components for COMET
    !+++++++++++++++++++++++++++++++++
!c  real(wp),pointer :: w_f        (:)
!c  real(wp),pointer :: w_f_coarse (:)
    !+++++++++++++++++++++++++++++++++
  end type t_datum
  !-------------------------------
  ! constant with default settings
  !-------------------------------
  type (t_datum) ,parameter :: inv_datum =                                    &
        t_datum('', rvind, rvind, -1._sp, 0._sp, 0._sp, 0._sp, 0._sp, -1._sp, &
             rvind, rvind, rvind, -1._sp,                                     &
             QC_MISS, -1, 0, -1, -1, 0, use_0, empty_set, SRC_NON, -1, 0, 0   &
!c          ,NULL(),NULL()                                                    &
                                                                              )
  !===========
  ! interfaces
  !===========
  interface print
    module procedure print_datum
  end interface print

contains
!===============================================================================
!----------------------------------------------
! set components of 't_datum' from BUFR-message
!----------------------------------------------
  subroutine set_datum (datum, value, corme)
  type (t_datum) ,intent(inout) :: datum
  real(sp)       ,intent(in)    :: value
  integer        ,intent(in)    :: corme
  !-------------------
  ! set observed datum
  !-------------------
    integer(i1) :: c
    if (value /= rvind) then
      datum% o   = value
      datum% qc  = iand (datum% qc, not(QC_MISS))
      if (datum% src == SRC_NON) datum% src = SRC_OBS
      c          = int (corme, i1)
      datum% src = ior (datum% src, min(SRC_CB3,c))
    endif
  end subroutine set_datum
!------------------------------------------------------------------------------
  subroutine set_qbits (datum, ivalue)
  type (t_datum) ,intent(inout) :: datum
  integer        ,intent(in)    :: ivalue
  !-----------------
  ! set quality bits
  !-----------------
! DWD table  0 02 201(former version; local definition)
! Quality control indication of following value Code figure 0 Good
! 0 no errors detected
! 1 slightly suspect
! 2 bad
! 3 corrected value
!
! WMO table  0 33 020
! Quality control indication of following value Code figure
! 0 Good
! 1 Inconsistent
! 2 Doubtful
! 3 Wrong
! 4 Not checked
! 5 Has been changed
! 6 Estimated
! 7 Missing value
!
    select case (qbit_conv)
    case (-1)
      !---------------
      ! dont set qbits
      !---------------
    case (0)
      !----------------------
      ! old Version (buggy ?)
      !----------------------
      select case (ivalue)
      case (0)                             ! OK
      case (1)                             ! KO
        datum% qc = ior (datum% qc, QC_DB)
      case (2)                             ! corrected
        datum% src = SRC_CDB
      case default                         ! should not happen
        datum% qc = ior (datum% qc, QC_DB)
      end select
      datum% src = ivalue
    case (1)
      !--------------------------------
      ! old Version (DWD), conservative
      !--------------------------------
      select case (ivalue)
      case (0)                             ! OK
      case (1,2)                           ! KO
        datum% qc  = ior (datum% qc, QC_DB)
      case (3)                             ! corrected
!       datum% src = SRC_CDB
        datum% qc  = ior (datum% qc, QC_DB)
      case default                         ! should not happen
        datum% qc  = ior (datum% qc, QC_DB)
      end select
    case (2)
      !------------------
      ! new Version (WMO)
      !------------------
      select case (ivalue)
      case (0) !  Good
      case (1) !  Inconsistent
        datum% qc = ior (datum% qc, QC_DB)
      case (2) !  Doubtful
        datum% qc = ior (datum% qc, QC_DB)
      case (3) !  Wrong
        datum% qc = ior (datum% qc, QC_DB)
      case (4) !  Not checked
      case (5) !  Has been changed
      case (6) !  Estimated
        datum% qc = ior (datum% qc, QC_DB)
      case (7) !  Missing value
      case default                         ! should not happen
        datum% qc = ior (datum% qc, QC_DB)
      end select
    end select
  end subroutine set_qbits
!------------------------------------------------------------------------------
  subroutine print_datum (d, comment)
  type (t_datum)   ,intent(in)           :: d
  character(len=*) ,intent(in) ,optional :: comment
  !----------------------------------------------
  ! printout of some components of type t_datum:
  !   mn  : mnemonic
  !   o   : observed value
  !   qc  : quality control flag
  !   u   : status (assimilate, monitor)
  !   src : source (observed, derived)
  !-----------------------------------------------
    if (present (comment)) then
      write (6,'(a,f15.5,2i4,1x,a)') d% mn, d% o, d% qc, d% src, comment
    else
      write (6,'(a,f15.5,2i4,1x,a)') d% mn, d% o, d% qc, d% src
    endif
  end subroutine print_datum
!------------------------------------------------------------------------------
  !-------------
  ! compare data
  !-------------
  subroutine cmp_datum (d1, d2, diff, m12, m21)
  type (t_datum) ,intent(in)    :: d1, d2 ! data to compare
  logical        ,intent(inout) :: diff   ! true if data differ
  logical        ,intent(inout) :: m12    ! d1 holds information, but d2 not
  logical        ,intent(inout) :: m21    ! d2 holds information, but d1 not
    logical :: ok1, ok2
    ok1 = (d1% src <= SRC_DER .and. d1% qc == QC_OK .and. d1% o /= rvind)
    ok2 = (d2% src <= SRC_DER .and. d2% qc == QC_OK .and. d2% o /= rvind)
    if (ok1.and.ok2) then
      diff = diff .or. d1% o /= d2% o
    else if (ok1) then
      m12  = .true.
    else if (ok2) then
      m21 = .true.
    endif
  end subroutine cmp_datum
!------------------------------------------------------------------------------
  !-----------
  ! merge data
  !-----------
  subroutine merge_datum (dst, src)
  type (t_datum) ,intent(inout) :: dst ! destination
  type (t_datum) ,intent(in)    :: src ! source
    character (len=8) :: mnem
    logical           :: diff, m12, m21
    integer           :: id, is
    integer(i2)       :: sig
    diff = .false.
    m12  = .false.
    m21  = .false.
    call cmp_datum (dst, src, diff, m12, m21)
    if (diff) then
      id = iand (dst% src, SRC_CB3)
      is = iand (src% src, SRC_CB3)
      if (is>id) then
        dst = src
      else if (is==id) then
        mnem    = dst% mn
        dst     = inv_datum
        dst% mn = mnem
        dst% qc = QC_INCON
      endif
    else
      sig = ior (dst% lev_sig, src% lev_sig)
      if (m21) then
        dst = src
      endif
      dst% lev_sig = sig
    endif
    if (dst% src >= SRC_DER) then
      mnem = dst% mn
      dst  = inv_datum
      dst% mn = mnem
    endif
  end subroutine merge_datum
!------------------------------------------------------------------------------
  subroutine count_datum (d1, d2, nsame, ndiff, only1, only2)
    type(t_datum) ,intent(in)    :: d1, d2 ! data to compare
    integer       ,intent(inout) :: nsame  ! count of equivalent data
    integer       ,intent(inout) :: ndiff  ! count of differing  data
    integer       ,intent(inout) :: only1  ! count of data only in set 1
    integer       ,intent(inout) :: only2  ! count of data only in set 2
    !-------------------
    ! count data content
    !-------------------
    logical :: ok1, ok2
    ok1 = (d1% src <= SRC_DER .and. d1% qc == QC_OK .and. d1% o /= rvind)
    ok2 = (d2% src <= SRC_DER .and. d2% qc == QC_OK .and. d2% o /= rvind)
    if (ok1.and.ok2) then
       if (d1% o == d2% o) then
          nsame = nsame + 1
       else
          ndiff = ndiff + 1
       end if
    else if (ok1) then
       only1 = only1 + 1
    else if (ok2) then
       only2 = only2 + 1
    endif
  end subroutine count_datum
!==============================================================================
end module mo_t_datum
