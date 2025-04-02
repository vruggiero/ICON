!
!+ WMO Satellite Instruments Identifier Table C-8 + RTTOV numbers
!
! $Id$
!
MODULE mo_instrid
!
! Description:
!   WMO Satellite Instruments Identifier Table C-8
!   and corresponding RTTOV instrument numbers.
!   (only instruments handled within 3D-Var so far)
!
! Current Code Owner: DWD, Andreas Rhodin
!    phone: +49 69 8062 2722
!    fax:   +49 69 8062 3721
!    email: andreas.rhodin@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_13        2011-11-01 Andreas Rhodin
!  WMO table C-8 sat id tables + RTTOV numbers
! V1_22        2013-02-13 Robin Faulwetter
!  Added ATMS and CrIS
! V1_23        2013-03-26 Robin Faulwetter
!  processing of CrIS data: cloud detection; Restructured MNW cloud detection
! V1_25        2013-05-06 Robin Faulwetter
!  Introduce HIRS on Metop-B in table instrid_t
! V1_35        2014-11-07 Robin Faulwetter
!  Enabled processing of chinese satellite radiances
! V1_47        2016-06-06 Robin Faulwetter
!  Enabled processing of SSMIS and AMSR-2, of GMI data
! V1_50        2017-01-09 Robin Faulwetter
!  Added features to rttov12
!
! Code Description:
! Language: Fortran 95.
! Software Standards:
!
!=================================================
implicit none

  !----------------
  ! Public entities
  !----------------
  private
  public :: instrid_t   ! instrument identifier Table
  public :: t_instrid   ! derived data type
  public :: instr_name  ! derive instrument identifier from name
  public :: name_instr  ! derive instrument name (8char mnemonic) from instrument id
  public :: rttov_name  ! derive RTTOV instrument identifier from name
  public :: name_rttov  ! derive instrument name from RTTOV instrument id
  public :: instr_rttov ! derive instrument identifier from RTTOV number
  public :: rttov_instr ! derive RTTOV number from instrument id
  public :: hss_instr   ! Evaluates, whether an instrumen is a hyperspectral sounder
  public :: mw_instr    ! Evaluates, whether an instrumen is a MW sounder
  public :: ir_instr    ! Evaluates, whether an instrumen is an IR sounder
  public :: vis_instr   ! Evaluates, whether an instrumen is an VIS sounder
  !-----------------------------------------------
  ! instrument identifier table entry derived type
  !-----------------------------------------------
  type t_instrid
    integer           :: code    ! WMO Code number
    integer           :: satid   ! WMO satellite id
    integer           :: rttov   ! RTTOV code number
    character(len=8)  :: mnem    ! Instrument short name
    character(len=8)  :: rtmnem  ! RTTOV      short name
    character(len=14) :: agency  ! Agency
    character(len=42) :: name    ! Instrument long name
  end type t_instrid

  !----------------------------
  ! instrument identifier table
  !----------------------------
  type(t_instrid) ,parameter :: instrid_t(28) = (/ &
    t_instrid(210,  -1, 61,'FCI     ','FCI     ','EUMETSAT      ','Flexible combined imager                  '),& ! MTG-I
    t_instrid(212,  -1, 57,'IRS     ','IRS     ','MTG           ','Infrared Sounder                          '),& ! MTG-IRS
    t_instrid(203,  -1, 15,'MHS     ','MHS     ','EUMETSAT      ','Microwave humidity sounder                '),& ! Metop, NOAA
    t_instrid(205,  -1, 20,'MVIRI   ','MVIRI   ','EUMETSAT      ','METEOSAT visible/infrared imager          '),& ! Meteosat first generation
    t_instrid(207,  -1, 21,'SEVIRI  ','SEVIRI  ','EUMETSAT      ','Spinning enhanced visible/infrared imager '),& ! Meteosat second generation
    t_instrid(221,  -1, 16,'IASI    ','IASI    ','CNES/EUMETSAT ','Infrared atmospheric sounding inferometer '),& ! Metop
    t_instrid(297,  -1, 56,'AHI     ','AHI     ','JMA           ','Advanced Himawar imager                   '),& ! Himawari
    t_instrid(478,  -1, 63,'AMSR2   ','AMSR2   ','JAXA          ','Advanced microwave scanning radiometer 2  '),& ! GCOM-W1
    t_instrid(519,  -1, 71,'GMI     ','GMI     ','NASA          ','GPM microwave imager                      '),& ! GPM-core
    t_instrid(570,  -1,  3,'AMSU-A  ','AMSU-A  ','NOAA          ','Advanced microwave sounding unit-A        '),&
    t_instrid(574,  -1,  4,'AMSU-B  ','AMSU-B  ','NOAA          ','Advanced microwave sounding unit-B        '),&
!    t_instrid(620,  -1, 27,'CRIS    ','CRIS    ','NOAA          ','Cross-track infrared sounder              '),& ! NPP
    t_instrid(620,  -1, 28,'CRIS-FSR','CRIS    ','NPP           ','Cross-track infrared sounder/full sp. res.'),& ! NPP
    t_instrid(621,  -1, 19,'ATMS    ','ATMS    ','NOAA          ','Advanced technology microwave sounder     '),& ! NPP
    t_instrid(606, 206,  0,'HIRS/3  ','HIRS    ','NOAA          ','High-resolution infrared sounder/3        '),& ! NOAA-15
    t_instrid(606, 207,  0,'HIRS/3  ','HIRS    ','NOAA          ','High-resolution infrared sounder/3        '),& ! NOAA-16
    t_instrid(606, 208,  0,'HIRS/3  ','HIRS    ','NOAA          ','High-resolution infrared sounder/3        '),& ! NOAA-17
    t_instrid(607, 209,  0,'HIRS/4  ','HIRS    ','NOAA          ','High-resolution infrared sounder/4        '),& ! NOAA-18
    t_instrid(607, 223,  0,'HIRS/4  ','HIRS    ','NOAA          ','High-resolution infrared sounder/4        '),& ! NOAA-19
    t_instrid(607,   4,  0,'HIRS/4  ','HIRS    ','NOAA          ','High-resolution infrared sounder/4        '),& ! Metop-A
    t_instrid(607,   3,  0,'HIRS/4  ','HIRS    ','NOAA          ','High-resolution infrared sounder/4        '),& ! Metop-A
    t_instrid(615,  -1, 22,'GOESIM  ','GOESIM  ','NOAA          ','GOES imager                               '),& ! GOES 12-15
    t_instrid(617,  -1, 44,'ABI     ','ABI     ','NOAA          ','Advanced baseline imager                  '),& ! GOES 16-19
    t_instrid(908,  -1, 10,'SSMI/S  ','SSMI/S  ','NOAA          ','Special sensor microwave imager sounder   '),& ! DMSP
    t_instrid(933,  -1, 42,'IRAS    ','IRAS    ','CMA           ','Infra Red Atmospheric Sounder             '),& ! Fen-Yun
    t_instrid(941,  -1, 34,'SAPHIR  ','SAPHIR  ','CNES          ','Saphir                                    '),& ! Megha-Tropiques
    t_instrid(954,  -1, 72,'MWTS-2  ','MWTS-2  ','CMA           ','Micro-Wave Temperature Sounder-2          '),& ! Fen-Yun
    t_instrid(953,  -1, 73,'MWHS-2  ','MWHS-2  ','CMA           ','Micro-Wave Humidity Sounder-2             '),& ! Megha-Tropiques
    t_instrid(966,  -1,132,'MWTS-3  ','MWTS-3  ','CMA           ','Micro-Wave Temperature Sounder-3          ')/) ! Fen-Yun


contains
!------------------------------------------------------------------------------
  elemental function instr_name (name)
  character(len=*) ,intent(in) :: name       ! mnemonic
  integer                      :: instr_name ! WMO instrument ID
  !---------------------------------------
  ! derive instrument identifier from name
  !---------------------------------------
    integer :: i
    instr_name = -1
    do i=1,size(instrid_t)
      if (name /= instrid_t(i)% mnem) cycle
      instr_name = instrid_t(i)% code
      exit
    end do
    if (name == '') instr_name = -1
  end function instr_name
!------------------------------------------------------------------------------
  elemental function name_instr (instr)
  character(len=8)    :: name_instr          ! mnemonic
  integer ,intent(in) :: instr               ! WMO instrument ID
  !-----------------------------------------------------------
  ! derive instrument name (8char mnemonic) from instrument id
  !-----------------------------------------------------------
    integer :: i
    name_instr = ''
    do i=1,size(instrid_t)
      if (instr /= instrid_t(i)% code) cycle
      name_instr = instrid_t(i)% mnem
      exit
    end do
  end function name_instr
!------------------------------------------------------------------------------
  elemental function rttov_name (name)
  character(len=*) ,intent(in) :: name       ! mnemonic
  integer                      :: rttov_name ! RTTOV instrument ID
  !---------------------------------------------
  ! derive RTTOV instrument identifier from name
  !---------------------------------------------
    integer :: i
    rttov_name = -1
    do i=1,size(instrid_t)
      if (name /= instrid_t(i)% mnem) cycle
      rttov_name = instrid_t(i)% rttov
      exit
    end do
    if (name == '') rttov_name = -1
  end function rttov_name
!------------------------------------------------------------------------------
  elemental function name_rttov (instr, satid)
  character(len=8)    :: name_rttov          ! mnemonic
  integer ,intent(in) :: instr               ! RTTOV instrument ID
  integer ,intent(in) :: satid               ! WMO   satellite  ID
  !------------------------------------------------
  ! derive instrument name from RTTOV instrument id
  !------------------------------------------------
    integer :: i
    name_rttov = ''
    do i=1,size(instrid_t)
      if (satid /= instrid_t(i)% satid .and. -1 < instrid_t(i)% satid) cycle
      if (instr /= instrid_t(i)% rttov) cycle
      name_rttov = instrid_t(i)% mnem
      exit
    end do
  end function name_rttov
!------------------------------------------------------------------------------
  elemental function instr_rttov (instr, satid)
  integer             :: instr_rttov         ! WMO   instrument ID
  integer ,intent(in) :: instr               ! RTTOV instrument ID
  integer ,intent(in) :: satid               ! WMO   satellite  ID
  !-----------------------------------------------
  ! derive instrument identifier from RTTOV number
  !-----------------------------------------------
    integer :: i
    instr_rttov = -1
    do i=1,size(instrid_t)
      if (satid /= instrid_t(i)% satid .and. -1 < instrid_t(i)% satid) cycle
      if (instr /= instrid_t(i)% rttov) cycle
      instr_rttov = instrid_t(i)% code
      exit
    end do
  end function instr_rttov
!------------------------------------------------------------------------------
  elemental function rttov_instr (instr, satid)
  integer                       :: rttov_instr ! RTTOV instrument ID
  integer, intent(in)           :: instr       ! WMO   instrument ID
  integer, intent(in), optional :: satid       ! WMO   satellite  ID
  !---------------------------------------
  ! derive RTTOV number from instrument id
  !---------------------------------------
    integer :: i
    rttov_instr = -1
    do i=1,size(instrid_t)
      if (instr /= instrid_t(i)% code) cycle
      if (present(satid) .and. instrid_t(i)%satid > 0) then
        if (satid /= instrid_t(i)% satid) cycle
      end if
      rttov_instr = instrid_t(i)% rttov
      exit
    end do
  end function rttov_instr
!------------------------------------------------------------------------------
  elemental logical function hss_instr(id, wmo)
    integer, intent(in)           :: id !< ID of the instrument, default: RTTOV-ID.
    logical, intent(in), optional :: wmo !< interpret ID as WMO-ID.
    !----------------------------------------------------------------
    ! Evaluates whether a given instrument is a hyperspectral sounder
    !----------------------------------------------------------------

    integer,  parameter :: hss_instr_rttov(5) = (/ 11, 16, 27, 28, 57/)
    logical :: l_wmo

    if (present(wmo)) then
      l_wmo = wmo
    else
      l_wmo = .false.
    end if

    if (l_wmo) then
      hss_instr = any(hss_instr_rttov(:) == rttov_instr(id))
    else
      hss_instr = any(hss_instr_rttov(:) == id)
    end if

  end function hss_instr

!------------------------------------------------------------------------------
  elemental logical function mw_instr(id, wmo)
    integer, intent(in)           :: id !< ID of the instrument, default: RTTOV-ID.
    logical, intent(in), optional :: wmo !< interpret ID as WMO-ID.
    !----------------------------------------------------------------
    ! Evaluates whether a given instrument is a microwave sounder
    !----------------------------------------------------------------

    integer,  parameter :: mw_instr_rttov(11) = (/ 3, 4, 10, 15, 19, 34, 63, 71, 72, 73, 132/)

    logical :: l_wmo

    if (present(wmo)) then
      l_wmo = wmo
    else
      l_wmo = .false.
    end if

    if (l_wmo) then
      mw_instr = any(mw_instr_rttov(:) == rttov_instr(id))
    else
      mw_instr = any(mw_instr_rttov(:) == id)
    end if

  end function mw_instr
!------------------------------------------------------------------------------
  elemental logical function ir_instr(id, wmo)
    integer, intent(in)           :: id !< ID of the instrument, default: RTTOV-ID.
    logical, intent(in), optional :: wmo !< interpret ID as WMO-ID.
    !----------------------------------------------------------------
    ! Evaluates whether a given instrument is an IR sounder
    !----------------------------------------------------------------

    integer,  parameter :: ir_instr_rttov(12) = (/ 0, 16, 20, 21, 22, 27, 28, 42, 44, 56, 57, 61/)

    logical :: lwmo

    if (present(wmo)) then
      lwmo = wmo
    else
      lwmo = .false.
    end if

    if (lwmo) then
      ir_instr = any(ir_instr_rttov(:) == rttov_instr(id))
    else
      ir_instr = any(ir_instr_rttov(:) == id)
    end if

  end function ir_instr
!------------------------------------------------------------------------------
  elemental logical function vis_instr(id, wmo)
    integer, intent(in)           :: id !< ID of the instrument, default: RTTOV-ID.
    logical, intent(in), optional :: wmo !< interpret ID as WMO-ID.
    !----------------------------------------------------------------
    ! Evaluates whether a given instrument is an IR sounder
    !----------------------------------------------------------------

    integer,  parameter :: vis_instr_rttov(6) = (/ 20, 21, 22, 44, 56, 61 /)

    logical :: lwmo

    if (present(wmo)) then
      lwmo = wmo
    else
      lwmo = .false.
    end if

    if (lwmo) then
      vis_instr = any(vis_instr_rttov(:) == rttov_instr(id))
    else
      vis_instr = any(vis_instr_rttov(:) == id)
    end if

  end function vis_instr
!------------------------------------------------------------------------------
end module mo_instrid
