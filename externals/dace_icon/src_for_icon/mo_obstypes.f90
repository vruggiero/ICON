!
!+ observation type tables: BUFR/ECMWF/DWD; tables, conversions
!
MODULE mo_obstypes
!
! Description:
!   This modules holds tables on:
!     WMO BUFR-type and (ECMWF) subtype
!     CMA obstype and codetype
!     mapping between DWD-database-id, ECMA-types and BUFR-types
!
!   functions to:
!     derive data types from CMA-codetyoe
!     derive data types from DWD-database-id (Datenbankkennzahl dbkz)
!
! Current Maintainer: DWD, Harald Anlauf, Alexander Cress
!    phone: +49 69 8062 4941
!    fax:   +49 69 8062 3721
!    email: harald.anlauf@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_1         2008/11/05 Andreas Rhodin
!  First operational 3D-Var release
! V1_2         2008/12/04 Andreas Rhodin
!  changes for SATEM from IASI
! V1_7         2009/08/24 Andreas Rhodin
!  add DBKZ=10385 (BUOY, new BUFR format)
! V1_8         2009/12/09 Harald Anlauf
!  Add "implicit none"
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Alexander Cress
!  changes for driftsondes, wind profilers
! V1_22        2013-02-13 Alexander Cress
!  changes for for wind profilers
! V1_29        2014/04/02 Andreas Rhodin
!  changes for ZTD processing, MODE-S data
! V1_31        2014-08-21 Andreas Rhodin
!  changes for ECMWF SYNOP (BUFR2)NetCDF input
! V1_35        2014-11-07 Alexander Cress
!  changes for TEMP BUFR reports (A.Cress)
! V1_36        2014-11-13 Gerhard Paul
!  Changes for new SYNOP BUFR reports
! V1_40        2015-02-27 Harald Anlauf
!  changes for new SHIP BUFR reports (A.Cress)
! V1_44        2015-09-30 Harald Anlauf
!  Preparations for Jason-2
! V1_45        2015-12-15 Harald Anlauf
!  Preparations for SARAL/Altika
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Andreas Rhodin  DWD  2004       original source
!                      2004-2007  updated
! Gerhard Paul    DWD  2008       DWD SKY data base
! Harald Anlauf   DWD  2008       mobile TEMP stations
!------------------------------------------------------------------------------
  use mo_wmo_tables, only: ec  => WMO0_ECMWF, &! generating center: ECMWF
                           dwd => WMO0_DWD,   &!                    DWD
                           eum => WMO0_EUMET, &!                    EUMETSAT
                           nca => WMO0_NCAR    !                    NCAR

  implicit none
  !----------------
  ! Public entities
  !----------------
  private
  public :: bufr         ! table on WMO BUFR-type and ECMWF subtype
  public :: t_BUFRtyp    ! derived type for 'bufr' table entry
  public :: cma          ! table on ECMA obstype and codetype
  public :: t_CMAtyp     ! derived type for 'cma' table entry
  public :: map          ! table with corresponding BUFR/CMA/WWD obs. types
  public :: t_obsid      ! derived type for 'map' table entry
  public :: obstype_dbkz ! function to derive BUFR/CMA types from DWD-dbkz
  public :: obstype_code ! function to derive BUFR/DWD types from CMA-codetype
  public :: obstype_bufr ! function to derive CMA/DWD types from BUFR type
  public :: dbkz_bufr    ! function to derive the DWD dbkz from BUFR type

  !----------------------------------------------------------------------------
  !
  ! WMO BUFR types and ECMWF subtypes
  !
  !----------------------------------------------------------------------------

  type t_BUFRtyp
    integer           :: type     ! BUFR type
    integer           :: subtype  ! BUFR subtype
    character(len=32) :: name     ! description
  end type t_BUFRtyp

  type (t_bufrtyp) ,parameter :: bufr (49) =         (/&
    t_bufrtyp(  0,  1,'Land SYNOP                   '),&! Land surface
    t_bufrtyp(  0,  3,'Automatic land SYNOP         '),&
    t_bufrtyp(  0, 14,'GNSS zenith total delay      '),&
!   t_bufrtyp(  0,140,'METAR                        '),&
    t_bufrtyp(  0,193,'METAR                        '),&
    !
    t_bufrtyp(  0, 20,'CMAN                         '),&! coastal stations
    !
    t_bufrtyp(  1,  9,'SHIP abbreviated             '),&! Sea surface
    t_bufrtyp(  1, 11,'SHIP 1                       '),&
    t_bufrtyp(  1, 13,'Automatic SHIP               '),&
    t_bufrtyp(  1, 19,'Reduced SHIP                 '),&
    t_bufrtyp(  1, 21,'DRIBU                        '),&
    t_bufrtyp(  1, 22,'BATHY                        '),&
    t_bufrtyp(  1, 23,'TESAC                        '),&
    !
    t_bufrtyp(  2, 91,'LAND PILOT                   '),&! Upper-air
    t_bufrtyp(  2, 92,'SHIP PILOT                   '),&! soundings
    t_bufrtyp(  2, 95,'Wind profiler (USA)          '),&
    t_bufrtyp(  2, 96,'Europ.&Japanese wind profiler'),&
    t_bufrtyp(  2,101,'Land TEMP                    '),&
    t_bufrtyp(  2,102,'SHIP TEMP                    '),&
    t_bufrtyp(  2,103,'DROP TEMP                    '),&
    t_bufrtyp(  2,106,'Mobile TEMP                  '),&
    !
    t_bufrtyp(  3, 51,'High-resolution TOVS 120km   '),&! Satellite
    t_bufrtyp(  3, 53,'RTOVS                        '),&! soundings
    t_bufrtyp(  3, 54,'TOVS-1B                      '),&
    t_bufrtyp(  3, 55,'ATOVS                        '),&
    t_bufrtyp(  3, 61,'Low-level temperature SATEM  '),&
    t_bufrtyp(  3, 62,'High-level temperature SATEM '),&
    t_bufrtyp(  3, 63,'PWC SATEM                    '),&
    t_bufrtyp(  3, 65,'Merged SATEM                 '),&
    t_bufrtyp(  3, 71,'Low-level temperature TOVS   '),&
    t_bufrtyp(  3, 72,'High-Level Temperature TOVS  '),&
    t_bufrtyp(  3, 73,'PWC TOVS                     '),&
    t_bufrtyp(  3, 75,'Merged TOVS                  '),&
    !
    t_bufrtyp(  4,142,'AIREP                        '),&! AIREP
    t_bufrtyp(  4,143,'COLBA                        '),&
    t_bufrtyp(  4,144,'AMDAR                        '),&
    t_bufrtyp(  4,145,'ACARS                        '),&
    t_bufrtyp(  4,146,'MODES                        '),&
    t_bufrtyp(  4,147,'TAMDAR                       '),&
    !
    t_bufrtyp(  5, 82,'Temperature and wind         '),&! SATOB
    t_bufrtyp(  5, 83,'Wind only                    '),&
    t_bufrtyp(  5, 84,'Surface temperature          '),&
    t_bufrtyp(  5, 85,'Clouds temperature           '),&
    t_bufrtyp(  5, 86,'High-resolution VIS wind     '),&
    t_bufrtyp(  5, 87,'AMV + quality control ?      '),&
    !
    t_bufrtyp(253,164,'PAOB                         '),&! PAOB
    !
    t_bufrtyp( -1, -1,'MODIS, GEOS, AIRS            '),&
    !
    t_bufrtyp( -1, 89,'Geostationary radiances      '),&
    t_bufrtyp( -1, 57,'Radiances                    '),&
    t_bufrtyp( -1,189,'Geostationary radiances      ')/)

  !----------------------------------------------------------------------------
  !
  ! CMA obstype and codetype
  !
  !----------------------------------------------------------------------------

  type t_CMAtyp
    integer            :: obstype  ! CMA observation type
    integer            :: codetype ! CMA code type
    character (len=48) :: name     ! description
  end type t_CMAtyp

  type (t_cmatyp) ,parameter :: cma (62) =                      (/&
    !
    ! SYNOP:
    !
    t_cmatyp(  1, 11,'Manual land station                      '),&
    t_cmatyp(  1, 14,'Automatic land station                   '),&
    t_cmatyp(  1, 15,'Road weather information system (SWIS)   '),&
    t_cmatyp(  1, 16,'Snow observations from SYNOP             '),&
    t_cmatyp(  1, 17,'Additional Land Surface Data / METAR USA '),&
    t_cmatyp(  1, 18,'Car data                                 '),&
    t_cmatyp(  1, 19,'Citizen Weather Station (Netatmo)        '),&
    t_cmatyp(  1, 20,'Coastal marine automated network (CMAN)  '),&
    t_cmatyp(  1, 21,'SHIP                                     '),&
    t_cmatyp(  1, 22,'SHIP abbreviated                         '),&
    t_cmatyp(  1, 23,'SHIP reduced (SHRED)                     '),&
    t_cmatyp(  1, 24,'Automatic SHIP                           '),&
    t_cmatyp(  1,140,'METAR                                    '),&
    !
    ! AIREP:
    !
    t_cmatyp(  2, 41,'CODAR                                    '),&
    t_cmatyp(  2,141,'Aircraft                                 '),&
    t_cmatyp(  2,142,'Simulated                                '),&
    t_cmatyp(  2,144,'AMDAR                                    '),&
    t_cmatyp(  2,145,'ACARS                                    '),&
    t_cmatyp(  2,146,'MODES                                    '),&
    t_cmatyp(  2,147,'TAMDAR                                   '),&
    t_cmatyp(  2,241,'COLBA                                    '),&
    !
    ! SATOB
    !
    t_cmatyp(  3, 88,'SATOB                                    '),&
    t_cmatyp(  3, 89,'High-resolution VIS wind                 '),&
    t_cmatyp(  3, 90,'AMV                                      '),&
    t_cmatyp(  3,188,'SST as DRIBU                             '),&
    !
    ! DRIBU:
    !
    t_cmatyp(  4, 63,'BATHY                                    '),&
    t_cmatyp(  4, 64,'TESAC                                    '),&
    t_cmatyp(  4,160,'ERS as DRIBU                             '),&
    t_cmatyp(  4,165,'DRIBU                                    '),&
    !
    ! TEMP:
    !
    t_cmatyp(  5, 35,'LAND                                     '),&
    t_cmatyp(  5, 36,'SHIP                                     '),&
    t_cmatyp(  5, 37,'Mobile                                   '),&
    t_cmatyp(  5, 39,'Land ROCOB                               '),&
    t_cmatyp(  5, 40,'SHIP ROCOB                               '),&
    t_cmatyp(  5,109,'LAND BUFR TEMP                           '),&
    t_cmatyp(  5,111,'SHIP BUFR TEMP                           '),&
    t_cmatyp(  5,135,'DROP                                     '),&
    t_cmatyp(  5,137,'Simulated                                '),&
    t_cmatyp(  5,230,'DROP BUFR TEMP                           '),&
    t_cmatyp(  5,231,'BUFR TEMP DESCENT                        '),&
    !
    ! PILOT:
    !
    t_cmatyp(  6, 32,'Land                                     '),&
    t_cmatyp(  6, 33,'SHIP                                     '),&
    t_cmatyp(  6, 34,'Wind profilers                           '),&
    !
    ! SATEM:
    !
    t_cmatyp(  7, 86,'GTS SATEM (500km)                        '),&
    t_cmatyp(  7,184,'High-res. simulated DWL TOVS             '),&
    t_cmatyp(  7,185,'High-res. simulated DWL SATEM            '),&
    t_cmatyp(  7,186,'High-res. SATEM (250km)                  '),&
    t_cmatyp(  7,200,'GTS BUFR SATEM 250km                     '),&
    t_cmatyp(  7,201,'GTS BUFR SATEM Clear Radiance            '),&
    t_cmatyp(  7,202,'GTS BUFR retr. profiles/clear radiances  '),&
    t_cmatyp(  7,210,'ATOVS                                    '),&
    t_cmatyp(  7,211,'RTOVS                                    '),&
    t_cmatyp(  7,212,'TOVS                                     '),&
    t_cmatyp(  7,215,'SSMI                                     '),&
    !
    ! PAOB:
    !
    t_cmatyp(  8,180,'PAOB                                     '),&
    !
    ! SCATT:
    !
    t_cmatyp(  9,  8,'Scatterometer 1                          '),&
    t_cmatyp(  9, 99,'Scatterometer 2                          '),&
    t_cmatyp(  9,122,'Scatterometer 3                          '),&
    t_cmatyp(  9,210,'Scatterometer 4                          '),&
    !
    ! GPSRO:
    !
    t_cmatyp( 11,250,'GPS Radio Occultation                    '),&
    !
    ! GPSGB:
    !
    t_cmatyp( 12,110,'Ground-based GPS (zenith delay)          '),&
    t_cmatyp( 12,251,'Ground-based GPS (slant delay)           ')/)

  !----------------------------------------------------------------------------
  !
  ! Mapping : DWD database id
  !         : BUFR type,subtype
  !         : ECMA obstype subtype
  !         : 3DVAR module
  !
  !----------------------------------------------------------------------------

  type t_obsid
    integer           :: dbkz      ! DWD   data-base id (Datenbankkennzahl)
    integer           :: bufrtype  ! WMO   BUFR type
    integer           :: subtype   !       BUFR subtype (depends on center)
    integer           :: centre    ! generating centre
    integer           :: obstype   !       CMA  obstype
    integer           :: codetype  !       CMA  codetype
    integer           :: modtype   ! 3DVAR module-type   0: not used in 3DVAR
    character(len=40) :: name      ! description
  end type t_obsid

  type (t_obsid) ,parameter :: map_1 (120) = (/&
  t_obsid(    0,  0,  1, ec, 1, 11,  1,'surface, SYNOP Sect.1-4 manual+PAST'),&
! t_obsid(    1,  0,140, ec, 1,140,  1,'surface, METAR'),&
  t_obsid(    1,  0,193,dwd, 1,140,  1,'surface, METAR'),&
  t_obsid(    2,  0,  1, ec, 1, 11,  1,'surface, SYNOP Sect.5'),&
  t_obsid(    3,  0, -1,dwd,-1, -1,  0,'climatol.data, SYNOP, DWD only'),&
  t_obsid(    4,  0, -1,dwd,-1, -1,  0,'surface, WEHI (DWD)'),&
  t_obsid(    5,  0,  1, ec, 1, 11,  1,'surface, SYNOP Sect.5, from 1.4.2001'),&
  t_obsid(    9,  0, -1,dwd,-1, -1,  0,'surface, PSEUDO-SYNOP (from TEMP)'),&
  t_obsid(   11,  0, -1,dwd,-1, -1,  0,'agr.met. AGRO'),&
  t_obsid(   12,  0, -1,dwd,-1, -1,  0,'agr.met. PHAEN'),&
  t_obsid(   13,  0, -1,dwd,-1, -1,  0,'agr.met. PHYTOPAT'),&
  t_obsid(   16, 10, -1,dwd,-1, -1,  0,'radioact. RADI'),&
  t_obsid(   17, 10, -1,dwd,-1, -1,  0,'radioact. RADA'),&
  t_obsid(   32,  0, -1,dwd,-1, -1,  0,'verification, SYNOP or SHIP obs-ana'),&
  t_obsid(   33,  0, -1,dwd,-1, -1,  0,'verification, Pseudo-SYNOP - ana'),&
  t_obsid(   64,  0, -1,dwd,-1, -1,  0,'climatol.data, CLIMAT'),&
  t_obsid(   65,  0, -1,dwd,-1, -1,  0,'climatol.data, EILKLIMA'),&
  t_obsid(   66,  0, -1,dwd,-1, -1,  0,'climatol.data, CLIMAT Sect. 1-4'),&
  t_obsid(   67,  0, -1,dwd,-1, -1,  0,'climatol.data, CLIMAT Sect. 5'),&
  t_obsid(   94,  0, 14,dwd,12,110,256,'Zenith Total Delay data'),&
  t_obsid(   95,  0, 14,dwd,12,251,256,'Slant Total Delay data'),&
  t_obsid(  128,  0,  3, ec, 1, 14,  1,'surface, SYNOP Sect 1-3 autom.'),&
  t_obsid(  131,  0,193,dwd, 1, 17,  1,'surface, METAR, USA'),&
  t_obsid(  170,  0,  3,dwd, 1, 15,  1,'surface, SWIS'),&
  t_obsid(  182,  0,255,dwd, 1, 16,  1,'surface, Snow observations'),&
!
  t_obsid(  256,  1,  9, ec, 1, 21,1,'SHIP manuell'),&
  t_obsid(  265,  1, 11, ec, 1, 21,1,'Pseudo-SHIP (from TEMP-SHIP)'),&
  t_obsid(  320,  1, -1,dwd,-1, -1,0,'climatol.data, CLIMAT SHIP'),&
  t_obsid(  384,  1, 13, ec, 1, 24,1,'SHIP, autom.'),&
  t_obsid(  385,  1, 21, ec, 4,165,1,'BUOY'),&
  t_obsid(  386,253,164, ec, 8,180,1,'PAOB'),&
  t_obsid(  400,  1, -1,dwd, 1, -1,0,'verification, DRIBU'),&
!
  t_obsid(  508,  2,211,dwd, 6, 32,2,'upper-air, PILOT Part A geopot'),&
  t_obsid(  509,  2,211,dwd, 6, 32,2,'upper-air, PILOT Part B geopot'),&
  t_obsid(  510,  2,211,dwd, 6, 32,2,'upper-air, PILOT Part C geopot'),&
  t_obsid(  511,  2,211,dwd, 6, 32,2,'upper-air, PILOT Part D geopot'),&
  t_obsid(  512,  2, 91, ec, 6, 32,2,'upper-air, PILOT Part A'),&
  t_obsid(  513,  2, 91, ec, 6, 32,2,'upper-air, PILOT Part B'),&
  t_obsid(  514,  2, 91, ec, 6, 32,2,'upper-air, PILOT Part C'),&
  t_obsid(  515,  2, 91, ec, 6, 32,2,'upper-air, PILOT Part D'),&
  t_obsid(  516,  2,215,dwd, 5, 37,2,'upper-air, TEMP  Part A mobil'),&
  t_obsid(  517,  2,215,dwd, 5, 37,2,'upper-air, TEMP  Part B mobil'),&
  t_obsid(  518,  2,215,dwd, 5, 37,2,'upper-air, TEMP  Part C mobil'),&
  t_obsid(  519,  2,215,dwd, 5, 37,2,'upper-air, TEMP  Part D mobil'),&
  t_obsid(  520,  2,101, ec, 5, 35,2,'upper-air, TEMP  Part A'),&
  t_obsid(  521,  2,101, ec, 5, 35,2,'upper-air, TEMP  Part B'),&
  t_obsid(  522,  2,101, ec, 5, 35,2,'upper-air, TEMP  Part C'),&
  t_obsid(  523,  2,101, ec, 5, 35,2,'upper-air, TEMP  Part D'),&
  t_obsid(  524,  2,215,dwd, 5, 37,2,'upper-air, TEMP  Part A mobil'),&
  t_obsid(  525,  2,215,dwd, 5, 37,2,'upper-air, TEMP  Part B mobil'),&
  t_obsid(  526,  2,215,dwd, 5, 37,2,'upper-air, TEMP  Part C mobil'),&
  t_obsid(  527,  2,215,dwd, 5, 37,2,'upper-air, TEMP  Part D mobil'),&
  t_obsid(10520,  2,101,dwd, 5, 35,2,'upper-air, TEMP  BUFR '),&
  t_obsid(10521,  2,101,dwd, 5, 35,2,'upper-air, TEMP  BUFR reduced'),&
  t_obsid(10526,  2,101,dwd, 5,109,2,'upper-air, TEMP  BUFR high res.'),&
  t_obsid(10527,  2,101,dwd, 5,109,2,'upper-air, TEMP  BUFR high res., reduced'),&
  t_obsid(10574,  2,101,dwd, 5,231,2,'upper-air, TEMP  DESCENT BUFR high res.'),&
  t_obsid(10516,  2,101,dwd, 5, 37,2,'upper-air, TEMP  mobile BUFR'),&
  t_obsid(10517,  2,101,dwd, 5, 37,2,'upper-air, TEMP  mobile BUFR, reduced'),&
  t_obsid(10570,  2,106,dwd, 5,231,2,'upper-air, TEMP  mobile DESCENT BUFR high res.'),&
!
! t_obsid(  528,  4, -1,dwd, 2, 41,2,'upper-air, CODAR'),&
  t_obsid(  528,  4,145, ec, 2,145,2,'upper-air, ACARS other'),&
  t_obsid(  529,  4,144, ec, 2,144,2,'upper-air, AMDAR'),&
! t_obsid(  530,  4,142, ec, 2, -1,2,'upper-air, AIREP'),&
  t_obsid(  530,  4,142,dwd, 2,141,2,'upper-air, AIREP'),&
  t_obsid(  531,  2, -1,dwd,-1, -1,0,'upper-air, TEMP from SATEM'),&
  t_obsid(  532,  4,145, ec, 2,145,2,'upper-air, ACARS Lufthansa'),&
  t_obsid(  533,  4,145, ec, 2,145,2,'upper-air, ACARS USA'),&
  t_obsid(  534,  4,145, ec, 2,145,2,'upper-air, ACARS Europe (Bracknell)'),&
  t_obsid(  535,  4,  4,dwd, 2,145,2,'upper-air, ACARS China'),&
  t_obsid(  538,  4,  4,dwd, 2,145,2,'upper-air, ACARS other'),&
  t_obsid(10532,  4,145,dwd, 2,145,2,'upper-air, ACARS (single level)'),&
  t_obsid(10533,  4,147,dwd, 2,147,2,'upper-air, TAMDAR and AFRIS global)'),&
  t_obsid(10534,  4,146,dwd, 2,146,2,'upper-air, MODES (Eurocontrol Europe)'),&
  t_obsid(  542,  4,146,dwd, 2,146,2,'upper-air, MODES (Eurocontrol Europe)'),&
  t_obsid(  536,  2,  2,dwd, 5, 35,2,'upper-air, merged TEMP, PILOT'),&
!
  t_obsid(  537,  2, -1,dwd,-1, -1,0,'Pseudo TEMPS, Pseudo TEMPS, land'),&
  t_obsid(  544,  2, -1,dwd,-1, -1,0,'verification, TEMP     obs-ana'),&
  t_obsid(  545,  2, -1,dwd,-1, -1,0,'verification, PILOT    obs-ana'),&
  t_obsid(  546,  2, -1,dwd,-1, -1,0,'verification, TEMP     obs-ana'),&
  t_obsid(  547,  4, -1,dwd,-1, -1,0,'verification, Aircraft obs-ana'),&
!
  t_obsid(10553,  1,  7,dwd, 6,132,2,'wind profiler, WINPROF'),&
  t_obsid(  548,  2,  0,dwd, 6,132,0,'wind profiler, USA'),&
  t_obsid(  549,  2,  0,dwd, 6,132,0,'wind profiler, USA'),&
  t_obsid(  550,  2,  0,dwd, 6,132,0,'wind profiler, USA'),&
  t_obsid(  551,  2,  1,dwd, 6,132,0,'wind profiler, USA'),&
!
  t_obsid(  552,  2,  7,dwd, 6,132,0,'wind profiler, USA'),&
  t_obsid(  553,  2,  1,dwd, 6,133,0,'wind profiler, Eu/Linddbg. u,v'),&
  t_obsid(  554,  2,  3,dwd, 6,134,0,'wind profiler, Eu/Linddbg. t,w'),&
  t_obsid(  555,  2,  0,dwd, 6,135,0,'wind profiler, Japan'),&
  t_obsid(  556,  2,  1,dwd, 6,136,0,'wind profiler, Rass Deutschl'),&
  t_obsid(  560,  2, -1,dwd,-1, -1,0,'verification, wind profiler obs-ana'),&
  t_obsid(  600,  6, -1,dwd,-1, -1,0,'radar, OPERA'),&
  t_obsid(  640,  6, -1,dwd, 6,137,0,'radar, wind'),&
  t_obsid(  650,  6, -1,dwd,-1, -1,0,'radar, DX-composit'),&
  t_obsid(  663,  6,  3,dwd, 6,138,0,'radar, radialwind'),&
  t_obsid(  584,  1,  7,dwd, 6,157,2,'SCADA, wind power turbines'),&
!
  t_obsid(  764,  2,212,dwd, 6, 33,2,'upper-air, PILOT SHIP A geopot'),&
  t_obsid(  765,  2,212,dwd, 6, 33,2,'upper-air, PILOT SHIP B geopot'),&
  t_obsid(  766,  2,212,dwd, 6, 33,2,'upper-air, PILOT SHIP C geopot'),&
  t_obsid(  767,  2,212,dwd, 6, 33,2,'upper-air, PILOT SHIP D geopot'),&
!
  t_obsid(  768,  2, 92, ec, 6, 33,2,'upper-air, PILOT SHIP A'),&
  t_obsid(  769,  2, 92, ec, 6, 33,2,'upper-air, PILOT SHIP B'),&
  t_obsid(  770,  2, 92, ec, 6, 33,2,'upper-air, PILOT SHIP C'),&
  t_obsid(  771,  2, 92, ec, 6, 33,2,'upper-air, PILOT SHIP D'),&
  t_obsid(  776,  2,102, ec, 5, 36,2,'upper-air, TEMP  SHIP A'),&
  t_obsid(  777,  2,102, ec, 5, 36,2,'upper-air, TEMP  SHIP B'),&
  t_obsid(  778,  2,102, ec, 5, 36,2,'upper-air, TEMP  SHIP C'),&
  t_obsid(  779,  2,102, ec, 5, 36,2,'upper-air, TEMP  SHIP D'),&
  t_obsid(  780,  2,213,dwd, 5,135,2,'upper-air, TEMP  DROP A'),&
  t_obsid(  781,  2,213,dwd, 5,135,2,'upper-air, TEMP  DROP B'),&
  t_obsid(  782,  2,213,dwd, 5,135,2,'upper-air, TEMP  DROP C'),&
  t_obsid(  783,  2,213,dwd, 5,135,2,'upper-air, TEMP  DROP D'),&
  t_obsid(10776,  2,102,dwd, 5, 36,2,'upper-air, TEMP SHIP BUFR'),&
  t_obsid(10777,  2,102,dwd, 5, 36,2,'upper-air, TEMP SHIP BUFR reduced'),&
  t_obsid(10780,  2,213,dwd, 5,230,2,'upper-air, TEMP DROP BUFR'),&
  t_obsid(10782,  2,102,dwd, 5,111,2,'upper-air, TEMP SHIP BUFR high res.'),&
  t_obsid(10783,  2,102,dwd, 5,111,2,'upper-air, TEMP SHIP BUFR high res.,red.'),&
  t_obsid(10785,  2,102,dwd, 5,231,2,'upper-air, TEMP SHIP BUFR high res.,descending'),&
!
  t_obsid(  792,  2,102, ec, 5, 36,2,'upper-air, merged TEMP/PILOT SHIP'),&
!
  t_obsid(   -1, -1, -1,dwd,-1, -1,0,'?'),&
  t_obsid(  800,  2, -1,dwd,-1, -1,0,'verification, TEMP SHIP obs-ana'),&
  t_obsid(  801, -1, -1,dwd,-1, -1,0,'verification, DROP')/)

type (t_obsid) ,parameter :: map_2 (50) = (/&
  t_obsid( 1664,  3, 61, ec, 7, -1,0,'satellite, SATEM A'),&
  t_obsid( 1666,  3, 62, ec, 7, -1,0,'satellite, SATEM B'),&
  t_obsid( 1672,  5, 83, ec, 3, 88,2,'satellite, SATOB Sect.2'),&
  t_obsid( 1673,  5, 83, ec, 3, 88,2,'satellite, SATOB Sect.3'),&
  t_obsid( 1674, 12, -1,dwd, 9, -1,1,'satellite, SATOB Sect.4'),&
  t_obsid( 1675,  5, -1,dwd, 3, -1,2,'satellite, SATOB Sect.5'),&
  t_obsid( 1677,  5, -1,dwd, 3, -1,2,'satellite, SATOB Sect.6'),&
  t_obsid( 1680,  3, -1,dwd,-1, -1,0,'verification, SATEM obs-ana'),&
  t_obsid( 1681,  5, -1,dwd,-1, -1,0,'verification, SATOB obs-ana'),&
  t_obsid( 1688,  3, -1,dwd,-1, -1,0,'satellite, proviles USA'),&
  t_obsid( 1689,  3, 55, ec,-1,210,4,'satellite, ATOVS HIRS'),&
  t_obsid( 1690,  3, 55, ec,-1,210,4,'satellite, ATOVS AMSU-A'),&
  t_obsid( 1691,  3, 55, ec,-1,210,4,'satellite, ATOVS AMSU-B'),&
  t_obsid( 1692,  3, -1,dwd,-1, -1,0,'satellite, AVHRR'),&
  t_obsid( 1696, 12, -1,dwd, 9, -1,1,'satellite, ERS'),&
  t_obsid( 1694,  3,255,dwd,11,250,8,'gnss radio occultations'),&
  t_obsid( 1694,  3,255,nca,11,250,8,'gnss radio occultations'),&
  t_obsid( 1695,  3,201,dwd,11,250,8,'METOP radio occultations thinned'),&
  t_obsid( 1697, 12, -1,dwd, 9,122,0,'satellite, QUICKSCAT'),&
  t_obsid( 1697,  1, -1,dwd, 9, 64,0,'satellite, QUICKSCATasDRIBU'),&
  t_obsid( 1698, 12, -1,dwd, 9,123,0,'satellite, ASCAT Eumetsat'),&
  t_obsid( 1699, 12, -1,dwd, 9,123,0,'satellite, ASCAT'),&
! t_obsid( 1699, 12,122,dwd,15,305,0,'satellite, ASCAT soil moisture'),&
  t_obsid( 1699,  1, -1,dwd, 9, 64,0,'satellite, ASCATasDRIBU'),&
! t_obsid( 1699, 12,122,dwd, 9, 64,0,'satellite, ASCAT level2'),&
! t_obsid( 1700, 12,255,dwd, 9, 64,0,'satellite, OSCAT level2'),&
  t_obsid( 1700,  1, -1,dwd, 9, 64,0,'satellite, OSCAT level2'),&
  t_obsid( 1770,  1, -1,dwd, 9, 64,0,'satellite, HSCAT level2'),&
  t_obsid( 1701,  1, -1,dwd, 9, 64,0,'satellite, JASON level2'),&
  t_obsid( 1702,  1, -1,dwd, 9, 64,0,'satellite, SARAL level2'),&
  t_obsid( 1780,  1, -1,dwd, 9, 64,0,'satellite, SENTINEL level2'),&
  t_obsid( 1704,  5, 87,dwd, 3, 90,2,'satellite, AMV, EUMETSAT'),&
  t_obsid( 1705,  5, 87,dwd, 3, 90,2,'satellite, AMV, GOES'),&
  t_obsid( 1706,  5,  0,dwd, 3, 90,2,'satellite, AMV, MODIS(ASCII)'),&
  t_obsid( 1707,  5,  1,dwd, 3, 90,2,'satellite, AMV, Chinese'),&
  t_obsid( 1708,  5,  0,dwd, 3, 90,2,'satellite, AMV, Japan'),&
  t_obsid( 1709,  5, 87,dwd, 3, 90,2,'satellite, AMV, NASA polar'),&
  t_obsid( 1710,  5, 87,dwd, 3, 90,2,'satellite, AMV, MODIS'),&
  t_obsid( 1720,  5, 87,dwd, 3, 90,2,'satellite, AMV, SENTINEL-3 A/B'),&
  t_obsid( 1794,  3,223,eum, 7, 86,128,'satellite, tq-retrieval, IASI'),&
  t_obsid( 1815, 23, -1,dwd,18, -1,0,'satellite wind lidar' ),&
! t_obsid(    0,  0,  1, ec, 1, 11,1,'surface, SYNOP Sect.1-4 manual+PAST'),&
  t_obsid(10000,  0,  1, ec, 1, 11,1,'surface, SYNOP manual new BUFR format'),&
! t_obsid(  128,  0,  3, ec, 1, 14,1,'surface, SYNOP Sect 1-3 autom.'),&
  t_obsid(10128,  0,  3, ec, 1, 14,1,'surface, SYNOP autom. new BUFR format'),&
  t_obsid(10015,  0,  1, ec, 1, 11,1,'surface, SYNOP manual, BUFR, WIGOS ID'),&
  t_obsid(10143,  0,  3, ec, 1, 14,1,'surface, SYNOP autom., BUFR, WIGOS ID'),&
! t_obsid(    5,  0,  1, ec, 1, 11,1,'surface, SYNOP Sect.5, from 1.4.2001'),&
  t_obsid(10005,  0,  1, ec, 1, 11,1,'surface, SYNOP Sect.5 new BUFR format'),&
  t_obsid(10158,  0, 20,dwd, 1, 20,1,'surface, CMAN coastal stations (USA) autom.'),&
  t_obsid(10170,  0,  3,dwd, 1, 15,1,'surface, SWIS'),&
  t_obsid(10256,  1,  9, ec, 1, 21,1,'SHIP, new BUFR format'),&
  t_obsid(10384,  1, 13, ec, 1, 24,1,'SHIP, new BUFR autom.'),&
  t_obsid(10385,  1, 21, ec, 4,165,1,'BUOY, new BUFR format'),&
  t_obsid(10600,  1,  7,dwd, 6,137,2,'radar wind, OPERA'    ),&
  t_obsid(99998, 22, -1,dwd,15, -1,0,'satellite radar'      )/) ! +++ temporary +++
!

!  type (t_obsid) ,parameter :: map (119) = (/map_1, map_2/)
  type (t_obsid) ,parameter :: map(size (map_1)+size (map_2)) = &
       (/map_1, map_2/)


!01031Code figure generating centre
!01031    00      WMO Secretariat
!01031 01 - 09    WMCs
!01031    01      Melbourne
!01031    07      US National Weather Service, National Centres for Environmental Prediction (NCEP)
!01031    08      US National Weather Service Telecommunications Gateway (NWSTG)
!01031 10 - 25    Centres in Region I
!01031    10      Cairo
!01031    12      Dakar
!01031    14      Nairobi
!01031    18      Tunis-Capaplnca (RSMC)
!01031    20      Las Palmas
!01031    21      Algiers
!01031    24      Pretoria
!01031    25      La R\202union
!01031 26 - 40    Centres in Region II
!01031    34      Tokyo (RSMC), Japan Meteorological Agency
!01031 41 - 50    Centres in Region III
!01031    46      Brazilian Space Agency - INPE
!01031 51 - 63    Centres in Region IV
!01031    51      Miami (RSMC/RAFC)
!01031    52      Miami (RSMC), National Hurricane Centre
!01031    53      Montreal (RSMC)
!01031    57      US Air Force - Air Force Global Weather Central
!01031    58      Fleet Numerical Meteorology and Oceanography Center, Monterey, CA, USA
!01031    59      The NOAA Forecast Systems Laboratory, Boulder, CO, USA
!01031 64 - 73    Centres in Region V
!01031 74 - 99    Centres in Region VI
!01031    74      UK Meteorological Office, Bracknell (RSMC)
!01031    78      Offenbach (RSMC)
!01031    85      Toulouse (RSMC)
!01031    97      European Space Agency (ESA)
!01031    98      European Centre for Medium Range Weather Forecasts (ECMWF) (RSMC)
!01031100 - 254   Reserved
!01031   255      Missing value
!01031NOTE:       See common Code table C-1 in the WMO Manuals on
!                 Code, Vol. I.2, Part C/c.

  type(t_obsid), parameter :: unknown = t_obsid (-1,-1,-1,-1,-1,-1,-1,'')

!------------------------------------------------------------------------------
contains
!------------------------------------------------------------------------------
  function obstype_dbkz (dbkz) result (obstype)
  type(t_obsid)             :: obstype
  integer       ,intent(in) :: dbkz
  !------------------------------------------------
  ! function to derive BUFR/CMA types from DWD-dbkz
  !------------------------------------------------
    integer            :: i

    do i = 1, size(map)
      if (      map(i)% centre == ec  &
          .and. map(i)% dbkz   == dbkz) then
        obstype = map(i)
        return
      endif
    end do

    do i = 1, size(map)
      if (      map(i)% centre == dwd &
          .and. map(i)% dbkz   == dbkz) then
        obstype = map(i)
        return
      endif
    end do

    do i = 1, size(map)
      if (      map(i)% centre == eum &
          .and. map(i)% dbkz   == dbkz) then
        obstype = map(i)
        return
      endif
    end do

    obstype       = unknown
    obstype% dbkz = dbkz

  end function obstype_dbkz
!------------------------------------------------------------------------------
  function obstype_code (obstype, code, centre) result (obsinfo)
  type(t_obsid)                 :: obsinfo
  integer ,intent(in)           :: obstype
  integer ,intent(in)           :: code
  integer ,intent(in) ,optional :: centre
  !----------------------------------------------------
  ! function to derive BUFR/DWD types from CMA-codetype
  !----------------------------------------------------
    integer            :: i
    integer            :: cntr
    !--------------------------------
    ! optional: centre;   default: ec
    !--------------------------------
    cntr = -1; if (present (centre)) cntr = centre
               if (cntr < 0)         cntr = ec
    !-----------------
    ! look up in table
    !-----------------
    do i = 1, size(map)
      if (      map(i)% obstype  == obstype &
          .and. map(i)% codetype == code    &
          .and. map(i)% centre   == cntr    ) then
        obsinfo = map(i)
        return
      endif
    end do
    !----------
    ! not found
    !----------
    obsinfo           = unknown
    obsinfo% obstype  = obstype
    obsinfo% codetype = code

  end function obstype_code
!------------------------------------------------------------------------------
  function obstype_bufr (bufrtype, subtype, centre) result (obstype)
  type(t_obsid)                       :: obstype
  integer       ,intent(in)           :: bufrtype, subtype
  integer       ,intent(in) ,optional :: centre
  !------------------------------------------------
  ! function to derive CMA/DWD types from BUFR type
  !------------------------------------------------
    integer :: i
    integer :: ctr

    ctr = ec; if (present (centre)) ctr = centre
    if (subtype > 0) then
      do i = 1, size(map)
        if (      map(i)% centre   == ctr      &
            .and. map(i)% bufrtype == bufrtype &
            .and. map(i)%  subtype ==  subtype ) then
          obstype = map(i)
          return
        endif
      end do
    else
      do i = 1, size(map)
        if (map(i)% bufrtype == bufrtype) then
          obstype          = map(i)
          obstype% subtype = 0
          return
        endif
      end do
    endif

    obstype           = unknown
    obstype% bufrtype = bufrtype
    obstype%  subtype = subtype
    if (present (centre)) obstype% centre = centre

  end function obstype_bufr
!------------------------------------------------------------------------------
  function dbkz_bufr (bufrtype, subtype, centre)
    integer                       :: dbkz_bufr
    integer, intent(in)           :: bufrtype, subtype, centre
    !---------------------------------------
    ! function to derive DBKZ from BUFR type
    !---------------------------------------
    integer :: i

    dbkz_bufr = -1
    do i = 1, size(map)
       if (      map(i)% centre   == centre   &
           .and. map(i)% bufrtype == bufrtype &
           .and. map(i)%  subtype ==  subtype ) then
          dbkz_bufr = map(i)% dbkz
          return
       endif
    end do
  end function dbkz_bufr
!------------------------------------------------------------------------------
end module mo_obstypes
