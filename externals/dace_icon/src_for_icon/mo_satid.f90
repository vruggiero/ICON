!
!+ WMO Satellite Identifier Table C-5
!
MODULE mo_satid
!
! Description:
!   WMO Satellite Identifiers from Common Code Table C-5.
!
! Current Maintainer: DWD, Robin Faulwetter, Alexander Cress, Harald Anlauf
!    phone: +49 69 8062 2746
!    fax:   +49 69 8062 3721
!    email: robin.faulwetter@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_1         2008/11/05 Andreas Rhodin
!  First operational 3D-Var release
! V1_7         2009/08/24 Andreas Rhodin
!  Add entry for NOAA-19
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Harald Anlauf
!  add new satellite ids
! V1_22        2013-02-13 Harald Anlauf
!  add Megha-Tropiques, CALIPSO, CLOUDSAT, SAC_D, KOMPSAT-5
! V1_35        2014-11-07 Robin Faulwetter
!  Enabled processing of chinese satellite radiances
! V1_37        2014-12-23 Alexander Cress
!  define satellite id METOP_12 = 852 (combined METOP 1 + METOP 2 product)
! V1_40        2015-02-27 Harald Anlauf
!  define satellite ids for HSCAT, RSCAT. (A.Cress)
! V1_43        2015-08-19 Harald Anlauf
!  mo_satid: implement satid_longname lookup via comment; add COSMIC-2 mission
! V1_45        2015-12-15 Alexander Cress
!  changes for Himawawari-8/9; Add satid for SARAL
! V1_47        2016-06-06 Harald Anlauf
!  Add WMO satids for AEOLUS, Jason-3, GPM-CORE, HY-2A, Sentinel 3A
!  Enabled processing of SSMIS and AMSR-2
! V1_48        2016-10-06 Harald Anlauf
!  fix satid for RSCAT to ISS (801)
!  Add CMA as processing center for GPSRO, handle FY-3C/GNOS
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Andreas Rhodin  DWD 2005-2008  original code, updates
! Gerhard Paul    DWD 2008       updates
!=================================================

#ifdef __RADIANCE__
  use mo_rad_util,     only: uppercase
#else
  use mo_dace_string,  only: uppercase => toupper
#endif

  implicit none

  !----------------
  ! Public entities
  !----------------
  private
  public :: satid_t          ! satellite identifier Table
  public :: t_satid          ! derived data type
  public :: satid            ! derive satellite identifier from name
  public :: satid_mnem       ! derive satellite identifier from mnemonic
  public :: satid_longname   ! derive satellite identifier from longname (comment)
  public :: satid_bufr2rttov ! derive RTTOV satellite identification
  public :: satname          ! derive satellite name (8char mnemonic) from satid
  public :: mnem             ! derive satellite name (4char mnemonic) from satid

  type t_satid
    integer           :: code
    character(len=8)  :: mnem
    character(len=16) :: comment
    character(len=4)  :: mnem4
    integer           :: rttov_id
    integer           :: rttov_platf
  end type t_satid

  type(t_satid) ,parameter :: satid_t(203) = [ &
    !   001-099             Europe
    t_satid(  1,'ERS_1   ','ERS 1       '    ,''    ,  1,  8),&
    t_satid(  2,'ERS_2   ','ERS 2       '    ,''    ,  2,  8),&
    t_satid(  3,'METOP-1 ','METOP-1     '    ,'MET1',  1, 10),&! Metop-B
    t_satid(  4,'METOP-2 ','METOP-2     '    ,'MET2',  2, 10),&! Metop-A
    t_satid(  5,'METOP-3 ','METOP-3     '    ,'MET3',  3, 10),&! Metop-C
    t_satid(852,'METOP-12','METOP-12    '    ,''    , -1, -1),&! Comb. of Metop-1 to Metop-3
    t_satid( 20,'SPOT_1  ','SPOT 1      '    ,''    , -1, -1),&
    t_satid( 21,'SPOT_2  ','SPOT 2      '    ,''    , -1, -1),&
    t_satid( 22,'SPOT_3  ','SPOT 3      '    ,''    , -1, -1),&
    t_satid( 23,'SPOT_4  ','SPOT 4      '    ,''    , -1, -1),&
    t_satid( 40,'OERSTED ','OERSTED     '    ,'OERS', -1, -1),&
    t_satid( 41,'CHAMP   ','CHAMP       '    ,'CHMP', -1, -1),&
    t_satid( 42,'TerraSAR','TerraSAR-X  '    ,'TSRX', -1, -1),&
    t_satid( 43,'TanDEM-X','TanDEM-X    '    ,'TDMX', -1, -1),&
    t_satid( 44,'PAZ     ','PAZ         '    ,'PAZ_', -1, -1),&
    t_satid( 46,'SMOS    ','SMOS        '    ,''    , -1, -1),&
    t_satid( 48,'AEOLUS  ','AEOLUS      '    ,''    , -1, -1),&
    t_satid( 50,'METEOS_3','METEOSAT 3  '    ,''    ,  3,  3),&
    t_satid( 51,'METEOS_4','METEOSAT 4  '    ,''    ,  4,  3),&
    t_satid( 52,'METEOS_5','METEOSAT 5  '    ,''    ,  5,  3),&
    t_satid( 53,'METEOS_6','METEOSAT 6  '    ,''    ,  6,  3),&
    t_satid( 54,'METEOS_7','METEOSAT 7  '    ,''    ,  7,  3),&
    t_satid( 55,'METEOS_8','METEOSAT 8  '    ,''    ,  1, 12),&
    t_satid( 56,'METEOS_9','METEOSAT 9  '    ,''    ,  2, 12),&
    t_satid( 57,'METEOS10','METEOSAT 10 '    ,''    ,  3, 12),&
    t_satid( 59,'METEOS_2','METEOSAT 2  '    ,''    ,  2,  3),&
    t_satid( 60,'ENVISAT ','ENVISAT     '    ,''    ,  1, 11),&
    t_satid( 61,'SENTI_3A','SENTINEL 3A '    ,''    ,  1, 19),&
    t_satid( 65,'SENTI_3B','SENTINEL 3B '    ,''    ,  1, 19),&
    t_satid( 66,'SENTI_6A','SENTINEL 6A '    ,'SE6A', -1, -1),&
    t_satid( 67,'SENTI_6B','SENTINEL 6B '    ,'SE6B', -1, -1),&
    t_satid( 70,'METEOS11','METEOSAT 11 '    ,''    ,  4, 12),&
    t_satid( 71,'MTG-1   ','MTG-1       '    ,''    ,  1, 32),&
    t_satid( 72,'MTG-2   ','MTG-2       '    ,''    ,  2, 32),&
    t_satid( 73,'MTG-3   ','MTG-3       '    ,''    ,  3, 32),&
    !   100-199             Japan
    t_satid(120,'ADEOS   ','ADEOS       '    ,''    ,  1, 14),&
    t_satid(121,'ADEOS_2 ','ADEOS 2     '    ,''    ,  2, 14),&
    t_satid(122,'GCOM-W1 ','GCOM-W1     '    ,''    ,  1, 29),&
    t_satid(150,'GMS_3   ','GMS 3       '    ,''    ,  3,  5),&
    t_satid(151,'GMS_4   ','GMS 4       '    ,''    ,  4,  5),&
    t_satid(152,'GMS_5   ','GMS 5       '    ,''    ,  5,  5),&
    t_satid(171,'MTSAT-1 ','MTSAT-1     '    ,''    ,  1, 15),&
    t_satid(172,'MTSAT-2 ','MTSAT-2     '    ,''    ,  2, 15),&
    t_satid(173,'HIMAWA-8','HIMAWARI-8  '    ,''    ,  8, 31),&
    t_satid(174,'HIMAWA-9','HIMAWARI-9  '    ,''    ,  9, 31),&
    !   200-299             USA
    t_satid(200,'NOAA-8  ','NOAA 8      '    ,''    ,  8,  1),&
    t_satid(201,'NOAA-9  ','NOAA 9      '    ,''    ,  9,  1),&
    t_satid(202,'NOAA-10 ','NOAA 10     '    ,''    , 10,  1),&
    t_satid(203,'NOAA-11 ','NOAA 11     '    ,''    , 11,  1),&
    t_satid(204,'NOAA-12 ','NOAA 12     '    ,''    , 12,  1),&
    t_satid(205,'NOAA-14 ','NOAA 14     '    ,''    , 14,  1),&
    t_satid(206,'NOAA-15 ','NOAA 15     '    ,''    , 15,  1),&
    t_satid(207,'NOAA-16 ','NOAA 16     '    ,''    , 16,  1),&
    t_satid(208,'NOAA-17 ','NOAA 17     '    ,''    , 17,  1),&
    t_satid(209,'NOAA-18 ','NOAA 18     '    ,''    , 18,  1),&
    t_satid(220,'LANDSA_5','LANDSAT 5   '    ,''    ,  5, 35),&
    t_satid(221,'LANDSA_4','LANDSAT 4   '    ,''    ,  4, 35),&
    t_satid(222,'LANDSA_7','LANDSAT 7   '    ,''    ,  7, 35),&
    t_satid(223,'NOAA-19 ','NOAA 19     '    ,''    , 19,  1),&
    t_satid(224,'NPP     ','NPP         '    ,'NPP_',  0, 17),&
    t_satid(225,'NOAA-20 ','NOAA 20     '    ,''    , 20,  1),&
    t_satid(226,'NOAA-21 ','NOAA 21     '    ,''    , 21,  1),&
    t_satid(240,'DMSP_7  ','DMSP 7      '    ,''    ,  7,  2),&
    t_satid(241,'DMSP_8  ','DMSP 8      '    ,''    ,  8,  2),&
    t_satid(242,'DMSP_9  ','DMSP 9      '    ,''    ,  9,  2),&
    t_satid(243,'DMSP_10 ','DMSP 10     '    ,''    , 10,  2),&
    t_satid(244,'DMSP_11 ','DMSP 11     '    ,''    , 11,  2),&
    t_satid(245,'DMSP_12 ','DMSP 12     '    ,''    , 12,  2),&
    t_satid(246,'DMSP_13 ','DMSP 13     '    ,''    , 13,  2),&
    t_satid(247,'DMSP_14 ','DMSP 14     '    ,''    , 14,  2),&
    t_satid(248,'DMSP_15 ','DMSP 15     '    ,''    , 15,  2),&
    t_satid(250,'GOES_6  ','GOES 6      '    ,''    ,  6,  4),&
    t_satid(251,'GOES_7  ','GOES 7      '    ,''    ,  7,  4),&
    t_satid(252,'GOES_8  ','GOES 8      '    ,''    ,  8,  4),&
    t_satid(253,'GOES_9  ','GOES 9      '    ,''    ,  9,  4),&
    t_satid(254,'GOES_10 ','GOES 10     '    ,''    , 10,  4),&
    t_satid(255,'GOES_11 ','GOES 11     '    ,''    , 11,  4),&
    t_satid(256,'GOES_12 ','GOES 12     '    ,''    , 12,  4),&
    t_satid(257,'GOES_13 ','GOES 13     '    ,''    , 13,  4),&
    t_satid(258,'GOES_14 ','GOES 14     '    ,''    , 14,  4),&
    t_satid(259,'GOES_15 ','GOES 15     '    ,''    , 15,  4),&
    t_satid(260,'JASON-1 ','JASON-1     '    ,''    ,  1, 36),&
    t_satid(261,'JASON-2 ','JASON-2     '    ,''    ,  2, 36),&
    t_satid(262,'JASON-3 ','JASON-3     '    ,''    ,  3, 36),&
    t_satid(265,'CICERO_1','CICERO OP1  '    ,'CIC1', -1, -1),&! GeoOptics
    t_satid(266,'CICERO_2','CICERO OP2  '    ,'CIC2', -1, -1),&! GeoOptics
    t_satid(267,'GNOMES-A','GNOMES-A    '    ,'GNMA', -1, -1),&! PlanetiQ
    t_satid(268,'GNOMES-B','GNOMES-B    '    ,'GNMB', -1, -1),&! PlanetiQ
    t_satid(269,'SPIRE_LM','SPIRE LEMUR '    ,'SPLM', -1, -1),&
    t_satid(270,'GOES_16 ','GOES 16     '    ,''    , 16,  4),&
    t_satid(271,'GOES_17 ','GOES 17     '    ,''    , 17,  4),&
    t_satid(272,'GOES_18 ','GOES 18     '    ,''    , 18,  4),&
    t_satid(273,'GOES_19 ','GOES 19     '    ,''    , 19,  4),&
    t_satid(281,'QUIKSCAT','QUIKSCAT    '    ,''    , -1, -1),&
    t_satid(282,'TRMM    ','TRMM        '    ,''    ,  1,  7),&
    t_satid(283,'CORIOLIS','CORIOLIS    '    ,''    ,  1, 16),&
    t_satid(285,'DMSP17  ','DMSP17      '    ,''    , 17,  2),&
    t_satid(286,'DMSP18  ','DMSP18      '    ,''    , 18,  2),&
    t_satid(287,'DMSP19  ','DMSP19      '    ,''    , 19,  2),&
    t_satid(288,'GPM     ','GPM         '    ,''    ,  1, 37),&
    !   300-399             Russian Federation
    t_satid(310,'GOMS_1  ','GOMS 1      '    ,''    , -1, -1),&
    t_satid(311,'GOMS_2  ','GOMS 2      '    ,''    , -1, -1),&
    t_satid(320,'MTRR2_21','METEOR 2-21 '    ,''    , -1, -1),&
    t_satid(321,'MTRR3-5 ','METEOR 3-5  '    ,''    , -1, -1),&
    t_satid(322,'MTRR3M-1','METEOR 3M-1 '    ,''    , -1, -1),&
    t_satid(323,'MTRR3M-2','METEOR 3M-2 '    ,''    , -1, -1),&
    t_satid(341,'RESURS  ','RESURS 01-4 '    ,''    , -1, -1),&
    !   400-499             India
    t_satid(421,'OSAT_2  ','OCEANSAT-2  '    ,'OCE2', -1, -1),&
    t_satid(422,'SCATSAT1','SCATSAT-1   '    ,''    , -1, -1),&
    t_satid(423,'OSAT_3  ','OCEANSAT-3  '    ,''    , -1, -1),&
    t_satid(430,'INSAT_1B','INSAT_1B    '    ,''    ,  2, 38),&
    t_satid(431,'INSAT_1C','INSAT_1C    '    ,''    ,  3, 38),&
    t_satid(432,'INSAT_1D','INSAT_1D    '    ,''    ,  4, 38),&
    t_satid(440,'MEGHA_T ','Megha-Tropiques' ,'MEGH',  1, 20),&
    t_satid(441,'SARAL   ','SARAL       '    ,''    ,  1, 33),&
    t_satid(450,'INSAT_2A','INSAT_2A    '    ,''    ,  1, 39),&
    t_satid(451,'INSAT_2B','INSAT_2B    '    ,''    ,  2, 39),&
    t_satid(452,'INSAT_2E','INSAT_2E    '    ,''    ,  5, 39),&
    t_satid(470,'INSAT_3A','INSAT_3A    '    ,''    ,  1, 40),&
    t_satid(471,'INSAT_3D','INSAT_3D    '    ,''    ,  4, 40),&
    t_satid(472,'INSAT_3E','INSAT_3E    '    ,''    ,  5, 40),&
    !   500-599             China
    t_satid(500,'FY-1C   ','FY-1C       '    ,''    ,  3, 13),&
    t_satid(501,'FY-1D   ','FY-1D       '    ,''    ,  4, 13),&
    t_satid(502,'HY-2A   ','HY-2A       '    ,''    , -1, -1),&
    t_satid(503,'HY-2B   ','HY-2B       '    ,''    , -1, -1),&
    t_satid(504,'HY-2C   ','HY-2C       '    ,''    , -1, -1),&
    t_satid(505,'HY-2D   ','HY-2D       '    ,''    , -1, -1),&
    t_satid(510,'FY-2    ','FY-2        '    ,''    ,  1,  6),&
    t_satid(512,'FY-2B   ','FY-2B       '    ,''    ,  2,  6),&
    t_satid(513,'FY-2C   ','FY-2C       '    ,''    ,  3,  6),&
    t_satid(514,'FY-2D   ','FY-2D       '    ,''    ,  4,  6),&
    t_satid(515,'FY-2E   ','FY-2E       '    ,''    ,  5,  6),&
    t_satid(516,'FY-2F   ','FY-2F       '    ,''    ,  6,  6),&
    t_satid(517,'FY-2G   ','FY-2G       '    ,''    ,  7,  6),&
    t_satid(520,'FY-3A   ','FY-3A       '    ,''    ,  1, 23),&
    t_satid(521,'FY-3B   ','FY-3B       '    ,''    ,  2, 23),&
    t_satid(522,'FY-3C   ','FY-3C       '    ,'FY3C',  3, 23),&
    t_satid(523,'FY-3D   ','FY-3D       '    ,'FY3D',  4, 23),&
    t_satid(524,'FY-3E   ','FY-3E       '    ,'FY3E',  5, 23),&
    !   600-699             Europe
    !   700-799             USA
    t_satid(700,'TIROS_M ','TIROS M (ITOS 1)',''    , -1, -1),&
    t_satid(701,'NOAA-1  ','NOAA 1      '    ,''    ,  1,  1),&
    t_satid(702,'NOAA-2  ','NOAA 2      '    ,''    ,  2,  1),&
    t_satid(703,'NOAA-3  ','NOAA 3      '    ,''    ,  3,  1),&
    t_satid(704,'NOAA-4  ','NOAA 4      '    ,''    ,  4,  1),&
    t_satid(705,'NOAA-5  ','NOAA 5      '    ,''    ,  5,  1),&
    t_satid(706,'NOAA-6  ','NOAA 6      '    ,''    ,  6,  1),&
    t_satid(707,'NOAA-7  ','NOAA 7      '    ,''    ,  7,  1),&
    t_satid(708,'TIROS_N ','TIROS N     '    ,''    , -1, -1),&
    t_satid(710,'SMS_1   ','GOES (SMS 1)'    ,''    , -1, -1),&
    t_satid(711,'SMS_2   ','GOES (SMS 2)'    ,''    , -1, -1),&
    t_satid(722,'GRACE_A ','GRACE A     '    ,'GRAA', -1, -1),&
    t_satid(723,'GRACE_B ','GRACE B     '    ,'GRAB', -1, -1),&
    t_satid(724,'COSM2_P1','COSMIC-2 P1 '    ,'C2P1', -1, -1),&! COSMIC-2
    t_satid(725,'COSM2_P2','COSMIC-2 P2 '    ,'C2P2', -1, -1),&! Polar component
    t_satid(726,'COSM2_P3','COSMIC-2 P3 '    ,'C2P3', -1, -1),&
    t_satid(727,'COSM2_P4','COSMIC-2 P4 '    ,'C2P4', -1, -1),&
    t_satid(728,'COSM2_P5','COSMIC-2 P5 '    ,'C2P5', -1, -1),&
    t_satid(729,'COSM2_P6','COSMIC-2 P6 '    ,'C2P6', -1, -1),&
    t_satid(731,'GOES_1  ','GOES 1      '    ,''    ,  1,  4),&
    t_satid(732,'GOES_2  ','GOES 2      '    ,''    ,  2,  4),&
    t_satid(733,'GOES_3  ','GOES 3      '    ,''    ,  3,  4),&
    t_satid(734,'GOES_4  ','GOES 4      '    ,''    ,  4,  4),&
    t_satid(735,'GOES_5  ','GOES 5      '    ,''    ,  5,  4),&
    t_satid(740,'COSMIC_1','COSMIC 1    '    ,'CO01', -1, -1),&! Formosat-3/COSMIC-1
    t_satid(741,'COSMIC_2','COSMIC 2    '    ,'CO02', -1, -1),&
    t_satid(742,'COSMIC_3','COSMIC 3    '    ,'CO03', -1, -1),&
    t_satid(743,'COSMIC_4','COSMIC 4    '    ,'CO04', -1, -1),&
    t_satid(744,'COSMIC_5','COSMIC 5    '    ,'CO05', -1, -1),&
    t_satid(745,'COSMIC_6','COSMIC 6    '    ,'CO06', -1, -1),&
    t_satid(750,'COSM2_E1','COSMIC-2 E1 '    ,'C2E1', -1, -1),&! COSMIC-2
    t_satid(751,'COSM2_E2','COSMIC-2 E2 '    ,'C2E2', -1, -1),&! Equatorial component
    t_satid(752,'COSM2_E3','COSMIC-2 E3 '    ,'C2E3', -1, -1),&
    t_satid(753,'COSM2_E4','COSMIC-2 E4 '    ,'C2E4', -1, -1),&
    t_satid(754,'COSM2_E5','COSMIC-2 E5 '    ,'C2E5', -1, -1),&
    t_satid(755,'COSM2_E6','COSMIC-2 E6 '    ,'C2E6', -1, -1),&
    t_satid(763,'NIMBUS_3','NIMBUS 3    '    ,''    ,  3, 30),&
    t_satid(764,'NIMBUS_4','NIMBUS 4    '    ,''    ,  4, 30),&
    t_satid(765,'NIMBUS_5','NIMBUS 5    '    ,''    ,  5, 30),&
    t_satid(766,'NIMBUS_6','NIMBUS 6    '    ,''    ,  6, 30),&
    t_satid(767,'NIMBUS_7','NIMBUS 7    '    ,''    ,  7, 30),&
    t_satid(780,'ERBS    ','ERBS        '    ,''    , -1, -1),&
    t_satid(781,'UARS    ','UARS        '    ,''    , -1, -1),&
    t_satid(782,'EARTH_PR','EARTH PROBE '    ,''    , -1, -1),&
    t_satid(783,'TERRA   ','TERRA       '    ,''    ,  1,  9),&
    t_satid(784,'AQUA    ','AQUA        '    ,''    ,  2,  9),&
    t_satid(785,'AURA    ','AURA        '    ,''    , -1, -1),&
    t_satid(786,'C/NOFS  ','C/NOFS      '    ,'CNFS', -1, -1),&
    t_satid(787,'CALIPSO ','CALIPSO     '    ,''    ,  1, 27),&
    t_satid(788,'CLOUDSAT','CLOUDSAT    '    ,''    , -1, -1),&
    !   800-998             Other satellite operators
    t_satid(800,'SUNSAT  ','SUNSAT      '    ,'SUNS', -1, -1),&
    t_satid(801,'RAPIDSCA','RAPIDSCAT   '    ,'RSCA', -1, -1),&! ISS (Int. Space Station)
    t_satid(802,'CFOSAT  ','CFOSAT      '    ,'SCAT', -1, -1),&
    t_satid(803,'GRACE_C ','GRACE C     '    ,'GRAC', -1, -1),&
    t_satid(804,'GRACE_D ','GRACE D     '    ,'GRAD', -1, -1),&
    t_satid(811,'GEOKO2A ','GEOKOMSAT-2A'    ,'GKOM', -1, -1),&
    t_satid(820,'SAC_C   ','SAC C       '    ,'SACC', -1, -1),&
    t_satid(821,'SAC_D   ','SAC D       '    ,'SACD', -1, -1),&
    t_satid(825,'KOMPSAT5','KOMPSAT-5   '    ,'KOM5', -1, -1),&
    t_satid(854,'LEOGEO  ','LEOGEO      '    ,'LEGE', -1, -1),&! LEOGEO AMVs
    t_satid(856,'SENTINEL','SENTINEL-3  '    ,'SENT', -1, -1),&
    t_satid(990,'HY2A    ','HY2A SCAT   '    ,'HSCA', -1, -1),&
    !  1000-                Other satellites (not in table C-5) or constellations
    t_satid(1001,'YUNYAO ','Yunyao      '    ,'YUNY', -1, -1),&
    t_satid(1002,'TIANMU ','Tianmu      '    ,'TIAN', -1, -1) ]

contains
!------------------------------------------------------------------------------
  elemental function satid (name)
  character(len=*) ,intent(in) :: name
  integer                      :: satid
  !------------------------------------------
  ! derive satellite identifier from mnemonic
  !------------------------------------------
    integer           :: i
    character (len=8) :: name_
    logical           :: l_repl
    if (name == '') then
      satid = -1                                     ! return -1 for empty name
      return
    endif
    satid = 0                                        ! return  0 if not found
    name_ = trim(name)
    if (name_(1:4) == 'MSG-') then
      select case(name_(5:5))
      case('1')
        name_ = 'METEOS_8'
      case('2')
        name_ = 'METEOS_9'
      case('3')
        name_ = 'METEOS10'
      case('4')
        name_ = 'METEOS11'
      end select
    end if
    do i=1,size(satid_t)
      if (name_ /= satid_t(i)% mnem) cycle
      satid = satid_t(i)% code
      exit
    end do
    if (satid == 0) then
      ! Handle confusion of '_' and '-'
      if (index(name_, '-') > 0) then
        l_repl = .true.
        i = index(name_, '-')
        name_(i:i) = '_'
      else if  (index(name_, '_') > 0) then
        l_repl = .true.
        i = index(name_, '_')
        name_(i:i) = '-'
      else
        l_repl = .false.
      end if
      if (l_repl) then
        do i=1,size(satid_t)
          if (name_ /= satid_t(i)% mnem) cycle
          satid = satid_t(i)% code
          exit
        end do
      end if
      if (satid == 0) then
        !
        name_ = trim(uppercase(name))
        if (trim(name_) == trim(name)) return
        do i=1,size(satid_t)
          if (name_ /= satid_t(i)% mnem) cycle
          satid = satid_t(i)% code
          exit
        end do
        if (satid == 0) then
          ! Handle confusion of '_' and '-'
          if (index(name_, '-') > 0) then
            l_repl = .true.
            i = index(name_, '-')
            name_(i:i) = '_'
          else if  (index(name_, '_') > 0) then
            l_repl = .true.
            i = index(name_, '_')
            name_(i:i) = '-'
          else
            l_repl = .false.
          end if
          if (l_repl) then
            do i=1,size(satid_t)
              if (name_ /= satid_t(i)% mnem) cycle
              satid = satid_t(i)% code
              exit
            end do
          end if
        end if
      end if
    end if
  end function satid
!------------------------------------------------------------------------------
  elemental function satid_mnem (mnem)
  character(len=4) ,intent(in) :: mnem
  integer                      :: satid_mnem
  !-------------------------------------------------
  ! derive satellite identifier from 4-char mnemonic
  !-------------------------------------------------
    integer :: i
    satid_mnem = 0
    do i=1,size(satid_t)
      if (mnem /= satid_t(i)% mnem4) cycle
      satid_mnem = satid_t(i)% code
      exit
    end do
    if (mnem == '') satid_mnem = -1
  end function satid_mnem
!------------------------------------------------------------------------------
  elemental function satid_longname (name)
  character(len=*) ,intent(in) :: name
  integer                      :: satid_longname
  !----------------------------------------------------
  ! derive satellite identifier from longname (comment)
  !----------------------------------------------------
    integer :: i
    satid_longname = 0
    do i=1,size(satid_t)
      if (name /= satid_t(i)% comment) cycle
      satid_longname = satid_t(i)% code
      exit
    end do
    if (name == '') satid_longname = -1
  end function satid_longname
!------------------------------------------------------------------------------
  pure subroutine satid_bufr2rttov (satid, satid_rttov, platform)
    !-------------------------------------------
    ! Conversion of Satellite ids: BUFR -> RTTOV
    !-------------------------------------------
    integer, intent(in)            :: satid   ! BUFR  satellite identifier
    integer, intent(out)           :: satid_rttov  ! RTTOV satellite identifier
    integer, intent(out), optional :: platform  ! RTTOV platform  identifier

    integer :: i, platform_
    !---------------------------------------
    ! set default values of return arguments
    !---------------------------------------
    satid_rttov = -1
    platform_   = -1
    !--------------------------------
    ! search for satellite identifier
    !--------------------------------
    do i=1, size(satid_t)
      if (satid /= satid_t(i)% code) cycle
      satid_rttov = satid_t(i)% rttov_id
      platform_   = satid_t(i)% rttov_platf
      exit
    end do

    if (present (platform)) platform = platform_

  end subroutine satid_bufr2rttov
!------------------------------------------------------------------------------
  elemental function satname (satid)
  character(len=8)    :: satname
  integer ,intent(in) :: satid
  !--------------------------------------------
  ! derive satellite name (mnemonic) from satid
  !--------------------------------------------
    integer :: i
    satname = ''
    do i=1,size(satid_t)
      if (satid /= satid_t(i)% code) cycle
      satname = satid_t(i)% mnem
      exit
    end do
  end function satname
!------------------------------------------------------------------------------
  function mnem (satid)
  character(len=4)    :: mnem
  integer ,intent(in) :: satid
  !--------------------------------------------------
  ! derive satellite name (4char mnemonic) from satid
  !--------------------------------------------------
    integer :: i
    mnem = '****'
    do i=1,size(satid_t)
      if (satid /= satid_t(i)% code) cycle
      mnem = satid_t(i)% mnem4
      if (mnem=='') write (mnem,'(i4.4)') satid
      exit
    end do
  end function mnem
!------------------------------------------------------------------------------
end module mo_satid
