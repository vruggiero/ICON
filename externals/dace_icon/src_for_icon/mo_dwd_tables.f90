!
!+ DWD Datenbankkennzahlen, constants, mnemonics
!
MODULE mo_dwd_tables
!
! Description:
!   Parameter definitions for DWD table entries (Datenbankkennzahlen).
!   Routine to derive mnemonics and descriptions.
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
!  New entry: DK_SATEM_IASI = 1794 ! Satellitenbeob, SATEM from IASI
! V1_7         2009/08/24 Andreas Rhodin
!  add dbkz=10385 (BUOY, new BUFR format)
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_11        2010/06/16 Harald Anlauf
!  Changes for GPSRO
! V1_13        2011/11/01 Alexander Cress
!  changes for SCATT, GPSRO, ACARS, Windprofilers
! V1_22        2013-02-13 Alexander Cress
!  changes for wind profilers
! V1_29        2014/04/02 Alexander Cress
!  Handle MODE-S data
! V1_35        2014-11-07 Alexander Cress
!  changes for TEMP BUFR reports (A.Cress)
! V1_36        2014-11-13 Gerhard Paul
!  changes for new SYNOP BUFR reports
! V1_40        2015-02-27 Harald Anlauf
!  changes for new SHIP BUFR reports (A.Cress)
! V1_44        2015-09-30 Harald Anlauf
!  Add DWD dbkz for altimeter data (Jason)
! V1_45        2015-12-15 Harald Anlauf
!  Preparations for SARAL/Altika
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Andreas Rhodin  DWD  2003-2007  initial version
! Gerhard Paul    DWD  2008       add SKY DBKZs
! Harald Anlauf   DWD  2008       add dropsondes, mobile TEMP stations
!------------------------------------------------------------------------------
  implicit none

  PUBLIC

  !---------------------
  ! Datenbank-Kennziffer
  !---------------------
  INTEGER, PARAMETER ::   &
    DK_SYNOP        =    0, & ! Bodenmeldungen, SYNOP, Sect. 1-4, manuell+PAST
    DK_METAR        =    1, & ! Bodenmeldungen, METAR
    DK_WEHI         =    4, & ! Bodenmeldungen, METAR, WEHI
    DK_SYNOP_5      =    5, & ! Bodenmeldungen, SYNOP, Section 5
    DK_P_SYNOP      =    9, & ! Bodenmeldungen, PSEUDO-SYNOP (aus TEMP)
    DK_SYNOP_A      =  128, & ! Bodenmeldungen, SYNOP, Sect. 1-3, autom.
    DK_METAR_USA    =  131, & ! Bodenmeldungen, erweiterte METAR USA
    DK_SEEMHMWS_1   =  156, & ! Bodenmeldungen, SYNOP, SEEMHMWS ECMWF
    DK_SEEMHMWS_2   =  166, & ! Bodenmeldungen, SYNOP, SEEMHMWS 10-Min
    DK_SWIS         =  170, & ! Bodenmeldungen, SWIS,  autom.
    DK_SSNOW        =  182, & ! Bodenmeldungen, SNOW
    DK_SHIP         =  256, & ! Bodenmeldungen, SHIP, manuell
    DK_P_SHIP       =  265, & ! Bodenmeldungen, PSEUDO-SHIP (aus TEMP-SHIP)
    DK_SHIP_AUTO    =  384, & ! Bodenmeldungen, SHIP, automatisch
    DK_BUOY         =  385, & ! Bodenmeldungen, BUOY
    DK_PAOB         =  386, & ! Bodenmeldungen, PAOB (Australien)
! new
    DK_PILOT_A_F    =  508, & ! In-Situ Beob.,  PILOT, PART A, geopotentielle Hoehe
    DK_PILOT_B_F    =  509, & ! In-Situ Beob.,  PILOT, PART B, geopotentielle Hoehe
    DK_PILOT_C_F    =  510, & ! In-Situ Beob.,  PILOT, PART C, geopotentielle Hoehe
    DK_PILOT_D_F    =  511, & ! In-Situ Beob.,  PILOT, PART D, geopotentielle Hoehe
! new
    DK_PILOT_A      =  512, & ! In-Situ Beob.,  PILOT, PART A
    DK_PILOT_B      =  513, & ! In-Situ Beob.,  PILOT, PART B
    DK_PILOT_C      =  514, & ! In-Situ Beob.,  PILOT, PART C
    DK_PILOT_D      =  515, & ! In-Situ Beob.,  PILOT, PART D
! old
    DK_TEMP_A_MOB   =  516, & ! In-Situ Beob.,  TEMP,  PART A, Mobil
    DK_TEMP_B_MOB   =  517, & ! In-Situ Beob.,  TEMP,  PART B, Mobil
    DK_TEMP_C_MOB   =  518, & ! In-Situ Beob.,  TEMP,  PART C, Mobil
    DK_TEMP_D_MOB   =  519, & ! In-Situ Beob.,  TEMP,  PART D, Mobil
!
    DK_TEMP_A       =  520, & ! In-Situ Beob.,  TEMP,  PART A
    DK_TEMP_B       =  521, & ! In-Situ Beob.,  TEMP,  PART B
    DK_TEMP_C       =  522, & ! In-Situ Beob.,  TEMP,  PART C
    DK_TEMP_D       =  523, & ! In-Situ Beob.,  TEMP,  PART D
! new
    DK_TEMP_A_MB    =  524, & ! In-Situ Beob.,  TEMP,  PART A, Mobil
    DK_TEMP_B_MB    =  525, & ! In-Situ Beob.,  TEMP,  PART B, Mobil
    DK_TEMP_C_MB    =  526, & ! In-Situ Beob.,  TEMP,  PART C, Mobil
    DK_TEMP_D_MB    =  527, & ! In-Situ Beob.,  TEMP,  PART D, Mobil
    DK_TEMP_MB_HR   = 10516, & ! In-Situ Beob.,  TEMP, Mobil, BUFR
    DK_TEMP_MB_HRR  = 10517, & ! In-Situ Beob.,  TEMP, Mobil, BUFR, reduced
    DK_TEMP_MB_DESC = 10570, & ! In-Situ Beob.,  TEMP, Mobil, descending BUFR
! new
    DK_TEMP_BUFR    = 10520, & ! In-Situ Beob.,  TEMP,  BUFR (whole report)
    DK_TEMP_BUFR_R  = 10521, & ! In-Situ Beob.,  TEMP,  BUFR (report till 100 hPa)
    DK_TEMP_BUFR_HR = 10526, & ! In-Situ Beob.,  TEMP, high res. BUFR (whole report)
    DK_TEMP_BUFR_HRR= 10527, & ! In-Situ Beob.,  TEMP, high res. BUFR (report till 100 hPa)
    DK_TEMP_DESC_HR = 10574, & ! In-Situ Beob.,  TEMP, high res. descending BUFR (whole report)
! new
!   DK_CODAR        =  528, & ! In-Situ Beob.,  CODAR
    DK_ACARS_OTH    =  528, & ! In-Situ Beob.,  ACARS-Daten sonstige
    DK_AMDAR        =  529, & ! In-Situ Beob.,  AMDAR
    DK_AIREP        =  530, & ! In-Situ Beob.,  AIREP
    DK_ACARS_LH     =  532, & ! In-Situ Beob.,  ACARS-Daten Lufthansa
    DK_ACARS_USA    =  533, & ! In-Situ Beob.,  ACARS-Daten USA
    DK_ACARS_EU     =  534, & ! In-Situ Beob.,  ACARS-Daten EUROPA (Bracknell)
    DK_ACARS_SINGLE =10532, & ! In-Situ Beob.,  ACARS-Daten, single level (unified format)
! old
!   DK_ACARS        =  535, & ! In-Situ Beob.,  ACARS-Daten sonstige
! new
    DK_ACARS_CH     =  535, & ! In-Situ Beob.,  ACARS-Daten China
    DK_ACARS        =  538, & ! In-Situ Beob.,  ACARS-Daten sonstige
    DK_TEMP_PILOT   =  536, & ! In-Situ Beob.,  zusammengefuegter TEMP/PILOT
! new
    DK_PSTEMP_L     =  537    ! Pseudotemps ueber Land

  !---------------------
  ! Datenbank-Kennziffer
  !---------------------
  INTEGER, PARAMETER ::   &
! new
    DK_WINDPROF_USA1  =  548, & ! Windprofiler u,v,(w) USA
    DK_WINDPROF_USA2  =  549, & ! Windprofiler u,v,(w) USA
    DK_WINDPROF_USA3  =  550, & ! Windprofiler u,v,(w) USA
    DK_WINDPROF_USA4  =  551, & ! Windprofiler u,v,(w) USA
    DK_WINDPROF_USA5  =  552, & ! Windprofiler u,v,(w) USA
    DK_WINDPROF_EUL   =  553, & ! Windprofiler u,v,(w) Europa . Lindenberg
    DK_WINDPROF_EUL_T =  554, & ! Windprofiler T,  (w) Europa . Lindenberg
    DK_WINDPROF_JAP   =  555, & ! Windprofiler u,v,(w) Japan
    DK_WINDPROF_RASS  =  556, & ! Windprofiler u,v,(w) RASS Deutschl.
    DK_RADAR_VAD      =  600, & ! Radar Vertical profile u,v
    DK_RADAR_RAD      =  663, & ! Radar Vertical profile radialwind component
!
    DK_PILOT_SHIP_A_F =  764, & ! In-Situ Beob.,  PILOT, SHIP, PART A, geoptentielle Hoehe
    DK_PILOT_SHIP_B_F =  765, & ! In-Situ Beob.,  PILOT, SHIP, PART B, geoptentielle Hoehe
    DK_PILOT_SHIP_C_F =  766, & ! In-Situ Beob.,  PILOT, SHIP, PART C, geoptentielle Hoehe
    DK_PILOT_SHIP_D_F =  767, & ! In-Situ Beob.,  PILOT, SHIP, PART D, geoptentielle Hoehe
! new
    DK_PILOT_SHIP_A   =  768, & ! In-Situ Beob.,  PILOT, SHIP, PART A
    DK_PILOT_SHIP_B   =  769, & ! In-Situ Beob.,  PILOT, SHIP, PART B
    DK_PILOT_SHIP_C   =  770, & ! In-Situ Beob.,  PILOT, SHIP, PART C
    DK_PILOT_SHIP_D   =  771, & ! In-Situ Beob.,  PILOT, SHIP, PART D
    DK_TEMP_SHIP_A    =  776, & ! In-Situ Beob.,  TEMP,  SHIP, PART A
    DK_TEMP_SHIP_B    =  777, & ! In-Situ Beob.,  TEMP,  SHIP, PART B
    DK_TEMP_SHIP_C    =  778, & ! In-Situ Beob.,  TEMP,  SHIP, PART C
    DK_TEMP_SHIP_D    =  779, & ! In-Situ Beob.,  TEMP,  SHIP, PART D
! new
    DK_TEMP_SHIP_BUFR = 10776, & ! In-Situ Beob.,  TEMP, SHIP BUFR (whole report)
    DK_TEMP_SHIP_BUFR_R=10777, & ! In-Situ Beob.,  TEMP, SHIP BUFR (report till 100 hPa)
    DK_TEMP_SHIP_HR   = 10782, & ! In-Situ Beob.,  TEMP, SHIP high res. BUFR (whole report)
    DK_TEMP_SHIP_HRR  = 10783, & ! In-Situ Beob.,  TEMP, SHIP high res. BUFR (report till 100 hPa)
    DK_TEMP_SHIP_DESC_HR=10785,& ! In-Situ Beob.,  TEMP, SHIP high res. descending BUFR (whole report)
! new
    DK_TEMP_DROP_A    =  780, & ! In-Situ Beob.,  TEMP,  DROP, PART A
    DK_TEMP_DROP_B    =  781, & ! In-Situ Beob.,  TEMP,  DROP, PART B
    DK_TEMP_DROP_C    =  782, & ! In-Situ Beob.,  TEMP,  DROP, PART C
    DK_TEMP_DROP_D    =  783, & ! In-Situ Beob.,  TEMP,  DROP, PART D
    DK_TEMP_PILOT_S   =  792, & ! In-Situ Beob.,zusammengefuegter TEMP/PILOT SHIP
    DK_PSTEMP_S       =  793, & ! Pseudotemps ueber See
    DK_TEMP_DROP_BUFR =10780, & ! In-Situ Beob.,  TEMP,  DROP BUFR
! new
    DK_SATEM_A        = 1664, & ! Satellitenbeob, SATEM, PART A
    DK_SATEM_C        = 1666, & ! Satellitenbeob, SATEM, PART C
    DK_SATOB_2        = 1672, & ! Satellitenbeob, SATOB, Section 2
    DK_SATOB_3        = 1673, & ! Satellitenbeob, SATOB, Section 3
    DK_GNSS_RO        = 1694, & ! Satellitenbeob, GNSS radio occultations
    DK_METOP_RO       = 1695, & ! Satellitenbeob, METOP radio occultations(thinned out)
    DK_SATEM_IASI     = 1794, & ! Satellitenbeob, SATEM from IASI
! new
    DK_BUOY_NBUFR     =10385, & ! Bodenmeldungen, BUOY, new BUFR format
    DK_MODES_OLD      =  542, & ! MODES data in AMDAR Format (Eurocontrol Europe)
    DK_TAMDAR         =10533, & ! TAMDAR and AFIRS Aircraft data global
    DK_MODES          =10534, & ! MODES data in AMDAR Format (Eurocontrol Europe)
    DK_WINDPROF       =10553, & ! Windprofiler u,v,(w) USA
    DK_SCADA          =  584, & ! Windprofiler u,v,(w) USA
    DK_VAD            =10600, & ! VAD wind profile u,v,(w) Canada
! new
    DK_SYNOP_NBUFR    =10000, & ! Bodenmeldungen, SYNOP, manual, new BUFR format (since 20141111 in DB)
    DK_SYNOP_5_NBUFR  =10005, & ! Bodenmeldungen, SYNOP, Section 5,  BUFR format
    DK_SYNOP_A_NBUFR  =10128, & ! Bodenmeldungen, SYNOP, autom., new BUFR format
    DK_SYNOP_WIGOS    =10015, & ! Bodenmeldungen, SYNOP, manual, WIGOS ID
    DK_SYNOP_A_WIGOS  =10143, & ! Bodenmeldungen, SYNOP, autom., WIGOS ID
    DK_SEEMHMWS_BUFR  =10150, & ! Bodenmeldungen, SYNOP, SEEMHMWS autom., new BUFR format
    DK_CMAN_A_NBUFR   =10158, & ! Bodenmeldungen, CMAN coastal stations (USA), autom., new BUFR format
    DK_SWIS_NBUFR     =10170, & ! Bodenmeldungen, SWIS,  autom., new BUFR format
    DK_SHIP_NBUFR     =10256, & ! Bodenmeldungen, SYNOP, manual, new BUFR format (since 20141111 in DB)
    DK_SHIP_A_NBUFR   =10384, & ! Bodenmeldungen, SYNOP, autom., new BUFR format
!
    DK_GPSGB_ZTD      =   94, & ! GNSS ground based, zenith delay
    DK_GPSGB_STD      =   95    ! GNSS ground based, slant delay

  !---------------------
  ! Datenbank-Kennziffer
  !---------------------
  INTEGER, PARAMETER ::   &
    DK_QUICKSCAT      = 1697, & ! Satellitenbeob, QUICKSCAT als DRIBU
! new
!   DK_QUICKSCAT      = 1697, & ! Satellitenbeob, QUICKSCAT  original(SKY)
    DK_ASCAT_EU       = 1698, & ! Satellitenbeob, ASCAT als DRIBU / original(SKY) level 2
    DK_ASCAT          = 1699, & ! Satellitenbeob, ASCAT als DRIBU / original(SKY) level 2
    DK_OSCAT          = 1700, & ! Satellitenbeob, OSCAT als DRIBU / original(SKY) level 2
    DK_ALTIM_JASON    = 1701, & ! Satellitenbeob, Altimetry (Jason)
    DK_ALTIM_SARAL    = 1702, & ! Satellitenbeob, Altimetry (SARAL)
    DK_HSCAT          = 1770, & ! Satellitenbeob, HSCAT
    DK_ALTIM_SENTINEL = 1780, & ! Satellitenbeob, Altimetry (SRAL)
! new
    DK_AMV_EUMETSAT   = 1704, & ! Satellitenbeob, AMV,   Eumetsat
    DK_AMV_GOES       = 1705, & ! Satellitenbeob, AMV,   GOES (NOAA/NESDIS)
! new
!   usage in global assimilation
!  < 200705        > 200705           > 200805
!    1706            1705               1710
!    ascii           bufr               bufr
!   pegasus          globus             sky
!   DK_AMV_MODIS_A   DK_AMV_MODIS_OLD   DK_AMV_MODIS
!
!   usage in parallel global assimilation
!   DK_AMV_MODIS     mld_file
!
!   result: modis are used (implicit in 1705)
!           but modis are thinned as GOES
    DK_AMV_MODIS_A    = 1706, & ! Satellitenbeob, AMV,   MODIS, ASCII Format archiv
    DK_AMV_MODIS_OLD  = 1705, & ! Satellitenbeob, AMV,   MODIS, BUFR  Format(>200705)
    DK_AMV_MODIS      = 1710    ! Satellitenbeob, AMV,   MODIS (NOAA/NESDIS  >200805)
! new
  INTEGER, PARAMETER ::   &
    DK_AMV_FY_X       = 1707, & ! Satellitenbeob, AMV,   Chinese
    DK_AMV_MTV        = 1708, & ! Satellitenbeob, AMV,   Japan
    DK_AMV_NOAA       = 1709, & ! Satellitenbeob, AMV,   NOAA experimental
    DK_AMV_SENTINEL   = 1720, & ! Satellitenbeob, AMV,   MODIS (EUMETSAT)
    DK_WLIDAR         = 1815    ! Satellitenbeob, HLOS Wind Lidar AEOLUS
contains

  subroutine dbkz_mnem (entry, mnemonic, comment)
  integer           ,intent(in)   :: entry
  character (len=*) ,intent(out)  :: mnemonic
  character (len=*) ,intent(out)  :: comment
  !--------------------------------------
  ! Derive mnemonic and comment from DBKZ
  !--------------------------------------
    select case (entry)
    case (DK_SYNOP)
      mnemonic = 'SYNOP'
      comment  = 'Bodenmeldungen, SYNOP, Sect. 1-4, manuell+PAST'
    case (DK_METAR)
      mnemonic = 'METAR'
      comment  = 'Bodenmeldungen, METAR'
    case (DK_WEHI)
      mnemonic = 'WEHI'
      comment  = 'Bodenmeldungen, METAR, WEHI'
    case (DK_SYNOP_5)
      mnemonic = 'SYNOP_5'
      comment  = 'Bodenmeldungen, SYNOP, Section 5'
    case (DK_P_SYNOP)
      mnemonic = 'P_SYNOP'
      comment  = 'Bodenmeldungen, PSEUDO-SYNOP (aus TEMP)'
    case (DK_ACARS_OTH)
      mnemonic = 'DK_ACARS_OTH'
      comment  = 'In-Situ Beob.,  ACARS-Daten sonstige'
    case (DK_SYNOP_A)
      mnemonic = 'SYNOP_A'
      comment  = 'Bodenmeldungen, SYNOP, Sect. 1-3, autom.'
    case (DK_SEEMHMWS_1)
       mnemonic = 'SEEMHMWS_1'
       comment  = 'Bodenmeldungen, SEEMHMWS ECMWF'
     case (DK_SEEMHMWS_2)
       mnemonic = 'SEEMHMWS_2'
       comment  = 'Bodenmeldungen, SEEMHMWS 10-Min'
    case (DK_METAR_USA)
      mnemonic = 'METAR_USA'
      comment  = 'Bodenmeldungen, METAR, USA '
!
    case (DK_SYNOP_NBUFR)
      mnemonic = 'SYNOP_NBUFR'
      comment  = 'Bodenmeldungen, SYNOP, manual, new BUFR format'
    case (DK_SYNOP_5_NBUFR)
      mnemonic = 'SYNOP_5_NBUFR'
      comment  = 'Bodenmeldungen, SYNOP, Sect.5, new BUFR format'
    case (DK_SYNOP_A_NBUFR)
      mnemonic = 'SYNOP_A_NBUFR'
      comment  = 'Bodenmeldungen, SYNOP, autom., new BUFR format'
    case (DK_SYNOP_WIGOS)
      mnemonic = 'SYNOP_WIGOS'
      comment  = 'Bodenmeldungen, SYNOP, manual, WIGOS ID'
    case (DK_SEEMHMWS_BUFR)
      mnemonic = 'SEEMHMWS_BUFR'
      comment  = 'Bodenmeldungen, SEEMHMWS Bufr format'
    case (DK_SYNOP_A_WIGOS)
      mnemonic = 'SYNOP_A_WIGOS'
      comment  = 'Bodenmeldungen, SYNOP, autom., WIGOS ID'
    case (DK_SWIS)
      mnemonic = 'SWIS'
      comment  = 'Bodenmeldungen, SWIS,  autom.'
    case (DK_SWIS_NBUFR)
      mnemonic = 'SWIS_NBUFR'
      comment  = 'Bodenmeldungen, SWIS,  autom., new BUFR format'
    case (DK_SSNOW)
      mnemonic = 'SSNOW'
      comment  = 'Bodenmeldungen, SNOW'
    case (DK_CMAN_A_NBUFR)
      mnemonic = 'CMAN_A_NBUFR'
      comment  = 'Bodenmeldungen, CMAN coastal stations (USA), autom.'
!
    case (DK_SHIP_NBUFR)
      mnemonic = 'SHIP_NBUFR'
      comment  = 'Bodenmeldungen, SHIP, manual, new BUFR format'
    case (DK_SHIP_A_NBUFR)
      mnemonic = 'SHIP_A_NBUFR'
      comment  = 'Bodenmeldungen, SHIP, autom., new BUFR format'
!
    case (DK_SHIP)
      mnemonic = 'SHIP'
      comment  = 'Bodenmeldungen, SHIP, manuell'
    case (DK_P_SHIP)
      mnemonic = 'P_SHIP'
      comment  = 'Bodenmeldungen, PSEUDO-SHIP (aus TEMP-SHIP)'
    case (DK_SHIP_AUTO)
      mnemonic = 'SHIP_AUTO'
      comment  = 'Bodenmeldungen, SHIP, automatisch'
    case (DK_BUOY)
      mnemonic = 'BUOY'
      comment  = 'Bodenmeldungen, BUOY'
    case (DK_BUOY_NBUFR)
      mnemonic = 'BUOY'
      comment  = 'Bodenmeldungen, BUOY (new BUFR format)'
    case (DK_PAOB)
      mnemonic = 'PAOB'
      comment  = 'Bodenmeldungen, PAOB (Australien)'
!
    case (DK_GPSGB_ZTD)
      mnemonic = 'GNSS ZTD'
      comment  = 'Zenit Delay an GNSS Station'
    case (DK_GPSGB_STD)
      mnemonic = 'GNSS STD'
      comment  = 'Slant Delay an GNSS Station'
!
    case (DK_PILOT_A_F)
      mnemonic = 'PILOT_A_F'
      comment  = 'In-Situ Beob.,  PILOT, PART A geopotentielle Hoehe'
    case (DK_PILOT_B_F)
      mnemonic = 'PILOT_B_F'
      comment  = 'In-Situ Beob.,  PILOT, PART B geopotentielle Hoehe'
    case (DK_PILOT_C_F)
      mnemonic = 'PILOT_C_F'
      comment  = 'In-Situ Beob.,  PILOT, PART C geopotentielle Hoehe'
    case (DK_PILOT_D_F)
      mnemonic = 'PILOT_D_F'
      comment  = 'In-Situ Beob.,  PILOT, PART D geopotentielle Hoehe'
!
    case (DK_PILOT_A)
      mnemonic = 'PILOT_A'
      comment  = 'In-Situ Beob.,  PILOT, PART A'
    case (DK_PILOT_B)
      mnemonic = 'PILOT_B'
      comment  = 'In-Situ Beob.,  PILOT, PART B'
    case (DK_PILOT_C)
      mnemonic = 'PILOT_C'
      comment  = 'In-Situ Beob.,  PILOT, PART C'
    case (DK_PILOT_D)
      mnemonic = 'PILOT_D'
      comment  = 'In-Situ Beob.,  PILOT, PART D'
!
    case (DK_TEMP_A_MOB)
      mnemonic = 'TEMP_A_MOB'
      comment  = 'In-Situ Beob.,  TEMP,  PART A Mobil'
    case (DK_TEMP_B_MOB)
      mnemonic = 'TEMP_B_MOB'
      comment  = 'In-Situ Beob.,  TEMP,  PART B Mobil'
    case (DK_TEMP_C_MOB)
      mnemonic = 'TEMP_C_MOB'
      comment  = 'In-Situ Beob.,  TEMP,  PART C Mobil'
    case (DK_TEMP_D_MOB)
      mnemonic = 'TEMP_D_MOB'
      comment  = 'In-Situ Beob.,  TEMP,  PART D Mobil'
!
    case (DK_TEMP_A)
      mnemonic = 'TEMP_A'
      comment  = 'In-Situ Beob.,  TEMP,  PART A'
    case (DK_TEMP_B)
      mnemonic = 'TEMP_B'
      comment  = 'In-Situ Beob.,  TEMP,  PART B'
    case (DK_TEMP_C)
      mnemonic = 'TEMP_C'
      comment  = 'In-Situ Beob.,  TEMP,  PART C'
    case (DK_TEMP_D)
      mnemonic = 'TEMP_D'
      comment  = 'In-Situ Beob.,  TEMP,  PART D'
!
    case (DK_TEMP_BUFR)
      mnemonic = 'TEMP_BUFR'
      comment  = 'In-Situ Beob.,  TEMP,  BUFR'
    case (DK_TEMP_BUFR_R)
      mnemonic = 'TEMP_BUFR'
      comment  = 'In-Situ Beob.,  TEMP,  BUFR reduced'
    case (DK_TEMP_BUFR_HR)
      mnemonic = 'TEMP_BUFR_HR'
      comment  = 'In-Situ Beob.,  TEMP,  BUFR high res.'
    case (DK_TEMP_BUFR_HRR)
      mnemonic = 'TEMP_BUFR_HR'
      comment  = 'In-Situ Beob.,  TEMP,  BUFR high res., reduced'
    case (DK_TEMP_DESC_HR)
      mnemonic = 'TEMP_DESC_BUFR_HR'
      comment  = 'In-Situ Beob.,  TEMP,  descending BUFR high res.'
!
    case (DK_TEMP_A_MB)
      mnemonic = 'TEMP_A_MB'
      comment  = 'In-Situ Beob.,  TEMP,  PART A Mobil'
    case (DK_TEMP_B_MB)
      mnemonic = 'TEMP_B_MB'
      comment  = 'In-Situ Beob.,  TEMP,  PART B Mobil'
    case (DK_TEMP_C_MB)
      mnemonic = 'TEMP_C_MB'
      comment  = 'In-Situ Beob.,  TEMP,  PART C Mobil'
    case (DK_TEMP_D_MB)
      mnemonic = 'TEMP_D_MB'
      comment  = 'In-Situ Beob.,  TEMP,  PART D Mobil'
    case (DK_TEMP_MB_HR)
      mnemonic = 'TEMP_MB_HR'
      comment  = 'In-Situ Beob.,  TEMP,  Mobil, BUFR'
    case (DK_TEMP_MB_HRR)
      mnemonic = 'TEMP_MB_HRR'
      comment  = 'In-Situ Beob.,  TEMP,  Mobil, BUFR, reduced'
    case (DK_TEMP_MB_DESC)
      mnemonic = 'TEMP_MB_DESC'
      comment  = 'In-Situ Beob.,  TEMP,  Mobil, descending BUFR'
!
    case (DK_AMDAR)
      mnemonic = 'AMDAR'
      comment  = 'In-Situ Beob.,  AMDAR'
    case (DK_AIREP)
      mnemonic = 'AIREP'
      comment  = 'In-Situ Beob.,  AIREP'
    case (DK_ACARS_LH)
      mnemonic = 'ACARS_LH'
      comment  = 'In-Situ Beob.,  ACARS-Daten Lufthansa'
    case (DK_ACARS_USA)
      mnemonic = 'ACARS_USA'
      comment  = 'In-Situ Beob.,  ACARS-Daten USA'
    case (DK_ACARS_EU)
      mnemonic = 'ACARS_EU'
      comment  = 'In-Situ Beob.,  ACARS-Daten EUROPA (Bracknell)'
!
    case (DK_ACARS_CH)
      mnemonic = 'ACARS_CH'
      comment  = 'In-Situ Beob.,  ACARS-Daten CHINA'
    case (DK_ACARS)
      mnemonic = 'ACARS'
      comment  = 'In-Situ Beob.,  ACARS-Daten others'
    case (DK_ACARS_SINGLE)
      mnemonic = 'DK_ACARS_SINGLE'
      comment  = 'In-Situ Beob.,  ACARS-Daten (single level)'
    case (DK_MODES_OLD)
      mnemonic = 'MODES_OLD'
      comment  = 'In-Situ Beob.,  MODES (Eurocontrol Europe)'
    case (DK_TAMDAR)
      mnemonic = 'TAMDAR'
      comment  = 'In-Situ Beob., TAMDAR and AFIRS-Daten (Global)'
    case (DK_MODES)
      mnemonic = 'MODES'
      comment  = 'In-Situ Beob.,  MODES (Eurocontrol Europe)'
!
    case (DK_TEMP_PILOT)
      mnemonic = 'TEMP_PILOT'
      comment  = 'zusammengefuegter TEMP/PILOT'
    case (DK_PSTEMP_L)
      mnemonic = 'PSTEMP_L'
      comment  = 'Pseudotemps ueber Land'
!
    case (DK_PILOT_SHIP_A_F)
      mnemonic = 'PILOT_SHIP_A_F'
      comment  = 'In-Situ Beob.,  PILOT, SHIP, PART A geopotentielle Hoehe'
    case (DK_PILOT_SHIP_B_F)
      mnemonic = 'PILOT_SHIP_B_F'
      comment  = 'In-Situ Beob.,  PILOT, SHIP, PART B geopotentielle Hoehe'
    case (DK_PILOT_SHIP_C_F)
      mnemonic = 'PILOT_SHIP_C_F'
      comment  = 'In-Situ Beob.,  PILOT, SHIP, PART C geopotentielle Hoehe'
    case (DK_PILOT_SHIP_D_F)
      mnemonic = 'PILOT_SHIP_D_F'
      comment  = 'In-Situ Beob.,  PILOT, SHIP, PART D geopotentielle Hoehe'
    case (DK_PILOT_SHIP_A)
      mnemonic = 'PILOT_SHIP_A'
      comment  = 'In-Situ Beob.,  PILOT, SHIP, PART A'
    case (DK_PILOT_SHIP_B)
      mnemonic = 'PILOT_SHIP_B'
      comment  = 'In-Situ Beob.,  PILOT, SHIP, PART B'
    case (DK_PILOT_SHIP_C)
      mnemonic = 'PILOT_SHIP_C'
      comment  = 'In-Situ Beob.,  PILOT, SHIP, PART C'
    case (DK_PILOT_SHIP_D)
      mnemonic = 'PILOT_SHIP_D'
      comment  = 'In-Situ Beob.,  PILOT, SHIP, PART D'
!
    case (DK_TEMP_SHIP_A)
      mnemonic = 'TEMP_SHIP_A'
      comment  = 'In-Situ Beob.,  TEMP,  SHIP, PART A'
    case (DK_TEMP_SHIP_B)
      mnemonic = 'TEMP_SHIP_B'
      comment  = 'In-Situ Beob.,  TEMP,  SHIP, PART B'
    case (DK_TEMP_SHIP_C)
      mnemonic = 'TEMP_SHIP_C'
      comment  = 'In-Situ Beob.,  TEMP,  SHIP, PART C'
    case (DK_TEMP_SHIP_D)
      mnemonic = 'TEMP_SHIP_D'
      comment  = 'In-Situ Beob.,  TEMP,  SHIP, PART D'
    case (DK_TEMP_SHIP_BUFR)
      mnemonic = 'TEMP_SHIP_BUFR'
      comment  = 'In-Situ Beob.,  TEMP,  SHIP, BUFR  '
    case (DK_TEMP_SHIP_BUFR_R)
      mnemonic = 'TEMP_SHIP_BUFR_R'
      comment  = 'In-Situ Beob.,  TEMP,  SHIP, BUFR reduced'
    case (DK_TEMP_SHIP_HR)
      mnemonic = 'TEMP_SHIP_BUFR_HR'
      comment  = 'In-Situ Beob.,  TEMP,  SHIP BUFR high res.'
    case (DK_TEMP_SHIP_HRR)
      mnemonic = 'TEMP_SHIP_BUFR_HR'
      comment  = 'In-Situ Beob.,  TEMP,  SHIP BUFR high res., reduced'
    case (DK_TEMP_SHIP_DESC_HR)
      mnemonic = 'TEMP_SHIP_DESC_HR'
      comment  = 'In-Situ Beob.,  TEMP,  SHIP descending BUFR high res.'
!
    case (DK_TEMP_DROP_A)
      mnemonic = 'TEMP_DROP_A'
      comment  = 'In-Situ Beob.,  TEMP,  DROP, PART A'
    case (DK_TEMP_DROP_B)
      mnemonic = 'TEMP_DROP_B'
      comment  = 'In-Situ Beob.,  TEMP,  DROP, PART B'
    case (DK_TEMP_DROP_C)
      mnemonic = 'TEMP_DROP_C'
      comment  = 'In-Situ Beob.,  TEMP,  DROP, PART C'
    case (DK_TEMP_DROP_D)
      mnemonic = 'TEMP_DROP_D'
      comment  = 'In-Situ Beob.,  TEMP,  DROP, PART D'
    case (DK_TEMP_DROP_BUFR)
      mnemonic = 'TEMP_DROP_BUFR'
      comment  = 'In-Situ Beob.,  TEMP,  DROP, BUFR'
!
    case (DK_TEMP_PILOT_S)
      mnemonic = 'TEMP_PILOT_S'
      comment  = 'zusammengefuegter TEMP/PILOT SHIP'
    case (DK_PSTEMP_S)
      mnemonic = 'PSTEMP_S'
      comment  = 'Pseudotemps ueber See'
!
    case (DK_WINDPROF_USA1)
      mnemonic = 'WINDPROF_USA1'
      comment  = 'Windprofiler u,v,(w) USA1'
    case (DK_WINDPROF_USA2)
      mnemonic = 'WINDPROF_USA2'
      comment  = 'Windprofiler u,v,(w) USA2'
    case (DK_WINDPROF_USA3)
      mnemonic = 'WINDPROF_USA3'
      comment  = 'Windprofiler u,v,(w) USA3'
    case (DK_WINDPROF_USA4)
      mnemonic = 'WINDPROF_USA4'
      comment  = 'Windprofiler u,v,(w) USA4'
    case (DK_WINDPROF_USA5)
      mnemonic = 'WINDPROF_USA5'
      comment  = 'Windprofiler u,v,(w) USA5'
    case (DK_WINDPROF_EUL)
      mnemonic = 'WINDPROF_EUL'
      comment  = 'Windprofiler u,v,(w) Europa . Lindenberg'
    case (DK_WINDPROF_EUL_T)
      mnemonic = 'WINDPROF_EUL_T'
      comment  = 'RASSprofiler T,  (w) Europa . Lindenberg'
    case (DK_WINDPROF_JAP)
      mnemonic = 'WINDPROF_JAP'
      comment  = 'Windprofiler u,v,(w) Japan'
    case (DK_WINDPROF_RASS)
      mnemonic = 'WINDPROF_RASS'
      comment  = 'Windprofiler u,v,(w) RASS Deutschl'
    case (DK_RADAR_VAD)
      mnemonic = 'DK_RADAR_VAD'
      comment  = 'Radar Vertical profile u,v'
    case (DK_RADAR_RAD)
      mnemonic = 'DK_RADAR_RAD'
      comment  = 'Radar Vertical profile radialwind component'
   case (DK_WINDPROF)
      mnemonic = 'WINDPROF'
      comment  = 'Windprofiler u,v new 10553'
   case (DK_SCADA)
      mnemonic = 'SCADA'
      comment  = 'Wind power turbines'
   case (DK_VAD)
      mnemonic = 'VAD PROF'
      comment  = 'Radar Vertical profile u,v new 10600'
!
    case (DK_SATEM_IASI)
      mnemonic = 'SATEM_IASI'
      comment  = 'Satellitenbeob, SATEM, IASI'
    case (DK_SATEM_A)
      mnemonic = 'SATEM_A'
      comment  = 'Satellitenbeob, SATEM, PART A'
    case (DK_SATEM_C)
      mnemonic = 'SATEM_C'
      comment  = 'Satellitenbeob, SATEM, PART C'
    case (DK_SATOB_2)
      mnemonic = 'SATOB_2'
      comment  = 'Satellitenbeob, SATOB, Section 2'
    case (DK_SATOB_3)
      mnemonic = 'SATOB_3'
      comment  = 'Satellitenbeob, SATOB, Section 3'
    case (DK_GNSS_RO)
      mnemonic = 'GNSS_RO'
      comment  = 'Satellitenbeob, GNSS radio occultations (CHAMP/GRACE,COSMIC,TSRX)'
    case (DK_METOP_RO)
      mnemonic = 'METOP_RO'
      comment  = 'Satellitenbeob, GNSS radio occultations (METOP/GRAS thinned)'
!
    case (DK_QUICKSCAT)
      mnemonic = 'QUICKSCAT'
      comment  = 'Satellitenbeob, QUICKSCAT (als DRIBU)'
    case (DK_ASCAT_EU)
      mnemonic = 'ASCAT EU'
      comment  = 'Satellitenbeob, ASCAT (Eumetsat)'
    case (DK_ASCAT)
      mnemonic = 'ASCAT'
      comment  = 'Satellitenbeob, ASCAT (als DRIBU)'
    case (DK_OSCAT)
      mnemonic = 'OSCAT'
      comment  = 'Satellitenbeob, OSCAT/RSCAT/SSCAT'
    case (DK_HSCAT)
      mnemonic = 'HSCAT'
      comment  = 'Satellitenbeob, HSCAT'
    case (DK_ALTIM_JASON)
      mnemonic = 'ALTIM_JASON'
      comment  = 'Satellitenbeob, Altimetry (Jason)'
    case (DK_ALTIM_SARAL)
      mnemonic = 'ALTIM_SARAL'
      comment  = 'Satellitenbeob, Altimetry (SARAL)'
    case (DK_ALTIM_SENTINEL)
      mnemonic = 'ALTIM_SENTINEL'
      comment  = 'Satellitenbeob, Altimetry (SRAL)'
    case (DK_AMV_EUMETSAT)
      mnemonic = 'AMV_EUMETSAT'
      comment  = 'Satellitenbeob, AMV,   Eumetsat'
    case (DK_AMV_GOES)
      mnemonic = 'AMV_GOES'
      comment  = 'Satellitenbeob, AMV,   GOES (NOAA/NESDIS); MODIS(>200705..200805)'
    case (DK_AMV_MODIS)
      mnemonic = 'AMV_MODIS'
      comment  = 'Satellitenbeob, AMV,   MODIS'
!   case (DK_AMV_MODIS_OLD)
!     mnemonic = 'AMV_MODIS_OLD'
!     comment  = 'Satellitenbeob, AMV,   MODIS, BUFR  Format >200705'
    case (DK_AMV_MODIS_A)
      mnemonic = 'AMV_MODIS_A'
      comment  = 'Satellitenbeob, AMV,   MODIS, ASCII Format archiv'
!
    case (DK_AMV_FY_X)
      mnemonic = 'AMV_FY_X'
      comment  = 'Satellitenbeob, AMV,   chinese'
    case (DK_AMV_MTV)
      mnemonic = 'AMV_MTV'
      comment  = 'Satellitenbeob, AMV,   japanese'
!
    case (DK_AMV_NOAA)
      mnemonic = 'AMV_NOAA'
      comment  = 'Satellitenbeob, AMV,   NOAA(polar) experimental'
!
    case (DK_WLIDAR)
      mnemonic = 'WLIDAR'
      comment  = 'Satellitenbeob, HLOS Wind Lidar, Aeolus'
!
    case (DK_AMV_SENTINEL)
      mnemonic = 'SENTINEL'
      comment  = 'Satellitenbeob, AMV, Sentinel-3 A/B'
!
    case default
      mnemonic = ''
      comment  = 'unknown'
    end select

  end subroutine dbkz_mnem

END MODULE mo_dwd_tables
