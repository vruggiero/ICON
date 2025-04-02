!
!+ Utility routines for Slant Total Delay operator
!
Module mo_std_coord
!
! Description:
!   Utility routines for Slant Total Delay operator
!
! Current Code Owner: DWD, Andreas Rhodin
!    phone: +49 69 8062 2722
!    fax:   +49 69 8062 3721
!    email: andreas.rhodin@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_22        2013-02-13 Michael Bender
!  utility routines for Slant Total Delay operator
! V1_23        2013-03-26 Andreas Rhodin
!  cleanup
! V1_26        2013/06/27 Michael Bender
!  Update
! V1_27        2013-11-08 Michael Bender
!  tl/adjoint routines
! V1_43        2015-08-19 Michael Bender
!  Raytracer added to STD operator; GNSS bias correction.
! V1_45        2015-12-15 Michael Bender
!  adjoint code fixed; Make R(dry)/R(vapor) consistent
! V1_46        2016-02-05 Michael Bender
!  some routines moved to mo_physics and mo_algorithms; Geoid2Geopot corrected
! V1_47        2016-06-06 Michael Bender
!  STD operator: old interpolation option removed; bugfix in BiLinearExp3D_ad
! V1_50        2017-01-09 Michael Bender
!  extended namelist STD_OBS; modules merged with COSMO code.
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
!==============================================================================
!
!> @file mo_std_coord.f90
!> Routines used by the  GNSS STD observation operator
!> (@see mo_std_operator.f90)
!>
!==============================================================================

#ifdef __COSMO__
use kind_parameters, only: wp,         &! working precision kind parameter
                           sp,         &! single  precision kind parameter
                           i2,i4        ! integer kind parameters
#else
use mo_kind,       only: wp,         &! working precision kind parameter
                         sp,         &! single  precision kind parameter
                         i2,i4        ! integer kind parameters
#endif

use mo_std_vector, only: line,       &! straight line derived type definition
                         linedef,    &! define straight line
                         linepos,    &!
                         linepos_tl, &!
                         linepos_ad, &!
                         quadeq,     &!
                         quadeq_tl,  &!
                         quadeq_ad    !

implicit none

!================
! Public entities
!================
private

public ::  deg2rad            ! degrees to radian
public ::  rad2deg            ! radian to degrees
public ::  pi05               ! pi/2
public ::  k1, k2, k3         ! constants for Thayer refrac.
public ::  invalid            ! invalid data, for initialization ...
public ::  invalsp            ! invalid data, for initialization ...
public ::  GNSSStation        ! type, GNSS station data
public ::  GNSSheader         ! BUFR/netCDF header information
public ::  Hchar              ! element of GNSSheader
public ::  Hint               ! element of GNSSheader
public ::  Hfloat             ! element of GNSSheader
public ::  SlantTotalDelay    ! type, STD observations
public ::  SlantHeader        ! type, header STD files
public ::  GPSStation         ! type,  GNSS station data
public ::  WGS84Param         ! type WGS84 reference ellips.
public ::  StrucDomes         ! type GNSS domes
public ::  ImportStakoSNXarr  ! read GNSS station coord.
public ::  ImportESTDarr      ! read STD data
public ::  ImportDomes        ! read GNSS station location
public ::  Ellips2Cart        ! ellips. to cartesian
public ::  Cart2Ellips        ! cartesian to ellips.
public ::  EllipsParam        ! derived type
public ::  LocalHorz2Cart     ! local to cartesian
public ::  LocalHorz2CartMat  ! LocalHorz2Cart with predefined matrix
public ::  Azimut2Horz        ! azimuth + elev. to local
public ::  CrossLineEllips    ! points on line and ellips.
!public ::  Geopot2Geoid       ! geopotential to height above geoid
!public ::  Geoid2Geopot       ! height above geoid to geopotential
public ::  NWein              ! Smith & Weintraub refractivity
public ::  NWein_tl           ! Smith & Weintraub refractivity, tl code
public ::  NWein_ad           ! Smith & Weintraub refractivity, adjoint code
public ::  NdWein             ! Smith & Weintraub dry refractivity
public ::  NwWein             ! Smith & Weintraub wet refractivity
public ::  ZHDsaas            ! Saastamoinen ZHD
public ::  BiLinear2D         ! bilinear interpolation between 4 points
public ::  BiLin2D            ! bilinear interpolation, regular grid
public ::  BiLinExp3D         ! interpolation
public ::  BiLinearExp3D      ! interpolation
public ::  BiLinearExp3D_tl   ! interpolation, tangent-linear code
public ::  BiLinearExp3D_ad   ! interpolation, adjoint code
public ::  BiLinearDiff2D     ! bilinear interpolation
public ::  BiLinearDiff2D_tl  ! bilinear interpolation, tangent-linear code
public ::  BiLinearDiff2D_ad  ! bilinear interpolation, adjoint code
public ::  Shepard2D          ! Shepard interpolation
public ::  Shepard2D_tl       ! Shepard interpolation, tangent-linear code
public ::  Shepard2D_ad       ! Shepard interpolation, adjoint code
public ::  ExpInt1D           ! vertical interpolation
public ::  ExpInt1D_tl        ! vertical interpolation, tangent-linear code
public ::  ExpInt1D_ad        ! vertical interpolation, adjoint code
public ::  ExpInt1DGrad       ! first derivatives with respect to "RefPt"
public ::  ExpInt1DGrad2      ! second derivatives with respect to "RefPt"
public ::  simpned            ! num. integration
public ::  UpdateK123         ! Change ki to namelist input
public ::  gmf                ! Global Mapping Function GMF
public ::  JulianDate         ! computes the Julian Date JD
public ::  GregorianDate      ! computes the Gregorian Date
public ::  JD2MJD             ! converts JD to MJD
public ::  MJD2JD             ! converts MJD to JD
public ::  GREG2DOY           ! computes the day od year (DOY)
public ::  DOY2GREG           ! computes the gregorian date from DOY
public ::  IntegPolyCube      ! numerical integration
public ::  IntegPolyCube_tl   ! numerical integration, tangent-linear code
public ::  IntegPolyCube_ad   ! numerical integration, adjoint code
public ::  LayerSearch        ! find vertical grid index k
public ::  SatPos             ! compute GNSS satellite position
public ::  Slant2Ellips       ! transform slant to ellipsoidal system
public ::  Slant2Ellips_tl    ! transform slant to ellipsoidal system, tl code
public ::  Slant2EllipsMat    ! Slant2Ellips with predefined matrix
public ::  Slant2EllipsMatGrad  ! first derivatives of "Slant2EllipsMat"
public ::  Slant2EllipsMatGrad2 ! second derivatives of "Slant2EllipsMat"
public ::  Slant2LocalHorzMat ! transform slant to local horizon system
public ::  GetFreeUnit        ! find the next free file unit
public ::  DistSphere         ! minimum distance between two points on a sphere
!-------------------
! Physical constants
!-------------------
public :: RDRD                ! R(dry)/R(vapor)
public :: EMRDRD              ! 1._wp - RDRD


! Umrechnung Winkel: Grad nach Radian - Winkel(rad) = 2*PI/360 * Winkel(Grad)
real (wp), parameter :: deg2rad = 3.14159265358979323_wp / 180.0_wp
! Umrechnung Winkel: Radian nach Grad - Winkel(Grad) = 360/(2*PI)*Winkel(rad)
real (wp), parameter :: rad2deg = 180.0_wp / 3.14159265358979323_wp

real (wp), parameter :: pi   = 3.14159265358979323_wp
real (wp), parameter :: pi2  = 6.28318530717958623_wp
real (wp), parameter :: pi05 = 1.57079632679489655_wp

! Some small real number
real (wp), parameter :: reps = 1.0E-10_wp

! Molare Masse von trockener Luft Md = 28,9644 g/mol = 0.0289644 kg/mol
real (wp), parameter :: Md = 0.0289644_wp
! Molare Masse des Wasserdampfes Mw = 18,01528 g/mol = 0.01801528 kg/mol
real (wp), parameter :: Mw = 0.01801528_wp

 real(wp), parameter :: RDRD   = Mw/Md          ! = R(dry)/R(vapor) ~ 0.62198
 real(wp), parameter :: EMRDRD = 1._wp - RDRD   !                   ~ 0.37802

!----------------------------------------------------------------------
! Definition einiger Konstanten fuer die Refraktivitaet der Atmosphaere
! nach Bevis, 1994
!----------------------------------------------------------------------
real (wp) :: k1 = 77.60_wp      ! k1 in K/hPa
real (wp) :: k2 = 70.40_wp      ! k2 in K/hPa
real (wp) :: k3 = 3.739E5_wp    ! k3 in K^2/hPa
! k22 = k2 - (Mw/Md)*k2 = 22.1 +- 2.2 K hPa^-1
!real (wp) :: k22 = 22.10_wp  ! K hPa^-1
!---------------------------------------------------------------------

! Grosse Halbachse eines GPS-Satelliten-Orbits (TLE) in Meter
! (Bezogen auf den Erdmittelpunkt)
real (wp), parameter :: GPSradius = 26561000.0_wp ! m

! Invalid data, for initialization ...
real(wp), parameter :: invalid = -999.0_wp
real(sp), parameter :: invalsp = -999.0_sp

!---------------------------------------------------------------------
! structure GNSSStation
!---------------------------------------------------------------------
!>
!> @brief Description of GNSS stations: Name, ID, coordinates, ...
!>
!> The derived type GNSSStation provides information about a GNSS station.
!> Three different sets of coordinates are allowed but not all of them
!> might be available.  When the station data are read first, the coordinates
!> available in this data set are written. These coordinates can later be
!> transformed to the other reference systems. \n
!> The station ID is taken from the provider and might not be available or
!> might not be unique if data from different providers are combined.
!>
!> @todo
!> Station ID: remove ID or define unique IDs within the DWD ?? \n
!> \n
!> StartMJD and EndMJD are meaningfull only if the station data are
!> read from a separate station list which provides coordinates for
!> a long period. In case of BUFR data each record should provide
!> the correct coordinates and no reference date is required.
!> Remove StartMJD and EndMJD ??
!>
!> <table>
!> <tr><th>variable  <th>description
!> <tr><td><b> Name </b>  <td>
!>    full name of GNSS station (most of the time not available)
!> <tr><td><b> SName </b>  <td>
!>    short name of GNSS station, should be an unique 4 character
!>    identifier of the station
!> <tr><td><b> ID </b>  <td>
!>    local station ID, defined while reading data, unique only
!             during program execution
!> <tr><td><b> CenterID </b>  <td>
!>    station ID defined by processing center and distributed
!             with the observations
!> <tr><td><b> CoordEll </b>  <td>
!>    ellipsoidal coordinates of the station: \n
!>             CoordEll(1) - longitude, radian \n
!>             CoordEll(2) - latitude, radian \n
!>             CoordEll(3) - height above ellipsoid, meter
!> <tr><td><b> CoordGeo </b>  <td>
!>    geographical coordinates of the station: \n
!>             CoordGeo(1) - longitude, radian \n
!>             CoordGeo(2) - latitude, radian  \n
!>             CoordGeo(3) - height above geoid = mean sea level, meter  \n
!>    The only difference between ellipsoidal coordinates and geographical
!>    coordinates is the height, latitude and longitude are identical. The
!>    height differs by the geoid undulation.
!> <tr><td><b> CoordCart </b>  <td>
!>    cartesian coordinates in a Earth centered Earth fixed
!>             frame of reference:   \n
!>             CoordCart(1) - X, m  \n
!>             CoordCart(2) - Y, m  \n
!>             CoordCart(3) - Z, m
!> <tr><td><b> GeoidCorr </b>  <td>
!>    geoid correction, meter, height difference between ellipsoid and
!>             geoid at the station: \n
!>             GeoidCorr =  CoordGeo(3) - CoordEll(3) \n
!>                       =  height above geoid - height above ellipsoid
!> <tr><td><b> StartMJD </b>  <td>
!>     station coordinates are valid from StartMJD (MJD - modified Julian Date)
!> <tr><td><b> EndMJD </b>  <td>
!>     station coordinates are valid until EndMJD
!> <tr><td><b> DataCategory </b>  <td>
!>     field "section1_data_category" in BUFR file
!> <tr><td><b> DataSubCategory  </b>  <td>
!>     field "section1_int_data_sub_category" in BUFR file
!> <tr><td><b> Center </b>  <td>
!>     field "section1_centre" in BUFR file
!> <tr><td><b> SubCenter  </b>  <td>
!>     field "section1_subcentre" in BUFR file
!> </table>
!>
!> The variables <b> DataCategory </b>, <b> DataSubCategory  </b>,
!> <b> Center </b> and <b> SubCenter  </b>
!> are not used by the program but need to
!> be passed from the BUFR file to the feedback file.
!> It is assumed that separate data sets are defined for each station
!> processd by a given processing center and for each data product.
!---------------------------------------------------------------------
type GNSSStation
   character (len=20)        :: Name       !< station name
   character (len=4)         :: SName      !< short name, 4 characters
   integer (i4)              :: ID         !< unique station ID
   integer (i4)              :: CenterID   !< ID defined by processing center
   character (len=20)        :: Country    !< country
   real (wp), dimension(1:3) :: CoordEll   !< ellipsoidal station coordinates
   real (wp), dimension(1:3) :: CoordGeo   !< geographical station coordinates
   real (wp), dimension(1:3) :: CoordCart  !< cartesian station coordinates
   real (wp)                 :: GeoidCorr  !< geoid correction
   real (wp)                 :: ModelSurf
   real (wp)                 :: StartMJD   !< coordinates valid from
   real (wp)                 :: EndMJD     !< coordinates valid until
   ! BUFR header
   integer (kind=2)          :: DataCategory
   integer (kind=2)          :: DataSubCategory
   integer (kind=2)          :: Center
   integer (kind=2)          :: SubCenter
end type GNSSStation



!---------------------------------------------------------------------
! structure StrucDomes
!---------------------------------------------------------------------
!
! Start of RoboDoc comments >>>>
!****T* GPSImport/StrucDomes
!
! Name
! StrucDomes
!
! Purpose
! Structure holding the information from one row of the file
! "smark_domes_SNX_RTT" which provides the ID, short name, full name,
! country and period (start date and end date).
! File location:
! /dsk/rtt17/m1/trop_base/ana/global/dat/smark_domes_SNX
! /dsk/igs2/dat/sta/smark_domes_SNX_RTT
!
! Parameters
! Name      - full name of GNSS station
! SName     - short name of GNSS station
! ID        - station ID
! Country   - station is located in this country
! domes     - unique station identifier from ITRF
! StartMJD  - station information is valid from StartMJD
! EndMJD    -                              until EndMJD
!---------------------------------------------------------------------
! Declaration
!---------------------------------------------------------------------
type StrucDomes
   character (len=23)     :: Name
   character (len=4)      :: SName
   integer (i4)           :: ID
   character (len=20)     :: Country
   character (len=9)      :: domes
   real (wp)              :: StartMJD
   real (wp)              :: EndMJD
end type StrucDomes
!****
! <<<<<<<< End of RoboDoc comments


!---------------------------------------------------------------------
! structure SlantTotalDelay
!---------------------------------------------------------------------
!
! Start of RoboDoc comments >>>>
!****T* DataLists/SlantTotalDelay
!
! Name
! SlantTotalDelay
!
! Purpose
! Stores data from slant-delays files, Ver1.0
!
! Parameters
! time      - Modified Julian date of data acquisiton
! station   - GPS station ID
! site      - site ID at GPS station
! satellite - satellite ID
! elevation - elevation in radian
! azimuth   - azimuth in radian
! slant     - slant delay in m
! zslant    - slant delay mapped into zenith direction in m
! --------------------------------------------------------------------
! simulated slant data
! STD       - Slant Total Delay (simulated)
! SZD       - Slant Zenith Delay (simulated)
! SDD       - Slant Dry Delay (simulated)
! SZDD      - Slant Zenith Dry Delay (simulated)
! SWD       - Slant Wet Delay (simulated)
! SZWD      - Slant Zenith Wet Delay (simulated)
! SWV       - Slant Water Vapour (simulated)
! SZWV      - Slant Zenith Water Vapour (simulated)
! --------------------------------------------------------------------
! Dry and wet delay computed from "slant"
! A model is used to estimate the dry delay based on the meteorological
! surface data defined below. The wet delay is assumed to be the difference
! swet = slant - sdry
! sdry      - dry slant delay
! swet      - slant delay due to water vapour
! --------------------------------------------------------------------
! Meteorological data from TRO file
! TROTOT    - total zenith delay (15 min. mean value)
! TROWET    - zenith wet delay (15 min. mean value)
! Press     - surface pressure at the GPS station in hPa
! Temp      - teperature in K
! Humi      - ?? humidity ??
!---------------------------------------------------------------------
! Declaration
!---------------------------------------------------------------------
type SlantTotalDelay
   real (wp)       :: time
   integer (i2)    :: station
   integer (i2)    :: site
   integer (i2)    :: satellite
   integer (i2)    :: GNSS
   real (wp)       :: elevation
   real (wp)       :: azimuth
   real (wp)       :: DryMap
   real (wp)       :: slant
   real (wp)       :: zslant
   !
   real (wp)       :: STDerr
   integer (i2)    :: STDacc
   real (wp)       :: SWD
   real (wp)       :: SIWV
   ! model equivalents of observations
   real (wp)       :: mSTD
   real (wp)       :: mSWD
   real (wp)       :: mSIWV
#ifdef __COSMO__
   ! state of observation => for feedback files  -  COSMO only
   integer         :: status      ! report body, flag single observation
   integer         :: check       ! report body, reason for rejecting obs.
   integer         :: icenter     ! GNSS processing center
   integer         :: processing  ! product of processing center
#endif
   ! for tests only
   integer (i2)    :: PE
   real (wp)       :: STDline

   ! simulated slant data
   !real (wp)       :: STD, SZD
   !real (wp)       :: SDD, SZDD
   !real (wp)       :: SWD, SZWD
   !real (wp)       :: SWV, SZWV
   ! dry and wet delay computed from "slant"
   !real (wp)       :: sdry
   !real (wp)       :: swet
   ! Meteorological data from TRO file
   !real (wp)       :: TROTOT
   !real (wp)       :: TROWET
   real (wp)       :: Press
   !real (wp)       :: Temp
   !real (wp)       :: Humi
   !
   real (wp), dimension(1:3) :: PosEll ! position on slant
                                       ! used as reference point

   character (len=8)  :: Name    ! COSMO test
end type SlantTotalDelay
!****
! <<<<<<<< End of RoboDoc comments


!---------------------------------------------------------------------
! structure SlantHeader
!---------------------------------------------------------------------
!
! Start of RoboDoc comments >>>>
!****T* GPSImport/SlantHeader
!
! Name
! SlantHeader
!
! Purpose
! Stores data from slant-delays files, Ver1.0 and Ver2.0
!
! Parameters
! Version  - Slant delay file format version    |
! Centre   - processing centre                  |
! Software - processing software                | information available
! Output   - data type, e.g. GPS slant delays   | since version 2.0
! Created  - date of data (re)processing        |
!            date string converted to MJD       |
! ..........................................
! NStat    - Number of stations in file
! NSat     - Number of satellites in file
! SRate    - sampling rate
! RefMJD   - Reference modified julien date
! ..........................................
! File     - File name of slant-delays files
! NSlants  - Number of slants in file
! First    - MJD of first slant (smallest MJD in file)
! Last     - MJD of last slant  (largest MJD in file)
! GPSWeek  - GPS week of data acquisition
! GPSWeekDay - Day within above GPS week (day = 1 (Sunday) - 7 (Saturday))
! Sort1    - Sorting criteria, primary key, Sort1 = 'stat', 'time', 'none'
!            Linked list of slants is sorted by station ID or time
! Sort2    - Secondary key, all data sets with identical primary key are
!            sorted by the secondary key:  Sort2 = 'stat', 'time', 'none'
! StatFirst - First station ID
! StatLast  - Last station ID
! ..........................................
! Optional information to simulated data
! ..........................................
! OrigFile  - Original slant delay file
! ExtraFile - Additional slant delay file, used if simulation period
!             covers more than one day, e.g. 23:00 - 1:00 h
! LMFile    - LM analysis file used for the simulation
! LMMJD     - Modified julian date of the LM analysis
! SimStart  - MJD start of simulation
! SimEnd    - MJD end of simulation
! NModel    - Model used to simulate refractivities N
!             e.g. Smich & Weintraub or Thayer
! ..........................................
! Status information
! ..........................................
! NSlantsRead  - number of slants copied to the linked list
! UpdateHeader - ExportESTD-Command: Compose new header even if an old header
!                could be copied
! ..........................................
! Additional comments:
!
! "First" and "Last" show the smallest and largest Modified
! Julian Date of all slants in the file. (This is independent of the order
! of slants in the file.) In contrast, "SimStart" and "SimEnd"
! show the beginning and the end of the simulation period as given in the
! parameter file. There might be no slants taken exactly at these dates, but
! there is at least one slant taken exactly at "First" and "Last".
!---------------------------------------------------------------------
! Declaration
!---------------------------------------------------------------------
type SlantHeader
   integer             :: Version
   character (len=256) :: Centre
   character (len=40)  :: Software
   character (len=40)  :: Output
   real (wp)           :: Created
   ! ...............................
   integer (i4)        :: NStat
   integer (i4)        :: NSat
   integer (i4)        :: SRate
   real (wp)           :: RefMJD
   ! ...............................
   character (len=256) :: File = ''
   integer (i4)        :: NSlants
   real (wp)           :: First
   real (wp)           :: Last
   integer (i4)        :: GPSWeek
   integer (i2)        :: GPSWeekDay
   character (len=4)   :: Sort1, Sort2
   integer (i4)        :: StatFirst, StatLast
   ! ...............................
   character (len=256) :: OrigFile = ''
   character (len=256) :: ExtraFile = ''
   character (len=256) :: LMFile = ''
   real (wp)           :: LMMJD
   real (wp)           :: SimStart
   real (wp)           :: SimEnd
   character (len=20)  :: NModel
   ! ...............................
   integer (i4)        :: NSlantsRead
   integer (i4)        :: NStatRead
   logical             :: UpdateHeader
end type SlantHeader
!****
! <<<<<<<< End of RoboDoc comments


!---------------------------------------------------------------------
! structure GPSStation
!---------------------------------------------------------------------
!
! Start of RoboDoc comments >>>>
!****T* GPSImport/GPSStation
!
! Name
! GPSStation
!
! Purpose
! GPS-Stationsdaten
!
! Parameters
! Name  - voller Name der GPS-Station
! SName - Kurzname der Station
! ID    - ID der Station
! Pos   - Position, kartesisch (X,Y,Z) aber in unterschiedlichen
!         Bezugssystemen, abh"angig vom Netzwerk
! GeoidCorr - geoid correction of station height
!---------------------------------------------------------------------
! Declaration
!---------------------------------------------------------------------
type GPSStation
   character (len=20)        :: Name
   character (len=4)         :: SName
   integer (i4)              :: ID
   real (wp), dimension(1:3) :: Pos
   real (wp)                 :: GeoidCorr
end type GPSStation
!****
! <<<<<<<< End of RoboDoc comments


!---------------------------------------------------------------------
! structure GNSSheader
!---------------------------------------------------------------------
!
! Name
! GNSSheader
!
! Purpose
! Structure containing the BUFR/netCDF header information.
!
! Parameters
!
!---------------------------------------------------------------------
  type Hchar
     character (len=30)                             :: name
     integer                                        :: ID
     character (len=1)                              :: fill
     character (len=1), dimension(:,:), allocatable :: vartmp
     character (len=20), dimension(:), allocatable  :: var
  end type Hchar

  type Hint
     character (len=30)                 :: name
     integer                            :: ID
     integer                            :: fill
     integer, dimension(:), allocatable :: var
  end type Hint

  type Hfloat
     character (len=30)              :: name
     integer                         :: ID
     real                            :: fill
     real, dimension(:), allocatable :: var
  end type Hfloat

  TYPE GNSSheader
     integer       :: icenter      ! processing center
     integer       :: processing   ! product
     integer       :: status       ! assimilation state of product
     integer       :: check        ! reason for rejecting observation
     TYPE (Hint)   :: center
     TYPE (Hint)   :: subcenter
     TYPE (Hint)   :: category
     TYPE (Hint)   :: subcategory
     TYPE (Hint)   :: year
     TYPE (Hint)   :: month
     TYPE (Hint)   :: day
     TYPE (Hint)   :: hour
     TYPE (Hint)   :: minute
     TYPE (Hchar)  :: station
     TYPE (Hfloat) :: lat
     TYPE (Hfloat) :: lon
     TYPE (Hint)   :: height
  END TYPE GNSSheader
!---------------------------------------------------------------------


!---------------------------------------------------------------------
! structure EllipsParam
!---------------------------------------------------------------------
!
! Start of RoboDoc comments >>>>
!****T* Coord/EllipsParam
!
! Name
! EllipsParam
!
! Struktur zur Aufnahme der Parameter verschiedener Referenz-Ellipsoide
!---------------------------------------------------------------------
! Declaration
!---------------------------------------------------------------------
type EllipsParam
   ! Name of parameter set
   character (len=15)   :: name
   ! Grosse Halbachse des Referenz-Ellipsoids a
   real (wp) :: a
   ! Kleine Halbachse des Referenz-Ellipsoids b
   real (wp) :: b
   ! Polkr"ummungsradius c
   real (wp) :: c
   ! Abplattung f
   real (wp) :: f
   ! Lineare Exzentrizit"at E
   real (wp) :: E
   ! Quadrat der ersten numerischen Exzentrizit"at EN1sup2
   real (wp) :: EN1sup2
   ! Quadrat der zweiten numerischen Exzentrizit"at EN2sup2
   real (wp) :: EN2sup2
   ! L"angenverh"altnis n
   real (wp) :: n
   ! Winkelgeschwindigkeit omega
   real (wp) :: omega
   !Geozentrische Gravitationskonstante (einschl. Atmosph"are) GM
   real (wp) :: GM
   ! 2. Zonale Harmonische C20
   real (wp) :: C20
end type EllipsParam
!****
! <<<<<<<< End of RoboDoc comments


!---------------------------------------------------------------------
! global variable WGS84Param
!---------------------------------------------------------------------
!
! Start of RoboDoc comments >>>>
!****V* Coord/WGS84Param
!
! Name
! WGS84Param
!
! Parameter der Referenz-Ellipsoide
! WGS84
!
! Parameter-Festlegung von Stefan Schl"ater, IDL-Programm
! ...........................
!  a = 6378137.000D
!  b = 6356752.314D
!  Ely.a    = a
!  Ely.b    = b
!  Ely.E   = SQRT(a^2-b^2)
!  Ely.f   = (a-b)/a
!  Ely.ea  = Ely.E/a
!  Ely.eb  = Ely.E/b
!  Ely.c   = a^2/b
!  Ely.n   = (a-b)/(a+b)
!  Ely.ea2 = Ely.ea^2
!  Ely.eb2 = Ely.eb^2
!  Ely.ea4 = Ely.ea^4
!  Ely.eb4 = Ely.eb^4
!  Ely.ea6 = Ely.ea^6
!  Ely.eb6 = Ely.eb^6
!  Ely.ea8 = Ely.ea^8
!  Ely.eb8 = Ely.eb^8
!---------------------------------------------------------------------
! Declaration
!---------------------------------------------------------------------
type (EllipsParam), parameter :: WGS84Param = EllipsParam &
 (                    &
   ! All numbers must be double precision, i. e. given with "_wp"
   ! The compiler will otherwise assume single precision, truncate the
   ! numbers and convert the truncatetd numbers to double precision.
   !
   ! name of parameter set
   'wgs84',           &
   ! Grosse Halbachse des Referenz-Ellipsoids a
   6378137.0_wp,       & ! m (Meter) 6378137 m
   ! Kleine Halbachse des Referenz-Ellipsoids b
   6356752.31424518_wp,       & ! m (Meter)
   ! Polkr"ummungsradius c
   6399593.626_wp,       & !
   ! Abplattung f
   1.0_wp/298.257223563_wp, & ! dimensionslos
   ! Lineare Exzentrizit"at E
   521854.0114151_wp,    & !
   ! Quadrat der ersten numerischen Exzentrizit"at EN1sup2
   !6.694379990141316E-003,      &
   6.6943799901411599688E-3_wp,     &
   ! Quadrat der zweiten numerischen Exzentrizit"at EN2sup2
   !6.739496742276434D-003,      &
   6.7394967422762758038E-3_wp,      &
   ! L"angenverh"altnis n
   0.0016792204_wp,      &
   ! Winkelgeschwindigkeit omega
   7.2921151467E-5_wp,   &
   !Geozentrische Gravitationskonstante (einschl. Atmosph"are) GM
   3.986005E14_wp,       & ! m^3 s^-2
   ! 2. Zonale Harmonische C20
   -484.16685E6_wp      &
  )
!****
! <<<<<<<< End of RoboDoc comments

TYPE :: DATE_TYPE
   INTEGER :: YEAR_J   ! year of end of Julian calendar
   INTEGER :: MONTH_J  ! month of end of Julian calendar
   INTEGER :: DAY_J    ! day of end of Julian calendar
   INTEGER :: YEAR_G   ! year of start of Gregorian calendar
   INTEGER :: MONTH_G  ! month of start of Gregorian calendar
   INTEGER :: DAY_G    ! day of start of Gregorian calendar
   INTEGER :: NDAYS    ! number of days dropped from calendar at switch
END TYPE DATE_TYPE

!=====================================================================
contains
!=====================================================================


!---------------------------------------------------------------------
! subroutine ImportStakoSNXarr
!---------------------------------------------------------------------
!
! Start of RoboDoc comments >>>>
!****S* GPSImport/ImportStakoSNXarr
!
! Name
! ImportStakoSNXarr
!
! Purpose
! Import data from StakoSNX files (e. g. stako_gps_SNX_G_07_HORI)
! Import station coordinates including the height above NN
!
! The last two columns of the HORI files give the start and end date for
! for each coordinate set. A new line is added if the stations position
! and/or name changes.
!
! HORI file:
! Col  1 : station's short name, 4 characters
! Col  2 : GFZ ID, 4 digit integer
! Col  3 : longitude, degrees, 0 - 360
! Col  4 : latitude, degrees, -90 - 0 - +90
! Col  5 : height above sea level, meter
! Col  6 : height above ellipsoid, meter
!  ...
!  ...   : format for synop, not read by this routine
!  ...
! Col 12 : data valid from MJD
! Col 13 : data valid until MJD
! Col 14 : station's short name, as in col. 1, not read
!
! Example (header and first row):
! *             Longitude    Latitude     Hei_NN     Hei_Ell
!  gsr1 1005   14.5437133   46.0481312    305.781    351.659
!
!  -------Format fuer Synop------  Start  End
! 1005  46.05  14.54   305m   305m 51000 99999 gsr1
!                                   MJD   MJD
!
! Call
! call ImportStakoSNXarr(File, STATlist, NStat)
! call ImportStakoSNXarr(File, STATlist, NStat, MinID, MaxID)
!
! Input
! File      - File name of StakoSNX file, optionally with path
! STATlist  - Array of station coordinates. Any old data will be deleted.
! MJD       - The stations coordinates will be returned for a specific MJD
!             (MJD - Modified Julian Date), optional
!             The MJD is optional only for backward compatibility, it must
!             be given if the correct coordinates at that day are required.
!             If the MJD is omitted, the most recent data are read (the
!             last line of coordinates of that station.)
! EndMJD    - If EndMJD is present, all stations with valid coordinates
!             between MJD (=StartMJD) and EndMJD will be regarded.
!             (EndMJD - Modified Julian Date), optional
!
! Output
! STATlist  - Array of of station coordinates, the array dimesion is
!             exactly the number of stations read, i.e. NStat.
!             Stucture "GNSSStation" is filled with the station data,
!             except the cartesian coordinates which cannot be computed here.
!             STATlist%CoordCart = 0.0
! NStat     - Number of stations read.
! MinID     - Min. station ID found in file (optional)
! MaxID     - Max. station ID found in file (optional)
!
! External References
! None
!
! Record of revisions
! Date        Programmer   Description of change
! ====        ==========   =====================
! 30.05.2012  M. Bender    new, copy of "ImportStakoSNX"
!                          linked lists replaced by arrays, new data
!                          structure "GNSSStation" used
! 08.06.2012  M. Bender    Lat/lon converted to radian, geoid correction
!                          computed, cartesian coordinates = 0
! 26.03.2014 M. Bender     HORI file entries with identical IDs are not
!                          copied to the station list => avoid wrong
!                          entries in HORI file
! 26.07.2017  M. Bender    Input arrays is automatically increased if number
!                          of station increases original size, support for
!                          new SINEX HORI format.
!---------------------------------------------------------------------
! Source
!---------------------------------------------------------------------
subroutine ImportStakoSNXarr(File, STATlist, NStat, MinID, MaxID, MJD, EndMJD)

! List of calling arguments:
character (len=*), intent(in)                           :: File
#ifdef __COSMO__
type (GNSSStation), dimension(:), allocatable           :: StatList
#else
type (GNSSStation), dimension(:), pointer               :: StatList
#endif
integer ,                         intent(out)           :: NStat
integer ,                         intent(out), optional :: MinID, MaxID
real(wp),                                      optional :: MJD
real(wp),                                      optional :: EndMJD

! List of local variables:
character (len=256)  :: row
integer              :: MinStatID, MaxStatID, OldID
integer              :: MJDi, EndMJDi

! >>> All stations in the HORI file must fit into this array,
!     the dimension must be reset if the file grows !!!!
type (GNSSStation), dimension(:), allocatable  :: TmpList, buf
type (GNSSStation)                             :: TmpRec

integer  :: unit
logical  :: offen
integer  :: iostat
logical  :: debug = .false.
!---------------------------------------------------------------------

!debug = .true.

if (debug) write(*,*) 'ImportStakoSNXarr> Start ...'

if (present(MJD)) then
   ! MJD was given, use the interger value
   MJDi = int(MJD)
else
   ! No MJD was given, read the last line with the most recent data
   ! The current coordinates are valid until MJD=99999
   MJDi = 99990
end if
if (present(EndMJD)) then
   ! MJD was given, use the interger value
   EndMJDi = int(EndMJD)
   if (MJDi .gt. EndMJDi) then
      ! MJD not given??? Set MJD to 51000, i.e. the oldest coordinates
      MJDi = 51000
   end if
else
   ! No EndMJD was given, don't read coordinates from a range but only for
   ! the given date => set EndMJDi to a negative number
   EndMJDi = -99
end if

! This routine doesn't provide cartesian coordinates of the stations:
! Initialize with zeros
TmpRec%CoordCart = 0.0_wp

unit = 40
offen = .true.
do while (offen)
   unit = unit + 1
   inquire(UNIT = unit, OPENED = offen )
end do

open( UNIT =   unit,                &
      FILE =   trim(File),          &
      STATUS = 'OLD',               &
      FORM =   'Formatted',         &
      ACTION = 'READ',              &
      IOSTAT = iostat                )

if (iostat .eq. 0) then
   ! File opened successfully

   if (debug) write(*,*) 'ImportStakoSNXarr> Opened file ', trim(File)

   ! >>> All stations in the HORI file must fit into this array,
   !     the dimension must be reset if the file grows !!!!
   allocate(TmpList(1:2000))

   MinStatID = 99999
   MaxStatID = 0
   OldID     = -99
   NStat = 0

   do
   ! Read ASCII file line by line until EOF or '-STA_COORDINATES'

      if (Nstat >= ubound(TmpList,1)) then
         ! increase input list size
         call move_alloc(TmpList, buf)
         allocate( TmpList(2*Nstat) )
         TmpList(1:Nstat) = buf
         deallocate( buf )
      end if

      read(UNIT=unit,IOSTAT=iostat,FMT='(a)') row

      if (iostat .ne. 0) exit      ! End of file ?
      if (index(row,'+STA_COORDINATES') .gt. 0) cycle ! begin of station block
      if (index(row(1:4),'*') .gt. 0) cycle  ! comment
      if (index(row(1:2),'%') .gt. 0) cycle  ! SINEX header/footer
      !write(*,*) row(1:10), index(row,'-STA_COORDINATES')
      if (index(row,'-STA_COORDINATES') .gt. 0) exit ! end of station block

      !if (debug) write(*,*) 'ImportStakoSNXarr> Record ', trim(row)

      read(row,'(tr1,a4,tr1,i4,tr1,f12.7,tr2,f12.7,tr2,f9.3,tr2,f9.2,tr35,  &
    &             f5.0,tr1,f5.0)')                                          &
          TmpRec%SName, TmpRec%ID,                                          &
          TmpRec%CoordGeo(1), TmpRec%CoordGeo(2), TmpRec%CoordGeo(3),       &
          TmpRec%CoordEll(3),                                               &
          TmpRec%StartMJD, TmpRec%EndMJD

      ! Longitude and latitude are degrees => convert to radian
      TmpRec%CoordGeo(1:2) = TmpRec%CoordGeo(1:2)*deg2rad

      ! Set Geoid correction for the station:
      TmpRec%GeoidCorr = TmpRec%CoordGeo(3) - TmpRec%CoordEll(3)

      ! Longitude and latitude are identical for geographical and
      ! elliptical coordinates: copy CoordGeo(1:2) to CoordEll(1:2)
      TmpRec%CoordEll(1:2) = TmpRec%CoordGeo(1:2)

      TmpRec%SName = to_upper(TmpRec%SName)
      !write(*,*) TmpRec%SName, TmpRec%ID, TmpRec%StartMJD, &
      !           TmpRec%EndMJD

      if (int(TmpRec%StartMJD) .le. MJDi .and.          &
          int(TmpRec%EndMJD)   .ge. MJDi       ) then
         ! The coordinates valid at MJDi have been found

         if ( TmpRec%ID .ne. OldID) then
            ! Only one entry for each station in station list
            NStat = NStat + 1
            TmpList(NStat) = TmpRec

            ! Find minimal and maximal station ID
            if (TmpRec%ID .gt. MaxStatID) MaxStatID = TmpRec%ID
            if (TmpRec%ID .lt. MinStatID) MinStatID = TmpRec%ID
            OldID = TmpRec%ID
         end if

      else if ( ( int(TmpRec%StartMJD) .ge.  MJDi    .and.               &
                  int(TmpRec%StartMJD) .le.  EndMJDi       ) .or.        &
                ( int(TmpRec%EndMJD)   .ge.  MJDi    .and.               &
                  int(TmpRec%EndMJD)   .le.  EndMJDi       )     ) then
         ! Station coordinates changed between MJDi and EndMJDi, use
         ! station

         if ( TmpRec%ID .ne. OldID) then
            ! Only one entry for each station in station list

            NStat = NStat + 1
            TmpList(NStat) = TmpRec

            ! Find minimal and maximal station ID
            if (TmpRec%ID .gt. MaxStatID) MaxStatID = TmpRec%ID
            if (TmpRec%ID .lt. MinStatID) MinStatID = TmpRec%ID
            OldID = TmpRec%ID

         end if

      else
         ! data set not corresponding to MJDi, skip

         !write(*,*) 'ImportStakoSNXarr> Skip station:'
         !write(*,*) TmpRec%SName, TmpRec%ID, TmpRec%StartMJD, &
         !           TmpRec%EndMJD
         !write(*,*) int(TmpRec%StartMJD), MJDi, int(TmpRec%EndMJD)


      end if

   end do  ! Read ASCII file line by line

   close(unit)

   if (present(MinID)) MinID = MinStatID
   if (present(MaxID)) MaxID = MaxStatID

   ! Copy temporary array TmpList to StatList
#ifdef __COSMO__
   if (allocated (StatList)) deallocate(StatList)
#else
   if (associated (StatList)) deallocate(StatList)
#endif
   allocate( StatList(1:NStat) )
   StatList = TmpList(1:NStat)
   deallocate(TmpList)

else
   write(*,*) 'ImportStakoSNXarr> Error opening file ', trim(File)
end if

if (debug) write(*,*) 'ImportStakoSNXarr> ... end'

end subroutine ImportStakoSNXarr
!****
! <<<<<<<< End of RoboDoc comments


!---------------------------------------------------------------------
! subroutine ImportESTDarr
!---------------------------------------------------------------------
!
! Start of RoboDoc comments >>>>
!****S* GPSImport/ImportESTDarr
!
! Name
! ImportESTDarr
!
! Purpose
! Import Extended Slant Total Delay Data (STD) from Slant delays file Ver. 1.0
! Reads file header and slant data in file
! Optionally a list of all stations and/or a list of all satellites can
! be provided. Slants from a given time window or a given Station may be
! selected.
!
! Call
! call ImportESTDarr(File, STDlist, Header, StatHash, StatList,  &
!                    NSat, SatList, MJDBegin, MJDEnd, StatID)
!
! Input
! File     - Filename slant data , including Path if required
! STDlist  - Array of STD data, old data will be deleted.
! MJDBegin - Modified Julian date (MJD), only slants after MJDBegin will be
!            copied to STDlist (optional)
! MJDEnd   - Modified Julian date (MJD), only slants before MJDEnd will be
!            copied to STDlist (optional)
!            Both parameters can be used to define a time window, but it is also
!            possible to give only one of them.
! StatID   - Station ID (4 digits), only slants from this station will be
!            copied to the linked list (optional).
!            StatID = 0 ca be used if only the header is to be read, the slant
!            section will be skiped for StatID < 1.
!            (So far either slants from all stations are read or slants
!             from only one station are read. It is not possible to provide
!             a station list.)
!
! Output
! STDlist  - Array of slant data.
!            STDlist is not associated if no data have been read, i.e.
!            if slant file doesn't exist, could not be opened or doesn't
!            contain any fitting data.
! Header   - Structure containing header information from slant file
! StatHash (optional) - Simple hash table providing the array index of a
!                       given station via the station ID
! StatList (optional) - List of stations in file (structure with station data)
!                       (1-dim array of station data)
! NSat (optional)     - Number of different satellites in file
! SatList (optional)  - List of satellites in file (IDs)
!                       1-dim array, array index is satellite id:
!                       SatList(i) = 1  if satellite with ID = i is present
!                       SatList(i) = 0  if there is no satellite with ID = i
!
! If "StatList" is used, "NStat" must not be omitted!
! If "SatList" is used, "NSat" must not be omitted!
!
! External References
! None
!
!
! Procedure
!
! Short description of slant delays file
! --------------------------------------
! Header:
! reference_mjd
! Modified Julian Date, reference for time data (first column, TIME).
! The data acquisition time of a specific slant is given as an offset to
! the reference time:
! reference_mjd + TIME = Modified Julian Date
! reference_mjd + TIME + 2400000.5 = Julian Date
!
! Body:
! TIME
! Modified Julian Date of data acquisition, offset to "reference_mjd" as
! given in the files header.
!
! STA
! Station ID (4 digits) + 2 digits indicating a specific site at the station
! Only the first 4 digits are required to identify a given station in
! the corresponding SINEX file
!
! SAT
! Satellite ID
!
! elev
! Elevation (degree) in the stations local horizon system
!
! azimu
! Azimuth (degree) in the stations local horizon system
!
! slant_delay(m)
! Total slant delay (STD) (m)
!
! slant_delay_norm(m)
! Total slant delay (STD) mapped into zenit direction (m)

! Record of revisions
! Date        Programmer   Description of change
! ====        ==========   =====================
! 17.01.2007  M. Bender    new
! 25.01.2007  M. Bender    Datum of first and last slant in file is saved
! 13.02.2007  M. Bender    Changed parameter list, use "SlantHeader"
! 28.02.2007  M. Bender    changes due to new slant format
! 09.03.2007  M. Bender    new slant format, station coordinates
! 21.03.2007  M. Bender    New optional arguments: MJDBegin, MJDEnd
!                          Read only Slants within this time window
! 29.03.2007  M. Bender    Allocatable arrays are deallocated before the
!                          first alocate statement - this is required if
!                          the routine is called multiple times
! 26.04.2007  M. Bender    New optional argument: StatID
!                          Read only slants from this station (so far only
!                          one station is supported)
! 16.05.2007  M. Bender    Skip slants section if "StatID" is present and < 1
!                          Save number of slants read in header.
! 22.05.2007  M. Bender    Hash table: So far a station list sorted by station
!                          ID has been assumed. Now usorted lists are supported.
! 23.05.2007  M. Bender    Wrong satellite list fixed: SatList(i)=1 if satellite
!                          with ID=i is present, SatList(i)=0 otherwise.
! 18.07.2007  M. Bender    The number of stations read from the file header is
!                          saved in haeder%NStatRead
! 01.04.2009  M. Bender    Current version copied to "ImportESTDOld", this
!                          version makes use of the linked list module
! 21.10.2009  M. Bender    Update to version 2.0 of slant file format
! 12.02.2010  M. Bender    SName convertet to capital letters
! 24.05.2012  M. Bender    Copy of "ImportESTD", linked lists removed,
!                          data are copied to arrays
! 14.11.2014  M. Bender    Reject rows with bad data (*** or NaN), simple
!                          quality check of data
!---------------------------------------------------------------------
! Source
!---------------------------------------------------------------------
subroutine ImportESTDarr(File, MaxObs, STDlist, Header, StatHash, StatList,  &
                         NSat, SatList, MJDBegin, MJDEnd, StatID)

! List of calling arguments:
character (len=*), intent(in)                              :: File
integer ,                         intent(in)               :: MaxObs
#ifdef __COSMO__
type (SlantTotalDelay), dimension(:), allocatable, intent(out) :: STDlist
#else
type (SlantTotalDelay), dimension(:), pointer              :: STDlist
#endif
type (SlantHeader), intent(out)                            :: Header
integer (i2), dimension(:),pointer,optional                :: StatHash
type (GPSStation), dimension(:), pointer, optional         :: StatList
integer (i4), intent(out), optional                        :: NSat
integer (i4),dimension(:), pointer, optional               :: SatList
real (wp), optional                                        :: MJDBegin, MJDEnd
integer (i4), intent(in), optional                         :: StatID

! List of local variables:
!type (LLSTDDataRecord), pointer    :: STD
character (len=256)                :: row
integer (i4)                       :: cstat, csat, NSlants, NStat, i
integer (i4)                       :: IDmin, IDmax, SatID
type (SlantTotalDelay), dimension(:), pointer          :: STDtmp
type (SlantTotalDelay)             :: lastSTD

integer (i4) :: unit
logical      :: offen
integer (i4) :: stat
integer      :: err
logical      :: debug, CopySlant, ReadSlants

integer (i4), parameter :: MaxSat  = 100
!---------------------------------------------------------------------

debug = .false.
!debug = .true.

if (debug) write(*,*) 'ImportESTDarr> Start reading file ', trim(File)
if (len_trim(File) .le. 0) then
   write(*,*) 'No input file given.'
   stop
end if

! The optional argument "StatID" can be used to read the header information
! only and to skip the (large) slant data part.
ReadSlants = .true.
if (present(StatID)) then
   ! Skip slant data if StatID < 1 (not a valid station ID)
   if (StatID .lt. 1) ReadSlants = .false.
end if

! Deallocate old STDlist array if it exists
#ifdef __COSMO__
if ( allocated(STDlist) ) deallocate(STDlist)
#else
if ( associated(STDlist) ) deallocate(STDlist)
#endif

! Open slant-delay file
unit = 40
offen = .true.
do while (offen)
   unit = unit + 1
   inquire(UNIT = unit, OPENED = offen )
end do

open( UNIT =   unit,                &
      FILE =   trim(File),          &
      STATUS = 'OLD',               &
      FORM =   'Formatted',         &
      ACTION = 'READ',              &
      IOSTAT = stat                )

if (stat .eq. 0) then
   ! Read slant-delays file

   Header%File = trim(File)
   Header%OrigFile = ' '
   Header%LMFile   = ' '
   Header%NModel   = ' '
   Header%Version  = -1

   ! Read file header
   row = ''
   do while (index(row,'-slant-delay_param') .eq. 0 .and.         &
             index(row,'-slant_delay_param') .eq. 0      )
      read(unit,'(a80)') row
      !write(*,*) '>>', trim(row), '<<'

      ! Read slant delay format version (supported: 1.0 and 2.0)
      if ( index(row,'Slant-delay') .ne. 0 ) then
         !read(row(27:27),'(i1)') Header%Version
         if ( index(row,'1.0') .ne. 0 ) Header%Version = 1
         if ( index(row,'2.0') .ne. 0 ) Header%Version = 2
         if (Header%Version .lt. 1 .or. Header%Version .gt. 2) then
            write(*,*) 'ImportESTDarr>>> ',                           &
                       'Slant file format version not supported:', &
                        Header%Version
            stop
         end if
         cycle
      end if

      ! Read "FILE/REFERENCE" section, available since version 2.0
      if (Header%Version .gt. 1) then

         if ( index(row,'Processing Centre') .ne. 0 ) then
            ! Read name of GPS processing centre providing the slant files
            Header%Centre = trim(row(20:))
            cycle
         end if
         if ( index(row,'SOFTWARE') .ne. 0 ) then
            ! Read name of GPS processing software
            Header%Software = trim(row(12:))
            cycle
         end if
         if ( index(row,'OUTPUT') .ne. 0 ) then
            ! Read data type written to this file
            Header%Output = trim(row(12:))
            cycle
         end if
         if ( index(row,'CREATED') .ne. 0 ) then
            ! Read (re)processing date of data
            !Header%Created =
            cycle
         end if

      end if
      ! "FILE/REFERENCE" section read

      if ( index(row,'number_of_stations') .ne. 0 ) then
         ! Read number of different stations which provide slants
         !write(*,*) '>>',row(3:15),'<<'
         read(row(21:23),'(i3)') Header%NStat
         NStat = Header%NStat
         if ( present(StatList) ) then
            if (associated(StatList)) then
               ! StatList already allocated
               !write(*,*) 'ImportESTD> Warning - Delete old station list'
               deallocate(StatList)
               allocate(StatList(1:NStat))
            else
               allocate(StatList(1:NStat))
            end if
            ! The station ID must allways have a defined state:
            ! ID = 0  - empty element
            ! ID > 0  - station data have been copied to this element
            StatList(:)%ID = 0
         end if
         cycle
      end if
      if ( index(row,'sampling_rate') .ne. 0 ) then
         ! Read time reference (Modified Julian Date)
         read(row(16:18),'(i3)') Header%SRate
         cycle
      end if
      if ( index(row,'number_of_satellites') .ne. 0 ) then
         ! Read number of satellites
         read(row(23:25),'(i3)') Header%NSat
         if (present(NSat)) then
            NSat = Header%NSat
         end if
         if ( present(SatList) ) then
            if (associated(SatList)) then
               ! SatList already allocated
               !write(*,*) 'ImportESTD> Warning - Delete old satellite list'
               deallocate(SatList)
               allocate(SatList(1:MaxSat))
            else
               allocate(SatList(1:MaxSat))
            end if
         end if
         cycle
      end if
      if ( index(row,'number_of_slants') .ne. 0 ) then
         ! Read number of slants
         read(row(19:27),'(i9)') Header%NSlants
         cycle
      end if
      if ( index(row(1:15),'reference_MJD') .ne. 0 ) then
         ! Read time reference (Modified Julian Date)
         read(row(16:27),'(f12.6)') Header%RefMJD
         cycle
      end if
      if ( index(row,'begin_MJD') .ne. 0 ) then
         ! Read Modified Julian Date of first slant
         read(row(12:23),'(f12.6)') Header%First
         cycle
      end if
      if ( index(row,'end_MJD') .ne. 0 ) then
         ! Read Modified Julian Date of last slant
         read(row(10:21),'(f12.6)') Header%Last
         cycle
      end if
      if ( index(row(1:10),'GPS_week') .ne. 0 ) then
         ! Read Modified Julian Date of last slant
         read(row(11:14),'(i4)') Header%GPSWeek
         cycle
      end if
      if ( index(row,'day_of_GPS_week') .ne. 0 ) then
         ! Read Modified Julian Date of last slant
         read(row(18:18),'(i1)') Header%GPSWeekDay
         cycle
      end if
      if ( index(row,'sort_key1') .ne. 0 ) then
         ! Read sorting criteria - primary key
         read(row(12:15),'(a4)') Header%Sort1
         cycle
      end if
      if ( index(row,'sort_key2') .ne. 0 ) then
         ! Read sorting criteria - secondary key
         read(row(12:15),'(a4)') Header%Sort2
      end if
   end do

   do while (index(row,'+sta_number_name_coordinate') .eq. 0)
      ! Look for simulated data
      read(unit,'(a256)') row
      !write(*,*) '2>>', trim(row), '<<'

      if ( index(row,'slant_orig_file') .ne. 0 ) then
         ! Read original file name
         !write(*,*) row(18:)
         read(row(18:),'(a)') Header%OrigFile
         !write(*,*) trim(Header%OrigFile)
      end if
      if ( index(row,'lm_analysis') .ne. 0 ) then
         ! Read name of LM analysis file
         read(row(14:),'(a)') Header%LMFile
      end if
      if ( index(row,'lm_mjd') .ne. 0 ) then
         ! Read Modified Julian Date of the LM analysis
         read(row(9:20),'(f12.6)') Header%LMMJD
      end if
      if ( index(row,'mjd_start') .ne. 0 ) then
         ! Read Modified Julian Date of first simulated slant
         read(row(12:23),'(f12.6)') Header%SimStart
      end if
      if ( index(row,'mjd_end') .ne. 0 ) then
         ! Read Modified Julian Date of last simulated slant
         read(row(10:21),'(f12.6)') Header%SimEnd
      end if
      if ( index(row,'n_model') .ne. 0 ) then
         ! Read model used to simulate atmospheric refractivity
         read(row(10:),'(a)') Header%NModel
      end if

   end do

   if (debug) write(*,*) 'ImportESTDarr> Header read, start reading stations'

   if ( present(StatList) ) then
      ! Read station ID's of all stations used in this file

!!$      do while (index(row,'+slant-delay_sta') .eq. 0)
!!$         ! Move foreward until "+slant-delay_sta" is found
!!$         read(unit,'(a80)') row
!!$      end do

      cstat = 0
      do
         read(unit,'(a80)') row
         !write(*,*) trim(row)
         if (index(row,'-sta_number_name_coordinate') .ne. 0) exit
         if (index(row(1:1),'*') .ne. 0) cycle
         cstat = cstat + 1
         read(row(2:5),'(i4)')      StatList(cstat)%ID
         read(row(7:10),'(a4)')     StatList(cstat)%SName
         read(row(12:25),'(f14.5)') StatList(cstat)%Pos(1)
         read(row(27:40),'(f14.5)') StatList(cstat)%Pos(2)
         read(row(42:55),'(f14.5)') StatList(cstat)%Pos(3)
         !write(*,*) trim(row), StatList(cstat)
         StatList(cstat)%SName = to_upper(StatList(cstat)%SName)
      end do

      if (Header%NStat .ne. cstat) then
         write(*,*) 'ImportESTDarr> Error - incorrect number of stations', &
                     NStat, cstat
      end if
      Header%NStatRead = cstat
      NStat = cstat

      if (present(StatHash)) then
         ! Allocate hash table which links the station ID to the
         ! array index of StatList.

          if (debug) write(*,*) 'ImportESTDarr> Creating Hash table ...'

         ! Find smallest and largest station ID
         IDmin = 999999
         IDmax = 0
         do i=1, cstat
            if (StatList(i)%ID .gt. IDmax) IDmax = StatList(i)%ID
            if (StatList(i)%ID .lt. IDmin) IDmin = StatList(i)%ID
         end do
         Header%StatFirst = IDmin
         Header%StatLast  = IDmax

         if (associated(StatHash)) then
            deallocate(StatHash)
            allocate( StatHash(IDmin:IDmax) )
         else
            allocate( StatHash(IDmin:IDmax) )
         end if

         ! Fill hash table
         StatHash = 0
         do i=1, NStat
            StatHash(StatList(i)%ID) = i
            !write(*,*) 'Hash : ', i, StatList(i)%ID
         end do

      end if

   end if

   if (debug)                                                             &
     write(*,*) 'ImportESTDarr> Stations read, start reading satellites'

   if ( present(SatList) ) then
      ! Read satellite ID's of all satellites used in this file

      SatList = 0

      do while (index(row,'+satellites_used') .eq. 0)
         ! Move foreward until "+slant-delay_sat" is found
         read(unit,'(a80)') row
      end do

      csat = 0
      do
         read(unit,'(a80)') row
         if  (index(row,'-satellites_used') .ne. 0) exit
         if (index(row(1:1),'*') .ne. 0) cycle

         read(row(2:5),'(i4)') SatID
         if (SatID .gt. 0 .and. SatID .le. MaxSat) then
            SatList(SatID) = 1
            csat = csat + 1
         else if (SatID .gt. MaxSat) then
            write(*,*) 'ImportESTDarr> Error - satellite ID is too lagre : ', &
                        SatID
         end if
         !write(*,*) trim(row), SatList(csat)
      end do

      if (Header%NSat .ne. csat) then
         write(*,*) 'ImportESTDarr> Error - incorrect number of satellites', &
                     NSat, csat
      end if

      NSat = csat

   end if   ! if ( present(SatList) ) then

   ! Skip rest of header
   do while (index(row,'--------') .eq. 0)
      read(unit,'(a80)') row
      !write(*,*) '>>', trim(row), '<<'
   end do

   if (debug)                                                             &
     write(*,*) 'ImportESTDarr> Header finished, start reading slants? ', &
                          ReadSlants

   ! Allocate temporary array for slant data. Maximum number of slants was
   ! given in the header: Header%NSlants
   if (Header%NSlants <= MaxObs) then
      ! Number of observations in file less than MaxObs
      allocate( STDtmp(1:Header%NSlants) )
   else
      ! Number of observations in file exceeds MaxObs, read only the
      ! first observations
      allocate( STDtmp(1:MaxObs) )

     WRITE(*,*)
     WRITE(*,*) '  *** WARNING: Number of STD observations exceeds limit ', &
                'set in namelist variable MaxSTDobs ***'
     WRITE(*,*) '      MaxSTDobs = ', MaxObs, ', number of STDs = ', &
                Header%NSlants
     WRITE(*,*) '      All STD observations above this limit will be skipped.'
     WRITE(*,*)

   end if

   Nslants = 0

   if (ReadSlants) then
      ! Read slant data

      do
         ! Read slant data

         if (Nslants >= MaxObs) exit   ! max. number of observations reached

         read(UNIT=unit,IOSTAT=stat,FMT='(a256)') row
         !read(UNIT=unit,IOSTAT=stat) row
         !write(*,*) '>>', trim(row), '<<'
         !write(*,*) '>>', len_trim(row)
         !write(*,*) 'read error : ', stat

         if (stat .ne. 0) exit                     ! end of file
         if  (index(row,'-slants') .ne. 0) exit    ! end of slants section
         !if (index(row(1:1),'*') .ne. 0) cycle     ! skip comment
         if (index(row,'*') .ne. 0) cycle          ! skip comment or bad data
         if (index(row,'NaN') .ne. 0) cycle        ! bad data: Not a Number

         Nslants = Nslants + 1

         ! copy current row to STD
         call ReadESTDRowarr(row, lastSTD, Header%Version, err)
         if (err .ne. 0) then
            write(*,*) 'ReadESTDRowarr, error : ', err
            cycle
         end if
         STDtmp(Nslants) = lastSTD

         ! Save Modified Julian Date
         STDtmp(Nslants)%time = STDtmp(Nslants)%time + Header%RefMJD

         CopySlant = .true.
         if (present(MJDBegin) .and. present(MJDEnd)) then
            ! Copy only slants from the given time window to the linked list
            if ( STDtmp(Nslants)%time .ge. MJDBegin .and.               &
                 STDtmp(Nslants)%time .le. MJDEnd         ) then
               CopySlant = .true.
            else
               CopySlant = .false.
            end if
         else if (present(MJDBegin)) then
            ! Copy only slants after MJDBegin to the linked list
            if (STDtmp(Nslants)%time .ge. MJDBegin) then
               CopySlant = .true.
            else
               CopySlant = .false.
            end if
         else if (present(MJDEnd)) then
            ! Copy only slants before MJDEnd to the linked list
            if (STDtmp(Nslants)%time .le. MJDEnd) then
               CopySlant = .true.
            else
               CopySlant = .false.
            end if
         end if
         if (present(StatID) .and. CopySlant) then
            ! Copy only slants from this single station to the linked list
            !write(*,*) 'Id :', STD%station, StatID
            if (STDtmp(Nslants)%station .eq. StatID) then
               CopySlant = .true.
            else
               CopySlant = .false.
            end if
         end if

         ! Quality check
         if (STDtmp(Nslants)%elevation < 0.0_wp .or.             &
             STDtmp(Nslants)%elevation > 90.0_wp     ) then
            CopySlant = .false.
         end if
         if (STDtmp(Nslants)%azimuth < -180.0_wp .or.            &
             STDtmp(Nslants)%azimuth > 360.0_wp       ) then
            CopySlant = .false.
         end if
         if (STDtmp(Nslants)%slant < 0.1 .or.                    &
             STDtmp(Nslants)%slant > 100.0    ) then
            CopySlant = .false.
         end if
         if (STDtmp(Nslants)%zslant < 0.1 .or.                   &
             STDtmp(Nslants)%zslant > 10.0     ) then
            CopySlant = .false.
         end if

         if (.not. CopySlant) then
            ! delete current record
            Nslants = Nslants - 1
         end if

      end do

      if (debug) write(*,*) 'ImportESTDarr> Slants read : ', Nslants

   end if   ! if (ReadSlants) then

   Header%NSlantsRead = Nslants
   !write(*,*) 'Nslants : ',  Nslants

   !if (NSlants .ne. Header%NSlants) then
   !   write(*,*) ' ImportESTD> Warning - wrong number of slants read ', &
   !                NSlants
   !end if

   close(unit)

   ! Copy temporary array to "STDlist"
   if (Nslants > 0) then
      allocate( STDlist(1:Nslants) )
      STDlist(1:Nslants) = STDtmp(1:Nslants)
   end if

else
   write(*,*) 'ImportESTDarr> Error opening file >>', trim(File), '<<'
   Header%NSlantsRead = 0
end if

if ( associated(STDtmp) ) deallocate(STDtmp)

if (debug) write(*,*) 'ImportESTDarr> ... end'

end subroutine ImportESTDarr
!****
! <<<<<<<< End of RoboDoc comments



!---------------------------------------------------------------------
! subroutine ReadESTDRowarr
!---------------------------------------------------------------------
!
! Start of RoboDoc comments >>>>
!****S* GPSImport/ReadESTDRowarr
!
! Name
! ReadESTDRowarr
!
! Purpose
! Extract extended slant data from one row of ASCII data from the slant file
! The extended data set provides simulated slant data in addition to the
! original data
!
! Call
! call ReadESTDRowarr(row, STD, err)
!
! Input
! row  - row of ASCII data
!
! Output
! STD  - Structure element with slant record
! err  - Returncode, err = 0 - no error, err <> 0 - error
!
! External References
! None
!
! Record of revisions
! Date        Programmer   Description of change
! ====        ==========   =====================
! 28.02.2007  M. Bender    new
! 26.03.2008  M. Bender    wrong substring read for "swet" -> fixed
! 01.04.2009  M. Bender    current version copied to "ReadESTDRowOld",
!                          type of pointer "STD" changed to "LLSTDDataRecord"
! 22.10.2009  M. Bender    update to slant data format ver. 2
! 25.05.2012  M. Bender    Copy of ReadESTDRow
!                          Pointer "STD" replaced by structure element
!---------------------------------------------------------------------
! Source
!---------------------------------------------------------------------
subroutine ReadESTDRowarr(row, STD, ver, err)

! List of calling arguments:
character (len=*), intent(in)       :: row
type (SlantTotalDelay), intent(out) :: STD
integer, intent(in)                 :: ver
integer, intent(out)                :: err

! List of local variables:
integer                            :: rowlen, star
!---------------------------------------------------------------------
!write(*,*) 'ReadESTDRow > start ', len_trim(row)

err = 0
rowlen = len_trim(row)
star = index(row,'*')
if (rowlen .lt. 56) then
   ! row too short, cannot contain slant data
   err = 1
else if (star .gt. 0) then
   ! row contains "****", i. e. invalid numbers
   err = 2
else
   ! Copy data in slant structure

   ! First copy the original data in the first part of the row
   call ReadSTDRowarr(row, STD, ver, err)
   !write(*,*) 'ReadSTDRow error : ', err
!!$   STD%sdry = 0.0
!!$   STD%swet = 0.0
!!$
!!$   if (err .eq. 0 .and. rowlen .gt. 80 .and. rowlen .lt. 100) then
!!$      ! read only sdry and swet
!!$
!!$      read(row(59:70),'(f12.4)')   STD%sdry
!!$      read(row(72:83),'(f12.4)')   STD%swet
!!$
!!$      !write(*,*) 'ReadESTDRow > ', STD%sdry, STD%swet
!!$
!!$      STD%STD = 0.0
!!$      STD%SZD = 0.0
!!$      STD%SDD = 0.0
!!$      STD%SZDD = 0.0
!!$      STD%SWD = 0.0
!!$      STD%SZWD = 0.0
!!$      STD%SWV = 0.0
!!$      STD%SZWV = 0.0
!!$
!!$   else if (err .eq. 0 .and. rowlen .gt. 100) then
!!$      ! Copy simulated data in the remainig part of the row
!!$
!!$      read(row(59:70),'(f12.4)')   STD%sdry
!!$      read(row(72:83),'(f12.4)')   STD%swet
!!$      read(row(85:96),'(f12.4)')   STD%STD
!!$      read(row(98:109),'(f12.4)')   STD%SZD
!!$      read(row(111:122),'(f12.4)')   STD%SDD
!!$      read(row(124:135),'(f12.4)')  STD%SZDD
!!$      read(row(137:148),'(f12.4)') STD%SWD
!!$      read(row(150:161),'(f12.4)') STD%SZWD
!!$      read(row(163:174),'(f12.4)') STD%SWV
!!$      read(row(176:189),'(f12.4)') STD%SZWV
!!$      !write(*,*) 'ReadESTDRow > ', STD%sdry, STD%swet
!!$
!!$   end if

end if

end subroutine ReadESTDRowarr
!****
! <<<<<<<< End of RoboDoc comments


!---------------------------------------------------------------------
! subroutine ReadSTDRowarr
!---------------------------------------------------------------------
!
! Start of RoboDoc comments >>>>
!****S* GPSImport/ReadSTDRowarr
!
! Name
! ReadSTDRowarr
!
! Purpose
! Extract slant data from one row of ASCII data from the slant file
!
! Call
! call ReadSTDRowarr(row, STD, err)
!
! Input
! row  - row of ASCII data
!
! Output
! STD  - Slant record
! err  - Returncode, err = 0 - no error, err <> 0 - error
!
! External References
! None
!
! Record of revisions
! Date        Programmer   Description of change
! ====        ==========   =====================
! 28.02.2007  M. Bender    new
! 01.04.2009  M. Bender    current version copied to "ReadSTDRowOld",
!                          type of pointer "STD" changed to "LLSTDDataRecord"
! 22.10.2009  M. Bender    update to slant data format ver. 2
! 15.08.2011  M. Bender    station and site: wrong format fixed (i6 => i4 + i2)
! 25.05.2012  M. Bender    Copy of ReadSTDRow,
!                          pointer STD replaced by structure "SlantTotalDelay"
!---------------------------------------------------------------------
! Source
!---------------------------------------------------------------------
subroutine ReadSTDRowarr(row, STD, ver, err)

! List of calling arguments:
character (len=*), intent(in)      :: row
type (SlantTotalDelay)             :: STD
integer, intent(in)                :: ver
integer, intent(out)               :: err

! List of local variables:

!---------------------------------------------------------------------

err = 0
if (len_trim(row) .lt. 56) then
   ! row too short, cannot contain slant data
   err = 1
else
   ! Copy data in slant structure

   if (ver .eq. 1) then
      ! Slant files version 1, EPOS 6, "old" format

      read(row(1:11),'(f11.6)') STD%time
      ! Read station ID (first 4 digits of "STA")
      !read(row(13:18),'(i6)') STD%station
      read(row(13:16),'(i4)') STD%station
      ! Read site ID (last 2 digits of "STA")
      read(row(17:18),'(i2)') STD%site
      read(row(20:22),'(i3)') STD%satellite
      !read(row(24:32),'(f9.4)') STD%elevation
      read(row(24:29),'(f6.2)') STD%elevation
      read(row(31:37),'(f7.2)') STD%azimuth
      read(row(39:47),'(f9.4)') STD%slant
      read(row(49:57),'(f9.4)') STD%zslant

   else
      ! Slant files version 2, EPOS 8, "new" format, "site" removed
      read(row,*)  STD%time, STD%station, STD%satellite, STD%elevation, &
                   STD%azimuth, STD%slant, STD%zslant
   end if
end if

end subroutine ReadSTDRowarr
!****
! <<<<<<<< End of RoboDoc comments

!---------------------------------------------------------------------
! function to_upper
!---------------------------------------------------------------------
!
! Purpose:
!         Convert a string or character to upper case
!          (valid for ASCII or EBCDIC processors)
! Input:
! string     - original string
!
! Output (return value):
! new_string - string with upper case letters only
!
! External References:
!
! Procedure:
!
! Record of revisions:
! Date        Programmer   Description of change
! ====        ==========   =====================
! 05.06.2007  Copy from Ed Akin's book "up_low_f.txt" with some changes
!---------------------------------------------------------------------
function  to_upper (string)  result (new_string) ! like C

! List of calling arguments:
character (len = *), intent(in) :: string      ! unknown length
character (len = len(string))   :: new_string  ! same length

! List of local variables:
character (len=26), parameter   ::               &
            UPPER = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ',  &
            lower = 'abcdefghijklmnopqrstuvwxyz'
integer :: k    ! loop counter
integer :: loc  ! position in alphabet
!---------------------------------------------------------------------

new_string = string       ! copy everything
do k = 1, len(string)     ! to change letters
   loc = index ( lower, string(k:k))   ! find letter in lower
   if (loc /= 0 ) new_string(k:k) = UPPER(loc:loc) ! change
end do ! over string characters

end function to_upper


!---------------------------------------------------------------------
! subroutine Ellips2Cart
!---------------------------------------------------------------------
!
! Start of RoboDoc comments >>>>
!****S* Coord/Ellips2Cart
!
! Name
! Ellips2Cart
!
! Call
! call Ellips2Cart( lambda, beta, height, X, Y, Z, Ellips)
!
! Purpose
! Umrechnung von ellipsoidischen in kartesische Koordinaten
!
! Die ellipsoidischen Koordinaten beziehen sich auf das Referenzellipsoid,
! das durch die Struktur "Ellips" festgelegt wird. Diese Routine kann
! damit f"ur beliebige Ellipsoide (nur Rotationsellipsoide ???) verwendet
! werden, wenn die Grundparameter des Ellipsoids in der Struktur
! "Ellips" gespeichert sind.
!
! R"ucktransformation mit der Routine "Cart2Ellips"
!
! Input
! lambda - ellipsoidische L"ange in rad
! beta   - ellipsoidische Breite in rad
! height - ellipsoidische H"ohe in m
! Ellips - Struktur vom Typ "EllipsParam" mit den Parametern des
!          Referenzellipsoids
!
! Output
! X, Y, Z - kartesische Koordinaten in m
!
! External References
!
! Examples
!
! Record of revisions
! Date        Programmer   Description of change
! ====        ==========   =====================
! 27. 1.2005  M. Bender    neu
! 19.04.2005  M. Bender    Reihenfolge der Parameter ver"andert
!---------------------------------------------------------------------
! Source
!---------------------------------------------------------------------
subroutine Ellips2Cart( lambda, beta, height, X, Y, Z, Ellips)

! List of calling arguments:
real (wp), intent(in) :: lambda, beta, height
real (wp), intent(out) :: X, Y, Z
type (EllipsParam) :: Ellips

! List of local variables:
real (wp) :: NK
!---------------------------------------------------------------------

! Querkr"ummungsradius des Ellipsoids
NK = Ellips%a/sqrt(1.0 - Ellips%EN1sup2*(sin(beta))**2)

! Berechnung der kartesischen Koordinaten
X = (NK+height)*cos(beta)*cos(lambda)
Y = (NK+height)*cos(beta)*sin(lambda)
Z = ((NK/(1+Ellips%EN2sup2)) + height)*sin(beta)

end subroutine Ellips2Cart
!****
! <<<<<<<< End of RoboDoc comments


!---------------------------------------------------------------------
! subroutine Cart2Ellips
!---------------------------------------------------------------------
!
! Start of RoboDoc comments >>>>
!****S* Coord/Cart2Ellips
!
! Name
! Cart2Ellips
!
! Call
! call Cart2Ellips(X, Y, Z, lambda, beta, height, Ellips)
!
! Purpose
! Umrechnung von kartesischen Koordinaten in ellipsoidische Koordinaten
! Direkte Rechnung ohne Iterationsverfahren
!
! Quelle: B. Hofmann-Wellenhof, H. Lichtenegger, J. Collins
!         Global Positioning System - Theory and Practice
!         Springer, Wien, New York, 1993, Seite 232
!         Die Transformationsformel fuer die Laenge wurde leicht abgewandelt.
!
!         Hier fehlt die Sonderbehandlung von X=Y=Z=0, Winkel = 0
!         und die separate Behandlung der Oktanden, so dass wirklich
!         Laengen zwischen 0 und 360 Grad und Breiten zwischen -90 und +90 Grad
!         herauskommen. Dies wurde analog zur Transformation in sphaerische
!         Koordinaten ergaenzt.
!
! Zum Testen dieser Transformation steht das Programm "TestWgs.f90" in
! ../ToolsTest/Coord zur Verfuegung. Dort werden bekannte Stationskoordinaten
! sowohl in kartesischen als auch in WGS84 Koordinaten eingelesen und dann
! mit den Ergebnissen der Transformationen verglichen.
!
! Input
! X, Y, Z - kartesische Koordinaten i m
! Ellips  - Struktur vom Typ "EllipsParam" mit den Daten des Referenzellipsoids
!
! Output
! lambda - ellipsoidische L"ange in rad
! beta   - ellipsoidische Breite in rad
! height - ellipsoidische H"ohe in m
!
! External References
!
! Examples
!
! Record of revisions
! Date        Programmer   Description of change
! ====        ==========   =====================
! 18 7.2005   M. Bender    neu, ersetzt die alte Routine gleichen Namens, die
!                          mit einem Iterations-Verfahren arbeitet
!                          (noch verf"ugbar als "Cart2EllipsIter")
! 19.08.2008  M. Bender    kleinere Aenderungen in den if-Abfragen
!---------------------------------------------------------------------
! Source
!---------------------------------------------------------------------
subroutine Cart2Ellips(X, Y, Z, lambda, beta, height, Ellips)

! List of calling arguments:
real (wp), intent(in) :: X, Y, Z
real (wp), intent(out) :: lambda, beta, height

type (EllipsParam) :: Ellips

! List of local variables:
real (wp) :: rho, phi, gamma, g2
real (wp) :: N, theta, h
!---------------------------------------------------------------------

! Einige Hilfsgr"ossen:
rho = sqrt(X**2+Y**2)

if (rho .gt. reps) then
   theta = atan((Z*Ellips%a)/(rho*Ellips%b))
else
   theta = pi05
end if

! Die ellipsoidische L"ange lambda:
if (abs(rho) > reps) then
   phi = 2.0_wp * atan(Y/(abs(X)+rho))
end if
if (X .ge. 0.0_wp) then
   if (abs(X) .lt. 1.0E-8_wp .and. abs(Y) .lt. 1.0E-8_wp) then
      ! if X=0 and Y=0
      lambda = 0.0_wp
   else if (Y .ge. 0.0_wp) then
      ! X>0 and Y>0
      lambda = phi
   else
      ! X>0 and y<0
      lambda = phi + 2.0_wp*pi
   end if
else  ! X < 0
   lambda = pi - phi
end if

! Die ellipsoidische Breite beta:
g2 = rho - Ellips%EN1sup2*Ellips%a*(cos(theta))**3  ! Nenner von gamma
if (abs(g2) > reps) then
   gamma = atan( (Z + Ellips%EN2sup2*Ellips%b*(sin(theta))**3) / g2 )
! else ???? +pi/2 oder -pi/2 ???
end if
if (abs(rho) .gt. 1.0E-8_wp) then
   ! rho > 0
   beta = gamma
else  ! rho = 0
   if (Z .gt. 0.0_wp) then
      beta = 0.500_wp*pi
   else if (abs(Z) < reps) then
      beta = 0.0_wp
   else if (Z .lt. 0.0_wp) then
      beta = -0.500_wp*pi
   end if
end if

! Bestimmung des Querkr"ummungsradius f"ur diese Breite:
N = Ellips%a/sqrt(1.0_wp - Ellips%EN1sup2*((sin(beta))**2))
!write(*,*) 'N = ', N, (sin(beta))**2

! Die ellipsoidische H"ohe h:
h = (rho/cos(beta)) - N
!write(*,*) 'h = ', h, cos(beta)

if (abs(rho) > reps) then
   height = h
else ! rho = 0
   if (abs(Z) > reps) then
      height = -Ellips%b
   else
      height = abs(Z) - Ellips%b
   end if
end if

end subroutine Cart2Ellips
!****
! <<<<<<<< End of RoboDoc comments


!---------------------------------------------------------------------
! subroutine LocalHorz2Cart
!---------------------------------------------------------------------
!
! Start of RoboDoc comments >>>>
!****S* Coord/LocalHorz2Cart
!
! Name
! LocalHorz2Cart
!
! Call
! call LocalHorz2Cart(x2, y2, z2, lambda, phi,   &
!                     ox, oy, oz, x1, y1, z1)
!
! Purpose
! Transformation from a local horizontal system to a global cartesian system
!
! Global system: right handed cartesian system
! Local horizontal system: left handed cartesian system with its origin
!                          at the position (lambda,phi,alt) of the
!                          observation point, e. g. the GPS station
!
! The transformation is identical for spheres and ellipsoids, see documentation
!
! Input
! x2,y2,z2  - cartesian coordinates in the local horizontal system
! lambda    - longitude of observation point, rad
! phi       - latitude of observation point, rad
! ox,oy,oz  - global cartesian coordinates of the observation point
!             (e. g. GPS station).
!             ox,oy,oz and lambda,phi are not independent values, but
!             are related by the choosen frame of reference
!             (ellipse, sphere...)
!
! Output
! x1,y1,z1  - global cartesian coordinates
!
! External References
! None
!
! Examples
!
! Record of revisions
! Date        Programmer   Description of change
! ====        ==========   =====================
! 22.01.2007  M. Bender    new
!---------------------------------------------------------------------
! Source
!---------------------------------------------------------------------
subroutine LocalHorz2Cart(x2, y2, z2, lambda, phi,   &
                          ox, oy, oz, x1, y1, z1)

! List of calling arguments:
real (wp), intent(in)  :: x2, y2, z2
real (wp), intent(in)  :: lambda, phi
real (wp), intent(in)  :: ox, oy, oz
real (wp), intent(out) :: x1, y1, z1

! List of local variables:
!---------------------------------------------------------------------

x1 = ox -sin(phi)*cos(lambda)*x2 -      &
         sin(lambda)*y2          +      &
         cos(phi)*cos(lambda)*z2
y1 = oy -sin(phi)*sin(lambda)*x2 +      &
         cos(lambda)*y2          +      &
         cos(phi)*sin(lambda)*z2
z1 = oz +cos(phi)*x2             +      &
         sin(phi)*z2

end subroutine LocalHorz2Cart
!****
! <<<<<<<< End of RoboDoc comments


!---------------------------------------------------------------------
! subroutine Azimut2Horz
!---------------------------------------------------------------------
!
! Start of RoboDoc comments >>>>
!****S* Coord/Azimut2Horz
!
! Name
! Azimut2Horz
!
! Call
! call Azimut2Horz(A, Z, dist, x2, y2, z2)
!
! Purpose
! Compute the cartesian coordinates in the local horizontal system of an
! object with given azimuth A, zenith distance Z and distance D
!
! Local horizontal system: left handed cartesian system with its origin
!                          at the position (lambda,phi) on the surface
!                          of the reference ellipsoid or sphere.
!
! Input
! A    -  azimuth of the object (angle in rad)
!         A = [0, ..., 2*pi] or
!         A = [-pi, ..., pi]
! Z    -  zenith distance (angle in rad)
!         Z = [0, ..., pi/2]
! dist - distance between the object and the origin of the local
!        horizontal system
!
! Output
! x2,y2,z2  - kartesian coordinates in the local horizontal system
!
! External References
! None
!
! Examples
!
! Record of revisions
! Date        Programmer   Description of change
! ====        ==========   =====================
! 12.05.2005  M. Bender    neu
!---------------------------------------------------------------------
! Source
!---------------------------------------------------------------------
subroutine Azimut2Horz(A, Z, dist, x2, y2, z2)

! List of calling arguments:
real (wp), intent(in)  :: A, Z, dist
real (wp), intent(out) :: x2, y2, z2

! List of local variables:
!---------------------------------------------------------------------

x2 = dist * cos(A) * sin(Z)
y2 = dist * sin(A) * sin(Z)
z2 = dist * cos(Z)

end subroutine Azimut2Horz
!****
! <<<<<<<< End of RoboDoc comments


!---------------------------------------------------------------------
! function NWein
!---------------------------------------------------------------------
!
! Start of RoboDoc comments >>>>
!****F* WVField/NWein
!
! Name
! NWein
!
! Call
! N = NWein( P, T, E )
!
! Purpose
! Refraktivitaet der Atmosphaere nach Smith und Weintraub:
!
! N = k1*(p/T) +(k2-k1)*(e/T) + k3*(e/(T*T))
!
! Brechungsindex n:  N = 10^6 (n-1)
!
! k1, k2, k3 - empirisch bestimmte Konstanten, definiert im Header
! p          - Druck in hPa
! T          - Temperatur in Kelvin
! e          - Dampfdruck des Wasserdampfes in hPa
! WARNUNG: Das Programm erwartet den Druck in Pa, nicht in hPa !!!
!
! Input
! p - Druck in Pa
! T - Temperatur in K
! e - Partialdruck des Wasserdampfes in Pa
!
! Output
! NWein - Refraktivitaet N der Atmosphaere
!
! External References
! None
!
! Record of revisions
! Date        Programmer   Description of change
! ====        ==========   =====================
! 29.01.2007  M. Bender    new
!---------------------------------------------------------------------
! Source
!---------------------------------------------------------------------
function NWein( P, T, E )

real (wp)   :: NWein

! List of calling arguments:
real (wp),  intent(in)   ::  T, P, E

! List of local variables:
real (wp)  ::  hP, hE
!---------------------------------------------------------------------

! Umrechnen der Druecke in hPa:
hP = P/100.0_wp
hE = E/100.0_wp

NWein = k1*(hP/T) + (k2 - k1)*(hE/T) + k3*(hE/(T*T))

end function NWein
!****
! <<<<<<<< End of RoboDoc comments


!---------------------------------------------------------------------
! function NdWein
!---------------------------------------------------------------------
!
! Start of RoboDoc comments >>>>
!****F* WVField/NdWein
!
! Name
! NdWein
!
! Call
! Nd = NdWein( Pd, T )
!
! Purpose
! Refraktivitaet der trockenen Atmosphaere nach Smith und Weintraub:
!
! Nd = k1*(pd/T)
!
! Brechungsindex n:  N = 10^6 (n-1)
!
! k1, k2, k3 - empirisch bestimmte Konstanten, definiert im Header
! pd         - Partialdruck der trockenen Luft in hPa
! T          - Temperatur in Kelvin
! WARNUNG: Das Programm erwartet den Druck in Pa, nicht in hPa !!!
!
! Input
! P  - Luftdruck in Pa
! T  - Temperatur in K
! e  - Partialdruck des Wasserdampfes in Pa
!
! Output
! NdWein - Refraktivitaet Nd der trockenen Atmosphaere
!
! External References
! None
!
! Record of revisions
! Date        Programmer   Description of change
! ====        ==========   =====================
! 01.02.2007  M. Bender    new
! 05.08.2015  M. Bender    parameters are now P, T, E
!---------------------------------------------------------------------
! Source
!---------------------------------------------------------------------
function NdWein( P, T, E )

real (wp)   :: NdWein

! List of calling arguments:
real (wp),  intent(in)   ::  P, T, E

! List of local variables:
real (wp)  ::  hP, hE
!---------------------------------------------------------------------

! convert pressure to hPa:
hP = P/100.0_wp
hE = E/100.0_wp

NdWein = k1*((hP-hE)/T)

end function NdWein
!****
! <<<<<<<< End of RoboDoc comments


!---------------------------------------------------------------------
! function NwWein
!---------------------------------------------------------------------
!
! Start of RoboDoc comments >>>>
!****F* WVField/NwWein
!
! Name
! NwWein
!
! Call
! Nw = NwWein( e, T )
!
! Purpose
! Refraktivitaet des Wasserdampfes nach Smith und Weintraub:
!
! Nw = k2*(e/T) + k3*(e/(T*T))
!
! Brechungsindex n:  N = 10^6 (n-1)
!
! k1, k2, k3 - empirisch bestimmte Konstanten, definiert im Header
! e          - Partialdruck des Wasserdampfes in hPa
! T          - Temperatur in Kelvin
! WARNUNG: Das Programm erwartet den Druck in Pa, nicht in hPa !!!
!
! Input
! e - Partialdruck des Wasserdampfes in Pa
! T - Temperatur in K
!
! Output
! NwWein - Refraktivitaet Nd des Wasserdampfes
!
! External References
! None
!
! Record of revisions
! Date        Programmer   Description of change
! ====        ==========   =====================
! 01.02.2007  M. Bender    new
!---------------------------------------------------------------------
! Source
!---------------------------------------------------------------------
function NwWein( e, T )

real (wp)   :: NwWein

! List of calling arguments:
real (wp),  intent(in)   ::  E, T

! List of local variables:
real (wp)  ::  hE
!---------------------------------------------------------------------

! convert pressure to hPa:
hE = E/100.0_wp

NwWein = k2*(hE/T) + k3*(hE/(T*T))

end function NwWein
!****
! <<<<<<<< End of RoboDoc comments


!---------------------------------------------------------------------
! function ZHDsaas
!---------------------------------------------------------------------
!
!> Compute the ZHD using the Saastamoinen formula
!>
!> <b> ZHD = ZHDsaas(P, lat, height) </b>
!>
!> J. Saastamoinen derived an approximation for the ZHD above a station
!> at a given latitude. The pressure at the station (surface prssure), the
!> height above ellipsoid of the station and its latitude are required to
!> estimate the ZHD. As the "station" might be at any height, the ZHD above
!> any height can be estimated, e.g. the ZHD above the top level of a
!> numerical weather model.
!>
!> References \n
!> Elgered, G.; Davis, J. L.; Herring, T. A. & Shapiro, I. I. \n
!> Geodesy by radio interferometry:
!> Water vapor radiometry for estimation of the wet delay \n
!> J. Geophys. Res., 1991, vol. 96, pp. 6541-6555 \n \n
!>
!> Bevis, M.; Businger, S.; Herring, T. A.; Rocken, C.; Anthes, R. A. &
!> Ware, R. H. \n
!> GPS Meteorology: Remote Sensing of Atmospheric Water Vapor
!> Using the Global Positioning System \n
!> J. Geophys. Res., 1992, vol. 97, pp. 15787-15801
!>
!> @param[in] P       surface pressure in Pa
!> @param[in] lat     latitude in radian
!> @param[in] height  station height above ellipsoid in m
!> @return    ZHDsaas ZHD above height in m
!
!---------------------------------------------------------------------
! Record of revisions
! Date        Programmer   Description of change
! ====        ==========   =====================
! 05.08.2015  M. Bender    new
!---------------------------------------------------------------------
function ZHDsaas( P, lat, height )

real (wp)   :: ZHDsaas

! List of calling arguments:
real (wp),  intent(in)   :: P, lat, height

! List of local variables:
real (wp)  ::  A
!---------------------------------------------------------------------

A = 1.0_wp - 0.00266_wp * cos(2*lat) - height*28.0E-7_wp

ZHDsaas = P*1.0E-2_wp*2.2779E-3_wp

end function ZHDsaas


!---------------------------------------------------------------------
! subroutine CrossLineEllips
!---------------------------------------------------------------------
!
! Start of RoboDoc comments >>>>
!****S* Grid/CrossLineEllips
!
! Name
! CrossLineEllips
!
! Call
! call CrossLineEllips ( RefEllips, Alt, SLine, CrossPt, lambda )
!
! Purpose
! Schnittpunkt einer Geraden mit einem Ellipsoid
!
! Hier wird ein Rotations-Ellipsoid vorausgesetzt, dessen grosse Halbachse(n)
! senkrecht zur Rotationsachse stehen:
! z - Rotationsachse
! x - Die x-Achse geht durch den Nullmeridian bei Greenwich
! y - Die y-Achse wird so erg"anzt, dass ein rechtshaendiges System entsteht
!
! Ellipse: x**2/a**2 + y**2/a**2 + z**2/b**2 = 1
!          a - grosse Halbachse
!          b - kleine Halbachse, parallel zur z-Achse
!
! Input
! RefEllips - Referenz-Ellipsoid
! Alt       - H"ohe "uber dem Referenz-Ellipsoid
! Line      - Gerade, deren Schnittpukt mit dem Ellipsoid berechnet werden soll
!
! Output
! CrossPt - Schnittpunkt Gerade - Ellipsoid in kartesischen Koordinaten
! lambda    - L"ange der Geraden, dabei wird angenommen, dass der Startpunkt
!             von "SLine" die GPS-Bodenstation ist. Lambda ist dann die L"ange
!             des Strahls bis zu einer Hoehe von "Alt".
! exist     - Gibt es einen Schnittpunkt zwischen der Geraden und dem
!             Ellipsoid?  exist = .true.  - Schnittpunkt existiert
!                         exist = .false. - Schnittpunkt existiert nicht
!
! External References
! None
!
! Record of revisions
! Date        Programmer   Description of change
! ====        ==========   =====================
! 19.02.2007  M. Bender    neu
! 31.01.2011  M. Bender    Optionalen Parameter "exist" eingef"uhrt
!---------------------------------------------------------------------
! Source
!---------------------------------------------------------------------
subroutine CrossLineEllips ( RefEllips, Alt, SLine, CrossPt, lambda, exist )

! List of calling arguments:
type (EllipsParam), intent(in)                    :: RefEllips
real (wp), intent(in)                  :: Alt
type (Line), intent(in)                           :: SLine
real (wp), dimension(1:3), intent(out) :: CrossPt
real (wp), intent(out)                 :: lambda
logical, optional,  intent(out)                   :: exist

! List of local variables:
real (wp) :: ha, hb
real (wp) :: a0, a1, a2
real (wp) :: lambda1, lambda2
logical              :: IsReal
!---------------------------------------------------------------------

! Hier wird nicht der Schnittpunkt der Geraden mit dem Referenz-Ellipsoid,
! sondern mit einem Ellipsoid in der H"ohe "Alt" gesucht. Die entsprechenden
! Halbachsen sind:
ha = RefEllips%a + Alt
hb = RefEllips%b + Alt

! Die beiden Schnittpunkte der Geraden mit der Ellipse ergeben sich als L"osung
! einer quadratischen Gleichung: a0 + a1*lambda + a2*lambda**2 = 0
a0 = (SLine%start(1)**2/ha**2) + (SLine%start(2)**2/ha**2) + &
     (SLine%start(3)**2/hb**2) - 1.0
a1 = 2*( (SLine%start(1)*SLine%unitvec(1)/ha**2) +  &
         (SLine%start(2)*SLine%unitvec(2)/ha**2) +  &
         (SLine%start(3)*SLine%unitvec(3)/hb**2)   )
a2 = (SLine%unitvec(1)/ha)**2 +  &
     (SLine%unitvec(2)/ha)**2 +  &
     (SLine%unitvec(3)/hb)**2

call QuadEq(a1/a2, a0/a2, lambda1, lambda2, IsReal)

if (IsReal) then
   ! Es gibt Schnittpunkte
   if (lambda1 .ge. lambda2) then
      ! Den n"aher gelegenen Punkt ausw"ahlen
      call LinePos(SLine, lambda1, CrossPt)
      lambda = lambda1
   else
      call LinePos(SLine, lambda2, CrossPt)
      lambda = lambda2
   end if
else
   ! Kein Schnittpunkt ??? Das sollte es nicht geben
   write(*,*) 'CrossLineEllips> Error: No intersecting point found'
   stop
end if

if (present(exist)) then
   exist = IsReal
end if

!write(*,*) IsReal, lambda1, lambda2

end subroutine CrossLineEllips
!****
! <<<<<<<< End of RoboDoc comments


!---------------------------------------------------------------------
! function BiLin2D
!---------------------------------------------------------------------
!
! Start of RoboDoc comments >>>>
!****F* Interpolation/BiLin2D
!
! Name
! BiLin2D
!
! Purpose
! 2D interpolation usig a bilinear interpoltion between four reference
! values.
!
!
!    y ^
!      |
!      |   node3                  node4             y3 = y4, x3 < x4
!      |
!      |          RfPt
!      |
!      |   node1                  node2             y1 = y2, x1 < x2
!      |
!      ---------------------------------------> x
!
!          x1 = x3                x2 = x4
!          y1 < y3                y2 < y4
!
! Test program: .../gpstomo/trunk/Tomo/TomoTest/Interpol/TestBiLin2D.f90
!
!
! Call
! val =  BiLin2D(RefPt, node)
!
! Input
! RefPt    - reference point, the field shall be interpolated at this
!            point, RefPoint(1) - longitude
!                   RefPoint(2) - latitude
!                   RefPoint(3) - altitude
! node     - array
!            node(i,:), i=1, ... , 4 : 4 nodes at the cell corner
!            node(i,j), j=1, 2, 3    : j = 1 - x
!                                      j = 2 - y
!                                      j = 3 - f_j = f(x,y)
!
! Output
! BiLin2D (return value) - interpolated value
!
! External References
! None
!
! Record of revisions
! Date        Programmer   Description of change
! ====        ==========   =====================
! 17.09.2009  M. Bender    new
!---------------------------------------------------------------------
! Source
!---------------------------------------------------------------------
function BiLin2D(RefPt, node)

! List of calling arguments:
real (wp)                                 :: BiLin2D
real (wp), dimension(1:2), intent(in)     :: RefPt
real (wp), dimension(1:4,1:3), intent(in) :: node

! List of local variables:
real (wp) :: x, y
!---------------------------------------------------------------------

! Normalized x and y coordinates:
x = (RefPt(1)-node(1,1)) / (node(2,1)-node(1,1))
y = (RefPt(2)-node(1,2)) / (node(3,2)-node(1,2))

!write(*,*) x, y

! Interpolation
BiLin2D  = node(1,3) * (1.0_wp-x)*(1.0_wp-y) +  &
           node(2,3) *     x    *(1.0_wp-y) +  &
           node(3,3) * (1.0_wp-x)*    y     +  &
           node(4,3) *     x    *    y

end function BiLin2D
!****
! <<<<<<<< End of RoboDoc comments


!---------------------------------------------------------------------
! function BiLinear2D
!---------------------------------------------------------------------
!
! Start of RoboDoc comments >>>>
!****F* Interpolation/BiLinear2D
!
! Name
! BiLinear2D
!
! Call
! val =  BiLinear2D(RefPt, node)
!
! Purpose
! 2D interpolation usig a bilinear interpoltion between four reference
! values. Any set of nodes is possible as long as (y2-y1) <> 0 and
! (x2-x1) <> 0 and (x4-x3) <> 0, i. e. 3 or 4 points are on a straoght line.
! If only 3 reference points ara available, the routine can be called with
! P4 = P2 (node 4 = node 2).
!
!
!    y ^
!      |
!      |   node3                  node4             y3 > y4, x3 < x4
!      |
!      |          RfPt
!      |
!      |   node1                  node2             y1 > y2, x1 < x2
!      |
!      ---------------------------------------> x
!
!          x1 < x2                x2 > x1
!          y1 < y3                y2 < y4
!
! Test program: .../gpstomo/trunk/Tomo/TomoTest/Interpol/TestBiLinear2D.f90
!
!

! Input
! RefPt    - reference point, the field shall be interpolated at this
!            point, RefPoint(1) - x, longitude
!                   RefPoint(2) - y, latitude
! node     - array
!            node(i,:), i=1, ... , 4 : 4 nodes at the cell corner
!            node(i,j), j=1, 2, 3    : j = 1 - x
!                                      j = 2 - y
!                                      j = 3 - f_j = f(x,y)
!
! Output
! BiLinear2D (return value) - interpolated value
!
! External References
! None
!
! Record of revisions
! Date        Programmer   Description of change
! ====        ==========   =====================
! 10.12.2012  M. Bender    new
!---------------------------------------------------------------------
! Source
!---------------------------------------------------------------------
function BiLinear2D(RefPt, node)

! List of calling arguments:
real (wp)                                 :: BiLinear2D
real (wp), dimension(1:2), intent(in)     :: RefPt
real (wp), dimension(1:4,1:3), intent(in) :: node

! List of local variables:
real (wp) :: c1, c2, d1, d2, d3, d4
real (wp) :: ya, yb, fa, fb, dy, a1, a2, a3
!---------------------------------------------------------------------

! c1 = x2 - x1
c1 = node(2,1) - node(1,1)
! c2 = x4 - x3
c2 = node(4,1) - node(3,1)
! d1 = (x - x1)/(x2-x1) =  (x - x1)/c1
d1 = (RefPt(1) - node(1,1))/c1
! d2 = (x2 - x)/(x2-x1) =  (x2 - x)/c1
d2 = (node(2,1) - RefPt(1))/c1
! d3 = (x - x3)/(x4-x3) =  (x - x3)/c2
d3 = (RefPt(1) - node(3,1))/c2
! d4 = (x4 - x)/(x4-x3) =  (x4 - x)/c2
d4 = (node(4,1) - RefPt(1))/c2
!write(*,*) 'BiLinear2D> c1, c2, d1, d2, d3, d4 : ', c1, c2, d1, d2, d3, d4

fa = d2*node(1,3) + d1*node(2,3)
fb = d4*node(3,3) + d3*node(4,3)
ya = d2*node(1,2) + d1*node(2,2)
yb = d4*node(3,2) + d3*node(4,2)
dy = yb - ya

if (abs(dy) .lt. 1.0E-9_wp .and. abs(fb-fa) .lt. 1.0E-7_wp) then
   ! x = x2, singularity, should happen only for 3 point interpolation
   !write(*,*) 'BiLinear2D> x = x2'

   a1 = - node(2,1)*node(3,2)*node(2,3)        &
        + node(1,1)*node(3,2)*node(2,3)        &
        - node(2,1)*node(2,2)*node(1,3)        &
        + node(3,1)*node(2,2)*node(1,3)        &
        + node(2,1)*node(1,2)*node(2,3)        &
        - node(3,1)*node(1,2)*node(2,3)        &
        + node(2,1)*node(2,2)*node(3,3)        &
        - node(1,1)*node(2,2)*node(3,3)

   a2 =  RefPt(2) * (   node(2,1)*node(1,3)    &
                      - node(3,1)*node(1,3)    &
                      + node(3,1)*node(2,3)    &
                      - node(2,1)*node(3,3)    &
                      + node(1,1)*node(3,3)    &
                      - node(1,1)*node(2,3)  )

   a3 =   node(2,1)*node(3,2)    &
        - node(1,1)*node(3,2)    &
        + node(1,1)*node(2,2)    &
        - node(2,1)*node(1,2)    &
        + node(3,1)*node(1,2)    &
        - node(3,1)*node(2,2)

   BiLinear2D = - (a1 + a2) / a3

else
   BiLinear2D = fa*(yb- RefPt(2))/dy + fb*(RefPt(2)-ya)/dy
end if

end function BiLinear2D
!****
! <<<<<<<< End of RoboDoc comments


!---------------------------------------------------------------------
! function BiLinearDiff2D
!---------------------------------------------------------------------
!
! Start of RoboDoc comments >>>>
!****F* Interpolation/BiLinearDiff2D
!
! Name
! BiLinearDiff2D
!
! Call
! val =  BiLinearDiff2D(RefPt, node)
!
! Purpose
! 2D interpolation usig a bilinear interpoltion between four reference
! values. Any set of nodes is possible as long as (y2-y1) <> 0 and
! (x2-x1) <> 0 and (x4-x3) <> 0, i. e. 3 or 4 points are on a straoght line.
! If only 3 reference points ara available, the routine can be called with
! P4 = P2 (node 4 = node 2).
!
! In case of 3 reference points there is a singularity at  x = x2
! i.e. RefPt(1) = x2, which is replaced by the mean value
! od two points close to the singularity which can be evaluated.
!
!
!    y ^
!      |
!      |   node3                  node4             y3 > y4, x3 < x4
!      |
!      |          RfPt
!      |
!      |   node1                  node2             y1 > y2, x1 < x2
!      |
!      ---------------------------------------> x
!
!          x1 < x2                x2 > x1
!          y1 < y3                y2 < y4
!
! Test program: .../gpstomo/trunk/Tomo/TomoTest/Interpol/TestBiLinear2D.f90
!
!

! Input
! RefPt    - reference point, the field shall be interpolated at this
!            point, RefPoint(1) - x, longitude
!                   RefPoint(2) - y, latitude
! node     - array
!            node(i,:), i=1, ... , 4 : 4 nodes at the cell corner
!            node(i,j), j=1, 2, 3    : j = 1 - x
!                                      j = 2 - y
!                                      j = 3 - f_j = f(x,y)
!
! Output
! BiLinearDiff2D (return value) - interpolated value
!
! External References
! None
!
! Record of revisions
! Date        Programmer   Description of change
! ====        ==========   =====================
! 17.05.2013  M. Bender    new
!---------------------------------------------------------------------
! Source
!---------------------------------------------------------------------
function BiLinearDiff2D(RefPt, nodein, Nnode)

! List of calling arguments:
real (wp)                                 :: BiLinearDiff2D
real (wp), dimension(1:2), intent(in)     :: RefPt
real (wp), dimension(1:4,1:3), intent(in) :: nodein
integer,  intent(in), optional            :: Nnode

! List of local variables:
real (wp), dimension(1:4,1:3) :: node
real (wp) :: c1, c2, d1, d2, d3, d4
real (wp) :: ya, yb, fa, fb, dy !, a1, a2, a3
!real (wp) :: epsilon, I1, I2
!---------------------------------------------------------------------

node = nodein

if (present(Nnode)) then
   if (Nnode == 3) then
      !if (abs(lon(1)-lon(2)) < 1.0E-6_wp) then
      if (abs(node(1,1)-node(2,1)) < 1.0E-6_wp) then
         ! change 2 and 3 to avoid singularity, copy 3 to 4
         node(4,:) = node(3,:)
         node(3,:) = node(2,:)
         node(2,:) = node(4,:)
      !else if (abs(lon(2)-lon(3)) < 1.0E-6_wp) then
      else if (abs(node(2,1)-node(3,1)) < 1.0E-6_wp) then
         ! copy 1 to 4 to avoid singularity
         node(4,:) = node(1,:)
      else
         node(4,:) = node(2,:)
      end if
   end if
end if

! c1 = x2 - x1
c1 = node(2,1) - node(1,1)
! c2 = x4 - x3
c2 = node(4,1) - node(3,1)
! d1 = (x - x1)/(x2-x1) =  (x - x1)/c1
d1 = (RefPt(1) - node(1,1))/c1
! d2 = (x2 - x)/(x2-x1) =  (x2 - x)/c1
d2 = (node(2,1) - RefPt(1))/c1
! d3 = (x - x3)/(x4-x3) =  (x - x3)/c2
d3 = (RefPt(1) - node(3,1))/c2
! d4 = (x4 - x)/(x4-x3) =  (x4 - x)/c2
d4 = (node(4,1) - RefPt(1))/c2
!write(*,*) 'BiLinearDiff2D> c1, c2, d1, d2, d3, d4 : ', c1, c2, d1, d2, d3, d4

fa = d2*node(1,3) + d1*node(2,3)
fb = d4*node(3,3) + d3*node(4,3)
ya = d2*node(1,2) + d1*node(2,2)
yb = d4*node(3,2) + d3*node(4,2)
dy = yb - ya

!!$if (abs(dy) .lt. 1.0E-7_wp .and. abs(fb-fa) .lt. 1.0E-7_wp) then
!!$   ! x = x2, singularity, should happen only for 3 point interpolation
!!$   ! The singularity is replaced by the mean value of two points close
!!$   ! to the singularity, i.e. computed for RefPt(1) + epsilon and
!!$   ! RefPt(1) - epsilon.
!!$   !write(*,*) 'BiLinearDiff2D> x = x2'
!!$
!!$   !write(*,*) 'Singularitaet: dy = ', dy
!!$
!!$   ! Assuming latitude and longitude given in radian, a difference
!!$   ! of 9*10^-5 rad is equivalent to a distance of ~ 10 m,
!!$   ! i. e. epsilon = 1.0D-5 should be very close to the singularity
!!$   ! and the error of the interpolation is < 0.1 %
!!$   epsilon = 1.0E-5_wp
!!$
!!$   ! x + epsilon (RefPt(1) + epsilon)
!!$   ! d1 = (x - x1)/(x2-x1) =  (x - x1)/c1
!!$   d1 = (RefPt(1)+epsilon - node(1,1))/c1
!!$   ! d2 = (x2 - x)/(x2-x1) =  (x2 - x)/c1
!!$   d2 = (node(2,1) - RefPt(1)-epsilon)/c1
!!$   ! d3 = (x - x3)/(x4-x3) =  (x - x3)/c2
!!$   d3 = (RefPt(1)+epsilon - node(3,1))/c2
!!$   ! d4 = (x4 - x)/(x4-x3) =  (x4 - x)/c2
!!$   d4 = (node(4,1) - RefPt(1)-epsilon)/c2
!!$
!!$   fa = d2*node(1,3) + d1*node(2,3)
!!$   fb = d4*node(3,3) + d3*node(4,3)
!!$   ya = d2*node(1,2) + d1*node(2,2)
!!$   yb = d4*node(3,2) + d3*node(4,2)
!!$   dy = yb - ya
!!$
!!$   I1 = fa*(yb- RefPt(2))/dy + fb*(RefPt(2)-ya)/dy
!!$
!!$   ! x - epsilon (RefPt(1) - epsilon)
!!$   ! d1 = (x - x1)/(x2-x1) =  (x - x1)/c1
!!$   d1 = (RefPt(1)-epsilon - node(1,1))/c1
!!$   ! d2 = (x2 - x)/(x2-x1) =  (x2 - x)/c1
!!$   d2 = (node(2,1) - RefPt(1)+epsilon)/c1
!!$   ! d3 = (x - x3)/(x4-x3) =  (x - x3)/c2
!!$   d3 = (RefPt(1)-epsilon - node(3,1))/c2
!!$   ! d4 = (x4 - x)/(x4-x3) =  (x4 - x)/c2
!!$   d4 = (node(4,1) - RefPt(1)+epsilon)/c2
!!$
!!$   fa = d2*node(1,3) + d1*node(2,3)
!!$   fb = d4*node(3,3) + d3*node(4,3)
!!$   ya = d2*node(1,2) + d1*node(2,2)
!!$   yb = d4*node(3,2) + d3*node(4,2)
!!$   dy = yb - ya
!!$
!!$   I2 = fa*(yb- RefPt(2))/dy + fb*(RefPt(2)-ya)/dy
!!$
!!$   BiLinearDiff2D = 0.50_wp * (I1 + I2)
if (abs(dy) .lt. 1.0E-8_wp) then
   ! slant%Nnghb = 3 : reference point is very close to one node
!!$   write(*,*) 'Singularitaet dy, c1, c2, d1, d2, d3, d4, fa, fb, ya, yb = ', &
!!$        dy, c1, c2, d1, d2, d3, d4,  fa, fb, ya, yb
!!$   if (abs(RefPt(2)-node(1,2)) < 1.0E-6_wp) then
!!$      BiLinearDiff2D = node(1,3)
!!$   else if (abs(RefPt(2)-node(2,2)) < 1.0E-6_wp) then
!!$      BiLinearDiff2D = node(2,3)
!!$   else if (abs(RefPt(2)-node(3,2)) < 1.0E-6_wp) then
!!$      BiLinearDiff2D = node(3,3)
!!$   else if (abs(RefPt(2)-node(4,2)) < 1.0E-6_wp) then
!!$      BiLinearDiff2D = node(4,3)
!!$   else
      BiLinearDiff2D = node(2,3)
   !end if
else
   ! No singularity, interpolate ...
   BiLinearDiff2D = fa*(yb- RefPt(2))/dy + fb*(RefPt(2)-ya)/dy
end if

end function BiLinearDiff2D


!---------------------------------------------------------------------
! function ExpInt1D
!---------------------------------------------------------------------
!
! Start of RoboDoc comments >>>>
!****F* Interpolation/ExpInt1D
!
! Name
! ExpInt1D
!
! Purpose
! One dimensional exponential interpolation using
! N(z) = N2 * (N1/N2)*exp((h2-z)/(h2-h1))
!
! To be used for, e. g., the vertical exponential interpolation of the
! atmospheric refractivity.
!
! Call
! val = ExpInt1D(RefPt, node)
!
! Input
!
! RefPt - scalar value, the values given in "node" will be interpolated
!         at this point.
! node  - array containing the coordinates and the values at these points
!            node(i,:), i=1, 2 : 2 nodes at the end of the interval
!            node(i,j), j=1, 2    : j = 1 - x, coordinate
!                                   j = 2 - f(x), value, observation, ...
!
! Output
! BiLin2D (return value) - interpolated value
!
! External References
! None
!
! Record of revisions
! Date        Programmer   Description of change
! ====        ==========   =====================
! 17.09.2009  M. Bender    new
!---------------------------------------------------------------------
! Source
!---------------------------------------------------------------------
function ExpInt1D(RefPt, node)

! List of calling arguments:
real (wp)                                    :: ExpInt1D
real (wp), intent(in)                        :: RefPt
real (wp), dimension(1:2,1:2), intent(inout) :: node

! List of local variables:
!real (wp) :: x, y
!---------------------------------------------------------------------

if (abs(node(1,2)) .lt. 1E-10_wp .and. abs(node(2,2)) .lt. 1E-10_wp ) then
   ! N1 = N2 = 0  => N(z) = 0

   ExpInt1D = 0.0_wp

else if (abs(node(2,1)-node(1,1)) < 1.0E-3_wp) then
   ! identical coordinates, interpolation not possible, return first value
   ExpInt1D = node(1,2)

else
   ! interpolate

   if (abs(node(1,2)) .lt. 1E-10_wp) then
      ! N1 = 0, N2 <> 0
      node(1,2) = node(2,2)/100.0_wp
   end if
   if (abs(node(2,2)) .lt. 1E-10_wp) then
      ! N2 = 0, N1 <> 0
      node(2,2) = node(1,2)/100.0_wp
   end if

   ExpInt1D = node(2,2) * (node(1,2)/node(2,2))**                     &
                          ((node(2,1)-RefPt)/(node(2,1)-node(1,1)))

end if

end function ExpInt1D
!****
! <<<<<<<< End of RoboDoc comments


!---------------------------------------------------------------------
! function BiLinExp3D
!---------------------------------------------------------------------
!
! Start of RoboDoc comments >>>>
!****F* Interpolation/BiLinExp3D
!
! Name
! BiLinExp3D
!
! Purpose
! 3D interpolation usig a bilinear interpoltion horizontally and an
! exponetial interpolation vertically.
!
! This function makes use of "ExpInt1D" and "BiLin2D".
!
!
! Call
! val =  BiLinExp3D(RefPt, node)
!
! Input
! RefPt    - reference point, the field shall be interpolated at this
!            point, RefPoint(1) - longitude
!                   RefPoint(2) - latitude
!                   RefPoint(3) - altitude
! node     - array containing the coordinates and the values to be interpolated
!              node(i,j,k) - i=1, ... , 4 : 4 nodes at the cell corner
!                            j=1, 2       : grid level, j=1 - lower level
!                                                       j=2 - upper level
!                            k=1, ... , 4 : k=1 - longitude
!                                           k=2 - latitude
!                                           k=3 - altitude
!                                           k=4 - value, refractivity, ...
! Output
! BiLinExp3D (return value) - interpolated value
!
! External References
! None
!
! Record of revisions
! Date        Programmer   Description of change
! ====        ==========   =====================
! 18.09.2009  M. Bender    new
!---------------------------------------------------------------------
! Source
!---------------------------------------------------------------------
function BiLinExp3D(RefPt, node)

! List of calling arguments:
real (wp)                                     :: BiLinExp3D
real (wp), dimension(1:3), intent(in)         :: RefPt
real (wp), dimension(1:4,1:2,1:4), intent(in) :: node

! List of local variables:
real (wp), dimension(1:2,1:2) :: Vnode
real (wp)                     :: VPt
real (wp), dimension(1:4,1:3) :: Hnode
real (wp), dimension(1:2)     :: HPt

integer :: i

!---------------------------------------------------------------------

! Vertical interpolation using "ExpInt1D"
VPt = RefPt(3)  ! interpolate to this altitude
do i=1, 4

   Vnode(1,1) = node(i,1,3)  ! altitude node i, lower level
   Vnode(2,1) = node(i,2,3)  ! altitude node i, upper level

   Vnode(1,2) = node(i,1,4)  ! refractivity node i, lower level
   Vnode(2,2) = node(i,2,4)  ! refractivity node i, upper level

   Hnode(i,3) = ExpInt1D(VPt, Vnode)  ! interpolated vertical refrac.

end do

! Horizontal interpolation using BiLin2D
HPt = RefPt(1:2)  ! interpolate to this (x,y)-coordinate

! x,y coordinates
Hnode(1,1) = node(1,1,1)   ! node1, x
Hnode(1,2) = node(1,1,2)   ! node1, y
Hnode(2,1) = node(2,1,1)   ! node2, x
Hnode(2,2) = node(2,1,2)   ! node2, y
Hnode(3,1) = node(3,1,1)   ! node3, x
Hnode(3,2) = node(3,1,2)   ! node3, y
Hnode(4,1) = node(4,1,1)   ! node4, x
Hnode(4,2) = node(4,1,2)   ! node4, y

BiLinExp3D = BiLin2D(HPt, Hnode)

end function BiLinExp3D
!****
! <<<<<<<< End of RoboDoc comments


!---------------------------------------------------------------------
! function BiLinearExp3D
!---------------------------------------------------------------------
!
! Start of RoboDoc comments >>>>
!****F* Interpolation/BiLinearExp3D
!
! Name
! BiLinearExp3D
!
! Purpose
! 3D interpolation usig a bilinear interpoltion horizontally and an
! exponetial interpolation vertically.
!
! This function makes use of "ExpInt1D" and "BiLinear2D".
!
!
! Call
! val =  BiLinearExp3D(RefPt, node)
! val =  BiLinearExp3D(RefPt, node, Npt)
!
! Input
! RefPt    - reference point, the field shall be interpolated at this
!            point, RefPoint(1) - longitude
!                   RefPoint(2) - latitude
!                   RefPoint(3) - altitude
! node     - array containing the coordinates and the values to be interpolated
!              node(i,j,k) - i=1, ... , 4 : 4 nodes at the cell corner
!                            j=1, 2       : grid level, j=1 - lower level
!                                                       j=2 - upper level
!                            k=1, ... , 4 : k=1 - longitude
!                                           k=2 - latitude
!                                           k=3 - altitude
!                                           k=4 - value, refractivity, ...
! Npt      - number of reference points, Npt = 3 or Npt = 4
!            optional
!
! Output
! BiLinearExp3D (return value) - interpolated value
!
! External References
! None
!
! Record of revisions
! Date        Programmer   Description of change
! ====        ==========   =====================
! 11.12.2012  M. Bender    new, copy of  BiLinExp3D
!---------------------------------------------------------------------
! Source
!---------------------------------------------------------------------
function BiLinearExp3D(RefPt, node, Npt)

! List of calling arguments:
real (wp)                                     :: BiLinearExp3D
real (wp), dimension(:), intent(in)           :: RefPt(1:3)
real (wp), dimension(1:4,1:2,1:4), intent(in) :: node
integer, optional                                        :: Npt

! List of local variables:
real (wp), dimension(1:2,1:2) :: Vnode
real (wp)                     :: VPt
real (wp), dimension(1:4,1:3) :: Hnode
real (wp), dimension(1:2)     :: HPt

integer :: p, i

!---------------------------------------------------------------------

! Number of reference points: 3 or 4
p = 4
if (present(Npt)) then
   if (Npt .eq. 3) p = 3
end if

! Vertical interpolation using "ExpInt1D"
VPt = RefPt(3)  ! interpolate to this altitude
do i=1, p

   Vnode(1,1) = node(i,1,3)  ! altitude node i, lower level
   Vnode(2,1) = node(i,2,3)  ! altitude node i, upper level

   Vnode(1,2) = node(i,1,4)  ! refractivity node i, lower level
   Vnode(2,2) = node(i,2,4)  ! refractivity node i, upper level

   Hnode(i,3) = ExpInt1D(VPt, Vnode)  ! interpolated vertical refrac.

end do

! Horizontal interpolation using BiLin2D
HPt = RefPt(1:2)  ! interpolate to this (x,y)-coordinate

! x,y coordinates
do i=1, p
   Hnode(i,1) = node(i,1,1)   ! node i, x
   Hnode(i,2) = node(i,1,2)   ! node i, y
end do
if (p .eq. 3) then
   ! Copy node 2 to node 4 as required by BiLinear2D
   Hnode(4,:) = Hnode(2,:)
end if

!BiLinearExp3D = BiLinear2D(HPt, Hnode)
!BiLinearExp3D = BiLinearDiff2D(HPt, Hnode, p)
BiLinearExp3D = Shepard2D(HPt, Hnode(1:p,:))

end function BiLinearExp3D
!****
! <<<<<<<< End of RoboDoc comments


subroutine simpned ( x, y, num, result )
!
!***********************************************************************
! Routine taken from the "intlib.f90" Library provided by John Burkardt
! http://orion.math.iastate.edu/burkardt/f_src/intlib/intlib.html
! Changed to double precision and fortran 90
!***********************************************************************
!
!! SIMPNE approximates the integral of unevenly spaced data.
!
!
!  Discussion:
!
!    The routine repeatedly interpolates a 3-point Lagrangian polynomial
!    to the data and integrates that exactly.
!
!  Reference:
!
!    Philip Davis and Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Blaisdell Publishing, 1967.
!
!  Modified:
!
!    30 October 2000
!
!  Parameters:
!
!    Input, real X(NUM), contains the X values of the data, in order.
!
!    Input, real Y(NUM), contains the Y values of the data.
!
!    Input, integer NUM, number of data points.  NUM must be at least 3.
!
!    Output, real RESULT.
!    RESULT is the approximate value of the integral.
!
  integer,   intent(in)  :: num
  real (wp), intent(in)  :: x(num)
  real (wp), intent(in)  :: y(num)
  real (wp), intent(out) :: result
!
  real (wp) :: del(3)
  real (wp) :: e
  real (wp) :: f
  real (wp) :: feints
  real (wp) :: g(3)
  integer       :: i
  integer       :: n
  real (wp) :: pi(3)
  real (wp) :: sum1
  real (wp) :: x1
  real (wp) :: x2
  real (wp) :: x3
!
  result = 0.0_wp

  if ( num <= 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SIMPNE - Fatal error!'
    write ( *, '(a)' ) '  NUM <= 2.'
    stop
  end if

  n = 1

  do

    x1 = x(n)
    x2 = x(n+1)
    x3 = x(n+2)
    e = x3*x3-x1*x1
    f = x3*x3*x3-x1*x1*x1
    feints = x3-x1
    del(1) = x3-x2
    del(2) = x1-x3
    del(3) = x2-x1
    g(1) = x2+x3
    g(2) = x1+x3
    g(3) = x1+x2
    pi(1) = x2*x3
    pi(2) = x1*x3
    pi(3) = x1*x2

    sum1 = 0.0_wp
    do i = 1, 3
      sum1 = sum1 + y(n-1+i)*del(i)*(f/3.0E+00-g(i)*0.50_wp+00*e+pi(i)*feints)
    end do
    result = result - sum1 / ( del(1) * del(2) * del(3) )

    n = n+2

    if ( n + 1 >= num ) then
      exit
    end if

  end do

  if ( mod(num,2) /= 0 ) then
    return
  end if

  n = num-2
  x3 = x(num)
  x2 = x(num-1)
  x1 = x(num-2)
  e = x3*x3-x2*x2
  f = x3*x3*x3-x2*x2*x2
  feints = x3-x2
  del(1) = x3-x2
  del(2) = x1-x3
  del(3) = x2-x1
  g(1) = x2+x3
  g(2) = x1+x3
  g(3) = x1+x2
  pi(1) = x2*x3
  pi(2) = x1*x3
  pi(3) = x1*x2

  sum1 = 0.0_wp
  do i = 1, 3
    sum1 = sum1 + y(n-1+i) * del(i) * &
      ( f / 3.0_wp+00 - g(i) * 0.50_wp+00 * e + pi(i) * feints )
  end do

  result = result - sum1 / ( del(1) * del(2) * del(3) )

  return
end subroutine simpned


!---------------------------------------------------------------------
! subroutine UpdateK123
!---------------------------------------------------------------------
!
! Start of RoboDoc comments >>>>
!****S* Coord/UpdateK123
!
! Name
! UpdateK123
!
! Call
! call  UpdateK123 (k1, k2, k3, err)
!
! Purpose
! Change the empirical constants k1, k2, k3 used to estimate
! the atmospheric refractivity with the Smith & Weintraub formula
! or the Thayer formula.
! The constants are checked, for invalid data the default set given
! by Bevis 1994 is used. If one of the constants is invalid none of
! them is changed to preserve a consistent set of constants.
! For k1 < 0 nothing is done but printing the default constants.
! (k1 < 0, e. g. if k1, k2, k3 are not given in the namelist.)
!
! Input
! nk1  - empirical constant k1 in K/hPa
! nk2  - empirical constant k2 in K/hPa
! nk3  - empirical constant k1 in K**2/hPa
!
! Output
! err - error code :  err = 0 - no errors
!                     err > 0 - invalid k1, k2, or k3
!
! External References
! None
!
! Record of revisions
! Date        Programmer   Description of change
! ====        ==========   =====================
! 09.08.2012  M. Bender    neu
!---------------------------------------------------------------------
! Source
!---------------------------------------------------------------------
subroutine UpdateK123 (nk1, nk2, nk3, err)

! List of calling arguments:
real (wp), intent(in) :: nk1, nk2, nk3
integer, intent(out)             :: err

! List of local variables:
real (wp), parameter :: mink1 = 77.0_wp, maxk1 = 78.0_wp
real (wp), parameter :: mink2 = 63.0_wp, maxk2 = 75.0_wp
real (wp), parameter :: mink3 = 3.70E5_wp, maxk3 = 3.80E5_wp
logical                         :: changek
!---------------------------------------------------------------------

err = 0
changek = .false.

if (nk1 > 0 .and. .not. (nk1 == k1 .and. nk2 == k2 .and. nk3 == k3)) then
   ! Check new set of constants

   changek = .true.
   if (nk1 .lt. mink1 .or. nk1 .gt. maxk1) then
      write(*,*) 'UpdateK123> Invalid k1: ', nk1
      write(*,*) 'UpdateK123> k1 must be within ', mink1,  &
                 ' and ', maxk1, ' K/hPa'
      changek = .false.
   end if
   if (nk2 .lt. mink2 .or. nk2 .gt. maxk2) then
      write(*,*) 'UpdateK123> Invalid k2: ', nk2
      write(*,*) 'UpdateK123> k2 must be within ', mink2,  &
                 ' and ', maxk2, ' K/hPa'
      changek = .false.
   end if
   if (nk3 .lt. mink3 .or. nk3 .gt. maxk3) then
      write(*,*) 'UpdateK123> Invalid k3: ', nk3
      write(*,*) 'UpdateK123> k3 must be within ', mink3,  &
                 ' and ', maxk3, ' K**2/hPa'
      changek = .false.
   end if
   if (changek) then
      ! Change k1, k2, k3
      k1 = nk1
      k2 = nk2
      k3 = nk3
!     k22 =  k2 * EMRDRD        ! = k2 - (Mw/Md)*k2
      write(*,*) 'UpdateK123> k1, k2, k3 was changed:'
      write(*,*) 'UpdateK123> k1 = ', k1
      write(*,*) 'UpdateK123> k2 = ', k2
      write(*,*) 'UpdateK123> k3 = ', k3
   else
      ! At least one of the constants was wrong, don't change
      ! anything to preserve a consistent set of constants
      write(*,*) 'UpdateK123> Error changing k1, k2, k3'
      err = 1
   end if

end if

if (.not. changek) then
   ! k1, k2, k3 will not be changed, use default
   write(*,*) 'UpdateK123> Use default empirical constants k1, k2, k3', &
                         ' given by Bevis, 1994:'
   write(*,*) 'UpdateK123> Default k1 = ', k1
   write(*,*) 'UpdateK123> Default k2 = ', k2
   write(*,*) 'UpdateK123> Default k3 = ', k3
end if

end subroutine UpdateK123
!****
! <<<<<<<< End of RoboDoc comments


!---------------------------------------------------------------------
! subroutine GMF
!---------------------------------------------------------------------
!
! Start of RoboDoc comments >>>>
!****S* Mapping/GMF
!
! Name
! gmf
!
! Call
! call gmf (dmjd,dlat,dlon,dhgt,zd,gmfh,gmfw)
!
! Purpose
! Global mapping function
!
! Reference:
! Boehm, J., A.E. Niell, P. Tregoning, H. Schuh (2006),
! Global Mapping Functions (GMF): A new empirical mapping function
! based on numerical weather model data,
! Geoph. Res. Letters, Vol. 33, L07304, doi:10.1029/2005GL025545.
!
! Download: http://ggosatm.hg.tuwien.ac.at/DELAY/
!
! Input
! dmjd -  modified julian date
! dlat -  ellipsoidal latitude in radians
! dlon -  longitude in radians
! dhgt -  height in m
! zd   -  zenith distance in radians
!
! Output
! gmfh -  hydrostatic mapping function
! gmfw -  wet mapping function
!
! External References
! None
!
! Record of revisions
! Date        Programmer   Description of change
! ====        ==========   =====================
! 31.07.2012  M. Bender    Copy of original GMF, some minor changes
!                          to fortran 95
!---------------------------------------------------------------------
! Source
!---------------------------------------------------------------------

      subroutine gmf (dmjd, dlat, dlon, dhgt, zd, gmfh, gmfw)

!     This subroutine determines the Global Mapping Functions GMF
!
!     Reference: Boehm, J., A.E. Niell, P. Tregoning, H. Schuh (2006),
!     Global Mapping Functions (GMF): A new empirical mapping
!     function based on numerical weather model data,
!     Geoph. Res. Letters, Vol. 33, L07304, doi:10.1029/2005GL025545.
!
!     input data
!     ----------
!     dmjd: modified julian date
!     dlat: ellipsoidal latitude in radians
!     dlon: longitude in radians
!     dhgt: height in m
!     zd:   zenith distance in radians
!
!     output data
!     -----------
!     gmfh: hydrostatic mapping function
!     gmfw: wet mapping function
!
!     Johannes Boehm, 2005 August 30
!     Rev. Boehm 21 July 2011: latitude -> ellipsoidal latitude
!
      ! implicit double precision (a-h,o-z)

      ! List of calling arguments:
      real (wp), intent(in) ::  dmjd    ! modified julian date
      real (wp), intent(in) ::  dlat    ! ellipsoidal latitude in radians
      real (wp), intent(in) ::  dlon    ! longitude in radians
      real (wp), intent(in) ::  dhgt    ! height in m
      real (wp), intent(in) ::  zd      ! zenith distance in radians

      real (wp), intent(out) ::  gmfh    ! hydrostatic mapping function
      real (wp), intent(out) ::  gmfw    ! wet mapping function

      ! List of local variables:
      real (wp), dimension(1:20)      :: dfac
      real (wp), dimension(1:10,1:10) :: P
      real (wp), dimension(1:55)      :: aP, bP, ah_mean, bh_mean
      real (wp), dimension(1:55)      :: ah_amp, bh_amp, aw_mean
      real (wp), dimension(1:55)      :: bw_mean, aw_amp, bw_amp
      !dimension dfac(20),P(10,10),aP(55),bP(55),                &
      !          ah_mean(55),bh_mean(55),ah_amp(55),bh_amp(55),  &
      !          aw_mean(55),bw_mean(55),aw_amp(55),bw_amp(55)

      integer       :: i, j, k, m, n, ir
      real (wp) :: doy, t, sum, bh, c0h, phh, c11h, c10h
      real (wp) :: ch, ahm, aha, ah, sine, beta, gamma
      real (wp) :: topcon, a_ht, b_ht, c_ht, hs_km
      real (wp) :: ht_corr_coef, ht_corr, bw, cw, awm
      real (wp) :: awa, aw
!---------------------------------------------------------------------

     ! pi = 3.14159265359d0

      data (ah_mean(i),i=1,55)/                                     &
      +1.2517d+02, +8.503d-01, +6.936d-02, -6.760d+00, +1.771d-01,  &
      +1.130d-02, +5.963d-01, +1.808d-02, +2.801d-03, -1.414d-03,   &
      -1.212d+00, +9.300d-02, +3.683d-03, +1.095d-03, +4.671d-05,   &
      +3.959d-01, -3.867d-02, +5.413d-03, -5.289d-04, +3.229d-04,   &
      +2.067d-05, +3.000d-01, +2.031d-02, +5.900d-03, +4.573d-04,   &
      -7.619d-05, +2.327d-06, +3.845d-06, +1.182d-01, +1.158d-02,   &
      +5.445d-03, +6.219d-05, +4.204d-06, -2.093d-06, +1.540d-07,   &
      -4.280d-08, -4.751d-01, -3.490d-02, +1.758d-03, +4.019d-04,   &
      -2.799d-06, -1.287d-06, +5.468d-07, +7.580d-08, -6.300d-09,   &
      -1.160d-01, +8.301d-03, +8.771d-04, +9.955d-05, -1.718d-06,   &
      -2.012d-06, +1.170d-08, +1.790d-08, -1.300d-09, +1.000d-10/

      data (bh_mean(i),i=1,55)/                                    &
      +0.000d+00, +0.000d+00, +3.249d-02, +0.000d+00, +3.324d-02,  &
      +1.850d-02, +0.000d+00, -1.115d-01, +2.519d-02, +4.923d-03,  &
      +0.000d+00, +2.737d-02, +1.595d-02, -7.332d-04, +1.933d-04,  &
      +0.000d+00, -4.796d-02, +6.381d-03, -1.599d-04, -3.685d-04,  &
      +1.815d-05, +0.000d+00, +7.033d-02, +2.426d-03, -1.111d-03,  &
      -1.357d-04, -7.828d-06, +2.547d-06, +0.000d+00, +5.779d-03,  &
      +3.133d-03, -5.312d-04, -2.028d-05, +2.323d-07, -9.100d-08,  &
      -1.650d-08, +0.000d+00, +3.688d-02, -8.638d-04, -8.514d-05,  &
      -2.828d-05, +5.403d-07, +4.390d-07, +1.350d-08, +1.800d-09,  &
      +0.000d+00, -2.736d-02, -2.977d-04, +8.113d-05, +2.329d-07,  &
      +8.451d-07, +4.490d-08, -8.100d-09, -1.500d-09, +2.000d-10/

      data (ah_amp(i),i=1,55)/                                      &
      -2.738d-01, -2.837d+00, +1.298d-02, -3.588d-01, +2.413d-02,   &
      +3.427d-02, -7.624d-01, +7.272d-02, +2.160d-02, -3.385d-03,   &
      +4.424d-01, +3.722d-02, +2.195d-02, -1.503d-03, +2.426d-04,   &
      +3.013d-01, +5.762d-02, +1.019d-02, -4.476d-04, +6.790d-05,   &
      +3.227d-05, +3.123d-01, -3.535d-02, +4.840d-03, +3.025d-06,   &
      -4.363d-05, +2.854d-07, -1.286d-06, -6.725d-01, -3.730d-02,   &
      +8.964d-04, +1.399d-04, -3.990d-06, +7.431d-06, -2.796d-07,   &
      -1.601d-07, +4.068d-02, -1.352d-02, +7.282d-04, +9.594d-05,   &
      +2.070d-06, -9.620d-08, -2.742d-07, -6.370d-08, -6.300d-09,   &
      +8.625d-02, -5.971d-03, +4.705d-04, +2.335d-05, +4.226d-06,   &
      +2.475d-07, -8.850d-08, -3.600d-08, -2.900d-09, +0.000d+00/

      data (bh_amp(i),i=1,55)/                                       &
      +0.000d+00, +0.000d+00, -1.136d-01, +0.000d+00, -1.868d-01,    &
      -1.399d-02, +0.000d+00, -1.043d-01, +1.175d-02, -2.240d-03,    &
      +0.000d+00, -3.222d-02, +1.333d-02, -2.647d-03, -2.316d-05,    &
      +0.000d+00, +5.339d-02, +1.107d-02, -3.116d-03, -1.079d-04,    &
      -1.299d-05, +0.000d+00, +4.861d-03, +8.891d-03, -6.448d-04,    &
      -1.279d-05, +6.358d-06, -1.417d-07, +0.000d+00, +3.041d-02,    &
      +1.150d-03, -8.743d-04, -2.781d-05, +6.367d-07, -1.140d-08,    &
      -4.200d-08, +0.000d+00, -2.982d-02, -3.000d-03, +1.394d-05,    &
      -3.290d-05, -1.705d-07, +7.440d-08, +2.720d-08, -6.600d-09,    &
      +0.000d+00, +1.236d-02, -9.981d-04, -3.792d-05, -1.355d-05,    &
      +1.162d-06, -1.789d-07, +1.470d-08, -2.400d-09, -4.000d-10/

      data (aw_mean(i),i=1,55)/                                      &
      +5.640d+01, +1.555d+00, -1.011d+00, -3.975d+00, +3.171d-02,    &
      +1.065d-01, +6.175d-01, +1.376d-01, +4.229d-02, +3.028d-03,    &
      +1.688d+00, -1.692d-01, +5.478d-02, +2.473d-02, +6.059d-04,    &
      +2.278d+00, +6.614d-03, -3.505d-04, -6.697d-03, +8.402d-04,    &
      +7.033d-04, -3.236d+00, +2.184d-01, -4.611d-02, -1.613d-02,    &
      -1.604d-03, +5.420d-05, +7.922d-05, -2.711d-01, -4.406d-01,    &
      -3.376d-02, -2.801d-03, -4.090d-04, -2.056d-05, +6.894d-06,    &
      +2.317d-06, +1.941d+00, -2.562d-01, +1.598d-02, +5.449d-03,    &
      +3.544d-04, +1.148d-05, +7.503d-06, -5.667d-07, -3.660d-08,    &
      +8.683d-01, -5.931d-02, -1.864d-03, -1.277d-04, +2.029d-04,    &
      +1.269d-05, +1.629d-06, +9.660d-08, -1.015d-07, -5.000d-10/

      data (bw_mean(i),i=1,55)/                                      &
      +0.000d+00, +0.000d+00, +2.592d-01, +0.000d+00, +2.974d-02,    &
      -5.471d-01, +0.000d+00, -5.926d-01, -1.030d-01, -1.567d-02,    &
      +0.000d+00, +1.710d-01, +9.025d-02, +2.689d-02, +2.243d-03,    &
      +0.000d+00, +3.439d-01, +2.402d-02, +5.410d-03, +1.601d-03,    &
      +9.669d-05, +0.000d+00, +9.502d-02, -3.063d-02, -1.055d-03,    &
      -1.067d-04, -1.130d-04, +2.124d-05, +0.000d+00, -3.129d-01,    &
      +8.463d-03, +2.253d-04, +7.413d-05, -9.376d-05, -1.606d-06,    &
      +2.060d-06, +0.000d+00, +2.739d-01, +1.167d-03, -2.246d-05,    &
      -1.287d-04, -2.438d-05, -7.561d-07, +1.158d-06, +4.950d-08,    &
      +0.000d+00, -1.344d-01, +5.342d-03, +3.775d-04, -6.756d-05,    &
      -1.686d-06, -1.184d-06, +2.768d-07, +2.730d-08, +5.700d-09/

      data (aw_amp(i),i=1,55)/                                       &
      +1.023d-01, -2.695d+00, +3.417d-01, -1.405d-01, +3.175d-01,    &
      +2.116d-01, +3.536d+00, -1.505d-01, -1.660d-02, +2.967d-02,    &
      +3.819d-01, -1.695d-01, -7.444d-02, +7.409d-03, -6.262d-03,    &
      -1.836d+00, -1.759d-02, -6.256d-02, -2.371d-03, +7.947d-04,    &
      +1.501d-04, -8.603d-01, -1.360d-01, -3.629d-02, -3.706d-03,    &
      -2.976d-04, +1.857d-05, +3.021d-05, +2.248d+00, -1.178d-01,    &
      +1.255d-02, +1.134d-03, -2.161d-04, -5.817d-06, +8.836d-07,    &
      -1.769d-07, +7.313d-01, -1.188d-01, +1.145d-02, +1.011d-03,    &
      +1.083d-04, +2.570d-06, -2.140d-06, -5.710d-08, +2.000d-08,    &
      -1.632d+00, -6.948d-03, -3.893d-03, +8.592d-04, +7.577d-05,    &
      +4.539d-06, -3.852d-07, -2.213d-07, -1.370d-08, +5.800d-09/

      data (bw_amp(i),i=1,55)/                                       &
      +0.000d+00, +0.000d+00, -8.865d-02, +0.000d+00, -4.309d-01,    &
      +6.340d-02, +0.000d+00, +1.162d-01, +6.176d-02, -4.234d-03,    &
      +0.000d+00, +2.530d-01, +4.017d-02, -6.204d-03, +4.977d-03,    &
      +0.000d+00, -1.737d-01, -5.638d-03, +1.488d-04, +4.857d-04,    &
      -1.809d-04, +0.000d+00, -1.514d-01, -1.685d-02, +5.333d-03,    &
      -7.611d-05, +2.394d-05, +8.195d-06, +0.000d+00, +9.326d-02,    &
      -1.275d-02, -3.071d-04, +5.374d-05, -3.391d-05, -7.436d-06,    &
      +6.747d-07, +0.000d+00, -8.637d-02, -3.807d-03, -6.833d-04,    &
      -3.861d-05, -2.268d-05, +1.454d-06, +3.860d-07, -1.068d-07,    &
      +0.000d+00, -2.658d-02, -1.947d-03, +7.131d-04, -3.506d-05,    &
      +1.885d-07, +5.792d-07, +3.990d-08, +2.000d-08, -5.700d-09/

!     reference day is 28 January
!     this is taken from Niell (1996) to be consistent
      doy = dmjd  - 44239.0_wp + 1 - 28

!     parameter t
      t = sin(dlat)

!     degree n and order m
      n = 9
      m = 9

! determine n!  (faktorielle)  moved by 1
      dfac(1) = 1
      do i = 1,(2*n + 1)
        dfac(i+1) = dfac(i)*i
      end do

!     determine Legendre functions (Heiskanen and Moritz,
!     Physical Geodesy, 1967, eq. 1-62)
      do i = 0,n
        do j = 0,min(i,m)
          ir = int((i - j)/2)
          sum = 0
          do k = 0,ir
            sum = sum + (-1)**k*dfac(2*i - 2*k + 1)/dfac(k + 1)/      &
              dfac(i - k + 1)/dfac(i - j - 2*k + 1)*t**(i - j - 2*k)
          end do
!         Legendre functions moved by 1
          P(i + 1,j + 1) = 1.0_wp/2**i*sqrt((1 - t**2)**(j))*sum
        end do
      end do

!     spherical harmonics
      i = 0
      do n = 0,9
        do m = 0,n
          i = i + 1
          aP(i) = P(n+1,m+1)*cos(m*dlon)
          bP(i) = P(n+1,m+1)*sin(m*dlon)
        end do
      end do

!     hydrostatic
      bh = 0.0029
      c0h = 0.062
      if (dlat.lt.0) then ! southern hemisphere
        phh  = pi
        c11h = 0.007
        c10h = 0.002
      else                ! northern hemisphere
        phh  = 0
        c11h = 0.005
        c10h = 0.001
      end if
      ch = c0h + ((cos(doy/365.25_wp*2*pi + phh)+1)*c11h/2 + c10h)*  &
                 (1-cos(dlat))

      ahm = 0.0_wp
      aha = 0.0_wp
      do i = 1,55
        ahm = ahm + (ah_mean(i)*aP(i) + bh_mean(i)*bP(i))*1d-5
        aha = aha + (ah_amp(i) *aP(i) + bh_amp(i) *bP(i))*1d-5
      end do
      ah  = ahm + aha*cos(doy/365.25_wp*2.0_wp*pi)

      sine   = sin(pi/2 - zd)
      beta   = bh/( sine + ch  )
      gamma  = ah/( sine + beta)
      topcon = (1.0_wp + ah/(1.0_wp + bh/(1.0_wp + ch)))
      gmfh   = topcon/(sine+gamma)

!     height correction for hydrostatic mapping function from Niell (1996)
      a_ht = 2.53d-5
      b_ht = 5.49d-3
      c_ht = 1.14d-3
      hs_km  = dhgt/1000.0_wp
!
      beta   = b_ht/( sine + c_ht )
      gamma  = a_ht/( sine + beta)
      topcon = (1.0_wp + a_ht/(1.0_wp + b_ht/(1.0_wp + c_ht)))
      ht_corr_coef = 1/sine - topcon/(sine + gamma)
      ht_corr      = ht_corr_coef * hs_km
      gmfh         = gmfh + ht_corr

!     wet
      bw = 0.00146
      cw = 0.04391

      awm = 0.0_wp
      awa = 0.0_wp
      do i = 1,55
        awm = awm + (aw_mean(i)*aP(i) + bw_mean(i)*bP(i))*1d-5
        awa = awa + (aw_amp(i) *aP(i) + bw_amp(i) *bP(i))*1d-5
      end do
      aw =  awm + awa*cos(doy/365.25_wp*2*pi)

      beta   = bw/( sine + cw )
      gamma  = aw/( sine + beta)
      topcon = (1.0_wp + aw/(1.0_wp + bw/(1.0_wp + cw)))
      gmfw   = topcon/(sine+gamma)

end subroutine gmf
!****
! <<<<<<<< End of RoboDoc comments


!---------------------------------------------------------------------
! function JulianDate
!---------------------------------------------------------------------
!
! Start of RoboDoc comments >>>>
!****S* GPSNav/JulianDate
!
! Name
! JulianDate
!
! Purpose
! Compute the Julian Date from a given Gregorian Date
! The time reference is UT (universal Time)
!
! WARNING: This routine is correct only for 1801 <= Year <= 2099
!
! Formulas and Fortran examples can be found in
! http://aa.usno.navy.mil/faq/docs/JD_Formula.html
!
! Procedure
!
!      INTEGER FUNCTION JD (YEAR,MONTH,DAY)
!C
!C---COMPUTES THE JULIAN DATE (JD) GIVEN A GREGORIAN CALENDAR
!C   DATE (YEAR,MONTH,DAY).
!C
!      INTEGER YEAR,MONTH,DAY,I,J,K
!C
!      I= YEAR
!      J= MONTH
!      K= DAY
!C
!      JD= K-32075+1461*(I+4800+(J-14)/12)/4+367*(J-2-(J-14)/12*12)
!     2    /12-3*((I+4900+(J-14)/12)/100)/4
!C
!      RETURN
!      END
!
! Record of revisions
! Date        Programmer   Description of change
! ====        ==========   =====================
! 25. 1.2005  M. Bender    new
! 29.03.2007  M. Bender    Code changed to the above formula
! 11.01.2010  M. Bender    Warning: output of year
!---------------------------------------------------------------------
! Source
!---------------------------------------------------------------------
function JulianDate(Year, Month, Day, Hour, Min, Sec) !result( JulianDate )

! List of calling arguments:
real (wp)                :: JulianDate
integer, intent(in)      :: Year, Month, Day, Hour, Min
real (wp), intent(in)    :: Sec

! List of local variables:
integer (i4)             :: JD

!integer :: Y, M
!---------------------------------------------------------------------

!!$if (Month .le. 2) then
!!$   Y = Year  -  1
!!$   M = Month + 12
!!$else
!!$   Y = Year
!!$   M = Month
!!$end if

!JulianDate = int(365.25*Y) + int(30.6001*(M+1)) + Day + Hour/24.0 + 1720981.5

if (year .lt. 1801 .or. year .gt. 2099) then
   write(*,*) 'JulianDate> Warning - year out of range, date might be wrong', &
              ' year = ', year
end if

JD = Day - 32075 + 1461*(Year+4800+(Month-14)/12)/4    &
                 + 367*(Month-2-(Month-14)/12*12)/12   &
                 - 3*((Year+4900+(Month-14)/12)/100)/4

! Hour, Min, Sec are fractions of one day
! 1 h = 1/24 = 0.041666
! 1 min = 1/(24*60) = 1/1440 = 6.9444e-4
! 1 sec = 1/(24*60*60) = 1/86400 = 1.1574e-5
JulianDate = real(JD) - 0.5_wp  + real(Hour)*(1.0_wp/24.0_wp)     &
                               + real(Min)*(1.0_wp/1440.0_wp)    &
                               + Sec*(1.0_wp/86400.0_wp)

! write(*,*) 'JulianDate >> ', JD, JulianDate, Hour

end function JulianDate
!****
! <<<<<<<< End of RoboDoc comments


subroutine GregorianDate( JulianDate, Year, Month, Day, Hour, Min, Sec)

! List of calling arguments:
real (wp), intent(in)     :: JulianDate
integer, intent(out) :: Year, Month, Day
integer, intent(out) :: Hour, Min
real (wp), intent(out)    :: Sec

! List of local variables:
integer (i4) :: JD, L, N, I, J, K
integer      :: f, g
real (wp)    :: s
!---------------------------------------------------------------------

!!$b = int(JulianDate + 0.5) + 1537
!!$c = int((b-122.1)/365.25)
!!$d = int(365.25*c)
!!$e = int((b-d)/30.6001)

! second of day (1 day = 24 h = 24*60*60 = 86400 s
s  = (JulianDate + 0.5_wp - int(JulianDate + 0.5_wp)) * 86400.0_wp
! write(*,*) s
Hour = int(s/3600.0)
f = Hour*3600
Min  = int((s-real(f))/60.0)
g = Min*60
Sec  = s-real(f)-real(g)
!write(*,*) 'Hour, Min, Sec = ', Hour, Min, Sec

!!$Day   = b - d - int(30.6001*e)
!!$Month = e - 1 - 12*int(3/14)
!!$if (Month .gt. 12) then
!!$   Month = Month - 12
!!$end if
!!$Year  = c - 4715 - int((7+Month)/10)

! The original code is referenced to 12:00 h am
JD = int(JulianDate+0.5_wp)

L  = JD+68569
N  = 4*L/146097
L  = L-(146097*N+3)/4
I  = 4000*(L+1)/1461001
L  = L-1461*I/4+31
J  = 80*L/2447
K  = L-2447*J/80
L  = J/11
J  = J+2-12*L
I  = 100*(N-49)+I+L

Year  = I
Month = J
Day   = K

if (year .lt. 1801 .or. year .gt. 2099) then
   write(*,*) 'GregorianDate> Warning - year out of range, ', &
              ' date might be wrong, year = ', year
end if

end subroutine GregorianDate


!---------------------------------------------------------------------
! subroutine JD2MJD
!---------------------------------------------------------------------
!
! Start of RoboDoc comments >>>>
!****S* GPSNav/JD2MJD
!
! Name
! JD2MJD
!
! Purpose
! Convert a Julian Date (JD) into a Modified Julian Date (MJD)
! ModifiedJulianDate = JulianDate - 2400000.5
!
! Input
! JulianDate - Julian Date
!
! Output
! JD2MJD - Modified Julian Date (return value)
!
! External References
! None
!
! Record of revisions
! Date        Programmer   Description of change
! ====        ==========   =====================
! 23. 5.2007  M. Bender    neu
!---------------------------------------------------------------------
! Source
!---------------------------------------------------------------------
function JD2MJD( JulianDate )  !result( JD2MJD )

! List of calling arguments:
real (wp) :: JD2MJD
real (wp) :: JulianDate
!---------------------------------------------------------------------

JD2MJD = JulianDate - 2400000.5_wp

end function JD2MJD
!****
! <<<<<<<< End of RoboDoc comments


!---------------------------------------------------------------------
! function MJD2JD
!---------------------------------------------------------------------
!
! Start of RoboDoc comments >>>>
!****F* GPSNav/MJD2JD
!
! Name
! MJD2JD
!
! Purpose
! Convert a Modified Julian Date (MJD) into a Julian Date (JD)
! JulianDate = ModifiedJulianDate + 2400000.5
!
! Input
!  ModifiedJulianDate -  Modified Julian Date
!
! Output
! MJD2JD - Julian Date (return value)
!
! External References
! None
!
! Record of revisions
! Date        Programmer   Description of change
! ====        ==========   =====================
! 23. 5.2007  M. Bender    neu
!---------------------------------------------------------------------
! Source
!---------------------------------------------------------------------
function MJD2JD(  ModifiedJulianDate )

! List of calling arguments:
real (wp) ::  ModifiedJulianDate
real (wp) :: MJD2JD
!---------------------------------------------------------------------

MJD2JD =  ModifiedJulianDate + 2400000.5_wp

end function MJD2JD
!****
! <<<<<<<< End of RoboDoc comments



! ****************************************************************************************
! The following section contains code obtained from  DAVID G. SIMPSON's                  *
! web page (NASA):                                                                       *
! http://www.davidgsimpson.com/software.html                                             *
!                                                                                        *
! Some date conversion routines are used here:                                           *
! GREG2DOY   -  converts a date on the Gregorian or Julian calendars to a day of year.   *
! DOY2GREG   -  converts a day of year to a date on the Gregorian or Julian calendars.   *
! LEAP       -  determines whether a given year on the Gregorian calendar is a leap year *
!                                                                                        *
! The original code and some other calendar routines are saved in the archive            *
! CodeSimpson.tgz                                                                        *
! ****************************************************************************************



!******************************************************************************
!
!                                                         G R E G 2 D O Y
!
!
!  Program:      GREG2DOY
!
!  Programmer:   David G. Simpson
!                NASA Goddard Space Flight Center
!                Greenbelt, Maryland  20771
!
!  Date:         November 20, 2001
!
!  Language:     ANSI Standard Fortran-90
!
!  Version:      1.00b  (October 25, 2004)
!
!  Description:  This program converts a date on the Gregorian or Julian
!                calendars to a day of year.
!
!  Note:         Array GREGORIAN_START defines the end dates of the Julian
!                calendar and start dates of the Gregorian calendar.
!                Set the parameter GREGORIAN_CHOICE to indicate the desired
!                start date of the Gregorian calendar, as listed in
!                array GREGORIAN_START.
!
!******************************************************************************

!******************************************************************************
!  Main program
!******************************************************************************

!      PROGRAM GREG2DOY
      subroutine GREG2DOY (Y, M, D, DOY)

!!$      TYPE :: DATE_TYPE
!!$         INTEGER :: YEAR_J                                   ! year of end of Julian calendar
!!$         INTEGER :: MONTH_J                                  ! month of end of Julian calendar
!!$         INTEGER :: DAY_J                                    ! day of end of Julian calendar
!!$         INTEGER :: YEAR_G                                   ! year of start of Gregorian calendar
!!$         INTEGER :: MONTH_G                                  ! month of start of Gregorian calendar
!!$         INTEGER :: DAY_G                                    ! day of start of Gregorian calendar
!!$         INTEGER :: NDAYS                                    ! number of days dropped from calendar at switch
!!$      END TYPE DATE_TYPE

      INTEGER :: D                                              ! day of month (+ fraction)
      INTEGER :: DOY
      INTEGER :: K
      INTEGER :: M                                              ! month (1-12)
      INTEGER :: Y                                              ! year
      LOGICAL :: GREGORIAN_FLAG                                 ! .TRUE. for Gregorian date, .FALSE. for Julian
      LOGICAL :: LEAP

      TYPE (DATE_TYPE), DIMENSION (3) :: GREGORIAN_START =   &
         (/ DATE_TYPE (1582, 10,  4, 1582, 10, 15, 10),      &  ! 1: Decree by Pope Gregory XIII
            DATE_TYPE (1752,  9,  2, 1752,  9, 14, 11),      &  ! 2: Great Britain
            DATE_TYPE (1918,  1, 31, 1918,  2, 14, 13)  /)      ! 3: Russia

      INTEGER, PARAMETER :: GREGORIAN_CHOICE = 1                ! set to 1 for 1582 date, 2 for 1752 date, etc.

      !LOGICAL :: GREGORIAN



!------------------------------------------------------------------------------
!  Main program code
!------------------------------------------------------------------------------

!!$      WRITE (UNIT=*, FMT='(A)', ADVANCE='NO') ' Enter month (1-12):  '  ! prompt for month
!!$      READ (UNIT=*, FMT=*) M
!!$
!!$      WRITE (UNIT=*, FMT='(A)', ADVANCE='NO') ' Enter day:  '           ! prompt for day of month
!!$      READ (UNIT=*, FMT=*) D
!!$
!!$      WRITE (UNIT=*, FMT='(A)', ADVANCE='NO') ' Enter year:  '          ! prompt for year
!!$      READ (UNIT=*, FMT=*) Y

      ! test for Gregorian calendar
      GREGORIAN_FLAG = GREGORIAN( Y, M, INT(D),                       &
                                  GREGORIAN_START(GREGORIAN_CHOICE) )

      LEAP = .FALSE.
      IF (MOD(Y,4) .EQ. 0) LEAP = .TRUE.

      IF (GREGORIAN_FLAG) THEN
         IF (MOD(Y,100) .EQ. 0) LEAP = .FALSE.
         IF (MOD(Y,400) .EQ. 0) LEAP = .TRUE.
      END IF

      IF (LEAP) THEN
         K = 1
      ELSE
         K = 2
      END IF

      DOY = ((275*M)/9) - K*((M+9)/12) + D - 30

      IF (GREGORIAN_FLAG .AND. (Y .EQ. GREGORIAN_START(GREGORIAN_CHOICE)%YEAR_G)) THEN
         DOY = DOY - GREGORIAN_START(GREGORIAN_CHOICE)%NDAYS
      END IF

      IF (.NOT. GREGORIAN_FLAG) THEN                                                ! print msg if Julian calendar in effect
         WRITE (UNIT=*, FMT='(/,A)') ' Julian calendar.'
      END IF

!      WRITE (UNIT=*, FMT='(/,A, I3)') ' Day of year = ', DOY                        ! print result

!      END PROGRAM GREG2DOY


      contains



!******************************************************************************
!  GREGORIAN
!
!  This function determines whether a given date is in the Gregorian calendar
!  (return value of .TRUE.) or on the Julian calendar (return .FALSE.).
!******************************************************************************

      FUNCTION GREGORIAN (YEAR, MONTH, DAY, GREG_START) RESULT (GREG_FLAG)

!!$      TYPE :: DATE_TYPE
!!$         INTEGER :: YEAR_J                            ! year of end of Julian calendar
!!$         INTEGER :: MONTH_J                           ! month of end of Julian calendar
!!$         INTEGER :: DAY_J                             ! day of end of Julian calendar
!!$         INTEGER :: YEAR_G                            ! year of start of Gregorian calendar
!!$         INTEGER :: MONTH_G                           ! month of start of Gregorian calendar
!!$         INTEGER :: DAY_G                             ! day of start of Gregorian calendar
!!$         INTEGER :: NDAYS                             ! number of days dropped from calendar at switch
!!$      END TYPE DATE_TYPE

      INTEGER, INTENT(IN) :: YEAR                        ! input year
      INTEGER, INTENT(IN) :: MONTH                       ! input month
      INTEGER, INTENT(IN) :: DAY                         ! input day of month
      TYPE (DATE_TYPE), INTENT(IN) :: GREG_START         ! contains Julian stop/Gregorian start dates

      LOGICAL :: GREG_FLAG                               ! result flag (.TRUE. for Gregorian)

      INTEGER :: CALTYPE = 0                             ! 0=unknown, 1=Julian, 2=Gregorian


      IF (YEAR .LT. GREG_START%YEAR_J) THEN              ! if year before end of Julian calendar..
         CALTYPE = 1                                     ! ..then this is a Julian date
      ELSE IF (YEAR .EQ. GREG_START%YEAR_J) THEN         ! if this is the last year of the Julian cal..
         IF (MONTH .LT. GREG_START%MONTH_J) THEN         ! ..then if this is before the ending month..
            CALTYPE = 1                                  ! ..then this is a Julian date
         ELSE IF (MONTH .EQ. GREG_START%MONTH_J) THEN    ! if this is the ending month..
            IF (DAY .LE. GREG_START%DAY_J) THEN          ! ..then if this is before/at the ending date..
               CALTYPE = 1                               ! ..then this is a Julian date
            END IF
         END IF
      END IF

      IF (YEAR .GT. GREG_START%YEAR_G) THEN              ! if year after start of Gregorian calendar..
         CALTYPE = 2                                     ! ..then this is a Gregorian date
      ELSE IF (YEAR .EQ. GREG_START%YEAR_G) THEN         ! if this is the first year of the Greg. cal..
         IF (MONTH .GT. GREG_START%MONTH_G) THEN         ! ..then if this is after the starting month..
            CALTYPE = 2                                  ! ..then this is a Gregorian date
         ELSE IF (MONTH .EQ. GREG_START%MONTH_G) THEN    ! if this is the starting month..
            IF (DAY .GE. GREG_START%DAY_G) THEN          ! ..then if this is at/after the starting date..
               CALTYPE = 2                               ! ..then this is a Gregorian date
            END IF
         END IF
      END IF

      SELECT CASE (CALTYPE)                              ! check calendar type
         CASE (0)                                        ! if unknown, we have an invalid date
            WRITE (UNIT=*, FMT='(A)') ' No such date.'   ! print error message
            STOP                                         ! stop program
         CASE (1)                                        ! if Julian date..
            GREG_FLAG = .FALSE.                          ! ..set return value to .false.
         CASE (2)                                        ! if Gregorian date..
            GREG_FLAG = .TRUE.                           ! ..set return value to .true.
      END SELECT

      END FUNCTION GREGORIAN


end subroutine GREG2DOY

!******************************************************************************
!
!                                                         D O Y 2 G R E G
!
!
!  Program:      DOY2GREG
!
!  Programmer:   David G. Simpson
!                NASA Goddard Space Flight Center
!                Greenbelt, Maryland  20771
!
!  Date:         November 20, 2001
!
!  Language:     ANSI Standard Fortran-90
!
!  Version:      1.00b (October 25, 2004)
!
!  Description:  This program converts a date on the Gregorian or Julian calendars to a day of year.
!
!  Note:         Array GREGORIAN_START defines the end dates of the Julian calendar and start dates of the Gregorian calendar.
!                Set the parameter GREGORIAN_CHOICE to indicate the desired start date of the Gregorian calendar, as listed in
!                array GREGORIAN_START.
!
!******************************************************************************

!******************************************************************************
!  Main program
!******************************************************************************

      !PROGRAM DOY2GREG
      subroutine DOY2GREG(DOY, Y, M, D)

      TYPE :: DATE_TYPE                                                       ! DATE_TYPE definition
         INTEGER :: YEAR_J                                                    !  year of end of Julian calendar
         INTEGER :: DOY_J                                                     !  day of year of end of Julian calendar
         INTEGER :: NDAYS                                                     !  num of days dropped from calendar at switch
         INTEGER :: TTLDAYS                                                   !  number of days in year of switch
      END TYPE DATE_TYPE

      INTEGER :: D                                                            ! day of month (+ fraction)
      INTEGER :: DOY
      INTEGER :: K
      INTEGER :: M                                                            ! month (1-12)
      INTEGER :: Y                                                            ! year
      LOGICAL :: GREGORIAN_FLAG                                               ! .TRUE. for Gregorian date, .FALSE. for Julian
      LOGICAL :: LEAP

      CHARACTER(LEN=9), DIMENSION(12), PARAMETER :: MONTH_NAME =            & ! month names
         (/ 'January  ', 'February ', 'March    ', 'April    ', 'May      ',&
            'June     ', 'July     ', 'August   ', 'September', 'October  ',&
            'November ', 'December ' /)

      TYPE (DATE_TYPE), DIMENSION (3) :: GREGORIAN_START =        &
         (/ DATE_TYPE (1582, 277, 10, 355),                       &           ! 1: Decree by Pope Gregory XIII
            DATE_TYPE (1752, 246, 11, 355),                       &           ! 2: Great Britain
            DATE_TYPE (1918,  31, 13, 352)  /)                                ! 3: Russia

      INTEGER, PARAMETER :: GREGORIAN_CHOICE = 1                              ! set to 1 for 1582 date, 2 for 1752 date, etc.

      !LOGICAL :: GREGORIAN


!------------------------------------------------------------------------------
!  Main program code
!------------------------------------------------------------------------------

!!$      WRITE (UNIT=*, FMT='(A)', ADVANCE='NO') ' Enter day of year:  '    ! prompt for day of year
!!$      READ (UNIT=*, FMT=*) DOY
!!$
!!$      WRITE (UNIT=*, FMT='(A)', ADVANCE='NO') ' Enter year:  '           ! prompt for year
!!$      READ (UNIT=*, FMT=*) Y

      GREGORIAN_FLAG = GREGORIAN(Y, DOY, GREGORIAN_START(GREGORIAN_CHOICE)) ! test for Gregorian calendar

      LEAP = .FALSE.                                                        ! test for leap year
      IF (MOD(Y,4) .EQ. 0) LEAP = .TRUE.

      IF (GREGORIAN_FLAG) THEN                                              ! additional Gregorian leap year tests
         IF (MOD(Y,100) .EQ. 0) LEAP = .FALSE.
         IF (MOD(Y,400) .EQ. 0) LEAP = .TRUE.
      END IF

      IF (LEAP) THEN                                                        ! set K based on calendar type
         K = 1
      ELSE
         K = 2
      END IF

      IF ((Y .EQ. GREGORIAN_START(GREGORIAN_CHOICE)%YEAR_J) .AND.   &       ! if this is the year we switch calendars..
          (DOY .GT. GREGORIAN_START(GREGORIAN_CHOICE)%DOY_J)) THEN          ! ..and we're on the Gregorian calendar..
         DOY = DOY + GREGORIAN_START(GREGORIAN_CHOICE)%NDAYS                ! ..then adjust for dropped days
      END IF

      M = INT(9.0_wp*(K+DOY)/275.0_wp + 0.98_wp)                               ! compute month
      IF (DOY .LT. 32) M = 1

      D = DOY - ((275*M)/9) + K*((M+9)/12) + 30                             ! compute day of month

!!$      IF (.NOT. GREGORIAN_FLAG) THEN                                     ! print msg if Julian calendar in effect
!!$         WRITE (UNIT=*, FMT='(/,A)') ' Julian calendar.'
!!$      END IF
!!$
!!$      IF (Y .GE. 1) THEN                                                 ! print results (AD)
!!$         WRITE (UNIT=*, FMT='(/,1X,A,1X,I2,", ",I7, " AD")')   &
!!$                TRIM(MONTH_NAME(M)), D, Y
!!$      ELSE                                                               ! print results (BC)
!!$         WRITE (UNIT=*, FMT='(/,1X,A,1X,I2,", ",I7, " BC")')   &
!!$                TRIM(MONTH_NAME(M)), D, -Y+1
!!$      END IF

      !END PROGRAM DOY2GREG


      contains



!******************************************************************************
!  GREGORIAN
!
!  This function determines whether a given date is in the Gregorian calendar (return value of .TRUE.) or on the Julian calendar
!  (return value of .FALSE.).
!******************************************************************************

      FUNCTION GREGORIAN (YEAR, DOY, GREG_START) RESULT (GREG_FLAG)

!!$      TYPE :: DATE_TYPE                                    ! DATE_TYPE definition
!!$         INTEGER :: YEAR_J                                 !  year of end of Julian calendar
!!$         INTEGER :: DOY_J                                  !  day of year of end of Julian calendar
!!$         INTEGER :: NDAYS                                  !  num of days dropped from calendar at switch
!!$         INTEGER :: TTLDAYS                                !  number of days in year of switch
!!$      END TYPE DATE_TYPE

      INTEGER, INTENT(IN) :: YEAR                             ! input year
      INTEGER, INTENT(IN) :: DOY                              ! input day of month
      TYPE (DATE_TYPE), INTENT(IN) :: GREG_START              ! contains Julian stop/Gregorian start dates

      LOGICAL :: GREG_FLAG                                    ! result flag (.TRUE. for Gregorian)

      INTEGER :: CALTYPE = 0                                  ! 0=unknown, 1=Julian, 2=Gregorian
      INTEGER :: TOTAL_DAYS
      LOGICAL :: LEAP_FLAG


      LEAP_FLAG = .FALSE.
      IF (MOD(YEAR,4) .EQ. 0) LEAP_FLAG = .TRUE.

      IF (YEAR .LT. GREG_START%YEAR_J) THEN                   ! if year before end of Julian calendar..
         CALTYPE = 1                                          ! ..then this is a Julian date
         IF (LEAP_FLAG) THEN
            TOTAL_DAYS = 366
         ELSE
            TOTAL_DAYS = 365
         END IF
      ELSE IF (YEAR .EQ. GREG_START%YEAR_J) THEN              ! if this is the last year of the Julian cal..
         IF (DOY .LE. GREG_START%DOY_J) THEN                  ! ..then if this is before the ending month..
            CALTYPE = 1                                       ! ..then this is a Julian date
         END IF
         TOTAL_DAYS = GREG_START%TTLDAYS
      END IF

      IF (YEAR .GT. GREG_START%YEAR_J) THEN                   ! if year after start of Gregorian calendar..
         CALTYPE = 2                                          ! ..then this is a Gregorian date
         IF (MOD(YEAR,100) .EQ. 0) LEAP_FLAG = .FALSE.
         IF (MOD(YEAR,400) .EQ. 0) LEAP_FLAG = .TRUE.
         IF (LEAP_FLAG) THEN
            TOTAL_DAYS = 366
         ELSE
            TOTAL_DAYS = 365
         END IF
      ELSE IF (YEAR .EQ. GREG_START%YEAR_J) THEN              ! if this is the first year of the Greg. cal..
         IF (DOY .GT. GREG_START%DOY_J) THEN                  ! ..then if this is after the starting month..
            CALTYPE = 2                                       ! ..then this is a Gregorian date
         END IF
         TOTAL_DAYS = GREG_START%TTLDAYS
      END IF

      IF (DOY .GT. TOTAL_DAYS) CALTYPE = 0

      SELECT CASE (CALTYPE)                                   ! check calendar type
         CASE (0)                                             ! if unknown, we have an invalid date
            WRITE (UNIT=*, FMT='(A)') ' No such day of year.' ! print error message
            STOP                                              ! stop program
         CASE (1)                                             ! if Julian date..
            GREG_FLAG = .FALSE.                               ! ..set return value to .false.
         CASE (2)                                             ! if Gregorian date..
            GREG_FLAG = .TRUE.                                ! ..set return value to .true.
      END SELECT

      END FUNCTION GREGORIAN


end subroutine DOY2GREG


!---------------------------------------------------------------------
! subroutine ImportDomes
!---------------------------------------------------------------------
!
! Start of RoboDoc comments >>>>
!****S* GPSImport/ImportDomes
!
! Name
! ImportDomes
!
! Call
! call ImportDomes
! call ImportDomes
!
! Purpose
! Reads the file "smark_domes_SNX_SNX" which provides the ID,
! short name, full name, country and period (start date and end date).
! /dsk/rtt17/m1/trop_base/ana/global/dat/smark_domes_SNX
!
! File is sorted by station ID (first column).
!
! There is no file header in "smark_domes_SNX_SNX", only the SINEX
! sections
! +SITE/DOMES  ...  -SITE/DOMES
! +SITE_ERR/DOMES ... -SITE_ERR/DOMES
!
! Two rows from the file:
! 9830 0400 A---- P      M    Tuebingen,  Germany    06300:00000 00000:00000
! 9831 0401 A---- P      M    Freiburg,  Germany     06293:00000 00000:00000
!
! Col.  1 : "9830"
! Station ID, the same ID is used in the HORI file
!
! Col.  2 : "0400"
! Station short name, 4 characters, the same short name is given
! in the HORI file.
!
! Col.  3 : "A----"
! PTSOLN ???
!
! Col.  4 : "P"
! T ???
!
! Col.  5 : "M"
! DOMES number, unique number of GNSS station, see below
!
! Col.  6 : "Tuebingen,  Germany"
! Full station name and country
!
! Col.  7 : "06300:00000"
! Start date - YYDDD:SSSSS, year (last two digits), doy and second of day
!
! Col.  8 : "00000:00000"
! End date - YYDDD:SSSSS, year (last two digits), doy and second of day
!
! http://itrf.ensg.ign.fr/domes_desc.php - The DOMES Numbering System
! Looking at the 10002M006 DOMES number:
!
! the first 3 digits indicate the area, usually the country (100=France)
! the next 2 digits designate the site number within the country (02=Grasse)
! the next letter indicates the tracking point ("M" is a monument point
!     such as pillar, pole, brass mark... ;  "S" indicates the reference
!     point of an instrument such as the intersection of axes of an SLR
!     telescope, reference point of a VLBI antenna, DORIS or GPS Antenna
!     Reference Point)
! the last 3 digits represent a sequential point number
!     (here 006=GPS Pillar/brass mark)
!
! Input
! File         - data file name, optionally including path
! StartMJD     - reference date: only data valid at this date are read
! EndMJD       - reference period: only data valid within the period
!                between StartMJD and EndMJD are read. If several entries
!                are valid, the last one is used.
!     EndMJD is optional, StartMJD is not!
!
! Output
! GPSdomes     - array of type "StrucDomes" containing the data
!                read from file "smark_domes_SNX_SNX"
!
! External References
! None
!
! Record of revisions
! Date        Programmer   Description of change
! ====        ==========   =====================
! 01.08.2012  M. Bender    new
!---------------------------------------------------------------------
! Source
!---------------------------------------------------------------------
subroutine ImportDomes(File, GPSdomes, StartMJD, EndMJD)

! List of calling arguments:
character (len=*), intent(in)                :: File
type (StrucDomes), dimension(:), pointer     :: GPSdomes
real (wp), intent(in)                    :: StartMJD
real (wp), optional, intent(in)          :: EndMJD

! List of local variables:
type (StrucDomes), dimension(:), pointer     :: domes
type (StrucDomes)                            :: drow
character (len=256)          :: row
integer :: N
real (wp) :: eMJD

character(len=5)   :: dummy5
character(len=1)   :: dummy1
character(len=11)  :: StartStr, EndStr
integer (i2)       :: year, doy
integer (i4)       :: daysec
integer            :: OldID

! Input buffer
integer (i4) :: unit
integer (i4) :: iostat

logical :: debug = .false.
!---------------------------------------------------------------------

!debug = .true.

if (debug) write(*,*) 'ImportDomes> Start ...'

if (present(EndMJD)) then
   ! MJD was given, use the interger value
   if (StartMJD .gt. EndMJD) then
      write(*,*) 'ImportDomes> Error - StartMJD > EndMJD ',   &
                  StartMJD, ' > ',  EndMJD
      stop
   end if
   eMJD = EndMJD
else
   ! No EndMJD was given, don't read coordinates from a range but only for
   ! the given date => set EndMJDi to a negative number
   eMJD = -99.0_wp
end if

OldID = -99

unit = GetFreeUnit()
open( UNIT =   unit,                &
      FILE =   trim(File),          &
      STATUS = 'OLD',               &
      FORM =   'Formatted',         &
      ACTION = 'READ',              &
      IOSTAT = iostat                )

if (iostat .eq. 0) then
   ! File opened successfully

   ! Buffer for the maximum number of stations
   allocate( domes(1:9999) )

   N = 0
   row = ' '

   ! Skip file header
   do while (index(row,'+SITE/DOMES') .eq. 0)
      read(UNIT=unit,IOSTAT=iostat,FMT='(a)') row
      if (iostat .ne. 0) exit
   end do

   do while (iostat .eq. 0)
      ! read data, skip comments (rows starting with "*" or "%")
      read(UNIT=unit,IOSTAT=iostat,FMT='(a)') row
      ! Exit if row contains data
      if (scan(row(1:3),'*%') .eq. 0) exit
   end do


   ReadDomes: do while (iostat .eq. 0)
      ! Read domes data until end of file is reached

      !read(UNIT=unit,IOSTAT=iostat,                                      &
      read(row,'(tr3,i4,tr1,a4,tr1,a5,tr1,a1,tr1,a9,tr1,a22,tr1,a11,tr1,a11)')  &
          drow%ID, drow%SName, dummy5, dummy1, drow%domes,          &
          drow%Name, StartStr, EndStr

      if (debug) write(*,*) trim(row)

!!$      write(*,*) drow%ID, ' | ',         &
!!$                 drow%SName, ' | ',         &
!!$                 dummy5, ' | ',         &
!!$                 dummy1, ' | ',         &
!!$                 drow%domes,     ' | ',         &
!!$                 drow%Name, ' | ',         &
!!$                 StartStr, ' | ',         &
!!$                 EndStr

      ! Only the last two digits of the year are given
      ! this is required by SINEXepoch2MJD !!!
      read(StartStr(1:2),'(i2)')  year
      read(StartStr(3:5),'(i3)') doy
      read(StartStr(7:11),'(i5)') daysec

      !write(*,*) 'year, doy, daysec ', year, doy, daysec
      call SINEXepoch2MJD (year, doy, daysec, drow%StartMJD)
      !write(*,*) 'drow%StartMJD ', drow%StartMJD

      if (index(EndStr,'00000:00000') .ne. 0) then
         ! For the last entries the expiration date is not yet known
         ! and EndStr = '00000:00000'. To be compatible with the HORI
         ! file, EndMJD is set to 99999.0 (i.e. 31.8.2132):
         drow%EndMJD = 99999.0_wp
      else
         read(EndStr(1:2),'(i2)')  year
         read(EndStr(3:5),'(i3)') doy
         read(EndStr(7:11),'(i5)') daysec
         call SINEXepoch2MJD (year, doy, daysec, drow%EndMJD)
      end if
      !write(*,*) 'drow%EndMJD ', drow%EndMJD

      drow%Country = adjustl(drow%Name(scan(drow%Name,',')+1:))
      drow%Name    = drow%Name(1:scan(drow%Name,',')-1)
      !write(*,*) 'drow%Country = ', drow%Country
      !write(*,*) 'drow%Name    = ', drow%Name

      if (debug) write(*,*) drow%ID, drow%SName, ' ', drow%domes, ' ', &
                            trim(drow%Name), ' ', trim(drow%Country),  &
                             ' ', drow%StartMJD, drow%EndMJD

      if (drow%StartMJD .le. StartMJD .and.          &
          drow%EndMJD   .ge. eMJD          ) then
         ! The coordinates valid at StartMJD have been found

         N = N + 1
         domes(N) = drow
         OldID = drow%ID

      else if ( ( drow%StartMJD .ge. StartMJD    .and.               &
                  drow%StartMJD .le.  eMJD             ) .or.        &
                ( drow%EndMJD   .ge.  StartMJD   .and.               &
                  drow%EndMJD   .le.  eMJD             )     ) then
         ! Station coordinates changed between MJDi and EndMJDi, use
         ! station

         if ( drow%ID .ne. OldID) then
            ! Only one entry for each station in station list:
            ! use last valid record
            N = N + 1
            domes(N) = drow
            OldID = drow%ID
         else
            ! Valid record with the same ID: Replace older record (same N)
            domes(N) = drow
         end if

      end if

       do while (iostat .eq. 0)
          ! read data, skip comments (rows starting with "*" or "%")
          read(UNIT=unit,IOSTAT=iostat,FMT='(a)') row
          ! Exit read loop at the end of the "SITE/DOMES" section
          if (index(row,'-SITE/DOMES') .ne. 0) exit ReadDomes
          ! Exit if row contains data
          if (scan(row(1:3),'*%') .eq. 0) exit
       end do

    end do ReadDomes

   close(unit)

   if( N .gt. 0) then
      ! Copy station data to an array of length N and deallocate the buffer

      ! Deallocate any old data from previous call to this routine
      if (associated(GPSdomes)) deallocate(GPSdomes)
      allocate( GPSdomes(1:N) )
      GPSdomes = domes(1:N)
      deallocate(domes)
      if (debug) write(*,*) 'ImportDomes> Number of stations read ', N
   else
      ! No valid data found, deallocate buffer
      deallocate(domes)
      if (debug) write(*,*) 'ImportDomes> No valid stations found'
   end if

else
   write(*,*) 'ImportDomes> Error opening file ', trim(File)
end if

if (debug) write(*,*) 'ImportDomes> ... end'

end subroutine ImportDomes
!****
! <<<<<<<< End of RoboDoc comments


!---------------------------------------------------------------------
! function GetFreeUnit
!---------------------------------------------------------------------
!
! Start of RoboDoc comments >>>>
!****F* Global/GetFreeUnit
!
! Name
! GetFreeUnit
!
! Call
! unit = GetFreeUnit
!
! Purpose
! Find the next free unit before opening file
!
! Input
! None
!
! Output
! GetFreeUnit - free unit, return value
!
! External References
! None
!
! Record of revisions
! Date        Programmer   Description of change
! ====        ==========   =====================
! 24.04.2009  M. Bender    new
!---------------------------------------------------------------------
! Source
!---------------------------------------------------------------------
function GetFreeUnit ()

! List of calling arguments:
integer :: GetFreeUnit

! List of local variables:
integer :: unit
logical :: opened
!---------------------------------------------------------------------

unit = 30
opened = .true.
do while (opened)
   unit = unit + 1
   inquire(UNIT = unit, OPENED = opened )
end do

GetFreeUnit = unit

end function GetFreeUnit
!****
! <<<<<<<< End of RoboDoc comments


!---------------------------------------------------------------------
! subroutine SINEXepoch2MJD
!---------------------------------------------------------------------
!
! Start of RoboDoc comments >>>>
!****S* GPSNav/SINEXepoch2MJD
!
! Name
! SINEXepoch2MJD
!
! Purpose
! Convert a SINEX format epoch into a Modified Julian Date
!
! SINEX Documentation:
! http://tau.fesg.tu-muenchen.de/~iers/web/sinex/sinex_v201_appendix1.pdf
!
! "The first header line and most blocks are related through epochs or time
! stamps in the following format:
! YY:DOY:SECOD  -- YY-year; DOY- day of year; SECOD -sec of day;
! E.g. the epoch
! 95:120:86399 denotes April 30, 1995 (23:59:59UT).
! The epochs 00:00:00000 are allowed in all blocks, except the first
! header line."
!
! YY:DDD:SSSSS. "UTC"
! YY = last 2 digits of the year,
! if YY <= 50 implies 21-st century,
! if YY > 50 implies 20-th century.
!
! Input
! year     - year               |
! doy      - day of year        | SINEX epoch: year:doy:secday
! secday   - seconds of day     |
!
! Output
! MJD      - Modified Julian Date
!
! External References
! None
!
! Record of revisions
! Date        Programmer   Description of change
! ====        ==========   =====================
! 26.06.2007  M. Bender    new
! 29.04.2009  M. Bender    "DOY2GREG" cannot be called with "intent(in)"
!                          variable "doy"
!---------------------------------------------------------------------
! Source
!---------------------------------------------------------------------
subroutine SINEXepoch2MJD (SINEXyear, doy, secday, MJD)

! List of calling arguments:
integer (i2), intent(in)     :: SINEXyear, doy
integer (i4), intent(in)     :: secday
real (wp), intent(out)       :: MJD

! List of local variables:
integer                  :: Year, Month, Day
integer                  :: Hour, Min
real (wp)                :: Sec, JD
integer                  :: M, D
integer                  :: doy2, year2
!---------------------------------------------------------------------

if (SINEXyear .le. 50) then
   ! if YY <= 50 implies 21-st century
   year = 2000 + SINEXyear
else
   ! if YY > 50 implies 20-th century
      year = 1900 + SINEXyear
end if
!write(*,*) year

!call DOY2GREG(int(doy), int(year), M, D)
doy2 = doy
year2 = year
call DOY2GREG(doy2, year2, M, D)

call secday2time (real(SecDay,wp), Hour, Min, Sec)
!write(*,*) Year, Month, Day,  Hour, Min, Sec
!write(*,*) Year, M, D,  Hour, Min, Sec
Month = M
Day = D
JD = JulianDate(Year, Month, Day, Hour, Min, Sec)
MJD = JD2MJD( JD )

end subroutine SINEXepoch2MJD


subroutine secday2time (SecDay, Hour, Min, Sec)

! List of calling arguments:
real (wp), intent(in)        :: SecDay
integer, intent(out)             :: Hour, Min
real (wp), intent(out)       :: Sec
!---------------------------------------------------------------------

Hour = int(SecDay/3600.0)
Min  = int((SecDay-Hour*3600)/60.0)
Sec  = mod(SecDay,60.0_wp)

end subroutine secday2time
!****
! <<<<<<<< End of RoboDoc comments


subroutine Ellips2Cart_tl( lambda, beta, height, X, Y, Z, Ellips,  &
                           lambda_tl, beta_tl, height_tl,          &
                           X_tl, Y_tl, Z_tl)

! List of calling arguments:
real (wp), intent(in) :: lambda, beta, height
real (wp), intent(in) :: lambda_tl, beta_tl, height_tl
real (wp), intent(out) :: X, Y, Z
real (wp), intent(out) :: X_tl, Y_tl, Z_tl

type (EllipsParam) :: Ellips

! List of local variables:
real (wp) :: NK, dNK
!---------------------------------------------------------------------

! Querkruemmungsradius des Ellipsoids
NK = Ellips%a/sqrt(1.0 - Ellips%EN1sup2*(sin(beta))**2)
dNK = beta_tl * (Ellips%a*Ellips%EN1sup2*sin(beta)*cos(beta)) /   &
                sqrt(1.0 - Ellips%EN1sup2*(sin(beta))**2)**3

! Berechnung der kartesischen Koordinaten
X = (NK+height)*cos(beta)*cos(lambda)
Y = (NK+height)*cos(beta)*sin(lambda)
Z = ((NK/(1+Ellips%EN2sup2)) + height)*sin(beta)

X_tl =   (dNK + height_tl) * cos(beta)*cos(lambda)      &
       - beta_tl *  (NK+height)*sin(beta)*cos(lambda)   &
       - lambda_tl *  (NK+height)*cos(beta)*sin(lambda)

Y_tl =   (dNK + height_tl) * cos(beta)*sin(lambda)      &
       - beta_tl * (NK+height)*sin(beta)*sin(lambda)    &
       + lambda_tl * (NK+height)*cos(beta)*cos(lambda)

Z_tl =   height_tl * sin(beta)     &
       + dNK *  sin(beta)/(1+Ellips%EN2sup2) &
       + beta_tl *  ((NK/(1+Ellips%EN2sup2)) + height)*cos(beta)

end subroutine Ellips2Cart_tl


subroutine Ellips2Cart_ad( lambda, beta, height, Ellips,  &
                           lambda_ad, beta_ad, height_ad,          &
                           X_ad, Y_ad, Z_ad)

! List of calling arguments:
real (wp), intent(in) :: lambda, beta, height
!real (wp), intent(out) :: X, Y, Z
type (EllipsParam) :: Ellips
real (wp), intent(inout) :: X_ad, Y_ad, Z_ad
real (wp), intent(out) :: lambda_ad, beta_ad, height_ad



! List of local variables:
real (wp) :: NK, dNK
!---------------------------------------------------------------------

! Querkruemmungsradius des Ellipsoids
NK = Ellips%a/sqrt(1.0 - Ellips%EN1sup2*(sin(beta))**2)

! Berechnung der kartesischen Koordinaten
!X = (NK+height)*cos(beta)*cos(lambda)
!Y = (NK+height)*cos(beta)*sin(lambda)
!Z = ((NK/(1+Ellips%EN2sup2)) + height)*sin(beta)

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++

! dZ/dh, ...
height_ad = Z_ad * sin(beta)
beta_ad   = Z_ad * ((NK/(1+Ellips%EN2sup2)) + height)*cos(beta)
dNK       = Z_ad * sin(beta)/(1+Ellips%EN2sup2)

! dY/dh, ...
lambda_ad =  Y_ad * (NK+height)*cos(beta)*cos(lambda)
beta_ad   = beta_ad -  (NK+height)*sin(beta)*sin(lambda) * Y_ad
height_ad = height_ad + cos(beta)*sin(lambda) * Y_ad
dNK       = dNK + cos(beta)*sin(lambda) * Y_ad

! dX/dh, ...
lambda_ad = lambda_ad - (NK+height)*cos(beta)*sin(lambda) * X_ad
beta_ad   = beta_ad - (NK+height)*sin(beta)*cos(lambda) * X_ad
height_ad = height_ad + cos(beta)*cos(lambda) * X_ad
dNK       = dNK + cos(beta)*cos(lambda) * X_ad

! dNK/d beta
beta_ad   = beta_ad +  dNK * (Ellips%a*Ellips%EN1sup2*sin(beta)*cos(beta)) /   &
                sqrt(1.0 - Ellips%EN1sup2*(sin(beta))**2)**3

X_ad = 0.0_wp
Y_ad = 0.0_wp
Z_ad = 0.0_wp

end subroutine Ellips2Cart_ad


subroutine Cart2Ellips_tl(X, Y, Z, lambda, beta, height, Ellips,  &
                          X_tl, Y_tl, Z_tl,                       &
                          lambda_tl, beta_tl, height_tl)

! List of calling arguments:
real (wp), intent(in) :: X, Y, Z
real (wp), intent(out) :: lambda, beta, height
real (wp), intent(in) :: X_tl, Y_tl, Z_tl
real (wp), intent(out) :: lambda_tl, beta_tl, height_tl

type (EllipsParam) :: Ellips

! List of local variables:
real (wp) :: rho, phi, gamma, g2
real (wp) :: N, theta, h
!integer :: i

real (wp) :: drho, dtheta, dphi, dg2, dgamma, dN, dh
real (wp) :: c1, c2
!---------------------------------------------------------------------

! Einige Hilfsgroessen:
rho = sqrt(X**2+Y**2)
drho = (X_tl * X +  Y_tl * Y) / rho

c1 = (Z*Ellips%a)/(rho*Ellips%b)
if (rho .gt. 1.0E-6_wp) then
   theta = atan(c1)
   dtheta = Z_tl * Ellips%a / (rho*Ellips%b*(1.0_wp+c1**2)) -         &
            drho * Z * Ellips%a /  (rho**2 * Ellips%b*(1.0_wp+c1**2))
else
   theta = pi05
   dtheta = 0.0_wp
end if

! Die ellipsoidische Laenge lambda:
if (abs(rho) > reps) then
   phi = 2.0_wp*atan(Y/(abs(X)+rho))
   dphi = ( -sign(2.0_wp*Y,X) * X_tl + 2.0_wp*(abs(X) + rho) * Y_tl   -    &
            2.0_wp*Y * drho ) / ((abs(X)+rho)**2 + Y**2)
end if
if (X .ge. 0.0_wp) then
   if (abs(X) .lt. 1.0E-8_wp .and. abs(Y) .lt. 1.0E-8_wp) then
      ! if X=0 and Y=0
      lambda = 0.0_wp
      lambda_tl = 0.0_wp
   else if (Y .ge. 0.0_wp) then
      ! X>0 and Y>0
      lambda = phi
      lambda_tl = dphi
   else
      ! X>0 and Y<0
      lambda = phi + 2.0_wp*pi
      lambda_tl = dphi
   end if
else  ! X < 0
   lambda = pi - phi
   lambda_tl = -dphi
end if

! Die ellipsoidische Breite beta:
g2 = rho - Ellips%EN1sup2*Ellips%a*(cos(theta))**3  ! Nenner von gamma
dg2 = drho + dtheta *                                                      &
             3.0_wp * Ellips%EN1sup2*Ellips%a * sin(theta) * cos(theta)**2
if (abs(g2) > reps) then
   gamma = atan( (Z + Ellips%EN2sup2*Ellips%b*(sin(theta))**3) / g2 )
   c2 = 1.0_wp/(1.0_wp +((Z+Ellips%EN2sup2*Ellips%b*(sin(theta))**3)/g2)**2)
   dgamma = c2 * (Z_tl/g2 -                                                 &
                  dg2 * (Z+Ellips%EN2sup2*Ellips%b*(sin(theta)**3))/g2**2 + &
                  dtheta * 3*Ellips%EN2sup2*Ellips%b*cos(theta) *           &
                           (sin(theta)**2)/g2 )

! else ???? +pi/2 oder -pi/2 ???
end if
if (abs(rho) .gt. 1.0E-8_wp) then
   ! rho > 0
   beta = gamma
   beta_tl = dgamma
else  ! rho = 0
   if (Z .gt. 0.0_wp) then
      beta = 0.50_wp*pi
      beta_tl = 0.0_wp
   else if (abs(Z) < reps) then
      beta = 0.0_wp
      beta_tl = 1.0_wp
   else if (Z .lt. 0.0_wp) then
      beta = -0.50_wp*pi
      beta_tl = 0.0_wp
   end if
end if

! Bestimmung des Querkruemmungsradius fuer diese Breite:
N = Ellips%a/sqrt(1.0_wp - Ellips%EN1sup2*((sin(beta))**2))
dN = beta_tl *  (Ellips%a*Ellips%EN1sup2*sin(beta)*cos(beta)) /     &
                 sqrt(1.0_wp - Ellips%EN1sup2*((sin(beta))**2))**3
!write(*,*) 'N = ', N, (sin(beta))**2

! Die ellipsoidische Hoehe h:
h = (rho/cos(beta)) - N
dh = -dN + drho/cos(beta) + beta_tl*rho*sin(beta)/(cos(beta)**2)
!write(*,*) 'h = ', h, cos(beta)

if (abs(rho) > reps) then
   height = h
   height_tl = dh
else ! rho = 0
   if (abs(Z) < reps) then
      height = -Ellips%b
   else
      height = abs(Z) - Ellips%b
   end if
end if

height_tl = dh ! ?????????

end subroutine Cart2Ellips_tl


subroutine Cart2Ellips_ad(X, Y, Z, lambda, beta, height, Ellips,   &
                          X_ad, Y_ad, Z_ad,                        &
                          lambda_ad, beta_ad, height_ad)

! List of calling arguments:
real (wp), intent(in) :: X, Y, Z
real (wp), intent(out) :: lambda, beta, height
real (wp), intent(out) :: X_ad, Y_ad, Z_ad
real (wp), intent(inout) :: lambda_ad, beta_ad, height_ad

type (EllipsParam) :: Ellips

! List of local variables:
real (wp) :: rho, phi, gamma, g2
real (wp) :: N, theta, h, c1, c2
!integer :: i

real (wp) :: dh, drho, dN, dgamma, dtheta, dg2, dphi
!---------------------------------------------------------------------

! Einige Hilfsgroessen:
rho = sqrt(X**2+Y**2)

if (rho .gt. 1.0E-6_wp) then
   theta = atan((Z*Ellips%a)/(rho*Ellips%b))
else
   theta = pi05
end if

! Die ellipsoidische Laenge lambda:
if (abs(rho) > reps) then
   phi = 2.0_wp*atan(Y/(abs(X)+rho))
end if
if (X .ge. 0.0_wp) then
   if (abs(X) .lt. 1.0E-8_wp .and. abs(Y) .lt. 1.0E-8_wp) then
      ! if X=0 and Y=0
      lambda = 0.0_wp
   else if (Y .ge. 0.0_wp) then
      ! X>0 and Y>0
      lambda = phi
   else
      ! X>0 and y<0
      lambda = phi + 2.0_wp*pi
   end if
else  ! X < 0
   lambda = pi - phi
end if

! Die ellipsoidische Breite beta:
g2 = rho - Ellips%EN1sup2*Ellips%a*(cos(theta))**3  ! Nenner von gamma
if (abs(g2) > reps) then
   gamma = atan( (Z + Ellips%EN2sup2*Ellips%b*(sin(theta))**3) / g2 )
! else ???? +pi/2 oder -pi/2 ???
end if
if (abs(rho) .gt. 1.0E-8_wp) then
   ! rho > 0
   beta = gamma
else  ! rho = 0
   if (Z .gt. 0.0_wp) then
      beta = 0.50_wp*pi
   else if (abs(Z) < reps) then
      beta = 0.0_wp
   else if (Z .lt. 0.0_wp) then
      beta = -0.50_wp*pi
   end if
end if

! Bestimmung des Querkruemmungsradius fuer diese Breite:
N = Ellips%a/sqrt(1.0_wp - Ellips%EN1sup2*((sin(beta))**2))
!write(*,*) 'N = ', N, (sin(beta))**2

! Die ellipsoidische Hoehe h:
h = (rho/cos(beta)) - N
!write(*,*) 'h = ', h, cos(beta)

if (abs(rho) > reps) then
   height = h
else ! rho = 0
   if (abs(Z) < reps) then
      height = -Ellips%b
   else
      height = abs(Z) - Ellips%b
   end if
end if

!=====================================================================
c1 = (Z*Ellips%a)/(rho*Ellips%b)

dh = 0.0_wp
if (abs(rho) > reps) then
   height = h
   dh =  height_ad
else ! rho = 0
   if (abs(Z) < reps) then
      height = -Ellips%b
   else
      height = abs(Z) - Ellips%b
      Z_ad = sign(1.0_wp,Z) * height_ad
   end if
end if

! Die ellipsoidische Hoehe h:
!h = (rho/cos(beta)) - N
drho = dh/cos(beta)
dN   = -dh
beta_ad =  beta_ad + dh * rho * sin(beta) / (cos(beta))**2

! Bestimmung des Querkruemmungsradius fuer diese Breite:
!N = Ellips%a/sqrt(1.0_wp - Ellips%EN1sup2*((sin(beta))**2))
beta_ad = beta_ad + dN * (Ellips%a*Ellips%EN1sup2*sin(beta)*cos(beta)) /     &
                 sqrt(1.0_wp - Ellips%EN1sup2*((sin(beta))**2))**3

dgamma = 0.0_wp
if (abs(rho) .gt. 1.0E-8_wp) then
   ! rho > 0
   beta = gamma
   dgamma = beta_ad
else  ! rho = 0
   if (Z .gt. 0.0_wp) then
      beta = 0.50_wp*pi
   else if (abs(Z) < reps) then
      beta = 0.0_wp
   else if (Z .lt. 0.0_wp) then
      beta = -0.50_wp*pi
   end if
end if

c2 = 1.0_wp/(1.0_wp +((Z+Ellips%EN2sup2*Ellips%b*(sin(theta))**3)/g2)**2)
if (abs(g2) > reps) then
   !gamma = atan( (Z + Ellips%EN2sup2*Ellips%b*(sin(theta))**3) / g2 )
   Z_ad = Z_ad + dgamma * c2/g2
   dtheta = dgamma * c2*3*Ellips%EN2sup2*Ellips%b*cos(theta) *           &
                           (sin(theta)**2)/g2
   dg2 =  -dgamma *  c2*(Z+Ellips%EN2sup2*Ellips%b*(sin(theta)**3))/g2**2
end if

!g2 = rho - Ellips%EN1sup2*Ellips%a*(cos(theta))**3  ! Nenner von gamma
drho = drho + dg2
dtheta = dtheta + dg2 *                                                 &
          3.0_wp * Ellips%EN1sup2*Ellips%a * sin(theta) * cos(theta)**2

dphi = 0.0_wp
if (X .ge. 0.0_wp) then
   if (abs(X) .lt. 1.0E-8_wp .and. abs(Y) .lt. 1.0E-8_wp) then
      ! if X=0 and Y=0
      lambda = 0.0_wp
   else if (Y .ge. 0.0_wp) then
      ! X>0 and Y>0
      lambda = phi
      dphi = lambda_ad
   else
      ! X>0 and y<0
      lambda = phi + 2.0_wp*pi
      dphi = lambda_ad
   end if
else  ! X < 0
   lambda = pi - phi
   dphi = -lambda_ad
end if

X_ad = 0.0_wp
if (abs(rho) > reps) then
   !phi = 2.0_wp*atan(Y/(abs(X)+rho))
   Y_ad = dphi *  2.0_wp*(abs(X) + rho) / ((abs(X)+rho)**2 + Y**2)
   X_ad = X_ad - dphi * (sign(2.0_wp*Y,X) / ((abs(X)+rho)**2 + Y**2))
   drho = drho  - dphi *  2.0_wp*Y  / ((abs(X)+rho)**2 + Y**2)
end if

if (rho .gt. 1.0E-6_wp) then
   !theta = atan((Z*Ellips%a)/(rho*Ellips%b))
   Z_ad = Z_ad + dtheta * Ellips%a / (rho*Ellips%b*(1.0_wp+c1**2))
   drho = drho - dtheta * Z * Ellips%a /  (rho**2 * Ellips%b*(1.0_wp+c1**2))
else
   theta = pi05
end if

!rho = sqrt(X**2+Y**2)
X_ad = X_ad + drho * X / rho
Y_ad = Y_ad + drho * Y / rho

!write(*,*) 'rho, drho : ', rho, drho

lambda_ad = 0.0_wp
beta_ad   = 0.0_wp
height_ad = 0.0_wp

end subroutine Cart2Ellips_ad


subroutine LocalHorz2Cart_tl(x2, y2, z2, lambda, phi,               &
                          ox, oy, oz, x1, y1, z1,                   &
                          x2_tl, y2_tl, z2_tl, lambda_tl, phi_tl,   &
                          ox_tl, oy_tl, oz_tl, x1_tl, y1_tl, z1_tl)

! List of calling arguments:
real (wp), intent(in)  :: x2, y2, z2
real (wp), intent(in)  :: lambda, phi
real (wp), intent(in)  :: ox, oy, oz
real (wp), intent(out) :: x1, y1, z1
! diff
real (wp), intent(in)  :: x2_tl, y2_tl, z2_tl
real (wp), intent(in)  :: lambda_tl, phi_tl
real (wp), intent(in)  :: ox_tl, oy_tl, oz_tl
real (wp), intent(out) :: x1_tl, y1_tl, z1_tl

! List of local variables:
real (wp) :: a, b, c, d
!---------------------------------------------------------------------

! Originalcode LocalHorz2Cart
x1 = ox -sin(phi)*cos(lambda)*x2 -      &
         sin(lambda)*y2          +      &
         cos(phi)*cos(lambda)*z2
y1 = oy -sin(phi)*sin(lambda)*x2 +      &
         cos(lambda)*y2          +      &
         cos(phi)*sin(lambda)*z2
z1 = oz +cos(phi)*x2             +      &
         sin(phi)*z2

! tangent linearer Code:
a = sin(phi) * cos(lambda)
b = cos(phi) * cos(lambda)
c = sin(phi) * sin(lambda)
d = cos(phi) * sin(lambda)

x1_tl = ox_tl - x2_tl * a - y2_tl * sin(lambda) + z2_tl * b  &
        - phi_tl * (x2*b + z2*a)                             &
        + lambda_tl * (x2*c - y2*cos(lambda) - z2*d)

y1_tl = oy_tl - x2_tl * c + y2_tl * cos(lambda) + z2_tl * d  &
        - phi_tl * (x2*d + z2*c)                       &
        - lambda_tl * (x2*a + y2*sin(lambda) - z2*b)

z1_tl = oz_tl + x2_tl*cos(phi) + z2_tl*sin(phi)  &
        - phi_tl * (x2*sin(phi) - z2*cos(phi))

end subroutine LocalHorz2Cart_tl


subroutine LocalHorz2Cart_ad(x2, y2, z2, lambda, phi, ox, oy, oz,    &
                             x1_ad, y1_ad, z1_ad,                    &
                             x2_ad, y2_ad, z2_ad, lambda_ad, phi_ad, &
                             ox_ad, oy_ad, oz_ad )

! List of calling arguments:
real (wp), intent(in)  :: x2, y2, z2
real (wp), intent(in)  :: lambda, phi
real (wp), intent(in)  :: ox, oy, oz
real (wp), intent(out) :: x2_ad, y2_ad, z2_ad
real (wp), intent(out) :: lambda_ad, phi_ad
real (wp), intent(out) :: ox_ad, oy_ad, oz_ad
!real (wp), intent(out) :: x1, y1, z1
real (wp), intent(inout) :: x1_ad, y1_ad, z1_ad

! List of local variables:
real (wp) :: a, b, c, d
!---------------------------------------------------------------------

!!$x1 = ox -sin(phi)*cos(lambda)*x2 -      &
!!$         sin(lambda)*y2          +      &
!!$         cos(phi)*cos(lambda)*z2
!!$y1 = oy -sin(phi)*sin(lambda)*x2 +      &
!!$         cos(lambda)*y2          +      &
!!$         cos(phi)*sin(lambda)*z2
!!$z1 = oz +cos(phi)*x2             +      &
!!$         sin(phi)*z2

a = sin(phi) * cos(lambda)
b = cos(phi) * cos(lambda)
c = sin(phi) * sin(lambda)
d = cos(phi) * sin(lambda)

!!$ x2_ad = 0.0_wp
!!$ y2_ad = 0.0_wp
!!$phi_ad = 0.0_wp
!!$lambda_ad = 0.0_wp

!!$! z1 = oz +cos(phi)*x2 + sin(phi)*z2
oz_ad = z1_ad
x2_ad = z1_ad * cos(phi)
z2_ad = z1_ad * sin(phi)
phi_ad = z1_ad * (-x2*sin(phi) + z2*cos(phi))

! y1=oy-sin(phi)*sin(lambda)*x2+cos(lambda)*y2+cos(phi)*sin(lambda)*z2
oy_ad = y1_ad
x2_ad = x2_ad - y1_ad * c
y2_ad = y1_ad * cos(lambda)
z2_ad = z2_ad + y1_ad * d
phi_ad = phi_ad - y1_ad * (x2*d + z2*c)
lambda_ad = -y1_ad * (x2*a + y2*sin(lambda) - z2*b)

! x1=ox-sin(phi)*cos(lambda)*x2-sin(lambda)*y2+cos(phi)*cos(lambda)*z2
ox_ad = x1_ad
x2_ad = x2_ad - x1_ad * a
y2_ad = y2_ad - x1_ad * sin(lambda)
z2_ad = z2_ad + x1_ad * b
phi_ad = phi_ad - x1_ad * (x2*b + z2*a)
lambda_ad = lambda_ad + x1_ad * (x2*c -y2*cos(lambda) - z2*d)

x1_ad = 0.0_wp
y1_ad = 0.0_wp
z1_ad = 0.0_wp

end subroutine LocalHorz2Cart_ad


subroutine Azimut2Horz_tl(A, Z, dist, x2, y2, z2,                    &
                          A_tl, Z_tl, dist_tl, x2_tl, y2_tl, z2_tl)

! List of calling arguments:
real (wp), intent(in)  :: A, Z, dist
real (wp), intent(in)  :: A_tl, Z_tl, dist_tl
real (wp), intent(out) :: x2, y2, z2
real (wp), intent(out) :: x2_tl, y2_tl, z2_tl

! List of local variables:
!---------------------------------------------------------------------

! Original code
x2 = dist * cos(A) * sin(Z)
y2 = dist * sin(A) * sin(Z)
z2 = dist * cos(Z)

! Tangent linear code
x2_tl = dist_tl * cos(A) * sin(Z) - A_tl * dist *sin(A) * sin(Z) + &
        Z_tl * dist * cos(A) * cos(Z)
y2_tl = dist_tl * sin(A) * sin(Z) + A_tl * dist * cos(A) * sin(Z) + &
        Z_tl * dist * sin(A) * cos(Z)
z2_tl = dist_tl * cos(Z) - Z_tl * dist * sin(Z)

end subroutine Azimut2Horz_tl


subroutine Azimut2Horz_ad(A, Z, dist, x2, y2, z2,     &
                          A_ad, Z_ad, dist_ad, x2_ad, y2_ad, z2_ad)

! List of calling arguments:
real (wp), intent(in)  :: A, Z, dist
real (wp), intent(out) :: A_ad, Z_ad, dist_ad
real (wp), intent(out) :: x2, y2, z2
real (wp), intent(inout) :: x2_ad, y2_ad, z2_ad

! List of local variables:
!---------------------------------------------------------------------

! z2 = dist * cos(Z)
dist_ad =  z2_ad * cos(Z)
Z_ad    = -z2_ad * dist * sin(Z)

! y2 = dist * sin(A) * sin(Z)
dist_ad = dist_ad + y2_ad * sin(A) * sin(Z)
A_ad    =           y2_ad * dist * sin(Z) * cos(A)
Z_ad    = Z_ad    + y2_ad * dist * sin(A) * cos(Z)

! x2 = dist * cos(A) * sin(Z)
dist_ad = dist_ad + x2_ad * cos(A) * sin(Z)
A_ad    = A_ad    - x2_ad * dist * sin(A) * sin(Z)
Z_ad    = Z_ad    + x2_ad * dist * cos(A) * cos(Z)

x2_ad = 0.0_wp
y2_ad = 0.0_wp
z2_ad = 0.0_wp

end subroutine Azimut2Horz_ad


!---------------------------------------------------------------------
! function IntegPolyCube
!---------------------------------------------------------------------
!
! Start of RoboDoc comments >>>>
!****F* Integration/IntegPolyCube
!
! Name
! IntegPolyCube
!
! Call
! integ = IntegPolyCube ( fx )
!
! Purpose
! Numerical integration of unequally spaced data
! A cubic interpolation polynomial with 4 supporting points is used to
! estimate the integral:
!
! P. E. Gill and G. F. Miller
! An algorithm for the integration of unequally spaced data
! The Computer Journal, Vol. 15, No. 1, p. 80-83, 1972
! http://comjnl.oxfordjournals.org/content/15/1/80.abstract
!
! The integral is always computed over the whole interval covered by the
! data. To obtain integrals over subintervals [ia, ib] call the function
! with the corresponding subinterval:
! integ = IntegPolyCube (fx(ia:ib))
!
! Warning:
! The interpolation polynomial doea work well only if the spacing of
! the data is is not too inhomogeneous. Single intervals with much larger
! spacings (> 100 times) can lead to very wrong results.
!
! Input
! fx(i,j) - i=1, n, index of nodes, n - number of array elements
!           j=1 - x, supporting points
!           j=2 - y, function y = f(x)
!
! Output
! IntegPolyCube (return value) - integral
!
!
! Record of revisions
! Date        Programmer   Description of change
! ====        ==========   =====================
! 04.03.2013  M. Bender    new
!---------------------------------------------------------------------
! Source
!---------------------------------------------------------------------
function IntegPolyCube ( fx )

! List of calling arguments:
real (wp)                              :: IntegPolyCube
real (wp), dimension(:,:), intent(in)  :: fx

!  List of local variables:
real (wp) :: sum, a, b
real (wp) :: Dx0, Dx1, Dx2, DDx1, DDx2, DDx3
real (wp) :: Df0, Df1, Df2, DDf1, DDf2, DDf3
integer              :: i, j, k,  n
!---------------------------------------------------------------------

j = lbound(fx,1)
k = ubound(fx,1)
n = k - j + 1

if (n .lt. 4) then
   ! At least 4 points are required to compute the integral
   write(*,*) 'IntegPolyCube> More than three points are required, n = ', n
   IntegPolyCube = 0.0_wp
else

   sum = 0.0_wp
   do i=j, k-3

      if (i .eq. j) then
         ! Compute first set of differences and the integral over the
         ! first interval:

         ! Delta x ( x => fx(i,1) )
         Dx0 = fx(i+1,1) - fx(i,1)
         Dx1 = fx(i+2,1) - fx(i+1,1)
         Dx2 = fx(i+3,1) - fx(i+2,1)
         DDx1 = Dx0 + Dx1
         DDx2 = Dx1 + Dx2
         DDx3 = DDx2 + Dx0
         ! Delta y / Delta x  ( y => fx(i,2) )
         Df0 = ( fx(i+1,2) - fx(i,2)   ) / Dx0
         Df1 = ( fx(i+2,2) - fx(i+1,2) ) / Dx1
         Df2 = ( fx(i+3,2) - fx(i+2,2) ) / Dx2
         DDf1 = ( Df1 - Df0 ) / DDx1
         DDf2 = ( Df2 - Df1 ) / DDx2
         DDf3 = ( DDf2 - DDf1 ) / DDx3

         ! Integral
         a = (Dx0 + 2.0_wp*Dx1) * DDf3 / 12.0_wp
         b = DDf1 / 6.0_wp
         sum = Dx0 * ( fx(i,2) + 0.50_wp*Dx0*Df0 +    &
                       (a - b) * Dx0**2           )

      end if

      ! Compute integral using 4 point formula

      ! Delta x ( x => fx(i,1) )
      !Dx0 = fx(i+1,1) - fx(i,1)
      !Dx1 = fx(i+2,1) - fx(i+1,1)
      Dx2 = fx(i+3,1) - fx(i+2,1)
      !DDx1 = Dx0 + Dx1
      DDx2 = Dx1 + Dx2
      DDx3 = DDx2 + Dx0
      ! Delta y / Delta x  ( y => fx(i,2) )
      !Df0 = ( fx(i+1,2) - fx(i,2)   ) / Dx0
      !Df1 = ( fx(i+2,2) - fx(i+1,2) ) / Dx1
      Df2 = ( fx(i+3,2) - fx(i+2,2) ) / Dx2
      !DDf1 = ( Df1 - Df0 ) / DDx1
      DDf2 = ( Df2 - Df1 ) / DDx2
      DDf3 = ( DDf2 - DDf1 ) / DDx3

      ! Integral
      a = 0.50_wp * ( DDf1 + DDf2 )
      b = 0.50_wp * ( Dx0 - Dx2 ) * DDf3
      sum = sum + Dx1 * ( 0.50_wp * ( fx(i+1,2) + fx(i+2,2) ) -       &
                          Dx1**2 * (a + b) / 6.0_wp )

      if (i .eq. k-3) then
         ! Compute  integral of the last interval using
         ! the same differences as above:

         ! Integral
         a = (Dx2 + 2.0_wp*Dx1) * DDf3 / 12.0_wp
         b = DDf2 / 6.0_wp
         sum = sum + Dx2 * ( fx(n,2) - 0.50_wp*Dx2*Df2 -    &
                             (a + b) * Dx2**2           )

         exit

      end if

      ! Reuse differences from last loop:
      Dx0  = Dx1
      Dx1  = Dx2
      DDx1 = DDx2
      Df0  = Df1
      Df1  = Df2
      DDf1 = DDf2

   end do

   IntegPolyCube = sum

end if

end function IntegPolyCube


!---------------------------------------------------------------------
! function IntegPolyCube_tl
!---------------------------------------------------------------------
!
! Start of RoboDoc comments >>>>
!****F* Integration/IntegPolyCube
!
! Name
! IntegPolyCube
!
! Call
! integ = IntegPolyCube ( fx )
!
! Purpose
! Numerical integration of unequally spaced data
! A cubic interpolation polynomial with 4 supporting points is used to
! estimate the integral:
!
! P. E. Gill and G. F. Miller
! An algorithm for the integration of unequally spaced data
! The Computer Journal, Vol. 15, No. 1, p. 80-83, 1972
! http://comjnl.oxfordjournals.org/content/15/1/80.abstract
!
! The integral is always computed over the whole interval covered by the
! data. To obtain integrals over subintervals [ia, ib] call the function
! with the corresponding subinterval:
! integ = IntegPolyCube (fx(ia:ib))
!
! Input
! fx(i,j) - i=1, n, index of nodes, n - number of array elements
!           j=1 - x, supporting points
!           j=2 - y, function y = f(x)
!
! Output
! IntegPolyCube (return value) - integral
!
!
! Record of revisions
! Date        Programmer   Description of change
! ====        ==========   =====================
! 07.03.2013  M. Bender    new
!---------------------------------------------------------------------
! Source
!---------------------------------------------------------------------
function IntegPolyCube_tl ( fx, fx_tl )

! List of calling arguments:
real (wp)                              :: IntegPolyCube_tl
real (wp), dimension (:,:), intent(in) :: fx, fx_tl

!  List of local variables:
real (wp) :: sum, a, b
real (wp) :: Dx0, Dx1, Dx2, DDx1, DDx2, DDx3
real (wp) :: Df0, Df1, Df2, DDf1, DDf2, DDf3
real (wp) :: tsum, ta, tb
real (wp) :: tDx0, tDx1, tDx2, tDDx1, tDDx2, tDDx3
real (wp) :: tDf0, tDf1, tDf2, tDDf1, tDDf2, tDDf3
integer   :: i, j, k,  n
!---------------------------------------------------------------------

j = lbound(fx,1)
k = ubound(fx,1)
n = k - j + 1

if (n .lt. 4) then
   ! At least 4 points are required to compute the integral
   write(*,*) 'IntegPolyCube> More than three points are required, n = ', n
   IntegPolyCube_tl = 0.0_wp
else

   sum = 0.0_wp
   do i=j, k-3

      if (i .eq. j) then
         ! Compute first set of differences and the integral over the
         ! first interval:

         ! Delta x ( x => fx(i,1) )
         Dx0 = fx(i+1,1) - fx(i,1)
         tDx0 =  fx_tl(i+1,1) - fx_tl(i,1)

         Dx1 = fx(i+2,1) - fx(i+1,1)
         tDx1 = fx_tl(i+2,1) - fx_tl(i+1,1)

         Dx2 = fx(i+3,1) - fx(i+2,1)
         tDx2 = fx_tl(i+3,1) - fx_tl(i+2,1)

         DDx1 = Dx0 + Dx1
         tDDx1 = tDx0 + tDx1

         DDx2 = Dx1 + Dx2
         tDDx2 = tDx1 + tDx2

         DDx3 = DDx2 + Dx0
         tDDx3 = tDDx2 + tDx0

         ! Delta y / Delta x  ( y => fx(i,2) )
         Df0 = ( fx(i+1,2) - fx(i,2)   ) / Dx0
         tDf0 = (fx_tl(i+1,2) - fx_tl(i,2))/Dx0 - &
                tDx0 * ( fx(i+1,2) - fx(i,2)   ) / Dx0**2

         Df1 = ( fx(i+2,2) - fx(i+1,2) ) / Dx1
         tDf1 = ( fx_tl(i+2,2) - fx_tl(i+1,2) ) / Dx1 -    &
                tDx1 * ( fx(i+2,2) - fx(i+1,2) ) / Dx1**2

         Df2 = ( fx(i+3,2) - fx(i+2,2) ) / Dx2
         tDf2 = ( fx_tl(i+3,2) - fx_tl(i+2,2) ) / Dx2 - &
                tDx2 * ( fx(i+3,2) - fx(i+2,2) ) / Dx2**2

         DDf1 = ( Df1 - Df0 ) / DDx1
         tDDf1 = ( tDf1 - tDf0 ) / DDx1 -          &
                  tDDx1 * ( Df1 - Df0 ) / DDx1**2

         DDf2 = ( Df2 - Df1 ) / DDx2
         tDDf2 = ( tDf2 - tDf1 ) / DDx2 -          &
                 tDDx2 * ( Df2 - Df1 ) / DDx2**2

         DDf3 = ( DDf2 - DDf1 ) / DDx3
         tDDf3 = ( tDDf2 - tDDf1 ) / DDx3 -           &
                 tDDx3 * ( DDf2 - DDf1 ) / DDx3**2

         ! Integral
         a = (Dx0 + 2.0_wp*Dx1) * DDf3 / 12.0_wp
         ta = (tDx0 * DDf3 + 2.0_wp*tDx1*DDf3 +    &
              (Dx0 + 2.0_wp*Dx1)*tDDf3) / 12.0_wp

         b = DDf1 / 6.0_wp
         tb = tDDf1 / 6.0_wp

         sum = Dx0 * ( fx(i,2) + 0.50_wp*Dx0*Df0 +    &
                       (a - b) * Dx0**2           )
         tsum = tDx0 * (fx(i,2) + Dx0*Df0 + 3.0_wp * (a - b) * Dx0**2) +  &
                fx_tl(i,2)*Dx0 + 0.50_wp*tDf0*Dx0**2 + (ta-tb)*Dx0**3

      end if

      ! Compute integral using 4 point formula

      ! Delta x ( x => fx(i,1) )
      !Dx0 = fx(i+1,1) - fx(i,1)
      !Dx1 = fx(i+2,1) - fx(i+1,1)
      Dx2 = fx(i+3,1) - fx(i+2,1)
      tDx2 = fx_tl(i+3,1) - fx_tl(i+2,1)

      !DDx1 = Dx0 + Dx1
      DDx2 = Dx1 + Dx2
      tDDx2 = tDx1 + tDx2

      DDx3 = DDx2 + Dx0
      tDDx3 = tDDx2 + tDx0

      ! Delta y / Delta x  ( y => fx(i,2) )
      !Df0 = ( fx(i+1,2) - fx(i,2)   ) / Dx0
      !Df1 = ( fx(i+2,2) - fx(i+1,2) ) / Dx1
      Df2 = ( fx(i+3,2) - fx(i+2,2) ) / Dx2
      tDf2 = ( fx_tl(i+3,2) - fx_tl(i+2,2) ) / Dx2 - &
               tDx2 * ( fx(i+3,2) - fx(i+2,2) ) / Dx2**2

      !DDf1 = ( Df1 - Df0 ) / DDx1
      DDf2 = ( Df2 - Df1 ) / DDx2
      tDDf2 = ( tDf2 - tDf1 ) / DDx2 -          &
                tDDx2 * ( Df2 - Df1 ) / DDx2**2

      DDf3 = ( DDf2 - DDf1 ) / DDx3
      tDDf3 = ( tDDf2 - tDDf1 ) / DDx3 -            &
                tDDx3 * ( DDf2 - DDf1 ) / DDx3**2

      ! Integral
      a = 0.50_wp * ( DDf1 + DDf2 )
      ta = 0.50_wp * ( tDDf1 + tDDf2 )

      b = 0.50_wp * ( Dx0 - Dx2 ) * DDf3
      tb = 0.50_wp * ( ( tDx0 - tDx2 ) * DDf3 +               &
                       tDDf3 * ( Dx0 - Dx2 )   )

      sum = sum + Dx1 * ( 0.50_wp * ( fx(i+1,2) + fx(i+2,2) ) -       &
                          Dx1**2 * (a + b) / 6.0_wp )
      tsum = tsum + tDx1 * ( 0.50_wp * (( fx(i+1,2) + fx(i+2,2) ) -       &
                             Dx1**2 * (a + b)) ) +                       &
                    0.50_wp*Dx1*(fx_tl(i+1,2)+fx_tl(i+2,2))       -       &
                    (tb + ta) * Dx1**3 / 6.0_wp

      if (i .eq. k-3) then
         ! Compute  integral of the last interval using
         ! the same differences as above:

         ! Integral
         a = (Dx2 + 2.0_wp*Dx1) * DDf3 / 12.0_wp
         ta = ( (tDx2 + 2.0_wp*tDx1) * DDf3 +                &
                (Dx2 + 2.0_wp*Dx1) * tDDf3   )  / 12.0_wp
         b = DDf2 / 6.0_wp
         tb = tDDf2 / 6.0_wp
         sum = sum + Dx2 * ( fx(n,2) - 0.50_wp*Dx2*Df2 -    &
                             (a + b) * Dx2**2           )
         tsum = tsum + Dx2*fx_tl(n,2) -  0.50_wp*tDf2*Dx2**2 -  &
                (ta + tb) * Dx2**3 +                           &
                ( fx(n,2) - Dx2*Df2 - 3.0_wp*(a + b) * Dx2**2 ) * tDx2

         exit

      end if

      ! Reuse differences from last loop:
      Dx0  = Dx1
      tDx0  = tDx1

      Dx1  = Dx2
      tDx1  = tDx2

      DDx1 = DDx2
      tDDx1 = tDDx2

      Df0  = Df1
      tDf0  = tDf1

      Df1  = Df2
      tDf1  = tDf2

      DDf1 = DDf2
      tDDf1 = tDDf2

   end do

   IntegPolyCube_tl = tsum

   if (IntegPolyCube_tl /= IntegPolyCube_tl) then
      write(*,*) 'IntegPolyCube_tl> IntegPolyCube_tl = NaN', char(10),  &
           'max fx_tl = ', maxval(fx_tl)
   end if

end if

end function IntegPolyCube_tl


subroutine IntegPolyCube_ad ( fx, fx_ad, IntegPolyCubeAdj )

! List of calling arguments:
real (wp), intent(inout)                :: IntegPolyCubeAdj
real (wp), dimension (:,:), intent(in)  :: fx
real (wp), dimension (:,:), intent(out) :: fx_ad

!  List of local variables:
real (wp) :: a, b
real (wp) :: Dx0, Dx1, Dx2, DDx1, DDx2, DDx3
real (wp) :: Df0, Df1, Df2, DDf1, DDf2, DDf3
real (wp) :: asum, aa, ab
real (wp) :: aDx0, aDx1, aDx2, aDDx1, aDDx2, aDDx3
real (wp) :: aDf0, aDf1, aDf2, aDDf1, aDDf2, aDDf3

integer              :: i, j, k,  n
!---------------------------------------------------------------------

j = lbound(fx,1)
k = ubound(fx,1)
n = k - j + 1

if (n .lt. 4) then
   ! At least 4 points are required to compute the integral
   write(*,*) 'IntegPolyCube> More than three points are required, n = ', n
   fx_ad = 0.0_wp
else

   ! Adjoint code
   fx_ad = 0.0_wp
   ! Delta x ( x => fx(i,1) )
   aDx0 = 0.0_wp
   aDx1 = 0.0_wp
   aDx2 = 0.0_wp
   aDDx1 = 0.0_wp
   aDDx2 = 0.0_wp
   aDDx3 = 0.0_wp
   ! Delta y / Delta x  ( y => fx(i,2) )
   aDf0 = 0.0_wp
   aDf1 = 0.0_wp
   aDf2 = 0.0_wp
   aDDf1 = 0.0_wp
   aDDf2 = 0.0_wp
   aDDf3 = 0.0_wp

   ! IntegPolyCube = sum
   asum = IntegPolyCubeAdj

   do i=k-3, j, -1

      ! Delta x ( x => fx(i,1) )
      Dx0 = fx(i+1,1) - fx(i,1)
      Dx1 = fx(i+2,1) - fx(i+1,1)
      Dx2 = fx(i+3,1) - fx(i+2,1)
      DDx1 = Dx0 + Dx1
      DDx2 = Dx1 + Dx2
      DDx3 = DDx2 + Dx0
      ! Delta y / Delta x  ( y => fx(i,2) )
      Df0 = ( fx(i+1,2) - fx(i,2)   ) / Dx0
      Df1 = ( fx(i+2,2) - fx(i+1,2) ) / Dx1
      Df2 = ( fx(i+3,2) - fx(i+2,2) ) / Dx2
      DDf1 = ( Df1 - Df0 ) / DDx1
      DDf2 = ( Df2 - Df1 ) / DDx2
      DDf3 = ( DDf2 - DDf1 ) / DDx3

      if (i .eq. k-3) then
         ! Compute  integral of the last interval using
         ! the same differences as above:

         ! Integral
         a = (Dx2 + 2.0_wp*Dx1) * DDf3 / 12.0_wp
         b = DDf2 / 6.0_wp

         ! sum = sum + Dx2 * ( fx(n,2) - 0.50D0*Dx2*Df2 -    &
         !                    (a + b) * Dx2**2           )
         aDx2 = aDx2 + ( fx(n,2)-Dx2*Df2-3.0_wp * (a+b)*Dx2**2 ) * asum
         fx_ad(n,2) = fx_ad(n,2) + Dx2 * asum
         aDf2 = aDf2 - asum * 0.50_wp * Dx2**2
         aa = -asum * Dx2**3
         ab = -asum * Dx2**3

         ! b = DDf2 / 6.0D0
         aDDf2 = aDDf2 + ab / 6.0_wp

         ! a = (Dx2 + 2.0D0*Dx1) * DDf3 / 12.0D0
         aDx2 =  aDx2 + (DDf3 / 12.0_wp) * aa
         aDx1 = aDx1 + (DDf3 / 6.0_wp) * aa
         aDDf3 = aDDf3 + ((Dx2 + 2.0_wp*Dx1)/12.0_wp) * aa

      end if

      ! Compute integral using 4 point formula
      a = 0.50_wp * ( DDf1 + DDf2 )
      b = 0.50_wp * ( Dx0 - Dx2 ) * DDf3

      ! Integral
      ! sum = sum + Dx1 * ( 0.50D0 * ( fx(i+1,2) + fx(i+2,2) ) -       &
      !                     Dx1**2 * (a + b) / 6.0D0 )
      aDx1 = aDx1 + 0.50_wp*( (fx(i+1,2)+fx(i+2,2)) - Dx1**2 * (a + b) ) * asum
      fx_ad(i+1,2) = fx_ad(i+1,2) + 0.50_wp*Dx1 * asum
      fx_ad(i+2,2) = fx_ad(i+2,2) + 0.50_wp*Dx1 * asum
      aa = -asum * Dx1**3  / 6.0_wp
      ab = -asum * Dx1**3  / 6.0_wp

      ! b = 0.50D0 * ( Dx0 - Dx2 ) * DDf3
      aDx0 = aDx0 + 0.50_wp * DDf3 * ab
      aDx2 = aDx2 - 0.50_wp * DDf3 * ab
      aDDf3 = aDDf3 + 0.50_wp * ( Dx0 - Dx2 ) * ab

      ! a = 0.50D0 * ( DDf1 + DDf2 )
      aDDf1 = aDDf1 + 0.50_wp * aa
      aDDf2 = aDDf2 + 0.50_wp * aa

      ! DDf3 = ( DDf2 - DDf1 ) / DDx3
      aDDf2 = aDDf2 + aDDf3 / DDx3
      aDDf1 = aDDf1 - aDDf3 / DDx3
      aDDx3 = aDDx3 - aDDf3 * ( DDf2 - DDf1 ) / DDx3**2

      ! DDf2 = ( Df2 - Df1 ) / DDx2
      aDf2 = aDf2 + aDDf2 / DDx2
      aDf1 = aDf1 - aDDf2 / DDx2
      aDDx2 = aDDx2 - aDDf2*( Df2 - Df1 ) / DDx2**2

      ! DDf1 = ( Df1 - Df0 ) / DDx1
      aDf1 = aDf1 + aDDf1 / DDx1
      aDf0 = aDf0 - aDDf1 / DDx1
      aDDx1 = aDDx1 - aDDf1*( Df1 - Df0 ) / DDx1**2

      ! Df2 = ( fx(i+3,2) - fx(i+2,2) ) / Dx2
      fx_ad(i+3,2) = fx_ad(i+3,2) + aDf2 / Dx2
      fx_ad(i+2,2) = fx_ad(i+2,2) - aDf2 / Dx2
      aDx2 = aDx2 - aDf2*( fx(i+3,2) - fx(i+2,2) ) / Dx2**2

      ! Df1 = ( fx(i+2,2) - fx(i+1,2) ) / Dx1
      fx_ad(i+2,2) = fx_ad(i+2,2) + aDf1 / Dx1
      fx_ad(i+1,2) = fx_ad(i+1,2) - aDf1 / Dx1
      aDx1 = aDx1 - aDf1*( fx(i+2,2) - fx(i+1,2) ) / Dx1**2

      ! Df0 = ( fx(i+1,2) - fx(i,2)   ) / Dx0
      fx_ad(i+1,2) = fx_ad(i+1,2) + aDf0 / Dx0
      fx_ad(i,2) = fx_ad(i,2) - aDf0 / Dx0
      aDx0 = aDx0 - aDf0*( fx(i+1,2) - fx(i,2)   ) / Dx0**2

      ! DDx3 = DDx2 + Dx0
      aDDx2 = aDDx2 + aDDx3
      aDx0 = aDx0 + aDDx3

      ! DDx2 = Dx1 + Dx2
      aDx1 = aDx1 + aDDx2
      aDx2 = aDx2  + aDDx2

      ! DDx1 = Dx0 + Dx1
      aDx0 = aDx0 + aDDx1
      aDx1 = aDx1 + aDDx1

      ! Dx2 = fx(i+3,1) - fx(i+2,1)
      fx_ad(i+3,1) = fx_ad(i+3,1) + aDx2
      fx_ad(i+2,1) = fx_ad(i+2,1) - aDx2

      ! Dx1 = fx(i+2,1) - fx(i+1,1)
      fx_ad(i+2,1) = fx_ad(i+2,1) + aDx1
      fx_ad(i+1,1) = fx_ad(i+1,1) - aDx1

      ! Dx0 = fx(i+1,1) - fx(i,1)
      fx_ad(i+1,1) = fx_ad(i+1,1) + aDx0
      fx_ad(i,1) = fx_ad(i,1) - aDx0

      if (i .eq. j) then
         ! Compute first set of differences and the integral over the
         ! first interval:

         ! Delta x ( x => fx(i,1) )
         aDx0 = 0.0_wp
         aDx1 = 0.0_wp
         aDx2 = 0.0_wp
         aDDx1 = 0.0_wp
         aDDx2 = 0.0_wp
         aDDx3 = 0.0_wp
         ! Delta y / Delta x  ( y => fx(i,2) )
         aDf0 = 0.0_wp
         aDf1 = 0.0_wp
         aDf2 = 0.0_wp
         aDDf1 = 0.0_wp
         aDDf2 = 0.0_wp
         aDDf3 = 0.0_wp

         a = (Dx0 + 2.0_wp*Dx1) * DDf3 / 12.0_wp
         b = DDf1 / 6.0_wp

         ! sum = Dx0 * ( fx(i,2) + 0.50D0*Dx0*Df0 +    &
         !               (a - b) * Dx0**2           )
         aDx0 = aDx0 + asum * ( fx(i,2) + Dx0*Df0 +       &
                                3.0D0*(a - b) * Dx0**2 )
         fx_ad(i,2) = fx_ad(i,2) + Dx0 * asum
         aDf0 = aDf0 + 0.50_wp * asum * Dx0**2
         aa = asum * Dx0**3
         ab = -asum * Dx0**3

         ! b = DDf1 / 6.0D0
         aDDf1 = aDDf1 + ab / 6.0_wp

         ! a = (Dx0 + 2.0D0*Dx1) * DDf3 / 12.0D0
         aDx0 = aDx0 + aa * DDf3 / 12.0_wp
         aDx1 = aDx1 + aa * DDf3 / 6.0_wp
         aDDf3 = aDDf3 + aa * (Dx0 + 2.0_wp*Dx1) / 12.0_wp

         ! DDf3 = ( DDf2 - DDf1 ) / DDx3
         aDDf2 = aDDf2 + aDDf3 / DDx3
         aDDf1 = aDDf1 - aDDf3 / DDx3
         aDDx3 = aDDx3 - aDDf3 * ( DDf2 - DDf1 ) / DDx3**2

         ! DDf2 = ( Df2 - Df1 ) / DDx2
         aDf2 = aDf2 + aDDf2 / DDx2
         aDf1 = aDf1 - aDDf2 / DDx2
         aDDx2 = aDDx2 - aDDf2*( Df2 - Df1 ) / DDx2**2

         ! DDf1 = ( Df1 - Df0 ) / DDx1
         aDf1 = aDf1 + aDDf1 / DDx1
         aDf0 = aDf0 - aDDf1 / DDx1
         aDDx1 = aDDx1 - aDDf1*( Df1 - Df0 ) / DDx1**2

         ! Df2 = ( fx(i+3,2) - fx(i+2,2) ) / Dx2
         fx_ad(i+3,2) = fx_ad(i+3,2) + aDf2 / Dx2
         fx_ad(i+2,2) = fx_ad(i+2,2) - aDf2 / Dx2
         aDx2 = aDx2 - aDf2*( fx(i+3,2) - fx(i+2,2) ) / Dx2**2

         ! Df1 = ( fx(i+2,2) - fx(i+1,2) ) / Dx1
         fx_ad(i+2,2) = fx_ad(i+2,2) + aDf1 / Dx1
         fx_ad(i+1,2) = fx_ad(i+1,2) - aDf1 / Dx1
         aDx1 = aDx1 - aDf1*( fx(i+2,2) - fx(i+1,2) ) / Dx1**2

         ! Df0 = ( fx(i+1,2) - fx(i,2)   ) / Dx0
         fx_ad(i+1,2) = fx_ad(i+1,2) + aDf0 / Dx0
         fx_ad(i,2) = fx_ad(i,2) - aDf0 / Dx0
         aDx0 = aDx0 - aDf0*( fx(i+1,2) - fx(i,2)   ) / Dx0**2

         ! DDx3 = DDx2 + Dx0
         aDDx2 = aDDx2 + aDDx3
         aDx0 = aDx0 + aDDx3

         ! DDx2 = Dx1 + Dx2
         aDx1 = aDx1 + aDDx2
         aDx2 = aDx2  + aDDx2

         ! DDx1 = Dx0 + Dx1
         aDx0 = aDx0 + aDDx1
         aDx1 = aDx1 + aDDx1

         ! Dx2 = fx(i+3,1) - fx(i+2,1)
         fx_ad(i+3,1) = fx_ad(i+3,1) + aDx2
         fx_ad(i+2,1) = fx_ad(i+2,1) - aDx2

         ! Dx1 = fx(i+2,1) - fx(i+1,1)
         fx_ad(i+2,1) = fx_ad(i+2,1) + aDx1
         fx_ad(i+1,1) = fx_ad(i+1,1) - aDx1

         ! Dx0 = fx(i+1,1) - fx(i,1)
         fx_ad(i+1,1) = fx_ad(i+1,1) + aDx0
         fx_ad(i,1) = fx_ad(i,1) - aDx0

      end if

      ! Delta x ( x => fx(i,1) )
      aDx0 = 0.0_wp
      aDx1 = 0.0_wp
      aDx2 = 0.0_wp
      aDDx1 = 0.0_wp
      aDDx2 = 0.0_wp
      aDDx3 = 0.0_wp
      ! Delta y / Delta x  ( y => fx(i,2) )
      aDf0 = 0.0_wp
      aDf1 = 0.0_wp
      aDf2 = 0.0_wp
      aDDf1 = 0.0_wp
      aDDf2 = 0.0_wp
      aDDf3 = 0.0_wp

   end do

end if

IntegPolyCubeAdj = 0.0_wp

end subroutine IntegPolyCube_ad


subroutine CrossLineEllips_tl ( RefEllips, Alt, SLine, CrossPt, lambda, exist, &
                                Alt_tl, SLine_tl, CrossPt_tl, lambda_tl )

! List of calling arguments:
type (EllipsParam), intent(in)                    :: RefEllips
real (wp), intent(in)                  :: Alt
type (Line), intent(in)                           :: SLine
real (wp), dimension(1:3), intent(out) :: CrossPt
real (wp), intent(out)                 :: lambda
logical, optional,  intent(out)                   :: exist

real (wp), intent(in)                  :: Alt_tl
type (Line), intent(in)                           :: SLine_tl
real (wp), dimension(1:3), intent(out) :: CrossPt_tl
real (wp), intent(out)                 :: lambda_tl

! List of local variables:
real (wp) :: ha, hb
real (wp) :: a0, a1, a2
real (wp) :: lambda1, lambda2
logical              :: IsReal
real (wp) :: Dha, Dhb
real (wp) :: Da0, Da1, Da2
real (wp) :: Dlambda1, Dlambda2
real (wp) :: Dp, Dq
!---------------------------------------------------------------------

! Hier wird nicht der Schnittpunkt der Geraden mit dem Referenz-Ellipsoid,
! sondern mit einem Ellipsoid in der Hoehe "Alt" gesucht. Die entsprechenden
! Halbachsen sind:
ha = RefEllips%a + Alt
hb = RefEllips%b + Alt

Dha = Alt_tl
Dhb = Alt_tl

! Die beiden Schnittpunkte der Geraden mit der Ellipse ergeben sich als Loesung
! einer quadratischen Gleichung: a0 + a1*lambda + a2*lambda**2 = 0
a0 = (SLine%start(1)**2/ha**2) + (SLine%start(2)**2/ha**2) + &
     (SLine%start(3)**2/hb**2) - 1.0
a1 = 2.0_wp*( (SLine%start(1)*SLine%unitvec(1)/ha**2) +  &
             (SLine%start(2)*SLine%unitvec(2)/ha**2) +  &
             (SLine%start(3)*SLine%unitvec(3)/hb**2)   )
a2 = (SLine%unitvec(1)/ha)**2 +  &
     (SLine%unitvec(2)/ha)**2 +  &
     (SLine%unitvec(3)/hb)**2

Da0 = 2.0_wp * (SLine%start(1)/ha**2) * SLine_tl%start(1) - &
      2.0_wp * (SLine%start(1)**2/ha**3) * Dha            + &
      2.0_wp * (SLine%start(2)/ha**2) * SLine_tl%start(2) - &
      2.0_wp * (SLine%start(2)**2/ha**3) * Dha            + &
      2.0_wp * (SLine%start(3)/hb**2) * SLine_tl%start(3) - &
      2.0_wp * (SLine%start(3)**2/hb**3) * Dhb

Da1 = (SLine%unitvec(1)/ha**2) * SLine_tl%start(1)          +   &
      (SLine%start(1)/ha**2) * SLine_tl%unitvec(1)          -   &
      2.0_wp * (SLine%start(1)*SLine%unitvec(1)/ha**3) * Dha +   &
      (SLine%unitvec(2)/ha**2) * SLine_tl%start(2)          +   &
      (SLine%start(2)/ha**2) * SLine_tl%unitvec(2)          -   &
      2.0_wp * (SLine%start(2)*SLine%unitvec(2)/ha**3) * Dha +   &
      (SLine%unitvec(3)/hb**2) * SLine_tl%start(3)          +   &
      (SLine%start(3)/hb**2) * SLine_tl%unitvec(3)          -   &
      2.0_wp * (SLine%start(3)*SLine%unitvec(3)/hb**3) * Dhb
Da1 = 2.0_wp * Da1

Da2 = 2.0_wp * (SLine%unitvec(1)/ha**2) * SLine_tl%unitvec(1) - &
      2.0_wp * (SLine%unitvec(1)**2/ha**3) * Dha              + &
      2.0_wp * (SLine%unitvec(2)/ha**2) * SLine_tl%unitvec(2) - &
      2.0_wp * (SLine%unitvec(2)**2/ha**3) * Dha              + &
      2.0_wp * (SLine%unitvec(3)/hb**2) * SLine_tl%unitvec(3) - &
      2.0_wp * (SLine%unitvec(3)**2/hb**3) * Dhb

! p = a1/a2
Dp = Da1/a2 - (a1 * Da2)/a2**2
! q = a0/a2
Dq =  Da0/a2 - (a0 * Da2)/a2**2

Dlambda1 = 0.0_wp
Dlambda2 = 0.0_wp
!call QuadEq(a1/a2, a0/a2, lambda1, lambda2, IsReal)
call QuadEq_tl(a1/a2, a0/a2, lambda1, lambda2, IsReal,    &
               Dp, Dq, Dlambda1, Dlambda2 )

if (IsReal) then
   ! Es gibt Schnittpunkte
   if (lambda1 .ge. lambda2) then
      ! Den naeher gelegenen Punkt auswaehlen
      CrossPt_tl = 0.0_wp
      call LinePos_tl(SLine, lambda1, CrossPt,           &
                      SLine_tl, Dlambda1, CrossPt_tl)
      lambda = lambda1
      lambda_tl = Dlambda1
   else
      CrossPt_tl = 0.0_wp
      call LinePos_tl(SLine, lambda2, CrossPt,           &
                      SLine_tl, Dlambda2, CrossPt_tl)
      lambda = lambda2
      lambda_tl = Dlambda2
   end if
else
   ! Kein Schnittpunkt: Die Gerade verlaeuft vollstaendig ausserhalb des
   ! Ellipsoids.
   write(*,*) 'CrossLineEllips> Warning: No intersecting point found'
   CrossPt_tl = 0.0_wp
   lambda_tl  = 0.0_wp
end if

if (present(exist)) then
   exist = IsReal
end if

!write(*,*) IsReal, lambda1, lambda2

end subroutine CrossLineEllips_tl


subroutine CrossLineEllips_ad( RefEllips, Alt, SLine, CrossPt, lambda, exist, &
                               Alt_ad, SLine_ad, CrossPt_ad, lambda_ad )

! List of calling arguments:
type (EllipsParam), intent(in)                    :: RefEllips
real (wp), intent(in)                  :: Alt
type (Line), intent(in)                           :: SLine
real (wp), dimension(1:3), intent(out) :: CrossPt
real (wp), intent(out)                 :: lambda
logical, optional,  intent(out)                   :: exist

real (wp), intent(out)                   :: Alt_ad
type (Line), intent(out)                            :: SLine_ad
real (wp), dimension(1:3), intent(inout) :: CrossPt_ad
real (wp), intent(inout)                 :: lambda_ad

! List of local variables:
real (wp) :: ha, hb
real (wp) :: a0, a1, a2
real (wp) :: lambda1, lambda2
logical              :: IsReal

real (wp) :: Dlambda1, Dlambda2
real (wp) :: Dp, Dq
real (wp) :: Da0, Da1, Da2
real (wp) :: Dha, Dhb
!---------------------------------------------------------------------

! Hier wird nicht der Schnittpunkt der Geraden mit dem Referenz-Ellipsoid,
! sondern mit einem Ellipsoid in der Hoehe "Alt" gesucht. Die entsprechenden
! Halbachsen sind:
ha = RefEllips%a + Alt
hb = RefEllips%b + Alt

! Die beiden Schnittpunkte der Geraden mit der Ellipse ergeben sich als Loesung
! einer quadratischen Gleichung: a0 + a1*lambda + a2*lambda**2 = 0
a0 = (SLine%start(1)**2/ha**2) + (SLine%start(2)**2/ha**2) + &
     (SLine%start(3)**2/hb**2) - 1.0
a1 = 2.0_wp*( (SLine%start(1)*SLine%unitvec(1)/ha**2) +  &
           (SLine%start(2)*SLine%unitvec(2)/ha**2) +  &
           (SLine%start(3)*SLine%unitvec(3)/hb**2)   )
a2 = (SLine%unitvec(1)/ha)**2 +  &
     (SLine%unitvec(2)/ha)**2 +  &
     (SLine%unitvec(3)/hb)**2

call QuadEq(a1/a2, a0/a2, lambda1, lambda2, IsReal)

if (IsReal) then
   ! Es gibt Schnittpunkte
   if (lambda1 .ge. lambda2) then
      ! Den naeher gelegenen Punkt auswaehlen
      call LinePos(SLine, lambda1, CrossPt)
      lambda = lambda1
   else
      call LinePos(SLine, lambda2, CrossPt)
      lambda = lambda2
   end if
else
   ! Kein Schnittpunkt: Die Gerade verlaeuft vollstaendig ausserhalb des
   ! Ellipsoids.
   write(*,*) 'CrossLineEllips> Warning: No intersecting point found'
end if

if (present(exist)) then
   exist = IsReal
end if

!write(*,*) IsReal, lambda1, lambda2

! LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL

if (IsReal) then
   ! Es gibt Schnittpunkte
   if (lambda1 .ge. lambda2) then
      ! Den naeher gelegenen Punkt auswaehlen
      Dlambda1 = lambda_ad
      call LinePos_ad(SLine, lambda1, CrossPt,           &
                      SLine_ad, Dlambda1, CrossPt_ad)
   Dlambda2 = 0.0_wp
   else
      Dlambda2 = lambda_ad
      call LinePos_ad(SLine, lambda2, CrossPt,           &
                      SLine_ad, Dlambda2, CrossPt_ad)
   Dlambda1 = 0.0_wp
   end if
else
   ! Kein Schnittpunkt: Die Gerade verlaeuft vollstaendig ausserhalb des
   ! Ellipsoids.
   Dlambda1 = 0.0_wp
   Dlambda2 = 0.0_wp
   SLine_ad%start   = 0.0_wp
   SLine_ad%unitvec = 0.0_wp
end if

! call QuadEq(a1/a2, a0/a2, lambda1, lambda2, IsReal)
call QuadEq_ad(a1/a2, a0/a2, IsReal, Dp, Dq, Dlambda1, Dlambda2)

! p = a1/a2
Da1 = Dp/a2
Da2 = - (a1 * Dp)/a2**2
! q = a0/a2
Da0 = Dq/a2
Da2 = Da2 - (a0 * Dq)/a2**2

! a2=(SLine%unitvec(1)/ha)**2+(SLine%unitvec(2)/ha)**2+(SLine%unitvec(3)/hb)**2
SLine_ad%unitvec(1) = SLine_ad%unitvec(1) + (2.0_wp*SLine%unitvec(1)/ha**2)*Da2
SLine_ad%unitvec(2) = SLine_ad%unitvec(2) + (2.0_wp*SLine%unitvec(2)/ha**2)*Da2
SLine_ad%unitvec(3) = SLine_ad%unitvec(3) + (2.0_wp*SLine%unitvec(3)/hb**2)*Da2
Dha = -2.0_wp * Da2 * SLine%unitvec(1)**2/ha**3
Dha = Dha - 2.0_wp * Da2 * SLine%unitvec(2)**2/ha**3
Dhb = -2.0_wp * Da2 * SLine%unitvec(3)**2/hb**3

!a1 = 2.0D0*( (SLine%start(1)*SLine%unitvec(1)/ha**2) +  &
!           (SLine%start(2)*SLine%unitvec(2)/ha**2) +  &
!           (SLine%start(3)*SLine%unitvec(3)/hb**2)   )
SLine_ad%start(1) = SLine_ad%start(1) +  2.0_wp*(SLine%unitvec(1)/ha**2) * Da1
SLine_ad%start(2) = SLine_ad%start(2) +  2.0_wp*(SLine%unitvec(2)/ha**2) * Da1
SLine_ad%start(3) = SLine_ad%start(3) +  2.0_wp*(SLine%unitvec(3)/hb**2) * Da1
Dha = Dha - 4.0_wp*(SLine%start(1)*SLine%unitvec(1)/ha**3) * Da1
Dha = Dha - 4.0_wp*(SLine%start(2)*SLine%unitvec(2)/ha**3) * Da1
Dhb = Dhb - 4.0_wp*(SLine%start(3)*SLine%unitvec(3)/hb**3) * Da1
SLine_ad%unitvec(1) = SLine_ad%unitvec(1) + 2.0_wp*(SLine%start(1)/ha**2) * Da1
SLine_ad%unitvec(2) = SLine_ad%unitvec(2) + 2.0_wp*(SLine%start(2)/ha**2) * Da1
SLine_ad%unitvec(3) = SLine_ad%unitvec(3) + 2.0_wp*(SLine%start(3)/hb**2) * Da1

!a0 = (SLine%start(1)**2/ha**2) + (SLine%start(2)**2/ha**2) + &
!     (SLine%start(3)**2/hb**2) - 1.0
SLine_ad%start(1) = SLine_ad%start(1) + 2.0_wp * (SLine%start(1)/ha**2) * Da0
SLine_ad%start(2) = SLine_ad%start(2) + 2.0_wp * (SLine%start(2)/ha**2) * Da0
SLine_ad%start(3) = SLine_ad%start(3) + 2.0_wp * (SLine%start(3)/hb**2) * Da0
Dha = Dha - 2.0_wp * ((SLine%start(1)**2)/ha**3) * Da0
Dha = Dha - 2.0_wp * ((SLine%start(2)**2)/ha**3) * Da0
Dhb = Dhb - 2.0_wp * ((SLine%start(3)**2)/hb**3) * Da0

! hb = RefEllips%b + Alt
Alt_ad = Dhb

! ha = RefEllips%a + Alt
Alt_ad = Alt_ad + Dha

lambda_ad  = 0.0_wp
CrossPt_ad = 0.0_wp

end subroutine CrossLineEllips_ad


function ExpInt1D_tl(RefPt, inode, RefPt_tl, inode_tl)

! List of calling arguments:
real (wp)                                 :: ExpInt1D_tl
real (wp), intent(in)                     :: RefPt
real (wp), dimension(1:2,1:2), intent(in) :: inode
real (wp), intent(in)                     :: RefPt_tl
real (wp), dimension(1:2,1:2), intent(in) :: inode_tl

! List of local variables:
!real (wp) :: x, y
real (wp), dimension(1:2,1:2) :: n, n_tl
real (wp) :: E, L, B, NE
!---------------------------------------------------------------------

if (abs(inode(1,2)) .lt. 1E-10_wp .and. abs(inode(2,2)) .lt. 1E-10_wp ) then
   ! N1 = N2 = 0  => N(z) = 0

   ExpInt1D_tl = 0.0_wp

else if (abs(inode(2,1)-inode(1,1)) < 1.0E-3_wp) then
   ! identical coordinates, interpolation not possible, return first value
   ExpInt1D_tl = inode_tl(1,2)

else
   ! interpolate

   n = inode
   n_tl = inode_tl

   if (abs(n(1,2)) .lt. 1E-10_wp) then
      ! N1 = 0, N2 <> 0
      n(1,2) = n(2,2)/100.0_wp
      n_tl(1,2) = inode_tl(1,2)/100.0_wp
   end if
   if (abs(n(2,2)) .lt. 1E-10_wp) then
      ! N2 = 0, N1 <> 0
      n(2,2) = n(1,2)/100.0_wp
      n_tl(2,2) = inode_tl(1,2)/100.0_wp
   end if


   !call SaveDiff('tl', n_tl(1,1), n_tl(1,2), n_tl(2,1),  &
   !                    n_tl(2,2), RefPt_tl)

   ! ExpInt1D = node(2,2) * (node(1,2)/node(2,2))**                     &
   !                        ((node(2,1)-RefPt)/(node(2,1)-node(1,1)))
   !          = n(2,2) * exp( ( (n(2,1)-RefPt)/(n(2,1)-n(1,1)) ) *      &
   !                            log(n(1,2)/n(2,2))                  )

   B = (n(2,1)-RefPt) / (n(2,1)-n(1,1))
   L = log( n(1,2)/n(2,2) )
   E = exp( B * L )
   NE = n(2,2) * E

   ! d I / d RefPt
   ExpInt1D_tl = -L/(n(2,1)-n(1,1)) * NE * RefPt_tl

   ! d I / d n(1,1)
   ExpInt1D_tl = ExpInt1D_tl + L * NE * (n(2,1)-RefPt)/(n(2,1)-n(1,1))**2  &
                                 * n_tl(1,1)

   ! d I / d n(1,2)
   ExpInt1D_tl = ExpInt1D_tl + B * NE / n(1,2) * n_tl(1,2)

   ! d I / d n(2,1)
   ExpInt1D_tl = ExpInt1D_tl + L * NE * (RefPt-n(1,1))/(n(2,1)-n(1,1))**2  &
                                 * n_tl(2,1)

   ! d I / d n(2,2)
   ExpInt1D_tl = ExpInt1D_tl + n_tl(2,2) *    &
                 (E - NE * B * (n(1,2)/n(2,2)**2) * (n(2,2)/n(1,2)))

end if

if (ExpInt1D_tl /= ExpInt1D_tl) then
   ! ExpInt1D_tl = NaN
   write(*,*) ' ExpInt1D_tl = NaN'

end if

end function ExpInt1D_tl

! For repeated calls the adjoint varaibles must be "inout" for summing up
! all increments !!!
subroutine ExpInt1D_ad(RefPt, inode, ExpInt_ad, RefPt_ad, inode_ad)

! List of calling arguments:
real (wp), intent(in)                      :: RefPt
real (wp), dimension(1:2,1:2), intent(in)  :: inode
real (wp), intent(inout)                   :: ExpInt_ad   ! input
real (wp), intent(out)                     :: RefPt_ad    ! output
real (wp), dimension(1:2,1:2), intent(inout) :: inode_ad  ! output


! List of local variables:
!real (wp) :: x, y
real (wp), dimension(1:2,1:2) :: n
real (wp) :: E, L, B, NE
!---------------------------------------------------------------------

if (abs(inode(1,2)) .lt. 1E-10_wp .and. abs(inode(2,2)) .lt. 1E-10_wp ) then
   ! N1 = N2 = 0  => N(z) = 0

   ExpInt_ad = 0.0_wp
   RefPt_ad  = 0.0_wp
   !inode_ad  = 0.0_wp

else if (abs(inode(2,1)-inode(1,1)) < 1.0E-3_wp) then
   ! identical coordinates, interpolation not possible, return first value
   ExpInt_ad =  inode_ad(1,2)
   RefPt_ad  = 0.0_wp

else
   ! interpolate

   n = inode

   B = (n(2,1)-RefPt) / (n(2,1)-n(1,1))
   L = log( n(1,2)/n(2,2) )
   E = exp( B * L )
   NE = n(2,2) * E

   RefPt_ad = -L/(n(2,1)-n(1,1)) * NE * ExpInt_ad

   inode_ad(1,1) = inode_ad(1,1) + &
                   L * NE * (n(2,1)-RefPt)/(n(2,1)-n(1,1))**2 * ExpInt_ad
   inode_ad(1,2) = inode_ad(1,2) + &
                   B * NE / n(1,2) * ExpInt_ad
   inode_ad(2,1) = inode_ad(2,1) + &
                   L * NE * (RefPt-n(1,1))/(n(2,1)-n(1,1))**2 * ExpInt_ad
   inode_ad(2,2) = inode_ad(2,2) +                                          &
                   (E - NE * B * (n(1,2)/n(2,2)**2) * (n(2,2)/n(1,2))) *    &
                   ExpInt_ad

   !call SaveDiff('ad', inode_ad(1,1), inode_ad(1,2), inode_ad(2,1),  &
   !                    inode_ad(2,2), RefPt_ad)

   if (abs(inode(1,2)) .lt. 1E-10_wp) then
      ! N1 = 0, N2 <> 0
      ! node(1,2) = node(2,2)/100.0D0
      inode_ad(2,2) = inode_ad(2,2) + inode_ad(1,2) / 100.0_wp
   end if
   if (abs(inode(2,2)) .lt. 1.0E-10_wp) then
      ! N2 = 0, N1 <> 0
      ! node(2,2) = node(1,2)/100.0D0
      inode_ad(1,2) = inode_ad(1,2) + inode_ad(2,2) / 100.0_wp
   end if

   ExpInt_ad = 0.0_wp

end if

end subroutine ExpInt1D_ad


function BiLinear2D_tl(RefPt, node, RefPt_tl, node_tl)

! List of calling arguments:
real (wp)                                 :: BiLinear2D_tl
real (wp), dimension(1:2), intent(in)     :: RefPt
real (wp), dimension(1:4,1:3), intent(in) :: node
real (wp), dimension(1:2), intent(in)     :: RefPt_tl
real (wp), dimension(1:4,1:3), intent(in) :: node_tl

! List of local variables:
real (wp) :: c1, c2, d1, d2, d3, d4
real (wp) :: ya, yb, fa, fb, dy      !, a1, a2, a3
real (wp) :: Dc1, Dc2, Dd1, Dd2, Dd3, Dd4
real (wp) :: Dya, Dyb, Dfa, Dfb, Ddy !, Da1, Da2, Da3
real (wp), dimension(1:2) :: Rpt
!---------------------------------------------------------------------

! c1 = x2 - x1
c1 = node(2,1) - node(1,1)
Dc1 = node_tl(2,1) - node_tl(1,1)

! c2 = x4 - x3
c2 = node(4,1) - node(3,1)
Dc2 = node_tl(4,1) - node_tl(3,1)

! d1 = (x - x1)/(x2-x1) =  (x - x1)/c1
d1 = (RefPt(1) - node(1,1))/c1
Dd1 = ((RefPt_tl(1) - node_tl(1,1))/c1) - ((RefPt(1) - node(1,1))/c1**2) * Dc1

! d2 = (x2 - x)/(x2-x1) =  (x2 - x)/c1
d2 = (node(2,1) - RefPt(1))/c1
Dd2 = ((node_tl(2,1) - RefPt_tl(1))/c1) - ((node(2,1) - RefPt(1))/c1**2) * Dc1

! d3 = (x - x3)/(x4-x3) =  (x - x3)/c2
d3 = (RefPt(1) - node(3,1))/c2
Dd3 = ((RefPt_tl(1) - node_tl(3,1))/c2) - ((RefPt(1) - node(3,1))/c2**2) * Dc2

! d4 = (x4 - x)/(x4-x3) =  (x4 - x)/c2
d4 = (node(4,1) - RefPt(1))/c2
Dd4 = ((node_tl(4,1) - RefPt_tl(1))/c2) - ((node(4,1) - RefPt(1))/c2**2) * Dc2

!write(*,*) 'BiLinear2D> c1, c2, d1, d2, d3, d4 : ', c1, c2, d1, d2, d3, d4

fa = d2*node(1,3) + d1*node(2,3)
Dfa = d2*node_tl(1,3) + Dd2*node(1,3) + d1*node_tl(2,3) + Dd1*node(2,3)

fb = d4*node(3,3) + d3*node(4,3)
Dfb = d4*node_tl(3,3) + Dd4*node(3,3) + d3*node_tl(4,3) + Dd3*node(4,3)

ya = d2*node(1,2) + d1*node(2,2)
Dya = d2*node_tl(1,2) + Dd2*node(1,2) + d1*node_tl(2,2) + Dd1*node(2,2)

yb = d4*node(3,2) + d3*node(4,2)
Dyb = d4*node_tl(3,2) + Dd4*node(3,2) + d3*node_tl(4,2) + Dd3*node(4,2)

dy = yb - ya
Ddy = Dyb - Dya
!call SaveDiff('tl', Dfa, Dfb, Dyb, Dya, RefPt_tl(2), Ddy)

if (abs(dy) .lt. 1.0E-9_wp .and. abs(fb-fa) .lt. 1.0E-7_wp) then
   ! x = x2, singularity, should happen only for 3 point interpolation
   write(*,*) 'BiLinear2D_tl> x = x2'

   ! Die Funktion ist hier nicht differenzierbar.
   ! Der Limes der originalen Funktion darf hier nicht einfach
   ! abgeleitet werden, sondern es muss der Limes der Ableitung
   ! berechnet werden !

   ! Notloesung:
   ! "RefPt" etwas verschieben, so dass keine Singularitaet mehr
   ! auftritt, dann alles nochmal rechnen:

   Rpt(1) = RefPt(1) + 1.0E-8_wp
   Rpt(2) = RefPt(2) + 1.0E-8_wp
   write(*,*) 'RefPt = ', RefPt
   write(*,*) 'RPt   = ', RPt

   ! Kopie von oben, RefPt => Rpt
   ! c1 = x2 - x1
   c1 = node(2,1) - node(1,1)
   Dc1 = node_tl(2,1) - node_tl(1,1)

   ! c2 = x4 - x3
   c2 = node(4,1) - node(3,1)
   Dc2 = node_tl(4,1) - node_tl(3,1)

   ! d1 = (x - x1)/(x2-x1) =  (x - x1)/c1
   d1 = (Rpt(1) - node(1,1))/c1
   Dd1 = (RefPt_tl(1) - node_tl(1,1))/c1 - (Rpt(1) - node(1,1))/c1**2 * Dc1

   ! d2 = (x2 - x)/(x2-x1) =  (x2 - x)/c1
   d2 = (node(2,1) - Rpt(1))/c1
   Dd2 = (node_tl(2,1) - RefPt_tl(1))/c1 - (node(2,1) - Rpt(1))/c1**2 * Dc1

   ! d3 = (x - x3)/(x4-x3) =  (x - x3)/c2
   d3 = (Rpt(1) - node(3,1))/c2
   Dd3 = (RefPt_tl(1) - node_tl(3,1))/c2 - (Rpt(1) - node(3,1))/c2**2 * Dc2

   ! d4 = (x4 - x)/(x4-x3) =  (x4 - x)/c2
   d4 = (node(4,1) - Rpt(1))/c2
   Dd4 = (node_tl(4,1) - RefPt_tl(1))/c2 - (node(4,1) - Rpt(1))/c2**2 * Dc2

   !write(*,*) 'BiLinear2D> c1, c2, d1, d2, d3, d4 : ', c1, c2, d1, d2, d3, d4

   fa = d2*node(1,3) + d1*node(2,3)
   Dfa = d2*node_tl(1,3) + Dd2*node(1,3) + d1*node_tl(2,3) + Dd1*node(2,3)

   fb = d4*node(3,3) + d3*node(4,3)
   Dfb = d4*node_tl(3,3) + Dd4*node(3,3) + d3*node_tl(4,3) + Dd3*node(4,3)

   ya = d2*node(1,2) + d1*node(2,2)
   Dya = d2*node_tl(1,2) + Dd2*node(1,2) + d1*node_tl(2,2) + Dd1*node(2,2)

   yb = d4*node(3,2) + d3*node(4,2)
   Dyb = d4*node_tl(3,2) + Dd4*node(3,2) + d3*node_tl(4,2) + Dd3*node(4,2)

   dy = yb - ya
   Ddy = Dyb - Dya

   BiLinear2D_tl = Dfa*(yb- Rpt(2))/dy + fa*(Dyb- RefPt_tl(2))/dy - &
                   fa*(yb- Rpt(2))/dy**2 * Ddy +                    &
                   Dfb*(Rpt(2)-ya)/dy + fb*(RefPt_tl(2)-Dya)/dy -   &
                   fb*(Rpt(2)-ya)/dy**2 * Ddy


else
   !BiLinear2D = fa*(yb- RefPt(2))/dy + fb*(RefPt(2)-ya)/dy
   BiLinear2D_tl = Dfa*(yb- RefPt(2))/dy + fa*(Dyb- RefPt_tl(2))/dy - &
                   fa*(yb- RefPt(2))/dy**2 * Ddy +                    &
                   Dfb*(RefPt(2)-ya)/dy + fb*(RefPt_tl(2)-Dya)/dy -   &
                   fb*(RefPt(2)-ya)/dy**2 * Ddy
end if


end function BiLinear2D_tl


subroutine BiLinear2D_ad(RefPt, node, RefPt_ad, node_ad, BiLinear2DVal_ad)

! List of calling arguments:
real (wp), intent(inout)                  :: BiLinear2DVal_ad
real (wp), dimension(1:2), intent(in)     :: RefPt
real (wp), dimension(1:4,1:3), intent(in) :: node
real (wp), dimension(1:2), intent(out)     :: RefPt_ad
real (wp), dimension(1:4,1:3), intent(out) :: node_ad

! List of local variables:
real (wp) :: BiLinear2D
real (wp) :: c1, c2, d1, d2, d3, d4
real (wp) :: ya, yb, fa, fb, dy, a1, a2, a3
real (wp) :: Dc1, Dc2, Dd1, Dd2, Dd3, Dd4
real (wp) :: Dya, Dyb, Dfa, Dfb, Ddy !, Da1, Da2, Da3
!---------------------------------------------------------------------

! c1 = x2 - x1
c1 = node(2,1) - node(1,1)
! c2 = x4 - x3
c2 = node(4,1) - node(3,1)
! d1 = (x - x1)/(x2-x1) =  (x - x1)/c1
d1 = (RefPt(1) - node(1,1))/c1
! d2 = (x2 - x)/(x2-x1) =  (x2 - x)/c1
d2 = (node(2,1) - RefPt(1))/c1
! d3 = (x - x3)/(x4-x3) =  (x - x3)/c2
d3 = (RefPt(1) - node(3,1))/c2
! d4 = (x4 - x)/(x4-x3) =  (x4 - x)/c2
d4 = (node(4,1) - RefPt(1))/c2
!write(*,*) 'BiLinear2D> c1, c2, d1, d2, d3, d4 : ', c1, c2, d1, d2, d3, d4

fa = d2*node(1,3) + d1*node(2,3)
fb = d4*node(3,3) + d3*node(4,3)
ya = d2*node(1,2) + d1*node(2,2)
yb = d4*node(3,2) + d3*node(4,2)
dy = yb - ya

if (abs(dy) .lt. 1.0E-9_wp .and. abs(fb-fa) .lt. 1.0E-7_wp) then
   ! x = x2, singularity, should happen only for 3 point interpolation
   !write(*,*) 'BiLinear2D> x = x2'

   a1 = - node(2,1)*node(3,2)*node(2,3)        &
        + node(1,1)*node(3,2)*node(2,3)        &
        - node(2,1)*node(2,2)*node(1,3)        &
        + node(3,1)*node(2,2)*node(1,3)        &
        + node(2,1)*node(1,2)*node(2,3)        &
        - node(3,1)*node(1,2)*node(2,3)        &
        + node(2,1)*node(2,2)*node(3,3)        &
        - node(1,1)*node(2,2)*node(3,3)

   a2 =  RefPt(2) * (   node(2,1)*node(1,3)    &
                      - node(3,1)*node(1,3)    &
                      + node(3,1)*node(2,3)    &
                      - node(2,1)*node(3,3)    &
                      + node(1,1)*node(3,3)    &
                      - node(1,1)*node(2,3)  )

   a3 =   node(2,1)*node(3,2)    &
        - node(1,1)*node(3,2)    &
        + node(1,1)*node(2,2)    &
        - node(2,1)*node(1,2)    &
        + node(3,1)*node(1,2)    &
        - node(3,1)*node(2,2)

   BiLinear2D = - (a1 + a2) / a3

else
   BiLinear2D = fa*(yb- RefPt(2))/dy + fb*(RefPt(2)-ya)/dy
end if

!---------------------------------------------------------------------

! BiLinear2D = fa*(yb- RefPt(2))/dy + fb*(RefPt(2)-ya)/dy
Dfa = ( (yb-RefPt(2))/dy ) * BiLinear2DVal_ad
Dfb = ( (RefPt(2)-ya)/dy ) * BiLinear2DVal_ad
Dyb = ( fa/dy ) * BiLinear2DVal_ad
Dya = -( fb/dy ) * BiLinear2DVal_ad
RefPt_ad(2) = ( (fb-fa)/dy ) * BiLinear2DVal_ad
Ddy =  -( fa*(yb- RefPt(2))/dy**2 + fb*(RefPt(2)-ya)/dy**2 ) * BiLinear2DVal_ad
! call SaveDiff('ad', Dfa, Dfb, Dyb, Dya, RefPt_ad(2), Ddy)

if (abs(dy) .lt. 1.0E-9_wp .and. abs(fb-fa) .lt. 1.0E-7_wp) then
   ! Hier waere eine Sonderbehandlung noetig, da die Funktion an diesem
   ! Punkt nicht stetig differenzierbar ist !!!
end if

! dy = yb - ya
Dyb = Dyb + Ddy
Dya = Dya - Ddy

! yb = d4*node(3,2) + d3*node(4,2)
Dd4 = node(3,2) * Dyb
Dd3 = node(4,2) * Dyb
node_ad(3,2) = d4 * Dyb
node_ad(4,2) = d3 * Dyb

!write(*,*) 'BiLinear2D_ad> BiLinear2DVal_ad = ', BiLinear2DVal_ad
!write(*,*) 'BiLinear2D_ad> Dyb              = ', Dyb


! ya = d2*node(1,2) + d1*node(2,2)
Dd2 = node(1,2) * Dya
Dd1 = node(2,2) * Dya
node_ad(1,2) = d2 * Dya
node_ad(2,2) = d1 * Dya

! fb = d4*node(3,3) + d3*node(4,3)
Dd4 = Dd4 + node(3,3) * Dfb
Dd3 = Dd3 + node(4,3) * Dfb
node_ad(3,3) = d4 * Dfb
node_ad(4,3) = d3 * Dfb

! fa = d2*node(1,3) + d1*node(2,3)
Dd2 = Dd2 + node(1,3) * Dfa
Dd1 = Dd1 + node(2,3) * Dfa
node_ad(1,3) = d2 * Dfa
node_ad(2,3) = d1 * Dfa

! d4 = (node(4,1) - RefPt(1))/c2
node_ad(4,1) = Dd4/c2
RefPt_ad(1) = -Dd4/c2
Dc2 = -( (node(4,1) - RefPt(1))/c2**2 ) * Dd4

! d3 = (RefPt(1) - node(3,1))/c2
RefPt_ad(1) = RefPt_ad(1) + Dd3/c2
node_ad(3,1) = -Dd3/c2
Dc2 = Dc2 - ( (RefPt(1) - node(3,1))/c2**2 ) * Dd3

! d2 = (node(2,1) - RefPt(1))/c1
node_ad(2,1) = Dd2/c1
RefPt_ad(1) = RefPt_ad(1) - Dd2/c1
Dc1 = -( (node(2,1) - RefPt(1))/c1**2  ) *Dd2

! d1 = (RefPt(1) - node(1,1))/c1
RefPt_ad(1) = RefPt_ad(1) + Dd1/c1
node_ad(1,1) = -Dd1/c1
Dc1 = Dc1 - ( (RefPt(1) - node(1,1))/c1**2  ) * Dd1

! c2 = node(4,1) - node(3,1)
node_ad(4,1) = node_ad(4,1) + Dc2
node_ad(3,1) = node_ad(3,1) - Dc2

! c1 = node(2,1) - node(1,1)
node_ad(2,1) = node_ad(2,1) + Dc1
node_ad(1,1) = node_ad(1,1) - Dc1

BiLinear2DVal_ad = 0.0_wp

end subroutine BiLinear2D_ad


!---------------------------------------------------------------------
! function BiLinearExp3D_tl
!---------------------------------------------------------------------
!
! Start of RoboDoc comments >>>>
!****F* Interpolation/BiLinearExp3D_tl
!
! Name
! BiLinearExp3D_tl
!
! Purpose
! 3D interpolation usig a bilinear interpoltion horizontally and an
! exponetial interpolation vertically.
!
! This function makes use of "ExpInt1D" and "BiLinear2D".
!
!
! Call
! val =  BiLinearExp3D_tl(RefPt, node)
! val =  BiLinearExp3D_tl(RefPt, node, Npt)
!
! Input
! RefPt    - reference point, the field shall be interpolated at this
!            point, RefPoint(1) - longitude
!                   RefPoint(2) - latitude
!                   RefPoint(3) - altitude
! node     - array containing the coordinates and the values to be interpolated
!              node(i,j,k) - i=1, ... , 4 : 4 nodes at the cell corner
!                            j=1, 2       : grid level, j=1 - lower level
!                                                       j=2 - upper level
!                            k=1, ... , 4 : k=1 - longitude
!                                           k=2 - latitude
!                                           k=3 - altitude
!                                           k=4 - value, refractivity, ...
! Npt      - number of reference points, Npt = 3 or Npt = 4
!            optional
!
! Output
! BiLinearExp3D (return value) - interpolated value
!
! External References
! None
!
! Record of revisions
! Date        Programmer   Description of change
! ====        ==========   =====================
! 01.03.2013  M. Bender    new
!---------------------------------------------------------------------
! Source
!---------------------------------------------------------------------
function BiLinearExp3D_tl(RefPt, node, RefPt_tl, node_tl, Npt)

! List of calling arguments:
real (wp)                                     :: BiLinearExp3D_tl
real (wp), dimension(:), intent(in)           :: RefPt(1:3)
real (wp), dimension(1:4,1:2,1:4), intent(in) :: node
real (wp), dimension(1:3), intent(in)         :: RefPt_tl
real (wp), dimension(1:4,1:2,1:4), intent(in) :: node_tl
integer, optional                                        :: Npt

! List of local variables:
real (wp), dimension(1:2,1:2) :: Vnode, DVnode
real (wp)                     :: VPt, DVPt
real (wp), dimension(1:4,1:3) :: Hnode, DHnode
real (wp), dimension(1:2)     :: HPt, DHPt

integer :: p, i
!---------------------------------------------------------------------

! Number of reference points: 3 or 4
p = 4
if (present(Npt)) then
   if (Npt .eq. 3) p = 3
end if

! Vertical interpolation using "ExpInt1D"
VPt = RefPt(3)  ! interpolate to this altitude
DVPt = RefPt_tl(3)

do i=1, p

   Vnode(1,1) = node(i,1,3)  ! altitude node i, lower level
   Vnode(2,1) = node(i,2,3)  ! altitude node i, upper level

   Vnode(1,2) = node(i,1,4)  ! refractivity node i, lower level
   Vnode(2,2) = node(i,2,4)  ! refractivity node i, upper level

   DVnode(1,1) = node_tl(i,1,3)
   DVnode(2,1) = node_tl(i,2,3)
   DVnode(1,2) = node_tl(i,1,4)
   DVnode(2,2) = node_tl(i,2,4)

   Hnode(i,3) = ExpInt1D(VPt, Vnode)  ! interpolated vertical refrac.

   DHnode(i,3) = ExpInt1D_tl(VPt, Vnode, DVPt, DVnode)

end do

! Horizontal interpolation using BiLin2D
HPt = RefPt(1:2)  ! interpolate to this (x,y)-coordinate
DHPt = RefPt_tl(1:2)

! x,y coordinates
do i=1, p
   Hnode(i,1) = node(i,1,1)   ! node i, x
   Hnode(i,2) = node(i,1,2)   ! node i, y

   DHnode(i,1) = node_tl(i,1,1)
   DHnode(i,2) = node_tl(i,1,2)

   ! Longitude and latitude are assumed to be identical for both
   ! levels, therefore
   ! node(i,2,1) = node(i,1,1) and
   ! node(i,2,2) = node(i,1,2) are never used and the corresponding
   ! node_tl elements are not used.
   !node_tl(i,2,1) = node_tl(i,1,1)
   !node_tl(i,2,2) = node_tl(i,1,2)

end do

if (p .eq. 3) then
   ! Copy node 2 to node 4 as required by BiLinear2D
   Hnode(4,:) = Hnode(2,:)
   DHnode(4,:) = DHnode(2,:)
end if

! BiLinearExp3D = BiLinear2D(HPt, Hnode)
!BiLinearExp3D_tl = BiLinearDiff2D_tl(HPt, Hnode, DHnode, p)
BiLinearExp3D_tl = Shepard2D_tl(HPt, Hnode(1:p,:), DHnode(1:p,:))

!if (BiLinearExp3D_tl /= BiLinearExp3D_tl) then
!   ! BiLinearExp3D_tl = NaN
!   write(*,*) 'BiLinearExp3D_tl = NaN'
!end if
!write(*,*) 'BiLinearExp3D_tl = ', BiLinearExp3D_tl

!!$write(*,*) 'BiLinearExp3D_tl', char(10), &
!!$     'max DHnode(i,3) = ', maxval(DHnode) , char(10), &
!!$     'BiLinearDiff2D_tl = ', BiLinearExp3D_tl

end function BiLinearExp3D_tl
!****
! <<<<<<<< End of RoboDoc comments


!---------------------------------------------------------------------
! subroutine BiLinearExp3D_ad
!---------------------------------------------------------------------
!
! Start of RoboDoc comments >>>>
!****F* Interpolation/BiLinearExp3D_ad
!
! Name
! BiLinearExp3D_ad
!
! Purpose
! 3D interpolation usig a bilinear interpoltion horizontally and an
! exponetial interpolation vertically.
!
! This function makes use of "ExpInt1D" and "BiLinear2D".
!
!
! Call
! val =  BiLinearExp3D_ad(RefPt, node)
! val =  BiLinearExp3D_ad(RefPt, node, Npt)
!
! Input
! RefPt    - reference point, the field shall be interpolated at this
!            point, RefPoint(1) - longitude
!                   RefPoint(2) - latitude
!                   RefPoint(3) - altitude
! node     - array containing the coordinates and the values to be interpolated
!              node(i,j,k) - i=1, ... , 4 : 4 nodes at the cell corner
!                            j=1, 2       : grid level, j=1 - lower level
!                                                       j=2 - upper level
!                            k=1, ... , 4 : k=1 - longitude
!                                           k=2 - latitude
!                                           k=3 - altitude
!                                           k=4 - value, refractivity, ...
! Npt      - number of reference points, Npt = 3 or Npt = 4
!            optional
!
! Output
! BiLinearExp3D (return value) - interpolated value
!
! External References
! None
!
! Record of revisions
! Date        Programmer   Description of change
! ====        ==========   =====================
! 01.03.2012  M. Bender    new
! 07.06.2012  M. Bende     first running version
!---------------------------------------------------------------------
! Source
!---------------------------------------------------------------------
subroutine BiLinearExp3D_ad(RefPt, node, RefPt_ad, node_ad,   &
                            ValBiLinearExp3D_ad, Npt)

! List of calling arguments:
real (wp), intent(inout)                       :: ValBiLinearExp3D_ad
real (wp), dimension(:), intent(in)            :: RefPt(1:3)
real (wp), dimension(1:4,1:2,1:4), intent(in)  :: node
real (wp), dimension(1:3), intent(out)         :: RefPt_ad
real (wp), dimension(1:4,1:2,1:4), intent(out) :: node_ad
integer, optional                              :: Npt

! List of local variables:
real (wp), dimension(1:2,1:2) :: Vnode, DVnode
real (wp)                     :: VPt, DVPt
real (wp), dimension(1:4,1:3) :: Hnode, DHnode
real (wp), dimension(1:2)     :: HPt !, DHPt

integer :: p, i

!---------------------------------------------------------------------

node_ad  = 0.0_wp
RefPt_ad = 0.0_wp

! Number of reference points: 3 or 4
p = 4
if (present(Npt)) then
   if (Npt .eq. 3) p = 3
end if

! Vertical interpolation using "ExpInt1D"
VPt = RefPt(3)  ! interpolate to this altitude
do i=1, p
   Vnode(1,1) = node(i,1,3)  ! altitude node i, lower level
   Vnode(2,1) = node(i,2,3)  ! altitude node i, upper level

   Vnode(1,2) = node(i,1,4)  ! refractivity node i, lower level
   Vnode(2,2) = node(i,2,4)  ! refractivity node i, upper level

   Hnode(i,3) = ExpInt1D(VPt, Vnode)  ! interpolated vertical refrac.
end do

! Horizontal interpolation using BiLin2D
HPt = RefPt(1:2)  ! interpolate to this (x,y)-coordinate

! x,y coordinates
do i=1, p
   Hnode(i,1) = node(i,1,1)   ! node i, x
   Hnode(i,2) = node(i,1,2)   ! node i, y
end do
if (p .eq. 3) then
   ! Copy node 2 to node 4 as required by BiLinear2D
   Hnode(4,:) = Hnode(2,:)
end if
!---------------- End of recomputation: HPt, Hnode ----------------------

! BiLinearExp3D = BiLinearDiff2D(HPt, Hnode)
!call BiLinearDiff2D_ad(HPt, Hnode, DHnode, ValBiLinearExp3D_ad, p)
call Shepard2D_ad(HPt, Hnode(1:p,:), DHnode, ValBiLinearExp3D_ad)

if (p .eq. 3) then
   ! Copy node 2 to node 4 as required by BiLinear2D
   ! Hnode(4,:) = Hnode(2,:)
   DHnode(2,:) = DHnode(2,:) + DHnode(4,:)
end if

do i=p, 1, -1
   ! Hnode(i,2) = node(i,1,2)
   node_ad(i,1,2) = DHnode(i,2)
   ! Hnode(i,1) = node(i,1,1)
   node_ad(i,1,1) = DHnode(i,1)

   ! Longitude and latitude are assumed to be identical for both
   ! levels, therefore
   ! node(i,2,1) = node(i,1,1) and
   ! node(i,2,2) = node(i,1,2) are never used and the corresponding
   ! node_tl elements are not used.
   !node_ad(i,2,2) = 0.0D0
   !node_ad(i,2,1) = 0.0D0
end do

! HPt = RefPt(1:2)
!RefPt_ad(1:2) = DHPt

do i=p, 1, -1
   ! Copy the required data to Vnode before calling ExpInt1D_ad
   Vnode(1,1) = node(i,1,3)  ! altitude node i, lower level
   Vnode(2,1) = node(i,2,3)  ! altitude node i, upper level
   Vnode(1,2) = node(i,1,4)  ! refractivity node i, lower level
   Vnode(2,2) = node(i,2,4)  ! refractivity node i, upper level
   DVnode = 0.0_wp           ! independent for each i, don't sum

   ! Hnode(i,3) = ExpInt1D(VPt, Vnode)
   call ExpInt1D_ad(VPt, Vnode, DHnode(i,3), DVPt, DVnode)

   ! Vnode(2,2) = node(i,2,4)
   node_ad(i,2,4) = DVnode(2,2)

   ! Vnode(1,2) = node(i,1,4)
   node_ad(i,1,4) = DVnode(1,2)

   ! Vnode(2,1) = node(i,2,3)
   node_ad(i,2,3) = DVnode(2,1)

   ! Vnode(1,1) = node(i,1,3)
   node_ad(i,1,3) = DVnode(1,1)

   ! VPt = RefPt(3)
   ! Needs to be inside the loop as each call of ExpInt1D_ad contributes
   ! to RefPt_ad(3) !!
   RefPt_ad(3) = RefPt_ad(3) + DVPt
end do

ValBiLinearExp3D_ad = 0.0_wp

end subroutine BiLinearExp3D_ad


!---------------------------------------------------------------------
! function BiLinearDiff2D_tl
!---------------------------------------------------------------------
!
! Start of RoboDoc comments >>>>
!****F* Interpolation/BiLinearDiff2D_tl
!
! Name
! BiLinearDiff2D_tl
!
! Call
! val =  BiLinearDiff2D_tl(RefPt, node)
!
! Purpose
!
! Tangent-linear code:
! Derivatives are computed with respect to the atmospheric quantities only.
! All coordinates are fixed and no derivatives with respect to coordinates
! will be regarded.
! RefPt_tl = 0            reference point
! nodein_tl(:,1) = 0      x-coordinate, longitude
! nodein_tl(:,2) = 0      y-coordinate, latitude
! nodein_tl(:,3) <> 0     p, T, q to be interpolated
!
! 2D interpolation usig a bilinear interpoltion between four reference
! values. Any set of nodes is possible as long as (y2-y1) <> 0 and
! (x2-x1) <> 0 and (x4-x3) <> 0, i. e. 3 or 4 points are on a straoght line.
! If only 3 reference points ara available, the routine can be called with
! P4 = P2 (node 4 = node 2).
!
! In case of 3 reference points there is a singularity at  x = x2
! i.e. RefPt(1) = x2, which is replaced by the mean value
! od two points close to the singularity which can be evaluated.
!
!
!    y ^
!      |
!      |   node3                  node4             y3 > y4, x3 < x4
!      |
!      |          RfPt
!      |
!      |   node1                  node2             y1 > y2, x1 < x2
!      |
!      ---------------------------------------> x
!
!          x1 < x2                x2 > x1
!          y1 < y3                y2 < y4
!
! Test program: .../gpstomo/trunk/Tomo/TomoTest/Interpol/TestBiLinear2D.f90
!
!

! Input
! RefPt    - reference point, the field shall be interpolated at this
!            point, RefPoint(1) - x, longitude
!                   RefPoint(2) - y, latitude
! node     - array
!            node(i,:), i=1, ... , 4 : 4 nodes at the cell corner
!            node(i,j), j=1, 2, 3    : j = 1 - x
!                                      j = 2 - y
!                                      j = 3 - f_j = f(x,y)
!
! Output
! BiLinearDiff2D (return value) - interpolated value
!
! External References
! None
!
! Record of revisions
! Date        Programmer   Description of change
! ====        ==========   =====================
! 17.05.2013  M. Bender    new
!---------------------------------------------------------------------
! Source
!---------------------------------------------------------------------
function BiLinearDiff2D_tl(RefPt, nodein, nodein_tl, Nnode)

! List of calling arguments:
real (wp)                                 :: BiLinearDiff2D_tl
real (wp), dimension(1:2), intent(in)     :: RefPt
real (wp), dimension(1:4,1:3), intent(in) :: nodein
real (wp), dimension(1:4,1:3), intent(in) :: nodein_tl
integer,  intent(in), optional            :: Nnode

! List of local variables:
real (wp), dimension(1:4,1:3) :: node, node_tl
real (wp) :: c1, c2, d1, d2, d3, d4
real (wp) :: ya, yb, fa, fb, dy !, a1, a2, a3
!real (wp) :: Dc1, Dc2, Dd1, Dd2, Dd3, Dd4
real (wp) :: Dfa, Dfb           !, Dya, Dyb, Ddy, Da1, Da2, Da3
!real (wp), dimension(1:2) :: Rpt
!real (wp) :: epsilon, I1, I2, DI1, DI2
!---------------------------------------------------------------------

node    = nodein
node_tl = nodein_tl

if (present(Nnode)) then
   if (Nnode == 3) then
      !if (abs(lon(1)-lon(2)) < 1.0E-6_wp) then
      if (abs(node(1,1)-node(2,1)) < 1.0E-6_wp) then
         ! change 2 and 3 to avoid singularity, copy 3 to 4
         node(4,:) = node(3,:)
         node(3,:) = node(2,:)
         node(2,:) = node(4,:)
         node_tl(4,:) = node_tl(3,:)
         node_tl(3,:) = node_tl(2,:)
         node_tl(2,:) = node_tl(4,:)
      !else if (abs(lon(2)-lon(3)) < 1.0E-6_wp) then
      else if (abs(node(2,1)-node(3,1)) < 1.0E-6_wp) then
         ! copy 1 to 4 to avoid singularity
         node(4,:) = node(1,:)
         node_tl(4,:) = node_tl(1,:)
      else
         node(4,:) = node(2,:)
         node_tl(4,:) = node_tl(2,:)
      end if
   end if
end if

! c1 = x2 - x1
c1 = node(2,1) - node(1,1)

! c2 = x4 - x3
c2 = node(4,1) - node(3,1)

! d1 = (x - x1)/(x2-x1) =  (x - x1)/c1
d1 = (RefPt(1) - node(1,1))/c1

! d2 = (x2 - x)/(x2-x1) =  (x2 - x)/c1
d2 = (node(2,1) - RefPt(1))/c1

! d3 = (x - x3)/(x4-x3) =  (x - x3)/c2
d3 = (RefPt(1) - node(3,1))/c2

! d4 = (x4 - x)/(x4-x3) =  (x4 - x)/c2
d4 = (node(4,1) - RefPt(1))/c2

!write(*,*) 'BiLinear2D> c1, c2, d1, d2, d3, d4 : ', c1, c2, d1, d2, d3, d4

fa = d2*node(1,3) + d1*node(2,3)
Dfa = d2*node_tl(1,3) + d1*node_tl(2,3)

fb = d4*node(3,3) + d3*node(4,3)
Dfb = d4*node_tl(3,3) + d3*node_tl(4,3)

ya = d2*node(1,2) + d1*node(2,2)

yb = d4*node(3,2) + d3*node(4,2)

dy = yb - ya
!call SaveDiff('tl', Dfa, Dfb, Dyb, Dya, RefPt_tl(2), Ddy)


if (abs(dy) .lt. 1.0E-8_wp) then
   ! slant%Nnghb = 3 : reference point is very close to point 2
!!$   write(*,*) 'Singularitaet: dy = ', dy
!!$   if (abs(RefPt(2)-node(1,2)) < 1.0E-6_wp) then
!!$      BiLinearDiff2D_tl = node_tl(1,3)
!!$   else if (abs(RefPt(2)-node(2,2)) < 1.0E-6_wp) then
!!$      BiLinearDiff2D_tl = node_tl(2,3)
!!$   else if (abs(RefPt(2)-node(3,2)) < 1.0E-6_wp) then
!!$      BiLinearDiff2D_tl = node_tl(3,3)
!!$   else if (abs(RefPt(2)-node(4,2)) < 1.0E-6_wp) then
!!$      BiLinearDiff2D_tl = node_tl(4,3)
!!$   else
      BiLinearDiff2D_tl = node_tl(2,3)
  ! end if
else
   !BiLinear2D = fa*(yb- RefPt(2))/dy + fb*(RefPt(2)-ya)/dy
   BiLinearDiff2D_tl = Dfa*(yb- RefPt(2))/dy + Dfb*(RefPt(2)-ya)/dy
end if

if (BiLinearDiff2D_tl /= BiLinearDiff2D_tl) then
   ! BiLinearDiff2D_tl = NaN
   write(*,*) ' BiLinearDiff2D_tl = NaN'
end if


end function BiLinearDiff2D_tl
!****
! <<<<<<<< End of RoboDoc comments


!---------------------------------------------------------------------
! function BiLinearDiff2D_full_tl
!---------------------------------------------------------------------
!
! Start of RoboDoc comments >>>>
!****F* Interpolation/BiLinearDiff2D_full_tl
!
! Name
! BiLinearDiff2D_full_tl
!
! Call
! val =  BiLinearDiff2D_full_tl(RefPt, node)
!
! Purpose
! 2D interpolation usig a bilinear interpoltion between four reference
! values. Any set of nodes is possible as long as (y2-y1) <> 0 and
! (x2-x1) <> 0 and (x4-x3) <> 0, i. e. 3 or 4 points are on a straoght line.
! If only 3 reference points ara available, the routine can be called with
! P4 = P2 (node 4 = node 2).
!
! In case of 3 reference points there is a singularity at  x = x2
! i.e. RefPt(1) = x2, which is replaced by the mean value
! od two points close to the singularity which can be evaluated.
!
!
!    y ^
!      |
!      |   node3                  node4             y3 > y4, x3 < x4
!      |
!      |          RfPt
!      |
!      |   node1                  node2             y1 > y2, x1 < x2
!      |
!      ---------------------------------------> x
!
!          x1 < x2                x2 > x1
!          y1 < y3                y2 < y4
!
! Test program: .../gpstomo/trunk/Tomo/TomoTest/Interpol/TestBiLinear2D.f90
!
!

! Input
! RefPt    - reference point, the field shall be interpolated at this
!            point, RefPoint(1) - x, longitude
!                   RefPoint(2) - y, latitude
! node     - array
!            node(i,:), i=1, ... , 4 : 4 nodes at the cell corner
!            node(i,j), j=1, 2, 3    : j = 1 - x
!                                      j = 2 - y
!                                      j = 3 - f_j = f(x,y)
!
! Output
! BiLinearDiff2D (return value) - interpolated value
!
! External References
! None
!
! Record of revisions
! Date        Programmer   Description of change
! ====        ==========   =====================
! 17.05.2013  M. Bender    new
!---------------------------------------------------------------------
! Source
!---------------------------------------------------------------------
function BiLinearDiff2D_full_tl(RefPt, nodein, RefPt_tl, nodein_tl, Nnode)

! List of calling arguments:
real (wp)                                 :: BiLinearDiff2D_full_tl
real (wp), dimension(1:2), intent(in)     :: RefPt
real (wp), dimension(1:4,1:3), intent(in) :: nodein
real (wp), dimension(1:2), intent(in)     :: RefPt_tl
real (wp), dimension(1:4,1:3), intent(in) :: nodein_tl
integer,  intent(in), optional            :: Nnode

! List of local variables:
real (wp), dimension(1:4,1:3) :: node, node_tl
real (wp) :: c1, c2, d1, d2, d3, d4
real (wp) :: ya, yb, fa, fb, dy      !, a1, a2, a3
real (wp) :: Dc1, Dc2, Dd1, Dd2, Dd3, Dd4
real (wp) :: Dya, Dyb, Dfa, Dfb, Ddy !, Da1, Da2, Da3
!real (wp), dimension(1:2) :: Rpt
!real (wp) :: epsilon, I1, I2, DI1, DI2
!---------------------------------------------------------------------

node    = nodein
node_tl = nodein_tl

if (present(Nnode)) then
   if (Nnode == 3) then
      !if (abs(lon(1)-lon(2)) < 1.0E-6_wp) then
      if (abs(node(1,1)-node(2,1)) < 1.0E-6_wp) then
         ! change 2 and 3 to avoid singularity, copy 3 to 4
         node(4,:) = node(3,:)
         node(3,:) = node(2,:)
         node(2,:) = node(4,:)
         node_tl(4,:) = node_tl(3,:)
         node_tl(3,:) = node_tl(2,:)
         node_tl(2,:) = node_tl(4,:)
      !else if (abs(lon(2)-lon(3)) < 1.0E-6_wp) then
      else if (abs(node(2,1)-node(3,1)) < 1.0E-6_wp) then
         ! copy 1 to 4 to avoid singularity
         node(4,:) = node(1,:)
         node_tl(4,:) = node_tl(1,:)
      else
         node(4,:) = node(2,:)
         node_tl(4,:) = node_tl(2,:)
      end if
   end if
end if

! c1 = x2 - x1
c1 = node(2,1) - node(1,1)
Dc1 = node_tl(2,1) - node_tl(1,1)

! c2 = x4 - x3
c2 = node(4,1) - node(3,1)
Dc2 = node_tl(4,1) - node_tl(3,1)

! d1 = (x - x1)/(x2-x1) =  (x - x1)/c1
d1 = (RefPt(1) - node(1,1))/c1
Dd1 = ((RefPt_tl(1) - node_tl(1,1))/c1) - ((RefPt(1) - node(1,1))/c1**2) * Dc1

! d2 = (x2 - x)/(x2-x1) =  (x2 - x)/c1
d2 = (node(2,1) - RefPt(1))/c1
Dd2 = ((node_tl(2,1) - RefPt_tl(1))/c1) - ((node(2,1) - RefPt(1))/c1**2) * Dc1

! d3 = (x - x3)/(x4-x3) =  (x - x3)/c2
d3 = (RefPt(1) - node(3,1))/c2
Dd3 = ((RefPt_tl(1) - node_tl(3,1))/c2) - ((RefPt(1) - node(3,1))/c2**2) * Dc2

! d4 = (x4 - x)/(x4-x3) =  (x4 - x)/c2
d4 = (node(4,1) - RefPt(1))/c2
Dd4 = ((node_tl(4,1) - RefPt_tl(1))/c2) - ((node(4,1) - RefPt(1))/c2**2) * Dc2

!write(*,*) 'BiLinear2D> c1, c2, d1, d2, d3, d4 : ', c1, c2, d1, d2, d3, d4

fa = d2*node(1,3) + d1*node(2,3)
Dfa = d2*node_tl(1,3) + Dd2*node(1,3) + d1*node_tl(2,3) + Dd1*node(2,3)

fb = d4*node(3,3) + d3*node(4,3)
Dfb = d4*node_tl(3,3) + Dd4*node(3,3) + d3*node_tl(4,3) + Dd3*node(4,3)

ya = d2*node(1,2) + d1*node(2,2)
Dya = d2*node_tl(1,2) + Dd2*node(1,2) + d1*node_tl(2,2) + Dd1*node(2,2)

yb = d4*node(3,2) + d3*node(4,2)
Dyb = d4*node_tl(3,2) + Dd4*node(3,2) + d3*node_tl(4,2) + Dd3*node(4,2)

dy = yb - ya
Ddy = Dyb - Dya
!call SaveDiff('tl', Dfa, Dfb, Dyb, Dya, RefPt_tl(2), Ddy)

!!$if (abs(dy) .lt. 1.0E-7_wp .and. abs(fb-fa) .lt. 1.0E-7_wp) then
!!$   ! x = x2, singularity, should happen only for 3 point interpolation
!!$   ! The singularity is replaced by the mean value of two points close
!!$   ! to the singularity, i.e. computed for RefPt(1) + epsilon and
!!$   ! RefPt(1) - epsilon.
!!$   !write(*,*) 'BiLinearDiff2D_tl> x = x2'
!!$
!!$   write(*,*) 'Singularitaet: dy = ', dy
!!$
!!$   ! Assuming latitude and longitude given in radian, a difference
!!$   ! of 9*10^-5 rad is equivalent to a distance of ~ 10 m,
!!$   ! i. e. epsilon = 1.0D-5 should be very close to the singularity
!!$   ! and the error of the interpolation is < 0.1 %
!!$   epsilon = 1.0E-5_wp
!!$
!!$   ! Interpolation 1
!!$   ! x + epsilon (RefPt(1) + epsilon)
!!$   !
!!$   ! d1 = (x - x1)/(x2-x1) =  (x - x1)/c1
!!$   d1 = (RefPt(1)+epsilon - node(1,1))/c1
!!$   Dd1 = ((RefPt_tl(1) - node_tl(1,1))/c1) -               &
!!$         ((RefPt(1) + epsilon - node(1,1))/c1**2) * Dc1
!!$
!!$   ! d2 = (x2 - x)/(x2-x1) =  (x2 - x)/c1
!!$   d2 = (node(2,1) - RefPt(1)-epsilon)/c1
!!$   Dd2 = ((node_tl(2,1) - RefPt_tl(1))/c1) -               &
!!$         ((node(2,1) - RefPt(1) - epsilon)/c1**2) * Dc1
!!$
!!$   ! d3 = (x - x3)/(x4-x3) =  (x - x3)/c2
!!$   d3 = (RefPt(1)+epsilon - node(3,1))/c2
!!$   Dd3 = ((RefPt_tl(1) - node_tl(3,1))/c2) -               &
!!$         ((RefPt(1) + epsilon - node(3,1))/c2**2) * Dc2
!!$
!!$   ! d4 = (x4 - x)/(x4-x3) =  (x4 - x)/c2
!!$   d4 = (node(4,1) - RefPt(1)-epsilon)/c2
!!$   Dd4 = ((node_tl(4,1) - RefPt_tl(1))/c2) -               &
!!$         ((node(4,1) - RefPt(1) - epsilon)/c2**2) * Dc2
!!$
!!$   fa = d2*node(1,3) + d1*node(2,3)
!!$   Dfa = d2*node_tl(1,3) + Dd2*node(1,3) + d1*node_tl(2,3) + Dd1*node(2,3)
!!$
!!$   fb = d4*node(3,3) + d3*node(4,3)
!!$   Dfb = d4*node_tl(3,3) + Dd4*node(3,3) + d3*node_tl(4,3) + Dd3*node(4,3)
!!$
!!$   ya = d2*node(1,2) + d1*node(2,2)
!!$   Dya = d2*node_tl(1,2) + Dd2*node(1,2) + d1*node_tl(2,2) + Dd1*node(2,2)
!!$
!!$   yb = d4*node(3,2) + d3*node(4,2)
!!$   Dyb = d4*node_tl(3,2) + Dd4*node(3,2) + d3*node_tl(4,2) + Dd3*node(4,2)
!!$
!!$   dy = yb - ya
!!$   Ddy = Dyb - Dya
!!$
!!$   I1 = fa*(yb- RefPt(2))/dy + fb*(RefPt(2)-ya)/dy
!!$   DI1 = Dfa*(yb- RefPt(2))/dy + fa*(Dyb- RefPt_tl(2))/dy - &
!!$         fa*(yb- RefPt(2))/dy**2 * Ddy +                    &
!!$         Dfb*(RefPt(2)-ya)/dy + fb*(RefPt_tl(2)-Dya)/dy -   &
!!$         fb*(RefPt(2)-ya)/dy**2 * Ddy
!!$
!!$   ! Interpolation 2
!!$   ! x - epsilon (RefPt(1) - epsilon)
!!$   !
!!$   ! d1 = (x - x1)/(x2-x1) =  (x - x1)/c1
!!$   d1 = (RefPt(1)-epsilon - node(1,1))/c1
!!$   Dd1 = ((RefPt_tl(1) - node_tl(1,1))/c1) -                 &
!!$         ((RefPt(1) - epsilon - node(1,1))/c1**2) * Dc1
!!$
!!$   ! d2 = (x2 - x)/(x2-x1) =  (x2 - x)/c1
!!$   d2 = (node(2,1) - RefPt(1)+epsilon)/c1
!!$   Dd2 = ((node_tl(2,1) - RefPt_tl(1))/c1) -                 &
!!$         ((node(2,1) - RefPt(1) + epsilon)/c1**2) * Dc1
!!$
!!$   ! d3 = (x - x3)/(x4-x3) =  (x - x3)/c2
!!$   d3 = (RefPt(1)-epsilon - node(3,1))/c2
!!$   Dd3 = ((RefPt_tl(1) - node_tl(3,1))/c2) -                 &
!!$         ((RefPt(1) - epsilon - node(3,1))/c2**2) * Dc2
!!$
!!$   ! d4 = (x4 - x)/(x4-x3) =  (x4 - x)/c2
!!$   d4 = (node(4,1) - RefPt(1)+epsilon)/c2
!!$   Dd4 = ((node_tl(4,1) - RefPt_tl(1))/c2) -                 &
!!$         ((node(4,1) - RefPt(1) + epsilon)/c2**2) * Dc2
!!$
!!$   fa = d2*node(1,3) + d1*node(2,3)
!!$   Dfa = d2*node_tl(1,3) + Dd2*node(1,3) + d1*node_tl(2,3) + Dd1*node(2,3)
!!$
!!$   fb = d4*node(3,3) + d3*node(4,3)
!!$   Dfb = d4*node_tl(3,3) + Dd4*node(3,3) + d3*node_tl(4,3) + Dd3*node(4,3)
!!$
!!$   ya = d2*node(1,2) + d1*node(2,2)
!!$   Dya = d2*node_tl(1,2) + Dd2*node(1,2) + d1*node_tl(2,2) + Dd1*node(2,2)
!!$
!!$   yb = d4*node(3,2) + d3*node(4,2)
!!$   Dyb = d4*node_tl(3,2) + Dd4*node(3,2) + d3*node_tl(4,2) + Dd3*node(4,2)
!!$
!!$   dy = yb - ya
!!$   Ddy = Dyb - Dya
!!$
!!$   I2 = fa*(yb- RefPt(2))/dy + fb*(RefPt(2)-ya)/dy
!!$   DI2 = Dfa*(yb- RefPt(2))/dy + fa*(Dyb- RefPt_tl(2))/dy - &
!!$         fa*(yb- RefPt(2))/dy**2 * Ddy +                    &
!!$         Dfb*(RefPt(2)-ya)/dy + fb*(RefPt_tl(2)-Dya)/dy -   &
!!$         fb*(RefPt(2)-ya)/dy**2 * Ddy
!!$
!!$   ! BiLinearDiff2D = 0.50D0 * (I1 + I2)
!!$   BiLinearDiff2D_tl = 0.50_wp * (DI1 + DI2)
if (abs(dy) .lt. 1.0E-8_wp) then
   ! slant%Nnghb = 3 : reference point is very close to point 2
   write(*,*) 'Singularitaet: dy = ', dy
   BiLinearDiff2D_full_tl = node_tl(2,3)
else
   !BiLinear2D = fa*(yb- RefPt(2))/dy + fb*(RefPt(2)-ya)/dy
   BiLinearDiff2D_full_tl = Dfa*(yb- RefPt(2))/dy + fa*(Dyb- RefPt_tl(2))/dy - &
                       fa*(yb- RefPt(2))/dy**2 * Ddy +                    &
                       Dfb*(RefPt(2)-ya)/dy + fb*(RefPt_tl(2)-Dya)/dy -   &
                       fb*(RefPt(2)-ya)/dy**2 * Ddy
end if

if (BiLinearDiff2D_full_tl /= BiLinearDiff2D_full_tl) then
   ! BiLinearDiff2D_tl = NaN
   write(*,*) ' BiLinearDiff2D_full_tl = NaN'
end if


end function BiLinearDiff2D_full_tl
!****
! <<<<<<<< End of RoboDoc comments


!---------------------------------------------------------------------
! function BiLinearDiff2D_ad
!---------------------------------------------------------------------
!
! Start of RoboDoc comments >>>>
!****F* Interpolation/BiLinearDiff2D_ad
!
! Name
! BiLinearDiff2D_ad
!
! Call
! val =  BiLinearDiff2D(RefPt, node)
!
! Purpose
!
! Adjoint code:
! Derivatives are computed with respect to the atmospheric quantities only.
! All coordinates are fixed and no derivatives with respect to coordinates
! will be regarded.
! node_ad(:,1) = 0      x-coordinate, longitude
! node_ad(:,2) = 0      y-coordinate, latitude
! node_ad(:,3) <> 0     p, T, q to be interpolated
! RefPt_ad = 0  - reference point is fixed
!
! 2D interpolation usig a bilinear interpoltion between four reference
! values. Any set of nodes is possible as long as (y2-y1) <> 0 and
! (x2-x1) <> 0 and (x4-x3) <> 0, i. e. 3 or 4 points are on a straoght line.
! If only 3 reference points ara available, the routine can be called with
! P4 = P2 (node 4 = node 2).
!
! In case of 3 reference points there is a singularity at  x = x2
! i.e. RefPt(1) = x2, which is replaced by the mean value
! od two points close to the singularity which can be evaluated.
!
!
!    y ^
!      |
!      |   node3                  node4             y3 > y4, x3 < x4
!      |
!      |          RfPt
!      |
!      |   node1                  node2             y1 > y2, x1 < x2
!      |
!      ---------------------------------------> x
!
!          x1 < x2                x2 > x1
!          y1 < y3                y2 < y4
!
! Test program: .../gpstomo/trunk/Tomo/TomoTest/Interpol/TestBiLinear2D.f90
!
!

! Input
! RefPt    - reference point, the field shall be interpolated at this
!            point, RefPoint(1) - x, longitude
!                   RefPoint(2) - y, latitude
! node     - array
!            node(i,:), i=1, ... , 4 : 4 nodes at the cell corner
!            node(i,j), j=1, 2, 3    : j = 1 - x
!                                      j = 2 - y
!                                      j = 3 - f_j = f(x,y)
!
! Output
! BiLinearDiff2D (return value) - interpolated value
!
! External References
! None
!
! Record of revisions
! Date        Programmer   Description of change
! ====        ==========   =====================
! 17.05.2013  M. Bender    new
!---------------------------------------------------------------------
! Source
!---------------------------------------------------------------------
subroutine BiLinearDiff2D_ad(RefPt, nodein, node_ad,   &
                             BiLinearDiff2DVal_ad, Nnode)

! List of calling arguments:
real (wp), intent(inout)                     :: BiLinearDiff2DVal_ad
real (wp), dimension(1:2), intent(in)        :: RefPt
real (wp), dimension(1:4,1:3), intent(in)    :: nodein
real (wp), dimension(1:4,1:3), intent(out)   :: node_ad
integer,  intent(in), optional               :: Nnode

! List of local variables:
real (wp), dimension(1:4,1:3) :: node
!real (wp) :: BiLinear2D
real (wp) :: c1, c2, d1, d2, d3, d4
real (wp) :: ya, yb, fa, fb, dy !, a1, a2, a3
!real (wp) :: Dc1, Dc2, Dd1, Dd2, Dd3, Dd4
real (wp) :: Dfa, Dfb !, Dya, Dyb, Ddy, Da1, Da2, Da3
!real (wp) :: epsilon, I1, I2, DI1, DI2
!---------------------------------------------------------------------

node = nodein

if (present(Nnode)) then
   if (Nnode == 3) then
      !if (abs(lon(1)-lon(2)) < 1.0E-6_wp) then
      if (abs(node(1,1)-node(2,1)) < 1.0E-6_wp) then
         ! change 2 and 3 to avoid singularity, copy 3 to 4
         node(4,:) = node(3,:)
         node(3,:) = node(2,:)
         node(2,:) = node(4,:)
      !else if (abs(lon(2)-lon(3)) < 1.0E-6_wp) then
      else if (abs(node(2,1)-node(3,1)) < 1.0E-6_wp) then
         ! copy 1 to 4 to avoid singularity
         node(4,:) = node(1,:)
      else
         node(4,:) = node(2,:)
      end if
   end if
end if

! c1 = x2 - x1
c1 = node(2,1) - node(1,1)
! c2 = x4 - x3
c2 = node(4,1) - node(3,1)
! d1 = (x - x1)/(x2-x1) =  (x - x1)/c1
d1 = (RefPt(1) - node(1,1))/c1
! d2 = (x2 - x)/(x2-x1) =  (x2 - x)/c1
d2 = (node(2,1) - RefPt(1))/c1
! d3 = (x - x3)/(x4-x3) =  (x - x3)/c2
d3 = (RefPt(1) - node(3,1))/c2
! d4 = (x4 - x)/(x4-x3) =  (x4 - x)/c2
d4 = (node(4,1) - RefPt(1))/c2
!write(*,*) 'BiLinear2D> c1, c2, d1, d2, d3, d4 : ', c1, c2, d1, d2, d3, d4

fa = d2*node(1,3) + d1*node(2,3)
fb = d4*node(3,3) + d3*node(4,3)
ya = d2*node(1,2) + d1*node(2,2)
yb = d4*node(3,2) + d3*node(4,2)
dy = yb - ya

! adjoint code
node_ad = 0.0_wp
node = nodein

if (abs(dy) .lt. 1.0E-8_wp) then
   ! slant%Nnghb = 3 : reference point is very close to one node
!!$   write(*,*) 'BiLinearDiff2D_ad> singularity ...'
!!$   if (abs(RefPt(2)-node(1,2)) < 1.0E-6_wp) then
!!$      node_ad(1,3) = BiLinearDiff2DVal_ad
!!$   else if (abs(RefPt(2)-node(2,2)) < 1.0E-6_wp) then
!!$      node_ad(2,3) = BiLinearDiff2DVal_ad
!!$   else if (abs(RefPt(2)-node(3,2)) < 1.0E-6_wp) then
!!$      node_ad(3,3) = BiLinearDiff2DVal_ad
!!$   else if (abs(RefPt(2)-node(4,2)) < 1.0E-6_wp) then
!!$      node_ad(4,3) = BiLinearDiff2DVal_ad
!!$   else
      node_ad(2,3) = BiLinearDiff2DVal_ad
   !end if
else
   ! No singularity, interpolate ...

   ! BiLinear2D = fa*(yb- RefPt(2))/dy + fb*(RefPt(2)-ya)/dy
   Dfa = ( (yb-RefPt(2))/dy ) * BiLinearDiff2DVal_ad
   Dfb = ( (RefPt(2)-ya)/dy ) * BiLinearDiff2DVal_ad

   ! dy = yb - ya
   !Dyb = Dyb + Ddy
   !Dya = Dya - Ddy

   ! yb = d4*node(3,2) + d3*node(4,2)
   !Dd4 = node(3,2) * Dyb
   !Dd3 = node(4,2) * Dyb
   !node_ad(3,2) = d4 * Dyb
   !node_ad(4,2) = d3 * Dyb

   !write(*,*) 'BiLinear2D_ad> BiLinear2DVal_ad = ', BiLinear2DVal_ad
   !write(*,*) 'BiLinear2D_ad> Dyb              = ', Dyb


   ! ya = d2*node(1,2) + d1*node(2,2)
   !Dd2 = node(1,2) * Dya
   !Dd1 = node(2,2) * Dya
   !node_ad(1,2) = d2 * Dya
   !node_ad(2,2) = d1 * Dya

   ! fb = d4*node(3,3) + d3*node(4,3)
   !Dd4 = Dd4 + node(3,3) * Dfb
   !Dd3 = Dd3 + node(4,3) * Dfb
   node_ad(3,3) = d4 * Dfb
   node_ad(4,3) = d3 * Dfb

   ! fa = d2*node(1,3) + d1*node(2,3)
   !Dd2 = Dd2 + node(1,3) * Dfa
   !Dd1 = Dd1 + node(2,3) * Dfa
   node_ad(1,3) = d2 * Dfa
   node_ad(2,3) = d1 * Dfa

   ! d4 = (node(4,1) - RefPt(1))/c2
   !node_ad(4,1) = Dd4/c2
   !RefPt_ad(1) = -Dd4/c2
   !Dc2 = -( (node(4,1) - RefPt(1))/c2**2 ) * Dd4

   ! d3 = (RefPt(1) - node(3,1))/c2
   !RefPt_ad(1) = RefPt_ad(1) + Dd3/c2
   !node_ad(3,1) = -Dd3/c2
   !Dc2 = Dc2 - ( (RefPt(1) - node(3,1))/c2**2 ) * Dd3

   ! d2 = (node(2,1) - RefPt(1))/c1
   !node_ad(2,1) = Dd2/c1
   !RefPt_ad(1) = RefPt_ad(1) - Dd2/c1
   !Dc1 = -( (node(2,1) - RefPt(1))/c1**2  ) *Dd2

   ! d1 = (RefPt(1) - node(1,1))/c1
   !RefPt_ad(1) = RefPt_ad(1) + Dd1/c1
   !node_ad(1,1) = -Dd1/c1
   !Dc1 = Dc1 - ( (RefPt(1) - node(1,1))/c1**2  ) * Dd1

   ! c2 = node(4,1) - node(3,1)
   !node_ad(4,1) = node_ad(4,1) + Dc2
   !node_ad(3,1) = node_ad(3,1) - Dc2

   ! c1 = node(2,1) - node(1,1)
   !node_ad(2,1) = node_ad(2,1) + Dc1
   !node_ad(1,1) = node_ad(1,1) - Dc1

end if

BiLinearDiff2DVal_ad = 0.0_wp

if (present(Nnode)) then
   if (Nnode == 3) then
      !if (abs(lon(1)-lon(2)) < 1.0E-6_wp) then
      if (abs(node(1,1)-node(2,1)) < 1.0E-6_wp) then
         ! change 2 and 3 to avoid singularity, copy 3 to 4
         node_ad(4,:) = node_ad(2,:)
         node_ad(2,:) = node_ad(3,:)
         node_ad(3,:) = node_ad(4,:)
      !else if (abs(lon(2)-lon(3)) < 1.0E-6_wp) then
      else if (abs(node(2,1)-node(3,1)) < 1.0E-6_wp) then
         ! copy 1 to 4 to avoid singularity
         node_ad(1,:) = node_ad(4,:)
      else
         node_ad(2,:) = node_ad(4,:)
      end if
   end if
end if

end subroutine BiLinearDiff2D_ad


!---------------------------------------------------------------------
! function BiLinearDiff2D_full_ad
!---------------------------------------------------------------------
!
! Start of RoboDoc comments >>>>
!****F* Interpolation/BiLinearDiff2D_full_ad
!
! Name
! BiLinearDiff2D_full_ad
!
! Call
! val =  BiLinearDiff2D_full_ad(RefPt, node)
!
! Purpose
! 2D interpolation usig a bilinear interpoltion between four reference
! values. Any set of nodes is possible as long as (y2-y1) <> 0 and
! (x2-x1) <> 0 and (x4-x3) <> 0, i. e. 3 or 4 points are on a straoght line.
! If only 3 reference points ara available, the routine can be called with
! P4 = P2 (node 4 = node 2).
!
! In case of 3 reference points there is a singularity at  x = x2
! i.e. RefPt(1) = x2, which is replaced by the mean value
! od two points close to the singularity which can be evaluated.
!
!
!    y ^
!      |
!      |   node3                  node4             y3 > y4, x3 < x4
!      |
!      |          RfPt
!      |
!      |   node1                  node2             y1 > y2, x1 < x2
!      |
!      ---------------------------------------> x
!
!          x1 < x2                x2 > x1
!          y1 < y3                y2 < y4
!
! Test program: .../gpstomo/trunk/Tomo/TomoTest/Interpol/TestBiLinear2D.f90
!
!

! Input
! RefPt    - reference point, the field shall be interpolated at this
!            point, RefPoint(1) - x, longitude
!                   RefPoint(2) - y, latitude
! node     - array
!            node(i,:), i=1, ... , 4 : 4 nodes at the cell corner
!            node(i,j), j=1, 2, 3    : j = 1 - x
!                                      j = 2 - y
!                                      j = 3 - f_j = f(x,y)
!
! Output
! BiLinearDiff2D (return value) - interpolated value
!
! External References
! None
!
! Record of revisions
! Date        Programmer   Description of change
! ====        ==========   =====================
! 17.05.2013  M. Bender    new
!---------------------------------------------------------------------
! Source
!---------------------------------------------------------------------
subroutine BiLinearDiff2D_full_ad(RefPt, nodein, RefPt_ad, node_ad,   &
                             BiLinearDiff2DVal_ad, Nnode)

! List of calling arguments:
real (wp), intent(inout)                     :: BiLinearDiff2DVal_ad
real (wp), dimension(1:2), intent(in)        :: RefPt
real (wp), dimension(1:4,1:3), intent(in)    :: nodein
real (wp), dimension(1:2), intent(out)        :: RefPt_ad
real (wp), dimension(1:4,1:3), intent(out)   :: node_ad
integer,  intent(in), optional               :: Nnode

! List of local variables:
real (wp), dimension(1:4,1:3) :: node
!real (wp) :: BiLinear2D
real (wp) :: c1, c2, d1, d2, d3, d4
real (wp) :: ya, yb, fa, fb, dy      !, a1, a2, a3
real (wp) :: Dc1, Dc2, Dd1, Dd2, Dd3, Dd4
real (wp) :: Dya, Dyb, Dfa, Dfb, Ddy !, Da1, Da2, Da3
!real (wp) :: epsilon, I1, I2, DI1, DI2
!---------------------------------------------------------------------

node = nodein

if (present(Nnode)) then
   if (Nnode == 3) then
      !if (abs(lon(1)-lon(2)) < 1.0E-6_wp) then
      if (abs(node(1,1)-node(2,1)) < 1.0E-6_wp) then
         ! change 2 and 3 to avoid singularity, copy 3 to 4
         node(4,:) = node(3,:)
         node(3,:) = node(2,:)
         node(2,:) = node(4,:)
      !else if (abs(lon(2)-lon(3)) < 1.0E-6_wp) then
      else if (abs(node(2,1)-node(3,1)) < 1.0E-6_wp) then
         ! copy 1 to 4 to avoid singularity
         node(4,:) = node(1,:)
      else
         node(4,:) = node(2,:)
      end if
   end if
end if

node_ad = 0.0_wp

! c1 = x2 - x1
c1 = node(2,1) - node(1,1)
! c2 = x4 - x3
c2 = node(4,1) - node(3,1)
! d1 = (x - x1)/(x2-x1) =  (x - x1)/c1
d1 = (RefPt(1) - node(1,1))/c1
! d2 = (x2 - x)/(x2-x1) =  (x2 - x)/c1
d2 = (node(2,1) - RefPt(1))/c1
! d3 = (x - x3)/(x4-x3) =  (x - x3)/c2
d3 = (RefPt(1) - node(3,1))/c2
! d4 = (x4 - x)/(x4-x3) =  (x4 - x)/c2
d4 = (node(4,1) - RefPt(1))/c2
!write(*,*) 'BiLinear2D> c1, c2, d1, d2, d3, d4 : ', c1, c2, d1, d2, d3, d4

fa = d2*node(1,3) + d1*node(2,3)
fb = d4*node(3,3) + d3*node(4,3)
ya = d2*node(1,2) + d1*node(2,2)
yb = d4*node(3,2) + d3*node(4,2)
dy = yb - ya

!!$if (abs(dy) .lt. 1.0E-7_wp .and. abs(fb-fa) .lt. 1.0E-7_wp) then
!!$   ! x = x2, singularity, should happen only for 3 point interpolation
!!$   ! The singularity is replaced by the mean value of two points close
!!$   ! to the singularity, i.e. computed for RefPt(1) + epsilon and
!!$   ! RefPt(1) - epsilon.
!!$
!!$   ! Assuming latitude and longitude given in radian, a difference
!!$   ! of 9*10^-5 rad is equivalent to a distance of ~ 10 m,
!!$   ! i. e. epsilon = 1.0D-5 should be very close to the singularity
!!$   ! and the error of the interpolation is < 0.1 %
!!$   epsilon = 1.0E-5_wp
!!$
!!$   ! Recompute second interpolation
!!$   ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!!$   ! x - epsilon (RefPt(1) - epsilon)
!!$   ! d1 = (x - x1)/(x2-x1) =  (x - x1)/c1
!!$   d1 = (RefPt(1)-epsilon - node(1,1))/c1
!!$   ! d2 = (x2 - x)/(x2-x1) =  (x2 - x)/c1
!!$   d2 = (node(2,1) - RefPt(1)+epsilon)/c1
!!$   ! d3 = (x - x3)/(x4-x3) =  (x - x3)/c2
!!$   d3 = (RefPt(1)-epsilon - node(3,1))/c2
!!$   ! d4 = (x4 - x)/(x4-x3) =  (x4 - x)/c2
!!$   d4 = (node(4,1) - RefPt(1)+epsilon)/c2
!!$
!!$   fa = d2*node(1,3) + d1*node(2,3)
!!$   fb = d4*node(3,3) + d3*node(4,3)
!!$   ya = d2*node(1,2) + d1*node(2,2)
!!$   yb = d4*node(3,2) + d3*node(4,2)
!!$   dy = yb - ya
!!$   ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!!$
!!$   ! BiLinearDiff2D = 0.50D0 * (I1 + I2)
!!$   DI1 =  0.50_wp * BiLinearDiff2DVal_ad
!!$   DI2 =  0.50_wp * BiLinearDiff2DVal_ad
!!$
!!$   ! I2 = fa*(yb- RefPt(2))/dy + fb*(RefPt(2)-ya)/dy
!!$   Dfa = ( (yb-RefPt(2))/dy ) * DI2
!!$   Dfb = ( (RefPt(2)-ya)/dy ) * DI2
!!$   Dyb = ( fa/dy ) * DI2
!!$   Dya = -( fb/dy ) * DI2
!!$   RefPt_ad(2) = ( (fb-fa)/dy ) * DI2
!!$   Ddy =  -( fa*(yb- RefPt(2))/dy**2 + fb*(RefPt(2)-ya)/dy**2 ) * DI2
!!$
!!$  ! dy = yb - ya
!!$   Dyb = Dyb + Ddy
!!$   Dya = Dya - Ddy
!!$
!!$   ! yb = d4*node(3,2) + d3*node(4,2)
!!$   Dd4 = node(3,2) * Dyb
!!$   Dd3 = node(4,2) * Dyb
!!$   node_ad(3,2) = d4 * Dyb
!!$   node_ad(4,2) = d3 * Dyb
!!$
!!$   ! ya = d2*node(1,2) + d1*node(2,2)
!!$   Dd2 = node(1,2) * Dya
!!$   Dd1 = node(2,2) * Dya
!!$   node_ad(1,2) = d2 * Dya
!!$   node_ad(2,2) = d1 * Dya
!!$
!!$   ! fb = d4*node(3,3) + d3*node(4,3)
!!$   Dd4 = Dd4 + node(3,3) * Dfb
!!$   Dd3 = Dd3 + node(4,3) * Dfb
!!$   node_ad(3,3) = d4 * Dfb
!!$   node_ad(4,3) = d3 * Dfb
!!$
!!$   ! fa = d2*node(1,3) + d1*node(2,3)
!!$   Dd2 = Dd2 + node(1,3) * Dfa
!!$   Dd1 = Dd1 + node(2,3) * Dfa
!!$   node_ad(1,3) = d2 * Dfa
!!$   node_ad(2,3) = d1 * Dfa
!!$
!!$   ! d4 = (node(4,1) - RefPt(1)+epsilon)/c2
!!$   node_ad(4,1) = Dd4/c2
!!$   RefPt_ad(1) = -Dd4/c2
!!$   Dc2 = -( (node(4,1) - RefPt(1) + epsilon)/c2**2 ) * Dd4
!!$
!!$   ! d3 = (RefPt(1)-epsilon - node(3,1))/c2
!!$   RefPt_ad(1) = RefPt_ad(1) + Dd3/c2
!!$   node_ad(3,1) = -Dd3/c2
!!$   Dc2 = Dc2 - ( (RefPt(1) - epsilon - node(3,1))/c2**2 ) * Dd3
!!$
!!$   ! d2 = (node(2,1) - RefPt(1)+epsilon)/c1
!!$   node_ad(2,1) = Dd2/c1
!!$   RefPt_ad(1) = RefPt_ad(1) - Dd2/c1
!!$   Dc1 = -( (node(2,1) - RefPt(1) + epsilon)/c1**2  ) *Dd2
!!$
!!$   !  d1 = (RefPt(1)-epsilon - node(1,1))/c1
!!$   RefPt_ad(1) = RefPt_ad(1) + Dd1/c1
!!$   node_ad(1,1) = -Dd1/c1
!!$   Dc1 = Dc1 - ( (RefPt(1) - epsilon - node(1,1))/c1**2  ) * Dd1
!!$
!!$   ! Recompute first interpolation
!!$   ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!!$   ! x + epsilon (RefPt(1) + epsilon)
!!$   ! d1 = (x - x1)/(x2-x1) =  (x - x1)/c1
!!$   d1 = (RefPt(1)+epsilon - node(1,1))/c1
!!$   ! d2 = (x2 - x)/(x2-x1) =  (x2 - x)/c1
!!$   d2 = (node(2,1) - RefPt(1)-epsilon)/c1
!!$   ! d3 = (x - x3)/(x4-x3) =  (x - x3)/c2
!!$   d3 = (RefPt(1)+epsilon - node(3,1))/c2
!!$   ! d4 = (x4 - x)/(x4-x3) =  (x4 - x)/c2
!!$   d4 = (node(4,1) - RefPt(1)-epsilon)/c2
!!$
!!$   fa = d2*node(1,3) + d1*node(2,3)
!!$   fb = d4*node(3,3) + d3*node(4,3)
!!$   ya = d2*node(1,2) + d1*node(2,2)
!!$   yb = d4*node(3,2) + d3*node(4,2)
!!$   dy = yb - ya
!!$   ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!!$
!!$   ! I1 = fa*(yb- RefPt(2))/dy + fb*(RefPt(2)-ya)/dy
!!$   Dfa = ( (yb-RefPt(2))/dy ) * DI1
!!$   Dfb = ( (RefPt(2)-ya)/dy ) * DI1
!!$   Dyb = ( fa/dy ) * DI1
!!$   Dya = -( fb/dy ) * DI1
!!$   RefPt_ad(2) = RefPt_ad(2) + ( (fb-fa)/dy ) * DI1
!!$   Ddy =  -( fa*(yb- RefPt(2))/dy**2 + fb*(RefPt(2)-ya)/dy**2 ) * DI1
!!$
!!$  ! dy = yb - ya
!!$   Dyb = Dyb + Ddy
!!$   Dya = Dya - Ddy
!!$
!!$   ! yb = d4*node(3,2) + d3*node(4,2)
!!$   Dd4 = node(3,2) * Dyb
!!$   Dd3 = node(4,2) * Dyb
!!$   node_ad(3,2) = node_ad(3,2) + d4 * Dyb
!!$   node_ad(4,2) = node_ad(4,2) + d3 * Dyb
!!$
!!$   ! ya = d2*node(1,2) + d1*node(2,2)
!!$   Dd2 = node(1,2) * Dya
!!$   Dd1 = node(2,2) * Dya
!!$   node_ad(1,2) = node_ad(1,2) + d2 * Dya
!!$   node_ad(2,2) = node_ad(2,2) + d1 * Dya
!!$
!!$   ! fb = d4*node(3,3) + d3*node(4,3)
!!$   Dd4 = Dd4 + node(3,3) * Dfb
!!$   Dd3 = Dd3 + node(4,3) * Dfb
!!$   node_ad(3,3) = node_ad(3,3) + d4 * Dfb
!!$   node_ad(4,3) = node_ad(4,3) + d3 * Dfb
!!$
!!$   ! fa = d2*node(1,3) + d1*node(2,3)
!!$   Dd2 = Dd2 + node(1,3) * Dfa
!!$   Dd1 = Dd1 + node(2,3) * Dfa
!!$   node_ad(1,3) = node_ad(1,3) + d2 * Dfa
!!$   node_ad(2,3) = node_ad(2,3) + d1 * Dfa
!!$
!!$   ! d4 = (node(4,1) - RefPt(1)-epsilon)/c2
!!$   node_ad(4,1) = node_ad(4,1) + Dd4/c2
!!$   RefPt_ad(1) = RefPt_ad(1) - Dd4/c2
!!$   Dc2 = -( (node(4,1) - RefPt(1) - epsilon)/c2**2 ) * Dd4
!!$
!!$   ! d3 = (RefPt(1)+epsilon - node(3,1))/c2
!!$   RefPt_ad(1) = RefPt_ad(1) + Dd3/c2
!!$   node_ad(3,1) = node_ad(3,1) - Dd3/c2
!!$   Dc2 = Dc2 - ( (RefPt(1) + epsilon - node(3,1))/c2**2 ) * Dd3
!!$
!!$   ! d2 = (node(2,1) - RefPt(1)-epsilon)/c1
!!$   node_ad(2,1) =  node_ad(2,1) + Dd2/c1
!!$   RefPt_ad(1) = RefPt_ad(1) - Dd2/c1
!!$   Dc1 = -( (node(2,1) - RefPt(1) - epsilon)/c1**2  ) *Dd2
!!$
!!$   ! d1 = d1 = (RefPt(1)+epsilon - node(1,1))/c1
!!$   RefPt_ad(1) = RefPt_ad(1) + Dd1/c1
!!$   node_ad(1,1) = node_ad(1,1) - Dd1/c1
!!$   Dc1 = Dc1 - ( (RefPt(1) + epsilon - node(1,1))/c1**2  ) * Dd1
!!$
!!$   ! c2 = node(4,1) - node(3,1)
!!$   node_ad(4,1) = node_ad(4,1) + Dc2
!!$   node_ad(3,1) = node_ad(3,1) - Dc2
!!$
!!$   ! c1 = node(2,1) - node(1,1)
!!$   node_ad(2,1) = node_ad(2,1) + Dc1
!!$   node_ad(1,1) = node_ad(1,1) - Dc1
if (abs(dy) .lt. 1.0E-8_wp) then
   ! slant%Nnghb = 3 : reference point is very close to point 2
   node_ad(2,3) = BiLinearDiff2DVal_ad
else
   ! No singularity, interpolate ...

   ! BiLinear2D = fa*(yb- RefPt(2))/dy + fb*(RefPt(2)-ya)/dy
   Dfa = ( (yb-RefPt(2))/dy ) * BiLinearDiff2DVal_ad
   Dfb = ( (RefPt(2)-ya)/dy ) * BiLinearDiff2DVal_ad
   Dyb = ( fa/dy ) * BiLinearDiff2DVal_ad
   Dya = -( fb/dy ) * BiLinearDiff2DVal_ad
   RefPt_ad(2) = ( (fb-fa)/dy ) * BiLinearDiff2DVal_ad
   Ddy =  -( fa*(yb- RefPt(2))/dy**2 + fb*(RefPt(2)-ya)/dy**2 ) *   &
         BiLinearDiff2DVal_ad

   ! dy = yb - ya
   Dyb = Dyb + Ddy
   Dya = Dya - Ddy

   ! yb = d4*node(3,2) + d3*node(4,2)
   Dd4 = node(3,2) * Dyb
   Dd3 = node(4,2) * Dyb
   node_ad(3,2) = d4 * Dyb
   node_ad(4,2) = d3 * Dyb

   !write(*,*) 'BiLinear2D_ad> BiLinear2DVal_ad = ', BiLinear2DVal_ad
   !write(*,*) 'BiLinear2D_ad> Dyb              = ', Dyb


   ! ya = d2*node(1,2) + d1*node(2,2)
   Dd2 = node(1,2) * Dya
   Dd1 = node(2,2) * Dya
   node_ad(1,2) = d2 * Dya
   node_ad(2,2) = d1 * Dya

   ! fb = d4*node(3,3) + d3*node(4,3)
   Dd4 = Dd4 + node(3,3) * Dfb
   Dd3 = Dd3 + node(4,3) * Dfb
   node_ad(3,3) = d4 * Dfb
   node_ad(4,3) = d3 * Dfb

   ! fa = d2*node(1,3) + d1*node(2,3)
   Dd2 = Dd2 + node(1,3) * Dfa
   Dd1 = Dd1 + node(2,3) * Dfa
   node_ad(1,3) = d2 * Dfa
   node_ad(2,3) = d1 * Dfa

   ! d4 = (node(4,1) - RefPt(1))/c2
   node_ad(4,1) = Dd4/c2
   RefPt_ad(1) = -Dd4/c2
   Dc2 = -( (node(4,1) - RefPt(1))/c2**2 ) * Dd4

   ! d3 = (RefPt(1) - node(3,1))/c2
   RefPt_ad(1) = RefPt_ad(1) + Dd3/c2
   node_ad(3,1) = -Dd3/c2
   Dc2 = Dc2 - ( (RefPt(1) - node(3,1))/c2**2 ) * Dd3

   ! d2 = (node(2,1) - RefPt(1))/c1
   node_ad(2,1) = Dd2/c1
   RefPt_ad(1) = RefPt_ad(1) - Dd2/c1
   Dc1 = -( (node(2,1) - RefPt(1))/c1**2  ) *Dd2

   ! d1 = (RefPt(1) - node(1,1))/c1
   RefPt_ad(1) = RefPt_ad(1) + Dd1/c1
   node_ad(1,1) = -Dd1/c1
   Dc1 = Dc1 - ( (RefPt(1) - node(1,1))/c1**2  ) * Dd1

   ! c2 = node(4,1) - node(3,1)
   node_ad(4,1) = node_ad(4,1) + Dc2
   node_ad(3,1) = node_ad(3,1) - Dc2

   ! c1 = node(2,1) - node(1,1)
   node_ad(2,1) = node_ad(2,1) + Dc1
   node_ad(1,1) = node_ad(1,1) - Dc1

end if

BiLinearDiff2DVal_ad = 0.0_wp

if (present(Nnode)) then
   if (Nnode == 3) then
      !if (abs(lon(1)-lon(2)) < 1.0E-6_wp) then
      if (abs(node(1,1)-node(2,1)) < 1.0E-6_wp) then
         ! change 2 and 3 to avoid singularity, copy 3 to 4
         node_ad(4,:) = node_ad(2,:)
         node_ad(2,:) = node_ad(3,:)
         node_ad(3,:) = node_ad(4,:)
      !else if (abs(lon(2)-lon(3)) < 1.0E-6_wp) then
      else if (abs(node(2,1)-node(3,1)) < 1.0E-6_wp) then
         ! copy 1 to 4 to avoid singularity
         node_ad(1,:) = node_ad(4,:)
      else
         node_ad(2,:) = node_ad(4,:)
      end if
   end if
end if

end subroutine BiLinearDiff2D_full_ad


function NWein_tl( P, T, E, P_tl, T_tl, E_tl )

real (wp)   :: NWein_tl

! List of calling arguments:
real (wp),  intent(in)   ::  T, P, E
real (wp),  intent(in)   ::  T_tl, P_tl, E_tl

! List of local variables:
real (wp)  ::  hP, hE
real (wp)  ::  DhP, DhE
!---------------------------------------------------------------------

! Umrechnen der Druecke in hPa:
hP = P/100.0
hE = E/100.0

DhP = P_tl/100.0_wp
DhE = E_tl/100.0_wp

!NWein = k1*(hP/T) + (k2 - k1)*(hE/T) + k3*(hE/(T*T))

NWein_tl = DhP * k1/T + DhE * ( (k2-k1)/T + k3/T**2 ) -  &
           T_tl * ( k1*(hP/T**2) + (k2 - k1)*(hE/T**2) + 2.0*k3*(hE/(T**3)) )

end function NWein_tl


subroutine NWein_ad( P, T, E, P_ad, T_ad, E_ad, N_ad )

!real (wp)   :: NWein

! List of calling arguments:
real (wp),  intent(in)    ::  T, P, E
real (wp),  intent(inout) ::  N_ad
real (wp),  intent(out)   ::  T_ad, P_ad, E_ad

! List of local variables:
real (wp)  ::  hP, hE
real (wp)  ::  DhP, DhE
!---------------------------------------------------------------------

! Umrechnen der Druecke in hPa:
hP = P/100.0_wp
hE = E/100.0_wp

! NWein = k1*(hP/T) + (k2 - k1)*(hE/T) + k3*(hE/(T*T))
DhP  = N_ad * k1/T
DhE  = N_ad * ( (k2 - k1)/T + k3/T**2 )
T_ad = -N_ad * ( k1*(hP/T**2) + (k2 - k1)*(hE/T**2) + 2.0*k3*(hE/(T**3)))

! hE = E/100.0
E_ad = DhE / 100.0_wp

! hP = P/100.0
P_ad = DhP / 100.0_wp

N_ad = 0.0_wp

end subroutine NWein_ad

function LayerSearch (a, value)

! List of calling arguments:
integer                  :: LayerSearch
real (wp), pointer            :: a(:)
real (wp), intent(in)         :: value

! List of local variables:
integer :: i, j, k
!---------------------------------------------------------------------

i = lbound(a,1)
j = ubound(a,1)

do while (i+1 .lt. j)
   k=(i+j)/2
   if (a(k) .le. value) then
      j = k
   else
      i = k
   end if
   !write(*,*) i, j, k
end do

LayerSearch = i+1

end function LayerSearch


subroutine SatPos (StationCoord, azimuth, elevation, GNSS, SatelliteCoord)

! List of calling arguments:
real (wp), dimension(1:3), intent(in)  :: StationCoord
real (wp), intent(in)                  :: azimuth, elevation
character (len=*),  intent(in)         :: GNSS
real (wp), dimension(1:3), intent(out) :: SatelliteCoord

! List of local variables:
real (wp)                 ::  px, py, pz,  gx, gy, gz
real (wp), dimension(1:3) :: CoordCart
type (Line)               :: ray
real (wp)                 :: a, b, c, lambda1, lambda2
logical                   :: IsReal
!---------------------------------------------------------------------

! Compute kartesian ECEF coordinates from the ellipsoidal station coordinates
call  Ellips2Cart( StationCoord(1),    & ! in: longitude
                   StationCoord(2),    & ! in: latitude
                   StationCoord(3),    & ! in: altitude
                   CoordCart(1),       & ! out: X
                   CoordCart(2),       & ! out: Y
                   CoordCart(3),       & ! out: Z
                   WGS84Param      )     ! in: reference ellipsoid

! Compute coordinates of an arbitrary point on the signal path in
! the local horizon system
call Azimut2Horz(azimuth, pi05-elevation, 5.0D4, px, py, pz)

! Transform these coordinates to global kartesian ECEF coordinates
call LocalHorz2Cart(px, py, pz,                 &
                    StationCoord(1),   &
                    StationCoord(2),   &
                    CoordCart(1),  &
                    CoordCart(2),  &
                    CoordCart(3),  &
                    gx, gy, gz)

! Define a straight line through the station and that point,
! i.e. an unit vector pointing from the station to the satellite
call LineDef(ray, CoordCart, (/gx,gy,gz/), 3)

! >>>>> Hier muessen die verschiedenen GNSS abgefragt wrden und der
!       jeweils richtige Radius eingesetzt werden !!!!
a = dot_product(CoordCart,CoordCart) - GPSradius**2

b = 2.0_wp * dot_product(CoordCart,ray%unitvec)
c = dot_product(ray%unitvec,ray%unitvec)

! Solve the quadratic equation for lambda
call QuadEq(b/c, a/c, lambda1, lambda2, IsReal)

!write(*,*) ' lambda1, lambda2 ',  lambda1, lambda2

if (IsReal) then
   ! There is a real solution to the quadratic equation
   ! Usually, there are two real solutions: The positive solution
   ! defines the satellite position above the station, the negative
   ! solution a satellite on the other side ofthe Earth, which cannot be
   ! seen at the given station.

   if (lambda1 .gt. 0.0) then
      SatelliteCoord = CoordCart + lambda1*ray%unitvec
   else if (lambda2 .gt. 0.0) then
      SatelliteCoord = CoordCart + lambda2*ray%unitvec
   else
      write(*,*) 'SatPos> Both lambdas are negative, the satellite'
      write(*,*) 'SatPos> would be on the other side of the Earth'
   end if

else
   write(*,*) 'SatPos> There is no real solution and the satellite position'
   write(*,*) 'SatPos> could not be computed'
end if

end subroutine SatPos


!---------------------------------------------------------------------
! subroutine Slant2Ellips
!---------------------------------------------------------------------
!
! Start of RoboDoc comments >>>>
!****S* Coord/Slant2Ellips
!
! Name
! Slant2Ellips
!
! Call
! call Slant2Ellips(Xslant, Xellips, azi, elev, RefPtEllips, RefPtECEF, Ellips)
!
! Purpose
! Transforms the coordinates of a vector from the slant system into
! ellipsoidal geographic coordinates like WGS84.
! Inverse transformation: Ellips2Slant
!
! This is a driver routine for three consecutive transformations:
! Slant2LocalHorz => LocalHorz2Cart => Cart2Ellips
!
! Ellipsoidal system: ellipsoidal geographic coordinates
!
! Slant system: right handed  cartesian system with its origin
!                          at the position (lambda,phi).
!               The x-axis points to the satellite with the given azimuth
!               and elevation, the y-axis is in the tangent plane and the
!               z-axis is the vertical axis if the elevation is zero.
!               For epsilon <> 0, the verticl component is mapped on the
!               z- and x-axis.
!
! Input
! Xslant - Coordinates of the slant system
! azi    - azimuth (angle in rad)
!          A = [0, ..., 2*pi] or
!          A = [-pi, ..., pi]
! elev   - elevation (angle in rad)
!          elev = [0, ..., pi/2]
! RefPtEllips - coordinates of reference point, e.g. station,
!               ellipsoidal geographic coordinates (lon, lat, height)
! RefPtECEF   - coordinates of reference point, e.g. station,
!               earth centered earth fixed cartesian coordinates (ECEF)
! Ellips      - reference ellipsoid, e.g. WGS84
!
! Output
! Xellips - geographic coordinates (longitude, latitude, height)
!
! External References
! None
!
! Record of revisions
! Date        Programmer   Description of change
! ====        ==========   =====================
! 18.07.2013  M. Bender    new
!---------------------------------------------------------------------
! Source
!---------------------------------------------------------------------
subroutine Slant2Ellips (Xslant, Xellips,                              &
                         azi, elev, RefPtEllips, RefPtECEF, Ellips)

! List of calling arguments:
real (wp), dimension(1:3), intent(in)  :: Xslant
real (wp), dimension(1:3), intent(out) :: Xellips
real (wp), intent(in)                  :: azi, elev
real (wp), dimension(1:3), intent(in)  :: RefPtECEF, RefPtEllips
type (EllipsParam), intent(in)         :: Ellips

! List of local variables:
real (wp), dimension(1:3) :: Xhorz, Xecef
!---------------------------------------------------------------------

! Transform coordinates from the slant system to the local horizon system
call Slant2LocalHorz (Xslant, azi, elev, Xhorz)

! ... local horizon system to earth centered earth fixed coordinates (ECEF)
call LocalHorz2Cart(Xhorz(1), Xhorz(2), Xhorz(3),               &
                    RefPtEllips(1), RefPtEllips(2),             &
                    RefPtECEF(1), RefPtECEF(2), RefPtECEF(3),   &
                    Xecef(1), Xecef(2), Xecef(3)              )

! ... ECEF to ellipsoidal coordinates
call Cart2Ellips(Xecef(1), Xecef(2), Xecef(3),    &
                 Xellips(1), Xellips(2), Xellips(3), Ellips)

end subroutine Slant2Ellips


!---------------------------------------------------------------------
! subroutine Slant2LocalHorz
!---------------------------------------------------------------------
!
! Start of RoboDoc comments >>>>
!****S* Coord/Slant2LocalHorz
!
! Name
! Slant2LocalHorz
!
! Call
! call Slant2LocalHorz (Xslant, A, elev, Xhorz)
!
! Purpose
! Transforms the coordinates of a vector from the slant system into
! the local horizon system.
! Inverse transformation: LocalHorz2Slant
!
! Local horizontal system: left handed cartesian system with its origin
!                          at the position (lambda,phi).
! Slant system: right handed  cartesian system with its origin
!                          at the position (lambda,phi).
!               The x-axis points to the satellite with the given azimuth
!               and elevation, the y-axis is in the tangent plane and the
!               z-axis is the vertical axis if the elevation is zero.
!               For epsilon <> 0, the verticl component is mapped on the
!               z- and x-axis.
!
!
! Input
! Xslant - Coordinates of the slant system
! A      - azimuth of the object (angle in rad)
!          A = [0, ..., 2*pi] or
!          A = [-pi, ..., pi]
! elev   - elevation (angle in rad)
!          elev = [0, ..., pi/2]
!
! Output
! Xhorz - Coordinates of the local horizon system
!
! External References
! None
!
! Examples
!
! Record of revisions
! Date        Programmer   Description of change
! ====        ==========   =====================
! 16.07.2013  M. Bender    neu
!---------------------------------------------------------------------
! Source
!---------------------------------------------------------------------
subroutine Slant2LocalHorz (Xslant, A, elev, Xhorz)

! List of calling arguments:
real (wp), dimension(1:3), intent(in)  :: Xslant
real (wp), intent(in)                  :: A, elev
real (wp), dimension(1:3), intent(out) :: Xhorz

! List of local variables:
!---------------------------------------------------------------------

Xhorz(1) = Xslant(1)*cos(elev)*cos(A) +    &
           Xslant(2)*sin(A)           -    &
           Xslant(3)*sin(elev)*cos(A)

Xhorz(2) =  Xslant(1)*cos(elev)*sin(A) -    &
            Xslant(2)*cos(A)           -    &
            Xslant(3)*sin(elev)*sin(A)

Xhorz(3) =  Xslant(1)*sin(elev) + Xslant(3)*cos(elev)

end subroutine Slant2LocalHorz


!---------------------------------------------------------------------
! subroutine Slant2Ellips_tl
!---------------------------------------------------------------------
!
! Start of RoboDoc comments >>>>
!****S* Coord/Slant2Ellips_tl
!
! Name
! Slant2Ellips_tl
!
! Call
! call Slant2Ellips_tl
!  (Xslant, Xellips, azi, elev, RefPtEllips, RefPtECEF, Ellips)
!
! Purpose
! Transforms the coordinates of a vector from the slant system into
! ellipsoidal geographic coordinates like WGS84.
! Inverse transformation: Ellips2Slant
!
! This is a driver routine for three consecutive transformations:
! Slant2LocalHorz => LocalHorz2Cart => Cart2Ellips
!
! Ellipsoidal system: ellipsoidal geographic coordinates
!
! Slant system: right handed  cartesian system with its origin
!                          at the position (lambda,phi).
!               The x-axis points to the satellite with the given azimuth
!               and elevation, the y-axis is in the tangent plane and the
!               z-axis is the vertical axis if the elevation is zero.
!               For epsilon <> 0, the verticl component is mapped on the
!               z- and x-axis.
!
! Input
! Xslant - Coordinates of the slant system
! azi    - azimuth (angle in rad)
!          A = [0, ..., 2*pi] or
!          A = [-pi, ..., pi]
! elev   - elevation (angle in rad)
!          elev = [0, ..., pi/2]
! RefPtEllips - coordinates of reference point, e.g. station,
!               ellipsoidal geographic coordinates (lon, lat, height)
! RefPtECEF   - coordinates of reference point, e.g. station,
!               earth centered earth fixed cartesian coordinates (ECEF)
! Ellips      - reference ellipsoid, e.g. WGS84
!
! Output
! Xellips - geographic coordinates (longitude, latitude, height)
!
! External References
! None
!
! Record of revisions
! Date        Programmer   Description of change
! ====        ==========   =====================
! 18.09.2013  M. Bender    new
!---------------------------------------------------------------------
! Source
!---------------------------------------------------------------------
subroutine Slant2Ellips_tl (Xslant, Xellips, Xslant_tl, Xellips_tl,     &
                            azi, elev, RefPtEllips, RefPtECEF, Ellips)

! List of calling arguments:
real (wp), dimension(1:3), intent(in)  :: Xslant, Xslant_tl
real (wp), dimension(1:3), intent(out) :: Xellips, Xellips_tl
real (wp), intent(in)                  :: azi, elev
real (wp), dimension(1:3), intent(in)  :: RefPtECEF, RefPtEllips
type (EllipsParam), intent(in)                    :: Ellips

! List of local variables:
real (wp), dimension(1:3) :: Xhorz, Xecef
real (wp), dimension(1:3) :: Xhorz_tl, Xecef_tl
real (wp)                 :: lambda_tl, phi_tl
real (wp)                 :: ox_tl, oy_tl, oz_tl

!---------------------------------------------------------------------

! Transform coordinates from the slant system to the local horizon system
call Slant2LocalHorz_tl(Xslant, Xslant_tl, azi, elev, Xhorz, Xhorz_tl)

! ... local horizon system to earth centered earth fixed coordinates (ECEF)
lambda_tl = 0.0_wp
phi_tl    = 0.0_wp
ox_tl     = 0.0_wp
oy_tl     = 0.0_wp
oz_tl     = 0.0_wp
call LocalHorz2Cart_tl(Xhorz(1), Xhorz(2), Xhorz(3),               &
                       RefPtEllips(1), RefPtEllips(2),             &
                       RefPtECEF(1), RefPtECEF(2), RefPtECEF(3),   &
                       Xecef(1), Xecef(2), Xecef(3),               &
                       Xhorz_tl(1), Xhorz_tl(2), Xhorz_tl(3),      &
                       lambda_tl, phi_tl,                          &
                       ox_tl, oy_tl, oz_tl,                        &
                       Xecef_tl(1), Xecef_tl(2), Xecef_tl(3)     )

! ... ECEF to ellipsoidal coordinates
call Cart2Ellips_tl(Xecef(1), Xecef(2), Xecef(3),    &
                    Xellips(1), Xellips(2), Xellips(3), Ellips ,  &
                    Xecef_tl(1), Xecef_tl(2), Xecef_tl(3),        &
                    Xellips_tl(1), Xellips_tl(2), Xellips_tl(3)  )

end subroutine Slant2Ellips_tl


!---------------------------------------------------------------------
! subroutine Slant2LocalHorz_tl
!---------------------------------------------------------------------
!
! Start of RoboDoc comments >>>>
!****S* Coord/Slant2LocalHorz_tl
!
! Name
! Slant2LocalHorz_tl
!
! Call
! call Slant2LocalHorz_tl (Xslant, A, elev, Xhorz)
!
! Purpose
! Transforms the coordinates of a vector from the slant system into
! the local horizon system.
! Inverse transformation: LocalHorz2Slant
!
! Local horizontal system: left handed cartesian system with its origin
!                          at the position (lambda,phi).
! Slant system: right handed  cartesian system with its origin
!                          at the position (lambda,phi).
!               The x-axis points to the satellite with the given azimuth
!               and elevation, the y-axis is in the tangent plane and the
!               z-axis is the vertical axis if the elevation is zero.
!               For epsilon <> 0, the verticl component is mapped on the
!               z- and x-axis.
!
!
! Input
! Xslant - Coordinates of the slant system
! A      - azimuth of the object (angle in rad)
!          A = [0, ..., 2*pi] or
!          A = [-pi, ..., pi]
! elev   - elevation (angle in rad)
!          elev = [0, ..., pi/2]
!
! Output
! Xhorz - Coordinates of the local horizon system
!
! External References
! None
!
! Examples
!
! Record of revisions
! Date        Programmer   Description of change
! ====        ==========   =====================
! 26.07.2013  M. Bender    neu
!---------------------------------------------------------------------
! Source
!---------------------------------------------------------------------
subroutine Slant2LocalHorz_tl ( Xslant, Xslant_tl, A, elev,     &
                                Xhorz, Xhorz_tl )

! List of calling arguments:
real (wp), dimension(1:3), intent(in)  :: Xslant, Xslant_tl
real (wp), intent(in)                  :: A, elev
real (wp), dimension(1:3), intent(out) :: Xhorz, Xhorz_tl

! List of local variables:
!---------------------------------------------------------------------

Xhorz(1) = Xslant(1)*cos(elev)*cos(A) +    &
           Xslant(2)*sin(A)           -    &
           Xslant(3)*sin(elev)*cos(A)

Xhorz(2) = Xslant(1)*cos(elev)*sin(A) -    &
           Xslant(2)*cos(A)           -    &
           Xslant(3)*sin(elev)*sin(A)

Xhorz(3) =  Xslant(1)*sin(elev) + Xslant(3)*cos(elev)


Xhorz_tl(1) = Xslant_tl(1)*cos(elev)*cos(A) +    &
              Xslant_tl(2)*sin(A)           -    &
              Xslant_tl(3)*sin(elev)*cos(A)

Xhorz_tl(2) = Xslant_tl(1)*cos(elev)*sin(A) -    &
              Xslant_tl(2)*cos(A)           -    &
              Xslant_tl(3)*sin(elev)*sin(A)

Xhorz_tl(3) =  Xslant_tl(1)*sin(elev) + Xslant_tl(3)*cos(elev)

end subroutine Slant2LocalHorz_tl


!---------------------------------------------------------------------
! subroutine Cart2EllipsGrad
!---------------------------------------------------------------------
!
! Start of RoboDoc comments >>>>
!****S* Coord/Cart2EllipsGrad
!
! Name
! Cart2EllipsGrad
!
! Call
! call Cart2Ellips(X, Y, Z, lambda, beta, height, Ellips)
!
! Purpose
! Umrechnung von kartesischen Koordinaten in ellipsoidische Koordinaten
! Direkte Rechnung ohne Iterationsverfahren
!
! Quelle: B. Hofmann-Wellenhof, H. Lichtenegger, J. Collins
!         Global Positioning System - Theory and Practice
!         Springer, Wien, New York, 1993, Seite 232
!         Die Transformationsformel fuer die Laenge wurde leicht abgewandelt.
!
!         Hier fehlt die Sonderbehandlung von X=Y=Z=0, Winkel = 0
!         und die separate Behandlung der Oktanden, so dass wirklich
!         Laengen zwischen 0° und 360° und Breiten zwischen -90° und +90°
!         herauskommen. Dies wurde analog zur Transformation in sphaerische
!         Koordinaten ergaenzt.
!
! Zum Testen dieser Transformation steht das Programm "TestWgs.f90" in
! ../ToolsTest/Coord zur Verfuegung. Dort werden bekannte Stationskoordinaten
! sowohl in kartesischen als auch in WGS84 Koordinaten eingelesen und dann
! mit den Ergebnissen der Transformationen verglichen.
!

! Cart2Ellips_dv => Cart2EllipsGrad
! Cc_dv = (1, 0, 0) ;  (0, 1, 0) ; (0, 0, 1)
! Ce_dv = d/d lambda , d/ d beta, d/d height



! Input
! X, Y, Z - kartesische Koordinaten i m
! Ellips  - Struktur vom Typ "EllipsParam" mit den Daten des Referenzellipsoids
!
! Output
! lambda - ellipsoidische Laenge in rad
! beta   - ellipsoidische Breite in rad
! height - ellipsoidische Hoehe in m
!
! External References
!
! Examples
!
! Record of revisions
! Date         Programmer   Description of change
! ====         ==========   =====================
! 20.11.2013   M. Bender    neu
!---------------------------------------------------------------------
! Source
!---------------------------------------------------------------------
!subroutine Cart2EllipsGrad(Cc, Ce, Ellips, Cc_dv, Ce_dv)
subroutine Cart2EllipsGrad(Cc, Ce, Ellips, CeGrad)

! List of calling arguments:
real (wp), dimension(1:3), intent(in)  :: Cc
real (wp), dimension(1:3), intent(out) :: Ce
!real (wp), dimension(:,:), intent(in)  :: Cc_dv
!real (wp), dimension(:,:), intent(out) :: Ce_dv
real (wp), dimension(1:3,1:3), intent(out) :: CeGrad

type (EllipsParam) :: Ellips

! List of local variables:
real (wp) :: lambda, beta, height
real (wp) :: rho, phi, gamma, g2
real (wp) :: N, theta, h
!integer   :: i

real (wp), dimension(1:3) :: drho, dtheta, dphi, dg2, dgamma, dN, dh
real (wp) :: c1, c2, f1, f2
!integer   :: multi
!---------------------------------------------------------------------

! Multidirectional differentiation: number of directions = 3

! Einige Hilfsgroessen:
rho = sqrt( Cc(1)**2 + Cc(2)**2 )
!do i=1, multi
!   drho(i) = Cc_dv(i,1)*Cc(1) + Cc_dv(i,2)*Cc(2)
!end do
drho(1) = Cc(1)
drho(2) = Cc(2)
drho(3) = 0.0_wp
drho = drho / rho
!write(*,*) 'Gr drho = ', drho

c1 = (Cc(3)*Ellips%a)/(rho*Ellips%b)
if (rho .gt. 1.0E-6_wp) then
   theta = atan(c1)
   f1 = Ellips%a / (rho*Ellips%b*(1.0_wp+c1**2))
   f2 =  Cc(3) * Ellips%a /  (rho**2 * Ellips%b*(1.0_wp+c1**2))
   !dtheta = Z_tl * Ellips%a / (rho*Ellips%b*(1.0D0+c1**2)) -         &
   !         drho * Z * Ellips%a /  (rho**2 * Ellips%b*(1.0D0+c1**2))
   !do i=1, multi
   !   dtheta(i) = Cc_dv(i,3) * f1 - drho(i) * f2
   !end do
   dtheta(1) = -drho(1) * f2
   dtheta(2) = -drho(2) * f2
   dtheta(3) = f1
else
   theta = pi05
   dtheta = 0.0_wp
end if
!write(*,*) 'Gr dtheta = ', dtheta

! Die ellipsoidische Laenge lambda:
if (abs(rho) > reps) then
   phi = 2.0_wp*atan(Cc(2)/(abs(Cc(1))+rho))
   f1 = -sign(2.0_wp*Cc(2),Cc(1))
   f2 =  2.0_wp*(abs(Cc(1)) + rho)
   !dphi = ( -sign(2.0D0*Y,X) * X_tl + 2.0D0*(abs(X) + rho) * Y_tl   -    &
   !         2.0D0*Y * drho ) / ((abs(X)+rho)**2 + Y**2)
   !do i=1, multi
   !   dphi(i) = Cc_dv(i,1) * f1 +  Cc_dv(i,2) * f2 -   &
   !             2.0_wp * Cc(2) * drho(i)
   !end do
   dphi(1) = f1 - 2.0_wp * Cc(2) * drho(1)
   dphi(2) = f2 - 2.0_wp * Cc(2) * drho(2)
   dphi(3) = 0.0_wp
   dphi = dphi /  ((abs(Cc(1))+rho)**2 + Cc(2)**2)
end if
if (Cc(1) .ge. 0.0_wp) then
   if (abs(Cc(1)) .lt. 1.0E-8_wp .and. abs(Cc(2)) .lt. 1.0E-8_wp) then
      ! if X=0 and Y=0
      lambda = 0.0_wp
      !lambda_tl = 0.0D0
      CeGrad(:,1) = 0.0_wp
   else if (Cc(2) .ge. 0.0_wp) then
      ! X>0 and Y>0
      lambda = phi
      !lambda_tl = dphi
      CeGrad(:,1) = dphi
   else
      ! X>0 and Y<0
      lambda = phi + 2.0_wp*pi
      !lambda_tl = dphi
      CeGrad(:,1) = dphi
   end if
else  ! X < 0
   lambda = pi - phi
   !lambda_tl = -dphi
   CeGrad(:,1) = -dphi
end if
!write(*,*) 'Gr dphi = ', dphi

! Die ellipsoidische Breite beta:
g2 = rho - Ellips%EN1sup2*Ellips%a*(cos(theta))**3  ! Nenner von gamma
!dg2 = drho + dtheta *                                                      &
!             3.0D0 * Ellips%EN1sup2*Ellips%a * sin(theta) * cos(theta)**2
f1 = 3.0_wp * Ellips%EN1sup2*Ellips%a * sin(theta) * cos(theta)**2
dg2 = drho + f1*dtheta

if (abs(g2) > reps) then
   gamma = atan( (Cc(3) + Ellips%EN2sup2*Ellips%b*(sin(theta))**3) / g2 )
   c2 = 1.0_wp/(1.0_wp +((Cc(3)+Ellips%EN2sup2*Ellips%b*(sin(theta))**3)/g2)**2)
   f1 = (Cc(3)+Ellips%EN2sup2*Ellips%b*(sin(theta)**3))/g2**2
   f2 = 3.0_wp*Ellips%EN2sup2*Ellips%b*cos(theta) * (sin(theta)**2)/g2
   !dgamma = c2 * (Z_tl/g2 -                                                 &
   !               dg2 * (Z+Ellips%EN2sup2*Ellips%b*(sin(theta)**3))/g2**2 + &
   !               dtheta * 3*Ellips%EN2sup2*Ellips%b*cos(theta) *           &
   !                        (sin(theta)**2)/g2 )
   !do i=1, multi
   !   dgamma(i) = Cc_dv(i,3)/g2 - dg2(i)*f1 + dtheta(i)*f2
   !end do
   dgamma(1) = -dg2(1)*f1 + dtheta(1)*f2
   dgamma(2) = -dg2(2)*f1 + dtheta(2)*f2
   dgamma(3) = 1.0_wp/g2 - dg2(3)*f1 + dtheta(3)*f2

   dgamma = dgamma * c2
! else ???? +pi/2 oder -pi/2 ???
end if
if (abs(rho) .gt. 1.0E-8_wp) then
   ! rho > 0
   beta = gamma
   !beta_tl = dgamma
   CeGrad(:,2) = dgamma
else  ! rho = 0
   if (Cc(3) .gt. 0.0_wp) then
      beta = 0.50_wp*pi
      !beta_tl = 0.0D0
      CeGrad(:,2) = 0.0_wp
   else if (abs(Cc(3)) < reps) then
      beta = 0.0_wp
      !beta_tl = 1.0D0
      CeGrad(:,2) = 1.0_wp
   else if (Cc(3) .lt. 0.0_wp) then
      beta = -0.5_wp*pi
      !beta_tl = 0.0_wp
      CeGrad(:,2) = 0.0_wp
   end if
end if
!write(*,*) 'Gr dgamma = ', dgamma

! Bestimmung des Querkruemmungsradius fuer diese Breite:
N = Ellips%a/sqrt(1.0_wp - Ellips%EN1sup2*((sin(beta))**2))
!dN = beta_tl *  (Ellips%a*Ellips%EN1sup2*sin(beta)*cos(beta)) /     &
!                 sqrt(1.0D0 - Ellips%EN1sup2*((sin(beta))**2))**3
f1 = (Ellips%a*Ellips%EN1sup2*sin(beta)*cos(beta)) /                &
                 sqrt(1.0_wp - Ellips%EN1sup2*((sin(beta))**2))**3
!dN = Ce_dv(:,2) * f1
!dN = 0.0_wp
dN = CeGrad(:,2) * f1


!write(*,*) 'Gr dN = ', dN
!write(*,*) 'N = ', N, (sin(beta))**2

! Die ellipsoidische Hoehe h:
h = (rho/cos(beta)) - N
!dh = -dN + drho/cos(beta) + Ce_dv(:,2)*rho*sin(beta)/(cos(beta)**2)
f1 = rho*sin(beta)/(cos(beta)**2)
dh = -dN  + drho/cos(beta) + CeGrad(:,2) * f1

!write(*,*) 'Gr dh = ', dh
!write(*,*) 'h = ', h, cos(beta)

if (abs(rho) > reps) then
   height = h
   !height_tl = dh
   CeGrad(:,3) = dh
else ! rho = 0
   if (abs(Cc(3)) < reps) then
      height = -Ellips%b
      CeGrad(:,3) = 0.0_wp
   else
      height = abs(Cc(3)) - Ellips%b
      CeGrad(:,3) = 0.0_wp
      CeGrad(3,3) = 1.0_wp
   end if
end if

!height_tl = dh ! ?????????
!CeGrad(:,3) = dh
Ce = (/lambda, beta, height/)

end subroutine Cart2EllipsGrad
!****
! <<<<<<<< End of RoboDoc comments

!---------------------------------------------------------------------
! subroutine Cart2EllipsGrad2
!---------------------------------------------------------------------
!
! Start of RoboDoc comments >>>>
!****S* Coord/Cart2EllipsGrad2
!
! Name
! Cart2EllipsGrad2
!
! Call
! call Cart2Ellips(X, Y, Z, lambda, beta, height, Ellips)
!
! Purpose
! Umrechnung von kartesischen Koordinaten in ellipsoidische Koordinaten
! Direkte Rechnung ohne Iterationsverfahren
!
! Quelle: B. Hofmann-Wellenhof, H. Lichtenegger, J. Collins
!         Global Positioning System - Theory and Practice
!         Springer, Wien, New York, 1993, Seite 232
!         Die Transformationsformel fuer die Laenge wurde leicht abgewandelt.
!
!         Hier fehlt die Sonderbehandlung von X=Y=Z=0, Winkel = 0
!         und die separate Behandlung der Oktanden, so dass wirklich
!         Laengen zwischen 0° und 360° und Breiten zwischen -90° und +90°
!         herauskommen. Dies wurde analog zur Transformation in sphaerische
!         Koordinaten ergaenzt.
!
! Zum Testen dieser Transformation steht das Programm "TestWgs.f90" in
! ../ToolsTest/Coord zur Verfuegung. Dort werden bekannte Stationskoordinaten
! sowohl in kartesischen als auch in WGS84 Koordinaten eingelesen und dann
! mit den Ergebnissen der Transformationen verglichen.
!

! Cart2Ellips_dv => Cart2EllipsGrad
! Cc_dv = (1, 0, 0) ;  (0, 1, 0) ; (0, 0, 1)
! Ce_dv = d/d lambda , d/ d beta, d/d height

! 2. derivative, i.e. derivative of CeGrad
!
! Ce - ellipsoidal coordinates
!      Ce(i) - i=1 lambda, longitude in rad
!              i=2 beta, latitude in rad
!              i=3 height above ellipsoid in m
!
! CeGrad - first derivatives
!          CeGrad(j,i) - i=1 lambda, longitude in rad
!                        i=2 beta, latitude in rad
!                        i=3 height above ellipsoid in m
!                      - j=1 d / dx
!                        j=2 d / dy
!                        j=3 d / dz
!
!               ( 11 12 13 )   ( lambda_x beta_x height_x )
! CeGrad(j,i) = ( 21 22 23 ) = ( lambda_y beta_y height_y )
!               ( 31 32 33 )   ( lambda_z beta_z height_z )
!
! Jacobian = (CeGrad)^T (transposed CeGrad)
!
!      ( 11 12 13 )   ( lambda_x lambda_y lambda_z )
!  J = ( 21 22 23 ) = ( beta_x   beta_y   beta_z   )
!      ( 31 32 33 )   ( height_x height_y height_z )
!
! CeGrad2 - second derivatives
!           CeGrad2(k,j,i) - i=1 lambda, longitude in rad
!                            i=2 beta, latitude in rad
!                            i=3 height above ellipsoid in m
!                          - j=1 d / dx
!                            j=2 d / dy
!                            j=3 d / dz
!                          - k=1 d / dx
!                            k=2 d / dy
!                            k=3 d / dz
!
! k = 1 - second derivatives with respect to x:
!
!                  ( 111 112 113 )   ( lambda_xx beta_xx height_xx )
! CeGrad2(1,j,i) = ( 121 122 123 ) = ( lambda_yx beta_yx height_yx )
!                  ( 131 132 133 )   ( lambda_zx beta_zx height_zx )
!
! k = 2 - second derivatives with respect to y:
!
!                  ( 211 212 213 )   ( lambda_xy beta_xy height_xy )
! CeGrad2(2,j,i) = ( 221 222 223 ) = ( lambda_yy beta_yy height_yy )
!                  ( 231 232 233 )   ( lambda_zy beta_zy height_zy )
!
! k = 3 - second derivatives with respect to z:
!
!                  ( 311 312 313 )   ( lambda_xz beta_xz height_xz )
! CeGrad2(3,j,i) = ( 321 322 323 ) = ( lambda_yz beta_yz height_yz )
!                  ( 331 332 333 )   ( lambda_zz beta_zz height_zz )
!
! beta_xz = d2 beta / dx dy

! first array index is highest derivative

! Input
! X, Y, Z - kartesische Koordinaten i m
! Ellips  - Struktur vom Typ "EllipsParam" mit den Daten des Referenzellipsoids
!
! Output
! lambda - ellipsoidische Laenge in rad
! beta   - ellipsoidische Breite in rad
! height - ellipsoidische Hoehe in m
!
! External References
!
! Examples
!
! Record of revisions
! Date         Programmer   Description of change
! ====         ==========   =====================
! 28.11.2013   M. Bender    neu
!---------------------------------------------------------------------
! Source
!---------------------------------------------------------------------
subroutine Cart2EllipsGrad2(Cc, Ce, Ellips, CeGrad, CeGrad2)

! List of calling arguments:
real (wp), dimension(1:3), intent(in)  :: Cc
real (wp), dimension(1:3), intent(out) :: Ce
!real (wp), dimension(:,:), intent(in)  :: Cc_dv
!real (wp), dimension(:,:), intent(out) :: Ce_dv
real (wp), dimension(1:3,1:3), intent(out) :: CeGrad
real (wp), dimension(1:3,1:3,1:3), intent(out) :: CeGrad2

type (EllipsParam) :: Ellips

! List of local variables:
real (wp) :: lambda, beta, height
real (wp) :: rho, phi, gamma, g2
real (wp) :: N, theta, h
integer   :: i, j, k

real (wp), dimension(1:3) :: drho, dtheta, dphi, dg2, dgamma, dN, dh
real (wp), dimension(1:3,1:3) :: d2rho, d2theta, d2phi, d2g2, d2gamma
real (wp), dimension(1:3,1:3) :: d2N, d2h
real (wp) :: c1, c2, f1, f2, e1, cc1, v1, a, b
real (wp), dimension(1:3) :: dc1, dc2, df1, df2, dv1, da, db
!integer   :: multi
!---------------------------------------------------------------------

! Einige Hilfsgroessen:
rho = sqrt( Cc(1)**2 + Cc(2)**2 )

drho(1) = Cc(1)
drho(2) = Cc(2)
drho(3) = 0.0_wp
drho = drho / rho

d2rho = 0.0_wp
e1 = rho**(-3)
d2rho(1,1) = e1 * Cc(2)**2  ! d2 rho / d x2
d2rho(2,2) = e1 * Cc(1)**2  ! d2 rho / d y2
d2rho(2,1) = - Cc(1)*Cc(2) * e1   ! d2 rho / dx d y
d2rho(1,2) = d2rho(2,1)           ! d2 rho / dy d x

c1 = (Cc(3)*Ellips%a)/(rho*Ellips%b)
cc1 = -c1/rho
dc1(1) = drho(1) * cc1              ! d c1 / d x
dc1(2) = drho(2) * cc1              ! d c1 / d y
dc1(3) = Ellips%a / (Ellips%b*rho)  ! d c1 / d z
if (rho .gt. 1.0E-6_wp) then
   theta = atan(c1)
   f1 = Ellips%a / (rho*Ellips%b*(1.0_wp+c1**2))
   f2 =  Cc(3) * Ellips%a /  (rho**2 * Ellips%b*(1.0_wp+c1**2))
   !dtheta = Z_tl * Ellips%a / (rho*Ellips%b*(1.0D0+c1**2)) -         &
   !         drho * Z * Ellips%a /  (rho**2 * Ellips%b*(1.0D0+c1**2))
   !do i=1, multi
   !   dtheta(i) = Cc_dv(i,3) * f1 - drho(i) * f2
   !end do
   df1(1) = - Ellips%a* Ellips%b*    &   ! d f1 / d x
              ((1.0_wp+c1**2)*drho(1) + 2.0_wp*rho*c1*dc1(1)) /   &
              (rho*Ellips%b*(1.0_wp+c1**2))**2
   df1(2) = - Ellips%a* Ellips%b*    &   ! d f1 / d y
              ((1.0_wp+c1**2)*drho(2) + 2.0_wp*rho*c1*dc1(2)) /   &
              (rho*Ellips%b*(1.0_wp+c1**2))**2
   df1(3) = - Ellips%a* Ellips%b*    &   ! d f1 / d z
              (2.0_wp*rho*c1*dc1(3)) /   &
              (rho*Ellips%b*(1.0_wp+c1**2))**2

   df2(1) = - 2.0_wp*rho* (Cc(3) * Ellips%a / Ellips%b) *     & ! d f2 / d x
              (rho*c1*dc1(1) + drho(1)*(1.0_wp+c1**2)) &
              / (rho**2 * (1.0_wp+c1**2))**2
   df2(2) = - 2.0_wp*rho* (Cc(3) * Ellips%a / Ellips%b) *  & ! d f2 / d y
              (rho*c1*dc1(2) + drho(2)*(1.0_wp+c1**2))     &
              / (rho**2 * (1.0_wp+c1**2))**2
   df2(3) =  (Ellips%a / Ellips%b) *                       & ! d f2 / d z
             ( (rho**2 * (1.0_wp+c1**2))**(-1) -  &
               (2.0_wp * Cc(3) * rho**2 * c1 * dc1(3))     &
               /  (rho**2 * (1.0_wp+c1**2))**2 )

   dtheta(1) = -drho(1) * f2
   dtheta(2) = -drho(2) * f2
   dtheta(3) = f1

   do j=1, 2
      ! first derivative: 1 - X, 2 - Y, 3 - Z
      do i=1, 3
         ! second derivative
         d2theta(i,j) = -(d2rho(i,j)*f2 + drho(j)*df2(i))
      end do
   end do
   d2theta(:,3) = df1
else
   theta = pi05
   dtheta = 0.0_wp
   d2theta = 0.0_wp
end if
!write(*,*) 'Gr dtheta = ', dtheta

! Die ellipsoidische Laenge lambda:
if (abs(rho) > reps) then
   phi = 2.0_wp*atan(Cc(2)/(abs(Cc(1))+rho))
   f1 = -sign(2.0_wp*Cc(2),Cc(1))
   f2 =  2.0_wp*(abs(Cc(1)) + rho)
   !dphi = ( -sign(2.0D0*Y,X) * X_tl + 2.0D0*(abs(X) + rho) * Y_tl   -    &
   !         2.0D0*Y * drho ) / ((abs(X)+rho)**2 + Y**2)
   !do i=1, multi
   !   dphi(i) = Cc_dv(i,1) * f1 +  Cc_dv(i,2) * f2 -   &
   !             2.0_wp * Cc(2) * drho(i)
   !end do
   dphi(1) = f1 - 2.0_wp * Cc(2) * drho(1)
   dphi(2) = f2 - 2.0_wp * Cc(2) * drho(2)
   dphi(3) = 0.0_wp
   v1 =  (abs(Cc(1))+rho)**2 + Cc(2)**2
   dphi = dphi / v1

   ! 2. derivatives
   df1    = 0.0_wp
   df1(2) = -sign(2.0_wp,Cc(1))

   df2(1) = 2.0_wp * ( sign(1.0_wp,Cc(1)) + drho(1) )
   df2(2) = 2.0_wp * drho(2)
   df2(3) = 0.0_wp

   dv1(1) = f2 * ( sign(1.0_wp,Cc(1)) + drho(1) )
   dv1(2) = f2 * drho(2) + 2.0_wp * Cc(2)
   dv1(3) = 0.0_wp

   d2phi = 0.0_wp
   d2phi(1,1) = - 2.0_wp * Cc(2) * d2rho(1,1)
   d2phi(2,1) = df1(2) - 2.0_wp * ( Cc(2) * d2rho(2,1) + drho(1) )
   d2phi(1,2) = df2(1) - 2.0_wp * Cc(2) * d2rho(1,2)
   d2phi(2,2) = df2(2) - 2.0_wp * ( Cc(2) * d2rho(2,2) + drho(2) )

   do j=1, 2
      do i=1, 2
         d2phi(j,i) = v1 * d2phi(j,i) - dphi(i) *v1 * dv1(j)
         d2phi(j,i) = d2phi(j,i) / v1**2
      end do
   end do

end if
if (Cc(1) .ge. 0.0_wp) then
   if (abs(Cc(1)) .lt. 1.0E-8_wp .and. abs(Cc(2)) .lt. 1.0E-8_wp) then
      ! if X=0 and Y=0
      lambda = 0.0_wp
      !lambda_tl = 0.0D0
      CeGrad(:,1) = 0.0_wp   ! first derivatives of lambda
      CeGrad2 = 0.0_wp       ! second derivatives of lambda
   else if (Cc(2) .ge. 0.0_wp) then
      ! X>0 and Y>0
      lambda = phi
      !lambda_tl = dphi
      CeGrad(:,1) = dphi              ! first derivatives of lambda
      do k=1, 3
         CeGrad2(k,:,1) = d2phi(k,:)  ! second derivatives of lambda
      end do
   else
      ! X>0 and Y<0
      lambda = phi + 2.0_wp*pi
      !lambda_tl = dphi
      CeGrad(:,1) = dphi              ! first derivatives of lambda
      do k=1, 3
         CeGrad2(k,:,1) = d2phi(k,:)  ! second derivatives of lambda
      end do
   end if
else  ! X < 0
   lambda = pi - phi
   !lambda_tl = -dphi
   CeGrad(:,1) = -dphi              ! first derivatives of lambda
   do k=1, 3
      CeGrad2(k,:,1) = -d2phi(k,:)  ! second derivatives of lambda
   end do
end if
!write(*,*) 'Gr dphi = ', dphi

! Die ellipsoidische Breite beta:
g2 = rho - Ellips%EN1sup2*Ellips%a*(cos(theta))**3  ! Nenner von gamma
!dg2 = drho + dtheta *                                                      &
!             3.0D0 * Ellips%EN1sup2*Ellips%a * sin(theta) * cos(theta)**2
f1 = 3.0_wp * Ellips%EN1sup2*Ellips%a * sin(theta) * cos(theta)**2
dg2 = drho + f1*dtheta

df1(:) = 3.0_wp * Ellips%EN1sup2*Ellips%a * cos(theta) * dtheta(:) *  &
         (cos(theta)**2 - 2.0_wp*sin(theta)**2)
do j=1, 3
   do i=1,3
      d2g2(j,i) = d2rho(j,i) + df1(j)*dtheta(i) + f1*d2theta(j,i)
   end do
end do

if (abs(g2) > reps) then
   gamma = atan( (Cc(3) + Ellips%EN2sup2*Ellips%b*(sin(theta))**3) / g2 )
   c2 = 1.0_wp/(1.0_wp +((Cc(3)+Ellips%EN2sup2*Ellips%b*(sin(theta))**3)/g2)**2)
   f1 = (Cc(3)+Ellips%EN2sup2*Ellips%b*(sin(theta)**3))/g2**2
   f2 = 3.0_wp*Ellips%EN2sup2*Ellips%b*cos(theta) * (sin(theta)**2) / g2
   !dgamma = c2 * (Z_tl/g2 -                                                 &
   !               dg2 * (Z+Ellips%EN2sup2*Ellips%b*(sin(theta)**3))/g2**2 + &
   !               dtheta * 3*Ellips%EN2sup2*Ellips%b*cos(theta) *           &
   !                        (sin(theta)**2)/g2 )
   !do i=1, multi
   !   dgamma(i) = Cc_dv(i,3)/g2 - dg2(i)*f1 + dtheta(i)*f2
   !end do
   dgamma(1) = -dg2(1)*f1 + dtheta(1)*f2
   dgamma(2) = -dg2(2)*f1 + dtheta(2)*f2
   dgamma(3) = 1.0_wp/g2 - dg2(3)*f1 + dtheta(3)*f2

   dgamma = dgamma * c2

   ! second derivatives of gamma
   dc2(1:2) = -2.0_wp * c2**2 *                                             &
              ( ( Cc(3)+Ellips%EN2sup2*Ellips%b*(sin(theta))**3 ) / g2 ) *  &
              ( 3.0_wp*Ellips%EN2sup2*Ellips%b*g2 *                         &
                cos(theta)*(sin(theta)**2)*dtheta(1:2) -                    &
              ((Cc(3)+Ellips%EN2sup2*Ellips%b*(sin(theta))**3)*dg2(1:2)))   &
              / g2**2
   dc2(3) =  -2.0_wp * c2**2 *                                                &
              ( ( Cc(3)+Ellips%EN2sup2*Ellips%b*(sin(theta))**3 ) / g2 ) *    &
              ( g2*(1.0_wp+3.0_wp*Ellips%EN2sup2*Ellips%b *                   &
                cos(theta)*(sin(theta)**2)*dtheta(3)) -                       &
              ((Cc(3)+Ellips%EN2sup2*Ellips%b*(sin(theta))**3)*dg2(3)))       &
              / g2**2

   df1(1:2) = ( 3.0_wp*Ellips%EN2sup2*Ellips%b *                            &
                g2*cos(theta)*(sin(theta)**2)*dtheta(1:2) - 2.0_wp *          &
                (Cc(3)+Ellips%EN2sup2*Ellips%b*(sin(theta))**3)*dg2(1:2) )  &
                / g2**3
   df1(3)   = ( g2*(1.0_wp + 3.0_wp*Ellips%EN2sup2*Ellips%b *               &
                    cos(theta)*(sin(theta)**2)*dtheta(3)) - 2.0_wp *        &
                (Cc(3)+Ellips%EN2sup2*Ellips%b*(sin(theta))**3)*dg2(3) )    &
                / g2**3

   df2(:) = ( g2*3.0_wp*Ellips%EN2sup2*Ellips%b *                    &
               ( 2.0_wp*sin(theta)*(cos(theta)**2)*dtheta(:) -       &
                 (sin(theta)**3)*dtheta(:)                    ) -    &
                 3.0_wp*Ellips%EN2sup2*Ellips%b*                     &
                 cos(theta)*(sin(theta)**2)*dg2(:)  ) / g2**2

   do j=1, 3
      do i=1, 2
         d2gamma(j,i) = ( -dg2(i)*df1(j) -f1*d2g2(j,i) + dtheta(i)*df2(j) + &
                        d2theta(j,i)*f2 ) * c2  +                           &
                        dc2(j) * dgamma(i)/c2
      end do
      d2gamma(j,3) = ( g2*dc2(j)-c2*dg2(j) ) / g2**2 +                     &
                     ( -dg2(3)*df1(j) -f1*d2g2(j,3) + dtheta(3)*df2(j) +   &
                        d2theta(j,3)*f2 ) * c2  +                          &
                        dc2(j) * (-dg2(3)*f1 + dtheta(3)*f2)
   end do

! else ???? +pi/2 oder -pi/2 ???
end if
if (abs(rho) .gt. 1.0E-8_wp) then
   ! rho > 0
   beta = gamma
   !beta_tl = dgamma
   CeGrad(:,2) = dgamma
   do k=1, 3
      CeGrad2(k,:,2) = d2gamma(k,:)  ! second derivatives of beta
   end do
else  ! rho = 0
   if (Cc(3) .gt. 0.0_wp) then
      beta = 0.50_wp*pi
      !beta_tl = 0.0D0
      CeGrad(:,2) = 0.0_wp
      CeGrad2(:,:,2) = 0.0_wp
   else if (abs(Cc(3)) < reps) then
      beta = 0.0_wp
      !beta_tl = 1.0D0
      CeGrad(:,2) = 1.0_wp
      CeGrad2(:,:,2) = 0.0_wp
   else if (Cc(3) .lt. 0.0_wp) then
      beta = -0.50_wp*pi
      !beta_tl = 0.0D0
      CeGrad(:,2) = 0.0_wp
      CeGrad2(:,:,2) = 0.0_wp
   end if
end if
!write(*,*) 'Gr dgamma = ', dgamma

! Bestimmung des Querkruemmungsradius fuer diese Breite:
N = Ellips%a/sqrt(1.0_wp - Ellips%EN1sup2*((sin(beta))**2))
!dN = beta_tl *  (Ellips%a*Ellips%EN1sup2*sin(beta)*cos(beta)) /     &
!                 sqrt(1.0D0 - Ellips%EN1sup2*((sin(beta))**2))**3
f1 = (Ellips%a*Ellips%EN1sup2*sin(beta)*cos(beta)) /                &
                 sqrt(1.0_wp - Ellips%EN1sup2*((sin(beta))**2))**3
!dN = Ce_dv(:,2) * f1
!dN = 0.0_wp
dN = CeGrad(:,2) * f1

! second derivatives
!a = (1.0_wp - Ellips%EN1sup2*((sin(beta))**2))**(-3/2)
a = (1.0_wp - Ellips%EN1sup2*((sin(beta))**2))**(-1.50_wp)
da(:) = CeGrad(:,2) * 3.0_wp * Ellips%EN1sup2 * sin(beta)*cos(beta) *  &
        (1.0_wp - Ellips%EN1sup2*((sin(beta))**2))**(-2.50_wp)
b = sin(beta)*cos(beta)
db(:) =  CeGrad(:,2) * ((cos(beta)**2) - (sin(beta)**2))
df1(:) = Ellips%a*Ellips%EN1sup2*(a*db(:) + da(:)*b)

do j=1, 3
   do i=1, 3
      d2N(j,i) = CeGrad2(j,i,2) * f1 +  CeGrad(i,2) * df1(j)
   end do
end do
!write(*,*) 'Gr dN = ', dN
!write(*,*) 'N = ', N, (sin(beta))**2

! Die ellipsoidische Hoehe h:
h = (rho/cos(beta)) - N
!dh = -dN + drho/cos(beta) + Ce_dv(:,2)*rho*sin(beta)/(cos(beta)**2)
f1 = rho*sin(beta)/(cos(beta)**2)
dh = -dN  + drho/cos(beta) + CeGrad(:,2) * f1

df1(:) = ( cos(beta)*(rho*cos(beta)*CeGrad(:,2) + drho(:)*sin(beta)) +  &
           2.0_wp*rho*(sin(beta)**2)*CeGrad(:,2) ) / (cos(beta)**3)
do j=1, 3
   do i=1, 3
      d2h(j,i) = -d2N(j,i) +                                              &
                 ( cos(beta)*d2rho(j,i) + sin(beta)*drho(i)*CeGrad(j,2) ) &
                 / (cos(beta)**2)  +                                      &
                 CeGrad(i,2)*df1(j) + CeGrad2(j,i,2)*f1
   end do
end do

!write(*,*) 'Gr dh = ', dh
!write(*,*) 'h = ', h, cos(beta)

if (abs(rho) > reps) then
   height = h
   !height_tl = dh
   CeGrad(:,3) = dh
   CeGrad2(:,:,3) = d2h
else ! rho = 0
   if (abs(Cc(3)) < reps) then
      height = -Ellips%b
      CeGrad(:,3) = 0.0_wp
      CeGrad2(:,:,3) = 0.0_wp
   else
      height = abs(Cc(3)) - Ellips%b
      CeGrad(:,3) = 0.0_wp
      CeGrad(3,3) = 1.0_wp
      CeGrad2(:,:,3) = 0.0_wp
   end if
end if

!height_tl = dh ! ?????????
!CeGrad(:,3) = dh
Ce = (/lambda, beta, height/)

end subroutine Cart2EllipsGrad2
!****
! <<<<<<<< End of RoboDoc comments


!---------------------------------------------------------------------
! subroutine Slant2EllipsMatGrad
!---------------------------------------------------------------------
!
! Start of RoboDoc comments >>>>
!****S* Coord/Slant2EllipsMatGrad
!
! Name
! Slant2EllipsMatGrad
!
! Call
! call Slant2Ellips_tl
!  (Xslant, Xellips, azi, elev, RefPtEllips, RefPtECEF, Ellips)
!
! Purpose
! Transforms the coordinates of a vector from the slant system into
! ellipsoidal geographic coordinates like WGS84.
! Inverse transformation: Ellips2Slant
!
! This is a driver routine for three consecutive transformations:
! Slant2LocalHorz => LocalHorz2Cart => Cart2Ellips
!
! Ellipsoidal system: ellipsoidal geographic coordinates
!
! Slant system: right handed  cartesian system with its origin
!                          at the position (lambda,phi).
!               The x-axis points to the satellite with the given azimuth
!               and elevation, the y-axis is in the tangent plane and the
!               z-axis is the vertical axis if the elevation is zero.
!               For epsilon <> 0, the verticl component is mapped on the
!               z- and x-axis.
!
! Input
! Xslant - Coordinates of the slant system
! azi    - azimuth (angle in rad)
!          A = [0, ..., 2*pi] or
!          A = [-pi, ..., pi]
! elev   - elevation (angle in rad)
!          elev = [0, ..., pi/2]
! RefPtEllips - coordinates of reference point, e.g. station,
!               ellipsoidal geographic coordinates (lon, lat, height)
! RefPtECEF   - coordinates of reference point, e.g. station,
!               earth centered earth fixed cartesian coordinates (ECEF)
! Ellips      - reference ellipsoid, e.g. WGS84
!
! Output
! Xellips - geographic coordinates (longitude, latitude, height)
!
! External References
! None
!
! Record of revisions
! Date        Programmer   Description of change
! ====        ==========   =====================
! 20.09.2013  M. Bender    new
!---------------------------------------------------------------------
! Source
!---------------------------------------------------------------------
subroutine Slant2EllipsMatGrad (Xslant, Xellips, XellipsGrad,     &
                               T1, T2, RefPtECEF, Ellips)

! List of calling arguments:
real (wp), dimension(1:3), intent(in)  :: Xslant
!real (wp), dimension(:,:), intent(in)  :: Xslant_dv
real (wp), dimension(1:3), intent(out) :: Xellips
!real (wp), dimension(:,:), intent(out) :: Xellips_dv
real (wp), dimension(1:3,1:3), intent(out) :: XellipsGrad
real (wp), dimension(:,:), intent(in)  :: T1, T2
real (wp), dimension(1:3), intent(in)  :: RefPtECEF
type (EllipsParam), intent(in)                    :: Ellips

! List of local variables:
real (wp), dimension(1:3)     :: Xhorz, Xecef
real (wp), dimension(1:3,1:3) :: Xhorz_dv, Xecef_dv
real (wp), dimension(1:3,1:3) :: XellipsGrad0 !, Xellips_dv
!real (wp)                 :: lambda_tl, phi_tl
!real (wp)                 :: ox_tl, oy_tl, oz_tl
integer                              :: i
!---------------------------------------------------------------------

!!$! Allocate local arrays of derivatives
!!$allocate( Xhorz_dv(lbound(Xslant_dv,1):ubound(Xslant_dv,1),    &
!!$                   lbound(Xslant_dv,1):ubound(Xslant_dv,1)) )
!!$allocate( Xecef_dv(lbound(Xslant_dv,1):ubound(Xslant_dv,1),    &
!!$                   lbound(Xslant_dv,1):ubound(Xslant_dv,1)) )

! Multidirectional differentiation: number of directions
!multi = ubound(Xslant_dv,1)

! Transform coordinates from the slant system to the local horizon system
!call Slant2LocalHorz_tl(Xslant, Xslant_tl, azi, elev, Xhorz, Xhorz_tl)
Xhorz    = matmul(T1, Xslant)
do i=1, 3
   Xhorz_dv(i,:) = T1(:,i)
end do
!!$do i=1, 3
!!$   write(*,*) 'Grad Xhorz_dv ', i, Xhorz_dv(i,:)
!!$end do
! ... local horizon system to earth centered earth fixed coordinates (ECEF)
!!$lambda_tl = 0.0_KindReal
!!$phi_tl    = 0.0_KindReal
!!$ox_tl     = 0.0_KindReal
!!$oy_tl     = 0.0_KindReal
!!$oz_tl     = 0.0_KindReal
!!$call LocalHorz2Cart_tl(Xhorz(1), Xhorz(2), Xhorz(3),               &
!!$                       RefPtEllips(1), RefPtEllips(2),             &
!!$                       RefPtECEF(1), RefPtECEF(2), RefPtECEF(3),   &
!!$                       Xecef(1), Xecef(2), Xecef(3),               &
!!$                       Xhorz_tl(1), Xhorz_tl(2), Xhorz_tl(3),      &
!!$                       lambda_tl, phi_tl,                          &
!!$                       ox_tl, oy_tl, oz_tl,                        &
!!$                       Xecef_tl(1), Xecef_tl(2), Xecef_tl(3)     )
Xecef = matmul(T2, Xhorz) + RefPtECEF
do i=1, 3
   Xecef_dv(i,:) = matmul(T2, Xhorz_dv(i,:))
end do

!!$do i=1, 3
!!$   write(*,*) 'Grad Xecef_dv ', i, Xecef_dv(i,:)
!!$end do

! ... ECEF to ellipsoidal coordinates
!!$call Cart2Ellips_tl(Xecef(1), Xecef(2), Xecef(3),    &
!!$                    Xellips(1), Xellips(2), Xellips(3), Ellips ,  &
!!$                    Xecef_tl(1), Xecef_tl(2), Xecef_tl(3),        &
!!$                    Xellips_tl(1), Xellips_tl(2), Xellips_tl(3)  )
!!$Xecef_dv(1,:) = (/1.0_wp, 0.0_wp, 0.0_wp/)
!!$Xecef_dv(2,:) = (/0.0_wp, 1.0_wp, 0.0_wp/)
!!$Xecef_dv(3,:) = (/0.0_wp, 0.0_wp, 1.0_wp/)
!write(*,*) 'Input Cart2Ellips_dv  Ce_dv(:,2) = ',  Ce_dv(:,2)
!call Cart2Ellips_dv(Xecef, Xellips, Ellips, Xecef_dv, Xellips_dv)
!write(*,*) 'DV Xellips = ', Xellips
call Cart2EllipsGrad(Xecef, Xellips, Ellips, XellipsGrad0)
!write(*,*) 'Gr Xellips = ', Xellips
XellipsGrad = matmul(Xecef_dv,XellipsGrad0)
!XellipsGrad = matmul(XellipsGrad0,Xecef_dv)
!!$XellipsGrad = 0.0_wp
!!$do i=1, 3
!!$   !XellipsGrad(i,i) =  Xecef_dv(i,i) * XellipsGrad0(i,i)
!!$   XellipsGrad(i,:) =  dot_product(Xecef_dv(i,:), XellipsGrad0(i,:))
!!$end do

!!$do i=1, 3
!!$   write(*,*) 'DV    ', Xellips_dv(i,:)
!!$   write(*,*) 'Grad0 ', XellipsGrad0(i,:)
!!$   write(*,*) 'Grad  ', XellipsGrad(i,:)
!!$end do

!!$deallocate( Xhorz_dv )
!!$deallocate( Xecef_dv )

end subroutine Slant2EllipsMatGrad
!****
! <<<<<<<< End of RoboDoc comments


!---------------------------------------------------------------------
! subroutine Slant2EllipsMatGrad2
!---------------------------------------------------------------------
!
! Start of RoboDoc comments >>>>
!****S* Coord/Slant2EllipsMatGrad2
!
! Name
! Slant2EllipsMatGrad2
!
! Call
! call Slant2Ellips_tl
!  (Xslant, Xellips, azi, elev, RefPtEllips, RefPtECEF, Ellips)
!
! Purpose
! Transforms the coordinates of a vector from the slant system into
! ellipsoidal geographic coordinates like WGS84.
! Inverse transformation: Ellips2Slant
!
! This is a driver routine for three consecutive transformations:
! Slant2LocalHorz => LocalHorz2Cart => Cart2Ellips
!
! Slant system (x,y,z) => Slant2LocalHorz => local horizon system (xh,yh,zh)
! local horizon system (xh,yh,zh) => LocalHorz2Cart => ECEF (Y, Y, Z)
! ECEF (Y, Y, Z) => Cart2Ellips => ellipsoidal coordinates (lambda,beta,height)
!
! Ellipsoidal system: ellipsoidal geographic coordinates
!
! Slant system: right handed  cartesian system with its origin
!                          at the position (lambda,phi).
!               The x-axis points to the satellite with the given azimuth
!               and elevation, the y-axis is in the tangent plane and the
!               z-axis is the vertical axis if the elevation is zero.
!               For epsilon <> 0, the verticl component is mapped on the
!               z- and x-axis.
!

! LocalHorz2Cart - local horizon system to ECEF system
! X = X(xh,yh,zh), Y = Y(xh,yh,zh), Z = Z(xh,yh,zh)
!
!         ( X_xh  X_yh  X_zh )
! Jecef = ( Y_xh  Y_yh  Y_zh )
!         ( Z_xh  Z_yh  Z_zh )
!
!           ( 11 12 13 )   ( X_xh  Y_xh  Z_xh )
! JtXecef = ( 21 22 23 ) = ( X_yh  Y_yh  Z_yh )
!           ( 31 32 33 )   ( X_zh  Y_zh  Z_zh )
!
! XellipsGrad2 => see Cart2EllipsGrad2, CeGrad2
!
! Input
! Xslant - Coordinates of the slant system
! azi    - azimuth (angle in rad)
!          A = [0, ..., 2*pi] or
!          A = [-pi, ..., pi]
! elev   - elevation (angle in rad)
!          elev = [0, ..., pi/2]
! RefPtEllips - coordinates of reference point, e.g. station,
!               ellipsoidal geographic coordinates (lon, lat, height)
! RefPtECEF   - coordinates of reference point, e.g. station,
!               earth centered earth fixed cartesian coordinates (ECEF)
! Ellips      - reference ellipsoid, e.g. WGS84
!
! Output
! Xellips - geographic coordinates (longitude, latitude, height)
!
! External References
! None
!
! Record of revisions
! Date        Programmer   Description of change
! ====        ==========   =====================
! 23.01.2014  M. Bender    new
!---------------------------------------------------------------------
! Source
!---------------------------------------------------------------------
subroutine Slant2EllipsMatGrad2 (Xslant, Xellips, XellipsGrad, XellipsGrad2, &
                                 T1, T2, RefPtECEF, Ellips)

! List of calling arguments:
real (wp), dimension(1:3), intent(in)  :: Xslant
!real (wp), dimension(:,:), intent(in)  :: Xslant_dv
real (wp), dimension(1:3), intent(out) :: Xellips
!real (wp), dimension(:,:), intent(out) :: Xellips_dv
real (wp), dimension(1:3,1:3), intent(out) :: XellipsGrad
real (wp), dimension(1:3,1:3,1:3), intent(out) :: XellipsGrad2

real (wp), dimension(:,:), intent(in)  :: T1, T2
real (wp), dimension(1:3), intent(in)  :: RefPtECEF
type (EllipsParam), intent(in)                    :: Ellips

! List of local variables:
real (wp), dimension(1:3)         :: Xhorz, Xecef
real (wp), dimension(1:3,1:3)     :: JtXecef !, Xhorz_dv
real (wp), dimension(1:3,1:3)     :: CartEllipsGrad
real (wp), dimension(1:3,1:3,1:3) :: CartEllipsGrad2
!real (wp)                 :: lambda_tl, phi_tl
!real (wp)                 :: ox_tl, oy_tl, oz_tl
integer                              :: i
!---------------------------------------------------------------------


! Transform coordinates from the slant system to the local horizon system
Xhorz    = matmul(T1, Xslant)
! linear transformation => T1 is the Jacobian containing the first derivatives,
! all second derivatives are 0

! ... local horizon system to earth centered earth fixed coordinates (ECEF)
Xecef = matmul(T2, Xhorz) + RefPtECEF
! linear transformation => T2 is the Jacobian containing the first derivatives,
! all second derivatives are 0
JtXecef = transpose( matmul(T2,T1) )

! ... ECEF to ellipsoidal coordinates
call Cart2EllipsGrad2(Xecef, Xellips, Ellips, CartEllipsGrad, CartEllipsGrad2)
! first derivatives
! "CartEllipsGrad" is the transposed Jacobian, therefore the transposed
! Jacobian "JtXecef" is required to compute the transposed Jacobian of all
! three consecitive transformations, i.e. "XellipsGrad":
XellipsGrad = matmul(JtXecef, CartEllipsGrad)

! second derivatives
do i=1, 3
   ! second derivative x x
   XellipsGrad2(1,1,i) = CartEllipsGrad2(1,1,i) * JtXecef(1,1)**2 +          &
                         CartEllipsGrad2(2,2,i) * JtXecef(1,2)**2 +          &
                         CartEllipsGrad2(3,3,i) * JtXecef(1,3)**2 +          &
                         2.0_wp * (                                          &
                          CartEllipsGrad2(2,1,i)*JtXecef(1,1)*JtXecef(1,2) + &
                          CartEllipsGrad2(3,1,i)*JtXecef(1,1)*JtXecef(1,3) + &
                          CartEllipsGrad2(3,2,i)*JtXecef(1,2)*JtXecef(1,3) )
   ! second derivative y y
   XellipsGrad2(2,2,i) = CartEllipsGrad2(1,1,i) * JtXecef(2,1)**2 +          &
                         CartEllipsGrad2(2,2,i) * JtXecef(2,2)**2 +          &
                         CartEllipsGrad2(3,3,i) * JtXecef(2,3)**2 +          &
                         2.0_wp * (                                          &
                          CartEllipsGrad2(2,1,i)*JtXecef(2,1)*JtXecef(2,2) + &
                          CartEllipsGrad2(3,1,i)*JtXecef(2,1)*JtXecef(2,3) + &
                          CartEllipsGrad2(3,2,i)*JtXecef(2,2)*JtXecef(2,3) )
   ! second derivative z z
   XellipsGrad2(3,3,i) = CartEllipsGrad2(1,1,i) * JtXecef(3,1)**2 +          &
                         CartEllipsGrad2(2,2,i) * JtXecef(3,2)**2 +          &
                         CartEllipsGrad2(3,3,i) * JtXecef(3,3)**2 +          &
                         2.0_wp * (                                          &
                          CartEllipsGrad2(2,1,i)*JtXecef(3,1)*JtXecef(3,2) + &
                          CartEllipsGrad2(3,1,i)*JtXecef(3,1)*JtXecef(3,3) + &
                          CartEllipsGrad2(3,2,i)*JtXecef(3,2)*JtXecef(3,3) )
   ! second derivative x y
   XellipsGrad2(2,1,i) = CartEllipsGrad2(1,1,i)*JtXecef(1,1)*JtXecef(2,1) +  &
                         CartEllipsGrad2(2,2,i)*JtXecef(1,2)*JtXecef(2,2) +  &
                         CartEllipsGrad2(3,3,i)*JtXecef(1,3)*JtXecef(2,3) +  &
                         CartEllipsGrad2(2,1,i)*(JtXecef(1,1)*JtXecef(2,2) + &
                                JtXecef(2,1)*JtXecef(1,2)) +                 &
                         CartEllipsGrad2(3,1,i)*(JtXecef(1,1)*JtXecef(2,3) + &
                                JtXecef(2,1)*JtXecef(1,3)) +                 &
                         CartEllipsGrad2(3,2,i)*(JtXecef(1,2)*JtXecef(2,3) + &
                                 JtXecef(2,2)*JtXecef(1,3))
   ! second derivative x z
   XellipsGrad2(3,1,i) = CartEllipsGrad2(1,1,i)*JtXecef(1,1)*JtXecef(3,1) +  &
                         CartEllipsGrad2(2,2,i)*JtXecef(1,2)*JtXecef(3,2) +  &
                         CartEllipsGrad2(3,3,i)*JtXecef(1,3)*JtXecef(3,3) +  &
                         CartEllipsGrad2(2,1,i)*(JtXecef(1,1)*JtXecef(3,2) + &
                                JtXecef(3,1)*JtXecef(1,2)) +                 &
                         CartEllipsGrad2(3,1,i)*(JtXecef(1,1)*JtXecef(3,3) + &
                                JtXecef(3,1)*JtXecef(1,3)) +                 &
                         CartEllipsGrad2(3,2,i)*(JtXecef(1,2)*JtXecef(3,3) + &
                                 JtXecef(3,2)*JtXecef(1,3))
   ! second derivative y x
   XellipsGrad2(1,2,i) = XellipsGrad2(2,1,i)
   ! second derivative y z
   XellipsGrad2(3,2,i) = CartEllipsGrad2(1,1,i)*JtXecef(2,1)*JtXecef(3,1) +  &
                         CartEllipsGrad2(2,2,i)*JtXecef(2,2)*JtXecef(3,2) +  &
                         CartEllipsGrad2(3,3,i)*JtXecef(2,3)*JtXecef(3,3) +  &
                         CartEllipsGrad2(2,1,i)*(JtXecef(2,1)*JtXecef(3,2) + &
                                JtXecef(3,1)*JtXecef(2,2)) +                 &
                         CartEllipsGrad2(3,1,i)*(JtXecef(2,1)*JtXecef(3,3) + &
                                JtXecef(3,1)*JtXecef(2,3)) +                 &
                         CartEllipsGrad2(3,2,i)*(JtXecef(2,2)*JtXecef(3,3) + &
                                 JtXecef(3,2)*JtXecef(2,3))
   ! second derivative z x
   XellipsGrad2(1,3,i) = XellipsGrad2(3,1,i)
   ! second derivative z y
   XellipsGrad2(2,3,i) = XellipsGrad2(3,2,i)
end do

end subroutine Slant2EllipsMatGrad2
!****
! <<<<<<<< End of RoboDoc comments



!---------------------------------------------------------------------
! subroutine LocalHorz2CartMat
!---------------------------------------------------------------------
!
! Start of RoboDoc comments >>>>
!****S* Coord/LocalHorz2CartMat
!
! Name
! LocalHorz2Cart
!
! Call
! call LocalHorz2Cart(x2, y2, z2, lambda, phi,   &
!                     ox, oy, oz, x1, y1, z1)
!
! Purpose
! Transformation from a local horizontal system to a global cartesian system
!
! Global system: right handed cartesian system
! Local horizontal system: left handed cartesian system with its origin
!                          at the position (lambda,phi,alt) of the
!                          observation point, e. g. the GPS station
!
! The transformation is identical for spheres and ellipsoids, see documentation
!
! Input
! x2,y2,z2  - cartesian coordinates in the local horizontal system
! lambda    - longitude of observation point, rad
! phi       - latitude of observation point, rad
! ox,oy,oz  - global cartesian coordinates of the observation point
!             (e. g. GPS station).
!             ox,oy,oz and lambda,phi are not independent values, but
!             are related by the choosen frame of reference
!             (ellipse, sphere...)
!
! Output
! x1,y1,z1  - global cartesian coordinates
!
! External References
! None
!
! Examples
!
! Record of revisions
! Date        Programmer   Description of change
! ====        ==========   =====================
! 18.11.2013  M. Bender    new
!---------------------------------------------------------------------
! Source
!---------------------------------------------------------------------
subroutine LocalHorz2CartMat(lambda, phi, T)

! List of calling arguments:
real (wp), intent(in)                      :: lambda, phi
real (wp), dimension(1:3,1:3), intent(out) :: T

! List of local variables:
!---------------------------------------------------------------------

!!$x1 = ox -sin(phi)*cos(lambda)*x2 -      &
!!$         sin(lambda)*y2          +      &
!!$         cos(phi)*cos(lambda)*z2
!!$y1 = oy -sin(phi)*sin(lambda)*x2 +      &
!!$         cos(lambda)*y2          +      &
!!$         cos(phi)*sin(lambda)*z2
!!$z1 = oz +cos(phi)*x2             +      &
!!$         sin(phi)*z2

! Matrix, row 1
T(1,1) =  -sin(phi)*cos(lambda)
T(1,2) = - sin(lambda)
T(1,3) =   cos(phi)*cos(lambda)

! Matrix, row 2
T(2,1) = -sin(phi)*sin(lambda)
T(2,2) =  cos(lambda)
T(2,3) =  cos(phi)*sin(lambda)

! Matrix, row 3
T(3,1) = cos(phi)
T(3,2) = 0.0_wp
T(3,3) = sin(phi)

end subroutine LocalHorz2CartMat
!****
! <<<<<<<< End of RoboDoc comments


!---------------------------------------------------------------------
! subroutine Slant2EllipsMat
!---------------------------------------------------------------------
!
! Start of RoboDoc comments >>>>
!****S* Coord/Slant2EllipsMat
!
! Name
! Slant2Ellips
!
! Call
! call Slant2Ellips(Xslant, Xellips, azi, elev, RefPtEllips, RefPtECEF, Ellips)
!
! Purpose
! Transforms the coordinates of a vector from the slant system into
! ellipsoidal geographic coordinates like WGS84.
! Inverse transformation: Ellips2Slant
!
! This is a driver routine for three consecutive transformations:
! Slant2LocalHorz => LocalHorz2Cart => Cart2Ellips
!
! Ellipsoidal system: ellipsoidal geographic coordinates
!
! Slant system: right handed  cartesian system with its origin
!                          at the position (lambda,phi).
!               The x-axis points to the satellite with the given azimuth
!               and elevation, the y-axis is in the tangent plane and the
!               z-axis is the vertical axis if the elevation is zero.
!               For epsilon <> 0, the verticl component is mapped on the
!               z- and x-axis.
!
! Input
! Xslant - Coordinates of the slant system
! azi    - azimuth (angle in rad)
!          A = [0, ..., 2*pi] or
!          A = [-pi, ..., pi]
! elev   - elevation (angle in rad)
!          elev = [0, ..., pi/2]
! RefPtEllips - coordinates of reference point, e.g. station,
!               ellipsoidal geographic coordinates (lon, lat, height)
! RefPtECEF   - coordinates of reference point, e.g. station,
!               earth centered earth fixed cartesian coordinates (ECEF)
! Ellips      - reference ellipsoid, e.g. WGS84
!
! Output
! Xellips - geographic coordinates (longitude, latitude, height)
!
! External References
! None
!
! Record of revisions
! Date        Programmer   Description of change
! ====        ==========   =====================
! 18.11.2013  M. Bender    new
!---------------------------------------------------------------------
! Source
!---------------------------------------------------------------------
subroutine Slant2EllipsMat (Xslant, Xellips, T1, T2, RefPtECEF, Ellips)

! List of calling arguments:
real (wp), dimension(1:3), intent(in)     :: Xslant
real (wp), dimension(1:3), intent(out)    :: Xellips
real (wp), dimension(1:3,1:3), intent(in) :: T1, T2
real (wp), dimension(1:3), intent(in)     :: RefPtECEF
type (EllipsParam), intent(in)            :: Ellips

! List of local variables:
real (wp), dimension(1:3) :: Xhorz, Xecef
!---------------------------------------------------------------------

! Transform coordinates from the slant system to the local horizon system
! call Slant2LocalHorz (Xslant, azi, elev, Xhorz)
Xhorz = matmul(T1, Xslant)

! ... local horizon system to earth centered earth fixed coordinates (ECEF)
!call LocalHorz2Cart(Xhorz(1), Xhorz(2), Xhorz(3),               &
!                    RefPtEllips(1), RefPtEllips(2),             &
!                    RefPtECEF(1), RefPtECEF(2), RefPtECEF(3),   &
!                    Xecef(1), Xecef(2), Xecef(3)              )
Xecef = matmul(T2, Xhorz) + RefPtECEF

! ... ECEF to ellipsoidal coordinates
call Cart2Ellips(Xecef(1), Xecef(2), Xecef(3),    &
                 Xellips(1), Xellips(2), Xellips(3), Ellips)

end subroutine Slant2EllipsMat
!****
! <<<<<<<< End of RoboDoc comments


!---------------------------------------------------------------------
! subroutine Slant2LocalHorzMat
!---------------------------------------------------------------------
!
! Start of RoboDoc comments >>>>
!****S* Coord/Slant2LocalHorz
!
! Name
! Slant2LocalHorz
!
! Call
! call Slant2LocalHorzMat (A, elev, T)
!
! Purpose
! Transformation from the slant system into ! the local horizon system:
! Computes the transformation matrix T for a given elevation and azimuth.
! In cases were the elevation and azimuth do not change the transformation
! can be computed with a simple matrix-vecor-product:
! Xhorz = matmul(T,Xslant)
!
! Inverse transformation: LocalHorz2SlantMat
!
! Transformation for A = const., elev = const.:
! Xhorz = T * Xslant
!
! Local horizontal system: left handed cartesian system with its origin
!                          at the position (lambda,phi).
! Slant system: right handed  cartesian system with its origin
!                          at the position (lambda,phi).
!               The x-axis points to the satellite with the given azimuth
!               and elevation, the y-axis is in the tangent plane and the
!               z-axis is the vertical axis if the elevation is zero.
!               For epsilon <> 0, the verticl component is mapped on the
!               z- and x-axis.
!
!
! Input
! A      - azimuth of the object (angle in rad)
!          A = [0, ..., 2*pi] or
!          A = [-pi, ..., pi]
! elev   - elevation (angle in rad)
!          elev = [0, ..., pi/2]
!
! Output
! T - transformation matrix, 3x3 matrix
!
! External References
! None
!
! Examples
!
! Record of revisions
! Date        Programmer   Description of change
! ====        ==========   =====================
! 18.11.2013  M. Bender    neu
!---------------------------------------------------------------------
! Source
!---------------------------------------------------------------------
subroutine Slant2LocalHorzMat (A, elev, T)

! List of calling arguments:
real (wp), intent(in)                      :: A, elev
real (wp), dimension(1:3,1:3), intent(out) :: T

! List of local variables:
!---------------------------------------------------------------------

! Slant2LocalHorz:
!!$Xhorz(1) = Xslant(1)*cos(elev)*cos(A) +    &
!!$           Xslant(2)*sin(A)           -    &
!!$           Xslant(3)*sin(elev)*cos(A)
!!$
!!$Xhorz(2) =  Xslant(1)*cos(elev)*sin(A) -    &
!!$            Xslant(2)*cos(A)           -    &
!!$            Xslant(3)*sin(elev)*sin(A)
!!$
!!$Xhorz(3) =  Xslant(1)*sin(elev) + Xslant(3)*cos(elev)

! Matrix, row 1
T(1,1) =  cos(elev)*cos(A)
T(1,2) =  sin(A)
T(1,3) = -sin(elev)*cos(A)

! Matrix, row 2
T(2,1) =  cos(elev)*sin(A)
T(2,2) = -cos(A)
T(2,3) = -sin(elev)*sin(A)

! Matrix, row 3
T(3,1) = sin(elev)
T(3,2) = 0.0_wp
T(3,3) = cos(elev)

end subroutine Slant2LocalHorzMat
!****
! <<<<<<<< End of RoboDoc comments


!---------------------------------------------------------------------
! subroutine ExpInt1DGrad
!---------------------------------------------------------------------
!
! Start of RoboDoc comments >>>>
!****S* Interpolation/ExpInt1DGrad
!
! Name
! ExpInt1DGrad
!
! Purpose
! One dimensional exponential interpolation using
! N(z) = N2 * (N1/N2)*exp((h2-z)/(h2-h1))
!
! To be used for, e. g., the vertical exponential interpolation of the
! atmospheric refractivity.
!
! Test of gradient with finite differences:
! .../gpstomo/trunk/Tomo/TomoTest/ModInterpolDiff/ExpInt1DGradTest.f90
!
! Call
! val = ExpInt1D_tl(RefPt, node)
!
! Input
!
! RefPt - scalar value, the values given in "node" will be interpolated
!         at this point.
! node  - array containing the coordinates and the values at these points
!            node(i,:), i=1, 2 : 2 nodes at the end of the interval
!            node(i,j), j=1, 2    : j = 1 - x, coordinate
!                                   j = 2 - f(x), value, observation, ...
!
! Output
! BiLin2D (return value) - interpolated value
!
! External References
! None
!
! Record of revisions
! Date        Programmer   Description of change
! ====        ==========   =====================
! 28.11.2013  M. Bender    new
!---------------------------------------------------------------------
! Source
!---------------------------------------------------------------------
subroutine ExpInt1DGrad(RefPt, inode, val, GradVal)

! List of calling arguments:
real (wp), intent(in)                     :: RefPt
real (wp), dimension(1:2,1:2), intent(in) :: inode
real (wp), intent(out)                    :: val
real (wp), intent(out)                    :: GradVal

! List of local variables:
real (wp), dimension(1:2,1:2) :: n
real (wp) :: E, L, B, NE
!---------------------------------------------------------------------

if (abs(inode(1,2)) .lt. 1E-10_wp .and. abs(inode(2,2)) .lt. 1E-10_wp) then
   ! N1 = N2 = 0  => N(z) = 0

   val     = 0.0_wp
   GradVal = 0.0_wp

else if (abs(inode(2,1)-inode(1,1)) < 1.0E-3_wp) then
   ! identical coordinates, interpolation not possible, return first value
   val = inode(1,2)
   GradVal = 0.0_wp

else
   ! interpolate

   n = inode

   if (abs(n(1,2)) .lt. 1E-10_wp) then
      ! N1 = 0, N2 <> 0
      n(1,2) = n(2,2)/100.0_wp
   end if
   if (abs(n(2,2)) .lt. 1E-10_wp) then
      ! N2 = 0, N1 <> 0
      n(2,2) = n(1,2)/100.0_wp
   end if

   val = n(2,2) * (n(1,2)/n(2,2))**((n(2,1)-RefPt)/(n(2,1)-n(1,1)))

   B = (n(2,1)-RefPt) / (n(2,1)-n(1,1))
   L = log( n(1,2)/n(2,2) )
   E = exp( B * L )
   NE = n(2,2) * E

   ! GradVal = d I / d RefPt
   GradVal = -L/(n(2,1)-n(1,1)) * NE

end if

end subroutine ExpInt1DGrad
!****
! <<<<<<<< End of RoboDoc comments


!---------------------------------------------------------------------
! subroutine ExpInt1DGrad2
!---------------------------------------------------------------------
!
! Start of RoboDoc comments >>>>
!****S* Interpolation/ExpInt1DGrad2
!
! Name
! ExpInt1DGrad2
!
! Purpose
! One dimensional exponential interpolation using
! N(z) = N2 * (N1/N2)*exp((h2-z)/(h2-h1))
!
! To be used for, e. g., the vertical exponential interpolation of the
! atmospheric refractivity.
!
! Test of gradient with finite differences:
! .../gpstomo/trunk/Tomo/TomoTest/ModInterpolDiff/ExpInt1DGradTest.f90
!
! Call
! val = ExpInt1D_tl(RefPt, node)
!
! Input
!
! RefPt - scalar value, the values given in "node" will be interpolated
!         at this point.
! node  - array containing the coordinates and the values at these points
!            node(i,:), i=1, 2 : 2 nodes at the end of the interval
!            node(i,j), j=1, 2    : j = 1 - x, coordinate
!                                   j = 2 - f(x), value, observation, ...
!
! Output
! BiLin2D (return value) - interpolated value
!
! External References
! None
!
! Record of revisions
! Date        Programmer   Description of change
! ====        ==========   =====================
! 22.01.2014  M. Bender    new
!---------------------------------------------------------------------
! Source
!---------------------------------------------------------------------
subroutine ExpInt1DGrad2(RefPt, inode, val, GradVal, GradVal2)

! List of calling arguments:
real (wp), intent(in)                     :: RefPt
real (wp), dimension(1:2,1:2), intent(in) :: inode
real (wp), intent(out)                    :: val
real (wp), intent(out)                    :: GradVal
real (wp), intent(out)                    :: GradVal2

! List of local variables:
real (wp), dimension(1:2,1:2) :: n
real (wp) :: exponent, E, L, B, NE, dB, dE
!---------------------------------------------------------------------

!!$integer :: i
!!$real(wp) :: bb(17), ee(17), r


!!$bb = (/1.1621608095682976, &
!!$      1.1621608095682976, &
!!$      1.1621608095682976, &
!!$      1.1621608095682976, &
!!$      1.1621608095682976, &
!!$      1.1621608095682976, &
!!$      1.1621608095682976, &
!!$      1.1621608095682976, &
!!$      1.1621608095682976, &
!!$      1.1621608095682976, &
!!$      1.1621608095682976, &
!!$      1.1621608095682976, &
!!$      1.1621608095682976, &
!!$      1.1621608095682976, &
!!$      1.1621608095682976, &
!!$      1.1621608095682976, &
!!$      1.1621608095682976/)
!!$
!!$
!!$ee = (/ -25.798663957450692, &
!!$       -28.324449502193534, &
!!$       -30.962485051019854, &
!!$       -33.712772906657833, &
!!$       -36.575314840237709, &
!!$       -39.550112187233445, &
!!$       -42.637165939029941, &
!!$       -45.83647682395452 , &
!!$       -49.148045375031792, &
!!$       -52.571871984021712, &
!!$       -56.107956942897118, &
!!$       -59.756300474507853, &
!!$       -63.516902754593168, &
!!$       -67.389763927120882, &
!!$       -71.374884114680356, &
!!$       -75.472263425383531, &
!!$       -79.681901957326275/)
!!$
!!$
!!$do i=1, size(bb)
!!$   r = bb(i)**ee(i)
!!$   !write(*,*) b(i), e(i), r
!!$end do
!!$
!!$
!!$
!!$
!!$write(*,*) 'ExpInt1DGrad2> RefPt = ', RefPt, ', inode = ', inode
!!$
!!$
!!$L = inode(1,2) / inode(2,2)
!!$do i= int(inode(2,1)-1000), 110000, 10
!!$
!!$   r = real(i, wp)
!!$   exponent = (inode(2,1)-r) / (inode(2,1)-inode(1,1))
!!$   E = L**exponent
!!$
!!$end do
!!$







if (abs(inode(1,2)) .lt. 1E-10_wp .and. abs(inode(2,2)) .lt. 1E-10_wp ) then
   ! N1 = N2 = 0  => N(z) = 0

   val     = 0.0_wp
   GradVal = 0.0_wp

else if (abs(inode(2,1)-inode(1,1)) < 1.0E-3_wp) then
   ! identical coordinates, interpolation not possible, return first value
   val     = inode(1,2)
   GradVal = 0.0_wp

else
   ! interpolate

   n = inode

   if (abs(n(1,2)) .lt. 1E-10_wp) then
      ! N1 = 0, N2 <> 0
      n(1,2) = n(2,2)/100.0_wp
   end if
   if (abs(n(2,2)) .lt. 1E-10_wp) then
      ! N2 = 0, N1 <> 0
      n(2,2) = n(1,2)/100.0_wp
   end if

!!$   !write(*,*) 'ExpInt1DGrad2> inode = ', n
!!$   write(*,*) 'ExpInt1DGrad2> RefPt = ', RefPt, ', inode = ', inode, char(10), &
!!$        '          n(2,1)-n(1,1) = ', n(2,1)-n(1,1), char(10), &
!!$        '          n(2,1)-RefPt  = ', n(2,1)-RefPt, char(10), &
!!$        '          n(1,2)/n(2,2) = ', n(1,2)/n(2,2), char(10), &
!!$        '          (n(2,1)-RefPt)/(n(2,1)-n(1,1))) = ', &
!!$                              (n(2,1)-RefPt) / (n(2,1)-n(1,1)), char(10)

!!$   exponent = (n(2,1)-RefPt) / (n(2,1)-n(1,1))
!!$   L = n(1,2) / n(2,2)
!!$   if (L < 0.0_wp) then
!!$      write(*,*) 'ExpInt1DGrad2> lllll L < 0 :', L
!!$   end if
!!$   E = abs(L)**exponent
!!$   !write(*,*) 'ExpInt1DGrad2>  L, exponent = ', L, exponent
!!$   val = n(2,2) * E
   !val = n(2,2) * (n(1,2)/n(2,2))**exponent

   val = n(2,2) * (n(1,2)/n(2,2))**((n(2,1)-RefPt)/(n(2,1)-n(1,1)))

   B = (n(2,1)-RefPt) / (n(2,1)-n(1,1))
   L = log( n(1,2)/n(2,2) )
   E = exp( B * L )
   NE = n(2,2) * E

   ! GradVal = d I / d RefPt  -  first derivative
   GradVal = -L/(n(2,1)-n(1,1)) * NE

   ! GradVal2 = d^2 I / d RefPt^2  -  second derivative
   dB = -1.0_wp / ( n(2,1)-n(1,1) )
   dE = L * dB * E
   GradVal2 = GradVal * dE/E

end if

end subroutine ExpInt1DGrad2
!****
! <<<<<<<< End of RoboDoc comments


!---------------------------------------------------------------------
! function Geopot2GeoidOld
!---------------------------------------------------------------------
!
!> Convert geopotential or geopotential height in height above geoid
!>
!> <b> AltGeoid = Geopot2Geoid( lat, hgeop ) </b>
!>
!> Function moved to oo-model/mo_physics.f90: geopot_geoid
!>
!> The geopotential is the potential energy due to gravity and therefore
!> related to the geoid or height above sea level.
!>
!> Geopotential: \f$ \Phi = \int_0^h g \, dh  \f$, g = g(h) - gravity \n
!> Geopotential height: \f$ h_{g} = \frac{1}{g_0} \int_0^h g \, dh \f$,
!>                      \f$g_0\f$ - standard gravity, \f$g_0 = 9.80665\f$
!>
!> Inverse function: ::Geoid2Geopot
!>
!> @param[in] lat
!>            latitude in radian \n
!> @param[in] geop
!>            geopotential or geopotential height in meter,
!>            depending on "geopot": \n
!>            geopot = .false. or missing: geop is geopotential height
!>                                                 in meter \n
!>            geopot = .true.            : geop is geopotential
!> @param[in] geopot
!>            logical giving the quantity in "geopot" \n
!>            optional
!> @return Geopot2Geoid - height above geoid in meter
!
!---------------------------------------------------------------------
! Record of revisions
! Date        Programmer   Description of change
! ====        ==========   =====================
! 16.12.2010  M. Bender    new, copy from Florian Zus
! 24.03.2015  M. Bender    extended to geopotential
!---------------------------------------------------------------------
function Geopot2GeoidOld(lat, geop, geopot)

! List of calling arguments:
real(wp)                       :: Geopot2GeoidOld
real(wp), intent(in)           :: geop
real(wp), intent(in)           :: lat
logical, optional,  intent(in) :: geopot

! List of local variables:
real(wp), parameter  :: g0 = 9.80665_wp
logical              :: IsPot

real(wp) :: rE
real(wp) :: gG
real(wp) :: SinLat2
!---------------------------------------------------------------------

! decide if "geop" is geopotential or geopotential height
if (present(geopot)) then
   IsPot = geopot
else
   ! assume "geop" is geopotential height
   IsPot = .false.
end if

SinLat2 = sin(lat)**2

! ellipsoidal radius of Earth, depending on latitude, radius in m
rE = 6378137.0_wp / (1.006803_wp - 0.006706_wp * SinLat2)

! local gravity on ellipsoid, depending on latitude
gG = 9.7803267714_wp *                            &
     (1.0_wp + 0.00193185138639_wp * SinLat2) /   &
     sqrt( 1.0_wp - 0.00669437999013_wp * SinLat2)

if (IsPot) then
   ! geopotential => height above geoid in m
   Geopot2GeoidOld = (rE * geop) / (gG * rE - geop)
else
   ! geopotential height in m => height above geoid in m
   Geopot2GeoidOld = (g0 * rE * geop) / (gG * rE - g0 * geop)
end if

end function Geopot2GeoidOld


!---------------------------------------------------------------------
! function Geoid2GeopotOld
!---------------------------------------------------------------------
!
!> Convert height above geoid to geopotential or geopotential height in m.
!>
!> <b> Hgp = Geoid2Geopot(lat, Hgeo) </b>
!>
!> Function moved to oo-model/mo_physics.f90: geoid_geopot
!>
!> The geopotential is the potential energy due to gravity and therefore
!> related to the geoid or height above sea level.
!>
!> Geopotential: \f$ \Phi = \int_0^h g \, dh  \f$, g = g(h) - gravity \n
!> Geopotential height: \f$ h_{g} = \frac{1}{g_0} \int_0^h g \, dh \f$,
!>                      \f$g_0\f$ - standard gravity, \f$g_0 = 9.80665\f$
!>
!> Inverse function: ::Geopot2Geoid
!>
!> \param[in] lat
!>            latitude in radian \n
!> \param[in] Hgeo
!>            height above geoid in meter
!> \param[in] geopot
!>            logical defining the output quantity "Geoid2Geopot" \n
!>            optional
!> \return    Geoid2Geopot -
!>            geopotential if geopot = .true. \n
!>            geopotential height in m if geopot = .false. or missing
!
!---------------------------------------------------------------------
! Record of revisions
! Date        Programmer   Description of change
! ====        ==========   =====================
! 24.03.2015  M. Bender    new, see "Geopot2Geoid"
!---------------------------------------------------------------------
function Geoid2GeopotOld(lat, Hgeo, geopot)

! List of calling arguments:
real(wp)                       :: Geoid2GeopotOld
real(wp), intent(in)           :: Hgeo
real(wp), intent(in)           :: lat
logical, optional,  intent(in) :: geopot

! List of local variables:
real(wp), parameter  :: g0 = 9.80665_wp
logical              :: IsPot

real(wp) :: rE
real(wp) :: gG
real(wp) :: SinLat2
!---------------------------------------------------------------------

! decide if "geop" is geopotential or geopotential height
if (present(geopot)) then
   IsPot = geopot
else
   ! assume "geop" is geopotential height
   IsPot = .false.
end if

SinLat2 = sin(lat)**2

! ellipsoidal radius of Earth, depending on latitude, radius in m
rE = 6378137.0_wp / (1.006803_wp - 0.006706_wp * SinLat2)

! local gravity on ellipsoid, depending on latitude
gG = 9.7803267714_wp *                            &
     (1.0_wp + 0.00193185138639_wp * SinLat2) /   &
     sqrt( 1.0_wp - 0.00669437999013_wp * SinLat2)

if (IsPot) then
   ! geopotential => height above geoid in m
   Geoid2GeopotOld =  (rE * Hgeo * gG) / (rE + Hgeo)
else
   ! geopotential height in m => height above geoid in m
   Geoid2GeopotOld = (rE * Hgeo * gG) / (g0*(rE + Hgeo))
end if

end function Geoid2GeopotOld


!---------------------------------------------------------------------
! function Shepard2D
!---------------------------------------------------------------------
!
!> 2D interpolation between 3 or more nodes using Shepard interpolation
!>
!> <b> val = Shepard2D(RefPt, node) </b>
!>
!> Some values, e.g. heights, refractivities, temperatures, are given at
!> the corners of a cell, e.g. a triangle or a square. These values are
!> interpolated at any given point within this cell (reference point).
!> For 2D interpolation each point has two coordiates and one value.
!> The Shepard interpolation is used, i.e. a weighted inverse distance
!> method.
!>
!> \param[in] RefPt
!>            coordinates of the reference point \n
!>            The value will be interpolated to this point. RefPt should
!>            be located within the points defined by <b>node</b>.
!> \param[in] node
!>            Grid nodes with coordinates and values. \n
!>            node(i,j)  \n
!>                 i=1, ... , 3 or 4 - node index describing 3 or 4 nodes \n
!>                 j=1, ... , 3  - coordinates an value at that point \n
!>                 j = 1 - longitude, j = 2 - latitude, j = 3 - value
!>
!> \return    Shepard2D
!>            interpolated value at RefPt
!
!---------------------------------------------------------------------
! Record of revisions
! Date        Programmer   Description of change
! ====        ==========   =====================
! 24.03.2015  M. Bender    new
!---------------------------------------------------------------------
function Shepard2D(RefPt, node)

! List of calling arguments:
real (wp)                             :: Shepard2D
real (wp), dimension(:), intent(in)   :: RefPt
real (wp), dimension(:,:), intent(in) :: node

! List of local variables:
real (wp), parameter :: ShepMu = 1.5_wp
real (wp), parameter :: inf = HUGE(Shepard2D)
real (wp), dimension(:), allocatable :: disti    ! inverse distances
real (wp) :: s1, s2
integer :: i
!---------------------------------------------------------------------

allocate( disti(lbound(node,1):ubound(node,1)) )

! compute inverse distances between all nodes and the reference point
disti = ( (RefPt(1)-node(:,1))**2 + (RefPt(2)-node(:,2))**2 )**(-ShepMu/2.0_wp)

if (any(disti > inf)) then
   ! reference point near one of the nodes, return value at that node
   do i=lbound(node,1), ubound(node,1)
      if (disti(i) > inf) then
         Shepard2D = node(i,3)
         exit
      end if
   end do
else
   ! interpolate, no singularity
   s1 = sum(disti)
   s2 = dot_product(node(:,3), disti)
   Shepard2D = s2 / s1
end if

end function Shepard2D


function Shepard2D_tl(RefPt, node, node_tl)

! List of calling arguments:
real (wp)                             :: Shepard2D_tl
real (wp), dimension(:), intent(in)   :: RefPt
real (wp), dimension(:,:), intent(in) :: node
real (wp), dimension(:,:), intent(in) :: node_tl

! List of local variables:
real (wp), parameter :: ShepMu = 1.5_wp
real (wp), parameter :: inf = HUGE(Shepard2D_tl)
real (wp), dimension(:), allocatable :: disti    ! inverse distances
real (wp) :: s1, s2
integer :: i
!---------------------------------------------------------------------

allocate( disti(lbound(node,1):ubound(node,1)) )

! compute inverse distances between all nodes and the reference point
disti = ( (RefPt(1)-node(:,1))**2 + (RefPt(2)-node(:,2))**2 )**(-ShepMu/2.0_wp)

if (any(disti > inf)) then
   ! reference point near one of the nodes, return value at that node
   do i=lbound(node,1), ubound(node,1)
      if (disti(i) > inf) then
         Shepard2D_tl = node_tl(i,3)
         exit
      end if
   end do
else
   ! interpolate, no singularity
   s1 = sum(disti)
   s2 = dot_product(node_tl(:,3), disti)
   Shepard2D_tl = s2 / s1
end if

end function Shepard2D_tl


subroutine Shepard2D_ad(RefPt, node, node_ad, Dipol)

! List of calling arguments:
real (wp), dimension(:), intent(in)    :: RefPt
real (wp), dimension(:,:), intent(in)  :: node
real (wp), dimension(:,:), intent(out) :: node_ad
real (wp),  intent(inout)              :: Dipol    ! input

! List of local variables:
real (wp), parameter :: ShepMu = 1.5_wp
real (wp), parameter :: inf = HUGE(Dipol)
real (wp), dimension(:), allocatable :: disti    ! inverse distances
real (wp) :: s1, Ds2
integer :: i
!---------------------------------------------------------------------

allocate( disti(lbound(node,1):ubound(node,1)) )

! compute inverse distances between all nodes and the reference point
disti = ( (RefPt(1)-node(:,1))**2 + (RefPt(2)-node(:,2))**2 )**(-ShepMu/2.0_wp)

! adjoint code
node_ad = 0.0_wp
if (any(disti > inf)) then
   ! reference point near one of the nodes, return value at that node
   do i=lbound(node,1), ubound(node,1)
      if (disti(i) > inf) then
         node_ad(i,3) = Dipol
         exit
      end if
   end do
else
   ! interpolate, no singularity
   s1 = sum(disti)
   Ds2 = Dipol/s1
   !node_ad(:,3) = disti(:) * Ds2
   node_ad(lbound(node,1):ubound(node,1),3) = disti(:) * Ds2
end if

end subroutine Shepard2D_ad

!---------------------------------------------------------------------
! function DistSphere
!---------------------------------------------------------------------
!
! Name
! DistSphere
!
! Call
! dist = DistSphere(pos1, pos2)
!
! Purpose
! Compute the minimum distance between two points on a sphere along
! the orthodrome.
! The angle between two points with longitude L1 and L2 and latitude
! F1 and F2 is given by
! phi = acos( sin(F1)*sin(F2) + cos(F1)*cos(F2)*cos(L2-L1) )
!       phi in radian, L1, L2, F1, F2 in radian
! phi = S/Re, Re = Earth's radius, S = distance between points
! S = phi*Re
!
! Input
! pos1   - position 1, longitude, latitude in rad
! pos2   - position 2, longitude, latitude in rad
!          The third array element pos(3) is usually the hight,
!          but will not be used by this routine.
!          pos1(1) = longitude in rad
!          pos1(2) = latitude in rad
!          pos1(3) = arbitrary value, not used by this function
!
! Output
! DistSphere - minimum distance between two points on a sphere in km
!              (return value)
!
!---------------------------------------------------------------------
function DistSphere (pos1, pos2)

! List of calling arguments:
real (wp)                 :: DistSphere
real (wp), dimension(1:3) :: pos1, pos2

! List of local variables:
real (wp)             :: r
real (wp), parameter  :: rE = 6378137.0_wp ! radius of Earth in m
!---------------------------------------------------------------------

r = sin(pos1(2)) * sin(pos2(2)) +     &
    cos(pos1(2)) *  cos(pos2(2)) * cos(pos2(1) - pos1(1))
r = acos(r)

! DistSphere = r * 6371.0        ! Distance in km
DistSphere = r * rE/1000.0_wp     ! Distance in km

end function DistSphere

!---------------------------------------------------------------------
End Module  mo_std_coord
!---------------------------------------------------------------------
