!
!+ Provide constants for WMO table entries used for GRIB encoding.
!
MODULE mo_wmo_tables
!
! Description:
!   Provide constants for WMO table entrues:
!     - Table 0 (C-1): center identifier
!     - Table 3      : level type
!     - Table 4      : unit of time
!     - Table 5      : range of time
!     - Table 6      : grid type
!     - Table 8      : scanning mode
!     - Table A      : BUFR data category
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
! V1_4         2009/03/26 Andreas Rhodin
!  define unit-of-time=13 (15 minutes)
! V1_7         2009/08/24 Andreas Rhodin
!  add DWD specific time range indicator values: 13=nudging, 15=EnKF analysis
! V1_8         2009/12/09 Andreas Rhodin
!  define: WMO0_CIMSS=176, WMO3_HEIGHT=103, DWD5_IFS_FC=14
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Harald Anlauf
!  add level type 102 (sea level), generating centers NCAR,COSMO
! V1_22        2013-02-13 Harald Anlauf
!  add generating centers DMI, KNMI, GPSRO; entries for ICON grid
! V1_23        2013-03-26 Andreas Rhodin
!  account for DWD time_range indicator DWD5_INIT_FC=11 (initialised forecast)
! V1_26        2013/06/27 Harald Anlauf
!  Add WMO definition for generalized vertical coordinate
!  Add WMO code for MeteoSwiss
! V1_27        2013-11-08 Andreas Rhodin
!  add GRIB1 codes for ECHAM
! V1_42        2015-06-08 Andreas Rhodin
!  define: DWD6_NONE = 999 ! no grid, just a collection of columns
! V1_43        2015-08-19 Harald Anlauf
!  Improve GRIB encoding of single-level/surface fields for flake
! V1_48        2016-10-06 Harald Anlauf
!  Add CMA as processing center for GPSRO, handle FY-3C/GNOS
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Andreas Rhodin  DWD  2003-2008
!------------------------------------------------------------------------------
  !-------------
  ! Modules used
  !-------------
  USE mo_exception, ONLY: finish
  IMPLICIT NONE

  PUBLIC

  !-------------------------------
  ! WMO table 0: center identifier
  ! (common code table C-1)
  !-------------------------------
  INTEGER, PARAMETER ::   &
    WMO0_NCEP     =   7,  & ! US    NCEP
    WMO0_NWSTG    =   8,  & ! US    NWS Telecommunications Gateway
    WMO0_NWS      =   9,  & ! US    NWS
    WMO0_IMD      =  28,  & ! IMD   New Delhi
    WMO0_JMA      =  34,  & ! JMA   Tokyo
    WMO0_CMA      =  38,  & ! CMA   Beijing
    WMO0_RSMC     =  39,  & ! RSMC  Beijing
    WMO0_KMA      =  40,  & ! KMA   Seoul
    WMO0_MSC      =  53,  & ! MSC   Montreal
    WMO0_MSC2     =  54,  & ! MSC   Montreal
    WMO0_ARINC    =  56,  & ! ARINC
    WMO0_NCAR     =  60,  & ! NCAR  US National Center for Atmospheric Research
    WMO0_UKMET    =  74,  & ! UK    Exeter
    WMO0_DWD      =  78,  & ! DWD   Offenbach
    WMO0_COMET    =  80,  & ! COMET Rome
    WMO0_DMI      =  94,  & ! DMI   Kopenhagen
    WMO0_ECMWF    =  98,  & ! ECMWF Reading
    WMO0_KNMI     =  99,  & ! KNMI  De Bilt
    WMO0_HKO      = 110,  & ! HKO   Hong-Kong
    WMO0_NOAA     = 160,  & ! NOAA  US NOAA/NESDIS
    WMO0_NASA     = 173,  & ! NASA  US NASA
    WMO0_UCAR     = 175,  & ! UCAR  University Corporation for Atm. Research
    WMO0_CIMSS    = 176,  & ! Cooperative Institute for Meteorol. Satellite
                            !   Studies, University of Wisconsin-Madison
    WMO0_SPIRE    = 178,  & ! Spire     (commercial)
    WMO0_GEOOPT   = 179,  & ! GeoOptics (commercial)
    WMO0_PLANET   = 180,  & ! PlanetiQ  (commercial)
    WMO0_MSWISS   = 215,  & ! MeteoSwiss Zurich
    WMO0_COSMO    = 250,  & ! COnsortium for Small scale MOdelling
    WMO0_MPIFM    = 252,  & ! Max Planck Institute for Meteorology
    WMO0_EUMET    = 254     ! EUMET EUMETSAT Operation Centre

  !------------------------
  ! WMO table 3: level type
  !------------------------
  INTEGER, PARAMETER ::   &
    WMO3_ISOBARIC  = 100, & ! isobaric level              [hPa]
    WMO3_BELOWSUR  = 111, & ! depth below surface         [cm]
    WMO3_SURFACE   =   1, & ! ground or water surface
    WMO3_CLD_BASE  =   2, & ! cloud base level
    WMO3_CLD_TOPS  =   3, & ! level of cloud tops
    WMO3_ZERO_ISO  =   4, & ! level of 0C isotherm
    WMO3_NOM_TOA   =   8, & ! nominal top of the atmosphere
    WMO3_ISOTHERM  =  20, & ! Isothermal level
    WMO3_SEALEVEL  = 102, & ! mean sea level
    WMO3_HEIGHT    = 103, & ! fixed height above sealevel [m]
    WMO3_ABOVESUR  = 105, & ! level above ground          [m]
    WMO3_LAYER     = 106, & ! layer b.levels above ground [m]
    WMO3_HYBRID    = 109, & ! hybrid levels
    WMO3_HYBRIDB   = 110, & ! level between 2 hybrid levels
    WMO3_HHYBRID   = 118, & ! height-based hybrid levels
    WMO3_GENV      = 150, & ! generalized vertical coordinate
    WMO3_LAKE_BOT  = 162, & ! lake or river bottom
    WMO3_SEDIM_BOT = 165, & ! bottom of sediment layer
    WMO3_MIX_LAYER = 166, & ! mixed layer
    WMO3_TILE_LAND = 181, & ! grid-tile land fraction as a model surface
    WMO3_TILE_SEA  = 182, & ! grid-tile water fraction as a model surface
    WMO3_TILE_SICE = 183, & ! grid-tile ice fraction on sea, lake or river
    WMO3_TILE_LICE = 184, & ! grid-tile glacier ice and inland ice fraction
    WMO3_SURF_HORZ = 208, & ! horizontal plane         (incl.orographic shading)
    WMO3_SURF_TANG = 209    ! tangent plane to terrain (incl.orographic shading)

  !--------------------------
  ! WMO table 4: unit of time
  !--------------------------
  INTEGER, PARAMETER ::       &
    WMO4_MINUTE     =   0,    &
    WMO4_HOUR       =   1,    &
    WMO4_DAY        =   2,    &
    WMO4_MONTH      =   3,    &
    WMO4_YEAR       =   4,    &
    WMO4_DECADE     =   5,    &
    WMO4_NORMAL     =   6,    &
    WMO4_CENTURY    =   7,    &
    WMO4_3_HOURS    =  10,    &
    WMO4_6_HOURS    =  11,    &
    WMO4_12_HOURS   =  12,    &
    DWD4_15_MINUTES =  13,    &
    WMO4_SECOND     = 254

  !---------------------------
  ! WMO table 5: range of time
  !---------------------------
  INTEGER, PARAMETER ::     &
    WMO5_FORECAST  =   0, &! Forecast valid for reference time + P1
    WMO5_INIT_ANA  =   1, &! Initialized analysis (P1=0).
    WMO5_RANGE     =   2, &! Range between ref.time+P1 and ref.time+P2
    WMO5_AVERAGE   =   3, &! Average (ref.time+P1 to ref.time+P2)
    WMO5_ACCU      =   4, &! Accumulation (ref.time+P1 to ref.time+P2)
    WMO5_DIFF      =   5, &! Difference (ref.time+P1 minus ref.time+P2)
    WMO5_VALID     =  10, &! Valid at ref.time+P1 (P1 occupies 2 octets)
    DWD5_INIT_FC   =  11, &! Initialized forecast
    DWD5_NUDGING   =  13, &! nudging analysis (DWD special)
    DWD5_IFS_FC    =  14, &! IFS forecast
    DWD5_ENKF_ANA  =  15, &! EnKF    analysis (DWD special)
    WMO5_FC_AV_1   = 113, &! Average of N forecasts.
    WMO5_FC_AC_1   = 114, &! Accumulation of N forecasts.
    WMO5_FC_AV_R   = 115, &! Average of N forecasts, with the same ref.time.
    WMO5_FC_AC_R   = 116, &! Accumulation of N forecasts, with same ref.time.
    WMO5_FC_AV_2   = 117, &! Average of N forecasts.
    WMO5_VAR_IANA  = 118, &! Temporal variance/covariance, of N init.analyses.
    WMO5_STDEV_FC  = 119, &! Standard deviation of N forecasts.
    WMO5_AV_IANA   = 123, &! Average of N uninitialized analyses.
    WMO5_AC_IANA   = 124   ! Accumulation of N uninitialized analyses.

  !------------------------------------------------
  ! WMO table 6: gridtype (GDS octet 6) code values
  !------------------------------------------------
  INTEGER, PARAMETER :: &
    WMO6_LATLON      =   0, & ! Latitude/Longitude Grid
    WMO6_MERCATOR    =   1, & ! Mercator projection grid
    WMO6_GNOMONIC    =   2, & ! Gnomonic projection grid
    WMO6_LAMBERT     =   3, & ! Normal Lambert Projection Grid
    WMO6_GAUSSIAN    =   4, & ! Gaussian latitude/longitude grid
    WMO6_POLAR       =   5, & ! Polar stereographic projection grid
    WMO6_ROTLL       =  10, & ! Rotated latitude/longitude grid
    WMO6_OBLIQUE     =  13, & ! Oblique Lambert Projection Grid
    WMO6_HARMONIC    =  50, & ! Spherical harmonic coefficients
    WMO6_SPACEVIEW   =  90, & ! Space view perspective or orthographic grid
    DWD6_ICOSAHEDRON = 192, & ! Icosahedral based triangular grid
    DWD6_ICON        = 193, & ! Temporary value for testing ICON
    DWD6_NONE        = 999    ! no grid, just a collection of columns

!--------------------------------------------------------------------------
!                   TABLE 7 - RESOLUTION AND COMPONENT FLAGS
!                            (GDS Octet 17)
!
!  Bit        Value        Meaning
!  1          0            Direction increments not given
!             1            Direction increments given
!  2          0            Earth assumed spherical with radius = 6367.47 km
!             1            Earth assumed oblate spheroid with size
!                          as determined by IAU in 1965:
!                          6378.160 km, 6356.775 km, f = 1/297.0
!  3-4                     reserved (set to 0)
!  5          0            u- and v-components of vector quantities resolved
!                          relative to easterly and northerly directions
!             1            u and v components of vector quantities resolved
!                          relative to the defined grid in the
!                          direction of increasing x and y (or i and j)
!                          coordinates respectively
!  6-8                     reserved (set to 0)
!---------------------------------------------------------------------------

  !-----------------------------------------------
  ! WMO table 8: scanning mode flag (GDS Octet 28)
  !-----------------------------------------------
  INTEGER, PARAMETER :: &
    WMO8_I_NEGATIVE    = 128, &! Bit 1: Points scan in -i direction
    WMO8_J_POSITIVE    =  64, &! Bit 2: Points scan in +j direction
    WMO8_J_CONSECUTIVE =  32   ! Bit 3: Adjacent points in j direction are
                               !        consecutive (FORTRAN: (J,I))
                               ! Note: i direction is defined as west to east
                               !       j direction is defined as south to north

!----------------------------------------------------------------------
!                 TABLE 9. SPECTRAL REPRESENTATION TYPE
!                            (GDS Octet 13)
!
!  VALUE       MEANING
!  1           Associated Legendre Polynomials of the First Kind with
!              normalization such that the integral equals 1
!
!                      TABLE 10. COEFFICIENT STORAGE MODE
!                             (GDS Octet 14)
!
!  VALUE        MEANING
!  1            The complex coefficients Xnm are stored for m > 0 as
!               pairs of real numbers Re(Xnm), Im(Xnm)
!               ordered with n increasing from m to N(m), first for m =
!               0 and then for m = 1, 2, 3,...M. The real part of the
!               (0,0) coefficient is stored in octets 12-15 of the BDS,
!               as a floating point number in the same manner
!               as the packing reference value, with units as in Table
!               2. The remaining coefficients, starting with the
!               imaginary part of the (0,0) coefficient, are packed
!               according to the GRIB packing algorithm, with units
!               as given in Table 5, in octets 16 and onward in the BDS.
!------------------------------------------------------------------------

  !----------------------------------
  ! WMO table A (BUFR): Data Category
  !----------------------------------
  INTEGER, PARAMETER :: &
    WMOA_SURF_LAND  =   0, &! Surface data - land
    WMOA_SURF_SEA   =   1, &! Surface data - sea
    WMOA_VERT_SOUND =   2, &! Vertical soundings (other than satellite)
    WMOA_VERT_S_SAT =   3, &! Vertical soundings (satellite)
    WMOA_SINGL_LEV  =   4, &! Single level upper-air data(other than satellite)
    WMOA_SINGL_SAT  =   5, &! Single level upper-air data (satellite)
    WMOA_RADAR      =   6, &! Radar data
    WMOA_SYNOPTIC   =   7, &! Synoptic data
    WMOA_CHEM       =   8, &! Physical/chemical constituents
    WMOA_DISP_TRANS =   9, &! Dispersal and transport
    WMOA_RADIO      =  10, &! Radiological data
    WMOA_REPL       =  11, &! BUFR tables, complete replacement or update
    WMOA_SURF_SAT   =  12, &! Surface data (satellite)
    WMOA_STATUS     =  20, &! Status information
    WMOA_RAD        =  21, &! Radiances (satellite)
    WMOA_RADAR_SAT  =  22, &! Radar (satellite) but not altimeter, scatterometer
    WMOA_LIDAR_SAT  =  23, &! Lidar (satellite)
    WMOA_SCATT      =  24, &! Scatterometry (satellite)
    WMOA_ALTIM      =  25, &! Altimetry (satellite)
    WMOA_OCEAN      =  31, &! Oceanographic data
    WMOA_IMAGE      = 101, &! Image data
    WMOA_LOCAL      = 255   ! Indicator for local use, with sub-category

contains

  subroutine wmo_mnem (table, entry, mnemonic, comment)
  character         ,intent(in)   :: table
  integer           ,intent(in)   :: entry
  character (len=*) ,intent(out)  :: mnemonic
  character (len=*) ,intent(out)  :: comment
    select case (table)
    case ('A')
      select case (entry)
      case (WMOA_SURF_LAND)
        mnemonic = 'SURF_LAND'
        comment  = 'Surface data - land'
      case (WMOA_SURF_SEA)
        mnemonic = 'SURF_SEA'
        comment  = 'Surface data - sea'
      case (WMOA_VERT_SOUND)
        mnemonic = 'VERT_SOUND'
        comment  = 'Vertical soundings (other than satellite)'
      case (WMOA_VERT_S_SAT)
        mnemonic = 'VERT_S_SAT'
        comment  = 'Vertical soundings (satellite)'
      case (WMOA_SINGL_LEV)
        mnemonic = 'SINGL_LEV'
        comment  = 'Single level upper-air data (other than satellite)'
      case (WMOA_SINGL_SAT)
        mnemonic = 'SINGL_SAT'
        comment  = 'Single level upper-air data (satellite)'
      case (WMOA_RADAR)
        mnemonic = 'RADAR'
        comment  = 'Radar data'
      case (WMOA_SYNOPTIC)
        mnemonic = 'SYNOPTIC'
        comment  = 'Synoptic data'
      case (WMOA_CHEM)
        mnemonic = 'CHEM'
        comment  = 'Physical/chemical constituents'
      case (WMOA_DISP_TRANS)
        mnemonic = 'DISP_TRANS'
        comment  = 'Dispersal and transport'
      case (WMOA_RADIO)
        mnemonic = 'RADIO'
        comment  = 'Radiological data'
      case (WMOA_REPL)
        mnemonic = 'REPL'
        comment  = 'BUFR tables, complete replacement or update'
      case (WMOA_SURF_SAT)
        mnemonic = 'SURF_SAT'
        comment  = 'Surface data (satellite)'
      case (WMOA_STATUS)
        mnemonic = 'STATUS'
        comment  = 'Status information'
      case (WMOA_RAD)
        mnemonic = 'RAD'
        comment  = 'Radiances'
      case (WMOA_RADAR_SAT)
        mnemonic = 'RADAR_SAT'
        comment  = 'Radar (satellite)'
      case (WMOA_LIDAR_SAT)
        mnemonic = 'LIDAR'
        comment  = 'Lidar (satellite)'
      case (WMOA_SCATT)
        mnemonic = 'SCATT'
        comment  = 'Scatterometry (satellite)'
      case (WMOA_ALTIM)
        mnemonic = 'ALTIM'
        comment  = 'Altimetry (satellite)'
      case (WMOA_OCEAN)
        mnemonic = 'OCEAN'
        comment  = 'Oceanographic data'
      case (WMOA_IMAGE)
        mnemonic = 'IMAGE'
        comment  = 'Image data'
      case (WMOA_LOCAL)
        mnemonic = 'LOCAL'
        comment  = 'Indicator for local use, with sub-category'
      case (240:254)
        mnemonic = 'EXPERIMENTAL'
        comment  = 'For experimental use'
      case default
        mnemonic = 'RESERVED'
        comment  = 'Reserved'
      end select
    case default
      call finish ('wmo_mnem','invalid table: '//table)
    end select
  end subroutine wmo_mnem

END MODULE mo_wmo_tables
