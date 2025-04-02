!
!+ core Slant Total Delay operator
!
module mo_std_operator
!
! Description:
!   core Slant Total Delay operator
!
! Current Maintainer: DWD, Michael Bender, Harald Anlauf
!    phone: +49 69 8062 4941
!    fax:   +49 69 8062 3721
!    email: harald.anlauf@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_22        2013-02-13 Michael Bender
!  core Slant Total Delay operator
! V1_26        2013/06/27 Andreas Rhodin
!  template for adjoint code
! V1_27        2013-11-08 Michael Bender
!  tl/adjoint routines
! V1_29        2014/04/02 Andreas Rhodin
!  changes for ZTD assimilation
! V1_42        2015-06-08 Andreas Rhodin
!  bias correction for ground based GNSS data
! V1_43        2015-08-19 Michael Bender
!  Raytracer added to STD operator; GNSS bias correction.
! V1_45        2015-12-15 Michael Bender
!  adjoint code fixed; Make R(dry)/R(vapor) consistent
! V1_46        2016-02-05 Michael Bender
!  some routines moved to mo_physics and mo_algorithms; disable test mode
! V1_47        2016-06-06 Michael Bender
!  old interpolation option removed; 'gpm' is geometric height, not geopot.
! V1_48        2016-10-06 Michael Bender
!  add recompiler commands for COSMO; define different verbosity levels
! V1_50        2017-01-09 Michael Bender
!  extended namelist STD_OBS; modules merged with COSMO code.
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
!==============================================================================
!
!> @file mo_std_operator.f90
!> Main module of the GNSS STD observation operator
!>
!> STD - Slant Total Delay
!> GNSS - Global Navigation Satellite System, like GPS, GLONASS, BeiDou, ...
!>
! COSMO
! Suche nach Gitterzellen optimieren
! Regulaeres Gitter fuer grobe Suche definieren, danach verfeinern
! siehe
! 3dvar/oo-model/mo_icon_grid.f90, subroutine set_search_grid
!
! Operator 2. Teil:
! Test einbauen, ob der Slant sich streckenweise unter der Modellorographie
! bewegt, wenn ja verwerfen.
!
!==============================================================================
! Modules used

#ifdef __COSMO__
  use kind_parameters, only: wp,              &! working precision kind
                             dp,              &! double  precision kind
                             sp,              &! single  precision kind
                             i2                ! 2-byte integer kind parameter
#else
  use mo_kind,         only: wp,              &! working precision kind
                             dp,              &! double  precision kind
                             sp,              &! single  precision kind
                             i2                ! 2-byte integer kind parameter
  use mo_mpi_dace,     only: dace              ! MPI group info
#endif

  use environment,     only: model_abort       ! abort in case of error
! use mo_dace_string,  only: char1,           &! conversion: integer -> char(1)
!                            char3,           &! conversion: integer -> char(3)
!                            split,           &! char string -> array
!                            eval_string       ! evaluate strings

#ifndef __COSMO__
  use mo_system,       only: flush             ! Flush I/O buffer
  use mo_run_params,   only: ana_time,        &! analysis date
                             flag_biasc_gpsgb  ! steering of bias correction
  use mo_time,         only: i_time            ! t_time->yyyy,mm,dd,hh,mi,ss
#endif

#ifndef __COSMO__
  use mo_atm_grid,    only:  &
       MO_UNKNOWN,  & ! model:
       MO_HRM,      & ! HRM
       MO_GME,      & ! GME
       MO_COSMO,    & ! COSMO
       MO_ICON,     & ! ICON
       MO_IFS         ! IFS
!
! For COSMO some parameters need to be defined (see below):
!  integer, parameter :: MO_UNKNOWN = 0  ! model:
!  integer, parameter :: MO_COSMO   = 3  ! COSMO
!  integer, parameter :: MO_ICON    = 4  ! ICON
#endif

  use utilities,       only: phi2phirot,        &! ...
                             rla2rlarot          ! ...
!                            rlarot2rla,        &! ...
!                            phirot2phi

#ifdef __COSMO__
  use mo_std_coord, only : &
       IntegPolyCube,      &! num. integration
       IntegPolyCube_tl,   &! num. integration, tl-code
       IntegPolyCube_ad     ! num. integration, adjoint code
#else
  use mo_algorithms, only : &
       IntegPolyCube,     &! num. integration
       IntegPolyCube_tl,  &! num. integration, tl-code
       IntegPolyCube_ad    ! num. integration, adjoint code
#endif

#ifndef __COSMO__
  use MSIS,            only: MSIS_Init,         &! Initialize MSIS module
                             MSIS_Refractivity   ! compute refractivity
  use Earth,           only: Geodetic            ! type, geod. coordinates
#endif

  use mo_std_coord,    only: deg2rad,           &! degrees to radian
                             rad2deg,           &! radian to degrees
                             pi05,              &! pi/2
                             k1, k2, k3,        &! constants for Thayer refrac.
                             RDRD,              &! R(dry)/R(vapor)
                             EMRDRD,            &! 1._wp - RDRD
                             invalid,           &! invalid data
                             GNSSStation,       &! type, GNSS station data
                             SlantTotalDelay,   &! type, STD observations
                             SlantHeader,       &! type, header STD files
                             GPSStation,        &! type,  GNSS station data
                             WGS84Param,        &! type WGS84 reference ellips.
                             StrucDomes,        &! type GNSS domes
                             ImportStakoSNXarr, &! read GNSS station coord.
                             ImportESTDarr,     &! read STD data
!                            ImportDomes,       &! read GNSS station location
                             Ellips2Cart,       &! ellips. to cartesian
                             Cart2Ellips,       &! cartesian to ellips.
                             LocalHorz2Cart,    &! local to cartesian
                             LocalHorz2CartMat, &! LocalHorz2Cart with matrix
                             Azimut2Horz,       &! azimuth + elev. to local
                             CrossLineEllips,   &! points on line and ellips.
!                            BiLinExp3D,        &! interpolation
                             BiLinearExp3D,     &! interpolation
                             BiLinearExp3D_tl,  &! interpolation, TL code
                             BiLinearExp3D_ad,  &! interpolation, adjoint code
                             Shepard2D,         &! Shepard interpolation
                             ExpInt1D,          &! vertical interpolation
                             ExpInt1D_tl,       &! vertical interpolation, TL
                             ExpInt1D_ad,       &! vertical interpolation, adj.
                             ExpInt1DGrad,      &! first derivatives of ExpInt1D
                             ExpInt1DGrad2,     &! second derivatives of ExpInt1D
!                            simpned,           &! num. integration
                             UpdateK123,        &! Change ki to namelist input
                             gmf,               &! Global Mapping Function GMF
                             JulianDate,        &! computes the Julian Date JD
                             GregorianDate,     &! computes the Gregorian Date
                             JD2MJD,            &! converts JD to MJD
                             MJD2JD,            &! converts MJD to JD
                             GREG2DOY,          &! computes day of year (DOY)
                             DOY2GREG,          &! comp. Grgorian date from DOY
                             NWein,             &! Refractivity
                             NWein_tl,          &! Refractivity, , TL code
                             NWein_ad,          &! Refractivity, , adjoint code
                             LayerSearch,       &! find vertical grid index k
                             SatPos,            &! comp. satellite position
                             Slant2Ellips,      &! transform slant to ellips.
!                            Slant2Ellips_tl,   &! transform slant to ellips.
                             Slant2EllipsMat,   &! Slant2Ellips with trans matrix
                             Slant2EllipsMatGrad,  &
                             Slant2EllipsMatGrad2, &
                             Slant2LocalHorzMat

  use mo_std_vector,   only: line,              &! straight line derived type
                             linedef,           &! define straight line
                             linepos,           &! coordinates of point on line
                             PointDist           ! distance between two points

#ifndef __COSMO__
  use mo_gnss_bc      ,only: biascor_mode,      &! mode used for bias correction
                             t_decay,           &! biasc. accum. decay time, days
                             n_required,        &! number of entries required
                                                 ! for bc
                             bc_fallback,       &! action if biasc-file not
                                                 ! present
                             BC_NOBC, BC_FG      ! values for 'biascor_mode'
#endif

  implicit none

#ifdef __COSMO__
  ! copy from mo_atm_grid:
  integer, parameter :: MO_UNKNOWN = 0  ! model:
  integer, parameter :: MO_COSMO   = 3  ! COSMO
  integer, parameter :: MO_ICON    = 4  ! ICON
#endif

  !================
  ! Public entities
  !================
  private
  !-------------------------
  ! derived type definitions
  !-------------------------
  public :: SlantData          ! slant information derived type
  public :: p_column           ! pointers to model columns
  !---------
  ! routines
  !---------
  public :: read_std_nml       ! read namelist
  public :: read_std_nml_cosmo ! read namelist within COSMO
  public :: read_gnss_stations ! read GNSS station information
  public :: scan_std_obs       ! read the observation file header
  public :: read_std_obs       ! read the observations
  public :: std_path_coord     ! compute supporting points on the slant path
  public :: std_delay          ! compute slant path delay
  public :: std_path_coord_line
  public :: std_delay_line
  public :: ztd_delay          ! optimized ZTD operator
! public :: std_validate       ! compare observations with model equivalents
  public :: destruct           ! deallocate components of type SlantData
  public :: GetFreeUnit        ! find free unit for opening file
  public :: Coord2Cell         ! find COSMO indices
  public :: RefracModel
  public :: RefracModelGrad
  public :: RefracModelGrad2
  public :: RefracModelOld
  public :: RefracModelGradOld
  public :: RefracModelGrad2Old
  !----------------------------------------
  ! namelist parameters (to be broadcasted)
  !----------------------------------------
  public :: mf                 ! string length in namelist
  public :: Nintervals         ! array length used in namelist
  public :: NStepVertMod       ! Number of vertical points inside the model
  public :: NStepVertTop       ! Number of vertical points above the model
  public :: NStepOpt           ! option for scaling points
  public :: Hmax               ! Maximum height for STD integration in m
  public :: HScaleP            ! air pressure scale height in model
  public :: HScaleP2           ! air pressure scale height above model
  public :: Hlevel             ! model levels per interval
  public :: Heights            ! height of model level at upper interval bound.
  public :: Hpoints            ! supporting points per interval
  public :: UseRaytracer       ! use raytracer for elev. below "UseRaytracer"
  public :: ZTDerror           ! ZTD observation error, m
  public :: Href               ! observation reference height (LETKF only)
  public :: ztd_col            ! use 1 model column for ZTDs

  public :: ZTDminUse          ! minimum ZTD used
  public :: ZTDmaxUse          ! maximum ZTD used
  public :: StatBelowSurface   ! reject stations below model surface
  public :: StatAboveSurface   ! reject stations above model surface
  public :: StatBelowColumn    ! reject stations below max. column

  public :: MaxSTDobs          ! max. number of STDs to read/process
  public :: MaxZTDobs          ! max. number of ZTDs to read/process

  public :: std_obs_file       ! STD observation input files
  public :: GeoidFile          ! Path + file name of geoid
  public :: HORIFile           ! GNSS station file
  public :: DomesFile          ! GNSS station locations
  ! ...
  public :: read_ascii
  public :: verbose            ! Verbosity level (0:quiet, >1:debug)
  public :: STDfile
  public :: StartMJD
  public :: EndMJD
  public :: SlantPath
  public :: pl_method          ! method to estimate vertical location
  !----------------------
  ! direct access to data
  !----------------------
  public :: STD                ! Collection of the slant information
  public :: STDobs             ! STDobs as part of STD
  public :: NSTDobs            ! number of STD observations in array STD
  public :: NSTDobsLen         ! Byte length of type "SlantTotalDelay"
  public :: CoordList          ! list of stations
  public :: GNSShash           ! station hash table (Id -> index)
  !-------------------------------------
  ! required by COSMO operator interface
  !-------------------------------------
  public :: MO_COSMO           ! COSMO model ID from mo_atm_grid

  !-----------
  ! interfaces
  !-----------
  interface destruct                     ! deallocate components of ..
    module procedure destruct_SlantData  ! .. derived type SlantData
  end interface destruct

  !interface std_path_coord_1
  !   module procedure std_path_coord, std_path_coord_line
  !end interface std_path_coord_1

  !interface std_delay_1
  !   module procedure std_delay_line, std_delay_ray3d
  !end interface std_delay_1

  !------------------------------------------
  ! Local data structures of the STD operator
  !------------------------------------------

  ! GNSS station coordinates
#ifdef __COSMO__
  type (GNSSStation), dimension(:), allocatable, target :: CoordList
#else
  type (GNSSStation), dimension(:), pointer  :: CoordList => NULL()
#endif
  integer, dimension(:), pointer             :: GNSShash => NULL()
  ! Hier muss eingetragen werden, fuer welchen Zeitraum die Koordinaten
  ! gelesen wurden. Wenn die STDs nicht zu diesem Zeitraum passen, muessen
  ! die Koordinaten neu eingelesen werden!
! real (wp), dimension(1:2)                  :: GnssStationValid


  !---------------------------------------------------------------------
  ! structure SlantData
  !---------------------------------------------------------------------
  !>
  !> @brief Collection of the slant information required by the STD operator.
  !>
  !> One array element of this type provides all GNSS data to carry out
  !> all computations for one single slant.
  !>
  ! --------------------------------------------------------------------
  !> @verbatim
  !> Slant       - slant delay data read from STD file
  !> Station     - GNSS station information, pointer to the corresponding
  !>               array element of "CoordList"
  !> SlantSteps      - distance to receiver in m of the supporting points,
  !>                   required for numerical integration along the slant
  !>                   path
  !>                   SlantSteps(1) = 0.0  -  receiver
  !>                   SlantSteps(i),   i = 1 ... Ntot
  !> SlantPointsCart - cartesian coordinates of the supporting points
  !>                   along the connecting line
  !>                   SlantPointsCart(i,1) - X (m)
  !>                   SlantPointsCart(i,2) - Y (m)
  !>                   SlantPointsCart(i,3) - Z (m)
  !>                   i - index of supporting point, i = 1 ... Ntot
  !> SlantPointsEll  - ellipsoidal coordinates of the supporting points
  !>                   along the connecting line
  !>                   path (lat/lon/alt)
  !>                   SlantPointsEll(i,1) - longitude (radian)
  !>                   SlantPointsEll(i,2) - latitude (radian)
  !>                   SlantPointsEll(i,3) - altitude (m) above ellipsoid
  !>                   i - index of supporting point, i = 1 ... Ntot
  !> SlantPointsGeo  - geoid coordinates of the supporting points
  !>                   along the connecting line
  !>                   => used by COSMO
  !>                   path (lat/lon/alt)
  !>                   SlantPointsGeo(i,1) - longitude (radian)
  !>                   SlantPointsGeo(i,2) - latitude (radian)
  !>                   SlantPointsGeo(i,3) - altitude (m) above geoid
  !>                                              (= above sea level)
  !>                   i - index of supporting point, i = 1 ... Ntot
  !> SlantPointsIdx  - index of the COSMO grid cell (lower left corner)
  !>                   containing the supporting point
  !>                   SlantPointsIdx(i,1-3) = (i,j,k)
  !> SlantPath       - true bended signal path estimated by the raytracer
  !>                   "SlantPath" provides the y and z coordinates, the x
  !>                   coordinates are in "SlantSteps".
  !>                   Signal path in the "slant system":
  !>                   (x,y,z)_i = ( SlantSteps(i), SlantPath(1,i),
  !>                                                SlantPath(2,i)  )
  !> SlantRefrac     - refractivity on the bended slant path
  !>                   (Same dimension and index as SlantSteps)
  !> LineRefrac      - refractivity at the supporting points on satellite-
  !>                   receiver axis (straight line)
  !>                   (Same dimension and index as SlantSteps)
  !>
  !> GeoidCorr       - Geoid undulation
  !>                   Identical to Station%GeoidCorr but pointer
  !>                   "Station" is not allocated while reading observations
  !> ExtraPolNode    - refractivity and height at the last two grid layers
  !>                   on the slant path, required for vertical extraplation
  !>                   of the refractivity profile
  !>
  !> Ntot        - total number of supporting points along the slant path
  !>               Ntot = Nmod + Nup
  !> Nmod        - number of supporting points inside the model grid, i.e.
  !>               below the model top layer
  !> Nup         - number of supporting points above the model top layer, i.e.
  !>               above the model top layer
  !>
  !> Nhor        - number of supporting points inside the model grid (<=Nmod)
  !> Naloc       - number of points used (either Nmod or 0 if Nhor<Nmod)
  !>
  !> htop        - model top, required for second call to std_path_coord
  !> dx          - horizontal grid spacing, required for second call ""
  !>
  !> idxi        - index to supporting model column in 3dvar derived type
  !>               't_cols'.
  !> wi          - horizontal interpolation weights
  !> assimilate  - assimilate this slant: .true./.false.
  !> STD         - slant total dely computed with raytracer along curved path
  !> STDline     - slant total dely computed along straight line
  !>               STDline = -1 if not computed
  !> @endverbatim
  ! --------------------------------------------------------------------
  type SlantData
     !> STD observation and meta data
     type (SlantTotalDelay)             :: Obs
     !> GNSS station data
     type (GNSSStation), pointer        :: Station => NULL()
     !> supporting points on connecting line, distance to GNSS receiver in m
     real (wp), dimension(:), pointer   :: SlantSteps => Null()
     !> cartesian coordinates of the supporting points
     real (wp), dimension(:,:), pointer :: SlantPointsCart => Null()
     !> ellipsoidal coordinates of the supporting points
     real (wp), dimension(:,:), pointer :: SlantPointsEll => Null()
     !> geoid coordinates of the supporting points
     real (wp), dimension(:,:), pointer :: SlantPointsGeo => Null()
     !> index of the COSMO grid cell containing the supporting point
     integer, dimension(:,:), pointer   :: SlantPointsIdx => Null()
     !> true bended signal path
     real (wp), dimension(:,:), pointer :: SlantPath => Null()
     !> refractivity on the bended slant path
     real (wp), dimension(:), pointer   :: SlantRefrac => Null()
     !> refractivity at the supporting points on connecting line
     real (wp), dimension(:), pointer   :: LineRefrac => Null()

     real (wp)                          :: GeoidCorr  !< Geoid undulation
     !> refractivity and height at the last two grid layers
     real(wp), dimension(1:2,1:2)       :: ExtraPolNode
     !> total number of supporting points along the slant path
     integer                            :: Ntot
     !> number of supporting points inside the model grid
     integer                            :: Nmod
     !> number of supporting points above the model top layer
     integer                            :: Nup
     !> ??? number of points in horizontal domain
     integer                            :: Nhor
     !> ??? number of points allocated / stored
     integer                            :: Naloc
     !> number of neighbours in model grid
     integer                            :: Nnghb
     real (wp)                          :: htop
     real (wp)                          :: dx
     !> indirect index array  (Nnghb,Naloc)
     integer,   dimension(:,:), pointer :: idxi => Null()
     !> interpolation weights (Nnghb,Naloc)
     real (wp), dimension(:,:), pointer :: wi   => Null()
     !> ??? not used ???
!     logical                            :: assimilate
     !> slant total dely computed with raytracer along curved path
     real(wp)                           :: STD
     !> slant total dely computed along straight line
     real(wp)                           :: STDline
  end type SlantData


  !---------------------------------------------------------------------
  ! structure p_column
  !---------------------------------------------------------------------
  !> \brief Derived type for column data required by the STD operator
  !>
  !> The type p_column provides the model columns surrounding the
  !> supporting points on the slant path.  \n
  !> The second part of the STD operator needs access to the atmospheric
  !> along the signal path. The derived type p_column provides one column
  !> of the required model fields, the coordinates of the grid nodes and
  !> geoid corrections.
  !>
  !> <table>
  !> <tr><th>variable  <th>description
  !> <tr><td><b> dlat </b>  <td>
  !>    latitude in degrees of all grid nodes in vertical column \n
  !>    COSMO: dlat = rlat * raddeg, COSMO field rlat is in radian \n
  !>    3D-Var: dlat = grid%rlat * raddeg, COSMO field rlat is in radian
  !> <tr><td><b> dlon </b>  <td>
  !>    longitude in degrees of all grid nodes in vertical column \n
  !>    COSMO: dlon = rlon * raddeg, COSMO field rlon is in radian \n
  !>    3D-Var: dlon = grid%rlon * raddeg, COSMO field rlon is in radian
  !> <tr><td><b> gpm </b>  <td>
  !>    geometric height = height above geoid = height above mean sea level
  !>    in meters at full model levels \n
  !>    gpm = 0.5 * ( hhl(i,j,k) + hhl(i,j,k+1) ) \n
  !>    (HHL - height at half levels => gpm - height at full levels)
  !> <tr><td><b> geoid </b>  <td>
  !>    geoid correction in meters at the position (dlon, dlat). All heights
  !>    in the array gpm need to be shifted by the geoid correction in order
  !>    to obtain heights above ellipsoid.
  !> <tr><td><b> p </b>  <td>
  !>    pressure in Pa at full model levels
  !> <tr><td><b> q </b>  <td>
  !>    specific humidity at full model levels
  !> <tr><td><b> t </b>  <td>
  !>    temperature in K at full model levels
  !> </table>
  !>
  !> Variable col of type p_column
  !>
  !> One single variable of type p_column provides one model column. For
  !> interpolation on the supporting points several neighbored columns are
  !> required (4 for COSMO, 3 for ICON). These data are required for all
  !> supporting points along the signal path. Therefore, a 2-dimensional array
  !> col(i,j) is used which provides all data required by the STD operator:
  !>
  !> @verbatim
  !> col  -  vertical grid columns
  !>         col(i,j), i=1, Nnghb, numberof neighbored columns
  !>                               Nnghb = 3 - GME
  !>                               Nnghb = 4 - COSMO
  !>                   j=1, slant%Nmod, number of columns along slant path
  !>                        Nmod - number of supporting points inside the model
  !>                        For each supporting point a pointer to the
  !>                        corresponding column is provided.
  !> col(i,j)%dlat - latitude, degrees, -90 <= dlat <= +90 degrees
  !> col(i,j)%dlon - longitude, degrees, -180 <= dlat <= +180 degrees
  !>                 Warning: The transformations require/provide the longitude
  !>                          in radian and 0 <= lon <= 2*Pi
  !> @endverbatim
  !---------------------------------------------------------------------
  type p_column
    real(wp) ,pointer :: p   (:) !< pressure          at full model levels (Pa)
    real(wp) ,pointer :: q   (:) !< specific humidity at full model levels
    real(wp) ,pointer :: t   (:) !< temperature       at full model levels  (K)
    real(wp) ,pointer :: gpm (:) !< geometric height  at full model levels  (m)
    real(wp)          :: geoid   !< geoid correction  at      model column  (m)
    real(wp)          :: dlat    !< latitude                           (degree)
    real(wp)          :: dlon    !< longitude                          (degree)
  end type p_column


  type (SlantData), dimension(:), pointer        :: STD => Null()
  integer                                        :: NSTDobs
  integer :: NSTDobsLen ! Byte length of one variable of type "SlantTotalDelay"
  ! NSTDobs number of STD observations in array STD:
  type (SlantTotalDelay), dimension(:), allocatable :: STDobs
! character (len=256)                            :: STDFile
  type (SlantHeader), save                       :: Header
  integer                                        :: NSat  !, NStat
  type (GPSStation), dimension(:), pointer       :: StatList => NULL()
  integer(i2),       dimension(:), pointer       :: StatHash => NULL()
  integer, dimension(:), pointer                 :: SatList => NULL()
  integer, parameter                             :: Nintervals = 10

  integer :: err, FileNr = 0


  !--- namelist parameter -----------------------------------------------------
  !>   namelist /STD_OBS/ - Parameter used by the STD operator
  !>
  !> @verbatim
  !> NStepVertMod - Number of supporting points on the signal path inside
  !>                the model. "NStepVertMod" is the number of points in
  !>                zenith direction. In case of STDs this number is 'mapped'
  !>                on the signal path using 1/sin(elevation).
  !> NStepVertTop - Number of supporting points on the signal path above the
  !>                model top. "NStepVertTop" is the number of points in
  !>                zenith direction. In case of STDs this number is 'mapped'
  !>                on the signal path using 1/sin(elevation).
  !>                Usually "NStepVertTop" can be much smaller than
  !>                "NStepVertMod" as the contribution to the ZTD/STD is rather
  !>                small and the information about the refractivity above the
  !>                model top is limited.
  !> Hmax         - Maximum height in m. It is assumed that the neutral
  !>                atmosphere reaches up to Hmax and the STD is integrated
  !>                from the station height to Hmax.
  !> HScaleP      - Scale hight in m of a hypothetical pressure profile. This
  !>                profile is used to scale the density of supporting points
  !>                with height. Smaller numbers lead to more points in the
  !>                lower atmosphere. HScaleP is used for the signal path
  !>                inside the model.
  !> HScaleP2     - Like "HScaleP" but used for the region above the model top.
  !>                "HScaleP2" should be much larger than "HScaleP" in order
  !>                to obtain a homogeneous distribution of supporting points
  !>                along the whole signal path from the receiver up to "Hmax".
  !> UseRaytracer - Maximum elevation in degrees for using the raytracer.
  !>                The raytracer is called only if the elevation of the STD
  !>                is below this elevation. For larger elevations it is
  !>                assumed that the signal path between the satellite and the
  !>                receiver is a straight line.
  !>                UseRaytracer = 90°  - The raytracer is called for all
  !>                                      ZTDs and STDs.
  !>                UseRaytracer = 0°   - The raytracer is never called and
  !>                                      a straight line is assumed for all
  !>                                      elevations.
  !> k1, k2, k3   - Refractivity coefficients used to compute the refractivity
  !>                with the Smith & Weintraub formula. If these coefficients
  !>                are not specified they are taken from Bevis, 1994.
  !>                Units: k1 - K / hPa
  !>                       k2 - K / hPa
  !>                       k3 - K^2 / hPa
  !>
  !> Parameter used by the assimilation system
  !>
  !> ZTDerror     - Assumed ZTD observation error in m.
  !>                For STDs "ZTDerror" is mapped
  !>                to the given elevation using the Global Mapping Function
  !>                (GME). The corresponding observation error is written to
  !>                the feedback files.
  !> Href         - Reference height in m above the GNSS station used by the
  !>                LETKF. In the feedback files the parameters plevel, dlat,
  !>                dlon give the position on the signal path where the
  !>                height above the station is Href. For Href = 0 the station
  !>                position and the pressure at the station are given.
  !> MaxSTDobs    - Maximum number of STDs to read/process
  !> MaxZTDobs    - Maximum number of ZTDs to read/process
  !>                The memory required to read and process the GNSS
  !>                observations grows with the number of observations and may
  !>                exceed the available memory. As the number of observations
  !>                is unpredictable an upper limit is set in order to
  !>                guarantee stable operation. The limit is applied to the read
  !>                buffer which starts allways at the beginning of the full
  !>                hour. The first MaxS/ZTDobs observations of the hour are
  !>                read and the number of selected observations might be
  !>                smaller or zero.
  !>
  !> Parameter governing the output messages
  !>
  !> verbose      - Verbosity level of program messages
  !>                At higher verbosity levels the output of all lower
  !>                levels will also be provided.
  !>                verbose = 0 (default) - output limited to some status
  !>                                        messages
  !>                verbose = 1 - messages necessary to follow the program flow
  !>                verbose = 2 - more messages and output of some important
  !>                              variables
  !>                verbose = 3 - debug mode, huge amount of output
  !>
  !> Parameter for bias correction (3D-VAR, LETKF)
  !>
  !> biascor_mode - ??
  !> t_decay      - ??
  !> n_required   - ??
  !> bc_fallback  - ??
  !>
  !> Paramter which define the location of external files
  !>
  !> GeoidFile    - Path and file name of the geoid file
  !>                The geoid file provides the geoid undulation required to
  !>                transform heights above geoid into heights above ellipsoid
  !>                and vice versa. So far only the "eigen-6c3stat" geoid is
  !>                supported.
  !>                This option is used by COSMO only, ICON EnVar uses the
  !>                EGM96 geoid. The STD operator in COSMO does not work
  !>                without a valid geoid and COSMO will be terminated if the
  !>                STD operator is started without a GeoidFile.
  !>
  !> Parameter for reading ASCII STD files
  !>           (BUFR data will be read automatically if available)
  !>
  !> read_ascii   - ? not used ?
  !> RefTimeStr   - Analysis reference time, string yyyymmddhh
  !>                The ASCII STD files containing data at this time plus/
  !>                minus STDperiod will be read.
  !> STDperiod    - Period in minutes around the  reference time
  !> SlantPath    - Path to the ASCII STD files
  !>                At the end of the path a directory indicating the year
  !>                will be appended automatically (e.g. .../y2015).
  !> std_obs_file - List of STD files which should be read. ? not active ?
  !> HORIFile     - Path and file name of GNSS station list.
  !> DomesFile    - Path and file name of Domes file.  ? not active ?
  !>
  !> @endverbatim
  !--- namelist parameter -----------------------------------------------------
! integer ::   NAMELIST_STD_OBS

  integer ,parameter :: mo           = 12     ! max number of observation
                                              ! input files
  integer ,parameter :: mf           = 256    ! max length of file names
  !-----------------------------
  ! To be read from the namelist
  !-----------------------------
  logical            :: read_ascii = .false.  ! read ascii file
  integer            :: verbose    = 0        ! Verbosity level
  character (len=10) :: RefTimeStr = ''       ! Reference time (analysis)
                                              ! yyyymmddhh
  integer            :: STDperiod    =  0     ! STD period, minutes
  integer            :: NStepVertMod = 60     ! Number of vertical points
                                              ! inside the model
  integer            :: NStepVertTop = 20     ! Number of vertical points
                                              ! above the model
  integer            :: NStepOpt = 2          ! option for scaling points
  real (wp)          :: Hmax    = 150000.     ! Maximum height for STD
                                              ! integration in m
  ! old settings
  !real               :: HScaleP = 6500.       ! air pressure scale height
  !real               :: HScaleP2 = 25000.     ! air pressure scale height

  real               :: HScaleP = -1.0       ! air pressure scale height
  real               :: HScaleP2 = -1.0     ! air pressure scale height

  ! Defaults are set depending on other parameters and the model
  integer, dimension(Nintervals) :: Hlevel  = -1.0  ! model levels per interval
  real, dimension(Nintervals)    :: Heights = -1.0  ! model level heights
  integer, dimension(Nintervals) :: Hpoints = -1.0  ! supporting points per int.

  integer            :: UseRaytracer = 50     ! use raytrace if elev. < UseRay.
  real               :: ZTDerror = 0.012      ! ZTD observation error, m
  real               :: Href = 0.0_wp         ! reference height in m, LETKF
  ! namelist parameters for checking/rejecting observations
  real               :: ZTDminUse = 0.0       ! minimum ZTD used
  real               :: ZTDmaxUse = 3.0       ! maximum ZTD used
  real               :: StatBelowSurface = 500.0  ! reject stat. below mod. surf
  real               :: StatAboveSurface = 1000.0 ! reject stat. above mod. surf
  real               :: StatBelowColumn  = 500.0  ! reject stat. below max. col
  integer            :: MaxSTDobs = 100000    ! max. number of STDs
  integer            :: MaxZTDobs =  50000    ! max. number of ZTDs
! real (wp)          :: k1 = 77.60_wp         ! refractivity coefficients
! real (wp)          :: k2 = 70.40_wp         ! refractivity coefficients
! real (wp)          :: k3 = 3.739E5_wp       ! refractivity coefficients
  character(len=mf)  :: GeoidFile = ''        ! Path + file name of geoid
  character(len=mf)  :: SlantPath = ''        ! Path to slant data
  character(len=mf)  :: std_obs_file (mo) ='' ! STD observation input files
  character(len=mf)  :: HORIFile = 'stako_gps_SNX_G_15_HORI'
                                              ! GNSS station file
  character(len=mf)  :: DomesFile = 'smark_domes_SNX_RTT'
                                              ! GNSS station locations
  integer            :: ztd_col   = 0         !   0: use all neighbours
                                              !   1: use nearest neighbour
                                              !   2: use lowest model surface
                                              !   3: interpolate
                                              ! < 0: duplicate columns
  integer            :: pl_method = 0         ! method to estimate vertical location
#ifdef __COSMO__
  ! The bias correction is not used in COSMO but the corresponding
  ! variables should exist (see USE statements)
  integer            :: biascor_mode = 0
  real(wp)           :: t_decay
  integer            :: n_required
  logical            :: bc_fallback
  integer            :: flag_biasc_gpsgb
  integer ,parameter :: BC_INVAL = -9  ! not set so far
  integer ,parameter :: BC_NOBC  =  0  ! no bias correction
  integer ,parameter :: BC_UP    =  1  ! only update bias corrrection file
  integer ,parameter :: BC_FG    =  2  ! apply  bias corr. from first guess
  integer ,parameter :: BC_AN    =  3  ! apply  bias corr. from analysis
  integer ,parameter :: BC_VARBC =  4  ! variational bias correction
#endif

  namelist /STD_OBS/ RefTimeStr, STDperiod, read_ascii, verbose,    &
                     NStepVertMod, NStepVertTop, Hmax, HScaleP,     &
                     HScaleP2, Hlevel, Heights, Hpoints,            &
                     k1, k2, k3, UseRaytracer, ZTDerror,            &
                     Href, MaxSTDobs, MaxZTDobs, GeoidFile,         &
                     ZTDminUse, ZTDmaxUse, StatBelowSurface,        &
                     StatAboveSurface, StatBelowColumn,             &
                     SlantPath, HORIFile, DomesFile, std_obs_file,  &
                     biascor_mode, t_decay, n_required, bc_fallback,&
                     ztd_col, pl_method

  !-------------------------------------------------
  ! Information derived from the namelist parameters
  !-------------------------------------------------
  character(len=mf)  :: STDfile               ! STD file with path
! character(len=mf)  :: StationFile           ! STD station file with path
! real (wp)          :: RefTimeMJD            ! Reference time (analysis)
                                              ! Modified Julian Date (MJD)
  ! Only slant data from the period [StartMJD, EndMJD] will be analysed:
  real (wp)          :: StartMJD              ! = RefTimeMJD - STDperiod
  real (wp)          :: EndMJD                ! = RefTimeMJD
  ! working copy of Hlevel, Heights, Hpoints
  integer, dimension(:), allocatable :: Plevel   ! model levels per interval
  real, dimension(:), allocatable    :: Pheights ! model level heights
  integer, dimension(:), allocatable :: Ppoints  ! supporting points per inter.

!==============================================================================
  contains
!==============================================================================

!  subroutine std_handler
!    !------------------------------
!    ! temporary routine for testing
!    !------------------------------
!    call read_std_nml
!    call read_gnss_stations
!    call read_std_obs(STDFile)
!!    call std_path_coord
!!   call std_redistribute data
!!   call std_call_operator
!  end subroutine std_handler
!
!==============================================================================
  subroutine read_std_nml (unit, initonly)
  !------------------------
  ! read namelist /STD_OBS/
  !------------------------
  integer, intent(in)           :: unit     ! Fortran unit number for namelist
  logical, intent(in), optional :: initonly ! initialize, do not read namelist

    integer            :: ierr, i, c
    integer            :: yyyy, mm, dd, hh, mi,ss
    integer            :: Year, Month, Day, Hour, Minit=0, DOY
    real (wp)          :: JD, Sec=0.0D0
    character (len=4)  :: YearStr
    character (len=3)  :: DoyStr
    character (len=2)  :: HourStr
    logical            :: FileExist
    logical ,save      :: first = .true.
    logical            :: debug = .false.
    character (len=2)  :: char2
    character (len=20) :: fmt1, fmt2
    logical            :: linitonly

    !------------------------
    ! read namelist only once
    !------------------------
    if (.not. first) return
    first = .false.

    linitonly = .false.; if (present (initonly)) linitonly = initonly
    !--------------------------------------------
    ! set defaults depending on 'ga3_biasc_airep'
    !--------------------------------------------
    select case (flag_biasc_gpsgb)
    case (-1)
      biascor_mode = BC_NOBC
      bc_fallback  = .false.
    case ( 0)
      biascor_mode = - BC_FG
      bc_fallback  = .true.
    case ( 1)
      biascor_mode = - BC_FG
      bc_fallback  = .false.
    end select

    if (debug) write(*,*) 'read_std_nml> start ..'
    if (.not. linitonly) then
    !--------------
    ! read namelist
    !--------------
#if defined(__ibm__)
    !-------------------
    ! catch error on IBM
    !-------------------
    read (unit ,nml=STD_OBS ,iostat=ierr)
    if (ierr/=0) call model_abort (-1,-1,                                &
                                   'GPSGB: ERROR in namelist /STD_OBS/', &
                                   'read_std_nml'                 )
#else
    read (unit ,nml=STD_OBS)
#endif
    end if

    !-----------------------------------------------------
    ! adjust 'biascor_mode' depending on 'ga3_biasc_airep'
    !-----------------------------------------------------
    if (biascor_mode < 0) then
      select case (flag_biasc_gpsgb)
      case (-1)
        biascor_mode = BC_NOBC
      case (0:1)
        biascor_mode = - biascor_mode
      end select
    endif

    !------------------------------------------------------------------------
    ! Scale heights and vertical ranges ():
    ! Set to default if no namelist entries are given
    !------------------------------------------------------------------------
    if (all(Hlevel < 0) .and. all(Heights < 0) .and. all(Hpoints < 0)  &
         .and. HScaleP > 0 .and. HScaleP2 > 0 ) then
       ! Use old profile only if HScaleP and HScaleP2 are given in the
       ! namelist and Hlevel, Hscale, Hpoints are not given
       ! The old defaults HScaleP  =  6500.
       !                  HScaleP2 = 25000.
       ! are no longer available!

       if (verbose >= 3) then
          write(*,*) 'read_std_nml> Use old options:'
          write(*,*) 'read_std_nml> HScaleP  = ', HScaleP
          write(*,*) 'read_std_nml> HScaleP2 = ', HScaleP2
          write(*,*) 'read_std_nml> NStepOpt = ', NStepOpt
       end if

    else
       ! use new settings in Hlevel, Hscale, Hpoints

       NStepOpt = 3

       if (all(Hlevel < 0) .and. all(Heights < 0) .and. all(Hpoints < 0)) then
          ! Use new default
          ! ICON defaults (read_std_nml is not called from COSMO)
          Hlevel(1:3)  = (/34, 16, 40/)
          Heights(1:3) = (/8821.524_sp, 15354.572_sp, 73681.773_sp/)
          Hpoints(1:3) = (/36, 18, 42/)
       end if

       ! Don't use NStepVertMod from namelist: Reset to sum of Hpoints
       c = count(Hpoints > 0)
       NStepVertMod = sum(Hpoints(1:c))

    end if

    !---------
    ! printout
    !---------
    write(char2,'(i2)') Nintervals       ! number of elements in array
    fmt1 = '(a,' // char2 // '(f10.3, tr2))'
    fmt2 = '(a,' // char2 // '(i3, tr2))'
    !write (6,*) fmt1
    !write (6,*) fmt2

    write (6,'(a)') repeat('-',79)
    write (6,'()')
    if (linitonly) then
      write (6,'(a)')    '  using defaults for /STD_OBS/:'
    else
      write (6,'(a)')    '  namelist /STD_OBS/:'
    end if
    write (6,'()')
    write (6,'(a,a)')    '    RefTimeStr      = ',trim (RefTimeStr)
    write (6,'(a,i6)')   '    STDperiod [min] = ',STDperiod
    write (6,'(a,i6)')   '    NStepOpt        = ',NStepOpt
    write (6,'(a,i6)')   '    NStepVertMod    = ',NStepVertMod
    write (6,'(a,i6)')   '    NStepVertTop    = ',NStepVertTop
    write (6,'(a,f7.0)') '    Hmax            = ',Hmax
    write (6,'(a,f7.0)') '    HScaleP         = ',HScaleP
    write (6,'(a,f7.0)') '    HScaleP2        = ',HScaleP2
    write (6,fmt2)       '    Hlevel          = ',pack (Hlevel,  Hlevel  >= 0)
    write (6,fmt1)       '    Heights         = ',pack (Heights, Heights >= 0._sp)
    write (6,fmt2)       '    Hpoints         = ',pack (Hpoints, Hpoints >= 0)
    write (6,'(a,f9.2)') '    Href  [m]       = ',Href
    write (6,'(a,f9.2)') '    ZTDminUse       = ',ZTDminUse
    write (6,'(a,f9.2)') '    ZTDmaxUse       = ',ZTDmaxUse
    write (6,'(a,f9.2)') '    StatBelowSurface= ',StatBelowSurface
    write (6,'(a,f9.2)') '    StatAboveSurface= ',StatAboveSurface
    write (6,'(a,f9.2)') '    StatBelowColumn = ',StatBelowColumn
    write (6,'(a,f10.3)')'    k1 [K/hPa]      = ',k1
    write (6,'(a,f10.3)')'    k2 [K/hPa]      = ',k2
    write (6,'(a,f8.1)') '    k3 [K^2/hPa]    = ',k3
    write (6,'(a,i6)')   '    UseRaytracer    = ',UseRaytracer
    write (6,'(a,i6)')   '    verbose         = ',verbose
    write (6,'(a,f6.4)') '    ZTDerror [m]    = ',ZTDerror
    write (6,'(a,a)')    '    SlantPath       = ',trim (SlantPath)
    write (6,'(a,a)')    '    HORIFile        = ',trim (HORIFile)
    write (6,'(a,a)')    '    DomesFile       = ',trim (DomesFile)
    write (6,'(a,a)')    '    std_obs_file    = ',trim (std_obs_file(1))
    write (6,'(a,i6)'  ) '    ztd_col         = ',ztd_col
    write (6,'(a,i6)'  ) '    pl_method       = ',pl_method
    write (6,'()')
    write (6,'(a,i6)'  ) '    biascor_mode    = ', biascor_mode
    write (6,'(a,l6)'  ) '    bc_fallback     = ', bc_fallback
    write (6,'(a,f6.1)') '    t_decay         = ', t_decay
    write (6,'(a,i6)'  ) '    n_required      = ', n_required
    write (6,'()')
    do i=2,mo
      if (std_obs_file(i)/='') &
      write (6,'(a,a)')  '                      ',trim (std_obs_file(i))
    end do
    write (6,'()')

    if (verbose > 1) debug = .true.

    if (Href < 0.0) then
       call model_abort (-1,-1, 'GPSGB: Href < 0 invalid', &
                                'read_std_nml'     )
    end if
    select case (pl_method)
    case (0, 20)
       ! Valid options, nothing to do
       continue
    case default
       call model_abort (-1,-1, &
                         'GPSGB: Invalid pl_method, valid are 0 and 20', &
                         'read_std_nml')
    end select

    STDFile = trim(std_obs_file(1))

#ifndef __COSMO__
    ! Use "yyyymmddhh_ana" from the RUN namelist if "RefTimeStr"
    ! is not given. The corresponding GPS slant data will be selected
    ! automatically.
    if (RefTimeStr == '') then
       call i_time (ana_time, yyyy, mm, dd, hh, mi, ss)
       if (yyyy .gt. 1900) then
          write(RefTimeStr(1:4),'(i4.4)') yyyy
          write(RefTimeStr(5:6),'(i2.2)') mm
          write(RefTimeStr(7:8),'(i2.2)') dd
          write(RefTimeStr(9:10),'(i2.2)') hh
       end if
    end if

    ! Reference Time is assumed to be
    ! YYYYMMDD or YYYYMMDDHH
    YearStr = ''
    DoyStr = ''
    HourStr = ''
    STDFile = ''
    if (RefTimeStr /= '' .and. SlantPath /= '') then
       read(RefTimeStr(1:4),'(i4.4)')  Year
       YearStr = RefTimeStr(1:4)
       read(RefTimeStr(5:6),'(i2.2)')  Month
       read(RefTimeStr(7:8),'(i2.2)')  Day
       if (len_trim(RefTimeStr) .gt. 8) then
          read(RefTimeStr(9:10),'(i2.2)')  Hour
          HourStr = RefTimeStr(9:10)
       else
          Hour = 0
       end if
       Minit = 0
       Sec   = 0.0
       ! (Analysis) reference time
       JD = JulianDate(Year, Month, Day, Hour, Minit, Sec)
       ! JD of first observations
       JD = JD - (0.5_wp*real(STDperiod,dp) / 1440.0_wp)
       StartMJD = JD2MJD(JD)
       ! MJD of last observations
       EndMJD   = StartMJD + (real(STDperiod,dp) / 1440.0_wp)
       i = 0
       do
          i = i + 1
          call GregorianDate( JD, Year, Month, Day, Hour, Minit, Sec)
          if (debug) write(*,*) 'JD = ', JD,  MJD2JD(EndMJD)
          if (debug) write(*,*) Year, Month, Day, Hour, Minit, Sec
          call GREG2DOY (int(Year), int(Month), int(Day), Doy)
          write(HourStr,'(i2.2)') Hour
          write(DoyStr,'(i3.3)') Doy
          write(YearStr,'(i4.4)') Year
          STDFile = 'slv_slant_'//YearStr//'_'//DoyStr//'_'//HourStr
          std_obs_file(i) = trim(SlantPath)//'/y'//YearStr//'/'//trim(STDFile)
          write(*,*) trim(std_obs_file(i))
          JD = JD + (60.0_wp / 1440.0_wp)
          if (JD > MJD2JD(EndMJD)) exit
       end do
    end if
    if (std_obs_file(1) /= "") then
       if (STDFile /= "" .and.                        &
            index(std_obs_file(1),trim(STDFile)) == 0 ) then
          write(*,*)
          write(*,*) 'read_std_nml> Warning: ',                       &
                     'Reference date and observation date do not match'
          write(*,*) 'STDFile         : ', trim(STDFile)
          write(*,*) 'std_obs_file(1) : ', trim(std_obs_file(1))
       end if
       STDFile = trim(std_obs_file(1))
       ! Auch hier StartMJD und EndMJD ausrechnen !!!
    else if (SlantPath /= "") then
       ! slant file name with path
       STDFile = trim(SlantPath)//'/y'//YearStr//'/'//trim(STDFile)
    end if
    if (STDFile /= "") then
       inquire(file=trim(STDFile), exist=FileExist)
       if (FileExist) then
          write (6,'(/a,a)')    '    Slant file : ', trim(STDFile)
       else
          write (6,'(/a,a,a)')  '    Slant file : ', trim(STDFile), &
                                ' does not exist'
          err = 1
       end if
    end if
    if (debug) write(*,*) 'read_std_nml> ... end'
    !std_obs_file(1) = trim(STDFile)
#endif


!!$
!!$    ! TEST TEST TEST TEST TEST TEST
!!$    std_obs_file(1) = trim(STDFile)
!!$    ! TEST TEST TEST TEST TEST TEST
!!$
!!$
!!$    ! Check for HORIFile in local directory and in SlantPath
!!$    inquire(file=trim(HORIFile), exist=FileExist)
!!$    if (FileExist) then
!!$       write (6,'(a,a)')    '    HORI  file : ', trim(HORIFile)
!!$    else
!!$       HORIFile = trim(SlantPath)//'/'//trim(HORIFile)
!!$       inquire(file=trim(HORIFile), exist=FileExist)
!!$       if (FileExist) then
!!$          write (6,'(a,a)')    '    HORI  file : ', trim(HORIFile)
!!$       else
!!$          write (6,'(a,a,a)')  '    HORI  file : ', trim(HORIFile), &
!!$                               ' does not exist'
!!$          err = 1
!!$       end if
!!$    end if
!!$
!!$    ! Check for DomesFile in local directory and in SlantPath
!!$    inquire(file=trim(DomesFile), exist=FileExist)
!!$    if (FileExist) then
!!$       write (6,'(a,a)')    '    Domes file : ', trim(DomesFile)
!!$    else
!!$       DomesFile = trim(SlantPath)//'/'//trim(DomesFile)
!!$       inquire(file=trim(DomesFile), exist=FileExist)
!!$       if (FileExist) then
!!$          write (6,'(a,a)')    '    Domes file : ', trim(DomesFile)
!!$       else
!!$          write (6,'(a,a,a)')  '    Domes file : ', trim(DomesFile), &
!!$                               ' does not exist'
!!$          err = 1
!!$       end if
!!$    end if

    !-------------------------------
    ! Update some internal constants
    !-------------------------------
    call UpdateK123 (k1, k2, k3, err)
    if (err .ne. 0) call model_abort (-1,-1,                           &
                                   'GPSGB: ERROR changing k1, k2, k3', &
                                   'read_std_nml'               )

  end subroutine read_std_nml

  !---------------------------------------------------------------------
  ! subroutine read_std_nml_cosmo
  !---------------------------------------------------------------------
  !
  !> @brief Reads the namelist of the STD operator from COSMO
  !>
  !> <b> call read_std_nml_cosmo (nuin, iz_err) </b>
  !>
  !> The namelist STD_OBS of the STD operator is read from an already open file.
  !> This routine does neither open nor close the file containing the
  !> STD namelist. \n
  !> The namelist variables are available in the mo_std_operator module.
  !>
  !> @param[in]  nuin   file unit, open file which contains the STD namelist
  !> @param[out] iz_err return code from the Fortran read statement
  !---------------------------------------------------------------------
  subroutine read_std_nml_cosmo (nuin, iz_err)

    ! Subroutine arguments:
    ! --------------------
    integer, intent(in)  :: nuin
    integer, intent(out) :: iz_err

    integer :: c
    !---------------------------------------------------------------------

    !----------------------
    ! Read namelist STD_OBS
    !----------------------
    READ (nuin, STD_OBS, IOSTAT=iz_err)

    !------------------------------------------------------------------------
    ! Scale heights and vertical ranges ():
    ! Set to default if no namelist entries are given
    !------------------------------------------------------------------------
    if (all(Hlevel < 0) .and. all(Heights < 0) .and. all(Hpoints < 0)  &
         .and. HScaleP > 0 .and. HScaleP2 > 0 ) then
       ! Use old profile only if HScaleP and HScaleP2 are given in the
       ! namelist and Hlevel, Heights, Hpoints are not given
       ! The old defaults HScaleP  =  6500.
       !                  HScaleP2 = 25000.
       ! are no longer available!

       if (verbose >= 3) then
          write(*,*) 'read_std_nml> Use old options:'
          write(*,*) 'read_std_nml> HScaleP  = ', HScaleP
          write(*,*) 'read_std_nml> HScaleP2 = ', HScaleP2
          write(*,*) 'read_std_nml> NStepOpt = ', NStepOpt
       end if

    else
       ! use new settings in Hlevel, Heights, Hpoints

       NStepOpt = 3

       if (all(Hlevel < 0) .and. all(Heights < 0) .and. all(Hpoints < 0)) then
          ! Use new default
          ! COSMO-DE defaults
          ! (2.8 km horz. grid spacing, 50 main levels up to 21500 m)
          Hlevel(1:4)  = (/20, 10, 10, 10/)
          Heights(1:4) = (/2829.47_sp, 6773.21_sp, 12877.68_sp, 21500.00_sp/)
          Hpoints(1:4) = (/21, 11, 11, 11/)
       end if

       ! NStepVertMod is not yet set: Set to sum of Hpoints in namelist
       c = count(Hpoints > 0)
       NStepVertMod = sum(Hpoints(1:c))


       ! Call ProfileSetup to print reference profile, only on PE=0,
       ! will later be called on all PEs to really setup the profile
       call ProfileSetup(sum(Hlevel(1:c)), .true.)

    end if

  end subroutine read_std_nml_cosmo

!==============================================================================

  subroutine read_gnss_stations(MJDstart, MJDend)
  !------------------------------------------------------
  ! Read gnss station information: name, coordinates, ...
  ! currently only read on PE 0
  ! keep results in derived type CoordList
  !------------------------------------------------------

    ! Subroutine arguments:
    ! --------------------
    real (wp), optional :: MJDstart, MJDend

    ! Local parameters:
    ! ----------------
! real (wp)                                :: StartDate, EndDate
  integer                                  :: NStat, i, n
  type (StrucDomes), dimension(:), pointer :: GPSdomes => NULL()
  integer, dimension(:), allocatable       :: DomI
  character (len=4)                        :: center
! integer                                  :: StatAlloc

  if (verbose >= 1) write(*,*) 'read_gnss_stations> Start ...'

  ! Use reference period given in subroutine arguments, if available
  if( present(MJDstart) ) StartMJD = MJDstart
  if( present(MJDend)   ) EndMJD   = MJDend

  if (verbose >= 1) write(*,*) 'read_gnss_stations> start ...'
  if (verbose >= 2) write(*,*) 'Read GNSS station coordinates :',      &
                     trim(HORIFile), '  MJD=', StartMJD, ' - ', EndMJD

  call ImportStakoSNXarr(HORIFile, CoordList, NStat,     &
                         MJD=StartMJD, EndMJD=EndMJD)

!!$  write(*,*) 'read_gnss_stations> ImportStakoSNXarr: NStat = ', &
!!$         NStat, ubound(CoordList,1)
!!$  do i=1, NStat
!!$     write(*,*) i, CoordList(i)
!!$  end do

  if (ubound(CoordList,1) .le. 0) then
     write(*,*) 'read_gnss_stations> Error - no stations read from ', &
                'file >>', trim(HORIFile), '<<'
     stop
  end if

!!$  call  ImportDomes(DomesFile, GPSdomes, StartMJD, EndMJD)
!!$  if (ubound(GPSDomes,1) .le. 0) then
!!$     write(*,*) 'read_gnss_stations> Error - no stations found in ', &
!!$                 trim(DomesFile)
!!$  end if

  !---------------------------------------------------------------------
  ! Cartesian station coordinates are not available in the
  ! HORI file: Start transformation ellipsoidal to cartesian coordinates
  ! GNSS stations need to be referenced by the station ID. Create simple
  ! Hash table with the station ID as array index.
  !---------------------------------------------------------------------
  center = 'GFS1'
  n = index(SlantPath, 'EPOS6')
  if (n > 0) center = 'GFS0'   ! GFZ slants processed with EPOS 6
  n = index(SlantPath, 'EPOS8')
  if (n > 0) center = 'GFS1'   ! GFZ slants processed with EPOS 8

  allocate( GNSShash(1:9999) )
  GNSShash = -99
  do i=1, NStat
     call  Ellips2Cart( CoordList(i)%CoordEll(1),    & ! in: longitude
                        CoordList(i)%CoordEll(2),    & ! in: latitude
                        CoordList(i)%CoordEll(3),    & ! in: altitude
                        CoordList(i)%CoordCart(1),   & ! out: X
                        CoordList(i)%CoordCart(2),   & ! out: Y
                        CoordList(i)%CoordCart(3),   & ! out: Z
                        WGS84Param)                 ! in: reference ellipsoid
     GNSShash(CoordList(i)%ID) = i
     !----------------------------------------------------------------------
     ! Create 8 character station name consisting of 4 char. GNSS short name
     ! and 4 char. processing center and product identifier following the
     ! E-GVAP convention for ZTD data (BUFR).
     ! GFS0: processing center is GFZ (=> GF)
     !       slant product (GPS only, EPOS 6) (=> S0)
     !----------------------------------------------------------------------
     CoordList(i)%Name = CoordList(i)%SName // center

     ! GeoidCorr and ModelSurf will be set later, initialize with -999
     CoordList(i)%GeoidCorr = -999.999_wp
     CoordList(i)%ModelSurf = -999.999_wp
     !CoordList(i)%Name      = ' '
     !
     ! Set some WMO identifiers which will later appear in the BUFR tables,
     ! will be copied to feedback files
     CoordList(i)%DataCategory    = 0   !=> surface data
     CoordList(i)%DataSubCategory = 14  !=> ground based GNSS data
     CoordList(i)%Center          = 23  !=> STDs provided by GFZ
     CoordList(i)%SubCenter       = 23  !=> STDs provided by GFZ
     !--------------------------------------------------------------------
     ! So far, the full station name is not available, copy the short name
     !--------------------------------------------------------------------
!    CoordList(i)%Name = CoordList(i)%SName

  end do

  if (verbose >= 2) write(*,*) 'Number of stations available : ', NStat

  if (ubound(CoordList,1) .ne. NStat) then
     !-----------------------------------------------
     ! Read error: Wrong number of GNSS stations read
     !-----------------------------------------------
     write(*,*) 'read_gnss_stations> Error reading ', trim(HORIFile)
     write(*,*) 'read_gnss_stations> Number of array elements :', &
                 ubound(CoordList,1)
     write(*,*) 'read_gnss_stations> NStat = ', NStat
  end if

  if (.false.) then
  !if (ubound(GPSDomes,1) .gt. 0) then
     ! Merge station information, i. e. copy full name and
     ! country to CoordList

     allocate( DomI(1:9999) )
     DomI = 0
     do i=1, ubound(GPSDomes,1)
        DomI(GPSDomes(i)%ID) = i
     end do
     do i=1, NStat
        !write(*,*) DomI(CoordList(i)%ID)
        if (DomI(CoordList(i)%ID) .gt. 0) then
           CoordList(i)%Name = GPSDomes(DomI(CoordList(i)%ID))%Name
           CoordList(i)%Country = GPSDomes(DomI(CoordList(i)%ID))%Country
        else
           CoordList(i)%Name = CoordList(i)%SName
           CoordList(i)%Country = 'Unknown'
        end if
     end do
     deallocate( DomI )

  end if
  if (associated(GPSDomes)) deallocate(GPSDomes)

  if (verbose >= 3) then
     ! print station list
     do i=1, NStat
        write(*,*) i, CoordList(i)%SName,  CoordList(i)%ID !,     &
                  ! CoordList(i)%Name, CoordList(i)%Country
     end do
  end if

  if (verbose >= 1) write(*,*) 'read_gnss_stations> ... end'

  end subroutine read_gnss_stations

!==============================================================================

 subroutine scan_std_obs (STDFile, Header)
 character(len=*)  ,intent(in)  :: STDFile
 type(SlantHeader) ,intent(out) :: Header
 !--------------------------------------------------------
 ! Read the observation data file header only (StatID < 0)
 ! the result is kept in derived type variable Header
 !--------------------------------------------------------
#ifdef __COSMO__
 type(SlantTotalDelay), allocatable  :: STDlist (:)  ! dummy
#else
 type(SlantTotalDelay), pointer  :: STDlist (:) =>NULL() ! dummy
#endif
 integer                         :: ID = -1

 call ImportESTDarr (STDFile, 10, STDlist, Header, StatID=ID)

end subroutine scan_std_obs

!==============================================================================

 subroutine read_std_obs (STDFile, DataStartMJD, DataEndMJD)

 character(len=*)  :: STDFile
 real (wp), optional, intent(out) ::  DataStartMJD, DataEndMJD

 !------------------------------------------------
 ! Read the observations, store the information
 ! currently ASCII data, called on PE 0 only
 ! the result is kept in derived type variable STD
 !-----------------------------------------------

 ! To Do
 ! _____
 ! Define reference date, STDperiod and the number of files to be read
 ! in a consistent way fpr 3D-Var and COSMO:
 ! 3D-Var: Get data from a period centered at the  reference date
 ! COSMO: Read all data in file(s)



 type (SlantTotalDelay), dimension(:),  pointer  :: STDlist =>NULL()
!type (SlantTotalDelay)                          :: OneSTD

#ifdef __COSMO__
 type STDFileData
    type (SlantTotalDelay), dimension(:),  allocatable  :: STDlist
 end type STDFileData
#else
 type STDFileData
    type (SlantTotalDelay), dimension(:),  pointer  :: STDlist =>NULL()
 end type STDFileData
#endif
 type (STDFileData), dimension(1:4) :: STDbuf

 !type (SlantHeader)                       :: Header
 !integer (i4)                             :: NStat, i, NSat, NSat2
 !type (GPSStation), dimension(:), pointer :: StatList=>NULL(), StatList2 =>NULL()
 !integer (i2), dimension(:),pointer   :: StatHash=>NULL(), StatHash2 =>NULL()
 !integer , dimension(:), pointer          :: SatList=>NULL(), SatList2 =>NULL()

 integer :: i, j, n, NSTD, Nlost, MaxBuf        !, n1, n2, Ntest
#ifdef __COSMO__
 integer :: Year, Doy, Hour, Month, Day, Minit
 real (wp) :: JD , Sec                               !, Sec, MJD
#endif
 logical :: test = .false.
 logical :: debug = .false.


 !debug = .true.
  if (verbose >= 1) write(*,*) 'read_std_obs> Start ...'

  if (debug) write(*,*) 'Read slants :', trim(STDfile)

  ! Nur wenige Testdaten prozessieren
  !test = .true.
  !Ntest = 20
  MaxBuf = 4     ! max. number of STD files to read

  ! STD file name convention:
  ! The time given in the file name indicates the beginning of the
  ! observations in that file, i.e. data from the period of time
  ! [t, t+1] hours can be found in that file.
  ! Exsample: slv_slant_2013_105_01 => t = 1:00 UTC
  ! begin_MJD 56397.041886 = 15.4.2013, 1:00 UTC
  ! end_MJD 56397.081816   = 15.4.2013, 1:58 UTC

  !STDFile = '/e/uhome/mbender/AssimCode/TestGNSS/Slants/slv_slant_2009_219_18'
  !call ImportESTDarr(STDFile, STDlist, Header, StatHash, StatList,  &
  !                   NSat, SatList)

!!$  i = index(STDfile,'slv_slant')
!!$  read(STDFile(i+10:i+13),'(i4)') Year
!!$  read(STDFile(i+15:i+17),'(i3)') Doy
!!$  read(STDFile(i+19:i+20),'(i2)') Hour
!!$  !write(*,*) 'Reaf : ', Year, Doy, Hour
!!$
!!$  call DOY2GREG(Doy, Year, Month, Day)
!!$  Minit = 0
!!$  Sec = 0.0D0
!!$  JD = JulianDate(Year, Month, Day, Hour, Minit, Sec)
!!$  !StartMJD = JD2MJD(JD) - 90.0D0 / 1440.0D0   ! minus 90 min
!!$  !EndMJD   = StartMJD + 180.0D0  / 1440.0D0   ! plus 180 min
!!$  StartMJD = JD2MJD(JD) - STDperiod / 1440.0D0   ! minus STDperiod min.
!!$  EndMJD   = JD2MJD(JD) + STDperiod / 1440.0D0   ! plus STDperiod min.

  if (debug) write(*,*) 'StartMJD, EndMJD ',   StartMJD, EndMJD

#ifdef __COSMO__
  ! COSMO: Compute StartMJD and EndMJD according to STD input file name
  i = index(STDfile,'slv_slant')
  read(STDFile(i+10:i+13),'(i4)') Year
  read(STDFile(i+15:i+17),'(i3)') Doy
  read(STDFile(i+19:i+20),'(i2)') Hour
  !write(*,*) 'Reaf : ', Year, Doy, Hour

  call DOY2GREG(Doy, Year, Month, Day)
  Minit = 0
  Sec = 0.0D0
  JD = JulianDate(Year, Month, Day, Hour, Minit, Sec)
  !StartMJD = JD2MJD(JD) - 90.0D0 / 1440.0D0   ! minus 90 min
  !EndMJD   = StartMJD + 180.0D0  / 1440.0D0   ! plus 180 min
  StartMJD = JD2MJD(JD) - STDperiod / 1440.0D0   ! minus STDperiod min.
  EndMJD   = JD2MJD(JD) + STDperiod / 1440.0D0   ! plus STDperiod min.
#endif

  if (present(DataStartMJD) .and. present(DataEndMJD)) then
     DataStartMJD = 500000.0   ! date in far future: 31.10.3227
     DataEndMJD   =      0.0
  end if

  MaxBuf = 1
  !JD = MJD2JD(StartMJD)
  do n=1, MaxBuf
!!$     call GregorianDate( JD, Year, Month, Day, Hour, Minit, Sec)
!!$     call GREG2DOY (int(Year), int(Month), int(Day), Doy)
!!$     write(STDFile(i+10:i+13),'(i4.4)') Year
!!$     write(STDFile(i+15:i+17),'(i3.3)') Doy
!!$     write(STDFile(i+19:i+20),'(i2.2)') Hour
!!$     write(*,*) n, JD, 'Read file ',  trim(STDfile)
     !if (len_trim(std_obs_file(n)) < 1) cycle
     !STDFile = trim(std_obs_file(n))
     call ImportESTDarr(STDFile, MaxSTDobs, STDbuf(n)%STDlist, Header,       &
                        StatHash, StatList, NSat, SatList, StartMJD, EndMJD)

     ! Find MJD of first and last observation
     if (present(DataStartMJD) .and. present(DataEndMJD)) then
        if(Header%First < DataStartMJD) DataStartMJD = Header%First
        if(Header%Last > DataEndMJD) DataEndMJD = Header%Last
     end if

     !JD = JD + 60.0/ 1440.0_dp  ! plus 1 h
  end do

  !---------------------------------------------------------
  ! Azimuth and elevation of the slants are given in degrees
  ! => convert to radian
  !---------------------------------------------------------
  !write(*,*) 'ImportESTDarr - Array'
  !write(*,*) 'Slant file version ', Header%Version
  !write(*,*) 'Number of slants read : ', Header%Nslants
  if (debug) write(*,*) 'Convert angles'
  NSTD = 0
  do n=1, MaxBuf
#ifdef __COSMO__
     if (allocated(STDbuf(n)%STDlist)) then
#else
     if (associated(STDbuf(n)%STDlist)) then
#endif
        NSTD = NSTD + ubound(STDbuf(n)%STDlist,1)
        do i=1, ubound(STDbuf(n)%STDlist,1)
           !write(*,*) STDlist(i)%time, STDlist(i)%azimuth,   &
           !            STDlist(i)%elevation, STDlist(i)%slant
           STDbuf(n)%STDlist(i)%elevation = STDbuf(n)%STDlist(i)%elevation  &
                                            * deg2rad
           STDbuf(n)%STDlist(i)%azimuth   = STDbuf(n)%STDlist(i)%azimuth    &
                                            * deg2rad
           !write(*,*) STDlist(i)%time, STDlist(i)%azimuth,   &
           !           STDlist(i)%elevation, STDlist(i)%slant
           !write(*,*) '---------------------------------'


           ! initialize variables which could not be read from file
           STDbuf(n)%STDlist(i)%GNSS = 0
           STDbuf(n)%STDlist(i)%DryMap = STDbuf(n)%STDlist(i)%slant /  &
                                         STDbuf(n)%STDlist(i)%zslant
           STDbuf(n)%STDlist(i)%STDerr = -99.99_wp
           STDbuf(n)%STDlist(i)%STDacc = -99
           STDbuf(n)%STDlist(i)%SWD    = -99.99_wp
           STDbuf(n)%STDlist(i)%SIWV   = -99.99_wp
           STDbuf(n)%STDlist(i)%mSTD   = -99.99_wp
           STDbuf(n)%STDlist(i)%mSWD   = -99.99_wp
           STDbuf(n)%STDlist(i)%mSIWV  = -99.99_wp
           STDbuf(n)%STDlist(i)%PE     = -99
           STDbuf(n)%STDlist(i)%STDline= -99.99_wp
           STDbuf(n)%STDlist(i)%Press  = -99.99_wp
           STDbuf(n)%STDlist(i)%PosEll = -999.99_wp
           STDbuf(n)%STDlist(i)%Name   = ''
#ifdef __COSMO__
           STDbuf(n)%STDlist(i)%status     = -99
           STDbuf(n)%STDlist(i)%check      = -99
           STDbuf(n)%STDlist(i)%icenter    = -99
           STDbuf(n)%STDlist(i)%processing = -99
#endif

        end do
     end if
  end do

  !if (ubound(STDlist,1) .gt. 1) then
  if (NSTD .gt. 1) then
     !-------------------------------------------------------------------------
     ! Copy the slant data to the "SlantData" structure and
     ! set the station pointer
     !-------------------------------------------------------------------------

    !write(*,*) 'Copy data'
     !allocate( STD(1:ubound(STDlist,1)) )
     allocate( STD(1:NSTD) )

     Nlost = 0
     j = 0
     do n=1, MaxBuf
        do i=1, ubound(STDbuf(n)%STDlist,1)
           j = j + 1   ! array index used by STD => STD(j)
           STD(j)%Obs = STDbuf(n)%STDlist(i)
           !STD(j)%assimilate = .false.   ! assimilate selected slants only
           STD(j)%Nmod = 0
           STD(j)%Nup  = 0
           STD(j)%Ntot = 0

           ! Test: site is used to save netCDF record for ZTDs,
           !       site = 0 => no record availabe
           STD(j)%Obs%site = 0

           !write(*,*) 'STDlist(i)%station = ', STDlist(i)%station
           !write(*,*) 'GNSShash(STDlist(i)%station) = ', GNSShash(STDlist(i)%station)
           !write(*,*) 'CoordList... & Name = ',  CoordList(GNSShash(STDlist(i)%station))%Name
           if (GNSShash(STDbuf(n)%STDlist(i)%station) .gt. 0) then
              !-----------------------------------
              ! Set pointer to station information
              !-----------------------------------
              STD(j)%Station =>  CoordList(GNSShash(  &
                   STDbuf(n)%STDlist(i)%station))
              STD(j)%Obs%Name = STD(j)%Station%Name
           else
              !-------------------------------------------
              ! No valid GNSS station informatin available
              ! >>> Error message
              !-------------------------------------------
              nullify(STD(j)%Station)
              write(*,*) 'Station not found: ', STDbuf(n)%STDlist(i)%station
              Nlost = Nlost + 1
           end if

        end do
     end do

     if (j .ne. NSTD) then
        write(*,*) 'read_std_obs> Array index j should be eqal to NSTD, but'
        write(*,*) 'read_std_obs> j = ', j, ' NSTD = ', NSTD
     else
        write(*,*) 'read_std_obs> Number of slants read: ', NSTD
     end if

     if (Nlost .gt. 0) then
        write(*,*) 'read_std_obs> There are slants with no station infos ', &
                                  Nlost
     end if

     if (associated(STDlist)) deallocate(STDlist)
     do n=1, 4
#ifdef __COSMO__
        if (allocated(STDbuf(n)%STDlist)) deallocate(STDbuf(n)%STDlist)
#else
        if (associated(STDbuf(n)%STDlist)) deallocate(STDbuf(n)%STDlist)
#endif
     end do

!!$     do i=1, 10
!!$        if( associated(STD(i)%Station) ) then
!!$           write(*,*) STD(i)%Obs%station,  STD(i)%Station%Name, STD(i)%Obs%slant
!!$        end if
!!$     end do

     if (NSTD .gt. 0) then
        Nlost = 0
        do i=1, ubound(STD, 1)
           if( .not. associated(STD(i)%Station) ) then
              Nlost = Nlost + 1
              write(*,*) 'read_std_obs>  Station not assoc. ',    &
                          STD(i)%Obs%station
           end if
        end do
        !if (Nlost .gt. 0)
        if (debug) write(*,*) 'read_std_obs> N Stations not assoc. ', Nlost
     end if

  else
     !----------------------
     ! No STD data available
     !----------------------
     write(*,*) 'read_std_obs> No STD data read, Nslants = ', Header%Nslants
  end if

  if (verbose >= 1) write(*,*) 'read_std_obs> ... end'

 end subroutine read_std_obs


 !---------------------------------------------------------------------
 ! subroutine std_path_coord
 !---------------------------------------------------------------------
 !
 !> @brief Defines the supporting points on the signal path
 !>
 !> <b> call std_path_coord(STD, htop, err) </b>
 !>
 !> Driver routine for std_path_coord_sat and std_path_coord_line which
 !> should be called from the model interface. Direct calls to
 !> std_path_coord_sat or std_path_coord_line are not recommended.
 !>
 !> @param[in,out] STD       derived type which contains all information
 !>                          about one STD
 !> @param[in]     htop      model height in meter
 !> @param[out]    err       return code, err = 0 - no errors
 !> @param[in]     model_id  optional, model identifier
 !> @param[in]     model_nz  optional, number of full levels in model
 !> @param[in]     model_dx  optional, estimated horizontal model grid spacing,
 !>                          degrees
 !
 !---------------------------------------------------------------------
 subroutine std_path_coord (STD, htop, err, model_id, model_nz, model_dx, opt)

 type (SlantData)  ,intent(inout) :: STD
 real(wp)          ,intent(in)    :: htop
 integer           ,intent(out)   :: err
 integer, optional ,intent(in)    :: model_id
 integer, optional ,intent(in)    :: model_nz
 real(wp),optional ,intent(in)    :: model_dx
 integer, optional ,intent(in)    :: opt

 integer  :: nz ! number of full levels in model
 real(wp) :: dx ! horizontal grid spacing, m
 integer  :: mid ! model id

 !write(*,*) 'std_path_coord> Htop = ', htop, ', ztd_col = ', ztd_col

 err = 0

 if (.not. associated(STD%Station)) then
    ! no station information available, return with error
    err = 1
    return
 end if

 ! Initialize some variables in STD
 STD%Obs%Press  = invalid
 STD%Obs%PosEll = invalid
 STD%Obs%Name   = ''


 if (ztd_col > 0) then
    !-------------------------------------------------------
    ! fast ZTD processing,
    ! no supporting points required, except station position
    !-------------------------------------------------------
!   replace with new routine for ZTD
!   call std_path_coord_line (STD, mid, htop, nz, dx, err)
    call ztd_path_coord      (STD, err)

 else
    !----------------------------------------------------------------------
    ! STD processing, neighbored columns for all supporting points required
    !----------------------------------------------------------------------
    if (present(model_id)) then
       mid = model_id
    else
       mid = MO_UNKNOWN
    end if

    ! set  mumber of full levels
    if (present(model_nz)) then
       nz = model_nz
    else
       if (present(model_id)) then
          if (model_id == MO_COSMO) then
             nz = 50
          else
             ! unknown model, unknown mumber of full levels
             nz = -1
          end if
       else
          if (htop < 23000.0) then
             ! assume COSMO-DE
             nz = 50
          else
             ! unknown model, unknown mumber of full levels
             nz = -1
          end if
       end if
    end if
    ! set horizontal grid spacing
    if (present(model_dx)) then
       dx = model_dx * 111200.0_wp  ! horizontal grid spacing in m
    else
       if (present(model_id)) then
          if (model_id == MO_COSMO) then
             dx = 2800
          else
             ! unknown model, unknown mumber of full levels
             dx = -1.0
          end if
       else
          if (htop < 23000.0) then
             ! assume COSMO-DE
             dx = 2800
          else
             ! unknown model, unknown mumber of full levels
             dx = -1.0
          end if
       end if
    end if

    !write(*,*) 'allocated(Plevel) .and. NStepOpt  = ', allocated(Plevel), NStepOpt
    if (.not. allocated(Plevel) .and. NStepOpt == 3) then
       call ProfileSetup(nz)

       if (verbose >= 3) then
          write(*,*) 'Use new namelist options:'
          write(*,*) 'HScaleP  = ', HScaleP
          write(*,*) 'HScaleP2 = ', HScaleP2
          write(*,*) 'NStepOpt = ', NStepOpt
          write(*,*) 'Hlevel = ', Hlevel
          write(*,*) 'Plevel = ', Plevel
          write(*,*) 'Heights = ', Heights
          write(*,*) 'Pheights = ', Pheights
          write(*,*) 'Hpoints = ', Hpoints
          write(*,*) 'Ppoints = ', Ppoints
       end if
    end if

    if (UseRaytracer >= 0) then
       ! compute supporting points up to the GNSS satellite
       !write(*,*) 'call std_path_coord_sat with mid, htop, nz, dx = ', &
       !     mid, htop, nz, dx
       call std_path_coord_sat (STD, mid, htop, nz, dx, err, opt)
       !call std_path_coord_sat_old (STD, htop, err)
    else
       ! for backward compatibility only: required by "std_delay_line"
       call std_path_coord_line (STD, mid, htop, nz, dx, err)
       !call std_path_coord_line_old (STD, htop, err)
    end if

 end if     ! fast ZTD processing

 end subroutine std_path_coord


 !---------------------------------------------------------------------
 ! subroutine ProfileSetup
 !---------------------------------------------------------------------
 !
 !> @brief Parameter setup for approximating model full levels
 !>
 !> <b> call ProfileSetup( model_nz ) </b>
 !>
 !> The supporting points along the signal path need to be defined without
 !> any access to the model fields, e.g. HHL. Nevertheless, the supporting
 !> points should have a vertical structure which is comparable to the
 !> model full levels. Three namelist variables can be used to adjust the
 !> vertical profile of supporting points to the model profile in a flexible
 !> way: \n
 !> <b> Hlevel </b>
 !> is an array of model points in a given vertical range. The sum of all
 !> points in Hlevel must be equal to the number of the model full levels.
 !> Array elements less than 0 are not defined and are not counted as levels.
 !> \n
 !> <b> Heights </b>
 !> is an array of model full level heights at the upper boundaries of the
 !> ranges defined by Hlevel. Thsese heights should be copied from the
 !> model documentation.
 !> \n
 !> <b> Hpoints </b>
 !> is an array with the number of supporting points in each vertical range
 !> defined in Hlevel and Heights. The number of array elements > 0 must be
 !> equal to the corresponding elements in Hlevel. Array elements < 0 are not
 !> defined and will not be used. Hpoints can be used to define an arbitrary
 !> number of supporting points in each vertical range which is independent
 !> from the model levels.
 !> \n
 !> Hlevel, Heights and Hpoints are used to approximate the model levels and
 !> to define the supporting points inside the model. For extending the
 !> supporting points up to Hmax the namelist variable
 !> <b> NStepVertTop </b> is used.
 !> \n
 !> This routine reads Hlevel, Heights, Hpoints and NStepVertTop and
 !> allocates three new arrays Plevel, Pheights and Ppoints which have
 !> the exact size of the vertical ranges, including the extended range
 !> above the model top. These arrays don't have any undefined elements < 0
 !> a will be used by the routine ProfileSquared to compute the heights of
 !> the supporting points. NStepVertTop is appended to Hlevel and Hpoints
 !> in order to create one set of parameters which can be used for the
 !> whole profile up to Hmax. \n
 !> The arrays Plevel, Pheights and Ppoints are used as the reference in
 !> zenith direction and must not be modified! For computing supporting
 !> points in slant direction a local copy of Ppoints with an increased
 !> number of points is used.
 !>
 !> @param[in]  model_nz  number of model full levels
 !> @param[in]  pp        print profile, optional
 !>
 !> Access to namelist parameters defined as module variables:
 !>
 !> @param[in]  Hlevel    array of model full levels in a given vertical range
 !> @param[in]  Heights   array of model full level heights in m
 !> @param[in]  Hpoints   array with the number of supporting points in
 !>                       each vertical range
 !> @param[in]  NStepVertTop  number of supporting points in the extended region
 !>                           above the model top
 !> @param[in]  NStepVertMod  number of supporting points inside the model
 !>
 !> Access to module variables:
 !>
 !> @param[out] Plevel   array of model full levels and levels above the model
 !>                      top for all vertical ranges
 !> @param[out] Pheights array of upper heights for all ranges defined
 !>                      in Plevel
 !> @param[out] Ppoints  array with the number of supporting points in
 !>                      each vertical range (in zenith direction)
 !
 !---------------------------------------------------------------------
 subroutine ProfileSetup( model_nz, pp )

 integer, intent(in)           :: model_nz
 logical, intent(in), optional :: pp          ! print profile

 integer                              :: Nint, i, j, n
!integer, dimension(:), allocatable   :: IdxProf, idx, ModIdx
 integer                              :: Ntot, shift
 real (wp)                            :: H0, Hs
!real (wp), dimension(:), allocatable :: H0int
 logical                              :: PProfile
 real (wp), dimension(:), allocatable :: Hprofile
 !----------------------------------------------------------

 if (verbose == 3) write(*,*) 'ProfileSetup> Start ...'

 ! Print reference profile to screen, default: don't print
 if ( present(pp) ) then
    PProfile = pp
 else
    PProfile = .false.
 end if

 ! get number of intervals within profile (elements > 0)
 Nint = count(Hlevel > 0)

 if (count(Heights > 0.0) /= Nint .or. count(Hpoints > 0) /= Nint) then
    call model_abort (-1,-1,                         &
       'GPSGB: inconsistent namelist settings in Hlevel, Hscale and Hpoints', &
       'ProfileSetup'                 )
 end if

 ! check if Hlevel has "model_nz" total levels
 if (sum(Hlevel(1:Nint)) /= model_nz) then
    call model_abort (-1,-1,                         &
       'GPSGB: Sum of levels in Hlevel not equal to model full levels', &
       'ProfileSetup'                 )
 end if

 ! working copy of namelist settings plus one interval above model top
 if (allocated(Plevel)) deallocate(Plevel)
 if (allocated(Pheights)) deallocate(Pheights)
 if (allocated(Ppoints)) deallocate(Ppoints)
 allocate( Plevel  (Nint+1) )
 allocate( Pheights(Nint+1) )
 allocate( Ppoints (Nint+1) )
 Plevel(1:Nint)   = Hlevel(1:Nint)
 Pheights(1:Nint) = Heights(1:Nint)
 Ppoints(1:Nint)  = Hpoints(1:Nint)
 Plevel(Nint+1)  = NStepVertTop
 Pheights(Nint+1)= Hmax
 Ppoints(Nint+1) = NStepVertTop

 if (PProfile) then
    ! print reference profile

!!$    write(*,*) 'ProfileSetup> Hlevel  = ', Hlevel
!!$    write(*,*) 'ProfileSetup> Plevel  = ', Plevel
!!$    write(*,*) 'ProfileSetup> Hscale  = ', Heights
!!$    write(*,*) 'ProfileSetup> Pscale  = ', Pheights
!!$    write(*,*) 'ProfileSetup> Hpoints = ', Hpoints
!!$    write(*,*) 'ProfileSetup> Ppoints = ', Ppoints
!!$    write(*,*) 'ProfileSetup> IdxProf = ', IdxProf
!!$    write(*,*) 'ProfileSetup> idx     = ', idx
!!$    write(*,*) 'ProfileSetup> ModIdx  = ', ModIdx

    write(*,*)
    write(*,*) 'Reference profile used to define supporting points ', &
               'on the GNSS signal path:'

    ! Compute profile: height of supporting points along the signal path
    Ntot = sum(Ppoints)
    allocate( Hprofile(1:Ntot) )

    n = 0
    do i=1, size(Ppoints)
       ! process intervals

       ! Compute H0 and Hs for the current profile
       if (i == 1) then
          shift = 0
          H0 = 0.0_wp                     ! first point at sea level
          Hs = Pheights(i) / (Ppoints(i)-1)**2
       else
          shift = 2
          H0 = Pheights(i-1)
          Hs = (Pheights(i)-Pheights(i-1)) / (Ppoints(i)+shift-1)**2
       end if

       ! Compute height of supporting points within current interval
       do j=shift, Ppoints(i)+shift-1
          n = n + 1
          Hprofile(n) = H0 + Hs * j**2
       end do

    end do

    do i=1, size(Hprofile)
       write(*,'(tr4,a,i3,tr2,f12.2)') 'i, H(i) = ', i, Hprofile(i)
    end do

 end if

 if (verbose == 3) write(*,*) 'ProfileSetup> ... end'

 end subroutine ProfileSetup


 !---------------------------------------------------------------------
 ! subroutine ProfileSquared
 !---------------------------------------------------------------------
 !
 !> @brief Computes the heights of the supporting points along the signal path
 !>
 !> <b> call ProfileSquared(STD, mod_htop, StatH, Hprofile) </b>
 !>
 !> The heights of the supporting points between the station and Htop are
 !> computed, i.e. the whole profile including the range above the model top.
 !> The number of points given in Ppoints is adjusted for STDs using
 !> STD%Nmod and STD%Nup. For ZTDs the number of points defined in the
 !> namelist is used. \n
 !> The reference profile defined by Plevel and Pheights is rescaled for the
 !> station height in order to follow the hybrid levels.
 !>
 !> tested with Hpoints = (/35, 17, 1/)  - works and ZTDs aren't too bad \n
 !>                     = (/ 5,  5, 1/)  - runs, ZTDs still surprisingly good\n
 !>                     = (/ 1,  1, 1/)  - crashes \n
 !>
 !>
 !> @param[in]  STD   derived type which contains all information about
 !>                   one STD, required to access STD%Nmod and STD%Nup
 !> @param[in]  mod_htop  model top full level
 !> @param[in]  StatH     station height = height of lowest supporting point
 !> @param[out] Hprofile  height profile of supporting points:
 !>                       First point at the GNSS station (h=StatH), last
 !>                       point at h=Hmax. (No point at the satellite!)
 !>
 !> Access to module variables:
 !>
 !> @param[in]  Plevel   array of model full levels and levels above the model
 !>                      top for all vertical ranges
 !> @param[in]  Pheights array of upper heights for all ranges defined
 !>                      in Plevel
 !> @param[in]  Ppoints  array with the number of supporting points in
 !>                      each vertical range (in zenith direction)
 !>
 !---------------------------------------------------------------------
 subroutine ProfileSquared(STD, mod_htop, StatH, Hprofile)

 type (SlantData) ,intent(in)                      :: STD       ! STD data
 real (wp), intent(in)                             :: mod_htop
 real (wp), intent(in)                             :: StatH
 real (wp), dimension(:), allocatable, intent(out) :: Hprofile

 integer   :: n, i, j, imax
 integer   :: Ntot, Nmod, shift
 real (wp) :: H0, Hs, Delta0, DeltaS, Dold, f

 real (wp), dimension(:), allocatable :: Lheights
 integer, dimension(:), allocatable   :: Lpoints
 !---------------------------------------------------------------------

 if (verbose == 3) write(*,*) 'ProfileSquared> Start ...'

 ! Adjust the number of supporting points along the signal path according
 ! to the elevation (=> STD%Nmod, STD%Nup)
 n = size(Ppoints)
 allocate( Lpoints(n) )
 Lpoints = Ppoints
 Nmod = sum(Ppoints(:n-1))
 !write(*,*) 'ProfileSquared> n, Ntot, STD%Nmod = ', n, Ntot, STD%Nmod

 if (Nmod /= STD%Nmod) then
    Lpoints(n) = STD%Nup
    do i=1, n-1
       Lpoints(i) = int(Ppoints(i) * STD%Nmod/Nmod)
    end do
    do i=1, STD%Nmod - sum(Lpoints(:n-1))
       Lpoints(i) =  Lpoints(i) + 1
    end do
    if ( sum(Lpoints(:n-1)) /= STD%Nmod) then
       call model_abort (-1,-1,                                              &
            'GPSGB: Sum of levels in Hlevel not equal to model full levels', &
            'ProfileSquared'                 )
    end if
 end if

 ! Rescale lower model levels to station height
 ! Assume that model levels above tropopause are constant, i.e. not
 ! terrain following
 ! Lheights is the local version of Pheights with rescaled lower levels
 allocate( Lheights(n) )
 Lheights = Pheights
 if (StatH > 5.0 .and. n > 1) then
    ! Rescale only if sation height is above sea level and if there are at
    ! least two levels used for profile approximation.

    ! Find first intervall that reaches above the tropopause (12000 m)
    do i=1, n
       if (Pheights(i) > 12000.0) then
          imax = i
          exit
       end if
    end do

    ! Rescale levels below the tropopause and preserve the relative
    ! layer thickness
    Delta0 = Pheights(imax)
    DeltaS = Delta0 - StatH
    f = DeltaS/Delta0
    do i=imax, 2, -1
       Dold =  Pheights(i) - Pheights(i-1)
       Lheights(i-1) = Lheights(i) - f*Dold
    end do

 end if

 ! Compute profile: height of supporting points along the signal path
 ! First point at the GNSS station (n=1, h=StatH), last point at h=Hmax
 ! (n = STD%Nmod + STD%Nup).
 ! (No point at the satellite!)
 Ntot = sum(Lpoints)
 allocate( Hprofile(1:Ntot) )
 Hprofile = invalid

 if (Ntot /= STD%Nmod + STD%Nup) then
    call model_abort (-1,-1,                   &
            'GPSGB: Wrong size of Hprofile',   &
            'ProfileSquared'                 )
 end if

 n = 0
 do i=1, size(Lpoints)
    ! process intervals

    ! Compute H0 and Hs for the current profile
    if (i == 1) then
       shift = 0
       H0 = StatH                     ! first point at station
       Hs = (Lheights(i)-StatH) / real((Lpoints(i)+shift-1)**2, wp)
    else
       shift = 2
       H0 = Lheights(i-1)
       Hs = (Lheights(i)-Lheights(i-1)) / real((Lpoints(i)+shift-1)**2, wp)
    end if
    !write(*,*) 'ProfileSquared> H0, Hs = ', H0, Hs
    ! Compute height of supporting points within current interval
    do j=shift, Lpoints(i)+shift-1
       n = n + 1
       Hprofile(n) = H0 + Hs * real(j**2, wp)
    end do

 end do

 !write(*,*) 'ProfileSquared> Pheights, Lheights, Hprofile = ', char(10), &
 !     Pheights, char(10), &
 !     Lheights, char(10), &
 !     Hprofile

 if (verbose == 3) write(*,*) 'ProfileSquared> ... end'

 end subroutine ProfileSquared


 !---------------------------------------------------------------------
 ! subroutine std_path_coord_sat
 !---------------------------------------------------------------------
 !
 !> @brief Defines the supporting points on the signal path for raytracing
 !>
 !> <b> call std_path_coord_sat(STD, htop, err) </b>
 !>
 !> This is the first part of the STD operator which computes the coordinates
 !> of the supporting points along the signal path. These coordinates are
 !> required by the interface to the weather model for collecting the model
 !> variables in the vicinity of the supporting points.  The second part
 !> of the STD operator (std_delay) gets access to this subset of model data,
 !> starts the raytracer, integrates alog the curved signal path and provides
 !> the STD. \n
 !> This  routine is called for each single slant, i.e. it processes only one
 !> single slant.
 !>
 !> @param[in,out] STD       derived type which contains all information about
 !>                          one STD
 !> @param[in]     mod_id    model identifier: MO_COSMO, MO_ICON, ...
 !> @param[in]     mod_htop  model top height in meter (full levels)
 !> @param[in]     mod_nz    number of vertical model levels (full levels)
 !> @param[in]     mod_dx    estimated horizontal grid spacing in meter
 !> @param[out]    err       return code, err = 0 - no errors
 !
 !---------------------------------------------------------------------
 subroutine std_path_coord_sat (STD, mod_id, mod_htop, mod_nz, mod_dx, &
                                err, opt)

 type (SlantData) ,intent(inout) :: STD       ! STD data
 integer          ,intent(in)    :: mod_id    ! model identifier
 real(wp)         ,intent(in)    :: mod_htop  ! height in m of model top level
 integer          ,intent(in)    :: mod_nz    ! number of vertical levels
 real(wp)         ,intent(in)    :: mod_dx    ! horizontal grid spacing, m
 integer          ,intent(out)   :: err       ! error code
 integer, optional,intent(in)    :: opt       ! option for re-computing

 ! mod_htop - altitude of model top level in meters  ??
 !            frame of reference ?? height above geoid? ellipsoid?

!real (wp) :: px, py, pz, gx, gy, gz
!real (wp) :: TestLat, TestLon, DeltaLat, DeltaLon

 integer :: a
!integer :: stat, unit
!logical :: opened
 !integer :: Nerr, Ncosmo, Nperiod, Nhorz, Nstd, Nperr, Nassim

 !integer                              :: inCOSMO

 integer                   :: StepOpt
 type (Line)               :: ray
 real (wp), dimension(1:3) :: CrossPt, point, wgspoint  !, point2,  wgspoint2
 real (wp), dimension(1:3) :: SatelliteCoord
 real (wp)                 :: RayLambda
 real (wp)                 :: DeltaP, P, lnP0, Ps, Pe, h, c
 real (wp), parameter      :: P0 = 1024.0_wp  ! reference pressure at sea level
 real (wp), parameter      :: rE = 6378137.0_wp ! radius of Earth in m
 real (wp)                 :: Hdist             ! horizontal distance
 character (len=10)        :: GNSS
 real (wp), dimension(:), allocatable :: Hprofile
 real (wp), dimension(:), allocatable :: SlantSteps
 !---------------------------------------------

 if (verbose >= 1) write(*,*) 'std_path_coord_sat> Start ...'

 StepOpt = NStepOpt
 if (present(opt)) then
    if (opt >= 10) StepOpt = opt
 end if

 if (verbose >= 3) then
    write(*,*) 'std_path_coord_sat> Recompute - StepOpt = ', StepOpt
    write(*,*) 'std_path_coord_sat> mod_id, mod_htop, mod_nz, mod_dx = ', &
                mod_id, mod_htop, mod_nz, mod_dx

    write(*,*) 'std_path_coord_sat> Station ', std% Station%SName, char(10),  &
            'Station ID = ', std% Station%ID, char(10),                 &
            'Geoid/Ellipsoid lon = ', std% Station%CoordGeo(1)*rad2deg, &
                       std% Station%CoordEll(1)*rad2deg, char(10),      &
            'Geoid/Ellipsoid lat = ', std% Station%CoordGeo(2)*rad2deg, &
                       std% Station%CoordEll(2)*rad2deg, char(10),      &
            'Geoid/Ellipsoid hei = ', std% Station%CoordGeo(3) ,        &
                       std% Station%CoordEll(3), char(10),              &
                       'Cartesian X/Y/Z m   = ',  std% Station%CoordCart

    write(*,*) 'std_path_coord_sat> Slant: Elevation ', &
         STD%Obs%elevation*rad2deg, &
          char(10), 'Azimuth ', STD%Obs%azimuth*rad2deg
 end if

 if (StepOpt == 10) then
    ! Re-compute new signal path using the old SlantSteps
    if (.not. (associated(std% SlantSteps)      .and.         &
               associated(std% SlantPointsCart) .and.         &
               associated(std% SlantPointsEll)       ) ) then
       ! Error: Re-computing of supporting points possible only if
       !        data from previous call are available
       call model_abort (-1,-1,                                           &
            'GPSGB: Re-computing of path not possible, old data missing', &
            'std_path_coord_sat'                 )
    end if

    ! Copy old SlantSteps
    allocate( SlantSteps(size(std% SlantSteps)) )
    SlantSteps(:) = std% SlantSteps(:)

    ! Delete old slant path
    deallocate ( std% SlantSteps      )
    deallocate ( std% SlantPointsCart )
    deallocate ( std% SlantPointsEll  )

 end if

 !-----------------------------------------
 ! Slant loop
 ! Slants are sorted by station and by time
 !-----------------------------------------
 if (.not. associated(std% Station) ) then
    write(*,*) 'std_path_coord_sat> No station data available => skip'
    !Nerr = Nerr + 1
    return
 end if

 if (std% Station%ID < 0) then
    write(*,*) 'std_path_coord_sat> Error : invalid station id = ', &
                std% Station%ID
    !Nerr = Nerr + 1
    return
 end if

 ! Compute satellite position and define a straight line from the receiver
 ! to the satellite
 GNSS = 'GPS'   ! add GLONASS, Galileo, Compass, ...
 call SatPos (STD%Station%CoordEll, STD%Obs%azimuth, &
                STD%Obs%elevation, GNSS, SatelliteCoord)
 call LineDef(ray,  std% Station%CoordCart, SatelliteCoord, 3)

 ! ??????????????
 ! mod_htop auch transformieren ???
 ! ??????????????
 call CrossLineEllips(WGS84Param, mod_htop, ray,                &
                      CrossPt, RayLambda)
 call Cart2Ellips(CrossPt(1),  CrossPt(2),  CrossPt(3),         &
      point(1), point(2), point(3), WGS84Param)
 !-----------------------------------------------------
 ! The signal path within the COSMO grid is now defined
 !-----------------------------------------------------

 !--------------------------------------------------------
 ! Scale the number of supporting points on the slant path
 ! according to the scaling option "NStepOpt".
 !--------------------------------------------------------
 select case (StepOpt)
    case (1)
       !--------------------------------------------------------
       ! Scale the number of supporting points on the slant path
       ! with 1/sin(elevation)
       ! For elevtions below 2 degrees, the number of supporting points is
       ! not further increased. (1/sin(2°) ~= 28)
       !--------------------------------------------------------
       if (verbose >= 2) &
            write(*,*) 'std_path_coord_sat> Start computing supporting points'
       if (STD%Obs%elevation .ge. 2.0*deg2rad) then
          c = sin(STD%Obs%elevation)
       else
          c = sin(2.0_wp * deg2rad)
       end if
       ! Empirical factor:
       ! ZTDs require a number of supporting points which is not below the
       ! number of vertical model layers. At lower elevations the number needs
       ! not to be increased that fast.
       c = 2.0*c
       std% Nmod = int(NStepVertMod / c)
       std% Nup  = int(NStepVertTop / c)
       std% Ntot = std% Nmod + std% Nup + 1

    case (2, 3)
       !--------------------------------------------------------
       !
       !--------------------------------------------------------

       ! horizontal distance between GNSS station and exit point of slant
       ! Hdist in meter (approx. sphere)
       Hdist = rE*( pi05 - STD%Obs%elevation -                           &
                    asin( cos(STD%Obs%elevation)*rE/(rE+mod_htop) ) )

       ! mod_dx = 2.8 km: COSMO grid spacing
       ! mod_dx =  13 km: ICON grid spacing
       std% Nmod = NStepVertMod + int(1.5 * Hdist/mod_dx)

       ! horizontal distance between GNSS station and exit point of slant
       ! Hdist in meter
       Hdist = rE*( pi05 - STD%Obs%elevation -                           &
                    asin( cos(STD%Obs%elevation)*rE/(rE+Hmax) ) )

       ! add one supporting point for 50 km
       std% Nup  = NStepVertTop + int(Hdist/50000.0_wp)

       if (NStepOpt < 3) then
          std% Ntot = std% Nmod + std% Nup + 1
       else
          std% Ntot = std% Nmod + std% Nup + 1
       end if
    case (10)
       ! Re-compute new signal path
       continue
    case default
       write(*,*) 'std_path_coord_sat> Unsupported option: NStepOpt = ', &
                  NStepOpt
       std% Nmod = NStepVertMod
       std% Nup  = NStepVertTop
       std% Ntot = std% Nmod + std% Nup + 1
 end select

 !--------------------------------------------------------------------
 ! Supporting points for the refractivity profile along the slant path
 !--------------------------------------------------------------------
 allocate ( std% SlantSteps(1:std% Ntot) )
 allocate ( std% SlantPointsCart(1:std% Ntot,1:3) )
 allocate ( std% SlantPointsEll(1:std% Ntot,1:3) )

 ! start at the GNSS receiver
 STD%SlantSteps(1)           = 0.0_wp  ! distance to receiver in m
 STD%SlantPointsCart(1,:)    = STD%Station%CoordCart
 STD%SlantPointsEll(1,:)     = STD%Station%CoordEll
 ! last point in the model grid
 STD%SlantSteps(STD%Nmod) = RayLambda ! distance to receiver in m
 STD%SlantPointsCart(STD%Nmod,:) = CrossPt
 ! last point at the satellite
 STD%SlantSteps(STD%Ntot) = PointDist(STD%Station%CoordCart, SatelliteCoord)
 STD%SlantPointsCart(STD%Ntot,:) = SatelliteCoord

 select case (StepOpt)
 case (1,2)
    ! Scale the density of supporting points with the air pressure
    ! P0 is the reference pressure at sea level, used to define
    ! a simple exponential profile: P(h) = P0*exp(-h/HscaleP)
    ! Use HscaleP for pressure profile within the model
    lnP0 = log( P0 )

    ! Pressure at the GPS station
    Ps = P0*exp(-STD%Station%CoordEll(3)/HscaleP)
    ! Pressure at the upper grid level
    Pe = P0*exp(-point(3)/HscaleP)

    DeltaP = (Ps - Pe) / real(STD%Nmod - 1)

    p =  Ps
    do a=2, STD%Nmod-1
       ! Compute supporting points at equidistant pressure levels
       p = p - DeltaP
       h =  -HscaleP*(log(p) - lnP0)
       call CrossLineEllips(WGS84Param, h, ray,         &
                            CrossPt, RayLambda)
       STD%SlantSteps(a) = RayLambda
       STD%SlantPointsCart(a,:) = CrossPt
    end do

    if ( STD%SlantSteps(STD%Nmod-1) .ge.                       &
         STD%SlantSteps(STD%Nmod)         ) then
       ! Check if supporting points are really below the upmost model level
       write(*,*) 'Error: Slant distance too large: ', &
            STD%SlantSteps(STD%Nmod-1), STD%SlantSteps(STD%Nmod)
    end if

    !------------------------------------------------------------
    ! Do the same for the slant path above the upmost model level
    !------------------------------------------------------------
    ! Last point with N > 0 at the height Hmax (N == 0 at the satellite)
    call CrossLineEllips(WGS84Param, Hmax, ray,         &
                         CrossPt, RayLambda)
    !call Cart2Ellips(CrossPt(1),  CrossPt(2),  CrossPt(3),         &
    !                 point(1), point(2), point(3), WGS84Param)
    STD%SlantSteps(STD%Ntot-1) = RayLambda ! last point on slant
    STD%SlantPointsCart(STD%Ntot-1,:) = CrossPt

    ! Pressure at the upper grid level
    ! Use HscaleP2 for pressure profile above the model top
    Ps = P0*exp(-point(3)/HscaleP2)
    ! Pressure at Hmax
    Pe = P0*exp(-Hmax/HscaleP2)

    DeltaP = (Ps - Pe) / real(STD%Nup)

    p =  Ps
    do a=STD%Nmod+1, STD%Ntot-2
       ! Compute supporting points at equidistant pressure levels
       p = p - DeltaP
       h =  -HscaleP2*(log(p) - lnP0)
       call CrossLineEllips(WGS84Param, h, ray,         &
                            CrossPt, RayLambda)
       STD%SlantSteps(a) = RayLambda
       STD%SlantPointsCart(a,:) = CrossPt
    end do

 case (3)
    ! use quadratic profile with several intervals

    ! Compute height above ellipsoid of all supporting points => Hprofile
    !write(*,*) 'std_path_coord_sat> call ProfileSquared'
    call ProfileSquared(std, mod_htop, std%Station%CoordEll(3), Hprofile)

    do a=2, STD%Ntot-1
       ! Compute supporting points on the ellipsoid at height "Hprofile"
       call CrossLineEllips(WGS84Param, Hprofile(a), ray, CrossPt, RayLambda)
       STD%SlantSteps(a) = RayLambda
       STD%SlantPointsCart(a,:) = CrossPt
    end do
 case (10)
    ! Re-compute new signal path using the old SlantSteps

    ! Copy old supporting points to new ray
    std% SlantSteps(2:STD%Ntot-1) = SlantSteps(2:STD%Ntot-1)

    ! Compute cartesian coordinates on ray
    do a=2, STD%Ntot-1
       call LinePos(ray, SlantSteps(a), STD%SlantPointsCart(a,:))
    end do

 end select

 ! So far, the distance to the GNSS receiver and the cartesian coordinates
 ! have been computed for all supporting ponts. Now compute the
 ! corresponding WGS84 coordinates.

 do a=2, std%Ntot
    !-----------------------------------------------------------------------
    ! SlantPointsCart - cartesian coordinates of a point on the slant path
    !                   with a distance = SlantSteps(a) to the GNSS receiver
    !-----------------------------------------------------------------------
    call Cart2Ellips( STD%SlantPointsCart(a,1),                &
                      STD%SlantPointsCart(a,2),                &
                      STD%SlantPointsCart(a,3),                &
                      wgspoint(1), wgspoint(2), wgspoint(3),   &
                      WGS84Param                             )
    !-----------------------------------------------------
    ! wgspoint - ellipsoidal coordinates of the same point
    !-----------------------------------------------------
    std% SlantPointsEll(a,:) = wgspoint
    !write(*,*) std%SlantPointsEll(a,1)*rad2deg,   &
    !           std%SlantPointsEll(a,2)*rad2deg, std%SlantPointsEll(a,3)

 end do

 if (verbose >= 3) then
    write(*,*) 'std_path_coord_sat> ',                  &
               'Coordinates of supporting points'
    do a=1, std%Ntot-1
       if (a==STD%Nmod) then
          write(*,*) 'i, SlantPointsEll = ', a,         &
                     std%SlantPointsEll(a,1:2)*rad2deg, &
                     std%SlantPointsEll(a,3),           &
                     ' <== model top'
       else
          write(*,*) 'i, SlantPointsEll = ', a,         &
                     std%SlantPointsEll(a,1:2)*rad2deg, &
                     std%SlantPointsEll(a,3)
       end if
    end do
    write(*,*) 'i, SlantPointsEll = ', a,         &
               std%SlantPointsEll(a,1:2)*rad2deg, &
               std%SlantPointsEll(a,3),           &
               ' <== GNSS satellite'
 end if

 if (verbose >= 1) write(*,*) 'std_path_coord_sat>  ... end'

end subroutine std_path_coord_sat


 !---------------------------------------------------------------------
 ! subroutine std_path_coord_sat_old
 !---------------------------------------------------------------------
 !
 !> @brief Defines the supporting points on the signal path for raytracing
 !>
 !> <b> call std_path_coord_sat_old(STD, htop, err) </b>
 !>
 !> This is the first part of the STD operator which computes the coordinates
 !> of the supporting points along the signal path. These coordinates are
 !> required by the interface to the weather model for collecting the model
 !> variables in the vicinity of the supporting points.  The second part
 !> of the STD operator (std_delay) gets access to this subset of model data,
 !> starts the raytracer, integrates alog the curved signal path and provides
 !> the STD. \n
 !> This  routine is called for each single slant, i.e. it processes only one
 !> single slant.
 !>
 !> @param[in,out] STD  derived type which contains all information about one STD
 !> @param[in]     htop model to height in meter
 !> @param[out]    err  return code, err = 0 - no errors
 !
 !---------------------------------------------------------------------
 subroutine std_path_coord_sat_old (STD, htop, err)

 type (SlantData) ,intent(inout) :: STD
 real(wp)         ,intent(in)    :: htop
 integer          ,intent(out)   :: err

 ! htop - altitude of model top level in meters  ??
 !        frame of reference ?? height above geoid? ellipsoid?

!real (wp) :: px, py, pz, gx, gy, gz
!real (wp) :: TestLat, TestLon, DeltaLat, DeltaLon

 integer :: a
!integer :: stat, unit
!logical :: opened
 !integer :: Nerr, Ncosmo, Nperiod, Nhorz, Nstd, Nperr, Nassim

 !integer                              :: inCOSMO

 type (Line)               :: ray
 real (wp), dimension(1:3) :: CrossPt, point, wgspoint  !, point2,  wgspoint2
 real (wp), dimension(1:3) :: SatelliteCoord
 real (wp)                 :: RayLambda
 real (wp)                 :: DeltaP, P, lnP0, Ps, Pe, h, c
 real (wp), parameter      :: P0 = 1024.0_wp  ! reference pressure at sea level
 real (wp), parameter      :: rE = 6378137.0_wp ! radius of Earth in m
 real (wp)                 :: Hdist             ! horizontal distance
 character (len=10)        :: GNSS

 !real (wp), dimension(:), allocatable :: SlantSteps
 !---------------------------------------------

 if (verbose >= 1) write(*,*) 'std_path_coord_sat> Start ...'

!!$  Nerr = 0
!!$  Ncosmo = 0
!!$  Nperiod = 0
!!$  Nhorz = 0
!!$  Nstd = 0
!!$  Nperr = 0

 if (verbose >= 3) then
    write(*,*) 'std_path_coord_sat> Station ', std% Station%SName, char(10),  &
            'Station ID = ', std% Station%ID, char(10),                 &
            'Geoid/Ellipsoid lon = ', std% Station%CoordGeo(1)*rad2deg, &
                       std% Station%CoordEll(1)*rad2deg, char(10),      &
            'Geoid/Ellipsoid lat = ', std% Station%CoordGeo(2)*rad2deg, &
                       std% Station%CoordEll(2)*rad2deg, char(10),      &
            'Geoid/Ellipsoid hei = ', std% Station%CoordGeo(3) ,        &
                       std% Station%CoordEll(3), char(10),              &
            'Cartesian X/Y/Z m   = ',  std% Station%CoordCart
 end if

 !-----------------------------------------
 ! Slant loop
 ! Slants are sorted by station and by time
 !-----------------------------------------
 if (.not. associated(std% Station) ) then
    write(*,*) 'std_path_coord_sat> No station data available => skip'
    !Nerr = Nerr + 1
    return
 end if

 if (std% Station%ID < 0) then
    write(*,*) 'std_path_coord_sat> Error : invalid station id = ', &
                std% Station%ID
    !Nerr = Nerr + 1
    return
 end if

 ! Compute satellite position and define a straight line from the receiver
 ! to the satellite
 GNSS = 'GPS'   ! add GLONASS, Galileo, Compass, ...
 call SatPos (STD%Station%CoordEll, STD%Obs%azimuth, &
                STD%Obs%elevation, GNSS, SatelliteCoord)
 call LineDef(ray,  std% Station%CoordCart, SatelliteCoord, 3)
 ! ??????????????
 ! htop auch transformieren ???
 ! ??????????????
 call CrossLineEllips(WGS84Param, htop, ray,                    &
                      CrossPt, RayLambda)
 call Cart2Ellips(CrossPt(1),  CrossPt(2),  CrossPt(3),         &
                  point(1), point(2), point(3), WGS84Param)
 !-----------------------------------------------------
 ! The signal path within the COSMO grid is now defined
 !-----------------------------------------------------

 !--------------------------------------------------------
 ! Scale the number of supporting points on the slant path
 ! according to the scaling option "NStepOpt".
 !--------------------------------------------------------
 select case (NStepOpt)
    case (1)
       !--------------------------------------------------------
       ! Scale the number of supporting points on the slant path
       ! with 1/sin(elevation)
       ! For elevtions below 2 degrees, the number of supporting points is
       ! not further increased. (1/sin(2°) ~= 28)
       !--------------------------------------------------------
       if (verbose >= 2) &
            write(*,*) 'std_path_coord_sat> Start computing supporting points'
       if (STD%Obs%elevation .ge. 2.0*deg2rad) then
          c = sin(STD%Obs%elevation)
       else
          c = sin(2.0_wp * deg2rad)
       end if
       ! Empirical factor:
       ! ZTDs require a number of supporting points which is not below the
       ! number of vertical model layers. At lower elevations the number needs
       ! not to be increased that fast.
       c = 2.0*c
       std% Nmod = int(NStepVertMod / c)
       std% Nup  = int(NStepVertTop / c)
       std% Ntot = std% Nmod + std% Nup + 1

    case (2)
       !--------------------------------------------------------
       !
       !--------------------------------------------------------

       ! horizontal distance between GNSS station and exit point of slant
       ! Hdist in meter
       Hdist = rE*( pi05 - STD%Obs%elevation -                           &
                    asin( cos(STD%Obs%elevation)*rE/(rE+htop) ) )

       ! 2.8 km: COSMO grid spacing
       ! => replace with grid spacing of the model!!!
       std% Nmod = NStepVertMod + int(1.5 * Hdist/2800.0_wp)

       !write(*,*) 'Hdist, Hdist/2800.0 ', Hdist, int(Hdist/2800.0_wp)

       ! horizontal distance between GNSS station and exit point of slant
       ! Hdist in meter
       Hdist = rE*( pi05 - STD%Obs%elevation -                           &
                    asin( cos(STD%Obs%elevation)*rE/(rE+Hmax) ) )

       ! add one supporting point for 50 km
       std% Nup  = NStepVertTop + int(Hdist/50000.0_wp)
       std% Ntot = std% Nmod + std% Nup + 1

    case default
       write(*,*) 'std_path_coord_sat> Unsupported option: NStepOpt = ', &
                  NStepOpt
       std% Nmod = NStepVertMod
       std% Nup  = NStepVertTop
       std% Ntot = std% Nmod + std% Nup + 1
 end select

 !--------------------------------------------------------------------
 ! Supporting points for the refractivity profile along the slant path
 !--------------------------------------------------------------------
 allocate ( std% SlantSteps(1:std% Ntot) )
 allocate ( std% SlantPointsCart(1:std% Ntot,1:3) )
 allocate ( std% SlantPointsEll(1:std% Ntot,1:3) )

 ! start at the GNSS receiver
 STD%SlantSteps(1)           = 0.0_wp  ! distance to receiver in m
 STD%SlantPointsCart(1,:)    = STD%Station%CoordCart
 ! last point in the model grid
 STD%SlantSteps(STD%Nmod) = RayLambda ! distance to receiver in m
 STD%SlantPointsCart(STD%Nmod,:) = CrossPt
 ! last point at the satellite
 STD%SlantSteps(STD%Ntot) = PointDist(STD%Station%CoordCart, SatelliteCoord)
 STD%SlantPointsCart(STD%Ntot,:) = SatelliteCoord

 ! Scale the density of supporting points with the air pressure
 ! P0 is the reference pressure at sea level, used to define
 ! a simple exponential profile: P(h) = P0*exp(-h/HscaleP)
 ! Use HscaleP for pressure profile within the model
 lnP0 = log( P0 )

 ! Pressure at the GPS station
 Ps = P0*exp(-STD%Station%CoordEll(3)/HscaleP)
 ! Pressure at the upper grid level
 Pe = P0*exp(-point(3)/HscaleP)

 DeltaP = (Ps - Pe) / real(STD%Nmod - 1)

 p =  Ps
 do a=2, STD%Nmod-1
    ! Compute supporting points at equidistant pressure levels
    p = p - DeltaP
    h =  -HscaleP*(log(p) - lnP0)
    call CrossLineEllips(WGS84Param, h, ray,         &
                         CrossPt, RayLambda)
    STD%SlantSteps(a) = RayLambda
    STD%SlantPointsCart(a,:) = CrossPt
 end do

 if ( STD%SlantSteps(STD%Nmod-1) .ge.                       &
      STD%SlantSteps(STD%Nmod)         ) then
    ! Check if supporting points are really below the upmost model level
    write(*,*) 'Error: Slant distance too large: ', &
         STD%SlantSteps(STD%Nmod-1), STD%SlantSteps(STD%Nmod)
 end if

 !------------------------------------------------------------
 ! Do the same for the slant path above the upmost model level
 !------------------------------------------------------------
 ! Last point with N > 0 at the height Hmax (N == 0 at the satellite)
 call CrossLineEllips(WGS84Param, Hmax, ray,         &
                      CrossPt, RayLambda)
 !call Cart2Ellips(CrossPt(1),  CrossPt(2),  CrossPt(3),         &
 !                 point(1), point(2), point(3), WGS84Param)
 STD%SlantSteps(STD%Ntot-1) = RayLambda ! last point on slant
 STD%SlantPointsCart(STD%Ntot-1,:) = CrossPt

 ! Pressure at the upper grid level
 ! Use HscaleP2 for pressure profile above the model top
 Ps = P0*exp(-point(3)/HscaleP2)
 ! Pressure at Hmax
 Pe = P0*exp(-Hmax/HscaleP2)

 DeltaP = (Ps - Pe) / real(STD%Nup)

 p =  Ps
 do a=STD%Nmod+1, STD%Ntot-2
    ! Compute supporting points at equidistant pressure levels
    p = p - DeltaP
    h =  -HscaleP2*(log(p) - lnP0)
    call CrossLineEllips(WGS84Param, h, ray,         &
                         CrossPt, RayLambda)
    STD%SlantSteps(a) = RayLambda
    STD%SlantPointsCart(a,:) = CrossPt
 end do

 if (verbose >= 2) write(*,*) 'std_path_coord_sat> ',             &
                              'Coordinates of supporting points'
 ! So far, the distance to the GNSS receiver and the cartesian coordinates
 ! have been computed for all supporting ponts. Now compute the
 ! corresponding WGS84 coordinates.

 do a=1, std%Ntot
    !-----------------------------------------------------------------------
    ! SlantPointsCart - cartesian coordinates of a point on the slant path
    !                   with a distance = SlantSteps(a) to the GNSS receiver
    !-----------------------------------------------------------------------
    call Cart2Ellips( STD%SlantPointsCart(a,1),                &
                      STD%SlantPointsCart(a,2),                &
                      STD%SlantPointsCart(a,3),                &
                      wgspoint(1), wgspoint(2), wgspoint(3),   &
                      WGS84Param                             )
    !-----------------------------------------------------
    ! wgspoint - ellipsoidal coordinates of the same point
    !-----------------------------------------------------
    std% SlantPointsEll(a,:) = wgspoint
    !write(*,*) std%SlantPointsEll(a,1)*rad2deg,   &
    !           std%SlantPointsEll(a,2)*rad2deg, std%SlantPointsEll(a,3)

 end do

 if (verbose >= 1) write(*,*) 'std_path_coord_sat>  ... end'

end subroutine std_path_coord_sat_old


 !---------------------------------------------------------------------
 ! subroutine ztd_path_coord
 !---------------------------------------------------------------------
 !
 !> @brief Defines one single supporting point with the
 !>        GNSS station coordinates
 !>
 !> <b> call ztd_path_coord(STD, err) </b>
 !>
 !> In case of ZTDs no predefined supporting points are required. The 3D-Var
 !> interface needs only the GNSS station coordinates and provides the model
 !> levels at that position. The ZTD operator will use any number of levels
 !> provided by  3D-Var interface.
 !>
 !>
 !> @param[in,out] STD       derived type which contains all information
 !>                          about one STD
 !> @param[out]    err       return code, err = 0 - no errors
 !
 !---------------------------------------------------------------------
 subroutine ztd_path_coord (STD, err)

 type (SlantData)  ,intent(inout) :: STD
 integer           ,intent(out)   :: err

 !-------------------------------------------------------------------
 ! Define just one supporting point at the GNSS station coordinates
 ! => will be re-defined in ztd_delay when the model levels are known
 !-------------------------------------------------------------------
 allocate ( std% SlantSteps(1) )
 allocate ( std% SlantPointsCart(1,1:3) )
 allocate ( std% SlantPointsEll(1,1:3) )

 STD%SlantSteps(1)        = 0.0_wp  ! distance to receiver in m
 STD%SlantPointsCart(1,:) = STD%Station%CoordCart
 STD%SlantPointsEll(1,:)  = STD%Station%CoordEll

 end subroutine ztd_path_coord


 !---------------------------------------------------------------------
 ! subroutine destruct_SlantData
 !---------------------------------------------------------------------
 !
 !> @brief Deallocates all internal arrays or pointers of the SlantData
 !>        derved type.
 !>
 !> <b> call destruct_SlantData(slant) </b>
 !>
 !> All data computed during the processing of one STD are stored in the
 !> SlantData derived type and several arrays are allocated. These arrays
 !> can be dealloceted in order to free memory.
 !>
 !> @param[in,out] STD  derived type which contains all information about one STD
 !
 !---------------------------------------------------------------------
 elemental subroutine destruct_SlantData (slant)

   type (SlantData) ,intent(inout) :: slant

   !----------------------------------------
   ! deallocate components of type SlantData
   !----------------------------------------
   if (associated (slant% Station))         nullify    (slant% Station)
   if (associated (slant% SlantSteps))      deallocate (slant% SlantSteps)
   if (associated (slant% SlantPointsCart)) deallocate (slant% SlantPointsCart)
   if (associated (slant% SlantPointsEll))  deallocate (slant% SlantPointsEll)
   if (associated (slant% SlantPointsGeo))  deallocate (slant% SlantPointsGeo)
   if (associated (slant% SlantPointsIdx))  deallocate (slant% SlantPointsIdx)
   if (associated (slant% SlantPath))       deallocate (slant% SlantPath)
   if (associated (slant% SlantRefrac))     deallocate (slant% SlantRefrac)
   if (associated (slant% LineRefrac))      deallocate (slant% LineRefrac)
   if (associated (slant% idxi))            deallocate (slant% idxi)
   if (associated (slant% wi))              deallocate (slant% wi)

 end subroutine destruct_SlantData


 !---------------------------------------------------------------------
 ! subroutine std_delay
 !---------------------------------------------------------------------
 !
 !> @brief Computes the slant total delay along the signal path.
 !>
 !> <b> call std_delay (slant, col, delay, ladj, col_ad, err) </b>
 !>
 !> Driver routine for std_delay_ray3d and std_delay_line which should
 !> be called from the model interface. Direct calls to std_delay_ray3d
 !> and std_delay_line are not recommended.
 !>
 !> Several checks are carried out within this routine: \n
 !> - Is the number of model columns equal to number of supporting points
 !>   on the signal path? \n
 !> - Are the supporting points within the grid cells defined by the
 !>   surrounding model columns, i.e. have the correct columns been
 !>   selected? \n
 !> - Is the station close to the model surface? \n
 !> - Are some of the supporting points on the signal path below the
 !>   model surface?  ??? Supporting points are always below the real
 !>                       signal path, doeas the raytracer converge ??? \n
 !> - Are the delays within a meanigful range,  e.g.  0 < ZTD < 3 m ? \n
 !>
 !>
 !> An error code is returned if one or more checks failed:
 !>
 !> <table>
 !> <tr><th>error  <th>description
 !> <tr><td><b> 0 </b>  <td>
 !>     no error
 !> <tr><td><b> 1 </b>  <td>
 !>     Wrong number of grid columns in variable <b>col</b>, doesn't
 !>     match the number of supporting points inside the model.
 !> <tr><td><b> 2 </b>  <td>
 !>     Some supporting points are outside the grid cells defind by the model
 !>     columns, i.e. the wrong model columns have been selected.
 !> <tr><td><b> 3 </b>  <td>
 !>     The station is located below the model surface, i.e. the station
 !>     is more than "StatBelowSurface" m below the model surface.
 !> <tr><td><b> 4 </b>  <td>
 !>     The station is located far above the model surface, i.e.
 !>     more than "StatAboveSurface" m above the model surface.
 !> <tr><td><b> 5 </b>  <td>
 !>     Some supporting points on the connecting line between the station
 !>     and the GNSS satellite are below the model surface. (Slants only)
 !> <tr><td><b> 6 </b>  <td>
 !>     The delay computed by the STD operator is meaningless and should
 !>     not be used, e.g.  not within 0 < ZTD < 3 m.
 !> </table>
 !>
 !> @param[in,out] slant  derived type which contains all information
 !>                       about one STD
 !> @param[in]     col    model columns: coordinates and model fields
 !>                       on the required grid nodes
 !> @param[in]     ladj   adjoint flag: Adjoint operator is called only if
 !>                       ladj = .true.
 !> @param[out]    delay  computed STD in meter \n
 !>                       The  same value is returned by the slant
 !>                       derived type: slant%STD = delay
 !> @param[out]    col_ad derivatives computed by the adjoint operator
 !>                       (only if ladj = .true.)
 !> @param[out]    err    return code, no errors if err = 0 \n
 !>                       see error code table
 !>
 !> @todo Height difference between lowest full level and model surface?
 !>       Currently 10 m hardcoded.
 !>
 !> @todo Station height check: Check only once per station and store result
 !>                             for following observations of the same station
 !
 !---------------------------------------------------------------------
 subroutine std_delay (slant, col, delay, ladj, col_ad, err)

 type (SlantData) ,intent(inout) :: slant        ! slant meta data
 type (p_column)  ,intent(in)    :: col (:,:)    ! model columns     (input)
 real(wp)         ,intent(out)   :: delay        ! Slant Total Delay (output)
 logical          ,intent(in)    :: ladj         ! flag to calculate adjoint
 type (p_column)  ,pointer       :: col_ad (:,:) ! adjoint
 integer          ,intent(out)   :: err          ! return code, error massege

 integer   :: err2        ! error codes from subroutine calls
 real(wp)  :: maxelev
 integer   :: i, j, k, Nlow
 real (wp) :: MaxLon, MinLon, MaxLat, MinLat, Pt(1:2)
 real (wp) :: ModelSurface, ModelMin !, delay2
 real (wp), dimension(:), allocatable   :: ModelColMins
 real (wp), dimension(:,:), allocatable :: node
 !integer, save :: Nfile = 20
 real(wp)      :: gmfh         ! GMF: hydrostatic mapping function
 real(wp)      :: gmfw         ! GMF: wet mapping function

 logical :: l1, l2

 logical       :: SkipOp       ! Skip operator, don't compute first guess
 !---------------------------------------------------------------------

 ! The input will be checked before the operator is called.
 ! If some errors or problems are found, the operator
 ! can be skipped:
 ! Skip operator in case of erroes: SkipOp = .true.
 ! Call the operator anyway:        SkipOp = .false. (default)
! SkipOp = .true.
 SkipOp = .false.

 if (verbose >= 3) then
    write(*,*) 'std_delay> Station ', slant% Station%SName, char(10),  &
            'Station ID = ', slant% Station%ID, char(10),                 &
            'Geoid/Ellipsoid lon = ', slant% Station%CoordGeo(1)*rad2deg, &
                       slant% Station%CoordEll(1)*rad2deg, char(10),      &
            'Geoid/Ellipsoid lat = ', slant% Station%CoordGeo(2)*rad2deg, &
                       slant% Station%CoordEll(2)*rad2deg, char(10),      &
            'Geoid/Ellipsoid hei = ', slant% Station%CoordGeo(3) ,        &
                       slant% Station%CoordEll(3), char(10),              &
            'Cartesian X/Y/Z m   = ',  slant% Station%CoordCart
 end if

 err = 0                 ! return code, error massege: err = 0 - no errors
 delay         = invalid  ! init. delay
 slant%STD     = invalid  ! init. delay
 slant%STDline = invalid  ! init. delay

 ! Check if the expected number of columns is available
 if (ztd_col <= 0 .and. slant%Nmod .ne. ubound(col,2)) then
    err = 1
    if (verbose > 0)                                                &
    write(*,*) 'std_delay> Wrong number of grid columns: point = ', &
         slant%Nmod, ' cols = ', ubound(col,2), ', err = ', err
 end if

 ! Check if supporting point is really surrounded by grid columns
 if (ztd_col == 0) then
    ! For ztd_col != 0 only one single column is provided and no check required
    do j=1, ubound(col,2)
       MaxLon = maxval(col(:,j)%dlon)
       MinLon = minval(col(:,j)%dlon)
       MaxLat = maxval(col(:,j)%dlat)
       MinLat = minval(col(:,j)%dlat)
       Pt = slant%SlantPointsEll(j,1:2)*rad2deg
       if (Pt(1) > 180.0_wp) Pt(1) = Pt(1) - 360.0_wp
       if ( Pt(1) >= MinLon .and.     &
            Pt(1) <= MaxLon .and.     &
            Pt(2) >= MinLat .and.     &
            Pt(2) <= MaxLat       ) then
          continue
       else
          err = 2
          if (verbose > 0)                                                  &
             write(*,*) 'std_delay> Supporting point outside grid columns', &
                   Pt,  MinLon, MaxLon,  MinLat, MaxLat, ', err = ', err
       end if
    end do
 end if

 ! Check if station height is close to the interpolated model surface
 if (ztd_col <= 0) then
    ! Interpolate model surface from surrounding columns
    allocate( node(slant%Nnghb,3) )
    k = ubound(col(1,1)%gpm, 1)   ! max. level index = lowest grid level
    do i=1, slant%Nnghb
       node(i,1) = col(i,1)%dlon    ! grid node longitude in degrees
       node(i,2) = col(i,1)%dlat    ! grid node latitude in degrees
       node(i,3) = col(i,1)%gpm(k)  ! grid node height above see level (geoid)
                                    ! lowest full level, not HSURF !!!
       !write(*,*) 'std_delay> i, node ', i, node(i,:)
    end do
    ! Interpolated model surface at the station coordinates
    ! Slant coordinares are given in radian: slant%Station%CoordGeo(1:2) = lon/lat
    ModelSurface = Shepard2D(slant%Station%CoordGeo(1:2)*rad2deg, node)
    ModelSurface = ModelSurface - 10.0_wp  ! - 10 m : model surface 10 m below
                                        !          lowest full level
    deallocate( node )

 else
    ! Only one column available, no interpolation required
    k = ubound(col(1,1)%gpm, 1)    ! max. level index = lowest grid level
    ModelSurface = col(1,1)%gpm(k) ! grid node height above see level (geoid)
    ModelSurface = ModelSurface - 10.0_wp
 end if

 if ((slant%Station%CoordGeo(3)-ModelSurface) < -StatBelowSurface) then
    ! station below model surface
    err = 3
    if (verbose > 0)                                              &
    write(*,*) 'std_delay> Station ', trim(slant% Station%SName), &
               ' below model surface : ',         &
         slant%Station%CoordGeo(3)-ModelSurface, ' m, err = ', err
 else if ((slant%Station%CoordGeo(3)-ModelSurface) > StatAboveSurface) then
    ! station far above model surface
    err = 4
    if (verbose > 0)                                              &
    write(*,*) 'std_delay> Station ', trim(slant% Station%SName), &
               ' above model surface : ',         &
         slant%Station%CoordGeo(3)-ModelSurface, ' m, err = ', err
 end if

 ! Checking the interpolated model surface at the station is not allways
 ! sufficient. In mountainous regions one model column might be far above
 ! the others (> 1000 m) and vertical extrapolation over a large height would
 ! be required even if the station is above the interpolated model surface.
 ! Therefore it is checked if the station is not too far below the maximum
 ! of the near by columns.
 if (ztd_col <= 0) then
    allocate( ModelColMins(slant%Nnghb) )
    do i=1, slant%Nnghb
       ModelColMins(i) = col(i,1)%gpm(k)
    end do
    ModelMin = maxval(ModelColMins)
    ModelMin = ModelMin - 10.0_wp  ! - 10 m : model surface 10 m below
                                   !          lowest full level
    deallocate( ModelColMins )
    if ((slant%Station%CoordGeo(3)-ModelMin) < -StatBelowColumn) then
       ! station below at least one model column => vertical extrapolation
       !                                            required
       err = 3
       if (verbose > 0)                                              &
            write(*,*) 'std_delay> Station ', trim(slant% Station%SName), &
            ' below one or more model columns : ',             &
            slant%Station%CoordGeo(3)-ModelMin, ' m, err = ', err
    end if
 end if

 !write(*,*) 'std_delay> Hstation, Hmodel, err = ', &
 !     slant%Station%CoordGeo(3), ModelSurface, err


 ! Check if supporting points along the signal path are below model surface
 ! Slants only, for ZTDs the station check is sufficient.
 ! WARNING:
 ! Ellipsoidal coordinates of the slant points (slant%SlantPointsEll) are
 ! compared against geoid coordinates of the model columns. This leads to
 ! an error in the order of the geoid undulation and slant%SlantPointsGeo
 ! should be used instead. However, SlantPointsGeo are not needed anywhere
 ! else and this is just a rough estimate, so transforming all supporting
 ! points to geoid coordinates seems to be a waste of resources.
 if (slant%Obs%elevation < pi05) then
    ! STD elevation < 90 deg.
    allocate( node(slant%Nnghb,3) )
    k = ubound(col(1,1)%gpm, 1)   ! max. level index = lowest grid level
    Nlow = 0
    do j=1, slant%Nmod
       do i=1, slant%Nnghb
          node(i,1) = col(i,j)%dlon    ! grid node longitude in degrees
          node(i,2) = col(i,j)%dlat    ! grid node latitude in degrees
          node(i,3) = col(i,j)%gpm(k)  ! grid node height ASL (geoid)
          !write(*,*) 'std_delay> i, node ', i, node(i,:)
       end do
       ! Interpolation of column surface
       ModelSurface = Shepard2D(slant%SlantPointsEll(j,1:2)*rad2deg, node)
       ModelSurface = ModelSurface - 10.0_wp  ! - 10 m: model surface 10 m
                                              !         below lowest full level
       if ((slant%Station%CoordEll(3)-ModelSurface) < -StatBelowSurface) then
          Nlow = Nlow + 1
       end if
    end do
    if (Nlow > 0) then
       err = 5
       if (verbose > 0)                                              &
       write(*,*) 'std_delay> Station ', trim(slant% Station%SName), &
                  ' Supporting points below model surface, ', &
            'Nlow = ', Nlow, ', err = ', err
    end if
 end if

 if (.false.) then
 !if (.true.) then
    ! write files with "true" and "estimated" profiles
    call PrintProfiles (slant, col, -1.0_wp, err)
 end if

 !---------------------------------------------------------------------
 ! The observation will be rejected if an error occured. However, these
 ! observations are not removed from the feedback files and will appear
 ! in the statistics. Usually this is the desired behaviour, otherwise
 ! the observations need to be DISMISSED.
 ! => The operator needs to run for these observations, even if no
 !    reliable result can be expected.
 !---------------------------------------------------------------------

 if (.not. (SkipOp .and. err > 0) ) then
 if (ztd_col <= 0) then
    ! surrounding model columns for all supporting points available
    if (UseRaytracer >= 0) then
       ! call raytracer and integrate along the curved signal path
       maxelev = real(UseRaytracer,wp) * deg2rad
       call std_delay_ray3d (slant, maxelev, col, delay, ladj, col_ad, err2)
    else
       ! for backward compatibility only: call old code with adjoint
       call std_delay_line (slant, col, delay, ladj, col_ad, err2)
       !      call ztd_delay (slant, col(1:1,1:1), delay2, ladj, col_ad, err)
       !      write(*,*) 'std_delay> old ZTD, new ZTD, Delta : ', &
       !           slant%Obs%slant, delay, delay2, delay-delay2,  &
       !           delay-slant%Obs%slant, delay2-slant%Obs%slant, slant%Nnghb
       !write(*,*) 'std_delay> ERROR std_delay_line'
    end if
 else
    ! Only one grid column for ZTD available: Call optimized ZTD operator
    call ztd_delay (slant, col, delay, ladj, col_ad, err2)
 end if
 end if

 ! Check if slant has a meaningful value, e.g.  0 < ZTD < 3 m
 if (slant%Obs%elevation < pi05) then
    ! STD, compute mapping function for mapping the maximum STD
    call gmf(slant%Obs%time,       & ! MJD
             slant% Station%CoordGeo(2),    & ! lat, rad
             slant% Station%CoordGeo(1),    & ! lon, rad
             slant% Station%CoordGeo(3),    & ! height, m
             pi05 - slant%Obs%elevation,    & ! zenith distance, rad
             gmfh, gmfw)                      ! => dry, wet mapping function
 else
    gmfh = 1.0_wp     ! "mapping" of ZTD
 end if

 !write(*,*) 'std_delay>  slant%STD, gmfh, gmfw ' , slant%STD, gmfh, gmfw

 if (delay < ZTDminUse .or. delay > ZTDmaxUse*gmfh) then
    err = 6
    if (verbose > 0)                                                     &
         write(*,*) 'std_delay> Strange delay found, ',                  &
         ' station = ', trim(slant%Station%SName),            &
         ' STD = ', delay,                                    &
         ' at an elevation of ', slant%Obs%elevation*rad2deg, &
         ' degrees, err = ', err
 end if

 if (err2 /= 0) then
    ! Error during iteration, slant path could not be estimated, no delay
    err = 6
    delay = 0.0_wp
    if (verbose > 0)   &
         write(*,*) 'std_delay> Error during iteration, reject slant'
 end if

 ! write detailled diagnostics
 !if (err /= 0 .and. verbose >= 4)  then
!!$ if (.true.) then
!!$    call PrintProfiles (slant, col, delay, err)
!!$ end if

 !call model_abort (-1,-1, 'Stop after first slant', 'std_delay')

 end subroutine std_delay


 subroutine PrintProfiles (slant, col, delay, err)

 type (SlantData) ,intent(in) :: slant        ! slant meta data
 type (p_column)  ,intent(in) :: col (:,:)    ! model columns     (input)
 real (wp),intent(in)         :: delay
 integer,intent(in)           :: err

 integer, save :: Nfile = 20
 integer :: k, j, i
 integer :: ounit

 ! write files with "true" and "estimated" profiles
 !-------------------------------------------------------------------
 Nfile = Nfile + 1

 write(*,*) 'Print profile to file fort ', dace% pe*100+Nfile

#ifndef __COSMO__
 ! Output unit, unit not open, write output files fort.unit, e.g. fort.1321
 ounit = dace% pe*100+Nfile

 write(ounit,*) '# Profiles'
 write(ounit,*) '#'
 write(ounit,*) '# Output written by >PrintProfiles<, module >mo_std_operator<'
 if (delay < 0.0) then
    write(ounit,*) '# >PrintProfiles< called before the STD operator'
 else
    write(ounit,*) '# >PrintProfiles< called after the STD operator'
 end if
 write(ounit,*)
 write(ounit,*) 'Delay , STD = ', slant%Obs%slant, ' m (obs)', char(10),       &
                '        STD = ', delay, ' m (mod)', char(10),                 &
               ' at an elevation of ', slant%Obs%elevation*rad2deg,' degrees', &
               char(10),                                                       &
               ', err = ', err
 write(ounit,*)
 if (err /= 0) then
    write(ounit,*) 'Error detected, observation will not be used.'
    write(ounit,*) 'err = 1 - Wrong number of grid columns'
    write(ounit,*) 'err = 2 - Supporting point outside grid columns'
    write(ounit,*) 'err = 3 - Station below model surface'
    write(ounit,*) 'err = 4 - Station  above model surface'
    write(ounit,*) 'err = 5 - Supporting points below model surface'
    write(ounit,*)
 end if

 write(ounit,*)                                                               &
      'Station ', slant% Station%SName,                              char(10),&
      'Station ID = ', slant% Station%ID,                            char(10),&
      'Geoid/Ellipsoid lon = ', slant% Station%CoordGeo(1)*rad2deg,           &
                                slant% Station%CoordEll(1)*rad2deg,  char(10),&
      'Geoid/Ellipsoid lat = ', slant% Station%CoordGeo(2)*rad2deg,           &
                                slant% Station%CoordEll(2)*rad2deg,  char(10),&
      'Geoid/Ellipsoid hei = ', slant% Station%CoordGeo(3) ,                  &
                                slant% Station%CoordEll(3), char(10),         &
      'Cartesian X/Y/Z m   = ', slant% Station%CoordCart
 write(ounit,*)

 write(ounit,*)                                               &
        'Slant info: ',                              char(10),&
        ' Station ',    slant%Station%SName,         char(10),& ! station short name
        ' ID ',         slant%Station%ID,            char(10),& ! station ID
        ' time ',       slant%Obs%time,              char(10),& ! observation time
        ' Sat ',        slant%Obs%satellite,         char(10),& ! observed satellite
        ' Elev ',       slant%Obs%elevation*rad2deg, char(10),& ! elevation
        ' Azi ',        slant%Obs%azimuth*rad2deg,   char(10),& ! azimuth
        ' ZTD_obs ',    slant%Obs%zslant,            char(10),& ! observed STD
        ' STD_obs ',    slant%Obs%slant,             char(10),& ! STD mapped to zenith
        ' STD_gerade ', slant%STDline,               char(10),& ! model STD, line
        ' STD_ray ',    slant%STD                               ! model STD, raytracer
 write(ounit,*)

 write(ounit,*) '#'
 write(ounit,*) '# Supporting points along the signal path'
 write(ounit,*) '# Last point in model: ', slant%Nmod
 write(ounit,*) '# Longitude              Latitude               Height (ellipsoid)'
 do i=1, size(slant%SlantPointsEll(:,3))
    write(ounit,*) i, slant%SlantPointsEll(i,1:2)*rad2deg, slant%SlantPointsEll(i,3)
 end do
 write(ounit,*)

 write(ounit,*) '#'
 write(ounit,*) '# Model profiles: gpm - p - t - q '
 write(ounit,*) '#'

 do i=1, slant%Nmod
     write(ounit,*) 'Supporting point ', i, ', coord = ', &
          slant%SlantPointsEll(i,1:2)*rad2deg, slant%SlantPointsEll(i,3)
     write(ounit,*) 'Surrounding model columns:'
     do j=1, size(col(:,1))
        write(ounit,*) 'Lon/lat = ', col(j,i)%dlon, col(j,i)%dlat
     end do
     write(ounit,*)
     do k=1, ubound(col(1,1)%p,1)
        if (ubound(col,2) > 2) then
           write(ounit,                                                           &
    '(i2,tr2,3(f9.2,tr2),tr5,3(f10.2,tr2),tr5,3(f6.2,tr2),tr5,3(f9.6,tr2))')   &
             k, col(1,i)%gpm(k), col(2,i)%gpm(k), col(3,i)%gpm(k), &
              col(1,i)%p(k), col(2,i)%p(k), col(3,i)%p(k), &
              col(1,i)%t(k), col(2,i)%t(k), col(3,i)%t(k), &
              col(1,i)%q(k), col(2,i)%q(k), col(3,i)%q(k)
        else
           write(ounit,       &
    '(i2,tr2,(f9.2,tr2),tr1,(f10.2,tr2),tr1,(f6.2,tr2),tr1,(f9.6,tr2))')   &
             k, col(1,i)%gpm(k), col(1,i)%p(k), col(1,i)%t(k),  &
             col(1,i)%q(k)
        end if
     end do
     write(ounit,*)
  end do

 write(ounit,*)
 write(ounit,*) '#'
 write(ounit,*) '# Refractivity profile along the signal path'
 write(ounit,*) '# No  SlantPos [m]  SlantRefrac  LineRefrac'
 if ( associated(Slant%SlantSteps)  .and.         &
      associated(Slant%SlantRefrac) .and.         &
      associated(Slant%LineRefrac)        ) then
    do j=1, Slant%Ntot
       write(dace% pe*100+Nfile,'(i4,tr2,f14.2,tr2,f17.12,tr2,f17.12)') &
            j, Slant%SlantSteps(j), Slant%SlantRefrac(j), Slant%LineRefrac(j)
    end do
 else if ( associated(Slant%SlantSteps)  .and.         &
           associated(Slant%SlantRefrac)       ) then
    do j=1, Slant%Ntot
       write(ounit,'(i4,tr2,f14.2,tr2,f17.12)') &
            j, Slant%SlantSteps(j), Slant%SlantRefrac(j)
    end do
 end if

!!$ write(dace% pe*100+Nfile,*) '#'
!!$ write(dace% pe*100+Nfile,*) '# true and estimated profiles'
!!$ write(dace% pe*100+Nfile,*) '# Hlevel = ', Hlevel
!!$ write(dace% pe*100+Nfile,*) '# Heights = ', Heights
!!$ write(dace% pe*100+Nfile,*) '# Hpoints = ', Hpoints
!!$ write(dace% pe*100+Nfile,*) '# '
!!$ write(dace% pe*100+Nfile,*) '# Plevel = ', Plevel
!!$ write(dace% pe*100+Nfile,*) '# Pheights = ', Pheights
!!$ write(dace% pe*100+Nfile,*) '# Ppoints = ', Ppoints
!!$ write(dace% pe*100+Nfile,*) '# '
!!$ write(dace% pe*100+Nfile,*) '# Profil Stuetzpunkte:'
!!$ write(dace% pe*100+Nfile,*) '# Station lon/lat/h ', &
!!$      slant%SlantPointsEll(1,1:2)*rad2deg,  slant%SlantPointsEll(1,3)
!!$ write(dace% pe*100+Nfile,*) slant%SlantPointsEll(:,3)
!!$ write(dace% pe*100+Nfile,*) '#'
!!$ write(dace% pe*100+Nfile,*) '# Modellprofil'
!!$ write(dace% pe*100+Nfile,*) '# latitude  = ', col(1,1)%dlat
!!$ write(dace% pe*100+Nfile,*) '# longitude = ', col(1,1)%dlon
!!$ write(dace% pe*100+Nfile,*) col(1,1)%gpm(:)

 CALL FLUSH(ounit)
#endif

 end subroutine PrintProfiles


 !---------------------------------------------------------------------
 ! subroutine ztd_delay
 !---------------------------------------------------------------------
 !
 !> @brief Computes the zenith total delay for one model column
 !>
 !> <b> call ztd_delay (slant, col, delay, ladj, col_ad, err) </b>
 !>
 !> This is an optimized ZTD operator which works with one single model
 !> column. The interpolation is done by the 3D-Var interface and the
 !> ZTD operator consists of three parts: \n
 !> 1) Computation of the refractivity on the model levels, \n
 !> 2) vertical extension of the profile above the model top and \n
 !> 3) vertical numerical integration which provides the ZTD.
 !>
 !> It is assumed that the GNSS station height as well as the geopotential
 !> height are heights above geoid, i.e. mean sea level, and no geoid
 !> correction and coordinate transformations are required. \n
 !> => use slant%Station%CoordGeo(3) as station height,
 !>        not CoordEll(3) as in the other delay routines \n
 !> => use slant%SlantPointsGeo to define the vertical profile,
 !>        not slant%SlantPointsEll
 !>
 !> @param[in,out] slant  derived type which contains all information
 !>                       about one STD or ZTD
 !> @param[in]     col    model column: coordinates and model fields
 !>                       on the given grid node
 !> @param[in]     ladj   adjoint flag: Adjoint operator is called only if
 !>                       ladj = .true.
 !> @param[out]    delay  computed ZTD in meter \n
 !>                       The computed delays are also written to the slant
 !>                       structure.
 !> @param[out]    col_ad derivatives computed by the adjoint operator
 !>                       (only if ladj = .true.)
 !> @param[out]    err    return code, no errors if err = 0
 !
 !---------------------------------------------------------------------
 subroutine ztd_delay (slant, col, delay, ladj, col_ad, err)

 type (SlantData) ,intent(inout) :: slant        ! slant meta data
 type (p_column)  ,intent(in)    :: col (:,:)    ! model columns     (input)
 real(wp)         ,intent(out)   :: delay        ! Slant Total Delay (output)
 logical          ,intent(in)    :: ladj         ! flag to calculate adjoint
 type (p_column)  ,pointer       :: col_ad (:,:) ! adjoint
 integer          ,intent(out)   :: err

 integer  :: i, j, k, kstat, n
 real(wp) :: Delta, Di, Hi
 real(wp) :: e, Refrac, e_ad, p_ad, t_ad !, Refrac_ad
 real(wp) :: STD_ad, RefPt_ad !, ExpInt_ad
 real(wp), dimension(:,:), allocatable :: RefracProfile, RefracProfile_ad
 real(wp), dimension(1:2,1:2) :: Anode, Anode_ad, Bnode, Bnode_ad
!real(wp) :: TotalDelay, LineDelay

 ! Test of tangent-linear and adjoint code
 logical                         :: testad = .false.  !.true.
 type (p_column)  ,pointer       :: col_tl (:,:) ! tangent_linear (input)
 type (p_column)  ,pointer       :: coltest_ad (:,:) ! adjoint (output)
 type (p_column)  ,pointer       :: Dcol (:,:) !
 real(wp)                        :: delay_tl     ! Slant Total Delay, TL
 real(wp)                        :: delay_ad     ! Slant Total Delay, adjoint
 real(wp)                        :: Ddelay, Diff, epsilon
 integer                         :: a, b
 real(wp)                        :: sumx, sumy !, delay2
 !-------------------------------------------------------------------------

 err = 0

 !----------------------------------------------------------------------------
 ! re-define supporting points according to "col"
 ! This is necessary as "ztd_path_coord" provides only the station coordinates
 !----------------------------------------------------------------------------

 ! Find column levels above the station and set the number of supporting points
 slant% Nup = NStepVertTop
 do k= ubound((col(1,1)%gpm),1), 1, -1
    if (slant%Station%CoordGeo(3) <= col(1,1)%gpm(k)) then
       ! GNSS station height below level k
       kstat = k
       ! set the number of supporting points inside the model to the number of
       ! col-levels above the GNSS station
       ! The namelist variable "NStepVertMod" is meaningless in this case and
       ! will be ignored!
       slant% Nmod = kstat + 1
       slant% Ntot = slant% Nmod + slant% Nup
       exit
    end if
 end do

 ! delete preliminary settings from "ztd_path_coord"
 if ( associated(slant%SlantSteps) ) deallocate(slant%SlantSteps)
 if ( associated(slant%SlantPointsCart) ) deallocate(slant%SlantPointsCart)
 if ( associated(slant%SlantPointsGeo) ) deallocate(slant%SlantPointsGeo)

 ! define supporting points as column levels above the station + NStepVertTop
 allocate ( slant% SlantSteps(1:slant% Ntot) )
 allocate ( slant% SlantPointsGeo(1:slant% Ntot,1:3) )

 slant%SlantSteps(1)           = 0.0_wp  ! distance to receiver in m
 ! use lon/lat (rad) of column for all supporting points
 forall (i=1:slant% Ntot)
    slant% SlantPointsGeo(i,1:2) = (/col(1,1)%dlon, col(1,1)%dlat/) * deg2rad
 end forall
 ! use GNSS station height in m as height of first point
 slant% SlantPointsGeo(1,3)    = slant%Station%CoordGeo(3)

 n = 1   ! array index of  slant%SlantSteps and SlantPointsGeo
 do k=kstat, 1, -1
    n = n + 1
    slant%SlantSteps(n)        = col(1,1)%gpm(k) - slant%Station%CoordGeo(3)
    slant% SlantPointsGeo(n,3) = col(1,1)%gpm(k)
 end do

 ! extended profile above model top
 ! heights increase by Di in each step
 Delta = col(1,1)%gpm(1) - col(1,1)%gpm(2)
 Di    = slant%Nup**(-2.0_wp) * (Hmax - col(1,1)%gpm(1) - slant%Nup*Delta)
 do i=n+1 , slant%Ntot
    j = i - n
    Hi = col(1,1)%gpm(1) + j*(Delta + j*Di)
    slant%SlantSteps(i)       = Hi - slant%Station%CoordGeo(3)
    slant%SlantPointsGeo(i,3) = Hi
 end do

 !------------------------------------------------------------------
 ! Estimate the zenith total delay (ZTD):
 ! Compute the refractivity N on the supporting points and integrate
 !------------------------------------------------------------------

 allocate( RefracProfile(1:slant%Ntot,1:2) )

 RefracProfile(:,1) = slant%SlantSteps
 RefracProfile(:,2) = 0.0_wp

 if (kstat < ubound((col(1,1)%gpm),1)) then
    ! station between two model levels:
    ! Compute refractivity also at level below station
    kstat = kstat + 1
    n = 0
 else
    ! station below lowest model level:
    ! extrapolate refractivity later
    n = 1
    RefracProfile(1,2) = invalid
 end if

 do k=kstat, 1, -1
    n = n + 1
    ! Partial pressure of water vapour e
    e = (col(1,1)%p(k) * col(1,1)%q(k)) / (RDRD +   &
               col(1,1)%q(k)*EMRDRD)
    ! Refractivity N
    Refrac = NWein(col(1,1)%p(k), col(1,1)%t(k), e)
    RefracProfile(n,2) = Refrac
 end do

 Anode(1,1) = col(1,1)%gpm(kstat)    ! height lower level
 Anode(2,1) = col(1,1)%gpm(kstat-1)  ! height upper level
 if (RefracProfile(1,2) > 0.0_wp) then
    ! Interpolate refractivity at GNSS station
    Anode(1,2) = RefracProfile(1,2)     ! refractivity lower level
    Anode(2,2) = RefracProfile(2,2)     ! refractivity upper level
 else
    ! Extrapolate refractivity below model to GNSS station
    Anode(1,2) = RefracProfile(2,2)     ! refractivity lower level
    Anode(2,2) = RefracProfile(3,2)     ! refractivity upper level
 end if
 RefracProfile(1,2) = ExpInt1D(slant%SlantPointsGeo(1,3), Anode)

 ! extrapolate refractivity above model top
 Bnode(1,1) = col(1,1)%gpm(2)               ! height lower level
 Bnode(2,1) = col(1,1)%gpm(1)               ! height upper level
 Bnode(1,2) = RefracProfile(Slant%Nmod-1,2) ! refractivity lower level
 Bnode(2,2) = RefracProfile(Slant%Nmod,2)   ! refractivity upper level
 do i=Slant%Nmod+1, Slant%Ntot
    RefracProfile(i,2) = ExpInt1D(slant%SlantPointsGeo(i,3), Bnode)
 end do

 ! Integrate refractivity from GNSS station to Hmax
 delay = 1.0e-6_wp * IntegPolyCube(RefracProfile)

 ! for backward compatibility ....
 slant%STD         = delay
 Slant%Obs%mSTD    = delay
 Slant%Obs%STDline = delay

 Slant%Obs%Press   = invalid

 if (pl_method > 0) then
    ! compute pressure at the GNSS station => plevel for LETKF
    ! find reference point on slant at a given height above the station
    ! => height above station: Href m
    do k=kstat, 1, -1
       if (col(1,1)%gpm(k)-col(1,1)%gpm(kstat) >= Href) then
          Slant%Obs%Press = col(1,1)%p(k)
          exit
       end if
    end do
 end if

 !========================================
 ! adjoint code follows (derives Jakobian)
 !========================================

 if (ladj) then
 !if (.true.) then

    allocate( RefracProfile_ad(1:slant%Ntot,1:2) )

    ! Integrate refractivity from GNSS station to Hmax
    delay_ad = 1._wp
    STD_ad = delay_ad * 1.0E-6_wp
    call IntegPolyCube_ad(RefracProfile, RefracProfile_ad, STD_ad)

    ! extrapolate refractivity above model top
    ! Bnode has not been reused and should be in the same state as above.
    Bnode_ad = 0.0_wp
    do i=Slant%Ntot, Slant%Nmod+1, -1
       call ExpInt1D_ad(slant%SlantPointsGeo(i,3), Bnode,           &
                        RefracProfile_ad(i,2), RefPt_ad, Bnode_ad)
    end do
    RefracProfile_ad(Slant%Nmod,2)   = RefracProfile_ad(Slant%Nmod,2) &
                                       + Bnode_ad(2,2)
    RefracProfile_ad(Slant%Nmod-1,2) = RefracProfile_ad(Slant%Nmod-1,2) &
                                       + Bnode_ad(1,2)

    ! Interpolate refractivity at GNSS station
    ! Anode has not been reused and should be in the same state as above.
    Anode_ad = 0.0_wp
    call ExpInt1D_ad(slant%SlantPointsGeo(1,3), Anode, &
                     RefracProfile_ad(1,2), RefPt_ad, Anode_ad)

    ! Reset kstat to original state
    kstat = slant%Nmod - 1
    if (kstat < ubound((col(1,1)%gpm),1)) then
       ! station between two model levels:
       ! Compute refractivity also at level below station
       kstat = kstat + 1
       n = 0
    else
       ! station below lowest model level:
       ! extrapolate refractivity later
       n = 1
       RefracProfile(1,2) = invalid
    end if

    if (RefracProfile(1,2) > 0.0_wp) then
       ! Interpolate refractivity at GNSS station
       RefracProfile_ad(1,2) = RefracProfile_ad(1,2) + Anode_ad(1,2)
       RefracProfile_ad(2,2) = RefracProfile_ad(2,2) + Anode_ad(2,2)
    else
       ! Extrapolate refractivity below model to GNSS station
       RefracProfile_ad(2,2) = RefracProfile_ad(2,2) + Anode_ad(1,2)
       RefracProfile_ad(3,2) = RefracProfile_ad(3,2) + Anode_ad(2,2)
    end if

    do k=kstat, 1, -1
       n = n + 1
       ! Re-compute partial pressure of water vapour e
       e = (col(1,1)%p(k) * col(1,1)%q(k)) / (RDRD +   &
               col(1,1)%q(k)*EMRDRD)

       ! Refractivity N
       call NWein_ad(col(1,1)%p(k), col(1,1)%t(k), e,           &
                     p_ad, t_ad, e_ad,  RefracProfile_ad(n,2))
       col_ad(1,1)%p(k) = col_ad(1,1)%p(k) + p_ad
       col_ad(1,1)%t(k) = col_ad(1,1)%t(k) + t_ad

       ! Partial pressure of water vapour e
       col_ad(1,1)%p(k) = col_ad(1,1)%p(k) +                     &
                          ( col(1,1)%q(k) / (RDRD +              &
                            col(1,1)%q(k)*EMRDRD) ) * e_ad
       col_ad(1,1)%q(k) = col_ad(1,1)%q(k) +                        &
                          ((( RDRD + col(1,1)%q(k) * EMRDRD ) *     &
                           col(1,1)%p(k) -                          &
                          col(1,1)%p(k) * col(1,1)%q(k) *  EMRDRD)  &
                           / (RDRD + col(1,1)%q(k)*EMRDRD)**2) *    &
                          e_ad
    end do
    !--------- end of adjoint code ----------------------------------------

    if (testad) then
       ! Test tangent-linear and adjoint code

       ! allocate  col_tl, coltest_ad and Dcol
       allocate( col_tl(1:1,1:1) )
       allocate( coltest_ad(1:1,1:1) )
       allocate( Dcol(1:1,1:1) )

       a = lbound(col(1,1)%p,1)
       b = ubound(col(1,1)%p,1)
       ! col_tl
       allocate( col_tl(1,1)%p(a:b) )
       allocate( col_tl(1,1)%q(a:b) )
       allocate( col_tl(1,1)%t(a:b) )
       col_tl(1,1)%p = 0.0_wp
       col_tl(1,1)%q = 0.0_wp
       col_tl(1,1)%t = 0.0_wp
       ! coltest_ad
       allocate( coltest_ad(1,1)%p(a:b) )
       allocate( coltest_ad(1,1)%q(a:b) )
       allocate( coltest_ad(1,1)%t(a:b) )
       coltest_ad(1,1)%p = 0.0_wp
       coltest_ad(1,1)%q = 0.0_wp
       coltest_ad(1,1)%t = 0.0_wp
       ! Dcol
       allocate( Dcol(1,1)%p(a:b) )
       allocate( Dcol(1,1)%q(a:b) )
       allocate( Dcol(1,1)%t(a:b) )
       Dcol(1,1)%gpm => col(1,1)%gpm
       Dcol(1,1)%geoid = col(1,1)%geoid
       Dcol(1,1)%dlat = col(1,1)%dlat
       Dcol(1,1)%dlon = col(1,1)%dlon

       ! Epsilon for finite differences
       epsilon = 5.0E-6_wp

       ! TL inkrement
       col_tl(1,1)%p(:) = 500.0_wp
       col_tl(1,1)%q(:) = 0.0010_wp
       col_tl(1,1)%t(:) = 0.10_wp

       ! Inkrement
       Dcol(1,1)%p(:) = col(1,1)%p(:) +  epsilon*col_tl(1,1)%p(:)
       Dcol(1,1)%q(:) = col(1,1)%q(:) +  epsilon*col_tl(1,1)%q(:)
       Dcol(1,1)%t(:) = col(1,1)%t(:) +  epsilon*col_tl(1,1)%t(:)

       ! Test of tangent-linear code using finite differences
       call  ztd_delay_tl (slant, Dcol, Ddelay, col_tl, delay_tl, err)
       ! The 2. call is necessary to obtain the correct delay_tl !!!!!!!
       call  ztd_delay_tl (slant, col, delay, col_tl, delay_tl, err)

       Diff = (Ddelay - delay) / epsilon

       write(*,*) 'Test tangent-linear code <=> finite differences', char(10), &
                  'delay, Ddelay : ', delay,  Ddelay,                char(10), &
                  'Diff, delay_tl, Delta :', Diff, delay_tl, Diff-delay_tl

       ! Test of adjoint code using array product test
       delay_ad = delay_tl
       call ztd_delay_ad (slant, col, delay, delay_ad, coltest_ad, err)

       sumy =  delay_tl**2
       sumx = 0.0_wp
       sumx = sumx + dot_product(coltest_ad(1,1)%p, col_tl(1,1)%p)
       sumx = sumx + dot_product(coltest_ad(1,1)%q, col_tl(1,1)%q)
       sumx = sumx + dot_product(coltest_ad(1,1)%t, col_tl(1,1)%t)

       write(*,*) 'Test adjoint code',  char(10), &
              'Norm X, Norm Y,  Norm X / Norm Y = ', sumx, sumy, sumx/sumy

    end if  ! if (testad) then

 end if   ! if (ladj) then

 end subroutine ztd_delay


 !---------------------------------------------------------------------
 ! subroutine ztd_delay_tl
 !---------------------------------------------------------------------
 !
 !> @brief Tangent-linear routine for ztd_delay
 !>
 !> <b> call ztd_delay_tl (slant, col, delay, col_tl, delay_tl, err) </b>
 !>
 !> The tangent-linear code contains the original code from ztd_delay and
 !> the derivatives.
 !>
 !> This is an optimized ZTD operator which works with one single model
 !> column. The interpolation is done by the 3D-Var interface and the
 !> ZTD operator consists of three parts: \n
 !> 1) Computation of the refractivity on the model levels, \n
 !> 2) vertical extension of the profile above the model top and \n
 !> 3) vertical numerical integration which provides the ZTD.
 !>
 !> It is assumed that the GNSS station height as well as the geopotential
 !> height are heights above geoid, i.e. mean sea level, and no geoid
 !> correction and coordinate transformations are required. \n
 !> => use slant%Station%CoordGeo(3) as station height,
 !>        not CoordEll(3) as in the other delay routines \n
 !> => use slant%SlantPointsGeo to define the vertical profile,
 !>        not slant%SlantPointsEll
 !>
 !> @param[in,out] slant  derived type which contains all information
 !>                       about one STD or ZTD
 !> @param[in]     col    model column: coordinates and model fields
 !>                       on the given grid node
 !> @param[out]    delay  computed ZTD in meter \n
 !>                       The computed delays are also written to the slant
 !>                       structure.
 !> @param[in]     col_tl    variation of the input parameters
 !> @param[out]    delay_tl  derivative of delay
 !> @param[out]    err       return code, no errors if err = 0
 !
 !---------------------------------------------------------------------
 subroutine ztd_delay_tl (slant, col, delay, col_tl, delay_tl, err)

 type (SlantData) ,intent(inout) :: slant        ! slant meta data
 type (p_column)  ,intent(in)    :: col (:,:)    ! model columns     (input)
 real(wp)         ,intent(out)   :: delay        ! Slant Total Delay (output)
 type (p_column)  ,pointer       :: col_tl (:,:) ! tangent_linear (input)
 real(wp)         ,intent(out)   :: delay_tl     ! Slant Total Delay, TL (out)
 integer          ,intent(out)   :: err          ! error code

 integer  :: i, j, k, kstat, n
 real(wp) :: Delta, Di, Hi
 real(wp) :: e, Refrac, e_tl, Refrac_tl
 real(wp), dimension(:,:), allocatable     :: RefracProfile, RefracProfile_tl
 real(wp), dimension(1:2,1:2) :: node, node_tl
!real(wp) :: TotalDelay, LineDelay
 !-------------------------------------------------------------------------

 err = 0

 !----------------------------------------------------------------------------
 ! re-define supporting points according to "col"
 ! This is necessary as "ztd_path_coord" provides only the station coordinates
 !----------------------------------------------------------------------------

 ! Find column levels above the station and set the number of supporting points
 slant% Nup = NStepVertTop
 do k= ubound((col(1,1)%gpm),1), 1, -1
    if (slant%Station%CoordGeo(3) <= col(1,1)%gpm(k)) then
       ! GNSS station height below level k
       kstat = k
       ! set the number of supporting points inside the model to the number of
       ! col-levels above the GNSS station
       ! The namelist variable "NStepVertMod" is meaningless in this case and
       ! will be ignored!
       slant% Nmod = kstat + 1
       slant% Ntot = slant% Nmod + slant% Nup
       exit
    end if
 end do

 ! delete preliminary settings from "ztd_path_coord"
 if ( associated(slant%SlantSteps) ) deallocate(slant%SlantSteps)
 if ( associated(slant%SlantPointsCart) ) deallocate(slant%SlantPointsCart)
 if ( associated(slant%SlantPointsGeo) ) deallocate(slant%SlantPointsGeo)

 ! define supporting points as column levels above the station + NStepVertTop
 allocate ( slant% SlantSteps(1:slant% Ntot) )
 allocate ( slant% SlantPointsGeo(1:slant% Ntot,1:3) )

 slant%SlantSteps(1)           = 0.0_wp  ! distance to receiver in m
 ! use lon/lat (rad) of column for all supporting points
 forall (i=1:slant% Ntot)
    slant% SlantPointsGeo(i,1:2) = (/col(1,1)%dlon, col(1,1)%dlat/) * deg2rad
 end forall
 ! use GNSS station height in m as height of first point
 slant% SlantPointsGeo(1,3)    = slant%Station%CoordGeo(3)

 n = 1   ! array index of  slant%SlantSteps and SlantPointsGeo
 do k=kstat, 1, -1
    n = n + 1
    slant%SlantSteps(n)        = col(1,1)%gpm(k) - slant%Station%CoordGeo(3)
    slant% SlantPointsGeo(n,3) = col(1,1)%gpm(k)
 end do

 ! extended profile above model top
 ! heights increase by Di in each step
 Delta = col(1,1)%gpm(1) - col(1,1)%gpm(2)
 Di    = slant%Nup**(-2.0_wp) * (Hmax - col(1,1)%gpm(1) - slant%Nup*Delta)
 do i=n+1 , slant%Ntot
    j = i - n
    Hi = col(1,1)%gpm(1) + j*(Delta + j*Di)
    slant%SlantSteps(i)       = Hi - slant%Station%CoordGeo(3)
    slant%SlantPointsGeo(i,3) = Hi
 end do

 !------------------------------------------------------------------
 ! Estimate the zenith total delay (ZTD):
 ! Compute the refractivity N on the supporting points and integrate
 !------------------------------------------------------------------

 allocate( RefracProfile   (1:slant%Ntot,1:2) )
 allocate( RefracProfile_tl(1:slant%Ntot,1:2) )

 RefracProfile(:,1) = slant%SlantSteps
 RefracProfile(:,2) = 0.0_wp
 RefracProfile_tl = 0.0_wp

 if (kstat < ubound((col(1,1)%gpm),1)) then
    ! station between two model levels:
    ! Compute refractivity also at level below station
    kstat = kstat + 1
    n = 0
 else
    ! station below lowest model level:
    ! extrapolate refractivity later
    n = 1
    RefracProfile(1,2) = invalid
 end if

 do k=kstat, 1, -1
    n = n + 1
    ! Partial pressure of water vapour e
    e = (col(1,1)%p(k) * col(1,1)%q(k)) / (RDRD +   &
               col(1,1)%q(k)*EMRDRD)
    e_tl = ((( RDRD + col(1,1)%q(k) * EMRDRD ) *       &
             col(1,1)%p(k) -                           &
             col(1,1)%p(k) * col(1,1)%q(k) *  EMRDRD)  &
             / (RDRD + col(1,1)%q(k)*EMRDRD)**2)       &
             * col_tl(1,1)%q(k) +                      &
               ((col(1,1)%q(k)) /                      &
                (RDRD + col(1,1)%q(k)*EMRDRD))         &
             * col_tl(1,1)%p(k)

    ! Refractivity N
    Refrac    = NWein(col(1,1)%p(k), col(1,1)%t(k), e)
    Refrac_tl = NWein_tl( col(1,1)%p(k), col(1,1)%t(k), e,      &
                          col_tl(1,1)%p(k), col_tl(1,1)%t(k),   &
                          e_tl )
    RefracProfile(n,2)    = Refrac
    RefracProfile_tl(n,2) = Refrac_tl
 end do

 node(1,1) = col(1,1)%gpm(kstat)    ! height lower level
 node(2,1) = col(1,1)%gpm(kstat-1)  ! height upper level
 node_tl = 0.0_wp
 if (RefracProfile(1,2) > 0.0_wp) then
    ! Interpolate refractivity at GNSS station
    node(1,2)    = RefracProfile(1,2)        ! refractivity lower level
    node(2,2)    = RefracProfile(2,2)        ! refractivity upper level
    node_tl(1,2) = RefracProfile_tl(1,2)     ! refractivity lower level
    node_tl(2,2) = RefracProfile_tl(2,2)     ! refractivity upper level
 else
    ! Extrapolate refractivity below model to GNSS station
    node(1,2)    = RefracProfile(2,2)        ! refractivity lower level
    node(2,2)    = RefracProfile(3,2)        ! refractivity upper level
    node_tl(1,2) = RefracProfile_tl(2,2)     ! refractivity lower level
    node_tl(2,2) = RefracProfile_tl(3,2)     ! refractivity upper level
 end if
 RefracProfile(1,2)    = ExpInt1D(slant%SlantPointsGeo(1,3), node)
 RefracProfile_tl(1,2) = ExpInt1D_tl(slant%SlantPointsGeo(1,3), node, &
                                     0.0_wp, node_tl)

 ! extrapolate refractivity above model top
 node(1,1) = col(1,1)%gpm(2)               ! height lower level
 node(2,1) = col(1,1)%gpm(1)               ! height upper level
 node(1,2) = RefracProfile(Slant%Nmod-1,2) ! refractivity lower level
 node(2,2) = RefracProfile(Slant%Nmod,2)   ! refractivity upper level
 node_tl = 0.0_wp
 node_tl(1,2) = RefracProfile_tl(Slant%Nmod-1,2) ! refractivity lower level
 node_tl(2,2) = RefracProfile_tl(Slant%Nmod,2)   ! refractivity upper level
 do i=Slant%Nmod+1, Slant%Ntot
    RefracProfile(i,2)    = ExpInt1D(slant%SlantPointsGeo(i,3), node)
    RefracProfile_tl(i,2) = ExpInt1D_tl(slant%SlantPointsGeo(i,3), node, &
                                        0.0_wp, node_tl)
 end do

 ! Integrate refractivity from GNSS station to Hmax
 delay    = 1.0e-6_wp * IntegPolyCube(RefracProfile)
 delay_tl = 1.0e-6_wp * IntegPolyCube_tl(RefracProfile, RefracProfile_tl)

 ! for backward compatibility ....
 slant%STD         = delay
 Slant%Obs%mSTD    = delay
 Slant%Obs%STDline = delay

 end subroutine ztd_delay_tl


 !---------------------------------------------------------------------
 ! subroutine ztd_delay_ad
 !---------------------------------------------------------------------
 !
 !> @brief Adjoint routine for ztd_delay
 !>
 !> <b> call ztd_delay_ad (slant, col, delay, delay_ad, col_ad, err) </b>
 !>
 !> The adjoint code contains the original code from ztd_delay and
 !> the derivatives.
 !>
 !> This is an optimized ZTD operator which works with one single model
 !> column. The interpolation is done by the 3D-Var interface and the
 !> ZTD operator consists of three parts: \n
 !> 1) Computation of the refractivity on the model levels, \n
 !> 2) vertical extension of the profile above the model top and \n
 !> 3) vertical numerical integration which provides the ZTD.
 !>
 !> It is assumed that the GNSS station height as well as the geopotential
 !> height are heights above geoid, i.e. mean sea level, and no geoid
 !> correction and coordinate transformations are required. \n
 !> => use slant%Station%CoordGeo(3) as station height,
 !>        not CoordEll(3) as in the other delay routines \n
 !> => use slant%SlantPointsGeo to define the vertical profile,
 !>        not slant%SlantPointsEll
 !>
 !> @param[in,out] slant  derived type which contains all information
 !>                       about one STD or ZTD
 !> @param[in]     col    model column: coordinates and model fields
 !>                       on the given grid node
 !> @param[out]    delay  computed ZTD in meter \n
 !>                       The computed delays are also written to the slant
 !>                       structure.
 !> @param[in]     delay_ad  variation of delay
 !> @param[out]    col_ad    derivatives computed by the adjoint operator
 !> @param[out]    err       return code, no errors if err = 0
 !
 !---------------------------------------------------------------------
 subroutine ztd_delay_ad (slant, col, delay, delay_ad, col_ad, err)

 type (SlantData) ,intent(inout) :: slant        ! slant meta data
 type (p_column)  ,intent(in)    :: col (:,:)    ! model columns     (input)
 real(wp)         ,intent(out)   :: delay        ! Slant Total Delay (output)
 real(wp)         ,intent(in)    :: delay_ad     ! variation of delay (input)
 type (p_column)  ,pointer       :: col_ad (:,:) ! adjoint (output)
 integer          ,intent(out)   :: err

 integer  :: i, j, k, kstat, n
 real(wp) :: Delta, Di, Hi
 real(wp) :: e, Refrac, e_ad, p_ad, t_ad !, Refrac_ad
 real(wp) :: STD_ad, RefPt_ad !, ExpInt_ad
 real(wp), dimension(:,:), allocatable     :: RefracProfile, RefracProfile_ad
 real(wp), dimension(1:2,1:2) :: Anode, Anode_ad, Bnode, Bnode_ad
!real(wp) :: TotalDelay, LineDelay
 !-------------------------------------------------------------------------

 err = 0

 !----------------------------------------------------------------------------
 ! re-define supporting points according to "col"
 ! This is necessary as "ztd_path_coord" provides only the station coordinates
 !----------------------------------------------------------------------------

 ! Find column levels above the station and set the number of supporting points
 slant% Nup = NStepVertTop
 do k= ubound((col(1,1)%gpm),1), 1, -1
    if (slant%Station%CoordGeo(3) <= col(1,1)%gpm(k)) then
       ! GNSS station height below level k
       kstat = k
       ! set the number of supporting points inside the model to the number of
       ! col-levels above the GNSS station
       ! The namelist variable "NStepVertMod" is meaningless in this case and
       ! will be ignored!
       slant% Nmod = kstat + 1
       slant% Ntot = slant% Nmod + slant% Nup
       exit
    end if
 end do

 ! delete preliminary settings from "ztd_path_coord"
 if ( associated(slant%SlantSteps) ) deallocate(slant%SlantSteps)
 if ( associated(slant%SlantPointsCart) ) deallocate(slant%SlantPointsCart)
 if ( associated(slant%SlantPointsGeo) ) deallocate(slant%SlantPointsGeo)

 ! define supporting points as column levels above the station + NStepVertTop
 allocate ( slant% SlantSteps(1:slant% Ntot) )
 allocate ( slant% SlantPointsGeo(1:slant% Ntot,1:3) )

 slant%SlantSteps(1)           = 0.0_wp  ! distance to receiver in m
 ! use lon/lat (rad) of column for all supporting points
 forall (i=1:slant% Ntot)
    slant% SlantPointsGeo(i,1:2) = (/col(1,1)%dlon, col(1,1)%dlat/) * deg2rad
 end forall
 ! use GNSS station height in m as height of first point
 slant% SlantPointsGeo(1,3)    = slant%Station%CoordGeo(3)

 n = 1   ! array index of  slant%SlantSteps and SlantPointsGeo
 do k=kstat, 1, -1
    n = n + 1
    slant%SlantSteps(n)        = col(1,1)%gpm(k) - slant%Station%CoordGeo(3)
    slant% SlantPointsGeo(n,3) = col(1,1)%gpm(k)
 end do

 ! extended profile above model top
 ! heights increase by Di in each step
 Delta = col(1,1)%gpm(1) - col(1,1)%gpm(2)
 Di    = slant%Nup**(-2.0_wp) * (Hmax - col(1,1)%gpm(1) - slant%Nup*Delta)
 do i=n+1 , slant%Ntot
    j = i - n
    Hi = col(1,1)%gpm(1) + j*(Delta + j*Di)
    slant%SlantSteps(i)       = Hi - slant%Station%CoordGeo(3)
    slant%SlantPointsGeo(i,3) = Hi
 end do

 !------------------------------------------------------------------
 ! Estimate the zenith total delay (ZTD):
 ! Compute the refractivity N on the supporting points and integrate
 !------------------------------------------------------------------

 allocate( RefracProfile(1:slant%Ntot,1:2) )

 RefracProfile(:,1) = slant%SlantSteps
 RefracProfile(:,2) = 0.0_wp

 if (kstat < ubound((col(1,1)%gpm),1)) then
    ! station between two model levels:
    ! Compute refractivity also at level below station
    kstat = kstat + 1
    n = 0
 else
    ! station below lowest model level:
    ! extrapolate refractivity later
    n = 1
    RefracProfile(1,2) = invalid
 end if

 do k=kstat, 1, -1
    n = n + 1
    ! Partial pressure of water vapour e
    e = (col(1,1)%p(k) * col(1,1)%q(k)) / (RDRD +   &
               col(1,1)%q(k)*EMRDRD)
    ! Refractivity N
    Refrac = NWein(col(1,1)%p(k), col(1,1)%t(k), e)
    RefracProfile(n,2) = Refrac
 end do

 Anode(1,1) = col(1,1)%gpm(kstat)    ! height lower level
 Anode(2,1) = col(1,1)%gpm(kstat-1)  ! height upper level
 if (RefracProfile(1,2) > 0.0_wp) then
    ! Interpolate refractivity at GNSS station
    Anode(1,2) = RefracProfile(1,2)     ! refractivity lower level
    Anode(2,2) = RefracProfile(2,2)     ! refractivity upper level
 else
    ! Extrapolate refractivity below model to GNSS station
    Anode(1,2) = RefracProfile(2,2)     ! refractivity lower level
    Anode(2,2) = RefracProfile(3,2)     ! refractivity upper level
 end if
 RefracProfile(1,2) = ExpInt1D(slant%SlantPointsGeo(1,3), Anode)

 ! extrapolate refractivity above model top
 Bnode(1,1) = col(1,1)%gpm(2)               ! height lower level
 Bnode(2,1) = col(1,1)%gpm(1)               ! height upper level
 Bnode(1,2) = RefracProfile(Slant%Nmod-1,2) ! refractivity lower level
 Bnode(2,2) = RefracProfile(Slant%Nmod,2)   ! refractivity upper level
 do i=Slant%Nmod+1, Slant%Ntot
    RefracProfile(i,2) = ExpInt1D(slant%SlantPointsGeo(i,3), Bnode)
 end do

 ! Integrate refractivity from GNSS station to Hmax
 delay = 1.0e-6_wp * IntegPolyCube(RefracProfile)

 ! for backward compatibility ....
 slant%STD         = delay
 Slant%Obs%mSTD    = delay
 Slant%Obs%STDline = delay

 ! ================== end of recomputation ==============================

 allocate( RefracProfile_ad(1:slant%Ntot,1:2) )

 ! Integrate refractivity from GNSS station to Hmax
 STD_ad = delay_ad * 1.0E-6_wp
 call IntegPolyCube_ad(RefracProfile, RefracProfile_ad, STD_ad)

 ! extrapolate refractivity above model top
 ! Bnode has not been reused and should be in the same state as above.
 Bnode_ad = 0.0_wp
 do i=Slant%Ntot, Slant%Nmod+1, -1
    call ExpInt1D_ad(slant%SlantPointsGeo(i,3), Bnode,           &
                     RefracProfile_ad(i,2), RefPt_ad, Bnode_ad)
 end do
 RefracProfile_ad(Slant%Nmod,2)   = RefracProfile_ad(Slant%Nmod,2) &
                                    + Bnode_ad(2,2)
 RefracProfile_ad(Slant%Nmod-1,2) = RefracProfile_ad(Slant%Nmod-1,2) &
                                    + Bnode_ad(1,2)

 ! Interpolate refractivity at GNSS station
 ! Anode has not been reused and should be in the same state as above.
 Anode_ad = 0.0_wp
 call ExpInt1D_ad(slant%SlantPointsGeo(1,3), Anode, &
                  RefracProfile_ad(1,2), RefPt_ad, Anode_ad)

 ! Reset kstat to original state
 kstat = slant%Nmod - 1
 if (kstat < ubound((col(1,1)%gpm),1)) then
    ! station between two model levels:
    ! Compute refractivity also at level below station
    kstat = kstat + 1
    n = 0
 else
    ! station below lowest model level:
    ! extrapolate refractivity later
    n = 1
    RefracProfile(1,2) = invalid
 end if

 if (RefracProfile(1,2) > 0.0_wp) then
    ! Interpolate refractivity at GNSS station
    RefracProfile_ad(1,2) = RefracProfile_ad(1,2) + Anode_ad(1,2)
    RefracProfile_ad(2,2) = RefracProfile_ad(2,2) + Anode_ad(2,2)
 else
    ! Extrapolate refractivity below model to GNSS station
    RefracProfile_ad(2,2) = RefracProfile_ad(2,2) + Anode_ad(1,2)
    RefracProfile_ad(3,2) = RefracProfile_ad(3,2) + Anode_ad(2,2)
 end if

 do k=kstat, 1, -1
    n = n + 1
    ! Re-compute partial pressure of water vapour e
    e = (col(1,1)%p(k) * col(1,1)%q(k)) / (RDRD +   &
               col(1,1)%q(k)*EMRDRD)

    ! Refractivity N
    call NWein_ad(col(1,1)%p(k), col(1,1)%t(k), e,           &
                  p_ad, t_ad, e_ad,  RefracProfile_ad(n,2))
    col_ad(1,1)%p(k) = col_ad(1,1)%p(k) + p_ad
    col_ad(1,1)%t(k) = col_ad(1,1)%t(k) + t_ad

    ! Partial pressure of water vapour e
    col_ad(1,1)%p(k) = col_ad(1,1)%p(k) +                     &
                       ( col(1,1)%q(k) / (RDRD +              &
                         col(1,1)%q(k)*EMRDRD) ) * e_ad
    col_ad(1,1)%q(k) = col_ad(1,1)%q(k) +                        &
                       ((( RDRD + col(1,1)%q(k) * EMRDRD ) *     &
                        col(1,1)%p(k) -                          &
                       col(1,1)%p(k) * col(1,1)%q(k) *  EMRDRD)  &
                        / (RDRD + col(1,1)%q(k)*EMRDRD)**2) *    &
                        e_ad
 end do

 end subroutine ztd_delay_ad


 !---------------------------------------------------------------------
 ! subroutine std_delay_ray3d
 !---------------------------------------------------------------------
 !
 !> @brief Computes the slant total delay along the signal path.
 !>
 !> <b> call std_delay_ray3d (slant, maxelev, col, delay, ladj, col_ad, err) </b>
 !>
 !> The slant total delay (STD) along the curved signal path is computed.
 !> The raytracer is called to estimate the curved signal path, the
 !> refractivity along this path is computed and the numerical integration
 !> provides the STD. \n
 !> The raytracer is called only if the elevation of the STD is less than
 !> or equal to the maximum elevation given in maxelev. For higher elevations
 !> it is assumed that the "true" signal path is very close to the
 !> connecting line and no raytracing is required.
 !>
 !> @param[in,out] slant  derived type which contains all information
 !>                       about one STD
 !> @param[in]    maxelev maximum elevation for raytracing, the raytracer is
 !>                       called only if slant\%Obs\%elevation <= maxelev
 !> @param[in]     col    model columns: coordinates and model fields
 !>                       on the required grid nodes
 !> @param[in]     ladj   adjoint flag: Adjoint operator is called only if
 !>                       ladj = .true.
 !> @param[out]    delay  computed STD in meter \n
 !>                       The computed delays are also written to the slant
 !>                       structure: \n
 !>                       slant%STD = delay is the "true" delay along the
 !>                                   curved signal path \n
 !>                       slant%STDline (optional) is the delay along the
 !>                       connecting line between satellie and GNSS receiver
 !> @param[out]    col_ad derivatives computed by the adjoint operator
 !>                       (only if ladj = .true.)
 !> @param[out]    err    return code, no errors if err = 0
 !
 !---------------------------------------------------------------------
 subroutine std_delay_ray3d (slant, maxelev, col, delay, ladj, col_ad, err)

 type (SlantData) ,intent(inout) :: slant        ! slant meta data
 real(wp)         ,intent(in)    :: maxelev      ! max. elev. for raytracing
 type (p_column)  ,intent(in)    :: col (:,:)    ! model columns     (input)
 real(wp)         ,intent(out)   :: delay        ! Slant Total Delay (output)
 logical          ,intent(in)    :: ladj         ! flag to calculate adjoint
 type (p_column)  ,pointer       :: col_ad (:,:) ! adjoint
 integer          ,intent(out)   :: err

 real(wp), dimension(:,:,:), pointer :: Hgeo => Null()

 integer  :: i, j
 real(wp) :: e, Refrc
 real(wp) :: TotalDelay, LineDelay

 logical :: printprofile = .false.
 logical :: checkcol = .false.
 !logical :: checkcol = .true.
 !-------------------------------------------------------------------------

 err = 0

 if (verbose >= 1) write(*,*)  'std_delay_ray3d> Start ...'
 if (verbose >= 2) then
    if (slant%Nnghb .eq. 3) write(*,*) 'GME'
    if (slant%Nnghb .eq. 4) write(*,*) 'COSMO'
end if

 if (checkcol) then
    ! check model columns
    call CheckModuleColumns(slant, col)
 end if

 ! Apply geoid corrections
 allocate(Hgeo( lbound(col,1):ubound(col,1),                       &
                lbound(col,2):ubound(col,2),                       &
                lbound(col(1,1)%gpm,1):ubound(col(1,1)%gpm,1) ) )

 do j=lbound(col,2), ubound(col,2)
    do i=lbound(col,1), ubound(col,1)
       Hgeo(i,j,:) = col(i,j)%gpm(:) + col(i,j)%geoid
    end do
 end do

 ! Use last column for extrapolation of the refractivity profile
 ! Here, the height and refractivity of the two topmost points are saved
 ! in "Slant%ExtraPolNod". The interpolation is done elsewhere using these
 ! two points.
 i = 1
 j = slant%Nmod
 ! k = 2 : Layer below top layer
 e = (col(i,j)%p(2) * col(i,j)%q(2)) / (RDRD +   &
      col(i,j)%q(2)*EMRDRD)
 Refrc =  NWein(col(i,j)%p(2), col(i,j)%t(2), e)
 Slant%ExtraPolNode(1,1) = Hgeo(i,j,2)
 Slant%ExtraPolNode(1,2) = Refrc
 ! k = 1 : top layer
 e = (col(i,j)%p(1) * col(i,j)%q(1)) / (RDRD +   &
               col(i,j)%q(1)*EMRDRD)
 Refrc =  NWein(col(i,j)%p(1), col(i,j)%t(1), e)
 Slant%ExtraPolNode(2,1) = Hgeo(i,j,1)
 Slant%ExtraPolNode(2,2) = Refrc

 LineDelay = invalid
 if (slant% Obs% elevation <= maxelev) then
    ! call raytracer to estimate the bended signal path through the atmosphere
    call SignalPath (slant, col, Hgeo, err)
    if (err /= 0) return

    ! interpolate the refractivity along the signal path
    call PathRefrac (slant, col, Hgeo)
    ! compute the STD
    call PathDelay (slant, TotalDelay)
    slant%STD = TotalDelay

    ! Test only
    call LineRefrac (Slant, col, Hgeo, LineDelay)
    slant%STDline = LineDelay
 else
    ! compute hypothetical delay along the straight line connecting
    ! satellite and receiver
    call LineRefrac (Slant, col, Hgeo, LineDelay)
    slant%STD     = LineDelay
    slant%STDline = LineDelay
 end if

 delay = slant%STD

 ! for COSMO ?? and ICON-LAM ???
 slant%Obs%mSTD    = slant%STD
 slant%Obs%STDline = slant%STDline

 Slant%Obs%Press   = invalid

 if (pl_method > 0) then
    ! compute pressure at the GNSS station => plevel for LETKF
    ! find reference point on slant at a given height above the station
    ! => height above station: Href m
    do i=2, ubound(Slant%SlantPointsEll,1)
       if (Slant%SlantPointsEll(i,3)-Slant%SlantPointsEll(1,3) > Href) then

          Slant%Obs%PosEll = Slant%SlantPointsEll(i,:)

          ! compute pressure at that point on the slant path
          call IpolModel (Slant, col, Hgeo, i, Slant%SlantPointsEll(i,:), &
                          Slant%Obs%Press)
          exit
       end if
    end do
 end if

 if (slant%STD > slant%STDline) then
    ! Error: Raytracer estimates path of minimum travel time => STD.
    !        If STD > STDline a wrong minimum was found!
    write(*,*)  'Slant info: ',                      &
        ' Station ', slant%Station%SName,            & ! station short name
        ' ID ', slant%Station%ID,                    & ! station ID
        ' time ', slant%Obs%time,                    & ! observation time
        ' Sat ', slant%Obs%satellite,                & ! observed satellite
        ' Elev ', slant%Obs%elevation*rad2deg,       & ! elevation
        ' Azi ', slant%Obs%azimuth*rad2deg,          & ! azimuth
        ' ZTD_obs ', slant%Obs%zslant,               & ! observed STD
        ' STD_obs ', slant%Obs%slant,                & ! STD mapped to zenith
        ' STD_line ', slant%STDline,                 & ! model STD, staright line
        ' STD_ray ', slant%STD,                      & ! model STD, raytracer
        'ray > line'
 else if (verbose >= 3) then
    write(*,*)  'Slant info: ',                      &
        ' Station ',    slant%Station%SName,         & ! station short name
        ' ID ',         slant%Station%ID,            & ! station ID
        ' time ',       slant%Obs%time,              & ! observation time
        ' Sat ',        slant%Obs%satellite,         & ! observed satellite
        ' Elev ',       slant%Obs%elevation*rad2deg, & ! elevation
        ' Azi ',        slant%Obs%azimuth*rad2deg,   & ! azimuth
        ' ZTD_obs ',    slant%Obs%zslant,            & ! observed STD
        ' STD_obs ',    slant%Obs%slant,             & ! STD mapped to zenith
        ' STD_gerade ', slant%STDline,               & ! model STD, line
        ' STD_ray ',    slant%STD                      ! model STD, raytracer
 end if

 if (printprofile) then
    call write_profile(slant, col)
 end if

 if (verbose >= 2) write(*,*)  'std_delay_ray3d> STD computation finished ...'
 if (verbose >= 2) write(*,*)  'std_delay_ray3d> deallocate arrays ...'

 if (associated(Hgeo)) deallocate(Hgeo)

 if (verbose >= 1) write(*,*)  'std_delay_ray3d> ... end'

 end subroutine std_delay_ray3d


  !---------------------------------------------------------------------
  ! subroutine write_profile
  !---------------------------------------------------------------------
  !
  !> @brief Saves slant profiles to disk
  !>
  !> <b> call write_profile (slant, col) </b>
  !>
  !> For each signal path several ASCII files are written which provide
  !> the coordinates and the refractivities along the signal path.
  !> This results in a large amount of output files and is very slow.
  !> Slant profiles shold be written only for testing or for plotting
  !> of a limited number of STDs.
  !>
  !> @param[in]  slant derived type which contains all information about one STD
  !> @param[in]  col   model columns: coordinates and model fields
  !>                       on the required grid nodes
  !---------------------------------------------------------------------
 subroutine write_profile (slant, col)

 type (SlantData) ,intent(in) :: slant        ! slant meta data
 type (p_column)  ,intent(in) :: col (:,:)    ! model columns     (input)

 integer  :: i, j
 integer :: stat, yunit, zunit, nunit, xunit, gunit, wunit
 character (len=6)  :: n6
 character (len=20) :: SlantFile
 integer, save      :: FileNr = 0

 real (wp), dimension(1:3) :: Cslant, Cwgs
  !---------------------------------------------------------------------

 ! Write slant path to file
 FileNr = FileNr + 1
 write(n6,'(i6.6)') FileNr
 SlantFile = 'SlantXY-'//n6//'.dat'
 yunit = GetFreeUnit()
 open(UNIT =   yunit                 , &
      FILE =   trim(SlantFile)       , &
      STATUS = 'UNKNOWN'             , &
      FORM =   'Formatted'           , &
      ACTION = 'WRITE'               , &
      IOSTAT = stat                   )

 SlantFile = 'SlantXZ-'//n6//'.dat'
 zunit = GetFreeUnit()
 open(UNIT =   zunit                 , &
      FILE =   trim(SlantFile)       , &
      STATUS = 'UNKNOWN'             , &
      FORM =   'Formatted'           , &
      ACTION = 'WRITE'               , &
      IOSTAT = stat                   )

 SlantFile = 'SlantNh-'//n6//'.dat'
 nunit = GetFreeUnit()
 open(UNIT =   nunit                 , &
      FILE =   trim(SlantFile)       , &
      STATUS = 'UNKNOWN'             , &
      FORM =   'Formatted'           , &
      ACTION = 'WRITE'               , &
      IOSTAT = stat                   )

 SlantFile = 'SlantNx-'//n6//'.dat'
 xunit = GetFreeUnit()
 open(UNIT =   xunit                 , &
      FILE =   trim(SlantFile)       , &
      STATUS = 'UNKNOWN'             , &
      FORM =   'Formatted'           , &
      ACTION = 'WRITE'               , &
      IOSTAT = stat                   )

 SlantFile = 'GeradeNx-'//n6//'.dat'
 gunit = GetFreeUnit()
 open(UNIT =   gunit                 , &
      FILE =   trim(SlantFile)       , &
      STATUS = 'UNKNOWN'             , &
      FORM =   'Formatted'           , &
      ACTION = 'WRITE'               , &
      IOSTAT = stat                   )

 SlantFile = 'SlantWGS-'//n6//'.dat'
 wunit = GetFreeUnit()
 open(UNIT =   wunit                 , &
      FILE =   trim(SlantFile)       , &
      STATUS = 'UNKNOWN'             , &
      FORM =   'Formatted'           , &
      ACTION = 'WRITE'               , &
      IOSTAT = stat                   )

 if (stat .eq. 0) then
    write(yunit,*) '# Slant ', Slant%Station%SName,    &
                           Slant%Obs%elevation*rad2deg, &
                           Slant%Obs%azimuth*rad2deg
    write(zunit,*) '# Slant ', Slant%Station%SName,    &
                           Slant%Obs%elevation*rad2deg, &
                           Slant%Obs%azimuth*rad2deg

    write(zunit,*)  '# ',  &
        ' Station ', slant%Station%SName,     & ! station short name
        ' ID ', slant%Station%ID,        & ! station ID
        ' time ', slant%Obs%time,          & ! observation time
        ' Sat ', slant%Obs%satellite,     & ! observed satellite
        ' Elev ', slant%Obs%elevation*rad2deg,     & ! elevation
        ' Azi ', slant%Obs%azimuth*rad2deg,       & ! azimuth
        ' ZTD_obs ', slant%Obs%zslant,         & ! observed STD
        ' STD_obs ', slant%Obs%slant,        & ! STD mapped to zenith
        ' STD_gerade ', slant%STDline,         & ! model STD, staright line
        ' STD_ray ',   slant%STD             ! model STD, raytracer

    write(nunit,*) '# Slant ', Slant%Station%SName,    &
                           Slant%Obs%elevation*rad2deg, &
                           Slant%Obs%azimuth*rad2deg
    write(xunit,*) '# Slant ', Slant%Station%SName,    &
                           Slant%Obs%elevation*rad2deg, &
                           Slant%Obs%azimuth*rad2deg
    write(gunit,*) '# Slant ', Slant%Station%SName,    &
                           Slant%Obs%elevation*rad2deg, &
                           Slant%Obs%azimuth*rad2deg
    !write(sunit,*) '# Hmax = ', Hmax, ' m'
    write(nunit,*) '# Ntot = ', Slant%Ntot
    write(nunit,*) '# Nmod = ', Slant%Nmod
    do j=1, Slant%Ntot
       write(yunit,'(f14.2,tr2,f12.7)') Slant%SlantSteps(j),        &
                                           Slant%SlantPath(1,j)
       write(zunit,'(f14.2,tr2,f12.7)') Slant%SlantSteps(j),        &
                                       Slant%SlantPath(2,j)
       write(nunit,'(f14.8,tr2,f14.2)') Slant%SlantRefrac(j),        &
                                       Slant%SlantPointsEll(j,3)
       write(xunit,'(f14.2,tr2,f12.7)') Slant%SlantSteps(j),        &
                                          Slant%SlantRefrac(j)
       write(gunit,'(f14.2,tr2,f12.7)') Slant%SlantSteps(j),        &
              Slant%LineRefrac(j)
    end do

    do i=1, Slant%Nmod
       !cslant = (/Slant%SlantSteps(i), Slant%SlantPath(1,i),   &
       !                                Slant%SlantPath(2,i)/ )
       cslant(1) = Slant%SlantSteps(i)
       cslant(2:3) = Slant%SlantPath(1:2,i)
       call Slant2Ellips(Cslant, Cwgs, Slant%Obs%azimuth,             &
                           Slant%Obs%elevation, Slant%Station%CoordEll, &
                           Slant%Station%CoordCart, WGS84Param)
       write(wunit,'(2(f9.4,tr2),f8.2,2(f9.4,tr2))')    &
              Cwgs(1)*rad2deg,  Cwgs(2)*rad2deg, Cwgs(3), &
              col(1,i)%dlon, col(1,i)%dlat
    end do

    close(yunit)
    close(zunit)
    close(nunit)
    close(xunit)
    close(gunit)
    close(wunit)

 end if

 end subroutine write_profile


!include 'valid.f90'

!!$ subroutine std_validate ()
!!$
!!$   implicit none
!!$
!!$   real (kind=8) ::  gmfh    ! hydrostatic mapping function
!!$   real (kind=8) ::  gmfw    ! wet mapping function
!!$
!!$   integer :: i
!!$
!!$   !---------------------------------------------------------------------
!!$
!!$
!!$   ! Map slants to zenith using the Global Mapping Function GMF
!!$   do i=1, ubound(STD,1)
!!$      if (STD(i)%assimilate) then
!!$
!!$         !call gmf (dmjd,dlat,dlon,dhgt,zd,gmfh,gmfw)
!!$         call gmf( STD(i)%Obs%time,             & ! MJD, observation time
!!$                   STD(i)%Station%CoordEll(2),  & ! station latitude
!!$                   STD(i)%Station%CoordEll(1),  & ! station longitude
!!$                   STD(i)%Station%CoordEll(3),  & ! station altitude, m
!!$                   pi05-STD(i)%Obs%elevation,   & ! zenith distance
!!$                   gmfh,gmfw)
!!$         STD(i)%Obs%SZD = STD(i)%Obs%STD / gmfh
!!$
!!$      end if
!!$   end do
!!$
!!$   ! Print ZTD differences
!!$   write(*,*)
!!$   write(*,*) 'ZTD differences: MJD, ZTDobs-ZTDmod in mm'
!!$   do i=1, ubound(STD,1)
!!$      if (STD(i)%assimilate) then
!!$         write(*,*) STD(i)%Obs%time,                    &
!!$                    1000.0D0 * (STD(i)%Obs%zslant - STD(i)%Obs%SZD)
!!$      end if
!!$   end do
!!$
!!$ end subroutine std_validate


!---------------------------------------------------------------------
! subroutine IpolModel
!---------------------------------------------------------------------
!
! Start of RoboDoc comments >>>>
!****S* Fermat/IpolModel
!
! Name
! IpolModel
!
! Call
! call IpolModel (slant, col, Hgeo, p, RefPtEll, Rfrc)
!
! Purpose
! Compute the refractivity at the position "RefPtEll" using temperature
! pressure and humidity from the weather model (provided by "col").
! The interpolation is carried ou in three steps:
!
! 1) Compute the refractivity N at the grid nodes of the cell containing
!    "RefPtEll".
! 2) Vertical interpolation to the required height (RefPtEll(3)).
! 3) Bilinear horizontal interpolation between these four points.
!
! Input
! Slant    - structure containing all required information about one slant
!            (output will also be written to this structure)
! col      - grid columns of the model field
! Hgeo     - height above ellipsoid computed for all columns "col"
!            Hgeo(i,j,:) = col(i,j)%gpm(:) + col(i,j)%geoid
! p        - index of the column to be used
! RefPtEll - ellipsoidal coordinates of the reference point, the refractivity
!            will be computed for this point
!            RefPtEll(1) - longitude, rad
!            RefPtEll(2) - latitude, rad
!            RefPtEll(3) - height above allipsoid, m
!
! Output
! Rfrc - refractivity N, N = 10^6*(n-1), n - refraction index
!
! External References
! None
!
! Record of revisions
! Date        Programmer   Description of change
! ====        ==========   =====================
! 10.09.2013  M. Bender    new
! 10.03.2014  M. Bender    longitude changed to -Pi <= lon <= +Pi
!---------------------------------------------------------------------
! Source
!---------------------------------------------------------------------
subroutine IpolModel (slant, col, Hgeo, p, RefPtEll, Rfrc)

implicit none

! List of calling arguments:
type (SlantData) ,intent(in)    :: slant        ! slant meta data
type (p_column)  ,intent(in)    :: col (:,:)    ! model columns     (input)
real(wp), dimension(:,:,:), pointer :: Hgeo     ! geoid coordinates (input)
integer, intent(in)                 :: p        ! index of point in "SlantPoints"
real(wp), dimension(1:3), intent(in) :: RefPtEll ! compute N here
real(wp)         ,intent(out)   :: Rfrc         ! refractivity (output)

! List of local variables:
real(wp), dimension(:), pointer :: Hcol => Null()
real(wp), dimension(1:3)        :: RefPtEllM

! vertical grid index of neighbored columns:
integer, dimension(:), pointer  :: k => Null()
!integer, dimension(:), pointer  :: k2 => Null()
!real(wp), dimension(:,:), pointer :: e => Null()
real(wp), dimension(:,:), pointer :: N => Null()
 real(wp), dimension(:), pointer :: Nvert => Null()
!real (wp), dimension(1:4,1:2,1:4)  :: node
!real(wp), dimension(1:2,1:2)    :: vnode
real(wp), dimension(1:4) :: Nh

real (wp) :: c1, c2, d1, d2, d3, d4
real (wp) :: ya, yb, fa, fb, dy     !, Rfrc2

integer :: i, m, h
!---------------------------------------------------------------------

allocate( k(1:slant%Nnghb) )
! allocate( k2(1:slant%Nnghb) )
!allocate( e(0:1,1:slant%Nnghb) )
allocate( N(0:1,1:slant%Nnghb) )
!!$if (IpolOpt .eq. 2) then
!!$   allocate( Nvert(1:slant%Nnghb) )
!!$end if

do i=1, slant%Nnghb
   ! Find the correct vertical indices for the next set of columns
   Hcol => Hgeo(i,p,:)
   !k(i) = LayerSearch (Hcol, slant%SlantPointsEll(p,3))
   !k(i) = LayerSearch (real(Hcol), real(RefPtEll(3)))
   k(i) = LayerSearch (Hcol, RefPtEll(3))
end do

! Convert longitude
! All transformations require longitudes 0° <= lon <= 360°
! but the model provides longitudes between  -180° <= lon <= 180°
! ????? Auch gleich in Grad umrechnen? Das spart die staendige Umrechnung ... ?????
RefPtEllM = RefPtEll
if (RefPtEll(1) .gt. 180.0*deg2rad) RefPtEllM(1) = RefPtEll(1)-360.0*deg2rad

! Compute the partial pressure of water vapour e and the refractivity
! on the corners of the grid cell, i.e 8 corners for COSMO
do i=1, slant%Nnghb
   ! neighbored columns ...

!!$       write(*,*) 'Pt ', p, ' Col. ', i, ' k = ', k(i), ' : ', &
!!$            Hgeo(i,p,k(i)), ' < ', slant%SlantPointsEll(p,3), &
!!$            ' < ', Hgeo(i,p,k(i)-1)
!!$       write(*,*) 'Druck ', col(i,p)%p(k(i)), ' > ', col(i,p)%p(k(i)-1)

   do m=0, 1
      ! upper (m=0) and lower (m=1) layer

      ! Vertical index:
      ! m=0 -  h = k(i)-1 - upper layer
      ! m=1 -  h = k(i)   - lower layer
      h = k(i) - 1 + m

      ! copy pressure values for interpolation
      N(m,i) = col(i,p)%p(h)

!!$      ! Partial pressure of water vapour e
!!$      e(m,i) = (col(i,p)%p(h) * col(i,p)%q(h)) / (0.62132_wp +   &
!!$                col(i,p)%q(h)*0.37868_wp)
!!$
!!$      ! Refractivity N
!!$      N(m,i) = NWein(col(i,p)%p(h), col(i,p)%t(h), e(m,i))

   end do
end do

!!$! Interpolation using function "BiLinearExp3D"
!!$if (slant%Nnghb .eq. 3) then
!!$   ! GME - interpolate between 3 points
!!$   node(4,:,:) = node(2,:,:)
!!$end if

!Rfrc = BiLinearExp3D(slant%SlantPointsEll(p,:), node, slant%Nnghb)
!Rfrc = BiLinearExp3D(RefPtEll(:), node, slant%Nnghb)

!write(*,*) 'N = ', N(0,1), N(1,1), ' ipol = ', RefracProfile(p,2)

! ExpInt1D
do i=1, slant%Nnghb
   ! vertical interpolation of the refractivity to the required height
   if (abs(N(1,i)) .lt. 1.0E-10_wp .and. abs(N(0,i)) .lt. 1.0E-10_wp ) then
      ! N1 = N2 = 0  => N(z) = 0
      Nh(i) = 0.0_wp
   else
      ! interpolate
      if (abs(N(1,i)) .lt. 1.0E-10_wp) then
         ! N1 = 0, N2 <> 0
         N(1,i) = N(0,i)/100.0_wp
      end if
      if (abs(N(0,i)) .lt. 1.0E-10_wp) then
         ! N2 = 0, N1 <> 0
         N(0,i) = N(1,i)/100.0_wp
      end if
      Nh(i) = N(0,i) * (N(1,i)/N(0,i))**                     &
           ((Hgeo(i,p,k(i)-1)-RefPtEllM(3))/(Hgeo(i,p,k(i)-1)-Hgeo(i,p,k(i))))
   end if
end do

! BiLinear2D
! Kopieren, Sonderbehandpung nur fuer Punkte innerhalb des Dreiecks
! => Eckpunkt kopieren

! c1 = x2 - x1
!c1 = node(2,1) - node(1,1)
c1 = (col(2,p)%dlon - col(1,p)%dlon) * deg2rad

! c2 = x4 - x3
!c2 = node(4,1) - node(3,1)
if (slant%Nnghb .eq. 4) then
   c2 = (col(4,p)%dlon - col(3,p)%dlon) * deg2rad
else
   ! slant%Nnghb = 3 : replace node 4 with node 2
   c2 = (col(2,p)%dlon - col(3,p)%dlon) * deg2rad
end if

! d1 = (x - x1)/(x2-x1) =  (x - x1)/c1
!d1 = (RefPt(1) - node(1,1))/c1
d1 = (RefPtEllM(1) - col(1,p)%dlon*deg2rad) / c1

! d2 = (x2 - x)/(x2-x1) =  (x2 - x)/c1
!d2 = (node(2,1) - RefPt(1))/c1
d2 = (col(2,p)%dlon*deg2rad - RefPtEllM(1)) / c1

! d3 = (x - x3)/(x4-x3) =  (x - x3)/c2
!d3 = (RefPt(1) - node(3,1))/c2
d3 = (RefPtEllM(1) - col(3,p)%dlon*deg2rad) / c2

! d4 = (x4 - x)/(x4-x3) =  (x4 - x)/c2
!d4 = (node(4,1) - RefPt(1)) / c2
if (slant%Nnghb .eq. 4) then
   d4 = (col(4,p)%dlon*deg2rad - RefPtEllM(1)) / c2
else
   ! slant%Nnghb = 3 : replace node 4 with node 2
   d4 = (col(2,p)%dlon*deg2rad - RefPtEllM(1)) / c2
end if

!write(*,*) 'BiLinearDiff2D> c1, c2, d1, d2, d3, d4 : ', c1, c2, d1, d2, d3, d4

!fa = d2*node(1,3) + d1*node(2,3)
!fb = d4*node(3,3) + d3*node(4,3)
!ya = d2*node(1,2) + d1*node(2,2)
!yb = d4*node(3,2) + d3*node(4,2)

!fa = d2*node(1,3) + d1*node(2,3)
fa = d2 * Nh(1) + d1 * Nh(2)

!fb = d4*node(3,3) + d3*node(4,3)
if (slant%Nnghb .eq. 4) then
   fb = d4*Nh(3) + d3*Nh(4)
else
   ! slant%Nnghb = 3 : replace node 4 with node 2
   fb = d4*Nh(3) + d3*Nh(2)
end if

!ya = d2*node(1,2) + d1*node(2,2)
ya = (d2*col(1,p)%dlat + d1*col(2,p)%dlat)*deg2rad

!yb = d4*node(3,2) + d3*node(4,2)
if (slant%Nnghb .eq. 4) then
   yb = (d4*col(3,p)%dlat + d3*col(4,p)%dlat)*deg2rad
else
   ! slant%Nnghb = 3 : replace node 4 with node 2
   yb = (d4*col(3,p)%dlat + d3*col(2,p)%dlat)*deg2rad
end if

dy = yb - ya

if (abs(dy) .lt. 1.0E-8_wp) then
   ! slant%Nnghb = 3 : reference point is very close to point 1
   Rfrc = Nh(1)
else
   ! No singularity, interpolate ...
   Rfrc = fa*(yb- RefPtEllM(2))/dy + fb*(RefPtEllM(2)-ya)/dy
end if

!    write(*,*) 'Ipol: 3dvar = ', Rfrc, ' bilin = ', RefracProfile(p,2), &
!               ' diff = ', Rfrc-RefracProfile(p,2), ' diff% = ', &
!                100.0*(Rfrc-RefracProfile(p,2))/Rfrc

!write(*,*) 'RefracModel> Rfrc - org, neu ', Rfrc, Rfrc2

deallocate( k )
!deallocate( e )
deallocate( N )
if (associated(Nvert)) deallocate( Nvert )

End subroutine IpolModel


!---------------------------------------------------------------------
! subroutine Coord2Cell
!---------------------------------------------------------------------
!
! Start of RoboDoc comments >>>>
!****S* Interpolation/Coord2Cell
!
! Name
! Coord2Cell
!
! Call
! call Coord2Cell( Coord, RefSys, LMHieader, LMCoords, CellIndex , inLM)
!
! Purpose
! Compute the grid indices of a LM cell containing a given point
!
! Input
! Coord    - position given in the coordinates of "RefSys"
! RefSys   - frame of reference used by "Coord"
! LMHeader - LM header information provided by the GRIB tables
!            Containes  the required information to compute the rotated
!            coordinates.

! polphi     3dvar: grid%dlatr
! pollam     3dvar: grid%dlonr

! RotLatLU   3dvar: grid%la1
! RotLonLU   3dvar: grid%lo1

! StepLat    3dvar: grid%dj
! StepLon    3dvar: grid%di



!
! Output
! Index    - LM grid indices (lower left corner) of the LM cell
!            containing the point "Coords".
!            Index returns in any case a valid LM index. If the point
!            is located outside the LM region, the index of the LM cell
!            closest to this point is given and "inLM" is set to a
!            value greater than 0.
! inLM     - flag indicating if the point is located inside the LM region
!            inLM = 0  -  point is inside the LM region
!            inLM = 1  -  point is in the west of the LM region
!            inLM = 2  -  point is in the north of the LM region
!            inLM = 3  -  point is in the east of the LM region
!            inLM = 4  -  point is in the south of the LM region
!            inLM = 5  -  point is above the LM region
!            inLM = 6  -  point is slightlybelow the surface layer of the LM
!                         interpolation should be possible
!                         altitude > MaxDepth
!            inLM = 7  -  point is far below the surface layer of the LM,
!                         interpolation might not be possible
!                         altitude < MaxDepth
!
! External References
! CoTrans
! lat2latrot
! lon2lonrot
!
! Record of revisions
! Date        Programmer   Description of change
! ====        ==========   =====================
! 08.09.2008  M. Bender    new
! 29.09.2008  M. Bender    check negative altitudes using MaxDepth, inLM = 7
!---------------------------------------------------------------------
! Source
!---------------------------------------------------------------------
subroutine Coord2Cell( wgspoint, CellIndex, inCOSMO,                 &
                       rlon, rlat, hflell,                           &
                       polphi, pollam, RotLatLU, RotLonLU, StepLat, StepLon)

implicit none

! List of calling arguments:
real (wp), dimension(1:3), intent(in)    :: wgspoint
integer, dimension(1:3), intent(out)     :: CellIndex
integer, intent(out)                     :: inCOSMO
real (wp), pointer                       :: hflell(:,:,:)
real (wp), pointer                       :: rlon(:,:,:,:)
real (wp), pointer                       :: rlat(:,:,:,:)
REAL (wp), INTENT (IN)                   :: polphi
REAL (wp), INTENT (IN)                   :: pollam

!REAL (wp), pointer                  :: polphi
!REAL (wp), pointer                  :: pollam

REAL (wp), INTENT (IN)                  :: RotLatLU
REAL (wp), INTENT (IN)                  :: RotLonLU

REAL (wp), INTENT (IN)                  :: StepLat
REAL (wp), INTENT (IN)                  :: StepLon



! List of local variables:
real (wp)                 :: phirot, lambdarot
!real (wp)                 :: lambda, phi, phirot2, lambdarot2
integer                   :: i0, j0, k0, k
!real (wp)                 :: RotLonLU, RotLatLU
real (wp), parameter      :: MaxDepth = -50.0

!real (wp)                :: MaxDeltaLon
!real (wp)                :: MaxDeltaLat
!real (wp)                :: TestLon, TestLat, d

logical :: check
logical :: debug = .false.
!---------------------------------------------------------------------

! debug = .true.
if (debug) then
  write(*,*) 'Coord2Cell> Start ...'
  write(*,*) 'Coord2Cell> polphi, pollam     : ', polphi, pollam
  write(*,*) 'Coord2Cell> RotLatLU, RotLonLU : ', RotLatLU, RotLonLU
  write(*,*) 'Coord2Cell> StepLat, StepLon   : ', StepLat, StepLon
end if




! The rotateted latitude/longitude grid of COSMO  is an equidistant grid
! with a constant grid spacing of LMHeader%StepLon in E-W-direction
! and LMHeader%StepLat in N-S-direction.
! Using the constant spacing it is easy to compute the indices of the
! COSMO gid cell from the rotated coordinates.

! Find the lower left corner of the LM cell containing the point
! All LM header angles are given in degree!
! Find the rotated coordinates of "point"
!phirot    = lat2latrot(wgspoint(2), wgspoint(1),                                &
!                       LMHeader%NPoleLat*deg2rad, LMHeader%NPoleLon*deg2rad)
!lambdarot = lon2lonrot(wgspoint(2), wgspoint(1),                                &
!                       LMHeader%NPoleLat*deg2rad, LMHeader%NPoleLon*deg2rad)

phirot    = phi2phirot ( wgspoint(2)*rad2deg, wgspoint(1)*rad2deg, polphi, pollam )
lambdarot = rla2rlarot ( wgspoint(2)*rad2deg, wgspoint(1)*rad2deg, polphi, pollam, 0.0_wp )

!TestLon = rlarot2rla (phirot, lambdarot, polphi, pollam, 0.0D0)
!TestLat = phirot2phi (phirot, lambdarot, polphi, pollam, 0.0D0)


!write(*,*) 'Lon ', wgspoint(1)*rad2deg, TestLon
!write(*,*) 'Lat ', wgspoint(2)*rad2deg, TestLat


! Compute the idices from the equidistant rotated grid
!i0 = 1 + floor(((lambdarot*rad22deg-LMHeader%RotLonLU)/LMHeader%StepLon))
!j0 = 1 + floor((phirot*rad22deg-LMHeader%RotLatLU)/LMHeader%StepLat)

!i0 = max(1 + floor(((lambdarot-RotLonLU)/StepLon)), 1)
!j0 = max(1 + floor((phirot-RotLatLU)/StepLat), 1)

i0 = 1 + floor(((lambdarot-RotLonLU)/StepLon))
j0 = 1 + floor((phirot-RotLatLU)/StepLat)

!if (i0 .eq. 0) i0 = 1
!if (j0 .eq. 0) j0 = 1

if (debug) then
write(*,*) 'Coord2Cell> lambdarot, phirot : ', lambdarot, phirot
write(*,*) 'Coord2Cell> lambdarot, RotLonLU, StepLon ', lambdarot, RotLonLU, StepLon, phirot-RotLatLU, &
                        (lambdarot-RotLonLU)/StepLon, floor(((lambdarot-RotLonLU)/StepLon))
write(*,*) 'Coord2Cell> i0, j0 ', i0, j0
write(*,*) 'Coord2Cell> wgspoint(1),  rlon(i0  , j0  , 1, 1) ', wgspoint(1)*rad2deg,  rlon(i0  , j0  , 1, 1)*rad2deg
write(*,*) 'Coord2Cell> wgspoint(1),  rlon(i0+1, j0  , 1, 1) ', wgspoint(1)*rad2deg,  rlon(i0+1, j0  , 1, 1)*rad2deg
write(*,*) 'Coord2Cell> wgspoint(1),  rlon(i0+1, j0+1, 1, 1) ', wgspoint(1)*rad2deg,  rlon(i0+1, j0+1, 1, 1)*rad2deg
write(*,*) 'Coord2Cell> wgspoint(2),  rlat(i0  , j0  , 1, 1) ', wgspoint(2)*rad2deg,  rlat(i0  , j0  , 1, 1)*rad2deg
write(*,*) 'Coord2Cell> wgspoint(2),  rlat(i0,   j0+1, 1, 1) ', wgspoint(2)*rad2deg,  rlat(i0  , j0+1, 1, 1)*rad2deg
write(*,*) 'Coord2Cell> wgspoint(2),  rlat(i0+1, j0+1, 1, 1) ', wgspoint(2)*rad2deg,  rlat(i0+1, j0+1, 1, 1)*rad2deg
end if
!write(*,*) LMCoords(i0  , j0  , 1), LMCoords(i0  , j0  , 2)
!write(*,*) 'Coord2Cell> i0, j0 ', i0, j0

! Some points are located in a neighboured cell due to rounding errors
!MaxDeltaLon = 1.4E-5_wp * deg2rad
!MaxDeltaLat = 9.0E-6_wp * deg2rad

!!$if (i0 .ge. 0 .and. j0 .ge. 0 .and.          &
!!$    i0 .lt. ubound(rlon,1)    .and.          &
!!$    j0 .lt. ubound(rlon,2)         ) then
!!$
!!$!if ( wgspoint(1) .lt. rlon(i0  , j0  , 1, 1) .and.          &
!!$!     wgspoint(2) .lt. rlat(i0  , j0  , 1, 1)       ) then
!!$if (abs(wgspoint(1)-rlon(i0+1, j0+1  , 1, 1)) .lt. MaxDeltaLon .and.       &
!!$    abs(wgspoint(2)-rlat(i0+1 , j0+1  , 1, 1)) .lt. MaxDeltaLat      ) then
!!$    !if (i0 .gt. 1) i0 = i0 - 1
!!$    !if (j0 .gt. 1) j0 = j0 - 1
!!$    if (i0 .lt. ubound(rlon,1)) i0 = i0 + 1
!!$    if (j0 .lt. ubound(rlat,2)) j0 = j0 + 1
!!$
!!$end if
!!$
!!$end if




!if (wgspoint(1) .lt. rlon(i0  , j0  , 1, 1)) then
!  if (abs(wgspoint(1)-rlon(i0  , j0  , 1, 1)) .lt. MaxDeltaLon) then
!    !if (i0 .gt. 1) i0 = i0 - 1
!    if (i0 .lt. ubound(rlon,1)) i0 = i0 + 1
!  end if
!end if
!if (wgspoint(1) .gt. rlon(i0+1, j0  , 1, 1)) then
!  if (i0 .lt. ubound(rlon,1)) i0 = i0 + 1
!end if
!if (wgspoint(2) .lt. rlat(i0  , j0  , 1, 1)) then
!  if (abs(wgspoint(2)-rlat(i0  , j0  , 1, 1)) .lt. MaxDeltaLat) then
!    if (j0 .gt. 1) j0 = j0 - 1
!  end if
!end if
!if (wgspoint(2) .gt. rlat(i0  , j0+1, 1, 1)) then
!  if (j0 .lt. ubound(rlat,2)) j0 = j0 + 1
!end if


! The point may not be inside the LM region: Check the resulting indices
! The indices returned by this routine should always be valid LM indices.
! If the point is outside the LM region return the indices of the LM cell
! which is closest to this point and set "inLM" accordingly.
! The point is located "inside" the LM if the cell enclosed by LM grid points:
! i0_max = LMHeader%ZonalPt - 1
! j0_max = LMHeader%MeridPt - 1
inCOSMO = 0
if (i0 .lt. 1) then
   ! Coord is located westwards to the LM region
   inCOSMO = 1
   i0 = 1
end if
if (j0 .lt. 1) then
   ! Coord is located southwards to the LM region
   inCOSMO = 4
   j0 = 1
end if
if (i0 .gt. ubound(rlon,1)-1) then
   ! Coord is located eastwards to the LM region
   inCOSMO = 3
   i0 = ubound(rlon,1)
end if
if (j0 .gt. ubound(rlon,2)-1) then
   ! Coord is located northwards to the LM region
   inCOSMO = 2
   j0 = ubound(rlon,2)
end if

! Find the correct COSMO layer
k0 = -1
do k=1, ubound(hflell,3)
   !write(*,*) 'NewPoint(3), LMCoords(i0,j0,k) : ', NewPoint(3), LMCoords(i0,j0,k)
   if ( wgspoint(3) .gt. hflell(i0,j0,k) ) then
      k0 = k
      exit
   end if
end do

! The point may be below the surface layer or above the LM region
if (k0 .eq. 1) then
   ! Point is above the LM region
   inCOSMO = 5
   !k0 = ubound(hflell,3)
else if (k0 .eq. -1) then
   ! The point is below the surface layer
   k0 = ubound(hflell,3)
   if (hflell(i0,j0,k0)-wgspoint(3) .lt. MaxDepth) then
      ! Point far below the surface, interpolation may not be possible
      inCOSMO = 7
   else
      ! Point sightly below the surface, interpolation should be possible
      inCOSMO = 6
   end if
end if

! Check
check = .false.
if (inCOSMO .eq. 0 .and. check) then
   check = .false.
   if ( wgspoint(1) .ge. rlon(i0,j0,1,1)    .and.           &
        wgspoint(1) .lt. rlon(i0+1,j0,1,1)  .and.           &
        wgspoint(2) .ge. rlat(i0,j0,1,1)    .and.           &
        wgspoint(2) .lt. rlat(i0,j0+1,1,1)  .and.           &
        wgspoint(3) .ge. hflell(i0,j0,k0)   .and.           &
        wgspoint(3) .lt. hflell(i0,j0,k0-1)        ) then

      check = .true.

   end if
   if (.not. check) then
      write(*,*) 'Coord2Cell> Error, wrong grid cell: i0=',  &
                  i0, ' j0=',j0, ' k0=', k0
      if ( .not. ( wgspoint(1) .ge. rlon(i0,j0,1,1) .and.           &
                   wgspoint(1) .lt. rlon(i0+1,j0,1,1)    ) ) then
         !d = (lambdarot-RotLonLU)/StepLon
         !if (abs(d-nint(d)) .gt. 0.08D0) then
         if ( .not. ( wgspoint(1) .ge. rlon(i0,j0+1,1,1) .and.           &
                   wgspoint(1) .lt. rlon(i0+1,j0+1,1,1)    ) ) then
         write(*,*) 'Coord2Cell> Laenge: ', rlon(i0,j0,1,1), ' < ',      &
                    wgspoint(1),  ' < ',  rlon(i0+1,j0,1,1)
         write(*,*) 'Coord2Cell> lambdarot-RotLonLU = ', lambdarot-RotLonLU, &
             (lambdarot-RotLonLU)/StepLon, floor(((lambdarot-RotLonLU)/StepLon))
         write(*,*) 'Coord2Cell> i0, j0+1 / i0+1, j0+1 ',   &
                    rlon(i0,j0+1,1,1), rlon(i0+1,j0+1,1,1)
         write(*,*) 'Coord2Cell> i0, j0 / i0+1, j0 ',   &
                    rlon(i0,j0,1,1), rlon(i0+1,j0,1,1)
         end if
      end if
      if ( .not. ( wgspoint(2) .ge. rlat(i0,j0,1,1) .and.           &
                   wgspoint(2) .lt. rlat(i0,j0+1,1,1)    ) ) then
         !d =  (phirot-RotLatLU)/StepLat
         !if (abs(d-nint(d)) .gt. 0.08D0) then
         if ( .not. ( wgspoint(2) .ge. rlat(i0+1,j0,1,1) .and.           &
                   wgspoint(2) .lt. rlat(i0+1,j0+1,1,1)    ) ) then
         write(*,*) 'Coord2Cell> Breite :', rlat(i0,j0,1,1), ' < ',      &
                    wgspoint(2),  ' < ',  rlat(i0,j0+1,1,1)
         write(*,*) 'Coord2Cell> phirot-RotLatLU = ', phirot-RotLatLU,  &
             (phirot-RotLatLU)/StepLat, floor((phirot-RotLatLU)/StepLat)
         write(*,*) 'Coord2Cell> i0, j0+1 / i0+1, j0+1 ',   &
                    rlat(i0,j0+1,1,1), rlat(i0+1,j0+1,1,1)
         write(*,*) 'Coord2Cell> i0, j0 / i0+1, j0 ',   &
                    rlat(i0,j0,1,1), rlat(i0+1,j0,1,1)
         end if
      end if
      if ( .not. ( wgspoint(3) .ge. hflell(i0,j0,k0)   .and.           &
                   wgspoint(3) .lt. hflell(i0,j0,k0-1)      ) ) then
      write(*,*) 'Coord2Cell> Hoehe:', hflell(i0,j0,k0), ' < ',      &
                  wgspoint(2),  ' < ',  hflell(i0,j0,k0-1)
   end if
   end if
end if

!write(*,*) i0, j0, k0
CellIndex = (/i0,j0,k0/)

if (debug) write(*,*) 'Coord2Cell>  ... end'

End Subroutine Coord2Cell
!****
! <<<<<<<< End of RoboDoc comments


!---------------------------------------------------------------------
! subroutine SignalPath
!---------------------------------------------------------------------
!
!> @brief Compute the bended slant path inside a given model field by calling
!> the raytracer and solving Fermat's principle.
!>
!> <b> call  SignalPath (Slant, col, Hgeo) </b>
!>
!>
!> @param[inout] Slant structure containing all required information about
!>                     one slant
!>                    (output will also be written to this structure)
!> @param[in] col  grid columns of the model field
!> @param[in] Hgeo height above ellipsoid computed for all columns "col"
!>                 Hgeo(i,j,:) = col(i,j)\%gpm(:) + col(i,j)%geoid
!> @param[out] err error, return and reject observation
!>
!> Output to the slant structure \n
!> The signal path as estimated by the raytracer is written to
!> Slant\%SlantPath. \n
!> The last column is extrapolated vertically in order to estimate
!> the refractivity field above the model domain. The two points used
!> for extrapolation are saved in Slant\%ExtraPolNode. \n
!> No horizontal extrapolation is applied.
!>
!> @verbatim
!>         Slant%ExtraPolNode((m,n): m - point, n value
!>                    m = 1 - lower point, k = 2
!>                    m = 2 - upper point, k = 1
!>                    n = 1 - height above ellipsoid in m
!>                    n = 2 - refractivity
!> @endverbatim
!
!---------------------------------------------------------------------
! Record of revisions
! Date        Programmer   Description of change
! ====        ==========   =====================
! 10.09.2013  M. Bender    new
! 19.06.2015  M. Bender    Copy data to ExtraPolNode now done in std_delay_ray3d
!---------------------------------------------------------------------
subroutine SignalPath (Slant, col, Hgeo, err)

implicit none

! List of calling arguments:
type (SlantData), intent(inout)      :: Slant
type (p_column)  ,intent(in)         :: col (:,:)  ! model columns     (input)
real(wp), dimension(:,:,:), pointer  :: Hgeo       ! geoid coordinates (input)
integer, intent(out)                 :: err

! List of local variables:
real (wp), dimension(:), pointer :: X0, X
real (wp)                        :: Delta   !, e, Refrc
!integer                          :: i, j
!---------------------------------------------------------------------

!!$! Use last column for extrapolation of the refractivity profile
!!$! Here, the height and refractivity of the two topmost points are saved
!!$! in "Slant%ExtraPolNod". The interpolation is done elsewhere using these
!!$! two points.
!!$i = 1
!!$j = slant%Nmod
!!$! k = 2 : Layer below top layer
!!$e = (col(i,j)%p(2) * col(i,j)%q(2)) / (RDRD +   &
!!$               col(i,j)%q(2)*EMRDRD)
!!$Refrc =  NWein(col(i,j)%p(2), col(i,j)%t(2), e)
!!$Slant%ExtraPolNode(1,1) = Hgeo(i,j,2)
!!$Slant%ExtraPolNode(1,2) = Refrc
!!$! k = 1 : top layer
!!$e = (col(i,j)%p(1) * col(i,j)%q(1)) / (RDRD +   &
!!$               col(i,j)%q(1)*EMRDRD)
!!$Refrc =  NWein(col(i,j)%p(1), col(i,j)%t(1), e)
!!$Slant%ExtraPolNode(2,1) = Hgeo(i,j,1)
!!$Slant%ExtraPolNode(2,2) = Refrc

!write(*,*) 'SignalPath> Node for vertical extrapolation', Slant%ExtraPolNode

err = 0

allocate(X0(1:2*Slant%Ntot-4))
allocate(X(1:2*Slant%Ntot-4))
X0 = 0.0_wp
if (associated(Slant%SlantPath)) deallocate(Slant%SlantPath)
allocate(Slant%SlantPath(1:2,1:Slant%Ntot))

!write(*,*) 'Aufruf von Newton03Ell, Punkte auf Slant: ', Slant%Ntot
Delta = 1.0e-5_wp
!call Newton02Ell (Slant, col, Hgeo, X0, Delta, X, 3, 3)
call Newton03Ell (Slant, col, Hgeo, X0, Delta, X, 3, 3, err)
!call Newton03EllTest (Slant, col, Hgeo, X0, Delta, X, 3, 3)

! Receiver position was fixed:
Slant%SlantPath(1,1) = 0.0_wp
Slant%SlantPath(2,1) = 0.0_wp
! Satellite position was fixed:
Slant%SlantPath(1,Slant%Ntot) = 0.0_wp
Slant%SlantPath(2,Slant%Ntot) = 0.0_wp

Slant%SlantPath(:,2:Slant%Ntot-1) = reshape(X,(/2,Slant%Ntot-2/))

deallocate(X0, X)

End subroutine SignalPath


!---------------------------------------------------------------------
! subroutine PathRefrac
!---------------------------------------------------------------------
!>
!> @brief Compute the refractivity N for all points on the curved signal path
!>
!> <b> call PathRefrac (Slant, col, Hgeo)  </b>
!>
!> The refractivity N is computed for each supporting point on the curved
!> signal path (inside the model and above the model). The resut is stored
!> in the array Slant\%SlantRefrac.
!>
!> @param[in,out] slant  derived type which contains all information
!>                       about one STD
!>                       (output will also be written to this structure)
!> @param[in]     col    model columns: coordinates and model fields
!>                       on the required grid nodes
!> @param[in]     Hgeo   height above ellipsoid computed for all columns "col"
!>                       Hgeo(i,j,:) = col(i,j)%gpm(:) + col(i,j)%geoid
!---------------------------------------------------------------------
! Record of revisions
! Date        Programmer   Description of change
! ====        ==========   =====================
! 10.09.2013  M. Bender    new
!---------------------------------------------------------------------
subroutine PathRefrac (Slant, col, Hgeo)

implicit none

! List of calling arguments:
type (SlantData), intent(inout) :: Slant
type (p_column)  ,intent(in)    :: col (:,:)    ! model columns     (input)
real(wp), dimension(:,:,:), pointer :: Hgeo     ! geoid coordinates (input)

! List of local variables:
real (wp), dimension(1:3) :: Cslant, Cwgs
!real (wp)                 :: a, b, c
integer                   :: i
!---------------------------------------------------------------------

! Allocate array for refractivities along the slant path
if (associated(Slant%SlantRefrac)) deallocate(Slant%SlantRefrac)
allocate (Slant%SlantRefrac(1:Slant%Ntot) )

! Compute refractivities N within the model
do i=1, Slant%Nmod
   cslant = (/Slant%SlantSteps(i), Slant%SlantPath(1,i), Slant%SlantPath(2,i)/)
   call Slant2Ellips(Cslant, Cwgs, Slant%Obs%azimuth, Slant%Obs%elevation,      &
                     Slant%Station%CoordEll, Slant%Station%CoordCart, WGS84Param)
   call RefracModel (slant, col, Hgeo, i, Cwgs, Slant%SlantRefrac(i))
end do

! Compute refractivities N above the model
do i=Slant%Nmod+1, Slant%Ntot-1
   cslant = (/Slant%SlantSteps(i), Slant%SlantPath(1,i), Slant%SlantPath(2,i)/)
   call Slant2Ellips(Cslant, Cwgs, Slant%Obs%azimuth, Slant%Obs%elevation,      &
                     Slant%Station%CoordEll, Slant%Station%CoordCart, WGS84Param)
   Slant%SlantRefrac(i) = ExpInt1D(Cwgs(3), Slant%ExtraPolNode)
end do

! N = 0 at the satellite position
Slant%SlantRefrac(Slant%Ntot) = 0.0_wp

End subroutine PathRefrac


!---------------------------------------------------------------------
! subroutine PathDelay
!---------------------------------------------------------------------
!
!> @brief Compute the path delay inside a given refractivity field
!>
!> <b> call PathDelay (Slant, TotalDelay) </b>
!>
!> The line integral \f$ STD = \int N(s) ds + (S-G) \f$
!> is solved for the curved signal path S of the given STD.
!>
!> @param[in]  Slant      derived type which contains all information
!>                        about one STD
!> @param[out] TotalDelay computed slant total delay in meter
!>
!---------------------------------------------------------------------
! Record of revisions
! Date        Programmer   Description of change
! ====        ==========   =====================
! 12.09.2013  M. Bender    new
!---------------------------------------------------------------------
subroutine PathDelay (Slant, TotalDelay)

implicit none

! List of calling arguments:
type (SlantData), intent(in)  :: Slant
real (wp),        intent(out) :: TotalDelay

! List of local variables:
real (wp), dimension(:,:), pointer :: profile
real (wp) :: DeltaX01, DeltaX02, DeltaX12
real (wp), dimension(1:3) :: lb1
real (wp) :: La1y, La1z
integer :: i, j
!---------------------------------------------------------------------

! Allocate profile for integration and copy x coordinates to first dim.
allocate( profile(1:Slant%Ntot-1,1:2) )
profile(:,1) = Slant%SlantSteps(1:Slant%Ntot-1)

do i=2, Slant%Ntot-2
   ! Compute the first derivatives dy/dy and dz/dx at all supporting
   ! points using Legendre poliomilas

   ! Lagrange basis polynomials for 3 supporting points:
   ! Compute first derivative at the center point
   DeltaX01 = Slant%SlantSteps(i-1) - Slant%SlantSteps(i)
   DeltaX02 = Slant%SlantSteps(i-1) - Slant%SlantSteps(i+1)
   DeltaX12 = Slant%SlantSteps(i) - Slant%SlantSteps(i+1)

   ! Compute polynomial for x = SlantX(i)
   ! lb1 - first derivative
   lb1(1) = DeltaX12 / (DeltaX01*DeltaX02)
   lb1(2) = (DeltaX01-DeltaX12) / (DeltaX01*DeltaX12)
   lb1(3) = -DeltaX01 / (DeltaX02*DeltaX12)

   ! Compute derivatives of Lagrange polynomials
   La1y = 0.0_wp
   La1z = 0.0_wp
   do j=0, 2
      La1y = La1y + lb1(j+1) * Slant%SlantPath(1,i+j-1)
      La1z = La1z + lb1(j+1) * Slant%SlantPath(2,i+j-1)
   end do

   ! Compute n(x,y,z) * sqrt( 1 + (dy/dx)^2 + (dz/dx)^2 )
   ! n = 1 + 1e-6 * N
   profile(i,2) = (1.0_wp + 1.0e-6_wp * Slant%SlantRefrac(i)) *  &
                  sqrt( 1.0_wp + La1y**2 + La1z**2 )
end do

! First and last point
! First point, derivative at the first point
DeltaX01 = Slant%SlantSteps(1) - Slant%SlantSteps(2)
DeltaX02 = Slant%SlantSteps(1) - Slant%SlantSteps(3)
DeltaX12 = Slant%SlantSteps(2) - Slant%SlantSteps(3)
! Compute polynomial for x = SlantX(1)
! lb1 - first derivative at first point, i.e. Slant%SlantSteps(1)
lb1(1) = (DeltaX01+DeltaX02) / (DeltaX01*DeltaX02)
lb1(2) = -(DeltaX02) / (DeltaX01*DeltaX12)
lb1(3) = DeltaX01 / (DeltaX02*DeltaX12)
! Compute derivatives of Lagrange polynomials
La1y = 0.0_wp
La1z = 0.0_wp
do j=1, 3
   La1y = La1y + lb1(j) * Slant%SlantPath(1,j)
   La1z = La1z + lb1(j) * Slant%SlantPath(2,j)
end do
profile(1,2) = (1.0_wp + 1.0e-6_wp * Slant%SlantRefrac(1)) *  &
                 sqrt( 1.0_wp + La1y**2 + La1z**2 )
!
! Last point in the neutral atmosphere, i.e. Slant%Ntot-1,
! derivative at that point
DeltaX01 = Slant%SlantSteps(Slant%Ntot-3) - Slant%SlantSteps(Slant%Ntot-2)
DeltaX02 = Slant%SlantSteps(Slant%Ntot-3) - Slant%SlantSteps(Slant%Ntot-1)
DeltaX12 = Slant%SlantSteps(Slant%Ntot-2) - Slant%SlantSteps(Slant%Ntot-1)
! Compute polynomial for x = SlantX(Slant%Ntot)
! lb1 - first derivative at last point, i.e. Slant%SlantSteps(Slant%Ntot)
lb1(1) = -(DeltaX12) / (DeltaX01*DeltaX02)
lb1(2) =  (DeltaX02) / (DeltaX01*DeltaX12)
lb1(3) = -(DeltaX12+DeltaX02) / (DeltaX02*DeltaX12)
La1y = 0.0_wp
La1z = 0.0_wp
do j=1, 3
   La1y = La1y + lb1(j) * Slant%SlantPath(1,Slant%Ntot+j-4)
   La1z = La1z + lb1(j) * Slant%SlantPath(2,Slant%Ntot+j-4)
end do
profile(Slant%Ntot-1,2) = ( 1.0_wp +                              &
                 1.0e-6_wp * Slant%SlantRefrac(Slant%Ntot-1) ) *  &
                 sqrt( 1.0_wp + La1y**2 + La1z**2 )

! Optical path length:
! WARNING: The numerical integration can lead to very wrong results
!          if the last step is much longer than all other steps!
!          The slant steps go up to Hmax, e.g. 150 km, the GNSS satellite
!          is at ~20000 km. The polynomial interpolation between these
!          two points can become very wrong and should be avoided. Therefore
!          integrate up to Slant%Ntot-1 and assume n=1.0 up to the
!          satellite.
!TotalDelay = IntegPolyCube(profile(1:Slant%Ntot-1,:))
!write(*,*) 'Int n-1 : ', TotalDelay
TotalDelay = IntegPolyCube(profile)
!write(*,*) 'Int n : ', TotalDelay
! Slant total delay = optical path length - geometric path length
TotalDelay =  TotalDelay -  Slant%SlantSteps(Slant%Ntot-1)
!write(*,*) 'Delay bi N-1 : ', TotalDelay

! Contribution of the last interval up to the satellite:
DeltaX01 = Slant%SlantSteps(Slant%Ntot) - Slant%SlantSteps(Slant%Ntot-1)
TotalDelay =  TotalDelay + sqrt(DeltaX01**2 +                         &
                                Slant%SlantPath(1,Slant%Ntot-1)**2 +  &
                                Slant%SlantPath(2,Slant%Ntot-1)**2  ) &
                         - DeltaX01

deallocate( profile )

End subroutine PathDelay


!---------------------------------------------------------------------
! subroutine LineRefrac
!---------------------------------------------------------------------
!>
!> @brief Computes the refractivity N for all supporting points
!>
!> <b>  call LineRefrac (Slant, col, Hgeo, LineDelay)  </b>
!>
!> Computes the refractivity N for all supporting points along the
!> satellite-receiver axis, i.e. for all coordinates in "SlantPointsEll".
!> The hypothetical slant total delay (STD) along this straight line is
!> also provided.
!>
!> @param[in,out] slant  derived type which contains all information
!>                       about one STD
!>                      (output will also be written to this structure)
!> @param[in]     col    model columns: coordinates and model fields
!>                       on the required grid nodes
!> @param[in]     Hgeo   height above ellipsoid computed for all columns "col"
!>                       Hgeo(i,j,:) = col(i,j)%gpm(:) + col(i,j)%geoid
!> @param[out]    LineDelay   STD along the satellite-receiver-axis \n
!>                            This is not the delay along the curved signal
!>                            path !!
!>
!---------------------------------------------------------------------
! Record of revisions
! Date        Programmer   Description of change
! ====        ==========   =====================
! 10.09.2013  M. Bender    new
!---------------------------------------------------------------------
subroutine LineRefrac (Slant, col, Hgeo, LineDelay)

implicit none

! List of calling arguments:
type (SlantData), intent(inout) :: Slant
type (p_column)  ,intent(in)    :: col (:,:)    ! model columns     (input)
real(wp), dimension(:,:,:), pointer :: Hgeo     ! geoid coordinates (input)
real (wp), intent(out)              :: LineDelay

! List of local variables:
!real (wp), dimension(1:3) :: Cslant, Cwgs
!real (wp)                 :: a, b, c
integer                   :: i
real (wp), dimension(:,:), pointer :: profile

!integer                   :: j
!real (wp)                        :: Delta, e, Refrc
!---------------------------------------------------------------------


!!$! => Nach std_delay_ray3d verschieben ???
!!$
!!$! Use last column for extrapolation of the refractivity profile
!!$! Here, the height and refractivity of the two topmost points are saved
!!$! in "Slant%ExtraPolNod". The interpolation is done elsewhere using these
!!$! two points.
!!$i = 1
!!$j = slant%Nmod
!!$! k = 2 : Layer below top layer
!!$e = (col(i,j)%p(2) * col(i,j)%q(2)) / (RDRD +   &
!!$               col(i,j)%q(2)*EMRDRD)
!!$Refrc =  NWein(col(i,j)%p(2), col(i,j)%t(2), e)
!!$Slant%ExtraPolNode(1,1) = Hgeo(i,j,2)
!!$Slant%ExtraPolNode(1,2) = Refrc
!!$! k = 1 : top layer
!!$e = (col(i,j)%p(1) * col(i,j)%q(1)) / (RDRD +   &
!!$               col(i,j)%q(1)*EMRDRD)
!!$Refrc =  NWein(col(i,j)%p(1), col(i,j)%t(1), e)
!!$Slant%ExtraPolNode(2,1) = Hgeo(i,j,1)
!!$Slant%ExtraPolNode(2,2) = Refrc


! Allocate array for refractivities along the slant axis
if (associated(Slant%LineRefrac)) deallocate(Slant%LineRefrac)
allocate (Slant%LineRefrac(1:Slant%Ntot) )

! Compute refractivities N within the model
do i=1, Slant%Nmod
   call RefracModel (slant, col, Hgeo, i,                         &
                     Slant%SlantPointsEll(i,:), Slant%LineRefrac(i) )
end do

! Compute refractivities N above the model
do i=Slant%Nmod+1, Slant%Ntot-1
   Slant%LineRefrac(i) = ExpInt1D(Slant%SlantPointsEll(i,3), Slant%ExtraPolNode)
end do

! N = 0 at the satellite position
Slant%LineRefrac(Slant%Ntot) = 0.0_wp

! Compute delay
allocate( profile(1:Slant%Ntot-1,1:2) )

profile(:,1) = Slant%SlantSteps(1:Slant%Ntot-1)
profile(:,2) = Slant%LineRefrac(1:Slant%Ntot-1)

LineDelay = IntegPolyCube(profile)
LineDelay = LineDelay * 1.0E-6_wp

deallocate( profile )

End subroutine LineRefrac


!---------------------------------------------------------------------
! subroutine BandSolv
!---------------------------------------------------------------------
!>
!> @brief solves a linear banded system: A * x = b
!>
!> <b> call BandSolv (A, b, ld, ud, code) </b>
!>
!> BandSolv solves a system of linear equations. It is assumed that A
!> is a band matrix with ud upper co-diagonals and ld lower co-diagonals.
!>
!>  Here A is a nonsingular n x (ld + ud + 1) matrix in condensed form, i.e.
!>  represented in a matrix with ld+ud+1 columns for its ld lower and ud upper
!>  co-diagonals. b denotes the right hand side of the system, and x
!>  is the solution.
!>
!>  The original code was taken from J.-P. Moreau
!>  http://jean-pierre.moreau.pagesperso-orange.fr/f_matrices.html
!>  subroutines banodec and banosol
!>
!>  Both subroutines were combined in BandSolv and modified:
!>  The array indices are starting!  with 1 and
!>  the matrix A is a  n x (ld + ud + 1) matrix.
!>
!>  Input matrix A in condensed form:
!>
!> @verbatim
!>     d  u1 u2 u3                           l3 l2 l1  d u1 u2 u3
!>     l1  d u1 u2 u3                        l3 l2 l1  d u1 u2 u3     ld = ud = 3
!>     l2 l1  d u1 u2 u3                     l3 l2 l1  d u1 u2 u3
!>     l3 l2 l1  d u1 u2 u3         => A =   l3 l2 l1  d u1 u2 u3
!>        l3 l2 l1  d u1 u2 u3               l3 l2 l1  d u1 u2 u3
!>           l3 l2 l1  d u1 u2 u3            l3 l2 l1  d u1 u2 u3
!>               .............               .. .. .. .. .. .. ..
!> @endverbatim
!>
!>  In condensed form the diagonal elements are stored as columns of the matrix A.
!>  The main diagonal is the column ld+1, the lower co-diagonals are stored in
!>  the columns 1, ..., ld, the upper co-diagonals are stored in  the columns
!>  ld+2, ..., ld+ud+1. The co-diagonals are ordered with respect to the main
!>  diagonal, i.e. the lower co-diagonals are counted from column ld to 1, ld
!>  beeing the nearest neighbor to the main diagonal (l1 in the example above).
!>  The column idices of the lower co-diagonals decrease from ld to 1
!>  while the column idices of the upper co-diagonals increase from ld+2
!>  to ld+ud+1.
!>
!>  Condensed LU decomposition, replaces the banded matrix A:
!>
!> @verbatim
!>       L3 L2 L1 U1 U2 U3 U4
!>       L3 L2 L1 U1 U2 U3 U4     The first ld columns contain L,
!>  A =  L3 L2 L1 U1 U2 U3 U4     the remaining ud+1 columns U.
!>       L3 L2 L1 U1 U2 U3 U4
!>       L3 L2 L1 U1 U2 U3 U4
!>       .. .. .. .. .. .. ..
!> @endverbatim
!>
!>  In general L and U consist of the main diagonal and ld and ud co-diagonals,
!>  i.e. would require ld+ud+2 columns. The main diagonal of L can be chosen to
!>  consist only of 1 and needs not to be stored in memory and ld+ud+1 columns
!>  are sufficient: \n
!>  L(i,i) = 1, i=1, ..., n
!>
!>  The true n x n matrices L and U are given by:
!>
!> @verbatim
!>        1                            U1 U2 U3 U4
!>       L1  1                            U1 U2 U3 U4
!>       L2 L1  1                            U1 U2 U3 U4
!>  L =  L3 L2 L1  1               U =          U1 U2 U3 U4
!>          L3 L2 L1  1                            U1 U2 U3 U4
!>             L3 L2 L1  1                            U1 U2 U3 U4
!>                .. .. .. ..                            .. .. .. ..
!> @endverbatim
!>
!>
!> @param[in,out] A     band matrix in compressed form: \n
!>                      Input: dimension(1:n,1:ld + ud + 1), i. e. the
!>                      number of columns matches the bandwidth and only the
!>                      minimum of memory is rquired. The array indices start
!>                      with 1, e.g. 1, ..., n rows \n
!>                      Output: LU decomposition of the original matrix A,
!>                              i.e. A is replaced by its LU decomposition.
!> @param[in,out] b     Input: vector conatining the right hand side of the
!>                             equation A * x = b, dimension(1:n) \n
!>                      Output: x, vector containing the result x of
!>                              A * x = b, b is replaced by x.
!> @param[in]     ld    number of lower co-diagonals
!> @param[in]     ud    number of upper co-diagonals
!> @param[out]    code  return code: \n
!>                      code = 0  -  no error, solution available \n
!>                      code = 1  -  invalid input parameter \n
!>                      code = 2  -  LU decompsition does not exist
!>
!---------------------------------------------------------------------
! Record of revisions
! Date        Programmer   Description of change
! ====        ==========   =====================
! 17.06.2013  M. Bender    new
!---------------------------------------------------------------------
Subroutine BandSolv (A, b, ld, ud, code)

implicit none

! List of calling arguments:
real (wp), dimension(:,:), intent(inout) :: A
real (wp), dimension(:), intent(inout)   :: b
integer, intent(in)                      :: ld, ud
integer, intent(out)                     :: code

! List of local variables:
integer   :: kend, kjend, jm, jk, k, j, i
integer   :: n, r
logical   :: LUexist
real (wp), parameter :: MACH_EPS = 1.0E-15_wp  ! some small value
!---------------------------------------------------------------------

! dimension of linear system: A <=> nxn matric, b <=> n-dim vector
n = size(b)
r = ubound(A,2)

if (n < 3 .or. ld < 0 .or. ud < 0 .or. ld+ud+1 > r) then
   ! Invalid input parameter
   code = 1
else
   ! LU decomposition of banded system

   LUexist = .true.

   do i = 1, n - 1                             ! loop over all rws
      kend  = min (ld + 1, n - i + 1)
      kjend = min (ud + 1, n - i + 1)

      if (abs(A(i,ld+1)) < MACH_EPS) then     ! LU decompsition does
         code = 2                                ! not exist
         LUexist = .false.
         exit
      end if

      do k=1, kend-1                             ! loop over all rows
         ! below row i
         A(k+i,ld-k+1) = A(k+i,ld-k+1) / A(i,ld+1)

         do j = 1, kjend-1
            jk = j + ld - k + 1
            jm = j + ld + 1
            A(k+i,jk) = A(k+i,jk) - A(k+i,ld-k+1) * A(i,jm)
         end do
      end do  ! k loop
   end do  ! i loop

   if (LUexist) then
      ! Solve linear banded system using the LU decomposition

      do i = 1, n - 1        ! forward substitution
         kend = min (ld + 1, n - i + 1)
         do k = 1, kend-1
            b(k+i) = b(k+i) - A(k+i,ld-k+1) * b(i)
         end do
      end do

      do i = n, 1, -1        ! back substitution
         kend = min (ud + 1, n - i + 1)
         do k = 1, kend-1
            b(i) = b(i) - A(i,k+ld+1) * b(i+k)
         end do
         b(i) = b(i) / A(i,ld+1)
      end do
      code = 0

   end if  ! if (LUexist) then

end if

end Subroutine BandSolv


!---------------------------------------------------------------------
! subroutine Newton02Ell
!---------------------------------------------------------------------
!>
!> @brief Solves a system of nonlinear equations using the Newton algorithm
!>
!> <b> call Newton02Ell (Slant, col, Hgeo, X0, Delta, X, ld, ud) </b>
!>
!> Newton02Ell solves a system of nonlinear equations using the Newton
!> algorithm and finite differences. The nonlinear equations are derived
!> for ellipsoidal coordinates. \n
!> Basic implementation, Jacobian is estimated using finite forward
!> differences, a band-shaped Jacobian is assumed.
!> The linear system of equations is solved with "BandSolv".
!> The subroutine FermatDgl  provides the nonlinear system of equations,
!> i.e. computes y = F(X).
!>
!> The algorithms can be found in the book
!>
!> Uwe Naumann
!> The Art of Differentiating Computer Programs
!> SIAM, 2012, 337 pages
!>
!> This routine is a straight forward implementation of algorithm 1.1,
!> page 3, with algorithm 1.2, page 4, for estimating the Jacobian.
!>
!> @param[in,out] slant  derived type which contains all information
!>                       about one STD
!> @param[in]     col    model columns: coordinates and model fields
!>                       on the required grid nodes
!> @param[in]     Hgeo   height above ellipsoid computed for all columns "col"
!>                       Hgeo(i,j,:) = col(i,j)%gpm(:) + col(i,j)%geoid
!> @param[in]     X0     Input: first guess of the solution \n
!>                       Output: solution, F(X) = 0
!> @param[in]     Delta  increment used for computing the forward
!>                       finite differences required by the Jacobian
!> @param[in]     ld     band matrix: number of lower co-diagonals
!> @param[in]     ud     band matrix: number of upper co-diagonals
!>
!---------------------------------------------------------------------
! Record of revisions
! Date        Programmer   Description of change
! ====        ==========   =====================
! 23.07.2013  M. Bender    new
!---------------------------------------------------------------------
subroutine Newton02Ell (Slant, col, Hgeo, X0, Delta, X, ld, ud)

!use Newton_Module, only: BandSolv

! List of calling arguments:
type (SlantData), intent(in)         :: Slant
type (p_column)  ,intent(in)         :: col (:,:)  ! model columns     (input)
real(wp), dimension(:,:,:), pointer  :: Hgeo       ! geoid coordinates (input)
real (wp), dimension(:), intent(in)  :: X0
real (wp), intent(in)                :: Delta
real (wp), dimension(:), intent(out) :: X
integer, intent(in)                  :: ld, ud

! List of local variables:
real (wp), dimension(:), pointer   :: yk, xk, dx, x2, y2
real (wp), dimension(:,:), pointer :: J
real (wp) :: absy, epsilon
integer              :: code
integer :: du, dl       !, dim
integer :: i, ii, k, n
!real    :: time1, time2, time3

! Maximum number of iterations:
integer, parameter :: MaxIter = 5
integer, parameter :: MinIter = 1
!---------------------------------------------------------------------

! array dimensions:
dl = lbound(X,1)
du = ubound(X,1)
!dim = du - dl + 1

allocate( yk(dl:du) )
allocate( y2(dl:du) )
allocate( xk(dl:du) )
allocate( x2(dl:du) )
allocate( dx(dl:du) )
allocate( J(dl:du,ld+ud+1) )

epsilon = 1.0E-9_wp
absy = invalid
k = 0

call FermatDgl(Slant, col, Hgeo, X0, yk)
xk = X0
do while ( (absy .gt. epsilon .and. k .le. MaxIter) .or. k .lt. MinIter)

   k = k + 1  ! iteration
   !write(*,*) 'Iteration ', k

   ! Compute the Jacobian F'(xk) and save it in condensed form
   ! The Jacobian is a band matrix with ld lower co-diagonals
   ! and ud upper co-diagonals
   ! Estimate derivatives using forward finite differences
   J = 0.0_wp
   do i=dl, du

      ! finite differences: ith column
      x2 = xk
      x2(i) = x2(i) + Delta
      call FermatDgl(Slant, col, Hgeo, x2, y2)

      n = 0
      do ii=i+ld, i-ud, -1
         ! copy the ith column to the Jacobian in condensed form.
         ! ii - rows which depend on the variables ld, ... , ud

         n = n + 1
         if (ii .lt. dl) cycle
         if (ii .gt. du) cycle

         !write(*,*) i, ii, y2(ii), yk(ii), (y2(ii)-yk(ii)) / Delta
         J(ii,n) = (y2(ii)-yk(ii)) / Delta

      end do
   end do

   ! Newton step: solve the linear system of equations for dx
   ! J dx = -yk
   yk = -yk
   call BandSolv (J, yk, ld, ud, code)
   dx = yk

   ! improved solution
   xk = xk + dx
   call FermatDgl(Slant, col, Hgeo, xk, yk)

   ! yk should be zero:
   absy = sqrt(dot_product(yk, yk))

end do

! return solution:
X = xk

deallocate( yk )
deallocate( y2 )
deallocate( xk )
deallocate( x2 )
deallocate( dx )
deallocate(  J )

end subroutine Newton02Ell


!---------------------------------------------------------------------
! subroutine Newton03Ell
!---------------------------------------------------------------------
!>
!> @brief Solves a system of nonlinear equations using the Newton algorithm
!>
!> <b> call Newton03Ell (Slant, col, Hgeo, X0, Delta, X, ld, ud) </b>
!>
!> Newton03Ell solves a system of nonlinear equations using the Newton
!> algorithm and the exact Jacobian. The nonlinear equations are derived
!> for ellipsoidal coordinates. \n
!> The exact Jacobian is computed by "FermatDglJacobi", no finite
!> differences required.
!> This is much faster and more accurate than "Newton02Ell" calling
!> "FermatDgl" for each
!> column of the Jacobian in order to estimate the derivatives with finite
!> differences as one single call to "FermatDglJacobi" privides the full
!> Jacobian with machine precision.
!> The linear system of equations is solved with "BandSolv".
!> The subroutine FermatDgl  provides the nonlinear system of equations,
!> i.e. computes y = F(X).
!>
!> The algorithms can be found in the book
!>
!> Uwe Naumann
!> The Art of Differentiating Computer Programs
!> SIAM, 2012, 337 pages
!>
!> This routine is a straight forward implementation of algorithm 1.1,
!> page 3, with algorithm 1.2, page 4, for estimating the Jacobian.
!>
!> @param[in,out] slant  derived type which contains all information
!>                       about one STD
!> @param[in]     col    model columns: coordinates and model fields
!>                       on the required grid nodes
!> @param[in]     Hgeo   height above ellipsoid computed for all columns "col"
!>                       Hgeo(i,j,:) = col(i,j)%gpm(:) + col(i,j)%geoid
!> @param[in]     X0     Input: first guess of the solution \n
!>                       Output: solution, F(X) = 0
!> @param[in]     Delta  increment used for computing the forward
!>                       finite differences required by the Jacobian
!> @param[in]     ld     band matrix: number of lower co-diagonals
!> @param[in]     ud     band matrix: number of upper co-diagonals
!> @param[out]    err    error during iteration, return and reject observation
!>
!---------------------------------------------------------------------
! Record of revisions
! Date        Programmer   Description of change
! ====        ==========   =====================
! 23.07.2013  M. Bender    new
!---------------------------------------------------------------------
subroutine Newton03Ell (Slant, col, Hgeo, X0, Delta, X, ld, ud, err)

! List of calling arguments:
type (SlantData), intent(in)         :: Slant
type (p_column)  ,intent(in)         :: col (:,:)  ! model columns     (input)
real(wp), dimension(:,:,:), pointer  :: Hgeo       ! geoid coordinates (input)
real (wp), dimension(:), intent(in)  :: X0
real (wp), intent(in)                :: Delta
real (wp), dimension(:), intent(out) :: X
integer, intent(in)                  :: ld, ud
integer, intent(out)                 :: err

! List of local variables:
real (wp), dimension(:), pointer   :: yk, xk, dx   !, x2, y2
real (wp), dimension(:,:), pointer :: J
real (wp) :: absy, epsilon
integer              :: code
integer :: du, dl       !, dim
integer :: k            !, i, ii, n
!real    :: time1, time2, time3

! Maximum number of iterations:
integer, parameter :: MaxIter = 10
integer, parameter :: MinIter = 2
!---------------------------------------------------------------------

if (verbose >= 1) write(*,*) 'Newton03Ell> Start ...'

err = 0

! array dimensions:
dl = lbound(X,1)
du = ubound(X,1)
!dim = du - dl + 1

allocate( yk(dl:du) )
!allocate( y2(dl:du) )
allocate( xk(dl:du) )
!allocate( x2(dl:du) )
allocate( dx(dl:du) )
allocate( J(dl:du,ld+ud+1) )

epsilon = 1.0E-9_wp
absy = invalid
k = 0

call FermatDglJacobi(Slant, col, Hgeo, X0, yk, J, err)
if (err /= 0) return
xk = X0

do while ( (absy .gt. epsilon .and. k .lt. MaxIter) .or. k .lt. MinIter)

   k = k + 1  ! iteration

   ! Newton step: solve the linear system of equations for dx
   ! J dx = -yk
   yk = -yk
   call BandSolv (J, yk, ld, ud, code)
   !dx = yk

   ! improved solution
   !xk = xk + dx
   xk = xk + yk
   call FermatDglJacobi(Slant, col, Hgeo, xk, yk, J, err)
   if (err /= 0) return

   ! yk should be zero:
   absy = sqrt(dot_product(yk, yk))

end do

! return solution:
X = xk

deallocate( yk )
deallocate( xk )
deallocate( dx )
deallocate(  J )

if (verbose >= 1) write(*,*) 'Newton03Ell> ... end'

end subroutine Newton03Ell


subroutine PrintSlantX(X)

real (wp), dimension(:), intent(in)  :: X

integer :: i, n

n = size(X) / 2

do i=1, n
   write(*,*) i, X(i), X(i+1)
end do

end subroutine PrintSlantX


!---------------------------------------------------------------------
! subroutine FermatDgl
!---------------------------------------------------------------------
!>
!> @brief Solves the difference equations defined by Fermat's principle
!>
!> <b> call FermatDgl (Slant, col, Hgeo, X, Y)  </b>
!>
!> Solves the nonlinear difference equations defined by
!> Fermat's principle, i.e. computes Y(X).
!>
!> The curved signal path is decribed by (x,y,z) in the "slant system".
!> Here, x are fixed parameters, y and z are the deviations from a straight
!> line. y and z need to be estimated in order to find the curved signal path.
!> The input vector X contains y_i and z_i for all supporting points:  \n
!> X = (y1, z1, y2, z2, y3, z3, ... , yN, zN) \n
!> The difference equations are a function of X => Y(X)
!>
!> @param[in] Slant structure containing all required information about
!>                     one slant
!>                    (output will also be written to this structure)
!> @param[in] col  grid columns of the model field
!> @param[in] Hgeo height above ellipsoid computed for all columns "col"
!>                 Hgeo(i,j,:) = col(i,j)\%gpm(:) + col(i,j)%geoid
!> @param[in]  X   X = (y1, z1, y2, z2, y3, z3, ... , yN, zN)
!> @param[out] Y   solution, Y(X)
!>
!---------------------------------------------------------------------
! Record of revisions
! Date        Programmer   Description of change
! ====        ==========   =====================
! 23.07.2013  M. Bender    new
!---------------------------------------------------------------------
subroutine FermatDgl (Slant, col, Hgeo, X, Y)

implicit none

! List of calling arguments:
type (SlantData), intent(in)         :: Slant
type (p_column)  ,intent(in)         :: col (:,:)  ! model columns     (input)
real(wp), dimension(:,:,:), pointer  :: Hgeo       ! geoid coordinates (input)
real (wp), dimension(:), intent(in)  :: X
real (wp), dimension(:), intent(out) :: Y

! List of local variables:
real (wp), dimension(:,:), pointer :: SlantYZ
real (wp), dimension(:), pointer   :: SlantX

integer :: Nvar, i, j
real (wp) :: DeltaX01, DeltaX02, DeltaX12
real (wp), dimension(1:3) :: lb1, lb2
real (wp) :: La1y, La2y, La1z, La2z
real (wp) :: RefIndex,  RefIndexX,  RefIndexY,  RefIndexZ
real (wp), dimension(1:3) :: Cslant, Cwgs
!real (wp), dimension(1:3) :: Cslant2
!real (wp), dimension(1:3) :: Cslant_tl, Cwgs_tl
real (wp), dimension(1:3,1:3) :: Cwgs_dv !, Cslant_dv
!real (wp), dimension(1:3) :: SxWgs, SyWgs, SzWgs
real (wp) :: Refrc, GradRefrcX, GradRefrcY, GradRefrcZ
!real (wp) :: Refrc2, GradRefrcX2, GradRefrcY2, GradRefrcZ2
!real (wp) :: Refrc3, GradRefrcX3, GradRefrcY3, GradRefrcZ3
real (wp) :: GradLambda, GradBeta, GradHeight
!real (wp) :: GradHeight2

real (wp), dimension(1:2,1:2) :: NodeCopy !, inode_tl

! Gradients by finite differences
!real (wp) :: DeltaL, DRefrc
!real (wp) :: FDRefrcX, FDRefrcY, FDRefrcZ
!logical :: GradientTest = .false.

real (wp), dimension(1:3,1:3) :: R1, R2

!---------------------------------------------------------------------
!GradientTest = .true.

if (verbose >= 1) write(*,*) 'FermatDgl> Start ...'

! Number of points on slant path, minus first and last point
Nvar = Slant%Ntot-2

if (verbose >= 3) then
   write(*,*) 'FermatDgl> Nref = ', Slant%Ntot,    &
                       ubound(Slant%SlantSteps,1)
   write(*,*) 'FermatDgl> Nvar = ', Nvar, ubound(X,1) / 2
end if

!inode_tl = 0.0_wp
NodeCopy = Slant%ExtraPolNode

! Set pointer on X coordinate of the slant system
allocate(SlantX(1:Slant%Ntot))
!SlantX => Slant%SlantSteps
SlantX = Slant%SlantSteps

! Map X on Y and Z coordinates of the slant system
! SlantYZ(1,i) = Y coordinates
! SlantYZ(2,i) = Z coordinates
allocate(SlantYZ(1:2,1:Nvar))
SlantYZ = reshape(X,(/2,Nvar/))

! Compute and save transformation slant system => local horizon system
call Slant2LocalHorzMat (Slant%Obs%azimuth, Slant%Obs%elevation, R1)
! Compute and save transformation local horizon system => cartesian ECEF system
call LocalHorz2CartMat(Slant%Station%CoordEll(1), Slant%Station%CoordEll(2), R2)

do i=1, Nvar
   ! Compute differential equations for y and z
   ! First and last point on slant are fixed: station and satellite
   ! => boundary value:
   ! SlantY(1) = SlantY(Nvar) = SlantZ(1) = SlantZ(Nvar) = 0

   ! Lagrange basis polynomials for 3 supporting points:
   ! Compute first and second derivatives
   DeltaX01 = SlantX(i) - SlantX(i+1)
   DeltaX02 = SlantX(i) - SlantX(i+2)
   DeltaX12 = SlantX(i+1) - SlantX(i+2)

   ! Compute polynomial for x = SlantX(i)
   ! lb1 - first derivative, lb2 - second derivative
   lb1(1) = DeltaX12 / (DeltaX01*DeltaX02)
   lb1(2) = (DeltaX01-DeltaX12) / (DeltaX01*DeltaX12)
   lb1(3) = -DeltaX01 / (DeltaX02*DeltaX12)

   ! lb2 - second derivative, does not depend on x
   lb2(1) = 2.0_wp / (DeltaX01*DeltaX02)
   lb2(2) = -2.0_wp / (DeltaX01*DeltaX12)
   lb2(3) = 2.0_wp / (DeltaX02*DeltaX12)

   ! Compute derivatives of Lagrange polynomials
   La1y = 0.0_wp
   La2y = 0.0_wp
   La1z = 0.0_wp
   La2z = 0.0_wp

   if (i .eq. 1) then
      ! Polynom on the first three points of the slant:
      ! y = z = 0 on the first point
      do j=2, 3
         La1y = La1y + lb1(j) * SlantYZ(1,i+j-2)
         La2y = La2y + lb2(j) * SlantYZ(1,i+j-2)
         La1z = La1z + lb1(j) * SlantYZ(2,i+j-2)
         La2z = La2z + lb2(j) * SlantYZ(2,i+j-2)
      end do
   else if (i .eq. Nvar) then
      ! Polynom on the last three points of the slant:
      ! y = z = 0 on the last point
      do j=1, 2
         La1y = La1y + lb1(j) * SlantYZ(1,i+j-2)
         La2y = La2y + lb2(j) * SlantYZ(1,i+j-2)
         La1z = La1z + lb1(j) * SlantYZ(2,i+j-2)
         La2z = La2z + lb2(j) * SlantYZ(2,i+j-2)
      end do
   else
      do j=1, 3
         La1y = La1y + lb1(j) * SlantYZ(1,i+j-2)
         La2y = La2y + lb2(j) * SlantYZ(1,i+j-2)
         La1z = La1z + lb1(j) * SlantYZ(2,i+j-2)
         La2z = La2z + lb2(j) * SlantYZ(2,i+j-2)
      end do
   end if

   !write(*,*) 'La1y, La2y, La1z, La2z = ', La1y, La2y, La1z, La2z

   ! Compute the refraction index n and its derivatives along
   ! the X, Y and Z axis of the slant system
   Cslant = (/SlantX(i+1),SlantYZ(1,i),SlantYZ(2,i)/)
   !call Slant2Ellips(Cslant, Cwgs,  Slant%Obs%azimuth, Slant%Obs%elevation,    &
   !                  Slant%Station%CoordEll, Slant%Station%CoordCart, WGS84Param)
   call Slant2EllipsMat (Cslant, Cwgs, R1, R2, Slant%Station%CoordCart, WGS84Param)

   ! Compute derivatives of N with respect to longitude, latitude and height
   if (i .le. Slant%Nmod-1) then
      ! interpolated refractivity inside the model
      call RefracModelGrad (slant, col, Hgeo, i+1, Cwgs, Refrc,        &
                               GradLambda, GradBeta, GradHeight)

   !write(*,*) 'FermatDgl>  Refrc, GradLambda, GradBeta, GradHeight : ', &
   !                    Refrc, GradLambda, GradBeta, GradHeight
   else
      ! extrapolated refractivity above the model
      call  ExpInt1Dgrad(Cwgs(3), NodeCopy, Refrc, GradHeight)
      GradLambda = 0.0
      GradBeta = 0.0
   end if

   call Slant2EllipsMatGrad(Cslant, Cwgs, Cwgs_dv,                           &
                            R1, R2, Slant%Station%CoordCart, WGS84Param)
   GradRefrcX = GradLambda*Cwgs_dv(1,1) + GradBeta*Cwgs_dv(1,2) +   &
                GradHeight*Cwgs_dv(1,3)
   GradRefrcY = GradLambda*Cwgs_dv(2,1) + GradBeta*Cwgs_dv(2,2) +   &
                GradHeight*Cwgs_dv(2,3)
   GradRefrcZ = GradLambda*Cwgs_dv(3,1) + GradBeta*Cwgs_dv(3,2) +   &
                GradHeight*Cwgs_dv(3,3)

   ! Covert refractivities to refraction indices
   RefIndex = 1.0D0 + Refrc * 1.0E-6_wp
   RefIndexX = GradRefrcX * 1.0E-6_wp
   RefIndexY = GradRefrcY * 1.0E-6_wp
   RefIndexZ = GradRefrcZ * 1.0E-6_wp

!!$   if (GradientTest .and. i .le. Slant%Nmod-1) then
!!$      ! Gradients by finite differences
!!$      write(*,*) 'Gradients : X ', GradRefrcX, FDRefrcX, GradRefrcX-FDRefrcX, &
!!$                            ' Y ', GradRefrcY, FDRefrcY, GradRefrcY-FDRefrcY, &
!!$                            ' Z ', GradRefrcZ, FDRefrcZ, GradRefrcZ-FDRefrcZ
!!$   end if

   !write(*,*) 'N, Nx, Ny, Nz = ', Refrc, GradRefrcX, GradRefrcY, GradRefrcZ
   !write(*,*) 'n, nx, ny, nz = ', RefIndex, RefIndexX, RefIndexY, RefIndexZ
   !write(*,*) 'Sx, WGS84   : ', SxWgs
   !write(*,*) 'Sy, WGS84   : ', SyWgs
   !write(*,*) 'Sz, WGS84   : ', SzWgs
   !write(*,*) 'Azimut, Elevation :', A, elev

   ! Difference equations obtained with Fermat's principle
   ! at the current point:
   Y(2*i-1) = La2y - ((RefIndexY/RefIndex) - (RefIndexX/RefIndex)*La1y)  &
              * (1.0_wp + La1y**2 + La1z**2)
   Y(2*i)   = La2z - ((RefIndexZ/RefIndex) - (RefIndexX/RefIndex)*La1z)  &
              * (1.0_wp + La1y**2 + La1z**2)

   !write(*,*) 'nx/n, ny/n, nz/n = ', Refrc, RefIndex, RefIndexX/RefIndex,   &
   !           RefIndexY/RefIndex, RefIndexZ/RefIndex

end do
deallocate(SlantX)
deallocate(SlantYZ)

if (verbose >= 1) write(*,*) 'FermatDgl> ... end'

End subroutine FermatDgl


!---------------------------------------------------------------------
! subroutine FermatDglJacobi
!---------------------------------------------------------------------
!>
!> @brief Solves the difference equations defined by Fermat's principle
!>
!> <b> call FermatDglJacobi (Slant, col, Hgeo, X, Y, Jc)  </b>
!>
!> Solves the nonlinear difference equations defined by
!> Fermat's principle, i.e. computes Y(X).
!>
!> The curved signal path is decribed by (x,y,z) in the "slant system".
!> Here, x are fixed parameters, y and z are the deviations from a straight
!> line. y and z need to be estimated in order to find the curved signal path.
!> The input vector X contains y_i and z_i for all supporting points:  \n
!> X = (y1, z1, y2, z2, y3, z3, ... , yN, zN) \n
!> The difference equations are a function of X => Y(X) \n
!> The Jacobian of Y is computed together with Y using the exact derivatives
!> of Y with respect to y_i and z_i.
!>
!> @param[in] Slant structure containing all required information about
!>                     one slant
!>                    (output will also be written to this structure)
!> @param[in] col  grid columns of the model field
!> @param[in] Hgeo height above ellipsoid computed for all columns "col"
!>                 Hgeo(i,j,:) = col(i,j)\%gpm(:) + col(i,j)%geoid
!> @param[in]  X   X = (y1, z1, y2, z2, y3, z3, ... , yN, zN)
!> @param[out] Y   solution, Y(X)
!> @param[out] Jc  Jacobian of Y in condensed form as required by "BandSolv"
!>                 The band width is 7: ld = ud = 3
!> @param[out] err error during iteration, return and reject observation
!
!>
!---------------------------------------------------------------------
! Record of revisions
! Date        Programmer   Description of change
! ====        ==========   =====================
! 23.07.2013  M. Bender    new
!---------------------------------------------------------------------
subroutine FermatDglJacobi (Slant, col, Hgeo, X, Y, Jc, err)

implicit none

! List of calling arguments:
type (SlantData), intent(in)         :: Slant
type (p_column)  ,intent(in)         :: col (:,:)  ! model columns     (input)
real(wp), dimension(:,:,:), pointer  :: Hgeo       ! geoid coordinates (input)
real (wp), dimension(:), intent(in)  :: X
real (wp), dimension(:), intent(out) :: Y
real (wp), dimension(:,:), pointer   :: Jc         ! Jacobian (output)
integer, intent(out)                 :: err

! List of local variables:
!real (wp), dimension(:,:), pointer :: SlantYZ
!real (wp), dimension(:), pointer   :: SlantX

real (wp), allocatable :: SlantYZ(:,:)
real (wp), allocatable :: SlantX(:)

integer :: Nvar, i, j
real (wp) :: DeltaX01, DeltaX02, DeltaX12
real (wp), dimension(1:3) :: lb1, lb2
real (wp) :: La1y, La2y, La1z, La2z
real (wp) :: RefIndex,  RefIndexX,  RefIndexY,  RefIndexZ
!real (wp) :: RefIndexXx,  RefIndexYx,  RefIndexZx
real (wp) :: RefIndexXy,  RefIndexYy,  RefIndexZy
real (wp) :: RefIndexXz,  RefIndexYz,  RefIndexZz
real (wp), dimension(1:3) :: Cslant, Cwgs
!real (wp), dimension(1:3) :: Cslant2
!real (wp), dimension(1:3) :: Cslant_tl, Cwgs_tl
real (wp), dimension(1:3,1:3) :: CwgsGrad !, Cslant_dv
real (wp), dimension(1:3,1:3,1:3) :: CwgsGrad2
!real (wp), dimension(1:3) :: SxWgs, SyWgs, SzWgs
real (wp) :: Refrc, GradRefrcX, GradRefrcY, GradRefrcZ
!real (wp) :: Refrc2, GradRefrcX2, GradRefrcY2, GradRefrcZ2
!real (wp) :: Refrc3, GradRefrcX3, GradRefrcY3, GradRefrcZ3
real (wp) :: GradLambda, GradBeta, GradHeight
real (wp) :: GradHeight2
real (wp) :: AY, AZ, B
real(wp), dimension(1:3,1:3) :: Grad2Rfrc

real (wp), dimension(1:2,1:2) :: NodeCopy !, inode_tl

! Gradients by finite differences
!real (wp) :: DeltaL, DRefrc
!real (wp) :: FDRefrcX, FDRefrcY, FDRefrcZ
!logical :: GradientTest = .false.

real (wp), dimension(1:3,1:3) :: R1, R2

integer :: Nbelow
real (wp) :: MinZ
!---------------------------------------------------------------------

if (verbose >=1) write(*,*) 'FermatDglJacobi> Start ...'

err = 0
!nullify(SlantX)
!nullify(SlantYZ)

! Number of points on slant path, minus first and last point
Nvar = Slant%Ntot-2

if (verbose >= 3) then
   write(*,*) 'FermatDgl> Nref = ', Slant%Ntot,    &
                       ubound(Slant%SlantSteps,1)
   write(*,*) 'FermatDgl> Nvar = ', Nvar, ubound(X,1) / 2
end if

!inode_tl = 0.0_wp
NodeCopy = Slant%ExtraPolNode

! Set pointer on X coordinate of the slant system
allocate(SlantX(1:Slant%Ntot))
!SlantX => Slant%SlantSteps
SlantX = Slant%SlantSteps

! Map X on Y and Z coordinates of the slant system
! SlantYZ(1,i) = Y coordinates
! SlantYZ(2,i) = Z coordinates
allocate(SlantYZ(1:2,1:Nvar))
SlantYZ = reshape(X,(/2,Nvar/))

! Check if the estimated signal path is meaningful
! => The signal path should be above the connecting line
! => If too many points are below the connecting line or
!    there are points far below the connecting line
!    then stop the iteration and reject the observation
Nbelow = count(SlantYZ(2,:) < -100.0)
MinZ = minval(SlantYZ(2,:))
if (Nbelow > 5 .or. MinZ < -500.0) then
   ! ZYX
   write(*,*) 'ZYX FermatDglJacobi Nbelow, MinZ = ', Nbelow, MinZ
   err = 2
   return
end if

! Allocate Jacobian, band matrix: band width = 7, ld = 3, ud = 3
if (associated(Jc)) deallocate(Jc)
allocate( Jc(1:2*Nvar,1:7) )
Jc = 0.0_wp

! Compute and save transformation slant system => local horizon system
call Slant2LocalHorzMat (Slant%Obs%azimuth, Slant%Obs%elevation, R1)
! Compute and save transformation local horizon system => cartesian ECEF system
call LocalHorz2CartMat(Slant%Station%CoordEll(1), Slant%Station%CoordEll(2), R2)

do i=1, Nvar
   ! Compute differential equations for y and z
   ! First and last point on slant are fixed: station and satellite
   ! => boundary value:
   ! SlantY(1) = SlantY(Nvar) = SlantZ(1) = SlantZ(Nvar) = 0

   ! Lagrange basis polynomials for 3 supporting points:
   ! Compute first and second derivatives
   DeltaX01 = SlantX(i) - SlantX(i+1)
   DeltaX02 = SlantX(i) - SlantX(i+2)
   DeltaX12 = SlantX(i+1) - SlantX(i+2)

   ! Compute polynomial for x = SlantX(i)
   ! lb1 - first derivative, lb2 - second derivative
   lb1(1) = DeltaX12 / (DeltaX01*DeltaX02)
   lb1(2) = (DeltaX01-DeltaX12) / (DeltaX01*DeltaX12)
   lb1(3) = -DeltaX01 / (DeltaX02*DeltaX12)

   ! lb2 - second derivative, does not depend on x
   lb2(1) = 2.0_wp / (DeltaX01*DeltaX02)
   lb2(2) = -2.0_wp / (DeltaX01*DeltaX12)
   lb2(3) = 2.0_wp / (DeltaX02*DeltaX12)

   ! Compute derivatives of Lagrange polynomials
   La1y = 0.0_wp
   La2y = 0.0_wp
   La1z = 0.0_wp
   La2z = 0.0_wp

   if (i .eq. 1) then
      ! Polynom on the first three points of the slant:
      ! y = z = 0 on the first point
      do j=2, 3
         La1y = La1y + lb1(j) * SlantYZ(1,i+j-2)
         La2y = La2y + lb2(j) * SlantYZ(1,i+j-2)
         La1z = La1z + lb1(j) * SlantYZ(2,i+j-2)
         La2z = La2z + lb2(j) * SlantYZ(2,i+j-2)
         ! Derivatives with respect to SlantYZ
         !DLa1y(j) = lb1(j)
      end do
   else if (i .eq. Nvar) then
      ! Polynom on the last three points of the slant:
      ! y = z = 0 on the last point
      do j=1, 2
         La1y = La1y + lb1(j) * SlantYZ(1,i+j-2)
         La2y = La2y + lb2(j) * SlantYZ(1,i+j-2)
         La1z = La1z + lb1(j) * SlantYZ(2,i+j-2)
         La2z = La2z + lb2(j) * SlantYZ(2,i+j-2)
      end do
   else
      do j=1, 3
         La1y = La1y + lb1(j) * SlantYZ(1,i+j-2)
         La2y = La2y + lb2(j) * SlantYZ(1,i+j-2)
         La1z = La1z + lb1(j) * SlantYZ(2,i+j-2)
         La2z = La2z + lb2(j) * SlantYZ(2,i+j-2)
      end do
   end if

   !write(*,*) 'La1y, La2y, La1z, La2z = ', La1y, La2y, La1z, La2z

   ! Compute the refraction index n and its derivatives along
   ! the X, Y and Z axis of the slant system
   Cslant = (/SlantX(i+1),SlantYZ(1,i),SlantYZ(2,i)/)

   call Slant2EllipsMat ( Cslant, Cwgs, R1, R2, Slant%Station%CoordCart, &
                          WGS84Param )

   ! Compute derivatives of N with respect to longitude, latitude and height
   if (i .le. Slant%Nmod-1) then

      !call RefracModelGrad (slant, col, Hgeo, i+1, Cwgs, Refrc,        &
      !                         GradLambda, GradBeta, GradHeight)
      call RefracModelGrad2 (slant, col, Hgeo, i+1, Cwgs, Refrc,        &
                             GradLambda, GradBeta, GradHeight,          &
                             Grad2Rfrc, err)
      if (Refrc .lt. 0.0_wp) &
          write(*,*) 'FermatDglJacobi> Modell i, Refrc, RefIndex = ',   &
                     i, Refrc, RefIndex

      if (err /= 0) then
         ! Error during minimization, stop iteration and reject observation
         return
      end if

   else

      !call  ExpInt1Dgrad(Cwgs(3), NodeCopy, Refrc, GradHeight)
      call  ExpInt1Dgrad2(Cwgs(3), NodeCopy, Refrc, GradHeight, GradHeight2)

      GradLambda = 0.0_wp
      GradBeta = 0.0_wp
      Grad2Rfrc = 0.0_wp
      Grad2Rfrc(3,3) = GradHeight2  ! N_height_height
      if (Refrc .lt. 0.0_wp) &
          write(*,*) 'FermatDglJacobi> Extrapol. i, Refrc, RefIndex = ', &
                     i, Refrc, RefIndex

   end if

   call Slant2EllipsMatGrad2(Cslant, Cwgs, CwgsGrad, CwgsGrad2,             &
                             R1, R2, Slant%Station%CoordCart, WGS84Param)

   ! First derivatives of the refractivity N
   GradRefrcX = GradLambda*CwgsGrad(1,1) + GradBeta*CwgsGrad(1,2) +   &
                GradHeight*CwgsGrad(1,3)
   GradRefrcY = GradLambda*CwgsGrad(2,1) + GradBeta*CwgsGrad(2,2) +   &
                GradHeight*CwgsGrad(2,3)
   GradRefrcZ = GradLambda*CwgsGrad(3,1) + GradBeta*CwgsGrad(3,2) +   &
                GradHeight*CwgsGrad(3,3)

   ! Second derivatives of the refractivity N
   ! n_xx - not required
   !RefIndexXx =
   ! n_xy
   RefIndexXy = GradLambda*CwgsGrad2(2,1,1) + GradBeta*CwgsGrad2(2,1,2) + &
                GradHeight*CwgsGrad2(2,1,3) +                             &
                Grad2Rfrc(1,1)*CwgsGrad(1,1)*CwgsGrad(2,1) +              &
                Grad2Rfrc(2,2)*CwgsGrad(1,2)*CwgsGrad(2,2) +              &
                Grad2Rfrc(3,3)*CwgsGrad(1,3)*CwgsGrad(2,3) +              &
                Grad2Rfrc(2,1)*( CwgsGrad(1,1)*CwgsGrad(2,2) +            &
                                 CwgsGrad(2,1)*CwgsGrad(1,2)   ) +        &
                Grad2Rfrc(3,1)*( CwgsGrad(1,1)*CwgsGrad(2,3) +            &
                                 CwgsGrad(2,1)*CwgsGrad(1,3)   ) +        &
                Grad2Rfrc(3,2)*( CwgsGrad(1,2)*CwgsGrad(2,3) +            &
                                 CwgsGrad(2,2)*CwgsGrad(1,3)   )
   ! n_xz
   RefIndexXz = GradLambda*CwgsGrad2(3,1,1) + GradBeta*CwgsGrad2(3,1,2) + &
                GradHeight*CwgsGrad2(3,1,3) +                             &
                Grad2Rfrc(1,1)*CwgsGrad(1,1)*CwgsGrad(3,1) +              &
                Grad2Rfrc(2,2)*CwgsGrad(1,2)*CwgsGrad(3,2) +              &
                Grad2Rfrc(3,3)*CwgsGrad(1,3)*CwgsGrad(3,3) +              &
                Grad2Rfrc(2,1)*( CwgsGrad(1,1)*CwgsGrad(3,2) +            &
                                 CwgsGrad(3,1)*CwgsGrad(1,2)   ) +        &
                Grad2Rfrc(3,1)*( CwgsGrad(1,1)*CwgsGrad(3,3) +            &
                                 CwgsGrad(3,1)*CwgsGrad(1,3)   ) +        &
                Grad2Rfrc(3,2)*( CwgsGrad(1,2)*CwgsGrad(3,3) +            &
                                 CwgsGrad(3,2)*CwgsGrad(1,3)   )
   ! n_yx - not required
!!$   RefIndexYx = GradLambda*CwgsGrad2(1,2,1) + GradBeta*CwgsGrad2(1,2,2) + &
!!$                GradHeight*CwgsGrad2(1,2,3) +                             &
!!$                Grad2Rfrc(1,1)*CwgsGrad(2,1)*CwgsGrad(1,1) +              &
!!$                Grad2Rfrc(2,2)*CwgsGrad(2,2)*CwgsGrad(1,2) +              &
!!$                Grad2Rfrc(3,3)*CwgsGrad(2,3)*CwgsGrad(1,3) +              &
!!$                Grad2Rfrc(2,1)*( CwgsGrad(2,1)*CwgsGrad(1,2) +            &
!!$                                 CwgsGrad(1,1)*CwgsGrad(2,2)   ) +        &
!!$                Grad2Rfrc(3,1)*( CwgsGrad(2,1)*CwgsGrad(1,3) +            &
!!$                                 CwgsGrad(1,1)*CwgsGrad(2,3)   ) +        &
!!$                Grad2Rfrc(3,2)*( CwgsGrad(2,2)*CwgsGrad(1,3) +            &
!!$                                 CwgsGrad(1,2)*CwgsGrad(2,3)   )

   ! n_yy
   RefIndexYy = GradLambda*CwgsGrad2(2,2,1) + GradBeta*CwgsGrad2(2,2,2) + &
                GradHeight*CwgsGrad2(2,2,3) +                             &
                Grad2Rfrc(1,1)*CwgsGrad(2,1)**2 +                         &
                Grad2Rfrc(2,2)*CwgsGrad(2,2)**2 +                         &
                Grad2Rfrc(3,3)*CwgsGrad(2,3)**2 +                         &
                2.0_wp * Grad2Rfrc(2,1) * CwgsGrad(2,1)*CwgsGrad(2,2) +   &
                2.0_wp * Grad2Rfrc(3,1) * CwgsGrad(2,1)*CwgsGrad(2,3) +   &
                2.0_wp * Grad2Rfrc(3,2) * CwgsGrad(2,2)*CwgsGrad(2,3)
   ! n_yz
   RefIndexYz = GradLambda*CwgsGrad2(3,2,1) + GradBeta*CwgsGrad2(3,2,2) + &
                GradHeight*CwgsGrad2(3,2,3) +                             &
                Grad2Rfrc(1,1)*CwgsGrad(2,1)*CwgsGrad(3,1) +              &
                Grad2Rfrc(2,2)*CwgsGrad(2,2)*CwgsGrad(3,2) +              &
                Grad2Rfrc(3,3)*CwgsGrad(2,3)*CwgsGrad(3,3) +              &
                Grad2Rfrc(2,1)*( CwgsGrad(2,1)*CwgsGrad(3,2) +            &
                                 CwgsGrad(3,1)*CwgsGrad(2,2)   ) +        &
                Grad2Rfrc(3,1)*( CwgsGrad(2,1)*CwgsGrad(3,3) +            &
                                 CwgsGrad(3,1)*CwgsGrad(2,3)   ) +        &
                Grad2Rfrc(3,2)*( CwgsGrad(2,2)*CwgsGrad(3,3) +            &
                                 CwgsGrad(3,2)*CwgsGrad(2,3)   )
   ! n_zx - not required
   !RefIndexZx =
   ! n_zz
   RefIndexZz = GradLambda*CwgsGrad2(3,3,1) + GradBeta*CwgsGrad2(3,3,2) + &
                GradHeight*CwgsGrad2(3,3,3) +                             &
                Grad2Rfrc(1,1) * (CwgsGrad(3,1)**2) +                     &
                Grad2Rfrc(2,2) * (CwgsGrad(3,2)**2) +                     &
                Grad2Rfrc(3,3) * (CwgsGrad(3,3)**2) +                     &
                2.0_wp * Grad2Rfrc(2,1) * CwgsGrad(3,1)*CwgsGrad(3,2) +   &
                2.0_wp * Grad2Rfrc(3,1) * CwgsGrad(3,1)*CwgsGrad(3,3) +   &
                2.0_wp * Grad2Rfrc(3,2) * CwgsGrad(3,2)*CwgsGrad(3,3)
   ! n_zy
   RefIndexZy = RefIndexYz

   ! Convert refractivities to refraction indices
   RefIndex = 1.0_wp + Refrc * 1.0E-6_wp
   RefIndexX = GradRefrcX * 1.0E-6_wp
   RefIndexY = GradRefrcY * 1.0E-6_wp
   RefIndexZ = GradRefrcZ * 1.0E-6_wp

   RefIndexXy = RefIndexXy * 1.0E-6_wp
   RefIndexXz = RefIndexXz * 1.0E-6_wp
   RefIndexYy = RefIndexYy * 1.0E-6_wp
   RefIndexYz = RefIndexYz * 1.0E-6_wp
   RefIndexZz = RefIndexZz * 1.0E-6_wp
   RefIndexZy = RefIndexZy * 1.0E-6_wp

   !write(*,*) 'FermatDglJacobi> i, Refrc, RefIndex = ',  i, Refrc, RefIndex

!   if (.false.) then
   ! Difference equations obtained with Fermat's principle
   ! at the current point:
   Y(2*i-1) = La2y - ((RefIndexY/RefIndex) - (RefIndexX/RefIndex)*La1y)  &
              * (1.0_wp + La1y**2 + La1z**2)
   Y(2*i)   = La2z - ((RefIndexZ/RefIndex) - (RefIndexX/RefIndex)*La1z)  &
              * (1.0_wp + La1y**2 + La1z**2)

   ! First derivatives of Y with respect to SlantYZ => Jacobian
   AY = (RefIndexY/RefIndex) - (RefIndexX/RefIndex)*La1y
   B  = 1.0_wp + La1y**2 + La1z**2
   AZ = (RefIndexZ/RefIndex) - (RefIndexX/RefIndex)*La1z

   ! Fy
   ! d Fy / d y_j-1
   Jc(2*i-1,2) = lb2(1) - ( AY * 2.0_wp*lb1(1)*La1y -             &
                             B  * (RefIndexX/RefIndex)*lb1(1) )
   ! d Fy / d z_j-1
   Jc(2*i-1,3) = -AY * 2.0_wp*lb1(1)*La1z

!!$   write(*,*) 'FermatDglJacobi> lb2(2) = ', lb2(2)
!!$   write(*,*) 'FermatDglJacobi> RefIndex**2 = ', RefIndex**2
!!$   write(*,*) 'FermatDglJacobi> 1 => ',   &
!!$        ((RefIndex*RefIndexYy-RefIndexY*RefIndexY)/RefIndex**2)
!!$   write(*,*) 'FermatDglJacobi> 2 => ',   &
!!$         (RefIndexX/RefIndex)*lb1(2)
!!$   write(*,*) 'FermatDglJacobi> 3 => ',   &
!!$        ((RefIndex*RefIndexXy-RefIndexY*RefIndexX)/RefIndex**2)

   ! d Fy / d y_j
   Jc(2*i-1,4) = lb2(2) - ( AY *  2.0_wp*lb1(2)*La1y +    &
                  B * ( ((RefIndex*RefIndexYy-RefIndexY*RefIndexY)/RefIndex**2) - &
                        (RefIndexX/RefIndex)*lb1(2) +  &
                        ((RefIndex*RefIndexXy-RefIndexY*RefIndexX)/RefIndex**2) &
                        *La1y ) )

   !write(*,*) 'FermatDglJacobiDWD> d Fy / d y_j = ', Jc(2*i-1,4)

   ! d Fy / d z_j
   Jc(2*i-1,5) = -AY * 2.0_wp*lb1(2)*La1z -    &
                 B * ( ((RefIndex*RefIndexYz-RefIndexY*RefIndexZ)/RefIndex**2) &
              - La1y*((RefIndex*RefIndexXz-RefIndexX*RefIndexZ)/RefIndex**2) )


    ! d Fy / d y_j+1
   Jc(2*i-1,6) = lb2(3) - ( AY * 2.0_wp*lb1(3)*La1y -             &
                             B  * (RefIndexX/RefIndex)*lb1(3) )
    ! d Fy / d z_j+1
   Jc(2*i-1,7) = -AY * 2.0_wp*lb1(3)*La1z

   ! Fz
   ! d Fz / d y_j-1
   Jc(2*i,1) =  -AZ * 2.0_wp*lb1(1)*La1y
   ! d Fz / d z_j-1
   Jc(2*i,2) = lb2(1) - ( AZ * 2.0_wp*lb1(1)*La1z -             &
                             B  * (RefIndexX/RefIndex)*lb1(1) )
   ! d Fz / d y_j
   Jc(2*i,3) =  -AZ * 2.0_wp*lb1(2)*La1y -    &
                 B * ( ((RefIndex*RefIndexZy-RefIndexZ*RefIndexY)/RefIndex**2) &
              - La1z*((RefIndex*RefIndexXy-RefIndexX*RefIndexY)/RefIndex**2) )

   ! d Fz / d z_j
   Jc(2*i,4) = lb2(2) - ( AZ *  2.0_wp*lb1(2)*La1z +    &
                  B * ( ((RefIndex*RefIndexZz-RefIndexZ*RefIndexZ)/RefIndex**2) - &
                        (RefIndexX/RefIndex)*lb1(2) +  &
                        ((RefIndex*RefIndexXz-RefIndexZ*RefIndexX)/RefIndex**2) &
                        *La1z ) )

   ! d Fz / d y_j+1
   Jc(2*i,5) = -AZ * 2.0_wp*lb1(3)*La1y
   ! d Fz / d z_j+1
   Jc(2*i,6) = lb2(3) - ( AZ * 2.0_wp*lb1(3)*La1z -             &
                             B  * (RefIndexX/RefIndex)*lb1(3) )

!   end if

end do

Jc(1,1:3) = 0.0_wp
Jc(2,1:2) = 0.0_wp
Jc(2*Nvar-1,6:7) = 0.0_wp
Jc(2*Nvar,5:7) = 0.0_wp

deallocate(SlantX)
deallocate(SlantYZ)

if (verbose >=1) write(*,*) 'FermatDglJacobi> ... end'

End subroutine FermatDglJacobi


!---------------------------------------------------------------------
! subroutine RefracModel
!---------------------------------------------------------------------
!
!> @brief Estimate the refractivity at some reference point by interpolating
!> between the model grid nodes.
!>
!> <b> call RefracModel (slant, col, Hgeo, p, RefPtEll, Rfrc) </b>
!>
!> Compute the refractivity at the position "RefPtEll" using temperature
!> pressure and humidity from the weather model (provided by "col").
!> The interpolation is carried out in three steps:
!>
!> 1) Compute the refractivity N at the grid nodes of the cell containing
!>    "RefPtEll".\n
!> 2) Vertical interpolation to the required height (RefPtEll(3)).\n
!> 3) Bilinear horizontal interpolation between the nodes.
!>
!> The horizontal bilinear interpolation works with 3 (ICON) or 4 (COSMO)
!> grid columns. In case of 3 columns one node is copied to the "missing"
!> fourth node and the interpolation is done in exactly the same way as with
!> four nodes. To avoid singularities some checks are made and the 3
!> nodes are reordered accordingly.
!>
!> ModuleTest.f90: Test ::RefracModel with TestRefracModel = .true.
!>
!> @param[in] slant  structure containing all required information about
!>                   one slant
!>                   (output will also be written to this structure)
!> @param[in] col    grid columns of the model field
!> @param[in] Hgeo   height above ellipsoid computed for all columns "col" \n
!>                   Hgeo(i,j,:) = col(i,j)%gpm(:) + col(i,j)%geoid
!> @param[in] p      index of the column to be used
!> @param[in] RefPtEll  ellipsoidal coordinates of the reference point,
!>                      the refractivity will be computed for this point \n
!>                      RefPtEll(1) - longitude, rad \n
!>                      RefPtEll(2) - latitude, rad \n
!>                      RefPtEll(3) - height above allipsoid, m
!>
!> @param[out] Rfrc  refractivity N, N = 10^6*(n-1), n - refraction index
!
!---------------------------------------------------------------------
! Record of revisions
! Date        Programmer   Description of change
! ====        ==========   =====================
! 10.09.2013  M. Bender    new
! 10.03.2014  M. Bender    longitude changed to -Pi <= lon <= +Pi
! 10.07.2015  M. Bender    major changes to avoid problems with two
!                          points having the same longitude.
!---------------------------------------------------------------------
subroutine RefracModel (slant, col, Hgeo, p, RefPtEll, Rfrc)

implicit none

! List of calling arguments:
type (SlantData) ,intent(in)    :: slant      ! slant meta data
type (p_column)  ,intent(in)    :: col (:,:)  ! model columns     (input)
real(wp), dimension(:,:,:), pointer :: Hgeo   ! geoid coordinates (input)
integer, intent(in)                 :: p      ! index of point in "SlantPoints"
real(wp), dimension(1:3), intent(in) :: RefPtEll ! compute N here
real(wp)         ,intent(out)   :: Rfrc         ! refractivity (output)

! List of local variables:
real(wp), dimension(:), pointer :: Hcol => Null()
real(wp), dimension(1:3)        :: RefPtEllM

! vertical grid index of neighbored columns:
integer, dimension(:), pointer  :: k => Null()
real(wp), dimension(:,:), pointer :: e => Null()
real(wp), dimension(:,:), pointer :: N => Null()
!real (wp), dimension(1:4,1:2,1:4)  :: node
!real(wp), dimension(1:2,1:2)    :: vnode
real(wp), dimension(1:4) :: lon, lat, Nh

real (wp) :: c1, c2, d1, d2, d3, d4
real (wp) :: ya, yb, fa, fb, dy     !, Rfrc2

integer :: i, m, h
!---------------------------------------------------------------------

allocate( k(1:slant%Nnghb) )
allocate( e(0:1,1:slant%Nnghb) )
allocate( N(0:1,1:slant%Nnghb) )

do i=1, slant%Nnghb
   ! Find the correct vertical indices for the next set of columns
   Hcol => Hgeo(i,p,:)
   k(i) = LayerSearch (Hcol, RefPtEll(3))
end do

! Convert longitude
! All transformations require longitudes 0° <= lon <= 360°
! but the model provides longitudes between  -180° <= lon <= 180°
! Convert RefPtEll to degrees (as in col%dlon,dlat) and to -180°<=lon<=180°
RefPtEllM(1:2) = RefPtEll(1:2) * rad2deg
RefPtEllM(3)   = RefPtEll(3)
if (RefPtEllM(1) > 180.0_wp) RefPtEllM(1) = RefPtEllM(1)-360.0_wp

! Compute the partial pressure of water vapour e and the refractivity
! on the corners of the grid cell, i.e 8 corners for COSMO
do i=1, slant%Nnghb
   ! neighbored columns ...

!!$       write(*,*) 'Pt ', p, ' Col. ', i, ' k = ', k(i), ' : ', &
!!$            Hgeo(i,p,k(i)), ' < ', slant%SlantPointsEll(p,3), &
!!$            ' < ', Hgeo(i,p,k(i)-1)
!!$       write(*,*) 'Druck ', col(i,p)%p(k(i)), ' > ', col(i,p)%p(k(i)-1)

   do m=0, 1
      ! upper (m=0) and lower (m=1) layer

      ! Vertical index:
      ! m=0 -  h = k(i)-1 - upper layer
      ! m=1 -  h = k(i)   - lower layer
      h = k(i) - 1 + m

      ! Partial pressure of water vapour e
      e(m,i) = (col(i,p)%p(h) * col(i,p)%q(h)) / (RDRD +   &
                col(i,p)%q(h)*EMRDRD)

      ! Refractivity N
      N(m,i) = NWein(col(i,p)%p(h), col(i,p)%t(h), e(m,i))

   end do
end do

do i=1, slant%Nnghb
   ! vertical interpolation of the refractivity to the required height
   if (abs(N(1,i)) .lt. 1.0E-10_wp .and. abs(N(0,i)) .lt. 1.0E-10_wp ) then
      ! N1 = N2 = 0  => N(z) = 0
      Nh(i) = 0.0_wp
   else
      ! interpolate
      if (abs(N(1,i)) .lt. 1.0E-10_wp) then
         ! N1 = 0, N2 <> 0
         N(1,i) = N(0,i)/100.0_wp
      end if
      if (abs(N(0,i)) .lt. 1.0E-10_wp) then
         ! N2 = 0, N1 <> 0
         N(0,i) = N(1,i)/100.0_wp
      end if
      Nh(i) = N(0,i) * (N(1,i)/N(0,i))**                     &
           ((Hgeo(i,p,k(i)-1)-RefPtEllM(3))/(Hgeo(i,p,k(i)-1)-Hgeo(i,p,k(i))))
   end if
   !write(*,*) 'Nh(i) = ', i, Nh(i)
end do

! Bilinear horizontal interpolation
lon(1:slant%Nnghb) = col(:,p)%dlon
lat(1:slant%Nnghb) = col(:,p)%dlat
if (abs(lon(1)) > 175.0) then
   ! avoid problems near plus/minus 180°
   ! convert longitude between +-180° to 0°-360°
   do i=1, slant%Nnghb
      if (lon(i) < 0.0) lon(i) = 360.0_wp - lon(i)
   end do
   if (RefPtEllM(1) < 0.0) RefPtEllM(1) = 360.0_wp - lon(i)
end if
if (slant%Nnghb .eq. 3) then
   ! ICON/GME - interpolate between 3 points
   if (abs(lon(1)-lon(2)) < 1.0E-6_wp) then
      ! change 2 and 3 to avoid singularity, copy 3 to 4
      !write(*,*) 'lon(1) = lon(2)'
      lon(4) = lon(3)
      lat(4) = lat(3)
      Nh(4)  = Nh(3)
      lon(3) = lon(2)
      lat(3) = lat(2)
      Nh(3)  = Nh(2)
      lon(2) = lon(4)
      lat(2) = lat(4)
      Nh(2)  = Nh(4)
   else if (abs(lon(2)-lon(3)) < 1.0E-6_wp) then
      ! copy 1 to 4 to avoid singularity
      !write(*,*) 'lon(2) = lon(3)'
      lon(4) = lon(1)
      lat(4) = lat(1)
      Nh(4)  = Nh(1)
   else
      lon(4) = lon(2)
      lat(4) = lat(2)
      Nh(4)  = Nh(2)
   end if
end if

! c1 = x2 - x1
c1 = lon(2) - lon(1)

! c2 = x4 - x3
c2 = lon(4) - lon(3)

! d1 = (x - x1)/(x2-x1) =  (x - x1)/c1
d1 = (RefPtEllM(1) - lon(1)) / c1

! d2 = (x2 - x)/(x2-x1) =  (x2 - x)/c1
d2 = (lon(2) - RefPtEllM(1)) / c1

! d3 = (x - x3)/(x4-x3) =  (x - x3)/c2
d3 = (RefPtEllM(1) - lon(3)) / c2

! d4 = (x4 - x)/(x4-x3) =  (x4 - x)/c2
d4 = (lon(4) - RefPtEllM(1)) / c2

!fa = d2*node(1,3) + d1*node(2,3)
fa = d2 * Nh(1) + d1 * Nh(2)

!fb = d4*node(3,3) + d3*node(4,3)
fb = d4*Nh(3) + d3*Nh(4)

!ya = d2*node(1,2) + d1*node(2,2)
ya = d2*lat(1) + d1*lat(2)

!yb = d4*node(3,2) + d3*node(4,2)
yb = d4*lat(3) + d3*lat(4)

dy = yb - ya

if (abs(dy) .lt. 1.0E-8_wp) then
   ! slant%Nnghb = 3 : reference point is very close to point 2
   Rfrc = Nh(2)
else
   ! No singularity, interpolate ...
   Rfrc = fa*(yb- RefPtEllM(2))/dy + fb*(RefPtEllM(2)-ya)/dy
end if

!    write(*,*) 'Ipol: 3dvar = ', Rfrc, ' bilin = ', RefracProfile(p,2), &
!               ' diff = ', Rfrc-RefracProfile(p,2), ' diff% = ', &
!                100.0*(Rfrc-RefracProfile(p,2))/Rfrc

!write(*,*) 'RefracModel> Rfrc - org, neu ', Rfrc, Rfrc2

deallocate( k )
deallocate( e )
deallocate( N )

End subroutine RefracModel


!---------------------------------------------------------------------
! subroutine RefracModelGrad
!---------------------------------------------------------------------
!
!> \brief Estimate the refractivity at some reference point by interpolating
!> between the model grid nodes and compute the derivatives with respect
!> to longitude, latitude and height.
!>
!> <b> call RefracModelGrad (slant, col, Hgeo, p, RefPtEll, Rfrc,   & \n
!>                       RfrcLambda, RfrcBeta, RfrcHeight)  </b>
!>
!>
!> Computes the refractivity at the position "RefPtEll" using temperature
!> pressure and humidity from the weather model (provided by "col") and the
!> first derivatives, i.e. additional calls to ::RefracModel are usually
!> not necessary.
!> The interpolation is carried out in three steps:
!>
!> 1) Compute the refractivity N at the grid nodes of the cell containing
!>    "RefPtEll".\n
!> 2) Vertical interpolation to the required height (RefPtEll(3)).\n
!> 3) Bilinear horizontal interpolation between the nodes.
!>
!> The horizontal bilinear interpolation works with 3 (ICON) or 4 (COSMO)
!> grid columns. In case of 3 columns one node is copied to the "missing"
!> fourth node and the interpolation is done in exactly the same way as with
!> four nodes. To avoid singularities some checks are made and the 3
!> nodes are reordered accordingly. \n
!> The derivatives of the refractivity N with respect to the
!> ellipsoidal coordinates longitude (lambda), latitude (beta) and height
!> are provided: \n
!> \f$ RfrcLambda = \frac{\partial N}{\partial \lambda} \f$,
!> \f$ RfrcBeta = \frac{\partial N}{\partial \beta} \f$ and
!> \f$ RfrcHeight = \frac{\partial N}{\partial h} \f$.
!>
!>
!> ModuleTest.f90: Test ::RefracModelGrad with TestRefracModel = .true. \n
!>                        Derivatives are tested with finite differences.
!>
!> @param[in] slant  structure containing all required information about
!>                   one slant
!>                   (output will also be written to this structure)
!> @param[in] col    grid columns of the model field
!> @param[in] Hgeo   height above ellipsoid computed for all columns "col" \n
!>                   Hgeo(i,j,:) = col(i,j)%gpm(:) + col(i,j)%geoid
!> @param[in] p      index of the column to be used
!> @param[in] RefPtEll  ellipsoidal coordinates of the reference point,
!>                      the refractivity will be computed for this point \n
!>                      RefPtEll(1) - longitude, rad \n
!>                      RefPtEll(2) - latitude, rad \n
!>                      RefPtEll(3) - height above allipsoid, m
!>
!> @param[out] Rfrc  refractivity N, N = 10^6*(n-1), n - refraction index
!> @param[out] RfrcLambda derivative \f$ \frac{\partial N}{\partial \lambda} \f$
!> @param[out] RfrcBeta   derivative \f$ \frac{\partial N}{\partial \beta} \f$
!> @param[out] RfrcHeight derivative \f$ \frac{\partial N}{\partial h} \f$
!
!---------------------------------------------------------------------
! Record of revisions
! Date        Programmer   Description of change
! ====        ==========   =====================
! 16.12.2013  M. Bender    new
! 10.07.2015  M. Bender    major changes to avoid problems with two
!                          points having the same longitude.
!---------------------------------------------------------------------
subroutine RefracModelGrad (slant, col, Hgeo, p, RefPtEll, Rfrc,   &
                            RfrcLambda, RfrcBeta, RfrcHeight)

implicit none

! List of calling arguments:
type (SlantData) ,intent(in)    :: slant        ! slant meta data
type (p_column)  ,intent(in)    :: col (:,:)    ! model columns     (input)
real(wp), dimension(:,:,:), pointer :: Hgeo     ! geoid coordinates (input)
integer, intent(in)                 :: p        ! index of pont in "SlantPoints"
real(wp), dimension(1:3), intent(in) :: RefPtEll ! compute N here
real(wp),intent(out)                :: Rfrc     ! refractivity (output)
real(wp),intent(out)                :: RfrcLambda, RfrcBeta, RfrcHeight
!real(wp),optional,intent(out)                :: A

! List of local variables:
real(wp), dimension(:), pointer :: Hcol => Null()
real(wp), dimension(1:3)        :: RefPtEllM

! vertical grid index of neighbored columns:
integer, dimension(:), pointer  :: k => Null()
real(wp), dimension(:,:), pointer :: e => Null()
real(wp), dimension(:,:), pointer :: N => Null()
!real(wp), dimension(1:3) :: RefPt_tl
real(wp), dimension(1:4) :: lon, lat, Nh, dNh

real (wp) :: c1, c2, d1, d2, d3, d4
real (wp) :: ya, yb, fa, fb, dy     !, Rfrc2
real (wp) :: B, L, EE
real (wp) :: Dd1lon,  Dd2lon,  Dd3lon,  Dd4lon
real(wp)  :: dfa1, dfa3, dfb1, dfb3
real (wp) :: Dyalon, Dyblon, Ddylon

!real(wp) :: RfrcLambda2, RfrcBeta2, RfrcHeight2, a1

integer :: i, m, h
!---------------------------------------------------------------------

allocate( k(1:slant%Nnghb) )
allocate( e(0:1,1:slant%Nnghb) )
allocate( N(0:1,1:slant%Nnghb) )

do i=1, slant%Nnghb
   ! Find the correct vertical indices for the next set of columns
   Hcol => Hgeo(i,p,:)
   k(i) = LayerSearch (Hcol, RefPtEll(3))
end do

! Convert longitude
! All transformations require longitudes 0° <= lon <= 360°
! but the model provides longitudes between  -180° <= lon <= 180°
! Convert RefPtEll to degrees (as in col%dlon,dlat) and to -180°<=lon<=180°
RefPtEllM(1:2) = RefPtEll(1:2) * rad2deg
RefPtEllM(3)   = RefPtEll(3)
if (RefPtEllM(1) > 180.0_wp) RefPtEllM(1) = RefPtEllM(1)-360.0_wp

! Compute the partial pressure of water vapour e and the refractivity
! on the corners of the grid cell, i.e 8 corners for COSMO
do i=1, slant%Nnghb
   ! neighbored columns ...

!!$       write(*,*) 'Pt ', p, ' Col. ', i, ' k = ', k(i), ' : ', &
!!$            Hgeo(i,p,k(i)), ' < ', slant%SlantPointsEll(p,3), &
!!$            ' < ', Hgeo(i,p,k(i)-1)
!!$       write(*,*) 'Druck ', col(i,p)%p(k(i)), ' > ', col(i,p)%p(k(i)-1)

   do m=0, 1
      ! upper (m=0) and lower (m=1) layer

      ! Vertical index:
      ! m=0 -  h = k(i)-1 - upper layer
      ! m=1 -  h = k(i)   - lower layer
      h = k(i) - 1 + m

      ! Partial pressure of water vapour e
      e(m,i) = (col(i,p)%p(h) * col(i,p)%q(h)) / (RDRD +   &
                col(i,p)%q(h)*EMRDRD)

      ! Refractivity N
      N(m,i) = NWein(col(i,p)%p(h), col(i,p)%t(h), e(m,i))

   end do
end do

do i=1, slant%Nnghb
   ! vertical interpolation of the refractivity to the required height
   if (abs(N(1,i)) .lt. 1.0E-10_wp .and. abs(N(0,i)) .lt. 1.0E-10_wp ) then
      ! N1 = N2 = 0  => N(z) = 0
      Nh(i) = 0.0_wp
   else
      ! interpolate
      if (abs(N(1,i)) .lt. 1.0E-10_wp) then
         ! N1 = 0, N2 <> 0
         N(1,i) = N(0,i)/100.0_wp
      end if
      if (abs(N(0,i)) .lt. 1.0E-10_wp) then
         ! N2 = 0, N1 <> 0
         N(0,i) = N(1,i)/100.0_wp
      end if
      Nh(i) = N(0,i) * (N(1,i)/N(0,i))**                     &
           ((Hgeo(i,p,k(i)-1)-RefPtEllM(3))/(Hgeo(i,p,k(i)-1)-Hgeo(i,p,k(i))))
      B = ( Hgeo(i,p,k(i)-1) - RefPtEllM(3) )   /             &
          ( Hgeo(i,p,k(i)-1)-Hgeo(i,p,k(i)) )
      L = log( N(1,i)/N(0,i) )
      EE = exp( B * L )
      dNh(i) = -L/( Hgeo(i,p,k(i)-1) - Hgeo(i,p,k(i)) ) * N(0,i) * EE
   end if
   !write(*,*) 'Nh(i) = ', i, Nh(i)
end do

! Bilinear horizontal interpolation
lon(1:slant%Nnghb) = col(:,p)%dlon
lat(1:slant%Nnghb) = col(:,p)%dlat
if (abs(lon(1)) > 175.0) then
   ! avoid problems near plus/minus 180°
   ! convert longitude between +-180° to 0°-360°
   do i=1, slant%Nnghb
      if (lon(i) < 0.0) lon(i) = 360.0_wp - lon(i)
   end do
   if (RefPtEllM(1) < 0.0) RefPtEllM(1) = 360.0_wp - lon(i)
end if
if (slant%Nnghb .eq. 3) then
   ! ICON/GME - interpolate between 3 points
   if (abs(lon(1)-lon(2)) < 1.0E-6_wp) then
      ! change 2 and 3 to avoid singularity, copy 3 to 4
      !write(*,*) 'lon(1) = lon(2)'
      lon(4) = lon(3)
      lat(4) = lat(3)
      Nh(4)  = Nh(3)
      dNh(4) = dNh(3)
      lon(3) = lon(2)
      lat(3) = lat(2)
      Nh(3)  = Nh(2)
      dNh(3) = dNh(2)
      lon(2) = lon(4)
      lat(2) = lat(4)
      Nh(2)  = Nh(4)
      dNh(2) = dNh(4)
   else if (abs(lon(2)-lon(3)) < 1.0E-6_wp) then
      ! copy 1 to 4 to avoid singularity
      !write(*,*) 'lon(2) = lon(3)'
      lon(4) = lon(1)
      lat(4) = lat(1)
      Nh(4)  = Nh(1)
      dNh(4) = dNh(1)
   else
      lon(4) = lon(2)
      lat(4) = lat(2)
      Nh(4)  = Nh(2)
      dNh(4) = dNh(2)
   end if
end if

! c1 = x2 - x1
c1 = lon(2) - lon(1)

! c2 = x4 - x3
c2 = lon(4) - lon(3)

! d1 = (x - x1)/(x2-x1) =  (x - x1)/c1
d1 = (RefPtEllM(1) - lon(1)) / c1
Dd1lon = rad2deg / c1

! d2 = (x2 - x)/(x2-x1) =  (x2 - x)/c1
d2 = (lon(2) - RefPtEllM(1)) / c1
Dd2lon = -rad2deg / c1

! d3 = (x - x3)/(x4-x3) =  (x - x3)/c2
d3 = (RefPtEllM(1) - lon(3)) / c2
Dd3lon = rad2deg / c2

! d4 = (x4 - x)/(x4-x3) =  (x4 - x)/c2
d4 = (lon(4) - RefPtEllM(1)) / c2
Dd4lon = -rad2deg / c2

!fa = d2*node(1,3) + d1*node(2,3)
fa = d2 * Nh(1) + d1 * Nh(2)
dfa1 = Dd2lon * Nh(1) + Dd1lon * Nh(2)   ! d fa / d lon
dfa3 = d2 * dNh(1) + d1 * dNh(2)    ! d fa / d height,  d fa / d lat = 0

!fb = d4*node(3,3) + d3*node(4,3)
fb = d4*Nh(3) + d3*Nh(4)
dfb1 = Dd4lon * Nh(3) + Dd3lon * Nh(4)   ! d fb / d lon
dfb3 = d4 * dNh(3) + d3 * dNh(4)    ! d fb / d height, d fb / d lat = 0

!ya = d2*node(1,2) + d1*node(2,2)
ya = d2*lat(1) + d1*lat(2)
Dyalon = Dd2lon*lat(1) + Dd1lon*lat(2)

!yb = d4*node(3,2) + d3*node(4,2)
yb = d4*lat(3) + d3*lat(4)
Dyblon = Dd4lon*lat(3) + Dd3lon*lat(4)

dy = yb - ya
Ddylon = Dyblon - Dyalon

if (abs(dy) .lt. 1.0E-8_wp) then
   ! slant%Nnghb = 3 : reference point is very close to point 2
   Rfrc = Nh(2)
   RfrcHeight = dNh(2)
   RfrcLambda = 0.0_wp
   RfrcBeta   = 0.0_wp
else
   ! No singularity, interpolate ...
   Rfrc = fa*(yb- RefPtEllM(2))/dy + fb*(RefPtEllM(2)-ya)/dy
   RfrcLambda = Dfa1*(yb - RefPtEllM(2))/dy + fa*Dyblon/dy -   &
                fa*Ddylon*(yb - RefPtEllM(2))/dy**2  +       &
                Dfb1*(RefPtEllM(2)-ya)/dy - fb*Dyalon/dy -   &
                fb*Ddylon*(RefPtEllM(2)-ya)/dy**2
   RfrcBeta   = rad2deg*(fb-fa)/dy
   RfrcHeight = dfa3*(yb-RefPtEllM(2))/dy + dfb3*(RefPtEllM(2)-ya)/dy
end if

deallocate( k )
deallocate( e )
deallocate( N )

End subroutine RefracModelGrad


!---------------------------------------------------------------------
! subroutine RefracModelGrad2
!---------------------------------------------------------------------
!>
!> @brief Second derivative of the refractivity N with respect to the
!>        ellipsoidal coordinates lambda, beta, height.
!>
!> <b> call RefracModelGrad2 (slant, col, Hgeo, p, RefPtEll, Rfrc,   &    \n
!>                           RfrcLambda, RfrcBeta, RfrcHeight, Grad2Rfrc) </b>
!>
!> Estimate the refractivity at some reference point by interpolating
!> between the model grid nodes and compute the first and second
!> derivatives with respect to longitude, latitude and height.
!>
!> Original routine without derivatives: @ref RefracModel \n
!> Routine with first derivatives:       @ref RefracModelGrad \n
!> \n
!> Computes the refractivity at the position "RefPtEll" using temperature
!> pressure and humidity from the weather model (provided by "col") and the
!> first and second derivatives, i.e. additional calls to ::RefracModel
!> or ::RefracModelGrad are usually  not necessary.
!> The interpolation is carried out in three steps:
!>
!> 1) Compute the refractivity N at the grid nodes of the cell containing
!>    "RefPtEll".\n
!> 2) Vertical interpolation to the required height (RefPtEll(3)).\n
!> 3) Bilinear horizontal interpolation between the nodes.
!>
!> The horizontal bilinear interpolation works with 3 (ICON) or 4 (COSMO)
!> grid columns. In case of 3 columns one node is copied to the "missing"
!> fourth node and the interpolation is done in exactly the same way as with
!> four nodes. To avoid singularities some checks are made and the 3
!> nodes are reordered accordingly. \n
!> The derivatives of the refractivity N with respect to the
!> ellipsoidal coordinates longitude (lambda), latitude (beta) and height
!> are provided: \n
!> \f$ RfrcLambda = \frac{\partial N}{\partial \lambda} \f$,
!> \f$ RfrcBeta = \frac{\partial N}{\partial \beta} \f$ and
!> \f$ RfrcHeight = \frac{\partial N}{\partial h} \f$.
!>
!> The second derivatives of N with respect to the ellipsoidal coordinates
!> are stored in the matrix Grad2Rfrc: \n
!> \verbatim
!>             Grad2Rfrc(j,i) :  i - first varaiable
!>                               j - second variable
!>
!> Grad2Rfrc(j,i) =
!>
!>           ( 11 12 13 )   ( N_lambda_lambda N_beta_lambda N_height_lambda )
!>           ( 21 22 23 ) = ( N_lambda_beta   N_beta_beta   N_height_beta   )
!>           ( 31 32 33 )   ( N_lambda_height N_beta_height N_height_height )
!> \endverbatim
!>
!> ModuleTest.f90: Test ::RefracModelGrad2 with TestRefracModel = .true. \n
!>                        Derivatives are tested with finite differences.
!>
!> @param[in] slant  structure containing all required information about
!>                   one slant
!>                   (output will also be written to this structure)
!> @param[in] col    grid columns of the model field
!> @param[in] Hgeo   height above ellipsoid computed for all columns "col" \n
!>                   Hgeo(i,j,:) = col(i,j)%gpm(:) + col(i,j)%geoid
!> @param[in] p      index of the column to be used
!> @param[in] RefPtEll  ellipsoidal coordinates of the reference point,
!>                      the refractivity will be computed for this point \n
!>                      RefPtEll(1) - longitude, rad \n
!>                      RefPtEll(2) - latitude, rad \n
!>                      RefPtEll(3) - height above allipsoid, m
!>
!> @param[out] Rfrc  refractivity N, N = 10^6*(n-1), n - refraction index
!> @param[out] RfrcLambda derivative \f$ \frac{\partial N}{\partial \lambda} \f$
!> @param[out] RfrcBeta   derivative \f$ \frac{\partial N}{\partial \beta} \f$
!> @param[out] RfrcHeight derivative \f$ \frac{\partial N}{\partial h} \f$
!> @param[out] Grad2Rfrc  matrix with second derivatives
!> @param[out] err        an error occured,
!>                        stop iteration and reject observation \n
!>                        err = 1  -  negative refractivity
!
!---------------------------------------------------------------------
! Record of revisions
! Date        Programmer   Description of change
! ====        ==========   =====================
! 16.12.2013  M. Bender    new
! 10.03.2014  M. Bender    longitude changed to -Pi <= lon <= +Pi
! 10.07.2015  M. Bender    major changes to avoid problems with two
!                          points having the same longitude.
!---------------------------------------------------------------------
subroutine RefracModelGrad2 (slant, col, Hgeo, p, RefPtEll, Rfrc,   &
                            RfrcLambda, RfrcBeta, RfrcHeight,       &
                            Grad2Rfrc, err)

implicit none

! List of calling arguments:
type (SlantData) ,intent(in)    :: slant        ! slant meta data
type (p_column)  ,intent(in)    :: col (:,:)    ! model columns     (input)
real(wp), dimension(:,:,:), pointer :: Hgeo     ! geoid coordinates (input)
integer, intent(in)                 :: p        ! index of pont in "SlantPoints"
real(wp), dimension(1:3), intent(in) :: RefPtEll ! compute N here
real(wp),intent(out)                :: Rfrc     ! refractivity (output)
real(wp),intent(out)                :: RfrcLambda, RfrcBeta, RfrcHeight
real(wp), dimension(1:3,1:3), intent(out) :: Grad2Rfrc
integer , intent(out)                     :: err

! List of local variables:
real(wp), dimension(:), pointer :: Hcol => Null()
real(wp), dimension(1:3)        :: RefPtEllM

! vertical grid index of neighbored columns:
integer, dimension(:), pointer  :: k => Null()
real(wp), dimension(:,:), pointer :: e => Null()
real(wp), dimension(:,:), pointer :: N => Null()
!real(wp), dimension(1:3) :: RefPt_tl
real(wp), dimension(1:4) ::  lon, lat, Nh, dNh, d2Nh

real (wp) :: c1, c2, d1, d2, d3, d4
real (wp) :: ya, yb, fa, fb, dy     !, Rfrc2
real (wp) :: B, L, EE, dB, dEE

real (wp) :: Dd1lon,  Dd2lon,  Dd3lon,  Dd4lon
real(wp)  :: dfa1, dfa3, dfb1, dfb3
real (wp) :: d2fa13, d2fa31, d2fa33, d2fb13, d2fb31, d2fb33
real (wp) :: Dyalon, Dyblon, Ddylon

!real(wp) :: RfrcLambda2, RfrcBeta2, RfrcHeight2, a1, a2

character (len=4) :: SName

integer :: i, m, h
!---------------------------------------------------------------------

err = 0

allocate( k(1:slant%Nnghb) )
allocate( e(0:1,1:slant%Nnghb) )
allocate( N(0:1,1:slant%Nnghb) )

do i=1, slant%Nnghb
   ! Find the correct vertical indices for the next set of columns
   Hcol => Hgeo(i,p,:)
   k(i) = LayerSearch (Hcol, RefPtEll(3))
end do

! Convert longitude
! All transformations require longitudes 0° <= lon <= 360°
! but the model provides longitudes between  -180° <= lon <= 180°
! Convert RefPtEll to degrees (as in col%dlon,dlat) and to -180°<=lon<=180°
RefPtEllM(1:2) = RefPtEll(1:2) * rad2deg
RefPtEllM(3)   = RefPtEll(3)
if (RefPtEllM(1) > 180.0_wp) RefPtEllM(1) = RefPtEllM(1)-360.0_wp

! Compute the partial pressure of water vapour e and the refractivity
! on the corners of the grid cell, i.e 8 corners for COSMO
do i=1, slant%Nnghb
   ! neighbored columns ...

!!$       write(*,*) 'Pt ', p, ' Col. ', i, ' k = ', k(i), ' : ', &
!!$            Hgeo(i,p,k(i)), ' < ', slant%SlantPointsEll(p,3), &
!!$            ' < ', Hgeo(i,p,k(i)-1)
!!$       write(*,*) 'Druck ', col(i,p)%p(k(i)), ' > ', col(i,p)%p(k(i)-1)

   do m=0, 1
      ! upper (m=0) and lower (m=1) layer

      ! Vertical index:
      ! m=0 -  h = k(i)-1 - upper layer
      ! m=1 -  h = k(i)   - lower layer
      h = k(i) - 1 + m

      ! Partial pressure of water vapour e
      e(m,i) = (col(i,p)%p(h) * col(i,p)%q(h)) / (RDRD +   &
                col(i,p)%q(h)*EMRDRD)

      ! Refractivity N
      N(m,i) = NWein(col(i,p)%p(h), col(i,p)%t(h), e(m,i))

   end do
end do

if (.false.) then
!if ( any(N(:,:) < 0.0) ) then
   write(*,*) 'RefracModelGrad2: Refractivity in columns'
   write(*,*) 'Index p = ', p
   write(*,*) 'Hoehe RefPt: ', RefPtEll(3)
   do i=1, slant%Nnghb
      write(*,*) i, k(i), col(i,p)%gpm(k(i)), N(1,i), ' -- ', &
                 i, k(i)-1, col(i,p)%gpm(k(i)-1),N(0,i)
   end do
end if

do i=1, slant%Nnghb
   ! vertical interpolation of the refractivity to the required height
   if (abs(N(1,i)) .lt. 1.0E-10_wp .and. abs(N(0,i)) .lt. 1.0E-10_wp ) then
      ! N1 = N2 = 0  => N(z) = 0
      Nh(i) = 0.0_wp
   else
      ! interpolate
      if (abs(N(1,i)) .lt. 1.0E-10_wp) then
         ! N1 = 0, N2 <> 0
         N(1,i) = N(0,i)/100.0_wp
      end if
      if (abs(N(0,i)) .lt. 1.0E-10_wp) then
         ! N2 = 0, N1 <> 0
         N(0,i) = N(1,i)/100.0_wp
      end if
      Nh(i) = N(0,i) * (N(1,i)/N(0,i))**                     &
           ((Hgeo(i,p,k(i)-1)-RefPtEllM(3))/(Hgeo(i,p,k(i)-1)-Hgeo(i,p,k(i))))
      B = ( Hgeo(i,p,k(i)-1) - RefPtEllM(3) )   /             &
          ( Hgeo(i,p,k(i)-1)-Hgeo(i,p,k(i)) )
      L = log( N(1,i)/N(0,i) )
      EE = exp( B * L )
      dNh(i) = -L/( Hgeo(i,p,k(i)-1) - Hgeo(i,p,k(i)) ) * N(0,i) * EE
      ! second derivatives:
      dB = -1.0_wp / ( Hgeo(i,p,k(i)-1)-Hgeo(i,p,k(i)) )
      dEE = L * dB * EE
      d2Nh(i) = dNh(i) * dEE/EE
   end if
   !write(*,*) 'Nh(i) = ', i, Nh(i)
end do

!write(*,*) 'RefracModelGrad2: Refractivity vert. interpol. in cols'
!write(*,*) Nh

! Bilinear horizontal interpolation
lon(1:slant%Nnghb) = col(:,p)%dlon
lat(1:slant%Nnghb) = col(:,p)%dlat
if (abs(lon(1)) > 175.0) then
   ! avoid problems near plus/minus 180°
   ! convert longitude between +-180° to 0°-360°
   do i=1, slant%Nnghb
      if (lon(i) < 0.0) lon(i) = 360.0_wp - lon(i)
   end do
   if (RefPtEllM(1) < 0.0) RefPtEllM(1) = 360.0_wp - lon(i)
end if
if (slant%Nnghb .eq. 3) then
   ! ICON/GME - interpolate between 3 points
   if (abs(lon(1)-lon(2)) < 1.0E-6_wp) then
      ! change 2 and 3 to avoid singularity, copy 3 to 4
      !write(*,*) 'lon(1) = lon(2)'
      lon(4) = lon(3)
      lat(4) = lat(3)
      Nh(4)  = Nh(3)
      dNh(4) = dNh(3)
      d2Nh(4) = d2Nh(3)
      lon(3) = lon(2)
      lat(3) = lat(2)
      Nh(3)  = Nh(2)
      dNh(3) = dNh(2)
      d2Nh(3) = d2Nh(2)
      lon(2) = lon(4)
      lat(2) = lat(4)
      Nh(2)  = Nh(4)
      dNh(2) = dNh(4)
      d2Nh(2) = d2Nh(4)
   else if (abs(lon(2)-lon(3)) < 1.0E-6_wp) then
      ! copy 1 to 4 to avoid singularity
      !write(*,*) 'lon(2) = lon(3)'
      lon(4) = lon(1)
      lat(4) = lat(1)
      Nh(4)  = Nh(1)
      dNh(4) = dNh(1)
      d2Nh(4) = d2Nh(1)
   else
      lon(4) = lon(2)
      lat(4) = lat(2)
      Nh(4)  = Nh(2)
      dNh(4) = dNh(2)
      d2Nh(4) = d2Nh(2)
   end if
end if

! c1 = x2 - x1
c1 = lon(2) - lon(1)

! c2 = x4 - x3
c2 = lon(4) - lon(3)

! d1 = (x - x1)/(x2-x1) =  (x - x1)/c1
d1 = (RefPtEllM(1) - lon(1)) / c1
Dd1lon = rad2deg / c1

! d2 = (x2 - x)/(x2-x1) =  (x2 - x)/c1
d2 = (lon(2) - RefPtEllM(1)) / c1
Dd2lon = -rad2deg / c1

! d3 = (x - x3)/(x4-x3) =  (x - x3)/c2
d3 = (RefPtEllM(1) - lon(3)) / c2
Dd3lon = rad2deg / c2

! d4 = (x4 - x)/(x4-x3) =  (x4 - x)/c2
d4 = (lon(4) - RefPtEllM(1)) / c2
Dd4lon = -rad2deg / c2

!fa = d2*node(1,3) + d1*node(2,3)
fa = d2 * Nh(1) + d1 * Nh(2)
! first derivatives of fa, d fa / d lat = 0
dfa1 = Dd2lon * Nh(1) + Dd1lon * Nh(2)   ! d fa / d lon
dfa3 = d2 * dNh(1) + d1 * dNh(2)    ! d fa / d height,  d fa / d lat = 0
! second derivatives of fa
! d dfa1 / d lon = 0,  d dfa1 / d lat = 0
d2fa13 = Dd2lon * dNh(1) + Dd1lon * dNh(2)  ! d dfa1 / d height
! d dfa3 / d lat = 0
d2fa31 = d2fa13                        ! d dfa3 / d lon
d2fa33 = d2 * d2Nh(1) + d1 * d2Nh(2)   ! d dfa3 / d height

!fb = d4*node(3,3) + d3*node(4,3)
fb = d4*Nh(3) + d3*Nh(4)
dfb1 = Dd4lon * Nh(3) + Dd3lon * Nh(4)   ! d fb / d lon
dfb3 = d4 * dNh(3) + d3 * dNh(4)    ! d fb / d height, d fb / d lat = 0
d2fb13 =  Dd4lon * dNh(3) + Dd3lon * dNh(4)
d2fb31 = d2fb13
d2fb33 = d4 * d2Nh(3) + d3 * d2Nh(4)   ! d dfb3 / d height

!ya = d2*node(1,2) + d1*node(2,2)
ya = d2*lat(1) + d1*lat(2)
Dyalon = Dd2lon*lat(1) + Dd1lon*lat(2)

!yb = d4*node(3,2) + d3*node(4,2)
yb = d4*lat(3) + d3*lat(4)
Dyblon = Dd4lon*lat(3) + Dd3lon*lat(4)

dy = yb - ya
Ddylon = Dyblon - Dyalon
! second derivatives are all 0

if (abs(dy) .lt. 1.0E-8_wp) then
   ! slant%Nnghb = 3 : reference point is very close to point 2
   Rfrc = Nh(2)
   RfrcHeight = dNh(2)
   RfrcLambda = 0.0_wp
   RfrcBeta   = 0.0_wp
   Grad2Rfrc(3,3) = d2Nh(1)
else
   ! No singularity, interpolate ...

   ! Refractivity N interpolated at the reference point (lon, lat alt)
   Rfrc = fa*(yb- RefPtEllM(2))/dy + fb*(RefPtEllM(2)-ya)/dy

   ! First derivative of N with respect to the coordinates of
   ! the reference point (lon, lat alt) = (lambda, beta, height)
   RfrcLambda = Dfa1*(yb - RefPtEllM(2))/dy + fa*Dyblon/dy -   &
                fa*Ddylon*(yb - RefPtEllM(2))/dy**2  +       &
                Dfb1*(RefPtEllM(2)-ya)/dy - fb*Dyalon/dy -   &
                fb*Ddylon*(RefPtEllM(2)-ya)/dy**2
   RfrcBeta   = rad2deg*(fb-fa)/dy
   RfrcHeight = dfa3*(yb-RefPtEllM(2))/dy + dfb3*(RefPtEllM(2)-ya)/dy

   ! Second derivative of N with respect to the coordinates of
   ! the reference point
   ! N_height_height
   Grad2Rfrc(3,3) = d2fa33*(yb-RefPtEllM(2))/dy + d2fb33*(RefPtEllM(2)-ya)/dy
   ! N_height_beta
   Grad2Rfrc(2,3) = rad2deg*(dfb3 - dfa3) / dy
   ! N_height_lambda
   Grad2Rfrc(1,3) = ( dy*(dfa3*Dyblon + d2fa31*yb - RefPtEllM(2)*d2fa31) -  &
                      dfa3*Ddylon*(yb-RefPtEllM(2)) +                        &
                      dy*(RefPtEllM(2)*d2fb31 - dfb3*Dyalon - d2fb31*ya ) - &
                      dfb3*Ddylon*(RefPtEllM(2)-ya)  ) / dy**2
   ! N_beta_height
   Grad2Rfrc(3,2) = Grad2Rfrc(2,3)
   ! N_beta_beta
   Grad2Rfrc(2,2) = 0.0_wp
   ! N_beta_lambda
   Grad2Rfrc(1,2) = ( dy*(dfb1-dfa1) - Ddylon*(fb-fa) ) / dy**2
   ! N_lambda_height
   Grad2Rfrc(3,1) =   Grad2Rfrc(1,3)
   ! N_lambda_beta
   Grad2Rfrc(2,1) = Grad2Rfrc(1,2)
   ! N_lambda_lambda
                    ! 1) d/dlambda ( Dfa1*(yb - RefPtEllM(2))/dy )
   Grad2Rfrc(1,1) = dfa1*(dy*Dyblon - (yb-RefPtEllM(2))*Ddylon) / dy**2 + &
                    ! 2) d/dlambda (  fa*Dyblon/dy )
                    Dyblon*(dy*dfa1-fa*Ddylon)/dy**2 - &
                    ! 3) d/dlambda ( fa*Ddylon*(yb - RefPtEllM(2))/dy**2 )
                    Ddylon*(yb-RefPtEllM(2)) *                      &
                           (dfa1*dy**2 - 2*dy*fa*Ddylon) / dy**4 - &
                           (fa*Ddylon*Dyblon) / dy**2    - &
                    ! 4) d/dlambda (  Dfb1*(RefPtEllM(2)-ya)/dy )
                    Dfb1*(dy*Dyalon + (RefPtEllM(2)-ya)*Ddylon)/dy**2  - &
                    ! 5) d/dlambda ( fb*Dyalon/dy )
                    Dyalon*(dy*dfb1 - fb*Ddylon) / dy**2  - &
                    ! 6) d/dlambda ( fb*Ddylon*(RefPtEllM(2)-ya)/dy**2 )
                    Ddylon*(RefPtEllM(2)-ya)*  &
                        (dy**2*dfb1 - 2.0_wp*fb*dy*Ddylon) / dy**4 + &
                        fb*Ddylon*Dyalon / dy**2

end if

!write(*,*) 'RefracModelGrad2: Refractivity horz. interpol. = ', Rfrc

if (Rfrc .lt. 0.0_wp) then
   ! This is a severe error and indicates that Newton's method does
   ! not converge.
   ! => Stop iteration before program crashes and reject observation.
   err = 1

   ! Print some debug messages
   if (associated(slant%Station)) then
      SName = slant%Station%SName
   else
      SName = '----'
   end if
   write(*,*) 'RefracModelGrad2> BiLin2D negative Refrakt. ', Rfrc, char(10), &
        'slant%Nnghb = ', slant%Nnghb,                             char(10), &
        'station short name = ', SName,                            char(10), &
        'azimuth, elevation = ', slant%Obs%azimuth*rad2deg,                  &
                                 slant%Obs%elevation*rad2deg,                &
                                                                   char(10), &
        'RefPt(1), Laenge   = ', RefPtEll(1), RefPtEll(1)*rad2deg, char(10), &
        'RefPt(2), Breite   = ', RefPtEll(2), RefPtEll(2)*rad2deg, char(10), &
        'RefPt(3), Hoehe    = ', RefPtEll(3),                      char(10), &
        'RefPtM(1), Laenge  = ', RefPtEllM(1),                     char(10), &
        'RefPtM(2), Breite  = ', RefPtEllM(2),                     char(10), &
        'RefPtM(3), Hoehe   = ', RefPtEllM(3),                     char(10), &
        'Eckpunkte,  Laenge = ', col(:,p)%dlon,                    char(10), &
        'Eckpunkte,  Breite = ', col(:,p)%dlat,                    char(10), &
!        'Eckpunkte,  Höhe   = ', col(:,p)%gpm,                     char(10), &
        'Nh = ', Nh,  char(10), &
        'c1 = ', c1,  char(10), &
        'c2 = ', c2,  char(10), &
        'd1 = ', d1,  char(10), &
        'd2 = ', d2,  char(10), &
        'd3 = ', d3,  char(10), &
        'd4 = ', d4,  char(10), &
        'fa = ', fa,  char(10), &
        'fb = ', fb,  char(10), &
        'ya = ', ya,  char(10), &
        'yb = ', yb,  char(10), &
        'dy = ', dy
end if

deallocate( k )
deallocate( e )
deallocate( N )

End subroutine RefracModelGrad2


!---------------------------------------------------------------------
! subroutine RefracModelOld
!---------------------------------------------------------------------
!
! Start of RoboDoc comments >>>>
!****S* Fermat/RefracModel
!
! Name
! RefracModel
!
! Call
! call RefracModel (slant, col, Hgeo, p, RefPtEll, Rfrc)
!
! Purpose
! Compute the refractivity at the position "RefPtEll" using temperature
! pressure and humidity from the weather model (provided by "col").
! The interpolation is carried out in three steps:
!
! 1) Compute the refractivity N at the grid nodes of the cell containing
!    "RefPtEll".
! 2) Vertical interpolation to the required height (RefPtEll(3)).
! 3) Bilinear horizontal interpolation between these four points.
!
! Input
! Slant    - structure containing all required information about one slant
!            (output will also be written to this structure)
! col      - grid columns of the model field
! Hgeo     - height above ellipsoid computed for all columns "col"
!            Hgeo(i,j,:) = col(i,j)%gpm(:) + col(i,j)%geoid
! p        - index of the column to be used
! RefPtEll - ellipsoidal coordinates of the reference point, the refractivity
!            will be computed for this point
!            RefPtEll(1) - longitude, rad
!            RefPtEll(2) - latitude, rad
!            RefPtEll(3) - height above allipsoid, m
!
! Output
! Rfrc - refractivity N, N = 10^6*(n-1), n - refraction index
!
! External References
! None
!
! Record of revisions
! Date        Programmer   Description of change
! ====        ==========   =====================
! 10.09.2013  M. Bender    new
! 10.03.2014  M. Bender    longitude changed to -Pi <= lon <= +Pi
!---------------------------------------------------------------------
! Source
!---------------------------------------------------------------------
subroutine RefracModelOld (slant, col, Hgeo, p, RefPtEll, Rfrc)

implicit none

! List of calling arguments:
type (SlantData) ,intent(in)    :: slant      ! slant meta data
type (p_column)  ,intent(in)    :: col (:,:)  ! model columns     (input)
real(wp), dimension(:,:,:), pointer :: Hgeo   ! geoid coordinates (input)
integer, intent(in)                 :: p      ! index of point in "SlantPoints"
real(wp), dimension(1:3), intent(in) :: RefPtEll ! compute N here
real(wp)         ,intent(out)   :: Rfrc         ! refractivity (output)

! List of local variables:
real(wp), dimension(:), pointer :: Hcol => Null()
real(wp), dimension(1:3)        :: RefPtEllM

! vertical grid index of neighbored columns:
integer, dimension(:), pointer  :: k => Null()
!integer, dimension(:), pointer  :: k2 => Null()
real(wp), dimension(:,:), pointer :: e => Null()
real(wp), dimension(:,:), pointer :: N => Null()
real (wp), dimension(1:4,1:2,1:4)  :: node
real(wp), dimension(1:4) :: Nh

real (wp) :: c1, c2, d1, d2, d3, d4
real (wp) :: ya, yb, fa, fb, dy     !, Rfrc2

integer :: i, m, h
!---------------------------------------------------------------------

allocate( k(1:slant%Nnghb) )
! allocate( k2(1:slant%Nnghb) )
allocate( e(0:1,1:slant%Nnghb) )
allocate( N(0:1,1:slant%Nnghb) )

do i=1, slant%Nnghb
   ! Find the correct vertical indices for the next set of columns
   Hcol => Hgeo(i,p,:)
   !k(i) = LayerSearch (Hcol, slant%SlantPointsEll(p,3))
   !k(i) = LayerSearch (real(Hcol), real(RefPtEll(3)))
   k(i) = LayerSearch (Hcol, RefPtEll(3))
end do

! Convert longitude
! All transformations require longitudes 0° <= lon <= 360°
! but the model provides longitudes between  -180° <= lon <= 180°
! ????? Auch gleich in Grad umrechnen? Das spart die staendige Umrechnung ... ?????
RefPtEllM = RefPtEll
if (RefPtEll(1) .gt. 180.0*deg2rad) RefPtEllM(1) = RefPtEll(1)-360.0*deg2rad

! Compute the partial pressure of water vapour e and the refractivity
! on the corners of the grid cell, i.e 8 corners for COSMO
do i=1, slant%Nnghb
   ! neighbored columns ...

!!$       write(*,*) 'Pt ', p, ' Col. ', i, ' k = ', k(i), ' : ', &
!!$            Hgeo(i,p,k(i)), ' < ', slant%SlantPointsEll(p,3), &
!!$            ' < ', Hgeo(i,p,k(i)-1)
!!$       write(*,*) 'Druck ', col(i,p)%p(k(i)), ' > ', col(i,p)%p(k(i)-1)

   do m=0, 1
      ! upper (m=0) and lower (m=1) layer

      ! Vertical index:
      ! m=0 -  h = k(i)-1 - upper layer
      ! m=1 -  h = k(i)   - lower layer
      h = k(i) - 1 + m

      ! Partial pressure of water vapour e
      e(m,i) = (col(i,p)%p(h) * col(i,p)%q(h)) / (RDRD +   &
                col(i,p)%q(h)*EMRDRD)

      ! Refractivity N
      N(m,i) = NWein(col(i,p)%p(h), col(i,p)%t(h), e(m,i))

      ! Interpolate the refractivity at the given point using the
      ! refractivities at the corners of the grid cell
      ! node(:,1,:) - lower level - k+1
      ! node(:,2,:) - upper level - k
      node(i,2-m,1) = col(i,p)%dlon*deg2rad
      node(i,2-m,2) = col(i,p)%dlat*deg2rad
      node(i,2-m,3) = Hgeo(i,p,h)
      node(i,2-m,4) = N(m,i)

   end do
end do

! Interpolation using function "BiLinearExp3D"
if (slant%Nnghb .eq. 3) then
   ! GME - interpolate between 3 points
   node(4,:,:) = node(2,:,:)
end if

!Rfrc = BiLinearExp3D(slant%SlantPointsEll(p,:), node, slant%Nnghb)
!Rfrc = BiLinearExp3D(RefPtEll(:), node, slant%Nnghb)

!write(*,*) 'N = ', N(0,1), N(1,1), ' ipol = ', RefracProfile(p,2)

! ExpInt1D
do i=1, slant%Nnghb
   ! vertical interpolation of the refractivity to the required height
   if (abs(N(1,i)) .lt. 1.0E-10_wp .and. abs(N(0,i)) .lt. 1.0E-10_wp ) then
      ! N1 = N2 = 0  => N(z) = 0
      Nh(i) = 0.0_wp
   else
      ! interpolate
      if (abs(N(1,i)) .lt. 1.0E-10_wp) then
         ! N1 = 0, N2 <> 0
         N(1,i) = N(0,i)/100.0_wp
      end if
      if (abs(N(0,i)) .lt. 1.0E-10_wp) then
         ! N2 = 0, N1 <> 0
         N(0,i) = N(1,i)/100.0_wp
      end if
      Nh(i) = N(0,i) * (N(1,i)/N(0,i))**                     &
           ((Hgeo(i,p,k(i)-1)-RefPtEllM(3))/(Hgeo(i,p,k(i)-1)-Hgeo(i,p,k(i))))
   end if
end do

! BiLinear2D
! Kopieren, Sonderbehandpung nur fuer Punkte innerhalb des Dreiecks
! => Eckpunkt kopieren

! c1 = x2 - x1
!c1 = node(2,1) - node(1,1)
c1 = (col(2,p)%dlon - col(1,p)%dlon) * deg2rad

! c2 = x4 - x3
!c2 = node(4,1) - node(3,1)
if (slant%Nnghb .eq. 4) then
   c2 = (col(4,p)%dlon - col(3,p)%dlon) * deg2rad
else
   ! slant%Nnghb = 3 : replace node 4 with node 2
   c2 = (col(2,p)%dlon - col(3,p)%dlon) * deg2rad
end if

! d1 = (x - x1)/(x2-x1) =  (x - x1)/c1
!d1 = (RefPt(1) - node(1,1))/c1
d1 = (RefPtEllM(1) - col(1,p)%dlon*deg2rad) / c1

! d2 = (x2 - x)/(x2-x1) =  (x2 - x)/c1
!d2 = (node(2,1) - RefPt(1))/c1
d2 = (col(2,p)%dlon*deg2rad - RefPtEllM(1)) / c1

! d3 = (x - x3)/(x4-x3) =  (x - x3)/c2
!d3 = (RefPt(1) - node(3,1))/c2
d3 = (RefPtEllM(1) - col(3,p)%dlon*deg2rad) / c2

! d4 = (x4 - x)/(x4-x3) =  (x4 - x)/c2
!d4 = (node(4,1) - RefPt(1)) / c2
if (slant%Nnghb .eq. 4) then
   d4 = (col(4,p)%dlon*deg2rad - RefPtEllM(1)) / c2
else
   ! slant%Nnghb = 3 : replace node 4 with node 2
   d4 = (col(2,p)%dlon*deg2rad - RefPtEllM(1)) / c2
end if

!write(*,*) 'BiLinearDiff2D> c1, c2, d1, d2, d3, d4 : ', c1, c2, d1, d2, d3, d4

!fa = d2*node(1,3) + d1*node(2,3)
!fb = d4*node(3,3) + d3*node(4,3)
!ya = d2*node(1,2) + d1*node(2,2)
!yb = d4*node(3,2) + d3*node(4,2)

!fa = d2*node(1,3) + d1*node(2,3)
fa = d2 * Nh(1) + d1 * Nh(2)

!fb = d4*node(3,3) + d3*node(4,3)
if (slant%Nnghb .eq. 4) then
   fb = d4*Nh(3) + d3*Nh(4)
else
   ! slant%Nnghb = 3 : replace node 4 with node 2
   fb = d4*Nh(3) + d3*Nh(2)
end if

!ya = d2*node(1,2) + d1*node(2,2)
ya = (d2*col(1,p)%dlat + d1*col(2,p)%dlat)*deg2rad

!yb = d4*node(3,2) + d3*node(4,2)
if (slant%Nnghb .eq. 4) then
   yb = (d4*col(3,p)%dlat + d3*col(4,p)%dlat)*deg2rad
else
   ! slant%Nnghb = 3 : replace node 4 with node 2
   yb = (d4*col(3,p)%dlat + d3*col(2,p)%dlat)*deg2rad
end if

dy = yb - ya

if (abs(dy) .lt. 1.0E-8_wp) then
   ! slant%Nnghb = 3 : reference point is very close to point 1
   Rfrc = Nh(1)
else
   ! No singularity, interpolate ...
   Rfrc = fa*(yb- RefPtEllM(2))/dy + fb*(RefPtEllM(2)-ya)/dy
end if

!    write(*,*) 'Ipol: 3dvar = ', Rfrc, ' bilin = ', RefracProfile(p,2), &
!               ' diff = ', Rfrc-RefracProfile(p,2), ' diff% = ', &
!                100.0*(Rfrc-RefracProfile(p,2))/Rfrc

!write(*,*) 'RefracModel> Rfrc - org, neu ', Rfrc, Rfrc2

deallocate( k )
deallocate( e )
deallocate( N )

End subroutine RefracModelOld
!****
! <<<<<<<< End of RoboDoc comments


!---------------------------------------------------------------------
! subroutine RefracModelGradOld
!---------------------------------------------------------------------
!
! Start of RoboDoc comments >>>>
!****S* Fermat/RefracModelGrad
!
! Name
! RefracModelGrad
!
! Call
! call RefracModelGrad (slant, col, Hgeo, p, RefPtEll, Rfrc,   &
!                       RfrcLambda, RfrcBeta, RfrcHeight)
!
! Purpose
! Compute the refractivity at the position "RefPtEll" using temperature
! pressure and humidity from the weather model (provided by "col").
! Compute the derivatives of "Rfrc" with respect to lambda, beta and height.
!
! 1) Compute the refractivity N at the grid nodes of the cell containing
!    "RefPtEll".
! 2) Vertical interpolation to the required height (RefPtEll(3)).
! 3) Bilinear horizontal interpolation between these four points.
! Second derivative of the refractivity N with respect to the
! ellipsoidal coordinates lambda, beta, height.
!
! Test with finite differences:
! .../gpstomo/trunk/Tomo/TomoProg/Fermatvergl/FermatVergeich.f90
!
! RefracModel      => N
! RefracModelGrad  => N, d N / d lambda,  d N / d beta,  d N / d height
! RefracModelGrad2 => N, dN, d2N
!
! RefracModel computes one single scalar function, i.e. the refractivity N
! at a certain point (lambda, beta, height).
! The Jacobian consists of 3 derivatives:
!
! J = (N_lambda, N_beta, N_height)
!   = (RfrcLambda, RfrcBeta, RfrcHeight)
!
! Input
! Slant    - structure containing all required information about one slant
!            (output will also be written to this structure)
! col      - grid columns of the model field
! Hgeo     - height above ellipsoid computed for all columns "col"
!            Hgeo(i,j,:) = col(i,j)%gpm(:) + col(i,j)%geoid
! p        - index of the column to be used
! RefPtEll - ellipsoidal coordinates of the reference point, the refractivity
!            will be computed for this point
!            RefPtEll(1) - longitude, rad
!            RefPtEll(2) - latitude, rad
!            RefPtEll(3) - height above allipsoid, m
!
! Output
! Rfrc       - refractivity N, N = 10^6*(n-1), n - refraction index
! RfrcLambda - d N / d lambda
! RfrcBeta   - d N / d beta
! RfrcHeight - d N / d height
!
! External References
! None
!
! Record of revisions
! Date        Programmer   Description of change
! ====        ==========   =====================
! 16.12.2013  M. Bender    new
!---------------------------------------------------------------------
! Source
!---------------------------------------------------------------------
subroutine RefracModelGradOld (slant, col, Hgeo, p, RefPtEll, Rfrc,   &
                            RfrcLambda, RfrcBeta, RfrcHeight)

implicit none

! List of calling arguments:
type (SlantData) ,intent(in)    :: slant        ! slant meta data
type (p_column)  ,intent(in)    :: col (:,:)    ! model columns     (input)
real(wp), dimension(:,:,:), pointer :: Hgeo     ! geoid coordinates (input)
integer, intent(in)                 :: p        ! index of pont in "SlantPoints"
real(wp), dimension(1:3), intent(in) :: RefPtEll ! compute N here
real(wp),intent(out)                :: Rfrc     ! refractivity (output)
real(wp),intent(out)                :: RfrcLambda, RfrcBeta, RfrcHeight
!real(wp),optional,intent(out)                :: A

! List of local variables:
real(wp), dimension(:), pointer     :: Hcol => Null()

! vertical grid index of neighbored columns:
integer, dimension(:), pointer  :: k => Null()
!integer, dimension(:), pointer  :: k2 => Null()
real(wp), dimension(:,:), pointer :: e => Null()
real(wp), dimension(:,:), pointer :: N => Null()
!real (wp), dimension(1:4,1:2,1:4)  :: node, node_tl
!real (wp), dimension(1:3) :: RefPt_tl
real(wp), dimension(1:4) :: Nh, dNh

real (wp) :: c1, c2, d1, d2, d3, d4
real (wp) :: ya, yb, fa, fb, dy     !, Rfrc2
real (wp) :: B, L, EE
real (wp) :: Dd1lon,  Dd2lon,  Dd3lon,  Dd4lon
real(wp)  :: dfa1, dfa3, dfb1, dfb3
real (wp) :: Dyalon, Dyblon, Ddylon

!real(wp) :: RfrcLambda2, RfrcBeta2, RfrcHeight2, a1

integer :: i, m, h
!---------------------------------------------------------------------

allocate( k(1:slant%Nnghb) )
! allocate( k2(1:slant%Nnghb) )
allocate( e(0:1,1:slant%Nnghb) )
allocate( N(0:1,1:slant%Nnghb) )

do i=1, slant%Nnghb
   ! Find the correct vertical indices for the next set of columns
   Hcol => Hgeo(i,p,:)
   !k(i) = LayerSearch (Hcol, slant%SlantPointsEll(p,3))
   !k(i) = LayerSearch (real(Hcol), real(RefPtEll(3)))
   k(i) = LayerSearch (Hcol, RefPtEll(3))
   !write(*,*) 'H k-1, k, k+1 ', Hcol(k(i)-1),  Hcol(k(i)),  Hcol(k(i)+1), RefPtEll(3)
end do

! Compute the partial pressure of water vapour e and the refractivity
! on the corners of the grid cell, i.e 8 corners for COSMO
do i=1, slant%Nnghb
   ! neighbored columns ...

!!$       write(*,*) 'Pt ', p, ' Col. ', i, ' k = ', k(i), ' : ', &
!!$            Hgeo(i,p,k(i)), ' < ', slant%SlantPointsEll(p,3), &
!!$            ' < ', Hgeo(i,p,k(i)-1)
!!$       write(*,*) 'Druck ', col(i,p)%p(k(i)), ' > ', col(i,p)%p(k(i)-1)

   do m=0, 1
      ! upper (m=0) and lower (m=1) layer

      ! Vertical index:
      ! m=0 -  h = k(i)-1 - upper layer
      ! m=1 -  h = k(i)   - lower layer
      h = k(i) - 1 + m

      ! Partial pressure of water vapour e
      e(m,i) = (col(i,p)%p(h) * col(i,p)%q(h)) / (RDRD +   &
                col(i,p)%q(h)*EMRDRD)

      ! Refractivity N
      N(m,i) = NWein(col(i,p)%p(h), col(i,p)%t(h), e(m,i))
      !write(*,*) 'e, N  =  ', e(m,i), N(m,i)

!!$      !if (IpolOpt .eq. 1) then
!!$         ! Interpolate the refractivity at the given point using the
!!$         ! refractivities at the corners of the grid cell
!!$         ! node(:,1,:) - lower level - k+1
!!$         ! node(:,2,:) - upper level - k
!!$         node(i,2-m,1) = col(i,p)%dlon*deg2rad
!!$         node(i,2-m,2) = col(i,p)%dlat*deg2rad
!!$         node(i,2-m,3) = Hgeo(i,p,h)
!!$         node(i,2-m,4) = N(m,i)
!!$         !write(*,*) 'RefracModelGrad> Node ', i, ' lon, lat, h, N = ', &
!!$         !            col(i,p)%dlon, col(i,p)%dlat, Hgeo(i,p,h),  N(m,i)
!!$      !end if

   end do
end do

!!$! Interpolation using function "BiLinearExp3D"
!!$if (slant%Nnghb .eq. 3) then
!!$   ! GME - interpolate between 3 points
!!$   node(4,:,:) = node(2,:,:)
!!$end if
!!$
!!$!Rfrc = BiLinearExp3D(slant%SlantPointsEll(p,:), node, slant%Nnghb)
!!$!write(*,*) 'RefracModelGrad> RefPtEll, slant%Nnghb ',    &
!!$!     RefPtEll(1)*rad22deg, RefPtEll(2)*rad22deg, RefPtEll(3), slant%Nnghb
!!$Rfrc2 = BiLinearExp3D(RefPtEll(:), node, slant%Nnghb)
!!$!write(*,*) 'RefracModelGrad>  Rfrc =',  Rfrc
!!$
!!$node_tl = 0.0_wp
!!$!Rfrc_tl = BiLinearExp3D_tl(RefPt, node, RefPt_tl, node_tl, Npt)
!!$RefPt_tl =  (/1.0_wp,0.0_wp,0.0_wp/)   ! dN / d lambda
!!$!RfrcLambda = BiLinearExp3D_tl(slant%SlantPointsEll(p,:), node,   &
!!$!                              RefPt_tl, node_tl, slant%Nnghb)
!!$RfrcLambda2 = BiLinearExp3D_tl(RefPtEll(:), node,                &
!!$     RefPt_tl, node_tl, slant%Nnghb)
!!$RefPt_tl =  (/0.0_wp,1.0_wp,0.0_wp/)   ! dN / d beta
!!$!RfrcBeta = BiLinearExp3D_tl(slant%SlantPointsEll(p,:), node,   &
!!$!                              RefPt_tl, node_tl, slant%Nnghb)
!!$RfrcBeta2 = BiLinearExp3D_tl(RefPtEll(:), node,   &
!!$     RefPt_tl, node_tl, slant%Nnghb)
!!$RefPt_tl =  (/0.0_wp,0.0_wp,1.0_wp/)   ! dN / d beta
!!$!RfrcHeight = BiLinearExp3D_tl(slant%SlantPointsEll(p,:), node,   &
!!$!                              RefPt_tl, node_tl, slant%Nnghb)
!!$RfrcHeight2 = BiLinearExp3D_tl(RefPtEll(:), node,   &
!!$     RefPt_tl, node_tl, slant%Nnghb)
!!$!write(*,*) 'N = ', N(0,1), N(1,1), ' ipol = ', RefracProfile(p,2)


! ExpInt1D
do i=1, slant%Nnghb
   ! vertical interpolation of the refractivity to the required height
   if (abs(N(1,i)) .lt. 1.0E-10_wp .and. abs(N(0,i)) .lt. 1.0E-10_wp ) then
      ! N1 = N2 = 0  => N(z) = 0
      Nh(i) = 0.0_wp
   else
      ! interpolate
      if (abs(N(1,i)) .lt. 1.0E-10_wp) then
         ! N1 = 0, N2 <> 0
         N(1,i) = N(0,i)/100.0_wp
      end if
      if (abs(N(0,i)) .lt. 1.0E-10_wp) then
         ! N2 = 0, N1 <> 0
         N(0,i) = N(1,i)/100.0_wp
      end if
      Nh(i) = N(0,i) * (N(1,i)/N(0,i))**                     &
           ((Hgeo(i,p,k(i)-1)-RefPtEll(3))/(Hgeo(i,p,k(i)-1)-Hgeo(i,p,k(i))))
      B = ( Hgeo(i,p,k(i)-1) - RefPtEll(3) ) / ( Hgeo(i,p,k(i)-1)-Hgeo(i,p,k(i)) )
      L = log( N(1,i)/N(0,i) )
      EE = exp( B * L )
      dNh(i) = -L/( Hgeo(i,p,k(i)-1) - Hgeo(i,p,k(i)) ) * N(0,i) * EE
   end if
end do
!A =  dNh(2)

! BiLinear2D
! Kopieren, Sonderbehandpung nur fuer Punkte innerhalb des Dreiecks
! => Eckpunkt kopieren

! c1 = x2 - x1
!c1 = node(2,1) - node(1,1)
c1 = (col(2,p)%dlon - col(1,p)%dlon) * deg2rad

! c2 = x4 - x3
!c2 = node(4,1) - node(3,1)
if (slant%Nnghb .eq. 4) then
   c2 = (col(4,p)%dlon - col(3,p)%dlon) * deg2rad
else
   ! slant%Nnghb = 3 : replace node 4 with node 2
   c2 = (col(2,p)%dlon - col(3,p)%dlon) * deg2rad
end if

! d1 = (x - x1)/(x2-x1) =  (x - x1)/c1
!d1 = (RefPt(1) - node(1,1))/c1
d1 = (RefPtEll(1) - col(1,p)%dlon*deg2rad) / c1
Dd1lon = 1.0_wp / c1

! d2 = (x2 - x)/(x2-x1) =  (x2 - x)/c1
!d2 = (node(2,1) - RefPt(1))/c1
d2 = (col(2,p)%dlon*deg2rad - RefPtEll(1)) / c1
Dd2lon = -1.0_wp / c1

! d3 = (x - x3)/(x4-x3) =  (x - x3)/c2
!d3 = (RefPt(1) - node(3,1))/c2
d3 = (RefPtEll(1) - col(3,p)%dlon*deg2rad) / c2
Dd3lon = 1.0_wp / c2

! d4 = (x4 - x)/(x4-x3) =  (x4 - x)/c2
!d4 = (node(4,1) - RefPt(1)) / c2
if (slant%Nnghb .eq. 4) then
   d4 = (col(4,p)%dlon*deg2rad - RefPtEll(1)) / c2
else
   ! slant%Nnghb = 3 : replace node 4 with node 2
   d4 = (col(2,p)%dlon*deg2rad - RefPtEll(1)) / c2
end if
Dd4lon = -1.0_wp / c2
!write(*,*) 'BiLinearDiff2D> c1, c2, d1, d2, d3, d4 : ', c1, c2, d1, d2, d3, d4

!fa = d2*node(1,3) + d1*node(2,3)
!fb = d4*node(3,3) + d3*node(4,3)
!ya = d2*node(1,2) + d1*node(2,2)
!yb = d4*node(3,2) + d3*node(4,2)

!fa = d2*node(1,3) + d1*node(2,3)
fa = d2 * Nh(1) + d1 * Nh(2)
dfa1 = Dd2lon * Nh(1) + Dd1lon * Nh(2)   ! d fa / d lon
dfa3 = d2 * dNh(1) + d1 * dNh(2)    ! d fa / d height,  d fa / d lat = 0

!fb = d4*node(3,3) + d3*node(4,3)
if (slant%Nnghb .eq. 4) then
   fb = d4*Nh(3) + d3*Nh(4)
   dfb1 = Dd4lon * Nh(3) + Dd3lon * Nh(4)   ! d fb / d lon
   dfb3 = d4 * dNh(3) + d3 * dNh(4)    ! d fb / d height, d fb / d lat = 0
else
   ! slant%Nnghb = 3 : replace node 4 with node 2
   fb = d4*Nh(3) + d3*Nh(2)
   dfb1 = Dd4lon * Nh(3) + Dd3lon * Nh(2)   ! d fb / d lon
   dfb3 = d4 * dNh(3) + d3 * dNh(2)     ! d fb / d height, d fb / d lat = 0
end if

!ya = d2*node(1,2) + d1*node(2,2)
ya = (d2*col(1,p)%dlat + d1*col(2,p)%dlat)*deg2rad
Dyalon = (Dd2lon*col(1,p)%dlat + Dd1lon*col(2,p)%dlat)*deg2rad

!yb = d4*node(3,2) + d3*node(4,2)
if (slant%Nnghb .eq. 4) then
   yb = (d4*col(3,p)%dlat + d3*col(4,p)%dlat)*deg2rad
   Dyblon = (Dd4lon*col(3,p)%dlat + Dd3lon*col(4,p)%dlat)*deg2rad
else
   ! slant%Nnghb = 3 : replace node 4 with node 2
   yb = (d4*col(3,p)%dlat + d3*col(2,p)%dlat)*deg2rad
   Dyblon = (Dd4lon*col(3,p)%dlat + Dd3lon*col(2,p)%dlat)*deg2rad
end if

dy = yb - ya
Ddylon = Dyblon - Dyalon

if (abs(dy) .lt. 1.0E-8_wp) then
   ! slant%Nnghb = 3 : reference point is very close to point 1
   Rfrc = Nh(1)
   RfrcHeight = dNh(1)
   RfrcLambda = 0.0_wp
   RfrcBeta   = 0.0_wp
else
   ! No singularity, interpolate ...
   Rfrc = fa*(yb-RefPtEll(2))/dy + fb*(RefPtEll(2)-ya)/dy
   RfrcLambda = Dfa1*(yb - RefPtEll(2))/dy + fa*Dyblon/dy -   &
                fa*Ddylon*(yb - RefPtEll(2))/dy**2  +       &
                Dfb1*(RefPtEll(2)-ya)/dy - fb*Dyalon/dy -   &
                fb*Ddylon*(RefPtEll(2)-ya)/dy**2
   RfrcBeta   = (fb-fa)/dy
   RfrcHeight = dfa3*(yb-RefPtEll(2))/dy + dfb3*(RefPtEll(2)-ya)/dy
end if

!write(*,*) 'RefracModelGrad> Ref.: Rfrc, Dlon, Dlat, Dh = ',  &
!     Rfrc2, RfrcLambda2, RfrcBeta2, RfrcHeight2
!write(*,*) 'RefracModelGrad> Neu.: Rfrc, Dlon, Dlat, Dh = ',  &
!     Rfrc, RfrcLambda, RfrcBeta, RfrcHeight

deallocate( k )
deallocate( e )
deallocate( N )

End subroutine RefracModelGradOld
!****
! <<<<<<<< End of RoboDoc comments


!---------------------------------------------------------------------
! subroutine RefracModelGrad2Old
!---------------------------------------------------------------------
!
! Start of RoboDoc comments >>>>
!****S* Fermat/RefracModelGrad2
!
! Name
! RefracModelGrad2
!
! Call
! call RefracModelGrad2 (slant, col, Hgeo, p, RefPtEll, Rfrc,   &
!                        RfrcLambda, RfrcBeta, RfrcHeight, Grad2Rfrc)
!
! Purpose

! Second derivative of the refractivity N with respect to the
! ellipsoidal coordinates lambda, beta, height.
!
!
! RefracModel      => M
! RefracModelGrad  => N, d N / d lambda,  d N / d beta,  d N / d height
! RefracModelGrad2 => N, dN, d2N
!
! RefracModel computes one single scalar function, i.e. the refractivity N
! at a certain point (lambda, beta, height).
! The Jacobian consists of 3 derivatives:
!
! J = (N_lambda, N_beta, N_height)
!   = (RfrcLambda, RfrcBeta, RfrcHeight)
!
! The Hessian consists of 9 derivatives:
!
!     (
! H = (
!     (
!
! Grad2Rfrc - second derivatives of N with respect to ellipsoidal coordinates
!             Grad2Rfrc(j,i) :  i - first varaiable
!                               j - second variable
!
! Grad2Rfrc(j,i) =
!
!           ( 11 12 13 )   ( N_lambda_lambda N_beta_lambda N_height_lambda )
!           ( 21 22 23 ) = ( N_lambda_beta   N_beta_beta   N_height_beta   )
!           ( 31 32 33 )   ( N_lambda_height N_beta_height N_height_height )
!
!
! Input
! X, Y, Z - cartesian coordinates
! Ellips  - type "EllipsParam" defining the reference ellipsoid
!
! Output
! lambda - ellipsoidal longitude in rad
! beta   - ellipsoidal latitude in rad
! height - ellipsoidal height in m
!
! External References
!
! Examples
!
! Record of revisions
! Date        Programmer   Description of change
! ====        ==========   =====================
! 16.12.2013  M. Bender    new
! 10.03.2014  M. Bender    longitude changed to -Pi <= lon <= +Pi
!---------------------------------------------------------------------
! Source
!---------------------------------------------------------------------
subroutine RefracModelGrad2Old (slant, col, Hgeo, p, RefPtEll, Rfrc,   &
                            RfrcLambda, RfrcBeta, RfrcHeight, Grad2Rfrc)

implicit none

! List of calling arguments:
type (SlantData) ,intent(in)    :: slant        ! slant meta data
type (p_column)  ,intent(in)    :: col (:,:)    ! model columns     (input)
real(wp), dimension(:,:,:), pointer :: Hgeo     ! geoid coordinates (input)
integer, intent(in)                 :: p        ! index of pont in "SlantPoints"
real(wp), dimension(1:3), intent(in) :: RefPtEll ! compute N here
real(wp),intent(out)                :: Rfrc     ! refractivity (output)
real(wp),intent(out)                :: RfrcLambda, RfrcBeta, RfrcHeight
real(wp), dimension(1:3,1:3), intent(out) :: Grad2Rfrc
!real(wp), optional, intent(out)                :: dA

! List of local variables:
real(wp), dimension(:), pointer :: Hcol => Null()
real(wp), dimension(1:3)        :: RefPtEllM

! vertical grid index of neighbored columns:
integer, dimension(:), pointer  :: k => Null()
!integer, dimension(:), pointer  :: k2 => Null()
real(wp), dimension(:,:), pointer :: e => Null()
real(wp), dimension(:,:), pointer :: N => Null()
!real (wp), dimension(1:4,1:2,1:4)  :: node, node_tl
!real (wp), dimension(1:3) :: RefPt_tl
real(wp), dimension(1:4) :: Nh, dNh, d2Nh

real (wp) :: c1, c2, d1, d2, d3, d4
real (wp) :: ya, yb, fa, fb, dy     !, Rfrc2
real (wp) :: B, L, EE, dB, dEE

real (wp) :: Dd1lon,  Dd2lon,  Dd3lon,  Dd4lon
real(wp)  :: dfa1, dfa3, dfb1, dfb3
real (wp) :: d2fa13, d2fa31, d2fa33, d2fb13, d2fb31, d2fb33
real (wp) :: Dyalon, Dyblon, Ddylon

!real(wp) :: RfrcLambda2, RfrcBeta2, RfrcHeight2, a1, a2

integer :: i, m, h
!---------------------------------------------------------------------

allocate( k(1:slant%Nnghb) )
! allocate( k2(1:slant%Nnghb) )
allocate( e(0:1,1:slant%Nnghb) )
allocate( N(0:1,1:slant%Nnghb) )

do i=1, slant%Nnghb
   ! Find the correct vertical indices for the next set of columns
   Hcol => Hgeo(i,p,:)
   !k(i) = LayerSearch (Hcol, slant%SlantPointsEll(p,3))
   !k(i) = LayerSearch (real(Hcol), real(RefPtEll(3)))
   k(i) = LayerSearch (Hcol, RefPtEll(3))
   !write(*,*) 'H k-1, k, k+1 ', Hcol(k(i)-1),  Hcol(k(i)),  Hcol(k(i)+1), RefPtEll(3)
end do

! Convert longitude
! All transformations require longitudes 0° <= lon <= 360°
! but the model provides longitudes between  -180° <= lon <= 180°
! ????? Auch gleich in Grad umrechnen? Das spart die staendige Umrechnung ... ?????
RefPtEllM = RefPtEll
if (RefPtEll(1) .gt. 180.0*deg2rad) RefPtEllM(1) = RefPtEll(1)-360.0*deg2rad

! Compute the partial pressure of water vapour e and the refractivity
! on the corners of the grid cell, i.e 8 corners for COSMO
do i=1, slant%Nnghb
   ! neighbored columns ...

!!$       write(*,*) 'Pt ', p, ' Col. ', i, ' k = ', k(i), ' : ', &
!!$            Hgeo(i,p,k(i)), ' < ', slant%SlantPointsEll(p,3), &
!!$            ' < ', Hgeo(i,p,k(i)-1)
!!$       write(*,*) 'Druck ', col(i,p)%p(k(i)), ' > ', col(i,p)%p(k(i)-1)

   do m=0, 1
      ! upper (m=0) and lower (m=1) layer

      ! Vertical index:
      ! m=0 -  h = k(i)-1 - upper layer
      ! m=1 -  h = k(i)   - lower layer
      h = k(i) - 1 + m

      ! Partial pressure of water vapour e
      e(m,i) = (col(i,p)%p(h) * col(i,p)%q(h)) / (RDRD +   &
                col(i,p)%q(h)*EMRDRD)

      ! Refractivity N
      N(m,i) = NWein(col(i,p)%p(h), col(i,p)%t(h), e(m,i))

      if ( N(m,i) .lt. 0.0_wp) then
         write(*,*) 'RefracModelGrad2> SmithW negative Refrakt. ', N(m,i)
      end if

      !write(*,*) 'e, N  =  ', e(m,i), N(m,i)

!!$      !if (IpolOpt .eq. 1) then
!!$         ! Interpolate the refractivity at the given point using the
!!$         ! refractivities at the corners of the grid cell
!!$         ! node(:,1,:) - lower level - k+1
!!$         ! node(:,2,:) - upper level - k
!!$         node(i,2-m,1) = col(i,p)%dlon*deg2rad
!!$         node(i,2-m,2) = col(i,p)%dlat*deg2rad
!!$         node(i,2-m,3) = Hgeo(i,p,h)
!!$         node(i,2-m,4) = N(m,i)
!!$         !write(*,*) 'RefracModelGrad> Node ', i, ' lon, lat, h, N = ', &
!!$         !            col(i,p)%dlon, col(i,p)%dlat, Hgeo(i,p,h),  N(m,i)
!!$      !end if

   end do
end do

! ExpInt1D
do i=1, slant%Nnghb
   ! vertical interpolation of the refractivity to the required height
   if (abs(N(1,i)) .lt. 1.0E-10_wp .and. abs(N(0,i)) .lt. 1.0E-10_wp ) then
      ! N1 = N2 = 0  => N(z) = 0
      Nh(i) = 0.0_wp
      dNh(i) = 0.0_wp
   else
      ! interpolate
      if (abs(N(1,i)) .lt. 1.0E-10_wp) then
         ! N1 = 0, N2 <> 0
         N(1,i) = N(0,i)/100.0_wp
      end if
      if (abs(N(0,i)) .lt. 1.0E-10_wp) then
         ! N2 = 0, N1 <> 0
         N(0,i) = N(1,i)/100.0_wp
      end if
      Nh(i) = N(0,i) * (N(1,i)/N(0,i))**                     &
           ((Hgeo(i,p,k(i)-1)-RefPtEllM(3))/(Hgeo(i,p,k(i)-1)-Hgeo(i,p,k(i))))
      B = ( Hgeo(i,p,k(i)-1) - RefPtEllM(3) ) / ( Hgeo(i,p,k(i)-1)-Hgeo(i,p,k(i)) )
      L = log( N(1,i)/N(0,i) )
      EE = exp( B * L )
      dNh(i) = (-L/( Hgeo(i,p,k(i)-1) - Hgeo(i,p,k(i)) )) * N(0,i) * EE
      ! second derivatives:
      dB = -1.0_wp / ( Hgeo(i,p,k(i)-1)-Hgeo(i,p,k(i)) )
      dEE = L * dB * EE
      d2Nh(i) = dNh(i) * dEE/EE
   end if

   if ( Nh(i) .lt. 0.0_wp) then
      write(*,*) 'RefracModelGrad2> ExpInt1D negative Refrakt. ', Nh(i)
   end if

end do

! BiLinear2D
! Kopieren, Sonderbehandpung nur fuer Punkte innerhalb des Dreiecks
! => Eckpunkt kopieren

! c1 = x2 - x1
!c1 = node(2,1) - node(1,1)
c1 = (col(2,p)%dlon - col(1,p)%dlon) * deg2rad

! c2 = x4 - x3
!c2 = node(4,1) - node(3,1)
if (slant%Nnghb .eq. 4) then
   c2 = (col(4,p)%dlon - col(3,p)%dlon) * deg2rad
else
   ! slant%Nnghb = 3 : replace node 4 with node 2
   c2 = (col(2,p)%dlon - col(3,p)%dlon) * deg2rad
end if

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! fix
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
if (c1 < 1.0E-8) c1 = 1.0E-7
if (c2 < 1.0E-8) c2 = 1.0E-7
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


! d1 = (x - x1)/(x2-x1) =  (x - x1)/c1
!d1 = (RefPt(1) - node(1,1))/c1
d1 = (RefPtEllM(1) - col(1,p)%dlon*deg2rad) / c1
Dd1lon = 1.0_wp / c1

! d2 = (x2 - x)/(x2-x1) =  (x2 - x)/c1
!d2 = (node(2,1) - RefPt(1))/c1
d2 = (col(2,p)%dlon*deg2rad - RefPtEllM(1)) / c1
Dd2lon = -1.0_wp / c1

! d3 = (x - x3)/(x4-x3) =  (x - x3)/c2
!d3 = (RefPt(1) - node(3,1))/c2
d3 = (RefPtEllM(1) - col(3,p)%dlon*deg2rad) / c2
Dd3lon = 1.0_wp / c2

! d4 = (x4 - x)/(x4-x3) =  (x4 - x)/c2
!d4 = (node(4,1) - RefPt(1)) / c2
if (slant%Nnghb .eq. 4) then
   d4 = (col(4,p)%dlon*deg2rad - RefPtEllM(1)) / c2
else
   ! slant%Nnghb = 3 : replace node 4 with node 2
   d4 = (col(2,p)%dlon*deg2rad - RefPtEllM(1)) / c2
end if
Dd4lon = -1.0_wp / c2
!write(*,*) 'BiLinearDiff2D> c1, c2, d1, d2, d3, d4 : ', c1, c2, d1, d2, d3, d4

! second derivatives of Dd1lon = Dd2lon = Dd3lon = Dd4lon = 0


!fa = d2*node(1,3) + d1*node(2,3)
!fb = d4*node(3,3) + d3*node(4,3)
!ya = d2*node(1,2) + d1*node(2,2)
!yb = d4*node(3,2) + d3*node(4,2)

!fa = d2*node(1,3) + d1*node(2,3)
fa = d2 * Nh(1) + d1 * Nh(2)
! first derivatives of fa, d fa / d lat = 0
dfa1 = Dd2lon * Nh(1) + Dd1lon * Nh(2)   ! d fa / d lon
dfa3 = d2 * dNh(1) + d1 * dNh(2)         ! d fa / d height
! second derivatives of fa
! d dfa1 / d lon = 0,  d dfa1 / d lat = 0
d2fa13 = Dd2lon * dNh(1) + Dd1lon * dNh(2)  ! d dfa1 / d height
! d dfa3 / d lat = 0
d2fa31 = d2fa13                        ! d dfa3 / d lon
d2fa33 = d2 * d2Nh(1) + d1 * d2Nh(2)   ! d dfa3 / d height

!fb = d4*node(3,3) + d3*node(4,3)
if (slant%Nnghb .eq. 4) then
   fb = d4*Nh(3) + d3*Nh(4)
   dfb1 = Dd4lon * Nh(3) + Dd3lon * Nh(4)   ! d fb / d lon
   dfb3 = d4 * dNh(3) + d3 * dNh(4)    ! d fb / d height, d fb / d lat = 0
   d2fb13 =  Dd4lon * dNh(3) + Dd3lon * dNh(4)
   d2fb31 = d2fb13
   d2fb33 = d4 * d2Nh(3) + d3 * d2Nh(4)   ! d dfb3 / d height
else
   ! slant%Nnghb = 3 : replace node 4 with node 2
   fb = d4*Nh(3) + d3*Nh(2)
   dfb1 = Dd4lon * Nh(3) + Dd3lon * Nh(2)   ! d fb / d lon
   dfb3 = d4 * dNh(3) + d3 * dNh(2)     ! d fb / d height, d fb / d lat = 0
   d2fb13 = Dd4lon * dNh(3) + Dd3lon * dNh(2)
   d2fb31 = d2fb13
   d2fb33 = d4 * d2Nh(3) + d3 * d2Nh(2)   ! d dfb3 / d height
end if

!ya = d2*node(1,2) + d1*node(2,2)
ya = (d2*col(1,p)%dlat + d1*col(2,p)%dlat)*deg2rad
Dyalon = (Dd2lon*col(1,p)%dlat + Dd1lon*col(2,p)%dlat)*deg2rad
!write(*,*) 'Dyalon ', Dd2lon, Dd1lon, col(1,p)%dlat, col(2,p)%dlat, Dyalon
!write(*,*) 'p = ', p

!yb = d4*node(3,2) + d3*node(4,2)
if (slant%Nnghb .eq. 4) then
   yb = (d4*col(3,p)%dlat + d3*col(4,p)%dlat)*deg2rad
   Dyblon = (Dd4lon*col(3,p)%dlat + Dd3lon*col(4,p)%dlat)*deg2rad
else
   ! slant%Nnghb = 3 : replace node 4 with node 2
   yb = (d4*col(3,p)%dlat + d3*col(2,p)%dlat)*deg2rad
   Dyblon = (Dd4lon*col(3,p)%dlat + Dd3lon*col(2,p)%dlat)*deg2rad
end if

dy = yb - ya
Ddylon = Dyblon - Dyalon
! second derivatives are all 0

Grad2Rfrc = 0.0_wp
if (abs(dy) .lt. 1.0E-8_wp) then
   ! slant%Nnghb = 3 : reference point is very close to point 1
   Rfrc = Nh(1)
   RfrcHeight = dNh(1)
   RfrcLambda = 0.0_wp
   RfrcBeta   = 0.0_wp
   Grad2Rfrc(3,3) = d2Nh(1)
else
   ! No singularity, interpolate ...

   ! Refractivity N interpolated at the reference point (lon, lat alt)
   Rfrc = fa*(yb- RefPtEllM(2))/dy + fb*(RefPtEllM(2)-ya)/dy

   ! First derivative of N with respect to the coordinates of
   ! the reference point (lon, lat alt) = (lambda, beta, height)
   RfrcLambda = Dfa1*(yb- RefPtEllM(2))/dy + fa*Dyblon/dy -   &
                fa*(yb- RefPtEllM(2))/dy**2 * Ddylon +       &
                Dfb1*(RefPtEllM(2)-ya)/dy - fb*Dyalon/dy -   &
                fb*(RefPtEllM(2)-ya)/dy**2 * Ddylon
   RfrcBeta   = (fb-fa)/dy
   RfrcHeight = dfa3*(yb-RefPtEllM(2))/dy + dfb3*(RefPtEllM(2)-ya)/dy

   ! Second derivative of N with respect to the coordinates of
   ! the reference point
   ! N_height_height
   Grad2Rfrc(3,3) = d2fa33*(yb-RefPtEllM(2))/dy + d2fb33*(RefPtEllM(2)-ya)/dy
   ! N_height_beta
   Grad2Rfrc(2,3) = (dfb3 - dfa3) / dy
   ! N_height_lambda
   Grad2Rfrc(1,3) = ( dy*(dfa3*Dyblon + d2fa31*yb - RefPtEllM(2)*d2fa31) -  &
                      dfa3*Ddylon*(yb-RefPtEllM(2)) +                        &
                      dy*(RefPtEllM(2)*d2fb31 - dfb3*Dyalon - d2fb31*ya ) - &
                      dfb3*Ddylon*(RefPtEllM(2)-ya)  ) / dy**2
   ! N_beta_height
   Grad2Rfrc(3,2) = Grad2Rfrc(2,3)
   ! N_beta_beta
   Grad2Rfrc(2,2) = 0.0_wp
   ! N_beta_lambda
   Grad2Rfrc(1,2) = ( dy*(dfb1-dfa1) - Ddylon*(fb-fa) ) / dy**2
   ! N_lambda_height
   Grad2Rfrc(3,1) = Grad2Rfrc(1,3)
   ! N_lambda_beta
   Grad2Rfrc(2,1) = Grad2Rfrc(1,2)
   ! N_lambda_lambda
                    ! 1) d/dlambda ( Dfa1*(yb - RefPtEllM(2))/dy )
   Grad2Rfrc(1,1) = dfa1*(dy*Dyblon - (yb-RefPtEllM(2))*Ddylon) / dy**2 + &
                    ! 2) d/dlambda (  fa*Dyblon/dy )
                    Dyblon*(dy*dfa1-fa*Ddylon)/dy**2 - &
                    ! 3) d/dlambda ( fa*Ddylon*(yb - RefPtEllM(2))/dy**2 )
                    Ddylon*(yb-RefPtEllM(2)) *                      &
                           (dfa1*dy**2 - 2*dy*fa*Ddylon) / dy**4 - &
                           (fa*Ddylon*Dyblon) / dy**2    - &
                    ! 4) d/dlambda (  Dfb1*(RefPtEllM(2)-ya)/dy )
                    Dfb1*(dy*Dyalon + (RefPtEllM(2)-ya)*Ddylon)/dy**2  - &
                    ! 5) d/dlambda ( fb*Dyalon/dy )
                    Dyalon*(dy*dfb1 - fb*Ddylon) / dy**2  - &
                    ! 6) d/dlambda ( fb*Ddylon*(RefPtEllM(2)-ya)/dy**2 )
                    Ddylon*(RefPtEllM(2)-ya)*  &
                        (dy**2*dfb1 - 2.0_wp*fb*dy*Ddylon) / dy**4 + &
                        fb*Ddylon*Dyalon / dy**2

end if

if (Rfrc .lt. 0.0_wp) then
   write(*,*) 'RefracModelGrad2> BiLin2D negative Refrakt. ', Rfrc, char(10), &
        'slant%Nnghb = ', slant%Nnghb,  char(10), &
        'RefPt(1), Laenge = ', RefPtEllM(1)*rad2deg, char(10), &
        'RefPt(2), Breite = ', RefPtEllM(2)*rad2deg, char(10), &
        'RefPt(3), Hoehe = ', RefPtEllM(3), char(10), &
        'Eckpunkte,  Laenge = ', col(:,p)%dlon,  char(10), &
        'Eckpunkte,  Breite = ', col(:,p)%dlat,  char(10), &
        'Nh = ', Nh,  char(10), &
        'c1 = ', c1,  char(10), &
        'c2 = ', c2,  char(10), &
        'd1 = ', d1,  char(10), &
        'd2 = ', d2,  char(10), &
        'd3 = ', d3,  char(10), &
        'd4 = ', d4,  char(10), &
        'fa = ', fa,  char(10), &
        'fb = ', fb,  char(10), &
        'ya = ', ya,  char(10), &
        'yb = ', yb,  char(10), &
        'dy = ', dy
end if


!write(*,*) dfa1, dy, Dyblon, yb, Ddylon
!write(*,*) dy*Dyblon, yb*Ddylon, dy*Dyblon - yb*Ddylon


!write(*,*) ' RefPtEll(2) = ',  RefPtEll(2),  RefPtEll(2)*rad22deg
!write(*,*) 'N1-4, ipol ', Nh, Rfrc

!write(*,*) 'RefracModelGrad> Ref.: Rfrc, Dlon, Dlat, Dh = ',  &
!     Rfrc2, RfrcLambda2, RfrcBeta2, RfrcHeight2
!write(*,*) 'RefracModelGrad> Neu.: Rfrc, Dlon, Dlat, Dh = ',  &
!     Rfrc, RfrcLambda, RfrcBeta, RfrcHeight

deallocate( k )
deallocate( e )
deallocate( N )

End subroutine RefracModelGrad2Old
!****
! <<<<<<<< End of RoboDoc comments


! Koordinaten der Punkte auf slant
! Koordinaten der Gittersaeulen
! Test: Liegt der Stuetzpunkt zwischen den Gittersaeulen?
! Test: Abstand zwischen den Gittersaeulen
! Test: Eckpunkte fuer Interpolation

subroutine CheckModuleColumns(slant, col)

 type (SlantData) ,intent(inout) :: slant        ! slant meta data
 type (p_column)  ,intent(in)    :: col (:,:)    ! model columns     (input)


 integer :: i, j, len
 character (len=512) :: row

 ! --------------------------------------------------------------------

 !do i=1, ubound(slant%SlantPointsEll,1)
 do i=1, ubound(col,2)
    row = ''
    write(row,*) 'Supporting point on transmitter - receiver axis : ', &
         slant%SlantPointsEll(i,1)*rad2deg,                            &
         slant%SlantPointsEll(i,2)*rad2deg,                            &
         slant%SlantPointsEll(i,3),  char(10)
    len = len_trim(row) + 2
    do j=1, slant%Nnghb
       write(row(len:),*) 'Grid column : ',   &
            j, col(j,i)%dlon, col(j,i)%dlat, char(10)
       len = len_trim(row) + 2
    end do
    write(*,*) trim(row)
 end do

end subroutine CheckModuleColumns





function GetFreeUnit ()

implicit none

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



!==============================================================================
! Old code without raytracer but with tangent-linear and adjoint
!==============================================================================

 subroutine std_path_coord_line (STD, mod_id, mod_htop, mod_nz, mod_dx, err2)
 !----------------------------------------------------------------------
 ! Compute the coordinates of the supporting points along the slant path
 ! => Part 1 of the STD operator
 ! as the original std_path_coord, but modified for the 3dvar environment:
 !   only calculate coordinates
 !   do not derive gridcells
 !   no check for out of domain conditions
 !   routine is called for each single slant, i.e. processes only one slant
 !----------------------------------------------------------------------
 type (SlantData) ,intent(inout) :: STD
 integer          ,intent(in)    :: mod_id    ! model identifier
 real(wp)         ,intent(in)    :: mod_htop  ! height in m of model top level
 integer          ,intent(in)    :: mod_nz    ! number of vertical levels
 real(wp)         ,intent(in)    :: mod_dx    ! horizontal grid spacing, m
 integer          ,intent(out)   :: err2

 ! htop - altitude of model top level in meters  ??
 !        frame of reference ?? height above geoid? ellipsoid?

 real (wp) :: px, py, pz, gx, gy, gz
!real (wp) :: TestLat, TestLon, DeltaLat, DeltaLon

 integer :: a      !, err, j, k,  ii, jj
!integer :: stat, unit
!logical :: opened
 !integer :: Nerr, Ncosmo, Nperiod, Nhorz, Nstd, Nperr, Nassim

 !integer                              :: inCOSMO

 type (Line)               :: ray
 real (wp), dimension(1:3) :: CrossPt, point, wgspoint,  point2,  wgspoint2
 real (wp)                 :: RayLambda
 real (wp)                 :: DeltaP, P, lnP0, Ps, Pe, h
 real (wp), parameter      :: P0 = 1024.0D0  ! reference pressure at sea level
 real (wp), parameter      :: rE = 6378137.0_wp ! radius of Earth in m
 real (wp)                 :: Hdist             ! horizontal distance
 real (wp), dimension(:), allocatable :: Hprofile
 !real (wp), dimension(:), allocatable :: SlantSteps

 logical :: debug = .false.
 !---------------------------------------------

!debug = .true.
 if (verbose > 1) debug = .true.
 if (debug) write(*,*) 'std_path_coord_line> Start ...'

!!$  Nerr = 0
!!$  Ncosmo = 0
!!$  Nperiod = 0
!!$  Nhorz = 0
!!$  Nstd = 0
!!$  Nperr = 0
!!$
!!$ write(*,*) 'std_path_coord_line> Geoid lon = ',                  &
!!$                   std% Station%CoordGeo(1)*rad2deg, char(10), &
!!$            '                  Geoid lat = ',                  &
!!$                   std% Station%CoordGeo(2)*rad2deg, char(10), &
!!$            '                  Geoid hei = ',                  &
!!$                   std% Station%CoordGeo(3), char(10),         &
!!$            '                  Ell  lon = ',                   &
!!$                   std% Station%CoordEll(1)*rad2deg, char(10), &
!!$            '                  Ell   lat = ',                  &
!!$                   std% Station%CoordEll(2)*rad2deg, char(10), &
!!$            '                  Ell   hei = ',                  &
!!$                   std% Station%CoordEll(3)

 if (debug) then
 write(*,*) 'std_path_coord_line> Station ', std% Station%SName, char(10), &
            'Geoid/Ellipsoid lon = ', std% Station%CoordGeo(1)*rad2deg, &
                       std% Station%CoordEll(1)*rad2deg, char(10),      &
            'Geoid/Ellipsoid lat = ', std% Station%CoordGeo(2)*rad2deg, &
                       std% Station%CoordEll(2)*rad2deg, char(10),      &
            'Geoid/Ellipsoid hei = ', std% Station%CoordGeo(3) ,        &
                       std% Station%CoordEll(3)
end if


 !-----------------------------------------
 ! Slant loop
 ! Slants are sorted by station and by time
 !-----------------------------------------
 if (.not. associated(std% Station) ) then
    write(*,*) 'No station data available => skip'
    !Nerr = Nerr + 1
    return
 end if

 if (debug) write(*,*) 'StatID, StatName ', std% Station%ID, std% Station%SName
 if (std% Station%ID < 0) then
    write(*,*) 'Error: std_path_coord_line, invalid station id =',  std% Station%ID
    !Nerr = Nerr + 1
    return
 end if
 !----------------------------------------------------------------
 ! Compute supporting points on the slant path
 !----------------------------------------------------------------
 ! Compute straight line from the receiver to the GNSS satellite:
 ! Transform azimuth and elevation in the local horizon system
 ! to ECEF cartesian coordinates
 !----------------------------------------------------------------
 if (debug) write(*,*) 'std_path_coord_line> Transform local => ECEF'
 call Azimut2Horz(std% Obs%azimuth, pi05-std% Obs%elevation, 5.0_wp,  &
                  px, py, pz)
 call LocalHorz2Cart(px, py, pz,                 &
                     std% Station%CoordEll(1),   &
                     std% Station%CoordEll(2),   &
                     std% Station%CoordCart(1),  &
                     std% Station%CoordCart(2),  &
                     std% Station%CoordCart(3),  &
                     gx, gy, gz)

 call LineDef(ray,  std% Station%CoordCart, (/gx,gy,gz/), 3)
 ! ??????????????
 ! htop auch transformieren ???
 ! ??????????????
 call CrossLineEllips(WGS84Param, mod_htop, ray,                    &
                      CrossPt, RayLambda)
 call Cart2Ellips(CrossPt(1),  CrossPt(2),  CrossPt(3),         &
                  point(1), point(2), point(3), WGS84Param)
 !-----------------------------------------------------
 ! The signal path within the COSMO grid is now defined
 !-----------------------------------------------------

 ! uuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuu
 call LinePos(ray, RayLambda, point2)
 call Cart2Ellips(point2(1),  point2(2),  point2(3),                     &
                  wgspoint2(1), wgspoint2(2), wgspoint2(3), WGS84Param)
 !write(*,*) 'Hoehe Kreuzungspunkt : ', point(3), mod_htop
 !write(*,*) 'Hoehe Geradenpunkt   : ',  wgspoint2(3)
 ! uuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuu

 !--------------------------------------------------------
 ! Scale the number of supporting points on the slant path
 ! according to the scaling option "NStepOpt".
 !--------------------------------------------------------
 select case (NStepOpt)
    case (1)
       !--------------------------------------------------------
       ! Scale the number of supporting points on the slant path
       ! with 1/sin(elevation)
       !--------------------------------------------------------
       if (debug) write(*,*)                                           &
            'std_path_coord_line> Start computing supporting points:', &
            NStepVertMod,NStepVertTop, std% Obs%elevation,             &
            sin(std% Obs%elevation)

       std% Nmod = int(NStepVertMod/sin(std% Obs%elevation))
       std% Nup  = int(NStepVertTop/sin(std% Obs%elevation))
       std% Ntot = std% Nmod + std% Nup
       if (debug) write(*,*) 'std_path_coord_line> Nmod, Nup, Ntot : ',          &
                                          std%Nmod, std%Nup, std%Ntot

    case (2, 3)
       !--------------------------------------------------------
       !
       !--------------------------------------------------------

       ! horizontal distance between GNSS station and exit point of slant
       ! Hdist in meter (approx. sphere)
       Hdist = rE*( pi05 - STD%Obs%elevation -                           &
                    asin( cos(STD%Obs%elevation)*rE/(rE+mod_htop) ) )

       ! mod_dx = 2.8 km: COSMO grid spacing
       ! mod_dx =  13 km: ICON grid spacing
       std% Nmod = NStepVertMod + int(1.5 * Hdist/mod_dx)

       ! horizontal distance between GNSS station and exit point of slant
       ! Hdist in meter
       Hdist = rE*( pi05 - STD%Obs%elevation -                           &
                    asin( cos(STD%Obs%elevation)*rE/(rE+Hmax) ) )

       ! add one supporting point for 50 km
       std% Nup  = NStepVertTop + int(Hdist/50000.0_wp)
       if (NStepOpt == 2) then
          !std% Ntot = std% Nmod + std% Nup + 1   ! ?????
          std% Ntot = std% Nmod + std% Nup
       else
          std% Ntot = std% Nmod + std% Nup
       end if

    case default
       write(*,*) 'std_path_coord_line> Unsupported option: NStepOpt = ', &
                  NStepOpt
       std% Nmod = NStepVertMod
       std% Nup  = NStepVertTop
       std% Ntot = std% Nmod + std% Nup + 1
 end select

 !--------------------------------------------------------------------
 ! Supporting points for the refractivity profile along the slant path
 !--------------------------------------------------------------------
 allocate ( std% SlantSteps(1:std% Ntot) )
 allocate ( std% SlantPointsCart(1:std% Ntot,1:3) )
 allocate ( std% SlantPointsEll(1:std% Ntot,1:3) )

 ! start at the GNSS receiver
 STD%SlantSteps(1)           = 0.0_wp  ! distance to receiver in m
 STD%SlantPointsCart(1,:)    = STD%Station%CoordCart
 ! last point in the model grid
 STD%SlantSteps(STD%Nmod) = RayLambda ! distance to receiver in m
 STD%SlantPointsCart(STD%Nmod,:) = CrossPt

 if (NStepOpt <= 2) then
    ! Scale the density of supporting points with the air pressure
    ! P0 is the reference pressure at sea level, used to define
    ! a simple exponential profile: P(h) = P0*exp(-h/HscaleP)
    ! Use HscaleP for pressure profile within the model

    lnP0 = log( P0 )

    ! Pressure at the GPS station
    Ps = P0*exp(-STD%Station%CoordEll(3)/HscaleP)
    ! Pressure at the upper grid level
    Pe = P0*exp(-point(3)/HscaleP)

    DeltaP = (Ps - Pe) / real(STD%Nmod - 1)

    p =  Ps
    do a=2, STD%Nmod-1
       ! Compute supporting points at equidistant pressure levels
       p = p - DeltaP
       h =  -HscaleP*(log(p) - lnP0)
       call CrossLineEllips(WGS84Param, h, ray,         &
                            CrossPt, RayLambda)
       STD%SlantSteps(a) = RayLambda
       STD%SlantPointsCart(a,:) = CrossPt
    end do

    if ( STD%SlantSteps(STD%Nmod-1) .ge.                       &
         STD%SlantSteps(STD%Nmod)         ) then
       ! Check if supporting points are really below the upmost model level
       write(*,*) 'Error: Slant distance too large: ', &
            STD%SlantSteps(STD%Nmod-1), STD%SlantSteps(STD%Nmod)
    end if

    !------------------------------------------------------------
    ! Do the same for the slant path above the upmost model level
    !------------------------------------------------------------
    call CrossLineEllips(WGS84Param, Hmax, ray,         &
                         CrossPt, RayLambda)
    !call Cart2Ellips(CrossPt(1),  CrossPt(2),  CrossPt(3),         &
    !                 point(1), point(2), point(3), WGS84Param)
    STD%SlantSteps(STD%Ntot) = RayLambda ! last point on slant
    STD%SlantPointsCart(STD%Ntot,:) = CrossPt

    ! Pressure at the upper grid level
    ! Use HscaleP2 for pressure profile above the model top
    Ps = P0*exp(-point(3)/HscaleP2)
    ! Pressure at Hmax
    Pe = P0*exp(-Hmax/HscaleP2)

    DeltaP = (Ps - Pe) / real(STD%Nup)

    p =  Ps
    do a=STD%Nmod+1, STD%Ntot-1
       ! Compute supporting points at equidistant pressure levels
       p = p - DeltaP
       h =  -HscaleP2*(log(p) - lnP0)
       call CrossLineEllips(WGS84Param, h, ray,         &
                            CrossPt, RayLambda)
       STD%SlantSteps(a) = RayLambda
       STD%SlantPointsCart(a,:) = CrossPt
    end do

 else   ! NStepOpt = 3
    ! use quadratic profile with several intervals

    ! Compute height above ellipsoid of all supporting points => Hprofile
    call ProfileSquared(std, mod_htop, std%Station%CoordEll(3), Hprofile)

    do a=2, STD%Ntot
       ! Compute supporting points on the ellipsoid at height "Hprofile"
       call CrossLineEllips(WGS84Param, Hprofile(a), ray, CrossPt, RayLambda)
       STD%SlantSteps(a) = RayLambda
       STD%SlantPointsCart(a,:) = CrossPt
    end do

 end if    ! ... NStepOpt

 do a=1, std%Ntot
    !-----------------------------------------------------------------------
    ! SlantPointsCart - cartesian coordinates of a point on the slant path
    !                   with a distance = SlantSteps(a) to the GNSS receiver
    !-----------------------------------------------------------------------
    call Cart2Ellips( STD%SlantPointsCart(a,1),                &
                      STD%SlantPointsCart(a,2),                &
                      STD%SlantPointsCart(a,3),                &
                      wgspoint(1), wgspoint(2), wgspoint(3),   &
                      WGS84Param                             )
    !-----------------------------------------------------
    ! wgspoint - ellipsoidal coordinates of the same point
    !-----------------------------------------------------
    std% SlantPointsEll(a,:) = wgspoint
    !write(*,*) std%SlantPointsEll(a,1)*rad2deg,   &
    !           std%SlantPointsEll(a,2)*rad2deg, std%SlantPointsEll(a,3)

 end do

 if (verbose > 2) then
    write(*,*) 'std_path_coord_line> Coordinates of supporting points'
    do a=1, STD%Ntot-1
       write(*,'(a,i4,tr2,2(f7.2,tr2),tr2,f9.2)') &
            'std_path_coord_line> supporting points lat/lon/height = ', a, &
            std%SlantPointsEll(a,1)*rad2deg,   &
            std%SlantPointsEll(a,2)*rad2deg,   &
            std%SlantPointsEll(a,3)
    end do
 end if

 if (debug) write(*,*) 'std_path_coord_line>  ... end'

end subroutine std_path_coord_line


subroutine std_path_coord_line_old (STD, htop, err2)
 !----------------------------------------------------------------------
 ! Compute the coordinates of the supporting points along the slant path
 ! => Part 1 of the STD operator
 ! as the original std_path_coord, but modified for the 3dvar environment:
 !   only calculate coordinates
 !   do not derive gridcells
 !   no check for out of domain conditions
 !   routine is called for each single slant, i.e. processes only one slant
 !----------------------------------------------------------------------
 type (SlantData) ,intent(inout) :: STD
 real(wp)         ,intent(in)    :: htop
 integer          ,intent(out)   :: err2

 ! htop - altitude of model top level in meters  ??
 !        frame of reference ?? height above geoid? ellipsoid?

 real (wp) :: px, py, pz, gx, gy, gz
!real (wp) :: TestLat, TestLon, DeltaLat, DeltaLon

 integer :: a      !, err, j, k,  ii, jj
!integer :: stat, unit
!logical :: opened
 !integer :: Nerr, Ncosmo, Nperiod, Nhorz, Nstd, Nperr, Nassim

 !integer                              :: inCOSMO

 type (Line)               :: ray
 real (wp), dimension(1:3) :: CrossPt, point, wgspoint,  point2,  wgspoint2
 real (wp)                 :: RayLambda
 real (wp)                 :: DeltaP, P, lnP0, Ps, Pe, h
 real (wp), parameter      :: P0 = 1024.0D0  ! reference pressure at sea level

 !real (wp), dimension(:), allocatable :: SlantSteps

 logical :: debug = .false.
 !---------------------------------------------

!debug = .true.
 if (verbose > 1) debug = .true.
 if (debug) write(*,*) 'std_path_coord_line> Start ...'

!!$  Nerr = 0
!!$  Ncosmo = 0
!!$  Nperiod = 0
!!$  Nhorz = 0
!!$  Nstd = 0
!!$  Nperr = 0
!!$
!!$ write(*,*) 'std_path_coord_line> Geoid lon = ',                  &
!!$                   std% Station%CoordGeo(1)*rad2deg, char(10), &
!!$            '                  Geoid lat = ',                  &
!!$                   std% Station%CoordGeo(2)*rad2deg, char(10), &
!!$            '                  Geoid hei = ',                  &
!!$                   std% Station%CoordGeo(3), char(10),         &
!!$            '                  Ell  lon = ',                   &
!!$                   std% Station%CoordEll(1)*rad2deg, char(10), &
!!$            '                  Ell   lat = ',                  &
!!$                   std% Station%CoordEll(2)*rad2deg, char(10), &
!!$            '                  Ell   hei = ',                  &
!!$                   std% Station%CoordEll(3)

 if (debug) then
 write(*,*) 'std_path_coord_line> Station ', std% Station%SName, char(10), &
            'Geoid/Ellipsoid lon = ', std% Station%CoordGeo(1)*rad2deg, &
                       std% Station%CoordEll(1)*rad2deg, char(10),      &
            'Geoid/Ellipsoid lat = ', std% Station%CoordGeo(2)*rad2deg, &
                       std% Station%CoordEll(2)*rad2deg, char(10),      &
            'Geoid/Ellipsoid hei = ', std% Station%CoordGeo(3) ,        &
                       std% Station%CoordEll(3)
end if


 !-----------------------------------------
 ! Slant loop
 ! Slants are sorted by station and by time
 !-----------------------------------------
 if (.not. associated(std% Station) ) then
    write(*,*) 'No station data available => skip'
    !Nerr = Nerr + 1
    return
 end if

 if (debug) write(*,*) 'StatID, StatName ', std% Station%ID, std% Station%SName
 if (std% Station%ID < 0) then
    write(*,*) 'Error: std_path_coord_line, invalid station id =',  std% Station%ID
    !Nerr = Nerr + 1
    return
 end if
 !----------------------------------------------------------------
 ! Compute supporting points on the slant path
 !----------------------------------------------------------------
 ! Compute straight line from the receiver to the GNSS satellite:
 ! Transform azimuth and elevation in the local horizon system
 ! to ECEF cartesian coordinates
 !----------------------------------------------------------------
 if (debug) write(*,*) 'std_path_coord_line> Transform local => ECEF'
 call Azimut2Horz(std% Obs%azimuth, pi05-std% Obs%elevation, 5.0_wp,  &
                  px, py, pz)
 call LocalHorz2Cart(px, py, pz,                 &
                     std% Station%CoordEll(1),   &
                     std% Station%CoordEll(2),   &
                     std% Station%CoordCart(1),  &
                     std% Station%CoordCart(2),  &
                     std% Station%CoordCart(3),  &
                     gx, gy, gz)

 call LineDef(ray,  std% Station%CoordCart, (/gx,gy,gz/), 3)
 ! ??????????????
 ! htop auch transformieren ???
 ! ??????????????
 call CrossLineEllips(WGS84Param, htop, ray,                    &
                      CrossPt, RayLambda)
 call Cart2Ellips(CrossPt(1),  CrossPt(2),  CrossPt(3),         &
                  point(1), point(2), point(3), WGS84Param)
 !-----------------------------------------------------
 ! The signal path within the COSMO grid is now defined
 !-----------------------------------------------------

 ! uuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuu
 call LinePos(ray, RayLambda, point2)
 call Cart2Ellips(point2(1),  point2(2),  point2(3),                     &
                  wgspoint2(1), wgspoint2(2), wgspoint2(3), WGS84Param)
 !write(*,*) 'Hoehe Kreuzungspunkt : ', point(3), htop
 !write(*,*) 'Hoehe Geradenpunkt   : ',  wgspoint2(3)
 ! uuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuu

 !--------------------------------------------------------
 ! Scale the number of supporting points on the slant path
 ! with 1/sin(elevation)
 !--------------------------------------------------------
 if (debug) write(*,*) 'std_path_coord_line> Start computing supporting points:',&
   NStepVertMod,NStepVertTop,std% Obs%elevation, sin(std% Obs%elevation)
 std% Nmod = int(NStepVertMod/sin(std% Obs%elevation))
 std% Nup  = int(NStepVertTop/sin(std% Obs%elevation))
 std% Ntot = std% Nmod + std% Nup
 if (debug) write(*,*) 'std_path_coord_line> Nmod, Nup, Ntot : ',          &
                                          std%Nmod, std%Nup, std%Ntot

 !--------------------------------------------------------------------
 ! Supporting points for the refractivity profile along the slant path
 !--------------------------------------------------------------------
 allocate ( std% SlantSteps(1:std% Ntot) )
 allocate ( std% SlantPointsCart(1:std% Ntot,1:3) )
 allocate ( std% SlantPointsEll(1:std% Ntot,1:3) )

 ! start at the GNSS receiver
 STD%SlantSteps(1)           = 0.0_wp  ! distance to receiver in m
 STD%SlantPointsCart(1,:)    = STD%Station%CoordCart
 ! last point in the model grid
 STD%SlantSteps(STD%Nmod) = RayLambda ! distance to receiver in m
 STD%SlantPointsCart(STD%Nmod,:) = CrossPt

 ! Scale the density of supporting points with the air pressure
 ! P0 is the reference pressure at sea level, used to define
 ! a simple exponential profile: P(h) = P0*exp(-h/HscaleP)
 ! Use HscaleP for pressure profile within the model
 lnP0 = log( P0 )

 ! Pressure at the GPS station
 Ps = P0*exp(-STD%Station%CoordEll(3)/HscaleP)
 ! Pressure at the upper grid level
 Pe = P0*exp(-point(3)/HscaleP)

 DeltaP = (Ps - Pe) / real(STD%Nmod - 1)

 p =  Ps
 do a=2, STD%Nmod-1
    ! Compute supporting points at equidistant pressure levels
    p = p - DeltaP
    h =  -HscaleP*(log(p) - lnP0)
    call CrossLineEllips(WGS84Param, h, ray,         &
                         CrossPt, RayLambda)
    STD%SlantSteps(a) = RayLambda
    STD%SlantPointsCart(a,:) = CrossPt
 end do

 if ( STD%SlantSteps(STD%Nmod-1) .ge.                       &
      STD%SlantSteps(STD%Nmod)         ) then
    ! Check if supporting points are really below the upmost model level
    write(*,*) 'Error: Slant distance too large: ', &
         STD%SlantSteps(STD%Nmod-1), STD%SlantSteps(STD%Nmod)
 end if

 !------------------------------------------------------------
 ! Do the same for the slant path above the upmost model level
 !------------------------------------------------------------
 call CrossLineEllips(WGS84Param, Hmax, ray,         &
                      CrossPt, RayLambda)
 !call Cart2Ellips(CrossPt(1),  CrossPt(2),  CrossPt(3),         &
 !                 point(1), point(2), point(3), WGS84Param)
 STD%SlantSteps(STD%Ntot) = RayLambda ! last point on slant
 STD%SlantPointsCart(STD%Ntot,:) = CrossPt

 ! Pressure at the upper grid level
 ! Use HscaleP2 for pressure profile above the model top
 Ps = P0*exp(-point(3)/HscaleP2)
 ! Pressure at Hmax
 Pe = P0*exp(-Hmax/HscaleP2)

 DeltaP = (Ps - Pe) / real(STD%Nup)

 p =  Ps
 do a=STD%Nmod+1, STD%Ntot-1
    ! Compute supporting points at equidistant pressure levels
    p = p - DeltaP
    h =  -HscaleP2*(log(p) - lnP0)
    call CrossLineEllips(WGS84Param, h, ray,         &
                         CrossPt, RayLambda)
    STD%SlantSteps(a) = RayLambda
    STD%SlantPointsCart(a,:) = CrossPt
 end do

 ! So far, the distance to the GNSS receiver and the cartesian coordinates
 ! have been computed for all supporting ponts. Now compute the
 ! corresponding WGS84 coordinates.
!write(*,*)
!write(*,*) std%station%sname, std% Obs%azimuth*rad2deg, &
!           std% Obs%elevation*rad2deg,  shape(std%SlantPointsEll)

 do a=1, std%Ntot
    !-----------------------------------------------------------------------
    ! SlantPointsCart - cartesian coordinates of a point on the slant path
    !                   with a distance = SlantSteps(a) to the GNSS receiver
    !-----------------------------------------------------------------------
    call Cart2Ellips( STD%SlantPointsCart(a,1),                &
                      STD%SlantPointsCart(a,2),                &
                      STD%SlantPointsCart(a,3),                &
                      wgspoint(1), wgspoint(2), wgspoint(3),   &
                      WGS84Param                             )
    !-----------------------------------------------------
    ! wgspoint - ellipsoidal coordinates of the same point
    !-----------------------------------------------------
    std% SlantPointsEll(a,:) = wgspoint
    !write(*,*) std%SlantPointsEll(a,1)*rad2deg,   &
    !           std%SlantPointsEll(a,2)*rad2deg, std%SlantPointsEll(a,3)

 end do

 if (verbose > 2) then
    write(*,*) 'std_path_coord_line> Coordinates of supporting points'
    do a=1, STD%Ntot-1
       write(*,'(a,i4,tr2,2(f7.2,tr2),tr2,f9.2)') &
            'std_path_coord_line> supporting points lat/lon/height = ', a, &
            std%SlantPointsEll(a,1)*rad2deg,   &
            std%SlantPointsEll(a,2)*rad2deg,   &
            std%SlantPointsEll(a,3)
    end do
 end if

 if (debug) write(*,*) 'std_path_coord_line>  ... end'

end subroutine std_path_coord_line_old


 subroutine std_delay_line (slant, col, delay, ladj, col_ad, err)
 !-------------------------------------------------------------------------
 ! Compute the slant total delay
 ! as the original std_delay but interface adapted to the 3dvar environment
 ! +++ work in progress +++
 !
 !
 ! col  -  vertical grid columns
 !         col(i,j), i=1, Nnghb, numberof neighbored columns
 !                               Nnghb = 3 - GME
 !                               Nnghb = 4 - COSMO
 !                   j=1, Nhor, number of slant points within grid
 !-------------------------------------------------------------------------
 type (SlantData) ,intent(inout) :: slant        ! slant meta data
 type (p_column)  ,intent(in)    :: col (:,:)    ! model columns     (input)
 real(wp)         ,intent(out)   :: delay        ! Slant Total Delay (output)
 logical          ,intent(in)    :: ladj         ! flag to calculate adjoint
 type (p_column)  ,pointer       :: col_ad (:,:) ! adjoint
 integer          ,intent(out)   :: err          ! return code

 real(wp), dimension(:,:,:), pointer :: Hgeo => Null()
 real(wp), dimension(:), pointer     :: Hcol => Null()
 real(wp), dimension(:,:), pointer   :: RefracProfile => Null()
!real(wp), dimension(:,:), pointer   :: RP => Null()
 integer  :: i, j, m, h, kk !, kmax, kmin

 ! vertical grid index of neighbored columns:
 integer, dimension(:), pointer  :: k => Null()
 !integer, dimension(:), pointer  :: k2 => Null()
 real(wp), dimension(:,:), pointer :: e => Null()
 real(wp), dimension(:,:), pointer :: N => Null()
 real (wp), dimension(1:4,1:2,1:4)  :: node

 ! Test of tangent-linear and adjoint code
 type (p_column)  ,pointer       :: col_tl (:,:) ! tangent_linear (input)
 type (p_column)  ,pointer       :: coltest_ad (:,:) ! adjoint (output)
 type (p_column)  ,pointer       :: Dcol (:,:) !
 real(wp)                        :: delay_tl     ! Slant Total Delay, TL
 real(wp)                        :: delay_ad     ! Slant Total Delay, adjoint
 real(wp)                        :: Ddelay, Diff, epsilon
 integer                         :: a, b
 real(wp)                        :: sumx, sumy  !, delay2
 ! Test of tangent-linear and adjoint code

 real(wp)                       :: STD2         !, STD
!real(wp)                       :: STD3, STD4
!real(wp)                       :: STDmod, STDtop
 real(dp)                       :: Rfrc

#ifndef __COSMO__
 type (Geodetic) :: EllCoord
#endif

!real(wp), dimension(:), pointer :: Nvert => Null()
!real(wp), dimension(1:2,1:2)    :: vnode

 real(wp), dimension(:,:), pointer :: e_ad => Null()
 real(wp), dimension(:,:), pointer :: N_ad => Null()
 real(wp), dimension(1:4,1:2,1:4)  :: node_ad
 real(wp), dimension(:,:), pointer :: RefracProfile_ad => Null()
!real(wp), dimension(1:2,1:2)      :: vnode_ad
!real(dp)                          :: Rfrc_ad, RPt_ad
 real(dp)                          :: STD2_ad
 real(wp), dimension(1:3)          :: RefPt_ad
 real(dp)                      :: p_ad, t_ad

 integer, save :: Nfile = 20

 logical :: testad = .false.
 logical :: debug = .false.
 ! .....

 !debug = .true.
 if (debug) write(*,*)  'std_delay_line> Start ...'
 if (debug .and. slant%Nnghb .eq. 3) write(*,*) 'GME'
 if (debug .and. slant%Nnghb .eq. 4) write(*,*) 'COSMO'

 if (debug) then
    write(*,*) 'Dimension Columns ', shape(col)
    write(*,*) 'Punkte auf Slant ', slant%Ntot
    write(*,*) 'Punkte im Gitter :', slant%Nmod, slant%Nhor
    write(*,*) 'Nachbarzellen    :', slant%Nnghb, slant%Naloc
 end if

 delay = invalid   ! init. delay

 ! Apply geoid corrections
 allocate(Hgeo( lbound(col,1):ubound(col,1),                       &
                lbound(col,2):ubound(col,2),                       &
                lbound(col(1,1)%gpm,1):ubound(col(1,1)%gpm,1) ) )
 !write(*,*) 'Dimension Hgeo ', shape(Hgeo)

 do j=lbound(col,2), ubound(col,2)
    do i=lbound(col,1), ubound(col,1)
       Hgeo(i,j,:) = col(i,j)%gpm(:) + col(i,j)%geoid
    end do
 end do

!!$ do j=lbound(col,2), ubound(col,2)
!!$    write(*,*) 'dlon, dlat ', col(1,j)%dlon, col(1,j)%dlat,  &
!!$         col(1,j)%dlon*rad2deg, col(1,j)%dlat*rad2deg
!!$ end do

!!$write(*,*)
!!$write(*,*) 'Dimension slant%SlantPointsEll ', shape(slant%SlantPointsEll)
!!$write(*,*) 'Dimension slant%SlantPointsEll ', lbound(slant%SlantPointsEll,1), &
!!$           ubound(slant%SlantPointsEll,1), lbound(slant%SlantPointsEll,2),    &
!!$           ubound(slant%SlantPointsEll,2)
!!$
!!$write(*,*) slant%station%sname, slant% Obs%azimuth*rad2deg, &
!!$           slant% Obs%elevation*rad2deg
!!$ do i=1, slant%Ntot
!!$
!!$
!!$    write(*,*) slant%SlantPointsEll(i,1)*rad2deg,   &
!!$               slant%SlantPointsEll(i,2)*rad2deg, slant%SlantPointsEll(i,3)
!!$
!!$ end do
!!$
!!$
!!$  write(*,*) 'H"ohenlevel:'
!!$ do i=1,  ubound(col(1,1)%gpm,1)
!!$    write(*,*) i, Hgeo(1,1,i), col(1,1)%gpm(i)
!!$ end do

 ! check if columns are equal
if (.false.) then
 do j=lbound(col,2), ubound(col,2)    ! j=1, slant%Nmod
    do i=lbound(col,1), ubound(col,1) ! i=1, Nnghb

       ! compare lon/lat
       if (abs(col(i,j)%dlon - col(1,1)%dlon) > 1.0e-8_wp) then
          write(*,*) 'std_delay_line> Col dlon differ'
       end if
       if (abs(col(i,j)%dlat - col(1,1)%dlat) > 1.0e-8_wp) then
          write(*,*) 'std_delay_line> Col dlat differ'
       end if

       do kk=1, ubound(col(1,1)%gpm,1)
          if (abs(col(i,j)%gpm(kk) - col(1,1)%gpm(kk)) > 1.0e-8_wp) then
             write(*,*) 'std_delay_line> Col gpm differ'
          end if
          if (abs(col(i,j)%p(kk) - col(1,1)%p(kk)) > 1.0e-8_wp) then
             write(*,*) 'std_delay_line> Col p differ'
          end if
          if (abs(col(i,j)%q(kk) - col(1,1)%q(kk)) > 1.0e-8_wp) then
             write(*,*) 'std_delay_line> Col q differ'
          end if
          if (abs(col(i,j)%t(kk) - col(1,1)%t(kk)) > 1.0e-8_wp) then
             write(*,*) 'std_delay_line> Col t differ'
          end if
       end do

    end do
 end do
end if

 allocate( RefracProfile(1:slant%Ntot,1:2) )
 allocate( k(1:slant%Nnghb) )
! allocate( k2(1:slant%Nnghb) )

 allocate( e(0:1,1:slant%Nnghb) )
 allocate( N(0:1,1:slant%Nnghb) )

 ! RefracProfile(i,j)
 !  i - supporting points on the slant path
 !  j - j = 1 distance on slant path (x coordinate)
 !            as given in slant%SlantSteps
 !      j = 2 refractivity at that point
 !RefracProfile(:,1) = slant%SlantSteps(1:slant%Nhor)
 RefracProfile(:,1) = slant%SlantSteps
 RefracProfile(:,2) = 0.0_wp

!write(*,*) 'k_0 = ', k, ' kmax = ', kmax, ' kmin = ', kmin
!!$
!!$ write(*,*) slant%SlantPointsEll(1,1), slant%SlantPointsEll(2,1), slant%SlantPointsEll(3,1)
!!$write(*,*) slant%SlantPointsEll(1,2), slant%SlantPointsEll(2,2), slant%SlantPointsEll(3,2)
!!$write(*,*) slant%SlantPointsEll(1,3), slant%SlantPointsEll(2,3), slant%SlantPointsEll(3,3)
!!$ write(*,*)
!!$ write(*,*)  'slant%SlantPointsEll(1,1:3)', slant%SlantPointsEll(1,1:3)
!!$ write(*,*)  'slant%SlantPointsEll(1:3,1)', slant%SlantPointsEll(1:3,1)

!!$

!!$ do i=1, slant%Nnghb
!!$    write(*,*) 'Col. ', i, ' k = ', k(i), ' : ', &
!!$                Hgeo(i,1,k(i)), ' < ', slant%SlantPointsEll(1,3), &
!!$                 ' < ', Hgeo(i,1,k(i)-1)
!!$ end do

 do j=1, slant%Nhor
    ! Process points on slant path. For each point one set of
    ! neighbored columns is provided.

    do i=1, slant%Nnghb
       ! Find the correct vertical indices for the next set of columns

       Hcol => Hgeo(i,j,:)
       k(i) = LayerSearch (Hcol, slant%SlantPointsEll(j,3))

!!$       if (k(i) .gt. kmin) then
!!$          write(*,*) 'Pt ', j, ' Col. ', i, ' k = ', k(i), ' : ', &
!!$               Hgeo(i,j,k(i)), ' < ', slant%SlantPointsEll(j,3), &
!!$               ' < ', Hgeo(i,j,k(i)-1)
!!$       else
!!$          write(*,*) 'Pt ', j, ' Col. ', i, ' k = ', k(i), ' : ', &
!!$               Hgeo(i,j,k(i))
!!$       end if

    end do  ! Find the correct vertical indices for the next set of columns

    ! Compute the partial pressure of water vapour e and the refractivity
    ! on the corners of the grid cell, i.e 8 corners for COSMO
    do i=1, slant%Nnghb
       ! neighbored columns ...

!!$       write(*,*) 'Pt ', j, ' Col. ', i, ' k = ', k(i), ' : ', &
!!$            Hgeo(i,j,k(i)), ' < ', slant%SlantPointsEll(j,3), &
!!$            ' < ', Hgeo(i,j,k(i)-1)
!!$       write(*,*) 'Druck ', col(i,j)%p(k(i)), ' > ', col(i,j)%p(k(i)-1)

       do m=0, 1
          ! upper (m=0) and lower (m=1) layer

          ! Vertical index:
          ! m=0 -  h = k(i)-1 - upper layer
          ! m=1 -  h = k(i)   - lower layer
          h = k(i) - 1 + m

          ! Partial pressure of water vapour e
          e(m,i) = (col(i,j)%p(h) * col(i,j)%q(h)) / (RDRD +   &
               col(i,j)%q(h)*EMRDRD)

          ! Refractivity N
          N(m,i) = NWein(col(i,j)%p(h), col(i,j)%t(h), e(m,i))
          if (debug) write(*,*) 'Refractivity N = ', m, i, N(m,i)

          ! Interpolate the refractivity at the given point using the
          ! refractivities at the corners of the grid cell
          ! node(:,1,:) - lower level - k+1
          ! node(:,2,:) - upper level - k
          node(i,2-m,1) = col(i,j)%dlon*deg2rad
          node(i,2-m,2) = col(i,j)%dlat*deg2rad
          node(i,2-m,3) = Hgeo(i,j,h)
          node(i,2-m,4) = N(m,i)

       end do
    end do

    ! Interpolation using function "BiLinearExp3D"
    if (slant%Nnghb .eq. 3) then
       ! GME - interpolate between 3 points
       node(4,:,:) = node(2,:,:)
    end if

    RefracProfile(j,2) = BiLinearExp3D(slant%SlantPointsEll(j,:),   &
                                        node, slant%Nnghb)
    !write(*,*) 'N = ', N(0,1), N(1,1), ' ipol = ', RefracProfile(j,2)



!    write(*,*) 'Ipol: 3dvar = ', Rfrc, ' bilin = ', RefracProfile(j,2), &
!               ' diff = ', Rfrc-RefracProfile(j,2), ' diff% = ', &
!                100.0*(Rfrc-RefracProfile(j,2))/Rfrc

 end do  ! do j=1, slant%Nhor

#ifndef __COSMO__
 do i=slant%Nhor+1, slant%Ntot

    ! .../3dvar/occ_lib/Earth.f90, used by MSIS
    ! Type Geodetic                ! Geodetic coordinates:
    !    Real(Double) :: H         ! Height above reference ellipsoid [km]
    !    Real(Double) :: Phi       ! Latitude from equator [degree]
    !    Real(Double) :: Lambda    ! Longitude [degree]
    ! End Type Geodetic

    ! Copy elliptical coordinates in "Geodetic" structure:
    EllCoord%Lambda = slant%SlantPointsEll(i,1) * rad2deg
    EllCoord%Phi    = slant%SlantPointsEll(i,2) * rad2deg
    EllCoord%H      = 0.0010_wp * slant%SlantPointsEll(i,3)


    !write(*,*) 'Geodetic ',  EllCoord%Lambda,  EllCoord%Phi, EllCoord%H
    ! WARNING: MSIS_Refractivity provides n-1, i. e. the refraction index
    !          minus one. To obtain the refractivity N it must be multiplied
    !          with 10^6 !!!
    call MSIS_Refractivity(EllCoord, Rfrc)
    RefracProfile(i,2) = Rfrc * 1.0E6_wp
    if (debug) write(*,*) 'MSIS Refractivity N = ', i, RefracProfile(i,2)

    !RefracProfile(i,2) = MSIS_Expand(EllCoord)
    !write(*,*) 'T, Tipol ', state%t(CellIndex(1)+1,   &
    !                                 CellIndex(2)+1,   &
    !                                 CellIndex(3),1),  &
    !                         RefracProfile(a,2)
    !write(*,*) 'top ', i,  RefracProfile(i,1), RefracProfile(i,2)

 end do
#endif


allocate (Slant%SlantRefrac(1:Slant%Ntot) )
Slant%SlantRefrac(:) = RefracProfile(:,2)

allocate (Slant%LineRefrac(1:Slant%Ntot) )
Slant%LineRefrac(:) = RefracProfile(:,2)

 ! Integrate along the slant path
!!$ ! 1) Use "simpned" to integrate the subpath inside and outside the model
!!$ call simpned (slant%SlantSteps(1:slant%Nhor), RefracProfile(1:slant%Nhor,2), &
!!$               slant%Nhor, STDmod )
!!$ call simpned (slant%SlantSteps(slant%Nhor:), RefracProfile(slant%Nhor:,2), &
!!$               slant%Ntot-slant%Nhor+1, STDtop )
!!$ delay =  (STDmod + STDtop)* 1.0D-6
!!$ ! 2) Use "simpned" to integrate the whole slant path
!!$ !    The sum of the subintegrals is in general not equal to the integral over
!!$ !    the  whole path !
!!$ call simpned (slant%SlantSteps, RefracProfile(:,2), &
!!$               slant%Ntot, delay2 )
!!$ delay2 = delay2*1.0D-6
!!$ ! 3) Use "IntegPolyCube" to integrate the subpath inside and outside the model
!!$ RP => RefracProfile(:slant%Nhor,:)
!!$ STD3 = IntegPolyCube(RP)
!!$ STD3 = STD3*1.0D-6
!!$ RP => RefracProfile(slant%Nhor:,:)
!!$ STD4 = IntegPolyCube(RP)
!!$ STD4 = STD4*1.0D-6
 ! 4) Use "simpned" to integrate the whole slant path

!!$write(*,*) 'RefracProfile(1,1) , RefracProfile(slant%Ntot,1), RefracProfile(1,2) , RefracProfile(slant%Ntot,2) = ', &
!!$     RefracProfile(1,1) , RefracProfile(slant%Ntot,1), &
!!$     RefracProfile(1,2) , RefracProfile(slant%Ntot,2)
!!$
!!$!write(*,*) 'GNSSFullProfile = ', RefracProfile
!!$
!!$do i=2, slant%Ntot
!!$   if (RefracProfile(i,1) < RefracProfile(i-1,1)) then
!!$      write(*,*) 'GNSSFullProfile: Invalid profile ', &
!!$           RefracProfile(i-1,1), RefracProfile(i,1)
!!$      exit
!!$   end if
!!$end do

 if (.false.) then
 !if (.true.) then
 ! write files with profiles
 !-------------------------------------------------------------------
 Nfile = Nfile + 1

#ifndef __COSMO__
 write(dace% pe*100+Nfile,*) '# Profile along signal path'
 write(dace% pe*100+Nfile,*) '# Coordinates on signal path :'
 write(dace% pe*100+Nfile,*) RefracProfile(:,1)
 write(dace% pe*100+Nfile,*) '# '
 write(dace% pe*100+Nfile,*) '# Refractivity on signal path :'
 write(dace% pe*100+Nfile,*) RefracProfile(:,2)
 CALL FLUSH(dace% pe*100+Nfile)
#endif

 !-------------------------------------------------------------------
 end if

 STD2 = IntegPolyCube(RefracProfile)
 !delay2 = IntegPolyCube(RefracProfile(1:slant%Nhor,:))
 STD2 = STD2*1.0E-6_wp
 !delay2 = delay2*1.0E-6_wp

 !write(*,*) 'std_delay_line> ZTD above model = ', STD2-delay2

 ! The delay is returned in two variables:
 delay     = STD2
 slant%STD = STD2

 Slant%Obs%mSTD    = slant%STD
 Slant%Obs%STDline = slant%STDline

! write(*,*) 'Slant Obs = ', slant%Obs%slant,              &
!            ' sim. = ', STDmod+STDtop,  &
!            ' Diff = ', slant%Obs%slant - STDmod - STDtop, &
!            ' model = ', STDmod,                          &
!            ' top = ', STDtop

if (debug) then
 write(*,*)  'std_delay_line: STDobs ',  &
        slant%Station%SName,     & ! station short name
        slant%Station%ID,        & ! station ID
        slant%Obs%time,          & ! observation time
        slant%Obs%satellite,     & ! observed satellite
        slant%Obs%elevation,     & ! elevation
        slant%Obs%azimuth,       & ! azimuth
        slant%Obs%zslant,        & ! observed STD
        slant%Obs%slant,         & ! STD mapped to zenith
        delay,                   & ! model STD
        slant%Obs%slant - delay  ! Diff obs - model
end if

!!$ write(*,*)
!!$ write(*,*) 'refraktivity profile (org):', slant%Nhor
!!$ do i=1, slant%Ntot
!!$    write(*,'(i4,tr2,f14.3,tr2,f13.8)') i, RefracProfile(i,:)
!!$ end do
!!$ write(*,*)

! error code needs to be defined
err = 0

if (debug) write(*,*)  'std_delay_line> STD computation finished ...'

!=======================================
! adjoint code follows (derives Jakobian
!=======================================

 !ladj = .false.
  if (ladj) then

     if (debug) write(*,*)  'std_delay_line> start of adjoint code ...'

     allocate( RefracProfile_ad(1:slant%Ntot,1:2) )
     RefracProfile_ad = 0.0_wp

     allocate( e_ad(0:1,1:slant%Nnghb) )
     allocate( N_ad(0:1,1:slant%Nnghb) )

     ! Adjoint code

     ! TEST TEST:
!    delay_ad = col_ad(1,1)%p(1)
!    col_ad(1,1)%p(1) = 0.0_wp
     delay_ad = 1.0_wp

     ! delay = STD2
     ! STD2 = STD2*1.0D-6
     STD2_ad = delay_ad * 1.0E-6_wp

     ! STD2 = IntegPolyCube(RefracProfile)
     call IntegPolyCube_ad (RefracProfile, RefracProfile_ad, STD2_ad)

     do j=slant%Nhor, 1, -1
        ! Process points on slant path. For each point one set of
        ! neighbored columns is provided.

        ! Recompute the required vertical indices
        do i=1, slant%Nnghb
           ! Find the correct vertical indices for the next set of columns
           Hcol => Hgeo(i,j,:)
           k(i) = LayerSearch (Hcol, slant%SlantPointsEll(j,3))
        end do

        ! Recompute refractivities
        ! Compute the partial pressure of water vapour e and the refractivity
        ! on the corners of the grid cell, i.e 8 corners for COSMO
        do i=1, slant%Nnghb
           ! neighbored columns ...

           do m=0, 1
              ! upper and lower layer

              ! Vertical index:
              ! m=0 -  h = k(i)-1 - upper layer
              ! m=1 -  h = k(i)   - lower layer
              h = k(i) - 1 + m

              ! Partial pressure of water vapour e
              e(m,i) = (col(i,j)%p(h) * col(i,j)%q(h)) /   &
                       (RDRD + col(i,j)%q(h)*EMRDRD)

              ! Refractivity N
              N(m,i) = NWein(col(i,j)%p(h), col(i,j)%t(h), e(m,i))

              ! Interpolate the refractivity at the given point using the
              ! refractivities at the corners of the grid cell
              ! node(:,1,:) - lower level - k+1
              ! node(:,2,:) - upper level - k
              node(i,2-m,1) = col(i,j)%dlon*deg2rad
              node(i,2-m,2) = col(i,j)%dlat*deg2rad
              node(i,2-m,3) = Hgeo(i,j,h)
              node(i,2-m,4) = N(m,i)

           end do
        end do

        N_ad = 0.0_wp

        ! Interpolation using function "BiLinearExp3D"

        ! recompute ...
        if (slant%Nnghb .eq. 3) then
           ! GME - interpolate between 3 points
           node(4,:,:) = node(2,:,:)
        end if

        ! RefracProfile(j,2) = BiLinearExp3D(slant%SlantPointsEll(j,:),   &
        !                                    node, slant%Nnghb)
        call BiLinearExp3D_ad(slant%SlantPointsEll(j,:), node, RefPt_ad,   &
                               node_ad, RefracProfile_ad(j,2), slant%Nnghb)

        if (slant%Nnghb .eq. 3) then
           ! GME - interpolate between 3 points
           ! node(4,:,:) = node(2,:,:)
           node_ad(2,:,:) = node_ad(2,:,:) + node_ad(4,:,:)
        end if

        ! Compute the partial pressure of water vapour e and the refractivity
        ! on the corners of the grid cell, i.e 8 corners for COSMO
        do i=slant%Nnghb, 1, -1
           ! neighbored columns ...

           do m=1, 0, -1
              ! upper and lower layer

              ! Vertical index:
              ! m=0 -  h = k(i)-1 - upper layer
              ! m=1 -  h = k(i)   - lower layer
              h = k(i) - 1 + m

              ! recompute e
              ! Partial pressure of water vapour e
              e(m,i) = (col(i,j)%p(h) * col(i,j)%q(h)) /    &
                   (RDRD + col(i,j)%q(h)*EMRDRD)


              ! Interpolate the refractivity at the given point using the
              ! refractivities at the corners of the grid cell
              ! node(:,1,:) - lower level - k+1
              ! node(:,2,:) - upper level - k
              !node(i,2-m,1) = col(i,j)%dlon*deg2rad
              !node(i,2-m,2) = col(i,j)%dlat*deg2rad
              !node(i,2-m,3) = Hgeo(i,j,h)

              ! node(i,2-m,4) = N(m,i)
              N_ad(m,i) = node_ad(i,2-m,4)

              ! Refractivity N
              ! N(m,i) = NWein(col(i,j)%p(k(i)+m), col(i,j)%t(k(i)+m), e(m,i))
              !call NWein_ad(col(i,j)%p(k(i)+m), col(i,j)%t(k(i)+m), e(m,i), &
              !     col_ad(i,j)%p(k(i)+m), col_ad(i,j)%t(k(i)+m), e_ad(m,i), &
              !     N_ad(m,i) )
              call NWein_ad(col(i,j)%p(h), col(i,j)%t(h), e(m,i), &
                            p_ad, t_ad, e_ad(m,i), N_ad(m,i) )
              col_ad(i,j)%p(h) = col_ad(i,j)%p(h) + p_ad
              col_ad(i,j)%t(h) = col_ad(i,j)%t(h) + t_ad

              ! Partial pressure of water vapour e
              !e(m,i) = (col(i,j)%p(k(i)+m) * col(i,j)%q(k(i)+m)) / (RDRD&
              !     col(i,j)%q(k(i)+m)*EMRDRD)
              col_ad(i,j)%p(h) =  col_ad(i,j)%p(h) +                &
                                       ( col(i,j)%q(h) / (RDRD +   &
                              col(i,j)%q(h)*EMRDRD) ) * e_ad(m,i)
              col_ad(i,j)%q(h) =  col_ad(i,j)%q(h) +                &
                ((( RDRD + col(i,j)%q(h) * EMRDRD ) *          &
                      col(i,j)%p(h) -                                    &
                      col(i,j)%p(h) * col(i,j)%q(h) *  EMRDRD)  &
                      / (RDRD + col(i,j)%q(h)*EMRDRD)**2) *    &
                      e_ad(m,i)

           end do
        end do

     end do  ! do j=slant%Nhor, 1, -1

     if (debug) write(*,*)  'std_delay_line> ... end of adjoint code'
  endif        ! adjoint




  !testad = .true.
  testad = .false.
  if (debug) write(*,*)  'std_delay_line> start test of adjoint?', testad

  if (testad) then
     ! Test tangent-linear and adjoint code

     if (debug) write(*,*)  'std_delay_line> start test of adjoint ...'

     ! allocate  col_tl and Dcol
     allocate( col_tl(lbound(col,1):ubound(col,1), lbound(col,2):ubound(col,2)) )
     allocate( coltest_ad(lbound(col,1):ubound(col,1),     &
                          lbound(col,2):ubound(col,2)) )
     allocate( Dcol(lbound(col,1):ubound(col,1), lbound(col,2):ubound(col,2)) )
     do i=lbound(col,1), ubound(col,1)
        do j=lbound(col,2), ubound(col,2)
           a = lbound(col(i,j)%p,1)
           b = ubound(col(i,j)%p,1)
           allocate( col_tl(i,j)%p(a:b) )
           allocate( col_tl(i,j)%q(a:b) )
           allocate( col_tl(i,j)%t(a:b) )
           col_tl(i,j)%p = 0.0_wp
           col_tl(i,j)%q = 0.0_wp
           col_tl(i,j)%t = 0.0_wp
           allocate( coltest_ad(i,j)%p(a:b) )
           allocate( coltest_ad(i,j)%q(a:b) )
           allocate( coltest_ad(i,j)%t(a:b) )
           coltest_ad(i,j)%p = 0.0_wp
           coltest_ad(i,j)%q = 0.0_wp
           coltest_ad(i,j)%t = 0.0_wp
           allocate( Dcol(i,j)%p(a:b) )
           allocate( Dcol(i,j)%q(a:b) )
           allocate( Dcol(i,j)%t(a:b) )
!!$           col_tl(i,j)%gpm => col(i,j)%gpm
!!$           col_tl(i,j)%geoid = col(i,j)%geoid
!!$           col_tl(i,j)%dlat = col(i,j)%dlat
!!$           col_tl(i,j)%dlon = col(i,j)%dlon
           Dcol(i,j)%gpm => col(i,j)%gpm
           Dcol(i,j)%geoid = col(i,j)%geoid
           Dcol(i,j)%dlat = col(i,j)%dlat
           Dcol(i,j)%dlon = col(i,j)%dlon
        end do
     end do

     ! Epsilon fuer den Differenzenquotienten
     epsilon = 5.0E-6_wp

     ! TL Inkremente
     do i=lbound(col,1), ubound(col,1)
        do j=lbound(col,2), ubound(col,2)
           col_tl(i,j)%p(:) = 500.0_wp
           col_tl(i,j)%q(:) = 0.0010_wp
           col_tl(i,j)%t(:) = 0.10_wp
        end do
     end do

     ! Inkrement
     do i=lbound(col,1), ubound(col,1)
        do j=lbound(col,2), ubound(col,2)
           Dcol(i,j)%p(:) = col(i,j)%p(:) +  epsilon*col_tl(i,j)%p(:)
           Dcol(i,j)%q(:) = col(i,j)%q(:) +  epsilon*col_tl(i,j)%q(:)
           Dcol(i,j)%t(:) = col(i,j)%t(:) +  epsilon*col_tl(i,j)%t(:)
        end do
     end do

     i = ubound(col,1) - 1
     j = lbound(col,2) + 1
     m = b-1
!!$     write(*,*) 'i, j, k lb, ub ', i, j, m, lbound(col(i,j)%p,1), ubound(col(i,j)%p,1)
!!$     write(*,*) 'p, p_tl, DeltaP ', col(i,j)%p(m), col_tl(i,j)%p(m),    &
!!$                                    Dcol(i,j)%p(m)
!!$     write(*,*) 'q, q_tl, DeltaQ ', col(i,j)%q(m), col_tl(i,j)%q(m),    &
!!$                                    Dcol(i,j)%q(m)
!!$     write(*,*) 'T, T_tl, DeltaT ', col(i,j)%t(m), col_tl(i,j)%t(m),    &
!!$                                    Dcol(i,j)%t(m)

     call  std_delay_line_tl (slant, Dcol, Ddelay, col_tl, delay_tl)
     ! The 2. call is necessary to obtain the correct delay_tl !!!!!!!
     call  std_delay_line_tl (slant, col, delay, col_tl, delay_tl)

     Diff = (Ddelay - delay) / epsilon

     !write(*,*) 'Test tangent-linear code <=> finite differences'
     !write(*,*) 'delay, Ddelay : ', delay,  Ddelay
     !write(*,*) 'Diff, delay_tl, Delta :', Diff, delay_tl, Diff-delay_tl

     !write(*,*) 'Test adjoint code'
     !write(*,*) 'Array product test:'
!!$     sumx = 0.0_wp
!!$     do i=lbound(col,1), ubound(col,1)
!!$        do j=lbound(col,2), ubound(col,2)
!!$           sumx = sumx + dot_product(coltest_ad(i,j)%p(:), col_tl(i,j)%p(:))
!!$           sumx = sumx + dot_product(coltest_ad(i,j)%q(:), col_tl(i,j)%q(:))
!!$           sumx = sumx + dot_product(coltest_ad(i,j)%t(:), col_tl(i,j)%t(:))
!!$        end do
!!$     end do
     !write(*,*) 'sumx = 0, sumx = ', sumx

!!$     write(*,*) 'coltest_ad(1,1)%p(:) = ',   coltest_ad(1,1)%p(:)
!!$     write(*,*) 'col_tl(1,1)%p(:) = ',   col_tl(1,1)%p(:)
!!$     write(*,*) ' dot_product =  ', dot_product(coltest_ad(1,1)%p(:), col_tl(1,1)%p(:))

     delay_ad = delay_tl
     call std_delay_line_ad (slant, col, delay_ad, coltest_ad)

     sumy =  delay_tl**2
     sumx = 0.0_wp
     do i=lbound(col,1), ubound(col,1)
        do j=lbound(col,2), ubound(col,2)
           a = lbound(col(i,j)%p,1)
           sumx = sumx + dot_product(coltest_ad(i,j)%p, col_tl(i,j)%p)
           sumx = sumx + dot_product(coltest_ad(i,j)%q, col_tl(i,j)%q)
           sumx = sumx + dot_product(coltest_ad(i,j)%t, col_tl(i,j)%t)
        end do
     end do

     !write(*,*) 'Norm X = ', sumx
     !write(*,*) 'Norm Y = ', sumy
     !write(*,*) 'Norm X, Norm Y,  Norm X / Norm Y = ', sumx, sumy, sumx/sumy
     if (abs(1.0_wp - sumx/sumy) .gt. 1E-8_wp) then
        write(*,*) 'Adjoint test failed !!!!', char(10), &
              'STD observation: ',  char(10), &
              slant%Station%SName,     & ! station short name
              slant%Station%ID,        & ! station ID
              slant%Obs%time,          & ! observation time
              slant%Obs%satellite,  char(10),     & ! observed satellite
              slant%Obs%elevation*rad2deg,     & ! elevation
              slant%Obs%azimuth*rad2deg,       & ! azimuth
              slant%Obs%zslant,        & ! observed STD
              slant%Obs%slant, char(10),   & ! STD mapped to zenith
              'Test tangent-linear code <=> finite differences',  char(10), &
              'delay, Ddelay : ', delay,  Ddelay,  char(10), &
              'Diff, delay_tl, Delta :', Diff, delay_tl, Diff-delay_tl, &
                                                              char(10), &
              'Test adjoint code',  char(10), &
              'Norm X, Norm Y,  Norm X / Norm Y = ', sumx, sumy, sumx/sumy, &
               char(10)
     else
        write(*,*) 'Adjoint increments',  char(10), &
             'Delay and its adjoint ', delay, delay_ad
     end if

     deallocate( col_tl )
     deallocate( coltest_ad )
     deallocate( Dcol   )

     if (debug) write(*,*)  'std_delay_line> ... test of adjoint finished'

  end if

 if (debug) write(*,*)  'std_delay_line> deallocate arrays ...'

 if (associated(Hgeo)) deallocate(Hgeo)
 if (associated(RefracProfile)) deallocate(RefracProfile)
 if (associated(k)) deallocate(k)
 if (associated(e)) deallocate(e)
 if (associated(N)) deallocate(N)

 if (associated(e_ad)) deallocate(e_ad)
 if (associated(N_ad)) deallocate(N_ad)
 if (associated(RefracProfile_ad)) deallocate(RefracProfile_ad)

 if (debug) write(*,*)  'std_delay_line> ... end'

 end subroutine std_delay_line


subroutine std_delay_line_tl (slant, col, delay, col_tl, delay_tl)
 !-------------------------------------------------------------------------
 ! Compute the slant total delay
 ! as the original std_delay but interface adapted to the 3dvar environment
 ! +++ work in progress +++
 !
 !
 ! col  -  vertical grid columns
 !         col(i,j), i=1, Nnghb, numberof neighbored columns
 !                               Nnghb = 3 - GME
 !                               Nnghb = 4 - COSMO
 !                   j=1, Nhor, number of slant points within grid
 !-------------------------------------------------------------------------
 type (SlantData) ,intent(in)    :: slant        ! slant meta data
 type (p_column)  ,intent(in)    :: col (:,:)    ! model columns     (input)
 real(wp)         ,intent(out)   :: delay        ! Slant Total Delay (output)

 type (p_column)  ,pointer       :: col_tl (:,:) ! tangent_linear (input)
 real(wp)         ,intent(out)   :: delay_tl     ! Slant Total Delay, TL


 real(wp), dimension(:,:,:), pointer :: Hgeo => Null()
 real(wp), dimension(:), pointer     :: Hcol => Null()
 real(wp), dimension(:,:), pointer   :: RefracProfile => Null()
!real(wp), dimension(:,:), pointer   :: RP => Null()
 integer  :: i, j, m, h !, kmax, kmin

 ! vertical grid index of neighbored columns:
 integer, dimension(:), pointer  :: k => Null()
 !integer, dimension(:), pointer  :: k2 => Null()
 real(wp), dimension(:,:), pointer :: e => Null()
 real(wp), dimension(:,:), pointer :: N => Null()
 real (wp), dimension(1:4,1:2,1:4)  :: node

 real(wp), dimension(:,:), pointer :: e_tl => Null()
 real(wp), dimension(:,:), pointer :: N_tl => Null()
 real(wp), dimension(1:4,1:2,1:4)  :: node_tl
 real(wp), dimension(:,:), pointer   :: RefracProfile_tl => Null()

!real(wp), dimension(1:2,1:2)    :: vnode_tl
!real(dp)                        :: Rfrc_tl, RPt_tl
 real(dp)                        :: STD2_tl
 real(wp), dimension(1:3)        :: RefPt_tl

 real(wp)                       :: STD2         !, STD, delay2
!real(wp)                       :: STD3, STD4
!real(wp)                       :: STDmod, STDtop
 real(dp)                       :: Rfrc

#ifndef __COSMO__
 type (Geodetic) :: EllCoord
#endif

!real(wp), dimension(1:2,1:2)    :: vnode

 logical :: debug = .false.
 ! .....

 !debug = .true.
 if (debug) write(*,*)  'std_delay_line_tl> Start ...'
 if (debug .and. slant%Nnghb .eq. 3) write(*,*) 'GME'
 if (debug .and. slant%Nnghb .eq. 4) write(*,*) 'COSMO'

 !write(*,*) 'Dimension Columns ', shape(col)
 !write(*,*) 'Punkte auf Slant ', slant%Ntot
 !write(*,*) 'Punkte im Gitter :', slant%Nmod, slant%Nhor
 !write(*,*) 'Nachbarzellen    :', slant%Nnghb, slant%Naloc, slant%Nnghb

!!$     i = ubound(col,1) - 1
!!$     j = lbound(col,2) + 1
!!$     m = ubound(col(i,j)%p,1)-1
!!$     write(*,*) 'i, j, k lb, ub ', i, j, m, lbound(col(i,j)%p,1), ubound(col(i,j)%p,1)
!!$     write(*,*) 'p, p_tl, DeltaP ', col(i,j)%p(m), col_tl(i,j)%p(m)
!!$     write(*,*) 'q, q_tl, DeltaQ ', col(i,j)%q(m), col_tl(i,j)%q(m)
!!$     write(*,*) 'T, T_tl, DeltaT ', col(i,j)%t(m), col_tl(i,j)%t(m)



 ! Apply geoid corrections
 allocate(Hgeo( lbound(col,1):ubound(col,1),                       &
                lbound(col,2):ubound(col,2),                       &
                lbound(col(1,1)%gpm,1):ubound(col(1,1)%gpm,1) ) )
 !write(*,*) 'Dimension Hgeo ', shape(Hgeo)

 do i=lbound(col,1), ubound(col,1)
    do j=lbound(col,2), ubound(col,2)
       Hgeo(i,j,:) = col(i,j)%gpm(:) + col(i,j)%geoid
    end do
 end do

 allocate( RefracProfile(1:slant%Ntot,1:2) )
 allocate( RefracProfile_tl(1:slant%Ntot,1:2) )
 allocate( k(1:slant%Nnghb) )

 allocate( e(0:1,1:slant%Nnghb) )
 allocate( e_tl(0:1,1:slant%Nnghb) )

 allocate( N(0:1,1:slant%Nnghb) )
 allocate( N_tl(0:1,1:slant%Nnghb) )

 ! RefracProfile(i,j)
 !  i - supporting points on the slant path
 !  j - j = 1 distance on slant path (x coordinate)
 !            as given in slant%SlantSteps
 !      j = 2 refractivity at that point
 !RefracProfile(:,1) = slant%SlantSteps(1:slant%Nhor)
 RefracProfile(:,1) = slant%SlantSteps
 RefracProfile(:,2) = 0.0_wp
 RefracProfile_tl(:,1) =  0.0_wp  ! slant%SlantSteps
 RefracProfile_tl(:,2) = 0.0_wp

 do j=1, slant%Nhor
    ! Process points on slant path. For each point one set of
    ! neighbored columns is provided.

    do i=1, slant%Nnghb
       ! Find the correct vertical indices for the next set of columns
       Hcol => Hgeo(i,j,:)
       k(i) = LayerSearch(Hcol, slant%SlantPointsEll(j,3))
    end do

    ! Compute the partial pressure of water vapour e and the refractivity
    ! on the corners of the grid cell, i.e 8 corners for COSMO
    do i=1, slant%Nnghb
       ! neighbored columns ...

       do m=0, 1
          ! upper and lower layer

          ! Vertical index:
          ! m=0 -  h = k(i)-1 - upper layer
          ! m=1 -  h = k(i)   - lower layer
          h = k(i) - 1 + m

          ! Partial pressure of water vapour e
          e(m,i) = (col(i,j)%p(h) * col(i,j)%q(h)) / (RDRD +   &
               col(i,j)%q(h)*EMRDRD)
          e_tl(m,i) = ((( RDRD + col(i,j)%q(h) * EMRDRD ) *     &
                      col(i,j)%p(h) -                                    &
                      col(i,j)%p(h) * col(i,j)%q(h) *  EMRDRD)  &
                      / (RDRD + col(i,j)%q(h)*EMRDRD)**2)      &
                      * col_tl(i,j)%q(h) +                               &
                      ((col(i,j)%q(h)) /                                 &
                      (RDRD + col(i,j)%q(h)*EMRDRD))           &
                      * col_tl(i,j)%p(h)

          ! Refractivity N
          N(m,i) = NWein(col(i,j)%p(h), col(i,j)%t(h), e(m,i))
          N_tl(m,i) = NWein_tl(col(i,j)%p(h), col(i,j)%t(h), e(m,i), &
                               col_tl(i,j)%p(h), col_tl(i,j)%t(h),   &
                               e_tl(m,i))

          !write(*,*) 'e, N ', e(m,i), N(m,i)

          ! Interpolate the refractivity at the given point using the
          ! refractivities at the corners of the grid cell
          ! node(:,1,:) - lower level - k+1
          ! node(:,2,:) - upper level - k
          node(i,2-m,1) = col(i,j)%dlon*deg2rad
          node(i,2-m,2) = col(i,j)%dlat*deg2rad
          node(i,2-m,3) = Hgeo(i,j,h)
          node(i,2-m,4) = N(m,i)

          node_tl(i,2-m,1) = 0.0_wp  ! col(i,j)%dlon*deg2rad
          node_tl(i,2-m,2) = 0.0_wp  ! col(i,j)%dlat*deg2rad
          node_tl(i,2-m,3) = 0.0_wp  ! Hgeo(i,j,k(i)+m)
          node_tl(i,2-m,4) = N_tl(m,i)

       end do
    end do

    ! Interpolation using function "BiLinearExp3D"
    if (slant%Nnghb .eq. 3) then
       ! GME - interpolate between 3 points
       node(4,:,:) = node(2,:,:)
       node_tl(4,:,:) = node_tl(2,:,:)
    end if

    RefracProfile(j,2) = BiLinearExp3D(slant%SlantPointsEll(j,:),   &
                                        node, slant%Nnghb)
    RefPt_tl = 0.0_wp
    RefracProfile_tl(j,2) = BiLinearExp3D_tl(slant%SlantPointsEll(j,:),   &
                                    node,  RefPt_tl, node_tl, slant%Nnghb)

    !write(*,*) 'N = ', N(0,1), N(1,1), ' ipol = ', RefracProfile(j,2)



!    write(*,*) 'Ipol: 3dvar = ', Rfrc, ' bilin = ', RefracProfile(j,2), &
!               ' diff = ', Rfrc-RefracProfile(j,2), ' diff% = ', &
!                100.0*(Rfrc-RefracProfile(j,2))/Rfrc

 end do  ! do j=1, slant%Nhor

#ifndef __COSMO__
 do i=slant%Nhor+1, slant%Ntot

    ! .../3dvar/occ_lib/Earth.f90, used by MSIS
    ! Type Geodetic                ! Geodetic coordinates:
    !    Real(Double) :: H         ! Height above reference ellipsoid [km]
    !    Real(Double) :: Phi       ! Latitude from equator [degree]
    !    Real(Double) :: Lambda    ! Longitude [degree]
    ! End Type Geodetic

    ! Copy elliptical coordinates in "Geodetic" structure:
    EllCoord%Lambda = slant%SlantPointsEll(i,1) * rad2deg
    EllCoord%Phi    = slant%SlantPointsEll(i,2) * rad2deg
    EllCoord%H      = 0.0010_wp * slant%SlantPointsEll(i,3)


    !write(*,*) 'Geodetic ',  EllCoord%Lambda,  EllCoord%Phi, EllCoord%H
    ! WARNING: MSIS_Refractivity provides n-1, i. e. the refraction index
    !          minus one. To obtain the refractivity N it must be multiplied
    !          with 10^6 !!!
    call MSIS_Refractivity(EllCoord, Rfrc)
    RefracProfile(i,2) = Rfrc * 1.0E6_wp
    !RefracProfile(i,2) = MSIS_Expand(EllCoord)
    !write(*,*) 'T, Tipol ', state%t(CellIndex(1)+1,   &
    !                                 CellIndex(2)+1,   &
    !                                 CellIndex(3),1),  &
    !                         RefracProfile(a,2)
    !write(*,*) 'top ', i,  RefracProfile(i,1), RefracProfile(i,2)

 end do
#endif

 ! 4) Use "simpned" to integrate the whole slant path
 STD2 = IntegPolyCube(RefracProfile)
 STD2 = STD2*1.0E-6_wp

 delay = STD2

 STD2_tl = IntegPolyCube_tl(RefracProfile, RefracProfile_tl)
 STD2_tl = STD2_tl*1.0E-6_wp
 delay_tl = STD2_tl

 i = 0
 do j=lbound(RefracProfile_tl,1), ubound(RefracProfile_tl,1)
    if (RefracProfile_tl(j,2) /= RefracProfile_tl(j,2)) then
       write(*,*) 'std_delay_line_tl, slant%Nhor, j, ', slant%Nhor, j, RefracProfile_tl(j,:)
       i = -999
       exit
    end if
 end do

!!$ write(*,*) 'std_delay_line_tl', char(10), &
!!$      'STD observation: ',  char(10), &
!!$      slant%Station%SName,     & ! station short name
!!$      slant%Station%ID,        & ! station ID
!!$      slant%Obs%time,          & ! observation time
!!$      slant%Obs%satellite,  char(10),     & ! observed satellite
!!$      slant%Obs%elevation*rad2deg,     & ! elevation
!!$      slant%Obs%azimuth*rad2deg,       & ! azimuth
!!$      slant%Obs%zslant,        & ! observed STD
!!$      slant%Obs%slant, char(10),   & ! STD mapped to zenith
!!$      'STD2_tl = ', STD2_tl, char(10),   &
!!$      'max RefracProfile_tl = ', maxval(RefracProfile_tl),  char(10),   &
!!$      'RefracProfile_tl contains NaN if i = -999 : i=', i


!!$ write(*,*)
!!$ write(*,*) 'refraktivity profile (tl):', slant%Nhor
!!$ do i=1, slant%Ntot
!!$    write(*,'(i4,tr2,f14.3,tr2,f13.8)') i, RefracProfile(i,:)
!!$ end do
!!$ write(*,*)

 if (associated(Hgeo)) deallocate(Hgeo)
 if (associated(RefracProfile)) deallocate(RefracProfile)
 if (associated(k)) deallocate(k)
 if (associated(e)) deallocate(e)
 if (associated(N)) deallocate(N)

 if (associated(e_tl)) deallocate(e_tl)
 if (associated(N_tl)) deallocate(N_tl)
 if (associated(RefracProfile_tl)) deallocate(RefracProfile_tl)

 if (debug) write(*,*)  'std_delay_line_tl> ... end'

 end subroutine std_delay_line_tl


subroutine std_delay_line_ad (slant, col, delay_ad, col_ad)
 !-------------------------------------------------------------------------
 ! Compute the slant total delay
 ! as the original std_delay but interface adapted to the 3dvar environment
 ! +++ work in progress +++
 !
 !
 ! col  -  vertical grid columns
 !         col(i,j), i=1, Nnghb, numberof neighbored columns
 !                               Nnghb = 3 - GME
 !                               Nnghb = 4 - COSMO
 !                   j=1, Nhor, number of slant points within grid
 !-------------------------------------------------------------------------
 type (SlantData) ,intent(in)    :: slant        ! slant meta data
 type (p_column)  ,intent(in)    :: col (:,:)    ! model columns     (input)

 type (p_column)  ,pointer       :: col_ad (:,:) ! adjoint model vars. (output)
 real(wp)         ,intent(in)    :: delay_ad     ! STD, adjoint (input)


 real(wp), dimension(:,:,:), pointer :: Hgeo => Null()
 real(wp), dimension(:), pointer     :: Hcol => Null()
 real(wp), dimension(:,:), pointer   :: RefracProfile => Null()
!real(wp), dimension(:,:), pointer   :: RP => Null()
 integer  :: i, j, m, h !, kmax, kmin

 ! vertical grid index of neighbored columns:
 integer, dimension(:), pointer  :: k => Null()
 real(wp), dimension(:,:), pointer :: e => Null()
 real(wp), dimension(:,:), pointer :: N => Null()
 real (wp), dimension(1:4,1:2,1:4)  :: node

 real(wp), dimension(:,:), pointer :: e_ad => Null()
 real(wp), dimension(:,:), pointer :: N_ad => Null()
 real(wp), dimension(1:4,1:2,1:4)  :: node_ad
 real(wp), dimension(:,:), pointer :: RefracProfile_ad => Null()

!real(wp), dimension(1:2,1:2)   :: vnode_ad
!real(dp)                       :: Rfrc_ad, RPt_ad
 real(dp)                      :: STD2_ad
 real(wp), dimension(1:3)          :: RefPt_ad
 real(dp)                      :: p_ad, t_ad

!real(wp)                       :: STD, STD2, delay2
!real(wp)                       :: STD3, STD4
!real(wp)                       :: STDmod, STDtop
 real(dp)                      :: Rfrc

#ifndef __COSMO__
 type (Geodetic) :: EllCoord
#endif

!real(wp), dimension(1:2,1:2)    :: vnode

 logical :: debug = .false.
 ! .....

 !debug = .true.
 if (debug) write(*,*)  'std_delay_line_ad> Start ...'
 if (debug .and. slant%Nnghb .eq. 3) write(*,*) 'GME'
 if (debug .and. slant%Nnghb .eq. 4) write(*,*) 'COSMO'

 ! Apply geoid corrections
 allocate(Hgeo( lbound(col,1):ubound(col,1),                       &
                lbound(col,2):ubound(col,2),                       &
                lbound(col(1,1)%gpm,1):ubound(col(1,1)%gpm,1) ) )

 do i=lbound(col,1), ubound(col,1)
    do j=lbound(col,2), ubound(col,2)
       Hgeo(i,j,:) = col(i,j)%gpm(:) + col(i,j)%geoid
    end do
 end do

 allocate( RefracProfile(1:slant%Ntot,1:2) )
 allocate( RefracProfile_ad(1:slant%Ntot,1:2) )

 allocate( k(1:slant%Nnghb) )

 allocate( e(0:1,1:slant%Nnghb) )
 allocate( e_ad(0:1,1:slant%Nnghb) )

 allocate( N(0:1,1:slant%Nnghb) )
 allocate( N_ad(0:1,1:slant%Nnghb) )

 ! RefracProfile(i,j)
 !  i - supporting points on the slant path
 !  j - j = 1 distance on slant path (x coordinate)
 !            as given in slant%SlantSteps
 !      j = 2 refractivity at that point
 RefracProfile(:,1) = slant%SlantSteps
 RefracProfile(:,2) = 0.0_wp
 RefracProfile_ad = 0.0_wp


 do j=1, slant%Nhor
    ! Process points on slant path. For each point one set of
    ! neighbored columns is provided.

    do i=1, slant%Nnghb
       ! Find the correct vertical indices for the next set of columns
       Hcol => Hgeo(i,j,:)
       k(i) = LayerSearch (Hcol, slant%SlantPointsEll(j,3))
    end do  ! Find the correct vertical indices for the next set of columns

    ! Compute the partial pressure of water vapour e and the refractivity
    ! on the corners of the grid cell, i.e 8 corners for COSMO
    do i=1, slant%Nnghb
       ! neighbored columns ...

       do m=0, 1
          ! upper and lower layer

          ! Vertical index:
          ! m=0 -  h = k(i)-1 - upper layer
          ! m=1 -  h = k(i)   - lower layer
          h = k(i) - 1 + m

          ! Partial pressure of water vapour e
          e(m,i) = (col(i,j)%p(h) * col(i,j)%q(h)) / (RDRD +   &
               col(i,j)%q(h)*EMRDRD)

          ! Refractivity N
          N(m,i) = NWein(col(i,j)%p(h), col(i,j)%t(h), e(m,i))

          ! Interpolate the refractivity at the given point using the
          ! refractivities at the corners of the grid cell
          ! node(:,1,:) - lower level - k+1
          ! node(:,2,:) - upper level - k
          node(i,2-m,1) = col(i,j)%dlon*deg2rad
          node(i,2-m,2) = col(i,j)%dlat*deg2rad
          node(i,2-m,3) = Hgeo(i,j,h)
          node(i,2-m,4) = N(m,i)

       end do
    end do

    ! Interpolation using function "BiLinearExp3D"
    if (slant%Nnghb .eq. 3) then
       ! GME - interpolate between 3 points
       node(4,:,:) = node(2,:,:)
    end if

    RefracProfile(j,2) = BiLinearExp3D(slant%SlantPointsEll(j,:),   &
                                        node, slant%Nnghb)
    !write(*,*) 'N = ', N(0,1), N(1,1), ' ipol = ', RefracProfile(j,2)


 end do  ! do j=1, slant%Nhor

#ifndef __COSMO__
 do i=slant%Nhor+1, slant%Ntot

    ! .../3dvar/occ_lib/Earth.f90, used by MSIS
    ! Type Geodetic                ! Geodetic coordinates:
    !    Real(Double) :: H         ! Height above reference ellipsoid [km]
    !    Real(Double) :: Phi       ! Latitude from equator [degree]
    !    Real(Double) :: Lambda    ! Longitude [degree]
    ! End Type Geodetic

    ! Copy elliptical coordinates in "Geodetic" structure:
    EllCoord%Lambda = slant%SlantPointsEll(i,1) * rad2deg
    EllCoord%Phi    = slant%SlantPointsEll(i,2) * rad2deg
    EllCoord%H      = 0.0010_wp * slant%SlantPointsEll(i,3)


    !write(*,*) 'Geodetic ',  EllCoord%Lambda,  EllCoord%Phi, EllCoord%H
    ! WARNING: MSIS_Refractivity provides n-1, i. e. the refraction index
    !          minus one. To obtain the refractivity N it must be multiplied
    !          with 10^6 !!!
    call MSIS_Refractivity(EllCoord, Rfrc)
    RefracProfile(i,2) = Rfrc * 1.0E6_wp

 end do
#endif

 ! ================== end of recomputation ==============================

 ! Adjoint code

 ! delay = STD2
 ! STD2 = STD2*1.0D-6
 STD2_ad = delay_ad * 1.0E-6_wp

 ! STD2 = IntegPolyCube(RefracProfile)
 call IntegPolyCube_ad (RefracProfile, RefracProfile_ad, STD2_ad)

 do j=slant%Nhor, 1, -1
    ! Process points on slant path. For each point one set of
    ! neighbored columns is provided.

    ! Recompute the required vertical indices
    do i=1, slant%Nnghb
       ! Find the correct vertical indices for the next set of columns
       Hcol => Hgeo(i,j,:)
       k(i) = LayerSearch (Hcol, slant%SlantPointsEll(j,3))
    end do  ! Find the correct vertical indices for the next set of columns

    ! Recompute refractivities
    ! Compute the partial pressure of water vapour e and the refractivity
    ! on the corners of the grid cell, i.e 8 corners for COSMO
    do i=1, slant%Nnghb
       ! neighbored columns ...

       do m=0, 1
          ! upper and lower layer

          ! Vertical index:
          ! m=0 -  h = k(i)-1 - upper layer
          ! m=1 -  h = k(i)   - lower layer
          h = k(i) - 1 + m

          ! Partial pressure of water vapour e
          e(m,i) = (col(i,j)%p(h) * col(i,j)%q(h)) / (RDRD +   &
               col(i,j)%q(h)*EMRDRD)

          ! Refractivity N
          N(m,i) = NWein(col(i,j)%p(h), col(i,j)%t(h), e(m,i))

          ! Interpolate the refractivity at the given point using the
          ! refractivities at the corners of the grid cell
          ! node(:,1,:) - lower level - k+1
          ! node(:,2,:) - upper level - k
          node(i,2-m,1) = col(i,j)%dlon*deg2rad
          node(i,2-m,2) = col(i,j)%dlat*deg2rad
          node(i,2-m,3) = Hgeo(i,j,h)
          node(i,2-m,4) = N(m,i)

       end do
    end do

    N_ad = 0.0_wp
    ! Interpolation using function "BiLinearExp3D"

    ! recompute ...
    if (slant%Nnghb .eq. 3) then
       ! GME - interpolate between 3 points
       node(4,:,:) = node(2,:,:)
    end if

    ! RefracProfile(j,2) = BiLinearExp3D(slant%SlantPointsEll(j,:),   &
    !                                    node, slant%Nnghb)
    call BiLinearExp3D_ad(slant%SlantPointsEll(j,:), node, RefPt_ad,   &
                             node_ad, RefracProfile_ad(j,2), slant%Nnghb)

    if (slant%Nnghb .eq. 3) then
       ! GME - interpolate between 3 points
       ! node(4,:,:) = node(2,:,:)
       node_ad(2,:,:) = node_ad(2,:,:) + node_ad(4,:,:)
    end if


    ! Compute the partial pressure of water vapour e and the refractivity
    ! on the corners of the grid cell, i.e 8 corners for COSMO
    do i=slant%Nnghb, 1, -1
       ! neighbored columns ...

       do m=1, 0, -1
          ! upper and lower layer

          ! Vertical index:
          ! m=0 -  h = k(i)-1 - upper layer
          ! m=1 -  h = k(i)   - lower layer
          h = k(i) - 1 + m

          ! recompute e
          ! Partial pressure of water vapour e
          e(m,i) = (col(i,j)%p(h) * col(i,j)%q(h)) / (RDRD +   &
               col(i,j)%q(h)*EMRDRD)

          ! Interpolate the refractivity at the given point using the
          ! refractivities at the corners of the grid cell
          ! node(:,1,:) - lower level - k+1
          ! node(:,2,:) - upper level - k
          !node(i,2-m,1) = col(i,j)%dlon*deg2rad
          !node(i,2-m,2) = col(i,j)%dlat*deg2rad
          !node(i,2-m,3) = Hgeo(i,j,h)

          ! node(i,2-m,4) = N(m,i)
          N_ad(m,i) = node_ad(i,2-m,4)

          ! Refractivity N
          ! N(m,i) = NWein(col(i,j)%p(k(i)+m), col(i,j)%t(k(i)+m), e(m,i))
          !call NWein_ad(col(i,j)%p(k(i)+m), col(i,j)%t(k(i)+m), e(m,i), & !P,T,E
          !     col_ad(i,j)%p(k(i)+m), col_ad(i,j)%t(k(i)+m), e_ad(m,i), & !..ad
          !     N_ad(m,i) )
          call NWein_ad(col(i,j)%p(h), col(i,j)%t(h), e(m,i), & !P,T,E
               p_ad, t_ad, e_ad(m,i), N_ad(m,i) )                         ! ..ad
          col_ad(i,j)%p(h) = col_ad(i,j)%p(h) + p_ad
          col_ad(i,j)%t(h) = col_ad(i,j)%t(h) + t_ad

          ! Partial pressure of water vapour e
          !e(m,i) = (col(i,j)%p(k(i)+m) * col(i,j)%q(k(i)+m)) / (RDRD +   &
          !     col(i,j)%q(k(i)+m)*EMRDRD)
          col_ad(i,j)%p(h) =  col_ad(i,j)%p(h) +                    &
                 ( col(i,j)%q(h) / (RDRD +                         &
                    col(i,j)%q(h)*EMRDRD) ) * e_ad(m,i)
          col_ad(i,j)%q(h) =  col_ad(i,j)%q(h) +                    &
                ((( RDRD + col(i,j)%q(h) * EMRDRD ) *          &
                      col(i,j)%p(h) -                                    &
                      col(i,j)%p(h) * col(i,j)%q(h) *  EMRDRD)  &
                      / (RDRD + col(i,j)%q(h)*EMRDRD)**2) *    &
                      e_ad(m,i)

       end do
    end do

 end do  ! do j=slant%Nhor, 1, -1

 if (associated(Hgeo)) deallocate(Hgeo)
 if (associated(RefracProfile)) deallocate(RefracProfile)
 if (associated(k)) deallocate(k)
 if (associated(e)) deallocate(e)
 if (associated(N)) deallocate(N)

 if (associated(e_ad)) deallocate(e_ad)
 if (associated(N_ad)) deallocate(N_ad)
 if (associated(RefracProfile_ad)) deallocate(RefracProfile_ad)

 if (debug) write(*,*)  'std_delay_line_ad> ... end'

 end subroutine std_delay_line_ad

!==============================================================================
! End of old code without raytracer but with tangent-linear and adjoint
!==============================================================================




!==============================================================================
end module mo_std_operator
