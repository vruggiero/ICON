!+ Utility routines for setting vertical grid and reference atmosphere(s)
!==============================================================================

MODULE vgrid_refatm_utils

!==============================================================================
!
! Description:
!  This module provides routines that have to deal with the vertical grid
!  and the reference atmosphere.
!  Especially when using GRIB2 and the general vertical coordinate, some
!  measures have to be taken to make sure that the same reference atmosphere
!  parameters and vertical grid specifications are taken in INT2LM and the
!  COSMO-Model.
!
!  It also contains the utilities to generate a new UUID and to compare
!  different UUIDs.
!
!  This module is shared between the COSMO model and DACE program packages (!),
!  execpt that a few lines following 'IF (vc_type%vert_coord(k) >= 0.0_wp)'
!  are commented out in DACE.
!
!  Routines (module procedures) currently contained:
!
!     - reference_atmosphere
!       Computes the reference atmosphere and vertical coordinate
!       parameters.
!
!     - reference_atmosphere_2
!       Computes the reference atmosphere and vertical coordinate
!       parameters for the new reference atmosphere
!
!     - reference_atmosphere_BVconst
!       Computes the reference atmosphere and vertical coordinate
!       parameters on assuming a constant Brunt-Vaisala frequency
!
!     - k_index_of_pressure_levels
!       Determination of the k-indices for those full levels which
!       correspond to about 850, 800, 500, 400 and 300 hPa.
!
!    - set_vcoordtype
!       to set a couple of default sets for the vertical coordinate parameters
!
!    - set_refatmtype
!       to set a couple of default sets for the reference atmosphere values
!
!    the following routines do access C-routines, which are located in a
!    (new) module c_utilities.c. It is used by Fortran2003 C-bindings
!
!    - uuid_create
!       to get a new uuid
!
!    - uuid_2char
!       to transform a uuid into a readable character array
!
!    - uuid_match
!       to test if two uuids match
!
! Current Code Owner (for COSMO and for DACE): 
!  DWD, Ulrich Schaettler
!  phone:  +49  69  8062 2739
!  fax:    +49  69  8062 3721
!  email:  ulrich.schaettler@dwd.de
!
! History:
! Version      Date       Name
! ----------   ---------- ----
! V4_28        2013/07/12 Ulrich Schaettler
!  Initial release for INT2LM and the COSMO-Model
! V4_29        2013/10/04 Ulrich Schaettler
!  Remove suffix _out from variables vcoord, refatm
! V5_1         2014-11-28 Ulrich Schaettler, Oliver Fuhrer
!  Enlarged khmax to 250
!  Initialize vertical coordinate information with -1.0
!  Modified all SR reference_atmosphere_xx:
!   - only compute ak, bk, and sigm_coord, if a new hhl has to be computed
!     (because only then the necessary vert_coord are given)
!   - introduced yerrmsg, ierror to all these routines, to check hhl consistency
!  Introduced new ivctype=4: modified SLEVE coordinate (DL)
!  Replaced ireals by wp (working precision) (OF)
! V5_2         2015-05-21 Ulrich Blahak, Ulrich Schaettler
!  Removed computation of vc_type%kflat from the reference_atmosphere_xxx() routines
!   in case of lnew_hhl=.FALSE., because the previous computation failed for
!   PE subdomains with entirely flat terrain (e.g. over sea) and the correct
!   computation involves a global communication. In this case, the computation has to be done
!   afterwards by the calling program (UB).
!  Modified SR dealloc_refatm_defaults, dealloc_vcoord_defaults to deallocate these
!   structures at end of the program  (US).
! V5_5a        2018-06-22 Ulrich Schaettler
!  Removed all iintegers (US)
!  Use kind_parameters instead of data_parameters (US)
!  Moved KIND parameters for the GRIB library from data_parameters to data_io
! upcoming based on V5_6  Christoph Schraff
!  - Set subroutines and variables used by other modules explicitly to public
!    (required by ICON Compiler directives when modules are ported to ICON);
!  - Some unused variables are removed (e.g. irealgrib) or commented out;
!  - Trailing blanks at eol removed.
!
! ==============================================================================
! DACE history:
! ------------
!  1_44       2015-09-30 Andreas Rhodin : take over from COSMO 5.03-beta
!  1_46       2016-02-05 Andreas Rhodin : remove debugging printout
!  1_51       2017-02-24 Andreas Rhodin : update shared modules from COSMO 5.04d
! ==============================================================================
! CAUTION: This module is used by both the DACE and COSMO model codes.       !!!
!!!        Therefore, anybody wanting to introduce a modification to this    !!!
!!!        module in the context of either of these programs must consult    !!!
!!!        the 'current code owner' of this module for the other program,    !!!
!!!        in order to allow for checking that the modification will comply  !!!
!!!        with both program packages. This must be done before the          !!!
!!!        modification is put into the Version Control System (VCS).        !!!
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================

! Modules used:
USE kind_parameters , ONLY :   &
  wp,        & ! KIND-type parameter for real variables
  i8           ! KIND-type parameter for integer*8 variables

!==============================================================================

IMPLICIT NONE

PRIVATE
          ! public subroutines
PUBLIC :: set_vcoord_defaults,  set_refatm_defaults,                         &
          reference_atmosphere,  reference_atmosphere_2,                     &
          dealloc_vcoord_structure,  dealloc_refatm_structure

          ! public data types and variables
PUBLIC :: vcoord_type, refatm_type, vcoord_defaults, refatm_defaults

          ! public subroutines and function used in the COSMO model code only
PUBLIC :: reference_atmosphere_BVconst,  k_index_of_pressure_levels,         &
          uuid_create,  uuid_2char,  uuid_match

          ! public data types and variables used in the COSMO model code only
PUBLIC :: khmax, lanalyt_calc_t0p0, lhhl_in_read, lhhl_hasbeenread,          &
          nfltvc, svc1, svc2, imax_vcoordtype, imax_refatmtype, rundefined,  &
          vcoord, refatm,                                                    &
          uuid_string_length, uuid_out, uuid_in_string, uuid_out_string

!==============================================================================

! parameters for vertical grid and reference atmosphere (from data_grid_lm)
! Global Parameters
  INTEGER, PARAMETER     ::                       &
    khmax = 250         ! max. number of vertical coordinate parameters
                        ! according to grib.

! to compute values p0, t0 of reference atmosphere analytically on full levels
  LOGICAL                     ::                                   &
    lanalyt_calc_t0p0,& ! necessary only for irefatm=1
    lhhl_in_read,     & ! to indicate, whether hhl for input data has been read
    lhhl_hasbeenread    ! to indicate, whether hhl for input data has been read
!   lhhl_lm_read        ! to indicate, whether hhl for COSMO data has been read

! to indicate, whether a new vertical grid HHL has to be written
! LOGICAL                     ::                                   &
!   lnewVGrid           ! necessary only for GRIB2

! Values for SLEVE coordinate
  INTEGER    ::                                   &
    nfltvc              ! number of applications of filter for
!   nfltvc_in           ! decomposition of topography in large-scale
                        ! and small-scale part for SLEVE coordinate

  REAL (KIND=wp)              ::                                   &
    svc1,   svc2        ! decay rates for large-scale and small-scale
!   svc1_in,svc2_in     ! topography parts for SLEVE coordinate

! variables for default vertical coordinate parameters
! REAL (KIND=wp)              ::                                   &
!   vcoord_d(khmax)     !  to read the namelist variables

! for allocating vcoord_type and refatm_type
INTEGER     ::    &
  imax_vcoordtype,    & ! dimension of vcoord-type
  imax_refatmtype       ! dimension of refatm-type

REAL (KIND=wp)               ::    &
  rundefined = -9999.0_wp ! a default value that indicates that a
                        !    parameter has not been set

! data type for defaults of vertical coordinates
TYPE vcoord_type
  INTEGER                      :: ivctype            ! vertical coordinate type
  INTEGER                      :: ivcoord_id
  INTEGER                      :: nlevels            ! number of half levels (ke+1)
  INTEGER                      :: kflat              ! index where levels become flat
  CHARACTER(LEN=1)             :: vc_uuid(16)        ! UUID of HHL file
  REAL (KIND=wp),     POINTER  :: vcflat             ! coordinate where levels become flat
  REAL (KIND=wp),     POINTER, DIMENSION(:) ::       &
                                vert_coord => NULL() ! height above sea level
  REAL (KIND=wp),     POINTER, DIMENSION(:) ::       &
                                sigm_coord => NULL() ! reference pressure normalized to [0,1]
END TYPE vcoord_type

TYPE (vcoord_type), ALLOCATABLE  :: vcoord_defaults(:)
TYPE (vcoord_type), SAVE         :: vcoord_in, vcoord

! data type for defaults of reference atmosphere settings
TYPE refatm_type
  INTEGER                      :: irefatm
  INTEGER                      :: irefatm_id
  REAL (KIND=wp),     POINTER  :: p0sl
  REAL (KIND=wp),     POINTER  :: t0sl
  REAL (KIND=wp),     POINTER  :: dt0lp
  REAL (KIND=wp),     POINTER  :: delta_t
  REAL (KIND=wp),     POINTER  :: h_scal
  REAL (KIND=wp),     POINTER  :: bvref
END TYPE refatm_type

TYPE (refatm_type), ALLOCATABLE  :: refatm_defaults(:)
TYPE (refatm_type)               :: refatm_in, refatm

!------------------------------------------------------------------------------

! UUID part:
INTEGER, PARAMETER :: uuid_string_length = 36

!CHARACTER(LEN=1) :: uuid_in(16)
CHARACTER(LEN=1) :: uuid_out(16)
CHARACTER(LEN=uuid_string_length) :: uuid_in_string, uuid_out_string

!==============================================================================

CONTAINS

!==============================================================================
!+ Set defaults for the vertical coordinate parameters
!------------------------------------------------------------------------------

SUBROUTINE set_vcoord_defaults

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine sets some defaults for the specification of the vertical
!   coordinate parameters, commonly used in the COSMO-Model. It is important
!   that the INT2LM and the COSMO-Model use the same specifications. While
!   these values are transferred by GRIB1 or NetCDF meta data, they cannot be
!   transferred with GRIB2 meta data from INT2LM to the COSMO-Model and
!   between COSMO-Nudging and COSMO-forecast runs.
!
! Method:
!   Defaults are set (which were set in src_namelists before)
!
!------------------------------------------------------------------------------

INTEGER       :: izerror, i

! End of header
!------------------------------------------------------------------------------

izerror         = 0

imax_vcoordtype = 4

ALLOCATE (vcoord_defaults(imax_vcoordtype), STAT=izerror)

! and allocate the different type entries
DO i =1, imax_vcoordtype
   ALLOCATE(vcoord_defaults(i)%vert_coord(khmax))
   ALLOCATE(vcoord_defaults(i)%sigm_coord(khmax))
   ALLOCATE(vcoord_defaults(i)%vcflat)
END DO

! some initializations for vcoord_in, vcoord
vcoord_in%ivcoord_id = -1
vcoord_in%kflat      = -1

ALLOCATE(vcoord_in%vcflat)
vcoord_in%vcflat     = -1.0_wp
IF (.NOT. ASSOCIATED(vcoord_in%vert_coord)) THEN
  ALLOCATE(vcoord_in%vert_coord(khmax))
  vcoord_in%vert_coord(:) = -1.0_wp
ENDIF
IF (.NOT. ASSOCIATED(vcoord_in%sigm_coord)) THEN
  ALLOCATE(vcoord_in%sigm_coord(khmax))
  vcoord_in%sigm_coord(:) = -1.0_wp
ENDIF

vcoord%ivcoord_id    = -1
vcoord%kflat         = -1

ALLOCATE(vcoord%vcflat)
vcoord%vcflat        = -1.0_wp
IF (.NOT. ASSOCIATED(vcoord%vert_coord)) THEN
  ALLOCATE(vcoord%vert_coord(khmax))
   vcoord%vert_coord = -1.0_wp
ENDIF
IF (.NOT. ASSOCIATED(vcoord%sigm_coord)) THEN
  ALLOCATE(vcoord%sigm_coord(khmax))
   vcoord%sigm_coord = -1.0_wp
ENDIF

! Defaults for COSMO-EU:
i = 1
vcoord_defaults(i)%ivctype          =  2
vcoord_defaults(i)%ivcoord_id       =  i
vcoord_defaults(i)%nlevels          = 41
vcoord_defaults(i)%kflat            = -1        ! will be set later
vcoord_defaults(i)%vc_uuid(:)       = 'x'       ! have to get it first
vcoord_defaults(i)%vcflat           = 11430.0_wp
vcoord_defaults(i)%vert_coord(1:vcoord_defaults(i)%nlevels) = (/             &
          22700.0000_wp, 20800.0000_wp,                              &
          19100.0000_wp, 17550.0000_wp, 16150.0000_wp,           &
          14900.0000_wp, 13800.0000_wp, 12785.0000_wp,           &
          11875.0000_wp, 11020.0000_wp, 10205.0000_wp,           &
           9440.0000_wp,  8710.0000_wp,  8015.0000_wp,           &
           7355.0000_wp,  6725.0000_wp,  6130.0000_wp,           &
           5565.0000_wp,  5035.0000_wp,  4530.0000_wp,           &
           4060.0000_wp,  3615.0000_wp,  3200.0000_wp,           &
           2815.0000_wp,  2455.0000_wp,  2125.0000_wp,           &
           1820.0000_wp,  1545.0000_wp,  1295.0000_wp,           &
           1070.0000_wp,   870.0000_wp,   695.0000_wp,           &
            542.0000_wp,   412.0000_wp,   303.0000_wp,           &
            214.0000_wp,   143.0000_wp,    89.0000_wp,           &
             49.0000_wp,    20.0000_wp,     0.0000_wp /)
vcoord_defaults(i)%sigm_coord(1:vcoord_defaults(i)%nlevels) = -1.0_wp


! Defaults for COSMO-DE:
i = 2
vcoord_defaults(i)%ivctype          =  2
vcoord_defaults(i)%ivcoord_id       =  i
vcoord_defaults(i)%nlevels          = 51
vcoord_defaults(i)%kflat            = -1        ! will be set later
vcoord_defaults(i)%vc_uuid(:)       = 'x'       ! have to get it first
vcoord_defaults(i)%vcflat           = 11357.0_wp
vcoord_defaults(i)%vert_coord(1:vcoord_defaults(i)%nlevels) = (/             &
                 22000.00_wp, 21000.00_wp, 20028.57_wp,          &
                 19085.36_wp, 18170.00_wp, 17282.14_wp,          &
                 16421.43_wp, 15587.50_wp, 14780.00_wp,          &
                 13998.57_wp, 13242.86_wp, 12512.50_wp,          &
                 11807.14_wp, 11126.43_wp, 10470.00_wp,          &
                  9837.50_wp,  9228.57_wp,  8642.86_wp,          &
                  8080.00_wp,  7539.64_wp,  7021.43_wp,          &
                  6525.00_wp,  6050.00_wp,  5596.07_wp,          &
                  5162.86_wp,  4750.00_wp,  4357.14_wp,          &
                  3983.93_wp,  3630.00_wp,  3295.00_wp,          &
                  2978.57_wp,  2680.36_wp,  2400.00_wp,          &
                  2137.14_wp,  1891.43_wp,  1662.50_wp,          &
                  1450.00_wp,  1253.57_wp,  1072.86_wp,          &
                   907.50_wp,   757.14_wp,   621.43_wp,          &
                   500.00_wp,   392.50_wp,   298.57_wp,          &
                   217.86_wp,   150.00_wp,    94.64_wp,          &
                    51.43_wp,    20.00_wp,     0.00_wp /)
vcoord_defaults(i)%sigm_coord(1:vcoord_defaults(i)%nlevels) = -1.0_wp

! some defaults for the SLEVE coordinate
i = 3
vcoord_defaults(i)%ivctype          =  3
vcoord_defaults(i)%ivcoord_id       =  i
vcoord_defaults(i)%nlevels          = 51
vcoord_defaults(i)%kflat            = -1        ! will be set later
vcoord_defaults(i)%vc_uuid(:)       = 'x'       ! have to get it first
vcoord_defaults(i)%vcflat           = 11430.0_wp
vcoord_defaults(i)%vert_coord(1:vcoord_defaults(i)%nlevels) = (/             &
                 22000.00_wp, 21000.00_wp, 20028.57_wp,          &
                 19085.36_wp, 18170.00_wp, 17282.14_wp,          &
                 16421.43_wp, 15587.50_wp, 14780.00_wp,          &
                 13998.57_wp, 13242.86_wp, 12512.50_wp,          &
                 11807.14_wp, 11126.43_wp, 10470.00_wp,          &
                  9837.50_wp,  9228.57_wp,  8642.86_wp,          &
                  8080.00_wp,  7539.64_wp,  7021.43_wp,          &
                  6525.00_wp,  6050.00_wp,  5596.07_wp,          &
                  5162.86_wp,  4750.00_wp,  4357.14_wp,          &
                  3983.93_wp,  3630.00_wp,  3295.00_wp,          &
                  2978.57_wp,  2680.36_wp,  2400.00_wp,          &
                  2137.14_wp,  1891.43_wp,  1662.50_wp,          &
                  1450.00_wp,  1253.57_wp,  1072.86_wp,          &
                   907.50_wp,   757.14_wp,   621.43_wp,          &
                   500.00_wp,   392.50_wp,   298.57_wp,          &
                   217.86_wp,   150.00_wp,    94.64_wp,          &
                    51.43_wp,    20.00_wp,     0.00_wp /)
vcoord_defaults(i)%sigm_coord(1:vcoord_defaults(i)%nlevels) = -1.0_wp

! some defaults for the SLEVE2 coordinate
i = 4
vcoord_defaults(i)%ivctype          =  4
vcoord_defaults(i)%ivcoord_id       =  i
vcoord_defaults(i)%nlevels          = 51
vcoord_defaults(i)%kflat            = -1        ! will be set later
vcoord_defaults(i)%vc_uuid(:)       = 'x'       ! have to get it first
vcoord_defaults(i)%vcflat           = 11430.0_wp
vcoord_defaults(i)%vert_coord(1:vcoord_defaults(i)%nlevels) = (/             &
                 22000.00_wp, 21000.00_wp, 20028.57_wp,          &
                 19085.36_wp, 18170.00_wp, 17282.14_wp,          &
                 16421.43_wp, 15587.50_wp, 14780.00_wp,          &
                 13998.57_wp, 13242.86_wp, 12512.50_wp,          &
                 11807.14_wp, 11126.43_wp, 10470.00_wp,          &
                  9837.50_wp,  9228.57_wp,  8642.86_wp,          &
                  8080.00_wp,  7539.64_wp,  7021.43_wp,          &
                  6525.00_wp,  6050.00_wp,  5596.07_wp,          &
                  5162.86_wp,  4750.00_wp,  4357.14_wp,          &
                  3983.93_wp,  3630.00_wp,  3295.00_wp,          &
                  2978.57_wp,  2680.36_wp,  2400.00_wp,          &
                  2137.14_wp,  1891.43_wp,  1662.50_wp,          &
                  1450.00_wp,  1253.57_wp,  1072.86_wp,          &
                   907.50_wp,   757.14_wp,   621.43_wp,          &
                   500.00_wp,   392.50_wp,   298.57_wp,          &
                   217.86_wp,   150.00_wp,    94.64_wp,          &
                    51.43_wp,    20.00_wp,     0.00_wp /)
vcoord_defaults(i)%sigm_coord(1:vcoord_defaults(i)%nlevels) = -1.0_wp

END SUBROUTINE set_vcoord_defaults

!==============================================================================
!==============================================================================
!+ Set defaults for the reference atmosphere parameters
!------------------------------------------------------------------------------

SUBROUTINE set_refatm_defaults

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine sets some defaults for the specification of the reference
!   atmosphere, commonly used in the COSMO-Model. It is important
!   that the INT2LM and the COSMO-Model use the same specifications. While
!   these values are transferred by GRIB1 or NetCDF meta data, they cannot be
!   transferred with GRIB2 meta data from INT2LM to the COSMO-Model and
!   between COSMO-Nudging and COSMO-forecast runs.
!
! Method:
!   Defaults are set (which were set in src_namelists before)
!
!------------------------------------------------------------------------------

INTEGER       :: izerror, i

! End of header
!------------------------------------------------------------------------------

izerror         = 0

imax_refatmtype = 2

ALLOCATE (refatm_defaults(imax_refatmtype), STAT=izerror)

DO i=1,imax_refatmtype
  ALLOCATE(refatm_defaults(i)%p0sl)
  ALLOCATE(refatm_defaults(i)%t0sl)
  ALLOCATE(refatm_defaults(i)%dt0lp)
  ALLOCATE(refatm_defaults(i)%delta_t)
  ALLOCATE(refatm_defaults(i)%h_scal)
  ALLOCATE(refatm_defaults(i)%bvref)
ENDDO

ALLOCATE(refatm_in%p0sl)
ALLOCATE(refatm_in%t0sl)
ALLOCATE(refatm_in%dt0lp)
ALLOCATE(refatm_in%delta_t)
ALLOCATE(refatm_in%h_scal)
ALLOCATE(refatm_in%bvref)

ALLOCATE(refatm%p0sl)
ALLOCATE(refatm%t0sl)
ALLOCATE(refatm%dt0lp)
ALLOCATE(refatm%delta_t)
ALLOCATE(refatm%h_scal)
ALLOCATE(refatm%bvref)

! Defaults for irefatm=1 (still used in COSMO-DE)
i = 1
refatm_defaults(i)%irefatm           =  1
refatm_defaults(i)%irefatm_id        =  i
refatm_defaults(i)%p0sl              = 100000.0_wp
refatm_defaults(i)%t0sl              =    288.15_wp
refatm_defaults(i)%dt0lp             =     42.0_wp
refatm_defaults(i)%delta_t           =  rundefined
refatm_defaults(i)%h_scal            =  rundefined
refatm_defaults(i)%bvref             =  rundefined

! Defaults for irefatm=2 (used in COSMO-EU)
i = 2
refatm_defaults(i)%irefatm           =  2
refatm_defaults(i)%irefatm_id        =  i
refatm_defaults(i)%p0sl              = 100000.0_wp
refatm_defaults(i)%t0sl              =    288.15_wp
refatm_defaults(i)%dt0lp             =  rundefined
refatm_defaults(i)%delta_t           =     75.0_wp
refatm_defaults(i)%h_scal            =  10000.0_wp
refatm_defaults(i)%bvref             =  rundefined


! initialize the reference atmosphere parameters for incoming data with defaults
! in case of llm2lm, they are updated with input data later on
refatm_in%irefatm    = 1
refatm_in%irefatm_id = 0
refatm_in%p0sl       = 100000.0_wp
refatm_in%t0sl       =    288.15_wp
refatm_in%dt0lp      =     42.0_wp
refatm_in%delta_t    =  rundefined
refatm_in%h_scal     =  rundefined
refatm_in%bvref      =  rundefined

refatm%irefatm       = 1
refatm%irefatm_id    = 0
refatm%p0sl          = 100000.0_wp
refatm%t0sl          =    288.15_wp
refatm%dt0lp         =     42.0_wp
refatm%delta_t       =  rundefined
refatm%h_scal        =  rundefined
refatm%bvref         =  rundefined

END SUBROUTINE set_refatm_defaults

!==============================================================================
!==============================================================================

!------------------------------------------------------------------------------

SUBROUTINE reference_atmosphere                                               &
 ( hhl, p0, p0hl, rho0, t0, t0hl, dp0, hsurf, hsurfs, ie, je, ke,             &
   ra_type, vc_type, svc1, svc2, r, g, lanalyt_calc_t0p0, lnew_hhl,           &
   yerrmsg, ierror)

!------------------------------------------------------------------------------
!
! Description:
!   This routine computes the constant fields of the reference atmosphere
!   (p0, t0, rho0 and dp0) from the height hhl of the model half levels.
!   Also, the height (vert_coord) and normalized base-state pressure values
!   (sigm_coord) referring to mean sea level (i.e. z = 0.0 and p0sl = 1000 hPa)
!   of the model half levels are computed.
!
! Method:
!   a) First, the reference pressure at half levels is obtained by
!      analytical integration of the hydrostatic equation from hhl and
!      the parameters p0sl (base state pressure at sea level, i.e. z=0),
!      t0sl (base state pressure at sea level) and dt0lp (a constant
!      rate of decrease of temperature with the logarithm of pressure,
!      dT/dlnp = const). The height of the lower boundary, hhl(ke+1), is
!      conformal to the terrain height.
!      For detail see Section 3.1.2 of the Scientific Documentation.
!   b) The (base-state) pressure thickness dp0 of the model layers is
!      calculated, dp0(K) = p0(K+1/2) - p0(k-1/2).
!   c) The base-state pressure at main levels is defined as the arithmetic
!      average of the corresponding half-level pressures,
!      p0(k) =  0.5*( po(k+1/2) - p0(k-1/2) )
!   d) The base-state densitiy rho0(k) at main levels is defined as a
!      layer average using the hydrostatic equation.
!
!-------------------------------------------------------------------------------

! Parameter list:

INTEGER                 , INTENT (IN)    ::    &
  ie, je, ke          ! dimensions of the fields

TYPE(vcoord_type),        INTENT (INOUT) ::    &
  vc_type             ! to specify vertical coordinate parameters

TYPE(refatm_type),        INTENT (IN)    ::    &
  ra_type             ! to specify reference atmosphere parameters

REAL (KIND=wp),     INTENT (INOUT)          ::    &
  hhl(ie,je,ke+1)     ! geometrical height of the model half levels

REAL (KIND=wp),     INTENT (IN)          ::    &
  hsurf (ie,je),    & ! height of surface topography
  hsurfs(ie,je,2)     ! height of splitted topography parts (SLEVE coordinate)

REAL (KIND=wp),     INTENT (OUT)         ::    &
  p0    (ie,je,ke), & ! base-state pressure at full levels
  p0hl  (ie,je,ke+1), & ! base-state pressure at half levels
  t0hl  (ie,je,ke+1), & ! base-state temperature at half levels
  dp0   (ie,je,ke), & ! base-state pressure thickness of model layers
  rho0  (ie,je,ke), & ! reference density at full levels
  t0    (ie,je,ke)    ! reference temperature at full levels

REAL (KIND=wp),     INTENT (IN)          ::    &
  svc1,             & ! vertical decay rate for large-scale topo part
  svc2,             & ! vertical decay rate for small-scale topo part
  r   ,             & ! gas constant for dry air
  g                   ! gravity acceleration

LOGICAL,             INTENT(IN)          ::    &
  lanalyt_calc_t0p0,& ! to calculate t0 analytically
  lnew_hhl            ! if a new hhl has to be computed for GRIB2

CHARACTER (LEN=*), INTENT(OUT)           ::    &
  yerrmsg

INTEGER,           INTENT(OUT)           ::    &
  ierror

!-------------------------------------------------------------------------------

! Local variables:

INTEGER         :: k
REAL (KIND=wp)  :: zn, zgdrt, ztdbe, zbetf, zpxx, zp0sl, zt0sl, zdt0lp
REAL (KIND=wp)  :: zak(ke+1), zbk(ke+1), zbk2(ke+1), zhfl(ie,je)

!-------------------------------------------------------------------------------

! Begin subroutine reference_atmosphere
  ierror   = 0
  yerrmsg  = ' '

  zp0sl  = ra_type%p0sl
  zt0sl  = ra_type%t0sl
  zdt0lp = ra_type%dt0lp

  zgdrt = g/r/zt0sl
  IF (zdt0lp /= 0.0_wp) THEN
    ztdbe = zt0sl/zdt0lp
  ELSE
    ztdbe = 0.0_wp
  ENDIF
  zbetf = 2.0_wp*zdt0lp*zgdrt/zt0sl

!------------------------------------------------------------------------------
! Section 1: Compute reference state (pressure, hhl) on HALF levels
!------------------------------------------------------------------------------

  ! Select the vertical coordinate type
  SELECT CASE (vc_type%ivctype)
  CASE ( 1 )
    ! Pressure-based hybrid vertical coordinate on input
    ! here hhl depends on the reference atmosphere!
    ! ivctype=1 is deprecated for GRIB2

    ! Calculate the inverse coordinate transformation, i.e. the zak's and zbk's
    vc_type%kflat = 0
    DO k = 1, ke+1
      IF( vc_type%sigm_coord(k) <= vc_type%vcflat ) THEN
        zak(k) = vc_type%sigm_coord(k)*zp0sl
        zbk(k) = 0.0_wp
        vc_type%kflat = k
      ELSE
        zak(k) = vc_type%vcflat * zp0sl *                                       &
                     (1.0_wp - vc_type%sigm_coord(k))/(1.0_wp - vc_type%vcflat)
        zbk(k) = (vc_type%sigm_coord(k) - vc_type%vcflat)/(1.0_wp - vc_type%vcflat)
      ENDIF
    ENDDO

    ! Compute the surface reference pressure from surface topography
    IF (lnew_hhl) THEN
      hhl(:,:,ke+1) = hsurf(:,:)
    ENDIF
    IF (zdt0lp == 0.0_wp) THEN
      p0hl (:,:,ke+1) = zp0sl*EXP ( - zgdrt*hhl(:,:,ke+1) )
    ELSE
      p0hl (:,:,ke+1) = zp0sl*EXP ( - ztdbe &
                      * (1.0_wp - SQRT(1.0_wp - zbetf*hhl(:,:,ke+1))) )
    ENDIF

    ! Compute the reference pressure at half levels from surface topography
    ! and vertical coordinate parameters ak and bk as well as the
    ! height of half levels from the hydrostatic equation
    DO  k = 1, ke
      p0hl(:,:,k) = zak(k) + zbk(k)*p0hl(:,:,ke+1)
      IF (lnew_hhl) THEN
        hhl (:,:,k) = (r/g)*LOG(zp0sl/p0hl(:,:,k)) &
                            *( zt0sl - 0.5_wp*zdt0lp*LOG(zp0sl/p0hl(:,:,k)) )
      ENDIF
    ENDDO

    DO  k = 1, ke
      zpxx = zak(k) + zbk(k)*zp0sl
      vc_type%vert_coord(k) = (r/g) * LOG(zp0sl/zpxx) &
                                    * (zt0sl - 0.5_wp*zdt0lp*LOG(zp0sl/zpxx))
    ENDDO
    vc_type%vert_coord(ke+1) = 0.0_wp

  CASE ( 2, 3, 4 )
    ! Height-based hybrid vertical coordinate on input
    ! Vertical grid specified in terms of hhl
    ! here hhl depends only on the zak, zbk and vcflat

    IF     (vc_type%ivctype == 2) THEN
      ! "standard" coordinate with zak, zbk

      IF (lnew_hhl) THEN
        ! Calculate the inverse coordinate transformation, i.e. the zak's and zbk's
        vc_type%kflat = 0
        DO k = 1, ke+1
          IF( vc_type%vert_coord(k) >= vc_type%vcflat ) THEN
            zak(k) = vc_type%vert_coord(k)
            zbk(k) = 0.0_wp
            vc_type%kflat = k
          ELSE
            zak(k) = vc_type%vert_coord(k)
            zbk(k) = (vc_type%vcflat - vc_type%vert_coord(k))/ vc_type%vcflat
          ENDIF
        ENDDO

        ! Compute the height of the model half-levels
        hhl(:,:,ke+1) = hsurf(:,:)
        DO  k = 1, ke
          hhl(:,:,k) = zak(k) + zbk(k)*hhl(:,:,ke+1)
        ENDDO
!     ELSE
        ! HERE, PREVIOUSLY THE COMPUTATION OF VC_TYPE%KFLAT WAS DONE. THIS HAS
        !  BEEN REMOVED BECAUSE IT INVOLVES A GLOBAL COMMUNICATION AND HAS NOW
        !  TO BE DONE AFTER THE CALL TO REFERENCE_ATMOSPHERE_XXX() SUBROUTINES!
      ENDIF

    ELSEIF (vc_type%ivctype == 3 .OR. vc_type%ivctype == 4) THEN

      IF (lnew_hhl) THEN

        ! switch between original (zn = 1) and revised (zn = 1.35) SLEVE coordinate:
        zn = 1.0_wp
        IF (vc_type%ivctype == 4) THEN
          zn = 1.35_wp
        ENDIF

        ! Calculate the inverse coordinate transformation, i.e. the zak's and zbk's
        vc_type%kflat = 0
        DO k = 1, ke+1
          IF (vc_type%vert_coord(k) >= vc_type%vcflat ) THEN
            zak (k) = vc_type%vert_coord(k)
            zbk (k) = 0.0_wp
            zbk2(k) = 0.0_wp
            vc_type%kflat   = k
          ELSE
            zak (k) = vc_type%vert_coord(k)
            zbk (k) = SINH( (vc_type%vcflat/svc1)**zn - (vc_type%vert_coord(k)/svc1)**zn ) / &
                      SINH( (vc_type%vcflat/svc1)**zn )
            zbk2(k) = SINH( (vc_type%vcflat/svc2)**zn - (vc_type%vert_coord(k)/svc2)**zn ) / &
                      SINH( (vc_type%vcflat/svc2)**zn )
          ENDIF
        ENDDO

        ! Compute the height of the model half-levels
        hhl(:,:,ke+1) = hsurf(:,:)
        DO  k = 1, ke
          hhl(:,:,k) = zak(k) + zbk(k)*hsurfs(:,:,1) + zbk2(k)*hsurfs(:,:,2)
        ENDDO
!     ELSE
        ! HERE, PREVIOUSLY THE COMPUTATION OF VC_TYPE%KFLAT WAS DONE. THIS HAS
        !  BEEN REMOVED BECAUSE IT INVOLVES A GLOBAL COMMUNICATION AND HAS NOW
        !  TO BE DONE AFTER THE CALL TO REFERENCE_ATMOSPHERE_XXX() SUBROUTINES!
      ENDIF

    ENDIF

    ! Compute the reference pressure at half levels
    ! (is the same for ivctype = 1/2)
    DO  k = 1, ke+1
      IF (zdt0lp == 0.0_wp) THEN
        p0hl              (:,:,k) = zp0sl * EXP ( - zgdrt*hhl(:,:,k) )
        !IF (vc_type%vert_coord(k) >= 0.0_wp) THEN
        !  vc_type%sigm_coord(    k) = EXP ( - zgdrt*vc_type%vert_coord(k) )
        !ENDIF
      ELSE
        p0hl                (:,:,k) = zp0sl * EXP ( - ztdbe*(1.0_wp        &
                        - SQRT(1.0_wp - zbetf*hhl(:,:,k))) )
        !IF (vc_type%vert_coord(k) >= 0.0_wp) THEN
        !  vc_type%sigm_coord(    k) = EXP ( - ztdbe*(1.0_wp                &
        !                - SQRT(1.0_wp - zbetf*vc_type%vert_coord(k))) )
        !ENDIF
      ENDIF
    ENDDO

  END SELECT

  ! Basic consistency check for hhl
  IF ( ANY( hhl(:,:,1:ke) - hhl(:,:,2:ke+1) <= 0.0_wp) ) THEN
    WRITE(*,*) ' *** ERROR: negative dz encountered for reference_atmosphere!'
    ierror = 1
    yerrmsg = ' *** ERROR: negative dz encountered for reference_atmosphere!'
  ENDIF

!------------------------------------------------------------------------------
! Section 2: Compute reference temperature on HALF levels: t0hl
!------------------------------------------------------------------------------

  ! Compute the base-state temperature of the reference state at HALF levels
  DO  k = 1, ke+1
    t0hl (:,:,k) = zt0sl * SQRT( 1.0_wp - zbetf * hhl(:,:,k) )
  ENDDO

!------------------------------------------------------------------------------
! Section 3: Compute corresponding fields on FULL levels
!------------------------------------------------------------------------------

  ! Compute pressure, pressure thickness, temperature,
  ! and density of the reference state at FULL levels
  ! starting from the arithmetic mean of geopotential height

  IF (lanalyt_calc_t0p0) THEN

    ! this is only done, if the new fast waves solver is used in combination
    ! with the first reference atmosphere (irefatm = 1)

    ! PRINT*, "reference_atmosphere: analytic calculation of t0, p0, ..."

    DO  k = 1, ke
      zhfl(:,:)   = 0.5_wp * ( hhl(:,:,k+1) + hhl(:,:,k) )

      t0  (:,:,k) = zt0sl * SQRT( 1.0_wp - zbetf * zhfl(:,:) )

      IF (zdt0lp == 0.0_wp) THEN
        p0  (:,:,k) = zp0sl * EXP ( - zgdrt*zhfl(:,:) )
      ELSE
        p0  (:,:,k) = zp0sl * EXP ( - ztdbe*(1.0_wp        &
          - SQRT(1.0_wp - zbetf*zhfl(:,:))) )
      ENDIF

      dp0 (:,:,k) =       ( p0hl(:,:,k+1) - p0hl(:,:,k) )
      rho0(:,:,k) = p0(:,:,k) / ( r * t0(:,:,k) )
    ENDDO

  ELSE

    ! PRINT *, "reference_atmosphere: calculate t0, p0, ... by averaging"

    DO  k = 1, ke
      p0  (:,:,k) = 0.5_wp * ( p0hl(:,:,k+1) + p0hl(:,:,k) )
      dp0 (:,:,k) =       ( p0hl(:,:,k+1) - p0hl(:,:,k) )
      rho0(:,:,k) = dp0(:,:,k) / ( hhl(:,:,k) - hhl(:,:,k+1) ) / g
      t0  (:,:,k) = p0(:,:,k) / ( r * rho0(:,:,k) )
    ENDDO

  ENDIF

END SUBROUTINE reference_atmosphere

!==============================================================================
!==============================================================================

SUBROUTINE reference_atmosphere_2                                             &
 ( hhl, p0, p0hl, rho0, t0, t0hl, dp0, hsurf, hsurfs, ie, je, ke,             &
   ra_type, vc_type, svc1, svc2, r, g, lnew_hhl, yerrmsg, ierror)

!-------------------------------------------------------------------------------
!
! Description:
!   This routine computes the constant fields of the reference atmosphere
!   (p0, t0, rho0 and dp0) from the height hhl of the model half levels.
!   Also, the height (vert_coord) and normalized base-state pressure values
!   (sigm_coord) referring to mean sea level (i.e. z = 0.0 and p0sl = 1000 hPa)
!   of the model half levels are computed.
!
!   The new reference atmosphere is based on the temperature profile
!   t0(z) = (t0sl-delta_t) + delta_t*exp(-z/h_scal),
!   where z = hhl(k) is the height of a model grid point.
!
! Method:
!   a) Depending on the vertical coordinate type, the procedure starts with
!      computing the height or pressure at half levels by analytical
!      integration of the hydrostatic equation, starting from p0sl (base
!      state pressure at sea level, i.e. z=0).
!   b) The (base-state) pressure thickness dp0 of the model layers is
!      calculated, dp0(K) = p0(K+1/2) - p0(k-1/2).
!   c) The base-state pressure at main levels is computed from the analytically
!      integrated hydrostatic equation, assuming that the height at full levels
!      (zhfl) is the arithmetic mean of the adjacent half levels.
!      The latter assumption is also made in the RK dynamical core.
!   d) The base-state density rho0(k) is computed from base-state pressure
!      and base-state temperature using the ideal gas equation.
!
!-------------------------------------------------------------------------------

! Parameter list:

INTEGER,           INTENT (IN)         ::    &
  ie, je, ke           ! dimensions of the fields

TYPE(vcoord_type), INTENT (INOUT)      ::    &
  vc_type              ! to specify vertical coordinate parameters

TYPE(refatm_type), INTENT (IN)         ::    &
  ra_type              ! to specify reference atmosphere parameters

REAL (KIND=wp),    INTENT (INOUT)      ::    &
  hhl(ie,je,ke+1)      ! geometrical height of the model half levels

REAL (KIND=wp),    INTENT (IN)         ::    &
  hsurf (ie,je),     & ! height of surface topography
  hsurfs(ie,je,2)      ! height of splitted topography parts (SLEVE coordinate)

REAL (KIND=wp),    INTENT (OUT)        ::    &
  p0    (ie,je,ke),  & ! base-state pressure at full levels
  p0hl  (ie,je,ke+1),& ! base-state pressure at half levels
  t0hl  (ie,je,ke+1),& ! base-state temperature at half levels
  dp0   (ie,je,ke),  & ! base-state pressure thickness of model layers
  rho0  (ie,je,ke),  & ! reference density at full levels
  t0    (ie,je,ke)     ! reference temperature at full levels

REAL (KIND=wp),    INTENT (IN)         ::    &
  svc1,              & ! vertical decay rate for large-scale topo part
  svc2,              & ! vertical decay rate for small-scale topo part
  r   ,              & ! gas constant for dry air
  g                    ! gravity acceleration

LOGICAL,           INTENT(IN)          ::    &
  lnew_hhl             ! if a new hhl has to be computed for GRIB2

CHARACTER (LEN=*), INTENT(OUT)         ::    &
  yerrmsg

INTEGER,           INTENT(OUT)         ::    &
  ierror

!-------------------------------------------------------------------------------

! Local variables:

INTEGER        :: k
REAL (KIND=wp) :: zpxx, zn, zt00, zp0sl, zt0sl, zdelta_t, zh_scal
REAL (KIND=wp) :: zak(ke+1), zbk(ke+1), zbk2(ke+1), zhfl(ie,je)

!-------------------------------------------------------------------------------

! Begin subroutine reference_atmosphere_2
  ierror   = 0
  yerrmsg  = ' '

  zp0sl    = ra_type%p0sl
  zt0sl    = ra_type%t0sl
  zdelta_t = ra_type%delta_t
  zh_scal  = ra_type%h_scal
  zt00     = zt0sl - zdelta_t

!------------------------------------------------------------------------------
! Section 1: Compute reference state (pressure, hhl) on HALF levels
!------------------------------------------------------------------------------

  ! Select the vertical coordinate type
  SELECT CASE (vc_type%ivctype)
  CASE ( 1 )
    ! Pressure-based hybrid vertical coordinate on input
    ! here hhl depends on the reference atmosphere!
    ! ivctype=1 is deprecated for GRIB2

    ! Calculate the inverse coordinate transformation, i.e. the zak's and zbk's
    vc_type%kflat = 0
    DO k = 1, ke+1
      IF( vc_type%sigm_coord(k) <= vc_type%vcflat ) THEN
        zak(k) = vc_type%sigm_coord(k)*zp0sl
        zbk(k) = 0.0_wp
        vc_type%kflat = k
      ELSE
        zak(k) =  vc_type%vcflat*zp0sl*(1.0_wp - vc_type%sigm_coord(k)) /     &
                                       (1.0_wp - vc_type%vcflat)
        zbk(k) = (vc_type%sigm_coord(k) - vc_type%vcflat)/(1.0_wp - vc_type%vcflat)
      ENDIF
    ENDDO

    ! Compute the surface reference pressure from surface topography
    IF (lnew_hhl) THEN
      hhl  (:,:,ke+1) = hsurf(:,:)
    ENDIF
    p0hl (:,:,ke+1) = zp0sl*EXP ( - g/r*zh_scal/zt00 * LOG( &
                       (EXP(hhl(:,:,ke+1)/zh_scal)*zt00 + zdelta_t)/(zt00 + zdelta_t)) )

    ! Compute the reference pressure at half levels from surface topography
    ! and vertical coordinate parameters zak and zbk as well as the
    ! height of half levels from the hydrostatic equation
    DO  k = 1, ke
      p0hl(:,:,k) = zak(k) + zbk(k)*p0hl(:,:,ke+1)

      IF (lnew_hhl) THEN
        hhl(:,:,k)  = zh_scal*LOG( (EXP( -r*zt00/(g*zh_scal) * &
                     LOG(p0hl(:,:,k)/zp0sl))*(zt00 + zdelta_t)-zdelta_t)/zt00)
      ENDIF
    ENDDO

    DO  k = 1, ke
      zpxx = zak(k) + zbk(k)*zp0sl
      vc_type%vert_coord(k) = zh_scal * LOG ( (EXP(-r*zt00/(g*zh_scal) * &
                                 LOG(zpxx/zp0sl))*(zt00 + zdelta_t)-zdelta_t)/zt00)
    ENDDO
    vc_type%vert_coord(ke+1) = 0.0_wp

  CASE ( 2, 3, 4 )
    ! Height-based hybrid vertical coordinate on input
    ! Vertical grid specified in terms of hhl
    ! here hhl depends only on the zak, zbk and vcflat

    IF     (vc_type%ivctype == 2) THEN
      ! "standard" coordinate with zak, zbk

      IF (lnew_hhl) THEN
        ! Calculate the inverse coordinate transformation, i.e. the zak's and zbk's
        vc_type%kflat = 0
        DO k = 1, ke+1
          IF( vc_type%vert_coord(k) >= vc_type%vcflat ) THEN
            zak(k) = vc_type%vert_coord(k)
            zbk(k) = 0.0_wp
            vc_type%kflat = k
          ELSE
            zak(k) = vc_type%vert_coord(k)
            zbk(k) = (vc_type%vcflat - vc_type%vert_coord(k))/ vc_type%vcflat
          ENDIF
        ENDDO

        ! Compute the height of the model half-levels
        hhl(:,:,ke+1) = hsurf(:,:)
        DO  k = 1, ke
          hhl(:,:,k) = zak(k) + zbk(k)*hhl(:,:,ke+1)
        ENDDO
!     ELSE
        ! HERE, PREVIOUSLY THE COMPUTATION OF VC_TYPE%KFLAT WAS DONE. THIS HAS
        !  BEEN REMOVED BECAUSE IT INVOLVES A GLOBAL COMMUNICATION AND HAS NOW
        !  TO BE DONE AFTER THE CALL TO REFERENCE_ATMOSPHERE_XXX() SUBROUTINES!
      ENDIF

    ELSEIF (vc_type%ivctype == 3 .OR. vc_type%ivctype == 4) THEN

      IF (lnew_hhl) THEN
        ! switch between original (zn = 1) and revised (zn = 1.35) SLEVE coordinate:
        zn = 1.0_wp
        IF (vc_type%ivctype == 4) THEN
          zn = 1.35_wp
        ENDIF

        ! Calculate the inverse coordinate transformation, i.e. the zak's and zbk's
        vc_type%kflat = 0
        DO k = 1, ke+1
          IF( vc_type%vert_coord(k) >= vc_type%vcflat ) THEN
            zak (k) = vc_type%vert_coord(k)
            zbk (k) = 0.0_wp
            zbk2(k) = 0.0_wp
            vc_type%kflat = k
          ELSE
            zak (k) = vc_type%vert_coord(k)
            zbk (k) = SINH( (vc_type%vcflat/svc1)**zn - (vc_type%vert_coord(k)/svc1)**zn ) / &
                      SINH( (vc_type%vcflat/svc1)**zn )
            zbk2(k) = SINH( (vc_type%vcflat/svc2)**zn - (vc_type%vert_coord(k)/svc2)**zn ) / &
                      SINH( (vc_type%vcflat/svc2)**zn )
          ENDIF
        ENDDO

        ! Compute the height of the model half-levels
        hhl(:,:,ke+1) = hsurf(:,:)
        DO  k = 1, ke
          hhl(:,:,k) = zak(k) + zbk(k)*hsurfs(:,:,1) + zbk2(k)*hsurfs(:,:,2)
        ENDDO
!     ELSE
        ! HERE, PREVIOUSLY THE COMPUTATION OF VC_TYPE%KFLAT WAS DONE. THIS HAS
        !  BEEN REMOVED BECAUSE IT INVOLVES A GLOBAL COMMUNICATION AND HAS NOW
        !  TO BE DONE AFTER THE CALL TO REFERENCE_ATMOSPHERE_XXX() SUBROUTINES!
      ENDIF

    ENDIF

    ! Compute the reference pressure at half levels
    ! (is the same for all ivctype)
    DO  k = 1, ke+1
      p0hl                (:,:,k) = zp0sl*EXP ( - g/r*zh_scal/zt00 *         &
                 LOG( (EXP(hhl(:,:,k)/zh_scal)*zt00 + zdelta_t)/(zt00 + zdelta_t)) )
      !IF (vc_type%vert_coord(k) >= 0.0_wp) THEN
      !  vc_type%sigm_coord(    k) =       EXP ( - g/r*zh_scal/zt00 *         &
      !           LOG( (EXP(vc_type%vert_coord(k)/zh_scal)*zt00 + zdelta_t)/(zt00 + zdelta_t)) )
      !ENDIF
    ENDDO

  END SELECT

  ! Basic consistency check for hhl
  IF ( ANY( hhl(:,:,1:ke) - hhl(:,:,2:ke+1) <= 0.0_wp) ) THEN
    WRITE(*,*) ' *** ERROR: negative dz encountered for reference_atmosphere_2!'
    ierror = 2
    yerrmsg = ' *** ERROR: negative dz encountered for reference_atmosphere_2!'
  ENDIF

!------------------------------------------------------------------------------
! Section 2: Compute reference temperature on HALF levels: t0hl
!------------------------------------------------------------------------------

  ! Compute the base-state temperature of the reference state at HALF levels
  DO  k = 1, ke+1
    t0hl  (:,:,k) = zt00 + zdelta_t * EXP( -hhl(:,:,k) / zh_scal )
  ENDDO

!------------------------------------------------------------------------------
! Section 3: Compute corresponding fields on FULL levels
!------------------------------------------------------------------------------

  ! Compute pressure, pressure thickness, temperature,
  ! and density of the reference state at FULL levels
  ! starting from the arithmetic mean of geopotential height

  DO  k = 1, ke
    dp0 (:,:,k) =       ( p0hl(:,:,k+1) - p0hl(:,:,k) )
    zhfl(:,:)   = 0.5_wp * ( hhl(:,:,k+1) + hhl(:,:,k) )
    p0  (:,:,k) = zp0sl*EXP ( - g/r*zh_scal/zt00 * LOG( &
                  (EXP(zhfl(:,:)/zh_scal)*zt00 + zdelta_t)/(zt00 + zdelta_t)) )
    t0  (:,:,k) = zt00 + zdelta_t * EXP( -zhfl(:,:) / zh_scal )
    rho0(:,:,k) = p0  (:,:,k)/(r*t0(:,:,k))
  ENDDO

END SUBROUTINE reference_atmosphere_2

!==============================================================================
!==============================================================================

SUBROUTINE reference_atmosphere_BVconst                                       &
 ( hhl, p0, p0hl, rho0, t0, t0hl, dp0, hsurf, hsurfs, ie, je, ke,             &
   ra_type, vc_type, svc1, svc2, r, g, lnew_hhl, yerrmsg, ierror)

!-------------------------------------------------------------------------------
!
! Description:
!   This routine computes the constant fields of the reference atmosphere
!   (p0, t0, rho0 and dp0) from the height hhl of the model half levels.
!   Also, the height (vert_coord) and normalized base-state pressure values
!   (sigm_coord) referring to mean sea level (i.e. z = 0.0 and p0sl = 1000 hPa)
!   of the model half levels are computed.
!
!   The reference atmosphere is based on assuming a constant Brunt-Vaisala
!   frequency, yielding the following reference temperature profile
!   t0(z) = thsl*EXP(N**2/g*z)*(p/p00)**(R/c_p)
!   where z = hhl(k) is the height of a model grid point.
!
! Method:
!   a) Depending on the vertical coordinate type, the procedure starts with
!      computing the height or pressure at half levels by analytical
!      integration of the hydrostatic equation, starting from p0sl (base
!      state pressure at sea level, i.e. z=0).
!   b) The (base-state) pressure thickness dp0 of the model layers is
!      calculated, dp0(K) = p0(K+1/2) - p0(k-1/2).
!   c) The base-state pressure at main levels is computed from the analytically
!      integrated hydrostatic equation, assuming that the height at full levels
!      (zhfl) is the arithmetic mean of the adjacent half levels.
!      The latter assumption is also made in the RK dynamical core.
!   d) The base-state density rho0(k) is computed from base-state pressure
!      and base-state temperature using the ideal gas equation.
!
!-------------------------------------------------------------------------------

! Parameter list:

INTEGER,             INTENT (IN)       ::    &
  ie, je, ke          ! dimensions of the fields

TYPE(vcoord_type), INTENT (INOUT)      ::    &
  vc_type              ! to specify vertical coordinate parameters

TYPE(refatm_type), INTENT (IN)         ::    &
  ra_type              ! to specify reference atmosphere parameters

REAL (KIND=wp),    INTENT (INOUT)      ::    &
  hhl(ie,je,ke+1)     ! geometrical height of the model half levels

REAL (KIND=wp),    INTENT (IN)         ::    &
  hsurf(ie,je),     & ! height of surface topography
  hsurfs(ie,je,2)     ! height of splitted topography parts (SLEVE coordinate)

REAL (KIND=wp),    INTENT (OUT)        ::    &
  p0    (ie,je,ke), & ! base-state pressure at full levels
  p0hl  (ie,je,ke+1), & ! base-state pressure at half levels
  t0hl  (ie,je,ke+1), & ! base-state temperature at half levels
  dp0   (ie,je,ke), & ! base-state pressure thickness of model layers
  rho0  (ie,je,ke), & ! reference density at full levels
  t0    (ie,je,ke)    ! reference temperature at full levels

REAL (KIND=wp),    INTENT (IN)         ::    &
  svc1,             & ! vertical decay rate for large-scale topo part
  svc2,             & ! vertical decay rate for small-scale topo part
  r   ,             & ! gas constant for dry air
  g                   ! gravity acceleration

LOGICAL,           INTENT(IN)          ::    &
  lnew_hhl             ! if a new hhl has to be computed for GRIB2

CHARACTER (LEN=*), INTENT(OUT)         ::    &
  yerrmsg

INTEGER,           INTENT(OUT)         ::    &
  ierror

!-------------------------------------------------------------------------------

! Local variables:

INTEGER         :: k
REAL (KIND=wp)  :: zn, zpxx, zp00, zc_p,   &
                   zp0sl, zt0sl, zthsl, zbvref, zg2_o_N2cp
REAL (KIND=wp)  :: zak(ke+1), zbk(ke+1), zbk2(ke+1), zhfl(ie,je)

!-------------------------------------------------------------------------------

! Begin subroutine reference_atmosphere_BVconst
  ierror   = 0
  yerrmsg  = ' '

  zc_p   = 3.5_wp * r ! specific heat capacity of air at constant pressure
  zp00   = 1.E5_wp    ! reference pressure (Pa) for potential temperature

  zp0sl  = ra_type%p0sl
  zt0sl  = ra_type%t0sl
  zbvref = ra_type%bvref
  zthsl  = ra_type%t0sl*(zp00/ra_type%p0sl)**(r/zc_p)

  zg2_o_N2cp = g**2 / ( ra_type%bvref**2 * zc_p )

!------------------------------------------------------------------------------
! Section 1: Compute reference state (pressure, hhl) on HALF levels
!------------------------------------------------------------------------------

  ! Select the vertical coordinate type
  SELECT CASE (vc_type%ivctype)
  CASE ( 1 )
    ! Pressure-based hybrid vertical coordinate on input
    ! here hhl depends on the reference atmosphere!

    ! Calculate the inverse coordinate transformation, i.e. the ak's and bk's
    vc_type%kflat = 0
    DO k = 1, ke+1
      IF( vc_type%sigm_coord(k) <= vc_type%vcflat ) THEN
        zak(k) = vc_type%sigm_coord(k)*zp0sl
        zbk(k) = 0.0_wp
        vc_type%kflat = k
      ELSE
        zak(k) = vc_type%vcflat*zp0sl*(1.0_wp - vc_type%sigm_coord(k)) /      &
                                      (1.0_wp - vc_type%vcflat)
        zbk(k) = (vc_type%sigm_coord(k) - vc_type%vcflat)/(1.0_wp - vc_type%vcflat)
      ENDIF
    ENDDO

    ! Compute the surface reference pressure from surface topography
    IF (lnew_hhl) THEN
      hhl  (:,:,ke+1) = hsurf(:,:)
    ENDIF
    p0hl (:,:,ke+1) = zp0sl*(1.0_wp + zg2_o_N2cp/zt0sl *                &
                      (EXP(-zbvref**2/g*hhl(:,:,ke+1))-1.0_wp))**(zc_p/r)

    ! Compute the reference pressure at half levels from surface topography
    ! and vertical coordinate parameters ak and bk as well as the
    ! height of half levels from the hydrostatic equation
    DO  k = 1, ke
      p0hl(:,:,k) = zak(k) + zbk(k)*p0hl(:,:,ke+1)

      IF (lnew_hhl) THEN
        hhl(:,:,k)  = -g/zbvref**2*LOG(1.0_wp+((p0hl(:,:,k)/zp00)**(r/zc_p)-&
                      (zp0sl/zp00)**(r/zc_p)) * zthsl/zg2_o_N2cp)
      ENDIF
    ENDDO

    DO  k = 1, ke
      zpxx = zak(k) + zbk(k)*zp0sl
      vc_type%vert_coord(k) = -g/zbvref**2*LOG(1.0_wp+((zpxx/zp00)**(r/zc_p)-&
                               (zp0sl/zp00)**(r/zc_p))* zthsl/zg2_o_N2cp)
    ENDDO
    vc_type%vert_coord(ke+1) = 0.0_wp

  CASE ( 2, 3, 4 )
    ! Height-based hybrid vertical coordinate on input
    ! Vertical grid specified in terms of hhl
    ! here hhl depends only on the zak, zbk and vcflat

    IF     (vc_type%ivctype == 2) THEN
      ! "standard" coordinate with zak, zbk

      IF (lnew_hhl) THEN
        ! Calculate the inverse coordinate transformation, i.e. the ak's and bk's
        vc_type%kflat = 0
        DO k = 1, ke+1
          IF( vc_type%vert_coord(k) >= vc_type%vcflat ) THEN
            zak(k) = vc_type%vert_coord(k)
            zbk(k) = 0.0_wp
            vc_type%kflat = k
          ELSE
            zak(k) = vc_type%vert_coord(k)
            zbk(k) = (vc_type%vcflat - vc_type%vert_coord(k))/ vc_type%vcflat
          ENDIF
        ENDDO

        ! Compute the height of the model half-levels
        hhl(:,:,ke+1) = hsurf(:,:)
        DO  k = 1, ke
          hhl(:,:,k) = zak(k) + zbk(k)*hhl(:,:,ke+1)
        ENDDO
!     ELSE
        ! HERE, PREVIOUSLY THE COMPUTATION OF VC_TYPE%KFLAT WAS DONE. THIS HAS
        !  BEEN REMOVED BECAUSE IT INVOLVES A GLOBAL COMMUNICATION AND HAS NOW
        !  TO BE DONE AFTER THE CALL TO REFERENCE_ATMOSPHERE_XXX() SUBROUTINES!
      ENDIF

    ELSEIF (vc_type%ivctype == 3 .OR. vc_type%ivctype == 4) THEN

      IF (lnew_hhl) THEN
        ! switch between original (zn = 1) and revised (zn = 1.35) SLEVE coordinate:
        zn = 1.0_wp
        IF (vc_type%ivctype == 4) THEN
          zn = 1.35_wp
        ENDIF

        ! Calculate the inverse coordinate transformation, i.e. the ak's and bk's
        vc_type%kflat = 0
        DO k = 1, ke+1
          IF (vc_type%vert_coord(k) >= vc_type%vcflat) THEN
            zak (k) = vc_type%vert_coord(k)
            zbk (k) = 0.0_wp
            zbk2(k) = 0.0_wp
            vc_type%kflat = k
          ELSE
            zak (k) = vc_type%vert_coord(k)
            zbk (k) = SINH( (vc_type%vcflat/svc1)**zn - (vc_type%vert_coord(k)/svc1)**zn ) / &
                      SINH( (vc_type%vcflat/svc1)**zn )
            zbk2(k) = SINH( (vc_type%vcflat/svc2)**zn - (vc_type%vert_coord(k)/svc2)**zn ) / &
                      SINH( (vc_type%vcflat/svc2)**zn )
          ENDIF
        ENDDO

        ! Compute the height of the model half-levels
        hhl(:,:,ke+1) = hsurf(:,:)
        DO  k = 1, ke
           hhl(:,:,k) = zak(k) + zbk(k)*hsurfs(:,:,1) + zbk2(k)*hsurfs(:,:,2)
        ENDDO
!     ELSE
        ! HERE, PREVIOUSLY THE COMPUTATION OF VC_TYPE%KFLAT WAS DONE. THIS HAS
        !  BEEN REMOVED BECAUSE IT INVOLVES A GLOBAL COMMUNICATION AND HAS NOW
        !  TO BE DONE AFTER THE CALL TO REFERENCE_ATMOSPHERE_XXX() SUBROUTINES!
      ENDIF

    ENDIF

    ! Compute the reference pressure at half levels
    ! (is the same for ivctype = 1/2)
    DO  k = 1, ke+1
      p0hl  (:,:,k)    = zp0sl * (1.0_wp + zg2_o_N2cp/zt0sl *                  &
                        (EXP(-zbvref**2/g*hhl(:,:,k))-1.0_wp))**(zc_p/r)

      IF (vc_type%vert_coord(k) >= 0.0_wp) THEN
        vc_type%sigm_coord (k) = (1.0_wp + zg2_o_N2cp/zt0sl*                   &
                        (EXP(-zbvref**2/g*vc_type%vert_coord(k))-1.0_wp))**(zc_p/r)
      ENDIF
    ENDDO

  END SELECT

  ! Basic consistency check for hhl
  IF ( ANY( hhl(:,:,1:ke) - hhl(:,:,2:ke+1) <= 0.0_wp) ) THEN
    WRITE(*,*) ' *** ERROR: negative dz encountered for reference_atmosphere_2!'
    ierror = 3
    yerrmsg = ' *** ERROR: negative dz encountered for reference_atmosphere_2!'
  ENDIF

!------------------------------------------------------------------------------
! Section 2: Compute reference temperature on HALF levels: t0hl
!------------------------------------------------------------------------------

  DO  k = 1, ke+1
    t0hl  (:,:,k) = zg2_o_N2cp + ( zt0sl - zg2_o_N2cp ) * EXP( zbvref**2 / g * hhl(:,:,k) )
  END DO

!------------------------------------------------------------------------------
! Section 3: Compute corresponding fields on FULL levels
!------------------------------------------------------------------------------

  ! Compute the (base-state) pressure thickness, the base-state density
  ! and the base-state pressure at FULL levels,
  DO  k = 1, ke
    dp0 (:,:,k) =       ( p0hl(:,:,k+1) - p0hl(:,:,k) )
    zhfl(:,:)   = 0.5_wp * ( hhl(:,:,k+1) + hhl(:,:,k) )
    p0  (:,:,k) = zp0sl * ( 1.0_wp + zg2_o_N2cp/zt0sl *                &
                  (EXP(-zbvref**2/g*zhfl(:,:))-1.0_wp) )**(zc_p/r)
    t0  (:,:,k) = zg2_o_N2cp                                               &
                + ( zt0sl - zg2_o_N2cp ) * EXP( zbvref**2 / g * zhfl(:,:) )
    rho0(:,:,k) = p0  (:,:,k)/(r*t0(:,:,k))
  ENDDO

END SUBROUTINE reference_atmosphere_BVconst

!==============================================================================
!==============================================================================

SUBROUTINE k_index_of_pressure_levels(  p0sl, sigmr, ke, llm,      &
  klv950, klv850, klv800, klv700, klv500, klv400, klv300 )

!-------------------------------------------------------------------------------
!
! Description:
! Determination of the k-indices for those full levels which correspond
! to about 850, 800, 500, 400 and 300 hPa.
!
!-------------------------------------------------------------------------------

  INTEGER,         INTENT (IN)              ::    &
    ke

  REAL (KIND=wp),  INTENT (IN)              ::    &
    p0sl,             & ! reference pressure at sea level
    sigmr (ke+1)        ! sigma-coordinate refering to PMSL

  LOGICAL,         INTENT (IN),  OPTIONAL   ::    &
    llm                 ! if .TRUE., running with lowered upper boundary


  INTEGER,         INTENT (OUT), OPTIONAL   ::    &
    klv950, klv850, klv800, klv700, klv500, klv400, klv300
  ! k-indices of the LM full levels corresponding to about xxx hPa


  INTEGER        :: k
  REAL (KIND=wp) :: zpnu, zpno


  DO k = 1,ke
    zpno = p0sl*sigmr(k  )
    zpnu = p0sl*sigmr(k+1)
    IF (PRESENT(klv950)) THEN
      IF ( (zpno <= 950.0E2_wp) .AND. (950.0E2_wp < zpnu) )    klv950 = k
    ENDIF
    IF (PRESENT(klv850)) THEN
      IF ( (zpno <= 850.0E2_wp) .AND. (850.0E2_wp < zpnu) )    klv850 = k
    ENDIF
    IF (PRESENT(klv800)) THEN
      IF ( (zpno <= 800.0E2_wp) .AND. (800.0E2_wp < zpnu) )    klv800 = k
    ENDIF
    IF (PRESENT(klv700)) THEN
      IF ( (zpno <= 700.0E2_wp) .AND. (700.0E2_wp < zpnu) )    klv700 = k
    ENDIF
    IF (PRESENT(klv500)) THEN
      IF ( (zpno <= 500.0E2_wp) .AND. (500.0E2_wp < zpnu) )    klv500 = k
    ENDIF
    IF (PRESENT(klv400)) THEN
      IF ( (zpno <= 400.0E2_wp) .AND. (400.0E2_wp < zpnu) )    klv400 = k
    ENDIF
    IF (PRESENT(klv300)) THEN
      IF ( (zpno <= 300.0E2_wp) .AND. (300.0E2_wp < zpnu) )    klv300 = k
    ENDIF
  ENDDO

  IF (PRESENT(llm)) THEN
    IF (llm) THEN
      ! set the k-indices for the upper hPa levels to 0
      klv300 = 0
      klv400 = 0
      klv500 = 0
    ENDIF
  ENDIF

END SUBROUTINE k_index_of_pressure_levels

!==============================================================================

! and the UUID Part

!==============================================================================
!==============================================================================

SUBROUTINE uuid_create (newuuid)

!------------------------------------------------------------------------------
!
! Description:
!  Creates a new pseudo-random universally unique identifier. The procedure is
!  based on the Fortran standard routine random_number, which is initialized
!  by random_seed. The values for random_seed are taken from the actual date
!  and time and are modified a bit.
!  random_number returns 16 positive real numbers in the range [0,1). These
!  numbers are multiplied with 256 and the INT-part is taken. The results are
!  16 integers in the range [0,256). These numbers are transferred to
!  characters, which then build the UUID.
!
!------------------------------------------------------------------------------

CHARACTER(LEN=1), INTENT(OUT) :: newuuid(16)

!------------------------------------------------------------------------------

! local variables
REAL             :: rtmp(16)  ! should be default real
INTEGER          :: itmp(16), i, ksize, idati(8), it1(2), it2(2), is1, is2, istat
INTEGER(i8)      :: itms1, itms2

INTEGER, ALLOCATABLE :: iseed(:)

!------------------------------------------------------------------------------

istat = 0

CALL random_seed
CALL random_seed (size=ksize)

ALLOCATE (iseed(ksize), STAT=istat)

CALL DATE_AND_TIME(values=idati)
itms1 = (idati(1) - 1970) * 365_i8 * 24 * 60 * 60 * 1000 &
       + idati(2) * 31_i8 * 24 * 60 * 60 * 1000 &
       + idati(3) * 24 * 60 * 60 * 60 * 1000 &
       + idati(5) * 60 * 60 * 1000 &
       + idati(6) * 60 * 1000      &
       + idati(7) * 1000 &
       + idati(8)
itms2 = ISHFTC (itms1, 16)

it1 = TRANSFER(itms1, it1)
is1 = IEOR (it1(1), it1(2))
it2 = TRANSFER(itms2, it2)
is2 = IEOR (it2(1), it2(2))

IF (ksize >= 6) THEN
  iseed(1) = it1(1) +  37211
  iseed(2) = it2(1) +  74421
  iseed(3) = is1    + 111631
  iseed(4) = it1(2) +  41113
  iseed(5) = it2(2) +  82225
  iseed(6) = is2    + 123341
  IF (ksize > 6) THEN
    DO i = 7, ksize
      iseed(i) = iseed(i-6) + 1313 * (i-6)
    ENDDO
  ENDIF
ELSE
  iseed(1) = it1(1)
  IF (ksize > 1) iseed(2) = it2(1)
  IF (ksize > 2) iseed(3) = is1
  IF (ksize > 3) iseed(4) = it2(2)
  IF (ksize > 4) iseed(5) = it1(2)
ENDIF

CALL random_seed (put=iseed)
CALL random_number(rtmp)

itmp(:)    = INT (rtmp(:) *256 )
newuuid(:) = ACHAR(itmp(:))

DEALLOCATE (iseed)

!------------------------------------------------------------------------------

END SUBROUTINE uuid_create

!==============================================================================
!==============================================================================

SUBROUTINE uuid_2char (uuid, uuid_string)

!------------------------------------------------------------------------------
!
! Description:
!  The character-array uuid is converted to a human-readable string.
!  The uuid is made up of bytes, and the contents of these bytes need not
!  be printable characters. But it could be converted to hexadecimal numbers,
!  which could be printed.
!
!------------------------------------------------------------------------------

  CHARACTER(LEN=1),                  INTENT(IN)  :: uuid(16)
  CHARACTER(LEN=uuid_string_length), INTENT(OUT) :: uuid_string

!------------------------------------------------------------------------------

! local variables
  CHARACTER(LEN=32)  :: tmp_string
  INTEGER            :: i

!------------------------------------------------------------------------------

! convert content of bytes to hexadecimal value
DO i = 1, 16
  tmp_string(2*i-1:2*i) = byte2hex(uuid(i))
ENDDO

! insert the "-" for the uuid format
uuid_string( 1: 8) = tmp_string( 1: 8)
uuid_string( 9: 9) = '-'
uuid_string(10:13) = tmp_string( 9:12)
uuid_string(14:14) = '-'
uuid_string(15:18) = tmp_string(13:16)
uuid_string(19:19) = '-'
uuid_string(20:23) = tmp_string(17:20)
uuid_string(24:24) = '-'
uuid_string(25:36) = tmp_string(21:32)

END SUBROUTINE uuid_2char

!==============================================================================
!==============================================================================

FUNCTION uuid_match (uuid1, uuid2)
  LOGICAL :: uuid_match
  CHARACTER(LEN=1), INTENT(IN)  :: uuid1(16)
  CHARACTER(LEN=1), INTENT(IN)  :: uuid2(16)

  uuid_match = .TRUE.
!NEC$ novector
  IF (ANY(uuid1(:) /= uuid2(:))) THEN
    uuid_match = .FALSE.
  ENDIF
END FUNCTION uuid_match

!==============================================================================
!==============================================================================
! Convert single byte to 'hexadecimal' string

PURE FUNCTION byte2hex (c) RESULT (hex)

CHARACTER(LEN=1), INTENT(IN) :: c
CHARACTER(LEN=2)             :: hex
INTEGER :: x

  x        = IACHAR (c)
  hex(1:1) = nibble (      x / 16)
  hex(2:2) = nibble (IAND (x,  15))

END FUNCTION byte2hex

!==============================================================================
!==============================================================================
! Convert single byte to 'hexadecimal' string
! Convert 'nibble' to 'hexadecimal'

PURE FUNCTION nibble (x)

INTEGER, INTENT(IN) :: x
CHARACTER           :: nibble

  SELECT CASE (x)
  CASE (0:9)
     nibble = ACHAR (IACHAR ('0') + x)
  CASE DEFAULT
     nibble = ACHAR (IACHAR ('a') - 10 + x)
  END SELECT

END FUNCTION nibble

!==============================================================================
!==============================================================================

SUBROUTINE dealloc_refatm_structure (ierror)
  INTEGER, INTENT(OUT) :: ierror

  INTEGER :: i

  ierror = 0

  DO i =1, imax_refatmtype
    DEALLOCATE(refatm_defaults(i)%p0sl,        STAT=ierror)
    DEALLOCATE(refatm_defaults(i)%t0sl,        STAT=ierror)
    DEALLOCATE(refatm_defaults(i)%dt0lp,       STAT=ierror)
    DEALLOCATE(refatm_defaults(i)%delta_t,     STAT=ierror)
    DEALLOCATE(refatm_defaults(i)%h_scal,      STAT=ierror)
    DEALLOCATE(refatm_defaults(i)%bvref,       STAT=ierror)
  ENDDO

  DEALLOCATE(refatm_defaults,                  STAT=ierror)

  DEALLOCATE(refatm_in%p0sl,                 STAT=ierror)
  DEALLOCATE(refatm_in%t0sl,                 STAT=ierror)
  DEALLOCATE(refatm_in%dt0lp,                STAT=ierror)
  DEALLOCATE(refatm_in%delta_t,              STAT=ierror)
  DEALLOCATE(refatm_in%h_scal,               STAT=ierror)
  DEALLOCATE(refatm_in%bvref,                STAT=ierror)

  DEALLOCATE(refatm%p0sl,                    STAT=ierror)
  DEALLOCATE(refatm%t0sl,                    STAT=ierror)
  DEALLOCATE(refatm%dt0lp,                   STAT=ierror)
  DEALLOCATE(refatm%delta_t,                 STAT=ierror)
  DEALLOCATE(refatm%h_scal,                  STAT=ierror)
  DEALLOCATE(refatm%bvref,                   STAT=ierror)

END SUBROUTINE dealloc_refatm_structure

!==============================================================================
!==============================================================================

SUBROUTINE dealloc_vcoord_structure (ierror)
  INTEGER, INTENT(OUT) :: ierror

  INTEGER :: i

  ierror = 0

  DO i =1, imax_vcoordtype
    DEALLOCATE(vcoord_defaults(i)%vert_coord, STAT=ierror)
    DEALLOCATE(vcoord_defaults(i)%sigm_coord, STAT=ierror)
    DEALLOCATE(vcoord_defaults(i)%vcflat,     STAT=ierror)
  ENDDO
  DEALLOCATE(vcoord_defaults,                 STAT=ierror)

  DEALLOCATE(vcoord_in%vert_coord,          STAT=ierror)
  DEALLOCATE(vcoord_in%sigm_coord,          STAT=ierror)
  DEALLOCATE(vcoord_in%vcflat,              STAT=ierror)

  DEALLOCATE(vcoord%vert_coord,             STAT=ierror)
  DEALLOCATE(vcoord%sigm_coord,             STAT=ierror)
  DEALLOCATE(vcoord%vcflat,                 STAT=ierror)

END SUBROUTINE dealloc_vcoord_structure

!==============================================================================

END MODULE vgrid_refatm_utils
