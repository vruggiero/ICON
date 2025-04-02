!+ Source module for COSMO observation operators for conventional observations
!-------------------------------------------------------------------------------

MODULE src_obs_operator_conv

!-------------------------------------------------------------------------------
!
! Description:
!   This module 'src_obs_operator_conv' contains the routines which perform the
!   COSMO forward observation operators for the conventional observations that
!   are also used by the nudging scheme in COSMO. This includes:
!    - multi-level reports (radiosondes, profilers, aircraft)
!    - upper-air single-level reports (aircraft, (satob))
!    - single-level surface reports (synop, ship, buoy)
!    - GPS zenith total path delay (ZPD) (and IWV) reports
!    - (prepared:) satellite retrievals
!
!   Note: This module belongs to a group of modules for the 'COSMO observation
!   ----  operators including quality control for conventional observations'.
!         It is shared between the COSMO model and the DACE program package !
!
! Method:
!   This module contains the following procedures:
!   - q2rh_col             : generalized relative humidity from specific humid.
!   - tvirt_col            : computes model columns of virtual temperature
!   - hhl_col              : computes model height on model half levels
!   - lhyphl_col           : computes LOG of hydrostatic pressure on half levels
!   - lapse_rate_levels    : searches model levels to calculate the lapse rate
!   - prep_vi_mo2ob_sg     : prepares interpol. from model levels to obs levels
!   - prep_vi_mo2ob        : prepares interpol. from model levels to obs levels
!   - sing_obs_operator    : forward obs operator for upper-air single-level obs
!   - surf_obs_operator    : forward obs operator for surface-level obs
!   - surf_obs_opr_1h      : forward obs operator for 1-hr sums (precip, rad.)
!   - cloud_obs_operator   : forward obs operator for cloud obs
!   - zpd_iwv_obs_operator : obs opr for GPS IWV and zenith path delay (ZPD)
!   - ps_obs_operator_sing : obs opr for surface pressure from a surface station
!   - ps_obs_operator_mult : obs opr for surf pressure from a multi-level report
!   - mult_obs_operator    : forward observation operator for multi-level obs
!   - mult_obs_2_modlev    : inverse observation operator for multi-level obs
!   - mult_vert_ipol_obs   : vertical interpolation of the obs to model level
!   - mult_obs_operator_z  : supplementary obs operator for multi-level T, z
!   - frac2octas           : conversion of fraction into octas
!   - conv_zpd_iwv         : conversion of ZPD into IWV (integr. water vapour)
!   - conv_iwv_2_zpd       : conversion of IWV into ZPD (zenith path delay)
!
!   All routines have to be called by external driving routines, except for:
!   - lapse_rate_levels   : called by ps_obs_operator_sing, ps_obs_operator_mult
!   - prep_vi_mo2ob_sg    : called by sing_obs_operator
!   - mult_vert_ipol_obs  : called by mult_obs_2_modlev
!
!   This module also contains elemental functions, formerly statement functions:
!   - fpvsw        : Magnus formula for water: saturation vapour pressure from T
!   - ftd          : inverse of Magnus formula: dewpoint T from vapour pressure
!   - fpv2q        : specific humidity from vapour pressure and air pressure
!   - fq2pv        : vapour pressure from specific humidity and air pressure
!
!   It uses from:
!    - utilities:    uv2df
!
!   Data modules used:
!    - kind_parameters
!    - data_obs_lib_cosmo
!    - data_obs_record
!    - mo_fdbk_tables
!
!   Note: This module does not contain any communication between PE's.
!
!   Note on the interface of this module:
!   -------------------------------------
!   All the variables used from external data-modules are parameters.
!   Only the the subroutine (input) arguments, have to be set to their correct
!   values in the calling routines outside this module.
!
! Current Code Owner (for COSMO and for DWD DACE):
!  DWD, Christoph Schraff
!  phone:  +49  69  8062 2725
!  fax:    +49  69  8062 3721
!  email:  christoph.schraff@dwd.de
!
! History:
! Version    Date       Name
! -------    ----       ----
! V4_28        2013/07/12 Christoph Schraff
!  Initial release, based on routines in modules 'src_mult_local.f90' and
!                   'src_sing_local.f90' of COSMO V4_22.
!  Previous major milestones in the developement of 'src_mult_local.f90' and
!  'src_sing_local.f90' related to the procedures in the current module:
!   - 1998 : Initial release.
!   - 2000 : Surface-level obs increments optionally derived by use of a
!            new surface layer parameterization.
!   - 2002 : Option for obs increments as differences of specific humidity.
!   - 2004 : Option to process satellite retrievals.
!   - 2012 : Introduction of observation operator of virtual temperature.
!  Changes with respect to V4_22:
!   - More modular code: Split-up of old large routine into new routines with
!     slim interfaces, e.g. using only data modules common to COSMO and 3DVAR-
!     LETKF package.
!   - Simulated observations instead of observation increments stored in SOR
!     array 'zsobdy', used for the feedback files. Minor bug fixes.
!   - The height observation increment at the lowest z-p-obs of a multi-level
!     report, that is written to the SOR (for LETKF), is now the one computed
!     in 'ps_obs_operator_mult' rather than 'mult_obs_operator(_z)'.
!   - Upper-air height obs, even from radiosondes, are now set to passive,
!     except for the lowest radiosonde height obs.
!     (Switch 'luse_mlz' prepares option for active use of upper-air height.)
!   - Statement functions replaced by elemental or intrinsic functions.
! V5_1         2014-11-28 Christoph Schraff, Oliver Fuhrer
!  Bug fix, to make if conditions containing 'epsy' more consistent.
!  Replaced ireals by wp (working precision) (OF)
! V5_3         2015-10-09 Christoph Schraff
!  Obs operator for wind speed, wind direction, and dewpoint temperature added
!  for surface observations.
! V5_4c        2016-10-06 Klaus Stephan
!  PATCH 20160922 to avoid division by zero in the very rare case of zpsvob = 0.
! V5_4d        2016-12-12 Christoph Schraff
!  - PATCH 20160922 removed, as zero surface pressure is excluded now already
!    when reading the obs (in module src_obs_cdfin_sing).
!  - Allow for IWV increments, if only IWV (and not ZTD) is observed
!    (modified interface to 'zpd_iwv_obs_operator').
! V5_4f        2017-09-01 Christoph Schraff
!  Update from DACE V1_51: Editorial changes only.
! V5_6b        2019-10-16 Christoph Schraff + Andreas Rhodin
!  - Obs operator 'surf_obs_opr_1h' added for 1-hour sums of precipitation and
!    radiation.
!  - Bug fix at computation of cloud base height (replacing cloud top height).
!  - Set subroutines and functions used by other modules explicitly to public
!    (required by ICON Compiler directives when modules are ported to ICON).
!  - All 'iintegers' removed and replaced by standard integers.
!  - Use of 'kind_parameters' instead of 'data_parameters'.
! @VERSION@    @DATE@     Sven Ulbrich
!  Computation of Td2m model equivalent added if RH2m is not observed.
! - 21.04.2022: Corrected loop range in 'prep_vi_mo2ob_sg' (for safety only).
! - 19.05.2022: Obs operator for ceiling added in 'cloud_obs_operator'.
! - 02.08.2022: Obs operator for visibility added in 'cloud_obs_operator'.
! - 24.04.2023: Routine 'shift_profile' added for processing of tower profiles.
! - 13.06.2023: 'lapse_rate_levels': if insufficient model levels, adjust choice
!               of levels and write warning message.
!
! ==============================================================================
! DACE history:
! ------------
!  1.51       2017-02-24 Andreas Rhodin : update shared modules from COSMO 5.04d
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
!===============================================================================
!
! Declarations:
!
! Modules used:
!
!-------------------------------------------------------------------------------

USE kind_parameters, ONLY :   &
    wp           ! KIND-type parameter for real variables

!-------------------------------------------------------------------------------

USE data_obs_lib_cosmo , ONLY :   &

! 1. General parameters
! ---------------------

    c0         ,& ! standard real constant 0.0
    c1         ,& ! standard real constant 1.0
    c2         ,& ! standard real constant 2.0
    c05        ,& ! standard real constant 0.5
!   c3600      ,& ! standard real constant 3600.0
    rmdi       ,& ! =-1.E31_wp : commonly used missing data indicator
    rmdich     ,& ! =-1.E30_wp : commonly used check value for miss data
    epsy       ,& ! = 1.E-8_wp : commonly used very small value > 0
    luse_mlz      ! if false then use multi-level T, not z, and set z-obs to
                  !          passive except for the lowest z-p-obs

! 4. I/O device numbers for obs processing / nudging and file names
! -----------------------------------------------------------------

!   nupr          ! unit number of file for all the remaining information

! end of data_obs_lib_cosmo

!-------------------------------------------------------------------------------

USE data_obs_record , ONLY :   &

! 1.1 ODR report header format
! ----------------------------

!   mxrhed       ,& ! header length of multi-level reports
!   nhilon       ,& ! longitude of observing station
!   nhjlat       ,& ! latitude  of observing station
!   nhalt        ,& ! station altitude [m]
!   nhtime       ,& ! time of observat. in forecast hours
!   nhzio        ,& ! latitude  of obs. station (or lowest datum) in g.pt units
!   nhzjo        ,& ! longitude of obs. station in grid pt. units
!   nhtvip       ,& ! observed multi-lev pressure interpol to lowest model level

!   mxrhdf       ,& ! header length of multi-level reports
!   mxghdf       ,& ! header length of GPS reports
!   mxthdf       ,& ! header length of satellite retrieval reports
!   nhio         ,& ! (local) x-coord. of grid pt. assigned to obs
!   nhjo         ,& ! (local) y-coord. of grid pt. assigned to obs
!   nhitot       ,& ! global x-coord. of grid pt. assigned to obs
!   nhjtot       ,& ! global y-coord. of grid pt. assigned to obs
!   nhobtp       ,& ! observation type
!   nhcode       ,& ! code type
!   nhschr       ,& ! station characteristics (packed as in VOF)
!   nhpass       ,& ! flag for report being set to 'passive' (as in VOF)
!   nhqcfw       ,& ! status of QC and of writing to feedobs files (+ p-QC flag)
!   nhnlev       ,& ! number of obs. levels (for multi-level reports)

! 1.3 ODR body format
! -------------------

!   maxrsl       ,& ! max. number of levels in multi-level ODR
    mxrbdy       ,& ! body length of multi-level reports
    nbtu         ,& ! u wind component [m/s]
    nbtv         ,& ! v wind component [m/s]
    nbtt         ,& ! temperature [K]
    nbtrh        ,& ! relative humidity [/]
    nbtp         ,& ! pressure [Pa]
    nbtz         ,& ! height [m]
    nbtlop       ,& ! LOG( pressure )
    mxrbdf       ,& ! body length of multi-level reports
    nbtflg       ,& ! main flag word          (bit pattern, see below: 'nb?flg')
    nbterr       ,& ! pre-processing status flags  (bit p., see below: 'nb?err')
    nbtqcf       ,& ! threshold quality control flags      (see below: 'nb?qcf')
    nbtlsg       ,& ! level id (bit pattern, as in NetCDF statistics file)
    nbtlid       ,& ! level identity          (bit pattern, see below: 'nb?lid')
    mxsbdy       ,& ! body length of single-level reports
    nbsu         ,& ! u wind component [m/s]
    nbsv         ,& ! v wind component [m/s]
    nbst         ,& ! temperature [K]
    nbsrh        ,& ! relative humidity [/]
    nbsp         ,& ! pressure [Pa]
    nbsz         ,& ! height [m]
    nbscbs       ,& ! (lowest) cloud base height CBH above surface       [m]
    nbscil       ,& ! ceiling (CBH of lowest cloud layer >= 5 octas)     [m]
    nbscl        ,& ! low       cloud cover        (BUFR Table 020011)   [octas]
    nbscm        ,& ! mid-level cloud cover        (BUFR Table 020011)   [octas]
    nbsch        ,& ! high      cloud cover        (BUFR Table 020011)   [octas]
    nbsct        ,& ! total     cloud cover        (BUFR Table 020011)   [octas]
    nbsvis       ,& ! visibility                                         [m]
    nbsff        ,& ! wind speed                                         [m/s]
    nbsdd        ,& ! wind direction                                     [deg]
    nbstd        ,& ! dewpoint temperature                               [K]
    nbsrr1       ,& ! precipitation amount over 1 hour                   [mm]
    nbsrad       ,& ! global  downward solar radiation, sum over 1 hour  [J/m2]
    nbsrdd       ,& ! diffuse downward solar radiation, sum over 1 hour  [J/m2]
    nbsrdl          ! long-wave downward radiation    , sum over 1 hour  [J/m2]

USE data_obs_record , ONLY :   &

! 1.4 Bit patterns for packed information in ODR body
! ---------------------------------------------------

!   nvfubp       ,& ! bit pos. for main flag on wind                  nb?flg
!   nvftbp       ,& ! bit pos. for main flag on temperature             "
!   nvfqbp       ,& ! bit pos. for main flag on humidity                "
    nvfzbp       ,& ! bit pos. for main flag on pressure / geopot.      "
!   nvfgbp       ,& ! bit pos. for main flag on IWV / ZPD               "
    nvfaoc       ,& ! no. of bits occ. by each main flag                "
    nvfbps       ,& ! bit pos. for main flags: 4: above 300 hPa         "
!   nvfboc       ,& ! no. of bits occ. for main flags                   "
    nvflbp       ,& ! bit pos. for level flag: level below surface      "
    nvlidp       ,& ! level id. bit pattern                           nb?lid
    nvlido       ,& ! no. bits occ. by each indicator in level id.      "

!       1.4.2  Other bit patt. for packed info in ODR (VOF) body, general words
!              ----------------------------------------------------------------
    nvru         ,& ! bit pos. for status/QC flags for horiz. wind nb?err/nb?qcf
    nvrt         ,& ! bit pos. for status/QC flags for temperature       "
    nvrq         ,& ! bit pos. for status/QC flags for humidity          "
    nvrz         ,& ! bit pos. for status/QC flags for pressure/height   "
!   nvriwv       ,& ! bit pos. for status/QC flags for IWV               "
!   nvrzpd       ,& ! bit pos. for status/QC flags for ZPD               "
!   nvrzbc       ,& ! bit pos. for temporary flag: QC ag. LBC pressure   "

! 1.5 Further quantities related to ODR
! -------------------------------------

    ystid        ,& ! obs. station identity to be printed
    fdoro        ,& ! scaling factor to vertical distances betw. model orography
                    ! and station height for (z-obs < z-mod)

! 2. Observation data records (ODR)
! ---------------------------------

!   omlbdy       ,& ! body   of multi-level ODR
!   omlhed       ,& ! header of multi-level ODR
!   momlbd       ,& ! body   of multi-level ODR

! 3. Simulated Observation Record (SOR)
! -------------------------------------

    mxsoml       ,& ! SOR body length for multi-level reports
    mxsosg       ,& ! SOR body length for single-level reports
    nso_u        ,& ! u wind component                               [m/s]
    nso_v        ,& ! v wind component                               [m/s]
    nso_t        ,& ! temperature                                    [K]
    nso_rh       ,& ! relative humidity                              [%]
    nso_p        ,& ! pressure (sfc.) | geopotential (upper-air)  [Pa] | [m2/s2]
    nso_ct       ,& ! total cloud cover                              [octas]
    nso_cl       ,& ! low cloud cover                                [octas]
    nso_cm       ,& ! mid-level cloud cover                          [octas]
    nso_ch       ,& ! high cloud cover                               [octas]
    nso_cbs      ,& ! cloud base height CBH above surface            [m]
    nso_cil      ,& ! ceiling (CBH of lowest cloud layer >= 5 octas) [m]
    nso_vis      ,& ! visibility                                     [m]
    nso_ff       ,& ! wind speed                                     [m/s]
    nso_dd       ,& ! wind direction                                 [deg]
    nso_td       ,& ! dewpoint temperature                           [K]
    nso_rr1      ,& ! precipitation                , sum over 1 hour [mm]
    nso_rad      ,& ! global downw. solar radiation, sum over 1 hour [J/m2]
    nso_rdd      ,& ! diffuse down. solar radiation, sum over 1 hour [J/m2]
    nso_rdl      ,& ! long-wave downward  radiation, sum over 1 hour [J/m2]
    nso_iq       ,& ! integrated water vapour (increment)            [mm]
    nso_zpd         ! zenith total path delay                        [mm]

! end of data_obs_record

!-------------------------------------------------------------------------------

USE mo_fdbk_tables, ONLY :   &

    OT_AIREP     ,& ! observation type for aircraft obs
    OT_TEMP      ,& ! observation type for TEMP radiosonde
    OT_PILOT     ,& ! observation type for PILOT (incl. wind profiler, RASS,...)
    OT_SATEM     ,& ! observation type for satellite retrievals
    LS_SURFACE      ! level significance: surface level

! end of mo_fdbk_tables

!-------------------------------------------------------------------------------

 USE utilities,                ONLY :  &
    uv2df           ! Converts the wind components to wind direction and speed

!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! Local declarations
!-------------------------------------------------------------------------------

IMPLICIT NONE
PRIVATE
          ! public subroutines and functions
PUBLIC :: q2rh_col,  tvirt_col,  hhl_col,  lhyphl_col,  shift_profile,         &
          sing_obs_operator,  surf_obs_operator,  surf_obs_opr_1h,             &
          cloud_obs_operator,  frac2octas,                                     &
          ps_obs_operator_sing,  ps_obs_operator_mult,  mult_obs_operator,     &
          mult_obs_2_modlev,  mult_obs_operator_z

! for nudging in COSMO model code only:
PUBLIC :: prep_vi_mo2ob,  zpd_iwv_obs_operator

!===============================================================================

!-------------------------------------------------------------------------------
! Public and Private Subroutines
!-------------------------------------------------------------------------------

CONTAINS



!-------------------------------------------------------------------------------
!+ Module procedure for computing relative humidity from specific humidity
!-------------------------------------------------------------------------------

SUBROUTINE q2rh_col ( ke, col_t, col_qv, col_qc, col_p, rdv, b1, b2w, b3, b4w  &
                    , nupr , col_rh )

!-------------------------------------------------------------------------------
! Description:
!   This routine computes (generalized) relative humdity from specific content
!   of water vapour (and cloud water).
!
! Method:
!   Use of Magnus formula over water.
!-------------------------------------------------------------------------------
! Modules used: These are declared in the module declaration section
! -------------

IMPLICIT NONE

! Subroutine arguments:
! --------------------

  INTEGER        , INTENT (IN)         ::       &
    ke             ,& ! number of vertical (main) levels in model column
    nupr              ! file unit number for control output
                      !   (if < 0 then no control output is written)

  REAL (KIND=wp) , INTENT (IN)         ::       &
                      ! model column (at the observation location) of:
    col_t   (ke)   ,& ! temperature                                    [  K  ]
    col_qv  (ke)   ,& ! specific water vapor content                   [kg/kg]
    col_qc  (ke)   ,& ! specific cloud water content                   [kg/kg]
    col_p   (ke)   ,& ! pressure (full value)        on main levels    [ Pa  ]
    rdv            ,& ! ratio of gas constant for dry air 'r_d'
                      !      and gas constant for water vapour 'r_v':  r_d/r_v
    b1             ,& ! \  constants for computing
    b2w            ,& !  \ the saturation vapour pressure
    b3             ,& !  / over water
    b4w               ! /  by the Magnus formula

  REAL (KIND=wp) , INTENT (OUT)        ::       &
    col_rh  (ke)      ! relative humidity                              [     ]

! Local parameters:
! ----------------

  REAL (KIND=wp) , PARAMETER  :: &
    c1  = 1.0_wp    ! standard real constant 1.0

! Local scalars:
! -------------

  LOGICAL        , SAVE ::  &
    lfirst = .TRUE.     ! true only if routine is called for the first time

  INTEGER        ::  &
    km                  ! vertical loop index over model levels

  REAL (KIND=wp) ::  &
    zpv       (ke)   ,& ! water vapour pressure
    ztd                 ! dewpoint temperature

!------------ End of header ----------------------------------------------------

!-------------------------------------------------------------------------------
! Begin Subroutine q2rh_col
!-------------------------------------------------------------------------------

  DO km = 1 , ke
    zpv (km) = fq2pv( col_qv(km) + col_qc(km) , col_p(km) , rdv )

    !   simulated relative humidity on model values
    col_rh (km) = zpv(km) / fpvsw( col_t(km), b1, b2w, b3, b4w )
  ENDDO

  IF ((lfirst) .AND. (nupr >= 0)) THEN
    WRITE( nupr,'("Example for dewpoint temperature:")' )
    DO km = 1 , ke
      ztd = ftd( MAX( zpv(km), 0.0001_wp), b1, b2w, b3, b4w )
      WRITE( nupr,'(3F9.1,F9.2,I4)' )                                          &
             col_t(km), ztd, col_p(km), col_rh(km), km
    ENDDO
    lfirst = .FALSE.
  ENDIF

!-------------------------------------------------------------------------------
! End of module procedure q2rh_col
!-------------------------------------------------------------------------------

END SUBROUTINE q2rh_col



!-------------------------------------------------------------------------------
!+ Module procedure for computing model columns of virtual temperature
!-------------------------------------------------------------------------------

SUBROUTINE tvirt_col ( ke, col_t, col_qv, col_qc, col_qrs, rdv ,  col_tv )

!-------------------------------------------------------------------------------
! Description:
!   This routine computes virtual temperature from temperature and specific
!   contents of water vapour and hydrometeors.
!
! Method:
!   Taking into account the water loading by hydrometeors.
!-------------------------------------------------------------------------------
! Modules used: These are declared in the module declaration section
! -------------

IMPLICIT NONE

! Subroutine arguments:
! --------------------

  INTEGER        , INTENT (IN)         ::       &
    ke                ! number of vertical (main) levels in model column

  REAL (KIND=wp) , INTENT (IN)         ::       &
                      ! model column (at the observation location) of:
    col_t   (ke)   ,& ! temperature                                    [  K  ]
    col_qv  (ke)   ,& ! specific water vapor content                   [kg/kg]
    col_qc  (ke)   ,& ! specific cloud water content                   [kg/kg]
    col_qrs (ke)   ,& ! spec. cont. of hydrometeors excl. cloud water  [kg/kg]
    rdv               ! ratio of gas constant for dry air 'r_d'
                      !      and gas constant for water vapour 'r_v':  r_d/r_v

  REAL (KIND=wp) , INTENT (OUT)        ::       &
    col_tv  (ke)      ! virtual temperature                            [  K  ]

! Local parameters:
! ----------------

  REAL (KIND=wp) , PARAMETER  :: &
    c1  = 1.0_wp    ! standard real constant 1.0

! Local scalars:
! -------------

  INTEGER        ::  &
    km                ! vertical loop index over model levels

  REAL (KIND=wp) ::  &
    rvd_m_o           ! r_v / r_d  -  1.0

!------------ End of header ----------------------------------------------------

!-------------------------------------------------------------------------------
! Begin Subroutine tvirt_col
!-------------------------------------------------------------------------------

  rvd_m_o  =  c1 / rdv  -  c1
  DO km = 1 , ke
    col_tv (km) = col_t(km) *(c1 + rvd_m_o *col_qv(km)                         &
                                 -          col_qc(km) - col_qrs(km))
!   !   simulated virtual (RASS) or dry bulb (other) temperature
!   IF (lvirt) THEN
!     zftv2t_m (km) = c1 / (c1 + rvd_m_o *col_qv(km) - col_qc(km) - col_qrs(km))
!     ztcd     (km) = col_t(km) / zftv2t_m(km)
!   ELSE
!     zftv2t_m (km) = c1
!     ztcd     (km) = col_t(km)
!   ENDIF
  ENDDO

!-------------------------------------------------------------------------------
! End of module procedure tvirt_col
!-------------------------------------------------------------------------------

END SUBROUTINE tvirt_col



!-------------------------------------------------------------------------------
!+ Module procedure for computing height on model half levels
!-------------------------------------------------------------------------------

SUBROUTINE hhl_col ( ke, col_z, hsurf , col_hhl )

!-------------------------------------------------------------------------------
! Description:
!   This routine computes LOG of hydrostatic pressure on half levels from
!   pressure and virtual temperature at main level and height at half levels.
!
! Method:
!   Integration of hydrostatic equation (approximate version).
!-------------------------------------------------------------------------------
! Modules used: These are declared in the module declaration section
! -------------

IMPLICIT NONE

! Subroutine arguments:
! --------------------

  INTEGER        , INTENT (IN)         ::       &
    ke                ! number of vertical (main) levels in model column

  REAL (KIND=wp) , INTENT (IN)         ::       &
                      ! model column (at the observation location) of:
    col_z   (ke)   ,& ! geometrical height           of main levels    [  m  ]
    hsurf             ! model orography                                [  m  ]

  REAL (KIND=wp) , INTENT (OUT)        ::       &
    col_hhl (ke+1)    ! geometrical height           of half levels    [  m  ]

! Local parameters: None
! ----------------

! Local scalars:
! -------------

  INTEGER        ::  &
    kk                ! vertical loop index over model levels

!------------ End of header ----------------------------------------------------

!-------------------------------------------------------------------------------
! Begin Subroutine hhl_col
!-------------------------------------------------------------------------------

  col_hhl (ke+1)  =  hsurf

  DO kk = ke , 1 , -1
    col_hhl (kk)  =  col_z(kk)  +  (col_z(kk) - col_hhl(kk+1))
  ENDDO

!-------------------------------------------------------------------------------
! End of module procedure hhl_col
!-------------------------------------------------------------------------------

END SUBROUTINE hhl_col


!-------------------------------------------------------------------------------
!+ Module procedure for computing LOG of hydrostatic pressure on half levels
!-------------------------------------------------------------------------------

SUBROUTINE lhyphl_col ( ke, col_p, col_tv, col_hhl, r_d, r_g , col_lhyphl )

!-------------------------------------------------------------------------------
! Description:
!   This routine computes LOG of hydrostatic pressure on half levels from
!   pressure and virtual temperature at main level and height at half levels.
!
! Method:
!   Integration of hydrostatic equation (approximate version).
!-------------------------------------------------------------------------------
! Modules used: These are declared in the module declaration section
! -------------

IMPLICIT NONE

! Subroutine arguments:
! --------------------

  INTEGER        , INTENT (IN)         ::       &
    ke                ! number of vertical (main) levels in model column

  REAL (KIND=wp) , INTENT (IN)         ::       &
                      ! model column (at the observation location) of:
    col_p   (ke)   ,& ! pressure (full value)        on main levels    [ Pa  ]
    col_tv  (ke)   ,& ! virtual temperature                            [  K  ]
    col_hhl (ke+1) ,& ! geometrical height           of half levels    [  m  ]
    r_d            ,& ! gas constant for dry air       (287.05)
    r_g               ! acceleration due to gravity    (  9.80665)

  REAL (KIND=wp) , INTENT (OUT)        ::       &
    col_lhyphl (ke+1)   ! LOG( hydrostatic pressure [pa] ) on half levels  [ - ]

! Local parameters: None
! ----------------

! Local scalars:
! -------------

  INTEGER        ::  &
    kk                ! vertical loop index over model levels

  REAL (KIND=wp) ::  &
    zhyphl (ke+1)  ,& ! hydrostatic pressure [pa] ) on half levels
    rdg               ! r_d / r_g

!------------ End of header ----------------------------------------------------

!-------------------------------------------------------------------------------
! Begin Subroutine lhyphl_col
!-------------------------------------------------------------------------------

  rdg   = r_d / r_g

  zhyphl (1)  =  col_p(1)  -  c05* (col_hhl(1) - col_hhl(2))                   &
                                  *col_p(1) / (rdg *col_tv(1))
  DO kk = 1 , ke
    zhyphl (kk+1)  =  zhyphl(kk) +  (col_hhl(kk) - col_hhl(kk+1))              &
                                   *col_p(kk) / (rdg *col_tv(kk))
  ENDDO

  DO kk = 1 , ke + 1
    col_lhyphl (kk)  =  LOG( zhyphl(kk) )
  ENDDO

!-------------------------------------------------------------------------------
! End of module procedure lhyphl_col
!-------------------------------------------------------------------------------

END SUBROUTINE lhyphl_col



!-------------------------------------------------------------------------------
!+ Module procedure for search for model levels used to calculate the lapse rate
!-------------------------------------------------------------------------------

SUBROUTINE lapse_rate_levels ( ke, col_z , klu, kll , hlapsu_in, hlapsl_in )

!-------------------------------------------------------------------------------
!
! Description:
!   The routine returns the indices of 2 model main levels which are close to
!   2 input heights above the model orography.
!
! Method:
!   Simple search by looping over model levels.
!
! Written by        :  Christoph Schraff, DWD  (original version: 18.09.97)
! Current Code Owner:  Christoph Schraff, DWD
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: European Standards for Writing and Documenting
!                       Exchangeable Fortran 90 Code + DWD Extensions.
!
!-------------------------------------------------------------------------------
!
! Declarations:
!
!-------------------------------------------------------------------------------
! Modules used: These are declared in the module declaration section
! -------------

IMPLICIT NONE

! Subroutine arguments:
! --------------------

  INTEGER        , INTENT (IN)         ::       &
    ke             ! number of vertical (main) levels in model column

  REAL (KIND=wp) , INTENT (IN)         ::       &
                   ! model column (at the observation location) of:
    col_z  (ke)    ! geometrical height           of main levels    [  m  ]

  INTEGER        , INTENT (OUT)        ::       &
    klu         ,& ! upper model main level used to compute lapse rate
    kll            ! lower model main level used to compute lapse rate

  REAL (KIND=wp) , INTENT (IN)    , OPTIONAL   ::       &
    hlapsu_in   ,& ! height [m] above model orography of upper level
    hlapsl_in      ! height [m] above model orography of lower level
                   !   used to compute lapse rate

! Local parameters: None
! ----------------

! Local scalars:
! -------------

  INTEGER        ::  &
    kk                ! loop incdex

  INTEGER , SAVE ::  &
    npr     =  0      ! loop incdex

  REAL (KIND=wp) ::  &
    hlapsu  =  1500._wp  ,& ! height [m] above model orography of upper
                                ! level used to compute lapse rate
    hlapsl  =   500._wp     ! height [m] above model orography of lower
                                ! level used to compute lapse rate

!------------ End of header ----------------------------------------------------


!-------------------------------------------------------------------------------
! Begin Subroutine lapse_rate_levels
!-------------------------------------------------------------------------------

  IF ((PRESENT( hlapsu_in )) .AND. (PRESENT( hlapsl_in ))) THEN
    hlapsu = hlapsu_in
    hlapsl = hlapsl_in
  ENDIF

  klu = 0
  get_lapse_levels: DO kk = ke-1 , 1 , -1
    IF (      (col_z(kk  ) - col_z(ke) >  hlapsl)                              &
        .AND. (col_z(kk+1) - col_z(ke) <= hlapsl))  kll = kk
    IF (      (col_z(kk  ) - col_z(ke) >  hlapsu)                              &
        .AND. (col_z(kk+1) - col_z(ke) <= hlapsu))  klu = kk
    IF (klu > 0)                                           EXIT get_lapse_levels
  ENDDO get_lapse_levels
  !   for safety, if insufficient model levels are input to MEC
  IF (klu == 0) THEN
    npr = npr + 1
    klu = 1
    kll = ke / 2
    IF (npr < 11)                                                              &
      WRITE(6,'("WARNING lapse_rate_levels: too few model levels",I4           &
              &," to compute proper lapse rate; use heights",F6.0,",",F6.0)' ) &
             ke, col_z(kll)-col_z(ke), col_z(klu)-col_z(ke)
  ENDIF

!-------------------------------------------------------------------------------
! End of module procedure lapse_rate_levels
!-------------------------------------------------------------------------------

END SUBROUTINE lapse_rate_levels



!-------------------------------------------------------------------------------
!+ Module procedure for search for model levels used to calculate the lapse rate
!-------------------------------------------------------------------------------

SUBROUTINE shift_profile ( zstalt, hsurf, ke, gacc, r_d, r_v, b1, b2w, b3, b4w &
                         , col_qx , col_lnp, col_p, col_z, col_t, col_qv )

!-------------------------------------------------------------------------------
!
! Description:
!   For particular observation reports (e.g. tower profiles) for which the model
!   equivalents should be computed at the same 'height of sensor above ground'
!   rather than the same geometric height or pressure as the observations,
!   the height and pressure levels of the model column are shifted vertically
!   to account for the difference between station height and model orography
!   as a first step of the forward (observation) operator.
!   This shifting procedure consists of the following steps:
!    - shift 'col_z' in the vertical such that height above model orography
!      becomes height above obs station
!    - apply a corresponding height correction to the temperature profile
!    - recompute pressure profile at new height levels
!    - correct specific water vapour content by preserving relative humidity
!
!   Note that this routine is not devised to be applied to reports with
!   observations (only) from the Prandtl layer, such as roadside sensors.
!   (Such observation types would require extrapolation and / or height
!    correction procedures instead of shifting of model profiles.)
!
!   Also note that 'hsurf' is not adjusted here but might be shifted too after
!   calling 'shift_profile'.
!
! Method:
!   See in-line description in the main section.
!
! Written by        :  Christoph Schraff, DWD  (original version: 31.07.20)
! Current Code Owner:  Christoph Schraff, DWD
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: European Standards for Writing and Documenting
!                       Exchangeable Fortran 90 Code + DWD Extensions.
!
!-------------------------------------------------------------------------------
!
! Declarations:
!
!-------------------------------------------------------------------------------
! Modules used: These are declared in the module declaration section
! -------------

IMPLICIT NONE

! Subroutine arguments:
! --------------------

  INTEGER        , INTENT (IN)         ::       &
    ke               ! number of vertical (main) levels in model column

  REAL (KIND=wp) , INTENT (IN)         ::       &
    zstalt        ,& ! station altitude
    hsurf         ,& ! model orography (at grid pt. assigned to obs)
    col_qx  (ke)  ,& ! specific content of (all) hydrometeors         [kg/kg]
    gacc          ,& ! acceleration due to gravity    (  9.80665)
    r_d           ,& ! gas constant for dry air       (287.05)
    r_v           ,& ! gas constant for water vapour  (461.51 J /(kg * K))
    b1            ,& ! \  constants for computing
    b2w           ,& !  \ the saturation vapour pressure
    b3            ,& !  / over water
    b4w              ! /  by the Magnus formula

  REAL (KIND=wp) , INTENT (INOUT)      ::       &
                     ! model column (at the observation location) of:
    col_p   (ke)  ,& ! pressure                    on main levels     [ Pa  ]
    col_lnp (ke)  ,& ! log( pressure )             on main levels     [ Pa  ]
    col_z   (ke)  ,& ! geometrical height          of main levels     [  m  ]
    col_t   (ke)  ,& ! temperature                                    [  K  ]
    col_qv  (ke)     ! specific water vapor content                   [kg/kg]

! Local parameters: None
! ----------------

! Local scalars:
! -------------

  INTEGER        ::  &
    klu              ,& ! upper model main level used to compute lapse rate
    kll              ,& ! lower model main level used to compute lapse rate
    kext             ,& ! top model level for which pressure has to extrapolated
    knew             ,& ! level index for shifted profile
    kk_sta           ,& ! level index where search for interpolat. levels starts
    kk                  ! loop index over levels

  REAL (KIND=wp) ::  &
                        ! shifted model column (at the observation location) of:
    col_p_new   (ke) ,& ! pressure                    on main levels     [ Pa  ]
    col_lnp_new (ke) ,& ! log( pressure )             on main levels     [ Pa  ]
    col_z_new   (ke) ,& ! geometrical height          of main levels     [  m  ]
    col_t_new   (ke) ,& ! temperature                                    [  K  ]
    col_qv_new  (ke) ,& ! specific water vapor content                   [kg/kg]
    zrh              ,& ! relative humidity                              [  1  ]
    zlapse           ,& ! lapse rate used for height correction of temperature
    zlapse_tv        ,& ! lapse rate of virtual temperature
    ztv         (2)  ,& ! virtual temperature at levels 'klu', 'kll'
    ztvke            ,& ! virtual temperature at lowest original level ('ke')
    zfrac            ,& ! interpolation weight
    z0      , z1     ,& ! auxilliary quantities used for extrapolation
    rdv              ,& ! ratio of gas constant for dry air 'r_d'
                        !      and gas constant for water vapour 'r_v':  r_d/r_v
    rvd_m_o             ! r_v / r_d  -  1.0

!------------ End of header ----------------------------------------------------


!-------------------------------------------------------------------------------
! Begin Subroutine shift_profile
!-------------------------------------------------------------------------------

  IF (ABS( zstalt-hsurf ) > 2._wp*epsy) THEN

    rdv      =  r_d / r_v
    rvd_m_o  =  c1 / rdv  -  c1

    ! correct (shift) height levels of model
    ! --------------------------------------
    !   col_z - hsurf: height above ground, to be added to station altitude
    col_z_new (:)  =  col_z(:) - hsurf + zstalt

    ! correct temperature profile according to height shift and lapse rate
    ! --------------------------------------------------------------------
    !   get lapse rate   (routine 'lapse_rate_levels' gets it from
    !                     approx. layer 500 - 1500 m above ground)
    CALL lapse_rate_levels ( ke, col_z(:) , klu, kll )
!   ======================
    ztv (1) = col_t(klu) * (1._wp + rvd_m_o*col_qv(klu) - col_qx(klu))
    ztv (2) = col_t(kll) * (1._wp + rvd_m_o*col_qv(kll) - col_qx(kll))
    zlapse_tv  =  (ztv(1)     - ztv(2)    ) / (col_z(klu) - col_z(kll))
    zlapse     =  (col_t(klu) - col_t(kll)) / (col_z(klu) - col_z(kll))
    !   temperature correction
    col_t_new (:)  =  col_t(:) + (zstalt - hsurf) *zlapse

    ! re-compute LOG( pressure ) levels
    ! ---------------------------------
    !   get target model levels [ke , kext] where extrapolation is needed
    !                           [kext+1, 1] where interpolation can be applied
    kext = ke + 1
    IF (zstalt < hsurf-epsy) THEN
      !   extrapolation is needed only for the (lowest) target level(s)
      !   which are below the lowest original level;
      !   this occurs only if zstalt < hsurf, i.e. if the levels are lowered
      !   (in the opposite case, when the levels are lifted, the interpolation
      !    formula is also used for extrapolating the top level(s) because those
      !    levels are not used anyway due to the absence of (e.g. tower) obs
      !    at these levels)
      DO knew = ke , 1 , -1
        IF (col_z_new(knew) < col_z(ke)-epsy) THEN
          kext = knew
        ELSE
          EXIT
        ENDIF
      ENDDO
    ENDIF
    !   carry out interpolation  (log(p) linear in z)
    kk = ke
    DO knew = kext-1 , 1 , -1
      kk_sta = kk
      find_orig_level: DO kk = kk_sta , 2 , -1
        IF ((col_z(kk-1) >= col_z_new(knew)-epsy) .OR. (kk == 2)) THEN
          zfrac = (col_z_new(knew) - col_z(kk)) / (col_z(kk-1) - col_z(kk))
          col_lnp_new (knew)  =  (1._wp-zfrac) * col_lnp(kk)                   &
                                +       zfrac  * col_lnp(kk-1)
          EXIT find_orig_level
        ENDIF
      ENDDO find_orig_level
    ENDDO
    !   carry out extrapolation, by integration of hydrostatic equation and 
    !                            assuming (non-zero!) constant lapse rate
    !                   (i.e. integrate from level ke:
    !                         d(ln(p)) = -g/R *1/(Tv(ke) + lapse*(z-z(ke)))*dz )
    !   (see also (German only!) wiki: Barometrische_Hoehenformel)
    IF (kext <= ke) THEN
      ztvke = col_t(ke) *(1._wp + rvd_m_o*col_qv(ke) - col_qx(ke))
      IF (ABS( zlapse_tv ) > 0.0001_wp) THEN
        z0    = - gacc /r_d /zlapse_tv
        z1    = zlapse_tv / ztvke
        DO knew = kext, ke
          col_lnp_new (knew)  =  col_lnp(ke)                                     &
                               + z0 *LOG( 1._wp + z1 *(  col_z_new(knew)         &
                                                       - col_z    (ke  )) )
        ENDDO
      ELSE
        !   near-zero lapse rate -> isothermal atm: d(ln(p)) = -g/R /Tv(ke) *dz
        z0    = - gacc /r_d /ztvke
        DO knew = kext, ke
          col_lnp_new (knew)  =  col_lnp(ke) + z0 *(col_z_new(knew) - col_z(ke))
        ENDDO
      ENDIF
    ENDIF

    ! correct specific water vapour content by preserving relative humidity
    ! ---------------------------------------------------------------------
    !   (qc is neglected here deliberately)
    DO kk = 1 , ke
      col_p_new  (kk) =  EXP( col_lnp_new(kk) )
      zrh             =  fq2pv( col_qv(kk) , col_p(kk) , rdv )                 &
                       / fpvsw( col_t (kk), b1, b2w, b3, b4w )
      col_qv_new (kk) =  fpv2q( zrh * fpvsw( col_t_new(kk), b1, b2w, b3, b4w ) &
                              , col_p_new(kk) , rdv )
    ENDDO

    ! assign new values to original variables
    ! and compute pressure
    ! ---------------------------------------
    DO kk = 1 , ke
      col_z   (kk)  =  col_z_new  (kk)
      col_t   (kk)  =  col_t_new  (kk)
      col_qv  (kk)  =  col_qv_new (kk)
      col_lnp (kk)  =  col_lnp_new(kk)
      col_p   (kk)  =  col_p_new  (kk)
    ENDDO
    ! note that 'hsurf' might be shifted too after calling 'shift_profile'
  ENDIF   !  (zstalt /= hsurf)

!-------------------------------------------------------------------------------
! End of module procedure shift_profile
!-------------------------------------------------------------------------------

END SUBROUTINE shift_profile



!-------------------------------------------------------------------------------
!+ Module procedure for preparing interpolation from model levels to obs levels
!-------------------------------------------------------------------------------

SUBROUTINE prep_vi_mo2ob_sg ( ke, col_p, zpob , kbz, zlpf )

!-------------------------------------------------------------------------------
! Description:
!   This routine prepares the vertical interpolation (linear in LOG( pressure ))
!   from model levels to 1 observation pressure level.
!
! Method:
!   Straightforward while search.
!-------------------------------------------------------------------------------
! Modules used: These are declared in the module declaration section
! -------------

IMPLICIT NONE

! Subroutine arguments:
! --------------------

  INTEGER        , INTENT (IN)         ::       &
    ke                ! number of vertical (main) levels in model column

  REAL (KIND=wp) , INTENT (IN)         ::       &
                      ! model column (at the observation location) of:
    col_p (ke)     ,& ! pressure (full value)    on main levels    [ Pa  ]
    zpob              ! observed pressure value

  INTEGER        , INTENT (OUT)        ::       &
    kbz               ! index of model level immediately below obs level
                      !   (if zero then obs is above top main model level)

  REAL (KIND=wp) , INTENT (OUT)        ::       &
    zlpf              ! weight factors for vertic. interpol. to obs levels
                      !   (weight of model level immediat. above obs lev.)

! Local parameters: None
! ----------------

! Local scalars:
! -------------

  INTEGER        ::  &
    kk             ,& ! vertical loop index over model levels
    ka                ! model level above and adjacent to the obs level

  REAL (KIND=wp) ::  &
    zlopob         ,& ! LOG( observation pressure )
    zlnp_b            ! log( pressure ) on model level 'kbz'

!------------ End of header ----------------------------------------------------

!-------------------------------------------------------------------------------
! Begin Subroutine prep_vi_mo2ob_sg
!-------------------------------------------------------------------------------

  !   do not extrapolate above top model main level
  IF (zpob <  col_p(1))  THEN
    kbz    = 0
    zlpf   = c0
  ELSE
    !   interpolate, or extrapolate below lowest model main level
    zlopob = LOG( zpob )
    kbz    = 2
    DO kk = 2 , ke-1
      IF (col_p(kk) <= zpob) kbz = kk + 1
    ENDDO
    ka = kbz - 1
    zlnp_b = LOG( col_p(kbz) )
    zlpf   = (zlnp_b - zlopob) / (zlnp_b - LOG( col_p(ka) ))
  ENDIF

!-------------------------------------------------------------------------------
! End of module procedure prep_vi_mo2ob_sg
!-------------------------------------------------------------------------------

END SUBROUTINE prep_vi_mo2ob_sg




!-------------------------------------------------------------------------------
!+ Module procedure for preparing interpolation from model levels to obs levels
!-------------------------------------------------------------------------------

SUBROUTINE prep_vi_mo2ob ( ke, col_p, col_lnp, nlev, zobbdy                    &
                         , kbotlev, ktoplev, kbz, zlpf )

!-------------------------------------------------------------------------------
! Description:
!   This routine prepares the vertical interpolation from model levels to
!   observation levels in 'mult_obs_operator_qc1' and 'mult_obs_qc_dz'
!   by computing auxilliary quantities (see routine intent(out) variables).
!
! Method:
!   Straightforward while search.
!-------------------------------------------------------------------------------
! Modules used: These are declared in the module declaration section
! -------------

IMPLICIT NONE

! Subroutine arguments:
! --------------------

  INTEGER        , INTENT (IN)         ::       &
    ke                   ,& ! number of vertical (main) levels in model column
    nlev                    ! number of vertical levels in report

  REAL (KIND=wp) , INTENT (IN)         ::       &
                            ! model column (at the observation location) of:
    col_p   (ke)         ,& ! pressure (full value)    on main levels    [ Pa  ]
    col_lnp (ke)         ,& ! log( pressure )          on main levels
    zobbdy  (nlev,mxrbdy)   ! multi-level obs. body (format: see 'omlbdy')

  INTEGER        , INTENT (OUT)        ::       &
    kbotlev              ,& ! number of obs levels below lowest    model level
    ktoplev              ,& ! number of obs levels above uppermost model level
    kbz     (nlev)          ! level indices for vertical interpol. to obs levels
                            !   (model level index immediately below obs level)

  REAL (KIND=wp) , INTENT (OUT)        ::       &
    zlpf   (nlev)           ! weight factors for vertic. interpol. to obs levels
                            !   (weight of model level immediat. above obs lev.)

! Local parameters:
! ----------------

  REAL (KIND=wp) , PARAMETER  :: &
    cdpext  =  101._wp  ! max. extent [Pa] of vertical extrapolation
                            ! of model values to obs. levels

! Local scalars:
! -------------

  INTEGER        ::  &
    ilev             ,& ! vertical loop index over observation levels
    ka                  ! model level above and adjacent to the obs level

!------------ End of header ----------------------------------------------------

!-------------------------------------------------------------------------------
! Begin Subroutine prep_vi_mo2ob
!-------------------------------------------------------------------------------

  !   obs levels below lowest model level
  ilev = 1
  DO WHILE (      (zobbdy(MIN( ilev,nlev ),nbtp)-cdpext > col_p(ke))           &
            .AND. (ilev <= nlev))
    kbz  (ilev)  =  ke
    ilev = ilev + 1
  ENDDO
  kbotlev = ilev - 1

  !   get model level indices
  ktoplev = 0
  ka = ke - 1

  DO ilev = 1+kbotlev , nlev
    DO WHILE ((col_lnp(ka) >= zobbdy(ilev,nbtlop)) .AND. (ka > 1))
      ka   = MAX( ka - 1 , 1 )
    ENDDO
    kbz  (ilev)  =  ka + 1
    IF (zobbdy(ilev,nbtp)+cdpext < col_p(1))  ktoplev = ktoplev + 1
  ENDDO

  !   all obs levels, compute interpolation weight factor:
  !                   linear interpolation in log(p) of u, v, T(v), RH
  DO ilev = 1 , nlev
    zlpf (ilev)  =   (col_lnp(kbz(ilev)) - zobbdy(ilev,nbtlop))                &
                   / (col_lnp(kbz(ilev)) - col_lnp(kbz(ilev)-1))
!                 mass-weighted linear interpolation (==> linear in p)
!   zlpf (ilev)  =   (col_p(kbz(ilev)) - zobbdy(ilev,nbtp))                    &
!                  / (col_p(kbz(ilev)) - col_p(kbz(ilev)-1))
  ENDDO

!-------------------------------------------------------------------------------
! End of module procedure prep_vi_mo2ob
!-------------------------------------------------------------------------------

END SUBROUTINE prep_vi_mo2ob



!-------------------------------------------------------------------------------
!+ Module procedure for observation operator for conventional single-level obs
!-------------------------------------------------------------------------------

SUBROUTINE sing_obs_operator ( ke, col_u, col_v, col_t, col_qv, col_qc, col_p  &
                             , zobbdy, lveridat, rdv, b1, b2w, b3, b4w         &
                             , zsobdy , vip, kb, zlopf )

!-------------------------------------------------------------------------------
!
! Description:
!   The routine applies the forward observation operator for (upper-air)
!   single-level reports    --------------------------------------------
!   -------------------- of horizontal wind, temperature and / or (specific)
!   humidity.
!   The forward operator consists here of a simple vertical interpolation
!   of model values to observation levels.
!   The interpolated values are also written to the Simulated Observation
!   Record SOR (for feedobs /feedback files).
!
!   Purpose: In COSMO, the forward observation operator is applied for
!            subsequent:
!              - filling NetCDF feedobs /feedback file for LETKF / verification
!              - quality control (when required)
!              - nudging scheme: spreading of obs. increments from obs. levels
!
! Method:
!   Vertical interpolation of model values to observation levels:
!     - interpolated variables :  wind components, (virtual) temperature,
!                                 relative humidity,  (height)
!     - interpolation method   :  linear in LOG( pressure )
!
!
! Written by        :  Christoph Schraff, DWD  (original version: 18.09.97)
! Current Code Owner:  Christoph Schraff, DWD
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: European Standards for Writing and Documenting
!                       Exchangeable Fortran 90 Code + DWD Extensions.
!
!-------------------------------------------------------------------------------
!
! Declarations:
!
!-------------------------------------------------------------------------------
! Modules used: These are declared in the module declaration section
! -------------

IMPLICIT NONE

! Subroutine arguments:
! --------------------

  INTEGER        , INTENT (IN)         ::       &
    ke                ! number of vertical (main) levels in model column

  LOGICAL        , INTENT (IN)         ::       &
    lveridat          ! fill simulated obs body (for writing to feedback file)

  REAL (KIND=wp) , INTENT (IN)         ::       &
                      ! model column (at the observation location) of:
    col_u   (ke)   ,& ! zonal wind speed             on Arakawa A grid [ m/s ]
    col_v   (ke)   ,& ! meridional wind speed        on Arakawa A grid [ m/s ]
    col_t   (ke)   ,& ! temperature                                    [  K  ]
    col_qv  (ke)   ,& ! specific water vapor content                   [kg/kg]
    col_qc  (ke)   ,& ! specific cloud water content                   [kg/kg]
!   col_rh  (ke)   ,& ! model relative humidity                        [     ]
    col_p   (ke)   ,& ! pressure (full value)        on main levels    [ Pa  ]
    rdv            ,& ! ratio of gas constant for dry air 'r_d'
                      !      and gas constant for water vapour 'r_v':  r_d/r_v
    b1             ,& ! \  constants for computing
    b2w            ,& !  \ the saturation vapour pressure
    b3             ,& !  / over water
    b4w               ! /  by the Magnus formula

  REAL (KIND=wp) , INTENT (IN)         ::       &
    zobbdy (mxsbdy)   ! single-level obs. body (format: see 'osgbdy')

  REAL (KIND=wp) , INTENT (INOUT)      ::       &
    zsobdy (mxsosg)   ! multi-level simulated obs. body (format: see SOR)

  REAL (KIND=wp) , INTENT (OUT)        ::       &
    vip    (4)     ,& ! model values interpolated to the obs level
    zlopf             ! weight factor for vertical interpol. to the obs level

  INTEGER        , INTENT (OUT)        ::       &
    kb                ! index of model level immediately below obs level
                      !   (if zero then obs is above top main model level)

! Local parameters: None
! ----------------

! Local scalars:
! -------------

  INTEGER        ::  &
    ka                ! upper model level adjacent to obs level

  REAL (KIND=wp) ::  &
    zrh   (2)         ! relat. humidity at model levels adjacent to obs level

!------------ End of header ----------------------------------------------------


!-------------------------------------------------------------------------------
! Begin Subroutine sing_obs_operator
!-------------------------------------------------------------------------------

  !   prepare vertical interpolation from model level to obs level
  !   by computing 'kb' and 'zlopf'

  CALL prep_vi_mo2ob_sg ( ke, col_p(:), zobbdy(nbsp) , kb, zlopf )
! =====================

  ! Linear vertical interpolation in log(p) of u, v, T, RH
  ! ------------------------------------------------------

  IF (kb > 0) THEN
    ka       =  kb - 1
!   zlopf    =  zlpf(ilev)
    zrh (1)  =    fq2pv( col_qv(ka) + col_qc(ka) , col_p(ka) , rdv )           &
                / fpvsw( col_t (ka), b1, b2w, b3, b4w )
    zrh (2)  =    fq2pv( col_qv(kb) + col_qc(kb) , col_p(kb) , rdv )           &
                / fpvsw( col_t (kb), b1, b2w, b3, b4w )
    vip (1)  =       (c1 -zlopf) * col_u (kb)   +   zlopf * col_u (ka)
    vip (2)  =       (c1 -zlopf) * col_v (kb)   +   zlopf * col_v (ka)
    vip (3)  =       (c1 -zlopf) * col_t (kb)   +   zlopf * col_t (ka)
    vip (4)  =  MIN( (c1 -zlopf) * zrh   (2)    +   zlopf * zrh   (1) , c1 )
  ENDIF

  ! write simulated observations to SOR 'zsobdy'
  ! --------------------------------------------

  IF ((lveridat) .AND. (kb > 0)) THEN
    IF (zobbdy(nbsu ) > rmdich)  zsobdy (nso_u )  =  vip(1)
    IF (zobbdy(nbsu ) > rmdich)  zsobdy (nso_v )  =  vip(2)
    IF (zobbdy(nbst ) > rmdich)  zsobdy (nso_t )  =  vip(3)
    IF (zobbdy(nbsrh) > rmdich)  zsobdy (nso_rh)  =  vip(4)
  ENDIF

!-------------------------------------------------------------------------------
! End of module procedure sing_obs_operator
!-------------------------------------------------------------------------------

END SUBROUTINE sing_obs_operator



!-------------------------------------------------------------------------------
!+ Module procedure for observation operator for conventional surface-level obs
!-------------------------------------------------------------------------------

SUBROUTINE surf_obs_operator ( ke, u_10m, v_10m, t_2m, td_2m, col_t, col_z     &
                             , ps, hsurf, zobbdy, zstalt, klu, kll, lveridat   &
                             , rdv, b1, b2w, b3, b4w                           &
                             , zsobdy , vip , loiqv2m_in , zpr )

!-------------------------------------------------------------------------------
!
! Description:
!   The routine applies the forward observation operator for surface-level
!   reports of horizontal   ----------------------------------------------
!   ------- wind, temperature and / or (specific) humidity.
!   The forward operator uses the 10-wind and 2-m quantities are computed by
!   the COSMO model parameterisation and applies a height correction for T-2m.
!   The simulated humidity observation is based on either the specific or the
!   relative humidity observation increment.
!   The interpolated values are also written to the Simulated Observation
!   Record SOR (for feedobs /feedback files).
!
!   Purpose: In COSMO, the forward observation operator is applied for
!            subsequent:
!              - filling NetCDF feedobs /feedback file for LETKF / verification
!              - quality control (when required)
!              - nudging scheme: spreading of obs. increments from obs. levels
!
! Method:
!   Vertical interpolation of model values to observation levels:
!     - interpolated variables :  wind components, (virtual) temperature,
!                                 relative humidity,  (height)
!     - interpolation method   :  linear in LOG( pressure )
!
!
! Written by        :  Christoph Schraff, DWD  (original version: 18.09.97)
! Current Code Owner:  Christoph Schraff, DWD
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: European Standards for Writing and Documenting
!                       Exchangeable Fortran 90 Code + DWD Extensions.
!
!-------------------------------------------------------------------------------
!
! Declarations:
!
!-------------------------------------------------------------------------------
! Modules used: These are declared in the module declaration section
! -------------

IMPLICIT NONE

! Subroutine arguments:
! --------------------

  INTEGER        , INTENT (IN)         ::       &
    ke             ,& ! number of vertical (main) levels in model column
    klu            ,& ! upper model level used for lapse rate
    kll               ! lower model level used for lapse rate

  LOGICAL        , INTENT (IN)         ::       &
    lveridat          ! fill simulated obs body (for writing to feedback file)

  REAL (KIND=wp) , INTENT (IN)         ::       &
                      ! model column (at the observation location) of:
    u_10m          ,& ! zonal wind speed in 10m      on Arakawa A grid ( m/s )
    v_10m          ,& ! meridional wind speed in 10m on Arakawa A grid ( m/s )
    t_2m           ,& ! temperature in 2m                              (  K  )
!   qv_2m          ,& ! specific water vapor content in 2m             (kg/kg)
    td_2m          ,& ! dew-point in 2m                                (  K  )
    col_t   (ke)   ,& ! temperature                                    [  K  ]
    col_z   (ke)   ,& ! geometrical height           of main levels    [  m  ]
    ps             ,& ! surface pressure (full value)                  [ Pa  ]
    hsurf          ,& ! model orography                                [  m  ]
    rdv            ,& ! ratio of gas constant for dry air 'r_d'
                      !      and gas constant for water vapour 'r_v':  r_d/r_v
    b1             ,& ! \  constants for computing
    b2w            ,& !  \ the saturation vapour pressure
    b3             ,& !  / over water
    b4w               ! /  by the Magnus formula

  REAL (KIND=wp) , INTENT (IN)         ::       &
    zstalt         ,& ! station altitude
    zobbdy (mxsbdy)   ! single-level obs. body (format: see 'osgbdy')

  REAL (KIND=wp) , INTENT (INOUT)      ::       &
    zsobdy (mxsosg)   ! single-level simulated obs. body (format: see SOR)

  REAL (KIND=wp) , INTENT (OUT)        ::       &
    vip    (4)        ! model values interpolated to the obs level

  LOGICAL        , INTENT (IN)    , OPTIONAL   ::       &
    loiqv2m_in        ! 2-m humidity obs increments based on differences of
                      !   specific humidity instead of relative humidity

  REAL (KIND=wp) , INTENT (OUT)   , OPTIONAL   ::       &
    zpr    (3)        ! auxilliary quantities for later diagnostic print-out

! Local parameters: None
! ----------------

! Local scalars:
! -------------

  REAL (KIND=wp) ::  &
    zlapse         ,& ! lapse rate used for the height correction of temperat.
    zpobz          ,& ! reference pressure assigned to observation level
    zqvvip         ,& ! model specific humidity extrapolated to obs. level
    zqvobs         ,& ! observed specific humidity at obs. level
!   viprh             ! relative humidity from model T-2m and Td-2m
    zff_ob         ,& !      observed          wind speed (not used further)
    zff_mo         ,& !      model (simulated) wind speed
    zdd_mo         ,& ! true model (simulated) wind direction
    zdd_mo_rot     ,& ! model (simulated) wind direction in rotated coordinates
    zdd_ob_rot     ,& ! observed          wind direction in rotated coordinates
    ewo            ,& ! vapour pressure
    eswo           ,& ! saturation vapour pressure over water
    ztd               ! dewpoint temperature

  LOGICAL        ::  &
    loiqv2m           ! 2-m humidity obs increments based on differences of
                      !   specific humidity instead of relative humidity

!------------ End of header ----------------------------------------------------


!-------------------------------------------------------------------------------
! Begin Subroutine surf_obs_operator
!-------------------------------------------------------------------------------

  vip (:) = rmdi

  ! set optional input parameter
  loiqv2m  = .FALSE.
  IF (PRESENT( loiqv2m_in ))  loiqv2m  =  loiqv2m_in

  ! initialise optional output parameter
  IF (PRESENT( zpr ))  zpr (:) = rmdi

  ! lapse rate for height correction for 2m temperature
  ! ---------------------------------------------------

! zlapse  =    (col_t (kml700) - col_t (kml850))                               &
!            / (col_z (kml700) - col_z (kml850))
  zlapse  =    (col_t (klu) - col_t (kll))                                     &
             / (col_z (klu) - col_z (kll))
  ! (alternative lapse rates : constant : -.0065 (e.g for height diff > 500m, or
  !  (Damrath): noon: -.0110 ; midnight : -.0045  --> with temporal interpolat.)
  !      ( for  Tmax: -.0130 ; for Tmin : -.0020 )

  ! Assigment of model values to observation level for u-10m, T-2m, RH-2m
  ! ---------------------------------------------------------------------

  vip (1)  =  u_10m
  vip (2)  =  v_10m
  !   height correction for 2m temperature
  vip (3)  =  t_2m  -  zlapse *(hsurf - zstalt)

  !   relative humidity from 2-m temperature and dewpoint temperature
  IF (zobbdy(nbsrh) > rmdich) THEN
    vip (4)  =    fpvsw( td_2m, b1, b2w, b3, b4w )                             &
                / fpvsw( t_2m , b1, b2w, b3, b4w )
    !   quality factor depending on temperature observation increment
    IF ((loiqv2m) .AND. (zobbdy(nbst) > rmdich)) THEN
      !  convert height corrected specific humidity increment into rel. humidity
      !   - for 'qv-vip': use model-RH-2m, height corrected model-T-2m, and
      !                       height corrected model-p-2m (set to model-ps)
      !   - for 'qv-obs': use observed RH-2m, T-2m, p-2m (set to model-ps)
      !     ==> pressure differences betw. obs and model (at station height) are
      !         not taken into account
      zpobz   =  ps
      zqvvip  =  fpv2q( vip(4)        *fpvsw( vip(3), b1, b2w, b3, b4w )       &
                      , zpobz , rdv )
      zqvobs  =  fpv2q( zobbdy(nbsrh) *fpvsw( zobbdy(nbst), b1, b2w, b3, b4w ) &
                      , zpobz , rdv )
      IF (PRESENT( zpr )) THEN
        zpr (1) = vip(4)
        zpr (2) = zqvvip
        zpr (3) = zqvobs
      ENDIF
      vip (4) =  zobbdy(nbsrh) +   SIGN ( c1 , zqvvip - zqvobs )               &
                                 * fq2pv( ABS( zqvvip - zqvobs ), zpobz, rdv ) &
                                 / fpvsw( zobbdy(nbst), b1, b2w, b3, b4w )
    ENDIF
  ENDIF

  ! !   enhance RH-2m obs error depending on temperature obs increment: factor
  ! IF ((lqfqv2m) .AND. (zobbdy(nbst) > rmdich)) THEN
  !   fqerr  =  exp( + (zobbdy(nbst) - vip(3))*(zobbdy(nbst) - vip(3))         &
  !                                  / (qqdts * qqdts) )


  ! compute derived simulated observations for writing to feedback file
  ! -------------------------------------------------------------------

  IF (lveridat) THEN
    zff_mo  =  rmdi
    zdd_mo  =  rmdi
    ztd     =  rmdi

    ! 10-m wind speed and direction:
    !    1.: compute model wind speed (possible from rotated wind components)
    !    2.: to avoid the need to rotate wind coordinates (requires variables
    !        such as 'rlat', 'pollat', etc.), compute:
    !        true model wind direction = true observed wind direction +
    !         (model dir. in rotated coord) - (observed dir. in rotated coord)
    IF ((zobbdy(nbsff) > rmdich) .AND. (zobbdy(nbsdd) > rmdich)) THEN

      CALL uv2df ( vip(1)      , vip(2)      , zdd_mo_rot, zff_mo )
    ! ==========
      CALL uv2df ( zobbdy(nbsu), zobbdy(nbsv), zdd_ob_rot, zff_ob )
    ! ==========
      zdd_mo  =  zobbdy(nbsdd)  +  (zdd_mo_rot - zdd_ob_rot)
    ENDIF

    ! 2-m dewpoint temperature:
    !   (re-)compute Td-2m from (height corrected) T-2m and derived RH-2m
    !   in order to attain consistency between simulated Td-2m and RH-2m
    !   (computation analogous to 'obs_rh2td' in src_obs_cdfin_util.f90)
    IF ((zobbdy(nbstd) > rmdich) .AND. (zobbdy(nbst) > rmdich)) THEN
      eswo  =  fpvsw( vip(3), b1, b2w, b3, b4w )
      IF (zobbdy(nbsrh) > rmdich) THEN
        ewo   =  eswo * MAX( MIN( vip(4) , c1 ) , 0.0001_wp )
      ELSE
        ! if observ. TD2M is given and RH2M is not, then use RH2M
        !   from original T2M and TD2M model values (irrespective of 'loiqv2m')
        !   Note: should appear only when processing some old fdbk-files
        !         (and does not occur when reading 'cdfin' files)
        ewo   =  eswo * MAX( MIN(   fpvsw( td_2m, b1, b2w, b3, b4w )           &
                                  / fpvsw( t_2m , b1, b2w, b3, b4w )           &
                                , c1 ) , 0.0001_wp )
      END IF
      !   (td, as e.g. stored in data files, is over water (WMO convention),
      !    thus Magnus formula over water instead of ice is used)
      IF (LOG(ewo/b1) < b2w-epsy)  ztd  =  ftd( ewo, b1, b2w, b3, b4w )
    ENDIF

  ! write simulated observations to SOR 'zsobdy'
  ! --------------------------------------------

    IF (zobbdy(nbsu ) > rmdich)  zsobdy (nso_u )  =  vip(1)
    IF (zobbdy(nbsu ) > rmdich)  zsobdy (nso_v )  =  vip(2)
    IF (zobbdy(nbst ) > rmdich)  zsobdy (nso_t )  =  vip(3)
    IF (zobbdy(nbsrh) > rmdich)  zsobdy (nso_rh)  =  vip(4)
    IF (zobbdy(nbsff) > rmdich)  zsobdy (nso_ff)  =  zff_mo
    IF (zobbdy(nbsdd) > rmdich)  zsobdy (nso_dd)  =  zdd_mo
    IF (zobbdy(nbstd) > rmdich)  zsobdy (nso_td)  =  ztd
  ENDIF

!-------------------------------------------------------------------------------
! End of module procedure surf_obs_operator
!-------------------------------------------------------------------------------

END SUBROUTINE surf_obs_operator



!-------------------------------------------------------------------------------
!+ Module procedure for observation operator for conventional cloud obs
!-------------------------------------------------------------------------------

SUBROUTINE surf_obs_opr_1h ( precip_1h, raddir, raddif, radlw, zobbdy , zsobdy )

!-------------------------------------------------------------------------------
!
! Description:
!   The routine applies the forward observation operator for some temporally
!   non-local observations, namely 1-hourly sums of precipitation and of
!   radiation.
!   The simulated values are written to the Simulated Observation Record SOR
!   (for feedobs /feedback files).
!
! Method:
!   Simple assignment, 1-hour sums or averaged (power) values have to be input.
!
! Written by        :  Christoph Schraff, DWD
! Current Code Owner:  Christoph Schraff, DWD
!-------------------------------------------------------------------------------
! Declarations:
!-------------------------------------------------------------------------------
! Modules used: These are declared in the module declaration section
! -------------

IMPLICIT NONE

! Subroutine arguments:
! --------------------

  REAL (KIND=wp) , INTENT (IN)         ::       &
    precip_1h      ,& ! 1-hr precipitation sum                           [mm]
    raddir         ,& ! direct (downward) solar radiation flux, 1-h mean [W/m2]
    raddif         ,& ! diffuse downward  solar radiation flux, 1-h mean [W/m2]
    radlw             ! long-wave downward      radiation flux, 1-h mean [W/m2]

  REAL (KIND=wp) , INTENT (IN)         ::       &
    zobbdy (mxsbdy)   ! single-level obs. body (format: see 'osgbdy')

  REAL (KIND=wp) , INTENT (INOUT)      ::       &
    zsobdy (mxsosg)   ! multi-level simulated obs. body (format: see SOR)

! Local parameters: None
! Local scalars: None

!------------ End of header ----------------------------------------------------

!-------------------------------------------------------------------------------
! Begin Subroutine surf_obs_opr_1h
!-------------------------------------------------------------------------------

  IF (zobbdy(nbsrr1) > rmdich)  zsobdy (nso_rr1) =  precip_1h
  IF (zobbdy(nbsrad) > rmdich)  zsobdy (nso_rad) = (raddir + raddif)*3600._wp
  IF (zobbdy(nbsrdd) > rmdich)  zsobdy (nso_rdd) =  raddif          *3600._wp
  IF (zobbdy(nbsrdl) > rmdich)  zsobdy (nso_rdl) =  radlw           *3600._wp

!-------------------------------------------------------------------------------
! End of module procedure surf_obs_opr_1h
!-------------------------------------------------------------------------------

END SUBROUTINE surf_obs_opr_1h



!-------------------------------------------------------------------------------
!+ Module procedure for observation operator for conventional cloud obs
!-------------------------------------------------------------------------------

SUBROUTINE cloud_obs_operator ( ke, clc, clct, clcl, clcm, clch, col_hhl       &
                              , zobbdy , zsobdy , ceil , vis )

!-------------------------------------------------------------------------------
!
! Description:
!   The routine applies the forward observation operator for cloud observations.
!                           ---------------------------------------------------
!   The forward operator uses the cloud fraction at model levels as well as
!   the integrated model cloud quantities 'total cloud', 'low cloud',
!   'midlevel cloud', and 'high cloud' and consist basically of a simple
!   assignment.
!   The simulated values are written to the Simulated Observation Record SOR
!   (for feedobs /feedback files).
!
! Method:
!   Simple assignment for cloud fraction; cloud base height as height
!   difference between orography and model layer with cloud fraction > 1% .
!
!
! Written by        :  Christoph Schraff, DWD
! Current Code Owner:  Christoph Schraff, DWD
!-------------------------------------------------------------------------------
! Declarations:
!-------------------------------------------------------------------------------
! Modules used: These are declared in the module declaration section
! -------------

IMPLICIT NONE

! Subroutine arguments:
! --------------------

  INTEGER        , INTENT (IN)         ::       &
    ke                ! number of vertical (main) levels in model column

! LOGICAL                 , INTENT (IN)         ::       &
!   lveridat          ! fill simulated obs body (for writing to feedback file)

  REAL (KIND=wp) , INTENT (IN)         ::       &
                      ! model column (at the observation location) of:
    clc     (ke)   ,& ! cloud fraction                                 [     ]
    col_hhl (ke+1) ,& ! geometrical height           of half levels    [  m  ]
    clct           ,& ! total cloud cover                              [     ]
    clcl           ,& ! low   cloud cover                              [     ]
    clcm           ,& ! mid-level cloud cover                          [     ]
    clch              ! high  cloud cover                              [     ]
!   td2m           ,& ! 2m temperature for CBH (laps rate method)      [  K  ]
!   t2m            ,& ! 2m dp temperature for CBH (laps rate method)   [  K  ]
!   col_p   (ke)   ,& ! pressure (full value)        on main levels    [ Pa  ]
!   col_t   (ke)   ,& ! temperature                                    [  K  ]
!   col_qv  (ke)   ,& ! specific water vapor content                   [kg/kg]
!   r_d            ,& ! gas constant for dry air       (287.05)
!   r_v            ,& ! gas constant for water vapour  (461.51 J /(kg * K))
!   t0_melt        ,& ! melting temperature of ice
!   b1             ,& ! \  constants for computing
!   b2w            ,& !  \ the saturation vapour pressure
!   b3             ,& !  / over water
!   b4w               ! /  by the Magnus formula

  REAL (KIND=wp) , INTENT (IN)    , OPTIONAL :: &
    ceil           ,& ! cloud ceiling (as computed by model diagnostic [  m  ]
    vis               ! visibility    (as computed by model diagnostic [  m  ]

  REAL (KIND=wp) , INTENT (IN)         ::       &
    zobbdy (mxsbdy)   ! single-level obs. body (format: see 'osgbdy')

  REAL (KIND=wp) , INTENT (INOUT)      ::       &
    zsobdy (mxsosg)   ! multi-level simulated obs. body (format: see SOR)

! Local parameters: None
! ----------------

! Local scalars:
! -------------

  LOGICAL        ::  &
    lceil             ! ceiling has already been computed

  INTEGER        ::  &
    kcl            ,& ! layer with cloud fraction > 1 %
    kk                ! loop incdex

! REAL (KIND=wp) ::  &
!   rdv               ! ratio of gas constant for dry air 'r_d'
                      !      and gas constant for water vapour 'r_v':  r_d/r_v

!------------ End of header ----------------------------------------------------

!-------------------------------------------------------------------------------
! Begin Subroutine cloud_obs_operator
!-------------------------------------------------------------------------------

  ! cloud amount
  IF (zobbdy(nbsct ) > rmdich)  zsobdy (nso_ct )  =  clct
  IF (zobbdy(nbscl ) > rmdich)  zsobdy (nso_cl )  =  clcl
  IF (zobbdy(nbscm ) > rmdich)  zsobdy (nso_cm )  =  clcm
  IF (zobbdy(nbsch ) > rmdich)  zsobdy (nso_ch )  =  clch

  ! visibility   (as diagnosed from model physics, if available)
  IF ((zobbdy(nbsvis) > rmdich) .AND. (PRESENT( vis ))) THEN
    IF (vis >= 0.0_wp)  zsobdy (nso_vis)  =  vis
  ENDIF

  ! ceiling
  lceil  = .FALSE.
  IF (zobbdy(nbscil) > rmdich) THEN
    IF (PRESENT( ceil )) THEN
      ! use ceiling as diagnosed from the model physics
      ! if it is available / has been read
      IF (ceil >= 0.0_wp) THEN
        ! 'ceil'  is defined ASL (above  sea   level) in ICON
        ! the obs is defined AGL (above ground level) in Synop
        zsobdy (nso_cil)  =  ceil - col_hhl(ke+1)
        lceil  = .TRUE.
      ENDIF
    ENDIF
    IF ((.NOT. lceil) .AND. (MAXVAL( clc ) >= 0.0_wp)) THEN
      ! use cloud cover profile: lowest level with cloud cover > 4 octas
      kcl = 0
      DO kk = ke , 1 , -1
        IF ((clc(kk) > 0.5625_wp) .AND. (kcl == 0))  kcl = kk
      ENDDO
      IF (kcl > 0)  zsobdy (nso_cil)  =  col_hhl(kcl+1) - col_hhl(ke+1)
    ENDIF
  ENDIF

  ! cloud base height
  IF (zobbdy(nbscbs) > rmdich) THEN
    IF (MAXVAL( clc ) >= 0.0_wp) THEN
      ! use cloud cover profile: lowest level with cloud cover > 4 octas
      kcl = 0
      DO kk = ke , 1 , -1
        IF ((clc(kk) > 0.01_wp  ) .AND. (kcl == 0))  kcl = kk
      ENDDO
      IF (kcl > 0)  zsobdy (nso_cbs)  =  col_hhl(kcl+1) - col_hhl(ke+1)
    ENDIF

    !  alternative estimations of cloud base height from model profile (Declair)
    !  - via lifted condensation level, epsy formula
!!  zsobdy (nso_cbs) = (t2m - td2m)* 400._wp *0.3048_wp
    !  - via lawrence formula
!!  rh2m  =  fpvsw( td2m ,b1, b2w, b3, b4w ) / fpvsw( t2m , b1, b2w, b3, b4w )
!!  zsobdy (nso_cbs) = (10._wp + 0.2_wp*(t2m - 273.15_wp)) * (100._wp - rh2m)
    !  - via relative humidity (water vapor pressure / saturation pressure)
!   kcl = 0
!   rdv =  r_d / r_v
!   DO kk = ke , 1 , -1
!     zrh  =  fq2pv( col_qv(kk) , col_p(kk) , rdv )                            &
!           / fpvsw( col_t(kk) , b1, b2w, b3, b4w )
!     IF ((zrh > 0.95_wp) .AND. (kcl == 0))  kcl = kk
!   ENDDO
!   IF (kcl > 0)  zsobdy (nso_cbs)  =  col_hhl(kcl+1) - col_hhl(ke+1)
  ENDIF

!-------------------------------------------------------------------------------
! End of module procedure cloud_obs_operator
!-------------------------------------------------------------------------------

END SUBROUTINE cloud_obs_operator


!-------------------------------------------------------------------------------
!+ Module procedure for conversion of fraction into octas
!-------------------------------------------------------------------------------

FUNCTION frac2octas ( frac )

!-------------------------------------------------------------------------------
!
! Description:
!   The routine convert an input quantity (typically cloud fraction) into octas
!   according to WMO descriptor table 020011:
!     0 octas   ---  0                           -->  0    - 0.01
!     1 octa    ---  1/10 or less, but not zero  -->  0.01 - 0.15
!     2 octas   ---  2/10 - 3/10                 -->  0.15 - 0.35
!     3 octas   ---  4/10                        -->  0.35 - 0.45
!     4 octas   ---  5/10                        -->  0.45 - 0.55
!     5 octas   ---  6/10                        -->  0.55 - 0.65
!     6 octas   ---  7/10 - 8/10                 -->  0.65 - 0.85
!     7 octa    ---  9/10 or more, but not 10/10 -->  0.85 - 0.99
!     8 octas   ---  10/10                       -->  0.99 - 1
!     other code figures for 'sky obscured' (e.g. by fog), 'scattered',
!     'broken', or 'few' are not set.
!
! Method:
!   Straightforward.
!
! Written by        :  Christoph Schraff, DWD
! Current Code Owner:  Christoph Schraff, DWD
!-------------------------------------------------------------------------------
! Declarations:
!-------------------------------------------------------------------------------
! Modules used: These are declared in the module declaration section
! -------------

IMPLICIT NONE

! Subroutine arguments:
! --------------------

  REAL (KIND=wp) , INTENT (IN)         ::       &
    frac          ! quantity as fraction

! Local parameters: None
! ----------------

! Local scalars:
! -------------

  REAL (KIND=wp) ::  &
    frac2octas    ! quantity as octa

  INTEGER        ::  &
    iocta         ! quantity as octa

!------------ End of header ----------------------------------------------------

!-------------------------------------------------------------------------------
! Begin Subroutine frac2octas
!-------------------------------------------------------------------------------

  IF (frac < 0.01_wp) THEN
    iocta = 0
  ELSEIF (frac > 0.99_wp) THEN
    iocta = 8
  ELSEIF (frac < 0.15_wp) THEN
    iocta = 1
  ELSEIF (frac > 0.85_wp) THEN
    iocta = 7
  ELSE
    iocta  =  MIN( 6 , MAX( 2 , NINT( frac *10._wp ) - 1 ) )
  ENDIF

  frac2octas  =  REAL( iocta , wp )

!-------------------------------------------------------------------------------
! End of module procedure frac2octas
!-------------------------------------------------------------------------------

END FUNCTION frac2octas



!-------------------------------------------------------------------------------
!+ Module procedure for observation operator for GPS ZPD (and IWV) obs
!-------------------------------------------------------------------------------

SUBROUTINE zpd_iwv_obs_operator ( ke, col_p, col_hhl, col_t, col_qv, col_qc    &
                                ,     col_qrs, r_g, r_d, r_v, t0_melt, degrad  &
                                ,     b1, b2w, b3, b4w, b2i, b4i, madj_hum     &
                                , zobzpd, zobiwv, zstalt, zstlat, lveridat     &
                                , zsobdy                                       &
                                , oiqint, oiqobs, oiqmod, oiqmic, oiqsat       &
                                , otqsat, oqvsat, orhmod, orhodz, mdabot       &
                                , ofracb, opint )

!-------------------------------------------------------------------------------
!
! Description:
!   The routine applies the forward observation operator for ground-based GPS
!   zenith path delay (ZPD) reports -----------------------------------------
!   ------------------------------- .
!   In addition,    the integrated water vapour (IWV) is computed both for the
!   for the observation ----------------------------- (from the observed ZPD)
!   ------------------- and for the simulated observation (from the model
!   humidity profile).  --------------------------------- This is a mixture of
!   forward and inverse observation operator.
!   Furthermore, model profiles and vertically integrated model quantities,
!   that are required for the quality control and/or for the nudging, are also
!   computed, as well as the index and fraction of the model level used to
!   compute the simulated IWV above the GPS station height. Note that all the
!   vertically integrated quantities refer to GPS station (antenna) height,
!   except for 'otqsat' (refers to model orography).
!
!   Purpose: In COSMO, the computed quantities are applied for subsequent:
!              - filling NetCDF feedobs /feedback file for LETKF / verification
!              - quality control (when required)
!              - nudging scheme: spreading of obs. increments from obs. levels
!
! Method:
!   By applying the algorithm of Bevis et al. (1994) and its inverse, using
!   station latitude and height together with model temperature and pressure
!   at station height as input.
!
! Written by        :  Christoph Schraff, DWD, 02.04.13
! Current Code Owner:  Christoph Schraff, DWD
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: European Standards for Writing and Documenting
!                       Exchangeable Fortran 90 Code + DWD Extensions.
!
!-------------------------------------------------------------------------------
!
! Declarations:
!
!-------------------------------------------------------------------------------
! Modules used: These are declared in the module declaration section
! -------------

IMPLICIT NONE

! Subroutine arguments:
! --------------------

  INTEGER        , INTENT (IN)         ::       &
    ke         ,& ! number of vertical (main) levels in model column
    madj_hum      ! if = 1 : adjust observed humidity (by ratio of saturation
                  !          vapour pressure over water to the one over ice,
                  !          to be applied if cloud ice is not a state variable)

  LOGICAL        , INTENT (IN)         ::       &
    lveridat          ! fill simulated obs body (for writing to feedback file)

  REAL (KIND=wp) , INTENT (IN)         ::       &
                      ! model column (at the observation location) of:
    col_p   (ke)   ,& ! pressure (full value)        on main levels    [ Pa  ]
    col_hhl (ke+1) ,& ! geometrical height           of half levels    [  m  ]
    col_t   (ke)   ,& ! temperature                                    [  K  ]
    col_qv  (ke)   ,& ! specific water vapor content                   [kg/kg]
    col_qc  (ke)   ,& ! specific cloud water content                   [kg/kg]
    col_qrs (ke)   ,& ! spec. cont. of hydrometeors excl. cloud water  [kg/kg]
!   col_rh  (ke)   ,& ! model relative humidity                        [     ]
    r_g            ,& ! acceleration due to gravity    (  9.80665)
    r_d            ,& ! gas constant for dry air       (287.05)
    r_v            ,& ! gas constant for water vapour  (461.51 J /(kg * K))
    t0_melt        ,& ! melting temperature of ice
    degrad         ,& ! factor for transforming degree to rad
    b1             ,& ! \  constants for computing
    b2w  , b2i     ,& !  \ the saturation vapour pressure
    b3             ,& !  / over water
    b4w  , b4i        ! /  by the Magnus formula

  REAL (KIND=wp) , INTENT (IN)         ::       &
    zstlat         ,& ! station latitude                            [deg]
    zstalt         ,& ! station altitude                            [ m ]
    zobzpd         ,& ! observed zenith total path delay ZPD = ZTD  [mm ]
    zobiwv            ! observed integrated water vapour     [kg/m2 ~ mm]
!   zobbdy (mxgbdy)   ! GPS IWV obs. body (format: see 'ogpbdy')

  REAL (KIND=wp) , INTENT (INOUT)      ::       &
    zsobdy (mxsosg)   ! multi-level simulated obs. body (format: see SOR)

  INTEGER        , INTENT (OUT)   , OPTIONAL     ::       &
    mdabot            ! index of lowest model level with GPS pseudo obs.

  REAL (KIND=wp) , INTENT (OUT)   , OPTIONAL     ::       &
    oiqint         ,& ! IWV derived from observed ZPD
    oiqobs         ,& ! obs.  integrated mass of water vapour q: obs.  IWV
    oiqmod         ,& ! model integrated mass of water vapour q: model IWV
    oiqmic         ,& ! IWV from water-to-ice corrected model IWV 'ziqmod'
    oiqsat         ,& ! IWV from model at saturation q-sat (from GPS height)
    otqsat         ,& ! as 'ziqsat', but from model orography
    oqvsat   (ke)  ,& ! model    specific humidity at saturation
    orhmod   (ke)  ,& ! model    relative humidity  (not generalized)
    orhodz   (ke)  ,& ! density of moisty air * thickness of model layer
    ofracb         ,& ! weight for contribution to model IWV from lowest layer
    opint             ! model pressure interpolated at GPS height [hPa]

! Local parameters: None
! ----------------

! Local scalars:
! -------------

  INTEGER        ::  &
    km             ,& ! vertical (loop) index of model level (== obs)
    kdabot         ,& ! index of lowest model level with GPS pseudo obs.
    kdatop            ! index of uppermost model level with GPS pseudo obs

  REAL (KIND=wp) ::  &
    ziqint         ,& ! IWV derived from observed ZPD
    ziqobs         ,& ! obs.  integrated mass of water vapour q: obs.  IWV
    ziqmod         ,& ! model integrated mass of water vapour q: model IWV
    ziqmic         ,& ! IWV from water-to-ice corrected model IWV 'ziqmod'
    ziqsat         ,& ! IWV from model at saturation q-sat (from GPS height)
    ztqsat         ,& ! as 'ziqsat', but from model orography
    zqvsat   (ke)  ,& ! model    specific humidity at saturation
    zrhmod   (ke)  ,& ! model    relative humidity  (not generalized)
    zrhodz   (ke)  ,& ! density of moisty air * thickness of model layer
    zfracb         ,& ! weight for contribution to model IWV from lowest layer
    zpdmod         ,& ! model derived zenith total delay
    zrho           ,& ! density of moisty air
    zpv            ,& ! water vapour pressure
    zpvsat         ,& ! saturation water vapour pressure
    pint           ,& ! model pressure interpolated at GPS height [hPa]
    tint           ,& ! model temperature interpolated at GPS height [K]
    hhke           ,& ! height of lowest full model level
    zzo            ,& ! height of upper model level
    zzu            ,& ! height of lower model level
    zfrac          ,& ! weight
    rdv            ,& ! ratio of gas constant for dry air 'r_d'
                      !      and gas constant for water vapour 'r_v':  r_d/r_v
    rvd_m_o           ! r_v / r_d  -  1.0

!------------ End of header ----------------------------------------------------


!-------------------------------------------------------------------------------
! Begin Subroutine zpd_iwv_obs_operator
!-------------------------------------------------------------------------------

  rdv      =  r_d / r_v
  rvd_m_o  =  c1 / rdv  -  c1

!-------------------------------------------------------------------------------
!  Section 1: Computation of the required model profiles
!             and vertically integrated model quantities,
!             including the simulated (obs-compatible) IWV 'ziqmic'
!-------------------------------------------------------------------------------

  ! determine (fraction of) model levels used to compute the model IWV
  ! (from the station height to the model top)
  ! ------------------------------------------------------------------

  !   if GPS antenna below surface take surface as lowest level
  kdabot = ke
  kdatop = 1

  !  determine fraction 'zfracb' of model layer that includes GPS antenna height
  !  Note: for GPS station below orography: zfracb > c1 ,i.e. constant extrapol.
  lowest_contribution: DO km = ke , 1 , -1
     kdabot = km
     zfracb =   ( col_hhl(km) - zstalt )                                       &
              / ( col_hhl(km) - col_hhl(km+1) )
     IF (zstalt-col_hhl(km) < c0)                       EXIT lowest_contribution
  ENDDO lowest_contribution

  ! compute required model profiles
  ! -------------------------------

  DO km = 1 , ke
    !--   Model Saturation specific humidity profile
    zpvsat       =  fpvsw( col_t(km) , b1, b2w, b3, b4w )
    zqvsat (km)  =  fpv2q( zpvsat , col_p(km) , rdv )
    !--   Model relative humidity profile
    zpv          =  fq2pv( col_qv(km) , col_p(km) , rdv )
    zrhmod (km)  =  zpv / zpvsat
    !--   Density of moist air
    zrho         =  col_p(km) /( r_d * col_t(km) *(  c1 + rvd_m_o*col_qv(km)   &
                                                   - col_qc(km) - col_qrs(km)) )
    !--   Mass per m^2 of moist air within the complete model layer
    zrhodz (km)  =  zrho *(col_hhl(km) - col_hhl(km+1))
  ENDDO

  ! compute required vertically integrated model quantities
  ! -------------------------------------------------------
  !   (and adjust 'zrhodz')

  ziqmod = c0
  ziqsat = c0
  ztqsat = c0

  DO km = ke , 1 , -1
    zfrac = c1
    IF ( km == kdabot)  zfrac = zfracb
    IF ((km >  kdabot) .OR. (km < kdatop))  zfrac = c0
    !--   Integrated mass of humidity in total model column assumed saturated
    ztqsat       =  ztqsat  +  zrhodz(km) * zqvsat(km)
    !--   Mass per m^2 of moist air within the (fraction of the) model layer
    zrhodz (km)  =  zrhodz(km) * zfrac
    !--   Integrated mass of model specific humidity (bottom layer with zfracb)
    !--   Integrated mass of model saturation specific humidity
    ziqmod       =  ziqmod  +  zrhodz(km) * col_qv(km)
    ziqsat       =  ziqsat  +  zrhodz(km) * zqvsat(km)
  ENDDO

  !   for model version without prognostic cloud ice:
  !   integrated mass of water-to-ice corrected model specific humidity
  ziqmic = ziqmod
  IF (madj_hum == 1) THEN
    ziqmic = c0
    DO km = kdabot , kdatop , -1
      IF (col_t(km) < t0_melt) THEN
        zpv    = fpvsw( col_t(km), b1, b2i, b3, b4i ) * zrhmod(km)
        ziqmic = ziqmic + zrhodz(km) * fpv2q( zpv , col_p(km) , rdv )
      ELSE
        ziqmic = ziqmic + zrhodz(km) * col_qv(km)
      ENDIF
    ENDDO
  ENDIF

!-------------------------------------------------------------------------------
!  Section 2:  Determination of the 'observed' IWV by computing it from
!              observed ZPD and interpolated pressure and temperature
!-------------------------------------------------------------------------------

  ! get model pressure and temperature interpolated to station height
  ! -----------------------------------------------------------------

  hhke  =  c05* (col_hhl(ke) + col_hhl(ke+1))

  IF (zstalt >= hhke) THEN
    pint = col_p(1)
    tint = col_t(1)
    zzo  = hhke
    DO  km = ke, 2, -1
      zzu = zzo
      zzo = c05* (col_hhl(km) + col_hhl(km-1))
      IF ((zzu <= zstalt) .AND. (zzo > zstalt)) THEN
        !   interpolate ln(p) linear in z
        zfrac  = ( zstalt - zzu ) / ( zzo - zzu )
        pint   =  EXP( (c1 - zfrac) * LOG( col_p(km  ) )                       &
                           + zfrac  * LOG( col_p(km-1) ) )
        !   interpolate T     linear in z
        tint   =       (c1 - zfrac) * col_t(km)                                &
                           + zfrac  * col_t(km-1)
    EXIT
      ENDIF
    ENDDO
  ELSE
    !   if station is below lowest model level (or even below model orography),
    !   then pressure is extrapolated from lowest model level, and temperature
    !   is used from lowest model level
    !   (note, temperature errors lead to small errors in derived IWV)
    pint  =  EXP( LOG( col_p(ke) ) + r_g/(r_d*col_t(ke)) *(hhke-zstalt) )
!   pint  =  col_p(ke) * EXP( r_g /(r_d*col_t(ke)) *(hhke-zstalt) )
    tint  =  col_t(ke)
  ENDIF

  pint = pint / 100_wp

  IF (zobzpd > rmdich) THEN

  ! convert ZPD (ZTD) into IWV
  ! --------------------------

    ziqint  =  conv_zpd_iwv ( zstlat, zstalt, pint, tint, zobzpd, r_v, degrad )
!              ============

  ! assign the observed IWV: prefer the converted ZPD over the reported IWV
  ! -----------------------------------------------------------------------

    ziqobs  =  ziqint
  ELSE
    ziqobs  =  zobiwv
    !   only for output
    ziqint  = -ziqobs
!   IF (zobbdy(nbgiwv) > rmdich)  ziqobs = zobbdy(nbgiwv)
  ENDIF

  ! adjust the observed IWV to make it compatible to the model IWV
  ! if the model does not know prognostic cloud ice
  ! --------------------------------------------------------------

  IF (madj_hum == 1)  ziqobs  =  ziqobs * ziqmod / ziqmic

  ! apply bias correction to GPS Integrated water vapour
  ! ----------------------------------------------------

! IF (lgpsbias)       ziqobs  =  ziqobs + ogpbdy(ngpob,nbgbia)

  ! set negative values to zero  (values < 2 mm will be not assimilated anyway!)
  ! ---------------------------

  IF (ziqobs <= c0)   ziqobs  =  0.001_wp

!-------------------------------------------------------------------------------
!  Section 3:  Write simulated observations to SOR 'zsobdy'.
!              - IWV: also written, if ODR does not contain an obs value;
!                     in COSMO/nudging, the obs value 'ziqobs' is written to
!                     the ODR after calling this routine
!              - ZPD: derived from model IWV 'ziqmic'
!                (if the model does not know progrostic cloud ice, model IWV
!                 'ziqmic' has been adjusted to make it observation compatible;
!                 but bias correction has not and should not be applied to the
!                 model value)
!              Also write optional variables, if required
!-------------------------------------------------------------------------------

  IF (lveridat) THEN

    zpdmod  =  conv_iwv_2_zpd ( zstlat, zstalt, pint, tint, ziqmic, r_v, degrad)
            !  ==============

                          zsobdy (nso_iq )  =  ziqmod
    IF (zobzpd > rmdich)  zsobdy (nso_zpd)  =  zpdmod
  ENDIF

  ! write optional variables, if present
  ! ------------------------------------

  IF (PRESENT( oiqint ))  oiqint  =  ziqint
  IF (PRESENT( oiqobs ))  oiqobs  =  ziqobs
  IF (PRESENT( oiqmod ))  oiqmod  =  ziqmod
  IF (PRESENT( oiqmic ))  oiqmic  =  ziqmic
  IF (PRESENT( oiqsat ))  oiqsat  =  ziqsat
  IF (PRESENT( otqsat ))  otqsat  =  ztqsat
  IF (PRESENT( oqvsat ))  oqvsat  =  zqvsat
  IF (PRESENT( orhmod ))  orhmod  =  zrhmod
  IF (PRESENT( orhodz ))  orhodz  =  zrhodz
  IF (PRESENT( mdabot ))  mdabot  =  kdabot
  IF (PRESENT( ofracb ))  ofracb  =  zfracb
  IF (PRESENT( opint  ))  opint   =  pint

!-------------------------------------------------------------------------------
! End of module procedure zpd_iwv_obs_operator
!-------------------------------------------------------------------------------

END SUBROUTINE zpd_iwv_obs_operator



!-------------------------------------------------------------------------------
!+ Module procedure for conversion of ZPD (ZTD) into IWV
!-------------------------------------------------------------------------------

FUNCTION conv_zpd_iwv ( zlat , zalt , press , ztt , zpd , r_v , degrad )

!-------------------------------------------------------------------------------
!
! Description :  Compute IWV from zenith total path delay (ZPD),
!                using station latitude and height together with (model)
!                pressure and temperature at station height in the algorithm of
!                Bevis et al. (1994).
!
! Written by        :  M. Tomassini, DWD, 09.01.2001
! Current Code Owner:  C. Schraff
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: European Standards for Writing and Documenting
!                       Exchangeable Fortran 90 Code + DWD Extensions.
!
!-------------------------------------------------------------------------------
!
! Declarations:
!
!-------------------------------------------------------------------------------

IMPLICIT NONE

! Function arguments:
! ------------------

  REAL (KIND=wp) , INTENT (IN)  ::      &
    zlat, zalt, press    ,& ! latitude [degree], altitude [m], pressure [hPa]
    ztt , zpd            ,& ! temperature [K], zenith path delay [mm]
    r_v                  ,& ! gas constant for water vapour (461.51 J /(kg * K))
    degrad                  ! factor for transforming degree to rad

! Local variables:
! ---------------

  REAL (KIND=wp) , PARAMETER  :: &
!   r_v  =  4.619E-3_wp ,& ! hPa * m^3 /(g * K)     =  461.51 J /(kg * K)
!   rhow =  1E6_wp      ,& ! water density  g/m^3
    ck2s = 22.1_wp      ,& ! +/- 2.2    K/hPa   (BEVIS,1994)
    ck3  =  3.739E5_wp     ! +/- 0.012  K^2/hPa

  REAL (KIND=wp) ::      &
    conv_zpd_iwv            ,& ! return value of function
    zwd, zhd, tm, pfak, flam   ! to derive IWV from ZPD (routine from Gendt,GFZ)

!------------ End of header ----------------------------------------------------

!-------------------------------------------------------------------------------
! Begin Function conv_zpd_iwv
!-------------------------------------------------------------------------------

! Zenith Hydrostatic Delay  (Runge et.al. )
  flam = c1 - 0.00266_wp *COS( c2 *zlat *degrad )                          &
            - 0.28E-6_wp *zalt
  zhd  = 2.2765_wp *press /flam

! Zenith Wet Delay
  zwd  = zpd - zhd

! Zenith Wet Delay conversion into IWV
  tm   = 70.2_wp + 0.72_wp *ztt
  pfak = 1E5_wp / (r_v *(ck2s + ck3/tm))

  conv_zpd_iwv = pfak * zwd

!-------------------------------------------------------------------------------
! End of function conv_zpd_iwv
!-------------------------------------------------------------------------------

END FUNCTION conv_zpd_iwv



!-------------------------------------------------------------------------------
!+ Module procedure to compute ZPD from IWV
!-------------------------------------------------------------------------------

FUNCTION conv_iwv_2_zpd ( zlat , zalt , press , ztt , ziwv , r_v , degrad )

!-------------------------------------------------------------------------------
!
! Description :  Compute zenith total path delay (ZPD) from integrated water
!                vapour (IWV) by inverting the algorithm of Bevis et al. (1994).
!                Input: station latitude and height, and
!                       (model) pressure and temperature at station height.
!                Output: ZPD [mm]
!
! Written by        :  C. Schraff, DWD, 27.03.2013
! Current Code Owner:  C. Schraff, DWD
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: European Standards for Writing and Documenting
!                       Exchangeable Fortran 90 Code + DWD Extensions.
!
!-------------------------------------------------------------------------------
!
! Declarations:
!
!-------------------------------------------------------------------------------

IMPLICIT NONE

! Function arguments:
! ------------------

  REAL (KIND=wp) , INTENT (IN)  ::      &
    zlat, zalt, press    ,& ! latitude [degree], altitude [m], pressure [hPa]
    ztt , ziwv           ,& ! temperature [K], integrated water vapour [mm]
    r_v                  ,& ! gas constant for water vapour (461.51 J /(kg * K))
    degrad                  ! factor for transforming degree to rad

! Local variables:
! ---------------

  REAL (KIND=wp) , PARAMETER  :: &
!   r_v  =  4.619E-3_wp ,& ! hPa * m^3 /(g * K)     =  461.51 J /(kg * K)
!   rhow =  1E6_wp      ,& ! water density  g/m^3
    ck2s = 22.1_wp      ,& ! +/- 2.2    K/hPa   (BEVIS,1994)
    ck3  =  3.739E5_wp     ! +/- 0.012  K^2/hPa

  REAL (KIND=wp) ::      &
    conv_iwv_2_zpd          ,& ! return value of function
    zwd, zhd, tm, pfak, flam   ! to derive IWV from ZPD (routine from Gendt,GFZ)

!------------ End of header ----------------------------------------------------

!-------------------------------------------------------------------------------
! Begin Function conv_iwv_2_zpd
!-------------------------------------------------------------------------------

! IWV conversion into Zenith Wet Delay
  tm   = 70.2_wp + 0.72_wp *ztt
  pfak = (r_v *(ck2s + ck3/tm)) * 1E-5_wp
  zwd  = ziwv * pfak

! Zenith Hydrostatic Delay  (Runge et.al. )
  flam = c1 - 0.00266_wp *COS( c2 *zlat *degrad )                          &
            - 0.28E-6_wp *zalt
  zhd  = 2.2765_wp *press /flam

! ZPD
  conv_iwv_2_zpd  =  zhd  +  zwd

!-------------------------------------------------------------------------------
! End of function conv_iwv_2_zpd
!-------------------------------------------------------------------------------

END FUNCTION conv_iwv_2_zpd



!-------------------------------------------------------------------------------
!+ Module procedure for obs operator for surface pressure from a surface station
!-------------------------------------------------------------------------------

SUBROUTINE ps_obs_operator_sing ( ke, np, col_p, col_z, col_t, col_qv, col_qc  &
                                , col_qrs, r_d, r_g, rdv, zobbdy, lveridat     &
                                , zsobdy , zpsvob, zpsdz, zpsvim )

!-------------------------------------------------------------------------------
!
! Description:
!   This module procedure extrapolates observed surface pressure to the lowest
!   model level (inverse observation operator) and model pressure to the level
!   of the pressure observation (forward observation operator) of a surface-
!   level report.
!
! Method:
!   IF the obs. station lies below the lowest model level 'ke':
!     vertical extrapolation of the obs. to 'ke' with a constant lapse rate and
!     the model temperature at the lowest model layer as starting point
!   IF the obs. station height 'zhsob' is greater than the height at level 'ke':
!     1): vertical interpolation of model pressure to height 'zhsob'
!         ==> obs. increment at height 'zhsob'
!     2): conveyance of this obs. incr. to the lowest model level by a pressure
!         correction, so that the height increment at 'zhsob' is unchanged.
!   This procedure assumes that a (finite) pressure obs exists.
!
! Written by        :  Christoph Schraff, DWD  (original version: 10.07.97)
! Current Code Owner:  Christoph Schraff, DWD
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: European Standards for Writing and Documenting
!                       Exchangeable Fortran 90 Code + DWD Extensions.
!
!-------------------------------------------------------------------------------
!
! Declarations:
!
!-------------------------------------------------------------------------------
! Modules used: These are declared in the module declaration section
! -------------

IMPLICIT NONE

! Subroutine arguments:
! --------------------

  INTEGER        , INTENT (IN)         ::       &
    ke               ,& ! number of vertical (main) levels in model column
    np                  ! number of pressure model fields
                        !   (e.g. 2: original field, boundary field)

  LOGICAL        , INTENT (IN)         ::       &
    lveridat            ! fill simulated obs body (for writing to feedback file)

  REAL (KIND=wp) , INTENT (IN)         ::       &
                        ! model column (at the observation location) of:
    col_p   (ke,np)  ,& ! pressure (full value)       on main levels     [ Pa  ]
    col_z   (ke)     ,& ! geometrical height          of main levels     [  m  ]
    col_t   (ke)     ,& ! temperature                                    [  K  ]
    col_qv  (ke)     ,& ! specific water vapor content                   [kg/kg]
    col_qc  (ke)     ,& ! specific cloud water content                   [kg/kg]
    col_qrs (ke)     ,& ! spec. cont. of hydrometeors excl. cloud water  [kg/kg]
    r_d              ,& ! gas constant for dry air       (287.05)
    r_g              ,& ! acceleration due to gravity    (  9.80665)
    rdv                 ! ratio of gas constant for dry air 'r_d'
                        !      and gas constant for water vapour 'r_v':  r_d/r_v

  REAL (KIND=wp) , INTENT (IN)         ::       &
    zobbdy  (mxsbdy)    ! single-level obs. body (format: see 'osgbdy')

  REAL (KIND=wp) , INTENT (INOUT)      ::       &
    zsobdy  (mxsosg)    ! (single-level) simulated obs. body (format: see SOR)

  REAL (KIND=wp) , INTENT (OUT)        ::  &
    zpsvob  (np)     ,& ! observed pressure interpolated to lowest model level
    zpsdz            ,& ! scaled height distance betw. obs sta. and level 'ke'
    zpsvim  (np)        ! model pressure interpolated to station height

! Local parameters:
! ----------------

  REAL (KIND=wp) , PARAMETER  ::  &
    dhsmin  =  0.5_wp   ! limit of difference between obs. level and model
                            ! orography below which the interpolated value is
                            ! set equal to the observed value

! Local scalars:
! -------------

  INTEGER        ::  &
    klu    , kll     ,& ! upper / lower model level used to calculate lapse rate
    ko               ,& ! model level index
    ip                  ! loop index (over model pressure states)

  REAL (KIND=wp) ::  &
    zhsob            ,& ! height at pressure observation
    zpsob            ,& ! observed (surface) pressure
    dhom             ,& ! zhsob - col_z(ke) ('surface' height difference)
    ztvke            ,& ! virtual temperature at lowest model level
    ztv  (2)         ,& ! virtual temperature at 2 levels used to get lapse rate
    zgamma           ,& ! lapse rate of the virtual temperature
!   ztdb             ,& ! auxiliary quantity (ztvke / dt0lp)
    rvd_m_o             ! r_v / r_d  -  1.0

! Local (automatic) arrays: None
! -------------------------

!
!------------ End of header ----------------------------------------------------


!-------------------------------------------------------------------------------
! Begin Subroutine ps_obs_operator_sing
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!  Section 1: Get height difference between observation and model levels
!-------------------------------------------------------------------------------

  zhsob  =  zobbdy(nbsz)
  zpsob  =  zobbdy(nbsp)
  dhom   =  zhsob - col_z(ke)

  !   very small height difference betw. obs. station and level 'ke':
  !   ---------------------------------------------------------------
  IF (ABS(dhom) <= dhsmin) THEN
    zpsdz  = c0
    DO ip = 1 , np
      zpsvob (ip) = zpsob
      zpsvim (ip) = col_p(ke,ip)
    ENDDO

!-------------------------------------------------------------------------------
!  Section 2: If obs. station lies below level 'ke': extrapolation
!-------------------------------------------------------------------------------

  ELSEIF (dhom < -dhsmin) THEN
    !  (for |dhom| = 100m, temp. error of 12K, the resulting p-error is ~0.5hPa)
!   zpsdz  = - fdoro(2) *dhm
    zpsdz  =             dhom

    rvd_m_o  =  c1 / rdv  -  c1

    CALL lapse_rate_levels ( ke, col_z(:) , klu, kll )
  ! ======================

    ! extrapolate pressure, assume - a constant lapse rate
    !                              - model temperature at level 'ke' is correct
    !   zpsvob =  zpsob * EXP( - g /r /zgamma *LOG( c1 - zgamma *dhom /ztvke ) )
    !   'zgamma': lapse rate of the virtual (!) temperature. This is computed
    !   from model values: temperature difference betw. levels 'kll', 'klu'
    !     (alternative: used temperature gradient as in base state
    !      (cf. reference_atmosphere' in meteo_utilities; Sci. Doc. Sec. 3.1.2):
    !  !  ztdb   = ztvke / dt0lp
    !  !  zpsvob =  zpsob                                                      &
    !  !          * EXP( ztdb * (c1 - SQRT( c1 - c2*g*dhom/(r_d*ztdb*ztvke) )) )

    ztvke   = col_t(ke)  * (c1 + rvd_m_o*col_qv(ke) - col_qc(ke) - col_qrs(ke) )
    ztv (1) = col_t(klu) * (c1 + rvd_m_o*col_qv(klu) -col_qc(klu) -col_qrs(klu))
    ztv (2) = col_t(kll) * (c1 + rvd_m_o*col_qv(kll) -col_qc(kll) -col_qrs(kll))
    zgamma     =  (ztv(1) - ztv(2)) / (col_z(klu) - col_z(kll))
    zpsvob (1) =  zpsob                                                        &
                * EXP( - r_g /r_d /zgamma *LOG( c1 - zgamma *dhom /ztvke ))

    !   model pressure extrapolated to station height (for threshold QC)
    zpsvim (1) = col_p(ke,1) * zpsob / zpsvob(1)

    !   for secondary model pressure fields
    DO ip = 2 , np
      zpsvob (ip) = zpsvob(1)
      zpsvim (ip) = col_p(ke,ip) * zpsob / zpsvob(1)
    ENDDO

!-------------------------------------------------------------------------------
!  Section 3: If obs. station lies above level 'ke':
!             1): interpol. of model-ln(p) to obs. station height ==> obs. incr.
!             2): height correction for conveying this obs. incr. to level 'ke'
!-------------------------------------------------------------------------------

  ELSEIF (dhom > dhsmin) THEN
    zpsdz   = dhom

    ! step 1: interpol. of model-ln(p) to obs station height ==> obs incr.
    ! --------------------------------------------------------------------
    ko    = ke
    DO WHILE ((col_z(ko-1) < zhsob) .AND. (ko > 2))
      ko = ko - 1
    ENDDO

    !   model pressure interpolated to station height (also for threshold QC)
    DO ip = 1 , np
      zpsvim (ip) = col_p(ko,ip) * EXP(   LOG( col_p(ko-1,ip) /col_p(ko,ip) )  &
                                       * (zhsob       - col_z(ko))             &
                                       / (col_z(ko-1) - col_z(ko))     )

      ! step 2): height correction for conveying this obs incr. to level 'ke'
      ! ---------------------------------------------------------------------
      zpsvob (ip) = col_p(ke,ip) + (zpsob - zpsvim(ip)) *( col_p(ke,ip)        &
                                                          /zpsvim(ip))
      ! !   obs incr. conveyed adiabatically
      ! zpsvob (ip) = col_p(ke,ip) + (zpsob - zpsvim(ip))
      !                             *EXP((c0-rdocp) *LOG(col_p(ke)/zpsvim))
    ENDDO
  ENDIF

!-------------------------------------------------------------------------------
!  Section 4: Write simulated observations to SOR 'zsobdy'
!-------------------------------------------------------------------------------

  IF ((lveridat) .AND. (zpsob > rmdich))  zsobdy (nso_p) = zpsvim(1)

!-------------------------------------------------------------------------------
! End of module procedure ps_obs_operator_sing
!-------------------------------------------------------------------------------

END SUBROUTINE ps_obs_operator_sing



!-------------------------------------------------------------------------------
!+ Module procedure for obs operator for surf pressure from a multi-level report
!-------------------------------------------------------------------------------

SUBROUTINE ps_obs_operator_mult ( ke, np, col_p, col_z, col_t, col_qv, col_qc  &
                                , col_qrs, r_d, r_g, rdv, doromx               &
                                , nlev, zobbdy, mzobbd, lpassiv, lveridat      &
                                , zsobdy , zpsvob, zpsdz, zpsvim, ilvvip )

!-------------------------------------------------------------------------------
!
! Description:
!   This module procedure interpolates observed pressure from multi-level
!   pressure / height data to the lowest model level (inverse observation
!   operator), and model pressure to the lowest level a with pressure / height
!   observation at or above the station height from the multi-level report
!   (forward observation operator).
!
! Method:
!   - IF z-p obs. are available below and above the lowest model level 'ke':
!      - Step 1: Inverse observation operator:
!                simple vertical interpolation (ln(p) linear in z) of observed.
!                pressure to level 'ke'     ==> obs. increment at level 'ke'
!      - Step 2: Forward observation operator:
!                extrapolation of the obs. increment at 'ke' to the lowest model
!                level by scaling the increment such that the height increment
!                at 'ke' remains unchanged  ==> obs. increment at height 'zhsob'
!   - IF the lowest z-p obs. (at z=zhsob) is above level 'ke':
!      - Step 1: Forward observation operator:
!                vertical interpolation (ln(p) linear in z) of model pressure
!                to height 'zhsob'          ==> obs. increment at height 'zhsob'
!      - Step 2: Inverse observation operator:
!                extrapolation of this obs. increment to the lowest model level
!                by scaling the increment such that the height increment at
!                'zhsob' remains unchanged  ==> obs. increment at level 'ke'
!
! Written by        :  Christoph Schraff, DWD  (original version: 07.07.97)
! Current Code Owner:  Christoph Schraff, DWD
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: European Standards for Writing and Documenting
!                       Exchangeable Fortran 90 Code + DWD Extensions.
!
!-------------------------------------------------------------------------------
!
! Declarations:
!
!-------------------------------------------------------------------------------
! Modules used: These are declared in the module declaration section
! -------------

IMPLICIT NONE

! Subroutine arguments:
! --------------------

  INTEGER        , INTENT (IN)         ::       &
    ke               ,& ! number of vertical (main) levels in model column
    np               ,& ! number of pressure model fields
                        !   (e.g. 2: original field, boundary field)
    nlev                ! number of vertical levels in report

  LOGICAL        , INTENT (IN)         ::       &
    lpassiv          ,& ! passive data are also used
    lveridat            ! fill simulated obs body (for writing to feedback file)

  REAL (KIND=wp) , INTENT (IN)         ::       &
                        ! model column (at the observation location) of:
    col_p   (ke,np)  ,& ! pressure (full value)       on main levels     [ Pa  ]
    col_z   (ke)     ,& ! geometrical height          of main levels     [  m  ]
    col_t   (ke)     ,& ! temperature                                    [  K  ]
    col_qv  (ke)     ,& ! specific water vapor content                   [kg/kg]
    col_qc  (ke)     ,& ! specific cloud water content                   [kg/kg]
    col_qrs (ke)     ,& ! spec. cont. of hydrometeors excl. cloud water  [kg/kg]
    r_d              ,& ! gas constant for dry air       (287.05)
    r_g              ,& ! acceleration due to gravity    (  9.80665)
    rdv              ,& ! ratio of gas constant for dry air 'r_d'
                        !      and gas constant for water vapour 'r_v':  r_d/r_v
    doromx              ! vertical cut-off radius for extrapolation of pressure
                        !   obs (interpolation range is scaled by fdoro)

  REAL (KIND=wp) , INTENT (IN)         ::       &
    zobbdy  (nlev,mxrbdy)    ! multi-level obs body (format: see 'omlbdy')

  INTEGER        , INTENT (INOUT)      ::       &
    mzobbd  (nlev,mxrbdf)    ! multi-level obs body (format: see 'momlbd')

  REAL (KIND=wp) , INTENT (INOUT)      ::       &
    zsobdy  (mxsoml,nlev)    ! multi-level simulated obs. body (format: see SOR)

  REAL (KIND=wp) , INTENT (OUT)        ::       &
    zpsvob  (np)     ,& ! observed pressure interpolated to lowest model level
    zpsdz            ,& ! scaled height distance level 'ke' and nearest obs.
                        !   (level 'ke' is the lowest model level)
    zpsvim  (np)        ! model pressure interpolated to obs level 'ilvvip(1)'

  INTEGER        , INTENT (OUT)        ::       &
    ilvvip  (2)         ! 1: index of lowest obs level with p-z obs
                        ! 2: index of (lower) obs level used for interpolation

! Local parameters:
! ----------------

  REAL (KIND=wp) , PARAMETER  ::  &
    dhsmin   =  0.5_wp   ! limit of difference between obs. level and model
                             ! orography below which the interpolated value is
                             ! set equal to the observed value

! Local scalars:
! -------------

  INTEGER        ::  &
    klu    , kll     ,& ! upper / lower model level used to calculate lapse rate
    ilev             ,& ! running obs. level index
    ilv              ,& ! index of lowest obs. level with p-z obs.
    ilva   , ilvb    ,& ! indices of neighbouring relevant obs. levels
    nzflg            ,& ! height flags
    ko               ,& ! model level index
    ip                  ! loop index (over model pressure states)

  REAL (KIND=wp) ::  &
    zdz              ,& ! negative height obs increment at lowest z-p obs level
    zhsob            ,& ! observed height   at lowest obs. point with z-p data
    zpsob            ,& ! observed pressure of lowest obs. point with z-p data
    dhom             ,& ! zhsob - col_z(ke) ('surface' height difference)
    dha    , dhb     ,& ! height differ. to neighbouring relevant obs. levels
    ztvke            ,& ! virtual temperature at lowest model level
    ztv  (2)         ,& ! virtual temperature at 2 levels used to get lapse rate
    zgamma           ,& ! lapse rate of the virtual temperature
!   zpmvi            ,& ! model pressure interpolated to relevant obs. level
    rvd_m_o             ! r_v / r_d  -  1.0

! Local (automatic) arrays: None
! -------------------------
!
!------------ End of header ----------------------------------------------------


!-------------------------------------------------------------------------------
! Begin Subroutine ps_obs_operator_mult
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!  Section 1: Get height difference between relevant observat. and model levels
!-------------------------------------------------------------------------------

  ! get surface obs level if available, or lowest obs level with p-z obs
  ! --------------------------------------------------------------------
  ilev = 0
  ilv  = 0
  DO WHILE (ilev < nlev)
    ilev  = ilev + 1
    nzflg = IBITS( mzobbd(ilev,nbtflg), nvfzbp, nvfaoc )
    !    if .not.lpassiv then obs has to be active,
    !                      or the only flag set is the height flag
    !                    (but not the flag for 'below surface')
    IF (      (zobbdy(ilev,nbtz) > rmdich)                                     &
        .AND. (     (      (     (BTEST( mzobbd(ilev,nbterr),nvrz   ))         &
                            .OR. (      (nzflg == ISHFT( 1, nvfbps(4) ))       &
                                  .AND. (.NOT. luse_mlz)))                     &
                     .AND. (.NOT. BTEST( mzobbd(ilev,nbtflg),nvflbp ))         &
                     .AND. (.NOT. BTEST( mzobbd(ilev,nbtqcf),nvrz   )))        &
               .OR. (lpassiv))) THEN
      IF (ilv   == 0) ilv = ilev
      !   level ID has surface-level bit set
      IF (BTEST( mzobbd(ilev,nbtlid), nvlidp(7) )) THEN
        ilv  = ilev
        ilev = nlev
      ENDIF
    ENDIF
  ENDDO
  ilvvip (1)    =  ilv

  IF (ilvvip(1) == 0)  dhom = rmdi
  IF (ilvvip(1) >= 1)  dhom = zobbdy(ilvvip(1),nbtz) - col_z(ke)

  !   very small or missing obs minus model 'surface' height difference
  !   -----------------------------------------------------------------
  IF (ABS(dhom) <= dhsmin) THEN
    zpsdz       =  c0
    zpsvob (:)  =  zobbdy(ilv,nbtp)
    zpsvim (:)  =  col_p(ke,:)
    ko          =  ke
  ELSEIF (ABS(dhom)*fdoro(2) > doromx) THEN
    zpsdz       =  c0
    zpsvob (:)  =  rmdi
    zpsvim (:)  =  rmdi
    dhom        =  c0
  ENDIF

!-------------------------------------------------------------------------------
!  Section 2: If the lowest p-z obs (usually at the obs. station) lies below the
!             lowest model level 'ke': interpolate multi-level p-z obs to 'ke'
!-------------------------------------------------------------------------------

  IF (dhom < -dhsmin) THEN
    !   find closest obs level below / above model level 'ke'
    ilvb = ilv
    ilev = ilv
    DO WHILE (ilev < nlev)
      ilev  = ilev + 1
      nzflg = IBITS( mzobbd(ilev,nbtflg), nvfzbp, nvfaoc )
      IF (      (zobbdy(ilev,nbtz) > rmdich)                                   &
          .AND. (     (      (     (BTEST( mzobbd(ilev,nbterr),nvrz ))         &
                              .OR. (      (nzflg == ishft( 1, nvfbps(4) ))     &
                                    .AND. (.NOT. luse_mlz)))                   &
                       .AND. (.NOT. BTEST( mzobbd(ilev,nbtqcf),nvrz )))        &
                 .OR. (lpassiv))) THEN
        IF (zobbdy(ilev,nbtz) < col_z(ke) ) THEN
          ilvb = ilev
        ELSE
          ilva = ilev
          ilev = nlev + 2
        ENDIF
      ENDIF
    ENDDO
    ilvvip (2) =  ilvb

    !   if both obs levels are found:  interpolate log( p ) linearly in z
    IF (ilev == nlev+2) THEN
      dhb         =  col_z(ke) - zobbdy(ilvb,nbtz)
      dha         =  zobbdy(ilva,nbtz) - col_z(ke)
      zpsdz       =  - MIN( dhb , dha )
      zpsvob (1)  =    zobbdy(ilvb,nbtp)                                       &
                     * EXP( LOG( zobbdy(ilva,nbtp)                             &
                                /zobbdy(ilvb,nbtp)) *dhb /(dhb + dha))
      IF (np >= 2)  zpsvob (2:np) = zpsvob(1)

    !   if no p-z obs. exists above the station height:
    !   extrapolate as for surface-level ps data, cf. 'ps_obs_operator_sing'
    ELSEIF (-dhom < doromx) THEN
      !  (for |dhom| = 100m, T- error of 12K, the resulting p-error is ~0.5hPa)
      zpsdz       =  dhom

      rvd_m_o     =  c1 / rdv  -  c1

      CALL lapse_rate_levels ( ke, col_z(:) , klu, kll )
    ! ======================

      ! extrapolate pressure, assume - a constant lapse rate
      !                              - model temperat. at level 'ke' is correct
      !   zpsvob =  zpsob * EXP( - g/r /zgamma *LOG( c1 - zgamma *dhom /ztvke ))
      !   'zgamma': lapse rate of the virtual (!) temperature. This is computed
      !   from model values: temperature difference betw. levels 'kll', 'klu'

      ztvke  = col_t(ke)  * (c1 + rvd_m_o*col_qv(ke) - col_qc(ke) -col_qrs(ke) )
      ztv(1) = col_t(klu) * (c1 + rvd_m_o*col_qv(klu) -col_qc(klu)-col_qrs(klu))
      ztv(2) = col_t(kll) * (c1 + rvd_m_o*col_qv(kll) -col_qc(kll)-col_qrs(kll))
      zgamma      =  (ztv(1) - ztv(2)) / (col_z(klu) - col_z(kll))
      zpsvob (1)  =  zobbdy(ilvvip(2),nbtp)                                    &
                   * EXP( - r_g /r_d /zgamma *LOG( c1 - zgamma *dhom /ztvke ))
      IF (np >= 2)  zpsvob (2:np) = zpsvob(1)

    !   if no p-z obs. exists within 'doromx(2)':
    ELSE
      zpsdz       =  c0
      zpsvob (:)  =  rmdi
      zpsvim (:)  =  rmdi
    ENDIF

    IF (zpsvob(1) > rmdich) THEN
      IF (np >= 2)  zpsvob (2:np) = zpsvob(1)
      !   for 'zpsvim', extrapolate the increment:
      !   zpsvim = (p(ke)-zpsvob) *zobbdy(ilv1,nbtp)/zpsvob + zobbdy(ilv1,nbtp)
      !          =  p(ke)         *zobbdy(ilv1,nbtp)/zpsvob
      zpsvim (:)  =  col_p(ke,:) * zobbdy(ilvvip(1),nbtp) / zpsvob(:)
      ko  =  ke
    ENDIF

!-------------------------------------------------------------------------------
!  Section 3: If the lowest p-z obs lies above the lowest model level 'ke':
!             1): interpol. of model-ln(p) to obs. station height ==> obs. incr.
!             2): height correction for conveying this obs. incr. to level 'ke'
!-------------------------------------------------------------------------------

  ELSEIF (dhom > dhsmin) THEN
    zpsob = zobbdy (ilvvip(1),nbtp)
    zhsob = zobbdy (ilvvip(1),nbtz)

    zpsdz   = dhom

    ! step 1: interpol. of model-ln(p) to obs station height ==> obs incr.
    ! --------------------------------------------------------------------
    ko    = ke
    DO WHILE ((col_z(ko-1) < zhsob) .AND. (ko > 2))
      ko = ko - 1
    ENDDO

    !   model pressure interpolated to station height (also for threshold QC)
    DO ip = 1 , np
      zpsvim (ip) = col_p(ko,ip) * EXP(   LOG( col_p(ko-1,ip) /col_p(ko,ip) )  &
                                       * (zhsob       - col_z(ko))             &
                                       / (col_z(ko-1) - col_z(ko))     )

      ! step 2): height correction for conveying this obs incr. to level 'ke'
      ! ---------------------------------------------------------------------
      zpsvob (ip) = col_p(ke,ip) + (zpsob - zpsvim(ip)) *( col_p(ke,ip)        &
                                                          /zpsvim(ip))
      ! !   obs incr. conveyed adiabatically
      ! zpsvob (ip) = col_p(ke,ip) + (zpsob - zpsvim(ip))
      !                             *EXP((c0-rdocp) *LOG(col_p(ke)/zpsvim))
    ENDDO
  ENDIF

!-------------------------------------------------------------------------------
!  Section 4: Write simulated observations to SOR 'zsobdy'
!-------------------------------------------------------------------------------

  IF ((lveridat) .AND. (zpsvim(1) > rmdich)) THEN
    !   convert simulated pressure into height increment (using approx. T)
    zdz  =  LOG( zpsvim(1) / zobbdy(ilvvip(1),nbtp) )  * r_d / r_g * col_t(ko)
    !   simulated height obs = height obs + height increment
    zsobdy (nso_p,ilvvip(1))  =  zobbdy(ilvvip(1),nbtz) + zdz
  ENDIF

!-------------------------------------------------------------------------------
! End of module procedure ps_obs_operator_mult
!-------------------------------------------------------------------------------

END SUBROUTINE ps_obs_operator_mult



!-------------------------------------------------------------------------------
!+ Module procedure for observation operator for conventional multi-level obs
!-------------------------------------------------------------------------------

SUBROUTINE mult_obs_operator ( ke, col_u,col_v, col_t_og, col_rh, col_p, col_z &
                             , nlev, zobbdy, mzobbd, kobtyp, lveridat, lqc     &
                             , nupr , zsobdy , vimtoob                         &
                             , lvirt, col_t , zt_o, kbotlev_out, ktoplev_out )

!-------------------------------------------------------------------------------
!
! Description:
!   The routine applies the forward observation operator for multi-level reports
!   which contain vertical  ----------------------------------------------------
!   profiles of horizontal wind, temperature, and / or relative humidity.
!   The forward operator consists here of a simple vertical interpolation
!   of model values to observation levels.
!   Additionally, this simple interpolation operator is also applied to multi-
!   level height data. Note that the resulting height observation increments
!   are not consistent with the temperature observation increments. Consistent
!   height increments are computed later on in the height and thickness quality
!   control check.
!
!   Purpose: In COSMO, the forward observation operator is applied for
!            subsequent:
!              - filling NetCDF feedobs /feedback file for LETKF / verification
!                  (for this purpose, simulated obs are written to the SOR)
!              - quality control (when required)
!              - nudging scheme: spreading of obs. increments from obs. levels,
!                                if required by specification from namelist
!              - nudging scheme: for inverse operator, which interpolates
!                                    obs data to model levels:
!                                a): for humidity (see 'mult_vert_ipol_obs')
!                                b): at the top / base of the obs. profile
!                                c): throughout the obs. profile if no
!                                    significant-level data exist
!
! Method:
!   Vertical interpolation of model values to observation levels:
!     - interpolated variables :  wind components, (virtual) temperature,
!                                 relative humidity,  (height)
!     - interpolation method   :  linear in LOG( pressure )
!
!
! Written by        :  Christoph Schraff, DWD  (original version: 18.09.97)
! Current Code Owner:  Christoph Schraff, DWD
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: European Standards for Writing and Documenting
!                       Exchangeable Fortran 90 Code + DWD Extensions.
!
!-------------------------------------------------------------------------------
!
! Declarations:
!
!-------------------------------------------------------------------------------
! Modules used: These are declared in the module declaration section
! -------------

IMPLICIT NONE

! Subroutine arguments:
! --------------------

  INTEGER        , INTENT (IN)         ::       &
    ke         ,& ! number of vertical (main) levels in model column
    nlev       ,& ! number of vertical levels in report
    kobtyp     ,& ! observation type
    nupr          ! file unit number for control output
                  !   (if < 0 then no control output is written)

  LOGICAL        , INTENT (IN)         ::       &
    lveridat   ,& ! fill simulated obs body (for writing to feedback file)
    lqc           ! threshold quality control to be applied

  REAL (KIND=wp) , INTENT (IN)         ::       &
                        ! model column (at the observation location) of:
    col_u     (ke)   ,& ! zonal wind speed             on Arakawa A grid [ m/s ]
    col_v     (ke)   ,& ! meridional wind speed        on Arakawa A grid [ m/s ]
    col_t_og  (ke)   ,& ! observed variable related to temperature
                        !   (temperature or virtual temperature)         [  K  ]
                        !   (used only for stability-dep. humidity QC thresholds
    col_rh    (ke)   ,& ! model relative humidity                        [     ]
    col_p     (ke)   ,& ! pressure (full value)        on main levels    [ Pa  ]
    col_z     (ke)      ! geometrical height           of main levels    [  m  ]
!   col_hhl   (ke+1) ,& ! geometrical height           of half levels    [  m  ]
!   col_dp0   (ke)   ,& ! reference pressure thickness of model layers   [ Pa  ]
!   col_rho0  (ke)   ,& ! reference density           (on main levels)   [ pa  ]
!   col_lhyphl(ke+1)    ! LOG( hydrostatic pressure [pa] ) on half levels  [ - ]

  REAL (KIND=wp) , INTENT (IN)         ::       &
    zobbdy (nlev,mxrbdy)    ! multi-level obs. body (format: see 'omlbdy')

  INTEGER        , INTENT (IN)         ::       &
    mzobbd (nlev,mxrbdf)    ! multi-level obs. body (format: see 'momlbd')

  REAL (KIND=wp) , INTENT (INOUT)      ::       &
    zsobdy (mxsoml,nlev)    ! multi-level simulated obs. body (format: see SOR)

  REAL (KIND=wp) , INTENT (OUT)        ::       &
    vimtoob (nlev,5)        ! model values interpolated to obs. levels

  LOGICAL        , INTENT (IN)    , OPTIONAL     ::       &
    lvirt                   ! virtual temperature observed, not temperature
                            !   ('lvirt' should exist, if 'zt_o' exists)

  REAL (KIND=wp) , INTENT (IN)    , OPTIONAL     ::       &
    col_t   (ke)            ! model column of temperature    [  K  ]

  REAL (KIND=wp) , INTENT (OUT)   , OPTIONAL     ::       &
    zt_o    (nlev)          ! observation-derived (dry bulb) temperature
                            !   ('zt_o' should exist, ONLY if quality control is
                            !    done subsequently by calling 'mult_obs_qc_fg')

  INTEGER        , INTENT (OUT)   , OPTIONAL     ::       &
    kbotlev_out          ,& ! number of obs levels below lowest    model level
    ktoplev_out             ! number of obs levels above uppermost model level

! Local parameters: None
! ----------------

! Local scalars:
! -------------

  INTEGER        ::  &
    kbotlev          ,& ! number of obs levels below lowest    model level
    ktoplev          ,& ! number of obs levels above uppermost model level
    ilev             ,& ! vertical loop index over observation levels
    km               ,& ! vertical loop index over model levels
    ka     , kb      ,& ! model levels adjacent to the obs.
    kma    , kmb     ,& ! model levels used for extrapol. below lowest main levl
    ivar                ! index of observed quantity: 1=(u,v); 3=T; 4=RH

  REAL (KIND=wp) ::  &
    zlopf               ! weight factor for vertical interpol. to the obs. point

! Local (automatic) arrays:
! -------------------------

  INTEGER        ::  &
    kbz      (nlev)    ! level indices for vertical interpol. to obs levels

  REAL (KIND=wp) ::  &
    zcol_lnp (ke)   ,& ! log( pressure )
    zlpf     (nlev)    ! weight factors for vertical interpolation to obs levels

!------------ End of header ----------------------------------------------------


!-------------------------------------------------------------------------------
! Begin Subroutine mult_obs_operator
!-------------------------------------------------------------------------------

  !   auxiliary quantities
  DO km = 1 , ke
    zcol_lnp (km) = LOG( col_p(km) )
  ENDDO

  !   prepare vertical interpolation from model level to obs level
  !   by computing 'kbotlev', 'ktoplev', 'kbz', and 'zlpf'

  CALL prep_vi_mo2ob ( ke, col_p, zcol_lnp, nlev, zobbdy                       &
                     , kbotlev, ktoplev, kbz, zlpf )
! ==================

  IF (nupr >= 0) WRITE( nupr,'("ilev  kn kf zlopf vm2o1 vm2o2 vm2o3 "          &
                             &," vimtoob4 RH(kb) RH(ka)    p")' )

  ! Linear vertical interpolation in log(p) of u, v, T, RH (and z)
  ! --------------------------------------------------------------
  ! (for model generated 'obs. profiles', store RH(qvc), else RH(qv))

  DO ilev = 1 + kbotlev , nlev
    kb               = kbz (ilev)
    ka               = kbz (ilev) - 1
    zlopf            = zlpf(ilev)
    vimtoob (ilev,1) = (c1 -zlopf  ) * col_u   (kb)   +   zlopf * col_u   (ka)
    vimtoob (ilev,2) = (c1 -zlopf  ) * col_v   (kb)   +   zlopf * col_v   (ka)
!   vimtoob (ilev,3) = (c1 -zlopf  ) * col_t   (kb)   +   zlopf * col_t   (ka)
    vimtoob (ilev,3) = (c1 -zlopf  ) * col_t_og(kb)   +   zlopf * col_t_og(ka)
    vimtoob (ilev,4) = (c1 -zlopf  ) * col_rh  (kb)   +   zlopf * col_rh  (ka)
    !   (non-hydrostatic) computation of height from model pressure
    vimtoob (ilev,5) = (c1 -zlopf  ) * col_z   (kb)   +   zlopf * col_z   (ka)
    !   hydrostatic computation of height from hydrostatic pressure
!   IF (zobbdy(ilev,nbtlop) >  col_lhyphl(kb) kmb = kb
!   IF (zobbdy(ilev,nbtlop) <= col_lhyphl(kb) kmb = ka
!   vimtoob (ilev,5) =  col_hhl(kmb+1)                                         &
!                     + (EXP(col_lhyphl(kmb+1)) - zobbdy(ilev,nbtp))           &
!                       *zrtgkml(icml,kmb) / col_p(kmb)
    IF (nupr >= 0) WRITE( nupr,'(2I4,I3,F5.2,2F6.1,F8.2,3F7.2,F8.0)')          &
                          ilev, kb, ka, zlopf, (vimtoob(ilev,km), km = 1,4)    &
                        , col_rh(kb), col_rh(ka), zobbdy(ilev,nbtp)
  ENDDO

  ! to prepare computation of simulated obs 'zsobdy'
  ! and quality control below lowest model level:
  ! extrapolate between levels (ke-3), ke
  ! -----------------------------------------------

  IF ((kobtyp /= OT_SATEM) .AND. (kbotlev >= 1)) THEN
    DO ilev = kbotlev , 1 , -1
      kmb = ke
      kma = ke - 3
      zlopf   =   (zcol_lnp(kmb) - zobbdy(ilev,nbtlop))                        &
                / (zcol_lnp(kmb) - zcol_lnp(kma))
      vimtoob (ilev,1) = (c1 -zlopf  ) * col_u   (kmb)  +  zlopf * col_u   (kma)
      vimtoob (ilev,2) = (c1 -zlopf  ) * col_v   (kmb)  +  zlopf * col_v   (kma)
!     vimtoob (ilev,3) = (c1 -zlopf  ) * col_t   (kmb)  +  zlopf * col_t   (kma)
      vimtoob (ilev,3) = (c1 -zlopf  ) * col_t_og(kmb)  +  zlopf * col_t_og(kma)
      vimtoob (ilev,4) = (c1 -zlopf  ) * col_rh  (kmb)  +  zlopf * col_rh  (kma)
      !   note that for 'z', this extrapolation assumes that T is constant over
      !   the extrapolation range (from 'kma' to ilev=1): not a good approxim.!!
      vimtoob (ilev,5) = (c1 -zlopf  ) * col_z   (kmb)  +  zlopf * col_z   (kma)
    ENDDO
  ENDIF

  ! write simulated observations to SOR 'zsobdy'
  ! --------------------------------------------

  IF (lveridat) THEN
    DO ilev = 1 , nlev
      IF ((zobbdy(ilev,nbtu ) > rmdich) .AND. (kobtyp /= OT_SATEM)) THEN
!       zsobdy (nso_u ,ilev) = zobbdy(ilev,nbtu ) - vimtoob(ilev,1)
        zsobdy (nso_u ,ilev) = vimtoob(ilev,1)
        zsobdy (nso_v ,ilev) = vimtoob(ilev,2)
      ENDIF
      IF (zobbdy(ilev,nbtt ) > rmdich)                                         &
        zsobdy (nso_t ,ilev) = vimtoob(ilev,3)
      IF (zobbdy(ilev,nbtrh) > rmdich)                                         &
        zsobdy (nso_rh,ilev) = MIN( c1, vimtoob(ilev,4) )
      !   height increments will be replaced in 'zsobdy' later on in the
      !   height and thickness check, (ONLY) IF that check is performed
      !   (i.e. if there are active temperature obs in the profile) !
      !   if (     luse_mlz) then pre-set height increments for real height obs
      !                 to avoid missing increments for active height obs in the
      !                 NetCDF feedobs file (except for the lowest z-p-obs which
      !                                     is approximated here by checking the
      !                                     surface level flag for convenience)
      !   if (.not.luse_mlz) then the only active height obs is the lowest
      !                 z-p-obs (usually the surface level) for which a height
      !                 increment is already set in 'ps_obs_operator_mult';
      !                 then set height increments only for passive real z-obs
      !   hence, pre-set height increments for true height observations to avoid
      !   missing increments for active height obs in the NetCDF feedobs file !!
      IF ((zobbdy(ilev,nbtz ) > rmdich) .AND. (kobtyp /= OT_SATEM)) THEN
        IF (      (      .NOT. BTEST( mzobbd(ilev,nbtflg),nvfzbp+nvfbps(6) ) ) &
            .AND. (     (.NOT. BTEST( mzobbd(ilev,nbtflg),nvfzbp+nvfbps(5) ))  &
                   .OR. (kobtyp /= OT_AIREP))                                  &
            .AND. (     (      (.NOT. luse_mlz)                                &
                         .AND. (.NOT. BTEST( mzobbd(ilev,nbterr),nvrz )     )) &
                   .OR. (      (      luse_mlz)                                &
                         .AND. (.NOT. BTEST( mzobbd(ilev,nbtlsg),LS_SURFACE))) &
                   .OR. (zsobdy(nso_p ,ilev) <= rmdich)))                      &
          zsobdy (nso_p ,ilev) = vimtoob(ilev,5)
      ENDIF
    ENDDO
  ENDIF

  ! re-set extrapolated model values to 'missing value'
  !   (if lqc =.true., then this must be postponed to
  !    after quality control of individual observations)
  ! ----------------------------------------------------

  IF (.NOT. lqc) THEN
    !   for obs levels below the lowest main model level,
    !   set extrapolated model values to 'missing value'
    !   (so that they are not used in the nudging or the height QC)
    DO ilev = 1 , kbotlev
      DO ivar = 1 , 5
        vimtoob (ilev,ivar) = rmdi
      ENDDO
    ENDDO

    !   for obs levels above the top main model level,
    !   set extrapolated model values to 'missing value'
    DO ilev = nlev - ktoplev + 1 , nlev
      DO ivar = 1 , 5
        vimtoob (ilev,ivar) = rmdi
      ENDDO
    ENDDO
  ENDIF

  ! set optional output variables  (which are required only if lqc =.true.)
  ! -----------------------------

  IF (PRESENT( kbotlev_out ))  kbotlev_out = kbotlev
  IF (PRESENT( ktoplev_out ))  ktoplev_out = ktoplev

  !   get (observed) (dry-bulb) temperature profile at observation levels
  !   to prepare for stability-dependent quality control in 'mult_obs_qc_fg'
  IF (PRESENT( zt_o )) THEN
    IF ((lqc) .AND. (PRESENT( lvirt ))) THEN
      DO ilev = 1 , nlev
        IF ((.NOT. lvirt) .OR. (zobbdy(ilev,nbtt) <= rmdich)) THEN
          zt_o (ilev) =   zobbdy(ilev,nbtt)
        ELSE
          zt_o (ilev) =   zobbdy(ilev,nbtt)                                    &
                        * ( (c1 -zlpf(ilev)) * col_t   (kbz(ilev)  )           &
                                             / col_t_og(kbz(ilev)  )           &
                           +     zlpf(ilev)  * col_t   (kbz(ilev)-1)           &
                                             / col_t_og(kbz(ilev)-1))
        ENDIF
      ENDDO
    ELSE
      zt_o (:)  =  zobbdy(:,nbtt)
    ENDIF
  ENDIF

!-------------------------------------------------------------------------------
! End of module procedure mult_obs_operator
!-------------------------------------------------------------------------------

END SUBROUTINE mult_obs_operator



!-------------------------------------------------------------------------------
!+ Module procedure for inverse observation operator for multi-level obs
!-------------------------------------------------------------------------------

SUBROUTINE mult_obs_2_modlev ( ke, col_u, col_v, col_t_og, col_rh, col_p       &
                             ,     col_lhyphl                                  &
                             , nlev, zobbdy, mzobbd, vimtoob, kobtyp, nupr     &
                             , iv1, iv2, lscadj                                &
                             , viobtom, mbotlv, mtoplv, lobinc_out, ilvpr )

!-------------------------------------------------------------------------------
!
! Description:
!   This routine performs the vertical interpolation of multi-level observations
!   to the model levels. This is the main part of the inverse observation
!   operator.
!   This is used in the nudging, but also (only temperature and humidity) for
!   the hydrostatic height / thickness quality control check. In this way, it
!   is also used in the LETKF.
!
! Method:   Vertical interpolation of multi-level observations to model levels
!           ------------------------------------------------------------------
!   - interpolated variables : wind components, (virtual) temperature,
!                              relative humidity
!   - interpolated quantities :
!      - observations   within TEMP / PILOT profiles with significant-level data
!                       in the lower half of the troposphere
!      - obs increments in profiles without significant-level data, including
!                       AIRCRAFT multi-level reports
!   - interpolation method (depending on namelist parameters and feasibility) :
!      - 'vertical scale adjustment': a weighted average of the continuous obs
!                       (increment) profile (piecewise linear in LOG(p))
!                       within the model layer;
!                Note:  This is approx. 'hydrostatically consistent':
!                       If the model temperature equals the resulting scale-
!                       adjusted observed temperature, the hydrostatic pressure
!                       interval between the model layer base and top
!                       (as defined by their fixed height) will then be the same
!                       as for the observed ('true') profile.
!                       (Assumption: the observed specific humidity is constant
!                        within the layer and equal to the model humidity.)
!                       ==> this allows to control the upper-air hydrostatic
!                           pressure by nudging of the 'surface' pressure and
!                           the temperature in this manner.
!      - linear interpolation in LOG( pressure ) to the model main level,
!                       done if: - (.NOT. lscadj), or
!                                - the model layer is not fully covered by the
!                                  continuous obs profile
!      - constant extrapolation of obs increment within a model layer at the
!                       base / top of the profile, if the spreading at that
!                       level is along model levels (e.g. stratosphere)
!
!
! Written by        :  Christoph Schraff, DWD, original version 18.09.97
! Current Code Owner:  Christoph Schraff, DWD
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: European Standards for Writing and Documenting
!                       Exchangeable Fortran 90 Code + DWD Extensions.
!
!-------------------------------------------------------------------------------
!
! Declarations:
!
!-------------------------------------------------------------------------------
! Modules used: These are declared in the module declaration section
! -------------

IMPLICIT NONE

! Subroutine arguments:
! --------------------

  INTEGER        , INTENT (IN)         ::       &
    ke               ,& ! number of vertical (main) levels in model column
    nlev             ,& ! number of vertical levels in report
    kobtyp           ,& ! observation type
                        ! --> Note on 'iv1','iv2': for temperature, output
                        !                          variables are always computed
    iv1              ,& ! must be either: 1: compute output also for wind
                        ! !!!         or  2: do not compute for wind
    iv2              ,& ! must be either: 3: compute output also for humidity
                        ! !!!         or  2: do not compute for humidity
                        ! --> If 'iv1','iv2' are not correctly specified then
                        !     array bounds may be violated !!!
    nupr                ! file unit number for control output
                        !   (if < 0 then no control output is written)

  REAL (KIND=wp) , INTENT (IN)         ::       &
                        ! model column (at the observation location) of:
    col_u     (ke)   ,& ! zonal wind speed             on Arakawa A grid [ m/s ]
    col_v     (ke)   ,& ! meridional wind speed        on Arakawa A grid [ m/s ]
    col_t_og  (ke)   ,& ! observed variable related to temperature
                        !   (temperature or virtual temperature)         [  K  ]
                        !   (used only for stability-dep. humidity QC thresholds
    col_rh    (ke)   ,& ! model relative humidity                        [     ]
    col_p     (ke)   ,& ! pressure (full value)        on main levels    [ Pa  ]
    col_lhyphl(ke+1)    ! LOG( hydrostatic pressure [pa] ) on half levels  [ - ]
                        !   (used for vertical correlations in dz-QC only)

  INTEGER        , INTENT (IN)         ::       &
    mzobbd  (nlev,mxrbdf)    ! multi-level obs. body (format: see 'momlbd')

  REAL (KIND=wp) , INTENT (IN)         ::       &
    zobbdy  (nlev,mxrbdy) ,& ! multi-level obs. body (format: see 'omlbdy')
    vimtoob (nlev,5)         ! model values interpolated to obs. levels

  LOGICAL        , INTENT (IN)         ::       &
    lscadj  (2*iv1-1:iv2+1)  ! method of vertical interpolation:
                             !   .true.  --> vertical scale adjustment
                             !   .false. --> linear interpolation (in log( p ))
                             !   (index: see 'viobtom')

  INTEGER        , INTENT (OUT)        ::       &
    mbotlv  (iv1:iv2)     ,& ! number of model levels above the surface
                             !   without observation increment
                             !   (index =1: wind; =2: temperature; =3: humidity)
    mtoplv  (iv1:iv2)        ! number of model levels at the top of the model
                             !   without observation increment
                             !   (this is used later on only in 'mult_obs_qc_dz'
                             !    and in the nudging)

  REAL (KIND=wp) , INTENT (OUT)        ::       &
    viobtom (ke            & ! obs. (incr.) interpolated to model levels
            ,2*iv1-1:iv2+1)  !   (2nd index = 1, 2 : interpolated wind increments;
                             !              = 3; 4 : full T; RH observed values)
                             !   (this is used later on only in 'mult_obs_qc_dz'
                             !    and in the nudging)

  LOGICAL        , INTENT (OUT)   , OPTIONAL   ::       &
    lobinc_out (iv1:iv2)     ! the quantities to be interpolated are observation
                             !   increments rather than observed values
                             !   (this is used later on only in the nudging)

  INTEGER        , INTENT (OUT)   , OPTIONAL   ::       &
    ilvpr      (ke)          ! obs level used for vertical interpolation
                             !   (used for control output only)

! Local parameters:
! ----------------

  REAL (KIND=wp) , PARAMETER  :: &
    c100r   =     0.01_wp     !

! Local scalars:
! -------------

  INTEGER        ::  &
    ivrs             ,& ! index of observed quantity: 1=(u,v); 2=T; 3=RH
    ivar             ,& ! index of observed quantity: 1=(u,v); 3=T; 4=RH
    ivio             ,& ! index of interpolated quantity: 3=T; 4=RH
    nlevm1           ,& ! number of obs. levels in report minus 1
    ilev             ,& ! vertical loop index over observation levels
    km     , mlev    ,& ! vertical loop index over model levels
    mbotlev          ,& ! number of lowest model levels without interpolated obs
    iallbot          ,& ! lowermost obs. level used for vertical interpolation
    ilevbot (3)      ,& ! lowermost obs. level used for vertical interpolation
    ilevtop (3)      ,& ! uppermost obs. level used for vertical interpolation
    ilvma            ,& ! index of nearest obs. level above the model main level
    ilvaa            ,& ! index of nearest obs. level above the upper half level
    ilvba            ,& ! index of nearest obs. level above the lower half level
    ilva   , ilvb    ,& ! indices of obs. levels containing data of type 'ivrs'
                        ! that are adjacent to the specified model level
    ilvxba , ilvxbb  ,& ! indices of obs. levels containing data of type 'ivrs'
                        ! that are adjacent to the lower model half level
    ilvxaa , ilvxab  ,& ! indices of obs. levels containing data of type 'ivrs'
                        ! that are adjacent to the upper model half level
    ilvs             ,& ! loop index over obs. levels within a model layer
    ilvsb            ,& ! obs. level index used for sum in 'vert. scale adjust.'
    nbtx             ,& ! ODR index for present variable
    nvrx             ,& ! bit pos. for status 'active' for variable
    nphpa            ,& ! pressure [hPa] of current obs. level
    nplast           ,& ! pressure [hPa] of 1 obs. level further below
    nplasi           ,& ! pressure of 1 significant-obs. level further below
    mic                 ! = 1, if interpolated quantity is obs. incr., else = 0

  REAL (KIND=wp) ::  &
    zotmax , zotmay  ,& ! obs. (incr.) (vector) interpolated to upper half level
    zotmbx , zotmby  ,& ! obs. (incr.) (vector) interpolated to lower half level
                        ! or to main level
    zviwp            ,& ! weight to the sum used for 'vertical scale adjustment'
    zobxa  , zobxb   ,& ! obs. (incr.) used for 'vertical scale adjustment'
    zhyphl (2)          ! EXP( col_lhypl ) at upper/lower limit of model layer

  LOGICAL        ::  &
    ldovip           ,& ! do the interpolation of obs. to current model levels
    linvip           ,& ! interpol. of obs. to model levels: linear in LOG( p )
    lsignif          ,& ! report has significant-level data
    lgap             ,& ! 'large' vertical gap betw. adjacent (sign.) obs levels
    lobincx          ,& ! interpolated quantity to model level: obs. incr.
    ltop             ,& ! top of (model or obs.) profile reached
    lobinc (iv1:iv2)    ! the quantities to be interpolated are observation
                        !   increments rather than observed values

! Local (automatic) arrays:
! -------------------------

  REAL (KIND=wp) ::  &
    zlopob   (nlev)  ,& ! LOG( pressure ) of observation levels
    zcol_lnp (ke)       ! log( pressure )

!------------ End of header ----------------------------------------------------


!-------------------------------------------------------------------------------
! Begin Subroutine mult_obs_2_modlev
!-------------------------------------------------------------------------------

  IF (PRESENT( ilvpr ))  ilvpr (:) = 0

  lobinc (iv1:iv2) = .TRUE.
  viobtom (:,:)    =  c0

  nlevm1  = nlev - 1

  !   auxiliary quantities
  DO ilev = 1 , nlev
    zlopob (ilev) = zobbdy(ilev,nbtlop)
  ENDDO
  DO km = 1 , ke
    zcol_lnp (km) = LOG( col_p(km) )
  ENDDO

!-------------------------------------------------------------------------------
!  Section 2: Preparation of the vertical interpolation
!-------------------------------------------------------------------------------

! Get the quantities which determine the mode of the vertical interpolation
! -------------------------------------------------------------------------
! indicator of the quantity ('obs' or 'obs increm.') to be interpolated. This is
! 'obs' if  - there are significant-level data
!      and  - there is no big gap (larger than 200 hPa in the troposphere, or
!             80 hPa around the tropopause and above) between general obs levels
!      and  - there is no big gap (larger than 150 hPa below ca. 700 hPa, or
!             300 hPa further above) between significant-level obs. levels
!      and  - the report is NOT an aircraft report or a satellite retrieval
!             (--> no significant-level data)
!             (for retrievals it makes more sense to interpolate increments
!              since the retrievals themselves contain the model profile as
!              background and therefore have similar fine-scale structures)

  iallbot = nlev
  DO ivrs = iv1 , iv2
    IF (ivrs == 1) nvrx   = nvru
    IF (ivrs == 2) nvrx   = nvrt
    IF (ivrs == 3) nvrx   = nvrq
    lsignif = .FALSE.
    lgap    = .FALSE.
    IF ((kobtyp == OT_TEMP) .OR. (kobtyp == OT_PILOT)) THEN
      nplast  = 0
      DO ilev = 1 , nlev
        IF (BTEST( mzobbd(ilev,nbterr),nvrx )) THEN
          nphpa = NINT( zobbdy(ilev,nbtp) * c100r )
          IF (nplast == 0) nplasi = nphpa
          IF (nplast > 0) nplast = nplast - nphpa
          IF (     (nplast > 200) .OR. (nplasi-nphpa > 300)                    &
              .OR. ((nphpa < 220) .AND. (nplast > 80))                         &
              .OR. ((nphpa > 650) .AND. (nplasi-nphpa > 150))) lgap = .TRUE.
          nplast = nphpa
          IF (    (MOD( nphpa,100 ) > 0) .AND. (nphpa /= 850)                  &
            .AND. (nphpa /= 925) .AND. (nphpa /= 250) .AND. (nphpa /= 150)) THEN
            lsignif = .TRUE.
            nplasi  = nphpa
          ENDIF
        ENDIF
      ENDDO
    ENDIF
    lobinc (ivrs) = (.NOT. lsignif) .OR. (lgap)

! Prepare the start of the vertical loop over model levels
! --------------------------------------------------------

    !   determine the range of obs. levels used for the interpolation
    !   (IF (lobinc) then interpolation of obs. incr, i.e. 'vimtoob' needed)
    ilevbot (ivrs) = nlev
    ilevtop (ivrs) = 0
    ivio = ivrs + MIN( ivrs-1 , 1 )
    DO ilev = 1 , nlev
      IF (      (      BTEST( mzobbd(ilev,nbterr),nvrx ))                      &
          .AND. (.NOT. BTEST( mzobbd(ilev,nbtqcf),nvrx ))                      &
          .AND. ((vimtoob(ilev,ivio) > rmdich) .OR. (.NOT. lobinc(ivrs)))) THEN
        IF (ilevbot(ivrs) == nlev) ilevbot (ivrs) = ilev
        ilevtop (ivrs) = ilev
      ENDIF
    ENDDO
    iallbot = MIN( iallbot , ilevbot(ivrs) )
  ENDDO

  !   indicate the number of model layers above the surface for which
  !   an observed value cannot be determined (limit: zlopob(1))
  km      = ke
  ltop    = .FALSE.
  DO WHILE ((zcol_lnp(km) > zlopob(iallbot)) .AND. (.NOT. ltop))
    DO ivar = 2*iv1 - 1 , iv2 + 1
      viobtom (km,ivar) = rmdi
    ENDDO
    ltop = (km == 1)
    km   = MAX( km - 1 , 1 )
  ENDDO
  mbotlev = ke - km
  DO ivrs = iv1 , iv2
    mbotlv (ivrs) = mbotlev
    mtoplv (ivrs) = 0
  ENDDO

  !   get obs. level just above the lowest model half level (pressure)
  ltop = .FALSE.
  ilev = iallbot
  DO WHILE ((zlopob(ilev) > col_lhyphl(km+1)) .AND. (.NOT. ltop))
    ltop = (ilev == nlev)
    ilev = MIN( ilev , nlevm1 ) + 1
  ENDDO

  ilvaa = ilev

! loop over model layers
! ----------------------

! ______________________________________________
  loop_over_model_layers:  DO mlev = km , 1 , -1
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ! get obs level just above current model level (if not found, skip interpol)
    !   - vertical scale adjustment          : specified level = half level
    !   - vertical interpol. linear in LOG(p): specified level = main level

    ilvba = ilvaa
    ilvma = ilvba
    DO WHILE ((zlopob(ilev) > col_lhyphl(mlev)) .AND. (.NOT. ltop))
      IF (     zlopob(ilev) > zcol_lnp(mlev)) ilvma = ilev + 1
      ltop = (ilev == nlev)
      ilev = MIN( ilev , nlevm1 ) + 1
    ENDDO
    ilvaa = ilev
    IF (ltop) ilvaa = nlev + 1

!-------------------------------------------------------------------------------
!  Section 3: Vertical interpolation of  -->  horizontal wind (u,v)
!             (used only in the nudging, not for the LETKF / feedback file)
!-------------------------------------------------------------------------------

    !   no vertical interpolation of wind, if only thickness QC to be done
    ldovip = (iv1 == 1)
    linvip = .FALSE.
    IF (ldovip)  linvip = (.NOT. lscadj(1))

! No data  ==>  no interpolation
! ------------------------------

    IF (ldovip) THEN
      IF (zcol_lnp(mlev) > zlopob(ilevbot(1))) THEN
        viobtom (mlev,1) = rmdi
        viobtom (mlev,2) = rmdi
        mbotlv (1) = mbotlv(1) + 1
        ldovip = .FALSE.

      !   if model level is above the top (u,v)-obs. (incr) of the obs. profile
      ELSEIF (     (ilevtop(1) == 0)                                           &
              .OR. (zcol_lnp(mlev) < zlopob(MAX( ilevtop(1),1 )))) THEN
        viobtom (mlev,1) = rmdi
        viobtom (mlev,2) = rmdi
        mtoplv (1) = mtoplv(1) + 1
        ldovip = .FALSE.
      ENDIF
    ENDIF
    IF ((ldovip) .AND. (.NOT. linvip)) THEN

! Interpolation by vertical scale adjustment:
! The interpolated value is weighted sum of 'arithmetic means between 2 adjacent
! observed values times the distances between the obs. points in LOG(p)-units)'.
! At the model layer boundaries, 'observed values' are obtained by interpolation
! linear in LOG(p).
! ------------------------------------------------------------------------------

      lobincx = lobinc(1)

      !   get interpolated values at lower / upper boundary of the model layer

      CALL mult_vert_ipol_obs ( ilvba , col_lhyphl(mlev+1) , 1 , lobincx       &
                              , nlev , zobbdy , mzobbd , vimtoob               &
                              , zotmbx , zotmby , ilva , ilvb )
!     =======================

      ilvxba = ilva
      ilvxbb = ilvb
      IF ((ilva > 0) .AND. (ilvb > 0))                                         &

        CALL mult_vert_ipol_obs ( ilvaa , col_lhyphl(mlev) , 1 , lobincx       &
                                , nlev , zobbdy , mzobbd , vimtoob             &
                                , zotmax , zotmay , ilva , ilvb )
!       =======================

      !   interpolate by vertical scale adjustment only if there are wind obs
      !   below and above the whole of the current model layer

      IF ((ilva > 0) .AND. (ilvb > 0)) THEN
        ilvxaa = ilva
        ilvxab = ilvb

        !   compute the weight to the sum
        zviwp = c1 / (col_lhyphl(mlev+1) - col_lhyphl(mlev))

        IF (ilvxba == ilvxaa) THEN
          !   without any wind obs within the model layer
          viobtom (mlev,1) = zviwp * c05 * (zotmax + zotmbx)                   &
                                   * (col_lhyphl(mlev+1) - col_lhyphl(mlev))
          viobtom (mlev,2) = zviwp * c05 * (zotmay + zotmby)                   &
                                   * (col_lhyphl(mlev+1) - col_lhyphl(mlev))
        ELSE
          !   with at least 1 wind obs within the model layer: compute the sum
          mic = 0
          IF (lobincx) mic = 1
          viobtom (mlev,1) =  (zotmbx             + zobbdy(ilvxba,nbtu)        &
                                                  - mic *vimtoob(ilvxba,1))    &
                             *(col_lhyphl(mlev+1) - zlopob(ilvxba))            &
                            + (zotmax             + zobbdy(ilvxab,nbtu)        &
                                                  - mic *vimtoob(ilvxab,1))    &
                             *(zlopob (ilvxab)    - col_lhyphl(mlev))
          viobtom (mlev,2) =  (zotmby             + zobbdy(ilvxba,nbtv)        &
                                                  - mic *vimtoob(ilvxba,2))    &
                             *(col_lhyphl(mlev+1) - zlopob(ilvxba))            &
                            + (zotmay             + zobbdy(ilvxab,nbtv)        &
                                                  - mic *vimtoob(ilvxab,2))    &
                             *(zlopob (ilvxab)    - col_lhyphl(mlev))
          ilvsb = ilvxba
          DO ilvs = ilvxba+1 , ilvxab
            IF (      (      BTEST( mzobbd(ilvs,nbterr),nvru ))                &
                .AND. (.NOT. BTEST( mzobbd(ilvs,nbtqcf),nvru ))) THEN
              viobtom (mlev,1) = viobtom(mlev,1)                               &
                                + (  zobbdy(ilvsb,nbtu)                        &
                                   - mic *vimtoob(ilvsb,1)                     &
                                   + zobbdy(ilvs ,nbtu)                        &
                                   - mic *vimtoob(ilvs ,1))                    &
                                  *(zlopob(ilvsb) - zlopob(ilvs))
              viobtom (mlev,2) = viobtom(mlev,2)                               &
                                + (  zobbdy(ilvsb,nbtv)                        &
                                   - mic *vimtoob(ilvsb,2)                     &
                                   + zobbdy(ilvs ,nbtv)                        &
                                   - mic *vimtoob(ilvs ,2))                    &
                                  *(zlopob(ilvsb) - zlopob(ilvs))
              ilvsb = ilvs
            ENDIF
          ENDDO
          !   weight the sum
          viobtom (mlev,1) = zviwp * c05 * viobtom(mlev,1)
          viobtom (mlev,2) = zviwp * c05 * viobtom(mlev,2)
        ENDIF
      ELSE
        linvip = .TRUE.
      ENDIF
    ENDIF
    IF ((ldovip) .AND. (linvip)) THEN

! Linear interpolation in LOG(p)
! ------------------------------

      lobincx = lobinc(1)

      CALL mult_vert_ipol_obs ( ilvma , zcol_lnp(mlev) , 1 , lobincx           &
                              , nlev , zobbdy , mzobbd , vimtoob               &
                              , zotmbx , zotmby , ilva , ilvb )
!     =======================

      IF ((ilva > 0) .AND. (ilvb > 0)) THEN
        viobtom (mlev,1) = zotmbx
        viobtom (mlev,2) = zotmby
      ELSE
        viobtom (mlev,1) = rmdi
        viobtom (mlev,2) = rmdi
        IF (ilvb == 0) THEN
          mbotlv (1) = mbotlv(1) + 1
        ELSEIF (ilva == 0) THEN
          mtoplv (1) = mtoplv(1) + 1
        ENDIF
      ENDIF
    ENDIF

! If the interpolated quantity is an observed value: compute the obs. increment
! -----------------------------------------------------------------------------

    IF (ldovip) THEN
      IF ((.NOT. lobincx) .AND. (viobtom(mlev,1) > rmdich)) THEN
        viobtom (mlev,1) = viobtom(mlev,1) - col_u(mlev)
        viobtom (mlev,2) = viobtom(mlev,2) - col_v(mlev)
      ENDIF
    ENDIF

!-------------------------------------------------------------------------------
!  Section 4: Vertical interpolation of  -->  (virtual) temperature T
!                                        -->  relative humidity RH
!             (used only in the nudging and the height / thickness QC check,
!              i.e. also used for the LETKF / feedback file)
!-------------------------------------------------------------------------------

    DO ivrs = 2 , iv2

      ivar = ivrs + 1
      ivio = ivar
      IF (ivrs == 2) THEN
        nbtx   = nbtt
        nvrx   = nvrt
      ELSE
        nbtx   = nbtrh
        nvrx   = nvrq
      ENDIF

      ldovip = .TRUE.
      linvip = (.NOT. lscadj(ivar))

! No data  ==>  no interpolation
! ------------------------------

      IF (zcol_lnp(mlev) > zlopob(ilevbot(ivrs))) THEN
        viobtom (mlev,ivio) = rmdi
        mbotlv (ivrs) = mbotlv(ivrs) + 1
        ldovip = .FALSE.

      !    if model level is above the top T / q-obs. (incr) of the obs. profile
      ELSEIF (     (ilevtop(ivrs) == 0)                                        &
              .OR. (zcol_lnp(mlev) < zlopob(MAX( ilevtop(ivrs),1 )))) THEN
        viobtom (mlev,ivio) = rmdi
        mtoplv (ivrs) = mtoplv(ivrs) + 1
        ldovip = .FALSE.
      ENDIF
      IF ((ldovip) .AND. (.NOT. linvip)) THEN

! Interpolation by vertical scale adjustment:
! The interpolated value is weighted sum of 'arithmetic means between 2 adjacent
! observed values times the distances between the obs. points in LOG(p)-units)'.
! At the model layer boundaries, 'observed values' are obtained by interpolation
! linear in LOG(p).
! ------------------------------------------------------------------------------

        lobincx = lobinc(ivrs)

        !   get interpolated values at lower / upper boundary of the model layer

        CALL mult_vert_ipol_obs ( ilvba , col_lhyphl(mlev+1) , ivrs , lobincx  &
                                , nlev , zobbdy , mzobbd , vimtoob             &
                                , zotmbx , zotmby , ilva , ilvb )
!       =======================

        ilvxba = ilva
        ilvxbb = ilvb
        IF ((ilva > 0) .AND. (ilvb > 0))                                       &

          CALL mult_vert_ipol_obs ( ilvaa , col_lhyphl(mlev) , ivrs , lobincx  &
                                  , nlev , zobbdy , mzobbd , vimtoob           &
                                  , zotmax , zotmay , ilva , ilvb )
!         =======================

        !   interpolate by vertical scale adjustment only if there are temper. /
        !   humidity data below and above the whole of the current model layer

        IF ((ilva > 0) .AND. (ilvb > 0)) THEN
          ilvxaa = ilva
          ilvxab = ilvb
          IF ((nupr >= 0) .AND. (mlev <= ke-2) .AND. (mlev >= ke-4))           &
            WRITE( nupr,'("after vertip_ob ",I4,2F8.2,4I4)' )                  &
                   mlev, zotmbx, zotmax, ilvxba, ilvxbb, ilvxaa, ilvxab

          !   compute the weight to the sum
          zviwp = c1 / (col_lhyphl(mlev+1) - col_lhyphl(mlev))

          IF (ilvxba == ilvxaa) THEN
            !   without any temperature / humidity obs. within the model layer:
            viobtom (mlev,ivio) = zviwp * c05 * (zotmax + zotmbx)              &
                                        *(col_lhyphl(mlev+1) - col_lhyphl(mlev))
          ELSE
            !   with at least 1 T/RH obs within the model layer: compute the sum
            zobxa = zobbdy(ilvxba,nbtx)
            zobxb = zobbdy(ilvxab,nbtx)
            IF ((ivrs == 3) .AND. (.NOT. lobincx)) THEN
              !   For saturated obs., take the interpolated model value of
              !   generalized relative humidity (including the cloud water)
              !   (if this is > 100%) in order to get the best estimate of the
              !   total moisture in the model layer and to prevent spurious
              !   drying (e.g. if the top of a cloud is within a model layer).
              !   Finally, the interpolated RH will be set <= 1.
              IF (zobxa >= c1-epsy) zobxa = MAX( c1 , vimtoob(ilvxba,ivio) )
              IF (zobxb >= c1-epsy) zobxb = MAX( c1 , vimtoob(ilvxab,ivio) )
            ELSEIF (ivrs == 3) THEN
              zobxa = zobxa - MIN( c1 , vimtoob(ilvxba,ivio) )
              zobxb = zobxb - MIN( c1 , vimtoob(ilvxab,ivio) )
            ELSEIF (lobincx) THEN
              zobxa = zobxa - vimtoob(ilvxba,ivio)
              zobxb = zobxb - vimtoob(ilvxab,ivio)
            ENDIF
            viobtom (mlev,ivio) =   (zotmbx + zobxb)                           &
                                   *(col_lhyphl(mlev+1) - zlopob(ilvxba))      &
                                  + (zotmax + zobxa)                           &
                                   *(zlopob (ilvxab) - col_lhyphl(mlev ))
            ilvsb = ilvxba
            DO ilvs = ilvxba+1 , ilvxab
              IF (      (      BTEST( mzobbd(ilvs,nbterr),nvrx ))              &
                  .AND. (.NOT. BTEST( mzobbd(ilvs,nbtqcf),nvrx ))) THEN
                zobxb = zobxa
                zobxa = zobbdy(ilvs,nbtx)
                IF ((ivrs == 3) .AND. (lobincx)) THEN
                  zobxa = zobxa - MIN( c1 , vimtoob(ilvs,ivio) )
                ELSEIF ((ivrs == 3) .AND. (zobxa >= c1-epsy)) THEN
                  zobxa = MAX( c1 , vimtoob(ilvs,ivio) )
                ELSEIF (lobincx) THEN
                  zobxa = zobxa - vimtoob(ilvs,ivio)
                ENDIF
                viobtom (mlev,ivio) =  viobtom(mlev,ivio)                      &
                                     + (zobxa + zobxb)                         &
                                      *(zlopob(ilvsb) - zlopob(ilvs))
                ilvsb = ilvs
              ENDIF
            ENDDO
            !   weight the sum
            viobtom (mlev,ivio) = zviwp * c05 * viobtom(mlev,ivio)
          ENDIF
        ELSE
          linvip = .TRUE.
        ENDIF
      ENDIF

! Linear interpolation in LOG(p)
! ------------------------------

      IF ((ldovip) .AND. (linvip)) THEN
        lobincx = lobinc(ivrs)

        CALL mult_vert_ipol_obs ( ilvma , zcol_lnp(mlev) , ivrs , lobincx      &
                                , nlev , zobbdy , mzobbd , vimtoob             &
                                , zotmbx , zotmby , ilva , ilvb )
!       =======================

        IF ((ilva > 0) .AND. (ilvb > 0)) THEN
          viobtom (mlev,ivio) = zotmbx
        ELSE
          viobtom (mlev,ivio) = rmdi
          IF (ilvb == 0) THEN
            mbotlv (ivrs) = mbotlv(ivrs) + 1
          ELSEIF (ilva == 0) THEN
            mtoplv (ivrs) = mtoplv(ivrs) + 1
          ENDIF
        ENDIF
      ENDIF

! If the interpolated quantity is an obs. increment: compute the observed value
! -----------------------------------------------------------------------------

      IF (ldovip) THEN
       IF ((lobincx) .AND. (viobtom(mlev,ivio) > rmdich)) THEN
        IF (ivrs == 2) viobtom (mlev,ivio) = viobtom(mlev,ivio) + col_t_og(mlev)
        IF (ivrs == 3) viobtom (mlev,ivio) = viobtom(mlev,ivio) + col_rh(mlev)
       ENDIF

       !   limit relative humidity to 100%
       IF (ivrs == 3) viobtom (mlev,ivio) = MIN( viobtom(mlev,ivio) , c1 )
      ENDIF

      !    printout for control
      IF ((ldovip) .AND. (nupr >= 0) .AND. (mlev <= ke-2)                      &
                                     .AND. (mlev >= ke-4)) THEN
        IF (linvip) ilvxbb = ilvb
        IF (linvip) ilvxaa = ilva
        DO ilvs = MAX( ilvxbb,1 ) , MIN( ilvxaa,nlev )
          WRITE( nupr,'("pressure / obs. val.: ", F7.0, F8.2)' )               &
                 zobbdy(ilvs,nbtp), zobbdy(ilvs,nbtx)
        ENDDO
        zhyphl(2) = EXP( col_lhyphl(mlev+1) )
        zhyphl(1) = EXP( col_lhyphl(mlev  ) )
        WRITE( nupr,'("lower/upper pressure / intpol obs. val.: ",2F7.0        &
                    &,F8.2)') zhyphl(2), zhyphl(1), viobtom(mlev,ivio)
      ENDIF

    ENDDO

    IF (PRESENT( ilvpr ))  ilvpr (mlev) = ilev

! ____________________________
  ENDDO loop_over_model_layers
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IF (PRESENT( lobinc_out ))  lobinc_out = lobinc

!-------------------------------------------------------------------------------
! End of module procedure mult_obs_2_modlev
!-------------------------------------------------------------------------------

END SUBROUTINE mult_obs_2_modlev



!-------------------------------------------------------------------------------
!+ Module procedure for vertical interpolation of multi-level obs to model level
!-------------------------------------------------------------------------------

SUBROUTINE mult_vert_ipol_obs ( ilev , targlop , ivrs , lobinc                 &
                              , nlev , zobbdy , mzobbd , vimtoob               &
                              , viobtp , viobtpv , ilva , ilvb )

!-------------------------------------------------------------------------------
!
! Description:
!   This routine performs the vertical interpolation of multi-level data of
!   variable 'ivrs' to a specified pressure level.
!
! Method:
!   Vertical interpolation is linear in LOG(pressure), if data of variable
!   'ivrs' exist above and below the specified pressure level.
!   The interpolated quantity is either the observed values or observation
!   increments, depending on the parameter input.
!   For interpolation of relative humidity between a saturated and a sub-
!   saturated observation, the saturated observation is replaced by the inter-
!   polated model value of generalized relative humidity (including the cloud
!   water), if this value exceeds 100%. This is done to prevent spurious
!   drying in such situations.
!
!
! Written by        :  Christoph Schraff, DWD  (25.09.97)
! Current Code Owner:  Christoph Schraff, DWD
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: European Standards for Writing and Documenting
!                       Exchangeable Fortran 90 Code + DWD Extensions.
!
!-------------------------------------------------------------------------------
!
! Declarations:
!
!-------------------------------------------------------------------------------
! Modules used: These are declared in the module declaration section
! -------------

IMPLICIT NONE

! Subroutine arguments:
! --------------------

  INTEGER        , INTENT (IN)         ::       &
    ilev  ,& ! nearest observation level just above the specified pressure level
             ! (this obs. level need not contain a datum of type 'ivrs'
    ivrs  ,& ! index of observed quantity: 1=(u,v); 2=T; 3=RH
    nlev     ! number of vertical levels in report

  REAL (KIND=wp) , INTENT (IN)         ::       &
    targlop  ! pressure (level) to which the obs. (incr.) are to be interpolated

  REAL (KIND=wp) , INTENT (IN)         ::       &
    zobbdy (nlev,mxrbdy) ,& ! multi-level obs. body (format: see 'omlbdy')
    vimtoob (nlev,5)        ! model values interpolated to obs. levels

  INTEGER        , INTENT (IN)         ::       &
    mzobbd (nlev,mxrbdf)    ! multi-level obs. body (format: see 'momlbd')

  LOGICAL        , INTENT (IN)         ::       &
    lobinc   ! the quantities to be interpolated are observation increments
             ! rather than observed values

  INTEGER        , INTENT (OUT)        ::       &
    ilva  ,& ! indices of observation levels used for the interpolation
    ilvb     ! (i.e. the nearest obs. levels above and below the specified tar-
             !  get pressure level 'targlop', which contain data of type 'ivrs')

  REAL (KIND=wp) , INTENT (OUT)        ::       &
    viobtp,& ! interpolated quantity
    viobtpv  !

! Local parameters: None
! ----------------

! Local scalars:
! -------------

  INTEGER        ::  &
    ivar             ,& ! index of observed quantity: 1=(u,v); 3=T; 4=RH
    ivay             ,& ! index of second comp. of obs. quantity: 2=v
    nbtx   , nbty    ,& ! ODR indices for present variable
    nvrx             ,& ! bit pos. for status 'active' for variable
    mic                 ! = 1 if quantit. to be interpol. are obs incr., else =0

  REAL (KIND=wp) ::  &
    zlopf            ,& ! weight factor for vertical interpol. to 'targlop'
    zvimrhb, zvimrha    ! model relative humidity interpolated to obs. levels
                        ! 'ilvb' resp. 'ilva'

  LOGICAL        ::  &
    ltop             ,& ! no data of type 'ivrs' available above 'targlop'
    lbot                ! no data of type 'ivrs' available below 'targlop'

! Local (automatic) arrays: None
! -------------------------
!
!------------ End of header ----------------------------------------------------


!-------------------------------------------------------------------------------
! Begin Subroutine mult_vert_ipol_obs
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!  Section 1: Setting of indices and parameters
!-------------------------------------------------------------------------------

  ivar  = ivrs + MIN( ivrs-1 , 1 )

  viobtpv = rmdi

  IF (ivrs == 1) THEN
    nbtx   = nbtu
    nvrx   = nvru
    nbty   = nbtv
    ivay   = ivar + 1
  ELSEIF (ivrs == 2) THEN
    nbtx   = nbtt
    nvrx   = nvrt
  ELSEIF (ivrs == 3) THEN
    nbtx   = nbtrh
    nvrx   = nvrq
  ENDIF

! Get the adjacent 'obs.' levels containing data of type 'ivrs'
! -------------------------------------------------------------

  !   adjacent level above the specified pressure
  ltop = (ilev > nlev)
  ilva = MIN( ilev , nlev )
  DO WHILE ((     (.NOT. BTEST( mzobbd(ilva,nbterr),nvrx ))                    &
             .OR. (      BTEST( mzobbd(ilva,nbtqcf),nvrx ))) .AND. (.NOT. ltop))
    ltop = (ilva == nlev)
    ilva = MIN( ilva + 1 , nlev )
  ENDDO

  !   adjacent level below the specified pressure
  lbot = (ilva <= 1)
  ilvb = MIN( MAX( ilev - 1 , 1 ) , nlev )
  DO WHILE ((     (.NOT. BTEST( mzobbd(ilvb,nbterr),nvrx ))                    &
             .OR. (      BTEST( mzobbd(ilvb,nbtqcf),nvrx ))) .AND. (.NOT. lbot))
    lbot = (ilvb == 1)
    ilvb = MAX( ilvb - 1 , 1 )
  ENDDO

! Indicate if there is no obs. level above or below
! -------------------------------------------------

! (Never set 'lobinc = .FALSE.' ! I.e., avoid the possibility of interpolating
!  obs. values over a large vertical interval)
! IF ((.NOT. lbot) .AND. (vimtoob(ilvb,ivar) < rmdich)) lobinc = .FALSE.
! IF ((.NOT. ltop) .AND. (vimtoob(ilva,ivar) < rmdich)) lobinc = .FALSE.
  IF ((lobinc) .AND. (vimtoob(ilvb,ivar) < rmdich)) lbot = .TRUE.
  IF ((lobinc) .AND. (vimtoob(ilva,ivar) < rmdich)) ltop = .TRUE.
  mic = 0
  IF (lobinc) mic = 1

!-------------------------------------------------------------------------------
!  Section 2: Perform the vertical interpolation (or extrapolation)
!-------------------------------------------------------------------------------

! Interpolation of 2 values or increments
! ---------------------------------------

  IF ((.NOT. ltop) .AND. (.NOT. lbot)) THEN
    zlopf  =  (zobbdy(ilvb,nbtlop) - targlop)                                  &
            / (zobbdy(ilvb,nbtlop) - zobbdy(ilva,nbtlop))

! Interpolation of (u,v) or T
! ---------------------------

    IF (ivrs <= 2)                                                             &
      viobtp  = (c1-zlopf)* (zobbdy(ilvb,nbtx) - mic*vimtoob(ilvb,ivar))       &
               +    zlopf * (zobbdy(ilva,nbtx) - mic*vimtoob(ilva,ivar))
    IF (ivrs == 1)                                                             &
      viobtpv = (c1-zlopf)* (zobbdy(ilvb,nbty) - mic*vimtoob(ilvb,ivay))       &
               +    zlopf * (zobbdy(ilva,nbty) - mic*vimtoob(ilva,ivay))

! Interpolation of RH:  If data point 'A' is saturated, and data pt. 'B' is not,
!                       then take the interpolated model-RH (if > 100%) for 'A'.
!                       (If (lobinc) then assume a zero increment at 'A'.)
! -------------------

    IF (ivrs == 3) THEN
      zvimrhb = MIN( c1 , vimtoob(ilvb,ivar) )
      zvimrha = MIN( c1 , vimtoob(ilva,ivar) )
      IF (      (zobbdy(ilvb,nbtx) >= c1-epsy)                                 &
          .AND. (zobbdy(ilva,nbtx) <  c1-epsy)                                 &
          .AND. (vimtoob(ilvb,ivar) > rmdich)) THEN
        viobtp =  (c1-zlopf) * (1-mic) * MAX( c1 , vimtoob(ilvb,ivar) )        &
                 +    zlopf  * (zobbdy(ilva,nbtx) - mic *zvimrha)
      ELSEIF (      (zobbdy(ilva,nbtx) >= c1-epsy)                             &
              .AND. (zobbdy(ilvb,nbtx) <  c1-epsy)                             &
              .AND. (vimtoob(ilva,ivar) > rmdich)) THEN
        viobtp =  (c1-zlopf) * (zobbdy(ilvb,nbtx) - mic *zvimrhb)              &
                 +    zlopf  * (1-mic) * MAX( c1 , vimtoob(ilva,ivar) )
      ELSE
        viobtp =  (c1-zlopf) * (zobbdy(ilvb,nbtx) - mic *zvimrhb)              &
                 +    zlopf  * (zobbdy(ilva,nbtx) - mic *zvimrha)
      ENDIF
    ENDIF

! Indicate if insufficient observations are available for interpolation
! ---------------------------------------------------------------------

  ELSE
    IF (lbot) ilvb   = 0
    IF (ltop) ilva   = 0
    viobtp = rmdi
  ENDIF

!-------------------------------------------------------------------------------
! End of module procedure mult_vert_ipol_obs
!-------------------------------------------------------------------------------

END SUBROUTINE mult_vert_ipol_obs



!-------------------------------------------------------------------------------
!+ Module procedure for supplementary observation operator for multi-level T, z
!-------------------------------------------------------------------------------

SUBROUTINE mult_obs_operator_z ( ke, col_t, col_p, col_z, col_tv, col_hhl      &
                               , nlev, zobbdy, mzobbd, zobdps, viobtom         &
                               , mbotlv, mtoplv, kobtyp, r_d, r_g, lveridat    &
                               , zsobdy , dzob, lqcdz , vcs_in, vcut_in )

!-------------------------------------------------------------------------------
!
! Description:
!   This routine computes height observation increments, which are consistent
!   with the temperature observation increments. The corresponding model-derived
!   height is then written to the Simulated Observation Record SOR.
!   Together with the necessary call of routine 'mult_obs_2_modlev', this
!   routine is a supplementary part to the forward observation operator for
!   multi-level temperature / height reports.
!
! Method:
!   1.   Compute on model levels vertical profile of temperature increments
!        (including increments given by the vertical weight function below and
!        above the profile),
!   2. - then convert it into a pressure increment profile using the exact
!        formulation of the hydrostatic equation in the COSMO model (and taking
!        into account a 'near-surface' pressure increment, if it exists),
!      - then convert it into a height increment profile,
!   3.   finally interpolate the height increments to the T-obs levels.
!
! Written by        :  Christoph Schraff, DWD  (original version: 1998)
! Current Code Owner:  Christoph Schraff, DWD
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: European Standards for Writing and Documenting
!                       Exchangeable Fortran 90 Code + DWD Extensions.
!
!-------------------------------------------------------------------------------
!
! Declarations:
!
!-------------------------------------------------------------------------------
! Modules used: These are declared in the module declaration section
! -------------

IMPLICIT NONE

! Subroutine arguments:
! --------------------

  INTEGER        , INTENT (IN)         ::       &
    ke               ,& ! number of vertical (main) levels in model column
    nlev             ,& ! number of vertical levels in report
    mbotlv           ,& ! number of model levels above the surface without
                        !   temperature observation increment
    mtoplv           ,& ! number of model levels at the top of the model without
                        !   temperature observation increment
    kobtyp              ! observation type

  LOGICAL        , INTENT (IN)         ::       &
    lveridat            ! fill Simulated Obs Record SOR (for feedback files)

  REAL (KIND=wp) , INTENT (IN)         ::       &
                        ! model column (at the observation location) of:
    col_t     (ke)   ,& ! temperature                                    [  K  ]
    col_p     (ke)   ,& ! pressure (full value)        on main levels    [ Pa  ]
    col_z     (ke)   ,& ! geometrical height           of main levels    [  m  ]
    col_hhl   (ke+1) ,& ! geometrical height           of half levels    [  m  ]
!   col_dp0   (ke)   ,& ! reference pressure thickness of model layers   [ Pa  ]
!   col_rho0  (ke)   ,& ! reference density           (on main levels)   [ pa  ]
    col_tv    (ke)   ,& ! virtual temperature                            [  K  ]
    r_d              ,& ! gas constant for dry air       (287.05)
    r_g                 ! acceleration due to gravity    (  9.80665)

  REAL (KIND=wp) , INTENT (IN)         ::       &
    zobdps               ,& ! pressure obs increment valid at lowest model level
    viobtom (ke)         ,& ! T-obs. interpolated to model levels
    zobbdy  (nlev,mxrbdy)   ! multi-level obs. body (format: see 'omlbdy')

  INTEGER        , INTENT (IN)         ::       &
    mzobbd  (nlev,mxrbdf)   ! multi-level obs. body (format: see 'momlbd')

  REAL (KIND=wp) , INTENT (INOUT)      ::       &
    zsobdy  (mxsoml,nlev)   ! multi-level simulated obs. body (format: see SOR)

  REAL (KIND=wp) , INTENT (OUT)   , OPTIONAL   ::       &
    dzob    (nlev)          ! profile of height observation increments

  LOGICAL        , INTENT (OUT)   , OPTIONAL   ::       &
    lqcdz                   ! height and thickness quality control to be done

  REAL (KIND=wp) , INTENT (IN)    , OPTIONAL   ::       &
    vcs_in  (2)          ,& ! vertical correlation scale for temperature [  m  ]
    vcut_in (2)             ! height range influenced by T-increments    [  m  ]
                            !   (index 1: below profile; index 2: above profile)
                            ! --> if 'vcs_in' or 'vcut_in' do not exist then
                            !     compute corresponding values (scale, range)
                            !     as deployed operationally in the nudging
                            !     scheme (using parameter 'f_vcs')

! Local parameters:
! ----------------

! REAL (KIND=wp) , PARAMETER  :: &
!   f_vcs  =  5.854190779_wp  ! factor [m/K] for vertical correlation scale,
                                  ! to be multiplied with virtual temperature
                                  !   (f_vcs = 0.2 *r_d/g)

! Local scalars:
! -------------

  INTEGER        ::  &
    kbotlev          ,& ! number of obs levels below lowest    model level
    ktoplev          ,& ! number of obs levels above uppermost model level
    ilev             ,& ! vertical loop index over observation levels
    mlev             ,& ! vertical loop index over model levels
    mdabot , mdatop  ,& ! index of lowest / top model level with interpol. obs.
    moibot , moitop     ! lowest / top model level used for thickness check

  REAL (KIND=wp) ::  &
!   zfbuoyt          ,& ! factor in the buoyancy term in the w-equation
!   zfbuoyb          ,& ! as 'zfbuoyt' but for a model level further below
    zvcs  (2)        ,& ! vertical correlation scale below / above profile
    zcut  (2)           ! lower / upper height limit influenced by T-increments

  LOGICAL        ::  &
    lex_vcs          ,& ! optional input variables 'vcs_in', 'vcut_in' exist
    ldodz               ! height and thickness quality control to be done

! Local (automatic) arrays:
! -------------------------

  INTEGER        ::  &
    kbz     (nlev)      ! level indices for vertical interpol. to obs levels

  REAL (KIND=wp) ::  &
    zlpf    (nlev)   ,& ! weight factors for vertical interpolation to obs levels
    zcol_lnp(ke)     ,& ! log( pressure )
    ztoi    (ke)     ,& ! temperature obs. increment profile at model layers
    zdz     (ke)     ,& ! hydrostatic upper-air height observation increments
    zdlnp   (ke+1)   ,& ! LOG( p ) increments at half levels
    zdlnpml (ke)     ,& ! LOG( p ) increments at main levels
!   zdz2    (ke)     ,& ! new height increments at main levels
!   zdzob2  (nlev)   ,& ! profile of height observation increments
!   zdzob3  (nlev)   ,& ! profile of height observation increments
    zdzob   (nlev)      ! profile of height observation increments

!------------ End of header ----------------------------------------------------


!-------------------------------------------------------------------------------
! Begin Subroutine mult_obs_operator_z
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!  Section 0: Preliminaries
!-------------------------------------------------------------------------------

  zdzob (:) = c0

! first decide, whether computation of height increments (and the
! height / thickness check) shall be done for the current report:
! skip this for profiles without (active !) temperature data or
! very short temperature profiles (vertical extent <= 50 hPa)
! ---------------------------------------------------------------

  ldodz  =  (mtoplv+mbotlv < ke)

  IF (ldodz) THEN
    mdatop  =  1 + mtoplv
    mdabot  = ke - mbotlv
    ldodz   = (col_p(mdabot)-col_p(mdatop) >= 5000.0_wp)
  ENDIF

  IF (ldodz) THEN

!-------------------------------------------------------------------------------
!  Section 1: Construction of an (interpolated) observed temperature profile;
!             besides the height and thickness check (for multi-level T),
!             this could be (but is not any more) used for the IWV check
!-------------------------------------------------------------------------------

    !   check for valid input values 'vcs_in' and 'vcut_in'
    lex_vcs  =  (PRESENT( vcs_in )) .AND. (PRESENT( vcut_in ))
    IF (lex_vcs) THEN ; lex_vcs =       (MIN( vcs_in(1), vcs_in(2) ) > -epsy)        &
                            .AND. (vcut_in(1) > -10000._wp)                &
                            .AND. (vcut_in(2) >   -400._wp)                &
                            .AND. (vcut_in(1) < vcut_in(2))
    ENDIF

    !   if no valid input values exist, then compute 'zvcs', 'zcut' by using
    !   the operational function from nudging scheme (given here by parameters,
    !   neglecting the adjustment of the vertical weight function for
    !   piecewise profiles)
    IF (lex_vcs) THEN
      zvcs (:)  =  vcs_in(:)
      zcut (:)  =  vcut_in(:)
    ELSE
      zvcs (1)  =  0.2_wp *r_d /r_g * col_tv(ke-mbotlv)
      zvcs (2)  =  0.2_wp *r_d /r_g * col_tv(1 +mtoplv)
      zcut (1)  =  col_z(ke-mbotlv)  -  zvcs(1)   !   (vcutof = 1)
      zcut (2)  =  col_z(1+ mtoplv)  +  zvcs(2)   !   (vcutof = 1)
         !   (vcutof: vertical cut-off radias in terms of scale of weight fn.)
    ENDIF

    DO mlev = ke , 1 , -1
      ztoi  (mlev) = c0
    ENDDO

    mdatop  =  1 + mtoplv
    mdabot  = ke - mbotlv
    DO mlev = mdatop , mdabot
      ztoi (mlev) = viobtom(mlev) - col_t(mlev)
    ENDDO
    moitop  = mdatop
    moibot  = mdabot
!   PRINT       '("sta. ",A,2X,": mbotoplv",6I4)', ystid                       &
!               , mdatop, mtoplv, moitop, mdabot, mbotlv, moibot

    !   below the observed profile, use the vertical weight function of nudging
    IF ((mdabot < ke) .AND. (mdabot >= mdatop)) THEN
      moibot  = mdabot
      DO mlev = mdabot + 1 , ke
        IF (col_z(mlev) >= zcut(1)) THEN
          ztoi (mlev)  = ztoi(mdabot) * EXP( -(  (col_z(mdabot) - col_z(mlev)) &
                                               / zvcs(1)) **2 )
          moibot       = mlev
        ENDIF
      ENDDO
      IF (moibot < ke)  ztoi (moibot+1)  =  c0
      moibot = MIN( moibot + 1 , ke )
    ENDIF

    !   above the observed profile, use the vertical weight function of nudging
    IF ((mdatop > 1) .AND. (mdatop <= mdabot)) THEN
      moitop  = mdatop
!     IF (mdatop > 1)  ztoi (mdatop-1)  =  c0
      DO mlev = mdatop - 1 , 1 , - 1
        IF (col_z(mlev) <= zcut(2)) THEN
          ztoi (mlev)  = ztoi(mdatop) * EXP( -(  (col_z(mlev) - col_z(mdatop)) &
                                               / zvcs(2)) **2 )
          moitop       = mlev
        ENDIF
      ENDDO
      IF (moitop > 1)  ztoi (moitop-1)  =  c0
!     moitop = MAX( MIN( moitop , mdatop - 1 ) , 1 )
      moitop = MAX( moitop - 1 , 1 )
    ENDIF

!   IF (ystid == '10739') PRINT '("IQ4a ",A,3I5,L3)',ystid, mtoplv, mbotlv     &
!!                              , ntstep, lqc

!-------------------------------------------------------------------------------
!  Section 2a: Computation of a profile of height obs increments at model levels
!              using the accurate COSMO formulation of the hydrostatic equation
!              (not deployed in order to avoid use of fields 'dp0', 'rho0')
!-------------------------------------------------------------------------------

! computation of a profile of pressure observation increments
! -----------------------------------------------------------

!   zdp (moibot) = c0
!!  IF (.NOT. lretv) THEN
!!    !   element 'nhtvip' is defined only if .not.lretv
!!    IF ((moibot == ke) .AND. (omlhed(nmlob,nhtvip) > c0))                    &
!!      zdp (moibot) = omlhed(nmlob,nhtvip) - col_p(ke)
!     IF (moibot == ke)  zdp (moibot) = zobdps
!!  ENDIF
!   zfbuoyt      = r_d * col_rho0(moibot) * col_t(moibot)
!   DO mlev = moibot-1 , moitop , -1
!     zfbuoyb    = zfbuoyt
!     zfbuoyt    = r_d * col_rho0(mlev) * col_t(mlev)
!     zdp (mlev) =          c1       / (c2 + col_dp0(mlev+1) /zfbuoyt)         &
!                 * (   zdp (mlev+1) * (c2 - col_dp0(mlev  ) /zfbuoyb)         &
!                    +  ztoi(mlev)   * col_dp0(mlev+1) /col_t(mlev  )          &
!                    +  ztoi(mlev+1) * col_dp0(mlev  ) /col_t(mlev+1) )
!   ENDDO
!   DO mlev = moitop-1 , 1 , -1
!     zfbuoyb    = zfbuoyt
!     zfbuoyt    = r_d * col_rho0(mlev) * col_t(mlev)
!     zdp (mlev) =          c1       / (c2 + col_dp0(mlev+1) /zfbuoyt)         &
!                 * (   zdp (mlev+1) * (c2 - col_dp0(mlev  ) /zfbuoyb) )
!   ENDDO

! conversion into a profile of height observation increments (at model levels)
! ----------------------------------------------------------------------------

!   DO mlev = ke , moibot+1 , -1
!     zdz2 (mlev) = c0
!   ENDDO
!!  DO mlev = moibot , moitop , -1
!   DO mlev = moibot , 1 , -1
!     zdz2 (mlev) = zdp(mlev) * r_d * col_t(mlev) / (r_g * col_p(mlev))
!   ENDDO

!-------------------------------------------------------------------------------
!  Section 2b: Computation of a profile of height obs increments at model levels
!              using a standard formulation of the hydrostatic equation
!-------------------------------------------------------------------------------

    zdlnp (moibot+1)  =  c0
    !   LOG( p ) increments at half levels
    !   - level ke-1: instead of 'col_p(ke)', surface pressure would be required
    !                 to obtain accurate increment of LOG( surface pressure)
    IF (moibot == ke)  zdlnp (ke+1) = LOG( (col_p(ke) + zobdps)/ col_p(ke) )

    !   - half levels aloft: use simple formulation of hydrostatic equation
    DO mlev = moibot , 1 , -1
      zdlnp (mlev) = zdlnp(mlev+1) - r_g/r_d *(  c1/(col_t(mlev)+ztoi(mlev))   &
                                               - c1/ col_t(mlev))              &
                                             *(col_hhl(mlev) - col_hhl(mlev+1))
    ENDDO

    !   LOG( p ) increments at main levels, by simple arithmetic averaging
    DO mlev = moibot , 1 , -1
      zdlnpml (mlev) = c05* (zdlnp(mlev+1) + zdlnp(mlev))
    ENDDO

    !   z increments (at main levels), using simple form. of hydrostatic eq.
    !   (note the sign: positive height increments for positive pressure incr.)
    zdz (:) = c0
    DO mlev = moibot , 1 , -1
      zdz (mlev) = zdlnpml(mlev) * r_d * (col_t(mlev) + c05*ztoi(mlev))/ r_g
    ENDDO

!-------------------------------------------------------------------------------
!  Section 3 : interpolation of height observation increments to obs levels
!-------------------------------------------------------------------------------

    !   prepare by computing neighbouring model levels 'kbz' and factors 'zlpf'
    DO mlev = 1 , ke
      zcol_lnp (mlev) = LOG( col_p(mlev) )
    ENDDO

    CALL prep_vi_mo2ob ( ke, col_p, zcol_lnp, nlev, zobbdy                     &
                       , kbotlev, ktoplev, kbz, zlpf )
!   ==================

    !   interpolate to obs levels
    DO ilev = 1 , nlev
      zdzob (ilev) = (c1-zlpf(ilev)) *zdz(kbz(ilev))                           &
                    +    zlpf(ilev)  *zdz(kbz(ilev)-1)
!!    zdzob2 (ilev) = (c1-zlpf(ilev)) *zdz2(kbz(ilev))                         &
!!                   +    zlpf(ilev)  *zdz2(kbz(ilev)-1)
!!    zdzob3 (ilev) =   -999._wp
!!    IF (zobbdy(ilev,nbtz) > rmdich)                                          &
!!      zdzob3 (ilev) =   zobbdy(ilev,nbtz)                                    &
!!                      - (c1-zlpf(ilev)) *col_z(kbz(ilev))                    &
!!                      -     zlpf(ilev)  *col_z(kbz(ilev)-1)
!!    IF (nlev > 20)                                                           &
!!      PRINT *,'ZZdz ', nlev, zobbdy(1,nbtz), ilev, zdzob(ilev), zdzob2(ilev) &
!!                                                 , zdzob3(ilev)
    ENDDO

!-------------------------------------------------------------------------------
!  Section 4 : Writing model-derived height to SOR, which is consistent with
!              the temperature observation increments
!-------------------------------------------------------------------------------

    IF (lveridat) THEN
      DO ilev = 1 , nlev
        !   note: at the lowest level with z-p obs, the value computed
        !         previously in 'ps_obs_operator_mult' is .NOT. replaced here
        !         (use same if condition as in 'mult_obs_operator'
        IF (zobbdy(ilev,nbtz) > rmdich) THEN
          IF (     (     .NOT. BTEST( mzobbd(ilev,nbtflg),nvfzbp+nvfbps(6) ) ) &
             .AND. (    (.NOT. BTEST( mzobbd(ilev,nbtflg),nvfzbp+nvfbps(5) ))  &
                   .OR. (kobtyp /= OT_AIREP))                                  &
             .AND. (    (      (.NOT. luse_mlz)                                &
                         .AND. (.NOT. BTEST( mzobbd(ilev,nbterr),nvrz )     )) &
                   .OR. (      (      luse_mlz)                                &
                         .AND. (.NOT. BTEST( mzobbd(ilev,nbtlsg),LS_SURFACE))) &
                   .OR. (zsobdy(nso_p ,ilev) <= rmdich)))                      &
            zsobdy (nso_p,ilev)  =  zobbdy(ilev,nbtz) - zdzob(ilev)
        ENDIF
      ENDDO
    ENDIF

  ENDIF

  !   set output variables, if present
  IF (PRESENT( dzob  ))  dzob (:)  =  zdzob(:)
  IF (PRESENT( lqcdz ))  lqcdz     =  ldodz

!-------------------------------------------------------------------------------
! End of module procedure mult_obs_operator_z
!-------------------------------------------------------------------------------

END SUBROUTINE mult_obs_operator_z


!===============================================================================

ELEMENTAL REAL (KIND=wp) FUNCTION fpvsw  ( zt, b1, b2w, b3, b4w )
  !----------------------------------------------------
  REAL    (KIND=wp)         , INTENT (IN)  ::  zt, b1, b2w, b3, b4w
  !---------------------------------------------------------------------------
  ! Magnus formula for water:  input  'zt'   : temperature
  !                            output 'fpvsw': saturation water vapour pressure
  ! Magnus formula for ice  :  if constants 'b2i', 'b4i' for ice are used for
  !                            'b2w', 'b4w' in the call
  !---------------------------------------------------------------------------
  !
  fpvsw  =  b1 * EXP( b2w*(zt-b3)/(zt-b4w) )
  !
END FUNCTION fpvsw

!-------------------------------------------------------------------------------

ELEMENTAL REAL (KIND=wp) FUNCTION ftd  ( zpv, b1, b2w, b3, b4w )
  !---------------------------------------------------
  REAL    (KIND=wp)         , INTENT (IN)  ::  zpv, b1, b2w, b3, b4w
  REAL    (KIND=wp)                        ::  zlogpv
  !---------------------------------------------------------------------------
  ! inverse of Magnus formula:  input  'zpv' :  water vapour pressure
  !                             output 'ftd' :  dewpoint temperature
  !---------------------------------------------------------------------------
  !
  zlogpv  =  LOG( zpv / b1 )
  ftd     =  (b3*b2w - zlogpv*b4w) / (b2w - zlogpv)
  !
END FUNCTION ftd

!-------------------------------------------------------------------------------

ELEMENTAL REAL (KIND=wp) FUNCTION fpv2q  ( zpv, zp, rdv )
  !--------------------------------------------
  REAL    (KIND=wp)         , INTENT (IN)  ::  zpv, zp, rdv
  !---------------------------------------------------------------------------
  ! specific humidity from water vapour pressure and air pressure
  !---------------------------------------------------------------------------
  !
  fpv2q  =  rdv * zpv / (zp - (c1-rdv)*zpv)
  !
END FUNCTION fpv2q

!-------------------------------------------------------------------------------

ELEMENTAL REAL (KIND=wp) FUNCTION fq2pv  ( zqv, zp, rdv )
  !--------------------------------------------
  REAL    (KIND=wp)         , INTENT (IN)  ::  zqv, zp, rdv
  !---------------------------------------------------------------------------
  ! water vapour pressure (at T = saturation vapour press. at Td)
  ! from specific humidity and air pressure
  !---------------------------------------------------------------------------
  !
  fq2pv  =  MAX( epsy , zqv ) * zp / (rdv + zqv*(c1-rdv))
  !
END FUNCTION fq2pv

!-------------------------------------------------------------------------------

END MODULE src_obs_operator_conv
