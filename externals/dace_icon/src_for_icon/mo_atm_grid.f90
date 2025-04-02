!
!+ Derived type definition to hold information on model grids
!
MODULE mo_atm_grid
!
! Description:
!   Defines data type 't_grid' to hold information on the model grid
!   and routines acting on this data type.
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
! V1_4         2009/03/26 Hendrik Reich
!  Changes for COSMO rotated lat-lon grid
! V1_7         2009/08/24 Harald Anlauf
!  construct_atm_grid: zero 'bkf' when using isobaric grid template
! V1_8         2009/12/09 Andreas Rhodin
!  subroutine print(grid): print correctly 3-d fields distributed over PEs
! V1_9         2010/04/20 Andreas Rhodin
!  bug fix for rev.1.65 (seriell/parallel printout)
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Andreas Rhodin
!  update meteo_utilities to COSMO revision V4_11
!  new component: grid% ptopf
! V1_20        2012-06-18 Harald Anlauf
!  construct_atm_grid: handle corner case for single-level field
!  print_atm_grid, gridinfo: Bugfix for single-level fields
! V1_22        2013-02-13 Harald Anlauf
!  initial support for ICON grid
! V1_23        2013-03-26 Harald Anlauf
!  implement ICON unstructured grid
! V1_26        2013/06/27 Harald Anlauf
!  Changes for GRIB2/GRIB_API/ICON
!  Current domain decomposition for ICON requires nproc2=1, enforce it.
! V1_27        2013-11-08 Harald Anlauf
!  Fixes for generalized vertical coordinate
! V1_28        2014/02/26 Harald Anlauf
!  Increase tolerances for detection of global latlon-grid,
!  GRIB2: properly read/pass/write uuidOfVGrid, numberOfVGridUsed
! V1_35        2014-11-07 Andreas Rhodin
!  vector version for destruct_atm_grid;
!  init_ctl_from_grid: add optional parameter 'tmpl'
! V1_42        2015-06-08 Andreas Rhodin
!  new component grid% d_gme; ICON local patch
! V1_43        2015-08-19 Andreas Rhodin
!  add 'fr_lake', 'depth_lk' to grid derived type
! V1_44        2015-09-30 Andreas Rhodin
!  construct_grid: set up 'xnglob' for grid DWD6_NONE
! V1_45        2015-12-15 Harald Anlauf
!  fixes for TR15581; print_atm_grid: adjust format for IFS w/ 137 levels
! V1_46        2016-02-05 Andreas Rhodin
!  revise setup of ivctype, refatm for GRIB2 usage
! V1_47        2016-06-06 Andreas Rhodin
!  fix ivctype lookup from defaults (for COSMO EU); add 'VCT_NO'
! V1_48        2016-10-06 Harald Anlauf
!  t_grid: domain decomposition info for dual grid
!  set_zlev_ref: estimate unperturbed model level heights for gen.vert.coor.
!  support COSMO ivctype=3,4 (sleeve coordinates)
!  determine mean lapse rate for SYNOP T2M correction
! V1_50        2017-01-09 Andreas Rhodin
!  t_grid: new components d_km,d_deg,htopf; adjust ni for ICON-Nest or ICON-LAM
! V1_51        2017-02-24 Andreas Rhodin
!  update shared modules from COSMO 5.04d
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Andreas Rhodin  MPIfM 1998-2002 original code
! Andreas Rhodin  DWD   2003-2008
! Harald Anlauf   DWD   2007-2008
!========================================================================
#include "tr15581.incf"
  !-------------
  ! Modules used
  !-------------
  USE mo_kind,       ONLY: sp, wp           ! working precision
  USE mo_algorithms, ONLY: hunt             ! index search routine
  USE mo_exception,  ONLY: finish           ! abort routine
  USE mo_memory,     ONLY: t_m,            &! data type holding 3D fields
                           destruct,       &! destruct t_m
                           allocate,       &!
                           deallocate,     &!
                           update,         &!
                           lgtd, lgti,     &! Legendre transformation, Inverse
                           ASSIGNMENT(=),  &! assign 3D fields
                           print            ! i/o routines
  USE mo_grads_atm,  ONLY: to_grads         ! write 3D fields to grads file
  USE mo_transform,  ONLY: gauaw            ! calculate Gaussian latitudes
  USE mo_grads,      ONLY: t_ctl,          &! data type holding grads .ctl info
                           init_ctl,       &! initialize t_ctl data type
                           destruct         ! destruct   t_ctl data type
  USE mo_mpi_dace,   ONLY: dace,           &! MPI communication info
                           p_max,          &! maximum value over PEs
                           p_min,          &! minimum value over PEs
                           p_sum,          &! sum     value over PEs
                           p_bcast,        &! overloaded broadcast routine
                           p_bcast_ptr      ! broadcast pointer + content
  USE mo_fortran_units,&
                     ONLY: get_unit_number,&!
                           return_unit_number
  USE mo_physics,    ONLY: gacc,           &! gravity acceleration
                           rearth,         &! radius of the earth
                           pi,             &! pi
                           d2r,            &! factor degree->radians (pi/180)
                           r2d,            &! factor radians->degree (180/pi)
                           R                ! gas constant of dry air[J/(kg*K)]
  USE mo_wmo_tables, ONLY: WMO6_LATLON,     &
                           WMO6_ROTLL,      &
                           WMO6_GAUSSIAN,   &
                           DWD6_ICOSAHEDRON,&
                           DWD6_ICON,       &
                           DWD6_NONE,       &
                           WMO6_HARMONIC,   &
                           WMO3_HYBRID,     &
                           WMO3_HYBRIDB,    &
                           WMO3_HHYBRID,    &
                           WMO3_ISOBARIC,   &
!                          WMO3_SEALEVEL,   &
!                          WMO3_SURFACE,    &
                           WMO3_GENV,       &
                           WMO8_J_POSITIVE
  USE mo_ico_grid,   ONLY: factorize_nir,  &!
                           global_coordinates
  USE mo_atm_decomp, ONLY: t_atm_dec,      &! parallel decomposition data type
                           setup_single,   &! setup single processor run
                           setup_parallel, &! setup multi  processor run
                           setup_decomp,   &! setup domain decomposition
                           setup_com_tri,  &! setup communication buffers
                           setup_com_reg,  &! setup communication buffers
                           setup_com_icon, &! dto. for unstructured grid
                           destruct         ! clean up t_atm_dec
  USE mo_atm_transp, ONLY: scatter_level,  &! scatter one horizontal grid level
                           gather_multi     ! gather  multi level field
  USE mo_run_params, ONLY: data,           &! path to constant data
                           path_file,      &! concatenate path/filename
                           cosmo_refatm,   &! fallback values for COSMO
                           cosmo_ivctype,  &! .. reference atmosphere and vct
                           cosmo_svc,      &! large/small-scale decay rate
                           cosmo_nfltvc,   &! number of filter applications
                           l_ke_in_gds,    &! explicit GDS entry for number of model levels
                           nproc1,         &! partitioning parameter
                           nproc2,         &! partitioning parameter
                           geoid_format,   &! file format of geoid
                           geoid_file       ! geoid file name
  use mo_t_obs,      only: ptop_lapse,     &! pressure bounds for derivation ..
                           pbot_lapse       ! of lapse rate (t2m extrapolation)
  use mo_usstd,      only: p_h_usstd,      &! Pressure from geopotential height
                           h_p_usstd        ! geopotential height from pressure
  use mo_icon_grid,  only: t_grid_icon,    &! ICON grid metadata
                           init_icongrid,  &! Read gridfile, init metadata
                           finish_icongrid,&! Cleanup of ICON grid metadata
                           compress_icon_grid ! Compress redundant grid info
  use mo_dace_string,only: char3,          &! integer to character(len=3)
                           byte2hex,       &! CHAR to 'hexadecimal'
                           tolower          ! Lowercase of string
  use environment,   only: model_abort      ! COSMO style abort routine
  use utilities,     only: sleve_split_oro  ! decomposes topography scales
  use vgrid_refatm_utils,                  &!
                     only: vcoord_type,    &! COSMO derived types
                           refatm_type,    &! "
                           vcoord_defaults,&! default vertical coordinates
                           refatm_defaults,&! default reference atmospheres
                           set_vcoord_defaults,     &!
                           dealloc_vcoord_structure,&!
                           set_refatm_defaults,     &!
                           dealloc_refatm_structure,&!
                           reference_atmosphere,    &! routines to calculate
                           reference_atmosphere_2    ! .. reference atmospheres

  !---------------------
  ! netCDF f90 interface
  !---------------------
  use netcdf,        only: nf90_Inquire_Attribute,&!
                           nf90_Inquire_Dimension,&!
                           nf90_inq_dimid,        &!
                           nf90_inq_varid,        &!
                           nf90_get_var,          &!
                           nf90_get_att,          &!
                           nf90_open,             &!
                           nf90_close ,           &!
                           NF90_NoWrite,          &!
                           NF90_GLOBAL

  implicit none
!==============================================================================
  !----------------
  ! Public entities
  !----------------
  PRIVATE
  !-----------------------------------------
  ! derived type definition and constructors
  !-----------------------------------------
  PUBLIC :: maxab                ! maximum number of vertical levels in t_grid
  PUBLIC :: maxls                ! maximum number of soil     layers in t_grid
  PUBLIC :: t_grid               ! Type definition for grid information
  PUBLIC :: t_vcoord             ! data type for COSMO vertical coordinates
  PUBLIC :: t_refatm             ! data type for COSMO reference atmosphere
  PUBLIC :: construct, destruct  ! Constructor and destructor routine
  PUBLIC :: allocate             ! Allocate t_grid pointer components
  PUBLIC :: deallocate           ! Deallocate t_grid pointer components
  PUBLIC :: print                ! Printout of grid information
  PUBLIC :: same_horizontal_grid ! .true. for same horizontal grid
  PUBLIC :: same_vertical_grid   ! .true. for same vertical grid
  !---------------------------------------------------
  ! DACE convention for vertical coordinate parameters
  !---------------------------------------------------
  PUBLIC :: VCT_NO               ! no vertical grid, unknown vct
  PUBLIC :: VCT_P_ISO            ! isobaric coordinate
  PUBLIC :: VCT_P_HYB            ! hybrid pressure coordinate      (GME)
  PUBLIC :: VCT_P_IFS            ! IFS hybrid pressure coordinate  (IFS)
  PUBLIC :: VCT_Z_HYB            ! hybrid z coordinate           (COSMO)
  PUBLIC :: VCT_Z_GEN            ! generalised z coordinate (ICON,COSMO)
  PUBLIC :: IVCTYPE_GME          ! (Internal) GME/HRM Vertical coordinate
  PUBLIC :: IVCTYPE_ICON         ! (Internal) Vertical coordinate type for ICON
  !------------------------------
  ! DACE convention for model Id.
  !------------------------------
  PUBLIC :: MO_UNKNOWN           ! unknown
  PUBLIC :: MO_HRM               ! HRM
  PUBLIC :: MO_GME               ! GME
  PUBLIC :: MO_COSMO             ! COSMO
  PUBLIC :: MO_ICON              ! ICON
  PUBLIC :: MO_IFS               ! IFS
  PUBLIC :: cmodel               ! return model in character representation
  !----------------------------------
  ! set specific components of t_grid
  !----------------------------------
  PUBLIC :: set_pointers         ! update pointers with memory table
  PUBLIC :: set_ptopf            ! set or estimate topmost full level pressure
  PUBLIC :: set_plev_indices     ! set level indices for lapse rate estimate
  PUBLIC :: set_zlev_ref         ! estimate unperturbed model level heights
  PUBLIC :: set_geoid            ! set geoid heights
  PUBLIC :: set_empty_grid       ! set default values for empty grid
  PUBLIC :: set_grid_fields      ! mnemonics and bounds for specific fields
  PUBLIC :: setup_global_coord   ! set cartesian coordinates of the mass points
  !--------------------------------------
  ! routines inherited/adapted from COSMO
  !--------------------------------------
  PUBLIC :: cosmo_ref_atm        ! set fields of COSMO reference atmosphere
  PUBLIC :: phirot2phi           ! convert latitude  from the rotated system
  PUBLIC :: phi2phirot           ! convert latitude  in   the rotated system
  PUBLIC :: rlarot2rla           ! convert longitude from the rotated system
  PUBLIC :: rla2rlarot           ! convert longitude in   the rotated system
  PUBLIC :: p_bcast              ! broadcast t_vcoord, t_refatm
  PUBLIC :: get_vertcoord_new    ! read COSMO vertical coordinate parameters
  PUBLIC :: cosmo_vcp            ! compute GDS vertical coordinate parameters
  !---------------------------
  ! routines adapted from ICON
  !---------------------------
  PUBLIC :: read_icon_metadata   ! Read gridfile, setup ICON grid metadata
  !-------------------------------------------
  ! write GRADS data sets (plot regular grids)
  !-------------------------------------------
  PUBLIC :: to_grads             ! dump t_grid to grads file
  PUBLIC :: init_ctl             ! initialise grads .ctl file from grid
  !----------------
  ! transformations
  !----------------
  PUBLIC :: lin_intpol           ! Linear interpolation of 2D field
  PUBLIC :: average              ! average to coarser grid
  PUBLIC :: lgtd                 ! perform direct Legendre transf.(on geosp)
  PUBLIC :: lgti                 ! perform inverse Legendre transf.
  PUBLIC :: coarsen              ! select strides of gridpoints
  !-------------------------------------------------------
  ! return nearest gridpoint or interpolation coefficients
  !-------------------------------------------------------
  PUBLIC :: hunt                 ! obtain indices
  PUBLIC :: w_intpol             ! obtain indices, weights for interpolation
  !--------------------------
  ! utility routine for LETKF
  !--------------------------
  PUBLIC :: reduced_grid         ! derive coarse grid
  PUBLIC :: fit_rotll            ! fit rotated lat-lon grid into given grid
  PUBLIC :: construct_local      ! construct local (rotated) lat-lon grid

!==============================================================================
  !-----------
  ! Interfaces
  !-----------
  interface p_bcast
    module procedure bcast_refatm   ! broadcast derived type t_refatm
    module procedure bcast_vcoord   ! broadcast derived type t_vcoord
  end interface p_bcast
!==============================================================================
  !------------------
  ! Public parameters
  !------------------
  integer, parameter :: IVCTYPE_GME  = 0 ! GME/HRM hybrid pressure coordinates
  integer, parameter :: IVCTYPE_ICON = 4 ! Vertical coordinate type for ICON
  !                     ivctype 1 and 4    only used internally in the DA code
  !                     ivctype 1 to  4    are used by COSMO

  integer, parameter :: VCT_NO     = 0  ! no vertical grid, unknown vct
  integer, parameter :: VCT_P_ISO  = 1  ! isobaric coordinate
  integer, parameter :: VCT_P_HYB  = 2  ! hybrid pressure coordinate      (GME)
  integer, parameter :: VCT_Z_HYB  = 3  ! hybrid z coordinate           (COSMO)
  integer, parameter :: VCT_Z_GEN  = 4  ! generalised z coordinate (ICON,COSMO)
  integer, parameter :: VCT_P_IFS  = 98 ! IFS hybrid pressure coordinate  (IFS)

  integer, parameter :: MO_UNKNOWN = 0  ! model:
  integer, parameter :: MO_HRM     = 1  ! HRM
  integer, parameter :: MO_GME     = 2  ! GME
  integer, parameter :: MO_COSMO   = 3  ! COSMO
  integer, parameter :: MO_ICON    = 4  ! ICON
  integer, parameter :: MO_IFS     = 5  ! IFS

  character(len=8), parameter :: cmodel (0:5) = ['UNKNOWN ', &! 0
                                                 'HRM     ', &! 1
                                                 'GME     ', &! 2
                                                 'COSMO   ', &! 3
                                                 'ICON    ', &! 4
                                                 'IFS     ']  ! 5
  !-----------------
  ! Local parameters
  !-----------------
  INTEGER ,PARAMETER :: maxab  = 250    ! maximum number of vertical levels
  INTEGER ,PARAMETER :: maxls  =   8    ! maximum number of soil     levels
!==============================================================================
  !======================
  ! Grid type definitions
  !======================
  !------------------------------------------------------------
  ! type t_vcoord and t_refatm reflect the COSMO derived types
  ! vcoord_type and refatm_type (defined in vgrid_refatm_utils)
  ! without unecessary pointer components
  !------------------------------------------------------------
  type t_vcoord
    !-----------------------------------------
    ! data type for COSMO vertical coordinates
    !-----------------------------------------
    integer           :: ivctype      =  0     ! vertical coordinate type
!   integer           :: ivcoord_id
    integer           :: nlevels      = -1     ! number of half levels (ke+1)
!   integer           :: kflat                 ! index where levels become flat
    character         :: vc_uuid(16)  =achar(0)! UUID of HHL file
    real(wp)          :: vcflat       = 0._wp  ! coord.where levels become flat
    real(wp) ,pointer :: vert_coord(:)=>NULL() ! height above sea level
    real(wp) ,pointer :: sigm_coord(:)=>NULL() ! reference pressure normalized
    !------------------
    ! sleve parameters:
    !------------------
    real(wp)          :: svc1         = 0._wp  ! decay rate large-scale topo
    real(wp)          :: svc2         = 0._wp  ! decay rate small-scale topo
    integer           :: nfltvc       = 0      ! number of filter applications
  end type t_vcoord

  type t_refatm
    !--------------------------------------------------
    ! data type for COSMO reference atmosphere settings
    !--------------------------------------------------
    INTEGER   :: irefatm    = -1
!   INTEGER   :: irefatm_id
    REAL(wp)  :: p0sl       = 0._wp
    REAL(wp)  :: t0sl       = 0._wp
    REAL(wp)  :: dt0lp      = 0._wp
    REAL(wp)  :: delta_t    = 0._wp
    REAL(wp)  :: h_scal     = 0._wp
!   REAL(wp)  :: bvref
  end type t_refatm

  TYPE t_grid
    !-----------------------
    ! parallel decomposition
    !-----------------------
    LOGICAL           :: ldec           ! true for decomposed domain
    TYPE (t_atm_dec)  :: dc             ! decomposition information
    !--------------------
    ! general information
    !--------------------
    REAL(wp)          :: a              ! radius of the planet
    REAL(wp)          :: g              ! gravity acceleration
    INTEGER           :: model          ! model
    INTEGER           :: gridtype       ! GDS octet 6 code values
    INTEGER           :: nx             ! number of grid points in x
    INTEGER           :: ny             ! number of grid points in y
    INTEGER           :: nz             ! number of levels
    INTEGER           :: nd             ! number of diamonds
    INTEGER           :: ns             ! number of soil layers
    INTEGER           :: nxny           ! nx * ny * nd
    INTEGER           :: size           ! nx * ny * nz * nd
    INTEGER           :: lbg   (4)      ! lower bound of global arrays
    INTEGER           :: ubg   (4)      ! upper bound of global arrays
    INTEGER           :: lb    (4)      ! lower bound of local  arrays
    INTEGER           :: ub    (4)      ! upper bound of local  arrays
    INTEGER           :: shape (4)      ! shape       of local  arrays
    INTEGER           :: scanmode       ! Orientation of quasi-regular grids
    INTEGER           :: d_gme (0:10)   ! ICON offset for GME diamond
    !--------------
    ! Gaussian grid
    !--------------
    INTEGER           :: ngl            ! number of gaussian latitudes
    !------------------------
    ! Spectral representation
    !------------------------
    INTEGER           :: nn             ! truncation
    !----------------
    ! Triangular grid
    !----------------
    INTEGER           :: ni
    INTEGER           :: ni2
    INTEGER           :: nir            ! root
    !--------------
    ! Regular grids
    !--------------
    REAL(wp)          :: la1            ! latitude of first grid point
    REAL(wp)          :: lo1            ! longitude of first grid point
    REAL(wp)          :: di             ! longitudinal increment
    REAL(wp)          :: dj             ! latitudinal increment
    REAL(wp)          :: d_deg          ! nominal resolution in degree
    REAL(wp)          :: d_km           ! nominal resolution in km
    LOGICAL           :: rot            ! rotated
    LOGICAL           :: cyc_x          ! cyclic boundary conditions in x
    LOGICAL           :: poly           ! poles at y boundaries
    LOGICAL           :: global         ! global grid
    REAL(wp) ,pointer :: dlon   (:)     ! (rotated) longitude [degree]
    REAL(wp) ,pointer :: dlat   (:)     ! (rotated) latitude  [degree]
    REAL(wp)          :: dlonr          ! longitude of N.pole of rotation
    REAL(wp)          :: dlatr          ! latitude  of N.pole of rotation
    CHARACTER(len=1)  :: arakawa        ! arakawa type
    !------------------------
    ! vertical discretization
    !------------------------
    INTEGER           :: levtyp         ! level type
    real(wp)          :: ptopf          ! model top pressure (full levels) [Pa]
    real(wp)          :: htopf          ! model top          (full levels) [gpm]
    integer           :: kmlbot         ! lower layer for lapse rate estimate
    integer           :: kmltop         ! upper layer for lapse rate estimate
    integer           :: vct            ! DACE unique vertical coordinate type
    !-----------------------------
    ! leveltype hybrid or isobaric
    !-----------------------------
    REAL(wp)          :: ak (maxab)     ! vertical coordinate parameter
    REAL(wp)          :: bk (maxab)     ! vertical coordinate parameter
    REAL(wp)          :: akf(maxab)     ! vertical coordinate p. (full levels)
    REAL(wp)          :: bkf(maxab)     ! vertical coordinate p. (full levels)
    !------------------------------
    ! COSMO vertical discretization
    !------------------------------
    type (t_vcoord)   :: vc             ! specification of vertical coordinate
    type (t_refatm)   :: ra             ! specification of reference atmosphere
    !-----------------------------
    ! soil vertical discretization
    !-----------------------------
    integer           :: dbs(maxls)     ! depth below surface [cm] (soil)
    !-------------------
    ! ICON grid metadata
    !-------------------
    type(t_grid_icon), pointer :: icongrid
    !------------------
    ! 2d surface fields
    !------------------
    TYPE(t_m)         :: m        (22)      ! 2D / 3D fields
    !-------------------------------------------
    ! pointers associated with 'grid% m(:)% ptr'
    !-------------------------------------------
    REAL(wp) _POINTER :: lsm      (:,:,:,:) ! land sea mask
    REAL(wp) _POINTER :: geosp    (:,:,:,:) ! surface geopotential       [m2/s2]
    REAL(wp) _POINTER :: rlon     (:,:,:,:) ! longitude                [radians]
    REAL(wp) _POINTER :: rlat     (:,:,:,:) ! latitude                 [radians]
    REAL(wp) _POINTER :: geo_sh   (:)       ! surface geopotential (spectral)
    REAL(wp) _POINTER :: geoid    (:,:,:,:) ! geoid heights                  [m]
    REAL(wp) _POINTER :: soiltyp  (:,:,:,:) ! soil type                   [1..9]
    REAL(wp) _POINTER :: fr_lake  (:,:,:,:) ! lake fraction
    REAL(wp) _POINTER :: depth_lk (:,:,:,:) ! lake depth                     [m]
    REAL(wp) _POINTER :: xnglob   (:,:,:,:) ! cartesian coord. on unit sphere
                                            ! fields for ICON/COSMO follow:
    REAL(wp) _POINTER :: hsurf    (:,:,:,:) ! geometrical height of surface
    REAL(wp) _POINTER :: hhl      (:,:,:,:) ! geometrical height of half levels
    REAL(wp) _POINTER :: p0       (:,:,:,:) ! base-state pressure
    REAL(wp) _POINTER :: dp0      (:,:,:,:) ! base-state pressure thickness
    REAL(wp) _POINTER :: rho0     (:,:,:,:) ! base-state density
    REAL(wp) _POINTER :: sso_stdh (:,:,:,:) ! SSO standard deviation         [m]
    REAL(wp) _POINTER :: sso_gamma(:,:,:,:) ! SSO
    REAL(wp) _POINTER :: sso_theta(:,:,:,:) ! SSO
    REAL(wp) _POINTER :: sso_sigma(:,:,:,:) ! SSO
    REAL(wp) _POINTER :: emis_rad (:,:,:,:) ! surface emissivity
    REAL(wp) _POINTER :: prs_min  (:,:,:,:) ! min. stomata resistance
    REAL(wp) _POINTER :: skc      (:,:,:,:) ! skin conductivity         [W/m2/K]
    !----------------------------------------------
    ! reference to 'original' indices and processor
    !----------------------------------------------
    INTEGER  _POINTER :: marr     (:,:,:,:) ! (4,i,j,d) index field: (pe,i,j,d)
    !-----------------------------------------
    ! Dual grid: decomposition, auxiliary data
    !-----------------------------------------
    TYPE (t_atm_dec)  :: dc_d           ! decomposition information
    INTEGER           :: nx_d           ! number of dual grid points (x)
    INTEGER           :: ny_d           ! number of dual grid points (y)
    INTEGER           :: nxny_d         ! number of dual grid points (total)
    INTEGER           :: lbg_d (4)      ! lower bound of global arrays
    INTEGER           :: ubg_d (4)      ! upper bound of global arrays
    INTEGER           :: lb_d  (4)      ! lower bound of local  arrays
    INTEGER           :: ub_d  (4)      ! upper bound of local  arrays
    INTEGER  _POINTER :: pe_d  (:,:,:)  ! pe holding dual grid point (i,j,d)
  END TYPE t_grid
!==============================================================================
enum, bind(c)
  enumerator :: LSM = 1, GEOSP, RLON, RLAT, GEO_S, GEOID, SOILTYP,           &
                FR_LAKE, DEPTH_LK, XNGLOB, HSURF, HHL, P0, DP0, RHO0,        &
                SSO_STDH, SSO_GAMMA, SSO_THETA, SSO_SIGMA, EMIS_RAD, PRS_MIN,&
                SKC
end enum
!==============================================================================
  !-----------
  ! Interfaces
  !-----------
  INTERFACE construct
    !------------------------------
    ! set up model grid description
    !------------------------------
    MODULE PROCEDURE construct_atm_grid   ! general case
    MODULE PROCEDURE construct_atm_grid_1 ! lon, lat provided as 1-d arrays
    MODULE PROCEDURE construct_atm_grid_2 ! lon, lat provided as 2-d arrays
  END INTERFACE construct
!------------------------------------------------------------------------------
  INTERFACE destruct
    MODULE PROCEDURE destruct_atm_grid
    MODULE PROCEDURE destruct_atm_grids
    MODULE PROCEDURE destruct_vcoord
  END INTERFACE destruct
!------------------------------------------------------------------------------
  INTERFACE PRINT
    MODULE PROCEDURE print_atm_grid
  END INTERFACE PRINT
!------------------------------------------------------------------------------
  INTERFACE ALLOCATE
    MODULE PROCEDURE allocate_atm_grid
  END INTERFACE ALLOCATE
!------------------------------------------------------------------------------
  INTERFACE DEALLOCATE
    MODULE PROCEDURE deallocate_atm_grid
  END INTERFACE DEALLOCATE
!------------------------------------------------------------------------------
  INTERFACE set_pointers
    MODULE PROCEDURE set_grid_pointers
  END INTERFACE set_pointers
!------------------------------------------------------------------------------
  INTERFACE lgtd
    MODULE PROCEDURE lgt_grid
  END INTERFACE lgtd
!------------------------------------------------------------------------------
  INTERFACE lgti
    MODULE PROCEDURE lgti_grid
  END INTERFACE lgti
!------------------------------------------------------------------------------
  INTERFACE average
    MODULE PROCEDURE average_atm
  END INTERFACE average
!------------------------------------------------------------------------------
  INTERFACE init_ctl
    MODULE PROCEDURE init_ctl_from_grid
  END INTERFACE init_ctl
!------------------------------------------------------------------------------
  INTERFACE to_grads
    MODULE PROCEDURE grid_to_grads
  END INTERFACE to_grads
!------------------------------------------------------------------------------
  INTERFACE hunt
    MODULE PROCEDURE hunt_grid
  END INTERFACE hunt
!------------------------------------------------------------------------------
  INTERFACE w_intpol
    MODULE PROCEDURE w_intpol
  END INTERFACE w_intpol
!------------------------------------------------------------------------------
  INTERFACE coarsen
    MODULE PROCEDURE coarsen_grid
  END INTERFACE coarsen
!==============================================================================
CONTAINS
!==============================================================================
!------------------------------------------------------------------------------
  SUBROUTINE construct_atm_grid_1 (grid, TEMPLATE,                    &
                                         gridtype,                    &
                                         nx, ny, ke, ngl, nn, ni, ns, &
                                         dlon, dlat,                  &
                                         lon, lat,                    &
                                         lo1, la1, di, dj,            &
                                         lor, lar,                    &
                                         levtyp,                      &
                                         ak, bk)
  !-------------------------------------------------------------------
  ! Initialize a grid.
  !
  ! lon, lat are provided as 1-dimensional arrays each.
  ! The precedence of the information in the optional parameters is as
  ! follows:
  !
  !  1) an empty grid is constructed.
  !  2) template:
  !     If a template grid is specified, an identical grid is
  !     constructed as long as no further parameters are given to
  !     overwrite this information:
  !  3) gridtype:
  !     GDS grid type
  !  4) nx, ny, ke:
  !     Number of longitudes, latitudes and levels
  !  5) dlon, dlat:
  !     longitudes and latitudes
  !  6) ak, bk:
  !     Sigma coordinate coefficients
  !-------------------------------------------------------------------
  TYPE (t_grid) ,INTENT(out)           :: grid
  TYPE (t_grid) ,INTENT(in)  ,OPTIONAL :: TEMPLATE
  INTEGER       ,INTENT(in)  ,OPTIONAL :: gridtype
  INTEGER       ,INTENT(in)  ,OPTIONAL :: nx, ny, ke, ngl, nn, ni, ns
  REAL(wp)      ,INTENT(in)  ,OPTIONAL :: dlon(:), dlat(:)
  REAL(wp)      ,INTENT(in)            :: lon(:), lat(:)
  REAL(wp)      ,INTENT(in)  ,OPTIONAL :: lo1, la1, di, dj
  REAL(wp)      ,INTENT(in)  ,OPTIONAL :: lor, lar
  INTEGER       ,INTENT(in)  ,OPTIONAL :: levtyp
  REAL(wp)      ,INTENT(in)  ,OPTIONAL :: ak(:), bk(:)
    !----------------
    ! local variables
    !----------------
    REAL(wp) :: lon2 (SIZE(lon),1,1)
    REAL(wp) :: lat2 (SIZE(lat),1,1)
    !-------------------------------
    ! wrapper for construct_atm_grid
    !-------------------------------
    lon2(:,1,1) = lon
    lat2(:,1,1) = lat
    CALL construct_atm_grid (grid, TEMPLATE,                    &
                                   gridtype,                    &
                                   nx, ny, ke, ngl, nn, ni, ns, &
                                   dlon, dlat,                  &
                                   lon2, lat2,                  &
                                   lo1, la1, di, dj,            &
                                   lor, lar,                    &
                          levtyp = levtyp,                      &
                              ak = ak,                          &
                              bk = bk                           )
  END SUBROUTINE construct_atm_grid_1
!------------------------------------------------------------------------------
  SUBROUTINE construct_atm_grid_2 (grid, TEMPLATE,                    &
                                         gridtype,                    &
                                         nx, ny, ke, ngl, nn, ni, ns, &
                                         dlon, dlat,                  &
                                         lon, lat,                    &
                                         lo1, la1, di, dj,            &
                                         lor, lar,                    &
                                         levtyp,                      &
                                         ak, bk)
  !-------------------------------------------------------------------
  ! Initialize a grid.
  !
  ! lon, lat are provided as 2-dimensional arrays each.
  ! The precedence of the information in the optional parameters is as
  ! follows:
  !
  !  1) an empty grid is constructed.
  !  2) template:
  !     If a template grid is specified, an identical grid is
  !     constructed as long as no further parameters are given to
  !     overwrite this information:
  !  3) gridtype:
  !     GDS grid type
  !  4) nx, ny, ke:
  !     Number of longitudes, latitudes and levels
  !  5) dlon, dlat:
  !     longitudes and latitudes
  !  6) ak, bk:
  !     Sigma coordinate coefficients
  !-------------------------------------------------------------------
  TYPE (t_grid) ,INTENT(out)           :: grid
  TYPE (t_grid) ,INTENT(in)  ,OPTIONAL :: TEMPLATE
  INTEGER       ,INTENT(in)  ,OPTIONAL :: gridtype
  INTEGER       ,INTENT(in)  ,OPTIONAL :: nx, ny, ke, ngl, nn, ni, ns
  REAL(wp)      ,INTENT(in)  ,OPTIONAL :: dlon(:), dlat(:)
  REAL(wp)      ,INTENT(in)            :: lon(:,:), lat(:,:)
  REAL(wp)      ,INTENT(in)  ,OPTIONAL :: lo1, la1, di, dj
  REAL(wp)      ,INTENT(in)  ,OPTIONAL :: lor, lar
  INTEGER       ,INTENT(in)  ,OPTIONAL :: levtyp
  REAL(wp)      ,INTENT(in)  ,OPTIONAL :: ak(:), bk(:)
    !----------------
    ! local variables
    !----------------
    REAL(wp) :: lon2 (SIZE(lon,1),SIZE(lon,2),1)
    REAL(wp) :: lat2 (SIZE(lat,1),SIZE(lat,2),1)
    !-------------------------------
    ! wrapper for construct_atm_grid
    !-------------------------------
    lon2(:,:,1) = lon
    lat2(:,:,1) = lat
    CALL construct_atm_grid (grid, TEMPLATE,                    &
                                   gridtype,                    &
                                   nx, ny, ke, ngl, nn, ni, ns, &
                                   dlon, dlat,                  &
                                   lon2, lat2,                  &
                                   lo1, la1, di, dj,            &
                                   lor, lar,                    &
                          levtyp = levtyp,                      &
                              ak = ak,                          &
                              bk = bk                           )
  END SUBROUTINE construct_atm_grid_2
!------------------------------------------------------------------------------
  SUBROUTINE construct_atm_grid (grid, TEMPLATE,                    &
                                       gridtype,                    &
                                       nx, ny, ke, ngl, nn, ni, ns, &
                                       dlon, dlat,                  &
                                       lon, lat,                    &
                                       lo1, la1, di, dj,            &
                                       lor, lar,                    &
                                       levtyp,                      &
                                       ivctype,                     &! COSMO
                                       vcoord,                      &! COSMO
                                       refatm,                      &! COSMO
                                       dbs,                         &
                                       ak, bk,                      &
                                       akf, bkf,                    &
                                       nproc1, nproc2, comm,        &
                                       arakawa,                     &
                                       scanmode,                    &
                                       uuid_v,                      &! ICON
                                       icongrid)
  !-------------------------------------------------------------------
  ! Initialize a grid.
  !
  ! The precedence of the information in the optional parameters is as
  ! follows:
  !
  !  1) an empty grid is constructed.
  !  2) template:
  !     If a template grid is specified, an identical grid is
  !     constructed as long as no further parameters are given to
  !     overwrite this information:
  !  3) gridtype:
  !     GDS grid type
  !  4) nx, ny, ke:
  !     Number of longitudes, latitudes and levels
  !  5) dlon, dlat:
  !     longitudes and latitudes
  !  6) ak, bk:
  !     Sigma coordinate coefficients
  !-------------------------------------------------------------------
  TYPE(t_grid)  ,INTENT(out)          :: grid            ! grid to set
  TYPE(t_grid)  ,INTENT(in) ,OPTIONAL :: template        ! template grid
  INTEGER       ,INTENT(in) ,OPTIONAL :: gridtype        ! horizontal gridtype
  INTEGER       ,INTENT(in) ,OPTIONAL :: nx, ny, ke      ! grid sizes: x,y,z
  INTEGER       ,INTENT(in) ,OPTIONAL :: ngl             ! no.Gaussian latitudes
  INTEGER       ,INTENT(in) ,OPTIONAL :: nn              ! spectral truncation
  INTEGER       ,INTENT(in) ,OPTIONAL :: ni              ! GME resolution
  INTEGER       ,INTENT(in) ,OPTIONAL :: ns              ! no.soil layers
  REAL(wp)      ,INTENT(in) ,OPTIONAL :: dlon(:), dlat(:)! lat/lon (degree)
  REAL(wp)      ,INTENT(in) ,OPTIONAL :: lon(:,:,:)      ! longitude  (rad)
  REAL(wp)      ,INTENT(in) ,OPTIONAL :: lat(:,:,:)      ! latitude   (rad)
  REAL(wp)      ,INTENT(in) ,OPTIONAL :: lo1, la1, di, dj! lat/lon start, delta
  REAL(wp)      ,INTENT(in) ,OPTIONAL :: lor, lar        ! pole of rotation
  INTEGER       ,INTENT(in) ,OPTIONAL :: levtyp          ! level type
  INTEGER       ,INTENT(in) ,OPTIONAL :: ivctype         ! COSMO
  TYPE(t_vcoord),INTENT(in) ,OPTIONAL :: vcoord          ! COSMO vert.coord.params.
  TYPE(t_refatm),INTENT(in) ,OPTIONAL :: refatm          ! COSMO reference atmosphere
  REAL(wp)      ,INTENT(in) ,OPTIONAL :: dbs(:)          ! depth below surface
  REAL(wp)      ,INTENT(in) ,OPTIONAL :: ak(:),  bk(:)   ! vert.coord.parameters
  REAL(wp)      ,INTENT(in) ,OPTIONAL :: akf(:), bkf(:)  ! v.c.p full levels
  INTEGER       ,INTENT(in) ,OPTIONAL :: nproc1, nproc2  ! no.processors/lon,lat
  INTEGER       ,INTENT(in) ,OPTIONAL :: comm            ! MPIcommunicator group
  CHARACTER     ,INTENT(in) ,OPTIONAL :: arakawa         ! grid type ('A','C')
  INTEGER       ,INTENT(in) ,OPTIONAL :: scanmode        ! gridpoint orientation
  CHARACTER     ,INTENT(in) ,OPTIONAL :: uuid_v(:)       ! uuidOfVGrid
  TYPE(t_grid_icon),POINTER ,OPTIONAL :: icongrid        ! ICON grid definition
    !----------------
    ! Local variables
    !----------------
    INTEGER               :: i,j
    INTEGER               :: ierr
    INTEGER               :: nz
    REAL(wp)              :: sinlat, coslat, sinlon, coslon
    TYPE (t_atm_dec)      :: dc ! temporary decomposition info
    logical               :: lico
    integer      ,pointer :: marr_(:,:,:,:) ! temporary index field

    !============================
    ! Set defaults for empty grid
    ! Nullify pointer components
    !============================
    call set_empty_grid (grid)
    !=====================================
    ! Specifications by optional arguments
    !=====================================
    CALL setup_single (dc)
    IF (PRESENT(TEMPLATE)) then
      grid = TEMPLATE
      dc   = TEMPLATE% dc
      if (associated (template% vc% vert_coord)) then
        allocate (grid% vc% vert_coord (size (template% vc% vert_coord)))
                  grid% vc% vert_coord =      template% vc% vert_coord
      endif
      if (associated (template% vc% sigm_coord)) then
        allocate (grid% vc% sigm_coord (size (template% vc% sigm_coord)))
                  grid% vc% sigm_coord =      template% vc% sigm_coord
      endif
      !------------------------------
      ! Nullify pointer components
      ! Don't associate with template
      !------------------------------
      nullify (grid% dlon)
      nullify (grid% dlat)
      nullify (grid% marr)
      nullify (grid% pe_d)
      DO i=1,SIZE(grid% m)
         NULLIFY (grid% m(i)% ptr)      ! With TR15581, ptr is ALLOCATABLE
      END DO
    endif
    IF (PRESENT(dlon))        grid% nx          = SIZE(dlon)
    IF (PRESENT(dlat))        grid% ny          = SIZE(dlat)
    IF (PRESENT(lon))         grid% nx          = SIZE(lon,1)
    IF (PRESENT(lat))         grid% ny          = SIZE(lat,2)
    IF (PRESENT(ak))          grid% nz          = SIZE(ak)-1
    IF (PRESENT(akf))         grid% nz          = SIZE(akf)
    IF (PRESENT(dbs))         grid% ns          = SIZE(dbs)
    IF (PRESENT(nx))          grid% nx          = nx
    IF (PRESENT(ny))          grid% ny          = ny
    IF (PRESENT(ngl))         grid% ngl         = ngl
    IF (PRESENT(ke))          grid% nz          = ke
    IF (PRESENT(nn))          grid% nn          = nn
    IF (PRESENT(ni))          grid% ni          = ni
    IF (PRESENT(ns))          grid% ns          = ns
    IF (PRESENT(di))          grid% di          = di
    IF (PRESENT(dj))          grid% dj          = dj
    IF (PRESENT(ivctype))     grid% vc% ivctype = ivctype
    IF (PRESENT(vcoord)) then
      grid% vc = vcoord
      if (associated (vcoord% vert_coord)) then
        allocate (grid% vc% vert_coord (size (vcoord% vert_coord)))
                  grid% vc% vert_coord =      vcoord% vert_coord
      endif
      if (associated (vcoord% sigm_coord)) then
        allocate (grid% vc% sigm_coord (size (vcoord% sigm_coord)))
                  grid% vc% sigm_coord =      vcoord% sigm_coord
      endif
    ENDIF
    IF (PRESENT(refatm))      grid% ra          = refatm
    IF (PRESENT(arakawa))     grid% arakawa     = arakawa
    IF (PRESENT(uuid_v))      grid% vc% vc_uuid = uuid_v
    IF (PRESENT(nproc1))      dc  % nproc1      = nproc1
    IF (PRESENT(nproc2))      dc  % nproc2      = nproc2
    IF (PRESENT(comm))        dc  % comm        = comm
    CALL setup_parallel (grid% dc, dc% nproc1, dc% nproc2, dc% comm)
    grid% ldec = (grid% dc% nproc1 * grid% dc% nproc2 /= 1)
    CALL setup_parallel (grid% dc_d, dc% nproc1, dc% nproc2, dc% comm)
    !-----------------------------
    ! check hard coded array sizes
    !-----------------------------
    if (grid% nz >= maxab) call finish ('construct_atm_grid','nz >= maxab')
    if (grid% ns >  maxls) call finish ('construct_atm_grid','ns >  maxls')
    !------------------
    ! estimate gridtype
    !------------------
    IF (PRESENT(nn)) then
       if (nn  > 0)        grid% gridtype    = WMO6_HARMONIC
    endif
    IF (PRESENT(ngl)) then
       if (ngl > 0)        grid% gridtype    = WMO6_GAUSSIAN
    end IF
    IF (PRESENT(ni)) then
       if (ni  > 0)        grid% gridtype    = DWD6_ICOSAHEDRON
    end IF
    IF (PRESENT(dj)) then
      select case(grid% gridtype)
      case(WMO6_LATLON, WMO6_ROTLL)
      case default
                           grid% gridtype    = WMO6_LATLON
      end select
    endif
    IF (PRESENT(gridtype)) grid% gridtype    = gridtype
    IF (PRESENT(levtyp)) then
                           grid% levtyp      = levtyp
      select case (levtyp)
      case (WMO3_ISOBARIC)
                           grid% vc% ivctype = 0
      case default
      end select
    end IF
    !-------------------------
    ! check for valid scanmode
    !-------------------------
    if (PRESENT(scanmode)) then
      select case (grid% gridtype)
      case (WMO6_GAUSSIAN)
        select case (scanmode)
        case (0,WMO8_J_POSITIVE)
          grid% scanmode = scanmode
        case default
          write(0,*)'construct_atm_grid, scanmode not implemented:',scanmode
          call finish('construct_atm_grid','scanmode not implemented')
        end select
      case default
        write(0,*)'construct_atm_grid, scanmode not implemented for gridtype',&
                   grid% gridtype
        call finish('construct_atm_grid',&
                   'scanmode not implemented for gridtype')
      end select
    endif
    !===============================
    ! calculate number of gridpoints
    !===============================
    SELECT CASE (grid% gridtype)
    !------------
    ! Latlon grid
    !------------
    CASE (WMO6_LATLON, WMO6_ROTLL)
      grid% lbg  = 1
      grid% ubg  = (/grid% nx, grid%  ny, grid% nz, grid% nd/)
      grid% lb   =   grid% lbg
      grid% ub   =   grid% ubg
      grid% nxny =   grid% nx * grid% ny * grid% nd

      CALL setup_decomp (grid% dc, grid% lb, grid% ub)

      CALL setup_com_reg (grid% dc,           &
                          grid% marr,         &
                          grid% lbg, grid% ubg)

    !------------------------------------
    ! Gaussian or spectral representation
    !------------------------------------
    CASE (WMO6_GAUSSIAN, WMO6_HARMONIC)
      IF (grid% ngl==0 .AND. grid% ny>0) grid% ngl = grid% ny
      !------------------------------------------------------
      ! estimate truncation from number of gaussian latitudes
      ! or vice versa
      !------------------------------------------------------
      IF (grid% nn == 0) THEN
        SELECT CASE (grid% ngl)
        CASE (512)
          grid% nn = 511 ! ECMWF linear grid
        CASE (480)
          grid% nn = 319
        CASE (320)
          grid% nn = 213
        CASE (160)
          grid% nn = 106
        CASE ( 96)
          grid% nn =  63
        CASE ( 94)
          grid% nn =  62
        CASE ( 64)
          grid% nn =  42
        CASE ( 48)
          grid% nn =  30
        CASE ( 32)
          grid% nn =  21
        CASE DEFAULT
          grid% nn = grid% ngl
        END SELECT
      ENDIF
      IF (grid% ngl==0) THEN
        SELECT CASE (grid% nn)
        CASE (511)
          grid% ngl = 512 ! ECMWF linear grid
        CASE (319)
          grid% ngl = 480
        CASE (213)
          grid% ngl = 320
        CASE (106)
          grid% ngl = 160
        CASE (63)
          grid% ngl =  96
        CASE (62)
          grid% ngl =  94
        CASE (42)
          grid% ngl =  64
        CASE (30)
          grid% ngl =  48
        CASE (21)
          grid% ngl =  32
        END SELECT
      ENDIF
      grid% ny = grid% ngl
      IF (grid% nx==0)                 grid% nx = grid% ny * 2
      grid% lbg  = 1
      grid% ubg  = (/grid% nx, grid%  ny, grid% nz, grid% nd/)
      grid% lb   =   grid% lbg
      grid% ub   =   grid% ubg
      grid% nxny =   grid% nx * grid% ny * grid% nd

      CALL setup_decomp (grid% dc, grid% lb, grid% ub)

      CALL setup_com_reg (grid% dc,           &
                          grid% marr,         &
                          grid% lbg, grid% ubg)

    !----------------
    ! triangular grid
    !----------------
    CASE (DWD6_ICOSAHEDRON)
      grid% nd     = 10
      grid% lbg    = (/      0 ,        1   ,       1 ,       1 /)
      grid% ubg    = (/grid% ni, grid%  ni+1, grid% nz, grid% nd/)
      grid% lb     =   grid% lbg
      grid% ub     =   grid% ubg
      grid% nx     =   grid% ni+1
      grid% ny     =   grid% ni+1
      grid% nxny   =   grid% nx * grid% ny * grid% nd
      grid% global = .true.

      CALL setup_decomp (grid% dc, grid% lb, grid% ub)

      CALL setup_com_tri (grid% dc,                             &
                          grid% marr,                           &
                          grid% lbg(1), grid% ubg(1),           &
                          grid% lbg(2), grid% ubg(2), grid% nd, &
                          grid% lb (1), grid% ub (1),           &
                          grid% lb (2), grid% ub (2), grid% ni  )

    !-----------------------
    ! ICON unstructured grid
    !-----------------------
    CASE (DWD6_ICON)
!print *, "construct_atm_grid(ICON): ni,nx,nz=",grid%ni,grid% nx,grid%nz
      if (grid%dc% nproc2 /= 1) then
         print *, "construct_atm_grid(ICON): nproc1,nproc2,npe =", &
              grid%dc% nproc1, grid%dc% nproc2, grid%dc% npe
         call finish ("construct_atm_grid","nproc2 /= 1 for ICON grid")
      end if
!     grid% nxny   =  20 * grid% ni ** 2    ! No nproma-type blocking (for now)
      grid% nxny   =   grid% nx
      grid% ny     =   grid% nx
      grid% nd     =   1
      if (grid% d_gme(1) == -1) then
        !---------------
        ! standard setup
        !---------------
        grid% lbg    = (/      1   ,      1   ,       1 ,       1 /)
        grid% ubg    = (/grid% nxny,      1   , grid% nz,       1 /)
        grid% lb     =   grid% lbg
        grid% ub     =   grid% ubg
        grid% global = .true.

        CALL setup_decomp (grid% dc, grid% lb, grid% ub)

        CALL setup_com_icon (grid% dc,           &
                             grid% marr,         &
                             grid% lbg, grid% ubg)
      else
        !--------------------------------------------
        ! non-standard setup, keep data from template
        !--------------------------------------------
        grid% ub (3) = grid% nz
        grid% ubg(3) = grid% nz
        allocate (grid% marr (4, grid% nxny, 1,1))
        grid% marr = template% marr
      endif
    !------------------------
    ! just a bunch of columns
    !------------------------
    CASE (DWD6_NONE)
      grid% nxny   = grid% nx * grid% ny * grid% nd
      grid% lbg    = 1
      grid% ubg    = (/grid% nx, grid%  ny, grid% nz, grid% nd/)
      grid% lb     = grid% lbg
      grid% ub     = grid% ubg
    !-------------
    ! unknown grid
    !-------------
    CASE default
      CALL finish ('construct_atm_grid','unknown grid type')
    END SELECT
    grid% size  = grid% nxny * max (1, grid% nz)
    grid% shape = grid% ub - grid% lb + 1
    !======================================================
    ! set mnemonics and bounds for the grid specific fields
    !======================================================
    call set_grid_fields (grid)

    IF (grid% nxny > 0) THEN
      !================
      ! allocate arrays
      !================
      ALLOCATE (grid% dlon (grid% nx)) ;grid% dlon = 0._wp
      ALLOCATE (grid% dlat (grid% ny)) ;grid% dlat = 0._wp
      CALL update (grid%m)
      CALL set_pointers (grid)
      !================
      ! set coordinates
      !================
      SELECT CASE (grid% gridtype)
      !------------------------
      ! just a bunch of columns
      !------------------------
      CASE (DWD6_NONE)
        if (present (lon) .and. present (lat)) then
          if (size (lon) /= size (lat)) &
             call finish ('construct_atm_grid','sizes of lon,lat must match')
          call setup_global_coord (grid% xnglob, grid% rlon, grid% rlat, &
                                   lon, lat, grid% lbg, grid% ubg, ierr  )
        endif
      !----------------
      ! triangular grid
      !----------------
      CASE (DWD6_ICOSAHEDRON)
        !----------
        ! factorize
        !----------
        IF (grid% ni < 4) CALL finish ('construct_atm_grid','ni must be >= 4')
        CALL factorize_nir (grid% ni, grid% ni2, grid% nir, ierr)
        !-----------------------------------------------------------------
        ! 2. Calculate the arrays xnglob, rlon, rlat, and corio for all 10
        ! diamonds
        !-----------------------------------------------------------------
        CALL global_coordinates (grid% xnglob, grid% rlon, grid% rlat,      &
          grid% ni, grid% ni2, grid% nir, grid% nd, grid%lbg, grid%ubg, ierr)
        IF (ierr/=0) CALL finish ('construct_atm_grid',&
                                  'Error in SUBROUTINE global_coordinates!')
      !-----------------------
      ! ICON unstructured grid
      !-----------------------
      CASE (DWD6_ICON)
        !-------------------------------------------------
        ! Check presence and validity of icongrid argument
        !-------------------------------------------------
                  lico = present (icongrid)
        if (lico) lico = associated (icongrid)
        !-----------------------------------------------------------------
        ! Set up arrays xnglob, rlon, rlat for given ICON grid
        !-----------------------------------------------------------------
        if (present (lon) .and. present (lat)) then
          if (size (lon) /= size (lat)) &
             call finish ('construct_atm_grid','sizes of lon,lat must match')
          call setup_global_coord (grid% xnglob, grid% rlon, grid% rlat, &
                                   lon, lat, grid% lbg, grid% ubg, ierr  )
        else if (lico) then
          if (grid% nxny /= icongrid% patch% n_patch_cells) then
             print *, "nxny         =", grid% nxny
             print *, "n_patch_cells=", icongrid% patch% n_patch_cells
             call finish ('construct_atm_grid','ICON grid sizes must match')
          end if
          grid% icongrid => icongrid
          grid% global   =  icongrid% global
          grid% dlon(:)  =  reshape (icongrid% patch% cells% center(:,:)% lon * r2d, &
                                     (/ grid% nxny /) )
          grid% dlat(:)  =  reshape (icongrid% patch% cells% center(:,:)% lat * r2d, &
                                     (/ grid% nxny /) )
          call setup_global_coord (grid% xnglob, grid% rlon, grid% rlat,         &
                                   reshape (grid% dlon * d2r,(/grid% nxny,1,1/)),&
                                   reshape (grid% dlat * d2r,(/grid% nxny,1,1/)),&
                                   grid% lbg, grid% ubg, ierr                    )
        else
          if (.not. present (template)) &
             call finish ('construct_atm_grid','ICON: lon and lat required')
          if (grid% gridtype /= template% gridtype) &
             call finish ('construct_atm_grid','ICON: gridtype must not change')
          if (grid% ni /= template% ni) &
             call finish ('construct_atm_grid','ICON: ni must not change')
          ! More checks (uuid etc.) needed...
          grid% xnglob = template% xnglob
          grid% rlon   = template% rlon
          grid% rlat   = template% rlat
        end if
        !-------------------------------------------
        ! Increment link count of ICON grid metadata
        !-------------------------------------------
        if (associated (grid% icongrid)) &
             grid% icongrid% refcount = grid% icongrid% refcount + 1
        !-----------------------------------------
        ! factorize to determine m and nn in RmBnn
        !-----------------------------------------
        if (grid% global) then
           CALL factorize_nir (grid% ni, grid% ni2, grid% nir, ierr)
        end if
        select case (grid% icongrid% grid_root)
        case (2)
           grid% ni2 = grid% icongrid% grid_level + 1
           grid% nir = 1
        case default
           grid% ni2 = grid% icongrid% grid_level
           grid% nir = grid% icongrid% grid_root
        end select
        if (dace% lpio) then
          if (grid% nir == 1) then
            print '(1x,A,       I2.2)', "ICON grid: R02B",            grid% ni2 - 1
          else
            print '(1x,A,I2.2,A,I2.2)', "ICON grid: R",grid% nir,"B", grid% ni2
          end if
        end if

        !---------------------------------
        ! Set up ICON dual grid (vertices)
        !---------------------------------
        if (grid% d_gme(1) /= -1) then
           !--------------------------------------------
           ! non-standard setup, keep data from template
           !--------------------------------------------
           allocate (grid% pe_d (grid% nxny_d, 1,1))
           grid% pe_d = template% pe_d
        else
           !---------------
           ! standard setup
           !---------------
           grid% nxny_d = grid% icongrid% patch% n_patch_verts_g
           grid% nx_d   = grid% nxny_d
           grid% ny_d   = 1
           grid% lbg_d  = [       1     , 1 ,       1 ,  1 ]
           grid% ubg_d  = [ grid% nxny_d, 1 , grid% nz,  1 ]
           grid% lb_d   = grid% lbg_d
           grid% ub_d   = grid% ubg_d

           CALL setup_decomp   (grid% dc_d,        grid% lb_d,  grid% ub_d)
           CALL setup_com_icon (grid% dc_d, marr_, grid% lbg_d, grid% ubg_d)

           allocate (grid% pe_d(grid% nxny_d, 1, 1))
           grid% pe_d(:,:,:) = marr_(1,:,:,:)
           deallocate (marr_)
           !---------------------------------------------
           ! Compress redundant grid information after
           ! having fixed the actual domain decomposition
           !---------------------------------------------
           call compress_icon_grid (grid% icongrid% patch,        &
                                    lbc=grid% lb,   ubc=grid% ub, &
                                    lbv=grid% lb_d, ubv=grid% ub_d)
        end if

      !------------
      ! other grids
      !------------
      CASE default
        IF (PRESENT(TEMPLATE)) THEN
          if(size(grid% dlon)==size(template% dlon)) &
            grid% dlon = template% dlon
          if(size(grid% dlat)==size(template% dlat)) &
            grid% dlat = template% dlat
        ENDIF
        IF (PRESENT(lo1))      grid% lo1 = lo1
        IF (PRESENT(di))       grid% di  = di
        grid% dlon =       (/ (grid% lo1    &
                           +i* grid% di,    &
                           i=0,grid% nx-1) /)
        IF (PRESENT(la1))      grid% la1 = la1
        IF (PRESENT(dj))       grid% dj  = dj
        IF (grid% dj/=0._wp)   grid% dlat = (/ (grid% la1    &
                           +i* grid% dj,    &
                           i=0,grid% ny-1) /)
        IF (PRESENT(dlon)) THEN
          grid% dlon = dlon
        ENDIF
        IF (PRESENT(dlat)) THEN
          grid% dlat = dlat
        ENDIF

        !---------------------
        ! GAUSSIAN or HARMONIC
        !---------------------
        IF (grid%ngl/=0) THEN
          IF (grid%ngl /= grid%ny) THEN
            WRITE (*,*) 'construct_atm_grid: ngl /= ny',grid%ngl, grid%ny
            CALL finish ('construct_atm_grid','ngl /= ny')
          ENDIF
          CALL gauaw (ga = grid% dlat)
          grid% dlat = r2d * ASIN (grid% dlat)
          if (grid% scanmode == WMO8_J_POSITIVE) grid% dlat = - grid% dlat
          grid% la1  = grid% dlat(1)
          grid% dj   = 0._wp
          grid% poly = .TRUE.
          IF(grid% di==0._wp .AND. grid% nx==2*grid%ngl) THEN
            grid% di = 360._wp / grid% nx
            grid% dlon =       (/ (grid% lo1     &
                                +i* grid% di,    &
                                i=0,grid% nx-1) /)
            grid% cyc_x = .TRUE.
          ENDIF
        ENDIF

#if defined(__ibm__)
        ! Work around xlf 9.1 bug with bounds-checking and spread:
        grid% rlon(:,:,1,1) = SPREAD (grid% dlon, 2, grid% ny) * d2r
        grid% rlat(:,:,1,1) = SPREAD (grid% dlat, 1, grid% nx) * d2r
#else
        grid% rlon(:,:,1,1) = SPREAD (grid% dlon * d2r, 2, grid% ny)
        grid% rlat(:,:,1,1) = SPREAD (grid% dlat * d2r, 1, grid% nx)
#endif
        IF (PRESENT(TEMPLATE)) THEN
          if (all(shape(grid% rlon)==shape(template% rlon))) THEN
            grid% rlon = template% rlon
            grid% rlat = template% rlat
          ENDIF
        ENDIF
        IF (PRESENT(lon)) grid% rlon(:,:,1,:) = lon * d2r
        IF (PRESENT(lat)) grid% rlat(:,:,1,:) = lat * d2r
        IF (PRESENT(nn) .OR. PRESENT(ngl))   grid% cyc_x = .TRUE.
        if (grid% gridtype == WMO6_GAUSSIAN) grid% cyc_x = .TRUE.
        if (grid% gridtype == WMO6_LATLON) then
          if (abs (grid% nx * grid% di  - 360._wp) < 1.e-3_wp      ) &
                                             grid% cyc_x = .TRUE.
          if (abs (grid% dlat(1)        +  90._wp) < 1.e-3_wp .and.  &
              abs (grid% dlat(grid% ny) -  90._wp) < 1.e-3_wp      ) &
                                             grid% poly  = .TRUE.
        endif
        grid% global = grid% poly .and. grid% cyc_x
        !---------------------------------------------------------
        ! transform rlon, rlat to geograph. coord. system (rad)
        ! determine lat,lon of rotated north pole and set rot to T
        !---------------------------------------------------------
        if (grid% gridtype == WMO6_ROTLL) then
           grid% rot = .true.
           if (present(lor)) grid% dlonr = lor-180._wp
           if (grid% dlonr < -180._wp) grid% dlonr = grid% dlonr + 360._wp
           if (present(lar)) grid% dlatr = -lar
!          print *, 'start trafo'
           do i=1,grid% nx
              do j=1,grid% ny
                 grid% rlon(i,j,1,1) = d2r*rlarot2rla(grid% dlat(j),grid% dlon(i),grid% dlatr,grid% dlonr,0._wp)
                 grid% rlat(i,j,1,1) = d2r*phirot2phi(grid% dlat(j),grid% dlon(i),grid% dlatr,grid% dlonr,0._wp)
              end do
           end do
!           print *, 'trafo finished'
!           print *, 'rlon(1,1), rlon(1,ny) = ', grid% rlon(1,1,1,1), grid% rlon(1,grid% ny,1,1)
!           print *, 'rlon(nx,ny), rlon(nx,1) = ', grid% rlon(grid% nx,grid% ny,1,1), grid% rlon(grid% nx,1,1,1)
!           print *, 'rlat(1,1), rlat(1,ny) = ', grid% rlat(1,1,1,1), grid% rlat(1,grid% ny,1,1)
!           print *, 'rlat(nx,ny), rlat(nx,1)', grid% rlat(grid% nx,grid% ny,1,1), grid% rlat(grid% nx,1,1,1)
!           print *, 'min(rlon),max(rlon) = ', minval(grid% rlon(:,:,1,1)), maxval(grid% rlon(:,:,1,1))
!           print *, 'min(rlat),max(rlat) = ', minval(grid% rlat(:,:,1,1)), maxval(grid% rlat(:,:,1,1))
!           print *, 'lat(rot SP), lon(rot SP)', lar, lor
!           print *, 'dlon(1), dlat(1) = ', grid% dlon(1), grid% dlat(1)
!           print *, 'dlon(nx), dlat(ny) = ', grid% dlon(grid% nx), grid% dlat(grid% ny)
        end if

        !-----------------------------
        ! set xnglob for regular grids
        !-----------------------------
        do j=1,grid% ny
          do i=1,grid% nx
            sinlat                 =  sin (grid% rlat(i,j,1,1))
            coslat                 =  cos (grid% rlat(i,j,1,1))
            sinlon                 =  sin (grid% rlon(i,j,1,1))
            coslon                 =  cos (grid% rlon(i,j,1,1))
            grid% xnglob (i,j,1,1) =  coslon * coslat ! x
            grid% xnglob (i,j,2,1) =  sinlon * coslat ! y
            grid% xnglob (i,j,3,1) =           sinlat ! z
          end do
        end do

      END SELECT
      !-------------------------------
      ! vertical coordinate parameters
      !-------------------------------
      nz = grid%nz
      if (nz > 0) then
        IF (PRESENT(ak))  grid% ak (:nz+1) = ak
        IF (PRESENT(bk))  grid% bk (:nz+1) = bk
        IF (PRESENT(ak) .or. PRESENT(bk)) THEN
          IF (grid% ak (nz+1) /= 0._wp .or. grid% bk (nz+1) /= 0._wp) then
                          grid% akf(:nz)   = 0.5_wp * &
                                             (grid% ak(1:nz) + grid% ak(2:nz+1))
                          grid% bkf(:nz)   = 0.5_wp * &
                                             (grid% bk(1:nz) + grid% bk(2:nz+1))
          else
                          grid% akf(:nz)   = grid% ak (:nz)
                          grid% ak         = 0._wp
                          grid% bkf(:nz)   = grid% bk (:nz)
                          grid% bk         = 0._wp
          endif
        END IF
        IF (PRESENT(akf)) then
                          grid% akf(:nz)   = akf
          if (grid% levtyp == WMO3_ISOBARIC) &
                          grid% bkf        = 0._wp
        end if
        IF (PRESENT(bkf)) grid% bkf(:nz)   = bkf
      else
        if (dace% lpio) &
          write(6,*) "construct_atm_grid: Warning: nz =", nz
      end if
      !--------------------------------------
      ! derive model top (moved to set_ptopf)
      !--------------------------------------
      grid% ptopf  = 0._wp
      grid% htopf  = 0._wp
      grid% kmlbot = 0
      grid% kmltop = 0
      !--------------------
      ! depth below surface
      !--------------------
      IF (PRESENT(dbs)) grid% dbs(:grid% ns) = nint (dbs(:grid% ns))
      !----------------------------------
      ! copy surface fields from template
      !----------------------------------
      if (present (template)) then
        if (    template% gridtype      == grid% gridtype       .and. &
            all(template% lb((/1,2,4/)) == grid% lb((/1,2,4/))) .and. &
            all(template% ub((/1,2,4/)) == grid% ub((/1,2,4/))) ) then
          do i=1,size(grid% m)
            select case (template% m(i)%i% name)
            case ('lsm','geosp','geoid','soiltyp','fr_lake','depth_lk', &
                  'sso_stdh','sso_gamma','sso_theta','sso_sigma',       &
                  'emis_rad','prs_min','skc')
              if (template% m(i)%i% alloc) then
                call allocate(grid% m(i))
                grid% m(i)% ptr = template% m(i)% ptr
              endif
            end select
          end do
          call set_pointers (grid)
        end if
      end if
    END IF ! nxny > 0
    !-------------------------------------------
    ! derive unique vertical grid representation
    !-------------------------------------------
    select case (grid% levtyp)
    case (WMO3_ISOBARIC)
      grid% vct   = VCT_P_ISO
    case (WMO3_HYBRID ,WMO3_HYBRIDB)
      select case (grid% vc% ivctype)
      case (0)
        grid% vct = VCT_P_HYB
      case (1:)
        grid% vct = VCT_Z_HYB
      case default
        write (6,*)  'construct_atm_grid:  invalid ivctype =',grid% vc% ivctype
        write (0,*)  'construct_atm_grid:  invalid ivctype =',grid% vc% ivctype
        call finish ('construct_atm_grid','invalid ivctype')
      end select
    case (WMO3_GENV)
      grid% vct   = VCT_Z_GEN
    case default
      if (dace% lpio) then
        write (6,*)  'construct_atm_grid:  invalid levtyp =',grid% levtyp
        write (6,*)  '                     vct set to        VCT_NO'
      endif
      grid% vct   = VCT_NO
    end select
    !---------------------------------------
    ! estimate model, set nominal resolution
    !---------------------------------------
    grid% model     = MO_UNKNOWN
    grid% d_km      = 0._wp
    grid% d_deg     = 0._wp
    select case (grid% gridtype)
    case (WMO6_GAUSSIAN)
      grid% model   = MO_IFS
      grid% d_deg   = grid% di
      grid% d_km    = d2r * grid% d_deg * 0.001_wp * grid% a
    case (DWD6_ICOSAHEDRON)
      grid% model   = MO_GME
      grid% d_km    = 7054._wp / grid% ni
      grid% d_deg   = r2d * grid% d_km / (0.001_wp * grid% a)
    case (DWD6_ICON)
      !------------------------------------
      ! adjust ni for ICON-Nest or ICON-LAM
      !------------------------------------
      if (associated (grid% icongrid)) then
        grid% ni = grid% icongrid% grid_root * 2 ** grid% icongrid% grid_level
      endif
      grid% model   = MO_ICON
      grid% d_km    = 7054._wp / (grid% ni * sqrt(2._wp))
      grid% d_deg   = r2d * grid% d_km / (0.001_wp * grid% a)
    case (WMO6_LATLON, WMO6_ROTLL)
      select case (grid% vct)
      case (VCT_Z_GEN, VCT_Z_HYB)
        grid% model = MO_COSMO
      case (VCT_P_HYB)
        if (grid% arakawa == "C") then
          grid% model = MO_HRM
        else
          grid% model = MO_IFS
        end if
      end select
      grid% d_deg   = grid% dj
      grid% d_km    = d2r * grid% d_deg * 0.001_wp * grid% a
    end select
  END SUBROUTINE construct_atm_grid
!==============================================================================
  subroutine set_ptopf (grid)
    !--------------------------------------------
    ! Set or estimate topmost full level pressure
    ! Set or estimate topmost height (gpm)
    !--------------------------------------------
    type (t_grid) ,INTENT(inout) :: grid

    real(wp) :: h, h2
    integer  :: i,j,l

    grid% ptopf = 0._wp
    grid% htopf = 0._wp
    select case (grid% levtyp)
    case (WMO3_HYBRID, WMO3_HYBRIDB, WMO3_ISOBARIC)
      if (grid% vc% ivctype == 0) then
        grid% ptopf = grid% akf(1)
      else if (associated (grid% p0)) then
        grid% ptopf = maxval (grid% p0 (:,:,1,:))
        grid% ptopf = p_max  (grid% ptopf, comm = grid% dc% comm)
      endif
      !--------------------------------------------
      ! approximate htopf from topmost hybrid level
      !--------------------------------------------
      grid% htopf = h_p_usstd (grid% ptopf)
    case (WMO3_HHYBRID, WMO3_GENV)
      !---------------------------------------
      ! Height-based hybrid levels
      ! Estimate ptopf from topmost full level
      !---------------------------------------
      if (associated (grid% hhl)) then
        h = huge(h)
        do l = grid% lb(4), grid% ub(4)
        do j = grid% lb(2), grid% ub(2)
        do i = grid% lb(1), grid% ub(1)
          h2 = (grid% hhl(i,j,1,l) + grid% hhl(i,j,2,l)) * 0.5_wp
          if (h2 /= 9999._wp) h = min (h, h2)
        end do
        end do
        end do
        h = p_min  (h, comm = grid% dc% comm)
      else
        if (grid% levtyp == WMO3_GENV) &
          call finish ('set_ptopf','HHL not associated!')
        h = grid% akf(1)
      end if
      grid% htopf =            h
      grid% ptopf = p_h_usstd (h)
    case default
!     write(0,*)   'set_ptopf:  unknown levtyp',grid% levtyp
!     call finish ('set_ptopf','unknown levtyp')
!   case (WMO3_SEALEVEL, WMO3_SURFACE)
      !---------------------------------------------------------------
      ! Don't do anything, but warn about missing vertical coordinates
      !---------------------------------------------------------------
      if (dace% lpio) then
        write(6,*)
        write(6,*) "  set_ptopf: Warning: levtyp =", grid% levtyp
        write(6,*)
      endif
    end select
  end subroutine set_ptopf
!==============================================================================
  subroutine set_plev_indices (grid)
  type (t_grid) ,intent(inout) :: grid
  !----------------------------------------------------
  ! set model level indices for typical pressure levels
  !----------------------------------------------------

    real(wp) :: z (grid% nz)   ! typical height   at model level
    real(wp) :: p (grid% nz)   ! typical pressure at model level
    real(wp) :: zh(grid% nz+1) ! typical model half-level height
    integer  :: k, ke, ix(3), i, j, d, pe

    !-----------------------------------------------
    ! set up typical pressure levels at model height
    !-----------------------------------------------
    grid% kmlbot = 0
    grid% kmltop = 0
    ke           = grid% nz

    select case (grid% vct)
    case (VCT_P_ISO)

      p = grid% akf (1:ke)

    case (VCT_P_HYB)

      p = grid% akf (1:ke) + grid% bkf (1:ke) * 100000._wp

    case (VCT_Z_HYB, VCT_Z_GEN)


      if (associated (grid% hsurf)) then
        ix = minloc (grid% hsurf(:,:,1,:), grid% hsurf(:,:,1,:) >= 0._wp)
      else if (associated (grid% geosp)) then
        ix = minloc (grid% geosp(:,:,1,:), grid% geosp(:,:,1,:) >= 0._wp)
      else
        goto 99
      endif

      if (.not.associated (grid% hhl  )) goto 99

      pe    = -huge (pe)
      zh    = -huge (zh)
      ix(1) = ix(1) + grid% lbg(1) - 1
      ix(2) = ix(2) + grid% lbg(2) - 1
      i     = ix(1)
      j     = ix(2)
      d     = ix(3)
      if (grid% lb(1) <= i .and. i <= grid% ub(1) .and. &
          grid% lb(2) <= j .and. j <= grid% ub(2) .and. &
          grid% lb(4) <= d .and. d <= grid% ub(4)       ) then
         pe = dace% pe
         zh = grid% hhl(i,j,1:ke+1,d)
      end if
      pe = p_max (pe)
      zh = p_max (zh)
      z  = (zh(1:ke) + zh(2:ke+1)) * 0.5_wp
      p  = p_h_usstd (z)
      call p_bcast (p, pe)
    case default
      goto 99
    end select
    !------------------
    ! determine indices
    !------------------
    grid% kmlbot = ke
    grid% kmltop = ke - 1
    do k = 1 , ke - 1
      if ((p(k)  <= pbot_lapse*100._wp) .and. &
          (p(k+1) > pbot_lapse*100._wp)       ) grid% kmlbot = k
      if ((p(k)  <= ptop_lapse*100._wp) .and. &
          (p(k+1) > ptop_lapse*100._wp)       ) grid% kmltop = k
    enddo
    grid% kmlbot = max( grid% kmlbot , grid% kmltop + 1 )

99  continue

  end subroutine set_plev_indices
!==============================================================================
  subroutine set_zlev_ref (grid)
    type(t_grid) ,intent(inout) :: grid
    !--------------------------------------------------------
    ! estimate unperturbed model level heights over sea level
    !--------------------------------------------------------
    real(wp) :: zh(grid% nz+1)          ! typical height at model level
    integer  :: ke, ix(3), i, j, d, pe

    if (.not. associated (grid% hhl  )) return
    if (.not. associated (grid% hsurf)) &
         call finish ("set_zlev_ref","hsurf not allocated")

    select case (grid% vct)
    case (VCT_Z_HYB, VCT_Z_GEN)
       ke    = grid% nz
       pe    = -huge (pe)
       zh    = -huge (zh)
       ix    = minloc (abs (grid% hsurf(:,:,1,:)))
       ix(1) = ix(1) + grid% lbg(1) - 1
       ix(2) = ix(2) + grid% lbg(2) - 1
       i     = ix(1)
       j     = ix(2)
       d     = ix(3)
       if (grid% lb(1) <= i .and. i <= grid% ub(1) .and. &
           grid% lb(2) <= j .and. j <= grid% ub(2) .and. &
           grid% lb(4) <= d .and. d <= grid% ub(4)       ) then
          pe = dace% pe
          zh = grid% hhl(i,j,1:ke+1,d)
       end if
       pe = p_max (pe)
       zh = p_max (zh)
       call p_bcast (zh, pe)
       grid% ak (1:ke+1) =  zh(1:ke+1)
       grid% akf(1:ke)   = (zh(1:ke) + zh(2:ke+1)) * 0.5_wp
    case default
       call finish ("set_zlev_ref","unsupported vct")
    end select
  end subroutine set_zlev_ref
!==============================================================================
  SUBROUTINE destruct_atm_grids (grid)
  !----------------------------------------
  ! deallocate pointer components of t_grid
  !----------------------------------------
  TYPE (t_grid) ,INTENT(inout) :: grid (:)
    integer :: i
    do i = 1, size(grid)
      call destruct_atm_grid (grid(i))
    end do
  END SUBROUTINE destruct_atm_grids
!------------------------------------------------------------------------------
  SUBROUTINE destruct_atm_grid (grid)
  !----------------------------------------
  ! deallocate pointer components of t_grid
  !----------------------------------------
  TYPE (t_grid) ,INTENT(inout) :: grid

    IF (grid% size > 0) THEN
      DEALLOCATE (grid%dlon)
      DEALLOCATE (grid%dlat)
      CALL destruct (grid% m)
      CALL destruct (grid% dc)
      CALL destruct (grid% dc_d)
      CALL destruct (grid% vc)
      IF (ASSOCIATED(grid% marr)) DEALLOCATE (grid% marr)
      IF (ASSOCIATED(grid% pe_d)) DEALLOCATE (grid% pe_d)
    ENDIF
    grid% size = 0

    IF (ASSOCIATED(grid% icongrid)) then
      !----------------------------------------------------
      ! Decrement and test link count of ICON grid metadata
      !----------------------------------------------------
      grid% icongrid% refcount = grid% icongrid% refcount - 1
      if (grid% icongrid% refcount <= 0) then
        call finish_icongrid (grid% icongrid)
        deallocate (grid% icongrid)
      end if
      nullify (grid% icongrid)
    END IF

  END SUBROUTINE destruct_atm_grid
!------------------------------------------------------------------------------
  SUBROUTINE destruct_vcoord (vc)
  TYPE (t_vcoord) ,INTENT(inout) :: vc
  !------------------------------------------
  ! deallocate pointer components of t_vcoord
  !------------------------------------------
    IF (ASSOCIATED (vc% vert_coord)) DEALLOCATE (vc% vert_coord)
    IF (ASSOCIATED (vc% sigm_coord)) DEALLOCATE (vc% sigm_coord)
  END SUBROUTINE destruct_vcoord
!==============================================================================
  subroutine set_empty_grid (grid)
    type(t_grid) TARGET ,INTENT(inout) :: grid
    integer :: i
    !============================
    ! Set defaults for empty grid
    !============================
    grid% gridtype    = -1
    grid% nx          = 0
    grid% ny          = 0
    grid% nz          = 1
    grid% nd          = 1
    grid% ns          = 0
    grid% nxny        = 0
    grid% di          = 0._wp
    grid% dj          = 0._wp
    grid% ngl         = 0
    grid% scanmode    = WMO8_J_POSITIVE    ! DWD default: +i,+j
    grid% nn          = 0
    grid% ni          = 0
    grid% d_gme       = -1
    grid% a           = rearth
    grid% g           = gacc
    grid% rot         = .FALSE.
    grid% cyc_x       = .FALSE.
    grid% poly        = .FALSE.
    grid% global      = .FALSE.
    grid% lo1         = 0._wp
    grid% la1         = 0._wp
    grid% dlonr       = 180._wp
    grid% dlatr       = 90._wp
    grid% ldec        = .FALSE.
    grid% levtyp      = WMO3_HYBRID
    grid% ak          = 0._wp
    grid% bk          = 0._wp
    grid% akf         = 0._wp
    grid% bkf         = 0._wp
    grid% dbs         = 0
    grid% arakawa     = 'A'
    grid% vct         = VCT_NO
    grid% icongrid    => NULL ()
    grid% nx_d        = 0
    grid% ny_d        = 0
    grid% nxny_d      = 0
    grid% lbg_d       = 1
    grid% lb_d        = 1
    grid% ubg_d       = 0
    grid% ub_d        = 0
    grid% ptopf       = 0._wp
    grid% htopf       = 0._wp
    grid% kmlbot      = 0
    grid% kmltop      = 0
    !-------------------------------
    ! Nullify all pointer components
    !-------------------------------
    nullify (grid% dlon)
    nullify (grid% dlat)
    nullify (grid% marr)
    nullify (grid% pe_d)
    DO i=1,SIZE(grid% m)
      NULLIFY (grid% m(i)% ptr)         ! With TR15581, ptr is ALLOCATABLE
    END DO
  end subroutine set_empty_grid
!==============================================================================
  subroutine set_grid_fields (grid)
    type(t_grid) TARGET ,INTENT(inout) :: grid
    !---------------------------------------------------------------
    ! set mnemonics and bounds for the grid specific fields
    ! set pointers to fields consistently with pointers in array 'm'
    !---------------------------------------------------------------
    integer :: i

    grid%m%i% name = [                                                           &
         'lsm      ','geosp    ','rlon     ','rlat     ','geo_s    ','geoid    ',&
         'soiltyp  ','fr_lake  ','depth_lk ','xnglob   ','hsurf    ','hhl      ',&
         'p0       ','dp0      ','rho0     ','sso_stdh ','sso_gamma','sso_theta',&
         'sso_sigma','emis_rad ','prs_min  ','skc      '                         ]
    !---------------------------------------------------
    ! default: single level fields, allocated on each PE
    !---------------------------------------------------
    DO i=1,SIZE(grid%m)
      grid%m(i)%i% lb    = grid% lbg
      grid%m(i)%i% ub    = grid% ubg
      grid%m(i)%i% ub(3) = 1
    END DO
    grid%m%i% nn = grid% nn
    grid%m(:)        %i% alloc = .FALSE.
    grid%m(RLON:RLAT)%i% alloc = .TRUE.        ! rlon, rlat
    grid%m(XNGLOB)   %i% alloc = .TRUE.        ! xnglob
    grid%m(GEO_S)    %i% rep   = 'sh'          ! geo_s
    grid%m(XNGLOB)   %i% ub(3) = 3             ! xnglob

    !-----------------------------------------
    ! multi-level fields, distributed over PEs
    !-----------------------------------------
    DO i = HHL, RHO0
      grid%m(i)%i% lb    = grid% lb       ! hhl, p0, dp0, rho0
      grid%m(i)%i% ub    = grid% ub
    END DO

    grid%m(HHL) %i% ub(3) = grid% ub(3)+1 ! hhl

    !------------------------------------------
    ! single-level fields, distributed over PEs
    !------------------------------------------
    DO i = SSO_STDH, SKC
      grid%m(i)%i% lb    = grid% lb       ! sso_*
      grid%m(i)%i% ub    = grid% ub
      grid%m(i)%i% lb(3) = 1
      grid%m(i)%i% ub(3) = 1
    END DO
    ! geoid decomposition
    grid%m(GEOID)%i% lb    = grid% lb
    grid%m(GEOID)%i% ub    = grid% ub
    grid%m(GEOID)%i% lb(3) = 1
    grid%m(GEOID)%i% ub(3) = 1
  end subroutine set_grid_fields
!==============================================================================
  SUBROUTINE set_grid_pointers (grid)
  !---------------------------------------------------------------
  ! set pointers to fields consistently with pointers in array 'm'
  !---------------------------------------------------------------
  TYPE (t_grid) TARGET ,INTENT(inout) :: grid
    grid% lsm       => grid%m(LSM      )% ptr
    grid% geosp     => grid%m(GEOSP    )% ptr
    grid% rlon      => grid%m(RLON     )% ptr
    grid% rlat      => grid%m(RLAT     )% ptr

    grid% geoid     => grid%m(GEOID    )% ptr
    grid% soiltyp   => grid%m(SOILTYP  )% ptr
    grid% fr_lake   => grid%m(FR_LAKE  )% ptr
    grid% depth_lk  => grid%m(DEPTH_LK )% ptr
    grid% xnglob    => grid%m(XNGLOB   )% ptr
    grid% hsurf     => grid%m(HSURF    )% ptr
    grid% hhl       => grid%m(HHL      )% ptr
    grid% p0        => grid%m(P0       )% ptr
    grid% dp0       => grid%m(DP0      )% ptr
    grid% rho0      => grid%m(RHO0     )% ptr
    grid% sso_stdh  => grid%m(SSO_STDH )% ptr
    grid% sso_gamma => grid%m(SSO_GAMMA)% ptr
    grid% sso_theta => grid%m(SSO_THETA)% ptr
    grid% sso_sigma => grid%m(SSO_SIGMA)% ptr
    grid% emis_rad  => grid%m(EMIS_RAD )% ptr
    grid% prs_min   => grid%m(PRS_MIN  )% ptr
    grid% skc       => grid%m(SKC      )% ptr
    nullify (grid% geo_sh)
    IF (ALLOCATED(grid%m(GEO_S)% ptr)) grid% geo_sh => grid%m(GEO_S)% ptr (:,1,1,1)
  END SUBROUTINE set_grid_pointers
!==============================================================================
  SUBROUTINE allocate_atm_grid (grid, field)
  !------------------------------------------
  ! allocate a specific field
  ! (array pointer components of type t_grid)
  !------------------------------------------
  TYPE (t_grid)     ,INTENT(inout) TARGET :: grid  ! grid variable to allocate
  CHARACTER (len=*) ,INTENT(in)           :: field ! mnemonic of field

    SELECT CASE (field)
    CASE ('lsm')
      CALL ALLOCATE (grid%m(LSM      )); grid% lsm       => grid%m(LSM      )%ptr
    CASE ('geosp')
      CALL ALLOCATE (grid%m(GEOSP    )); grid% geosp     => grid%m(GEOSP    )%ptr
    CASE ('geoid')
      CALL ALLOCATE (grid%m(GEOID    )); grid% geoid     => grid%m(GEOID    )%ptr
    CASE ('soiltyp')
      CALL ALLOCATE (grid%m(SOILTYP  )); grid% soiltyp   => grid%m(SOILTYP  )%ptr
    CASE ('fr_lake')
      CALL ALLOCATE (grid%m(FR_LAKE  )); grid% fr_lake   => grid%m(FR_LAKE  )%ptr
    CASE ('depth_lk')
      CALL ALLOCATE (grid%m(DEPTH_LK )); grid% depth_lk  => grid%m(DEPTH_LK )%ptr
    CASE ('hsurf')
      CALL ALLOCATE (grid%m(HSURF    )); grid% hsurf     => grid%m(HSURF    )%ptr
    CASE ('hhl')
      CALL ALLOCATE (grid%m(HHL      )); grid% hhl       => grid%m(HHL      )%ptr
    CASE ('p0')
      CALL ALLOCATE (grid%m(P0       )); grid% p0        => grid%m(P0       )%ptr
    CASE ('dp0')
      CALL ALLOCATE (grid%m(DP0      )); grid% dp0       => grid%m(DP0      )%ptr
    CASE ('rho0')
      CALL ALLOCATE (grid%m(RHO0     )); grid% rho0      => grid%m(RHO0     )%ptr
    CASE ('sso_stdh')
      CALL ALLOCATE (grid%m(SSO_STDH )); grid% sso_stdh  => grid%m(SSO_STDH )%ptr
    CASE ('sso_gamma')
      CALL ALLOCATE (grid%m(SSO_GAMMA)); grid% sso_gamma => grid%m(SSO_GAMMA)%ptr
    CASE ('sso_theta')
      CALL ALLOCATE (grid%m(SSO_THETA)); grid% sso_theta => grid%m(SSO_THETA)%ptr
    CASE ('sso_sigma')
      CALL ALLOCATE (grid%m(SSO_SIGMA)); grid% sso_sigma => grid%m(SSO_SIGMA)%ptr
    CASE ('emis_rad')
      CALL ALLOCATE (grid%m(EMIS_RAD )); grid% emis_rad  => grid%m(EMIS_RAD )%ptr
    CASE ('prs_min')
      CALL ALLOCATE (grid%m(PRS_MIN  )); grid% prs_min   => grid%m(PRS_MIN  )%ptr
    CASE ('skc')
      CALL ALLOCATE (grid%m(SKC      )); grid% skc       => grid%m(SKC      )%ptr
    CASE DEFAULT
      call finish ('allocate_atm_grid','unknown field: '//field)
    END SELECT
  END SUBROUTINE allocate_atm_grid
!==============================================================================
  SUBROUTINE deallocate_atm_grid (grid, field)
  !------------------------------------------
  ! deallocate a specific field
  ! (array pointer components of type t_grid)
  !------------------------------------------
  TYPE (t_grid)     ,INTENT(inout) TARGET :: grid  ! grid variable to dealloc
  CHARACTER (len=*) ,INTENT(in)           :: field ! mnemonic of field

    SELECT CASE (field)
    CASE ('lsm')
      CALL DEALLOCATE (grid%m(LSM      )); nullify (grid% lsm)
    CASE ('geosp')
      CALL DEALLOCATE (grid%m(GEOSP    )); nullify (grid% geosp)
    CASE ('geoid')
      CALL DEALLOCATE (grid%m(GEOID    )); nullify (grid% geoid)
    CASE ('soiltyp')
      CALL DEALLOCATE (grid%m(SOILTYP  )); nullify (grid% soiltyp)
    CASE ('fr_lake')
      CALL DEALLOCATE (grid%m(FR_LAKE  )); nullify (grid% fr_lake)
    CASE ('depth_lk')
      CALL DEALLOCATE (grid%m(DEPTH_LK )); nullify (grid% depth_lk)
    CASE ('hsurf')
      CALL DEALLOCATE (grid%m(HSURF    )); nullify (grid% hsurf)
    CASE ('hhl')
      CALL DEALLOCATE (grid%m(HHL      )); nullify (grid% hhl)
    CASE ('p0')
      CALL DEALLOCATE (grid%m(P0       )); nullify (grid% p0)
    CASE ('dp0')
      CALL DEALLOCATE (grid%m(DP0      )); nullify (grid% dp0)
    CASE ('rho0')
      CALL DEALLOCATE (grid%m(RHO0     )); nullify (grid% rho0)
    CASE ('sso_stdh')
      CALL DEALLOCATE (grid%m(SSO_STDH )); nullify (grid% sso_stdh)
    CASE ('sso_gamma')
      CALL DEALLOCATE (grid%m(SSO_GAMMA)); nullify (grid% sso_gamma)
    CASE ('sso_theta')
      CALL DEALLOCATE (grid%m(SSO_THETA)); nullify (grid% sso_theta)
    CASE ('sso_sigma')
      CALL DEALLOCATE (grid%m(SSO_SIGMA)); nullify (grid% sso_sigma)
    CASE ('emis_rad')
      CALL DEALLOCATE (grid%m(EMIS_RAD )); nullify (grid% emis_rad)
    CASE ('prs_min')
      CALL DEALLOCATE (grid%m(PRS_MIN  )); nullify (grid% prs_min)
    CASE ('skc')
      CALL DEALLOCATE (grid%m(SKC      )); nullify (grid% skc)
    CASE DEFAULT
      call finish ('deallocate_atm_grid','unknown field: '//field)
    END SELECT
  END SUBROUTINE deallocate_atm_grid
!==============================================================================
  SUBROUTINE init_ctl_from_grid (ctl, file, grid, title,         &
                                      tdefn, tdefi, tdefd, undef,&
                                      zrev, tmpl, comment        )
  !-----------------------------------------------------------------------
  ! set the GRADS meta data (type t_ctl) from the grid metat data (t_grid)
  !-----------------------------------------------------------------------
  TYPE (t_ctl)     ,INTENT(out)          :: ctl        ! GRADS meta data
  CHARACTER(len=*) ,INTENT(in) ,OPTIONAL :: file       ! GRADS file name
  TYPE (t_grid)    ,INTENT(in)           :: grid       ! grid meta data
  character(len=*) ,INTENT(in) ,OPTIONAL :: title      ! title
  INTEGER          ,INTENT(in) ,OPTIONAL :: tdefn      ! # of time slices
  CHARACTER(len=*) ,INTENT(in) ,OPTIONAL :: tdefi      ! first time
  CHARACTER(len=*) ,INTENT(in) ,OPTIONAL :: tdefd      ! time increments
  REAL(wp)         ,INTENT(in) ,OPTIONAL :: undef      ! GRADS undefined value
  LOGICAL          ,INTENT(in) ,OPTIONAL :: zrev       ! reversed z coordinate
  TYPE (t_ctl)     ,INTENT(in) ,OPTIONAL :: tmpl       ! template
  CHARACTER(len=*) ,INTENT(in) ,OPTIONAL :: comment(:) ! comments for .ctl-file

    logical                        :: lpio
    real(wp)                       :: di, dj
    integer                        :: nc
    character(len=78) ,allocatable :: c (:)
    lpio = .true.
    if (grid% dc% nproc1 * grid% dc% nproc2 /= 1) lpio = (dace% lpio)

    if (lpio) then
      !----------------------
      ! prepare comment lines
      !----------------------
      nc = 0
      if (present (comment)) nc = size (comment)
      allocate (c (nc+2))
      if (present (comment)) c (1:nc) = comment
      c (nc+1) = 'gridtype = '//char3(grid% gridtype)
      c (nc+2) = 'levtyp   = '//char3(grid% levtyp)
      !----------------------------------
      ! handle different horizontal grids
      !----------------------------------
      select case (grid% gridtype)
      case (DWD6_ICOSAHEDRON)
        di = 63.486_wp / grid% ni
        dj = di
      case default
        di = grid% di
        dj = grid% dj
      end select
      !--------------------------------
      ! handle different vertical grids
      !--------------------------------
      select case (grid% levtyp)
      case (WMO3_ISOBARIC)
        CALL init_ctl (ctl                                             ,&
                       file    =       file                            ,&
                       title   =       title                           ,&
                       nx      = grid% nx                              ,&
                       ny      = grid% ny                              ,&
                       lo1     = grid% lo1                             ,&
                       la1     = grid% la1                             ,&
                       di      =       di                              ,&
                       dj      =       dj                              ,&
                       ngl     = grid% ngl                             ,&
                       ke      = grid% nz                              ,&
                       zlev    = grid% akf / 100.,    &! levels in hPa
                       tdefn   =       tdefn                           ,&
                       tdefi   =       tdefi                           ,&
                       tdefd   =       tdefd                           ,&
                       undef   =       undef                           ,&
                       comment =       c                               ,&
                       zrev    =       zrev                            ,&
                       yrev    = grid% dlat(1) > grid% dlat(grid% ny)  ,&
                       tmpl    =       tmpl                             )
      case default
        CALL init_ctl (ctl                                             ,&
                       file    =       file                            ,&
                       nx      = grid% nx                              ,&
                       ny      = grid% ny                              ,&
                       lo1     = grid% lo1                             ,&
                       la1     = grid% la1                             ,&
                       di      =       di                              ,&
                       dj      =       dj                              ,&
                       ngl     = grid% ngl                             ,&
                       ke      = grid% nz                              ,&
                       tdefn   =       tdefn                           ,&
                       tdefi   =       tdefi                           ,&
                       tdefd   =       tdefd                           ,&
                       undef   =       undef                           ,&
                       comment =       c                               ,&
                       zrev    =       zrev                            ,&
                       yrev    = grid% dlat(1) > grid% dlat(grid% ny)  ,&
                       tmpl    =       tmpl                             )
      end select
    endif
  END SUBROUTINE init_ctl_from_grid
!------------------------------------------------------------------------------
  SUBROUTINE grid_to_grads (ctl, grid)
  TYPE (t_ctl)  ,INTENT(inout) :: ctl  ! GRADS meta data to use and extend
  TYPE (t_grid) ,INTENT(in)    :: grid ! model grid data
  !------------------------------------------------------------------
  ! write the invariant fields of derived type t_grid to a GRADS file
  !------------------------------------------------------------------
    CALL to_grads (ctl, grid% m, grid% dc,                  &
                   yrev=grid% dlat(1) > grid% dlat(grid% ny))
  END SUBROUTINE grid_to_grads
!==============================================================================
  SUBROUTINE print_atm_grid (grid, iunit, intend, verbose)
  !--------------------------------------------------------------
  ! Produce a printout of the grid meta data.
  !   The 'intend' string is written at the start of each line
  !     for printout of contained derived type components.
  !   If verbose is .true. min/max values of invariant fields are
  !     printed level by level.
  !--------------------------------------------------------------
  TYPE (t_grid)    ,INTENT(in)           :: grid    ! grid description to print
  INTEGER          ,INTENT(in) ,OPTIONAL :: iunit   ! output unit   (default=6)
  CHARACTER(len=*) ,INTENT(in) ,OPTIONAL :: intend  ! intend string (default'')
  LOGICAL          ,INTENT(in) ,OPTIONAL :: verbose ! verbosity flag(default=F)

    INTEGER           :: iu, ni, k
    CHARACTER(len=32) :: ci
    LOGICAL           :: v
    REAL(wp) ,pointer :: la1(:), lo1(:)
    real(wp)          :: ps_ref         ! Surface pressure for US standard atm.
    real(wp)          :: pf
    CHARACTER(len=32) :: uuid
    !--------------------
    ! optional parameters
    !--------------------
    ni=0
    ci=''
    IF (PRESENT(intend)) THEN
      ci = intend
      ni = LEN_TRIM(ci)
    ENDIF
    !-------------------------------------------
    ! print globally allocated fields on PE p_io
    !-------------------------------------------
    IF (dace% lpio) THEN
      iu = 6       ;IF (PRESENT(iunit))   iu = iunit
      v  = .FALSE. ;IF (PRESENT(verbose)) v  = verbose
      !------
      ! print
      !------
      nullify(la1); IF (ASSOCIATED(grid% rlat)) la1 => grid% rlat (1,:,1,1)
      nullify(lo1); IF (ASSOCIATED(grid% rlon)) lo1 => grid% rlon (:,1,1,1)
      IF(v) WRITE(iu,"(a,a        )") ci(:ni),' (t_grid) :'
            WRITE(iu,"(a,a,i10    )") ci(:ni),' grid   : ',grid% gridtype
            WRITE(iu,"(a,a,i10    )") ci(:ni),' nx     : ',grid% nx
            WRITE(iu,"(a,a,i10    )") ci(:ni),' ny     : ',grid% ny
            WRITE(iu,"(a,a,i10    )") ci(:ni),' ngl    : ',grid% ngl
            WRITE(iu,"(a,a,i10    )") ci(:ni),' ni     : ',grid% ni
         if (grid% d_gme(0) < 0) then
            WRITE(iu,"(a,a,4x,11i6)") ci(:ni),' d_gme  : ',grid% d_gme(0)
         else
            WRITE(iu,"(a,a,4x,11i6)") ci(:ni),' d_gme  : ',grid% d_gme(0:)
         end if
!           WRITE(iu,"(a,a,4x,11i6)") ci(:ni),' lb     : ',grid% lb
!           WRITE(iu,"(a,a,4x,11i6)") ci(:ni),' ub     : ',grid% ub
            WRITE(iu,"(a,a,i10,3i6)") ci(:ni),' lbg    : ',grid% lbg
            WRITE(iu,"(a,a,i10,3i6)") ci(:ni),' ubg    : ',grid% ubg
            if (grid% gridtype == DWD6_ICON) then
              WRITE(iu,"(a,a,i10)") ci(:ni),' gridnum: ',&
                                             grid% icongrid% grid_num
              uuid = byte2hex (grid% icongrid% uuid)
              WRITE(iu,"(a,a,a  )") ci(:ni),' uuid(h): ',uuid
            end if
            WRITE(iu,"(a,a,i10  )") ci(:ni),' levtyp : ',grid% levtyp
            WRITE(iu,"(a,a,i10  )") ci(:ni),' ivctype: ',grid% vc% ivctype
            if (grid% levtyp == WMO3_GENV) then
               uuid = byte2hex (grid% vc% vc_uuid)
               WRITE(iu,"(a,a,a )") ci(:ni),' uuid(v): ',uuid
            end if
            !-----------------------------------------------------------
            ! COSMO reference atmosphere parameters (hybrid levels only)
            !-----------------------------------------------------------
!           else if (grid% ivctype /= 0 .and. grid% ivctype /= IVCTYPE_ICON)
            WRITE(iu,"(a,a,i10  )") ci(:ni),' irefatm: ',grid% ra% irefatm
         if (grid% ra% irefatm >= 0) then
            WRITE(iu,"(a,a,f16.5)") ci(:ni),' p0sl   : ',grid% ra% p0sl
            WRITE(iu,"(a,a,f16.5)") ci(:ni),' t0sl   : ',grid% ra% t0sl
            WRITE(iu,"(a,a,f16.5)") ci(:ni),' dt0lp  : ',grid% ra% dt0lp
            WRITE(iu,"(a,a,f16.5)") ci(:ni),' vcflat : ',grid% vc% vcflat
         end if

            if (grid% gridtype == WMO6_GAUSSIAN) &
            WRITE(iu,"(a,a,i10  )") ci(:ni),' scan   : ',grid% scanmode
            WRITE(iu,"(a,a,i10  )") ci(:ni),' nn     : ',grid% nn
            WRITE(iu,"(a,a,i10  )") ci(:ni),' nz     : ',grid% nz
            WRITE(iu,"(a,a,i10  )") ci(:ni),' ns     : ',grid% ns
            WRITE(iu,"(a,a,l1   )") ci(:ni),' rot    : ',grid% rot
            WRITE(iu,"(a,a,l1   )") ci(:ni),' cyc_x  : ',grid% cyc_x
            WRITE(iu,"(a,a,l1   )") ci(:ni),' poly   : ',grid% poly
            WRITE(iu,"(a,a,l1   )") ci(:ni),' global : ',grid% global
            WRITE(iu,"(a,a,a    )") ci(:ni),' arakaw : ',grid% arakawa
            WRITE(iu,"(a,a,f16.5)") ci(:ni),' d_km   : ',grid% d_km
            WRITE(iu,"(a,a,f16.5)") ci(:ni),' d_deg  : ',grid% d_deg
      if (grid% rot) then
            WRITE(iu,"(a,a,f16.5)") ci(:ni),' dlonr  : ',grid% dlonr
            WRITE(iu,"(a,a,f16.5)") ci(:ni),' dlatr  : ',grid% dlatr
      end if
      IF(v) THEN
            WRITE(iu,"(a,a,f16.5)") ci(:ni),' la1    : ',grid% la1
            WRITE(iu,"(a,a,f16.5)") ci(:ni),' lo1    : ',grid% lo1
            WRITE(iu,"(a,a,f16.5)") ci(:ni),' di     : ',grid% di
            WRITE(iu,"(a,a,f16.5)") ci(:ni),' dj     : ',grid% dj
            WRITE(iu,"(a,a,f16.5)") ci(:ni),' a      : ',grid% a
            WRITE(iu,"(a,a,f16.5)") ci(:ni),' g      : ',grid% g
            WRITE(iu,"(a,a,f16.5)") ci(:ni),' ptopf  : ',grid% ptopf
            WRITE(iu,"(a,a,f16.5)") ci(:ni),' htopf  : ',grid% htopf
            WRITE(iu,"(a,a,i10  )") ci(:ni),' kmlbot : ',grid% kmlbot
            WRITE(iu,"(a,a,i10  )") ci(:ni),' kmltop : ',grid% kmltop
            WRITE(iu,"(a,a,i10  )") ci(:ni),' size   : ',grid% size
            CALL prfield1n (grid% dlon,    ' dlon     ')
            CALL prfield1n (grid% dlat,    ' dlat     ')
            CALL prfield1n (      lo1,     ' lon(1)   ')
            CALL prfield1n (      la1,     ' lat(1)   ')
      ENDIF
            CALL prfield  (grid% rlon,     ' rlon     ',fak=57.296_wp,dim='[degree]')
            CALL prfield  (grid% rlat,     ' rlat     ',fak=57.296_wp,dim='[degree]')
            CALL prfield  (grid% lsm,      ' lsm      ')
            CALL prfield  (grid% geosp,    ' geosp    ',fak=0.10197_wp,dim='[m]')
            CALL prfield  (grid% geoid,    ' geoid    ')
            CALL prfield  (grid% soiltyp,  ' soiltyp  ')
            CALL prfield  (grid% fr_lake,  ' fr_lake  ')
            CALL prfield  (grid% depth_lk, ' depth_lk ')
      IF(v) THEN
            CALL prfield  (grid% xnglob,   ' xnglob   ')
            CALL prfield1 (grid% geo_sh,   ' geo_sh   ')
            CALL prfield  (grid% hsurf,    ' hsurf    ')
            CALL prfield  (grid% sso_stdh, ' sso_stdh ')
            CALL prfield  (grid% sso_gamma,' sso_gamma')
            CALL prfield  (grid% sso_theta,' sso_theta')
            CALL prfield  (grid% sso_sigma,' sso_sigma')
            CALL prfield  (grid% emis_rad, ' emis_rad ')
            CALL prfield  (grid% prs_min,  ' prs_min  ')
            CALL prfield  (grid% skc,      ' skc      ')
!           if (grid% ivctype == IVCTYPE_ICON) &
            CALL prfield  (grid% hhl  ,    ' hhl      ')
!           CALL prfield  (grid% p0   ,    ' p0       ')
!           CALL prfield  (grid% dp0  ,    ' dp0      ')
!           CALL prfield  (grid% rho0 ,    ' rho0     ')
            if (grid% levtyp  /= WMO3_GENV) then
            !---------------------------------------
            ! Print estimate of model level pressure
            !---------------------------------------
            ps_ref = p_h_usstd (0.0_wp)
            DO k=1, grid% nz
              select case (grid% levtyp)
              case default
                pf = grid% akf(k) + grid% bkf(k)*ps_ref
              case (WMO3_HHYBRID)       ! Estimate height-based hybrid levels
                pf = p_h_usstd (grid% akf(k))
              end select
              WRITE(iu,"(a,1x,'k=',i3,' ak=',f11.3, ' bk=',f9.6,&
                                    &' akf=',f11.3,' bkf=',f9.6,' pf=',f9.1)")&
                 ci(:ni), k, grid% ak(k), grid% bk(k), grid% akf(k), &
                 grid% bkf(k), pf
            END DO
            if (grid% nz > 0) then
              k = grid% nz + 1
              WRITE(iu,"(a,1x,'k=',i3,' ak=',f11.3, ' bk=',f9.6)")&
                 ci(:ni), k,  grid%ak(k), grid% bk(k)
            end if
            end if
            DO k=1, grid% ns
              WRITE(iu,"(a,1x,'dbs    : k=',i2,i8)") ci(:ni), k, grid% dbs(k)
            END DO
      ENDIF
    ENDIF
    !---------------------------------------
    ! print 3-d fields, distributed over PEs
    !---------------------------------------
    call print (grid% m(HHL:), iunit=iu, intend=ci(:ni),verbose=verbose)
  CONTAINS
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !---------------------------------------------
    ! routines to print minimum and maximum values
    ! of array components
    !---------------------------------------------
    SUBROUTINE prfield (field,string,fak,dim)
    REAL(wp)          ,POINTER              :: field(:,:,:,:)
    CHARACTER (len=*) ,INTENT(in)           :: string
    REAL(wp)          ,INTENT(in) ,OPTIONAL :: fak
    CHARACTER (len=*) ,INTENT(in) ,OPTIONAL :: dim
    !----------------------------------
    ! print minimum,maximum of 4d array
    !----------------------------------
      CHARACTER(len=8) :: d
      real(wp)         :: minv, maxv
      d=''; IF(PRESENT(dim))d=dim
      IF (ASSOCIATED(field)) THEN
        minv = MINVAL(field)
        maxv = MAXVAL(field)
        IF(PRESENT(fak)) THEN
          WRITE (iu,"(a,a,a,2g10.3,' | ',2g10.3,1x,a)") ci(:ni), string ,&
                ': min,max=',minv,maxv,fak*minv,fak*maxv,trim(d)
        ELSE
          WRITE (iu,"(a,a,a,2g10.3,1x,a)") ci(:ni), string ,&
                ': min,max=',minv,maxv,trim(d)
        ENDIF
      ELSE
        WRITE(iu,"(a,a,a)")   ci(:ni), string ,': not associated!'
      ENDIF
    END SUBROUTINE prfield
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    SUBROUTINE prfield1 (field,string)
    REAL(wp), POINTER :: field(:)
    CHARACTER (len=*) :: string
    !----------------------------------
    ! print minimum,maximum of 1d array
    !----------------------------------
      IF (ASSOCIATED(field)) THEN
        WRITE (iu,"(a,a,a,2g10.3)") ci(:ni), string ,&
              ': min,max=',MINVAL(field),MAXVAL(field)
      ELSE
        WRITE(iu,"(a,a,a)")   ci(:ni), string ,': not associated!'
      ENDIF
    END SUBROUTINE prfield1
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    SUBROUTINE prfield1n (field,string)
    REAL(wp), POINTER :: field(:)
    CHARACTER (len=*) :: string
    !-----------------------------------------
    ! print first and last element of 1d array
    !-----------------------------------------
      IF (ASSOCIATED(field)) THEN
        WRITE (iu,"(a,a,a,2g10.3)") ci(:ni), string ,&
              ': (1),(n)=',field(1),field(SIZE(field))
      ELSE
        WRITE(iu,"(a,a,a)")   ci(:ni), string ,': not associated!'
      ENDIF
    END SUBROUTINE prfield1n
  END SUBROUTINE print_atm_grid
!==============================================================================
  SUBROUTINE lin_intpol (ingrid, ogrid, in, out)
  !------------------------------------------------------
  ! perform linear interpolation of 2D (horizontal) field
  ! (Gaussian or lat-lon non rotated grids so far)
  !------------------------------------------------------
  TYPE (t_grid) ,INTENT(in)  :: ingrid    ! grid description of input field
  TYPE (t_grid) ,INTENT(in)  ::  ogrid    ! grid description of output field
  REAL(wp)      ,INTENT(in)  :: in (:,:)  ! input field
  REAL(wp)      ,INTENT(out) :: out(:,:)  ! output field

    INTEGER  :: ix(ogrid%nx), iy(ogrid%ny)
    REAL(wp) :: rx(ogrid%nx), ry(ogrid%ny)
    INTEGER  :: i,j
    !----------------------------
    ! check for non-rotated grids
    !----------------------------
    IF(ogrid% rot .OR. ingrid% rot) &
      CALL finish('lin_intpol','cannot handle rotated grids')
    !-----------------------
    ! determine coefficients
    !-----------------------
    ix = 1
    CALL hunt (ingrid%dlon, ogrid%dlon, ix)
    iy = 1
    CALL hunt (ingrid%dlat, ogrid%dlat, iy)
    ix = MAX(1, MIN(ingrid%nx, ix))
    iy = MAX(1, MIN(ingrid%ny, iy))
    rx = (ogrid%dlon        - ingrid%dlon(ix)) &
       / (ingrid%dlon(ix+1) - ingrid%dlon(ix))
    ry = (ogrid%dlat        - ingrid%dlat(iy)) &
       / (ingrid%dlat(iy+1) - ingrid%dlat(iy))
    !------------
    ! interpolate
    !------------
    DO j=1,ogrid% ny
      DO i=1,ogrid% nx
        out(i,j) = (1._wp - rx(i)) * (1._wp - ry(j)) * in(ix(i)  ,iy(j))   &
                 +          rx(i)  * (1._wp - ry(j)) * in(ix(i)+1,iy(j))   &
                 + (1._wp - rx(i)) *          ry(j)  * in(ix(i)  ,iy(j)+1) &
                 +          rx(i)  *          ry(j)  * in(ix(i)+1,iy(j)+1)
      END DO
    END DO
  END SUBROUTINE lin_intpol
!------------------------------------------------------------------------------
  SUBROUTINE average_atm (ingrid, ogrid, in, out, n, ixy, jxy, nmin, invalid)
  !-------------------------------------------------------
  ! average field to coarser grid (nearest neighbour)
  ! (Gaussian and non rotated lat-lon grids so far)
  !-------------------------------------------------------
  TYPE (t_grid)      ,INTENT(in)  :: ingrid   ! grid of source field
  TYPE (t_grid)      ,INTENT(in)  ::  ogrid   ! grid of resulting field
  REAL(wp)           ,INTENT(in)  :: in (:,:) ! source field
  REAL(wp) ,OPTIONAL ,INTENT(out) :: out(:,:) ! resulting field
  INTEGER  ,OPTIONAL ,INTENT(out) :: n  (:,:) ! entries
  INTEGER  ,OPTIONAL ,INTENT(out) :: ixy(:,:) ! 1st index
  INTEGER  ,OPTIONAL ,INTENT(out) :: jxy(:,:) ! 2nd index
  INTEGER  ,OPTIONAL ,INTENT(in)  :: nmin     ! min number of entries (def=1)
  REAL(wp) ,OPTIONAL ,INTENT(in)  :: invalid  ! invalid value for n<nmin (0.)

    INTEGER  :: ix   (ingrid%nx),   iy  (ingrid%ny)
    REAL(wp) :: rx   (ingrid%nx),   ry  (ingrid%ny)
    REAL(wp) :: dlon ( ogrid%nx+2), dlat( ogrid%ny+2)
    REAL(wp) :: inv
    REAL(wp) :: r
    INTEGER  :: m(ogrid%nx,ogrid%ny)
    REAL(wp) :: o(ogrid%nx,ogrid%ny)
    INTEGER  :: i, j, nx, ny, nm, loc(2)
    LOGICAL  :: rot
    inv = -HUGE(inv); IF(PRESENT(invalid)) inv = invalid
    nm  = 1;          IF(PRESENT(nmin))    nm  = nmin
    rot = ingrid% rot .OR. ogrid% rot
    IF(PRESENT(ixy)) ixy = 0
    IF(PRESENT(jxy)) jxy = 0
    !----------------------------------------
    ! determine indices for non rotated grids
    !----------------------------------------
    IF(.NOT. ogrid% rot) THEN
      !-----------------
      ! 1 gridpoint halo
      !-----------------
      nx = ogrid%nx
      ny = ogrid%ny
      dlon (1)      = 2*ogrid% dlon(1)    - ogrid%dlon(2)
      dlon (2:nx+1) =   ogrid% dlon
      dlon (  nx+2) = 2*ogrid% dlon(  nx) - ogrid%dlon(nx-1)
      dlat (1)      = 2*ogrid% dlat(1)    - ogrid%dlat(2)
      dlat (2:ny+1) =   ogrid% dlat
      dlat (  ny+2) = 2*ogrid% dlat(  ny) - ogrid%dlat(ny-1)
    ENDIF
    IF(.NOT. rot) THEN
      !-----------------------
      ! determine coefficients
      !-----------------------
      ix = 1
      iy = 1
      CALL hunt (dlon, ingrid%dlon, ix)
      CALL hunt (dlat, ingrid%dlat, iy)
      WHERE (ix == nx+2) ix = 0
      WHERE (iy == ny+2) iy = 0
      rx = 0._wp
      ry = 0._wp
      WHERE (ix /= 0)
        rx = (ingrid%dlon - dlon(ix)) / (dlon(ix+1) - dlon(ix))
        WHERE (rx>0.5) ix=ix+1
      endwhere
      WHERE (iy /= 0)
        ry = (ingrid%dlat - dlat(iy)) / (dlat(iy+1) - dlat(iy))
        WHERE (ry>0.5) iy=iy+1
      endwhere
      !----------------------------------------------------------
      ! neglect halo, consider cyclic longitudes for global model
      !----------------------------------------------------------
      WHERE (iy == 1 .OR. iy == ny+2) iy = 0
      IF (ogrid%ngl /= 0) THEN
        WHERE (ix == 1)    ix = nx+1
        WHERE (ix == nx+2) ix = 2
      ELSE
        WHERE (ix == 1 .OR. ix == nx+2) ix = 0
      ENDIF
      ix = ix - 1
      iy = iy - 1
    ENDIF
    !--------
    ! average
    !--------
    o   = 0._wp
    m   = 0
    DO j=1,ingrid% ny
      DO i=1,ingrid% nx
        !------------------------------------
        ! determine indices for rotated grids
        ! no treatment of boundaries so far
        ! no diamonds so far
        !------------------------------------
        IF (ogrid% rot) THEN
          loc = MINLOC(SQRT( (ingrid%rlon(i,j,1,1)-ogrid%rlon(:,:,1,1))**2 &
                            +(ingrid%rlat(i,j,1,1)-ogrid%rlat(:,:,1,1))**2))
          ix(i)=loc(1)
          iy(j)=loc(2)
        ELSE IF (ingrid% rot) THEN
          CALL hunt (dlon, ingrid%rlon(i,j,1,1)*r2d, ix(i))
          CALL hunt (dlat, ingrid%rlat(i,j,1,1)*r2d, iy(j))
          IF (ix(i) == nx+2) ix(i) = 0
          IF (iy(j) == ny+2) iy(j) = 0
          r = 0._wp
          IF (ix(i) /= 0) THEN
            r = (ingrid%rlon(i,j,1,1)*r2d - dlon(ix(i))) / &
                (dlon(ix(i)+1)-dlon(ix(i)))
            IF (r>0.5) ix(i)=ix(i)+1
          ENDIF
          r = 0._wp
          IF (iy(j) /= 0) THEN
            r = (ingrid%rlat(i,j,1,1)*r2d - dlat(iy(j))) / &
                (dlat(iy(j)+1)-dlat(iy(j)))
            IF (r>0.5) iy(j)=iy(j)+1
          ENDIF
          !----------------------------------------------------------
          ! neglect halo, consider cyclic longitudes for global model
          !----------------------------------------------------------
          IF (iy(j) == 1 .OR. iy(j) == ny+2) iy(j) = 0
          IF (ogrid%ngl /= 0) THEN
            IF (ix(i) == 1)    ix(i) = nx+1
            IF (ix(i) == nx+2) ix(i) = 2
          ELSE
            IF (ix(i) == 1 .OR. ix(i) == nx+2) ix(i) = 0
          ENDIF
          ix(i) = ix(i) - 1
          iy(j) = iy(j) - 1
         ENDIF
         !------------
         ! now average
         !------------
         IF(PRESENT(ixy)) ixy(i,j) = ix(i)
         IF(PRESENT(jxy)) jxy(i,j) = iy(j)
         IF(ix(i)>0 .AND. iy(j)>0 .AND. in(i,j)/=inv ) THEN
           o  (ix(i),iy(j)) = o  (ix(i),iy(j)) + in(i,j)
           m  (ix(i),iy(j)) = m  (ix(i),iy(j)) + 1
         ENDIF
      END DO
    END DO
    !---------------
    ! invalid values
    !---------------
    WHERE (m>=nm)
      o   = o   / m
    ELSEWHERE
      o   = inv
    endwhere
    !---------------------------
    ! optional output parameters
    !---------------------------
    IF(PRESENT(n))     n = m
    IF(PRESENT(out)) out = o
  END SUBROUTINE average_atm
!==============================================================================
  SUBROUTINE lgt_grid (grid, nn)
  !-----------------------------------------
  ! perform Legendre transformation on geosp
  !-----------------------------------------
  TYPE (t_grid) ,INTENT(inout) TARGET    :: grid
  INTEGER       ,INTENT(in)    ,OPTIONAL :: nn
    grid% m(GEO_S)       = grid% m(GEOSP)
    grid% m(GEO_S)%i% nn = grid% nn
    CALL lgtd (grid%m(GEO_S), nn=nn)
    grid% geo_sh => grid%m(GEO_S)% ptr (:,1,1,1)
  END SUBROUTINE lgt_grid
!------------------------------------------------------------------------------
  SUBROUTINE lgti_grid (grid, nx, ny)
  !-----------------------------------------
  ! perform inverse Legendre transformation
  ! interpolate to grid point space
  !-----------------------------------------
  TYPE (t_grid) ,INTENT(inout) TARGET    :: grid
  INTEGER       ,INTENT(in)    ,OPTIONAL :: nx, ny
    INTEGER :: i
    LOGICAL :: flip
    flip = grid% dlat(1) < grid% dlat(grid% ny)
    IF(PRESENT(nx)) grid% nx = nx
    IF(PRESENT(ny)) grid% ny = ny
    grid% di  = 360._wp / nx
    grid% dj  = 180._wp / ny
    grid% ngl = ny
    DEALLOCATE (grid% dlon); ALLOCATE (grid% dlon (grid% nx))
    DEALLOCATE (grid% dlat); ALLOCATE (grid% dlat (grid% ny))
    CALL gauaw (ga = grid% dlat)
    grid% dlat = r2d * ASIN (grid% dlat)
    grid% dlon = (/(grid% di*i, i=0, grid% nx-1)/)
    IF(flip) grid% dlat = - grid% dlat
    grid% la1 = grid% dlat (1)
    grid% lo1 = grid% dlon (1)
    grid%m%i% lb(1) = 1
    grid%m%i% lb(2) = 1
    grid%m%i% ub(1) = grid% nx
    grid%m%i% ub(2) = grid% ny
    CALL update (grid% m)
    CALL set_pointers(grid)
#if defined(__ibm__)
    ! Work around xlf 9.1 bug with bounds-checking and spread:
    grid% rlon(:,:,1,1) = SPREAD (grid% dlon, 2, grid% ny) * d2r
    grid% rlat(:,:,1,1) = SPREAD (grid% dlat, 1, grid% nx) * d2r
#else
    grid% rlon(:,:,1,1) = SPREAD (grid% dlon * d2r, 2, grid% ny)
    grid% rlat(:,:,1,1) = SPREAD (grid% dlat * d2r, 1, grid% nx)
#endif
    grid% m(GEOSP)  = grid% m(GEO_S)
    CALL lgti (grid%m(GEOSP))
    grid% geosp => grid%m(GEOSP)% ptr
  END SUBROUTINE lgti_grid
!==============================================================================
  SUBROUTINE hunt_grid (grid, x, y, ix, iy)
  !------------------------------------------------------------
  ! return indices to grid
  ! for regular lat-lon or Gaussian grids
  !   ix, iy on entry: first guess
  !          on exit:  indices of gridpoint east-south to (x,y)
  !                    <0: not found
  !------------------------------------------------------------
  TYPE(t_grid) ,INTENT(in)    :: grid ! grid information
  REAL(wp)     ,INTENT(in)    :: x    ! x coordinate value
  REAL(wp)     ,INTENT(in)    :: y    ! y coordinate value
  INTEGER      ,INTENT(inout) :: ix   ! x direction index
  INTEGER      ,INTENT(inout) :: iy   ! y direction index

    INTEGER  :: nx
    INTEGER  :: ny
    REAL(wp) :: lon, lon1, lon2
    nx = grid% nx
    ny = grid% ny
    IF (grid% rot) THEN
      WRITE (0,*) 'hunt_grid: cannot handle rotated grid!'
      CALL finish( 'hunt_grid','cannot handle rotated grid!')
    ELSE
      !--------------------------
      ! care about cyclic bc in x
      !--------------------------
      lon2 = (grid% dlon(nx) + grid% dlon(1) + 360._wp) * 0.5_wp
      lon1 = lon2 - 360._wp
      lon  = x
      IF (lon < lon1) lon = lon + 360._wp
      IF (lon > lon2) lon = lon - 360._wp
      !------------
      ! first guess
      !------------
      IF (ix < 1 .OR. ix > nx) ix = 1 + int ((lon-grid%dlon(1)) &
                                 / (grid%dlon(nx)-grid%dlon(1)) )
      IF (iy < 1 .OR. iy > ny) iy = 1 + int ((  y-grid%dlat(1)) &
                                 / (grid%dlat(ny)-grid%dlat(1)) )
      !----------
      ! get index
      !----------
      CALL hunt (grid% dlon,lon,ix)
      CALL hunt (grid% dlat,  y,iy)
    ENDIF
  END SUBROUTINE hunt_grid
!------------------------------------------------------------------------------
  SUBROUTINE w_intpol (grid, x, y, n, ix, iy, id, w)
  !-----------------------------------------------------------------
  ! return indices and interpolation weights to grid
  !   ix, iy on entry: first guess
  !          on exit:  indices of gridpoint east-south to (x,y)
  !   n      on exit:  >0 OK   : number of gridpoints required
  !                    <0 ERROR: number required < array sizes
  !                     0 ERROR: interpolation impossible, at bounds
  ! linear interpolation so far
  !-----------------------------------------------------------------
  TYPE(t_grid) ,INTENT(in)    :: grid  ! grid information
  REAL(wp)     ,INTENT(in)    :: x     ! x coordinate value
  REAL(wp)     ,INTENT(in)    :: y     ! y coordinate value
  INTEGER      ,INTENT(out)   :: n     ! number of coefficients required
  INTEGER      ,INTENT(inout) :: ix(:) ! x direction indices
  INTEGER      ,INTENT(inout) :: iy(:) ! y direction indices
  INTEGER      ,INTENT(inout) :: id(:) ! diamond indices
  REAL(wp)     ,INTENT(out)   :: w (:) ! weights
    !----------------
    ! local variables
    !----------------
    INTEGER  :: nx   ! array bound in x direction
    INTEGER  :: ny   ! array bound in y direction
    INTEGER  :: i, j ! indices
    INTEGER  :: il   ! excess gridpoints in E-S direction
    INTEGER  :: iu   ! excess gridpoints in N-W direction
    REAL(wp) :: lon  ! longitude
    REAL(wp) :: dx,dy! grid spacing
    REAL(wp) :: w12 ,w21, w13, w31
    id = 1 ! ++ diamond index
    nx = grid% nx
    ny = grid% ny
    i = ix(1)
    j = iy(1)
    CALL hunt (grid, x, y, i, j)
    !---------------------------------------------
    ! coefficients for linear interpolation so far
    !---------------------------------------------
    il = 0
    iu = 1
    IF (j-il<1 .OR.  j+iu>ny) THEN
      !-----------------
      ! skip pole so far
      !-----------------
      n = 0
    ELSE
      n = 4
      IF (SIZE(ix)<n) THEN
        !--------------------------
        ! array arguments too small
        !--------------------------
        n = -n
      ELSE
        !------------
        ! get indices
        !------------
        ix(1:n) = (/i, i+1, i  , i+1/)
        iy(1:n) = (/j, j  , j+1, j+1/)
        IF (i-il<1 .OR. i+iu>nx) THEN
          !-----------------------------
          ! correct for east-west bounds
          !-----------------------------
          IF(.NOT. grid% cyc_x) THEN
            n = 0
          ELSE
            WHERE (ix>nx) ix = ix - nx
            WHERE (ix<1 ) ix = ix + nx
          ENDIF
        ENDIF
        IF (n == 4) THEN
          !---------------------------------
          ! weights for linear interpolation
          !---------------------------------
          dy  = grid% dlat(iy(3)) - grid% dlat(iy(1))
          dx  = grid% dlon(ix(2)) - grid% dlon(ix(1))
          IF (dx < 0._wp) dx =   dx + 360._wp ! cyclic boundary condition
          lon = x
          IF(lon<grid% dlon(ix(1))) lon = lon + 360._wp
          w21  = (lon - grid% dlon(ix(1))) / dx
          w12  = 1._wp - w21
          w31  = (y   - grid% dlat(iy(1))) / dy
          w13  = 1._wp - w31
          w(1) = w12 * w13
          w(2) = w21 * w13
          w(3) = w12 * w31
          w(4) = w21 * w31
        ENDIF
      ENDIF
    ENDIF
  END SUBROUTINE w_intpol
!==============================================================================
  SUBROUTINE set_geoid (grid)
  TYPE (t_grid) ,INTENT(inout) :: grid

  character (len=256) :: geoid_path_file
  logical             :: exist

  geoid_path_file = path_file(data, trim(geoid_file))

  ! Check if "geoid_file" exists
  if (dace% lpio) then
     inquire(FILE=geoid_path_file, EXIST=exist)
     if (.not. exist) then
        call finish('set_geoid',                                   &
                    'geoid file not found: '//trim(geoid_path_file))
     end if
  end if

  ! Call input routine for given geoid format
  geoid_format = tolower(geoid_format)
  select case (geoid_format)
  case ('ascii')
     call set_geoid_EGM96 (grid, geoid_path_file)
  case ('nc', 'netcdf')
     call set_geoid_6c3stat (grid, geoid_path_file)
  case default
     call finish('set_geoid',                                 &
                 'geoid format not supported: '//trim(geoid_format))
  end select

  END SUBROUTINE set_geoid
!==============================================================================
  SUBROUTINE set_geoid_EGM96 (grid, geoid_path_file)
  TYPE (t_grid) ,INTENT(inout)  :: grid
  CHARACTER (len=*), INTENT(in) :: geoid_path_file
  !----------------------------------------------------------------------------
  !    read geoid heights
  !
  ! grid file structure:
  !
  !  This FORTRAN program will interpolate values from a WORLD-WIDE grid file.
  !  The grid file must have a header as record 1 with the following format:
  !
  !       south    north    west    east    spacing N-S   spacing W-E
  !
  !  The grid data are in subsequent records arranged from North to South, West
  !  to East, i.e. record 2 would be the northern most latitude band of data.
  !
  !   -90.000000   90.000000     .000000  360.000000     .250000     .250000
  !
  !    13.606   13.606   13.606   13.606   13.606   13.606   13.606   13.606
  !    13.606   13.606   13.606   13.606   13.606   13.606   13.606   13.606
  !    13.606   13.606   13.606   13.606   13.606   13.606   13.606   13.606
  !----------------------------------------------------------------------------

    REAL(wp) :: lonf, lonl, latf, latl, lond, latd

    INTEGER, PARAMETER    :: nlat =  721
    INTEGER, PARAMETER    :: nlon = 1441
    INTEGER               :: i,j,l, ig,jg
    INTEGER               :: iu, ios
    REAL(wp) ,allocatable :: geoid (:,:)
    REAL(wp) ,POINTER     :: g (:,:,:,:), rlon (:,:,:,:), rlat (:,:,:,:)
    REAL(wp)              :: wi, wj, lon

    IF (ASSOCIATED(grid%geoid)) RETURN

    allocate (geoid (nlon,nlat))
    !--------------------------------
    ! Read grid data on I/O processor
    !--------------------------------
    ios = -1
    if (dace% lpio) then
       iu = get_unit_number ()
       OPEN (iu, file=geoid_path_file,    &
                 status='old', action='read', iostat=ios)
       if (ios == 0) then
          READ (iu,*) latf, latl, lonf, lonl, latd, lond

          DO j=1,nlat
             READ (iu,*) geoid (:,nlat-j+1)
          END DO

          CLOSE (iu)
          CALL return_unit_number (iu)
       end if
    end if
    call p_bcast (ios,   dace% pio)
    if (ios /= 0) call finish('set_geoid_EGM96',                           &
                              'cannot open '//trim(geoid_path_file))

    !--------------------------------------
    ! Broadcast to other processor elements
    !--------------------------------------
    call p_bcast (latf,  dace% pio)
    call p_bcast (latl,  dace% pio)
    call p_bcast (lonf,  dace% pio)
    call p_bcast (lonl,  dace% pio)
    call p_bcast (latd,  dace% pio)
    call p_bcast (lond,  dace% pio)
    call p_bcast (geoid, dace% pio)

    CALL ALLOCATE (grid, 'geoid')

    g    => grid% geoid
    rlon => grid% rlon
    rlat => grid% rlat

    latd =  latd * pi / 180._wp
    lond =  lond * pi / 180._wp
    latf =  latf * pi / 180._wp
    lonf =  lonf * pi / 180._wp

    DO l     = LBOUND(g,4),UBOUND(g,4)
      DO j   = LBOUND(g,2),UBOUND(g,2)
        DO i = LBOUND(g,1),UBOUND(g,1)
          lon = rlon(i,j,1,l)
          IF(lon<0._wp) lon=lon+2._wp*pi
          wi = (lon-lonf) / lond
          ig = INT(wi)
          wi = wi - ig
          ig = ig + 1
          IF (ig==nlon) ig = 1
          wj = (rlat(i,j,1,l)-latf) / latd
          jg = INT(wj)
          wj = wj - jg
          jg = jg + 1
          IF (jg==nlat) jg = 1
          g(i,j,1,l) = (1._wp-wi) * (1._wp-wj) * geoid (ig  ,jg)   &
                     +        wi  * (1._wp-wj) * geoid (ig+1,jg)   &
                     + (1._wp-wi) *        wj  * geoid (ig  ,jg+1) &
                     +        wi  *        wj  * geoid (ig+1,jg+1)
        END DO
      END DO
   END DO

   if (dace% lpio) then
      write(*,'(a)') 'set_geoid: EGM96 geoid read'
      !call print_geoid (grid, 'set_geoid_EGM96')
   endif

  END SUBROUTINE set_geoid_EGM96
!------------------------------------------------------------------------------
  SUBROUTINE set_geoid_6c3stat (grid, geoid_path_file)
  TYPE (t_grid) ,INTENT(inout)  :: grid
  CHARACTER (len=*), INTENT(in) :: geoid_path_file
  !----------------------------------------------------------------------
  ! Read geoid undulation from netCDF file
  !
  ! The eigen-6c3stat is originally provided as a large ASCII file.
  ! Use "geoid2netcdf" to convert the ASCII file to netCDF and to select
  ! the desired region. Reading the ASCII file is not supported.
  !----------------------------------------------------------------------

  INTEGER             :: i, j, l, ig, jg
  REAL(wp)            :: wi, wj, lon
  REAL(wp)            :: lonf, latf, lond, latd
  REAL(wp) ,POINTER   :: g (:,:,:,:), rlon (:,:,:,:), rlat (:,:,:,:)

  ! netCDF
  INTEGER :: ncid, ncret, ierr         ! netCDF ID, return code, error code
  INTEGER :: LonDimID, LatDimID        ! netCDF dimension IDs
  INTEGER :: LonVarID, LatVarID        ! netCDF variable IDs
  INTEGER :: GeoidVarID

  CHARACTER (len=40)      :: DimName   ! name of netCDF dimension
  INTEGER, DIMENSION(1:2) :: ncstart   ! NF90_Get_Var start indices
  INTEGER, DIMENSION(1:2) :: nccount   ! NF90_Get_Var number of variables

  CHARACTER (len=40)      :: GeoidName          ! netCDF attributes "modelname"
  REAL (kind=wp)          :: LatNorth, LatSouth ! netCDF attributes
  REAL (kind=wp)          :: LonWest, LonEast   ! netCDF attributes
  REAL (kind=wp)          :: GridStep           ! netCDF attributes
  INTEGER                 :: Nlon, Nlat, Npoint ! netCDF attributes
  REAL (kind=sp), dimension(:), pointer :: &
                   RefLon,                 & ! geoid longitudes, from netCDF
                   RefLat                    ! geoid latitudes, from netCDF
  REAL (kind=sp), dimension(:,:), allocatable :: &
                   ncgeoid                   ! geoid undulation from netCDF
  REAL (kind=wp), dimension(:,:,:), pointer :: &
                   tmpgeoid                  ! geoid on full grid,
                                             ! for interpolation

  INTEGER :: MPIopt = 1   ! MPIopt = 1 - broadcast netCDF geoid and
                          !              interpolate on each processor
                          ! MPIopt = 2 - interpolate on I/O processor and
                          !              scatter sub-domains to processors

  !---------------------------------
  ! Read geoid data on I/O processor
  !---------------------------------
  if (dace% lpio) then

     ! Open netCDF file read only
     ncret = NF90_Open(TRIM(geoid_path_file), NF90_NoWrite, ncid)
     IF (ncret /= 0) THEN
        ! Error opening geoid file
        call finish('set_geoid_6c3stat',                     &
                    'cannot open '//trim(geoid_path_file))
     END IF

   ! Get global attributes
   ncret = nf90_get_att(ncid, NF90_GLOBAL, 'modelname', GeoidName)
   ncret = nf90_get_att(ncid, NF90_GLOBAL, 'latlimit_nort', LatNorth)
   ncret = nf90_get_att(ncid, NF90_GLOBAL, 'latlimit_south', LatSouth)
   ncret = nf90_get_att(ncid, NF90_GLOBAL, 'longlimit_west', LonWest)
   ncret = nf90_get_att(ncid, NF90_GLOBAL, 'longlimit_east', LonEast)
   ncret = nf90_get_att(ncid, NF90_GLOBAL, 'gridstep', GridStep)
   ncret = nf90_get_att(ncid, NF90_GLOBAL, 'latitude_parallels', Nlat)
   ncret = nf90_get_att(ncid, NF90_GLOBAL, 'longitude_parallels', Nlon)
   ncret = nf90_get_att(ncid, NF90_GLOBAL, 'number_of_gridpoints', Npoint)

   ! inquire dimensions
   ncret = NF90_Inq_DimID(ncid, 'londim', LonDimID)
   ncret = NF90_Inq_DimID(ncid, 'latdim', LatDimID)
   ncret = NF90_Inquire_Dimension(ncid, LonDimID, DimName, Nlon)
   ncret = NF90_Inquire_Dimension(ncid, LatDimID, DimName, Nlat)

   ! inquire variable IDs
   ncret = NF90_Inq_VarID(ncid, 'geoid', GeoidVarID)
   ncret = NF90_Inq_VarID(ncid, 'lon', LonVarID)
   ncret = NF90_Inq_VarID(ncid, 'lat', LatVarID)

   ! Read longitudes and latitudes
   allocate( RefLon(1:Nlon) )
   allocate( RefLat(1:Nlat) )
   ncret = NF90_Get_Var(ncid, LonVarID, RefLon)
   ncret = NF90_Get_Var(ncid, LatVarID, RefLat)

   ! Read whole geoid
   allocate( ncgeoid(1:Nlon,1:Nlat) )
   ncstart = (/1, 1/)
   nccount = (/Nlon, Nlat/)
   ncret = NF90_Get_Var(ncid, GeoidVarID, ncgeoid,     &
                        start=ncstart, count=nccount)

   ! Close the netCDF file.
   ncret = nf90_close(ncid)

  end if   !   if (dace% lpio) then

  !-------------------------------
  ! Distribute geoid and meta-data
  !-------------------------------
  call p_bcast (ncret,   dace% pio)
  if (ncret /= 0) call finish('set_geoid_6c3stat',                    &
                       'cannot open '//trim(geoid_path_file))

  call p_bcast (LatSouth, dace% pio)
  call p_bcast (LatNorth, dace% pio)
  call p_bcast (LonWest,  dace% pio)
  call p_bcast (LonEast,  dace% pio)
  call p_bcast (GridStep, dace% pio)
  call p_bcast (Nlat,     dace% pio)
  call p_bcast (Nlon,     dace% pio)

  if (MPIopt == 1) then
     ! Broadcast entire netCDF geoid and interpolate on each processor

     if (.not. dace% lpio) then
        allocate( ncgeoid(1:Nlon,1:Nlat) )
     end if

     call p_bcast (ncgeoid,  dace% pio)

     !---------------------------
     ! Interpolate geoid on grid
     !---------------------------
     CALL ALLOCATE (grid, 'geoid')

     grid% geoid = -999.0

     g    => grid% geoid
     rlon => grid% rlon
     rlat => grid% rlat

     latd =  GridStep * d2r
     lond =  latd
     latf =  LatSouth * d2r
     lonf =  LonWest  * d2r

     DO l     = LBOUND(g,4),UBOUND(g,4)
        DO j   = LBOUND(g,2),UBOUND(g,2)
           DO i = LBOUND(g,1),UBOUND(g,1)
              lon = rlon(i,j,1,l)
              IF(lon<0._wp) lon=lon+2._wp*pi
              wi = (lon-lonf) / lond
              ig = INT(wi)
              wi = wi - ig
              ig = ig + 1
              IF (ig==nlon) ig = 1
              wj = (rlat(i,j,1,l)-latf) / latd
              jg = INT(wj)
              wj = wj - jg
              jg = Nlat - jg !+ 1
              !IF (jg==nlat) jg = 1
              IF (jg==0) jg = 1
              IF (jg>=Nlat) jg = Nlat - 1
              g(i,j,1,l) = (1._wp-wi) * (1._wp-wj) * ncgeoid (ig  ,jg)   &
                         +        wi  * (1._wp-wj) * ncgeoid (ig+1,jg)   &
                         + (1._wp-wi) *        wj  * ncgeoid (ig  ,jg+1) &
                         +        wi  *        wj  * ncgeoid (ig+1,jg+1)
           END DO
        END DO
     END DO

     IF (ANY(g == -999.0)) THEN
        call finish('set_geoid_6c3stat',                            &
                    'geoid file does not cover entire model region')
     END IF
  end if   ! if (MPIopt == 1) then

  if (MPIopt == 2) then
     ! Interpolate on I/O processor and scatter sub-domains to processors

     allocate( tmpgeoid(grid%lbg(1):grid%ubg(1), grid%lbg(2):grid%ubg(2),   &
                        grid%lbg(4):grid%ubg(4))  )

     tmpgeoid = -999.0

     if (dace% lpio) then
        rlon => grid% rlon
        rlat => grid% rlat

        latd =  GridStep * d2r
        lond =  latd
        latf =  LatSouth * d2r
        lonf =  LonWest  * d2r

        DO l     = grid%lbg(4), grid%ubg(4)
           DO j   = grid%lbg(2), grid%ubg(2)
              DO i = grid%lbg(1), grid%ubg(1)

                 lon = rlon(i,j,1,l)
                 IF(lon<0._wp) lon=lon+2._wp*pi
                 wi = (lon-lonf) / lond
                 ig = INT(wi)
                 wi = wi - ig
                 ig = ig + 1
                 IF (ig==nlon) ig = 1
                 wj = (rlat(i,j,1,l)-latf) / latd
                 jg = INT(wj)
                 wj = wj - jg
                 jg = Nlat - jg !+ 1
                 !IF (jg==nlat) jg = 1
                 IF (jg==0) jg = 1
                 IF (jg>=Nlat) jg = Nlat - 1
                 tmpgeoid(i,j,l) =                                    &
                        (1._wp-wi) * (1._wp-wj) * ncgeoid (ig  ,jg)   &
                      +        wi  * (1._wp-wj) * ncgeoid (ig+1,jg)   &
                      + (1._wp-wi) *        wj  * ncgeoid (ig  ,jg+1) &
                      +        wi  *        wj  * ncgeoid (ig+1,jg+1)
              END DO
           END DO
        END DO

        IF (ANY(tmpgeoid == -999.0)) THEN
           call finish('set_geoid_6c3stat',                            &
                       'geoid file does not cover entire model region')
        END IF

        ! Deallocate netCDF fields before geoid is distributed
        deallocate( RefLon )
        deallocate( RefLat )
        deallocate( ncgeoid )

     end if     ! if (dace% lpio) then

     ! Allocate geoid with domain decomposition
     CALL ALLOCATE (grid, 'geoid')

     ! Copy "tmpgeoid" as long domain decomposition is not used for geoid
     call p_bcast (tmpgeoid,  dace% pio)
     DO l     = grid%lbg(4), grid%ubg(4)
        DO j   = grid%lbg(2), grid%ubg(2)
           DO i = grid%lbg(1), grid%ubg(1)
              grid%geoid(i,j,1,l) = tmpgeoid(i,j,l)
           END DO
        END DO
     END DO

     ! Distribute "tmpgeoid" regarding the domain decomposition
     !call scatter_level (tmpgeoid, grid%geoid(:,:,1,:), grid%dc, dace%pio)

  end if   ! if (MPIopt == 2) then

  if (dace% lpio) then
      write(*,'(a)') 'set_geoid: ' // trim(GeoidName) // ' geoid read'
     !call print_geoid (grid, 'set_geoid_6c3stat')
  end if

  END SUBROUTINE set_geoid_6c3stat
!------------------------------------------------------------------------------
  SUBROUTINE print_geoid (grid, msg)
    TYPE (t_grid) ,INTENT(inout) :: grid
    CHARACTER (*)                :: msg

    INTEGER                      :: iu, ios
    INTEGER                      :: i,j,l
    REAL (kind=wp), target, allocatable :: glogeoid(:,:,:,:)  ! global geoid
    REAL (kind=wp), pointer      :: gg(:,:,:,:)


    gg => grid%geoid
    if (size(grid%geoid) < size(grid%rlon)) then
       if (dace% lpio) then
          allocate( glogeoid(grid%lbg(1):grid%ubg(1), grid%lbg(2):grid%ubg(2), &
                                       1:1          , grid%lbg(4):grid%ubg(4)) )
       end if
       call gather_multi (grid%geoid, glogeoid, grid%dc, dace%pio)
       gg => glogeoid
    end if

  if (dace% lpio) then
     iu = get_unit_number ()
     OPEN (iu, file='geoid_grid.dat',    &
           status='unknown', action='write', iostat=ios)
     if (ios == 0) then
        write(iu,'(a)') '# Geoidundulation, interpolated on grid% geoid'
        write(iu,'(a)') '# ' // trim(msg)
        write(iu,'(a,4(i10,tr2))') '# shape(geoid) = ', shape(gg)
        write(iu,'(a)') '# '
        write(iu,'(a)') '# lon (deg.)  lat (deg.) geoid (m)'

        DO l     = LBOUND(gg,4),UBOUND(gg,4)
           DO j   = LBOUND(gg,2),UBOUND(gg,2)
              DO i = LBOUND(gg,1),UBOUND(gg,1)

                 write(iu,'(3(f9.3,tr2))') &
                      grid%rlon(i,j,1,l)*r2d, grid%rlat(i,j,1,l)*r2d, &
                      gg(i,j,1,l)

              END DO
           END DO
        END DO

        CLOSE (iu)
        CALL return_unit_number (iu)
     end if
  end if

  END SUBROUTINE print_geoid
!------------------------------------------------------------------------------
  SUBROUTINE coarsen_grid (y, x, stride)
  !----------------------------------------------------------------------------
  ! derive a grid of coarser resolution by taking only each stride.th gridpoint
  ! loop over invariant fields, pick up selected gridpoints
  !----------------------------------------------------------------------------
  TYPE (t_grid) ,INTENT(inout) :: y       ! output grid
  TYPE (t_grid) ,INTENT(in)    :: x       ! input grid
  INTEGER       ,INTENT(in)    :: stride  ! stride

    !----------------
    ! local variables
    !----------------
    INTEGER :: i
    !--------------------------------------
    ! loop over variables, get local bounds
    !--------------------------------------
    DO i=1, SIZE(x%m)
      IF (x%m(i)%i% alloc .AND. y%m(i)%i% alloc) THEN
        !--------
        ! coarsen
        !--------
        y%m(i)% ptr(:,:,:,:) = x%m(i)% ptr(::stride,::stride,:,:)
      ELSE
        !---------------------------
        ! mark variables not present
        !---------------------------
        y%m(i)%i% alloc = .FALSE.
      ENDIF
    END DO
    !--------
    ! cleanup
    !--------
    CALL update (y% m)
    CALL set_pointers (y)
  END SUBROUTINE coarsen_grid
!------------------------------------------------------------------------------
  subroutine cosmo_ref_atm (grid, lnew_hhl)
  type (t_grid) ,intent (inout)           :: grid
  logical       ,intent (in)    ,optional :: lnew_hhl ! flag to calculate hhl
  !-------------------------------------------
  ! set fields of COSMO reference atmosphere :
  !   hhl  : geometrical height of half levels
  !   p0   : base-state pressure
  !   rho0 : base-state density
  !   dp0  : base-state pressure thickness
  !-------------------------------------------
    !----------------
    ! local variables
    !----------------
    integer               :: ie,je,ke,iu,il,ju,jl !local indices
    integer               :: i, j, n
    integer               :: ierror
    real(wp), allocatable :: hsurfs(:,:,:)
    real(wp), allocatable :: t0    (:,:,:)
    real(wp), allocatable :: p0hl  (:,:,:)
    real(wp), allocatable :: t0hl  (:,:,:)
    logical               :: lanalyt_calc_t0p0, lnewhhl
    character (len=80)    :: yerrmsg
    type (refatm_type)    :: refatm
    type (vcoord_type)    :: vcoord

    lnewhhl = .true.; if (present (lnew_hhl)) lnewhhl = lnew_hhl

    iu = grid% ub(1)
    il = grid% lb(1)
    ju = grid% ub(2)
    jl = grid% lb(2)

    ie = iu - il + 1
    je = ju - jl + 1
    ke = grid% nz
    !-----------------------------------------------------------
    ! fallback to default vertical coordinates if not set so far
    !-----------------------------------------------------------
    if (grid% vc% ivctype < 1) then
      call set_vcoord_defaults
      i = -1
      do j = 1, size (vcoord_defaults)
        if (vcoord_defaults(j)% ivctype == cosmo_ivctype .and. &
            vcoord_defaults(j)% nlevels == ke+1                ) then
          i = j
          exit
        endif
      end do
      if (i < 1) then
        !--------------------------------------------------------
        ! no matching ivctype found, abort if new hhl is required
        !--------------------------------------------------------
        if (lnewhhl .or. cosmo_ivctype == 1) then
          write(0,*) "  grid% nz + 1  =", ke+1
          write(0,*) "  cosmo_ivctype =", cosmo_ivctype
          write(0,*) "vcoord_defaults: ivctype, =", &
                      vcoord_defaults (:)% ivctype
          write(0,*) "vcoord_defaults: levels   =", &
                      vcoord_defaults (:)% nlevels
          write(0,*) "lnew_hhl                  =",lnewhhl
          call finish ('cosmo_ref_atm','unknown ivctype')
        endif
      else
        !-------------------------------------
        ! matching ivctype found, copy entries
        !-------------------------------------
        n = vcoord_defaults (i)% nlevels
        allocate (grid% vc% vert_coord (n))
        allocate (grid% vc% sigm_coord (n))
        grid% vc% ivctype    = vcoord_defaults (i)% ivctype
        grid% vc% nlevels    = vcoord_defaults (i)% nlevels
        ! Do not overwrite given uuid!
!!      grid% vc% vc_uuid    = vcoord_defaults (i)% vc_uuid
        grid% vc% vcflat     = vcoord_defaults (i)% vcflat
        grid% vc% vert_coord = vcoord_defaults (i)% vert_coord (1:n)
        grid% vc% sigm_coord = vcoord_defaults (i)% sigm_coord (1:n)
        grid% vc% svc1       = cosmo_svc(1)
        grid% vc% svc2       = cosmo_svc(2)
        grid% vc% nfltvc     = cosmo_nfltvc
        call dealloc_vcoord_structure (ierror)
      endif
      !---------------------------------------------------------
      ! no setup of new hhl required, just copy namelist entries
      !---------------------------------------------------------
      if (i < 1) then
        grid% vc% ivctype    = cosmo_ivctype
        grid% vc% svc1       = cosmo_svc(1)
        grid% vc% svc2       = cosmo_svc(2)
        grid% vc% nfltvc     = cosmo_nfltvc
      endif
    endif
    !--------------------------------------------------------------------
    ! copy vertical coordinate parameters from DACE to COSMO derived type
    !--------------------------------------------------------------------
    allocate (vcoord% vcflat)
    vcoord% ivctype    = grid% vc% ivctype
    vcoord% nlevels    = grid% vc% nlevels
    vcoord% vc_uuid    = grid% vc% vc_uuid
    vcoord% vcflat     = grid% vc% vcflat
    if (lnewhhl) then
       allocate (vcoord% sigm_coord (size (grid% vc% sigm_coord)))
       allocate (vcoord% vert_coord (size (grid% vc% vert_coord)))
       vcoord% vert_coord = grid% vc% vert_coord
       vcoord% sigm_coord = grid% vc% sigm_coord
    end if
    !---------------------------------------------------
    ! allocate temporary reference atmosphere components
    !---------------------------------------------------
    allocate (refatm% p0sl)
    allocate (refatm% t0sl)
    allocate (refatm% dt0lp)
    allocate (refatm% delta_t)
    allocate (refatm% h_scal)
    allocate (refatm% bvref)
    !-----------------------------------------------------------
    ! fallback to default reference atmosphere if not set so far
    !------------------------------------------------------------
    if (grid% ra% irefatm < 1) then
      i = cosmo_refatm
      if (i < 1) call finish ('cosmo_ref_atm','unknown refatm')
      call set_refatm_defaults
      grid% ra% irefatm = i
      grid% ra% p0sl    = refatm_defaults (i)% p0sl
      grid% ra% t0sl    = refatm_defaults (i)% t0sl
      grid% ra% dt0lp   = refatm_defaults (i)% dt0lp
      grid% ra% delta_t = refatm_defaults (i)% delta_t
      grid% ra% h_scal  = refatm_defaults (i)% h_scal
      call dealloc_refatm_structure (ierror)
    endif
    !---------------------------------------------------------------------
    ! copy reference atmosphere parameters from DACE to COSMO derived type
    !---------------------------------------------------------------------
    refatm% irefatm = grid% ra% irefatm
    refatm% p0sl    = grid% ra% p0sl
    refatm% t0sl    = grid% ra% t0sl
    refatm% dt0lp   = grid% ra% dt0lp
    refatm% delta_t = grid% ra% delta_t
    refatm% h_scal  = grid% ra% h_scal
    !-----------------------------------------
    ! set fields of COSMO reference atmosphere
    !-----------------------------------------
    call allocate(grid, 'hhl')
    call allocate(grid, 'p0')
    call allocate(grid, 'rho0')
    call allocate(grid, 'dp0')
    !---------------------
    !allocate local arrays
    !---------------------
    allocate (t0    (ie,je,ke  ))
    allocate (p0hl  (ie,je,ke+1))
    allocate (t0hl  (ie,je,ke+1))
    lanalyt_calc_t0p0 = .false.
    !-----------------------------------------------
    ! split topography scales for sleeve coordinates
    !-----------------------------------------------
    if (lnew_hhl .and. (vcoord% ivctype == 3 &
                  .or.  vcoord% ivctype == 4)) then
      allocate (hsurfs(grid% lbg(1):grid% ubg(1),   &
                       grid% lbg(2):grid% ubg(2), 2)); hsurfs = 0
      call sleve_split_oro (   &!
        grid% hsurf (:,:,1,1), &! <- height of surface topography
              hsurfs,          &! -> split orography
        grid% nx,              &! <- dimensions of hsurf
        grid% ny,              &! <- dimensions of hsurf
        grid% vc% nfltvc,      &! <- number of filter applications
              0,               &! <- nextralines ?
        grid% vc% svc1,        &! <- decay rate for large-scale part
        grid% vc% svc2,        &! <- decay rate for small-scale to
        grid% vc% vcflat,      &! <- coord.where levels become flat
              6,               &! <- output unit
              dace% pe,            &! <- MPI PE
              ierror,          &! -> error return code
              yerrmsg          )! -> error return message
      if (ierror /= 0) call finish('sleve_split_oro',yerrmsg)
    else
      allocate (hsurfs(grid% lb (1):grid% ub (1),   &
                       grid% lb (2):grid% ub (2), 2)); hsurfs = 0
    endif
    !----------------------------------------------------
    ! calculate reference atmosphere (and optionally hhl)
    !----------------------------------------------------
    select case (refatm% irefatm)
    case (1)
!     if (dace% lpio) write(6,*) 'SR cosmo_ref_atm: call reference_atmosphere'
      call reference_atmosphere (     &!
        grid% hhl  (:,:,:,1),         &! <-> geometrical height of half levels
        grid% p0   (:,:,:,1),         &!  -> base-state  full level pressure
              p0hl,                   &!  -> base-state pressure at half levels
        grid% rho0 (:,:,:,1),         &!  -> at full level reference density
               t0,                    &!  -> full level reference temp.
               t0hl,                  &!  -> base-state half level temp.
        grid% dp0  (:,:,:,1),         &!  -> base-state layer thickness
        grid% hsurf (il:iu,jl:ju,1,1),&! <-  height of surface topography
              hsurfs(il:iu,jl:ju,:),  &! <-  split orography, sleve coordinates
              ie,                     &! <-  dimensions of the fields
              je,                     &! <-  dimensions of the fields
              grid% nz,               &! <-  dimensions of the fields
              refatm,                 &! <-  (contains dt0lp, p0sl, t0sl)
              vcoord,                 &! <-  -> only kflat returned
        grid% vc% svc1,               &! <-  decay rate large-scale topo part
        grid% vc% svc2,               &! <-  decay rate small-scale topo part
              R,                      &! <-  r, gas constant for dry air
              gacc,                   &! <-  g, gravity acceleration
              lanalyt_calc_t0p0,      &! <-  calculate t0 analytically
              lnewhhl,                &! <-  compute new hhl for GRIB2
              yerrmsg,                &!  ->
              ierror                  )!  ->
      if (ierror /= 0) call finish('reference_atmosphere',yerrmsg)
    case (2)
!     if (dace% lpio) &
!       write(6,*) 'SR cosmo_ref_atm: call reference_atmosphere 2'
      call reference_atmosphere_2 (   &!
        grid% hhl  (:,:,:,1),         &! <-> geometrical height of half levels
        grid% p0   (:,:,:,1),         &!  -> base-state pressure at full levels
              p0hl,                   &!  -> base-state pressure at half levels
        grid% rho0 (:,:,:,1),         &!  -> reference density at full levels
              t0,                     &!  -> full level reference temp.
              t0hl,                   &!  -> base-state half level temp.
        grid% dp0  (:,:,:,1),         &!  -> base-state layer thickness
        grid% hsurf (il:iu,jl:ju,1,1),&! <-  height of surface topography
              hsurfs(il:iu,jl:ju,:),  &! <-  split orography, sleve coordinates
              ie,                     &! <-  dimensions of the fields
              je,                     &! <-  dimensions of the fields
        grid% nz,                     &! <-  dimensions of the fields
              refatm,                 &! <-  (contains dt0lp, p0sl, t0sl)
              vcoord,                 &! <-  -> only kflat returned
        grid% vc% svc1,               &! <-  decay rate large-scale topo part
        grid% vc% svc2,               &! <-  decay rate small-scale topo part
              R,                      &! <-  r, gas constant for dry air
              gacc,                   &! <-  g, gravity acceleration
              lnewhhl,                &! <-  compute new hhl for GRIB2
              yerrmsg,                &!  ->
              ierror                  )!  ->
      if (ierror /= 0) call finish('reference_atmosphere_2',yerrmsg)
    end select
    !--------------------------
    ! derive model top pressure
    !--------------------------
    grid% ptopf = maxval (grid% p0 (:,:,1,:))
    grid% ptopf = p_max  (grid% ptopf, comm = grid% dc% comm)
!   !-----------------------------------
!   ! set hhl if required
!   ! taken from Cosmo 5.2 src_input
!   ! (so far kflat is not used in DACE)
!   !-----------------------------------
!   IF (lnewhhl) THEN
!      !---------------------------------------------------------
!      ! 1) set a missing value for initialization of the search:
!      !---------------------------------------------------------
!      vcoord%kflat = -1
!      !---------------------------------------------------------------------
!      ! 2) search the level from top to bottom, where for the first time the
!      !    minimum and maximum of the local hhl (local PE) are different.
!      !---------------------------------------------------------------------
!      kloop1: DO k = 2, ke+1
!         IF (MAXVAL(grid% hhl(:,:,k,1)) /= MINVAL(grid% hhl(:,:,k,1))) THEN
!            vcoord%kflat = k-1
!            EXIT kloop1
!         ENDIF
!      ENDDO kloop1
!      !-----------------------------------------------------------------
!      ! 3) In the global grid, the global kflat must be the maximum over
!      !    the local kflat's, because for PEs with only sea points
!      !    (=flat orography), the
!      !    above search did not find the correct kflat yet.
!      !-----------------------------------------------------------------
!      vcoord%kflat = p_max (vcoord% kflat)
!      !------------------------------------------------------------------
!      ! 4) If now the global kflat is still -1, this means that the whole
!      !    domain only contains sea points.
!      !    In this case, kflat has no real meaning but we set it to ke-1
!      !------------------------------------------------------------------
!      IF (vcoord% kflat == -1) THEN
!       vcoord%kflat = ke-1
!     END IF
!   END IF
    !---------------------------------------------------
    ! deallocate temporary COSMO derived type components
    !---------------------------------------------------
    deallocate (refatm% p0sl)
    deallocate (refatm% t0sl)
    deallocate (refatm% dt0lp)
    deallocate (refatm% delta_t)
    deallocate (refatm% h_scal)
    deallocate (refatm% bvref)

    if (lnewhhl) then
       deallocate (vcoord% vert_coord)
       deallocate (vcoord% sigm_coord)
    end if
    deallocate (vcoord% vcflat)

  END SUBROUTINE cosmo_ref_atm
!==============================================================================
!------------------------------------------------------------------------------
! Functions from COSMO-code for trafo geogr.<--> rotated coord.system
!------------------------------------------------------------------------------

  FUNCTION  phirot2phi ( phirot, rlarot, polphi, pollam, polgam )
  !--------------------------------------------------------------------------
  ! This function converts phi (latitude) from one rotated system to another
  ! system. If the optional argument polgam is present, the other system
  ! can also be a rotated one, where polgam is the angle between the two
  ! north poles.
  ! If polgam is not present, the other system is the real geographical
  ! system.
  !--------------------------------------------------------------------------
  REAL (KIND=wp), INTENT (IN)      ::        &
       polphi,   & ! latitude of the rotated north pole
       pollam,   & ! longitude of the rotated north pole
       phirot,   & ! latitude in the rotated system
       rlarot      ! longitude in the rotated system

  REAL (KIND=wp), INTENT (IN)      ::        &
       polgam      ! angle between the north poles of the systems

  REAL (KIND=wp)                   ::        &
       phirot2phi  ! latitude in the geographical system

    ! Local variables

    REAL (KIND=wp)                   ::        &
         zsinpol, zcospol, zphis, zrlas, zarg, zgam

    REAL (KIND=wp), PARAMETER        ::        &
         zrpi18 = r2d,                  &
         zpir18 = d2r
    ! replace by values from mo_physics

    !--------------------------------------------------------------------------

    ! Begin function phirot2phi

    zsinpol     = SIN (zpir18 * polphi)
    zcospol     = COS (zpir18 * polphi)

    zphis       = zpir18 * phirot
    IF (rlarot > 180.0_wp) THEN
       zrlas = rlarot - 360.0_wp
    ELSE
       zrlas = rlarot
    ENDIF
    zrlas       = zpir18 * zrlas

    IF (polgam /= 0.0_wp) THEN
       zgam  = zpir18 * polgam
       zarg  = zsinpol*SIN (zphis) +                                           &
            zcospol*COS(zphis) * ( COS(zrlas)*COS(zgam) - SIN(zgam)*SIN(zrlas) )
    ELSE
       zarg  = zcospol * COS (zphis) * COS (zrlas) + zsinpol * SIN (zphis)
    ENDIF

    phirot2phi  = zrpi18 * ASIN (zarg)

  END FUNCTION phirot2phi
!------------------------------------------------------------------------------

  FUNCTION  phi2phirot ( phi, rla, polphi, pollam )
  !--------------------------------------------------------------------------
  ! This routine converts phi (latitude) from the real geographical system
  ! to phi in the rotated system.
  !--------------------------------------------------------------------------
  REAL (KIND=wp), INTENT (IN)      ::        &
       polphi,  & ! latitude of the rotated north pole
       pollam,  & ! longitude of the rotated north pole
       phi,     & ! latitude in the geographical system
       rla        ! longitude in the geographical system

  REAL (KIND=wp)                   ::        &
       phi2phirot ! longitude in the rotated system

    ! Local variables

    REAL (KIND=wp)                       ::    &
         zsinpol, zcospol, zlampol, zphi, zrla, zarg1, zarg2, zrla1

    REAL (KIND=wp), PARAMETER            ::    &
         zrpi18 = r2d,                  & !
         zpir18 = d2r

    !--------------------------------------------------------------------------

    ! Begin function phi2phirot

    zsinpol  = SIN (zpir18 * polphi)
    zcospol  = COS (zpir18 * polphi)
    zlampol  =      zpir18 * pollam
    zphi     =      zpir18 * phi
    IF (rla > 180.0_wp) THEN
       zrla1  = rla - 360.0_wp
    ELSE
       zrla1  = rla
    ENDIF
    zrla     = zpir18 * zrla1

    zarg1    = SIN (zphi) * zsinpol
    zarg2    = COS (zphi) * zcospol * COS (zrla - zlampol)

    phi2phirot = zrpi18 * ASIN (zarg1 + zarg2)

  END FUNCTION phi2phirot
!------------------------------------------------------------------------------

  FUNCTION  rlarot2rla (phirot, rlarot, polphi, pollam, polgam)
  !----------------------------------------------------------------------------
  ! This function converts lambda (longitude) from one rotated system to
  ! another system. If the optional argument polgam is present, the other
  ! system can also be a rotated one, where polgam is the angle between the two
  ! north poles.
  ! If polgam is not present, the other system is the real geographical
  ! system.
  !----------------------------------------------------------------------------
  REAL (KIND=wp), INTENT (IN)      ::        &
       polphi,   & ! latitude of the rotated north pole
       pollam,   & ! longitude of the rotated north pole
       phirot,   & ! latitude in the rotated system
       rlarot      ! longitude in the rotated system

  REAL (KIND=wp), INTENT (IN)      ::        &
       polgam      ! angle between the north poles of the systems

  REAL (KIND=wp)                   ::        &
       rlarot2rla  ! latitude in the geographical system

    ! Local variables

    REAL (KIND=wp)                   ::        &
         zsinpol, zcospol, zlampol, zphis, zrlas, zarg1, zarg2, zgam

    REAL (KIND=wp), PARAMETER        ::        &
         zrpi18 = r2d,                  & !
         zpir18 = d2r

    !--------------------------------------------------------------------------

    ! Begin function rlarot2rla

    zsinpol = SIN (zpir18 * polphi)
    zcospol = COS (zpir18 * polphi)

    zlampol = zpir18 * pollam
    zphis   = zpir18 * phirot
    IF (rlarot > 180.0_wp) THEN
       zrlas = rlarot - 360.0_wp
    ELSE
       zrlas = rlarot
    ENDIF
    zrlas   = zpir18 * zrlas

    IF (polgam /= 0.0_wp) THEN
       zgam    = zpir18 * polgam
       zarg1   = SIN (zlampol) *                                                &
            (- zsinpol*COS(zphis) * (COS(zrlas)*COS(zgam) - SIN(zrlas)*SIN(zgam))  &
            + zcospol * SIN(zphis))                                               &
            - COS (zlampol)*COS(zphis) * (SIN(zrlas)*COS(zgam) + COS(zrlas)*SIN(zgam))

       zarg2   = COS (zlampol) *                                                &
            (- zsinpol*COS(zphis) * (COS(zrlas)*COS(zgam) - SIN(zrlas)*SIN(zgam))  &
            + zcospol * SIN(zphis))                                               &
            + SIN (zlampol)*COS(zphis) * (SIN(zrlas)*COS(zgam) + COS(zrlas)*SIN(zgam))
    ELSE
       zarg1   = SIN (zlampol) * (-zsinpol * COS(zrlas) * COS(zphis)  +    &
            zcospol *              SIN(zphis)) -    &
            COS (zlampol) *             SIN(zrlas) * COS(zphis)
       zarg2   = COS (zlampol) * (-zsinpol * COS(zrlas) * COS(zphis)  +    &
            zcospol *              SIN(zphis)) +   &
            SIN (zlampol) *             SIN(zrlas) * COS(zphis)
    ENDIF

    IF (zarg2 == 0.0) zarg2 = 1.0E-20_wp

    rlarot2rla = zrpi18 * ATAN2(zarg1,zarg2)

  END FUNCTION rlarot2rla
!------------------------------------------------------------------------------

  FUNCTION  rla2rlarot ( phi, rla, polphi, pollam, polgam )
  !---------------------------------------------------------------------------
  ! This routine converts lambda (longitude) from the real geographical system
  !  to lambda in the rotated system.
  !---------------------------------------------------------------------------
  REAL (KIND=wp), INTENT (IN)      ::        &
       polphi,  & ! latitude of the rotated north pole
       pollam,  & ! longitude of the rotated north pole
       phi,     & ! latitude in geographical system
       rla        ! longitude in geographical system

  REAL (KIND=wp), INTENT (IN)      ::        &
       polgam      ! angle between the north poles of the systems

  REAL (KIND=wp)                   ::        &
       rla2rlarot ! latitude in the the rotated system

    ! Local variables

    REAL (KIND=wp)                       ::    &
         zsinpol, zcospol, zlampol, zphi, zrla, zarg1, zarg2, zrla1

    REAL (KIND=wp), PARAMETER            ::    &
         zrpi18 = r2d,                  & !
         zpir18 = d2r

    !--------------------------------------------------------------------------

    ! Begin function rla2rlarot

    zsinpol  = SIN (zpir18 * polphi)
    zcospol  = COS (zpir18 * polphi)
    zlampol  =      zpir18 * pollam
    zphi     =      zpir18 * phi
    IF (rla > 180.0_wp) THEN
       zrla1  = rla - 360.0_wp
    ELSE
       zrla1  = rla
    ENDIF
    zrla     = zpir18 * zrla1

    zarg1    = - SIN (zrla-zlampol) * COS(zphi)
    zarg2    = - zsinpol * COS(zphi) * COS(zrla-zlampol) + zcospol * SIN(zphi)

    IF (zarg2 == 0.0) zarg2 = 1.0E-20_wp

    rla2rlarot = zrpi18 * ATAN2 (zarg1,zarg2)

    IF (polgam /= 0.0_wp ) THEN
       rla2rlarot = polgam + rla2rlarot
       IF (rla2rlarot > 180._wp) rla2rlarot = rla2rlarot -360._wp
    ENDIF

  END FUNCTION rla2rlarot
!==============================================================================
  function same_horizontal_grid (g1, g2, reordered) result (same)
  !----------------------------------------------------------------------
  ! Determine if two grids are the 'same' (check the relevant meta data).
  ! If the parameter 'reordered' is present two grids are assumed to be
  ! the same even if the grid-points are re-ordered (in order to fit with
  ! the meshes of a coarser grid used for the LETKF calculations.
  !----------------------------------------------------------------------
  type(t_grid) ,intent(in)            :: g1
  type(t_grid) ,intent(in)            :: g2
  logical      ,intent(out) ,optional :: reordered
  logical                             :: same

    same = .false.

    if        (g1% gridtype /= g2% gridtype)  return
    if        (g1% nx       /= g2% nx      )  return
    if        (g1% ny       /= g2% ny      )  return
    if        (g1% nd       /= g2% nd      )  return
    if        (g1% nn       /= g2% nn      )  return
    if        (g1% ni       /= g2% ni      )  return
    if        (g1% ngl      /= g2% ngl     )  return
    if        (g1% di       /= g2% di      )  return
    if        (g1% dj       /= g2% dj      )  return
    if        (g1% scanmode /= g2% scanmode)  return
    if        (g1% lo1      /= g2% lo1     )  return
    if        (g1% la1      /= g2% la1     )  return
    if        (g1% dlonr    /= g2% dlonr   )  return
    if        (g1% dlatr    /= g2% dlatr   )  return
    if (present (reordered)) then
      reordered = any (g1% d_gme /= g2% d_gme)
    else
      if (any (g1% d_gme    /= g2% d_gme   )) return
    endif
    same = .true.

  end function same_horizontal_grid
!==============================================================================
  function same_vertical_grid (g1, g2) result (same)
    !---------------------------------------------------------
    ! Compare vertical discretization parameters of two grids.
    !---------------------------------------------------------
    type(t_grid) ,intent(in) :: g1
    type(t_grid) ,intent(in) :: g2
    logical                  :: same

    same = .false.

    if      (g1% model       /= g2% model      )  return
    if      (g1% nz          /= g2% nz         )  return
    if      (g1% ns          /= g2% ns         )  return
    if      (g1% levtyp      /= g2% levtyp     )  return
    if      (g1% htopf       /= g2% htopf      )  return
    if      (g1% vct         /= g2% vct        )  return
    if      (g1% vc% ivctype /= g2% vc% ivctype)  return
    if      (g1% vc% nlevels /= g2% vc% nlevels)  return
    if (any (g1% vc% vc_uuid /= g2% vc% vc_uuid)) return
    if      (g1% vc% vcflat  /= g2% vc% vcflat )  return
    if (any (g1% ak          /= g2% ak         )) return
    if (any (g1% bk          /= g2% bk         )) return
    if (any (g1% akf         /= g2% akf        )) return
    if (any (g1% bkf         /= g2% bkf        )) return
    same = .true.

  end function same_vertical_grid
!==============================================================================
  SUBROUTINE setup_global_coord (xn, rlon, rlat, lon, lat, lb, ub, kerror)
    !------------------------------------------------------------------
    ! Calculate Cartesian coordinates of the mass points (cell centers)
    ! of the given unstructured grid on the unit sphere.
    !------------------------------------------------------------------
    integer  ,intent(in)  :: lb(4), ub(4)
    real(wp) ,intent(out) :: xn   (lb(1):,lb(2):,lb(3):,lb(4):)
    real(wp) ,intent(out) :: rlon (lb(1):,lb(2):,lb(3):,lb(4):)
    real(wp) ,intent(out) :: rlat (lb(1):,lb(2):,lb(3):,lb(4):)
    real(wp) ,intent(in)  :: lon  (:,:,:)   ! Gridpoint longitudes [radian]
    real(wp) ,intent(in)  :: lat  (:,:,:)   ! Gridpoint latitudes  [radian]
    integer  ,intent(out) :: kerror

    INTEGER  :: j1, j2, jd      ! Loop indices
    integer  :: k
    real(wp) :: sinlat, coslat  ! sin(lat), cos(lat)
    real(wp) :: sinlon, coslon  ! sin(lon), cos(lon)

    kerror = 0

    if (size (rlon) /= size (lon) .or. size (rlat) /= size (lat)) then
       print *, "size (rlon,lon,rlat,lat):", &
            size (rlon), size (lon), size (rlat), size (lat)
       call finish ('setup_global_coord',&
                    'sizes of lon, rlon, lat, and rlat must match')
    end if

    rlon = reshape (lon, shape (rlon))
    rlat = reshape (lat, shape (rlat))

    k = lb(3)
    do jd = lb(4), ub(4)
       do j2 = lb(2), ub(2)
          do j1 = lb(1), ub(1)
             sinlat = sin (rlat(j1,j2,k,jd))
             coslat = cos (rlat(j1,j2,k,jd))
             sinlon = sin (rlon(j1,j2,k,jd))
             coslon = cos (rlon(j1,j2,k,jd))
             xn(j1,j2,1,jd) = coslon * coslat ! x
             xn(j1,j2,2,jd) = sinlon * coslat ! y
             xn(j1,j2,3,jd) =          sinlat ! z
          end do
       end do
    end do

  end SUBROUTINE setup_global_coord
!==============================================================================
  subroutine read_icon_metadata (grid, gridfile, pio)
    type(t_grid)      ,intent(inout) :: grid       ! grid variable
    character(len=*)  ,intent(in)    :: gridfile   ! ICON grid metadata file
    integer, optional ,intent(in)    :: pio        ! pe to use for reading

    if (grid% gridtype == DWD6_ICON) then
       allocate (grid% icongrid)
       call init_icongrid (grid% icongrid, gridfile, pio=pio)
       if (grid% nxny /= grid% icongrid% patch% n_patch_cells_g) then
          if (dace% lpio) then
             write(0,*) "grid size provided:", grid% nxny
             write(0,*) "grid size required:", grid% icongrid% patch% n_patch_cells_g
          end if
          call finish ("read_icon_metadata", "grid sizes do not match")
       end if
       grid% icongrid% refcount = 1
    else
       if (dace% lpio) write(0,*) "gridtype =", grid% gridtype
       call finish ("read_icon_metadata", "gridtype /= ICON")
    end if
  end subroutine read_icon_metadata
!==============================================================================
!--------------------------------------------------------
! broadcast routines for derived types t_refatm, t_vcoord
!--------------------------------------------------------
#define DERIVED type(t_refatm)
#define p_bcast_DERIVED bcast_refatm
#undef  MPI_TYPE
#include "p_bcast.incf"
!------------------------------------------------------------------------------
#undef DERIVED
#undef p_bcast_DERIVED
#define DERIVED type(t_vcoord)
#define p_bcast_DERIVED bcast_vc
#undef  MPI_TYPE
#include "p_bcast.incf"
!------------------------------------------------------------------------------
  subroutine bcast_vcoord (vc, source)
  type(t_vcoord)    ,intent(inout) :: vc
  integer           ,intent(in)    :: source

    if (dace% pe /= source) then
      if (associated (vc% vert_coord)) deallocate (vc% vert_coord)
      if (associated (vc% sigm_coord)) deallocate (vc% sigm_coord)
    endif
    call   bcast_vc  (vc,             source)
    if (dace% pe /= source) then
      nullify (vc% vert_coord)
      nullify (vc% sigm_coord)
    endif
    call p_bcast_ptr (vc% vert_coord, source)
    call p_bcast_ptr (vc% sigm_coord, source)

  end subroutine bcast_vcoord
!==============================================================================

  subroutine reduced_grid (rf, rni, nzr, enkf_grid, grid_coarse)
  !----------------------------------------------------------------------------
  ! Derives reduced (coarser) grid for LETKF weight calculations:
  !   global   GME   or ICON : coarse grid is GME grid.
  !   regional COSMO or ICON : coarse grid is rotated lat-lon grid.
  ! In case of a global grid 'rni' gives the resolution of the coarse GME grid.
  ! In case of a regional grid 'rf' is the coarsening factor
  !----------------------------------------------------------------------------
  integer        ,intent(in)   :: rf         ! horizontal coarsening factor
  integer        ,intent(in)   :: rni        ! reduced ni-resolution
  integer        ,intent(in)   :: nzr        ! number of height levels
  type (t_grid)  ,intent(in)   :: enkf_grid  ! grid of fc. and an.ensemble
  type (t_grid)  ,pointer      :: grid_coarse! reduced grid

    !----------------
    ! local variables
    !----------------
    integer             :: nxr             ! reduced nx
    integer             :: nyr             ! reduced ny
!   integer             :: gridtype        ! gridtype of enkf_grid
    real(wp)            :: delta_x         ! delta for index in x-direction
    real(wp)            :: delta_y         ! delta for index in y-direction
    real(wp) ,parameter :: eps = 1.e-6_wp  ! fix rounding errors
    real(wp)            :: dlat,  dlon     ! S pole of axis of rotation
    real(wp)            :: dlatr, dlonr    ! N pole of axis of rotation
    real(wp)            :: dlat1, dlon1    ! latitude  bounds of rot.grid
    real(wp)            :: dlatn, dlonn    ! longitude bounds of rot.grid
    real(wp)            :: di, dj          ! coarse grid spacing

    !----------------------------------------------------------
    ! consistency checks of horizontal interpolation parameters
    !----------------------------------------------------------
    if (enkf_grid% global) then
      if (rf/=1)                                                       &
        call finish('reduced_grid',                                    &
                    'gridtype and reduced grid parameters inconsistent')
    else
      if (rni/=1)                                                      &
        call finish('reduced_grid',                                    &
                    'gridtype and reduced grid parameters inconsistent')
    endif

    if (rf /= 1 .or. rni/=1) then
      allocate (grid_coarse)
      select case(enkf_grid% gridtype)
      !-----------
      ! COSMO grid
      !-----------
      case(WMO6_LATLON, WMO6_ROTLL)
        !-----------------------------------------
        ! derive number of grid-points, increments
        !-----------------------------------------
        nxr     = enkf_grid% nx / rf
        nyr     = enkf_grid% ny / rf
        delta_x = (enkf_grid% nx - 1 + 2 * eps) / (nxr - 1)
        delta_y = (enkf_grid% ny - 1 + 2 * eps) / (nyr - 1)
        !---------------------
        !construct coarse grid
        !---------------------
!       gridtype = enkf_grid% gridtype
        call construct (grid_coarse,                                    &
                        template = enkf_grid,                           &
                        nx       = nxr,                                 &
                        ny       = nyr,                                 &
                        ke       = nzr,                                 &
                        di       = delta_x*enkf_grid% di,               &
                        dj       = delta_y*enkf_grid% dj,               &
                        lo1      = enkf_grid% lo1 - eps * enkf_grid% di,&
                        la1      = enkf_grid% la1 - eps * enkf_grid% dj )
        call print (grid_coarse)
        if (dace% lpio) write(6,'(a,i4)') ' coarse grid (COSMO) chosen',grid_coarse% gridtype
      !----------------
      ! GME / ICON grid
      !----------------
      case(DWD6_ICOSAHEDRON, DWD6_ICON, WMO6_GAUSSIAN)
        !--------------
        ! global domain
        !--------------
        if (enkf_grid% global) then
          if (dace% lpio) write(6,*) 'before construct grid: gridtype, rni = ',&
                                      enkf_grid% gridtype, rni
          call construct (grid_coarse,                 &
                          gridtype = DWD6_ICOSAHEDRON, &
                          ni       = rni,              &
                          ke       = nzr,              &
                            nproc1 = nproc1,           &
                            nproc2 = nproc2,           &
                              comm = dace% comm        )
          if (dace% lpio) write(6,'(a)') ' coarse grid (GME) chosen'
          call print (grid_coarse)
        !-------------------------
        ! regional grid (ICON-LAM)
        !-------------------------
        else
          !----------------------------------------------
          ! determine parameters for rotated lat-lon grid:
          ! mid-point of the domain, bounds
          !----------------------------------------------
          call fit_rotll (enkf_grid, dlonr, dlatr, dlon1, dlat1, dlonn, dlatn)
          dlat = - dlatr
          dlon =   dlonr - 180._wp; if (dlon < -180._wp) dlon = dlonr + 180._wp
          !-------------
          ! grid spacing
          !-------------
          di  = enkf_grid% d_deg * rf
          nxr = ceiling ( (dlonn - dlon1) / di) + 1
          nyr = ceiling ( (dlatn - dlat1) / di) + 1
          di  = (dlonn - dlon1) / (nxr - 1)
          dj  = (dlatn - dlat1) / (nyr - 1)
          !----------------------------------
          !construct coarse grid for ICON-LAM
          !----------------------------------
          call construct (grid_coarse,              &
                          gridtype = WMO6_ROTLL,    &
                          nx       = nxr,           &
                          ny       = nyr,           &
                          ke       = nzr,           &
                          di       = di,            &
                          dj       = dj,            &
                          lor      = dlon,          &
                          lar      = dlat,          &
                          lo1      = dlon1,         &
                          la1      = dlat1,         &
                          nproc1   = nproc1,        &
                          nproc2   = nproc2,        &
                          comm     = dace% comm     )
          if (dace% lpio) write(6,'(a)') ' coarse grid (rotated lat-lon) chosen'
          call print (grid_coarse)
        endif
      case default
        call finish('reduced_grid: grid not supported')
      end select
    end if

  end subroutine reduced_grid
!------------------------------------------------------------------------------
  subroutine fit_rotll (grid, dlonr, dlatr, dlon1, dlat1, dlonn, dlatn)
    type(t_grid) ,intent(in)  :: grid           ! Limited area grid
    real(wp)     ,intent(out) :: dlonr, dlatr   ! N pole of axis of rotation
    real(wp)     ,intent(out) :: dlon1, dlonn   ! longitude bounds of rot.grid
    real(wp)     ,intent(out) :: dlat1, dlatn   ! latitude  bounds of rot.grid
    !---------------------------------------------------------------------
    ! determine parameters of rotated lat-lon grid fitting into given grid
    !---------------------------------------------------------------------
    integer             :: i,j,d           ! grid indices
    real(wp)            :: e(3), f(3)      ! unit vector
    real(wp)            :: dlat,  dlon     ! S pole of axis of rotation
    real(wp)            :: dla,   dlo      ! coordinates in rotated grid
    real(wp)            :: sinlat, coslat
    real(wp)            :: coslon, sinlon
    real(wp) ,parameter :: eps = 1.e-6_wp  ! fix rounding errors

    if (grid% global) call finish('fit_rotll', 'grid must not be global')

    !----------------------------------------------
    ! determine parameters for rotated lat-lon grid:
    ! mid-point of the domain
    !----------------------------------------------
    e = 0._wp
    do d = grid% lb(4), grid% ub(4)
    do j = grid% lb(2), grid% ub(2)
    do i = grid% lb(1), grid% ub(1)
      e = e + grid% xnglob (i,j,:,d)
    end do
    end do
    end do
    e = p_sum (e)
    e = e / sqrt (sum (e*e))
    !------------------------------
    ! poles of the axis of rotation
    !------------------------------
    f(1) =  e(3)
    f(2) =  e(2)
    f(3) = -e(1)
    sinlat =       f (3)
    coslat = sqrt (f (1) ** 2 + f (2) ** 2)
    coslon =       f (1) / coslat
    sinlon =       f (2) / coslat
    dlat   = r2d * atan2 (sinlat, coslat)
    dlon   = r2d * atan2 (sinlon, coslon)
    dlatr  = - dlat
    dlonr  = dlon - 180._wp; if (dlonr < -180._wp) dlonr = dlonr + 360._wp
    !-------------------------------------------
    ! bounds of the domain in the rotated system
    !-------------------------------------------
    dlat1 =  huge(1._wp)
    dlon1 =  huge(1._wp)
    dlatn = -huge(1._wp)
    dlonn = -huge(1._wp)
    do d = grid% lb(4), grid% ub(4)
    do j = grid% lb(2), grid% ub(2)
    do i = grid% lb(1), grid% ub(1)
      dla = phi2phirot ( r2d * grid% rlat(i,j,1,d), &
                         r2d * grid% rlon(i,j,1,d), &
                         dlatr, dlonr )
      dlo = rla2rlarot ( r2d * grid% rlat(i,j,1,d), &
                         r2d * grid% rlon(i,j,1,d), &
                         dlatr, dlonr, 0._wp)
      dlat1 = min (dlat1, dla)
      dlon1 = min (dlon1, dlo)
      dlatn = max (dlatn, dla)
      dlonn = max (dlonn, dlo)
    end do
    end do
    end do
    dlat1 = p_min (dlat1) - eps
    dlon1 = p_min (dlon1) - eps
    dlatn = p_max (dlatn) + eps
    dlonn = p_max (dlonn) + eps

  end subroutine fit_rotll
!------------------------------------------------------------------------------
  subroutine construct_local (grid, ref_grid, d_km, ke, verbose)
    type(t_grid)      ,intent(out) :: grid     ! local grid (lat-lon or rot-ll)
    type(t_grid)      ,intent(in)  :: ref_grid ! reference grid (limited area)
    real(wp),optional ,intent(in)  :: d_km     ! horizontal resolution
    integer ,optional ,intent(in)  :: ke       ! number of levels
    logical ,optional ,intent(in)  :: verbose  ! print resulting grid?
    !----------------------------------------------------------------------
    ! Derive local grid with given resolution from given limited-area grid.
    !----------------------------------------------------------------------
    integer             :: nx, ny          ! grid size
    real(wp)            :: di, dj          ! grid spacing
    real(wp)            :: lo1, la1        ! lat/lon start
    real(wp)            :: dlat,  dlon     ! S pole of axis of rotation
    real(wp)            :: dlatr, dlonr    ! N pole of axis of rotation
    real(wp)            :: dlat1, dlon1    ! latitude  bounds of rot.grid
    real(wp)            :: dlatn, dlonn    ! longitude bounds of rot.grid
    real(wp)            :: d               ! resolution (km)
    real(wp)            :: rf              ! horizontal coarsening factor
    logical             :: verb            ! print grid?
    real(wp) ,parameter :: eps = 1.e-6_wp  ! fix rounding errors

    if (ref_grid% global) &
         call finish ('construct_local', 'ref_grid cannot be global')

    d = 0._wp
    if (present (d_km)) d = d_km
    verb = .false.
    if (present (verbose)) verb = verbose

    select case(ref_grid% gridtype)
    !-----------
    ! COSMO grid
    !-----------
    case(WMO6_LATLON, WMO6_ROTLL)
       if (d <= 0._wp) then
          nx  = ref_grid% nx
          ny  = ref_grid% ny
          di  = ref_grid% di
          dj  = ref_grid% dj
          lo1 = ref_grid% lo1
          la1 = ref_grid% la1
       else
          !-----------------------------------------
          ! derive number of grid-points, increments
          !-----------------------------------------
          rf  = d / ref_grid% d_km
          nx  = int (ref_grid% nx / rf)
          ny  = int (ref_grid% ny / rf)
          di  = (ref_grid% nx - 1 + 2 * eps) / (nx - 1) * ref_grid% di
          dj  = (ref_grid% ny - 1 + 2 * eps) / (ny - 1) * ref_grid% dj
          lo1 = ref_grid% lo1 - eps * ref_grid% di
          la1 = ref_grid% la1 - eps * ref_grid% dj
       end if
       !----------------------
       ! construct coarse grid
       !----------------------
       call construct (grid,                &
                       template = ref_grid, &
                       nx       = nx,       &
                       ny       = ny,       &
                       ke       = ke,       &
                       di       = di,       &
                       dj       = dj,       &
                       lo1      = lo1,      &
                       la1      = la1  )
    !--------------
    ! ICON-LAM grid
    !--------------
    case(DWD6_ICON)
       !-------------------------------------
       ! Use equivalent resolution by default
       !-------------------------------------
       if (d <= 0._wp) d = ref_grid% d_km
       !----------------------------------------------
       ! determine parameters for rotated lat-lon grid:
       ! mid-point of the domain, bounds
       !----------------------------------------------
       call fit_rotll (ref_grid, dlonr, dlatr, dlon1, dlat1, dlonn, dlatn)
       dlat = - dlatr
       dlon =   dlonr - 180._wp; if (dlon < -180._wp) dlon = dlonr + 180._wp
       !-------------
       ! grid spacing
       !-------------
       rf = d / ref_grid% d_km
       di = ref_grid% d_deg * rf
       nx = ceiling ( (dlonn - dlon1) / di) + 1
       ny = ceiling ( (dlatn - dlat1) / di) + 1
       di = (dlonn - dlon1) / (nx - 1)
       dj = (dlatn - dlat1) / (ny - 1)
       !--------------------------------------------
       ! construct rotated lat-lon grid for ICON-LAM
       !--------------------------------------------
       call construct (grid,                  &
                       gridtype = WMO6_ROTLL, &
                       nx       = nx,         &
                       ny       = ny,         &
                       ke       = ke,         &
                       di       = di,         &
                       dj       = dj,         &
                       lor      = dlon,       &
                       lar      = dlat,       &
                       lo1      = dlon1,      &
                       la1      = dlat1,      &
                       nproc1   = nproc1,     &
                       nproc2   = nproc2,     &
                       comm     = dace% comm  )
    case default
       write(0,*) "Unsupported gridtype:", ref_grid% gridtype
       call finish('construct_local: grid not supported')
    end select
    if (verb) call print (grid)

  end subroutine construct_local

!==============================================================================
  subroutine cosmo_vcp (vcp, nvcp, grid)
  real(wp)     ,intent(out) :: vcp (:)  ! vertical coordinate parameters
  integer      ,intent(out) :: nvcp     ! number of vcp
  type(t_grid) ,intent(in)  :: grid     ! grid derived type
  !-----------------------------------------------
  ! compute the GDS vertical coordinate parameters
  ! for the COSMO hybrid z-coordinate.
  ! (old coding style)
  !-----------------------------------------------
    integer :: i

    IF (.NOT. l_ke_in_gds) THEN      !old style of coding the vertical coordinate parameters
      select case (grid% vc% ivctype)
      case (1:2)

        nvcp = grid% nz+1 + 4
        if (nvcp>size(vcp)) call finish('cosmo_vcp','size of vcp is too small')
        vcp (1) = grid% ra% p0sl
        vcp (2) = grid% ra% t0sl
        vcp (3) = grid% ra% dt0lp
        vcp (4) = grid% vc% vcflat
        do i = 1, grid% nz + 1
          vcp (i+4) = grid% vc% vert_coord (i)
        end do
        if (grid% ra% irefatm == 2) then ! Write parameters for new reference atmosphere
          nvcp = grid% nz+1 + 4 + 6
          if (nvcp>size(vcp)) call finish('cosmo_vcp','size of vcp is too small')
          vcp (4+grid% nz+1+1) = REAL(grid% vc% ivctype+100, wp)
          vcp (4+grid% nz+1+2) = 0.0_wp
          vcp (4+grid% nz+1+3) = 0.0_wp
          vcp (4+grid% nz+1+4) = 0.0_wp
          vcp (4+grid% nz+1+5) = grid% ra% delta_t
          vcp (4+grid% nz+1+6) = grid% ra% h_scal
        endif
      case default
        call finish('cosmo_vcp','invalid ivctype')
      end select
    ELSE           !new style of coding the vertical coordinate parameters
      select case (grid% vc% ivctype)
      case (1:2)

        nvcp = grid% nz+1 + 6
        if (nvcp>size(vcp)) call finish('cosmo_vcp','size of vcp is too small')
        if (grid% ra% irefatm == 2) then ! Write parameters for new reference atmosphere
          vcp (1) = REAL(grid% vc% ivctype+100, wp)
        else
          vcp (1) = REAL(grid% vc% ivctype, wp)
        endif
        vcp (2) = REAL (grid% nz,  wp)
        vcp (3) = grid% ra% p0sl
        vcp (4) = grid% ra% t0sl
        vcp (5) = grid% ra% dt0lp
        vcp (6) = grid% vc% vcflat
        do i = 1, grid% nz + 1
          vcp (i+6) = grid% vc% vert_coord (i)
        end do
        if (grid% ra% irefatm == 2) then ! Write parameters for new reference atmosphere
          nvcp = grid% nz+1 + 11
          if (nvcp>size(vcp)) call finish('cosmo_vcp','size of vcp is too small')
          vcp (6+grid% nz+1+1) = 0.0_wp
          vcp (6+grid% nz+1+2) = 0.0_wp
          vcp (6+grid% nz+1+3) = 0.0_wp
          vcp (6+grid% nz+1+4) = grid% ra% delta_t
          vcp (6+grid% nz+1+5) = grid% ra% h_scal
        endif
      case default
        call finish('cosmo_vcp','invalid ivctype')
      end select

    ENDIF
  end subroutine cosmo_vcp

!------------------------------------------------------------------------------
  SUBROUTINE get_vertcoord_new (inv, pv, vcoord, refatm)
  USE data_parameters ,ONLY: iintegers ! KIND parameter for standard integer
  !----------------------------------------------------------------------------
  ! Description:
  !   Read the vertical coordinate parameters (p0sl, t0sl, dt0lp, vcflat and
  !   coordinate parameters) from the grid description section of a record.
  !
  ! Method:
  !   The values are de-gribed with the routine REFSTF from the grib library.
  !----------------------------------------------------------------------------
  INTEGER,            INTENT(IN)    :: inv
  REAL    (KIND=wp),  INTENT(IN)    :: pv(inv)!
! TYPE (vcoord_type), INTENT(INOUT) :: vcoord ! original COSMO derived types,
! TYPE (refatm_type), INTENT(INOUT) :: refatm ! .. not used here
  TYPE (t_vcoord),    INTENT(OUT)   :: vcoord ! reference coordinate parameters
  TYPE (t_refatm),    INTENT(OUT)   :: refatm ! vertical coordinate parameters


  ! Local Variables
  INTEGER (KIND=iintegers)    :: k, kec, idummy, ierrstat
  CHARACTER (LEN=25) yroutine
  CHARACTER (LEN=80) yerrmsg

  integer :: my_cart_id = 0
  integer :: ke1

  yroutine = 'get_vertcoord_new'
  idummy = NINT (pv( 1), iintegers)

  IF ((idummy >= 1) .AND. (idummy <= 200)) THEN
    !------------------------------------------------------------------------
    ! This is the new grib GDS coding style introduced with LM 3.18
    !
    ! new reference atmosphere has vertical coordinate types 101-103 in Grib
    ! (which correspond to types 1-3 in old reference atmosphere)
    !------------------------------------------------------------------------
    IF (idummy <= 100) THEN
      vcoord%ivctype = idummy
      refatm%irefatm = 1
    ELSE
      !--------------------------------------------------------------------
      ! if the reference atmosphere is determined, vertical coordinate type
      ! is set to values from 1 - 3
      !--------------------------------------------------------------------
      vcoord%ivctype = idummy - 100
      refatm%irefatm = 2
    ENDIF
!   !-------------------------------------
!   ! Check the number of vertical levels:
!   !-------------------------------------
!   IF (NINT(pv(2), iintegers) /= ke) THEN
!     WRITE (yerrmsg,'(A)')                                       &
!      ' ERROR *** Number of vertical levels of input data is not correct *** '
!     PRINT *,                                                                &
!    ' ERROR *** Number of vertical levels of input data is not correct *** ',&
!      NINT(pv(2)), ke
!     ierrstat = 2008
!     CALL model_abort (my_cart_id, ierrstat, yerrmsg, yroutine)
!   ENDIF

    ke1=NINT(pv(2))+1

    refatm% p0sl    = pv( 3)
    refatm% t0sl    = pv( 4)
    refatm% dt0lp   = pv( 5)
    vcoord% vcflat  = pv( 6)
    vcoord% nlevels = ke1

    ALLOCATE (vcoord% vert_coord(ke1))
    ALLOCATE (vcoord% sigm_coord(ke1))

    IF     (vcoord%ivctype == 1) THEN
      DO k = 1,ke1
        vcoord%sigm_coord(k)  = pv(6+k)
      ENDDO
    ELSEIF ( ANY( vcoord%ivctype == (/2,3,4/) ) ) THEN
      DO k = 1,ke1
        vcoord%vert_coord(k)  = pv(6+k)
      ENDDO
    ENDIF

    IF ( ANY( vcoord%ivctype == (/3,4/) ) ) THEN
      !---------------------------------
      ! read three more SLEVE parameters
      !---------------------------------
      vcoord% svc1   = pv(6 + ke1 + 1)
      vcoord% svc2   = pv(6 + ke1 + 2)
      vcoord% nfltvc = NINT (pv(6 + ke1 + 3), iintegers)
    ELSE
      !-------------------------------
      ! just to have an initialization
      !-------------------------------
      vcoord% svc1   = 0.0_wp
      vcoord% svc2   = 0.0_wp
      vcoord% nfltvc = 0_iintegers
    ENDIF

    IF (refatm%irefatm == 2) THEN
      !--------------------------------------------------------
      ! read additional parameters for new reference atmosphere
      !--------------------------------------------------------
      refatm%delta_t = pv(6 + ke1 + 4)
      refatm%h_scal  = pv(6 + ke1 + 5)
    ENDIF

  ELSE
    !-------------------------------------------------------------------
    ! this is the old style of coding the vertical coordinate parameters
    !-------------------------------------------------------------------
    refatm%p0sl   = pv( 1)
    refatm%t0sl   = pv( 2)
    refatm%dt0lp  = pv( 3)
    vcoord%vcflat = pv( 4)

    !--------------------------------------------------
    ! check, how many vertical levels are in input data
    ! check whether pv is descending or ascending
    !--------------------------------------------------
    IF     (pv(5) < pv(6)) THEN
      !------------------------------------------------------------------------
      ! ascending (pressure based): count number of values until 1.0 is reached
      !------------------------------------------------------------------------
      kec = 2
      checkp: DO k = 6, inv
        IF ( (pv(k) > pv(k-1)) .AND. (pv(k) <= 0.9999_wp) ) THEN
          kec = kec + 1
        ELSE
          EXIT checkp
        ENDIF
      ENDDO checkp
    ELSEIF (pv(5) > pv(6)) THEN
      !-----------------------------------------------------------------------
      ! descending (height based): count number of values until 0.0 is reached
      !-----------------------------------------------------------------------
      kec = 2
      checkh: DO k = 6, inv
        IF ( (pv(k) < pv(k-1)) .AND. (pv(k) >= 0.5_wp) ) THEN
          kec = kec + 1
        ELSE
          EXIT checkh
        ENDIF
      ENDDO checkh
    ENDIF

    ke1 = kec

    IF (kec /= ke1) THEN
!     WRITE (yerrmsg,'(A)')                                       &
!      ' ERROR *** Number of vertical levels of input data is not correct *** '
!     ierrstat = 2008
!     CALL model_abort (my_cart_id, ierrstat, yerrmsg, yroutine)
    ELSE
      vcoord%nlevels = ke1
    ENDIF

    !-------------------------------------------------------------------------
    ! Get the type of the vertical coordinate and the reference atmosphere
    ! As this is not coded properly in the old coding style, we have to use
    ! a rather crude method here: see, how many vertical coordinate parameters
    ! are coded:
    !-------------------------------------------------------------------------
    SELECT CASE (inv - ke1 - 4)
    CASE (0, 1)
      ! case=1 can occur because between INT2LM and COSMO we wrote
      ! ivctype as another meta data (but not between COSMO-Nudging and COSMO)
      refatm%irefatm = 1   ! for irefatm > 1, more coord. parameters are coded
      IF     (pv(5) < pv(6)) THEN
        ! ascending (pressure based):
        vcoord%ivctype = 1
      ELSEIF (pv(5) > pv(6)) THEN
        ! descending (height based):
        vcoord%ivctype = 2
      ENDIF

    CASE (4)
      refatm%irefatm = 1   ! for irefatm > 1, more coordinate parameters are coded
      ! NOTE: here we allow only ivctype=3, since for SLEVE2 we impose coding using
      !       l_ke_in_gds=.TRUE. (both in int2lm and cosmo), so this part of the
      !       code should never be called with ivctype=4
      vcoord%ivctype = 3   ! Sleve coordinate

      ! read three more SLEVE parameters
      vcoord% svc1   = pv(4 + ke1 + 2)
      vcoord% svc2   = pv(4 + ke1 + 3)
      vcoord% nfltvc = NINT (pv(4 + ke1 + 4), iintegers)

    CASE (5)
      refatm%irefatm = 3   ! for irefatm > 1, more coord. parameters are coded
      vcoord%ivctype = NINT (pv(4 + ke1 + 1), iintegers)

    CASE (6)
      refatm%irefatm = 2   ! for irefatm > 1, more coord. parameters are coded
      vcoord%ivctype = NINT (pv(4 + ke1 + 1), iintegers)
      IF (vcoord%ivctype > 100) THEN
        vcoord%ivctype = vcoord%ivctype - 100
      ENDIF

    CASE DEFAULT
      yerrmsg = ' ERROR *** ivctype and irefatm for input data not available ***'
      ierrstat = 2008
      CALL model_abort (my_cart_id, ierrstat, yerrmsg, yroutine)
    END SELECT

    ALLOCATE (vcoord%vert_coord(ke1))
    ALLOCATE (vcoord%sigm_coord(ke1))

    IF     (vcoord%ivctype == 1) THEN
      DO k = 1,ke1
        vcoord%sigm_coord(k)  = pv(4+k)
      ENDDO
    ELSEIF ( ANY( vcoord%ivctype == (/2,3,4/) ) ) THEN
      DO k = 1,ke1
        vcoord%vert_coord(k)  = pv(4+k)
      ENDDO
    ENDIF
    IF ( ANY( vcoord%ivctype == (/3,4/) ) ) THEN
      ! read three more SLEVE parameters
      vcoord% svc1   = pv(4 + ke1 + 2)
      vcoord% svc2   = pv(4 + ke1 + 3)
      vcoord% nfltvc = NINT (pv(4 + ke1 + 4), iintegers)

      ! Check for meaningful values of svc1, svc2 and nfltvc
      IF ((vcoord%svc1 > vcoord%vert_coord(1)) .OR. (vcoord%svc1 < 0.0_wp)) THEN
         yerrmsg = ' ERROR *** svc1 not in allowed range for ivctype = 3/4 ***'
         ierrstat = 2008
         CALL model_abort (my_cart_id, ierrstat, yerrmsg, yroutine)
      ENDIF

      IF ((vcoord%svc2 > vcoord%vert_coord(1)) .OR. (vcoord%svc2 < 0.0_wp)) THEN
        yerrmsg = ' ERROR *** svc2 not in allowed range for ivctype = 3/4 ***'
        ierrstat = 2008
        CALL model_abort (my_cart_id, ierrstat, yerrmsg, yroutine)
      ENDIF

      IF (vcoord%nfltvc <= 0) THEN
        yerrmsg = ' ERROR *** nfltvc must be greater than or equal to zero ***'
        ierrstat = 2008
        CALL model_abort (my_cart_id, ierrstat, yerrmsg, yroutine)
      ENDIF

    ELSE
       ! just to have an initialization
       vcoord% svc1   = 0.0_wp
       vcoord% svc2   = 0.0_wp
       vcoord% nfltvc = 0_iintegers
    ENDIF

    IF (refatm%irefatm == 2) THEN
      ! read additional parameters for new reference atmosphere
      refatm%delta_t = pv(4 + ke1 + 5)
      refatm%h_scal  = pv(4 + ke1 + 6)
    ENDIF

   ENDIF

  END SUBROUTINE get_vertcoord_new

!==============================================================================
END MODULE mo_atm_grid
