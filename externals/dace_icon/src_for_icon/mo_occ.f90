!
!+ GNSS occultation observation raytracing-operator derived type and routines.
!
MODULE mo_occ
!
! Description:
! Definition of data type holding GNSS occultation observations
! and operations thereon.
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
!  Cleanup
! V1_5         2009/05/25 Detlef Pingel
!  read quality flags; distinguish rising and setting occultations (D.Pingel)
! V1_6         2009/06/10 Harald Anlauf
!  fix precision of constants
!  process_occ_3d  : workaround for bug in sxf90 rev.360;
!  read_gpsro_bufr : handle invalid BUFR values in MSEC
! V1_7         2009/08/24 Harald Anlauf
!  check_store_occ : add check for bad level ordering
!                    fix memory leak
!  read_gpsro_bufr : fix double rounding issue
! V1_8         2009/12/09 Andreas Rhodin
!  new: subroutine read_gpsro_netcdf (read files from BUFRX2NetCDF)
! V1_9         2010/04/20 Harald Anlauf
!  Various extensions, optimisations and bug fixes
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_11        2010/06/16 Harald Anlauf
!  check_store_occ: bugfix for corner case in ray selection
!  /GPS_RO/: implement sgm_mono for more flexible monotony check
! V1_13        2011/11/01 Harald Anlauf
!  improved checks on GPSRO fg and obs
! V1_15        2011/12/06 Harald Anlauf
!  check_mono_occ: handle observation status symmetrically w.r.t. obs<->fg
!                  fix off-by-one error for level used in check
! V1_17        2011/12/21 Harald Anlauf
!  Modify workaround for bug in SunStudio 12.0
! V1_19        2012-04-16 Andreas Rhodin
!  option to specify separate input directory for observations
! V1_20        2012-06-18 Harald Anlauf
!  set up full occultation geometry
! V1_22        2013-02-13 Harald Anlauf
!  Add vertical coordinate type for ICON
!  Add workarounds for GRAS data from wave-optics processor
!  provide positions and velocities in Earth frame for GPSRO ray-tracer
! V1_27        2013-11-08 Harald Anlauf
!  horizontal interpolation, remove ak,bk from t_global, account for waterload
! V1_28        2014/02/26 Andreas Rhodin
!  changed interface to new_int
! V1_31        2014-08-21 Harald Anlauf
!  read_occ_ropp: suppress warning for missing ROPP inputs files
! V1_35        2014-11-07 Harald Anlauf
!  set_occ_data: add check for ph
! V1_42        2015-06-08 Harald Anlauf
!  horint_mode; adaptions to MEC, feedback file I/O
! V1_43        2015-08-19 Harald Anlauf
!  enhanced obs.error model: correct polar tropopause parameters,
!  namelist /OBSERR_GNSSRO/
! V1_44        2015-09-30 Harald Anlauf
!  occ_obserr: shape options for blob near cold-point TP
! V1_45        2015-12-15 Andreas Rhodin
!  detect/flag suspicious profiles after FG scan, in LETKF, in MEC
! V1_46        2016-02-05 Andreas Rhodin
!  base decisions on new flag 'vct', not 'ivctype'
! V1_47        2016-06-06 Andreas Rhodin
!  add obstype as selection criterium in namelist /rules/
! V1_51        2017-02-24 Andreas Rhodin
!  use generalised humidity for GPSRO, namelist flag for bug-compatibility
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Andreas Rhodin  MPIfM/DWD  2000-2008  original code
! Maria Tomassini DWD        2004       input format changed to NetCDF
! Oliver Schmid   DWD        2005       new obs data type
! Andreas Rhodin  DWD        2006       read ROPP data format, use 1d-operator
! Andreas Rhodin  DWD        2007       read BUFR data format
! Detlef Pingel   DWD        2006-2008  various changes for exps and operation
!==============================================================================

!-------------------------------------------------
! uncomment to exit in case of unknown BUFR codes:
!
#define CHECKCODES
!-------------------------------------------------

  !=============
  ! modules used
  !=============
  !------------------------
  ! general purpose modules
  !------------------------
  use mo_kind,          only: wp,                 &! working precision
                              sp,                 &! single precision
                              dp,                 &! double precision
                              i8                   ! 64bit integer
  use mo_exception,     only: finish, message      ! abort routine
  use mo_mpi_dace,      only: dace,               &! MPI group info
                              p_bcast,            &! generic broadcast routine
                              p_max,              &! MPI generic max routine
                              p_sum                ! MPI generic sum routine
! use mo_fortran_units, only: get_unit_number,    &! obtain a free unit number
!                             return_unit_number   ! release the unit number
  use mo_namelist,      only: position_nml,       &! position namelist
                              POSITIONED,         &! return value
                              nnml                 ! namelist file unit number
  use mo_constants,     only: omcor                ! (2*pi/day_len)
  use mo_physics,       only: gacc,               &! gravity acceleration
                              Rgas=>R,            &! gas constant
                              d2r,                &! factor degree -> radians
                              r2d,                &! factor radians -> degree
                              rh_q,               &! relative <- specific hum.
                              q_rh,               &! relative -> specific hum.
                              q_rh_adj             ! adjoint routine
  use mo_run_params,    only: obsinput,           &! BUFR file input path
                              path_file            ! concatenate path/filename
  use mo_bufr_dwd,     only: t_bufr,              &! BUFR record data type
                             inv_bufr,            &! invalid value indicator
                             ty_loop_cnt,         &! loop counter  indicator
                             bufr_get_entry_texts,&!
                             bufr_get_entry_units  !
  use mo_p_output,     only: oline,               &! output line buffer
                             iol,                 &! indexof next line to write
                       nl => nextline,            &! increment line number
                             flush_buf_pio,       &! write buffer on I/O PE
                             flush_buf             ! write buffer
  use mo_cntrlvar,     only: gh_rh,               &! generalised humidity from rh
                             rh_gh                 ! relative    humidity from gh
  !---------------------------------------
  ! atmospheric state data type definition
  !---------------------------------------
  use mo_atm_state,     only: t_atm              ! atmospheric state data type
  use mo_atm_grid,      only: t_grid,           &! atmospheric grid  data type
                              VCT_P_ISO,        &! isobaric coordinate
                              VCT_P_HYB,        &! GME/HRM/IFS hybrid vertical coordinate
                              VCT_Z_HYB,        &! hybrid z coordinate           (COSMO)
                              VCT_Z_GEN          ! generalised z coordinate (ICON,COSMO)
  use mo_time,          only: t_time,           &! date&time data type
                              init_time,        &! initialise time data type
                              iyyyymmdd,ihhmmss,&! conversion routines
                              cyyyymmdd,chhmmss,&!
                              i_time             ! t_time->yyyy,mm,dd,hh,mi,ss
  !-----------------
  ! matrix data type
  !-----------------
  use mo_dec_matrix,    only: t_vector_segm      ! segment data type
  !----------------------
  ! observation data type
  !----------------------
  use mo_obs_tables,    only: decr_rpt_use,  &! degrade status of report
                              idb_dbk,       &! index in table rept_stat
                              check_report_0,&! basic checks on reports
                              check_report_1,&! basic checks on reports
                              rept_use        ! status flag table
  use mo_t_use,         only: STAT_DISMISS,  &! status flag: dismiss report
                              STAT_OBS_ONLY, &!            : no model equivalent
                              STAT_REJECTED, &!            : rejected
                              STAT_PASSIVE,  &!            : passive
                              STAT_PAS_REJ,  &!            : passive rejected
                              STAT_ACTIVE,   &!            : active
                              STAT_ACTIVE_0I,&!            : active
                              STAT_ACCEPTED, &!            : accepted
                              CHK_FG,        &!            : first guess check
                              CHK_NOIMPL,    &!            : not implemented
                              CHK_INSDAT,    &!            : insufficient data
                              CHK_BLACKLIST, &!            : invalid sat.id.
                              CHK_NONE,      &!            : no check fired
                              CHK_NOTUSED,   &!            : not used
                              CHK_CONSIST,   &!            : inconsis. profiles
                              CHK_QI,        &!            : quality index
                              CHK_HEIGHT,    &!            : height range
                              CHK_RULE,      &! specific rule
                              t_use,         &! data type to hold state
                              decr_use        ! decrease state of datum
  use mo_t_datum,       only: t_datum,       &! date+time derived type
                              rvind,         &! missing value (real)
                              inv_datum       ! default t_datum variable
  use mo_vqc,           only: svqc,          &! default var.qual.cntrl stdev.
                              vqc_form        ! formulation of obs-costfunction
  use mo_obs_set,       only: t_obs_set,     &! obs data type
                              t_obs_block     ! component of t_obs_set
  use mo_t_obs,         only: t_obs,         &! observation data type
                              t_spot,        &! observation meta data type
                              t_head,        &! component of t_spot
                              t_coord,       &! coordinates
                              new_spot,      &! reserve memory for meta data
                              unitvector,    &! set unit vector on the sphere
                              set_xuv,       &! set unit vector, solar zenith
                              new_obs,       &! reserve memory for observations
                              new_int,       &! reserve memory for interp.space
                              new_par,       &! reserve memory for parameters
                              shrink_report, &! remove passive observations
                              OBS_TV,        &! interp. observation type: Tv
                              OBS_RH,        &! interp. observation type: rh
                              OBS_HS,        &! interp. observation type: geop.
                              ITY_MCOLS,     &! interpolation type: column
                              source,        &! list of Report source files
                              n_source,      &! number of Report source files
                              add_source,    &! add source files to list
                              FT_ROPIC,      &! NetCDF file type id.
                              tsk_name,      &! name of task
  !---------------------
  ! predefined constants
  !---------------------
                              GPSRO,          &! GNSS occultation type id
                              TSK_INIT,       &! task value: initialise module
                              TSK_READ,       &! read observations
                              TSK_SET_CHR,    &! set observation characteristics
                              TSK_SHRINK,     &! release unused obs. in report
                              TSK_SETUP_COLS ,&! setup model columns required
                              TSK_SETUP_FUL0 ,&! setup interpolation space
                              TSK_SETUP_FULL ,&! final setup of PSAS-space
                              TSK_R,          &! setup observational error
                              TSK_Y,          &! run forward operator
                              TSK_K,          &! evaluate linear operator
                              CHR_NONL,       &! characteristics: nonlinear
                              CHR_EXP          !                  expensive
  use mo_fdbk_tables,   only: VN_BENDANG,     &! bending angle observation type
                              OT_GPSRO         ! report type value for GPSRO
  use mo_t_col,         only: t_cols,         &! model columns data type
                              t_col,          &! 1 column of the model
                              COL_T,          &! req. observations: temperature
                              COL_Q,          &! specific humidity
                              COL_X,          &! water + ice load
#ifdef DEBUG_OCC
                              COL_GEOH,       &! geopotential (half levels)
#endif
                              COL_PH           ! pressure     (half levels)
  use mo_wmo_tables,    only: WMO0_EUMET,     &! WMO table 0: center identifier
                              WMO6_LATLON,    &! WMO table 6: gridtypes
                              WMO6_GAUSSIAN,  &!
                              DWD6_ICON,      &!
                              DWD6_ICOSAHEDRON !
  use mo_satid,         only: mnem, satid_mnem !  derive: mnemonic <-> satid
  use mo_obstypes,      only: t_obsid,        &! observation id table entry
                              obstype_dbkz     ! derive obsids from dbkz
  use mo_obs_rules,     only: get_rule,       &! routine to get a rule
                              t_set,          &! result data type
                              iud,            &! undefined integer value
                              rud              ! undefined real value
  !---------------------------------------------------
  ! Interface to Michael Gorbunovs Raytracing operator
  !---------------------------------------------------
  use Earth,            only: geodetic,      &! Geodetic coordinates data type
                              Geod_from_Cart,&! Convert Cartesian to geodetic
                              R_Earth,       &! Average Earth radius [km]
                              GAST            ! Greenwich Apparent Sidereal Time Angle
  use Coordinates,      only: cartesian,     &! Cartesian coordinates data type
                              Rotate          ! Rotation of vector
  use Occ_Coordinates,  only: Perigee         ! Derive perigee point
  use Errors,           only: error_status,  &! status data type definition
                              enter_callee,  &!
              error_callee => error,         &!
                              display_status,&!
                              clear_status    !
  use Defaults,         only: wpr => workpr   !
  use ICO_grid,         only: t_global,      &! Global fields data type
                              gf              !
  use ECHAM_fields,     only: ng1,           &! # grid points for lat/lon interpol.
                              FST_NULL,      &! Field status: not allocated
                              ECHAM_cleanup, &! deallocate module variables
                              ECHAM_init      ! initialize global fields
  use ECHAM_fields_adj, only: ECHAM_init_adj,&!
                              ECHAM_cleanup_adj
  use ECHAM_rays_adj,   only: ECHAM_rays_adj_init,&
                              ECHAM_rays_adj_cleanup,&
                              ECHAM_Refraction,   &!
                              RAStat,             &! Dynamical field status
                              zn_t,               &! d(ZN(i))/d(T(j))
                              zn_q,               &! d(ZN(i))/d(Q(j))
                              zn_p                 ! d(ZN(i))/d(Psur(j))
  use mo_grid_intpol,   only: ECHAM_rays_idx_init,&!
                              echam_1d_idx_init    !
  !-------------------
  ! NetCDF header data
  !-------------------
  use mo_head_netcdf, only:           &!
                 cdfinid => ncid,     &! NetCDF file id
                            s1date,   &! nominal (synoptic) date
                            s1cent,   &! data centre
                            s1cents,  &! data sub centre
                            s1cat,    &! data category
                            s1catls,  &! data sub category
                            stime,    &! header observation time (section1)
                            db_time,  &! data bank time
                            mgg,      &! hour
                            ngg,      &! minute
                            s2ikz,    &! DWD-internal classifier
                            istidn,   &! numeric station number (here: satid)
                            obs_time, &! actual time
                            mlah,     &! latitude
                            mloh       ! longitude
  !---------------------
  ! netCDF f90 interface
  !---------------------
  use netcdf,      only: nf90_open,             &!
                         nf90_close,            &!
                         nf90_Inquire_Dimension,&!
                         nf90_Inquire_Variable, &!
                         nf90_inq_dimid,        &!
                         nf90_inq_varid,        &!
                         nf90_get_var,          &!
                         nf90_get_att,          &!
                         nf90_strerror,         &!
                         NF90_NOWRITE,          &! mode flag to open a dataset
                         NF90_NOERR              ! status return value: no error

  implicit none
!==============================================================================
  !----------------
  ! public entities
  !----------------
  private
  public :: read_occ_nml        ! read namelist /GPS_RO/
  public :: process_occ_3d      ! routine to process GNSS occultation data
  public :: setup_occ_msis      ! Derive pressure level for transfer to MSIS
  public :: operator            ! namelist parameter: 1=1d ,3=3d, 2=hybrid
  public :: t_occ, t_ray        ! GNSS occultation data type
  public :: load_occ, store_occ ! load, store t_occ in t_obs
  public :: read_quality, t_pcd ! read occ quality flag
  public :: invalid             ! invalid value
  public :: occ_col2xi          ! convert t,q ,gp  to  t,rh,gh
  public :: set_occ             ! set occ atmospheric data structures
  public :: destruct_occ        ! deallocate occ atmospheric data structures
  public :: const_occ_data      ! test: set atmosphere to homogeneous state
  public :: restore_occ_data    ! test: restore non-homogeneous atmosphere
  public :: skip_invalid        ! skip invalid rays in subsequent iteration
  public :: skip_lowest         ! skip lowest ray
  public :: chk_1d_strict       ! Strict checks on FG in 1d operator
  public :: rm_obs_only         ! Dismiss rays with status OBS_ONLY
  public :: z_oro_min           ! Min. height above surface, 1d op.
  public :: dz_duct             ! Min. distance above ducting layer
  public :: dxdz_min            ! Lower bound on dx/dz, 1d operator
  public :: verbose             ! 0:silent; >=3 verbose in 1d-operator
  public :: write_profile       ! write N,t,q profile (1=1d, 3=3d)
  public :: institution         ! NetCDF attribute
  public :: model               ! NetCDF attribute
  public :: p_send              ! MPI send    type t_occ
  public :: p_recv              ! MPI receive type t_occ
  public :: read_occ_ropp       ! read occultation data (ROPP format)
  public :: read_gpsro_bufr     ! read GPSRO observation from BUFR record
  public :: read_gpsro_netcdf   ! read GPSRO observation from NetCDF file
  public :: read_fdbk_occ       ! restore t_occ, t_ray after feedback file read
  public :: ztop_f              ! upper bound feedback info
  public :: dz_f                ! increment   feedback info
  public :: check_mono_occ      ! check occ. obs. profile for consistency
  public :: use_waterload       ! Account for condensates (water/ice)
  public :: use_gh              ! use generalised humidity
  public :: jac_bugfix          ! use proper state dependence of Jacobian
  public :: horint_mode         ! Horizontal interpolation mode
  public :: refract_model       ! Refractivity model for forward operator
  public :: check_susp_occ      ! Check for suspicious profiles after FG scan

!==============================================================================
  real(wp), parameter :: invalid     = -999._wp

  integer,  parameter :: nsat        = 16             ! Max. # of satellites
  integer,  parameter :: pcc_invalid = inv_datum% pcc ! Per Cent Confidence

  !-----------------------------------
  ! private GNSS occultation data type
  !-----------------------------------

  type(geodetic) ,parameter :: geo_inv  = &! invalid variable of type geodetic
       geodetic (0._wp, invalid, invalid)

  type(cartesian),parameter :: cart_inv = &! invalid variable of type cartesian
       cartesian (invalid)

  integer, parameter :: PCD_INVALID = 65535 ! 16 bits set

  type t_pcd
     integer             :: PCD          = PCD_INVALID  ! Product quality flags
     Integer             :: quality      = 0
     Integer             :: product      = 0
     Integer             :: occult_type  = 0
     Integer             :: exsphs_proc  = 0
     Integer             :: bangle_procc = 0
     Integer             :: refrac_procc = 0
     Integer             :: meteor_procc = 0
     Integer             :: open_loop    = 0
     Integer             :: reflection   = 0
     Integer             :: l2_signal    = 0
     Integer             :: dummy4       = 0
     Integer             :: dummy5       = 0
     Integer             :: dummy6       = 0
     Integer             :: bckgrnd_prof = 0
     Integer             :: retriev_prof = 0
     Integer             :: missing      = 0
  end type t_pcd


  type t_ray                              ! obs for one bending angle
    real(wp)        :: p       = invalid  ! impact parameter               [km]
    real(wp)        :: eps     = invalid  ! bending angle                 [rad]
    real(wp)        :: var     = invalid  ! bending angle variance      [rad^2]
    real(wp)        :: var2    = invalid  ! variance interpolated-original
    type (geodetic) :: geo     = geo_inv  ! geodetic coordinates          [deg]
    real(wp)        :: bearing = invalid  ! GNSS->LEO line of sight  [deg.true]
    integer         :: ix(4)   = 0        ! indices to imcol
    type(cartesian) :: rleo    = cart_inv ! LEO position
    type(cartesian) :: vleo    = cart_inv ! LEO speed
    type(cartesian) :: rgps    = cart_inv ! GPS position
    type(cartesian) :: vgps    = cart_inv ! GPS speed
    integer         :: pcc     = pcc_invalid ! Per Cent Confidence
  end type t_ray

  type t_occ          !  GNSS occultation data type (ray tracer)
    character(len=64) :: file     = ''      ! Input file name
    type (t_time)     :: time               ! Occultation date
    integer           :: occid    = 0       ! Occultation ID
    character(len=4)  :: leoid    = ''      ! Leo         ID
    integer           :: gnssid   = 0       ! GNSS series ID
    integer           :: prn      = 0       ! GNSS PRN
    real(wp)          :: timeinc  = invalid ! Time increment since start [s]
    type (geodetic)   :: gp       = geo_inv ! lat,lon (georeferencing) [deg]
    type (geodetic)   :: gp1d     = geo_inv ! lat,lon (1d-op.)         [deg]
    type (geodetic)   :: gp3d     = geo_inv ! lat,lon (3d-op.)         [deg]
    type (cartesian)  :: rleo     = cart_inv! LEO position              [km]
    type (cartesian)  :: vleo     = cart_inv! LEO speed               [km/s]
    type (cartesian)  :: rgps     = cart_inv! GPS position              [km]
    type (cartesian)  :: vgps     = cart_inv! GPS speed               [km/s]
    type (cartesian)  :: xlc      = cart_inv! Local curvature center
    real(wp)          :: rlc      = invalid ! Local curvature radius    [km]
    real(wp)          :: bearing  = invalid ! GNSS->LEO line of sight  [deg]
    real(wp)          :: wf       = invalid ! Filter width for RO data  [km]
    real(wp)          :: w        = invalid ! Filter width for BP/CT    [km]
    real(wp)          :: wi       = invalid ! " ionosph. correction     [km]
    integer           :: nray     =  0      ! number of samples
    logical           :: checked  = .false. ! rays are all valid
    logical           :: critical = .false. ! rays in critical domain of BA op.
    type (t_pcd)      :: pcd
  end type t_occ

  !+++++++++++++++++++++++++++++
  ! GNSS RO data type, temporary
  !+++++++++++++++++++++++++++++
  type t_prf_k
    real (wp) :: t
    real (wp) :: q
  end type t_prf_k

  type t_srf_k
    integer           :: i
    integer           :: j
    integer           :: l ! ++ diamond index
    integer           :: im
    real(wp)          :: p
  end type t_srf_k

  type t_occ_k
    integer                 :: n     = 0       ! number of profiles
    real(wp)                :: e     = invalid ! reference bending angle
    type (t_prf_k) ,pointer :: p(:,:)=>NULL()  ! (levels,profiles)
    type (t_srf_k) ,pointer :: s  (:)=>NULL()  !        (profiles)
  end type t_occ_k

!------------------------------------------------------------------------------

  !---------------
  ! NetCDF indices
  !---------------
  integer            :: ncid        ! NetCDF file id
  character(len=128) :: ncfile      ! NetCDF file name
  integer            :: status      ! NetCDF return variable from functions
  integer            :: stanc  (3)  ! NetCDF start  parameter
  integer            :: counc  (3)  ! NetCDF count  parameter
  integer            :: strnc  (3)  ! NetCDF stride parameter

  !=========
  ! namelist
  !=========
  integer  :: operator      = 3       ! 1: 1d operator (Abel transform)
                                      ! 2: multiple 1d operator(at perigee pts)
                                      ! 3: 3d raytracing operator
                                      ! 4: hybrid 1d+3d operator (not implem.)
  integer  :: verbose       = 3       ! even:silent; >=3 verbose in 1d-operator
  logical  :: homogeneous   = .false. ! test: homogeneous atmospheric fields
  logical  :: use_waterload = .false. ! Account for condensates (water/ice)
  logical  :: use_gh        = .true.  ! use generalised humidity (for bug-compatibility)
  logical  :: jac_bugfix    = .true.  ! use proper state dependence of Jacobian
  logical  :: use_gnssid    = .false. ! use gnss series ident as part of statid
  logical  :: effective_geo = .true.  ! store effective(T)/nominal(F) geolocation
  !----------------------------------
  ! input, impact parameter selection
  !----------------------------------
  character(len=64) :: file(20) = ''  ! NetCDF input input file names
  integer  :: sat_ids(nsat) = -1      ! valid satellite ids:
                                      !   3: METOP-B
                                      !   4: METOP-A
                                      !  41: CHAMP
                                      !  42: TerraSAR-X
                                      !  43: TanDEM-X
                                      ! 421: Oceansat-2
                                      ! 722: GRACE A
                                      ! 723: GRACE B
                                      ! 740: COSMIC 1
                                      ! 741: COSMIC 2
                                      ! 742: COSMIC 3
                                      ! 743: COSMIC 4
                                      ! 744: COSMIC 5
                                      ! 745: COSMIC 6
                                      ! 786: C/NOFS
                                      ! 820: SAC C
  real(wp) :: ray_incr (10) = 0._wp   ! increments between rays to pick up
  real(wp) :: ray_levs (10) = 0._wp   ! levels to define increments
  real(wp) :: ray_dist_ref  = 10._wp  ! Max. distance from georef. pt. [degree]
  real(wp) :: occ_start_max = 20._wp  ! Occultation must start below [km]
  logical  :: smooth_eps    = .true.  ! smooth observed  bending angles
  logical  :: keep_lowest   = .false. ! keep lowest ray (but don't smooth it)
  logical  :: skip_lowest   = .false. ! skip lowest ray
  logical  :: skip_invalid  = .true.  ! skip invalid rays in subsequent iter.
  logical  :: chk_1d_strict = .false. ! Strict checks on FG in 1d operator
  logical  :: rm_obs_only   = .false. ! Dismiss rays with status OBS_ONLY
  real(wp) :: h_mean        = 20._wp  ! eval. 1dop. at mean perigee below [km]
  real(wp) :: z_1d          = 0._wp   ! eval. 1dop. at perigee closest to [km]
  real(wp) :: z_oro_min     = 0._wp   ! Min. height above surface, 1d op. [km]
  real(wp) :: dz_duct       =-1._wp   ! Min. distance above ducting layer [km]
  real(wp) :: dxdz_min      = 0._wp   ! Lower bound on dx/dz, 1d operator
  real(wp) :: obs_err_sigma = 0._wp   ! correlation length for correlated R
  integer  :: chk_mono      = 1       ! monotony check: 1=obs, 2=fg, 3=both
  real(wp) :: level_mono    = 0._wp   ! top border for monotony check of bda
  real(sp) :: sgm_mono      = 0._sp   ! obs.err.scale for monotony check of bda
  integer  :: occult_type   = 1       ! 1: setting, 2: rising occs, 3: both
  integer  :: qimin         =-1       ! Minimum quality index (pcc)
  integer  :: horint_mode   = 1       ! Horizontal interpolation mode:
                                      ! 0: nearest model column
                                      ! 1: interpolation by forward operator
  integer  :: refract_model = 1       ! Refractivity model for forward operator
                                      ! 1: Rüeger (2002)
                                      ! 2: Aparicio & Laroche (JGR 116, 2011)
  real(wp) :: refract_scale = 1._wp   ! Refractivity scaling factor
  integer  :: chk_susp      = 0       ! Check for suspicious profiles using FG
                                      ! 0: don't check
                                      ! 1: print only
                                      ! 2: print and reject
  real(sp) :: thresh_susp   = 0.25_sp ! Min. fraction of FG rejected for susp.
  integer  :: max_gaps      = 0       ! Max. number of allowed gaps
  real(wp) :: gap_tol       = 0._wp   ! Threshold for gap detection
  !--------------------
  ! observational error
  !--------------------
  integer  :: modify_obs_err= 1       ! prescribe obs. error:
                                      !   0  : never
                                      !   1,4: for missing obs. error in data
                                      !   2,5: if larger than obs.error in data
                                      !   3,6: always
                                      !   7  : enhanced obs. error model
  real(wp) :: obs_err_0km   = 0.1_wp  ! rel. obs. error at 0 km altitude
  real(wp) :: obs_err_10km  = 0.01_wp ! rel. obs. error at 10 km altitude
  real(wp) :: obs_err_abs   = 6.e-6_wp! minimal obs. error (absolute value)
  real(wp) :: obs_err_height(8)       ! Piecewise linear parameterization
  real(wp) :: obs_err_value (8)       !           of observation error.


  !---------------------------------------------------------
  ! Satellite specific parameters of observation error model
  !---------------------------------------------------------
  type t_occ_obserr
     integer  :: satid(8)     ! LEO  satellite ID
     integer  :: gnssid       ! GNSS series ID
     integer  :: center       ! generating center ID
!    integer  :: subcenter    !         subcenter ID
     real(wp) :: s_strat_t    ! relative error in mid stratosphere, tropics
     real(wp) :: s_strat_p    ! relative error in mid stratosphere, poles
     real(wp) :: h0_tr        ! transition height in tropics [km]
     real(wp) :: h0_pol       ! transition height near poles [km]
     real(wp) :: s0_tr        ! relative error in tropics at 0km
     real(wp) :: s0_pol       ! relative error near poles at 0km
     real(wp) :: s_iono_t     ! absolute ionospheric residual, tropics [rad]
     real(wp) :: s_iono_p     ! absolute ionospheric residual, poles   [rad]
     real(wp) :: h_cpt_t      ! height of tropical cold-point TP "blob" [km]
     real(wp) :: s_cpt_t      ! relative error amplitude of "blob"
     real(wp) :: dh_cpt_t     ! vertical   half width of "blob" [km]
     real(wp) :: dl_cpt_t     ! meridional half width of "blob" [deg]
     real(wp) :: h_cpt_p      ! height of polar cold-point TP "blob" [km]
     real(wp) :: s_cpt_p      ! relative error amplitude of "blob"
     real(wp) :: dh_cpt_p     ! vertical   half width of "blob" [km]
     integer  :: cpt_h_shape  ! latitudinal dependence of CPT "blob" height
  end type t_occ_obserr

  integer                     :: nerr = 0                 ! Number of entries
  type(t_occ_obserr), pointer :: occ_obserr(:) => NULL ()
  !-----------------------------
  ! GPSRO specific feedback file
  !-----------------------------
  integer           :: write_profile = 0      ! write N,t,q profile (1=1d, 3=3d)
  character(len=32) :: institution   = 'DWD'  ! NetCDF attribute
  character(len=32) :: model         = 'GME'  ! NetCDF attribute
  real(wp)          :: ztop_f        = 30._wp ! upper bound feedback info
  real(wp)          :: dz_f          = 0.2_wp ! increment   feedback info

  namelist /GPS_RO/ ray_levs, ray_incr, ray_dist_ref, occ_start_max, &
                    qimin, smooth_eps, keep_lowest, skip_lowest,     &
                    skip_invalid, write_profile, operator,           &
                    file, verbose, h_mean, z_1d, homogeneous,        &
                    modify_obs_err, obs_err_0km, obs_err_10km,       &
                    obs_err_abs, obs_err_height, obs_err_value,      &
                    institution, model, ztop_f, dz_f,                &
                    sat_ids, obs_err_sigma, level_mono, sgm_mono,    &
                    chk_mono, occult_type, chk_1d_strict, z_oro_min, &
                    dz_duct, dxdz_min, use_waterload, horint_mode,   &
                    refract_model, chk_susp, thresh_susp,            &
                    max_gaps, gap_tol, use_gh, jac_bugfix,           &
                    use_gnssid, effective_geo, rm_obs_only,          &
                    refract_scale


  !================================
  ! private constants and variables
  !================================
  integer   ,parameter :: mray = 20000      ! max. no. of rays in I/O
  integer        ,save :: occ_int_size  = 0 ! length of data type t_occ
  integer        ,save :: ray_int_size  = 0 ! length of data type t_ray
  integer        ,save :: occ_byte_size = 0 ! length of data type t_occ
! integer        ,save :: ray_byte_size = 0 ! length of data type t_ray
  type(t_global) ,save :: gf_back           ! backup used by const_occ_data
  real(wp)       ,save :: p0_msis       = 0 ! MSIS base pressure level [Pa]
  real(wp)  ,parameter :: mtokm = 1._wp/1000._wp ! m -> km
  real(wp),       save :: obs_err_hgt(0:9)  ! Helper arrays for
  real(wp),       save :: obs_err_val(0:9)  ! obs.error parameterization
  integer              :: obs_err_nlev  = 0 ! No.levels for obs.err.param.

  !===========
  ! interfaces
  !===========
  interface p_send
    module procedure p_send_occ
    module procedure p_send_occs
  end interface p_send

  interface p_recv
    module procedure p_recv_occ
    module procedure p_recv_occs
  end interface p_recv

!==============================================================================
contains
!==============================================================================
  subroutine read_occ_ropp (obs)
  type (t_obs)     ,intent(inout)         :: obs    ! observations data type
  !------------------------------------------------------------------------
  ! GNSS Radio Occultation data are read from a NetCDF file in ROPP format.
  ! Approximately 30 rays are selected and stored in the observation
  ! data type OBS.
  !------------------------------------------------------------------------
    integer, parameter       ::  mxvob = 360 ! max. no. rays per occultation
    !----------------
    ! local variables
    !----------------
    type(t_spot)             :: spt       ! report meta data variable
    type(t_spot)             :: empty     ! default initialised
    type(t_head)             :: head      ! data usually stored in BUFR header
    type(t_use)              :: use       ! state of the report
    type(t_occ)              :: occ       ! occultation meta data
    type(t_ray) ,pointer     :: rays(:)   ! occultation data (rays) to store
    integer                  :: n         ! number of rays read
    integer                  :: i         ! index variable
    integer                  :: j         ! index variable
    integer                  :: i_file    ! index variable
    logical                  :: lkeep

    !-------------------------------
    ! extra variables for necdf file
    !-------------------------------

    integer         :: noccs       ! occultation events in input file
    integer         :: yyyy,mo,dd  ! occultation date YYYYMMDD
    integer         :: hh,mi,ss,ms ! occultation time HHMMSS
    integer         :: ioc         ! occultation id
    integer         :: iob         ! obs below 35 km
    real(wp)        :: erl(3)      ! local curvature center coord
    real(wp)        :: re          ! local curvature radius km

    real(wp)        :: pn(mxvob)   ! Impact parameter  [m]
    real(wp)        :: en(mxvob)   ! Refraction angles [rad]
    real(wp)        :: vn(mxvob)   ! Refraction angle variance [rad^2]

    real(wp)        :: xl(mxvob,3) ! LEO coordinates in Earth frame
    real(wp)        :: vl(mxvob,3) ! LEO velocities m/s
    real(wp)        :: xg(mxvob,3) ! GPS coordinates in Earth frame
    real(wp)        :: vg(mxvob,3) ! GPS velocities m/s
    real(wp)        :: la_tp(mxvob)! lat deg perigee
    real(wp)        :: lo_tp(mxvob)! lon deg perigee

    real(wp)        :: la          ! lat deg
    real(wp)        :: lo          ! lon degf

    integer         :: buf_type = 240
    character(len=8):: statid
    character(len=5):: leo_id, gns_id

    nullify (rays)
    !--------------------------------------
    ! check for usage of Radio occultations
    !--------------------------------------
    if (rept_use(OT_GPSRO)% use(CHK_NONE) <= STAT_DISMISS) return

    !-----------------------------------
    ! read from files (from one PE only)
    !-----------------------------------
    if (all(file=='')) then
#ifdef ROPP
      if (dace% lpio) &
           call message ('read_occ_ropp','no input file given')
#endif
      return
    endif
    do i_file=1,size(file)
      if (file(i_file) /=' ') then
        call add_source (obsinput, file(i_file),&
                         filetype=FT_ROPIC, obstype=OT_GPSRO)
      endif
    end do
    if (.not.dace% lpio) return
    do i_file=1,size(file)
      if (file(i_file) ==' ') cycle

      write(6,'(a,a)') ' reading ',file(i_file)

      !-----------------
      ! open NetCDF file
      !-----------------
      ncfile = path_file(obsinput, file(i_file))
      status = nf90_open (ncfile, NF90_NOWRITE, ncid)
      if (status /= NF90_NOERR) call error ('nf90_open')

      !----------------------------------
      ! read number of occultation events
      !  and number of levels
      !----------------------------------
      call get_dim     ( noccs , 'dim_unlim')
      print * ,' occultations read: ' , noccs
      noccs = min (noccs, rept_use(OT_GPSRO)% max_proc)
      print * ,' occultations used: ' , noccs
      call get_dim     ( iob   , 'dim_lev1b')
      print * ,' level1b levels   : ' , iob
      if (iob > mxvob) call finish('read_occ_ropp','iob > mxvob')

      !-------------------------
      ! read netcdf file, header
      !-------------------------
      ! loop over occultation events
      do j=1,noccs

      ! read header info of occultation event nocc
        stanc  = j          ! read from occultation event j
        counc  = 1
        strnc  = 1
        call get_var_char   ( leo_id  ,'leo_id')
        call get_var_char   ( gns_id  ,'gns_id')
        call get_var_int    ( yyyy    ,'year')
        call get_var_int    ( mo      ,'month')
        call get_var_int    ( dd      ,'day')
        call get_var_int    ( hh      ,'hour')
        call get_var_int    ( mi      ,'minute')
        call get_var_int    ( ss      ,'second')
        call get_var_int    ( ms      ,'msec')
        call get_var_real   ( re      ,'roc')
        call get_var_real   ( la      ,'lat')
        call get_var_real   ( lo      ,'lon')

        stanc  = (/1,j,0/)
        counc  = (/3,1,0/)
        strnc  = (/1,1,0/)
        call get_var_real_1 (erl      ,'r_coc')

        !--------------------
        ! derive station name
        !--------------------
        ioc = hh * 100 + mi
        write(statid,'(a4,i4.4)') leo_id, ioc

        !-----------------------
        ! store in occ data type
        !-----------------------
        occ%      file   = file(i_file)
        occ%      nray   = iob
        occ%      occid  = ioc
        occ%      leoid  = leo_id
!       occ%      gnsid  = gns_id
        occ% gp%  phi    = la
        occ% gp%  lambda = lo
        occ% xlc% x      = erl * mtokm
        occ%      rlc    = re  * mtokm
        occ% gp1d        = occ% gp
        call init_time (occ% time, yyyy, mo, dd, hh, mi, ss)

        !---------------------------------------------------
        ! store header data, insert DWD dbkz, perform checks
        !---------------------------------------------------
        head% modtype  = GPSRO
        head% obstype  = OT_GPSRO
        head% buf_type = buf_type
        head% source   = n_source  ! source file number
        head% record   = j         ! record in file
        head% time     = occ% time
        head% dbkz     = 1694
        call check_report_0 (use, head, 1)
        if (use% state <= STAT_DISMISS) cycle

        !------------------------------
        ! read level-1b data (profiles)
        !------------------------------
        n = iob
        if(associated(rays)) deallocate(rays)
        allocate(rays(n))
        stanc  = (/1,j,0/)    ! read from occultation event j
        counc  = (/n,1,0/)    ! n data
        strnc  = (/1,1,0/)    !

        call get_var_real_1 ( pn,    'impact' )
        call get_var_real_1 ( en,    'bangle' )
        call get_var_real_1 ( vn,    'bangle_sigma' )
        call get_var_real_1 ( la_tp, 'lat_tp')
        call get_var_real_1 ( lo_tp, 'lon_tp')

        stanc  = (/1,1,j/)    ! read from occultation event j
        counc  = (/n,3,1/)    ! n data
        strnc  = (/1,1,1/)    !

        call get_var_real_2 (  xl(1:n,:), 'r_leo_1b' )
        call get_var_real_2 (  vl(1:n,:), 'v_leo_1b' )
        call get_var_real_2 (  xg(1:n,:), 'r_gns_1b' )
        call get_var_real_2 (  vg(1:n,:), 'v_leo_1b' )

        rays%         p    = pn   (1:n) * mtokm
        rays%         eps  = en   (1:n)
        rays%         var  = vn   (1:n) ** 2
        rays% geo% Lambda  =      lo_tp(1:n)
        rays% geo% Phi     =      la_tp(1:n)
        do i=1,3
          rays% rleo% x(i) = xl (1:n,i) * mtokm
          rays% vleo% x(i) = vl (1:n,i)
          rays% rgps% x(i) = xg (1:n,i) * mtokm
          rays% vgps% x(i) = vg (1:n,i)
        end do

        !-------------------------------------------
        ! Final preparation of observation type data
        ! Storage into components of 'obs'
        !-------------------------------------------
        spt = empty
        spt% use                 = use
        spt% hd                  = head
        spt% hd% satid           = satid_mnem (leo_id)
        spt% ident               = spt% hd% satid
        spt% statid              = statid

        call check_store_occ (occ, rays, spt, obs, lkeep)
        call flush_buf_pio
        deallocate (rays)
      end do ! read occultation events

      !------------------
      ! close NetCDF file
      !------------------
      status = nf90_close (ncid)
      if (status /= NF90_NOERR) call error ('nf90_close')
    enddo

  end subroutine read_occ_ropp
!============================================================================
! Auxiliary routines to read NetCDF files follow
!----------------------------------------------------------------------------
  subroutine get_dim (dim ,name, stat)
  !-------------------------------
  ! get dimension from Netcdf file
  !-------------------------------
  integer       ,intent(out) :: dim
  character(len=*)  ,intent(in)  :: name
  integer, optional ,intent(out) :: stat
  integer :: dimid
  status = nf90_inq_dimid   (ncid, name, dimid)
  if (present (stat)) stat = status
  if (status /= NF90_NOERR) then
    if (present (stat)) return
    call error ('nf90_inq_dimid ('//name//')')
  endif
  status = nf90_Inquire_Dimension (ncid, dimid, len=dim)
  if (status /= NF90_NOERR) call error ('nf90_inq_dimlen ('//name//')')
  end subroutine get_dim
!----------------------------------------------------------------------------
  subroutine get_var_int (int ,name)
  !---------------------
  ! get 0-D int variable
  !---------------------
  integer      ,intent(out) :: int
  character(len=*) ,intent(in)  :: name
  integer :: varid
  status = nf90_inq_varid (ncid, name, varid)
  if (status /= NF90_NOERR) call error ('nf90_inq_varid ('//name//')')
  status = nf90_get_var (ncid, varid, int, stanc)
  if (status /= NF90_NOERR) call error ('nf90_get_var::int ('//name//')')
  end subroutine get_var_int
!----------------------------------------------------------------------------
  subroutine get_var_char (c ,name)
  !---------------------
  ! get character string
  !---------------------
  character(len=*) ,intent(out) :: c
  character(len=*) ,intent(in)  :: name
  integer :: varid
  status = nf90_inq_varid (ncid, name, varid)
  if (status /= NF90_NOERR) call error ('nf90_inq_varid ('//name//')')
  status = nf90_get_var  (ncid, varid, c,       &
    (/1,stanc(1)/), (/len(c),counc(1)/), (/1,1/))
  if (status /= NF90_NOERR) call error ('nf90_get_var::char ('//name//')')
  end subroutine get_var_char
!----------------------------------------------------------------------------
  subroutine get_var_real (x ,name)
  !---------------------
  ! get 0-D real variable
  !---------------------
  real(dp)         ,intent(out) :: x
  character(len=*) ,intent(in)  :: name
  integer :: varid
  status = nf90_inq_varid (ncid, name, varid)
  if (status /= NF90_NOERR) call error ('nf90_inq_varid ('//name//')')
  status = nf90_get_var (ncid, varid, x, stanc)
  if (status /= NF90_NOERR) call error ('nf90_get_var::double ('//name//')')
  end subroutine get_var_real
!----------------------------------------------------------------------------
  subroutine get_var_real_1 (x ,name)
  !----------------------
  ! get 1-D real variable
  !----------------------
  real(dp)         ,intent(out) :: x (:)
  character(len=*) ,intent(in)  :: name
  integer  :: varid
  status = nf90_inq_varid (ncid, name, varid)    ! get varid
  if (status /= NF90_NOERR) call error ('nf90_inq_varid ('//name//')')
  status = nf90_get_var (ncid, varid, x, stanc, counc, strnc)
  if (status /= NF90_NOERR) call error ('nf90_get_var::double ('//name//')')
  end subroutine get_var_real_1
!----------------------------------------------------------------------------
  subroutine get_var_real_2 (x ,name)
  !----------------------
  ! get 1-D real variable
  !----------------------
  real(dp)         ,intent(out) :: x (:,:)
  character(len=*) ,intent(in)  :: name
  integer  :: varid
  status = nf90_inq_varid (ncid, name, varid)    ! get varid
  if (status /= NF90_NOERR) call error ('nf90_inq_varid ('//name//')')
  status = nf90_get_var (ncid, varid, x, stanc, counc, strnc)
  if (status /= NF90_NOERR) call error ('nf90_get_var::double ('//name//')')
  end subroutine get_var_real_2
!----------------------------------------------------------------------------
  subroutine error (string)
  !-------------------------
  ! abort on error condition
  !-------------------------
  character(len=*) string
  character(len=80) :: str
  str = nf90_strerror (status)
  write (0,*) 'read_occ_ropp: '//trim(string)//', file='//trim(ncfile)
  write (0,*) trim(str)
  write (0,*) 'start  =',stanc
  write (0,*) 'count  =',counc
  write (0,*) 'stride =',strnc
  call finish ('read_occ_ropp',trim(string)//', file='//trim(ncfile))
  end subroutine error
!==============================================================================
  subroutine process_occ_3d (task, spot, obs, atm, cols, xi, y, Jo, Jo_atm, &
                             state)
  integer            ,intent(in)             :: task    ! what to do
  type(t_spot)       ,intent(inout),optional :: spot    ! SPOT observations
  type(t_obs_block)  ,intent(inout),optional :: obs     ! observation data type
  type(t_atm)        ,intent(in)             :: atm     ! atmospheric state
  type(t_cols)       ,intent(in)   ,optional :: cols    ! model columns
  type(t_vector_segm),intent(in)   ,optional :: xi      ! interpolated values
  type(t_vector_segm),intent(inout),optional :: y       ! observed quantity
  real(wp)           ,intent(inout),optional :: Jo      ! obs. cost funct. Jo
  type(t_atm)        ,intent(inout),optional :: Jo_atm  ! gradient:d Jo/d atm
  integer            ,intent(in)   ,optional :: state   ! status flag

  !----------------------------------------------------------------------------
  ! Main observation processing routine
  ! similar for all observation operators
  !----------------------------------------------------------------------------
    !----------------
    ! local variables
    !----------------
    type(t_occ)          :: o           ! occultation data type
    type(t_ray) ,pointer :: rs(:)       ! rays data type
    type(t_ray) ,pointer :: r           ! ray data type pointer
    integer              :: stat        ! return status
    real(wp)             :: em          ! model refraction angle
    real(wp)             :: fe_zn(6)    ! d(e)/d(zn)
    integer              :: k, kd       ! ray index
    integer              :: i, j        ! lon, lat index
    integer              :: id          ! diamond index
    integer              :: l, m, mm    ! horizontal coordinate index
    integer              :: ix          ! column index
    integer              :: n, nn       ! index
    real(wp)             :: Jo_o        ! d Jo / d y
    real(wp)             :: Jk          ! J contribution from level
    integer              :: tsk         ! task, local copy
    integer(i8)          :: iatm        ! model columns required
    integer              :: natm        ! # of model columns required
    logical ,allocatable :: mask(:)     ! determine number of profiles used
    real(wp),allocatable :: Hnew(:,:)
    type(t_grid),pointer :: grid
    real(wp),allocatable :: rh_a(:), q_a(:), t_a(:), p_a(:) ! adjoint variables
    real(wp),allocatable :: rh  (:), q  (:), t  (:), p  (:) !         variables
    real(wp),allocatable :: gh_a(:), gh (:)                 ! gen.humidity
!   real(wp),allocatable :: qx  (:)     ! water load (condensates)
    real(wp)             :: z_a         !, z
    integer              :: ke          ! number of model levels
    logical     ,pointer :: msk(:)      ! true for observations kept
    logical              :: change      ! true for observations changed
    real(wp)             :: corr_obs_ij ! dummy for setup correlated R
    real(wp)             :: d_ij        ! dummy for setup correlated R
    real(wp)             :: bda         ! dummy for setup correlated R
    real(wp)             :: fc          ! dummy for setup correlated R
    real(wp),allocatable :: Rnew(:,:)   ! correlated R (test purpose)
    integer              :: iii         !+++ workaround for NEC sxf90 bug
    type (t_set)         :: set         ! set of characteristics for channels
    logical              :: nneighb     ! nearest neighbour mode
    type(t_occ_k),pointer:: occ_k(:) =>NULL()

    ! TSK_INIT       =     1 ! initialize modules
    ! TSK_READ       =     2 ! read observations
    ! TSK_SET_CHR    =     4 ! set observ. characteristics
    ! TSK_SHRINK     =     8 ! release unused obs. in report
    ! TSK_SETUP_COLS =    16 ! setup columns
    ! TSK_SETUP_FUL0 =    32 ! setup interpolation space
    ! TSK_SETUP_FULL =    64 ! setup description of PSAS-space
    ! TSK_R          =   128 ! setup observational error
    ! TSK_Y          =   256 ! run forward operator
    ! TSK_YH         =   512 ! run linear or forward operator
    ! TSK_H          =  1024 ! run tangent linear operator
    ! TSK_K          =  2048 ! evaluate linear operator

    if (rept_use(OT_GPSRO)% use(CHK_NONE) <= STAT_DISMISS) return
    !==============================
    ! observation non_specific part
    !==============================
    !----------------
    ! tsk == TSK_INIT
    ! read namelist
    !----------------
    grid => atm% grid
    tsk = task
    if (iand (TSK_INIT,tsk) /= 0) then
      call read_occ_nml
      tsk=tsk-TSK_INIT
    endif
    if (tsk==0) return
    !-----------------
    ! tsk == TSK_READ:
    ! read data
    !-----------------
    if (iand (TSK_READ,tsk) /= 0) then
      call read_occ_ropp (obs% o)
      tsk=tsk-TSK_READ
    endif
    if (tsk==0) return
    !--------------------------------
    ! set observation characteristics
    !--------------------------------
    if (iand (TSK_SET_CHR,tsk) /= 0) then
      tsk=tsk-TSK_SET_CHR
    endif
    if (tsk==0) return

    if (verbose >= 2 .and. dace% lpio) print *, "process_occ_3d: ", tsk_name (task)
    !==========================
    ! observation specific part
    !==========================
    if (.not.present(spot)) call finish('process_occ_3d','spot not present')
    if (.not.present(obs )) call finish('process_occ_3d','obs  not present')
    if (spot% hd% modtype /= GPSRO) return
    call load_occ (obs% o, spot, o, rs)

    !--------------------
    ! setup model columns
    !--------------------
    if (iand (TSK_SETUP_COLS,tsk) /= 0) then
      select case (grid% vct)
      case (VCT_P_ISO)
         call finish ("process_occ_3d","isobaric model grid not supported")
      end select

      iatm = COL_T + COL_Q + COL_PH         ! parameters required (T+Q)
      natm = 3                              ! number of parameters required (3)
#ifdef DEBUG_OCC
      iatm = iatm + COL_GEOH
      natm = natm + 1
#endif
      if (use_waterload) then
         iatm = iatm + COL_X                ! condensates (qcl,qci,qr,qs,qg)
         natm = natm + 1
      end if
      !-----------------------------------
      ! mean tangent point for 1d operator
      !-----------------------------------
      nneighb = (horint_mode == 0)
      call ECHAM_1d_idx_init (spot% col,    &!
                              obs% o,       &!
                              iatm,         &! <-- atmospheric parameters
                              natm,         &! <-- no. atmospheric parameters
                              grid,         &! <-- grid description variable
                              spot% i_time, &! <-- time slot
                              nneighb,      &! <-- nearest neighbour mode
                              spot% imcol)   ! --> model columns required
      !-------------------------------------------------------------
      ! determine model columns required by the ray-tracing operator
      !-------------------------------------------------------------
!     kd = max(1, size(rs)/10)
      kd = 1
      k = 1
      do
        call ECHAM_rays_idx_init &
          (rs(k)% rgps,  &! <-- Transmitter position
           rs(k)% rleo,  &! <-- Receiver position
           o% XLC,       &! <-- Curvature center
           grid,         &! <-- grid description variable
           iatm,         &! <-- atmospheric parameters
           natm,         &! <-- no. atmospheric parameters
           spot% i_time, &! <-- time slot
           obs% o,       &! <-> Observation data type
           spot% imcol)   ! --> List of model columns required
        if (k==size(rs)) exit
        k = min (k + kd, size(rs))
      end do
      spot% n_spt = size (spot% imcol)
      spot% mke   = grid% nz
      tsk=tsk-TSK_SETUP_COLS

!!$if (verbose >= 3) then
!!$print *, "ECHAM_rays_idx_init: size(imcol)=", size (spot% imcol)
!!$do k = 1, size (spot% imcol)
!!$print '(a,i3,a,2f9.3)', " imcol(",k," ): lat,lon =", spot%imcol(k)%c%dlat, spot%imcol(k)%c%dlon
!!$!print *, "pe,idx,iatm =", spot% imcol(k)%pe, spot% imcol(k)%idx, spot% imcol(k)%iatm
!!$end do
!!$end if

    endif
    if(tsk==0) goto 999
    !----------------------------------------
    ! final setup of PSAS interpolation space
    !----------------------------------------
    if (iand (TSK_SETUP_FUL0,tsk) /= 0) then
      if (dace% pe==obs% o% pe) then
        m = 2* atm% grid% nz + 1
        call new_int (obs% o, spot, m * size(spot% imcol))
      endif
      tsk=tsk-TSK_SETUP_FUL0
    endif

    !========================================================
    ! set interpolation space: observed quantities and levels
    !========================================================
    if (iand (TSK_SETUP_FULL,tsk) /= 0) then
      if (dace% pe==obs% o% pe) then
        m = 2* atm% grid% nz + 1
        do k=0,size(spot% imcol)-1
          i = spot% imcol(k+1)% imc(1)
          n = k*m
          obs% o% t_int (spot%i%i+n+1:spot%i%i+n+m-1:2) = OBS_TV
          obs% o% t_int (spot%i%i+n+2:spot%i%i+n+m-1:2) = OBS_RH
          obs% o% t_int (             spot%i%i+n+m    ) = OBS_HS
          obs% o% lev   (spot%i%i+n+1:spot%i%i+n+m-1:2) =    (cols% col(i)% p)
          obs% o% lev   (spot%i%i+n+2:spot%i%i+n+m-1:2) =    (cols% col(i)% p)
          obs% o% lev   (             spot%i%i+n+m    ) = log(cols% col(i)% &
                                                                         s% ps)
        end do
      endif
      tsk=tsk-TSK_SETUP_FULL
    endif
    if(tsk==0) goto 999

    !===============================================
    ! tsk == TSK_R
    ! set up R (observation error covariance matrix)
    !===============================================
    if (iand (TSK_R,tsk) /= 0) then
      if (dace% pe==obs% o% pe) then
        !------------------------
        ! set observational error
        !------------------------
         k = obs% R% ia (spot% o% i+1)
         do i=1,spot% o% n
!           obs% R% ia (spot% o% i + i) = k
            iii = spot% o% i + i   !+++ workaround for bug in sxf90 rev.360
            obs% R% ia   (iii) = k
            obs% R% ja     (k) = iii
            obs% R% packed (k) = rs(i)% var
            k = k + 1
         end do
         obs% R% ia (spot% o% i + spot% o% n + 1) = k

         !=================================================
         ! for test purpose:  set up correlated Rnew
         ! (correlation function according to Gaspari-Cohn)
         !=================================================
         if (obs_err_sigma > 0._wp) then
            allocate (Rnew(size(obs% o% olev),size(obs% o% olev)))
            Rnew=0._wp
            do i=1,obs% R% m
               do j=i,obs% R% n
                  fc   = sqrt(10._wp/3._wp)
                  d_ij = abs (obs% o% olev (i) - obs% o% olev (j))
                  bda = d_ij / fc
                  if (d_ij <= fc) then
                     corr_obs_ij =                  &
                          - 1._wp/4._wp*bda**5      &
                          + 1._wp/2._wp*bda**4      &
                          + 5._wp/8._wp*bda**3      &
                          - 5._wp/3._wp*bda**2      &
                          + 1._wp
                  else
                     if (d_ij <= 2*fc) then
                        corr_obs_ij =               &
                               1._wp/12._wp*bda**5  &
                             - 1._wp/2._wp *bda**4  &
                             + 5._wp/8._wp *bda**3  &
                             + 5._wp/3._wp *bda**2  &
                             - 5._wp       *bda     &
                             + 4._wp                &
                             - 2._wp/3._wp *bda**(-1)
                     else
                        corr_obs_ij = 0._wp
                     endif
                  endif
                  Rnew(i,j) = corr_obs_ij * sqrt (rs(i)% var * rs(j)% var)
                  Rnew(j,i) = Rnew(i,j)
               enddo
            enddo
         endif

        !------------------------------------------
        ! setup variational quality control bounds
        !------------------------------------------
        if (.not. associated (obs% o% s_vqc)) then
          allocate (obs% o% s_vqc (obs% o% n_obs))
          obs% o% s_vqc = svqc
        endif
        if (.not. associated (obs% o% f_vqc)) then
          allocate (obs% o% f_vqc (obs% o% n_obs))
          obs% o% f_vqc = vqc_form
        endif
        call get_rule (type     = spot% hd% modtype,  &! <- module      type
                       obstype  = spot% hd% obstype,  &! <- observation type
                       codetype = spot% hd% codetype, &! <- code type
                       bf_type  = iud,                &! <- no BUFR     type
                       bf_subt  = iud,                &! <- no BUFR  subtype
                       db_kz    = iud,                &! <- no Datenbankkennzahl
                       stname   = '',                 &! <- no Station Name
                       lat      = spot% col% c% dlat, &! <- latitude
                       lon      = spot% col% c% dlon, &! <- longitude
                       o        = set            )     ! -> channel information
        if ( set% sgm_vq /= rud) then
          obs% o% s_vqc (spot% o% i+1 : spot% o% i + spot% o% n) = set% sgm_vq
        endif
        if (set% frm_vq /= iud) then
          obs% o% f_vqc (spot% o% i+1 : spot% o% i + spot% o% n) = set% frm_vq
        endif
      endif
      tsk = tsk - TSK_R
    endif
    if(tsk==0) goto 999

    !--------------------------------------------------
    ! setup K matrix: allocate memory for Jakobi matrix
    !--------------------------------------------------
    if (iand (TSK_K+TSK_Y,tsk) /= 0) then
    if (spot% pe_eval == dace% pe)   then
    if (iand (TSK_K      ,tsk) /= 0) then
      !+++++++++++++++++++++++++++++++++++++++++++
      ! temporary: deallocate occ_k, allocate mask
      !+++++++++++++++++++++++++++++++++++++++++++
      if (associated (occ_k)) then
        call finish ("process_occ_3d","occ_k must not be associated")
!       do k=0,ubound (occ_k,1)
!         if (occ_k(k)% n /= 0) then
!           deallocate (occ_k(k)% p)
!           deallocate (occ_k(k)% s)
!         endif
!       end do
!       deallocate (occ_k)
      endif
      allocate (occ_k (0:o% nray))
      allocate (mask(size(spot% imcol)))
      mask = .false.
      !
      ! final
      !
      k = obs% H% ia (spot% i% i + 1)
      do j=1,spot% i% n
        obs% H% ia (spot% i% i +j) = k
      end do
      obs% H% ia (spot% i% i + spot% i% n + 1) = k
!     obs% xi% x (spot% i% i+1:spot% i% i + spot% i% n) = 0._wp
      obs% yi% x (spot% o% i+1:spot% o% i + spot% o% n) = invalid
    endif

    !------------------------------------------------
    ! copy atmospheric state to echam_fields module
    ! Initialization of echam_constituents_adj module
    !------------------------------------------------
!   call set_occ (gf, grid, cols, obs% o, xi)
    call const_occ_data (gf, spot)
    call echam_init_adj
    !-------------------------------
    ! loop over observational levels
    !-------------------------------
    do k = 1, o% nray
      r => rs(k)
      if(r% eps == invalid) then
        !=============
        ! ray not used
        !=============
        if (verbose > 1) write(0,                                        &
          '("process_occ: pe=",i3,", ray=",i3,", z=",f7.3,", skipped")') &
          dace% pe,k,r% p-o% rlc
      else
        !====================
        ! run nonlinear model
        !====================
        !--------------------------
        ! call observation operator
        !--------------------------
        call echam_rays_adj_init
        if (verbose > 1) write(0,                    &
          '("process_occ: pe=",i3,", ray=",i3,", z=",f7.3)') &
          dace% pe,k,r% p-o% rlc
        Call ECHAM_Refraction &
          (r% p,     & ! <-- RO impact parameter [km]
           o% xlc,   & ! <-- Curvature center in Earth frame
           r% rleo,  & ! <-- LEO coordinates in Earth frame
           r% vleo,  & ! <-- LEO velocities in Earth frame
           r% rgps,  & ! <-- GPS coordinates in Earth frame
           r% vgps,  & ! <-- GPS velocities in Earth frame
           EM,       & ! --> Model refraction angle
           FE_ZN,    & ! --> D(E)/D(ZN)
           stat)       ! --> error status
!!$print *, "After ECHAM_Refraction:"
!!$print *, "EM    =", EM
!!$print *, "FE_ZN =", FE_ZN
        if(stat/=0) then
          if (verbose > 1) write(0,                                        &
            '("process_occ: pe=",i3,", ray=",i3,", z=",f7.3,", stat=",i1)')&
            dace% pe,k,r% p-o% rlc,stat
        else
          if (verbose > 1) write(0,&
            '("process_occ: pe=",i3,", ray=",i3,", z=",f7.3,", stat=",i1, &
            &", Em-o=",3f11.8)') dace% pe,k,r% p-o% rlc,stat,em - r % eps,&
            sqrt(r % var2), em
        endif
        if(o% checked) then
          if(stat/=0) then
!           call finish ('process_occ','invalid ray')
            write(0,&
            '("process_occ: pe=",i3,", ray=",i3,", z=",f7.3,", SKIP RAY !!")')&
              dace% pe,k,r% p-o% rlc
            if (skip_invalid) r% eps = invalid
          endif
        else
          if(stat/=0) then
            r% eps = invalid
          else
            if (skip_lowest) r% eps = invalid ! skip lowest ray
            o% checked = .true.
          endif
        endif
        !----------------------------------
        ! setup K matrix: set bending angle
        !----------------------------------
        if (iand (TSK_K,tsk) /= 0) then
          l = 0
          occ_k(k)% e              = invalid  ! temporary
          obs% yi% x (spot% o%i+k) = invalid  ! final
          if(r% eps /= invalid .and. stat==0) then
            do m = 1,size (spot% imcol)
              ix = spot% imcol(m)% imc(1)
              i  = cols% col(ix)% i
              j  = cols% col(ix)% j
              id = cols% col(ix)% l
              if (RAStat(i,j,id) /= FST_NULL) l = l + 1
            end do
            occ_k(k)% e              = em ! temporary
            obs% yi% x (spot% o%i+k) = em ! final
          endif
          occ_k(k)% n = l
          if (l > 0) then
            allocate (occ_k(k)% p (grid% nz,l))
            allocate (occ_k(k)% s          (l))
          endif
          if (verbose > 1) write(0,*) ' number of profiles used:',l
        endif
        !--------------------------
        ! bending angle observation
        !--------------------------
        if (present(y)) y% x (spot% o% i+k) = invalid
        if(r% eps /= invalid .and. stat==0) then
          !--------------------------------------
          ! calculate observational cost function
          !--------------------------------------
          Jk   = 0.5_wp * (em - r % eps)**2 / r % var
          Jo_o =          (em - r % eps)    / r % var
          if (present(Jo)) Jo = Jo + Jk
          !---------------------------
          ! store analysis or forecast
          !---------------------------
          if (present(y)) y% x (spot% o% i+k) = em
          !-------------------
          ! calculate gradient
          !-------------------
          m = 0
          do l=1,size(spot% imcol)
            ix = spot% imcol(l)% imc(1)
            i  = cols% col(ix)% i
            j  = cols% col(ix)% j
            id = cols% col(ix)% l
            if (RAStat(i,j,id) /= FST_NULL) then
              if (present (Jo_atm)) then
                Jo_atm% t (i,j,:,id) = Jo_atm% t (i,j,:,id)  &
                  + Jo_o * matmul ( FE_ZN , zn_t(i,j,id) %m)
                Jo_atm% q (i,j,:,id) = Jo_atm% q (i,j,:,id)  &
                  + Jo_o * matmul ( FE_ZN , zn_q(i,j,id) %m)
                Jo_atm% ps(i,j,1,id) = Jo_atm% ps(i,j,1,id)  &
                  + Jo_o * sum    ( FE_ZN * zn_p(i,j,id) %p)
              endif
              !---------------
              ! setup K matrix
              !---------------
              m = m + 1
              if (iand (TSK_K,tsk) /= 0) then
                occ_k(k)% s  (m)% im = l
                occ_k(k)% s  (m)% i  = i
                occ_k(k)% s  (m)% j  = j
                occ_k(k)% s  (m)% l  = id ! ++ diamond index
                occ_k(k)% s  (m)% p  = sum   (FE_ZN*zn_p(i,j,id) %p)
                occ_k(k)% p(:,m)% t  = matmul(FE_ZN,zn_t(i,j,id) %m)
                occ_k(k)% p(:,m)% q  = matmul(FE_ZN,zn_q(i,j,id) %m)
                mask (l) = .true.
              endif
            endif
          end do
        endif
        call echam_rays_adj_cleanup
      endif
    end do
    !--------
    ! cleanup
    !--------
    call restore_occ_data (gf)
!   call destruct_occ
    call echam_cleanup_adj
    !---------------
    ! setup K matrix
    !---------------
    if (iand (TSK_K,tsk) /= 0) then
      i = count (mask)
      allocate (occ_k(0)% p (grid% nz,i))
      allocate (occ_k(0)% s          (i))
      occ_k(0)% n = i
      occ_k(0)% e = 0._wp
      if (verbose > 1)print *,' total number of profiles used:',i
      ke = cols% ke
      allocate (rh_a(ke), q_a(ke), t_a(ke), p_a(ke), gh_a(ke))
      allocate (rh  (ke), q  (ke), t  (ke), p  (ke), gh  (ke))
!     obs% xi% x (spot% i% i + 1 : spot% i% i + spot% i% n) = 0._wp
      k = obs% H% ia (spot% i% i + 1)
      do j=1,spot% i% n
        obs% H% ia (spot% i% i +j) = k
      end do
      obs% H% ia (spot% i% i + spot% i% n + 1) = k
      m = 0
      do l=1,size(spot% imcol)
        ix = spot% imcol(l)% imc(1)
        i  = cols% col(ix)% i
        j  = cols% col(ix)% j
        id = cols% col(ix)% l
        nn = 2 * ke + 1
        n  = (l-1) * nn

        if (jac_bugfix) then
           !-------------------------------------
           ! use same atm. profile as occ-library
           !-------------------------------------
           q  = gf% q (:,ix)
           t  = gf% t (:,ix)
           p  = gf% p (:,ix)
        else
           q  =     cols% col(ix)% q
           t  =     cols% col(ix)% t
           p  = exp(cols% col(ix)% p)
        end if
!       if (use_waterload) qx = cols% col(ix)% x
!       z  =     cols% col(ix)% s% geosp / gacc
        rh = rh_q (q, t, p)

        if (mask(l)) then
          !----------------------
          ! temporary: set t,q,ps
          !----------------------
          m  = m + 1
          q_a = 1._wp; rh_a=0._wp; t_a=0._wp; p_a=0._wp
          call q_rh_adj (q_a, rh_a, t_a, p_a, rh, t, p)
          !---------------------------------
          ! account for generalised humidity
          !---------------------------------
          if (use_gh) then
            gh   =  gh_rh (rh, .false.)
            call rh_gh (rh, gh, gh_a)
            rh_a = rh_a * gh_a
          endif
          ! Careful! (TODO: see related code in mo_occ_1d, occ-library)
          z_a = cols% col(ix)% s% ps * gacc / (cols% col(ix)% t (ke) * Rgas)
          occ_k(0)% s  (m)% i  = i
          occ_k(0)% s  (m)% j  = j
          occ_k(0)% s  (m)% l  = id
          occ_k(0)% s  (m)% im = l
          occ_k(0)% s  (m)% p  = cols% col(ix)% s% ps
          occ_k(0)% p(:,m)% t  = t
          occ_k(0)% p(:,m)% q  = q
          !---------------
          ! Jakobi matrix:
          !---------------
          allocate (Hnew(spot% o% n, nn       ))
          Hnew=0._wp
          do k=1,o% nray
          do mm=1,occ_k(k) % n
          if (occ_k(k) % s(mm) % im == l) then
            Hnew (k,1:2*ke:2) = occ_k(k) % p(:,mm)% t +    &
                                occ_k(k) % p(:,mm)% q * t_a
            Hnew (k,2:2*ke:2) = occ_k(k) % p(:,mm)% q * rh_a
            Hnew (k,  nn    ) = occ_k(k) % s  (mm)% p * z_a
            exit
          end if
          end do
          end do
          k = obs% H% ia (spot% i% i + n + 1)
          do j=1,nn
            obs% H% ia (spot% i% i + n + j) = k
            do i=1,spot% o% n
              if (Hnew (i,j) /= 0._wp) then
                obs% H% packed (k) = Hnew (i,j)
                obs% H% ja     (k) = spot% o% i + i
                k = k + 1
              endif
            end do
          end do
          obs% H% ia (spot% i% i + n + nn + 1) = k
          deallocate (Hnew)
        else
          k = obs% H% ia (spot% i% i + n + 1)
          do j=1,nn
            obs% H% ia (spot% i% i + n + j) = k
          end do
          obs% H% ia (spot% i% i + n + nn + 1) = k
        endif
!       !------------------
!       ! final: set t,rh,z
!       !------------------
!       obs% xi% x (spot%o%i+1:spot%o%i+2*ke:2) = t    ! temperature
!       obs% xi% x (spot%o%i+2:spot%o%i+2*ke:2) = rh   ! relative humidity
!       obs% xi% x (           spot%o%i+nn    ) = z    ! geopotential
      end do
      deallocate (gh_a, rh_a, q_a, t_a, p_a, gh, rh, q, t, p)
      deallocate (mask)
      !++++++++++++++++++++++++++++
      ! temporary: deallocate occ_k
      !++++++++++++++++++++++++++++
      if (associated (occ_k)) then
        do k=0,ubound (occ_k,1)
          if (occ_k(k)% n /= 0) then
            deallocate (occ_k(k)% p)
            deallocate (occ_k(k)% s)
          endif
        end do
        deallocate (occ_k)
        nullify    (occ_k)
      endif
    endif ! TSK_K
    endif ! pe_eval == dace% pe
    endif ! TSK_K+TSK_Y
    if (iand (TSK_K,tsk) /= 0) tsk = tsk - TSK_K
    if (iand (TSK_Y,tsk) /= 0) tsk = tsk - TSK_Y
    !--------------------
    ! write feedback file
    !--------------------
!   if (iand(tsk,TSK_FEEDBACK)/=0) then
!     if (spot% pe_eval == dace% pe) call write_fdb (spot, obs% o, o, rs)
!   endif

    !------------------------------------------------------------------
    ! TSK_SHRINK:
    ! release unused observations in the report
    !------------------------------------------------------------------
    if (iand (TSK_SHRINK,tsk) /= 0) then
      call shrink_report (spot, obs%o, state, change, mask=msk)
      if (change) then
        call load_occ (obs% o, spot, o, rs)
        j = 0
        do i=1,size(msk)
          if (msk(i)) then
            j      = j + 1
            rs (j) = rs (i)
          endif
        end do
        o%         nray = j
        spot% col% nlev = j
        call store_occ (obs% o, spot, o, rs(1:j))
        deallocate (msk)
      endif
      tsk = tsk - TSK_SHRINK
      if (tsk == 0) return
    endif

    !==========================
    ! abort if any task is left
    !==========================
    if (tsk /= 0) then
      write(0,*)  'process_occ_3d:  unknown task =',tsk
      call finish('process_occ_3d','unknown task')
    endif

    !--------
    ! cleanup
    !--------
999 continue
    call store_occ (obs% o, spot, o, rs)
    deallocate (rs)

  end subroutine process_occ_3d
!==============================================================================
#if 0
  !---------------
  ! write feedback
  !---------------
  subroutine write_fdb (spot, obs, occ, rays)
  type (t_spot) ,intent(in) :: spot    ! general information on observation
  type (t_obs)  ,intent(in) :: obs     ! observation data type
  type (t_occ)  ,intent(in) :: occ     ! meta information on gnss occ.
  type (t_ray)  ,intent(in) :: rays(:) ! information on rays

    integer               :: iu        ! Fortran unit number
    character             :: ns        ! hemisphere indicator ('N' or 'S')
    character             :: ew        ! east/west  indicator ('E' or 'W')
    character             :: co = '#'  ! trailing characters (# for header)
    integer               :: nray      ! number of rays used
    integer               :: k         ! loop index
    real(wp) ,allocatable :: fc  (:)   ! forecasted bending angles
    real(wp) ,allocatable :: ana (:)   ! analysed bending angles
!   integer               :: is, ie    ! indices to access obs% fc, obs% ana
    !----------
    ! open file
    !----------
    iu = get_unit_number ()
    open (iu, file=trim(occ% file)//'.fdb')
    !----------------------------------------
    ! write header information from .inf-file
    !----------------------------------------
    ns   = 'S'; if (occ% gp% phi    >= 0) ns = 'N'
    ew   = 'W'; if (occ% gp% lambda >= 0) ew = 'E'
!   ms   = occ% sec * 1000._wp
    nray = count (rays(:)% eps /= invalid)
    write (iu,'(a)') co
!   write (iu,'(a,i4,2("/",i2.2),t35,a)')&
!     co, occ% yyyy, occ% mo, occ% dd,'! Occultation date (yyyy/mm/dd)'
!   write (iu,'(a,2(i2.2,":"),i2.2,".",i3.3,t35,a)')&
!     co, occ% hh, occ% mi, ms/1000, mod(ms,1000),  &
!                                     '! Occultation UTC (hh:mm:ss.sss)'
    write (iu,'(a,i4.4,t35,a)')        &
      co, spot% ident,                '! Occultation ID'
    write (iu,'(a,f4.1,a1,1x,f5.1,a1,t35,a)')             &
      co, abs(occ% gp% phi), ns, abs(occ% gp% lambda), ew,&
                                      '! Geodetic latitude and longitude [deg]'
    write (iu,'(a,3(f9.3,1x),t35,a)')  &
      co, occ% xlc,                   '! Local curvature center (cartesian)'
    write (iu,'(a,f9.3,t35,a)')        &
      co, occ% rlc,                   '! Local curvature radius [km]'
    write (iu,'(a,f9.3,t35,a)')          &
      co, occ% wf,                    '! GPS data smoothing window [points]'
    write (iu,'(a,f9.3,t35,a)')          &
      co, occ% w,                     '! BP data smoothing window [points]'
    write (iu,'(a,f9.3,t35,a)')          &
      co, occ% wi,                    '! Ionospheric smoothing window [points]'
    !------------------------------------
    ! feedback information on occultation
    !------------------------------------
    write (iu,'(a,t35,a)') co,        '!'
    write (iu,'(a,i9,t35,a)')          &
      co, nray,                       '! Number of rays used'
!   write (iu,'(a,i9,t35,a)')          &
!     co, spot% qcf,                  '! Quality control flag'
    !-----------------------------
    ! feedback information on rays
    !-----------------------------
    write (iu,'(a,t35,a)') co,        '!'
    write (iu,'(a,a)') co,' Bending angles:'
    write (iu,'(a,a)') co,             &
'        z imp.param    observed  forecasted    analysed      fc-obs     ana-obs     obs_err'
    write (iu,'(a)') co
    allocate (fc  (spot% col% nlev))
    allocate (ana (spot% col% nlev))
!   is  = spot% o% i + 1
!   ie  = spot% o% i + spot% o% n
    fc  = invalid;!if(associated(obs% fc )) fc  = obs% fc  (is:ie)
    do k=1, spot% col% nlev
      if (rays(k)% eps == invalid) cycle
      write (iu,'(2f10.3,6f12.8)') rays(k)% p-occ% rlc, rays(k)% p,   &
        rays(k)% eps, fc(k), ana(k),                                  &
        fc(k) - rays(k)% eps, ana(k) - rays(k)% eps, sqrt(rays(k)% var)
    end do
    !---------------------------
    ! close file, release memory
    !---------------------------
    close (iu)
    call return_unit_number (iu)
    deallocate (fc)
    deallocate (ana)
  end subroutine write_fdb
#endif
!==============================================================================
  !---------------------------
  ! Private auxiliary routines
  !---------------------------
!------------------------------------------------------------------------------
  subroutine set_size
  !-------------------------------------------------------
  ! store sizes of derived data types T_OCC and T_RAY
  ! (in termes of size of component OBS% PAR) into
  ! private module variables OCC_INT_SIZE and RAY_INT_SIZE
  !-------------------------------------------------------
    type (t_occ)  :: occ
    type (t_ray)  :: ray
    type (t_obs)  :: obs
    if (occ_int_size == 0) then
      occ_int_size  = size (transfer (occ, obs% par))
      ray_int_size  = size (transfer (ray, obs% par))
      occ_byte_size = size (transfer (occ, (/' '/) ))
!     ray_byte_size = size (transfer (ray, (/' '/) ))
    endif
  end subroutine set_size
!------------------------------------------------------------------------------
  subroutine store_occ (obs, spot, occ, rays)
  type (t_obs)  ,intent(inout) :: obs      ! data of all observations
  type (t_spot) ,intent(inout) :: spot     ! meta data of this observation
  type (t_occ)  ,intent(in)    :: occ      ! occultation meta data
  type (t_ray)  ,intent(in)    :: rays (:) ! occultation data
  !-----------------------------------------------------------------------
  ! Store the data from variables OCC and RAYS in the component PAR of
  ! OBS at position provided by SPOT. Allocate memory for PAR if required.
  !-----------------------------------------------------------------------
    integer ,pointer :: par (:)
    integer          :: n, m
    if (occ_int_size == 0) call set_size
    n = spot%col%nlev
    m = occ_int_size + n * ray_int_size
    if (spot% p% i < 0) call new_par (obs, m, spot=spot)
    if (m < spot% p% n) spot% p% n = m
    if (m > spot% p% n) call finish('store_occ','m > spot% p% n')
    par => obs % par (spot% p% i+1 : spot% p% i + spot% p% n)
    par (1 : occ_int_size) = transfer(occ ,par)
    par (occ_int_size+1 :) = transfer(rays,par)
  end subroutine store_occ
!------------------------------------------------------------------------------
  subroutine load_occ (obs, spot, occ, rays)
  type (t_obs)  ,intent(in)  :: obs      ! data of all observations
  type (t_spot) ,intent(in)  :: spot     ! meta data of this observation
  type (t_occ)  ,intent(out) :: occ      ! occultation meta data
  type (t_ray)  ,pointer     :: rays (:) ! occultation data
  !------------------------------------------------------------------
  ! Load the data from component PAR of OBS from position provided by
  ! SPOT. Store into OCC and RAYS. allocate RAYS with size required.
  !------------------------------------------------------------------

    integer :: i1, in, i

    allocate (rays (spot% col% nlev))
    if (occ_int_size == 0) call set_size
    occ  = transfer (obs% par (spot% p% i+1 : spot% p% i + occ_int_size), occ)

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! work around NEC SX compiler bug
!
!   rays = transfer (obs% par (spot% p% i+1 + occ_int_size : &
!                              spot% p% i   + spot% p% n), rays)
    i1=spot% p% i + occ_int_size
    do i=1,size      (rays)
      in=i1+ray_int_size
      rays(i)=transfer (obs% par (i1+1:in), rays(i))
      i1=in
    end do
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  end subroutine load_occ
!------------------------------------------------------------------------------
  subroutine p_send_occ (occ, dest)
  type (t_occ) ,intent(in) :: occ
  integer      ,intent(in) :: dest
    if (occ_byte_size == 0) call set_size
    call p_send_derivedtype (occ, occ_byte_size, dest, 1, (dace% comm))
  end subroutine p_send_occ
!------------------------------------------------------------------------------
  subroutine p_send_occs (occ, dest)
  type (t_occ) ,intent(in) :: occ(:)
  integer      ,intent(in) :: dest
    if (occ_byte_size == 0) call set_size
    call p_send_derivedtype2 (occ, size(occ)*occ_byte_size, dest,     &
                                                       1, (dace% comm))
  end subroutine p_send_occs
!------------------------------------------------------------------------------
  subroutine p_recv_occ (occ, src)
  type (t_occ) ,intent(out) :: occ
  integer      ,intent(in)  :: src
    if (occ_byte_size == 0) call set_size
    call p_recv_derivedtype (occ, occ_byte_size, src, 1, (dace% comm))
  end subroutine p_recv_occ
!------------------------------------------------------------------------------
  subroutine p_recv_occs (occ, src)
  type (t_occ) ,intent(out) :: occ(:)
  integer      ,intent(in)  :: src
    if (occ_byte_size == 0) call set_size
    call p_recv_derivedtype2 (occ, size(occ)*occ_byte_size, src,     &
                                                      1, (dace% comm))
  end subroutine p_recv_occs
!------------------------------------------------------------------------------
!!$  subroutine check_par (obs, spot)
!!$  type (t_obs)  ,intent(in)  :: obs
!!$  type (t_spot) ,intent(in)  :: spot
!!$  !--------------------------------------------------------
!!$  ! build checksum over segment of OBS% PAR associated with
!!$  ! observation SPOT.
!!$  !--------------------------------------------------------
!!$    integer :: ix, i
!!$    ix = 0
!!$    do i=spot% p% i+1 , spot% p% i + spot% p% n
!!$      ix = ieor (ix, obs% par(i))
!!$    end do
!!$  end subroutine check_par
!==============================================================================

  subroutine read_occ_nml
  !-----------------------
  ! read namelist /GPS_RO/
  !-----------------------
    integer :: ierr
    integer :: i
    !-------------
    ! set defaults
    !-------------
    ray_levs ( : ) = 0._wp
    ray_incr ( : ) = 0._wp
    ray_dist_ref   = 10._wp   ! Max. distance from georef. pt. [degree]
    occ_start_max  = 20._wp   ! Occultation must start below [km]
    smooth_eps     = .true.
    keep_lowest    = .false.  ! keep lowest ray (but don't smooth it)
    skip_lowest    = .false.
    skip_invalid   = .true.
    chk_1d_strict  = .false.  ! Strict checks on FG in 1d operator
    rm_obs_only    = .false.  ! Dismiss rays with status OBS_ONLY
    operator       = 3        ! 1=1d ,3=3d, 2=hybrid
    file           = ''
    verbose        = 3
    sat_ids        = -1       ! valid satellite ids
    h_mean         = 20._wp   ! evaluate 1d-operator at mean perigee below [km]
    z_1d           =  0._wp   ! evaluate 1d-operator at perigee closest to [km]
    z_oro_min      =  0._wp   ! Min. height above surface, 1d op. [km]
    dz_duct        = -1._wp   ! Min. distance above ducting layer [km]
    dxdz_min       =  0._wp   ! Lower bound on dx/dz, 1d operator
    homogeneous    = .false.  ! test: homogeneous atmospheric fields
    use_waterload  = .false.  ! Account for condensates (water/ice)
    use_gh         = .true.   ! use generalised humidity
    jac_bugfix     = .true.   ! use proper state dependence of Jacobian
    use_gnssid     = .false.  ! use gnssid as part of statid
    effective_geo  = .true.   ! store effective(T)/nominal(F) geolocation
    modify_obs_err = 1        ! modify obs. error:
                              ! (1-3: Healy model; 4-6: piecewise linear)
                              !   0  : never
                              !   1,4: for missing obs. error in data
                              !   2,5: if larger than obs.error in data
                              !   3,6: always
                              !   7  : enhanced obs. error model
    obs_err_0km    = 0.1_wp   ! rel. obs. error at 0 km altitude
    obs_err_10km   = 0.01_wp  ! rel. obs. error at 10 km altitude
    obs_err_abs    = 6.e-6_wp ! minimal obs. error (absolute value)
    obs_err_height = HUGE(0._wp)
    obs_err_value  = -1._wp   ! Piecewise linear obs.err.model
    write_profile  = 0        ! write N,t,q profile (1=1d, 3=3d)
    institution    = 'DWD'    ! NetCDF attribute
    model          = 'GME'    ! NetCDF attribute
    ztop_f         = 30._wp   ! upper bound (30 km)
    dz_f           = 0.2_wp   ! increment   (200 m)
    chk_mono       = 1        ! monotony check: 1= obs, 2=fg, 3=both
    level_mono     = 0._wp    ! top border for monotony check of bda [m]
    sgm_mono       = 0._sp    ! obs.err.scale for monotony check of bda
    occult_type    = 1        ! 1: setting, 2: rising occs, 3: both
    obs_err_sigma  = 0._wp    ! correlation length for correlated R
    qimin          =-1        ! Minimum quality index (pcc)
    horint_mode    = 1        ! Horizontal interpolation mode
    refract_model  = 1        ! Refractivity model for forward operator
    refract_scale  = 1._wp    ! Refractivity scaling factor
    chk_susp       = 0        ! Check for suspicious profiles using FG
    thresh_susp    = 0.25_sp  ! Min. fraction of FG rejected for suspicious
    max_gaps       = 0        ! Max. number of allowed gaps
    gap_tol        = 0._wp    ! Threshold for gap detection

    !--------------
    ! read namelist
    !--------------
    if (dace% lpio) then
      call position_nml ('GPS_RO', status=ierr)
      select case (ierr)
      case (POSITIONED)
#if defined(__ibm__)
        read (nnml ,nml=GPS_RO, iostat=ierr)
        if (ierr/=0) call finish ('read_occ_nml','ERROR in namelist /GPS_RO/')
#else
        read (nnml ,nml=GPS_RO)
#endif
      end select
      if (all(ray_levs==0._wp .and. ray_incr==0._wp)) then
        ray_levs (1:4) = (/2.0_wp,5.0_wp,10._wp,20._wp/)
        ray_incr (1:4) = (/0.3_wp,0.5_wp, 1._wp, 0._wp/)
      endif
      write(6,'(a)') repeat('-',79)
      write(6,'(a)')
      write(6,'(a)')    ' Namelist /GPS_RO/:'
      write(6,'(a)')
      write(6,'(a,i3)')       ' operator       =', operator
      write(6,'(a,i3)')       ' horint_mode    =', horint_mode
      write(6,'(a,i3)')       ' occult_type    =', occult_type
      write(6,'(a,i3)')       ' refract_model  =', refract_model
      write(6,'(a,f10.6)')    ' refract_scale  =', refract_scale
      write(6,'(a,i3)')       ' modify_obs_err =', modify_obs_err
      write(6,'(a,i3)')       ' qimin          =', qimin
      write(6,'(a,i3)')       ' verbose        =', verbose
      write(6,'(a,i3)')       ' write_profile  =', write_profile
      write(6,'(a,l3)')       ' use_waterload  =', use_waterload
      write(6,'(a,l3)')       ' use_gh         =', use_gh
      write(6,'(a,l3)')       ' jac_bugfix     =', jac_bugfix
      write(6,'(a,l3)')       ' use_gnssid     =', use_gnssid
      write(6,'(a,l3)')       ' effective_geo  =', effective_geo
      write(6,'(a,l3)')       ' smooth_eps     =', smooth_eps
      write(6,'(a,l3)')       ' homogeneous    =', homogeneous
      write(6,'(a,l3)')       ' keep_lowest    =', keep_lowest
      write(6,'(a,l3)')       ' skip_lowest    =', skip_lowest
      write(6,'(a,l3)')       ' skip_invalid   =', skip_invalid
      write(6,'(a,l3)')       ' chk_1d_strict  =', chk_1d_strict
      write(6,'(a,l3)')       ' rm_obs_only    =', rm_obs_only
      write(6,'(a,es10.2,a)') ' occ_start_max  =', occ_start_max, " km"
      write(6,'(a,es10.2,a)') ' h_mean         =', h_mean, " km"
      write(6,'(a,es10.2,a)') ' z_1d           =', z_1d, " km"
      write(6,'(a,es10.2,a)') ' z_oro_min      =', z_oro_min, " km"
      write(6,'(a,es10.2,a)') ' dz_duct        =', dz_duct, " km"
      write(6,'(a,es10.2)')   ' dxdz_min       =', dxdz_min
      write(6,'(a,i7)')       ' chk_mono       =', chk_mono
      write(6,'(a,es10.2,a)') ' level_mono     =', level_mono, " m"
      write(6,'(a,es10.2)')   ' sgm_mono       =', sgm_mono
      write(6,'(a,i7)')       ' chk_susp       =', chk_susp
      write(6,'(a,f10.2)')    ' thresh_susp    =', thresh_susp
      write(6,'(a,i7)')       ' max_gaps       =', max_gaps
      write(6,'(a,f10.2,a)')  ' gap_tol        =', gap_tol, " km"
      write(6,'(a,es10.2)')   ' obs_err_0km    =', obs_err_0km
      write(6,'(a,es10.2)')   ' obs_err_10km   =', obs_err_10km
      write(6,'(a,es10.2)')   ' obs_err_abs    =', obs_err_abs
      write(6,'(a,9f10.2)')   ' obs_err_height =', pack (obs_err_height, &
                                                   obs_err_height < HUGE(0._wp))
      write(6,'(a,9es10.2)')  ' obs_err_value  =', pack (obs_err_value, &
                                                         obs_err_value >= 0._wp)
      write(6,'(a,es10.2,a)') ' ray_dist_ref   =', ray_dist_ref, " degree"
      write(6,'(a)')
      write(6,'(a)') ' ray levels increment [km]'
      do i=1,size(ray_levs)
        write(6,'(1x,2f10.3)') ray_levs(i), ray_incr(i)
        if (ray_incr(i)==0._wp) exit
      end do
      write(6,'(a)')
    endif
    call p_bcast (ray_levs      ,dace% pio)
    call p_bcast (ray_incr      ,dace% pio)
    call p_bcast (ray_dist_ref  ,dace% pio)
    call p_bcast (smooth_eps    ,dace% pio)
    call p_bcast (keep_lowest   ,dace% pio)
    call p_bcast (skip_lowest   ,dace% pio)
    call p_bcast (skip_invalid  ,dace% pio)
    call p_bcast (chk_1d_strict ,dace% pio)
    call p_bcast (rm_obs_only   ,dace% pio)
    call p_bcast (write_profile ,dace% pio)
    call p_bcast (operator      ,dace% pio)
    call p_bcast (file          ,dace% pio)
    call p_bcast (sat_ids       ,dace% pio)
    call p_bcast (verbose       ,dace% pio)
    call p_bcast (h_mean        ,dace% pio)
    call p_bcast (z_1d          ,dace% pio)
    call p_bcast (z_oro_min     ,dace% pio)
    call p_bcast (dz_duct       ,dace% pio)
    call p_bcast (dxdz_min      ,dace% pio)
    call p_bcast (homogeneous   ,dace% pio)
    call p_bcast (use_waterload ,dace% pio)
    call p_bcast (use_gh        ,dace% pio)
    call p_bcast (jac_bugfix    ,dace% pio)
    call p_bcast (use_gnssid    ,dace% pio)
    call p_bcast (effective_geo ,dace% pio)
    call p_bcast (modify_obs_err,dace% pio)
    call p_bcast (obs_err_0km,   dace% pio)
    call p_bcast (obs_err_10km,  dace% pio)
    call p_bcast (obs_err_abs,   dace% pio)
    call p_bcast (obs_err_height,dace% pio)
    call p_bcast (obs_err_value ,dace% pio)
    call p_bcast (institution,   dace% pio)
    call p_bcast (model,         dace% pio)
    call p_bcast (ztop_f,        dace% pio)
    call p_bcast (dz_f,          dace% pio)
    call p_bcast (level_mono,    dace% pio)
    call p_bcast (sgm_mono,      dace% pio)
    call p_bcast (chk_mono,      dace% pio)
    call p_bcast (chk_susp,      dace% pio)
    call p_bcast (thresh_susp,   dace% pio)
    call p_bcast (max_gaps,      dace% pio)
    call p_bcast (gap_tol,       dace% pio)
    call p_bcast (occult_type,   dace% pio)
    call p_bcast (obs_err_sigma, dace% pio)
    call p_bcast (qimin,         dace% pio)
    call p_bcast (horint_mode,   dace% pio)
    call p_bcast (refract_model, dace% pio)
    call p_bcast (refract_scale, dace% pio)

    select case (modify_obs_err)
    case (1:3,7)
       obs_err_height(1)  = 0
       obs_err_height(2)  = 10
       obs_err_height(3:) = HUGE(0._wp)
       obs_err_value (1)  = obs_err_0km
       obs_err_value (2)  = obs_err_10km
    case (:-1,8:)
       call finish ("read_occ_nml","modify_obs_err out of range")
    end select
    if (modify_obs_err > 0) then
       obs_err_nlev = count (obs_err_height /= HUGE(0._wp))
       if (obs_err_nlev == 0) &
            call finish ("read_occ_nml","obs_err_height not set")
       obs_err_hgt(1:obs_err_nlev) = obs_err_height(1:obs_err_nlev)
       obs_err_hgt(0)              = - HUGE(0._wp) / 2
       obs_err_hgt(obs_err_nlev+1) =   HUGE(0._wp) / 2
       obs_err_val(1:obs_err_nlev) = obs_err_value (1:obs_err_nlev)
       obs_err_val(0)              = obs_err_value (1)
       obs_err_val(obs_err_nlev+1) = obs_err_value (  obs_err_nlev)
    end if
    if (modify_obs_err >= 7) then
       call read_occ_obserr ()
    end if
    if (refract_model < 1 .or. refract_model > 2) &
       call finish ("read_occ_nml","refract_model out of range (1,2)")
    if (abs (refract_scale - 1._wp) > 1.e-2_wp) &
       call finish ("read_occ_nml","refract_scale out of range (0.99,1.01)")
  end subroutine read_occ_nml
!==============================================================================
  subroutine read_occ_obserr ()
    !-------------------------------------------------------
    ! Set parameters for enhanced RO observation error model
    !-------------------------------------------------------
    integer  :: satid(8)     ! LEO  satellite ID
    integer  :: gnssid       ! GNSS series ID
    integer  :: center       ! generating center ID
!   integer  :: subcenter    !         subcenter ID
    real(wp) :: s_strat_t    ! relative error in mid stratosphere, tropics
    real(wp) :: s_strat_p    ! relative error in mid stratosphere, poles
    real(wp) :: h0_tr        ! transition height in tropics [km]
    real(wp) :: h0_pol       ! transition height near poles [km]
    real(wp) :: s0_tr        ! relative error in tropics at 0km
    real(wp) :: s0_pol       ! relative error near poles at 0km
    real(wp) :: s_iono_t     ! absolute ionospheric residual, tropics [rad]
    real(wp) :: s_iono_p     ! absolute ionospheric residual, poles   [rad]
    real(wp) :: h_cpt_t      ! height of tropical cold-point TP "blob" [km]
    real(wp) :: s_cpt_t      ! relative error amplitude of "blob"
    real(wp) :: dh_cpt_t     ! vertical   half width of "blob" [km]
    real(wp) :: dl_cpt_t     ! meridional half width of "blob" [deg]
    real(wp) :: h_cpt_p      ! height of polar cold-point TP "blob" [km]
    real(wp) :: s_cpt_p      ! relative error amplitude of "blob"
    real(wp) :: dh_cpt_p     ! vertical   half width of "blob" [km]
    integer  :: cpt_h_shape  ! latitudinal dependence of CPT blob height (1..3)
    !---------------------------------------------------------------
    namelist /OBSERR_GNSSRO/ satid, gnssid, center,                &
                             s_strat_t, s_strat_p,                 &
                             h0_tr, h0_pol, s0_tr, s0_pol,         &
                             s_iono_t, s_iono_p,                   &
                             h_cpt_t, s_cpt_t, dh_cpt_t, dl_cpt_t, &
                             h_cpt_p, s_cpt_p, dh_cpt_p, cpt_h_shape
    !---------------------------------------------------------------
    integer                     :: ierr, i
    logical                     :: first
    type(t_occ_obserr), pointer :: e, p(:)
#if defined(__ibm__)
    integer                     :: ios
#endif

    if (dace% lpio) then
      write(6,'(a)') repeat('-',79)
      write(6,'()')
      write(6,'(a)') ' Reading namelist /OBSERR_GNSSRO/'
      write(6,'()')
    endif
    nerr  = 0
    first = .true.
    do
      !-------------
      ! set defaults
      !-------------
      satid       = -1          ! Any satellite
      gnssid      = -1          ! Any GNSS series
      center      = -1          ! Any generating center
      s_strat_t   =  0.010_wp   ! rel.error
      s_strat_p   =  6.0e-3_wp  ! rel.error
      h0_tr       = 12.0_wp     ! km
      h0_pol      =  9.0_wp     ! km
      s0_tr       =  0.18_wp    ! rel.error
      s0_pol      =  0.09_wp    ! rel.error
      s_iono_t    =  2.25e-6_wp ! rad
      s_iono_p    =  2.0e-6_wp  ! rad
      s_cpt_t     =  0.020_wp   ! rel.error
      h_cpt_t     = 18.0_wp     ! km
      dh_cpt_t    =  2.5_wp     ! km
      dl_cpt_t    = 15.0_wp     ! deg
      s_cpt_p     =  0.010_wp   ! rel.error
      h_cpt_p     = 10.0_wp     ! km
      dh_cpt_p    =  1.5_wp     ! km
      cpt_h_shape =  1          ! mode (1..3)
      !------------------------------
      ! read namelist /OBSERR_GNSSRO/
      !------------------------------
      if (dace% lpio) then
        call position_nml ('OBSERR_GNSSRO' ,lrewind=first ,status=ierr)
        first = .false.
        select case (ierr)
        case (POSITIONED)
#if defined(__ibm__)
          read (nnml ,nml=OBSERR_GNSSRO, iostat=ios)
          if(ios/=0)call finish('read_occ_obserr',&
                                'ERROR in namelist /OBSERR_GNSSRO/')
#else
          read (nnml ,nml=OBSERR_GNSSRO)
#endif
        end select
      endif
      !-------------------------------------------
      ! exit if no further namelist group is found
      !-------------------------------------------
      call p_bcast (ierr, dace% pio)
      if (ierr /= POSITIONED) exit
      !------------------------
      ! broadcast to other PE's
      !------------------------
      call p_bcast (satid      , dace% pio)
      call p_bcast (gnssid     , dace% pio)
      call p_bcast (center     , dace% pio)
      call p_bcast (s_strat_t  , dace% pio)
      call p_bcast (s_strat_p  , dace% pio)
      call p_bcast (h0_tr      , dace% pio)
      call p_bcast (h0_pol     , dace% pio)
      call p_bcast (s0_tr      , dace% pio)
      call p_bcast (s0_pol     , dace% pio)
      call p_bcast (s_iono_t   , dace% pio)
      call p_bcast (s_iono_p   , dace% pio)
      call p_bcast (s_cpt_t    , dace% pio)
      call p_bcast (h_cpt_t    , dace% pio)
      call p_bcast (dh_cpt_t   , dace% pio)
      call p_bcast (dl_cpt_t   , dace% pio)
      call p_bcast (s_cpt_p    , dace% pio)
      call p_bcast (h_cpt_p    , dace% pio)
      call p_bcast (dh_cpt_p   , dace% pio)
      call p_bcast (cpt_h_shape, dace% pio)

      nerr =  nerr + 1
      p    => occ_obserr
      allocate (occ_obserr(nerr))
      if (nerr > 1) then
         occ_obserr(:nerr-1) = p
         deallocate (p)
      end if

      e => occ_obserr(nerr)
      e% satid       = satid
      e% gnssid      = gnssid
      e% center      = center
      e% s_strat_t   = s_strat_t
      e% s_strat_p   = s_strat_p
      e% h0_tr       = h0_tr
      e% h0_pol      = h0_pol
      e% s0_tr       = s0_tr
      e% s0_pol      = s0_pol
      e% s_iono_t    = s_iono_t
      e% s_iono_p    = s_iono_p
      e% s_cpt_t     = s_cpt_t
      e% h_cpt_t     = h_cpt_t
      e% dh_cpt_t    = dh_cpt_t
      e% dl_cpt_t    = dl_cpt_t
      e% s_cpt_p     = s_cpt_p
      e% h_cpt_p     = h_cpt_p
      e% dh_cpt_p    = dh_cpt_p
      e% cpt_h_shape = cpt_h_shape
    end do

    if (nerr == 0) call finish ("read_occ_obserr","namelist not found")

    if (dace% lpio) then
      write(6,'(a)')   ' GNSS-RO observation errors:'
      write(6,'()')
      do i = 1, nerr
         e => occ_obserr(i)
         write(6,'(a,9i6)')       ' satid          =', pack (e% satid, e% satid > 0)
         write(6,'(a,i6)')        ' gnssid         =', e% gnssid
         write(6,'(a,i6)')        ' center         =', e% center
         write(6,'(a,es10.2)')    ' s_strat_t      =', e% s_strat_t
         write(6,'(a,es10.2)')    ' s_strat_p      =', e% s_strat_p
         write(6,'(a,f6.2,4x,a)') ' h0_tr          =', e% h0_tr,  " km"
         write(6,'(a,f6.2,4x,a)') ' h0_pol         =', e% h0_pol, " km"
         write(6,'(a,es10.2)')    ' s0_tr          =', e% s0_tr
         write(6,'(a,es10.2)')    ' s0_pol         =', e% s0_pol
         write(6,'(a,es10.2,a)')  ' s_iono_t       =', e% s_iono_t, " rad"
         write(6,'(a,es10.2,a)')  ' s_iono_p       =', e% s_iono_p, " rad"
         write(6,'(a,es10.2)')    ' s_cpt_t        =', e% s_cpt_t
         write(6,'(a,f6.2,4x,a)') ' h_cpt_t        =', e% h_cpt_t,  " km"
         write(6,'(a,f6.2,4x,a)') ' dh_cpt_t       =', e% dh_cpt_t, " km"
         write(6,'(a,f6.2,4x,a)') ' dl_cpt_t       =', e% dl_cpt_t, " deg"
         write(6,'(a,es10.2)')    ' s_cpt_p        =', e% s_cpt_p
         write(6,'(a,f6.2,4x,a)') ' h_cpt_p        =', e% h_cpt_p,  " km"
         write(6,'(a,f6.2,4x,a)') ' dh_cpt_p       =', e% dh_cpt_p, " km"
         write(6,'(a,i6)')        ' cpt_h_shape    =', e% cpt_h_shape
         write(6,'(/)')
      end do
      write(6,'(a)') repeat('-',79)
      write(6,'()')
    end if

  end subroutine read_occ_obserr
!==============================================================================
! Routines to handle the data structure used by the occultation
! operators to store atmospheric fields and adjoints:
!
!
! Variables used:
!
!   gf    type(t_global)      data used by the occultation operator
!   grid  type(t_grid)        model grid meta data
!   cols  type(t_cols)        atmospheric data, gathered at observation PEs
!   xi    type(t_vector_segm) atmospheric data, interpolated
!   obs   type(t_obs_block)   observation meta data
!   spot  type(t_spot)        observation meta data, individual report
!
!
! Subroutines defined below:
!
!   occ_col2xi                            convert t,q ,gp  to  t,rh,gh
!   occ_2xi2col                           convert t,rh,gh  to  t,q ,gp
!
!        set_occ      (gf, grid, cols)    set data used by occultation operator
!   destruct_occ      (gf)                deallocate data by the occ. operator
!        set_occ_grid (gf, grid)          set grid meta data
!        set_occ_data (gf, cols)          set atmospheric data: cols -> gf
!        set_occ_data (gf, cols, obs, xi) set atmospheric data: xi   -> gf
!        get_occ_data (gf, obs, xi)       set atmospheric data: xi   <- gf
!      const_occ_data (gf, spot)          patch gf for horiz.homogeneous data
!    restore_occ_data (gf, spot)          restore gf for horiz.inhomog.  data
!==============================================================================
  subroutine occ_col2xi (t, q, geop, p, xi, fg)
  !------------------------------------------------------------------------
  ! converts columns of t,q  and surface geopotential (from different PEs)
  !                  to t,rh and surface geop.height  (interpolation space)
  ! called by subroutine interpolate, module psasutil
  ! called by get_occ_data (currently not used)
  !------------------------------------------------------------------------
    real(wp) ,intent(in)  :: t  (:)  ! temperature         (K)
    real(wp) ,intent(in)  :: q  (:)  ! specific humidity   (kg/kg)
    real(wp) ,intent(in)  :: geop    ! geopotential        (m2/s2)
    real(wp) ,intent(in)  :: p  (:)  ! pressure            (Pa)
    real(wp) ,intent(out) :: xi (:)  ! output: t(K),gh(1),gh(gpm)
    logical  ,intent(in)  :: fg      ! if true: account for rh0fg

      integer :: ke

      ke = size(t)
      xi(1:2*ke  :2) = t              ! temperature          (K)
    if (use_gh) then
      xi(2:2*ke  :2) = gh_rh(        &! generalised humidity (1)
                       rh_q (q, t, p)&! relative humidity
                     , fg, 200._wp, p )
    else
      xi(2:2*ke  :2) = rh_q (q, t, p) ! relative humidity   (1)
    endif
      xi(  2*ke+1  ) = geop / gacc    ! geopotential height (m)

  end subroutine occ_col2xi
!------------------------------------------------------------------------------
  subroutine occ_xi2col (t, q, geop, p, xi)
  !------------------------------------------------------------------------
  ! converts columns of t,rh and surface geop.height  (interpolation space)
  !                  to t,rh and surface geop.height
  ! called by set_occ_data (to set GNSS specific data type)
  !------------------------------------------------------------------------
    real(wpr) ,intent(out) :: t  (:)  ! temperature         (K)
    real(wpr) ,intent(out) :: q  (:)  ! specific humidity   (kg/kg)
    real(wpr) ,intent(out) :: geop    ! geopotential        (m2/s2)
    real(wpr) ,intent(in)  :: p  (:)  ! pressure            (Pa)
    real(wp)  ,intent(in)  :: xi (:)  ! virtual temperature (K)
                                      ! generalised humidity(1)
                                      ! geopotential height (m)
      integer  :: ke
      real(wp) :: rh (size(q))

      ke   = size(t)
    if (use_gh) then
      call rh_gh (rh, xi(2:2*ke  :2))
    else
      rh   =          xi(2:2*ke  :2)
    endif
      t    =          xi(1:2*ke  :2)         ! temperature         (K)
      q    = q_rh    (     rh       ,&       ! specific humidity   (kg/kg)
                      real(t,wp)    ,&
                      real(p,wp)     )
      geop =          xi(  2*ke+1  ) * gacc  ! geopotential        (m2/s2)

  end subroutine occ_xi2col
!------------------------------------------------------------------------------
  subroutine set_occ (gf, grid, cols, obs, xi)
  type (t_global)      ,intent(inout)        :: gf
  type (t_grid)        ,pointer              :: grid
  type (t_cols)        ,intent(in)           :: cols
  type (t_obs)         ,intent(in) ,optional :: obs
  type (t_vector_segm) ,intent(in) ,optional :: xi
  !--------------------------------------------------------------------------
  ! set model data (GNSS specific data type) used by the occultation operator
  !--------------------------------------------------------------------------
    type(error_status) ,pointer :: errstat ! error status variable

    if (rept_use(OT_GPSRO)% use(CHK_NONE) <= STAT_DISMISS) return
    nullify (errstat)
    call enter_callee ('set_occ',errstat)
    call echam_cleanup
    call set_occ_grid (gf, grid)
    call set_occ_data (gf, cols, obs, xi)
    call echam_init (errstat)
    if (error_callee (errstat)) then
      call display_status (errstat)
      call finish ('set_occ','error in echam_init')
    endif
!   call echam_init_adj
    call clear_status (errstat)
    if (associated (errstat)) deallocate (errstat)
  end subroutine set_occ
!------------------------------------------------------------------------------
  subroutine destruct_occ ! (gf)
! type (t_global) ,intent (inout) :: gf
  !-------------------------------------------------
  ! deallocate data used by the occultation operator
  !-------------------------------------------------
    call echam_cleanup
!   call echam_cleanup_adj
  end subroutine destruct_occ
!------------------------------------------------------------------------------
  subroutine set_occ_data (gf, cols, obs, xi)
  type (t_global)      ,intent(inout)        :: gf
  type (t_cols)        ,intent(in)           :: cols
  type (t_obs)         ,intent(in) ,optional :: obs
  type (t_vector_segm) ,intent(in) ,optional :: xi

    integer                :: i, j, k
    integer                :: is
    integer                :: ke21
    type (t_spot) ,pointer :: si
    type (t_col)  ,pointer :: c

    gf% n         = cols% ncol
    gf% yyyymmdd  = iyyyymmdd (cols% time)
    gf% hhmmss    = ihhmmss   (cols% time)

    allocate (gf% s (          gf% n))
    allocate (gf% t (  gf% nz, gf% n))
    allocate (gf% q (  gf% nz, gf% n))
    allocate (gf% qx(  gf% nz, gf% n))
    allocate (gf% p (  gf% nz, gf% n))
    allocate (gf% ph(0:gf% nz, gf% n))
    allocate (gf% i (gf% lbg(1):gf% ubg(1),&
                     gf% lbg(2):gf% ubg(2),&
                     gf% lbg(3):gf% ubg(3)))
    !-------------------------------------------------------------------
    ! Height-based (hybrid) vertical coordinates: heights of half-levels
    !-------------------------------------------------------------------
    if (gf% vctype /= VCT_P_HYB) then
      allocate (gf% z (gf% nz+1, gf% n))
      gf% z = 0  ! Not yet implemented...
    end if

    gf% i = 0
    do j=1, gf% n
      c => cols% col(j)
      if (associated (c% t) .and. &
          associated (c% q) .and. &
          associated (c% p) .and. &
          associated (c% ph)      ) then
        gf% i   (                c% i,         &
                                 c% j,         &
                                 c% l ) = j       ! index array
        gf% s   (j)% i1    =     c% i             ! horizontal index 1
        gf% s   (j)% i2    =     c% j             ! horizontal index 2
        gf% s   (j)% id    =     c% l             ! diamond index
        gf% t (:,j)        =     c% t             ! temperature      [K]
        gf% q (:,j)        =     c% q             ! spec.humid.  [kg/kg]
        if (use_waterload) then
            gf% qx(:,j)    =     c% x             ! water load   [kg/kg]
        else
            gf% qx(:,j)    =     0._wp
        end if
        gf% p (:,j)        = exp(c% p)            ! full lev. pres. [Pa]
        gf% ph(:,j)        = exp(c% ph)           ! half lev. pres. [Pa]
        if (gf% vctype == VCT_P_HYB) &
            gf% ph(0,j)    =     0._wp            ! top level for hybrid p
#ifdef DEBUG_OCC
        if (gf% vctype /= VCT_P_HYB) &
        gf% z (:,j)        =     c% geoh          ! geop.height    [gpm]
#endif
        gf% s   (j)% Hsur  =     c% s% geosp      ! surface geop.  [gpm]
        gf% s   (j)% Gundu =     c% s% geoid      ! geoid undulation [m]
        gf% s   (j)% Psur  =     c% s% ps         ! surface press.  [Pa]
        gf% s   (j)% XLon  =     c% c% dlon       ! longitude      [deg]
        gf% s   (j)% XLat  =     c% c% dlat       ! latitude       [deg]
      else
        gf% s   (j)% XLat  = 0._wp
      end if
    end do

    if (present (xi)) then
      if (.not.present(obs)) &
        call finish ('set_occ_data','xi is present, obs not')
      do is = 1, obs% n_spot
        si   => obs% spot(is)
        ke21 = 2 * si% mke + 1
        if (si% hd% obstype /= OT_GPSRO) cycle
        if (si% use% state <= STAT_DISMISS) cycle
        if (si% n_spt * ke21 /= si% i% n) &
          call finish ('set_occ_data','n_spt*(2*ke+1)/=n')
        i = si% i% i
        do j = 1, si% n_spt
          k = si% imcol(j)% imc(1)
          call occ_xi2col (gf% t(:,k),      &
                           gf% q(:,k),      &
                           gf% s  (k)% Hsur,&
                           gf% p(:,k),      &
                           xi% x(i+1:i+ke21))
          i = i + ke21
        end do
      end do
    endif

  end subroutine set_occ_data
!------------------------------------------------------------------------------
! subroutine get_occ_data (gf, obs, xi)
! type (t_global)      ,intent(in)  :: gf
! type (t_obs)         ,intent(in)  :: obs
! type (t_vector_segm) ,intent(out) :: xi
! !-----------------------------------
! ! this routine is currently not used
! !-----------------------------------
!   integer                :: i, j, k
!   integer                :: is
!   integer                :: ke21
!   type (t_spot) ,pointer :: si
!
!   do is = 1, obs% n_spot
!     si   => obs% spot(is)
!     ke21 = 2 * si% mke + 1
!     if (si% hd% obstype /= OT_GPSRO) cycle
!     if (si% n_spt * ke21 /= si% i% n) &
!       call finish ('set_occ_data','n_spt*(2*ke+1)/=n')
!     i = si% i% i
!     do j = 1, si% n_spt
!         k = si% imcol(j)% imc(1)
!         call occ_col2xi (real(gf% t(:,k),      wp), &
!                          real(gf% q(:,k),      wp), &
!                          real(gf% s  (k)% Hsur,wp), &
!                          real(gf% p(:,k),      wp), &
!                          xi% x(i+1:i+ke21))
!       i = i + ke21
!     end do
!   end d
!
! end subroutine get_occ_data
!------------------------------------------------------------------------------
  subroutine set_occ_grid (gf, grid)
  type (t_global) ,intent (inout) :: gf
  type (t_grid)   ,pointer        :: grid
    !-------------
    ! set metadata
    !-------------
    gf% Grid_Type = grid% gridtype
    gf% vctype    = grid% vct           ! VCT_P_HYB (GME/IFS/HRM), else ICON
    gf% hint_mode = horint_mode         ! Horizontal interpolation mode
    gf% ref_model = refract_model       ! Refractivity model ident.
    gf% ref_scale = refract_scale       ! Refractivity scaling factor
    gf% nz        = grid% nz
    gf% nd        = grid% nd
    gf% ni        = grid% ni
    gf% ni2       = grid% ni2
    gf% nir       = grid% nir
    gf% lbg       = grid% lbg( (/1,2,4/) )
    gf% ubg       = grid% ubg( (/1,2,4/) )
    !-----------------------------------
    ! set subgrid size for interpolation
    !-----------------------------------
    select case (gf% Grid_Type)
    case (WMO6_GAUSSIAN, WMO6_LATLON)
      gf% ngp = ng1**2
    case (DWD6_ICOSAHEDRON, DWD6_ICON)
      gf% ngp = 3
    end select
    !-----------------------------------------------
    ! number of columns per profile in interpolation
    !-----------------------------------------------
    select case (horint_mode)
    case (0)
      gf% ncol = 1                      ! nearest neighbour mode
    case default
      gf% ncol = gf% ngp                ! same as subgrid
    end select
    !-------------------------------------
    ! set vertical coordinate coefficients
    !-------------------------------------
!   allocate (gf% a    (0:grid% nz))           ! a (0:nz) !!
!   allocate (gf% b    (0:grid% nz))           ! b (0:nz) !!
    allocate (gf% b_ad (0:grid% nz))
!   gf% a                  = grid% ak (1:grid% nz+1)
!   gf% b                  = grid% bk (1:grid% nz+1)
    gf% b_ad (:grid% nz-1) = 0._wp
    gf% b_ad ( grid% nz  ) = 1._wp
    gf% p0_msis            = p0_msis           ! MSIS base level pressure
    if (p0_msis <= 0._wp) then
       write(0,*) dace% pe, ": set_occ_grid: p0_msis =", p0_msis
       call finish ("set_occ_grid","p0_msis <= 0")
    end if
    !--------------------------
    ! set longitudes, latitudes
    !--------------------------
    select case (grid% gridtype)
    case (WMO6_LATLON, WMO6_GAUSSIAN)
      allocate (gf% glon (grid% nx))
      allocate (gf% glat (grid% ny))
      gf% glon = grid% dlon
      gf% glat = grid% dlat
    case (DWD6_ICOSAHEDRON, DWD6_ICON)
      gf% xnglob => grid% xnglob
      gf% marr   => grid% marr
      gf% grid   => grid
    end select
  end subroutine set_occ_grid
!------------------------------------------------------------------------------
  subroutine const_occ_data (gf, spot)
  type (t_global) ,intent (inout) :: gf
  type (t_spot)   ,intent (in)    :: spot
  !------------------------------------------
  ! test: set atmosphere to homogeneous state
  !------------------------------------------
    integer :: ix, i

    if (.not. homogeneous) return

    ix = sum (maxloc (spot% col% h% w))
    ix = spot% col% h% imc(ix,1)

    allocate (gf_back% t (size(gf% t ,1), size(gf% t ,2)))
    allocate (gf_back% q (size(gf% q ,1), size(gf% q ,2)))
    allocate (gf_back% qx(size(gf% qx,1), size(gf% qx,2)))
    allocate (gf_back% p (size(gf% p ,1), size(gf% p ,2)))
    allocate (gf_back% ph(size(gf% ph,1), size(gf% ph,2)))
    allocate (gf_back% s (size(gf% s )))
    gf_back%   t  = gf% t
    gf_back%   q  = gf% q
    gf_back%   qx = gf% qx
    gf_back%   p  = gf% p
    gf_back%   ph = gf% ph
    gf_back%   s  = gf% s
    if (gf% vctype /= VCT_P_HYB) then
      allocate (gf_back% z (size(gf% z ,1), size(gf% z ,2)))
      gf_back% z  = gf% z
    end if

    do i=1,size(spot% imcol)
      gf%   t (:,i) = gf% t (:,ix)
      gf%   q (:,i) = gf% q (:,ix)
      gf%   qx(:,i) = gf% qx(:,ix)
      gf%   p (:,i) = gf% p (:,ix)
      gf%   ph(:,i) = gf% ph(:,ix)
      gf%   s   (i) = gf% s   (ix)
      if (gf% vctype /= VCT_P_HYB) &
        gf% z (:,i) = gf% z (:,ix)
    end do

  end subroutine const_occ_data
!------------------------------------------------------------------------------
  subroutine restore_occ_data (gf)
  type (t_global) ,intent (inout) :: gf
  !-----------------------------------
  ! restore non-homogeneous atmosphere
  !-----------------------------------
    if (.not. homogeneous) return
    deallocate (gf% t, gf% q, gf% qx, gf% p, gf% ph, gf% s)
    gf% t   => gf_back% t
    gf% q   => gf_back% q
    gf% qx  => gf_back% qx
    gf% p   => gf_back% p
    gf% ph  => gf_back% ph
    gf% s   => gf_back% s
    if (gf% vctype /= VCT_P_HYB) then
      deallocate (gf% z)
      gf% z => gf_back% z
    end if
    nullify (gf_back% t, gf_back% q, gf_back% qx, gf_back% p, gf_back% ph, &
             gf_back% s, gf_back% z)
  end subroutine restore_occ_data
!==============================================================================
  !-------------------------------------------------------
  ! Derive pressure level for blending of MSIS climatology
  ! The transfer is done at the topmost level/layer.
  !-------------------------------------------------------
  subroutine setup_occ_msis (atm)
  type(t_atm) ,intent(in) :: atm      ! atmospheric state

    real(wp)                :: p0
    !--------------------------------------
    ! check for usage of Radio occultations
    !--------------------------------------
    if (rept_use(OT_GPSRO)% use(CHK_NONE) <= STAT_DISMISS) return

    select case (atm% grid% vct)
    case (VCT_P_ISO)
       p0 = atm% grid% akf(1)
    case (VCT_P_HYB)
       p0 = atm% grid% ak (2)
    case (VCT_Z_HYB, VCT_Z_GEN)
       if (associated (atm% ph)) then
         p0 = maxval (atm% ph(:,:,2,:))
       else if (associated (atm% pf)) then
         p0 = maxval (atm% pf(:,:,2,:))
       else
         call finish('setup_occ_msis','ph,pf not associated!')
       endif
       p0 = p_max  (p0)
    case default
       call finish('setup_occ_msis','invalid levtyp')
    end select
    if (p0 <= 0) call finish('setup_occ_msis','p0 <= 0')
    !-------------------------------------------
    ! In case of repeated calls keep the maximum
    !-------------------------------------------
    p0_msis = max (p0_msis, p0)
  end subroutine setup_occ_msis
!==============================================================================
! Processing of BUFR reports
!---------------------------
  subroutine read_gpsro_bufr (bufr, spt, obs, lkeep)

  type (t_bufr) ,intent(inout)        :: bufr  ! BUFR record to decode
  type (t_spot) ,intent(inout)        :: spt   ! meta information to set
  type (t_obs)  ,intent(inout)        :: obs   ! observations data type to set
  logical       ,intent(out)          :: lkeep ! accept observation ?

    !----------------
    ! index variables
    !----------------
    integer          :: is           ! sub-set    index
    integer          :: ie           ! entry in sub-set
    integer          :: id           ! descriptor index
    !-----------------------------------------
    ! quantities derived from the BUFR message
    !-----------------------------------------
    integer          :: ival  ! value decoded  from BUFR (integer)
    real(sp)         :: rval  ! value decoded  from BUFR (real)
    real(dp)         :: scale ! scaling factor from BUFR (double precision)
    character(len=8) :: ymnem ! value decoded from BUFR (string)
    integer          :: itype ! type of BUFR message

    integer          :: unused = 0

    integer          :: yyyy,mo,dd  ! occultation date YYYYMMDD
    integer          :: hh,mi,ss    ! occultation time HHMMSS
    integer          :: ifreq       ! frequency, 0 for ionospheric correction
    !----------------------------
    ! GPS RO data (derived types)
    !----------------------------
    type (t_occ)          :: occ     ! occultation meta data type
    type (t_ray) ,pointer :: r  (:)  ! ray meta data type
    integer               :: nlev    ! number of rays (step 1b levels)
    !-------------------------------------------
    ! estimate number of levels, allocate memory
    !-------------------------------------------
    lkeep = .false.
    if (bufr% sec3% num_subsets == 0) then
      call decr_rpt_use (spt, CHK_NOIMPL, &
                         comment='read_gpsro_bufr: cannot handle zero subsets')
      return
    endif
    if (bufr% sec3% num_subsets > 1) then
       print*,'!!!WARNING!!!   read_gpsro_bufr: read only 1st of several subsets'
       return
    endif
    is   = 1
    nlev = 0
    levelcount: do ie = 1, bufr% nbufdat(is)
      id     = bufr% idescidx(is,ie)
      itype  = bufr% itype(id)
      if (itype == ty_loop_cnt) then
        select case (bufr% ymnem(id))
        case ('Lcnt000')
          nlev = nlev + 1
        case ('Lcnt002')
          exit levelcount
        end select
      endif
    end do levelcount

    if(nlev < 1) then
      call decr_rpt_use (spt, CHK_INSDAT, comment='no. levels < 1')
      return
    endif

    allocate (r(nlev))
    !---------------------
    ! loop over data, copy
    !---------------------
    is    = 1
    nlev  = 0
    lkeep = .true.

    ss = 0
    readdata: do ie=1,bufr% nbufdat(is)
      !-------------
      ! decode datum
      !-------------
      ival  = bufr% ibufdat (is,ie)
      !----------------------------------------
      ! if valid, copy datum from BUFR record
      !----------------------------------------
!     if(ival /= inv_bufr) then
      if(.true.          ) then
        !---------------------------------------
        ! decode real/character datum, mnemonics
        !---------------------------------------
        id    = bufr% idescidx(is,ie)
        itype = bufr% itype(id)
        ymnem = bufr% ymnem(id)
        rval  = rvind
        scale = 1._dp
        if (bufr% is_char(id) == 0) then
          IF (ival /= inv_bufr) then
!           rval = ival * bufr% scale(id)
            ! Work around problem with double-rounding (before/after multiply)
            ! when ival > 2^24 and bufr%scale and rval are single-precision.
            rval = real (ival,wp) * bufr% scale(id)
            ! Higher accuracy for impact parameter etc. .
            ! Gives same results as in BUFR2NetCDF input files.
            ! Least sign. digit (10cm) of impact parameter still truncated
            ! due to rval declared as real(bp) instead of real(dp)
            scale = 10._dp**(-bufr% iscale(id))
            rval  = ival * scale
          endif
        endif

        select case (ymnem)
        !-------------------------
        ! Radio Occultation header
        !-------------------------
        case ('310026')  ! Sequence: SATELLITE RADIO OCCULTATION D[]
        case ('310022')  ! Sequence: SATEL.VERTIC.SOUND.DATA: TECH[]
        case ('MI1I2')   ! SATELLITE IDENTIFIER                   [CODE_TABLE]
          spt% ident      = ival
          spt% hd% satid  = ival
        case ('MSAIN')   ! SATELLITE INSTRUMENTS                  [CODE_TABLE]
          spt% sttyp  = ival
        case ('MMIOGC')  ! IDENTIFICATION OF ORIG./GENER. CENTRE  [CODE_TABLE]
          spt% hd% center = ival
        case ('MPTAG')   ! PRODUCT TYPE FOR RETRIEVED ATMOSH. GASE[CODE_TABLE]
          ! 2 = limb sounding
        case ('MSWID')   ! SOFTWARE IDENTIFICATION                [NUMERIC]
          ! uniquely defines processing algorithm used by org. centre
        !--------------------------
        ! Time of occultation start
        !--------------------------
        case ('MTISI')   ! TIME SIGNIFICANCE                      [CODE_TABLE]
          ! 17 = start of phenomenon
        case ('301011')  ! Sequence: DATE                         []
        case ('MJJJ')    ! YEAR                                   [YEAR]
          yyyy = ival
        case ('MMM')     ! MONTH                                  [MONTH]
          mo   = ival
        case ('MYY')     ! DAY                                    [DAY]
          dd   = ival
        case ('301012')  ! Sequence: HOUR/MINUTE                  []
        case ('MGG')     ! HOUR                                   [HOUR]
          hh   = ival
        case ('NGG')     ! MINUTE                                 [MINUTE]
          mi   = ival
        case ('MSEC')    ! SECOND                                 [SECOND]
          if (ival /= inv_bufr) ss = max (min (rval, 60._sp), 0._sp)
        !-------------------------------
        ! RO summary quality information
        !-------------------------------
        case ('NQROD')   ! QUALITY FLAGS FOR RADIO OCCULTATION DAT[FLAG_TABLE]
           occ% pcd% PCD = pcd_invalid
           if (ival /= inv_bufr) occ% pcd% PCD  = ival
        case ('MPCCO')   ! PER CENT CONFIDENCE                    [%]
        !--------
        ! LEO POD
        !--------
        case ('304030')  ! Sequence: LOCATION OF PLATFORM         []
        case ('MDFECL')  ! DIST. FROM EARTHS CENTRE, 0 DEGREES LON[M]
!         occ% rleo % x(1) = rval
        case ('MDFECE')  ! DIST. FROM EARTHS CENTRE, 90 DEGR. EAST[M]
!         occ% rleo % x(2) = rval
        case ('MDFECN')  ! DIST. FROM EARTHS CENTRE, NORTH POLE   [M]
!         occ% rleo % x(3) = rval
        case ('304031')  ! Sequence: SPEED OF PLATFORM            []
        case ('MAPVFC')  ! ABS. PLATFORM VELOCITY - 1. COMPONENT  [M/S]
!         occ% vleo % x(1) = rval
        case ('MAPVSC')  ! ABS. PLATFORM VELOCITY - 2. COMPONENT  [M/S]
!         occ% vleo % x(2) = rval
        case ('MAPVTC')  ! ABS. PLATFORM VELOCITY - 3. COMPONENT  [M/S]
!         occ% vleo % x(3) = rval
        !---------
        ! GNSS POD
        !---------
        case ('MSACL')   ! SATELLITE CLASSIFICATION               [CODE_TABLE]
          occ% gnssid      = ival
        case ('MPTID')   ! PLATFORM TRANSMITTER ID NUMBER         [NUMERIC]
          occ% prn         = ival
!       case ('304030')  ! Sequence: LOCATION OF PLATFORM         []
        case ('MDFECL0') ! DIST. FROM EARTHS CENTRE, 0 DEGREES LON[M]

        case ('MDFECE0') ! DIST. FROM EARTHS CENTRE, 90 DEGR. EAST[M]
        case ('MDFECN0') ! DIST. FROM EARTHS CENTRE, NORTH POLE   [M]
!       case ('304031')  ! Sequence: SPEED OF PLATFORM            []
        case ('MAPVFC0') ! ABS. PLATFORM VELOCITY - 1. COMPONENT  [M/S]
        case ('MAPVSC0') ! ABS. PLATFORM VELOCITY - 2. COMPONENT  [M/S]
        case ('MAPVTC0') ! ABS. PLATFORM VELOCITY - 3. COMPONENT  [M/S]
        !-----------------------
        ! Local Earth parameters
        !-----------------------
        case ('MSETI')   ! TIME INCREMENT                         [SECOND]
        case ('301021')  ! Sequence: LATITUDE/LONGITUDE           []
        case ('MLAH')    ! LATITUDE (HIGH ACCURACY)               [DEGREE]
           occ% gp%  phi    = invalid
           if (ival /= inv_bufr) occ% gp%  phi = rval
        case ('MLOH')    ! LONGITUDE (HIGH ACCURACY)              [DEGREE]
           occ% gp%  lambda = invalid
           if (ival /= inv_bufr) occ% gp%  lambda = rval
        case ('MDFECL1') ! DIST. FROM EARTHS CENTRE, 0 DEGREES LON[M]
           occ% xlc% x(1)  = invalid
           if (ival /= inv_bufr) occ% xlc% x(1) = rval * mtokm ! ->  [KM]
        case ('MDFECE1') ! DIST. FROM EARTHS CENTRE, 90 DEGR. EAST[M]
           occ% xlc% x(2)    = invalid
           if (ival /= inv_bufr) occ% xlc% x(2) = rval * mtokm ! ->  [KM]
        case ('MDFECN1') ! DIST. FROM EARTHS CENTRE, NORTH POLE   [M]
            occ% xlc% x(3)   = invalid
           if (ival /= inv_bufr) occ% xlc% x(3) = rval * mtokm ! ->  [KM]
        case ('MELRC')   ! EARTH'S LOCAL RADIUS OF CURVATURE      [M]
           occ% rlc        = invalid
           if (ival /= inv_bufr) occ% rlc = rval * mtokm       ! ->  [KM]
        case ('MDA')     ! BEARING OR AZIMUTH                     [DEGREE_TRUE]
           occ% bearing        = invalid
           if (ival /= inv_bufr) occ% bearing     = rval
        case ('NGEUN')   ! GEOID UNDULATION                       [M]
        !----------------
        ! RO Step 1b data
        !----------------
        case ('Loop000') ! Start of Loop - 113000                 []
        case ('MEDRE')   ! EXTENDED DELAYED DESCRIPT.REPLIC.FACTOR[NUMERIC]
        case ('Lcnt000') ! Loop Counter                           []
          nlev = nlev + 1
          if (nlev > size(r)) call finish('read_gpsro_bufr','nlev > size(r)')
        case ('MLAH0')   ! LATITUDE (HIGH ACCURACY)               [DEGREE]
             if (ival /= inv_bufr) r(nlev)% geo% phi    = rval
        case ('MLOH0')   ! LONGITUDE (HIGH ACCURACY)              [DEGREE]
             if (ival /= inv_bufr) r(nlev)% geo% lambda = rval
        case ('MDA0')    ! BEARING OR AZIMUTH                     [DEGREE_TRUE]
          r(nlev)% bearing           = rval
        case ('Loop001') ! Start of Loop - 108000                 []
          ifreq = - 999
        case ('MDREP')   ! DELAYED DESCRIPTOR REPLICATION FACTOR  [NUMERIC]
        case ('Lcnt001') ! Loop Counter                           []
        case ('MMFREQ')  ! MEAN FREQUENCY                         [1/S]
          ifreq = ival
        case ('MIMPA')   ! IMPACT PARAMETER                       [M]
          if (ifreq==0) then
             r(nlev)% p = invalid
             if (ival /= inv_bufr) r(nlev)% p   = rval * mtokm   ! ->  [KM]
          endif
        case ('NBEAN')   ! BENDING ANGLE                          [RADIANS]
          if (ifreq==0) then
             r(nlev)% eps = 0._wp
             if (ival /= inv_bufr) r(nlev)% eps = rval
          endif
        case ('MFOST')   ! FIRST ORDER STATISTICS                 [CODE_TABLE]
        case ('NBEAN0')  ! BENDING ANGLE                          [RADIANS]
          if (ifreq==0) then
             r(nlev)% var = 0._wp
             if (ival /= inv_bufr) r(nlev)% var = rval * rval
          endif
        case ('MFOST0')  ! FIRST ORDER STATISTICS                 [CODE_TABLE]
        case ('MPCCO0')  ! PER CENT CONFIDENCE                    [%]
          r(nlev)% pcc = pcc_invalid
          if (ival /= inv_bufr) r(nlev)% pcc = ival
        !----------------
        ! RO Step 2a data
        !----------------
        case ('Loop002') ! Start of Loop - 108000                 []
#ifndef CHECKCODES
        case default
          !-------------------------
          ! ignore subsequent fields
          !-------------------------
          exit readdata  ! ignore subsequent fields
#else
        case ('MEDRE0')  ! EXTENDED DELAYED DESCRIPT.REPLIC.FACTOR[NUMERIC]
        case ('Lcnt002') ! Loop Counter                           []
        case ('MH')      ! HEIGHT                                 [M]
        case ('NATRE')   ! ATMOSPHERIC REFRACTIVITY               [N_UNITS]
        case ('MFOST1')  ! FIRST ORDER STATISTICS                 [CODE_TABLE]
        case ('NATRE0')  ! ATMOSPHERIC REFRACTIVITY               [N_UNITS]
        case ('MFOST2')  ! FIRST ORDER STATISTICS                 [CODE_TABLE]
        case ('MPCCO1')  ! PER CENT CONFIDENCE                    [%]
        !----------------
        ! RO Step 2b data
        !----------------
        case ('Loop003') ! Start of Loop - 116000                 []
        case ('MEDRE1')  ! EXTENDED DELAYED DESCRIPT.REPLIC.FACTOR[NUMERIC]
        case ('Lcnt003') ! Loop Counter                           []
        case ('NHHH')    ! GEOPOTENTIAL HEIGHT                    [GPM]
        case ('MPPP')    ! PRESSURE                               [PA]
        case ('MTN')     ! TEMPERATURE/DRY BULB TEMPERATURE       [K]
        case ('MSPHU')   ! SPECIFIC HUMIDITY                      [KG/KG]
        case ('MFOST3')  ! FIRST ORDER STATISTICS                 [CODE_TABLE]
        case ('MPPP0')   ! PRESSURE                               [PA]
        case ('MTN0')    ! TEMPERATURE/DRY BULB TEMPERATURE       [K]
        case ('MSPHU0')  ! SPECIFIC HUMIDITY                      [KG/KG]
        case ('MFOST4')  ! FIRST ORDER STATISTICS                 [CODE_TABLE]
        case ('MPCCO2')  ! PER CENT CONFIDENCE                    [%]
        !----------------
        ! RO Step 2c data
        !----------------
        case ('MVTSA')   ! VERTICAL SIGNIFICANCE (SATELL.OBSERV.) [CODE_TABLE]
        case ('NHHH0')   ! GEOPOTENTIAL HEIGHT                    [GPM]
        case ('MPPP1')   ! PRESSURE                               [PA]
        case ('MFOST5')  ! FIRST ORDER STATISTICS                 [CODE_TABLE]
        case ('MPPP2')   ! PRESSURE                               [PA]
        case ('MFOST6')  ! FIRST ORDER STATISTICS                 [CODE_TABLE]
        case ('MPCCO3')  ! PER CENT CONFIDENCE                    [%]
        !---------------
        ! codes ignored:
        !---------------
        case ('End Seq') ! End Sequence
        case ('EndLoop') ! End Loop
        !------------------------
        ! check for unknown codes
        !------------------------
        case default
          call bufr_get_entry_texts (bufr)
          call bufr_get_entry_units (bufr)
          write(0,'(a,a,a)') ymnem, bufr% ytext(id),      &
                        '['//trim(bufr% yunit(id))//']'
          unused = unused + 1
          if (unused > 50) exit readdata
#endif
        end select
      end if
    end do readdata
    !----------------------------------------------------------
    ! abort if unkown codes are found or wrong number of levels
    !----------------------------------------------------------
    if (unused /= 0)   call finish('read_gpsro_bufr','code(s) not implemented')
    if (nlev/=size(r)) call finish('read_gpsro_bufr','nlev /= size(r)')

    !--------------------------------------------
    ! store information in occ specific data type
    !--------------------------------------------
    call init_time (occ% time, yyyy, mo, dd, hh, mi, ss)
    occ% file      = source(spt% hd% source)% file
    occ% occid     = hh * 100 + mi
    occ% leoid     = mnem (spt% ident)
    occ% gp1d      = occ% gp
    occ% nray      = nlev
    spt% col% nlev = nlev
    spt% phase     = occ% pcd% PCD
    write(spt% statid,'(a4,i4.4)') occ% leoid, occ% occid

    !-------------------------------------------
    ! Final preparation of observation type data
    ! Storage into components of 'obs'
    !-------------------------------------------
    call check_store_occ (occ, r, spt, obs, lkeep)

    deallocate (r)

  end subroutine read_gpsro_bufr
!==============================================================================
  subroutine read_gpsro_netcdf (ifile, i_source, obs, lkeep, nkeep)
  !----------------------------------------------------------
  ! Read GPSRO observations from netCDF (converted BUFR data)
  !----------------------------------------------------------
  integer      ,intent(in)    :: ifile    ! Number of netCDF file to read
  integer      ,intent(inout) :: i_source ! record number in source-file
  type (t_obs) ,intent(inout) :: obs      ! observations to set
  logical      ,intent(out)   :: lkeep    ! any observation accepted ?
  integer      ,intent(out)   :: nkeep    ! number of accepted observations

    !================
    ! local variables
    !================
    logical             :: lk                  ! observation accepted ?
    type (t_obsid)      :: obsid               ! observation id table entry
    !-----------
    ! dimensions
    !-----------
    integer       :: nrep         ! number of reports
    integer       :: mlev         ! max. number of levels
    integer       :: nfreq        ! number of frequencies
    integer       :: nlev         ! number of rays
    integer       :: i, j, k, i0  ! loop indices
    !-----------------
    ! temporary arrays
    !-----------------
    type (t_spot) ,allocatable :: spt      (:) ! 3dvar report meta data
    type (t_occ)  ,allocatable :: occ      (:) ! occultation meta data type
    type (t_ray)  ,allocatable :: rays   (:,:) ! ray meta data type
    type (t_ray)  ,pointer     :: r      (:)   ! ray meta data type
    integer       ,allocatable :: ifreq  (:,:) ! index of L1+L2 'frequency'
    real(wp)      ,allocatable :: tmp  (:,:,:) ! temporary

    !--------------------------------------------------------
    ! GNSS identification for GPS, GLONASS, Galileo, BDS, ...
    ! using satellite codes from SINEX TRO V2.00
    !--------------------------------------------------------
!   character(1),  parameter   :: gnss_ident(401:*) = ["G","R","E","C"]
    ! workaround for AMD flang and old PGI:
    character(1),  parameter   :: gnss_ident(401:404) = ["G","R","E","C"]

    !----------------------------------------
    ! get dimensions of fields in NetCDF file
    !----------------------------------------
    ncid  =  cdfinid
    nfreq =  dimlen ('MIMPA',1) ! number of frequencies
    mlev  =  dimlen ('MIMPA',2) ! number of levels
    nrep  =  size   (s1date)    ! number of reports

    !---------------------
    ! allocate temporaries
    !---------------------
    allocate (spt                (nrep))
    allocate (occ                (nrep))
    allocate (rays         (mlev ,nrep))
    allocate (ifreq        (mlev ,nrep))
    allocate (tmp   (nfreq ,mlev ,nrep))

    !-----------
    ! set header
    !-----------
    spt% hd% dbkz        = s2ikz         ! DWD-internal classifier
    spt% hd% satid       = istidn        ! satellite identifier
    spt% hd% obstype     = OT_GPSRO      ! observation type
!   spt% hd% codetype    =
    spt% hd% modtype     = GPSRO         ! module type
    spt% hd% buf_type    = s1cat         ! data category
    spt% hd% buf_subtype = s1catls       ! local data sub category
    spt% hd% time        = stime         ! header observation time (section1)
    spt% hd% db_time     = db_time       ! data bank time
    spt% hd% source      = ifile         ! file index
    i0                   = sum (source(1:ifile-1)% entries)
    spt% hd% id          = (/(i,i=i0+1,i0+nrep)/)
    spt% hd% record      = (/(i,i=   1,   nrep)/)
    spt% hd% center      = s1cent        ! data centre
    spt% hd% subcenter   = s1cents       ! data subcentre

    !---------------------------------
    ! derive information from DWD dbkz
    !---------------------------------
    if (s2ikz(1) >= 0) then
       obsid             = obstype_dbkz (s2ikz(1))
       spt% hd% codetype = obsid% codetype
    end if

    spt% col% c% dlat    = mlah          ! latitude
    spt% col% c% dlon    = mloh          ! longitude
    spt% actual_time     = obs_time      ! time of occultation start
    spt% ident           = istidn        ! satellite identifier
    spt% col% nlev       = get_int ('MEDRE' ,0 ,nrep) ! number of levels
    spt% sttyp           = get_int ('MSAIN' ,-1,nrep) ! satellite instruments
    spt% center_id       = get_int ('MMIOGC',-1,nrep) ! processing center

    occ% time       = obs_time            ! time of occultation start
    occ% file       = source(ifile)% file ! file name
    occ% occid      = mgg * 100 + ngg     ! hhmm
    occ% nray       = spt% col% nlev

    !-------------------------------------------------------------------------
    ! not used:
    ! MPTAG (2 = limb sounding)        PRODUCT TYPE FOR RETRIEVED ATMOSH. GASE
    ! MSWID                            SOFTWARE IDENTIFICATION
    ! MTISI (17 = start of phenomenon) TIME SIGNIFICANCE
    ! MPCCO                            PER CENT CONFIDENCE [%]
    !-------------------------------------------------------------------------
    !-------------------------------
    ! RO summary quality information
    !-------------------------------
    occ% pcd% PCD = get_int ('NQROD', pcd_invalid, nrep)
    !-----------------------------
    ! LEO POD (currently not used)
    !-----------------------------
    ! MDFECL -> occ% rleo % x(1)  DIST. FROM EARTHS CENTRE, 0 DEGREES LON[M]   -> km
    ! MDFECE -> occ% rleo % x(2)  DIST. FROM EARTHS CENTRE, 90 DEGR. EAST[M]   -> km
    ! MDFECN -> occ% rleo % x(3)  DIST. FROM EARTHS CENTRE, NORTH POLE   [M]   -> km
    ! MAPVFC -> occ% vleo % x(1)  ABS. PLATFORM VELOCITY - 1. COMPONENT  [M/S] -> km/s
    ! MAPVSC -> occ% vleo % x(2)  ABS. PLATFORM VELOCITY - 2. COMPONENT  [M/S] -> km/s
    ! MAPVTC -> occ% vleo % x(3)  ABS. PLATFORM VELOCITY - 3. COMPONENT  [M/S] -> km/s
    !
    occ% rleo% x(1)  = get_real ('MDFECL' ,invalid ,mtokm ,nrep)
    occ% rleo% x(2)  = get_real ('MDFECE' ,invalid ,mtokm ,nrep)
    occ% rleo% x(3)  = get_real ('MDFECN' ,invalid ,mtokm ,nrep)
    occ% vleo% x(1)  = get_real ('MAPVFC' ,invalid ,mtokm ,nrep)
    occ% vleo% x(2)  = get_real ('MAPVSC' ,invalid ,mtokm ,nrep)
    occ% vleo% x(3)  = get_real ('MAPVTC' ,invalid ,mtokm ,nrep)

    !---------
    ! GNSS POD
    !---------
    !     used: MDFECL0 DIST. FROM EARTHS CENTRE, 0 DEGREES LON[M]   -> km
    !           MDFECE0 DIST. FROM EARTHS CENTRE, 90 DEGR. EAST[M]   -> km
    !           MDFECN0 DIST. FROM EARTHS CENTRE, NORTH POLE   [M]   -> km
    !           MAPVFC0 ABS. PLATFORM VELOCITY - 1. COMPONENT  [M/S] -> km/s
    !           MAPVSC0 ABS. PLATFORM VELOCITY - 2. COMPONENT  [M/S] -> km/s
    !           MAPVTC0 ABS. PLATFORM VELOCITY - 3. COMPONENT  [M/S] -> km/s
    !           MSACL   SATELLITE CLASSIFICATION
    !           MPTID   PLATFORM TRANSMITTER ID NUMBER
    !
    occ% rgps% x(1)  = get_real ('MDFECL0' ,invalid ,mtokm ,nrep)
    occ% rgps% x(2)  = get_real ('MDFECE0' ,invalid ,mtokm ,nrep)
    occ% rgps% x(3)  = get_real ('MDFECN0' ,invalid ,mtokm ,nrep)
    occ% vgps% x(1)  = get_real ('MAPVFC0' ,invalid ,mtokm ,nrep)
    occ% vgps% x(2)  = get_real ('MAPVSC0' ,invalid ,mtokm ,nrep)
    occ% vgps% x(3)  = get_real ('MAPVTC0' ,invalid ,mtokm ,nrep)

    occ% gnssid      = get_int  ('MSACL', -1 ,nrep)
    occ% prn         = get_int  ('MPTID', -1 ,nrep)
    spt% tracking    = occ% gnssid
    spt% sender_id   = occ% prn

    !-----------------------
    ! Local Earth parameters
    !-----------------------
    !     used: MSETI    TIME INCREMENT                         [SECOND]
    !           MLAH     LATITUDE (HIGH ACCURACY)               [DEGREE]
    !           MLOH     LONGITUDE (HIGH ACCURACY)              [DEGREE]
    !           MDFECL1  DIST. FROM EARTHS CENTRE, 0 DEGREES LON[M] -> km
    !           MDFECE1  DIST. FROM EARTHS CENTRE, 90 DEGR. EAST[M] -> km
    !           MDFECN1  DIST. FROM EARTHS CENTRE, NORTH POLE   [M] -> km
    !           MELRC    EARTH'S LOCAL RADIUS OF CURVATURE      [M] -> km
    !           MDA      BEARING OR AZIMUTH                     [DEGREE_TRUE]
    ! not used: NGEUN    GEOID UNDULATION                       [M]
    !
    occ% timeinc     = get_real ('MSETI'   ,invalid ,1._wp ,nrep)
    occ% gp%  phi    = mlah
    occ% gp%  lambda = mloh
    occ% gp1d        = occ% gp
    occ% xlc% x(1)   = get_real ('MDFECL1' ,invalid ,mtokm ,nrep)
    occ% xlc% x(2)   = get_real ('MDFECE1' ,invalid ,mtokm ,nrep)
    occ% xlc% x(3)   = get_real ('MDFECN1' ,invalid ,mtokm ,nrep)
    occ% rlc         = get_real ('MELRC'   ,invalid ,mtokm ,nrep)
    occ% bearing     = get_real ('MDA'     ,invalid ,1._wp ,nrep)

    !============================
    ! Read radio occultation body
    !============================
    !----------------
    ! RO Step 1b data
    !----------------
    !
    ! not used: 'MDREP'  DELAYED DESCRIPTOR REPLICATION FACTOR  [NUMERIC]
    !           'MFOST'  FIRST ORDER STATISTICS                 [CODE_TABLE]
    !           'MFOST0' FIRST ORDER STATISTICS                 [CODE_TABLE]
    !           'MPCCO0' PER CENT CONFIDENCE                    [%]
    !
    !     used: 'MEDRE'  EXTENDED DELAYED DESCRIPT.REPLIC.FACTOR[NUMERIC]
    !           'MLAH0'  LATITUDE (HIGH ACCURACY)               [DEGREE]
    !           'MLOH0'  LONGITUDE (HIGH ACCURACY)              [DEGREE]
    !           'MDA0'   BEARING OR AZIMUTH                     [DEGREE_TRUE]
    !           'MMFREQ' MEAN FREQUENCY                         [1/S]
    !           'MIMPA'  IMPACT PARAMETER                       [M] -> km
    !           'NBEAN'  BENDING ANGLE                          [RADIANS]
    !           'NBEAN0' BENDING ANGLE error                    [RADIANS]
    !
    occ %      nray    = spt% col% nlev
    rays% geo% phi     = get_real2 ('MLAH0' ,invalid ,1._wp ,mlev ,nrep)
    rays% geo% lambda  = get_real2 ('MLOH0' ,invalid ,1._wp ,mlev ,nrep)
    rays%      bearing = get_real2 ('MDA0'  ,invalid ,1._wp ,mlev ,nrep)
    !-----------------------------------------
    ! pick up frequency '0' (combined L1 + L2)
    !-----------------------------------------
    tmp   = get_real3 ('MMFREQ' ,0._wp ,1._wp ,nfreq ,mlev ,nrep)
    ifreq = 0
    do j = 1, nfreq
      where (tmp (j,:,:) == 0._wp) ifreq = j
    end do
    !--------------------------------
    ! combined L1+L2 impact parameter
    !--------------------------------
    tmp   = get_real3 ('MIMPA' ,invalid ,mtokm ,nfreq ,mlev ,nrep)
    rays% p = invalid
    do j = 1, nrep
      do i = 1, mlev
        k = ifreq (i,j)
        if (k>0) rays(i,j)% p = tmp (k,i,j)
      end do
    end do
    !-----------------------------
    ! combined L1+L2 bending angle
    !-----------------------------
    tmp   = get_real3 ('NBEAN' ,0._wp ,1._wp ,nfreq ,mlev ,nrep)
    rays% eps = 0._wp
    do j = 1, nrep
      do i = 1, mlev
        k = ifreq (i,j)
        if (k>0) rays(i,j)% eps = tmp (k,i,j)
      end do
    end do
    !----------------------------------------------
    ! combined L1+L2 bending angle error (variance)
    !----------------------------------------------
    tmp   = get_real3 ('NBEAN0' ,0._wp ,1._wp ,nfreq ,mlev ,nrep)
    rays% var = 0._wp
    do j = 1, nrep
      do i = 1, mlev
        k = ifreq (i,j)
        if (k>0) rays(i,j)% var = tmp (k,i,j) ** 2
      end do
    end do
    !--------------------
    ! Per Cent Confidence
    !--------------------
    rays% pcc = get_int2 ('MPCCO0', int (pcc_invalid) ,1 ,mlev ,nrep)
    !---------------------------
    ! RO Step 2a data (not used)
    !---------------------------
    ! 'MEDRE0' EXTENDED DELAYED DESCRIPT.REPLIC.FACTOR[NUMERIC]
    ! 'MH'     HEIGHT                                 [M]
    ! 'NATRE'  ATMOSPHERIC REFRACTIVITY               [N_UNITS]
    ! 'MFOST1' FIRST ORDER STATISTICS                 [CODE_TABLE]
    ! 'NATRE0' ATMOSPHERIC REFRACTIVITY               [N_UNITS]
    ! 'MFOST2' FIRST ORDER STATISTICS                 [CODE_TABLE]
    ! 'MPCCO1' PER CENT CONFIDENCE                    [%]
    !---------------------------
    ! RO Step 2b data (not used)
    !---------------------------
    ! 'MEDRE1' EXTENDED DELAYED DESCRIPT.REPLIC.FACTOR[NUMERIC]
    ! 'NHHH'   GEOPOTENTIAL HEIGHT                    [GPM]
    ! 'MPPP'   PRESSURE                               [PA]
    ! 'MTN'    TEMPERATURE/DRY BULB TEMPERATURE       [K]
    ! 'MSPHU'  SPECIFIC HUMIDITY                      [KG/KG]
    ! 'MFOST3' FIRST ORDER STATISTICS                 [CODE_TABLE]
    ! 'MPPP0'  PRESSURE                               [PA]
    ! 'MTN0'   TEMPERATURE/DRY BULB TEMPERATURE       [K]
    ! 'MSPHU0' SPECIFIC HUMIDITY                      [KG/KG]
    ! 'MFOST4' FIRST ORDER STATISTICS                 [CODE_TABLE]
    ! 'MPCCO2' PER CENT CONFIDENCE                    [%]
    !---------------------------
    ! RO Step 2c data (not used)
    !---------------------------
    ! 'MVTSA'  VERTICAL SIGNIFICANCE (SATELL.OBSERV.) [CODE_TABLE]
    ! 'NHHH0'  GEOPOTENTIAL HEIGHT                    [GPM]
    ! 'MPPP1'  PRESSURE                               [PA]
    ! 'MFOST5' FIRST ORDER STATISTICS                 [CODE_TABLE]
    ! 'MPPP2'  PRESSURE                               [PA]
    ! 'MFOST6' FIRST ORDER STATISTICS                 [CODE_TABLE]
    ! 'MPCCO3' PER CENT CONFIDENCE                    [%]

    !-----------------------------------
    ! final loop over individual reports
    !-----------------------------------
    nkeep = 0
    do i = 1, nrep
      !--------------------------------
      ! set remaining meta data entries
      !--------------------------------
      spt(i)% hd% idbk  = idb_dbk (spt(i)% hd% dbkz, OT_GPSRO)
      if(occ(i)% pcd% PCD < 32768) then
        spt(i)%   phase = occ(i)% pcd% PCD
      else
        spt(i)%   phase = occ(i)% pcd% PCD - 65536
      endif
      occ(i)%     leoid = mnem (spt(i)% ident)
      write(spt(i)% statid,'(a4,i4.4)') occ(i)% leoid, occ(i)% occid
      if (use_gnssid) then
         ! Annotate all but not GPS occultations
         if (occ(i)% gnssid >  lbound (gnss_ident,1) .and. &
             occ(i)% gnssid <= ubound (gnss_ident,1)       ) then
            spt(i)% statid(9:10) = "_" // gnss_ident(occ(i)% gnssid)
         end if
      end if
      !--------------------------------------------
      ! perform simple generic check on report type
      !--------------------------------------------
      call check_report_0 (spt(i)% use, spt(i)% hd, 1)
      if (spt(i)% use% state <= STAT_DISMISS) cycle

#if defined (__SUNPRO_F90) || defined (__SUNPRO_F95)
      ! Workaround for optimization(?) bug in Sun f95 (SunStudio 12.0)
      if (spt(i)% use% state < 0) &
      write(*,*) "read_gpsro_netcdf:check_report_0: use =", spt(i)% use% state
#endif
      !---------------------
      ! check number of rays
      !---------------------
      nlev = spt(i)% col% nlev
      if (nlev < 1) then
        call decr_rpt_use (spt(i), CHK_INSDAT,                      &
                           comment='read_gpsro_netcdf: no. rays < 1')
        cycle
      endif

      !----------------
      ! standard checks
      !----------------
      call check_report_1 (spt(i))
      if (spt(i)% use% state <= STAT_DISMISS) cycle

#if defined (__SUNPRO_F90) || defined (__SUNPRO_F95)
      ! Diagnose bug in Sun f95 (SunStudio 12.0)
      if (spt(i)% use% state < 1) &
      write(*,*) "read_gpsro_netcdf:check_report_1: use =", spt(i)% use% state
#endif
      !-------------------------------------------
      ! Final preparation of observation type data
      ! Storage into components of 'obs'
      !-------------------------------------------
      allocate (r (nlev))
      r = rays (1:nlev,i)
      call check_store_occ (occ(i), r, spt(i), obs, lk)
      deallocate (r)
      if (lk) nkeep = nkeep + 1
    end do

    lkeep = nkeep > 0

  contains
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function get_int (name, missing, n) result (values)
    !----------------------------
    ! decode integers from header
    !----------------------------
    integer          ,intent(in) :: n           ! number of reports
    character(len=*) ,intent(in) :: name        ! name of variable
    integer          ,intent(in) :: missing     ! replace fillvalues with
    integer                      :: values (n)  ! decoded values with

      integer :: fillvalue
      integer :: status
      integer :: varid

      status = nf90_inq_varid (ncid, name, varid)
      if (status /= NF90_NOERR) then
        values = missing
      else
        status = nf90_get_att (ncid, varid, '_FillValue', fillvalue)
        status = nf90_get_var (ncid, varid, values)
        if    (status /= NF90_NOERR) then
          values = missing
        else
          where (values == fillvalue) values = missing
        endif
      endif

    end function get_int
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function get_int2 (name, missing, factor, l, n) result (values)
    !--------------------------------
    ! decode integer values from body
    !--------------------------------
    integer          ,intent(in) :: l           ! number of levels
    integer          ,intent(in) :: n           ! number of reports
    character(len=*) ,intent(in) :: name        ! name of variable
    integer          ,intent(in) :: missing     ! replace fillvalues with
    integer          ,intent(in) :: factor      ! multiply valid entries with
    integer                      :: values (l,n)! decoded values

      integer  :: fillvalue
      integer  :: status
      integer  :: varid

      status = nf90_inq_varid (ncid, name, varid)
      if (status /= NF90_NOERR) then
        values = missing
      else
        status = nf90_get_att (ncid, varid, '_FillValue', fillvalue)
        status = nf90_get_var (ncid, varid, values)
        if    (status /= NF90_NOERR) then
          values = missing
        else
          where (values == fillvalue)
            values = missing
          elsewhere
            values = values * factor
          endwhere
        endif
      endif

    end function get_int2
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function get_real (name, missing, factor, n) result (values)
    !-------------------------------
    ! decode real values from header
    !-------------------------------
    integer          ,intent(in) :: n           ! number of reports
    character(len=*) ,intent(in) :: name        ! name of variable
    real(wp)         ,intent(in) :: missing     ! replace fillvalues with
    real(wp)         ,intent(in) :: factor      ! multiply valid entries with
    real(wp)                     :: values (n)  ! decoded values

      real(wp) :: fillvalue
      integer  :: status
      integer  :: varid

      status = nf90_inq_varid (ncid, name, varid)
      if (status /= NF90_NOERR) then
        values = missing
      else
        status = nf90_get_att (ncid, varid, '_FillValue', fillvalue)
        status = nf90_get_var (ncid, varid, values)
        if    (status /= NF90_NOERR) then
          values = missing
        else
          where (values == fillvalue)
            values = missing
          elsewhere
            values = values * factor
          endwhere
        endif
      endif

    end function get_real
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function get_real2 (name, missing, factor, l, n) result (values)
    !-----------------------------
    ! decode real values from body
    !-----------------------------
    integer          ,intent(in) :: l           ! number of levels
    integer          ,intent(in) :: n           ! number of reports
    character(len=*) ,intent(in) :: name        ! name of variable
    real(wp)         ,intent(in) :: missing     ! replace fillvalues with
    real(wp)         ,intent(in) :: factor      ! multiply valid entries with
    real(wp)                     :: values (l,n)! decoded values

      real(wp) :: fillvalue
      integer  :: status
      integer  :: varid

      status = nf90_inq_varid (ncid, name, varid)
      if (status /= NF90_NOERR) then
        values = missing
      else
        status = nf90_get_att (ncid, varid, '_FillValue', fillvalue)
        status = nf90_get_var (ncid, varid, values)
        if    (status /= NF90_NOERR) then
          values = missing
        else
          where (values == fillvalue)
            values = missing
          elsewhere
            values = values * factor
          endwhere
        endif
      endif

    end function get_real2
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function get_real3 (name, missing, factor, m, l, n) result (values)
    !-----------------------------------------
    ! decode real values (3D-fields) from body
    !-----------------------------------------
    integer          ,intent(in) :: m             ! innermost dimension len
    integer          ,intent(in) :: l             ! number of levels
    integer          ,intent(in) :: n             ! number of reports
    character(len=*) ,intent(in) :: name          ! name of variable
    real(wp)         ,intent(in) :: missing       ! replace fillvalues with
    real(wp)         ,intent(in) :: factor        ! multiply valid entries with
    real(wp)                     :: values (m,l,n)! decoded values

      real(wp) :: fillvalue
      integer  :: status
      integer  :: varid

      status = nf90_inq_varid (ncid, name, varid)
      if (status /= NF90_NOERR) then
        values = missing
      else
        status = nf90_get_att (ncid, varid, '_FillValue', fillvalue)
        status = nf90_get_var (ncid, varid, values)
        if    (status /= NF90_NOERR) then
          values = missing
        else
          where (values == fillvalue)
            values = missing
          elsewhere
            values = values * factor
          endwhere
        endif
      endif

    end function get_real3
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function dimlen (name, dim)
    !--------------------------------------------
    ! get dimension size of fields in NetCDF file
    ! return -1 in case of error
    !--------------------------------------------
    character(len=*) ,intent(in) :: name   ! name of variable
    integer          ,intent(in) :: dim    ! dimension index
    integer                      :: dimlen ! len of dimension

      integer :: status
      integer :: varid
      integer :: dimids (10)
      integer :: ndims

      dimlen = -1
      status = nf90_inq_varid         (ncid, name,  varid)
      if (status /= NF90_NOERR) return
      status = nf90_Inquire_Variable  (ncid, varid, ndims=ndims, dimids=dimids)
      if (status /= NF90_NOERR) return
      status = nf90_Inquire_Dimension (ncid, dimids(dim), len=dimlen)
      if (status /= NF90_NOERR) dimlen = -1

    end function dimlen
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine read_gpsro_netcdf
!==============================================================================
  subroutine check_occ_geom (occ, spot, rays, status)
    !---------------------------------------------
    ! Validate occultation geometry:
    ! Check metadata and fix values where possible
    ! Return value of status:
    ! 0: ok or fixes successfully applied
    ! 1: inconsistent data
    ! 2: insufficient or missing data
    !---------------------------------------------
    type(t_occ)  ,intent(inout)  :: occ
    type(t_spot) ,intent(in   )  :: spot
    type(t_ray)  ,intent(in   )  :: rays(:)
    integer      ,intent(  out)  :: status

    integer                :: yyyy,mm,dd,hh,mi,ss
    real(wp)               :: tmp           ! temporary
    real(wp), dimension(3) :: v1, v2        ! temporary 3-vectors
    real(wp), dimension(3) :: x_occ         ! occultation point
    real(wp), dimension(3) :: los           ! GPS -> LEO line of sight
    real(wp)               :: phi           ! Rotation angle ECF->ECI
    real(wp)               :: timeinc       ! Time increment since start
!   real(wp)               :: a1, a2        ! projections of sat. positions
    type(t_coord)          :: b             ! Local basis at occ. point
    type(cartesian)        :: rgps_e        ! GPS coordinates in Earth frame
    type(cartesian)        :: rleo_e        ! LEO coordinates in Earth frame
    type(cartesian)        :: vgps_e        ! GPS velocity    in Earth frame
    type(cartesian)        :: vleo_e        ! LEO velocity    in Earth frame
    type(cartesian)        :: xlc_e         ! Loc. of curvature, Earth frame
    type(cartesian)        :: xp, xps       ! Approximate perigee point
    logical                :: debug = .true.
    logical                :: l_eci2ecf     ! Transform from ECI->ECF

    type(Cartesian), parameter :: PA = Cartesian((/0,0,1/))  ! Polar axis

    debug = (verbose > 1)
    !-------------------------------------
    ! Check presence of essential metadata
    !-------------------------------------
    status = 2
    !---------------------------------------------------------------
    ! Header (georeferencing): Centre of curvature coordinates (ECF)
    !---------------------------------------------------------------
    if (all (occ% xlc % x == cart_inv% x)      ) return

    !---------------------------------------------
    ! Level 1a: satellite positions and velocities
    !---------------------------------------------
    if (all (occ% rgps% x == cart_inv% x) .or. &
        all (occ% rleo% x == cart_inv% x) .or. &
        all (occ% vgps% x == cart_inv% x) .or. &
        all (occ% vleo% x == cart_inv% x)      ) return

    if (all (occ% rgps% x == occ% rleo% x)) then
       print *, "*** Not funny!  (r_gps == r_leo)"
       print *, occ% rgps% x, cart_inv% x
       return
    end if

    status = 0

    !-------------------------------------------------------------
    ! Construct local basis of tangent space at Occultation Point:
    ! u = e_r, v ~ e_z .x. e_r, w = u .x. v
    !-------------------------------------------------------------
    b% dlon = occ% gp% Lambda
    b% dlat = occ% gp% Phi
    call set_xuv (b)

    !-------------------------------------------------
    ! Determine Greenwich Apparent Sidereal Time Angle
    ! for transformation between ECF and ECI frames.
    !-------------------------------------------------
    call i_time (occ% time, yyyy, mm, dd, hh, mi, ss)
    timeinc = occ% timeinc
    if (timeinc == invalid) then
       print *, "check_occ_geom: ", spot% statid, " : time increment not set!"
       timeinc = 0
    end if
    phi = GAST (yyyy, mm, dd, hh, mi, real (ss,wp), timeinc)

    if (debug) then
       print *
       print *, "check_occ_geom: ", spot% statid, "  ", &
            merge ("(Setting)", "(Rising) ", occ% PCD% occult_type == 0)
    end if

    rgps_e = occ% rgps
    rleo_e = occ% rleo
    xlc_e  = occ% xlc

    !------------------------------------------------
    ! Hack for MetOp-A data from EUMETSAT (not KNMI),
    ! generated by old geometric optics processing:
    ! apply transformation from ECI to ECF
    !------------------------------------------------
    l_eci2ecf = .false.
    if (spot% hd% center == WMO0_EUMET) then
       if (occ% pcd% open_loop == 0) l_eci2ecf = .true.
    end if

    if (l_eci2ecf) then
!      print *, "GAST =", phi
       rgps_e = Rotate (occ% rgps, PA, -phi)
       rleo_e = Rotate (occ% rleo, PA, -phi)
!      xlc_e  = Rotate (occ% xlc , PA, -phi) ! Unrotated agrees with DMI!
    end if

    !--------------------------------------------------
    ! Line-of-Sight GPS->LEO should be in tangent space
    ! (probably not true for some deep occultations!?)
    !--------------------------------------------------
    los = rgps_e% x - rleo_e% x
    tmp = cos_angle (los, b% x)
    if (abs (tmp) > 0.02) then
       print *, spot% statid, ": LOS check 1 failed, LOS not tangent (orthogonal)!?"
       print *, "Scalar product 0:", tmp
!      status = 1
    end if

    !-----------------------------------
    ! Check consistency of azimuth angle
    ! Construct vector in tangent space,
    ! rotated by angle MDA (clockwise)
    !-----------------------------------
    tmp = occ% bearing * d2r
    v1  = (-sin (tmp)) * b% du + (-cos (tmp)) * b% dv
    tmp = cos_angle (los, v1)
    if (tmp < 0.999) then
       print *, spot% statid, ": Azimuth check failed:", tmp
       status = 1
    end if

    if (.not. all (rays(:)% bearing == invalid)) then
       if (any (rays(:)% bearing == invalid)) then
          print *, spot% statid,": Some azimuth angles are undefined!"
       else
          tmp = maxval (abs (modulo (rays(:)% bearing + 180._wp &
                                     -   occ% bearing,  360._wp ) - 180._wp))
          if (tmp > 5) then
             print *,  spot% statid, ": Consistency check of azimuth angles:"
             print *, "Max. difference to reference value:", tmp
             status = 1
          end if
       end if
    end if

    !-----------------------------------------
    ! Consistency checks for occultation point
    !-----------------------------------------
!   x_occ     = occ% rlc * b% x + xlc_e% x      ! Estimate of tangent point
    x_occ     = R_Earth  * b% x                 ! Estimate of tangent point
    xp        = Perigee (rgps_e, rleo_e, 0._wp) ! Perigee from straight line
    occ% gp3d = Geod_from_Cart (xp)

!   tmp    = sqrt (sum (x_occ**2) / sum (xp% x**2))
    tmp    = R_Earth / sqrt (sum (xp% x**2))
    xps% x = tmp * xp% x                        ! Perigee, rescaled

    tmp = cos_angle (x_occ, xps% x)
    if (tmp < 0.999) then
       if (debug) then
          print *, "Scalar product 1:", tmp
!         print *, "xlc_0 =", occ% xlc
          print *, "xlc   =", xlc_e% x
          print *, "x_occ =", x_occ
          print *, "xp~   =", xps% x
!         print *, "|xp~|,  |x_occ| =", sqrt (sum (xps% x**2)), sqrt (sum (x_occ**2))
          print *, "|xp~  -  x_occ| =", sqrt (sum ((x_occ - xps% x)**2))
          print *, "|x_occ-xlc|,rlc =", sqrt (sum ((x_occ - xlc_e% x)**2)),occ% rlc
       end if
    end if

    !------------------------------------------
    ! Check collinearity of GPS->OP and OP->LEO
    !------------------------------------------
    v1  = x_occ - rleo_e% x
    v2  = rgps_e% x - x_occ
    tmp = cos_angle (v1, v2)

    if (tmp < 0.99) then
       print *, spot% statid, ": LOS check 2 failed, GPS->OP and OP->LEO not collinear!?"
       print *, "Scalar product 2:", tmp
!      status = 1
    end if

    !-------------------------------------------------
    ! Project satellite positions onto line through OP
    !-------------------------------------------------
!!$    v2 = xp% x / sqrt (sum (xp% x**2))
!!$    a1 = sum (b% x * rgps_e% x)
!!$    a2 = sum (b% x * rleo_e% x)
!!$    if (abs (a1 / occ% rlc - 1._wp) > 0.01_wp .or. &
!!$        abs (a2 / occ% rlc - 1._wp) > 0.01_wp ) then
!!$       if (debug) then
!!$          print *, "Local radius of curvature [km]:", occ% rlc
!!$          print *, "Projected impact parameter GPS:", a1, sum (v2 * rgps_e% x)
!!$          print *, "Projected impact parameter LEO:", a2, sum (v2 * rleo_e% x)
!!$       end if
!!$    end if

    !----------------------------------------------
    ! Transform satellite velocities to Earth frame
    !----------------------------------------------
    vgps_e    = Rotate (occ% vgps, PA, -phi)
    v1(1)     =   omcor * rgps_e% x(2)
    v1(2)     = - omcor * rgps_e% x(1)
    v1(3)     = 0
    vgps_e% x = vgps_e% x + v1

    vleo_e    = Rotate (occ% vleo, PA, -phi)
    v2(1)     =   omcor * rleo_e% x(2)
    v2(2)     = - omcor * rleo_e% x(1)
    v2(3)     = 0
    vleo_e% x = vleo_e% x + v2

    if (debug) then
       print *, "|vgps(inertial frame)| =", sqrt (sum (occ% vgps% x **2))
       print *, "|vgps(earth    frame)| =", sqrt (sum (vgps_e% x **2))
!      print *, "|vgps(frame rotation)| =", sqrt (sum (v1 **2))
       tmp = cos_angle (vgps_e% x, rgps_e% x)
       print *, " cos_angle (vgps,rgps) =", tmp
       print *, "|vleo(inertial frame)| =", sqrt (sum (occ% vleo% x **2))
       print *, "|vleo(earth    frame)| =", sqrt (sum (vleo_e% x **2))
!      print *, "|vleo(frame rotation)| =", sqrt (sum (v2 **2))
       tmp = cos_angle (vleo_e% x, rleo_e% x)
       print *, " cos_angle (vleo,rleo) =", tmp
    end if

    !--------------------------------------------------------
    ! Project relative velocity onto GPS -> LEO line of sight
    !--------------------------------------------------------
    if (debug) then
       tmp = cos_angle (los, vleo_e% x)
       print *, "Scalar product 3:", tmp, merge ("Setting", "Rising ", tmp < 0)
    end if
    tmp = cos_angle (los, vgps_e% x - vleo_e% x)
    if (debug) then
       print *, "Scalar product 4:", tmp, merge ("Setting", "Rising ", tmp > 0)
    end if

    if (occ% PCD% occult_type == 0) then
       !-----------------------------------
       ! Setting occultation: d must be > 0
       !-----------------------------------
       if (tmp < 0) then
          call message ('check_occ_geom',spot% statid//' rising occultation?')
          write(0,'(a,4f11.3)') 'rleo-rgps :', rleo_e% x - rgps_e% x
          write(0,'(a,4f11.3)') 'vleo-vgps :', vleo_e% x - vgps_e% x
          write(0,'(a,4f11.3)') 'cos(angle):', tmp
          status = 1
       end if
    else
       !-----------------------------------
       ! Rising occultation: d must be < 0
       !-----------------------------------
       if (tmp > 0) then
          call message ('check_occ_geom',spot% statid//' setting occultation?')
          write(0,'(a,4f11.3)') 'rleo-rgps :', rleo_e% x - rgps_e% x
          write(0,'(a,4f11.3)') 'vleo-vgps :', vleo_e% x - vgps_e% x
          write(0,'(a,f11.6)')  'cos(angle):', tmp
          status = 1
       end if
    end if

    !---------------------------------------------------
    ! Write back positions and velocities in Earth frame
    !---------------------------------------------------
    occ% rgps = rgps_e
    occ% rleo = rleo_e
!   occ% xlc  = xlc_e

    occ% vgps = vgps_e
    occ% vleo = vleo_e

  contains

    pure function cos_angle (r1, r2)
      real(wp)             :: cos_angle      ! Cosine of angle(r1,r2)
      real(wp), intent(in) :: r1(3), r2(3)

      cos_angle = sum (r1*r2) / sqrt (sum (r1*r1) * sum (r2*r2))
    end function cos_angle
  end subroutine check_occ_geom
!==============================================================================
  !-------------------------------------------------
  ! Called from routine reading BUFR or NetCDF files
  !   Final preparation of observation type data
  !   Storage into components of 'obs'
  !-------------------------------------------------
  subroutine check_store_occ (occ, rays, spot, obs, lkeep)
  type(t_occ)  ,intent(inout)  :: occ
  type(t_ray)  ,pointer        :: rays (:)
  type(t_spot) ,intent(inout)  :: spot
  type(t_obs)  ,intent(inout)  :: obs
  logical      ,intent(out)    :: lkeep ! keep or reject observation

!   integer   ,parameter :: mxvob = 5000 ! max. no. rays per occultation
    integer              :: ind(occ%nray)! indices to pick up
    integer              :: k, i         ! indices
    real(wp)             :: p_low        ! lowermost impact parameter
    integer              :: n_mean       ! no.levels used for mean perigeepoint
    real(wp)             :: p_1d         ! impact parameter closest to z_1d
    integer              :: nray         ! number of rays used
    integer              :: nold         ! number of rays select for smoothing
    real(wp)             :: z            ! height
    real(wp)             :: h            ! height [km]
    type(t_ray) ,pointer :: old (:)      ! temporary
    type(t_spot),pointer :: spt          ! temporary
!   real(wp)             :: p            ! temporary p (for interpolation)
    real(wp)             :: w, ww        ! temporary weight           ('')
    real(wp)             :: eps          ! temporary bending angle    ('')
    real(wp)             :: var          ! temporary variance
    real(wp)             :: eps2         ! temporary bending angle**2 ('')
    real(wp)             :: dlon         ! temporary longitude
    real(wp)             :: dlat         ! temporary latitude
    integer              :: pcc          ! temporary per cent confidence
    integer              :: k1,k2,k3     ! indices
    integer              :: id           ! observation id
    type(t_datum)        :: bod          ! body derived type
    logical              :: up           ! GPSRO-profile sorted bot->top
    real(wp)             :: x    (3)     ! Zenith vector, temporary
    real(wp)             :: x_ref(3)     ! Zenith vector of georeference point
    real(wp)             :: cut          ! Cosine of max.distance to georef.pt.
    real(wp)             :: htop         ! Highest ray level
    integer              :: u            ! use flag
    integer              :: status       ! occultation geometry status
    integer              :: ndup         ! number of duplicates

    lkeep = .true.

    !---------------------
    ! read RO quality flag
    !---------------------
    call read_quality(occ)
    if (occ% PCD% pcd /= pcd_invalid) then
       if (occult_type < 3) then
          if (occult_type /= 2**occ% PCD% occult_type) then
             call decr_rpt_use (spot, CHK_NOTUSED, STAT_PASSIVE, &
                                      comment='occult_type')
          endif
       endif
       if (occ% PCD% exsphs_proc  /= 0 .or. &
           occ% PCD% bangle_procc /= 0) then
          call decr_rpt_use (spot, CHK_QI)
       end if
    endif

    !-----------------------------------
    ! exit if essential metadata missing
    !-----------------------------------
    if (occ% rlc == invalid) then
      lkeep = .false.
      call decr_rpt_use (spot, CHK_INSDAT, comment='MELRC')
      return
    endif
    !------------------------------
    ! Select rays with valid
    ! values of bangle, error.
    ! Check consistency of metadata
    !------------------------------
    x_ref = unitvector (dlat=occ% gp% Phi, dlon=occ% gp% Lambda)
    cut   = cos (d2r * ray_dist_ref)
    nray = 0
    do k = 1, occ% nray
       dlat = rays(k)% geo% Phi
       dlon = rays(k)% geo% Lambda
       if ((  rays(k)% eps <= 0._wp  )   .or. &
            ( rays(k)% var <  0._wp  )   .or. &
            ( rays(k)% p   == invalid)   .or. &
            ( dlat         == invalid)   .or. &
            ( dlon         == invalid)   .or. &
            ( rays(k)% p   <= occ% rlc)  .or. &
            ( abs(dlat)    >   90._wp)   .or. &
            ( abs(dlon)    >  360._wp) )      &
            cycle
       if ( qimin > 0                   .and. &
            rays(k)% pcc /= pcc_invalid .and. &
            rays(k)% pcc <  qimin             ) cycle
       if ((modify_obs_err == 0)      .and. &
           (rays(k)% var   == 0._wp)) cycle           ! discard obs with var=0
       x = unitvector (dlat, dlon)
       if (sum (x * x_ref) < cut)     cycle
!      if (nray > 0) then
!         if (rays(k)% p == rays(ind(nray))% p) cycle ! duplicate tp heights
!      end if
       nray = nray + 1
       ind(nray) = k
    enddo

    !--------------------------------------
    ! exit if no valid data in occ profile
    !--------------------------------------
    if (nray == 0) then
      lkeep = .false.
      call decr_rpt_use (spot, CHK_INSDAT, comment='nray <= 0')
      return
    endif

    old  => rays
    allocate (rays (nray))
    occ% nray = nray
    rays = old (ind(1:nray))
    deallocate (old)
    nullify (old)

    !-------------------------------
    ! sort occ. profile if necessary
    ! (bottom -> top)
    !-------------------------------
    up = (rays(occ% nray)% p >= rays(1)% p)
    if (.not. up) then
       old => rays
       allocate (rays (nray))
       rays = old (occ% nray:1:-1)
       deallocate (old)
       nullify (old)
    endif

    !---------------------------------------------------------------------
    ! Beware: BUFR resolution (10cm) may be too limited for high-res data,
    ! so allow for "same" impact parameters
    !---------------------------------------------------------------------
    if (any (rays(2:occ% nray)% p < rays(1:occ% nray-1)% p)) then
       write(*,*) "check_store_occ: impact parameters not sorted: rays% p ="
       write(*,*) rays% p
!      call finish ("check_store_occ","impact parameters not sorted!")
       call message ("check_store_occ","impact parameters not sorted!")
       lkeep = .false.
       call decr_rpt_use (spot, CHK_CONSIST, comment='bad level ordering')
       return
    end if

    ndup = count (rays(2:occ% nray)% p == rays(1:occ% nray-1)% p)
    if (ndup > 0) then
       write(*,*) "check_store_occ: duplicate impact parameters for ", &
            spot% statid, ndup
    end if

    !----------------------------------------------------
    ! Check occultation geometry and metadata consistency
    !----------------------------------------------------
    if (operator >= 3) then
       call check_occ_geom (occ, spot, rays, status)
       select case (status)
       case default
          ! OK
       case (1)
          call decr_rpt_use (spot, CHK_CONSIST, comment='occ. geometry')
          lkeep = .false.
          return
       case (2)
          call decr_rpt_use (spot, CHK_INSDAT , comment='occ. geometry')
          lkeep = .false.
          return
       end select
    end if

    !----------------------------------------------
    ! Check if occultation starts below given bound
    !----------------------------------------------
    p_low = minval (rays(:)% p)
    if (occ_start_max > 0._wp) then
       z = p_low - occ% rlc
       if (z > occ_start_max) then
          call decr_rpt_use (spot, CHK_HEIGHT, STAT_PASSIVE, &
                                   comment='occultation_start')
       end if
    end if

    !-----------------------------------------------------
    ! average lowermost h_mean [km] of perigee coordinates
    !-----------------------------------------------------
    occ% gp1d  = occ% gp        ! default: 1d same georeferencing as 3d
    if (h_mean /= 0._wp) then
      !----------------------------------
      ! Derive average via zenith vectors
      !----------------------------------
      n_mean = 0
      x      = 0._wp
      do k = 1, occ% nray
        if (rays(k)% p  >  p_low + h_mean) cycle
        x      = x + unitvector (rays(k)% geo% Phi, rays(k)% geo% Lambda)
        n_mean = n_mean + 1
      enddo

      if (n_mean>0) then
         occ% gp1d%  Phi    = r2d * asin  (x(3) / sqrt (sum (x**2)))
         occ% gp1d%  Lambda = r2d * atan2 (x(2), x(1))
      else
        occ% gp1d  = occ% gp
      endif

    endif

    !--------------------------------------
    ! find perigee coordinates closest
    ! to height z_1d [km] of bangle profile
    !--------------------------------------
    if (z_1d /= 0._wp) then
       p_low=minval(rays(:)% p)
       p_1d = p_low + z_1d
       i = sum(minloc(abs(rays(:)% p - p_1d)))
       occ% gp1d%  phi    = rays(i)% geo% phi
       occ% gp1d%  lambda = rays(i)% geo% lambda
    endif

    !----------------------
    ! Set observation error
    !----------------------
    call set_occ_obserr (spot, occ, rays)

    !----------------------
    ! pick up rays
    ! default:
    !    < 5km : dh =  250m
    !    <10km : dh =  500m
    !    <20km : dh = 1000m
    !----------------------
    ind  = 0
    z    = rays(1)% p - occ% rlc        ! Level of first ray
    k    = maxloc (ray_levs, dim=1)
    i    = count  (ray_levs(1:k) <= z)
    i    = max (i, 1)
    h    = ray_levs(i)                  ! start at lowest level
    htop = ray_levs(k)                  ! Highest ray level
    nray = 0

#if defined (__GFORTRAN__)
    ! Barrier for optimizer to work around a gfortran 4.3 problem
    if (occ% nray < 0) then
       write (0,*) "ray_levs(1:4) =", ray_levs(1:4)
    end if
#endif

l1: do k = 1, occ% nray
      z = rays(k)% p - occ% rlc         ! loop over impact params in dataset
      if (z > htop) exit l1
      if (z >= h) then
        if (nray >= size(ind)) exit l1
        nray      = nray + 1            ! chose ray
        ind(nray) = k
        do
          if (z < ray_levs(i+1)) exit   ! chose new increment ?
          i=i+1
          if (i==size(ray_levs)) exit l1
        end do
        if (ray_incr(i) <= 0)    exit l1
        do while(z >= h)
          h = h + ray_incr(i)           ! chose new lower bound
        end do
      endif
    end do l1

    !-------------------
    ! keep selected rays
    !-------------------
    nold = nray
    i    = 0
    if (smooth_eps) then
       if (keep_lowest) then
          nray = nray - 1               ! Drop only highest ray of selected
          i    = 0
       else
          nray = nray - 2               ! Drop lowest and highest ray
          i    = 1
       end if
    endif
    if (nray <= 0) then
      lkeep = .false.
      call decr_rpt_use (spot, CHK_INSDAT, comment='nray <= 0')
      return
    endif

    old => rays
    allocate (rays (nray))
    occ% nray = nray
    rays = old (ind(1+i:occ% nray+i))

    !------------------------------------
    ! smooth observational bending angles
    !------------------------------------
    if (smooth_eps) then
      if (keep_lowest) then
        !-----------------------------------------
        ! Don't smooth bending angle of lowest ray
        !-----------------------------------------
        rays(1)% var2 = 0._wp
!       k=1
!       z = rays(k)% p - occ% rlc
!       write(0,'(a,i4,4ES24.15)') "k,z,eps,sqrt(var) =", &
!            k,z,rays(k)% eps,sqrt(rays(k)%var)
      end if
      do i = 2, nold-1
        ww   = 0._wp
        eps  = 0._wp
        eps2 = 0._wp
        var  = 0._wp
        x    = 0._wp
!       p    = 0._wp
        z    = 0._wp
        pcc  = 100
        k1 = ind (i-1)
        k2 = ind (i)
        k3 = ind (i+1)
        do k = k1+1, k2
          w    = (old(k)%p - old(k1)%p) / (old(k2)%p - old(k1)%p)
          ww   = ww   + w
          x    = x    + w * unitvector (old(k)% geo% Phi, old(k)% geo% Lambda)
!         p    = p    + w * old(k)% p
          z    = z    + w *(old(k)% p - occ% rlc) ! Omit offset in averaging
          eps  = eps  + w * log(old(k)% eps)
          eps2 = eps2 + w * log(old(k)% eps**2)
          var  = var  + w * log(old(k)% var)
          pcc  = min (pcc,  old(k)% pcc)
        end do
        do k = k2 + 1, k3-1
          w    = (old(k)%p - old(k3)%p) / (old(k2)%p - old(k3)%p)
          ww   = ww   + w
          x    = x    + w * unitvector (old(k)% geo% Phi, old(k)% geo% Lambda)
!         p    = p    + w * old(k)% p
          z    = z    + w *(old(k)% p - occ% rlc) ! Omit offset in averaging
          eps  = eps  + w * log(old(k)% eps)
          eps2 = eps2 + w * log(old(k)% eps**2)
          var  = var  + w * log(old(k)% var)
          pcc  = min (pcc,  old(k)% pcc)
        end do
        if (keep_lowest) then
           k = i
        else
           k = i-1
        end if
        rays(k)% geo% Phi    = r2d * asin  (x(3) / sqrt (sum (x**2)))
        rays(k)% geo% Lambda = r2d * atan2 (x(2), x(1))
!       rays(k)% p           =      p    / ww
        rays(k)% p           =      z    / ww + occ% rlc
        rays(k)% eps         = exp (eps  / ww)
        rays(k)% var         = exp (var  / ww)
        rays(k)% var2        = max (0._wp, exp(eps2/ww) - rays(k)% eps **2)
        rays(k)% pcc         =      pcc
!z = rays(k)% p - occ% rlc
!write(0,'(a,2i4,4ES24.15)') "i,k,z,eps,sqrt(var) =", &
!  i,k,z,rays(k)% eps,sqrt(rays(k)%var)
      end do
    else
      rays% var2 = 0._wp
    endif
    deallocate (old)
    nullify (old)

    !---------
    ! printout
    !---------
    if (verbose > 0) then

      call nl
      call nl; write (oline(iol),'(a,a)')     ' statid: ',spot% statid
      call nl; write (oline(iol),'(a,a)')     ' date  : ',cyyyymmdd(occ% time)
      call nl; write (oline(iol),'(a,a)')     ' time  : ',chhmmss  (occ% time)
      call nl; write (oline(iol),'(a,2f10.3)')' gp    : ',occ% gp%   phi,  &
                                                          occ% gp%   lambda
      call nl; write (oline(iol),'(a,2f10.3)')' gp1d  : ',occ% gp1d% phi,  &
                                                          occ% gp1d% lambda
      if (operator >= 3) then
       call nl;write (oline(iol),'(a,2f10.3)')' gp3d  : ',occ% gp3d% phi,  &
                                                          occ% gp3d% lambda
      end if
      if ( any (occ% rleo% x /= cart_inv% x) .or. &
           any (occ% vleo% x /= cart_inv% x) .or. &
           any (occ% rgps% x /= cart_inv% x) .or. &
           any (occ% vgps% x /= cart_inv% x) .or. &
           any (occ% xlc % x /= cart_inv% x) ) then
       call nl;write (oline(iol),'(a,4f10.3)')' rleo  : ',occ% rleo !, sqrt(sum(occ%rleo%x**2))
       call nl;write (oline(iol),'(a,4f10.3)')' vleo  : ',occ% vleo !, sqrt(sum(occ%vleo%x**2))
       call nl;write (oline(iol),'(a,4f10.3)')' rgps  : ',occ% rgps !, sqrt(sum(occ%rgps%x**2))
       call nl;write (oline(iol),'(a,4f10.3)')' vgps  : ',occ% vgps !, sqrt(sum(occ%vgps%x**2))
       call nl;write (oline(iol),'(a,3f10.3)')' xlc   : ',occ% xlc
      end if
      call nl; write (oline(iol),'(a,3f10.3)')' rlc   : ',occ% rlc
      if (occ% pcd% pcd /= pcd_invalid) then
        ! Clear dummy bits 3..5 out of bit 0..15:
        k = iand (occ% pcd% pcd, (2**16-1) - (8+16+32))
        call nl;write(oline(iol),'(a,b16.16)')' pcd   : ',k
      else
        call nl;write(oline(iol),'(a)')       ' pcd   : (invalid)'
      end if
      call nl
      call nl; write (oline(iol),'(a)') &
                      ' selected impact parameters and bending angles:'
      call nl; write (oline(iol),'(a)') &
 '   k         p       lat       lon         eps        deps       deps2  pcc'
      do k = 1, occ% nray
        call nl; write (oline(iol),'(i4,3f10.3,3f12.8,2x,i3)') k,       &
          rays(k)%p-occ%rlc, rays(k)% geo% Phi, rays(k)% geo% Lambda,   &
          rays(k)%eps, sqrt(rays(k)%var), sqrt(rays(k)%var2), rays(k)%pcc
      end do
      ! Write satellite positions only when set to non-fill value in BUFR
      if (any (rays(1)% rleo% x /= cart_inv% x) .or. &
          any (rays(1)% rgps% x /= cart_inv% x)      ) then
        call nl
        call nl; write (oline(iol),'(a,a)') &
                   '   k     rcurv       x_leo       y_leo',&
                   '       z_leo       x_gps       y_gps       z_gps'
        do k = 1, occ% nray
          if (k > 1 .and. all (rays(k)% rleo% x == rays(1)% rleo% x) &
                    .and. all (rays(k)% rgps% x == rays(1)% rgps% x)) cycle
          call nl; write (oline(iol),'(i4,f10.3,6f12.3)') &
               k, occ% rlc, rays(k)% rleo% x, rays(k)% rgps% x
        end do
      end if
      call nl
      call nl; write (oline(iol),'(a)') repeat('-',79)
    endif

    !====================================
    ! store data in observation data type
    !====================================
    if (occ_int_size == 0) call set_size
    !----------
    ! meta data
    !----------
    spot%         int_type    = ITY_MCOLS
!   spot%         ident       = occ% occid
    if (operator >= 3) then
     spot%col% c% dlat        = occ% gp3d% phi
     spot%col% c% dlon        = occ% gp3d% lambda
    else
     spot%col% c% dlat        = occ% gp1d% phi
     spot%col% c% dlon        = occ% gp1d% lambda
    end if
    spot% z                   = occ% rlc * 1000._wp
    spot% col%    nlev        = occ% nray
    spot%         cost        = occ% nray * 100._wp
    spot%         char        = CHR_NONL+CHR_EXP
    spot%         actual_time = occ% time
    call set_xuv  (spot)
    !--------------------------
    ! report selection (filter)
    !--------------------------
    call check_report_1 (spot)
    if (sat_ids(1) /= -1 .and. all(sat_ids /= spot% hd% satid)) &
      call decr_rpt_use (spot, CHK_BLACKLIST)

    if (spot% use% state > STAT_DISMISS) then
      call new_spot (obs,1, set_id=.true.)
      spt => obs% spot (obs% n_spot)
      id      = spt% id
      spt     = spot
      spt% id = id
      !------------------------------------------------------------------------
      ! Set satellite positions&velocities for all rays from report as fallback
      !------------------------------------------------------------------------
      if (all (rays(1)% rleo% x == cart_inv% x)) rays(:)% rleo    = occ% rleo
      if (all (rays(1)% vleo% x == cart_inv% x)) rays(:)% vleo    = occ% vleo
      if (all (rays(1)% rgps% x == cart_inv% x)) rays(:)% rgps    = occ% rgps
      if (all (rays(1)% vgps% x == cart_inv% x)) rays(:)% vgps    = occ% vgps
      where   (rays(:)% bearing == invalid     ) rays(:)% bearing = occ% bearing
      !-----------------
      ! observation body
      !-----------------
      bod % mn          = 'bendangl'
      bod % use % state = STAT_ACTIVE
      bod % use % check = CHK_NONE
      call new_obs (obs, spt% col% nlev, spot=spt)
      spt% nr                                                = spt% col% nlev
      obs % varno (spt% o% i+1 : spt% o% i + spt% o% n)      = VN_BENDANG
      obs %  olev (spt% o% i+1 : spt% o% i + spt% o% n)      = &
                                            1000._wp * (rays% p - occ% rlc)
      obs %  body (spt% o% i+1 : spt% o% i + spt% o% n)      = bod
      obs %  body (spt% o% i+1 : spt% o% i + spt% o% n)% o   = rays% eps
      obs %  body (spt% o% i+1 : spt% o% i + spt% o% n)% pcc = rays% pcc
      obs %  body (spt% o% i+1 : spt% o% i + spt% o% n)% eo  = sqrt (rays% var)
      if (operator == 1 .and. effective_geo) then
        ! 1d operator without tangent point drift using effective geolocation
        obs% body (spt% o% i+1 : spt% o% i + spt% o% n)% lat = occ% gp1d% phi
        obs% body (spt% o% i+1 : spt% o% i + spt% o% n)% lon = occ% gp1d% lambda
      else
        obs% body (spt% o% i+1 : spt% o% i + spt% o% n)% lat = rays% geo% phi
        obs% body (spt% o% i+1 : spt% o% i + spt% o% n)% lon = rays% geo% lambda
      end if
      obs %  body (spt% o% i+1 : spt% o% i + spt% o% n)% obs_par(1)    &
                                                             = rays% bearing
      !-------------------------------
      ! Check for GPSRO specific rules
      !-------------------------------
      do i = spt% o% i+1, spt% o% i + spt% o% n
        u = STAT_ACCEPTED
        call get_rule (type     =      spt% hd% modtype,    &! <- module      type
                       obstype  =      spt% hd% obstype,    &! <- observation type
                       codetype =      spt% hd% codetype,   &! <- code type
                       bf_type  =      spt% hd% buf_type,   &! <- BUFR    type
                       bf_subt  =      spt% hd% buf_subtype,&! <- BUFR subtype
                       db_kz    =      spt% hd% dbkz,       &! <- Datenbankkennzahl
                       center   =      spt% center_id,      &! <- Processing center
                       satid    = int (spt% hd% satid),     &! <- Satellite id.
                       stname   =      spt% statid,         &! <- Station Name
                       lat      =      spt% col% c% dlat,   &! <- Latitude
                       lon      =      spt% col% c% dlon,   &! <- Longitude
                       zobs     =      obs% olev (i),       &! <- z
                       phase    =      occ% pcd% PCD,       &! <- PCD flag
                       use      =      u)                    ! -> use flag
        u = min (u, int (obs% body(i)% use% state))
        call decr_use (obs% body(i)% use, u, check=CHK_RULE)
      end do
      !------------------------------------------------
      ! store information read in observation data type
      !------------------------------------------------
      call store_occ (obs, spt, occ, rays)
    else
      lkeep = .false.
    endif

  end subroutine check_store_occ
!==============================================================================
  subroutine set_occ_obserr (spot, occ, rays)
    !------------------------------------------------
    ! set observation error variance for GNSS RO data
    !------------------------------------------------
    type(t_spot),intent(in)  :: spot
    type(t_occ) ,intent(in)  :: occ
    type(t_ray) ,pointer     :: rays (:)
    !----------------
    ! Local variables
    !----------------
    integer  :: k, i, j       ! indices
    real(wp) :: z             ! height [km]
    real(wp) :: r             ! relative error
    real(wp) :: w, w2, ws, wt ! temporary weights
    real(wp) :: new_var       ! observation error variance
    real(wp) :: obs2          ! observed bending angle ** 2
    real(wp) :: ht            ! transition height (troposphere <-> stratosphere)
    real(wp) :: lat, lat_     ! temporary: latitude, restricted latitude
    real(wp) :: coslat        ! temporary: cos(latitude)
    real(wp) :: h_cpt         ! temporary: effective CP tropopause height
    real(wp) :: dh_cpt        ! temporary: effective CP tropopause width
    real(wp) :: s_cpt         ! temporary: effective TP representativeness error
    !---------------------------------------------
    ! Tuning parameters of observation error model
    !---------------------------------------------
    integer  :: order = 4    ! shape parameter of latitudinal dependence
    type(t_occ_obserr), pointer :: e

    real(wp), parameter :: ln2 = log (2._wp)

    if (modify_obs_err >= 7) then
       do j = nerr, 1, -1
          e => occ_obserr(j)
          if (e% satid(1) > 0 .and. all (e% satid  /= spot% hd% satid)) cycle
          if (e% center   > 0 .and.      e% center /= spot% hd% center) cycle
          if (e% gnssid   > 0 .and.      e% gnssid /= occ% gnssid     ) cycle
          exit
       end do
       if (j < 1) call finish ("set_occ_obserr",                           &
                               "no obserr specification for "//spot% statid)
    else
       e => NULL()
    end if

    if (modify_obs_err > 0) then
       do k=1, occ% nray
          z = rays(k)% p - occ% rlc
          obs2 = rays(k)% eps ** 2
          if (modify_obs_err < 7) then
             !------------------------------------
             ! modify observation errors (S.Healy)
             !------------------------------------
!!$          if (z > 10._wp)  then
!!$             new_var = (max(rays(k)% eps, 0._wp) * obs_err_10km) ** 2
!!$          else
!!$             new_var = (max(rays(k)% eps, 0._wp) *                     &
!!$                  (obs_err_10km + (10._wp - z) * &
!!$                  (obs_err_0km - obs_err_10km) / 10._wp)) ** 2
!!$          endif
             !-----------------------------------------------
             ! Piecewise linear parameterization of obs.error
             !-----------------------------------------------
             do i = 1, obs_err_nlev
                if (z < obs_err_hgt(i)) exit
             end do
             w = (obs_err_hgt(i) - z) / (obs_err_hgt(i) - obs_err_hgt(i-1))
             r = w * obs_err_val(i-1) + (1._wp - w) * obs_err_val(i)
             new_var = max (0._wp,   obs2 * r ** 2)
             new_var = max (new_var, obs_err_abs ** 2)
          else
             lat     = rays(k)% geo% phi
             coslat  = abs (cos (lat * d2r))
             w       = coslat ** order
             w2      = coslat ** 2
             ws      = exp (-ln2*(lat / e% dl_cpt_t)**2)
!            ws      = 1 / (1 + (lat / e% dl_cpt_t)**2)
             select case (e% cpt_h_shape)
             case (1)
               h_cpt  = e% h_cpt_p      * (1-w ) + e% h_cpt_t      * w
               dh_cpt = e% dh_cpt_p     * (1-w ) + e% dh_cpt_t     * w
             case (2)
               h_cpt  = e% h_cpt_p      * (1-w2) + e% h_cpt_t      * w2
               dh_cpt = e% dh_cpt_p     * (1-w2) + e% dh_cpt_t     * w2
             case default ! 3
               lat_   = max (min (abs (lat), 60._wp), 30._wp)
               wt     = abs (cos ((lat_ - 30._wp) * 3 * d2r)) ** 2
               h_cpt  = e% h_cpt_p      * (1-wt) + e% h_cpt_t      * wt
               dh_cpt = e% dh_cpt_p     * (1-wt) + e% dh_cpt_t     * wt
             end select
             s_cpt   =  e% s_cpt_p      * (1-ws) + e% s_cpt_t      * ws
             ht      =  e% h0_pol       * (1-w ) + e% h0_tr        * w
             z       =  max (z, epsilon (z))
             r       = (e% s_strat_p**2 * (1-w ) + e% s_strat_t**2 * w)        &
                     + (e% s0_pol   **2 * (1-w ) + e% s0_tr    **2 * w)        &
                       * max (0._wp, 1._wp - z / ht) ** 2                      &
                     + s_cpt**2 * exp (-ln2*(log (z/h_cpt) / (dh_cpt/h_cpt))**2)
             new_var = (e% s_iono_p**2 * (1-w) + e% s_iono_t**2 * w) + r * obs2
          end if
          select case (modify_obs_err)
          case (1,4)
             if (rays(k)% var == 0._wp ) rays(k)% var = new_var
          case (2,5)
             if (rays(k)% var < new_var) rays(k)% var = new_var
          case (3,6:)
             rays(k)% var = new_var
          end select
!write(0,'(a,i4,3ES24.15)') "k, z, eps, var =", k,z,rays(k)% eps,rays(k)%var
       enddo
    endif
  end subroutine set_occ_obserr
!==============================================================================
  subroutine check_mono_occ (obs)
  !-----------------------------------------------------------------
  ! check occultation observation profile for consistency (monotony)
  ! check             first guess profile as well
  !-----------------------------------------------------------------
    type (t_obs_set) ,intent(inout) :: obs ! observation data

    integer               :: ib, is, i ! indices
    type(t_spot) ,pointer :: si        ! pointer to report meta data
    type(t_obs)  ,pointer :: oi        ! pointer to observation box
    logical               :: crop_o
    logical               :: crop_fg
    real(sp)              :: ocut_o    ! effective cut for bending angle (Obs)
    real(wp)              :: ocut_fg   ! effective cut for bending angle (FG)
    real(sp)              :: tol       ! tolerance for non-monotony
    integer               :: n_o, n_fg ! count of removed rays for reason
    logical               :: l_o, l_fg ! check obs, fg
    logical               :: p_o, p_fg ! print obs, fg
    integer               :: n_oo      ! count of rays with status OBS_ONLY

    l_o  = iand (chk_mono, 1) /= 0
    l_fg = iand (chk_mono, 2) /= 0
    p_o  = iand (chk_mono, 4) /= 0
    p_fg = iand (chk_mono, 8) /= 0
    l_o  = l_o  .or. p_o
    l_fg = l_fg .or. p_fg
    if (.not. (l_o .or. l_fg .or. rm_obs_only)) return
    n_o  = 0
    n_fg = 0
    n_oo = 0
    do ib = 1, size(obs% o)
      if (obs% o(ib)% pe /= dace% pe) cycle
      oi => obs% o(ib)
      do is = 1, oi% n_spot
        si      => oi% spot(is)
        if (si% hd% obstype /= OT_GPSRO)     cycle
        if (si% use% state  <= STAT_DISMISS) cycle
        crop_o  = .false.
        crop_fg = .false.
        ocut_o  = 0._sp
        ocut_fg = 0._wp

        if (rm_obs_only) then
           do i = si%o%i+1, si%o%i+si%o%n
              if (oi% body(i)% use% state == STAT_OBS_ONLY) then
                 call decr_use (oi% body(i)% use, STAT_DISMISS, check=CHK_FG)
                 n_oo = n_oo + 1
              end if
           end do
        end if

        do i =si%o%i+si%o%n, si%o%i+2, -1
          if (oi% varno(i) == VN_BENDANG) then

            tol     = sgm_mono * oi% body(i)% eo
            ocut_o  = max (ocut_o,  oi% body(i)% o         - tol)
            ocut_fg = max (ocut_fg, obs%b% yb% s(ib)% x(i) - tol)

            if   (oi% olev(i-1) <= level_mono) then
              if (l_o) then
                if (oi% body(i-1)% o < ocut_o .or. crop_o) then
                  call decr_use (oi% body(i-1)% use, check=CHK_FG)
                  crop_o = .true.
                  n_o    = n_o + 1
                  if (p_o) write(*,*) "check_mono_occ: ", si% statid, &
                       ": o @ z =", real (oi% olev(i-1))
                endif
              endif
              if (l_fg .and. oi% body(i-1)% use% state > STAT_OBS_ONLY) then
                if (obs%b% yb% s(ib)% x(i-1) < ocut_fg .or. crop_fg) then
                  call decr_use (oi% body(i-1)% use, check=CHK_FG)
                  crop_fg = .true.
                  n_fg    = n_fg + 1
                  if (p_fg) write(*,*) "check_mono_occ: ", si% statid, &
                       ": fg @ z =", real (oi% olev(i-1))
                endif
              endif
            endif

          endif
        end do
      end do
    end do

    n_fg = p_sum (n_fg)
    n_o  = p_sum (n_o)
    n_oo = p_sum (n_oo)

    if (dace% lpio .and. n_fg + n_o + n_oo > 0) then
      write(6,'(/,a,/)') repeat('-',79)
      if (n_fg + n_o > 0) then
        write(6,'(a,/)')   "  Occultations with strong departure from monotony"
        write(6,'(a,i7)')  "  first-guess: n_fg =", n_fg
        write(6,'(a,i7)')  "  observation: n_o  =", n_o
        write(6,'()')
      end if
      if (n_oo > 0) then
        write(6,'(a,/)')   "  Rays removed (OBS_ONLY)"
        write(6,'(a,i7)')  "  observation: n_oo =", n_oo
        write(6,'()')
      end if
    end if

  end subroutine check_mono_occ
!==============================================================================
  subroutine check_susp_occ (obs)
    type (t_obs_set) ,intent(inout) :: obs ! observation data
    !------------------------------------------------------
    ! Check for suspicious occultation observation profiles
    ! after first guess check has been applied
    !------------------------------------------------------
    integer               :: ib, is, i ! indices
    type(t_spot) ,pointer :: si        ! pointer to report meta data
    type(t_obs)  ,pointer :: oi        ! pointer to observation box
    integer               :: nacc      ! active/passive ray count
    integer               :: nrej      ! fg rejected rays
    integer               :: np        ! # processed  profiles
    integer               :: ns        ! # suspicious profiles
    integer               :: na        ! # accepted   profiles
    integer               :: ng        ! # profiles with gaps
    integer               :: nb        ! # fg rejected profiles
    real(sp)              :: tol       ! tolerance for rejected rays
    type(t_occ)           :: o         ! occultation data type
    type(t_ray)  ,pointer :: rs(:)     ! rays data type
    integer               :: ngap      ! # gaps in profile
    real(wp)              :: delta     ! impact parameter difference
    real(wp)              :: gapmax    ! max. gap size
    real(wp)              :: hgap      ! impact height of largest gap
    logical               :: lrej      ! reject profile?

    if (chk_susp    <= 0) return
    if (thresh_susp <= 0) return
    tol = min (thresh_susp, 1.0)

    if (dace% lpio) then
      write(6,'(a)') repeat('-',79)
      write(6,'()')
      write(6,'(a)') "  Check for suspicious occultation profiles"
      write(6,'()')
    end if

    np = 0
    na = 0
    ns = 0
    ng = 0
    nb = 0
    do ib = 1, size(obs% o)
      if (obs% o(ib)% pe /= dace% pe) cycle
      oi => obs% o(ib)
      !------------------------
      ! loop over GPSRO reports
      !------------------------
      do is = 1, oi% n_spot
        if  (oi% spot(is)% hd% obstype /= OT_GPSRO)     cycle
        if  (oi% spot(is)% use% state  <= STAT_DISMISS) cycle
        si   => oi% spot(is)
        np   = np + 1
        nacc = 0
        nrej = 0
        ngap = 0
        lrej = .false.
        !----------------------------------------
        ! Detect data gaps (for information only)
        !----------------------------------------
        if (gap_tol > 0._wp) then
           gapmax = 0._wp
           hgap   = 0._wp
           call load_occ (oi, si, o, rs)
           do i = 1, size (rs)-1
              delta = rs(i+1)% p - rs(i)% p
              if (delta > gap_tol) then
                 ngap = ngap + 1
                 if (delta > gapmax) then
                    gapmax = delta
                    hgap   = rs(i)% p - o% rlc
                 end if
              end if
           end do
           if (ngap > 0) then
              call nl
              write (oline(iol),'(2x,a,a,2f10.3,3x,a,i6,2f6.1)') &
                   'occ lat lon: ', si% statid,                  &
                   si% col% c% dlat, si% col% c% dlon,           &
                   'ngap max  z:', ngap, gapmax, hgap
              ng = ng + 1
              if (ngap > max_gaps) lrej = .true.
           end if
           deallocate (rs)
        end if
        !-------------------------------
        ! Count rejected rays in profile
        !-------------------------------
        do i = si%o%i+1, si%o%i+si%o%n
           if (oi% varno(i) == VN_BENDANG) then
              select case (oi% body(i) % use% state)
              case (STAT_PASSIVE,STAT_ACTIVE_0I:STAT_ACCEPTED)
                 nacc = nacc + 1
              case (STAT_PAS_REJ,STAT_REJECTED,:STAT_OBS_ONLY)
                 nrej = nrej + 1
              end select
           end if
        end do
        if (nrej > 0 .and. nrej >= (nacc+nrej)*tol) then
           call nl
           write (oline(iol),'(2x,a,a,2f10.3,3x,a,2i6,f6.1)') &
                'occ lat lon: ', si% statid,                  &
                si% col% c% dlat, si% col% c% dlon,           &
                'acc. rej. %:', nacc, nrej, 100.*nrej/(nacc+nrej)
           nb   = nb + 1
           lrej = .true.
        end if
        if (lrej) then
           ns = ns + 1
           if (chk_susp >= 2) &
                call decr_rpt_use (si, CHK_CONSIST, STAT_REJECTED, &
                                   comment="suspicious")
        else
           na = na + 1
        end if
      end do ! is
    end do   ! ib

    np = p_sum (np)
    na = p_sum (na)
    ns = p_sum (ns)
    ng = p_sum (ng)
    nb = p_sum (nb)

    call flush_buf ()

    if (dace% lpio) then
      write(6,'()')
      write(6,'(a,i7)')     '  processed    profiles = ',np
      write(6,'(a,i7)')     '  accepted     profiles = ',na
      write(6,'(a,i7)')     '  suspicious   profiles = ',ns
      write(6,'(a,i7)')     '  gap-affected profiles = ',ng
      write(6,'(a,i7)')     '  fg rejected  profiles = ',nb
      write(6,'()')
    end if

  end subroutine check_susp_occ
!==============================================================================
  subroutine read_quality(occ)

    type(t_occ)  ,intent(inout)  :: occ

    !-------------------------------------------------------------------------------------------------------
    !                       !                                                                              !
    !    product quality    !                 value of product quality flag                                !
    !          flag         !                                                                              !
    !                       !         value 0                    !              value 1                    !
    !                       !                                    !                                         !
    !-------------------------------------------------------------------------------------------------------
    !    pcd_quality        !  nominal quality                   !   non-nominal quality                   !
    !    pcd_product        !  NRT product                       !   offline product                       !
    !    pcd_occult_type    !  setting occultation               !   rising occultation                    !
    !    pcd_exsphs_proc    !  excess phase processing nominal   !   excess phase processing non-nominal   !
    !    pcd_bangle_procc   !  bending angle processing nominal  !   bending angle processing non-nominal  !
    !    pcd_refrac_procc   !  refractivity  processing nominal  !   refractivity  processing non-nominal  !
    !    pcd_meteor_procc   !  meteorological processing nominal !   meteorological processing non-nominal !
    !    pcd_open_loop      !  open loop tracking not used       !   open loop tracking                    !
    !    pcd_reflection     !  no surface reflections detected   !   surface reflections detected          !
    !    pcd_l2_signal      !  L2P GPS signal used               !   L2C GPS signal used                   !
    !    pcd_bckgrnd_prof   !  background profile nominal        !   background profile non-nominal        !
    !    pcd_retriev_prof   !  retrieved profile                 !   background profile                    !
    !    pcd_missing        !  pcd bits 1-15 valid               !   pcd missing: bits 1-15 invalid        !
    !-------------------------------------------------------------------------------------------------------
    !    See also GRAS SAF ROPP User Guide, IO module, Version 6.0, Table 2.5                              !
    !-------------------------------------------------------------------------------------------------------

    !------------------
    ! set quality flags
    !------------------
    occ% pcd% quality      = 0
    occ% pcd% product      = 0
    occ% pcd% occult_type  = 0
    occ% pcd% exsphs_proc  = 0
    occ% pcd% bangle_procc = 0
    occ% pcd% refrac_procc = 0
    occ% pcd% meteor_procc = 0
    occ% pcd% open_loop    = 0
    occ% pcd% reflection   = 0
    occ% pcd% l2_signal    = 0
    occ% pcd% dummy4       = 0
    occ% pcd% dummy5       = 0
    occ% pcd% dummy6       = 0
    occ% pcd% bckgrnd_prof = 0
    occ% pcd% retriev_prof = 0
    occ% pcd% missing      = 0

    if (btest(occ% pcd% PCD,15)) occ% pcd% quality      = 1
    if (btest(occ% pcd% PCD,14)) occ% pcd% product      = 1
    if (btest(occ% pcd% PCD,13)) occ% pcd% occult_type  = 1
    if (btest(occ% pcd% PCD,12)) occ% pcd% exsphs_proc  = 1
    if (btest(occ% pcd% PCD,11)) occ% pcd% bangle_procc = 1
    if (btest(occ% pcd% PCD,10)) occ% pcd% refrac_procc = 1
    if (btest(occ% pcd% PCD, 9)) occ% pcd% meteor_procc = 1
    if (btest(occ% pcd% PCD, 8)) occ% pcd% open_loop    = 1
    if (btest(occ% pcd% PCD, 7)) occ% pcd% reflection   = 1
    if (btest(occ% pcd% PCD, 6)) occ% pcd% l2_signal    = 1
    if (btest(occ% pcd% PCD, 5)) occ% pcd% dummy4       = 1
    if (btest(occ% pcd% PCD, 4)) occ% pcd% dummy5       = 1
    if (btest(occ% pcd% PCD, 3)) occ% pcd% dummy6       = 1
    if (btest(occ% pcd% PCD, 2)) occ% pcd% bckgrnd_prof = 1
    if (btest(occ% pcd% PCD, 1)) occ% pcd% retriev_prof = 1
    if (btest(occ% pcd% PCD, 0)) occ% pcd% missing      = 1

  end subroutine read_quality

!==============================================================================
  subroutine read_fdbk_occ (o)
  type (t_obs) ,intent(inout) :: o          ! observation data type variable
  !-----------------------------------------------
  ! restore t_occ, t_ray from t_spot, body
  ! after observation was read from feedback file
  !-----------------------------------------------
    type (t_spot) ,pointer     :: s      ! pointer to report
    type (t_occ)               :: occ    ! GPSRO specific type to restore
    type (t_ray)  ,allocatable :: r(:)   ! GPSRO specific type to restore
    integer                    :: n      ! number of rays
    integer                    :: na     ! number of rays allocated
    integer                    :: i      ! report

    na = 100
    allocate (r (na))
    !------------------
    ! loop over reports
    !------------------
    do i = 1, o% n_spot
      if  (o% spot(i)% hd% obstype /= OT_GPSRO) cycle
      s => o% spot(i)
      n =  s% o% n
      !-------------------------
      ! GPSRO specific meta data
      !-------------------------
      s% int_type     = ITY_MCOLS
      s% char         = CHR_NONL+CHR_EXP
      s% nr           = n                       ! size of R matrix block
      !----------------------
      ! re-allocate temporary
      !----------------------
      if (n > na) then
        na = n
        deallocate (r)
        allocate   (r(na))
      endif
      !--------------------------
      ! restore t_occ from t_spot
      !--------------------------
      occ% rlc        = s% z / 1000._wp
      occ% nray       = s% col%    nlev
      occ% time       = s% actual_time
      occ% gp% phi    = s% col% c% dlat
      occ% gp% lambda = s% col% c% dlon
      occ% gp1d       = occ% gp
      occ% gp3d       = occ% gp
      occ% pcd% pcd   = s% phase
      occ% gnssid     = s% tracking
      occ% prn        = s% sender_id
      !-------------------------------------------------------------
      ! for MEC: observations from feedback file should be 'checked'
      !-------------------------------------------------------------
      occ% checked    = .true.
      !------------------------
      ! restore t_ray from body
      !------------------------
      r(:n)% p           = (o% olev (s%o%i+1 : s%o%i+n) + s% z) / 1000._wp
      r(:n)% eps         =  o% body (s%o%i+1 : s%o%i+n)% o
      r(:n)% pcc         =  o% body (s%o%i+1 : s%o%i+n)% pcc
      r(:n)% var         =  o% body (s%o%i+1 : s%o%i+n)% eo ** 2
      r(:n)% var2        =  0._wp
      r(:n)% geo% phi    =  o% body (s%o%i+1 : s%o%i+n)% lat
      r(:n)% geo% lambda =  o% body (s%o%i+1 : s%o%i+n)% lon
      r(:n)% bearing     =  o% body (s%o%i+1 : s%o%i+n)% obs_par(1)
      !------------------------------
      ! store in component t_obs% par
      !------------------------------
      call store_occ (o, s, occ, r(1:n))
    end do

  end subroutine read_fdbk_occ
!======================================================================
end module mo_occ
