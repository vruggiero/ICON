!
!+ Atmospheric state derived type declarations and operator definitions
!
MODULE mo_atm_state
!
! Description:
!   Defines data type 't_atm' holding the variables of an atmospheric state
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
! V1_2         2008/12/04 Andreas Rhodin
!  add prognostic variables: rain, snow, pressure deviation
! V1_4         2009/03/26 Andreas Rhodin
!  Changes for COSMO LETKF
! V1_5         2009/05/25 Andreas Rhodin
!  read runtype,runclass,expid from GRIB-file
! V1_6         2009/06/10 Andreas Rhodin
!  fix vertical interpolation (half level -> any level)
! V1_7         2009/08/24 Andreas Rhodin
!  add fields to pass through for COSMO
! V1_8         2009/12/09 Andreas Rhodin
!  call print(grib) from parallel section
!  account for flag geof=T in case of COSMO vertical coordinate
! V1_9         2010/04/20 Hendrik Reich
!  select_params0: pass u_c, v_c (for COSMO C-grid)
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_11        2010/06/16 Andreas Rhodin
!  account for prognostic rain,snow,graupel in height calculations
! V1_13        2011/11/01 Andreas Rhodin
!  hold ensemble_id in meta data
! V1_19        2012-04-16 Harald Anlauf
!  optimize uvrot2uv_vec, uv2rotuv_vec
!  vectorise ensemble mean and spread calculation for LETKF
!  prepare fields to pass through for GME LETKF: h_ice,t_ice,h_snow,tqr,tqs
! V1_20        2012-06-18 Andreas Rhodin
!  move routine derive_params from mo_letkf to mo_atm_state
!  new routine set_p1, atm_over_i, new option 'pp' in derive_params
! V1_22        2013-02-13 Andreas Rhodin
!  changes for handling of soil/surface variables in LETKF
!  changes for ICON vertical grid
! V1_23        2013-03-26 Andreas Rhodin
!  add components for LETKF: f_inflat, p_energy
! V1_26        2013/06/27 Jaison Ambadan
!  subroutine (ens_corr) for calculating correlations
!  new subroutine from_grads: read atmospheric state from GRADS file
! V1_27        2013-11-08 Harald Anlauf
!  Fixes for generalized vertical coordinate
! V1_28        2014/02/26 Andreas Rhodin
!  vector version for set_geo, set_pf, set_tv, new: set_rh
!  add infraread surface emissivity principal component coefficients
!  fixes for the case that nproc1/2 exceeds the respective # of grid-points
! V1_29        2014/04/02 Andreas Rhodin
!  new specific function: integer * atmospheric_state
! V1_31        2014-08-21 Andreas Rhodin
!  vector versions for set_ps, set_ph, C_to_A_grid, A_to_C_grid
!  LETKF: interpolate analysis increment to deterministic resolution
! V1_35        2014-11-07 Andreas Rhodin
!  new vector versions for:p_sum_atm, assign, destruct
!  new function levels (select subset of levels for isobaric vertical coord.)
!  add 'pmsl' (mean sea level pressure) to atmospheric state vector
! V1_37        2014-12-23 Andreas Rhodin
!  set_p for ICON: make sure ps if present before calling set_p
! V1_42        2015-06-08 Harald Anlauf
!  OpenMP parallelization
! V1_43        2015-08-19 Andreas Rhodin
!  add variables for 'flake' to derived type 't_atm'
! V1_44        2015-09-30 Andreas Rhodin
!  adapt interpolate_soil_flake to tail approach
! V1_45        2015-12-15 Andreas Rhodin
!  LETKF GRIB2 compatibility; revise soil/surface interpolation
!  TR15581/F2003 compatibility
! V1_46        2016-02-05 Michael Bender
!  New routine set_tqv; t_atm: new component 'alb_dif'
! V1_47        2016-06-06 Andreas Rhodin
!  add variables aswdir_s, aswdifd_s, alwd_s, qh,nccloud,...,nchail, dpsdt
!  include 'ps' (call set_ps) in 'derive_params'
!  set_rh2m: derive 2m relative humidity from dewpoint
! V1_48        2016-10-06 Harald Anlauf
!  t_atm: ICON, vorticity lives on dual grid
!  set_ps, set_ph: fixes for COSMO
! V1_50        2017-01-09 Andreas Rhodin
!  derive_params: new options 'geof', 'tv', 'den'
! V1_51        2017-02-24 Andreas Rhodin
!  add fields sso_stdh,sso_gamma,sso_theta,sso_sigma to atmospheric state
!  new subroutines p_min_atm,p_max_atm
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Andreas Rhodin  MPIfM      1998-2002
! Andreas Rhodin  DWD        2003-2008
! Harald Anlauf   DWD        2007-2008
! Le Duc,         Uni Hanoi  2005      Arakawa C-grid handling in vert.interp.
! J.-O.Beismann   NEC        2008      Vectorisation of vintp_spline
! Michael Bender  DWD        2016
!==============================================================================

#include "tr15581.incf"

  !-------------
  ! used modules
  !-------------
  use mo_exception,     only : finish               ! abort routine
  use mo_kind,          only : wp, sp, i8           ! working precision
  use mo_fortran_units, only : get_unit_number,    &! get free Fortran unit
                               return_unit_number   ! release  Fortran unit
  use mo_atm_grid,      only : t_grid, print,      &! grid information
                               VCT_P_ISO,          &! isobaric coordinate
                               VCT_P_HYB,          &! hybrid p coordinate
                               VCT_P_IFS,          &! IFS hybrid p coordinate
                               VCT_Z_HYB,          &! hybrid z coordinate
                               VCT_Z_GEN,          &! generalised z coordinate
                               MO_COSMO,           &! COSMO model
                               MO_ICON,            &! ICON  model
                               MO_GME,             &! GME   model
                               MO_HRM,             &! HRM   model
                               MO_IFS               ! IFS  model
  use mo_time,          only : t_time,             &! time derived type
                               print,              &! print time information
                               operator (+),       &! operations on ..
                               operator (*),       &! .. type t_time
                               operator (-),       &! ..
                               hours                ! extract hours from t_time
  use mo_dace_string,   only : split                ! string -> array
  use mo_memory,        only : t_m,                &! table of atm.fields
                               destruct,           &! memory management
                               allocate,           &!
                               update,             &!
                               deallocate,         &!
                               assignment (=),     &! table = real
                               assign_m_reals,     &! PRELIMINARY
                               plus, minus,        &!
                               times, over,        &!
                               power, sum, sqrt_m, &!
                               max_m,              &!
                               lgtd,               &! Legendre transf.
                               lgti,               &! inverse ''
                               lgtd_ad, lgti_ad,   &! adjoint ''
                               operator (/=),      &!
                               print, dump,        &!
                               retrieve,           &! i/o routines
                               nwords,             &! allocated reals
                               cut,                &! cut regions
                               p_min_m,            &! parallel min
                               p_max_m,            &! parallel max
                               p_sum_m              ! parallel sum
  use spherepack95,     only : sp95_init,          &!
                               sp95_cleanup,       &!
                               sp95_sfvp,          &! u,v->psi,chi
                               sp95_sfvp_gg,       &! " Gauss grid
                               sp95_sfvp_vrtdiv_gg  ! ... from vrt,div
  use mo_dft,           only : poisson_solve_dct    ! Poisson solver (DCT)
  use mo_math_oper,     only : hor_curl,           &! Curl of hor. vector field
                               hor_div              ! Divergence of  "    field
  use mo_physics,       only : pi,                 &! 3,1415...
                               q_rh,               &! spec.hum.from rh
                               rh_q,               &! rel.hum. from  q
                               tv_t_q,             &! virt.temp. <- t,q
                               t_tv_q,             &! t from tv, q
!                              tq_tvrh,            &! t,q from tv, rh
                               tq_tvrh_noadj,      &!+++ workaround
                               tq_tvrh_noxfail,    &!+++ workaround
                               t_tq_tvrh,          &! debug info
                               p_rho_thetav,       &! p <- rho,theta_v
                               exner,              &! Exner function
                               esw_t,              &! sat.vapor pres. / water
                               R,                  &! gas constant
                               RD,                 &! gas constant of water vapor
                               RDRD,               &! R/RD
                               EMRDRD,             &! 1 - RDRD
                               rddrm1,             &! r_v/r_d - 1
                               r2d, d2r,           &! degree<->radians
                               lapse_cl,           &! Clim. lapse rate
                               gacc,               &! gravity accel.
                               geopot_geoid         ! geopotential->height
  use mo_grads,         only : t_ctl                ! GRADS metadata
  use mo_grads_atm,     only : to_grads,           &! write state to GRADS file
                               from_grads           ! read      from GRADS file
  use mo_mpi_dace,      only : dace,               &! MPI communication info
                               p_isend,            &! unblocked _ISEND
                               p_recv,             &!        MPI_RECV
                               p_wait,             &!        MPI_WAIT
                               p_sum                !
  use mo_atm_transp,    only : scatter_level,      &!
                               gather_level
  use mo_algorithms,    only : init_spline2,       &! 2nd derivatives for ...
                               vintp_spline,       &! ... spline interpolation
                               integpolycube        ! numerical integration
  use mo_wmo_tables,    only : WMO6_GAUSSIAN,      &!
                               WMO6_LATLON,        &!
                               WMO6_ROTLL,         &!
                               DWD6_ICOSAHEDRON,   &!
                               DWD6_ICON            !
  use mo_grid_intpol,   only : hor_intp,           &! hor.interpolation
                               hor_intp_coeff,     &! coef.for hor_intp
                               t_h_idx              ! type ind.,weights
  use meteo_utilities,  only : calps                ! compute surf.pres
  use mo_soil,          only : n_st,               &! number of soil types
                               ind_wso,            &! index from soil moisture
                               wso_ind,            &! soil moisture from index
                               ST_ICE, ST_ROCK, ST_SEAWATER, ST_SEAICE
#ifdef __ICON__
  use sfc_flake,        only : flake_init           ! Flake consistency check
#else
  use mo_flake,         only : flake_init           ! Flake consistency check
#endif
#ifndef TR15581
  use mo_allocate,      only : call_level           ! calllevel counter
#endif
  implicit none
  !----------------
  ! public entities
  !----------------
  private
  !-------------------------------------------------
  ! derived type definition, constructor, destructor
  !-------------------------------------------------
  public :: t_atm              ! atmospheric state data-type
  public :: construct, destruct! initialize, free atmospheric state
  public :: allocate           ! allocate specific components of t_atm
  public :: deallocate         ! deallocate speci. components of t_atm
  public :: print              ! formatted printout of atmospheric state
  public :: nm                 ! max. number of fields in t_atm
  !-------------------------------------------
  ! enumerator positions of fields in in t_atm
  !-------------------------------------------
  public :: T,U,V,Q,T_SO,W_SO,RH, &!
            T2M,VMAX_10M,TMIN_2M, &!
            TMAX_2M,TOT_PREC,RH2M,&!
            T2M_LAND,RH2M_LAND,   &!
            ASWDIR_S,ASWDIFD_S     !
  !-------------------------------------
  ! assignment and arithmetic operations
  !-------------------------------------
  public :: assignment (=)     ! assign atmospheric state
  public :: operator (+)
  public :: operator (-)
  public :: operator (*)
  public :: operator (/)
  public :: operator (**)
  public :: sqrt               ! square root of atmospheric fields
  public :: max                ! maximum value
  public :: sum                ! sum of field elements
  public :: lgtd, lgti         ! Legendre transformations: direct, inverse
  public :: lgtd_ad, lgti_ad   ! Legendre transformations: adjoints
  !--------------------------------
  ! set certain components of t_atm
  !--------------------------------
  public :: set_p, set_p_adj   ! set pressure levels
  public :: set_pf, set_ph     ! set pressure levels (full, half)
  public :: set_ps             ! set surface pressure (COSMO)
  public :: set_geo            ! set geopotential
  public :: set_t              ! set temperature from virt.pot.temperature
  public :: set_tv             ! set virtual temperature
  public :: set_tv_geo         ! set virtual temperature from geopotential
  public :: set_tq_tv          ! set temperature from virtual temp. + moisture
  public :: set_q              ! set specific humidity from rh, t, p
  public :: set_rh             ! set relative humidity from q, t, p
  public :: set_rh2m           ! set 2m relative humidity from t2m, td2m
  public :: set_rh2m_land      ! set 2m relative humidity from t2m, td2m (land)
  public :: set_psi_chi        ! set streamfunction and velocity potential
  public :: set_psi_chi_u_v    ! set streamfunction, velocity potential, u,v
  public :: set_vrt_div        ! set vorticity and divergence from u,v
  public :: set_index          ! set soil moisture index from soil moisture
  public :: set_wso            ! set soil moisture from soil moisture index
  public :: mean_stdev         ! calculate mean & std.deviation of fields
  public :: set_tqv            ! set integrated water vapour
  public :: set_pointers       ! update pointers with memory table
  public :: ens_corr           ! calculate the ensemble correlation of fields
  !----------------------------------
  ! maintain lists of fields in t_atm
  !----------------------------------
  public :: select_params      ! gather fields from different sources (states)
  public :: print_sel_par      ! print out parameters selected by select_params
  public :: derive_params      ! transform fields (eg. q <-> rh etc.)
  public :: keep_params        ! keep only parameters given in list
  public :: merge              ! merge states
  !-----------------------------------------
  ! change representation of wind components
  !-----------------------------------------
  public :: A_to_C_grid        ! transform u,v,u_10m,v_10m from A-grid to C-grid
  public :: C_to_A_grid        ! transform u,v,u_10m,v_10m from C-grid to A-grid
  public :: uvtrafo            ! trafo of (u,v) rot <--> geo coord. system
  public :: rotate_pole_global ! Rotate the winds at the poles (global coord.)
  public :: rotate_pole_local  ! Rotate the winds at the poles (local  coord.)
  public :: scatter_wind_poles ! Scatter the winds at the poles over diamonds
  public :: average_pole_scalar! average poles from neigbours (on icosahedron)
  public :: average_pole_wind  ! average poles from neigbours (on icosahedron)
  !-----------------
  ! support plotting
  !-----------------
  public :: to_grads           ! dump t_atm to   grads file
  public :: from_grads         ! read t_atm from grads file
  public :: plot_t             ! plot time series
  public :: plot_v             ! plot profile
  !----------------
  ! transformations
  !----------------
  public :: vert_intp          ! vertical   interpolation
  public :: hor_intp           ! horizontal interpolation
  public :: lin_int_2          ! linear interpolation of 2 states
  public :: coarsen            ! select strides of gridpoints
  public :: hor_intp_coef_atm  ! set horizontal interpolation coefficients
  public :: hydro_balanc       ! hydrostatic balancing of p and t profile
  !--------------
  ! miscellaneous
  !--------------
  public :: levels             ! select subset of levels on isobaric coordinate
  public :: dump               ! dump to binary file
  public :: retrieve           ! retrieve from binary file
  public :: nwords             ! return size of all fields
  public :: real_wp            ! return rank 1 array holding all fields
  public :: cut                ! cutout region
  public :: p_sum_atm          ! sum over states on different processors
  public :: p_min_atm          ! min over states on different processors
  public :: p_max_atm          ! max over states on different processors
  public :: keep_diamond       ! keep only the desired diamond
!==============================================================================
  !---------------------
  ! data type definition
  !---------------------
  integer ,parameter :: nm = 120 ! number of fields in data type
  type t_atm
#ifndef TR15581
    !----------------------------
    ! function call depth counter
    !----------------------------
    integer               :: allocation_level
#endif
    !---------------------
    ! run type information
    !---------------------
    character(len=16)     :: name        = ''        ! arbitrary name
    character(len=8)      :: runtype     = ''        ! forecast, analysis, ..
    integer               :: runclass    = -1        ! 0..3 haupt,vor,ass,test
    integer               :: expid       = -1        ! experiment id
    integer               :: member      = -1        ! ensemble member
    integer               :: members     = -1        ! ensemble size
    integer               :: ensemble_id = -1        ! ensemble id
    !-----------------
    ! grid information
    !-----------------
    type(t_time)          :: time                    ! verification time
    type(t_time)          :: ref_time                ! reference time
    type(t_grid) ,pointer :: grid      =>NULL()      ! pointer to grid
    integer               :: lb (4)    = [1,1,1,1]   ! lower bounds
    integer               :: ub (4)    = [0,1,1,1]   ! upper bounds
    integer               :: lbd(4)    = [1,1,1,1]   ! lower bounds (dual grid)
    integer               :: ubd(4)    = [0,1,1,1]   ! upper bounds (dual grid)
    type(t_m)             :: m  (nm)                 ! atmospheric fields
    integer(i8)           :: size      = 0           ! allocated reals
    character             :: gp_wind   = ' '         ! current wind allocation
    !---------------------
    ! prognostic variables
    !---------------------
    real(wp)     _POINTER :: ps       (:,:,:,:) =>NULL() ! surface pressure
    real(wp)     _POINTER :: psr      (:,:,:,:) =>NULL() ! reference pressure
    real(wp)     _POINTER :: pp       (:,:,:,:) =>NULL() ! pressure deviation
    real(wp)     _POINTER :: dpsdt    (:,:,:,:) =>NULL() ! pressure tendency
    real(wp)     _POINTER :: t        (:,:,:,:) =>NULL() ! temperature
    real(wp)     _POINTER :: u        (:,:,:,:) =>NULL() ! longitudinal wind
    real(wp)     _POINTER :: v        (:,:,:,:) =>NULL() ! meridional   wind
    real(wp)     _POINTER :: w        (:,:,:,:) =>NULL() ! vertical     wind
    real(wp)     _POINTER :: q        (:,:,:,:) =>NULL() ! specific humidity
    real(wp)     _POINTER :: qcl      (:,:,:,:) =>NULL() ! spec. cl. liq.water
    real(wp)     _POINTER :: qci      (:,:,:,:) =>NULL() ! specific cloud ice
    real(wp)     _POINTER :: qr       (:,:,:,:) =>NULL() ! prognostic rain
    real(wp)     _POINTER :: qs       (:,:,:,:) =>NULL() ! prognostic snow
    real(wp)     _POINTER :: qg       (:,:,:,:) =>NULL() ! prognostic graupel
    real(wp)     _POINTER :: qh       (:,:,:,:) =>NULL() ! prognostic hail
    real(wp)     _POINTER :: tke      (:,:,:,:) =>NULL() ! prognostic tke
    !----------------------
    ! number concentrations
    !----------------------
    real(wp)     _POINTER :: nccloud  (:,:,:,:) =>NULL() ! number of cloud droplets per unit mass of air
    real(wp)     _POINTER :: ncrain   (:,:,:,:) =>NULL() ! specific number concentration of rain
    real(wp)     _POINTER :: ncsnow   (:,:,:,:) =>NULL() ! specific number concentration of snow
    real(wp)     _POINTER :: ncice    (:,:,:,:) =>NULL() ! number of cloud ice particles per unit mass of air
    real(wp)     _POINTER :: ncgraupel(:,:,:,:) =>NULL() ! specific number concentration of graupel
    real(wp)     _POINTER :: nchail   (:,:,:,:) =>NULL() ! specific number concentration of hail
    !-----
    ! soil
    !-----
    real(wp)     _POINTER :: t_so     (:,:,:,:) =>NULL() ! soil temperature
    real(wp)     _POINTER :: w_so     (:,:,:,:) =>NULL() ! soil humidity
    real(wp)     _POINTER :: w_so_ice (:,:,:,:) =>NULL() ! soil ice
    real(wp)     _POINTER :: smi      (:,:,:,:) =>NULL() ! soil moisture index
    !--------
    ! surface
    !--------
    real(wp)     _POINTER :: t_s      (:,:,:,:) =>NULL() ! surface temperature
    real(wp)     _POINTER :: qv_s     (:,:,:,:) =>NULL() ! surface spec.humid.
    real(wp)     _POINTER :: skt      (:,:,:,:) =>NULL() ! skin temperature
    real(wp)     _POINTER :: t_snow   (:,:,:,:) =>NULL() ! snow temperature
    real(wp)     _POINTER :: freshsnw (:,:,:,:) =>NULL() ! snow age
    real(wp)     _POINTER :: rho_snow (:,:,:,:) =>NULL() ! snow density
    real(wp)     _POINTER :: w_snow   (:,:,:,:) =>NULL() ! snow water content
    real(wp)     _POINTER :: snowc    (:,:,:,:) =>NULL() ! snow cover
    real(wp)     _POINTER :: w_i      (:,:,:,:) =>NULL() ! interception reserv.
    real(wp)     _POINTER :: iremispc (:,:,:,:) =>NULL() ! IR surf.emis.PCcoef.
    !------
    ! Flake
    !------
    real(wp)     _POINTER :: t_mnw_lk (:,:,:,:) =>NULL() ! Mean water column t.
    real(wp)     _POINTER :: t_wml_lk (:,:,:,:) =>NULL() ! Mixed-layer temp.
    real(wp)     _POINTER :: h_ml_lk  (:,:,:,:) =>NULL() ! Mixed-layer depth
    real(wp)     _POINTER :: t_bot_lk (:,:,:,:) =>NULL() ! Bottom temperature
    real(wp)     _POINTER :: c_t_lk   (:,:,:,:) =>NULL() ! Shape factor
!   real(wp)     _POINTER :: t_b1_lk  (:,:,:,:) =>NULL() !
!   real(wp)     _POINTER :: h_b1_lk  (:,:,:,:) =>NULL() !
    !------------------
    ! derived variables
    !------------------
    real(wp)     _POINTER :: clc      (:,:,:,:) =>NULL() ! cloud cover
    real(wp)     _POINTER :: ph       (:,:,:,:) =>NULL() ! half level pressure
    real(wp)     _POINTER :: pf       (:,:,:,:) =>NULL() ! full level pressure
    real(wp)     _POINTER :: rh       (:,:,:,:) =>NULL() ! relative humidity
    real(wp)     _POINTER :: t2m      (:,:,:,:) =>NULL() ! 2m temperature
    real(wp)     _POINTER :: rh2m     (:,:,:,:) =>NULL() ! 2m rel. humidity
    real(wp)     _POINTER :: geoh     (:,:,:,:) =>NULL() ! geopotential (half)
    real(wp)     _POINTER :: geof     (:,:,:,:) =>NULL() ! geop. (full levels)
    real(wp)     _POINTER :: tsurf    (:,:,:,:) =>NULL() ! surface temperature
    real(wp)     _POINTER :: fr_ice   (:,:,:,:) =>NULL() ! sea ice fraction
    real(wp)     _POINTER :: td2m     (:,:,:,:) =>NULL() ! 2m dewpoint temp.
    real(wp)     _POINTER :: tv       (:,:,:,:) =>NULL() ! virtual temperature
    real(wp)     _POINTER :: psi      (:,:,:,:) =>NULL() ! streamfunction
    real(wp)     _POINTER :: chi      (:,:,:,:) =>NULL() ! velocity potential
    real(wp)     _POINTER :: vrt      (:,:,:,:) =>NULL() ! rel. vorticity
    real(wp)     _POINTER :: div      (:,:,:,:) =>NULL() ! divergence
    !--------------------------------------------
    ! land-tile averaged diagnostic fields (ICON)
    !--------------------------------------------
    real(wp)     _POINTER :: t2m_land (:,:,:,:) =>NULL() ! 2m temperature   (land)
    real(wp)     _POINTER :: rh2m_land(:,:,:,:) =>NULL() ! 2m rel. humidity (land)
    real(wp)     _POINTER :: td2m_land(:,:,:,:) =>NULL() ! 2m dewpoint temp.(land)
    !-------------------------------------------------------
    ! averaged / accumulated / time-range-related quantities
    !-------------------------------------------------------
    real(wp)     _POINTER :: vmax_10m (:,:,:,:) =>NULL() ! max. 10m wind speed
    real(wp)     _POINTER :: tmin_2m  (:,:,:,:) =>NULL() ! min.  2m temperature
    real(wp)     _POINTER :: tmax_2m  (:,:,:,:) =>NULL() ! max.  2m temperature
    real(wp)     _POINTER :: tot_prec (:,:,:,:) =>NULL() ! total precipitation
    real(wp)     _POINTER :: aswdir_s (:,:,:,:) =>NULL() ! downw.direct SW rad
    real(wp)     _POINTER :: aswdifd_s(:,:,:,:) =>NULL() ! down. diffusive r.
    real(wp)     _POINTER :: alwd_s   (:,:,:,:) =>NULL() ! down. long wave rad.
    !---------------------
    ! diagnostic variables
    !---------------------
    real(wp)     _POINTER :: clct     (:,:,:,:) =>NULL() ! total cloud cover[%]
    real(wp)     _POINTER :: clcl     (:,:,:,:) =>NULL() ! low cloud cover  [%]
    real(wp)     _POINTER :: clcm     (:,:,:,:) =>NULL() ! med.cloud cover  [%]
    real(wp)     _POINTER :: clch     (:,:,:,:) =>NULL() ! high cloud cover [%]
    real(wp)     _POINTER :: u_10m    (:,:,:,:) =>NULL() ! 10m wind component
    real(wp)     _POINTER :: v_10m    (:,:,:,:) =>NULL() ! 10m wind component
    real(wp)     _POINTER :: pmsl     (:,:,:,:) =>NULL() ! sea level pressure
    real(wp)     _POINTER :: qv_dia   (:,:,:,:) =>NULL() ! specific humidity
    real(wp)     _POINTER :: qc_dia   (:,:,:,:) =>NULL() ! spec. cloud water
    real(wp)     _POINTER :: qi_dia   (:,:,:,:) =>NULL() ! spec. cloud ice
    real(wp)     _POINTER :: vis      (:,:,:,:) =>NULL() ! hor. visibility
    real(wp)     _POINTER :: ceiling  (:,:,:,:) =>NULL() ! ceiling
    real(wp)     _POINTER :: clc_rad  (:,:,:,:) =>NULL() ! cloud cover (mod.)
    !-----------------------------------------
    ! backup for wind components on the C-grid
    !-----------------------------------------
    real(wp)     _POINTER :: u_10m_c  (:,:,:,:) =>NULL() ! 10m wind component
    real(wp)     _POINTER :: v_10m_c  (:,:,:,:) =>NULL() ! 10m wind component
    real(wp)     _POINTER :: u_c      (:,:,:,:) =>NULL() ! u component
    real(wp)     _POINTER :: v_c      (:,:,:,:) =>NULL() ! v component
    !-------------------------------------
    ! fields to pass through (COSMO LETKF)
    !-------------------------------------
    real(wp)     _POINTER :: z0       (:,:,:,:) =>NULL() ! Rauhigkeitslaenge
    real(wp)     _POINTER :: plcov    (:,:,:,:) =>NULL() ! plant coverage
    real(wp)     _POINTER :: lai      (:,:,:,:) =>NULL() ! leaf area index
    real(wp)     _POINTER :: rootdp   (:,:,:,:) =>NULL() ! root depth
    real(wp)     _POINTER :: vio3     (:,:,:,:) =>NULL() ! vertical integr. O3
    real(wp)     _POINTER :: hmo3     (:,:,:,:) =>NULL() ! height of O3 maximum
    real(wp)     _POINTER :: for_e    (:,:,:,:) =>NULL() ! Bedeckung Nadelwald
    real(wp)     _POINTER :: for_d    (:,:,:,:) =>NULL() ! Bedeckung Laubwald
    real(wp)     _POINTER :: alb_dif  (:,:,:,:) =>NULL() ! diffusive solar albedo
    real(wp)     _POINTER :: aer_so4  (:,:,:,:) =>NULL() !
    real(wp)     _POINTER :: aer_dust (:,:,:,:) =>NULL() !
    real(wp)     _POINTER :: aer_org  (:,:,:,:) =>NULL() !
    real(wp)     _POINTER :: aer_bc   (:,:,:,:) =>NULL() !
    real(wp)     _POINTER :: aer_ss   (:,:,:,:) =>NULL() !
!   real(wp)     _POINTER :: sso_stdh (:,:,:,:) =>NULL() ! ~> moved to t_grid!
!   real(wp)     _POINTER :: sso_gamma(:,:,:,:) =>NULL() !
!   real(wp)     _POINTER :: sso_theta(:,:,:,:) =>NULL() !
!   real(wp)     _POINTER :: sso_sigma(:,:,:,:) =>NULL() !
    !-----------------------------------
    ! fields to pass through (GME LETKF)
    !-----------------------------------
    real(wp)     _POINTER :: h_ice    (:,:,:,:) =>NULL() ! ice thickness
    real(wp)     _POINTER :: t_ice    (:,:,:,:) =>NULL() ! ice surface temp.
    real(wp)     _POINTER :: h_snow   (:,:,:,:) =>NULL() ! snow depth       (m)
    real(wp)     _POINTER :: tqr      (:,:,:,:) =>NULL() ! integrated rainwater
    real(wp)     _POINTER :: tqs      (:,:,:,:) =>NULL() ! integrated snowwater
    real(wp)     _POINTER :: tqv      (:,:,:,:) =>NULL() ! integrated water vapor
    !--------------------------
    ! ICON prognostic variables
    !--------------------------
    real(wp)     _POINTER :: den      (:,:,:,:) =>NULL() ! air density
    real(wp)     _POINTER :: vptmp    (:,:,:,:) =>NULL() ! virt.pot.temperature
    !--------------------------
    ! ICON diagnostic variables
    !--------------------------
    real(wp)     _POINTER :: reff_qc  (:,:,:,:) =>NULL() ! effect. radius of cloud water
    real(wp)     _POINTER :: reff_qi  (:,:,:,:) =>NULL() ! " " of cloud ice
    real(wp)     _POINTER :: reff_qr  (:,:,:,:) =>NULL() ! " " of rain droplets
    real(wp)     _POINTER :: reff_qg  (:,:,:,:) =>NULL() ! " " of graupel
    real(wp)     _POINTER :: reff_qs  (:,:,:,:) =>NULL() ! " " of cloud snow
    real(wp)     _POINTER :: reff_qh  (:,:,:,:) =>NULL() ! " " of hail
    !----------------------
    ! LETKF specific fields
    !----------------------
    real(wp)     _POINTER :: f_inflat (:,:,:,:) =>NULL() ! inflation factor
    real(wp)     _POINTER :: p_energy (:,:,:,:) =>NULL() ! perturbation energy
    !-------------------------
    ! ICON ART specific fields
    !-------------------------
    real(wp)     _POINTER :: dusta    (:,:,:,:) =>NULL() ! mineral dust (fine)
    real(wp)     _POINTER :: dustb    (:,:,:,:) =>NULL() ! mineral dust (medium)
    real(wp)     _POINTER :: dustc    (:,:,:,:) =>NULL() ! mineral dust (coarse)
    !--------------------------------
    ! Optional input fields for RTTOV
    !--------------------------------
    real(wp)     _POINTER :: o3       (:,:,:,:) =>NULL() ! ozone concentration kg/kg
    real(wp)     _POINTER :: co2      (:,:,:,:) =>NULL() ! ozone concentration kg/kg
    ! allsky
    !----------------------------------------------------------
    ! inputs for allsky reflectances for cloud scheme TL and AL
    !----------------------------------------------------------
    real(wp)     _POINTER :: rcld     (:,:,:,:) =>NULL() ! stddev. satur. deficit
  end type t_atm
!==============================================================================
enum, bind(c)
  enumerator :: PS = 1   ,PSR      ,PP       ,DPSDT    ,&! 1..4
                T        ,U        ,V        ,W        ,&
                Q        ,QCL      ,QCI      ,QR       ,&
                QS       ,QG       ,QH       ,NCCLOUD  ,&
                NCRAIN   ,NCSNOW   ,NCICE    ,NCGRAUPEL,&
                NCHAIL   ,T_SO     ,W_SO     ,W_SO_ICE ,&! 21..24
                T_S      ,QV_S     ,T_SNOW   ,FRESHSNW ,&
                RHO_SNOW ,W_SNOW   ,W_I      ,IREMISPC ,&
                T_MNW_LK ,T_WML_LK ,H_ML_LK  ,T_BOT_LK ,&
                C_T_LK   ,CLC      ,PH       ,PF       ,&
                RH       ,T2M      ,RH2M     ,GEOH     ,&! 41..44
                GEOF     ,TSURF    ,FR_ICE   ,TD2M     ,&
                TV       ,PSI      ,CHI      ,VRT      ,&
                DIV      ,VMAX_10M ,TMIN_2M  ,TMAX_2M  ,&
                TOT_PREC ,ASWDIR_S ,ASWDIFD_S,ALWD_S   ,&
                CLCT     ,CLCL     ,CLCM     ,CLCH     ,&! 61..64
                U_10M    ,V_10M    ,PMSL     ,U_10M_C  ,&
                V_10M_C  ,U_C      ,V_C      ,Z0       ,&
                PLCOV    ,LAI      ,ROOTDP   ,VIO3     ,&
                HMO3     ,FOR_E    ,FOR_D    ,ALB_DIF  ,&
                AER_SO4  ,AER_DUST ,AER_ORG  ,AER_BC   ,&! 81..84
!               SSO_STDH ,SSO_GAMMA,SSO_THETA,SSO_SIGMA,&! ~> moved to t_grid!
                AER_SS   ,H_ICE    ,T_ICE    ,H_SNOW   ,&
                TQR      ,TQS      ,TQV      ,DEN      ,&
                VPTMP    ,F_INFLAT ,P_ENERGY ,TKE      ,&! 93..96
                DUSTA    ,DUSTB    ,DUSTC    ,QV_DIA   ,&
                QC_DIA   ,QI_DIA   ,O3       ,REFF_QC  ,&
                REFF_QI  ,REFF_QR  ,REFF_QG  ,REFF_QS  ,&!105..108
                REFF_QH  ,SKT      ,CO2      ,VIS      ,&
                CEIL     ,SNOWC    ,RCLD     ,SMI      ,&
                T2M_LAND ,TD2M_LAND,RH2M_LAND,CLC_RAD    !117..120
end enum
!==============================================================================
  !-----------
  ! interfaces
  !-----------
  interface construct
    module procedure construct_state
    module procedure construct_states
    module procedure construct_states2
  end interface construct
!------------------------------------------------------------------------------
  interface destruct
    module procedure destruct_state
    module procedure destruct_states
    module procedure destruct_states2
  end interface destruct
!------------------------------------------------------------------------------
  interface allocate
    module procedure allocate_state
    module procedure allocate_state_trange
  end interface allocate
!------------------------------------------------------------------------------
  interface deallocate
    module procedure deallocate_state
  end interface deallocate
!------------------------------------------------------------------------------
  interface print
    module procedure print_state
  end interface print
!------------------------------------------------------------------------------
  interface dump
    module procedure dump_state
  end interface dump
!------------------------------------------------------------------------------
  interface retrieve
    module procedure retrieve_state
  end interface retrieve
!------------------------------------------------------------------------------
  interface assignment (=)
    module procedure assign_atm_to_atm
    module procedure assign_real_to_atm
    module procedure assign_real_to_atms
    module procedure assign_real_to_atms2
    module procedure assign_array_to_atm
  end interface ! assignment (=)
!------------------------------------------------------------------------------
  interface operator (+)
    module procedure atm_plus_atm
  end interface ! operator (+)
!------------------------------------------------------------------------------
  interface operator (-)
    module procedure atm_minus_atm
  end interface ! operator (-)
!------------------------------------------------------------------------------
  interface operator (*)
    module procedure atm_times_r
    module procedure atm_times_atm
    module procedure   r_times_atm
    module procedure   i_times_atm
  end interface ! operator (*)
!------------------------------------------------------------------------------
  interface operator (/)
    module procedure atm_over_atm  ! t_atm / t_atm
    module procedure atm_over_r    ! t_atm / real(wp)
    module procedure atm_over_i    ! t_atm / integer
  end interface ! operator (/)
!------------------------------------------------------------------------------
  interface operator (**)
    module procedure atm_power_i
  end interface ! operator (**)
!------------------------------------------------------------------------------
  interface p_sum_atm
    module procedure p_sum_atm
    module procedure p_sum_atm1
    module procedure p_sum_atm2
  end interface p_sum_atm
!------------------------------------------------------------------------------
  interface max
    module procedure max_r_atm
  end interface max
!------------------------------------------------------------------------------
  interface sqrt
    module procedure sqrt_atm
  end interface ! sqrt
!------------------------------------------------------------------------------
  interface sum
    module procedure sum_atm
  end interface sum
!------------------------------------------------------------------------------
  interface lin_int_2
    module procedure lin_int_2_atm
  end  interface lin_int_2
!------------------------------------------------------------------------------
  interface lgtd
    module procedure lgt_atm
  end interface lgtd
!------------------------------------------------------------------------------
  interface lgti
    module procedure lgti_atm
  end interface lgti
!------------------------------------------------------------------------------
  interface lgtd_ad
    module procedure lgt_ad_atm
  end interface lgtd_ad
!------------------------------------------------------------------------------
  interface lgti_ad
    module procedure lgti_ad_atm
  end interface lgti_ad
!------------------------------------------------------------------------------
  interface plot_t
    module procedure plot_atm_t
  end  interface plot_t
!------------------------------------------------------------------------------
  interface plot_v
    module procedure plot_atm_v
  end  interface plot_v
!------------------------------------------------------------------------------
  interface nwords
    module procedure nwords_atm
  end  interface nwords
!------------------------------------------------------------------------------
  interface real_wp
    module procedure real_atm
  end interface real_wp
!------------------------------------------------------------------------------
#ifndef TR15581
  interface
    pure subroutine delete_m (x,i)
      use mo_memory
      integer    ,intent(in) :: i
      type (t_m) ,intent(in) :: x(i)
    end subroutine delete_m
  end interface
#endif
!------------------------------------------------------------------------------
  interface set_pointers
    module procedure set_atm_pointers
  end interface set_pointers
!------------------------------------------------------------------------------
  interface set_p
    module procedure set_p
    module procedure set_p1
  end interface set_p
!------------------------------------------------------------------------------
  interface set_geo
    module procedure set_geo
    module procedure set_geo1
  end interface set_geo
!------------------------------------------------------------------------------
  interface set_ps
    module procedure set_ps
    module procedure set_ps1
  end interface set_ps
!------------------------------------------------------------------------------
  interface set_pf
    module procedure set_pf
    module procedure set_pf1
  end interface set_pf
!------------------------------------------------------------------------------
  interface set_ph
    module procedure set_ph
    module procedure set_ph1
  end interface set_ph
!------------------------------------------------------------------------------
  interface set_t
    module procedure set_t
    module procedure set_t1
  end interface set_t
!------------------------------------------------------------------------------
  interface set_tv
    module procedure set_tv
    module procedure set_tv1
  end interface set_tv
!------------------------------------------------------------------------------
  interface set_rh
    module procedure set_rh
    module procedure set_rh1
  end interface set_rh
!------------------------------------------------------------------------------
  interface set_rh2m
    module procedure set_rh2m
    module procedure set_rh2m1
  end interface set_rh2m
!------------------------------------------------------------------------------
  interface set_rh2m_land
    module procedure set_rh2m_land
    module procedure set_rh2m_land1
  end interface set_rh2m_land
!------------------------------------------------------------------------------
  interface set_q
    module procedure set_q
    module procedure set_q1
  end interface set_q
!------------------------------------------------------------------------------
  interface to_grads
    module procedure atm_to_grads
  end interface to_grads
!------------------------------------------------------------------------------
  interface from_grads
    module procedure atm_from_grads
  end interface from_grads
!------------------------------------------------------------------------------
  interface cut
    module procedure cut_atm
  end interface cut
!------------------------------------------------------------------------------
  interface merge
    module procedure merge_atm
  end interface
!------------------------------------------------------------------------------
  interface coarsen
    module procedure coarsen_atm
  end interface coarsen
!------------------------------------------------------------------------------
  interface hor_intp
    module procedure hor_intp_atm
  end interface hor_intp
!------------------------------------------------------------------------------
  interface derive_params           ! transform fields (eg. q <-> rh etc.)
    module procedure derive_params1 ! vector version
    module procedure derive_params0 ! scalar version
  end interface derive_params
!------------------------------------------------------------------------------
  interface keep_params
    module procedure keep_params
    module procedure keep_params_name
  end interface keep_params
!------------------------------------------------------------------------------
  interface A_to_C_grid           ! transform u,v,u_10m,v_10m from A- to C-grid
    module procedure A_to_C_grid  ! scalar version
    module procedure A_to_C_grid1 ! vector version
  end interface A_to_C_grid
!------------------------------------------------------------------------------
  interface C_to_A_grid           ! transform u,v,u_10m,v_10m from C- to A-grid
    module procedure C_to_A_grid  ! scalar version
    module procedure C_to_A_grid1 ! vector version
  end interface C_to_A_grid
!------------------------------------------------------------------------------
  interface rotate_pole_global           ! Rotate winds at poles (global coord.)
    module procedure rotate_pole_global  ! scalar version
    module procedure rotate_pole_global1 ! vector version
  end interface rotate_pole_global
!------------------------------------------------------------------------------
  interface rotate_pole_local            ! Rotate winds at poles (local  coord.)
    module procedure rotate_pole_local   ! scalar version
    module procedure rotate_pole_local1  ! vector version
  end interface rotate_pole_local
!------------------------------------------------------------------------------
  interface uvtrafo               ! trafo of (u,v) rot <--> geo coord. system
    module procedure uvtrafo      ! scalar version
    module procedure uvtrafo1     ! vector version
  end interface uvtrafo
!==============================================================================
contains
!==============================================================================


!==============================================================================
! Subroutines derived from COSMO-code:
!   (u,v)-trafo geogr.<--> rotated coord.system
!==============================================================================
elemental SUBROUTINE uvrot2uv (urot, vrot, rlat, rlon, pollat, pollon, u, v)
!------------------------------------------------------------------------------
!
! Description:
!   This routine converts the wind components u and v from the rotated system
!   to the real geographical system.
!
! Method:
!   Transformation formulas for converting between these two systems.
!
!------------------------------------------------------------------------------

! Parameter list:
  REAL (KIND=wp), INTENT (IN)          ::    &
       urot, vrot,     & ! wind components in the rotated grid
       rlat, rlon,     & ! latitude and longitude in the true geographical system
       pollat, pollon    ! latitude and longitude of the north pole of the
! rotated grid

  REAL (KIND=wp), INTENT (OUT)         ::    &
       u, v              ! wind components in the true geographical system

! Local variables

  REAL (KIND=wp)                       ::    &
       zsinpol, zcospol, zlonp, zlat, zarg1, zarg2, znorm

  REAL (KIND=wp), parameter            ::    &
       zpir18 = d2r

!------------------------------------------------------------------------------
! Begin subroutine uvrot2uv
!------------------------------------------------------------------------------

! Converting from degree to radians
  zsinpol = SIN(pollat * zpir18)
  zcospol = COS(pollat * zpir18)
  zlonp   = (pollon-rlon) * zpir18
  zlat    =         rlat  * zpir18

  zarg1   = zcospol*SIN(zlonp)
  zarg2   = zsinpol*COS(zlat) - zcospol*SIN(zlat)*COS(zlonp)
  znorm   = 1._wp/SQRT(zarg1**2 + zarg2**2)

! Convert the u- and v-components
  u       =   urot*zarg2*znorm + vrot*zarg1*znorm
  v       = - urot*zarg1*znorm + vrot*zarg2*znorm

END SUBROUTINE uvrot2uv

!------------------------------------------------------------------------------

  subroutine uvrot2uv_vec(u, v, rlat, rlon, pollat, pollon, idim)
  !----------------------------------------------------------------------------
  !
  ! Description:
  !   This routine converts the wind components u and v from the rotated
  !   system to the real geographical system. This is the vectorized form
  !   of the routine above, i.e. the computation is for a whole 2D field.
  !
  ! Method:
  !   Transformation formulas for converting between these two systems.
  !
  !----------------------------------------------------------------------------
  integer  ,intent(in)    :: idim        ! dimension of the field
  real(wp) ,intent(inout) :: u   (idim)  ! wind components
  real(wp) ,intent(inout) :: v   (idim)  !
  real(wp), intent(in)    :: rlat(idim)  ! coordinates in the true
  real(wp), intent(in)    :: rlon(idim)  !   geographical system (radians)
  real(wp), intent(in)    :: pollat      ! latitude (degree) and
  real(wp), intent(in)    :: pollon      ! longitude of the north pole of the
                                         ! rotated grid

    real(wp) :: zsinpol, zcospol, zlonp, zlat, zarg1, zarg2, znorm, zugeo, zvgeo
    integer  :: i

    ! Converting from degree to radians
    zsinpol = sin(pollat * d2r)
    zcospol = cos(pollat * d2r)

    do i = 1, idim
      zlonp   = d2r * pollon - rlon(i)
      zlat    =                rlat(i)

      zarg1   = zcospol*sin(zlonp)
      zarg2   = zsinpol*cos(zlat) - zcospol*sin(zlat)*cos(zlonp)
      znorm   = 1._wp/sqrt(zarg1**2 + zarg2**2)

      ! Convert the u- and v-components
      zugeo   = ( u(i)*zarg2 + v(i)*zarg1)*znorm
      zvgeo   = (-u(i)*zarg1 + v(i)*zarg2)*znorm
      u(i)    = zugeo
      v(i)    = zvgeo
    enddo

  end subroutine uvrot2uv_vec

!------------------------------------------------------------------------------

elemental SUBROUTINE uv2uvrot(u, v, rlat, rlon, pollat, pollon, urot, vrot)

!------------------------------------------------------------------------------
!
! Description:
!   This routine converts the wind components u and v from the real
!   geographical system to the rotated system.
!
! Method:
!   Transformation formulas for converting between these two systems.
!
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Parameter list:
REAL (KIND=wp), INTENT (IN)          ::    &
  u   , v   ,     & ! wind components in the true geographical system
  rlat, rlon,     & ! coordinates in the true geographical system
  pollat, pollon    ! latitude and longitude of the north pole of the
                    ! rotated grid

REAL (KIND=wp), INTENT (OUT)         ::    &
  urot, vrot        ! wind components in the rotated grid

! Local variables

REAL (KIND=wp)                       ::    &
  zsinpol, zcospol, zlonp, zlat, zarg1, zarg2, znorm

REAL (KIND=wp), parameter            ::    &
  zpir18 = d2r

!------------------------------------------------------------------------------
! Begin Subroutine uv2uvrot
!------------------------------------------------------------------------------

  zsinpol = SIN(pollat * zpir18)
  zcospol = COS(pollat * zpir18)
  zlonp   = (pollon-rlon) * zpir18
  zlat    =         rlat  * zpir18

  zarg1   = zcospol*SIN(zlonp)
  zarg2   = zsinpol*COS(zlat) - zcospol*SIN(zlat)*COS(zlonp)
  znorm   = 1.0_wp/SQRT( zarg1**2 + zarg2**2 )

! Transform the u and v wind components
  urot   =  u*zarg2*znorm - v*zarg1*znorm
  vrot   =  u*zarg1*znorm + v*zarg2*znorm

END SUBROUTINE uv2uvrot

!------------------------------------------------------------------------------
  subroutine uv2uvrot_vec (u, v, rlat, rlon, pollat, pollon, idim)
  !----------------------------------------------------------------------------
  !
  ! Description:
  !   This routine converts the wind components u and v from the real
  !   geographical system to the rotated system. This is the vectorized form
  !   of the routine above, i.e. the computation is for a whole array.
  !
  ! Method:
  !   Transformation formulas for converting between these two systems.
  !
  !----------------------------------------------------------------------------
  integer,  intent(in)    :: idim       ! dimensions of the field
  real(wp), intent(inout) :: u   (idim) ! wind components in the true
  real(wp), intent(inout) :: v   (idim) !   geographical system
  real(wp), intent(in)    :: rlat(idim) ! coordinates in the true geographical
  real(wp), intent(in)    :: rlon(idim) !   system (radians)
  real(wp), intent(in)    :: pollat     ! latitude (degree) and
  real(wp), intent(in)    :: pollon     ! longitude of the north pole of the
                                        ! rotated grid

    real(wp) :: zsinpol, zcospol, zlonp, zlat, zarg1, zarg2, znorm, zurot, zvrot
    integer  :: i

    zsinpol = sin ( pollat * d2r )
    zcospol = cos ( pollat * d2r )

    do i = 1, idim
      zlonp   =  d2r * pollon - rlon(i)
      zlat    =                 rlat(i)

      zarg1   = zcospol*sin(zlonp)
      zarg2   = zsinpol*cos(zlat) - zcospol*sin(zlat)*cos(zlonp)
      znorm   = 1.0_wp/sqrt( zarg1**2 + zarg2**2 )

      ! Transform the u and v wind components
      zurot = (u(i)*zarg2 - v(i)*zarg1)*znorm
      zvrot = (u(i)*zarg1 + v(i)*zarg2)*znorm
      u(i)  = zurot
      v(i)  = zvrot
    enddo

  end subroutine uv2uvrot_vec
!-------------------------------------------------------------------------------
  subroutine uvtrafo1 (atms,trafo)
    !-------------------------------------
    ! Vector version of subroutine uvtrafo
    !-------------------------------------
    type(t_atm) ,intent(inout) :: atms (:)
    integer     ,intent(in)    :: trafo
    integer                    :: i        ! loop index
    do i = 1, size (atms)
      call uvtrafo (atms (i), trafo)
    end do
  end subroutine uvtrafo1
!------------------------------------------------------------------------------
  subroutine uvtrafo (atm,trafo)
  !----------------------------------------------------------------
  ! Description:
  ! This routine converts the wind components u and v from the real
  ! geographical system to the rotated system for 3d
  !
  ! Method:
  ! uses the above subroutines from the COSMO-code
  !-----------------------------------------------------------------
  type(t_atm) ,intent(inout) :: atm    ! atm. state
  integer     ,intent(in)    :: trafo  ! flag,0: geo-->rot,1:rot--> geo

    integer :: k, lx, ly, il, iu, jl, ju

    if (.not. atm% grid% rot) return

    iu = atm% grid% ub(1)
    il = atm% grid% lb(1)
    ju = atm% grid% ub(2)
    jl = atm% grid% lb(2)

    lx = iu - il + 1
    ly = ju - jl + 1

    ! loop over vertical levels
    select case(trafo)
    case(0)
      do k=1, atm% grid% nz
        call uv2uvrot_vec (atm%       u     (  :  ,  :  ,k,1), &
                           atm%       v     (  :  ,  :  ,k,1), &
                           atm% grid% rlat  (il:iu,jl:ju,1,1), &
                           atm% grid% rlon  (il:iu,jl:ju,1,1), &
                           atm% grid% dlatr,                   &
                           atm% grid% dlonr,                   &
                                      lx*ly                    )
      end do
    case(1)
      do k=1, atm% grid% nz
        call uvrot2uv_vec (atm%       u     (  :  ,  :  ,k,1), &
                           atm%       v     (  :  ,  :  ,k,1), &
                           atm% grid% rlat  (il:iu,jl:ju,1,1), &
                           atm% grid% rlon  (il:iu,jl:ju,1,1), &
                           atm% grid% dlatr,                   &
                           atm% grid% dlonr,                   &
                           lx*ly                               )
       end do
     end select
  end subroutine uvtrafo
!==============================================================================

  subroutine construct_states2 (states, grid, lb, ub, nz, alloc, time, &
                                ref_time, template, name               )
  !--------------------------------------------------------------
  ! Description:
  ! Initialise array of derived type 't_atm' before first use
  !
  ! Method:
  ! Call scalar version for each array element in turn
  !--------------------------------------------------------------
  type (t_atm)      ,intent(out)          :: states(:,:) ! variable to set
  type (t_grid)     ,pointer    ,optional :: grid        ! grid information
  integer           ,intent(in) ,optional :: lb(:)       ! lower bound
  integer           ,intent(in) ,optional :: ub(:)       ! upper bound
  integer           ,intent(in) ,optional :: nz          ! number of levels
  character (len=*) ,intent(in) ,optional :: alloc       ! fields to allocate
  type (t_time)     ,intent(in) ,optional :: time        ! verification time
  type (t_time)     ,intent(in) ,optional :: ref_time    ! reference time
  type (t_atm)      ,intent(in) ,optional :: template    ! template
  character (len=*) ,intent(in) ,optional :: name        ! arbitrary name
    integer i
    do i=1,size(states,2)
      call construct (states(:,i), grid,                               &
                      lb, ub, nz, alloc, time, ref_time, template, name)
    end do
  end subroutine construct_states2
!------------------------------------------------------------------------------
  subroutine construct_states (states, grid, lb, ub, nz, alloc, time, &
                               ref_time, template, name               )
  !--------------------------------------------------------------
  ! Description:
  ! Initialise array of derived type 't_atm' before first use
  !
  ! Method:
  ! Call scalar version for each array element in turn
  !--------------------------------------------------------------
  type (t_atm)      ,intent(out)          :: states(:) ! variable to initialize
  type (t_grid)     ,pointer    ,optional :: grid      ! grid information
  integer           ,intent(in) ,optional :: lb(:)     ! lower bound
  integer           ,intent(in) ,optional :: ub(:)     ! upper bound
  integer           ,intent(in) ,optional :: nz        ! number of levels
  character (len=*) ,intent(in) ,optional :: alloc     ! fields to allocate
  type (t_time)     ,intent(in) ,optional :: time      ! verification time
  type (t_time)     ,intent(in) ,optional :: ref_time  ! reference time
  type (t_atm)      ,intent(in) ,optional :: template  ! template
  character (len=*) ,intent(in) ,optional :: name      ! arbitrary name
    integer i
    do i=1,size(states)
      call construct (states(i),                                             &
                      grid, lb, ub, nz, alloc, time, ref_time, template, name)
    end do
  end subroutine construct_states
!------------------------------------------------------------------------------
  subroutine construct_state (state, grid, lb, ub, nz, alloc, time, ref_time, &
                              template, name                                  )
  !---------------------------------------------------------------------------
  ! Description:
  ! Initialise variable of derived type 't_atm' before first use
  !
  ! Set pointer to grid information variable 'grid'
  ! Allocate arrays
  ! Use optional parameters 'lb','ub'to overwrite array bound values
  ! Use optional parameter  'alloc' to specify components to allocate
  ! Use optional parameters 'time','ref_time' to set the respective components
  ! Use optional 'template' to derive above parameters from another state
  ! Either parameter 'grid' or 'template' is mandatory
  !---------------------------------------------------------------------------
  type (t_atm)     ,intent(out)          :: state    ! variable to initialize
  type (t_grid)    ,pointer    ,optional :: grid     ! grid information
  integer          ,intent(in) ,optional :: lb(:)    ! lower bound
  integer          ,intent(in) ,optional :: ub(:)    ! upper bound
  integer          ,intent(in) ,optional :: nz       ! number of levels
  character(len=*) ,intent(in) ,optional :: alloc    ! fields to allocate
  type (t_time)    ,intent(in) ,optional :: time     ! verification time
  type (t_time)    ,intent(in) ,optional :: ref_time ! reference time
  type (t_atm)     ,intent(in) ,optional :: template ! template
  character(len=*) ,intent(in) ,optional :: name     ! arbitrary name

    integer            :: i
    character(len=512) :: al
    integer            :: ll(4), uu(4)
    integer            :: ll_d(4), uu_d(4)      ! Bounds of dual grid
    ALLOCATIONLEVEL (state)
    !---------------------
    ! name, time, and grid
    !---------------------
    nullify (state% grid)
    if (present (name))     state% name     = name
    if (present (time))     state% ref_time = time
    if (present (template)) then
      state% time     =  template% time
      state% ref_time =  template% ref_time
      state% grid     => template% grid
      state% gp_wind  =  template% gp_wind
      state% member   =  template% member
      state% members  =  template% members
      ll              =  template% lb
      uu              =  template% ub
      ll_d            =  template% lbd          ! Bounds of dual grid
      uu_d            =  template% ubd          ! (if applicable)
    endif
    if (present (time))     state% time     = time
    if (present (ref_time)) state% ref_time = ref_time
    if (present (grid)) then
      state% grid => grid
      ll          =  grid% lb
      uu          =  grid% ub
      ll_d        =  grid% lb_d
      uu_d        =  grid% ub_d
    endif
    if (.not.associated(state% grid)) &
      call finish('construct_state','no grid present')
    !---------------
    ! determine size
    !---------------
    if (present(ub)) uu(1:size(ub)) = ub
    if (present(lb)) ll(1:size(lb)) = lb
    if (present(nz)) uu(3)          = nz
    state% lb   =  ll
    state% ub   =  uu
    if (uu_d(1) < ll_d(1)) then
       ll_d     =  ll
       uu_d     =  uu
    else
       uu_d(3)  =  uu(3)        ! Align number of levels only
    end if
    state% lbd  =  ll_d
    state% ubd  =  uu_d
    !------------------
    ! initialize memory
    !------------------
    state%m%i% name = (/'ps       ','psr      ','pp       ','dpsdt    ',&
                        't        ','u        ','v        ','w        ',&
                        'q        ','qcl      ','qci      ','qr       ',&
                        'qs       ','qg       ','qh       ','nccloud  ',&
                        'ncrain   ','ncsnow   ','ncice    ','ncgraupel',&
                        'nchail   ','t_so     ','w_so     ','w_so_ice ',&
                        't_s      ','qv_s     ','t_snow   ','freshsnw ',&
                        'rho_snow ','w_snow   ','w_i      ','iremispc ',&
                        't_mnw_lk ','t_wml_lk ','h_ml_lk  ','t_bot_lk ',&
                        'c_t_lk   ','clc      ','ph       ','pf       ',&
                        'rh       ','t2m      ','rh2m     ','geoh     ',&
                        'geof     ','tsurf    ','fr_ice   ','td2m     ',&
                        'tv       ','psi      ','chi      ','vrt      ',&
                        'div      ','vmax_10m ','tmin_2m  ','tmax_2m  ',&
                        'tot_prec ','aswdir_s ','aswdifd_s','alwd_s   ',&
                        'clct     ','clcl     ','clcm     ','clch     ',&
                        'u_10m    ','v_10m    ','pmsl     ','u_10m_c  ',&
                        'v_10m_c  ','u_c      ','v_c      ','z0       ',&
                        'plcov    ','lai      ','rootdp   ','vio3     ',&
                        'hmo3     ','for_e    ','for_d    ','alb_dif  ',&
                        'aer_so4  ','aer_dust ','aer_org  ','aer_bc   ',&
!                       'sso_stdh ','sso_gamma','sso_theta','sso_sigma',&
                        'aer_ss   ','h_ice    ','t_ice    ','h_snow   ',&
                        'tqr      ','tqs      ','tqv      ','den      ',&
                        'vptmp    ','f_inflat ','p_energy ','tke      ',&
                        'dusta    ','dustb    ','dustc    ','qv_dia   ',&
                        'qc_dia   ','qi_dia   ','o3       ','reff_qc  ',&
                        'reff_qi  ','reff_qr  ','reff_qg  ','reff_qs  ',&
                        'reff_qh  ','skt      ','co2      ','vis      ',&
                        'ceiling  ','snowc    ','rcld     ','smi      ',&
                        't2m_land ','td2m_land','rh2m_land','clc_rad  '/)

    do i=1,size(state%m)
      state%m(i)%i% lb    = ll
      state%m(i)%i% ub    = uu
    end do
    state%m%i%            nn    = state% grid% nn
    state%m(VRT     )%i%  lb    = ll_d       ! vorticity on dual grid
    state%m(VRT     )%i%  ub    = uu_d       ! may have different bounds
    state%m(PS      )%i%  lb(3) = 1          ! ps:      (:,:,1)
    state%m(PS      )%i%  ub(3) = 1
    state%m(PSR     )%i%  lb(3) = 1          ! psr:     (:,:,1)
    state%m(PSR     )%i%  ub(3) = 1
    state%m(DPSDT   )%i%  lb(3) = 1          ! dpsdt:   (:,:,1)
    state%m(DPSDT   )%i%  ub(3) = 1
    state%m(W       )%i%  ub(3) = uu(3) + 1  ! w:       (:,:,:ke+1)
    state%m(TKE     )%i%  ub(3) = uu(3) + 1  ! TKE:     (:,:,:ke+1)
    state%m(T_SO    )%i%  lb(3) = 0          ! t_so     (:,:,0:8)
    state%m(T_SO    )%i%  ub(3) = state% grid% ns
    state%m(W_SO    )%i%  lb(3) = 1          ! w_so     (:,:,1:8)
    state%m(W_SO    )%i%  ub(3) = state% grid% ns
    state%m(W_SO_ICE)%i%  lb(3) = 1          ! w_so_ice (:,:,1:8)
    state%m(W_SO_ICE)%i%  ub(3) = state% grid% ns
    state%m(SMI     )%i%  lb(3) = 1          ! smi      (:,:,1:8)
    state%m(SMI     )%i%  ub(3) = state% grid% ns
    state%m(T_S     )%i%  lb(3) = 1
    state%m(T_S     )%i%  ub(3) = 1          ! t_s
    state%m(QV_S    )%i%  lb(3) = 1
    state%m(QV_S    )%i%  ub(3) = 1          ! qv_s
    state%m(SKT     )%i%  lb(3) = 1
    state%m(SKT     )%i%  ub(3) = 1          ! skt:     (:,:,1)
    state%m(VIS     )%i%  lb(3) = 1
    state%m(VIS     )%i%  ub(3) = 1          ! vis:     (:,:,1)
    state%m(CEIL    )%i%  lb(3) = 1
    state%m(CEIL    )%i%  ub(3) = 1          ! ceiling: (:,:,1)
    state%m(T_SNOW  )%i%  lb(3) = 1
    state%m(T_SNOW  )%i%  ub(3) = 1          ! t_snow
    state%m(FRESHSNW)%i%  lb(3) = 1
    state%m(FRESHSNW)%i%  ub(3) = 1          ! freshsnw
    state%m(RHO_SNOW)%i%  lb(3) = 1
    state%m(RHO_SNOW)%i%  ub(3) = 1          ! rho_snow
    state%m(W_SNOW  )%i%  lb(3) = 1
    state%m(W_SNOW  )%i%  ub(3) = 1          ! w_snow
    state%m(SNOWC   )%i%  lb(3) = 1
    state%m(SNOWC   )%i%  ub(3) = 1          ! snowc
    state%m(W_I     )%i%  lb(3) = 1
    state%m(W_I     )%i%  ub(3) = 1          ! w_i
    state%m(IREMISPC)%i%  lb(3) = 0
    state%m(IREMISPC)%i%  ub(3) = 6          ! iremispc
    do i = T_MNW_LK, C_T_LK
      state%m(i)%i% lb(3) = 1                ! Flake
      state%m(i)%i% ub(3) = 1
    end do
    state%m(PH      )%i%  ub(3) = uu(3) + 1  ! ph:      (:,:,:ke+1)
    state%m(GEOH    )%i%  ub(3) = uu(3) + 1  ! geoh:    (:,:,:ke+1)
    state%m(RCLD    )%i%  ub(3) = uu(3) + 1  ! rcld:    (:,:,:ke+1) ! allsky: rcld is on 1/2 levels
    state%m(TSURF   )%i%  lb(3) = 1          ! tsurf:   (:,:,1)
    state%m(TSURF   )%i%  ub(3) = 1
    state%m(FR_ICE  )%i%  lb(3) = 1          ! fr_ice:  (:,:,1)
    state%m(FR_ICE  )%i%  ub(3) = 1
    state%m(TD2M    )%i%  lb(3) = 1          ! td2m:    (:,:,1)
    state%m(TD2M    )%i%  ub(3) = 1
    do i = T2M, RH2M
      state%m(i)%i% lb(3) = 1          ! t2m, rh2m
      state%m(i)%i% ub(3) = 1
    end do
    do i = T2M_LAND, RH2M_LAND
      state%m(i)%i% lb(3) = 1          ! t2m_land, td2m_land, rh2m_land
      state%m(i)%i% ub(3) = 1
    end do
    do i = VMAX_10M, V_10M_C
      state%m(i)%i% lb(3) = 1          ! vmax_10m,tmin_2m,tmax_2m, clct,clcl,..
      state%m(i)%i% ub(3) = 1          ! ..clcm,clch, tot_prec, u_10m,v_10m,._c
    end do                             ! ..pmsl
    do i = Z0, TQV
      state%m(i)%i% lb(3) = 1          ! fields to pass through (single level)
      state%m(i)%i% ub(3) = 1          !
    end do
    !++++++++++++++++++++++++++++++++++++++++++++++
    ! Temporary hack for incomplete ICON-ART fields
    !++++++++++++++++++++++++++++++++++++++++++++++
    do i = DUSTA, DUSTC
      state%m(i)%i% lb(3) = 25         ! only levels >= 25 available
    end do

    select case (state% grid% vct)
    case (VCT_P_HYB)                     ! GME, IFS, HRM :
      state%m(PSR )%i%  ref   = .true.   !   psr
      state%m(PH  )%i%  ref   = .true.   !   ph
      state%m(PF  )%i%  ref   = .true.   !   pf
    case (VCT_Z_HYB, VCT_Z_GEN)          ! COSMO, ICON :
      state%m(GEOH)%i%  ref   = .true.   !   geoh
      state%m(GEOF)%i%  ref   = .true.   !   geof
    end select
    !-----------------------------
    ! determine fields to allocate
    !-----------------------------
    if(present(template)) then
      state%m%i% alloc = template%m%i% alloc
!     state%m%i% lb(3) = template%m%i% lb(3)
!     state%m%i% ub(3) = template%m%i% ub(3)
      state%m%i% code  = template%m%i% code
      state%m%i% table = template%m%i% table
    endif
    if(present(alloc)) then
      if ( len_trim(alloc) > len(al) - 2 ) &
        call finish ('construct_state','len(alloc)>len(al)-2')
      al(1:1) = ' '
      al(2: ) = trim(alloc)
      state%m%i% alloc = .false.
      if (index(al,      ' all ')/=0) state%m           %i% alloc = .true.
      if (index(al,       ' ps ')/=0) state%m(       PS)%i% alloc = .true.
      if (index(al,      ' psr ')/=0) state%m(      PSR)%i% alloc = .true.
      if (index(al,       ' pp ')/=0) state%m(       PP)%i% alloc = .true.
      if (index(al,    ' dpsdt ')/=0) state%m(    DPSDT)%i% alloc = .true.
      if (index(al,        ' t ')/=0) state%m(        T)%i% alloc = .true.
      if (index(al,        ' u ')/=0) state%m(        U)%i% alloc = .true.
      if (index(al,        ' v ')/=0) state%m(        V)%i% alloc = .true.
      if (index(al,        ' w ')/=0) state%m(        W)%i% alloc = .true.
      if (index(al,        ' q ')/=0) state%m(        Q)%i% alloc = .true.
      if (index(al,      ' qcl ')/=0) state%m(      QCL)%i% alloc = .true.
      if (index(al,      ' qci ')/=0) state%m(      QCI)%i% alloc = .true.
      if (index(al,       ' qr ')/=0) state%m(       QR)%i% alloc = .true.
      if (index(al,       ' qs ')/=0) state%m(       QS)%i% alloc = .true.
      if (index(al,       ' qg ')/=0) state%m(       QG)%i% alloc = .true.
      if (index(al,       ' qh ')/=0) state%m(       QH)%i% alloc = .true.
      if (index(al,  ' nccloud ')/=0) state%m(  NCCLOUD)%i% alloc = .true.
      if (index(al,   ' ncrain ')/=0) state%m(   NCRAIN)%i% alloc = .true.
      if (index(al,   ' ncsnow ')/=0) state%m(   NCSNOW)%i% alloc = .true.
      if (index(al,    ' ncice ')/=0) state%m(    NCICE)%i% alloc = .true.
      if (index(al,' ncgraupel ')/=0) state%m(NCGRAUPEL)%i% alloc = .true.
      if (index(al,   ' nchail ')/=0) state%m(   NCHAIL)%i% alloc = .true.
      if (index(al,     ' t_so ')/=0) state%m(     T_SO)%i% alloc = .true.
      if (index(al,     ' w_so ')/=0) state%m(     W_SO)%i% alloc = .true.
      if (index(al, ' w_so_ice ')/=0) state%m( W_SO_ICE)%i% alloc = .true.
      if (index(al,      ' t_s ')/=0) state%m(      T_S)%i% alloc = .true.
      if (index(al,     ' qv_s ')/=0) state%m(     QV_S)%i% alloc = .true.
      if (index(al,   ' t_snow ')/=0) state%m(   T_SNOW)%i% alloc = .true.
      if (index(al, ' freshsnw ')/=0) state%m( FRESHSNW)%i% alloc = .true.
      if (index(al, ' rho_snow ')/=0) state%m( RHO_SNOW)%i% alloc = .true.
      if (index(al,   ' w_snow ')/=0) state%m(   W_SNOW)%i% alloc = .true.
      if (index(al,    ' snowc ')/=0) state%m(    SNOWC)%i% alloc = .true.
      if (index(al,      ' w_i ')/=0) state%m(      W_I)%i% alloc = .true.
      if (index(al, ' iremispc ')/=0) state%m( IREMISPC)%i% alloc = .true.
      if (index(al, ' t_mnw_lk ')/=0) state%m( T_MNW_LK)%i% alloc = .true.
      if (index(al, ' t_wml_lk ')/=0) state%m( T_WML_LK)%i% alloc = .true.
      if (index(al,  ' h_ml_lk ')/=0) state%m(  H_ML_LK)%i% alloc = .true.
      if (index(al, ' t_bot_lk ')/=0) state%m( T_BOT_LK)%i% alloc = .true.
      if (index(al,   ' c_t_lk ')/=0) state%m(   C_T_LK)%i% alloc = .true.
      if (index(al,      ' clc ')/=0) state%m(      CLC)%i% alloc = .true.
      if (index(al,       ' ph ')/=0) state%m(       PH)%i% alloc = .true.
      if (index(al,       ' pf ')/=0) state%m(       PF)%i% alloc = .true.
      if (index(al,       ' rh ')/=0) state%m(       RH)%i% alloc = .true.
      if (index(al,      ' t2m ')/=0) state%m(      T2M)%i% alloc = .true.
      if (index(al,     ' rh2m ')/=0) state%m(     RH2M)%i% alloc = .true.
      if (index(al,     ' geoh ')/=0) state%m(     GEOH)%i% alloc = .true.
      if (index(al,     ' geof ')/=0) state%m(     GEOF)%i% alloc = .true.
      if (index(al,    ' tsurf ')/=0) state%m(    TSURF)%i% alloc = .true.
      if (index(al,   ' fr_ice ')/=0) state%m(   FR_ICE)%i% alloc = .true.
      if (index(al,     ' td2m ')/=0) state%m(     TD2M)%i% alloc = .true.
      if (index(al,       ' tv ')/=0) state%m(       TV)%i% alloc = .true.
      if (index(al,      ' psi ')/=0) state%m(      PSI)%i% alloc = .true.
      if (index(al,      ' chi ')/=0) state%m(      CHI)%i% alloc = .true.
      if (index(al,      ' vrt ')/=0) state%m(      VRT)%i% alloc = .true.
      if (index(al,      ' div ')/=0) state%m(      DIV)%i% alloc = .true.
      if (index(al, ' vmax_10m ')/=0) state%m( VMAX_10M)%i% alloc = .true.
      if (index(al,  ' tmin_2m ')/=0) state%m(  TMIN_2M)%i% alloc = .true.
      if (index(al,  ' tmax_2m ')/=0) state%m(  TMAX_2M)%i% alloc = .true.
      if (index(al, ' tot_prec ')/=0) state%m( TOT_PREC)%i% alloc = .true.
      if (index(al, ' aswdir_s ')/=0) state%m( ASWDIR_S)%i% alloc = .true.
      if (index(al,' aswdifd_s ')/=0) state%m(ASWDIFD_S)%i% alloc = .true.
      if (index(al,   ' alwd_s ')/=0) state%m(   ALWD_S)%i% alloc = .true.
      if (index(al,     ' clct ')/=0) state%m(     CLCT)%i% alloc = .true.
      if (index(al,     ' clcl ')/=0) state%m(     CLCL)%i% alloc = .true.
      if (index(al,     ' clcm ')/=0) state%m(     CLCM)%i% alloc = .true.
      if (index(al,     ' clch ')/=0) state%m(     CLCH)%i% alloc = .true.
      if (index(al,    ' u_10m ')/=0) state%m(    U_10M)%i% alloc = .true.
      if (index(al,    ' v_10m ')/=0) state%m(    V_10M)%i% alloc = .true.
      if (index(al,     ' pmsl ')/=0) state%m(     PMSL)%i% alloc = .true.
      if (index(al,  ' u_10m_c ')/=0) state%m(  U_10M_C)%i% alloc = .true.
      if (index(al,  ' v_10m_c ')/=0) state%m(  V_10M_C)%i% alloc = .true.
      if (index(al,      ' u_c ')/=0) state%m(      U_C)%i% alloc = .true.
      if (index(al,      ' v_c ')/=0) state%m(      V_C)%i% alloc = .true.
      if (index(al,       ' z0 ')/=0) state%m(       Z0)%i% alloc = .true.
      if (index(al,    ' plcov ')/=0) state%m(    PLCOV)%i% alloc = .true.
      if (index(al,      ' lai ')/=0) state%m(      LAI)%i% alloc = .true.
      if (index(al,   ' rootdp ')/=0) state%m(   ROOTDP)%i% alloc = .true.
      if (index(al,     ' vio3 ')/=0) state%m(     VIO3)%i% alloc = .true.
      if (index(al,     ' hmo3 ')/=0) state%m(     HMO3)%i% alloc = .true.
      if (index(al,    ' for_e ')/=0) state%m(    FOR_E)%i% alloc = .true.
      if (index(al,    ' for_d ')/=0) state%m(    FOR_D)%i% alloc = .true.
      if (index(al,  ' alb_dif ')/=0) state%m(  ALB_DIF)%i% alloc = .true.
      if (index(al,  ' aer_so4 ')/=0) state%m(  AER_SO4)%i% alloc = .true.
      if (index(al, ' aer_dust ')/=0) state%m( AER_DUST)%i% alloc = .true.
      if (index(al,  ' aer_org ')/=0) state%m(  AER_ORG)%i% alloc = .true.
      if (index(al,   ' aer_bc ')/=0) state%m(   AER_BC)%i% alloc = .true.
      if (index(al,   ' aer_ss ')/=0) state%m(   AER_SS)%i% alloc = .true.
!     if (index(al, ' sso_stdh ')/=0) state%m( SSO_STDH)%i% alloc = .true.
!     if (index(al,' sso_gamma ')/=0) state%m(SSO_GAMMA)%i% alloc = .true.
!     if (index(al,' sso_theta ')/=0) state%m(SSO_THETA)%i% alloc = .true.
!     if (index(al,' sso_sigma ')/=0) state%m(SSO_SIGMA)%i% alloc = .true.
      if (index(al,    ' h_ice ')/=0) state%m(    H_ICE)%i% alloc = .true.
      if (index(al,    ' t_ice ')/=0) state%m(    T_ICE)%i% alloc = .true.
      if (index(al,   ' h_snow ')/=0) state%m(   H_SNOW)%i% alloc = .true.
      if (index(al,      ' tqr ')/=0) state%m(      TQR)%i% alloc = .true.
      if (index(al,      ' tqs ')/=0) state%m(      TQS)%i% alloc = .true.
      if (index(al,      ' tqv ')/=0) state%m(      TQV)%i% alloc = .true.
      if (index(al,      ' den ')/=0) state%m(      DEN)%i% alloc = .true.
      if (index(al,    ' vptmp ')/=0) state%m(    VPTMP)%i% alloc = .true.
      if (index(al, ' f_inflat ')/=0) state%m( F_INFLAT)%i% alloc = .true.
      if (index(al, ' p_energy ')/=0) state%m( P_ENERGY)%i% alloc = .true.
      if (index(al,      ' tke ')/=0) state%m(      TKE)%i% alloc = .true.
      if (index(al,    ' dusta ')/=0) state%m(    DUSTA)%i% alloc = .true.
      if (index(al,    ' dustb ')/=0) state%m(    DUSTB)%i% alloc = .true.
      if (index(al,    ' dustc ')/=0) state%m(    DUSTC)%i% alloc = .true.
      if (index(al,   ' qv_dia ')/=0) state%m(   QV_DIA)%i% alloc = .true.
      if (index(al,   ' qc_dia ')/=0) state%m(   QC_DIA)%i% alloc = .true.
      if (index(al,   ' qi_dia ')/=0) state%m(   QI_DIA)%i% alloc = .true.
      if (index(al,       ' o3 ')/=0) state%m(       O3)%i% alloc = .true.
      if (index(al,  ' reff_qc ')/=0) state%m(  REFF_QC)%i% alloc = .true.
      if (index(al,  ' reff_qi ')/=0) state%m(  REFF_QI)%i% alloc = .true.
      if (index(al,  ' reff_qr ')/=0) state%m(  REFF_QR)%i% alloc = .true.
      if (index(al,  ' reff_qg ')/=0) state%m(  REFF_QG)%i% alloc = .true.
      if (index(al,  ' reff_qs ')/=0) state%m(  REFF_QS)%i% alloc = .true.
      if (index(al,  ' reff_qh ')/=0) state%m(  REFF_QH)%i% alloc = .true.
      if (index(al,      ' skt ')/=0) state%m(      SKT)%i% alloc = .true.
      if (index(al,      ' co2 ')/=0) state%m(      CO2)%i% alloc = .true.
      if (index(al,      ' vis ')/=0) state%m(      VIS)%i% alloc = .true.
      if (index(al,  ' ceiling ')/=0) state%m(     CEIL)%i% alloc = .true.
      if (index(al,     ' rcld ')/=0) state%m(     RCLD)%i% alloc = .true.
      if (index(al,      ' smi ')/=0) state%m(      SMI)%i% alloc = .true.
      if (index(al, ' t2m_land ')/=0) state%m( T2M_LAND)%i% alloc = .true.
      if (index(al,' td2m_land ')/=0) state%m(TD2M_LAND)%i% alloc = .true.
      if (index(al,' rh2m_land ')/=0) state%m(RH2M_LAND)%i% alloc = .true.
      if (index(al,  ' clc_rad ')/=0) state%m(  CLC_RAD)%i% alloc = .true.
    endif
    call update (state%m)
    call set_pointers (state)
  end subroutine construct_state
!------------------------------------------------------------------------------
  subroutine allocate_state_trange (a, field, trange, nlev, code, table, &
                                    idx, ierr                            )
  type (t_atm)     TARGET,intent(inout) :: a         ! atm_state-variable
  character (len=*)      ,intent(in)    :: field     ! mnemonic of field
  integer                ,intent(in)    :: trange(:) ! time range
  integer      ,optional ,intent(in)    :: nlev      ! number of levels
  integer      ,optional ,intent(in)    :: code      ! GRIB code  no.
  integer      ,optional ,intent(in)    :: table     ! GRIB table no.
  integer      ,optional ,intent(out)   :: idx       ! index
  integer      ,optional ,intent(out)   :: ierr      ! error return code
  !--------------------------------------------
  ! allocate single array pointer components of
  ! of variable of type t_atm
  ! for variables with time range indicator
  !--------------------------------------------
    integer :: i, n

    n = size (trange)
    if (all (trange <= 0)) then
      call allocate_state (a, field, code=code, nlev=nlev,         &
                                     table=table, idx=i, ierr=ierr )
      if (i > 0) a%m(i)% i% tr = 0
    else
      call allocate_state (a, field, code=code, nlev=n,            &
                                     table=table, idx=i, ierr=ierr )
      if (i > 0) a%m(i)% i% tr(1:n) = trange
    endif
    if (present (idx)) idx = i

  end subroutine allocate_state_trange
!------------------------------------------------------------------------------
  elemental subroutine allocate_state (a, field, nlev, code, table, idx, ierr)
  type (t_atm)     TARGET,intent(inout) :: a     ! atm_state-variable
  character (len=*)      ,intent(in)    :: field ! mnemonic of field
  integer      ,optional ,intent(in)    :: nlev  ! number of levels to allocate
  integer      ,optional ,intent(in)    :: code  ! GRIB code  no.
  integer      ,optional ,intent(in)    :: table ! GRIB table no.
  integer      ,optional ,intent(out)   :: idx   ! index
  integer      ,optional ,intent(out)   :: ierr  ! error return code
  !--------------------------------------------
  ! allocate single array pointer components of
  ! of variable of type t_atm
  !--------------------------------------------
    integer           :: i, j, n, ierror
    character(len=16) :: pars(nm)

    call split (pars, field, n)
    ierror =  0
    i      = -1
    do j = 1, n
      select case (pars(j))
      case ('ps')
        i =  PS        ;call allocate (a%m(i),ptr=a% ps                 )
      case ('psr')
        i =  PSR       ;call allocate (a%m(i),ptr=a% psr                )
      case ('pp')
        i =  PP        ;call allocate (a%m(i),ptr=a% pp       ,nlev=nlev)
      case ('dpsdt')
        i =  DPSDT     ;call allocate (a%m(i),ptr=a% dpsdt              )
      case ('t')
        i =  T         ;call allocate (a%m(i),ptr=a% t        ,nlev=nlev)
      case ('u')
        i =  U         ;call allocate (a%m(i),ptr=a% u        ,nlev=nlev)
      case ('v')
        i =  V         ;call allocate (a%m(i),ptr=a% v        ,nlev=nlev)
      case ('w')
        i =  W         ;call allocate (a%m(i),ptr=a% w        ,nlev=nlev)
      case ('q')
        i =  Q         ;call allocate (a%m(i),ptr=a% q        ,nlev=nlev)
      case ('qcl')
        i =  QCL       ;call allocate (a%m(i),ptr=a% qcl      ,nlev=nlev)
      case ('qci')
        i =  QCI       ;call allocate (a%m(i),ptr=a% qci      ,nlev=nlev)
      case ('qr')
        i =  QR        ;call allocate (a%m(i),ptr=a% qr       ,nlev=nlev)
      case ('qs')
        i =  QS        ;call allocate (a%m(i),ptr=a% qs       ,nlev=nlev)
      case ('qg')
        i =  QG        ;call allocate (a%m(i),ptr=a% qg       ,nlev=nlev)
      case ('qh')
        i =  QH        ;call allocate (a%m(i),ptr=a% qh       ,nlev=nlev)
      case ('nccloud')
        i =  NCCLOUD   ;call allocate (a%m(i),ptr=a% nccloud  ,nlev=nlev)
      case ('ncrain')
        i =  NCRAIN    ;call allocate (a%m(i),ptr=a% ncrain   ,nlev=nlev)
      case ('ncsnow')
        i =  NCSNOW    ;call allocate (a%m(i),ptr=a% ncsnow   ,nlev=nlev)
      case ('ncice')
        i =  NCICE     ;call allocate (a%m(i),ptr=a% ncice    ,nlev=nlev)
      case ('ncgraupel')
        i =  NCGRAUPEL ;call allocate (a%m(i),ptr=a% ncgraupel,nlev=nlev)
      case ('nchail')
        i =  NCHAIL    ;call allocate (a%m(i),ptr=a% nchail   ,nlev=nlev)
      case ('t_so')
        i =  T_SO      ;call allocate (a%m(i),ptr=a% t_so     ,nlev=nlev)
      case ('w_so')
        i =  W_SO      ;call allocate (a%m(i),ptr=a% w_so     ,nlev=nlev)
      case ('w_so_ice')
        i =  W_SO_ICE  ;call allocate (a%m(i),ptr=a% w_so_ice ,nlev=nlev)
      case ('t_s')
        i =  T_S       ;call allocate (a%m(i),ptr=a% t_s                )
      case ('qv_s')
        i =  QV_S      ;call allocate (a%m(i),ptr=a% qv_s               )
      case ('t_snow')
        i =  T_SNOW    ;call allocate (a%m(i),ptr=a% t_snow             )
      case ('freshsnw')
        i =  FRESHSNW  ;call allocate (a%m(i),ptr=a% freshsnw           )
      case ('rho_snow')
        i =  RHO_SNOW  ;call allocate (a%m(i),ptr=a% rho_snow           )
      case ('w_snow')
        i =  W_SNOW    ;call allocate (a%m(i),ptr=a% w_snow             )
      case ('snowc')
        i =  SNOWC     ;call allocate (a%m(i),ptr=a% snowc              )
      case ('w_i')
        i =  W_I       ;call allocate (a%m(i),ptr=a% w_i                )
      case ('iremispc')
        i =  IREMISPC  ;call allocate (a%m(i),ptr=a% iremispc           )
      case ('t_mnw_lk')
        i =  T_MNW_LK  ;call allocate (a%m(i),ptr=a% t_mnw_lk           )
      case ('t_wml_lk')
        i =  T_WML_LK  ;call allocate (a%m(i),ptr=a% t_wml_lk           )
      case ('h_ml_lk' )
        i =  H_ML_LK   ;call allocate (a%m(i),ptr=a% h_ml_lk            )
      case ('t_bot_lk')
        i =  T_BOT_LK  ;call allocate (a%m(i),ptr=a% t_bot_lk           )
      case ('c_t_lk'  )
        i =  C_T_LK    ;call allocate (a%m(i),ptr=a% c_t_lk             )
      case ('clc')
        i =  CLC       ;call allocate (a%m(i),ptr=a% clc      ,nlev=nlev)
      case ('ph')
        i =  PH        ;call allocate (a%m(i),ptr=a% ph       ,nlev=nlev)
      case ('pf')
        i =  PF        ;call allocate (a%m(i),ptr=a% pf       ,nlev=nlev)
      case ('rh')
        i =  RH        ;call allocate (a%m(i),ptr=a% rh       ,nlev=nlev)
      case ('t2m')
        i =  T2M       ;call allocate (a%m(i),ptr=a% t2m                )
      case ('rh2m')
        i =  RH2M      ;call allocate (a%m(i),ptr=a% rh2m               )
      case ('geoh')
        i =  GEOH      ;call allocate (a%m(i),ptr=a% geoh     ,nlev=nlev)
      case ('geof')
        i =  GEOF      ;call allocate (a%m(i),ptr=a% geof     ,nlev=nlev)
      case ('tsurf')
        i =  TSURF     ;call allocate (a%m(i),ptr=a% tsurf              )
      case ('fr_ice')
        i =  FR_ICE    ;call allocate (a%m(i),ptr=a% fr_ice             )
      case ('td2m')
        i =  TD2M      ;call allocate (a%m(i),ptr=a% td2m               )
      case ('tv')
        i =  TV        ;call allocate (a%m(i),ptr=a% tv       ,nlev=nlev)
      case ('psi')
        i =  PSI       ;call allocate (a%m(i),ptr=a% psi      ,nlev=nlev)
      case ('chi')
        i =  CHI       ;call allocate (a%m(i),ptr=a% chi      ,nlev=nlev)
      case ('vrt')
        i =  VRT       ;call allocate (a%m(i),ptr=a% vrt      ,nlev=nlev)
      case ('div')
        i =  DIV       ;call allocate (a%m(i),ptr=a% div      ,nlev=nlev)
      case ('vmax_10m')
        i =  VMAX_10M  ;call allocate (a%m(i),ptr=a% vmax_10m ,nlev=nlev)
      case ('tmin_2m ')
        i =  TMIN_2M   ;call allocate (a%m(i),ptr=a% tmin_2m  ,nlev=nlev)
      case ('tmax_2m ')
        i =  TMAX_2M   ;call allocate (a%m(i),ptr=a% tmax_2m  ,nlev=nlev)
      case ('tot_prec')
        i =  TOT_PREC  ;call allocate (a%m(i),ptr=a% tot_prec ,nlev=nlev)
      case ('aswdir_s')
        i =  ASWDIR_S  ;call allocate (a%m(i),ptr=a% aswdir_s ,nlev=nlev)
      case ('aswdifd_s')
        i =  aswdifd_s ;call allocate (a%m(i),ptr=a% aswdifd_s,nlev=nlev)
      case ('alwd_s')
        i =  ALWD_S    ;call allocate (a%m(i),ptr=a% alwd_s   ,nlev=nlev)
      case ('clct')
        i =  CLCT      ;call allocate (a%m(i),ptr=a% clct               )
      case ('clcl')
        i =  CLCL      ;call allocate (a%m(i),ptr=a% clcl               )
      case ('clcm')
        i =  CLCM      ;call allocate (a%m(i),ptr=a% clcm               )
      case ('clch')
        i =  CLCH      ;call allocate (a%m(i),ptr=a% clch               )
      case ('u_10m')
        i =  U_10M     ;call allocate (a%m(i),ptr=a% u_10m              )
      case ('v_10m')
        i =  V_10M     ;call allocate (a%m(i),ptr=a% v_10m              )
      case ('pmsl')
        i =  PMSL      ;call allocate (a%m(i),ptr=a% pmsl               )
      case ('u_10m_c')
        i =  U_10M_C   ;call allocate (a%m(i),ptr=a% u_10m_c            )
      case ('v_10m_c')
        i =  V_10M_C   ;call allocate (a%m(i),ptr=a% v_10m_c            )
      case ('u_c')
        i =  U_C       ;call allocate (a%m(i),ptr=a% u_c                )
      case ('v_c')
        i =  V_C       ;call allocate (a%m(i),ptr=a% v_c                )
      case ('z0')
        i =  Z0        ;call allocate (a%m(i),ptr=a% z0                 )
      case ('plcov')
        i =  PLCOV     ;call allocate (a%m(i),ptr=a% plcov              )
      case ('lai')
        i =  LAI       ;call allocate (a%m(i),ptr=a% lai                )
      case ('rootdp')
        i =  ROOTDP    ;call allocate (a%m(i),ptr=a% rootdp             )
      case ('vio3')
        i =  VIO3      ;call allocate (a%m(i),ptr=a% vio3               )
      case ('hmo3')
        i =  HMO3      ;call allocate (a%m(i),ptr=a% hmo3               )
      case ('for_e')
        i =  FOR_E     ;call allocate (a%m(i),ptr=a% for_e              )
      case ('for_d')
        i =  FOR_D     ;call allocate (a%m(i),ptr=a% for_d              )
      case ('alb_dif')
        i =  ALB_DIF   ;call allocate (a%m(i),ptr=a% alb_dif            )
      case ('aer_so4')
        i =  AER_SO4   ;call allocate (a%m(i),ptr=a% aer_so4            )
      case ('aer_dust')
        i =  AER_DUST  ;call allocate (a%m(i),ptr=a% aer_dust           )
      case ('aer_org')
        i =  AER_ORG   ;call allocate (a%m(i),ptr=a% aer_org            )
      case ('aer_bc')
        i =  AER_BC    ;call allocate (a%m(i),ptr=a% aer_bc             )
      case ('aer_ss')
        i =  AER_SS    ;call allocate (a%m(i),ptr=a% aer_ss             )
!     case ('sso_stdh')
!       i =  SSO_STDH ;call allocate (a%m(i),ptr=a% sso_stdh            )
!     case ('sso_gamma')
!       i =  SSO_GAMMA;call allocate (a%m(i),ptr=a% sso_gamma           )
!     case ('sso_theta')
!       i =  SSO_THETA;call allocate (a%m(i),ptr=a% sso_theta           )
!     case ('sso_sigma')
!       i =  SSO_SIGMA;call allocate (a%m(i),ptr=a% sso_sigma           )
      case ('h_ice')
        i =  H_ICE     ;call allocate (a%m(i),ptr=a% h_ice              )
      case ('t_ice')
        i =  T_ICE     ;call allocate (a%m(i),ptr=a% t_ice              )
      case ('h_snow')
        i =  H_SNOW    ;call allocate (a%m(i),ptr=a% h_snow             )
      case ('tqr')
        i =  TQR       ;call allocate (a%m(i),ptr=a% tqr                )
      case ('tqs')
        i =  TQS       ;call allocate (a%m(i),ptr=a% tqs                )
      case ('tqv')
        i =  TQV       ;call allocate (a%m(i),ptr=a% tqv                )
      case ('den')
        i =  DEN       ;call allocate (a%m(i),ptr=a% den                )
      case ('vptmp')
        i =  VPTMP     ;call allocate (a%m(i),ptr=a% vptmp              )
      case ('f_inflat')
        i =  F_INFLAT  ;call allocate (a%m(i),ptr=a% f_inflat           )
      case ('p_energy')
        i =  P_ENERGY  ;call allocate (a%m(i),ptr=a% p_energy           )
      case ('tke')
        i =  TKE       ;call allocate (a%m(i),ptr=a% tke      ,nlev=nlev)
      case ('dusta')
        i =  DUSTA     ;call allocate (a%m(i),ptr=a% dusta              )
      case ('dustb')
        i =  DUSTB     ;call allocate (a%m(i),ptr=a% dustb              )
      case ('dustc')
        i =  DUSTC     ;call allocate (a%m(i),ptr=a% dustc              )
      case ('qv_dia')
        i =  QV_DIA    ;call allocate (a%m(i),ptr=a% qv_dia             )
      case ('qc_dia')
        i =  QC_DIA    ;call allocate (a%m(i),ptr=a% qc_dia             )
      case ('qi_dia')
        i =  QI_DIA    ;call allocate (a%m(i),ptr=a% qi_dia             )
      case ('o3')
        i =  O3        ;call allocate (a%m(i),ptr=a% o3                 )
      case ('reff_qc')
        i =  REFF_QC   ;call allocate (a%m(i),ptr=a% reff_qc  ,nlev=nlev)
       case ('reff_qi')
        i =  REFF_QI   ;call allocate (a%m(i),ptr=a% reff_qi  ,nlev=nlev)
      case ('reff_qr')
        i =  REFF_QR   ;call allocate (a%m(i),ptr=a% reff_qr  ,nlev=nlev)
      case ('reff_qg')
        i =  REFF_QG   ;call allocate (a%m(i),ptr=a% reff_qg  ,nlev=nlev)
      case ('reff_qs')
        i =  REFF_QS   ;call allocate (a%m(i),ptr=a% reff_qs  ,nlev=nlev)
      case ('reff_qh')
        i =  REFF_QH   ;call allocate (a%m(i),ptr=a% reff_qh  ,nlev=nlev)
      case ('skt')
        i =  SKT       ;call allocate (a%m(i),ptr=a% skt                )
      case ('co2')
        i =  CO2       ;call allocate (a%m(i),ptr=a% co2                )
      case ('vis')
        i =  VIS       ;call allocate (a%m(i),ptr=a% vis                )
      case ('ceiling')
        i =  CEIL      ;call allocate (a%m(i),ptr=a% ceiling            )
      case ('rcld')
        i =  RCLD      ;call allocate (a%m(i),ptr=a% rcld     ,nlev=nlev)
      case ('smi')
        i =  SMI       ;call allocate (a%m(i),ptr=a% smi      ,nlev=nlev)
      case ('t2m_land')
        i =  T2M_LAND  ;call allocate (a%m(i),ptr=a% t2m_land           )
      case ('td2m_land')
        i =  TD2M_LAND ;call allocate (a%m(i),ptr=a% td2m_land          )
      case ('rh2m_land')
        i =  RH2M_LAND ;call allocate (a%m(i),ptr=a% rh2m_land          )
      case ('clc_rad')
        i =  CLC_RAD   ;call allocate (a%m(i),ptr=a% clc_rad  ,nlev=nlev)
      case default
        ierror = 1
      end select
      if (ierror == 0) then
        if (present (code )) a%m(i)%i% code  = code
        if (present (table)) a%m(i)%i% table = table
      endif
    end do
    a% size = nwords(a% m)
    if (present (ierr )) ierr = ierror
    if (present (idx  )) idx  = i
  end subroutine allocate_state
!------------------------------------------------------------------------------
  elemental subroutine deallocate_state (state, fields)
  type (t_atm)      ,intent(inout) :: state  ! atm_state-variable
  character (len=*) ,intent(in)    :: fields ! mnemonic of fields to dealloc.
  !--------------------------------------------
  ! deallocate single array pointer components of
  ! of variable of type t_atm
  !--------------------------------------------
    character(len=16) :: pars(nm)
    integer           :: i, n

    call split (pars, fields, n)
    do i = 1, n
      select case (pars(i))
      case ('ps')
        call deallocate (state%m(PS       )); nullify (state% ps)
      case ('psr')
        call deallocate (state%m(PSR      )); nullify (state% psr)
      case ('pp')
        call deallocate (state%m(PP       )); nullify (state% pp)
      case ('dpsdt')
        call deallocate (state%m(DPSDT    )); nullify (state% dpsdt)
      case ('t')
        call deallocate (state%m(T        )); nullify (state% t)
      case ('u')
        call deallocate (state%m(U        )); nullify (state% u)
      case ('v')
        call deallocate (state%m(V        )); nullify (state% v)
      case ('w')
        call deallocate (state%m(W        )); nullify (state% w)
      case ('q')
        call deallocate (state%m(Q        )); nullify (state% q)
      case ('qcl')
        call deallocate (state%m(QCL      )); nullify (state% qcl)
      case ('qci')
        call deallocate (state%m(QCI      )); nullify (state% qci)
      case ('qr')
        call deallocate (state%m(QR       )); nullify (state% qr)
      case ('qs')
        call deallocate (state%m(QS       )); nullify (state% qs)
      case ('qg')
        call deallocate (state%m(QG       )); nullify (state% qg )
      case ('qh')
        call deallocate (state%m(QH       )); nullify (state% qh )
      case ('nccloud')
        call deallocate (state%m(NCCLOUD  )); nullify (state% nccloud )
      case ('ncrain')
        call deallocate (state%m(NCRAIN   )); nullify (state% ncrain )
      case ('ncsnow')
        call deallocate (state%m(NCSNOW   )); nullify (state% ncsnow )
      case ('ncice')
        call deallocate (state%m(NCICE    )); nullify (state% ncice )
      case ('ncgraupel')
        call deallocate (state%m(NCGRAUPEL)); nullify (state% ncgraupel )
      case ('nchail')
        call deallocate (state%m(NCHAIL   )); nullify (state% nchail )
      case ('t_so')
        call deallocate (state%m(T_SO     )); nullify (state% t_so )
      case ('w_so')
        call deallocate (state%m(W_SO     )); nullify (state% w_so )
      case ('w_so_ice')
        call deallocate (state%m(W_SO_ICE )); nullify (state% w_so_ice )
      case ('t_s')
        call deallocate (state%m(T_S      )); nullify (state% t_s )
      case ('qv_s')
        call deallocate (state%m(QV_S     )); nullify (state% qv_s )
      case ('t_snow')
        call deallocate (state%m(T_SNOW   )); nullify (state% t_snow )
      case ('freshsnw')
        call deallocate (state%m(FRESHSNW )); nullify (state% freshsnw )
      case ('rho_snow')
        call deallocate (state%m(RHO_SNOW )); nullify (state% rho_snow )
      case ('w_snow')
        call deallocate (state%m(W_SNOW   )); nullify (state% w_snow )
      case ('snowc')
        call deallocate (state%m(SNOWC    )); nullify (state% snowc )
      case ('w_i')
        call deallocate (state%m(W_I      )); nullify (state% w_i )
      case ('iremispc')
        call deallocate (state%m(IREMISPC )); nullify (state% iremispc )
      case ('t_mnw_lk')
        call deallocate (state%m(T_MNW_LK )); nullify (state% t_mnw_lk )
      case ('t_wml_lk')
        call deallocate (state%m(T_WML_LK )); nullify (state% t_wml_lk )
      case ('h_ml_lk')
        call deallocate (state%m(H_ML_LK  )); nullify (state% h_ml_lk )
      case ('t_bot_lk')
        call deallocate (state%m(T_BOT_LK )); nullify (state% t_bot_lk )
      case ('c_t_lk')
        call deallocate (state%m(C_T_LK   )); nullify (state% c_t_lk )
      case ('clc')
        call deallocate (state%m(CLC      )); nullify (state% clc)
      case ('ph')
        call deallocate (state%m(PH       )); nullify (state% ph)
      case ('pf')
        call deallocate (state%m(PF       )); nullify (state% pf)
      case ('rh')
        call deallocate (state%m(RH       )); nullify (state% rh)
      case ('t2m')
        call deallocate (state%m(T2M      )); nullify (state% t2m)
      case ('rh2m')
        call deallocate (state%m(RH2M     )); nullify (state% rh2m)
      case ('geoh')
        call deallocate (state%m(GEOH     )); nullify (state% geoh)
      case ('geof')
        call deallocate (state%m(GEOF     )); nullify (state% geof)
      case ('tsurf')
        call deallocate (state%m(TSURF    )); nullify (state% tsurf)
      case ('fr_ice')
        call deallocate (state%m(FR_ICE   )); nullify (state% fr_ice)
      case ('td2m')
        call deallocate (state%m(TD2M     )); nullify (state% td2m)
      case ('tv')
        call deallocate (state%m(TV       )); nullify (state% tv)
      case ('psi')
        call deallocate (state%m(PSI      )); nullify (state% psi)
      case ('chi')
        call deallocate (state%m(CHI      )); nullify (state% chi)
      case ('vrt')
        call deallocate (state%m(VRT      )); nullify (state% vrt)
      case ('div')
        call deallocate (state%m(DIV      )); nullify (state% div)
      case ('vmax_10m')
        call deallocate (state%m(VMAX_10M )); nullify (state% vmax_10m)
      case ('tmin_2m')
        call deallocate (state%m(TMIN_2M  )); nullify (state% tmin_2m)
      case ('tmax_2m')
        call deallocate (state%m(TMAX_2M  )); nullify (state% tmax_2m)
      case ('tot_prec')
        call deallocate (state%m(TOT_PREC )); nullify (state% tot_prec)
      case ('aswdir_s')
        call deallocate (state%m(ASWDIR_S )); nullify (state% aswdir_s)
      case ('aswdifd_s')
        call deallocate (state%m(ASWDIFD_S)); nullify (state% aswdifd_s)
      case ('alwd_s')
        call deallocate (state%m(ALWD_S   )); nullify (state% alwd_s)
      case ('clct')
        call deallocate (state%m(CLCT     )); nullify (state% clct)
      case ('clcl')
        call deallocate (state%m(CLCL     )); nullify (state% clcl)
      case ('clcm')
        call deallocate (state%m(CLCM     )); nullify (state% clcm)
      case ('clch')
        call deallocate (state%m(CLCH     )); nullify (state% clch)
      case ('u_10m')
        call deallocate (state%m(U_10M    )); nullify (state% u_10m)
      case ('v_10m')
        call deallocate (state%m(V_10M    )); nullify (state% v_10m)
      case ('pmsl')
        call deallocate (state%m(PMSL     )); nullify (state% pmsl)
      case ('u_10m_c')
        call deallocate (state%m(U_10M_C  )); nullify (state% u_10m_c)
      case ('v_10m_c')
        call deallocate (state%m(V_10M_C  )); nullify (state% v_10m_c)
      case ('u_c')
        call deallocate (state%m(U_C      )); nullify (state% u_c)
      case ('v_c')
        call deallocate (state%m(V_C      )); nullify (state% v_c)
      case ('z0')
        call deallocate (state%m(Z0       )); nullify (state% z0)
      case ('plcov')
        call deallocate (state%m(PLCOV    )); nullify (state% plcov)
      case ('lai')
        call deallocate (state%m(LAI      )); nullify (state% lai)
      case ('rootdp')
        call deallocate (state%m(ROOTDP   )); nullify (state% rootdp)
      case ('vio3')
        call deallocate (state%m(VIO3     )); nullify (state% vio3)
      case ('hmo3')
        call deallocate (state%m(HMO3     )); nullify (state% hmo3)
      case ('for_e')
        call deallocate (state%m(FOR_E    )); nullify (state% for_e)
      case ('for_d')
        call deallocate (state%m(FOR_D    )); nullify (state% for_d)
      case ('alb_dif')
        call deallocate (state%m(ALB_DIF  )); nullify (state% alb_dif)
      case ('aer_so4')
        call deallocate (state%m(AER_SO4  )); nullify (state% aer_so4)
      case ('aer_dust')
        call deallocate (state%m(AER_DUST )); nullify (state% aer_dust)
      case ('aer_org')
        call deallocate (state%m(AER_ORG  )); nullify (state% aer_org)
      case ('aer_bc')
        call deallocate (state%m(AER_BC   )); nullify (state% aer_bc)
      case ('aer_ss')
        call deallocate (state%m(AER_SS   )); nullify (state% aer_ss)
!     case ('sso_stdh')
!       call deallocate (state%m(SSO_STDH )); nullify (state% sso_stdh)
!     case ('sso_gamma')
!       call deallocate (state%m(SSO_GAMMA)); nullify (state% sso_gamma)
!     case ('sso_theta')
!       call deallocate (state%m(SSO_THETA)); nullify (state% sso_theta)
!     case ('sso_sigma')
!       call deallocate (state%m(SSO_SIGMA)); nullify (state% sso_sigma)
      case ('h_ice')
        call deallocate (state%m(H_ICE    )); nullify (state% h_ice)
      case ('t_ice')
        call deallocate (state%m(T_ICE    )); nullify (state% t_ice)
      case ('h_snow')
        call deallocate (state%m(H_SNOW   )); nullify (state% h_snow)
      case ('tqr')
        call deallocate (state%m(TQR      )); nullify (state% tqr)
      case ('tqs')
        call deallocate (state%m(TQS      )); nullify (state% tqs)
      case ('tqv')
        call deallocate (state%m(TQV      )); nullify (state% tqv)
      case ('den')
        call deallocate (state%m(DEN      )); nullify (state% den)
      case ('vptmp')
        call deallocate (state%m(VPTMP    )); nullify (state% vptmp)
      case ('f_inflat')
        call deallocate (state%m(F_INFLAT )); nullify (state% f_inflat)
      case ('p_energy')
        call deallocate (state%m(P_ENERGY )); nullify (state% p_energy)
      case ('tke')
        call deallocate (state%m(TKE      )); nullify (state% tke )
      case ('dusta')
        call deallocate (state%m(DUSTA    )); nullify (state% dusta)
      case ('dustb')
        call deallocate (state%m(DUSTB    )); nullify (state% dustb)
      case ('dustc')
        call deallocate (state%m(DUSTC    )); nullify (state% dustc)
      case ('qv_dia')
        call deallocate (state%m(QV_DIA   )); nullify (state% qv_dia)
      case ('qc_dia')
        call deallocate (state%m(QC_DIA   )); nullify (state% qc_dia)
      case ('qi_dia')
        call deallocate (state%m(QI_DIA   )); nullify (state% qi_dia)
      case ('o3')
        call deallocate (state%m(O3       )); nullify (state% o3)
      case ('reff_qc')
        call deallocate (state%m(REFF_QC  )); nullify (state% reff_qc)
      case ('reff_qi')
        call deallocate (state%m(REFF_QI  )); nullify (state% reff_qi)
      case ('reff_qr')
        call deallocate (state%m(REFF_QR  )); nullify (state% reff_qr)
      case ('reff_qg')
        call deallocate (state%m(REFF_QG  )); nullify (state% reff_qg)
      case ('reff_qs')
        call deallocate (state%m(REFF_QS  )); nullify (state% reff_qs)
      case ('reff_qh')
        call deallocate (state%m(REFF_QH  )); nullify (state% reff_qh)
      case ('skt')
        call deallocate (state%m(SKT      )); nullify (state% skt)
      case ('co2')
        call deallocate (state%m(CO2      )); nullify (state% co2)
      case ('vis')
        call deallocate (state%m(VIS      )); nullify (state% vis)
      case ('ceiling')
        call deallocate (state%m(CEIL     )); nullify (state% ceiling)
      case ('rcld')
        call deallocate (state%m(RCLD     )); nullify (state% rcld)
      case ('smi')
        call deallocate (state%m(SMI      )); nullify (state% smi)
      case ('t2m_land')
        call deallocate (state%m(T2M_LAND )); nullify (state% t2m_land)
      case ('td2m_land')
        call deallocate (state%m(TD2M_LAND)); nullify (state% td2m_land)
      case ('rh2m_land')
        call deallocate (state%m(RH2M_LAND)); nullify (state% rh2m_land)
      case ('clc_rad')
        call deallocate (state%m(CLC_RAD  )); nullify (state% clc_rad)
      end select
    end do
    state% size = nwords(state% m)
  end subroutine deallocate_state
!------------------------------------------------------------------------------
  subroutine destruct_states (states)
  type (t_atm) ,intent(inout) :: states(:)
    integer i
    do i=1,size(states)
      call destruct (states(i))
    end do
  end subroutine destruct_states
!------------------------------------------------------------------------------
  subroutine destruct_states2 (states)
  type (t_atm) ,intent(inout) :: states(:,:)
    integer i, j
    do j=1,size(states,2)
    do i=1,size(states,1)
      call destruct (states(i,j))
    end do
    end do
  end subroutine destruct_states2
!------------------------------------------------------------------------------
  subroutine destruct_state (state)
  type (t_atm) ,intent(inout) :: state
  !---------------------------------------
  ! deallocate pointer components of t_atm
  !---------------------------------------
    nullify (state% grid)
    state% lb = 1
    state% ub = 0
    call destruct (state% m)
  end subroutine destruct_state
!------------------------------------------------------------------------------
#ifndef TR15581
  elemental subroutine delete (state)
  type (t_atm) ,intent(in) :: state
  !---------------------------------------
  ! deallocate pointer components of t_atm
  !---------------------------------------
    call delete_m (state% m, size(state% m))
  end subroutine delete
#endif
!------------------------------------------------------------------------------
  subroutine copy_head (y, x)
  type (t_atm) ,intent(inout) :: y
  type (t_atm) ,intent(in)    :: x
    y% time     =  x% time
    y% ref_time =  x% ref_time
    y% grid     => x% grid
    y% lb       =  x% lb
    y% ub       =  x% ub
    y% gp_wind  =  x% gp_wind
    y% member   =  x% member
    y% members  =  x% members
  end subroutine copy_head
!------------------------------------------------------------------------------
  subroutine assign_atm_to_atm (y, x)
  type (t_atm) ,intent(inout) :: y
  type (t_atm) ,intent(in)    :: x
    ENTERFUNCTION
    if (any(y%lb/=x%lb.or.y%ub/=x%ub)) call destruct (y%m)
    call copy_head (y, x)
    y% m        =  x% m
    call set_pointers (y)
    DELETESTORAGE (x)
    LEAVEFUNCTION
  end subroutine assign_atm_to_atm
!------------------------------------------------------------------------------
!
!NAG f95 beta (edit 112) bug
!Implicit pointer assignment from X
!illegal in ELEMENTAL procedure ASSIGN_REAL_TO_ATM
!
#ifdef TR15581
#  define ASSIGN(Y,X) Y = X
#else
#  define ASSIGN(Y,X) call assign_m_reals (Y, X)
#endif
  subroutine assign_real_to_atm (y, x)
  type (t_atm) ,intent(inout) :: y
  real(wp)     ,intent(in)    :: x
    ASSIGN(y%m, x)
  end subroutine assign_real_to_atm
!------------------------------------------------------------------------------
  subroutine assign_real_to_atms (y, x)
  type (t_atm) ,intent(inout) :: y(:)
  real(wp)     ,intent(in)    :: x
    integer :: i
    do i=1,size(y)
      ASSIGN(y(i)%m, x)
    end do
  end subroutine assign_real_to_atms
!------------------------------------------------------------------------------
  subroutine assign_real_to_atms2 (y, x)
  type (t_atm) ,intent(inout) :: y(:,:)
  real(wp)     ,intent(in)    :: x
    integer :: i, j
    do j=1,size(y,2)
    do i=1,size(y,1)
      ASSIGN(y(i,j)%m, x)
    end do
    end do
  end subroutine assign_real_to_atms2
!------------------------------------------------------------------------------
  subroutine binop_head (y, x1, x2)
  type (t_atm) ,intent(inout) :: y
  type (t_atm) ,intent(in)    :: x1
  type (t_atm) ,intent(in)    :: x2
    y% time     =  x1% time
    y% ref_time =  x1% ref_time
    y% grid     => x1% grid
    y% lb       =  x1% lb
    y% ub       =  x1% ub
    y% members  = max (x1% members,  x2% members)
    if (x1% member == x2% member) then
      y% member = x1% member
    else if (min(x1% member, x2% member) <= 0) then
      y% member = max (x1% member, x2% member)
    else
      y% member = -1
    endif
  end subroutine binop_head
!------------------------------------------------------------------------------
  function atm_plus_atm (x1, x2) result (y)
  type (t_atm) ,intent(in)    :: x1
  type (t_atm) ,intent(in)    :: x2
  type (t_atm)                :: y
    ENTERFUNCTION
    ALLOCATIONLEVEL (y)
    call binop_head (y, x1, x2)
!    y% m        =  x1% m + x2% m
    call plus (y% m, x1% m, x2% m)
    call set_pointers (y)
    DELETESTORAGE (x1)
    DELETESTORAGE (x2)
    LEAVEFUNCTION
  end function atm_plus_atm
!------------------------------------------------------------------------------
  function atm_minus_atm (x1, x2) result (y)
  type (t_atm) ,intent(in)    :: x1
  type (t_atm) ,intent(in)    :: x2
  type (t_atm)                :: y
    ENTERFUNCTION
    ALLOCATIONLEVEL (y)
    call binop_head (y, x1, x2)
!    y%m         =  x1%m - x2%m
    call minus (y% m, x1% m, x2% m)
    call set_pointers (y)
    DELETESTORAGE (x1)
    DELETESTORAGE (x2)
    LEAVEFUNCTION
  end function atm_minus_atm
!------------------------------------------------------------------------------
  function atm_times_atm (x1, x2) result (y)
  type (t_atm) ,intent(in)    :: x1
  type (t_atm) ,intent(in)    :: x2
  type (t_atm)                :: y
    ENTERFUNCTION
    ALLOCATIONLEVEL (y)
    call binop_head (y, x1, x2)
!    y%m         =  x1%m * x2%m
    call times (y% m, x1% m, x2% m)
    call set_pointers (y)
    DELETESTORAGE (x1)
    DELETESTORAGE (x2)
    LEAVEFUNCTION
  end function atm_times_atm
!------------------------------------------------------------------------------
  function atm_times_r (x1, x2) result (y)
  type (t_atm) ,intent(in)    :: x1
  real(wp)     ,intent(in)    :: x2
  type (t_atm)                :: y
    ENTERFUNCTION
    ALLOCATIONLEVEL (y)
    call copy_head (y, x1)
!    y%m         =  x1%m * x2
    call times (y% m, x1% m, x2)
    call set_pointers (y)
    DELETESTORAGE (x1)
    LEAVEFUNCTION
  end function atm_times_r
!------------------------------------------------------------------------------
  function r_times_atm (x1, x2) result (y)
  real(wp)     ,intent(in)    :: x1
  type (t_atm) ,intent(in)    :: x2
  type (t_atm)                :: y
    ENTERFUNCTION
    ALLOCATIONLEVEL (y)
    call copy_head (y, x2)
!    y%m         =  x2%m * x1
    call times (y% m, x1, x2% m)
    call set_pointers (y)
    DELETESTORAGE (x2)
    LEAVEFUNCTION
  end function r_times_atm
!------------------------------------------------------------------------------
  function i_times_atm (x1, x2) result (y)
  integer      ,intent(in)    :: x1
  type (t_atm) ,intent(in)    :: x2
  type (t_atm)                :: y
    ENTERFUNCTION
    ALLOCATIONLEVEL (y)
    call copy_head (y, x2)
!    y%m         =  x2%m * x1
    call times (y% m, x1, x2% m)
    call set_pointers (y)
    DELETESTORAGE (x2)
    LEAVEFUNCTION
  end function i_times_atm
!------------------------------------------------------------------------------
  function atm_over_atm (x1, x2) result (y)
  type (t_atm) ,intent(in)    :: x1
  type (t_atm) ,intent(in)    :: x2
  type (t_atm)                :: y
    ENTERFUNCTION
    ALLOCATIONLEVEL (y)
    call binop_head (y, x1, x2)
!    y%m         =  x1%m / x2%m
    call over (y% m, x1% m, x2% m)
    call set_pointers (y)
    DELETESTORAGE (x1)
    DELETESTORAGE (x2)
    LEAVEFUNCTION
 end function atm_over_atm
!------------------------------------------------------------------------------
  function atm_over_r (x1, x2) result (y)
  type (t_atm) ,intent(in)    :: x1
  real(wp)     ,intent(in)    :: x2
  type (t_atm)                :: y
    ENTERFUNCTION
    ALLOCATIONLEVEL (y)
    call copy_head (y, x1)
!    y%m         =  x1%m / x2
    call times (y% m, x1% m, 1._wp/x2)
    call set_pointers (y)
    DELETESTORAGE (x1)
    LEAVEFUNCTION
  end function atm_over_r
!------------------------------------------------------------------------------
  function atm_over_i (x1, x2) result (y)
  type (t_atm) ,intent(in)    :: x1
  integer      ,intent(in)    :: x2
  type (t_atm)                :: y
    ENTERFUNCTION
    ALLOCATIONLEVEL (y)
    call copy_head (y, x1)
!    y%m         =  x1%m / x2
    call times (y% m, x1% m, 1._wp/x2)
    call set_pointers (y)
    DELETESTORAGE (x1)
    LEAVEFUNCTION
  end function atm_over_i
!------------------------------------------------------------------------------
  function atm_power_i (x1, x2) result (y)
  type (t_atm) ,intent(in)    :: x1
  integer      ,intent(in)    :: x2
  type (t_atm)                :: y
    ENTERFUNCTION
    ALLOCATIONLEVEL (y)
    call copy_head (y, x1)
!    y%m         =  x1%m ** x2
    call power (y% m, x1% m, x2)
    call set_pointers (y)
    DELETESTORAGE (x1)
    LEAVEFUNCTION
  end function atm_power_i
!------------------------------------------------------------------------------
  function sqrt_atm (x) result (y)
  type (t_atm) ,intent(in)    :: x
  type (t_atm)                :: y
    ENTERFUNCTION
    ALLOCATIONLEVEL (y)
    call copy_head (y, x)
    call sqrt_m (y% m, x% m)
    call set_pointers (y)
    DELETESTORAGE (x)
    LEAVEFUNCTION
  end function sqrt_atm
!------------------------------------------------------------------------------
  function max_r_atm (x1, x2) result (y)
  real(wp)     ,intent(in)    :: x1
  type (t_atm) ,intent(in)    :: x2
  type (t_atm)                :: y
    ENTERFUNCTION
    ALLOCATIONLEVEL (y)
    call copy_head (y, x2)
!    y%m         =  max (x2%m, x1)
    call max_m (y% m, x1, x2% m)
    call set_pointers (y)
    DELETESTORAGE (x2)
    LEAVEFUNCTION
  end function max_r_atm
!------------------------------------------------------------------------------
  function levels (x, g) result (y)
  type (t_atm)  ,intent(in)   :: x
  type (t_grid) ,pointer      :: g
  type (t_atm)                :: y
  !---------------------------------------------------------
  ! select subset of levels for isobaric vertical coordinate
  !---------------------------------------------------------
    integer :: i, j, k
    ENTERFUNCTION
    if (x%grid% vct /= VCT_P_ISO) call finish('levels','x not ISOBARIC')
    if (g%      vct /= VCT_P_ISO) call finish('levels','g not ISOBARIC')
    call construct (y, g, template = x)
    do i = 1, size (x% m)
      if (.not. x% m(i)% i% alloc) cycle
      if (size(y% m(i)% ptr, 3) ==  size(x% m(i)% ptr, 3) ) then
        y% m(i)% ptr = x% m(i)% ptr
      else
        if (size(y% m(i)% ptr, 3) /= y% grid% nz) &
          call finish ('levels','size(y,3)/=nz')
        if (size(x% m(i)% ptr, 3) /= x% grid% nz) &
          call finish ('levels','size(x,3)/=nz')
        k=1
        do j = 1, x% grid% nz
          if (y% grid% akf(k) == x% grid% akf(j)) then
            y% m(i)% ptr(:,:,k,:) = x% m(i)% ptr(:,:,j,:)
            k=k+1
          endif
        enddo
        if (k-1 /= y% grid% nz) call finish ('levels','k-1 /= nz')
      endif
    enddo
    DELETESTORAGE (x)
    LEAVEFUNCTION
  end function levels
!------------------------------------------------------------------------------
  function sum_atm (x) result (y)
  type (t_atm) ,intent(in)    :: x
  real(wp)                    :: y
    y = sum(x%m)
    ENTERFUNCTION
    DELETESTORAGE (x)
    LEAVEFUNCTION
  end function sum_atm
!------------------------------------------------------------------------------
  subroutine p_min_atm (x)
  type (t_atm) ,intent(inout) :: x
    call p_min_m (x%m)
  end subroutine p_min_atm
!------------------------------------------------------------------------------
  subroutine p_max_atm (x)
  type (t_atm) ,intent(inout) :: x
    call p_max_m (x%m)
  end subroutine p_max_atm
!------------------------------------------------------------------------------
  subroutine p_sum_atm (x)
  type (t_atm) ,intent(inout) :: x
    call p_sum_m (x%m)
  end subroutine p_sum_atm
!------------------------------------------------------------------------------
  subroutine p_sum_atm1 (x)
  type (t_atm) ,intent(inout) :: x (:)
    integer :: i
    do i = 1, size (x)
      call p_sum_m (x(i)%m)
    end do
  end subroutine p_sum_atm1
!------------------------------------------------------------------------------
  subroutine p_sum_atm2 (x)
  type (t_atm) ,intent(inout) :: x (:,:)
    integer :: i, j
    do j = 1, size (x,2)
    do i = 1, size (x,1)
      call p_sum_m (x(i,j)%m)
    end do
    end do
  end subroutine p_sum_atm2
!------------------------------------------------------------------------------
  subroutine lin_int_2_atm (y, weight, x1, x2)
  type (t_atm) ,intent(inout) :: y
  real(wp)     ,intent(in)    :: weight
  type (t_atm) ,intent(in)    :: x1, x2
  !---------------------------------------------------------
  ! linear interpolation: y = weight * x1 + (1.-weight) * x2
  !---------------------------------------------------------
    integer  :: i
    ENTERFUNCTION
! SGI BUG WORKAROUND
!   if (any(x1%m%i /= x2%m%i)) call finish('lin_int_2_atm')
    do i=1,size(x1%m)
      if (x1%m(i)%i /= x2%m(i)%i) call finish('lin_int_2_atm')
    end do
! END SGI BUG WORKAROUND
    call binop_head (y, x1, x2)
    y% time     =  weight * x1% time + (1._wp-weight) * x2% time
    y% ref_time =  x1% ref_time
    call update (y% m)
    call set_pointers (y)
    do i=1,size(y% m)
      if (y%m(i)%i% ref) then
        y = x1
      else if (y%m(i)%i% alloc) then
        y%m(i)% ptr = weight * x1%m(i)% ptr + (1._wp-weight) * x2%m(i)% ptr
      endif
    end do
    DELETESTORAGE (x1)
    DELETESTORAGE (x2)
    LEAVEFUNCTION
  end subroutine lin_int_2_atm

!==============================================================================
  subroutine mean_stdev (mean, stdev, x, field)
  !-------------------------------------------
  ! calculate mean & std.deviation of fields
  ! optimised for SX9/LETKF/COSMO (nproc2 = 1)
  !-------------------------------------------
  type (t_atm) ,intent(inout) :: mean
  type (t_atm) ,intent(inout) :: stdev
  type (t_atm) ,intent(in)    :: x (:)
  character(*), intent(in)    :: field ! (Re)calculate mean/stdev
  optional                    :: field
    integer :: i   ! 1st horizontal index (~ 1..15)
    integer :: j   ! 2nd horizontal index (~ 1..500)
    integer :: k   ! vertical index       (~ 1..50)
    integer :: d   ! diamond index        (~ 1)
    integer :: m   ! model field index
    integer :: l   ! ensemble index       (~ 1..40)
    integer :: n   ! ensemble size        (~ 40)
    integer :: lb(4), ub(4)             ! Copy of grid local bounds
    real(wp), allocatable :: sum_(:)    ! Temporary sum
    real(wp), allocatable :: var_(:)    ! Temporary sum of squares
    ENTERFUNCTION

    n  = size(x)
    lb = mean% grid% lb
    ub = mean% grid% ub

    allocate (sum_(lb(1):ub(1)))
    allocate (var_(lb(1):ub(1)))

!#ifndef _CRAYFTN /* work around bug in Crayftn <= 8.3.5 */
!$omp parallel do private(i,j,k,d,l,m,sum_,var_) schedule(dynamic)
!#endif
!NEC$ novector
    do m=1, size(mean% m)
      if (       mean% m(m)% i% alloc) then
        if (present (field)) then
          if (field /= mean% m(m)% i% name) cycle
          if (.not. x(1)% m(m)% i% alloc) &
            call finish ("mean_stdev","not allocated: " // trim (field))
        end if
        do d=              lb(4), ub(4)
        do k=mean% m(m)%i% lb(3), mean% m(m)%i% ub(3)
        do j=              lb(2), ub(2)

        sum_ = 0._wp
        var_ = 0._wp
!DIR$ IVDEP
!NEC$ ivdep
        do i=              lb(1), ub(1)
!NEC$ nounroll
          do l = 1, n
            sum_(i) = sum_(i) +  x(l)% m(m)% ptr(i,j,k,d)
          end do
          !-------------------------------------------------------
          ! Compute mean first to avoid rounding issues with
          ! ("invariant") fields that are the same for all members
          !-------------------------------------------------------
          sum_(i) = sum_(i) * (1._wp/n)
!NEC$ nounroll
          do l = 1, n
            var_(i) = var_(i) + (x(l)% m(m)% ptr(i,j,k,d) - sum_(i)) ** 2
          end do
          mean % m(m)% ptr(i,j,k,d) =       sum_(i)
          stdev% m(m)% ptr(i,j,k,d) = sqrt (var_(i) * (1._wp/(n-1)))
        end do
        end do
        end do
        end do
      endif
    end do
!#ifndef _CRAYFTN /* work around bug in Crayftn <= 8.3.5 */
!$omp end parallel do
!#endif

    do l=1, n
      DELETESTORAGE (x(l))
    end do
    LEAVEFUNCTION

  end subroutine mean_stdev
!==============================================================================
  subroutine ens_corr (corr, mean, stdev, x, vnam)
  !-------------------------------------------
  ! calculate ensemble correlation of fields
  !-------------------------------------------
  type (t_atm)     ,intent(inout) :: corr  ! -> derived correlation
  type (t_atm)     ,intent(inout) :: mean  ! <- mean
  type (t_atm)     ,intent(inout) :: stdev ! <- standard deviation
  type (t_atm)     ,intent(in)    :: x (:) ! <- ensemble state
  character(len=*) ,intent(in)    :: vnam  ! variable to be correlated with others

    integer :: my, ky ! indices for 'vnam'
    integer :: i   ! 1st horizontal index (~ 1..15)
    integer :: j   ! 2nd horizontal index (~ 1..500)
    integer :: k   ! vertical index       (~ 1..50)
    integer :: d   ! diamond index        (~ 1)
    integer :: m   ! model field index
    integer :: l   ! ensemble index       (~ 1..40)
    integer :: n   ! ensemble size        (~ 40)

    do m=1, size(corr% m)
      if ( corr% m(m)% i% alloc ) then
        if ( trim(corr% m(m)% i% name) == trim(vnam) ) then
          my = m
          ky = corr% m(m)%i% lb(3)
          ! print *, 'vnam% idx = ', vnam, "%", idx
        endif
      endif
    enddo

    ENTERFUNCTION
    n = size(x)
!NEC$ novector
    do m=1, size(corr% m)
     if ( corr% m(m)% i% alloc) then
        corr% m(m)% ptr(:,:,:,:) = 0._wp
        do d=corr% grid%   lb(4), corr% grid%   ub(4)
          do k=corr% m(m)%i% lb(3), corr% m(m)%i% ub(3)
            do j=corr% grid%   lb(2), corr% grid%   ub(2)
              do i=corr% grid%   lb(1), corr% grid%   ub(1)
!NEC$ nounroll
                 do l = 1, n
                    corr% m(m)% ptr(i,j,k,d) =  corr% m(m )% ptr(i,j,k ,d)                              &
                                            + ((x(l)% m(m )% ptr(i,j,k ,d) - mean% m(m )% ptr(i,j,k ,d))&
                                            *  (x(l)% m(my)% ptr(i,j,ky,d) - mean% m(my)% ptr(i,j,ky,d)))
                 end do
                 corr% m(m)% ptr(i,j,k,d) = corr% m(m)% ptr(i,j,k,d) * (1._wp/(n-1))
                 if ((stdev% m(m)% ptr(i,j,k,d) * stdev% m(my)% ptr(i,j,ky,d)) == 0._wp ) then
                     corr% m(m)% ptr(i,j,k,d) = 0._wp
                 else
                     corr% m(m)% ptr(i,j,k,d) = corr% m(m )% ptr(i,j,k ,d) &
                                            / (stdev% m(m )% ptr(i,j,k ,d) &
                                            *  stdev% m(my)% ptr(i,j,ky,d) )
                 endif
!NEC$ nounroll
              end do
            end do
          end do
        end do
     endif
    end do

    do l=1, n
      DELETESTORAGE (x(l))
    end do
    DELETESTORAGE (mean)
    DELETESTORAGE (stdev)
    LEAVEFUNCTION

  end subroutine ens_corr
!==============================================================================
  subroutine lgt_atm (state, nn, norm, only)
  type (t_atm)               ,intent(inout) :: state
  integer          ,optional ,intent(in)    :: nn
  integer          ,optional ,intent(in)    :: norm
  character(len=*) ,optional ,intent(in)    :: only
    integer           :: i
    character(len=64) :: onl
    if (present(only)) then
      onl = only
      do i=1,size(state%m)
        if (index(onl,trim(state%m(i)%i%name)//' ')/=0) &
          call lgtd (state%m(i), nn=nn, norm=norm)
      end do
    else
      do i=1,size(state%m)
        call lgtd (state%m(i), nn=nn, norm=norm)
      end do
    endif
    call set_pointers (state)
  end subroutine lgt_atm
!------------------------------------------------------------------------------
  subroutine lgti_atm (state, ni, nj, norm, only)
  type (t_atm)               ,intent(inout) :: state
  integer          ,optional ,intent(in)    :: ni
  integer          ,optional ,intent(in)    :: nj
  integer          ,optional ,intent(in)    :: norm
  character(len=*) ,optional ,intent(in)    :: only
    integer           :: i
    character(len=64) :: onl
    if (present(only)) then
      onl = only
      do i=1,size(state%m)
        if (index(onl,trim(state%m(i)%i%name)//' ')/=0) then
          if(present(ni)) then
            state%m(i)%i% lb(1) = 1
            state%m(i)%i% ub(1) = ni
          endif
          if(present(nj)) then
            state%m(i)%i% lb(2) = 1
            state%m(i)%i% ub(2) = nj
          endif
          call lgti (state%m(i), norm=norm)
        endif
      end do
    else
      if(present(ni)) then
        state%m%i% lb(1) = 1
        state%m%i% ub(1) = ni
      endif
      if(present(nj)) then
        state%m%i% lb(2) = 1
        state%m%i% ub(2) = nj
      endif
      do i=1,size(state%m)
        call lgti (state%m(i), norm=norm)
      end do
    endif
    call set_pointers (state)
  end subroutine lgti_atm
!------------------------------------------------------------------------------
  subroutine lgti_ad_atm (state, nn, norm, only)
  type (t_atm)               ,intent(inout) :: state
  integer          ,optional ,intent(in)    :: nn
  integer          ,optional ,intent(in)    :: norm
  character(len=*) ,optional ,intent(in)    :: only
    integer           :: i
    character(len=64) :: onl
    if (present(only)) then
      onl = only
      do i=1,size(state%m)
        if (index(onl,trim(state%m(i)%i%name)//' ')/=0) &
          call lgti_ad (state%m(i), nn=nn, norm=norm)
      end do
    else
      do i=1,size(state%m)
        call lgti_ad (state%m(i), nn=nn, norm=norm)
      end do
    endif
    call set_pointers (state)
  end subroutine lgti_ad_atm
!------------------------------------------------------------------------------
  subroutine lgt_ad_atm (state, norm, only)
  type (t_atm)               ,intent(inout) :: state
  integer          ,optional ,intent(in)    :: norm
  character(len=*) ,optional ,intent(in)    :: only
    integer           :: i
    character(len=64) :: onl
    if (present(only)) then
      onl = only
      do i=1,size(state%m)
        if (index(onl,trim(state%m(i)%i%name)//' ')/=0) &
          call lgtd_ad (state%m(i), norm=norm)
      end do
    else
      do i=1,size(state%m)
        call lgtd_ad (state%m(i), norm=norm)
      end do
    endif
    call set_pointers (state)
  end subroutine lgt_ad_atm
!==============================================================================
  elemental subroutine set_atm_pointers (state)
  type (t_atm) TARGET,intent(inout) :: state
  !---------------------------------
  ! set pointers to allocated memory
  !---------------------------------
      state% ps        => state%m(PS       )% ptr
      state% psr       => state%m(PSR      )% ptr
      state% pp        => state%m(PP       )% ptr
      state% dpsdt     => state%m(DPSDT    )% ptr
      state% t         => state%m(T        )% ptr
      state% u         => state%m(U        )% ptr
      state% v         => state%m(V        )% ptr
      state% w         => state%m(W        )% ptr
      state% q         => state%m(Q        )% ptr
      state% qcl       => state%m(QCL      )% ptr
      state% qci       => state%m(QCI      )% ptr
      state% qr        => state%m(QR       )% ptr
      state% qs        => state%m(QS       )% ptr
      state% qg        => state%m(QG       )% ptr
      state% qh        => state%m(QH       )% ptr
      state% nccloud   => state%m(NCCLOUD  )% ptr
      state% ncrain    => state%m(NCRAIN   )% ptr
      state% ncsnow    => state%m(NCSNOW   )% ptr
      state% ncice     => state%m(NCICE    )% ptr
      state% ncgraupel => state%m(NCGRAUPEL)% ptr
      state% nchail    => state%m(NCHAIL   )% ptr
      state% t_so      => state%m(T_SO     )% ptr
      state% w_so      => state%m(W_SO     )% ptr
      state% w_so_ice  => state%m(W_SO_ICE )% ptr
      state% t_s       => state%m(T_S      )% ptr
      state% qv_s      => state%m(QV_S     )% ptr
      state% t_snow    => state%m(T_SNOW   )% ptr
      state% freshsnw  => state%m(FRESHSNW )% ptr
      state% rho_snow  => state%m(RHO_SNOW )% ptr
      state% w_snow    => state%m(W_SNOW   )% ptr
      state% snowc     => state%m(SNOWC    )% ptr
      state% w_i       => state%m(W_I      )% ptr
      state% iremispc  => state%m(IREMISPC )% ptr
      state% t_mnw_lk  => state%m(T_MNW_LK )% ptr
      state% t_wml_lk  => state%m(T_WML_LK )% ptr
      state% h_ml_lk   => state%m(H_ML_LK  )% ptr
      state% t_bot_lk  => state%m(T_BOT_LK )% ptr
      state% c_t_lk    => state%m(C_T_LK   )% ptr
      state% clc       => state%m(CLC      )% ptr
      state% ph        => state%m(PH       )% ptr
      state% pf        => state%m(PF       )% ptr
      state% rh        => state%m(RH       )% ptr
      state% t2m       => state%m(T2M      )% ptr
      state% rh2m      => state%m(RH2M     )% ptr
      state% geoh      => state%m(GEOH     )% ptr
      state% geof      => state%m(GEOF     )% ptr
      state% tsurf     => state%m(TSURF    )% ptr
      state% fr_ice    => state%m(FR_ICE   )% ptr
      state% td2m      => state%m(TD2M     )% ptr
      state% tv        => state%m(TV       )% ptr
      state% psi       => state%m(PSI      )% ptr
      state% chi       => state%m(CHI      )% ptr
      state% vrt       => state%m(VRT      )% ptr
      state% div       => state%m(DIV      )% ptr
      state% vmax_10m  => state%m(VMAX_10M )% ptr
      state% tmin_2m   => state%m(TMIN_2M  )% ptr
      state% tmax_2m   => state%m(TMAX_2M  )% ptr
      state% tot_prec  => state%m(TOT_PREC )% ptr
      state% aswdir_s  => state%m(ASWDIR_S )% ptr
      state% aswdifd_s => state%m(ASWDIFD_S)% ptr
      state% alwd_s    => state%m(ALWD_S   )% ptr
      state% clct      => state%m(CLCT     )% ptr
      state% clcl      => state%m(CLCL     )% ptr
      state% clcm      => state%m(CLCM     )% ptr
      state% clch      => state%m(CLCH     )% ptr
      state% u_10m     => state%m(U_10M    )% ptr
      state% v_10m     => state%m(V_10M    )% ptr
      state% pmsl      => state%m(PMSL     )% ptr
      state% u_10m_c   => state%m(U_10M_C  )% ptr
      state% v_10m_c   => state%m(V_10M_C  )% ptr
      state% u_c       => state%m(U_C      )% ptr
      state% v_c       => state%m(V_C      )% ptr
      state% z0        => state%m(Z0       )% ptr
      state% plcov     => state%m(PLCOV    )% ptr
      state% lai       => state%m(LAI      )% ptr
      state% rootdp    => state%m(ROOTDP   )% ptr
      state% vio3      => state%m(VIO3     )% ptr
      state% hmo3      => state%m(HMO3     )% ptr
      state% for_e     => state%m(FOR_E    )% ptr
      state% for_d     => state%m(FOR_D    )% ptr
      state% alb_dif   => state%m(ALB_DIF  )% ptr
      state% aer_so4   => state%m(AER_SO4  )% ptr
      state% aer_dust  => state%m(AER_DUST )% ptr
      state% aer_org   => state%m(AER_ORG  )% ptr
      state% aer_bc    => state%m(AER_BC   )% ptr
      state% aer_ss    => state%m(AER_SS   )% ptr
!     state% sso_stdh  => state%m(SSO_STDH )% ptr
!     state% sso_gamma => state%m(SSO_GAMMA)% ptr
!     state% sso_theta => state%m(SSO_THETA)% ptr
!     state% sso_sigma => state%m(SSO_SIGMA)% ptr
      state% h_ice     => state%m(H_ICE    )% ptr
      state% t_ice     => state%m(T_ICE    )% ptr
      state% h_snow    => state%m(H_SNOW   )% ptr
      state% tqr       => state%m(TQR      )% ptr
      state% tqs       => state%m(TQS      )% ptr
      state% tqv       => state%m(TQV      )% ptr
      state% den       => state%m(DEN      )% ptr
      state% vptmp     => state%m(VPTMP    )% ptr
      state% f_inflat  => state%m(F_INFLAT )% ptr
      state% p_energy  => state%m(P_ENERGY )% ptr
      state% tke       => state%m(TKE      )% ptr
      state% dusta     => state%m(DUSTA    )% ptr
      state% dustb     => state%m(DUSTB    )% ptr
      state% dustc     => state%m(DUSTC    )% ptr
      state% qv_dia    => state%m(QV_DIA   )% ptr
      state% qc_dia    => state%m(QC_DIA   )% ptr
      state% qi_dia    => state%m(QI_DIA   )% ptr
      state% o3        => state%m(O3       )% ptr
      state% reff_qc   => state%m(REFF_QC  )% ptr
      state% reff_qi   => state%m(REFF_QI  )% ptr
      state% reff_qr   => state%m(REFF_QR  )% ptr
      state% reff_qg   => state%m(REFF_QG  )% ptr
      state% reff_qs   => state%m(REFF_QS  )% ptr
      state% reff_qh   => state%m(REFF_QH  )% ptr
      state% skt       => state%m(SKT      )% ptr
      state% co2       => state%m(CO2      )% ptr
      state% vis       => state%m(VIS      )% ptr
      state% ceiling   => state%m(CEIL     )% ptr
      state% rcld      => state%m(RCLD     )% ptr
      state% smi       => state%m(SMI      )% ptr
      state% t2m_land  => state%m(T2M_LAND )% ptr
      state% td2m_land => state%m(TD2M_LAND)% ptr
      state% rh2m_land => state%m(RH2M_LAND)% ptr
      state% clc_rad   => state%m(CLC_RAD  )% ptr
      state% size      =  nwords(state%m)
  end subroutine set_atm_pointers
!------------------------------------------------------------------------------
  subroutine set_ps (state)
  !-----------------------------------------------
  ! set surface pressure (ps).
  ! for the COSMO hybrid z-coordinate derive ps
  !   from pp (pressure deviation at full levels);
  ! for the ICON hybrid z-coordinate
  !   extrapolate ps from lowest full level.
  !-----------------------------------------------
  type (t_atm) ,intent(inout) :: state
    integer               :: ie, je, ke
    real(wp) ,allocatable :: qrs(:,:)   ! precipitation water (rain,snow,...)
    real(wp) ,allocatable :: qcl(:,:)   ! specific cloud water content
    logical               :: ltv
    integer               :: lb1, ub1
    type(t_grid) ,pointer :: grid
    integer               :: i, j, k, d ! grid/loop indices
    real(wp)              :: fis        ! surface geopotential
    real(wp)              :: fi0, p0    ! reference values for interpolation

    !----------------------------------------------
    ! Generalised vertical coordinate (ICON, COSMO)
    !----------------------------------------------
    if ( state% grid% model == MO_ICON  .or.                              &
        (state% grid% model == MO_COSMO .and. .not.associated (state% pp))) then
      if (.not.associated (state% grid% hhl)) &
        call finish('set_ps','hhl not present')
      if (.not.associated (state% pf)) &
        call finish('set_ps','pf not present')
      !----------------------------------------------
      ! Extrapolate from lowest full level to surface
      !----------------------------------------------
      ke  = size (state% pf,3)
      ltv = .not.associated (state% tv)
      if (ltv) then
        if(associated(state% q)) then
          call set_tv (state)
        else
          call allocate (state, 'tv')
          state% tv(:,:,ke,:) = state% t(:,:,ke,:)
        endif
      endif
      call allocate (state, 'ps')
      lb1 = state% grid% lb(1)
      ub1 = state% grid% ub(1)
      state% ps(:,:,1,:) &
        = state% pf(:,:,ke,:) *                                     &
          exp (gacc/R * 0.5_wp*(state% grid% hhl(lb1:ub1,:,ke  ,:)- &
          state% grid% hhl(lb1:ub1,:,ke+1,:))                       &
          / state% tv       (   :   ,:,ke  ,:)  )
      if (ltv) call deallocate (state, 'tv')
    !-----------------------
    ! COSMO old z-coordinate
    !-----------------------
    else if (state% grid% model == MO_COSMO) then
      !--------------------------------------
      ! check for presence of required fields
      !--------------------------------------
      grid => state% grid
      if (.not.associated(state% pp )) call finish('set_ps','pp   not present')
      if (.not.associated(state% t  )) call finish('set_ps','t    not present')
      if (.not.associated(state% q  )) call finish('set_ps','q    not present')
      if (.not.associated(grid% rho0)) call finish('set_ps','rho0 not present')
      if (.not.associated(grid% p0  )) call finish('set_ps','p0   not present')
      if (.not.associated(grid% dp0 )) call finish('set_ps','dp0  not present')
      !------------------
      ! derive ps from pf
      !------------------
      ie = size(state% pp,1)
      je = size(state% pp,2)
      ke = size(state% pp,3)
      allocate (qrs (ie,je)) ;qrs = 0._wp
      allocate (qcl (ie,je)) ;qcl = 0._wp
      if (associated(state% qr )) qrs = qrs + state% qr (:,:,ke,1) ! rain
      if (associated(state% qs )) qrs = qrs + state% qs (:,:,ke,1) ! snow
      if (associated(state% qg )) qrs = qrs + state% qg (:,:,ke,1) ! graupel
      if (associated(state% qh )) qrs = qrs + state% qh (:,:,ke,1) ! hail
      if (associated(state% qci)) qrs = qrs + state% qci(:,:,ke,1) ! ice
      if (associated(state% qcl)) qcl = qcl + state% qcl(:,:,ke,1) ! cloud water
      call allocate (state, 'ps')
      call calps (state% ps  (:,:,1 ,1) ,&! surface pressure
                  state% pp  (:,:,ke,1) ,&! perturbation pressure
                  state% t   (:,:,ke,1) ,&! temperature
                  state% q   (:,:,ke,1) ,&! specific water vapor content
                         qcl            ,&! specific cloud water content
                         qrs            ,&! precipitation water
                  grid%  rho0(:,:,ke,1) ,&! reference dens., full levels
                  grid%  p0  (:,:,ke,1) ,&! reference press.,full levels
                  grid%  dp0 (:,:,ke,1) ,&! fulllevel pressure thickness
                  ie, je                ,&! array sizes
                  rddrm1                ,&! r_v/r_d - 1
                  R                     ,&! gas constant for dry air
                  1, ie, 1, je          ) ! index range for computation
    !----------------------------
    ! IFS/GME hybrid p-coordinate
    !----------------------------
    else if ((state% grid% model == MO_IFS    .or. &
              state% grid% model == MO_HRM    .or. &
              state% grid% model == MO_GME         ) .and. &
             (state% grid% vct   == VCT_P_HYB .or.  &
              state% grid% vct   == VCT_P_IFS      )       ) then
      if (.not.associated (state% ps)) &
        call finish('set_ps','cannot derive ps for IFS/GME')
    !-----------------------------------------------------
    ! Isobaric grid ("other models"): derive ps from pmsl.
    !-----------------------------------------------------
    else if (state% grid% vct == VCT_P_ISO) then
      grid => state% grid
      if (.not. associated (grid% geosp)) &
           call finish ('set_ps','geosp not associated')
      if (.not. associated (state% geof)) &
           call finish ('set_ps','geof not present')
      if (.not. associated (state% pmsl)) &
           call finish ('set_ps','pmsl not present')

      ke = size (state% geof,3)
      if (ke < 2) call finish ('set_ps','geof must have at least two levels')
      call allocate (state, 'ps')

      do     d = grid% lb(4), grid% ub(4)
        do   j = grid% lb(2), grid% ub(2)
          do i = grid% lb(1), grid% ub(1)
             fis = grid% geosp(i,j,1,d)
             if (fis >= state% geof(i,j,ke,d)) then
                !------------------------------------------------------------
                ! Interpolate where model surface above lowest isobaric level
                !------------------------------------------------------------
                do k = 2, ke
                   if (fis >= state% geof(i,j,k,d)) exit
                end do
                k   = k - 1
                fi0 = state% geof(i,j,k+1,d)
                p0  = grid%  akf (k+1)
             else
                !------------------------------
                ! Extrapolate to mean sea level
                !------------------------------
                k   = ke
                fi0 = 0._wp
                p0  = state% pmsl(i,j,1  ,d)
                !-------------------------------------------------------
                ! Beware: use tolerance due to (grib) rounding;
                ! next-to-lowest level sometimes numerically more stable
                !-------------------------------------------------------
                if (abs (state% geof(i,j,ke,d)) < 1.e-1_wp) then
                   k = ke-1
                end if
             end if
             state% ps(i,j,1,d) = p0 * exp (log (grid%  akf(k) / p0)         &
                                              * (       fis           - fi0) &
                                              / (state% geof(i,j,k,d) - fi0) )
          end do
        end do
      end do
    !-----------------------------------
    ! no generalised vertical coordinate
    !-----------------------------------
    else
      write (0,*) "set_ps: model,vct =", state% grid% model,  state% grid% vct
      call finish('set_ps','no COSMO or ICON z-coordinate')
    end if

  end subroutine set_ps
!------------------------------------------------------------------------------
  subroutine set_p1 (states)
  !----------------------
  ! set pressure level(s)
  !----------------------
  type (t_atm) ,intent(inout) :: states (:)
    integer :: i
    do i = 1, size (states)
      call set_p (states (i))
    end do
  end subroutine set_p1
!------------------------------------------------------------------------------
  subroutine set_ps1 (states)
  !----------------------
  ! set pressure level(s)
  !----------------------
  type (t_atm) ,intent(inout) :: states (:)
    integer :: i
    do i = 1, size (states)
      call set_ps (states (i))
    end do
  end subroutine set_ps1
!------------------------------------------------------------------------------
  subroutine set_pf1 (states)
  !----------------------
  ! set pressure level(s)
  !----------------------
  type (t_atm) ,intent(inout) :: states (:)
    integer :: i
    do i = 1, size (states)
      call set_pf (states (i))
    end do
  end subroutine set_pf1
!------------------------------------------------------------------------------
  subroutine set_ph1 (states)
  !----------------------
  ! set pressure level(s)
  !----------------------
  type (t_atm) ,intent(inout) :: states (:)
    integer :: i
    do i = 1, size (states)
      call set_ph (states (i))
    end do
  end subroutine set_ph1
!------------------------------------------------------------------------------
  subroutine set_p (state)
  !----------------------
  ! set pressure level(s)
  !----------------------
  type (t_atm) ,intent(inout) :: state
    select case (state% grid% vct)
    case (VCT_P_ISO)
      !--------------
      ! isobaric grid
      !--------------
      call set_pf (state)
    case (VCT_Z_HYB, VCT_Z_GEN)
      select case (state% grid% model)
      case (MO_ICON, MO_COSMO)
      !------------------------------------------------------------------
      ! COSMO: derive pressure from pressure deviation pp if present
      ! ICON: derive pressure from density, virtual potential temperature
      !------------------------------------------------------------------
        if (.not. associated (state% pf)) call set_pf (state)
        if (.not. associated (state% ps)) call set_ps (state)
        if (.not. associated (state% ph)) call set_ph (state)
      case default
        if (associated (state% pf)) return
        call finish('set_p','invalid model')
      end select
    case (VCT_P_HYB)
      !----------
      ! GME / HRM
      !----------
      call set_ph (state)
      call set_pf (state)
    case default
      call finish('set_p','invalid vct')
    end select
  end subroutine set_p
!------------------------------------------------------------------------------
  subroutine set_ph (state)
  !----------------------
  ! set pressure level(s)
  !----------------------
  type (t_atm) ,intent(inout) :: state
    integer   :: k, ke
    character :: cvct
    ke = state% ub(3)
    select case (state% grid% vct)
    case default
      write (cvct,'(i1)') state% grid% vct
      call finish ('set_ph','unknown vct '//cvct)
    case (VCT_P_HYB)
      !--------------------------------------------
      ! hybrid pressure coordinates: set ph from ak ,bk ,psr
      !                              set pf from akf,bkf,psr
      !--------------------------------------------
      call allocate (state,'ph')
      if (associated(state% psr)) then
        do k = 1, ke + 1
          state% ph(:,:,k,:) = state%grid% ak(k) + state%grid% bk(k) &
                             * state% psr(:,:,1,:)
        end do
      else
        if (.not.associated(state% ps)) &
          call finish('set_ph','neither ps nor psr is allocated')
        do k = 1, ke + 1
          state% ph(:,:,k,:) = state%grid% ak(k) + state%grid% bk(k) &
                             * state% ps(:,:,1,:)
        end do
      endif
    case (VCT_Z_HYB, VCT_Z_GEN)
      !--------------------------
      ! not implemented for COSMO
      !--------------------------
!     if (state% grid% model == MO_COSMO) &
!       call finish ('set_ph','model == COSMO')
      !------------------------
      ! ICON (preliminary):
      ! ph(k)*ph(k+1) = pf(k)^2
      !------------------------
      if (.not.associated (state% ps)) call finish('set_ph',          &
                                                   'ps not associated')
      if (.not.associated (state% pf)) call finish('set_ph',          &
                                                   'pf not associated')
      call allocate (state,'ph')
      state% ph(:,:,ke+1,:) = state% ps(:,:,1,:)
      do k = ke, 1, -1
         state% ph(:,:,k,:) = state% pf(:,:,k,:)**2 / state% ph(:,:,k+1,:)
      end do
    end select
  end subroutine set_ph
!------------------------------------------------------------------------------
  subroutine set_pf (state)
  !----------------------
  ! set pressure level(s)
  !----------------------
  type (t_atm) ,intent(inout) :: state

    integer :: k, ke

    ke = state% ub(3)
    call allocate (state,'pf')
    select case (state% grid% vct)
    case default
      call finish ('set_pf','unknown vct')
    case (VCT_P_ISO)
      !--------------------------------------
      ! isobaric coordinates: set pf from akf
      !--------------------------------------
      do k = 1, ke
        state% pf(:,:,k,:) = state%grid% akf(k)
      end do
    case (VCT_P_HYB)
      !-----------------------------------------------------
      ! hybrid pressure coordinates: set ph from ak ,bk ,psr
      !                     set pf from akf,bkf,psr
      !-----------------------------------------------------
      if (associated(state% psr)) then
        do k = 1, ke
          state% pf(:,:,k,:) = state%grid% akf(k) + state%grid% bkf(k) &
                             * state% psr(:,:,1,:)
        end do
      else
        if (.not.associated(state% ps)) &
          call finish ('set_pf','neither ps nor psr is allocated')
        do k = 1, ke
          state% pf(:,:,k,:) = state%grid% akf(k) + state%grid% bkf(k) &
                             * state% ps(:,:,1,:)
        end do
      endif
    case (VCT_Z_HYB, VCT_Z_GEN)
      select case (state% grid% model)
      !----------------------------------------------------------
      ! COSMO: derive from pp and reference atmosphere if present
      !----------------------------------------------------------
      case (MO_COSMO)
        if (associated (state% pp) .and. associated (state% grid% p0)) then
          state% pf = state% pp + state% grid% p0
          return
        end if
      end select
      !--------------------------------------------
      ! ICON, COSMO: either derive from half levels
      ! pf(k)^2 = ph(k)*ph(k+1)
      !--------------------------------------------
      if (associated (state% ph)) then
        do k = 1, ke
          state% pf(:,:,k,:) = sqrt (state% ph(:,:,k+1,:)*state% ph(:,:,k,:))
        end do
      !---------------------------------------------------------
      ! or derive from density and virtual potential temperature
      !---------------------------------------------------------
      else if (associated (state% den) .and. associated (state% vptmp)) then
        state% pf = p_rho_thetav (rho=state% den, thetav=state% vptmp)
      else
        call finish('set_pf','neither ph nor (den,vptmp) is associated')
      end if
    end select
  end subroutine set_pf
!------------------------------------------------------------------------------
  subroutine set_p_adj (sta_a, state)
  type (t_atm) ,intent(inout) :: sta_a
  type (t_atm) ,intent(in)    :: state
!---------------
! nonlinear code
!---------------
!   integer :: k, kep1
!   if (.not. associated (state% ps)) call finish ('set_p','ps not associated')
!   call allocate (state,'ph')
!   call allocate (state,'pf')
!   kep1 = state% ub(3) + 1
!   do k = 1,kep1
!     state% ph(:,:,k) = state%grid% ak(k) + state%grid% bk(k) * state% ps
!   end do
!   state% pf = 0.5_wp * (state% ph(:,:,:kep1-1) + state% ph(:,:,2:))
!--------------------
! tangent linear code
!--------------------
!   do k = 1,kep1
!     sta_a% ph(:,:,k) = state%grid% bk(k) * sta_a% ps
!   end do
!   sta_a% pf = 0.5_wp * (sta_a% ph(:,:,:kep1-1) + sta_a% ph(:,:,2:))
!-------------
! adjoint code
!-------------
    integer :: k, ke
    ke = state% ub(3)
    select case (state% grid% vct)
    case default
      !--------------------------------
      ! not implemented for COSMO, ICON
      !--------------------------------
      call finish ('set_p_adj','unknown vct')
    case (VCT_P_ISO)
      !--------------------------------------
      ! isobaric coordinates: set pf from akf
      !--------------------------------------
      sta_a% pf(:,:,:,:)     = 0._wp
    case (VCT_P_HYB)
      !-----------------------------------------------------
      ! hybrid pressure coordinates: set ph from ak ,bk ,psr
      !                     set pf from akf,bkf,psr
      !-----------------------------------------------------
      do k = 1,ke
        sta_a% psr(:,:,1,:)  = sta_a% psr(:,:,1,:) &
                             + sta_a% pf (:,:,k,:) * state%grid% bkf(k)
        sta_a% pf(:,:,k,:)   = 0._wp
      end do
      do k = 1,ke+1
        sta_a% psr(:,:,1,:)  = sta_a% psr(:,:,1,:) &
                             + sta_a% ph (:,:,k,:) * state%grid% bk(k)
        sta_a% ph(:,:,k,:)   = 0._wp
      end do
    end select
  end subroutine set_p_adj
!------------------------------------------------------------------------------
  subroutine set_geo1 (states, geof)
  !------------------------
  ! set geopotantial height
  !------------------------
  type(t_atm) ,intent(inout)        :: states (:)
  logical     ,intent(in) ,optional :: geof ! flag for geop. on full levels
    integer :: i
    do i = 1, size (states)
      call set_geo (states (i), geof)
    end do
  end subroutine set_geo1
!------------------------------------------------------------------------------
  subroutine set_geo (a, geof)
  !-----------------------
  ! calculate geopotential
  !-----------------------
  type (t_atm) ,intent(inout)        :: a    ! atmospheric state
  logical      ,intent(in) ,optional :: geof ! flag for geop. on full levels
    !-----------------------------------------------------------------------
    ! local variables: indices, virtual temp., cloud ice+water , delta log p
    !-----------------------------------------------------------------------
    integer               :: k, kep
    logical               :: lgeof
    integer               :: lb1, ub1
    real(wp), allocatable :: tv(:,:,:,:), x(:,:,:,:), dlnp(:,:,:,:)

    lgeof = .false.; if(present(geof)) lgeof = geof
    call allocate (a,'geoh')
    kep = a% ub(3)
    select case (a% grid% vct)
    !----------
    ! GME / HRM
    !----------
    case (VCT_P_HYB)
      !-----------------------------------------------------
      ! sanity checks, allocate components in state variable
      !-----------------------------------------------------
      if (.not.associated (a% ph)) call finish('set_geo','ph not associated')
      if (.not.associated (a% t )) call finish('set_geo','t  not associated')
      if (.not.associated (a% q )) call finish('set_geo','q  not associated')
      allocate (tv   (size(a% t,1),size(a% t,2),size(a% t,3),size(a% t,4)))
      allocate (x    (size(a% t,1),size(a% t,2),size(a% t,3),size(a% t,4)))
      allocate (dlnp (size(a% t,1),size(a% t,2),size(a% t,3),size(a% t,4)))
      !---------------------------
      ! derive virtual temperature
      !---------------------------
      x = 0._wp
      if (associated(a% qcl)) x = x + a% qcl
      if (associated(a% qci)) x = x + a% qci
      if (associated(a% qr )) x = x + a% qr
      if (associated(a% qs )) x = x + a% qs
      if (associated(a% qg )) x = x + a% qg
      if (associated(a% qh )) x = x + a% qh
      tv = tv_t_q (a% t, a% q, x)
      !-------------------------------
      ! integrate hydrostatic equation
      !-------------------------------
      a% geoh(:,:,kep+1,:) = a% grid% geosp(a% lb(1):a% ub(1),&
                                            a% lb(2):a% ub(2),&
                                                    1        ,&
                                            a% lb(4):a% ub(4))
      do k = kep, 2, -1
        dlnp   (:,:,k,:) = R * log (a% ph(:,:,k+1,:) / a% ph(:,:,k,:))
        a% geoh(:,:,k,:) = a% geoh(:,:,k+1,:) + dlnp(:,:,k,:) * tv(:,:,k,:)
      end do
      !----------
      ! top layer
      !----------
!     a% geoh(:,:,1,:) = huge(1._sp)
      a% geoh(:,:,1,:) = a% geoh(:,:,2,:) + R * log (2._wp)* tv(:,:,1,:)
      !-------------------------------
      ! derive full level geopotential
      !-------------------------------
!     a% geof(:,:,1,:) = 0.5_wp * a% geoh(:,:,1,:)
!     do k = 2,kep
!       a% geof(:,:,k,:) = 0.5_wp * (a% geoh(:,:,k,:) + a% geoh(:,:,k+1,:))
!     end do
      if (lgeof) then
        if (.not.associated (a% pf)) call finish('set_geo','pf not associated')
        call allocate (a,'geof')
        do k = 1, kep
#ifdef _CRAYFTN
!DIR$ IVDEP
#endif
          a% geof(:,:,k,:) = a% geoh(:,:,k+1,:) + tv(:,:,k,:) &
                           * R * log (a% ph(:,:,k+1,:) / a% pf(:,:,k,:))
        end do
      endif
    case (VCT_Z_HYB, VCT_Z_GEN)
    !------------
    ! COSMO, ICON
    !------------
      if (.not.associated (a% grid% hhl)) &
           call finish('set_geo','hhl not associated')
      lb1 = a% grid% lb(1)
      ub1 = a% grid% ub(1)
#ifdef _CRAYFTN
!DIR$ IVDEP
#endif
      a% geoh = a% grid% hhl(lb1:ub1,:,:,:) * gacc
      !-------------------------------
      ! derive full level geopotential
      !-------------------------------
      if (lgeof) then
        call allocate (a,'geof')
        do k = 1, kep
#ifdef _CRAYFTN
!DIR$ IVDEP
#endif
          a% geof(:,:,k,:) = 0.5_wp * (a% geoh(:,:,k+1,:) + a% geoh(:,:,k,:))
        end do
      endif
    case (VCT_P_ISO)
    !-------------------------------
    ! Isobaric grid ("other models")
    !-------------------------------
      if (.not. associated (a% geof)) &
           call finish ('set_geo','geof not associated, but isobaric grid!')
      call deallocate (a,'geoh')
    case default
      call finish('set_geo','invalid vct')
    end select
  end subroutine set_geo
!------------------------------------------------------------------------------
  subroutine set_t1 (a)
    !------------------------
    ! vector version of set_t
    !------------------------
    type(t_atm), intent(inout) :: a(:)
    integer :: i
    do i = 1, size (a)
      call set_t (a(i))
    end do
  end subroutine set_t1
!------------------------------------------------------------------------------
  subroutine set_t (a)
    !---------------------------------------------------------
    ! calculate temperature from virtual potential temperature
    !---------------------------------------------------------
    type (t_atm) ,intent(inout)        :: a    ! atmospheric state
    !----------------
    ! local variables
    !----------------
    !-----------------------------------------------------
    ! sanity checks, allocate components in state variable
    !-----------------------------------------------------
    if (.not.associated (a% vptmp)) call finish('set_t','vptmp not associated')
    if (.not.associated (a% pf))    call finish('set_t','pf    not associated')
    if (.not.associated (a% q ))    call finish('set_t','q     not associated')
    call allocate (a,'tv')
    !--------------------------------------------------------------
    ! Derive virtual temperature from virtual potential temperature
    !--------------------------------------------------------------
    a% tv = a% vptmp * exner (a% pf)

    call allocate (a,'t')
    if (associated (a% qcl) .and. associated (a% qci)) then
      a% t = t_tv_q (a% tv, a% q, a% qcl, a% qci)
    else if (associated (a% qcl)) then
      a% t = t_tv_q (a% tv, a% q, a% qcl)
    else
      a% t = t_tv_q (a% tv, a% q)
    end if
  end subroutine set_t
!------------------------------------------------------------------------------
  subroutine set_tv_geo (a)
  !-------------------------------------------------------
  ! calculate virtual temperature from geopotential height
  ! (inverse of set_geo)
  !-------------------------------------------------------
  type (t_atm) ,intent(inout)        :: a    ! atmospheric state
    !----------------
    ! local variables
    !----------------
    integer  :: k, kep
    real(wp) :: dlnp (size(a% geoh,1),size(a% geoh,2),&
                      size(a% geoh,3),size(a% geoh,4))
    !-----------------------------------------------------
    ! sanity checks, allocate components in state variable
    !-----------------------------------------------------
    if(.not.associated(a%ph))   call finish('set_tv_geo','ph   not associated')
    if(.not.associated(a%geoh)) call finish('set_tv_geo','geoh not associated')
    call allocate (a,'tv')
    !-----------------------------------
    ! differentiate hydrostatic equation
    !-----------------------------------
    kep = a% ub(3)
    do k = kep, 2, -1
      dlnp (:,:,k,:) = R * log (a% ph(:,:,k+1,:) / a% ph(:,:,k,:))
      a% tv(:,:,k,:) = (a% geoh(:,:,k,:) - a% geoh(:,:,k+1,:)) / dlnp(:,:,k,:)
    end do
    a% tv(:,:,1,:) = (a% geoh(:,:,1,:) - a% geoh(:,:,2,:)) / (R * log (2._wp))
  end subroutine set_tv_geo
!------------------------------------------------------------------------------
  subroutine set_tq_tv (a)
  !------------------------------------------------------------
  ! calculate temperature from virtual temperature and moisture
  ! if specific humidity is not allocated, it is allocated and
  ! derived from relative humidity
  !------------------------------------------------------------
  type (t_atm) ,intent(inout)        :: a    ! atmospheric state
    !----------------
    ! local variables
    !----------------
    integer  :: irhq
    real(wp) :: x     (size(a% tv,1),size(a% tv,2),size(a% tv,3),size(a% tv,4))
    integer  :: i_fail(size(a% tv,1),size(a% tv,2),size(a% tv,3),size(a% tv,4))
    !---------------------------------
    ! local variables for error report
    !---------------------------------
    integer                        :: n_fail, k, ijd (3), n, i, if
    real(wp), dimension(dace% npe) :: pf, rh, tv, xx, dlat, dlon
    type (t_tq_tvrh)               :: xf
    real(wp)                       :: tt, qq
    !-----------------------------------------------------
    ! sanity checks, allocate components in state variable
    !-----------------------------------------------------
    irhq = 0
    if (associated (a% rh)) irhq = 1
    if (associated (a% q )) irhq = 2
    if (irhq == 0) call finish('set_t_tv','neither rh nor q is present')
    if (.not.associated (a% tv)) call finish('set_t_tv','tv not associated')
    call allocate (a,'t')
    !-------------------------------------
    ! sum liquid and ice content of clouds
    !-------------------------------------
    x = 0._wp
    if (associated(a% qcl)) x = x + a% qcl
    if (associated(a% qci)) x = x + a% qci
    if (associated(a% qr )) x = x + a% qr
    if (associated(a% qs )) x = x + a% qs
    if (associated(a% qg )) x = x + a% qg
    if (associated(a% qh )) x = x + a% qh
    !-------------------
    ! derive temperature
    !-------------------
    select case (irhq)
    case (1)
      !---------------
      ! t, q <- tv, rh
      !---------------
      if (.not.associated (a% pf)) call finish('set_t_tv','pf not associated')
      call allocate (a,'q')
      call tq_tvrh_noxfail (a% t, a% q, a% tv, a% rh, a% pf, x, i_fail=i_fail)
      n_fail = p_sum (count (i_fail < 0))
      if (n_fail/=0) then
        !-------------
        ! error report
        !-------------
        if (dace% lpio) then
          write (6,*)
          write (6,*) 'set_t_tv: tq_tvrh failed',n_fail,'times'
          write (0,*) 'set_t_tv: tq_tvrh failed',n_fail,'times'
        endif
        do k = 1, size (a% tv,3)
          !-----------------
          ! loop over levels
          !-----------------
          n      = count (i_fail(:,:,k,:) < 0)
          n_fail = p_sum (n)
          if (n_fail/=0) then
            write (6,*) 'set_t_tv: tq_tvrh failed',n_fail,'times at level',k
            !--------------------------------
            ! explore at most one fail per PE
            !--------------------------------
            pf           = 0._wp
            rh           = 0._wp
            tv           = 0._wp
            xx           = 0._wp
            dlat         = 0._wp
            dlon         = 0._wp
            if (n/=0) then
              ijd = minloc (i_fail(:,:,k,:))
              pf  (dace% pe+1) = a%       pf   (ijd(1), ijd(2), k, ijd(3))
              rh  (dace% pe+1) = a%       rh   (ijd(1), ijd(2), k, ijd(3))
              tv  (dace% pe+1) = a%       tv   (ijd(1), ijd(2), k, ijd(3))
              xx  (dace% pe+1) =          x    (ijd(1), ijd(2), k, ijd(3))
              dlat(dace% pe+1) = a% grid% rlat (ijd(1), ijd(2), 1, ijd(3)) * r2d
              dlon(dace% pe+1) = a% grid% rlat (ijd(1), ijd(2), 1, ijd(3)) * r2d
            endif
            pf           = p_sum (pf)
            rh           = p_sum (rh)
            tv           = p_sum (tv)
            xx           = p_sum (xx)
            dlat         = p_sum (dlat)
            dlon         = p_sum (dlon)
            if (dace% lpio) then
              do n=1,dace% npe
                if (pf(n) /= 0._wp) then
                  !--------------
                  ! rerun tq_tvrh
                  !--------------
                  write (6,*) 'pe,i,j,d=',n-1,ijd,                        &
                  ', pf=',pf(n),', rh=',rh(n),', tv=',tv(n), ', x=',xx(n),&
                  ', dlat=',dlat(n),', dlon=',dlon(n)
                  call tq_tvrh_noadj (tt, qq, tv(n), rh(n), pf(n), xx(n), &
                                      i_fail=if, x_fail=xf)
                  write (6,*) '      tmin,tmax:',xf%tmin, xf%tmax,&
                                  ', qmin,qmax:',xf%qmin, xf%qmax
                  do i = 1, abs(if)
                    write (6,*) '      iter:',i,&
                      't,q,tv,rh:',xf%t(i),xf%q(i),xf%tv(i),xf%rh(i)
                  end do
                endif
              end do
            endif
          endif
        end do
        call finish('set_t_tv','tq_tvrh failed')
      endif
    case (2)
      !-----------
      ! t <- tv, q
      !-----------
      a% t = t_tv_q (a% tv, a% q, x)
    end select
  end subroutine set_tq_tv
!------------------------------------------------------------------------------
  subroutine set_tv1 (states)
  !-------------------------------------------
  ! set virtual temperature ensemble of states
  !-------------------------------------------
  !----------------------
  ! set pressure level(s)
  !----------------------
  type (t_atm) ,intent(inout) :: states (:)
    integer :: i
    do i = 1, size (states)
      call set_tv (states (i))
    end do
  end subroutine set_tv1
!------------------------------------------------------------------------------
  subroutine set_tv (a)
  !------------------------------
  ! calculate virtual temperature
  !------------------------------
  type (t_atm) ,intent(inout) :: a
    !-----------------------------------------------------------------------
    ! local variables: indices, virtual temp., cloud ice+water , log delta p
    !-----------------------------------------------------------------------
    real(wp) ,allocatable :: x (:,:,:,:)
    if(.not.associated(a% t)) call finish('set_tv','t is not associated')
    if(.not.associated(a% q)) call finish('set_tv','q is not associated')
    allocate (x (size(a% t,1),size(a% t,2),size(a% t,3),size(a% t,4)))
    !-------------------------------------
    ! allocate component in state variable
    !-------------------------------------
    call allocate (a,'tv')
    !---------------------------
    ! derive virtual temperature
    !---------------------------
    x = 0._wp
    if (associated(a% qcl)) x = x + a% qcl
    if (associated(a% qci)) x = x + a% qci
    if (associated(a% qr )) x = x + a% qr
    if (associated(a% qs )) x = x + a% qs
    if (associated(a% qg )) x = x + a% qg
    if (associated(a% qh )) x = x + a% qh
    a% tv = tv_t_q (a% t, a% q, x)
    deallocate (x)
  end subroutine set_tv
!------------------------------------------------------------------------------
  subroutine set_den (a)
  !--------------------------------------------------------
  ! calculate density from virtual temperature and pressure
  !--------------------------------------------------------
  type (t_atm) ,intent(inout) :: a
    !-----------------------------------------------------------------------
    ! local variables: indices, virtual temp., cloud ice+water , log delta p
    !-----------------------------------------------------------------------
    if(.not.associated(a% tv)) call finish('set_den','tv is not associated')
    if(.not.associated(a% pf)) call finish('set_den','pf is not associated')
    !-------------------------------------
    ! allocate component in state variable
    !-------------------------------------
    call allocate (a,'den')
    !---------------
    ! derive density
    !---------------
    a% den = a% pf / (R * a% tv)
  end subroutine set_den
!------------------------------------------------------------------------------
  subroutine set_q1 (states)
  !---------------------------------------------
  ! set specific humidity for ensemble of states
  !---------------------------------------------
  type (t_atm) ,intent(inout) :: states (:)
    integer :: i
    do i = 1, size (states)
      call set_q (states (i))
    end do
  end subroutine set_q1
!------------------------------------------------------------------------------
  subroutine set_q (a)
  !----------------------------
  ! calculate specific humidity
  !----------------------------
  type (t_atm) ,intent(inout) :: a
    !--------------------------------------
    ! check for presence of required fields
    !--------------------------------------
    if(.not.associated(a% t )) call finish('set_q','t  is not associated')
    if(.not.associated(a% rh)) call finish('set_q','rh is not associated')
    if(.not.associated(a% pf)) call finish('set_q','pf is not associated')
    !-------------------------------------
    ! allocate component in state variable
    !-------------------------------------
    call allocate (a,'q')
    !-------------------------
    ! derive specific humidity
    !-------------------------
    a% q = q_rh (a% rh, a% t, a% pf)
  end subroutine set_q
!------------------------------------------------------------------------------
  subroutine set_rh1 (states)
  !---------------------------------------------
  ! set relative humidity for ensemble of states
  !---------------------------------------------
  type (t_atm) ,intent(inout) :: states (:)
    integer :: i
    do i = 1, size (states)
      call set_rh (states (i))
    end do
  end subroutine set_rh1
!------------------------------------------------------------------------------
  subroutine set_rh (a)
  !----------------------------
  ! calculate relative humidity
  !----------------------------
  type (t_atm) ,intent(inout) :: a
    !--------------------------------------
    ! check for presence of required fields
    !--------------------------------------
    if(.not.associated(a% t )) call finish('set_rh','t  is not associated')
    if(.not.associated(a% q )) call finish('set_rh','q  is not associated')
    if(.not.associated(a% pf)) call finish('set_rh','pf is not associated')
    !-------------------------------------
    ! allocate component in state variable
    !-------------------------------------
    call allocate (a,'rh')
    !-------------------------
    ! derive relative humidity
    !-------------------------
!NEC$ inline
    a% rh = rh_q (a% q, a% t, a% pf)
    !------------------------------------------
    ! derive 2m relative humidity from dewpoint
    !------------------------------------------
    if (associated (a% t2m) .and. associated (a% td2m)) call set_rh2m (a)
  end subroutine set_rh
!------------------------------------------------------------------------------
  subroutine set_rh2m1 (states, rhmax)
    type(t_atm) ,intent(inout)          :: states(:)
    real(wp)    ,intent(in)   ,optional :: rhmax    ! Upper limit on rel.hum.
    !-----------------------------------------------------------
    ! derive 2m relative humidity from dewpoint (vector version)
    !-----------------------------------------------------------
    integer :: i
    do i = 1, size (states)
      call set_rh2m (states(i), rhmax)
    end do
  end subroutine set_rh2m1
!------------------------------------------------------------------------------
  subroutine set_rh2m (a, rhmax)
    type(t_atm) ,intent(inout)          :: a
    real(wp)    ,intent(in)   ,optional :: rhmax    ! Upper limit on rel.hum.
    !------------------------------------------
    ! derive 2m relative humidity from dewpoint
    !------------------------------------------
    if (.not.associated (a% t2m )) call finish('set_rh2m','t2m  not associated')
    if (.not.associated (a% td2m)) call finish('set_rh2m','td2m not associated')
    call allocate (a,'rh2m')
    a% rh2m = esw_t (a% td2m) / esw_t (a% t2m)
    if (present (rhmax)) then
      a% rh2m = min (a% rh2m, rhmax)
    end if
  end subroutine set_rh2m
!------------------------------------------------------------------------------
  subroutine set_rh2m_land (a, rhmax)
    type(t_atm) ,intent(inout)          :: a
    real(wp)    ,intent(in)   ,optional :: rhmax    ! Upper limit on rel.hum.
    !------------------------------------------
    ! derive 2m relative humidity from dewpoint
    ! for land tile averages (done separately)
    !------------------------------------------
    if (.not.associated (a% t2m_land )) call finish('set_rh2m_land',          &
                                                    't2m_land  not associated')
    if (.not.associated (a% td2m_land)) call finish('set_rh2m_land',          &
                                                    'td2m_land not associated')
    call allocate (a,'rh2m_land')
    a% rh2m_land = esw_t (a% td2m_land) / esw_t (a% t2m_land)
    if (present (rhmax)) then
      a% rh2m_land = min (a% rh2m_land, rhmax)
    end if
  end subroutine set_rh2m_land
!------------------------------------------------------------------------------
  subroutine set_rh2m_land1 (states, rhmax)
    type(t_atm) ,intent(inout)          :: states(:)
    real(wp)    ,intent(in)   ,optional :: rhmax    ! Upper limit on rel.hum.
    !-----------------------------------------------------------
    ! derive 2m relative humidity from dewpoint (vector version)
    !-----------------------------------------------------------
    integer :: i
    do i = 1, size (states)
      call set_rh2m_land (states(i), rhmax)
    end do
  end subroutine set_rh2m_land1
!------------------------------------------------------------------------------
  subroutine set_wso (a, mode)
    type(t_atm) ,intent(inout) :: a
    integer     ,intent(in)    :: mode
    !----------------------------------------------
    ! derive soil moisture from soil moisture index
    !----------------------------------------------
      integer :: l
      if (.not.associated (a% w_so))          call finish('set_wso','w_so not associated')
      if (.not.associated (a% grid% soiltyp)) call finish('set_wso','soiltyp not present')
      do l = 1, a% grid% ns
        a% w_so (:,:,l,:) = wso_ind (a% w_so          (:,:,l,:) ,                      &
                                int (a% grid% soiltyp (a% grid% lb(1):a% grid% ub(1),  &
                                                       a% grid% lb(2):a% grid% ub(2),  &
                                                     1,a% grid% lb(4):a% grid% ub(4))),&
                                                        l, mode                        )
      end do
  end subroutine set_wso
!------------------------------------------------------------------------------
  subroutine set_index (a, mode)
    type(t_atm) ,intent(inout) :: a
    integer     ,intent(in)    :: mode
    !----------------------------------------------
    ! derive soil moisture index from soil moisture
    !----------------------------------------------
      integer :: l
      if (.not.associated(a% w_so))          call finish('set_index','w_so not associated')
      if (.not.associated(a% grid% soiltyp)) call finish('set_index','soiltyp not present')
      do l = 1, a% grid% ns
        a% w_so (:,:,l,:) = ind_wso (a% w_so          (:,:,l,:) ,                      &
                                int (a% grid% soiltyp (a% grid% lb(1):a% grid% ub(1),  &
                                                       a% grid% lb(2):a% grid% ub(2),  &
                                                     1,a% grid% lb(4):a% grid% ub(4))),&
                                                        l, mode                        )
      end do
  end subroutine set_index
!------------------------------------------------------------------------------
  subroutine set_psi_chi (a)
    type (t_atm) ,intent(inout) :: a
    !----------------------------------------------------------------
    ! calculate streamfunction and velocity potential
    ! - Global fields, lat-lon grid: requires u,v as input
    !   (based on spherepack, no parallelization)
    ! - Limited area (rotated) lat-lon: requires vorticity,divergence
    !----------------------------------------------------------------
    integer  :: ig
    real(wp) :: dx, dy
    !-------------------------------------
    ! allocate component in state variable
    !-------------------------------------
    call allocate (a,'psi')
    call allocate (a,'chi')
    !---------------------------------------------
    ! derive streamfunction and velocity potential
    !---------------------------------------------
    if (a% grid% global) then
      call sp95_init
      ig = 0; if (a%grid% dlat(2) < a%grid% dlat(1)) ig = 1  ! 1 for N->S
      select case (a% grid% gridtype)
      case (WMO6_GAUSSIAN)
        call sp95_sfvp_gg (a% u  (:,:,:,1), a% v  (:,:,:,1),            &
                           a% psi(:,:,:,1), a% chi(:,:,:,1), a% grid% a,&
                           ig                                           )
      case (WMO6_LATLON)
        if (ig/=0) call finish('set_psi_chi',                            &
                               'sp95_sfvp: N->S ordering not implemented')
        call sp95_sfvp    (a% u  (:,:,:,1), a% v  (:,:,:,1),            &
                           a% psi(:,:,:,1), a% chi(:,:,:,1), a% grid% a )
      case default
        call finish ('set_psi_chi','grid is not Gaussian or LatLon')
      end select
      call sp95_cleanup
    else
      if (.not. associated (a% vrt)) then
         call finish ("set_psi_chi","vorticity not associated!")
      end if
      if (.not. associated (a% div)) then
         call finish ("set_psi_chi","divergence not associated!")
      end if
      select case (a% grid% gridtype)
      case (WMO6_LATLON, WMO6_ROTLL)
        dx = a% grid% a * d2r * a% grid% di
        dy = a% grid% a * d2r * a% grid% dj
        call poisson_solve_dct (a% psi(:,:,:,1), a% vrt(:,:,:,1), dx, dy)
        call poisson_solve_dct (a% chi(:,:,:,1), a% div(:,:,:,1), dx, dy)
      case default
        call finish ('set_psi_chi','grid is not LatLon')
      end select
    end if
  end subroutine set_psi_chi
!------------------------------------------------------------------------------
  subroutine set_psi_chi_u_v (a)
  !----------------------------------------------------------------
  ! Calculate streamfunction, velocity potential, and hor. wind u,v
  ! (using Spherepack, no parallelization)
  !----------------------------------------------------------------
    type (t_atm) ,intent(inout) :: a
    integer :: ig
    if (.not. associated (a% vrt)) then
       call finish ("set_psi_chi_u_v","vorticity not associated!")
    end if
    if (.not. associated (a% div)) then
       call finish ("set_psi_chi_u_v","divergence not associated!")
    end if
    !--------------------------------------
    ! Allocate components in state variable
    !--------------------------------------
    call allocate (a,'psi')
    call allocate (a,'chi')
    call allocate (a,'u')
    call allocate (a,'v')
    !---------------------------------------------------
    ! Derive streamfunction, velocity potential, u and v
    !---------------------------------------------------
    call sp95_init ()
    ig = 0; if (a%grid% dlat(2) < a%grid% dlat(1)) ig = 1  ! 1 for N->S
    select case (a% grid% gridtype)
    case (WMO6_GAUSSIAN)
      call sp95_sfvp_vrtdiv_gg (a% vrt(:,:,:,1), a% div(:,:,:,1), &
                                a% psi(:,:,:,1), a% chi(:,:,:,1), &
                                a% u  (:,:,:,1), a% v  (:,:,:,1), &
                                rsphere=a% grid% a, ig=ig)
    case default
      call finish ('set_psi_chi_u_v','grid is not Gaussian')
    end select
    call sp95_cleanup ()
  end subroutine set_psi_chi_u_v
!------------------------------------------------------------------------------
  subroutine set_vrt_div (a)
    type (t_atm) ,intent(inout) :: a
    !--------------------------------------------
    ! Calculate vorticity and divergence from u,v
    !--------------------------------------------
    if (.not. associated (a% u)) call finish ("set_vrt_div","u not associated!")
    if (.not. associated (a% v)) call finish ("set_vrt_div","v not associated!")
    !-------------------------------------
    ! allocate component in state variable
    !-------------------------------------
    call allocate (a,'vrt')
    call allocate (a,'div')
    select case (a% grid% gridtype)
    case (WMO6_LATLON, WMO6_ROTLL, DWD6_ICON)
       call hor_curl (a% u, a% v, a% vrt, a% grid)
       call hor_div  (a% u, a% v, a% div, a% grid)
    case default
       call finish ('set_vrt_div','grid is not LatLon')
    end select
  end subroutine set_vrt_div
!------------------------------------------------------------------------------
  subroutine set_tqv (state)
  !---------------------------
  ! set integrated water vpour
  !---------------------------

    type (t_atm) ,intent(inout) :: state

    real (wp), pointer :: rhowprofile(:,:)
    real (wp) :: e
    integer   :: Nx1, Nx2, Ny1, Ny2, Nz1, Nz2, Nr1, Nr2
    integer   :: i, j, k, r

    call allocate (state, 'tqv')

    select case (state% grid% vct)
    !case (IVCTYPE_GME)
    case (VCT_P_HYB)
       !----
       ! GME
       !----
       if (.not.associated (state% pf))   call set_p   (state)
       if (.not.associated (state% geof)) call set_geo (state, geof=.true.)
    case (VCT_Z_HYB, VCT_Z_GEN)
       select case (state% grid% model)
       !case (IVCTYPE_ICON)
       case (MO_ICON)
          !-----
          ! ICON
          !-----
          !call set_ps  (atm)
          !call set_geo (atm, geof=.true.)
       case (MO_COSMO)
          !------
          ! COSMO
          !------
          if (.not.associated (state% pf)) call set_p (state)
       end select
    case default
       !--------------------------------
       ! try so set all necessary fields
       !--------------------------------
       if (.not.associated (state% pf)) call set_p (state)
       if (.not.associated (state% geof)) call set_geo (state, geof=.true.)
    end select

    if (.not.associated (state% pf))                                    &
                call finish('set_tqv','pf not associated')
    if (.not.associated (state% q))                                     &
                call finish('set_tqv','q not associated')
    if (.not.associated (state% t))                                     &
                call finish('set_tqv','t not associated')
    if (.not.associated (state% tqv))                                   &
                call finish('set_tqv','tqv not associated')

    ! Number of points in each dimension:
    Nx1 = lbound(state%t,1)
    Nx2 = ubound(state%t,1)
    Ny1 = lbound(state%t,2)
    Ny2 = ubound(state%t,2)
    Nz1 = lbound(state%t,3)
    Nz2 = ubound(state%t,3)

    ! GME: diamond index
    Nr1 = lbound(state%t,4)
    Nr2 = ubound(state%t,4)

    ! Array for one absolute humidity profile:
    ! rhowprofile(i,j) - i - height index
    !                    j = 1 - altitude, m
    !                    j = 2 - absolute humidity
    allocate( rhowprofile(Nz1:Nz2,1:2) )

    do r=Nr1, Nr2
       do j=Ny1, Ny2
          do i=Nx1, Nx2
             do k=Nz1, Nz2
                ! Partial pressure of water vapour e
                !e = (state%pf(i,j,k,r) * state%q(i,j,k,r)) /         &
                !     (0.62132_wp + state%q(i,j,k,r)*0.37868_wp)
                e = (state%pf(i,j,k,r) * state%q(i,j,k,r)) /         &
                     (RDRD + state%q(i,j,k,r)*EMRDRD)
                rhowprofile(k,2) =  e / (RD * state%t(i,j,k,r))

                select case (state% grid% vct)
                case (VCT_P_HYB)
                   ! height above geoid, estimated from geopotential
                   rhowprofile(k,1) = geopot_geoid( state%grid%rlat(i,j,1,r), &
                                               state%geof(i,j,k,r), .true. )
                case default
                   ! Field HHL is available in ICON and COSMO
                   ! HHL - half levels, up to Nz2+1 layers
                   !        => interpolate to full levels
                   rhowprofile(k,1) = 0.50_wp * ( state%grid%hhl(i,j,k,1) +   &
                        state%grid%hhl(i,j,k+1,1) )
                end select
             end do
             ! Integrate absolute humidity from top to bottom (=> minus)
             state%tqv(i,j,1,r) = - integpolycube(rhowprofile)
          end do
       end do
    end do

    if (associated(rhowprofile)) deallocate(rhowprofile)

  end subroutine set_tqv
!------------------------------------------------------------------------------
  subroutine print_state (state, iunit, intend, verbose, comment, grid)
  type (t_atm)     ,intent(in)           :: state  ! variable to print
  integer          ,intent(in) ,optional :: iunit  ! unit   (default=6)
  character(len=*) ,intent(in) ,optional :: intend ! intend string ('')
  logical          ,intent(in) ,optional :: verbose
  character(len=*) ,intent(in) ,optional :: comment
  logical          ,intent(in) ,optional :: grid   ! print grid?
  !-------------------------
  ! print t_grid (formatted)
  ! (for debugging)
  !-------------------------
    integer           :: iu, n
    character(len=32) :: c
!   logical           :: v
    logical           :: lgrid
    ENTERFUNCTION
    !--------------------
    ! optional parameters
    !--------------------
    n=0
    c=''
    if (present(intend)) then
      c = intend
      n = len(c)
    endif
    iu = 6; if (present(iunit)) iu = iunit
!   v = .false.; if (present(verbose)) v = verbose
    lgrid = .true. ; if (present (grid)) lgrid = grid
    !------
    ! print
    !------
    if (dace% lpio) then
      if (present(comment)) then
        write(iu,'(a)') repeat('_',79)
        write(iu,'()')
        write(iu,'(4x,a)') trim(comment)
        write(iu,'()')
      endif
      write(iu,"(a,' (t_atm)')") c(:n)
      write(iu,"(a,' name           :',a)"  ) c(:n), state% name
      call print(state% time,     iunit=iu, intend=c(:n)//' time%')
      call print(state% ref_time, iunit=iu, intend=c(:n)//' ref_time%')
    endif
    if(associated(state% grid)) then
      if (lgrid) call print(state% grid, iunit=iu, &
                              intend=c(:n)//' grid%',verbose=verbose)
    else
      write(iu,"(a,' grid: NOT ASSOCIATED!')") c(:n)
    endif
    if (dace% lpio) then
      write(iu,"(a,' run type       :' ,a8)" ) c(:n), state% runtype
      write(iu,"(a,' run class      :' ,i8)" ) c(:n), state% runclass
      write(iu,"(a,' experiment Id  :' ,i8)" ) c(:n), state% expid
      write(iu,"(a,' ensemble member:' ,i8)" ) c(:n), state% member
      write(iu,"(a,' ensemble size  :' ,i8)" ) c(:n), state% members
      write(iu,"(a,' lb             :' ,4i8)") c(:n), state% lb
      write(iu,"(a,' ub             :' ,4i8)") c(:n), state% ub
      write(iu,"(a,' size           :' ,i12)") c(:n), state% size
      write(iu,"(a,' gp_wind        : ',a)"  ) c(:n), state% gp_wind
    endif
    call print(state% m, iunit=iu, intend=c(:n),verbose=verbose)
    if (dace% lpio) then
      write(iu,"(a,' size =',i12)")           c(:n), nwords(state% m)
    endif
    DELETESTORAGE (state)
    LEAVEFUNCTION
  end subroutine print_state
!------------------------------------------------------------------------------
  subroutine atm_to_grads (ctl, atm, t, comment, iostat)
  type (t_ctl)               ,intent(inout) :: ctl     ! GRADS metadata
  type (t_atm)               ,intent(in)    :: atm     ! atmospheric state
  integer          ,optional ,intent(in)    :: t       ! time slice
  character(len=*) ,optional ,intent(in)    :: comment ! comment to write
  integer          ,optional ,intent(out)   :: iostat  ! I/O error status
  !--------------------------------------
  ! write atmospheric state to GRADS file
  !--------------------------------------
    ENTERFUNCTION
    call to_grads (ctl, atm% m, atm% grid% dc,                               &
                   yrev    = atm%grid%dlat(1) > atm%grid%dlat(atm%grid% ny), &
                   t       = t,                                              &
                   comment = comment,                                        &
                   iostat  = iostat                                          )
    DELETESTORAGE (atm)
    LEAVEFUNCTION
  end subroutine atm_to_grads
!------------------------------------------------------------------------------
  subroutine atm_from_grads (ctl, atm, t)
  type (t_ctl)               ,intent(in)    :: ctl ! GRADS metadata
  type (t_atm)               ,intent(inout) :: atm ! atmospheric state
  integer          ,optional ,intent(in)    :: t   ! time slice
  !---------------------------------------
  ! read atmospheric state from GRADS file
  !---------------------------------------
    ENTERFUNCTION
    call from_grads (ctl, atm% m, atm% grid% dc,                            &
                     yrev = atm%grid%dlat(1) > atm%grid%dlat(atm%grid% ny), &
                     t    = t                                               )
    DELETESTORAGE (atm)
    LEAVEFUNCTION
  end subroutine atm_from_grads
!------------------------------------------------------------------------------
  subroutine dump_state (state, unit, destr)
  type (t_atm) ,intent(inout) :: state
  integer      ,intent(in)    :: unit
  logical      ,intent(in)    :: destr
    character(len=11) :: formatted
    inquire(unit,form=formatted)
    if (formatted=='FORMATTED') then
      write (unit,*) state% time,   state% ref_time, state% lb, state% ub, &
                     state% member, state% members
    else
      write (unit)   state% time,   state% ref_time, state% lb, state% ub, &
                     state% member, state% members
    endif
    call dump (state%m, unit)
    if (destr) call destruct (state)
  end subroutine dump_state
!------------------------------------------------------------------------------
  subroutine retrieve_state (state, unit, grid)
  type (t_atm)   ,intent(out) :: state
  integer        ,intent(in)  :: unit
  type (t_grid)  ,pointer     :: grid
    character(len=11) :: formatted
    inquire(unit,form=formatted)
    if (formatted=='FORMATTED') then
      read (unit,*) state% time,   state% ref_time, state% lb, state% ub, &
                    state% member, state% members
    else
      read (unit)   state% time,   state% ref_time, state% lb, state% ub, &
                    state% member, state% members
    endif
    state% grid => grid
    call retrieve (state%m, unit)
    call set_pointers (state)
  end subroutine retrieve_state
!------------------------------------------------------------------------------
  subroutine plot_atm_t (x,base,i,j)
  type (t_atm) ,intent(in)           :: x (:)
  character (len=*)     ,intent(in)           :: base
  integer               ,intent(in) ,optional :: i,j
    integer  :: ii, jj, k, l, lb, ub, iu, ios
    real(wp) :: h0
    ENTERFUNCTION
    h0 = hours (x(1)%time)
    iu = get_unit_number ()
    do k = 1,size(x(1)%m)
      if (x(1)%m(k)%i% alloc) then
        open(iu,file=trim(base)//x(1)%m(k)%i% name,iostat=ios)
        if (ios/=0) cycle
        ii=x(1)%m(k)%i%lb(1); if(present(i)) ii=i
        jj=x(1)%m(k)%i%lb(2); if(present(j)) jj=j
        lb = x(1)%m(k)%i%lb(3)
        ub = x(1)%m(k)%i%ub(3)
        do l = 1,size(x)
          write(iu,*) l,hours (x(l)%time)-h0, x(l)%m(k)%ptr(ii,jj,ub:lb:-1,1)
        end do
        close(iu)
      endif
    end do
    call return_unit_number (iu)
#   ifndef TR15581
      do l=1,size(x)
        DELETESTORAGE (x(l))
      end do
#   endif
    LEAVEFUNCTION
  end subroutine plot_atm_t
!------------------------------------------------------------------------------
  subroutine plot_atm_v (x,base,i,j)
  type (t_atm) ,intent(in)           :: x (:)
  character (len=*)     ,intent(in)           :: base
  integer               ,intent(in) ,optional :: i,j
    integer  :: ii, jj, k, l, m, lb, ub, iu, ios
    real(sp) :: pf, y(size(x))
    integer  :: id = 1 ! ++ diamond
    ENTERFUNCTION
    iu = get_unit_number ()
    do k = 1,size(x(1)%m)
      if (x(1)%m(k)%i% alloc) then
        open(iu,file=trim(base)//x(1)%m(k)%i% name,iostat=ios)
        if (ios/=0) cycle
        ii=x(1)%m(k)%i%lb(1); if(present(i)) ii=i
        jj=x(1)%m(k)%i%lb(2); if(present(j)) jj=j
        lb = x(1)%m(k)%i%lb(3)
        ub = x(1)%m(k)%i%ub(3)
        do l = ub,lb,-1
          if (ub==ubound(x(1)%pf,3)) then
            pf = x(1)%pf(ii,jj,l,id)
          else if (ub==ubound(x(1)%ph,3)) then
            pf = x(1)%ph(ii,jj,l,id)
          else
            pf = l
          endif
          y  = (/(x(m)%m(k)%ptr(ii,jj,l,1),m=1,size(x))/)
          if (associated(x(1)%geof)) then
            write(iu,*) l,pf,x(1)%geof(ii,jj,l,id),y
          else
            write(iu,'(i4,200g12.3)') l,pf,y
          endif
        end do
        close(iu)
      endif
    end do
    call return_unit_number (iu)
#   ifndef TR15581
      do m=1,size(x)
        DELETESTORAGE (x(m))
      end do
#   endif
    LEAVEFUNCTION
  end subroutine plot_atm_v
!------------------------------------------------------------------------------
  function  nwords_atm (atm) result (n)
  type (t_atm) ,intent(in) :: atm
  integer(i8)              :: n
  ENTERFUNCTION
  n = nwords (atm%m)
  DELETESTORAGE (atm)
  LEAVEFUNCTION
  end function nwords_atm
!------------------------------------------------------------------------------
  function real_atm (atm) result (r)
  type (t_atm) ,intent(in) :: atm
  real(wp)                 :: r (atm%size)
    integer     :: i
    integer(i8) :: l,n
    ENTERFUNCTION
    l=0
    do i=1,size(atm%m)
      if (atm%m(i)%i% alloc) then
        n = size(atm%m(i)%ptr)
        r(l+1:l+n) = reshape (atm%m(i)% ptr, (/n/))
        l = l+n
      endif
    end do
    DELETESTORAGE (atm)
    LEAVEFUNCTION
  end function real_atm
!------------------------------------------------------------------------------
  subroutine assign_array_to_atm (atm, r)
  type (t_atm) ,intent(inout) :: atm
  real(wp)     ,intent(in)    :: r (:)
    integer     :: i
    integer(i8) :: l,n
    l=0
    if (size(r)/=nwords (atm%m)) &
      call finish ('assign_array_to_atm','size(r)/=nwords(atm%m)')
    do i=1,size(atm%m)
      if (atm%m(i)%i% alloc) then
        n = size(atm%m(i)%ptr)
        atm%m(i)% ptr = reshape (r(l+1:l+n), shape(atm%m(i)% ptr))
        l = l+n
      endif
    end do
  end subroutine assign_array_to_atm
!==============================================================================
  subroutine cut_atm (y, x, lb, ub)
  type (t_atm) ,intent(inout) :: y
  type (t_atm) ,intent(in)    :: x
  integer      ,intent(in)    :: lb(2)
  integer      ,intent(in)    :: ub(2)
    ENTERFUNCTION
    if (any(y%lb(1:2)/=lb.or.y%ub(1:2)/=ub)) call destruct (y%m)
    y% time        =  x% time
    y% ref_time    =  x% ref_time
    y% grid        => x% grid
    y% lb(3)       =  x% lb(3);  y% lb(1:2) = lb
    y% ub(3)       =  x% ub(3);  y% ub(1:2) = ub
    y% member      =  x% member
    y% members     =  x% members
    call cut (y% m, x% m, lb, ub)
    call set_pointers (y)
    DELETESTORAGE (x)
    LEAVEFUNCTION
  end subroutine cut_atm
!------------------------------------------------------------------------------
  subroutine merge_atm (dest, src, par)
    type(t_atm),  intent(inout) :: dest         ! destination state
    type(t_atm),  intent(in)    :: src          ! source state
    character(*), intent(in)    :: par          ! list of variables to merge
    !------------------------------------------------
    ! Merge fields from source into destination state
    !------------------------------------------------
    character(len=16) :: pars(nm) ! names of variables to merge
    character(len=16) :: name     ! name of variable
    integer           :: n        ! actual number of variables
    integer           :: m        ! loop index
    ENTERFUNCTION
    call split (pars, par, n)
    do m = 1, size (dest% m)
       name = src% m(m)% i% name
       if (any (pars(:n) == name)) then
          if (.not. src % m(m)% i% alloc) &
               call finish ("merge_atm","not allocated: "//trim(name))
          if (.not. dest% m(m)% i% alloc) call allocate (dest, name)
          dest% m(m)% ptr = src% m(m)% ptr
       end if
    end do
    DELETESTORAGE (src)
    LEAVEFUNCTION
  end subroutine merge_atm
!==============================================================================
  subroutine coarsen_atm (y, x, stride)
  type (t_atm) ,intent(inout) :: y
  type (t_atm) ,intent(in)    :: x
  integer      ,intent(in)    :: stride
    !----------------
    ! local variables
    !----------------
    real(wp), pointer :: xg  (:,:,:) ! temporary global array
    real(wp), pointer :: yg  (:,:,:) ! temporary global array
    integer           :: lgx (4)     ! lower bound of source      (global)
    integer           :: ugx (4)     ! upper bound of source      (global)
    integer           :: lgy (4)     ! lower bound of destination (global)
    integer           :: ugy (4)     ! upper bound of destination (global)
    integer           :: lx  (4)     ! lower bound of source      (local)
    integer           :: ux  (4)     ! upper bound of source      (local)
!   integer           :: ly  (4)     ! lower bound of destination (local)
!   integer           :: uy  (4)     ! upper bound of destination (local)
    integer           :: i           ! field index
    integer           :: k           ! level index
    !----------------------
    ! executable statements
    !----------------------
    ENTERFUNCTION
    y% time     =  x% time
    y% ref_time =  x% ref_time
    y% member   =  x% member
    y% members  =  x% members
    !------------------------
    ! get global array bounds
    !------------------------
    lgx = x% grid% lbg
    ugx = x% grid% ubg
    lgy = y% grid% lbg
    ugy = y% grid% ubg
    !--------------------------
    ! allocate temporary arrays
    !--------------------------
    if (dace% lpio) then
      allocate (xg (lgx(1):ugx(1),lgx(2):ugx(2),lgx(4):ugx(4)))
      allocate (yg (lgy(1):ugy(1),lgy(2):ugy(2),lgy(4):ugy(4)))
    else
      nullify (xg, yg)
    endif
    !--------------------------------------
    ! loop over variables, get local bounds
    !--------------------------------------
    do i=1, size(x%m)
      if (x%m(i)%i% alloc .and. y%m(i)%i% alloc) then
        lx = x%m(i)%i% lb
        ux = x%m(i)%i% ub
!       ly = y%m(i)%i% lb
!       uy = y%m(i)%i% ub
        !-------------------------------------------
        ! loop over levels: gather, coarsen, scatter
        !-------------------------------------------
        do k = lx(3), ux(3)
          call gather_level  (xg, x%m(i)% ptr(:,:,k,:), x%grid%dc, dace% pio)
          if (dace% lpio) yg = xg(::stride,::stride,:)
          call scatter_level (yg, y%m(i)% ptr(:,:,k,:), y%grid%dc, dace% pio)
        end do
      else
        !---------------------------
        ! mark variables not present
        !---------------------------
        y%m(i)%i% alloc = .false.
      endif
    end do
    !--------
    ! cleanup
    !--------
    if (dace% lpio) deallocate (xg, yg)
    call update (y% m)
    call set_pointers (y)
    DELETESTORAGE (x)
    LEAVEFUNCTION
  end subroutine coarsen_atm
!==============================================================================
  subroutine hor_intp_coef_atm (gin, gout, l_soil, l_land, l_sea, ind_atm,    &
                    ind_soil, ind_land, ind_sea, fr_lake, fr_land, fr_land_sst)
  type (t_grid) ,intent(inout)        :: gin
  type (t_grid) ,intent(inout)        :: gout
  logical       ,intent(in)           :: l_soil           ! interpolate soil ?
  logical       ,intent(in)           :: l_land           !          surface ?
  logical       ,intent(in)           :: l_sea            !          surface ?
  type (t_h_idx),pointer              :: ind_atm  (:,:,:) ! interp.coef. atm.
  type (t_h_idx),pointer    ,optional :: ind_soil (:,:,:) ! interp.coef. soil
  type (t_h_idx),pointer    ,optional :: ind_land (:,:,:) ! interp.coef. surf.
  type (t_h_idx),pointer    ,optional :: ind_sea  (:,:,:) ! interp.coef. surf.
  real(wp)      ,intent(in) ,optional :: fr_lake(2) ! lake fraction bound(in/out)
  real(wp)      ,intent(in) ,optional :: fr_land(2) ! land fraction bound(in/out)
  real(wp)      ,intent(in) ,optional :: fr_land_sst(2)
  !--------------------------------------------------------------------------
  ! Calculate source grid indices and interpolation coefficients
  ! for horizontal interpolation.
  ! Optionally nearest neighbour interpolation for surface and soil variables
  !--------------------------------------------------------------------------

!   type (t_h_idx) ,pointer :: ind      (:,:,:)     ! interpolation coeffs.
    logical                 :: same_lsm ( 3  , 3  ) ! same land/sea list
    logical                 :: same_sty (n_st,n_st) ! same soiltyp  list
    logical                 :: same_ris (n_st,n_st) ! same characteristics
    integer                 :: i, j                 ! soiltyp indices
    integer    ,allocatable :: stin  (:,:,:)        ! input  'soil type'
    integer    ,allocatable :: stout (:,:,:)        ! output 'soil type'
    real(wp)   ,allocatable :: di    (:,:,:)        ! input  'soil type'
    real(wp)   ,allocatable :: do    (:,:,:)        ! output 'soil type'
    real(wp)                :: frl   (2)            ! lake fraction bounds
    real(wp)                :: fls   (2)            ! land/sea fraction bounds
    real(wp)                :: flss  (2)            ! land/sea fraction bounds
    logical    ,allocatable :: lsi   (:,:,:)        ! input  land-sea mask
    logical    ,allocatable :: lso   (:,:,:)        ! output land-sea mask

    !---------------------------------------
    ! set up cross-references for soil types
    !---------------------------------------
    frl  = 0.05_wp; if (present(fr_lake    )) frl  = fr_lake
    fls  = 0.05_wp; if (present(fr_land    )) fls  = fr_land
    flss = 0.05_wp; if (present(fr_land_sst)) flss = fr_land_sst
    same_lsm = .false.
    same_sty = .false.
    same_ris = .false.

    do i = 1, 3
      same_lsm (i,i) = .true.
    end do

    do j = 1, n_st
    do i = 1, n_st
      !---------------
      ! same soil type
      !---------------
      if (i==j) same_sty (i,j) = .true.
      !-----------------------------------------------
      ! same characteristics: sea - rock - ice - other
      !-----------------------------------------------
      select case (j)
      case (ST_SEAWATER, ST_SEAICE)
        select case (i)
        case (ST_SEAWATER, ST_SEAICE)
          same_ris (i,j) = .true.
        end select
      case (ST_ICE)
        select case (i)
        case (ST_ICE)
          same_ris (i,j) = .true.
        end select
      case (ST_ROCK)
        select case (i)
        case (ST_ROCK)
          same_ris (i,j) = .true.
        end select
      case default
        select case (i)
        case (ST_SEAWATER, ST_SEAICE, ST_ICE, ST_ROCK)
        case default
          same_ris (i,j) = .true.
        end select
      end select
    end do
    end do

    nullify                         (ind_atm)
    if (present (ind_soil)) nullify (ind_soil)
    if (present (ind_sea )) nullify (ind_sea )
    if (present (ind_land)) nullify (ind_land)
    allocate (stin  (gin % lbg(1) : gin % ubg(1),&
                     gin % lbg(2) : gin % ubg(2),&
                     gin % lbg(4) : gin % ubg(4)))
    allocate (stout (gout% lbg(1) : gout% ubg(1),&
                     gout% lbg(2) : gout% ubg(2),&
                     gout% lbg(4) : gout% ubg(4)))
    allocate (di    (gin % lbg(1) : gin % ubg(1),&
                     gin % lbg(2) : gin % ubg(2),&
                     gin % lbg(4) : gin % ubg(4)))
    allocate (do    (gout% lbg(1) : gout% ubg(1),&
                     gout% lbg(2) : gout% ubg(2),&
                     gout% lbg(4) : gout% ubg(4)))
    allocate (lsi   (gin % lbg(1) : gin % ubg(1),&
                     gin % lbg(2) : gin % ubg(2),&
                     gin % lbg(4) : gin % ubg(4)))
    allocate (lso   (gout% lbg(1) : gout% ubg(1),&
                     gout% lbg(2) : gout% ubg(2),&
                     gout% lbg(4) : gout% ubg(4)))

    !-----------------------------------------
    ! interpolation indices for the atmosphere
    ! (linear interpolation)
    !-----------------------------------------
    call hor_intp_coeff (gout, &
                         gin , &
                         ind_atm       )

    !--------------------------------------
    ! calculate interpolation coefficients
    !--------------------------------------
    if (l_soil) then
      !---------------------------
      ! interpolation indices soil
      !---------------------------
      if (.not. associated (gout% soiltyp))                 &
        call finish('hor_intp_coef_atm','destination soiltyp not present')
      if (.not. associated (gin % soiltyp))            &
        call finish('hor_intp_coef_atm','source soiltyp not present')
      if (.not. associated (ind_soil)) then
        stin  = nint (gin % soiltyp(:,:,1,:))
        stout = nint (gout% soiltyp(:,:,1,:))
        call hor_intp_coeff (gout, &
                             gin , &
                             ind_soil,     &
                     idxin = ind_atm,      &
                       sti = stin,         &
                       sto = stout,        &
                      stc1 = same_sty,     &
                      stc2 = same_sty      )
!!!                   stc2 = same_ris
      endif
    endif
    if (l_land .or. l_sea) then
      !------------------------------
      ! interpolation indices surface
      !------------------------------
      if (.not. associated (gout% lsm))                              &
        call finish('hor_intp_coef_atm','destination lsm not present')
      if (.not. associated (gin % lsm))                         &
        call finish('hor_intp_coef_atm','source lsm not present')
    endif

    if (l_land) then

      if (.not. associated (gout% fr_lake))                              &
        call finish('hor_intp_coef_atm','destination fr_lake not present')
      if (.not. associated (gin % fr_lake))                         &
        call finish('hor_intp_coef_atm','source fr_lake not present')
      if (.not. associated (gout% depth_lk))                              &
        call finish('hor_intp_coef_atm','destination depth_lk not present')
      if (.not. associated (gin % depth_lk))                         &
        call finish('hor_intp_coef_atm','source depth_lk not present')

      if (.not. associated (ind_land)) then
        if (associated (gin% fr_lake)) then
          lsi = gin% lsm(:,:,1,:) + gin% fr_lake(:,:,1,:) < fls(1)
        else
          lsi = gin% lsm(:,:,1,:) < fls(1)
        endif
        where (lsi)
          stin  = 1
          di    = 0._wp
        elsewhere     (gin % fr_lake (:,:,1,:) >= frl(1)&
                 .and. gin % depth_lk(:,:,1,:) >  0._wp )
          stin  = 2
          di    = log (gin % depth_lk(:,:,1,:))
        elsewhere
          stin  = 3
          di    = 0._wp
        endwhere

        if (associated (gout% fr_lake)) then
          lso = gout% lsm(:,:,1,:) + gout% fr_lake(:,:,1,:) < fls(2)
        else
          lso = gout% lsm(:,:,1,:) < fls(2)
        endif
        where (lso)
          stout = 1
          do    = 0._wp
        elsewhere     (gout% fr_lake (:,:,1,:) >= frl(2)&
                 .and. gout% depth_lk(:,:,1,:) >  0._wp )
          stout = 2
          do    = log (gout% depth_lk(:,:,1,:))
        elsewhere
          stout = 3
          do    = 0._wp
        endwhere
        call hor_intp_coeff (gout, &
                             gin , &
                             ind_land,     &
                     idxin = ind_atm,      &
                       sti = stin,         &
                       sto = stout,        &
                        di = di,           &
                        do = do,           &
                      stc1 = same_lsm,     &
                      stc2 = same_lsm      )
      endif
    endif

    if (l_sea) then
      if (.not. associated (ind_sea)) then
        if (associated (gin% fr_lake)) then
          lsi = gin% lsm(:,:,1,:) + gin% fr_lake(:,:,1,:) < flss(1)
        else
          lsi = gin% lsm(:,:,1,:) < flss(1)
        endif
         where (lsi)
          stin  = 1
          di    = 0._wp
!       elsewhere     (gin % fr_lake (:,:,1,:) >= frl(1)&
!                .and. gin % depth_lk(:,:,1,:) >  0._wp )
!         stin  = 2
!         di    = log (gin % depth_lk(:,:,1,:))
        elsewhere
          stin  = 3
          di    = 0._wp
        endwhere

        if (associated (gout% fr_lake)) then
          lso = gout% lsm(:,:,1,:) + gout% fr_lake(:,:,1,:) < flss(2)
        else
          lso = gout% lsm(:,:,1,:) < flss(2)
        endif
        where (lso)
          stout = 1
          do    = 0._wp
!       elsewhere     (gout% fr_lake (:,:,1,:) >= frl(2)&
!                .and. gout% depth_lk(:,:,1,:) >  0._wp )
!         stout = 2
!         do    = log (gout% depth_lk(:,:,1,:))
        elsewhere
          stout = 3
          do    = 0._wp
        endwhere
        call hor_intp_coeff (gout, &
                             gin , &
                             ind_sea,      &
                     idxin = ind_atm,      &
                       sti = stin,         &
                       sto = stout,        &
                        di = di,           &
                        do = do,           &
                      stc1 = same_lsm,     &
                      stc2 = same_lsm      )
      endif
    endif

  end subroutine hor_intp_coef_atm
!==============================================================================
  subroutine hor_intp_atm (atmin, atmout, fr_lake, fr_land)
  type (t_atm) ,intent(inout)        :: atmin
  type (t_atm) ,intent(inout)        :: atmout
  real(wp)     ,intent(in) ,optional :: fr_lake(2) ! lake fraction bound(in/out)
  real(wp)     ,intent(in) ,optional :: fr_land(2) ! land fraction bound(in/out)
  !-------------------------------------------
  ! horizontally interpolate atmospheric state
  !-------------------------------------------
    type (t_h_idx) ,pointer :: ind_atm  (:,:,:)     ! interp.coeffs. for atm.
    type (t_h_idx) ,pointer :: ind_soil (:,:,:)     ! interp.coeffs. for soil
    type (t_h_idx) ,pointer :: ind_surf (:,:,:)     ! interp.coeffs. for surf.
    type (t_h_idx) ,pointer :: ind      (:,:,:)     ! interpolation coeffs.
    logical                 :: l_soil               ! interpolate soil ?
    logical                 :: l_surf               !          surface ?
    real(wp)   ,allocatable :: x      (:,:,:,:)     ! temporary soil moisture
    integer                 :: m                    ! atmospheric field index
    integer                 :: k                    ! level index

    character(len=8), parameter :: list_soil (3) = &
      ['t_so    ','w_so    ','w_so_ice']
    character(len=8), parameter :: list_surf(19) = &
      ['freshsnw','t_snow  ','rho_snow','w_snow  ','h_snow  ','z0      ', &
       'tsurf   ','qv_s    ','w_i     ','t_ice   ','h_ice   ','fr_ice  ', &
       't_mnw_lk','t_wml_lk','h_ml_lk ','t_bot_lk','c_t_lk  ','skt     ', &
       'snowc   ' ]

    !--------------------------------------
    ! loop over fields in atmospheric state
    ! calculate interpolation coefficients
    !--------------------------------------
    l_soil = .false.
    l_surf = .false.
    do m = 1, size (atmout% m)
      !------------------------------------------------------
      ! look for presence of fields in source and destination
      !------------------------------------------------------
      if (.not. atmout% m(m)%i% alloc) cycle
      if (.not. atmin % m(m)%i% alloc) then
        call deallocate (atmout% m(m))
        cycle
      endif
      !---------------------------------------
      ! calculate interpolation indices,
      !---------------------------------------
      if      (any (list_soil == atmout% m(m)%i% name)) then
        l_soil = .true.
      else if (any (list_surf == atmout% m(m)%i% name)) then
        l_surf = .true.
      endif
    end do

    call hor_intp_coef_atm                                         &
      (atmin% grid, atmout% grid, l_soil, l_surf, .false., ind_atm,&
       ind_soil, ind_surf, fr_lake=fr_lake, fr_land=fr_land        )

    !--------------------------------------
    ! loop over fields in atmospheric state
    ! interpolate
    !--------------------------------------
    do m = 1, size (atmout% m)
      if (.not. atmout% m(m)%i% alloc) cycle

      if      (any (list_soil == atmout% m(m)%i% name)) then
        ind => ind_soil
      else if (any (list_surf == atmout% m(m)%i% name)) then
        ind => ind_surf
      else
        ind => ind_atm
      endif

      select case (atmout% m(m)%i% name)
      case ('w_so')
        !---------------------------------------------------------
        ! interpolate soil moisture index instead of soil moisture
        !---------------------------------------------------------
        allocate (x (atmin% grid% lbg(1):atmin% grid% ubg(1),&
                     atmin% grid% lbg(2):atmin% grid% ubg(2),&
                     1,                                      &
                     atmin% grid% lbg(4):atmin% grid% ubg(4)))
        do k = 1, atmin% grid% ns
          x = ind_wso (atmin% w_so (:,:,k:k,:), int(atmin% grid% soiltyp), k, 1)
          call hor_intp (atmout% m(m)% ptr (:,:,k:k,:), &
                         x,                             &
                         atmout% grid,                  &
                         atmin % grid,                  &
                         ind                            )
          atmout% m(m)% ptr (:,:,k:k,:) = wso_ind (atmout% m(m)% ptr (:,:,k:k,:), &
                                               int(atmout% grid% soiltyp), k, 1   )
        end do
        deallocate (x)
      case ('w_so_ice')
        !------------------------------------------------------
        ! interpolate soil ice normalised by soil water content
        !------------------------------------------------------
        if (.not.associated (atmout% w_so))                       &
          call finish('hor_intp_atm','atmout% w_so not associated')
        if (.not.associated (atmin % w_so))                       &
          call finish('hor_intp_atm','atmin % w_so not associated')
        allocate (x (atmin% grid% lbg(1):atmin% grid% ubg(1),&
                     atmin% grid% lbg(2):atmin% grid% ubg(2),&
                     1,                                      &
                     atmin% grid% lbg(4):atmin% grid% ubg(4)))
        do k = 1, atmin% grid% ns
          x = atmin% w_so_ice (:,:,k:k,:) / max (atmin% w_so (:,:,k:k,:), 1.e-10_wp)
          call hor_intp (atmout% m(m)% ptr (:,:,k:k,:), &
                         x,                             &
                         atmout% grid,                  &
                         atmin % grid,                  &
                         ind                            )
          atmout% m(m)% ptr (:,:,k:k,:) = atmout% m(m)% ptr  (:,:,k:k,:) &
                                        * atmout%       w_so (:,:,k:k,:)
        end do
        deallocate (x)
      case default
        !---------------------------------
        ! interpolate all other quantities
        !---------------------------------
        call hor_intp (atmout% m(m)% ptr, &
                       atmin % m(m)% ptr, &
                       atmout% grid,      &
                       atmin % grid,      &
                       ind                )
      end select
    end do

    if (associated(ind_atm )) deallocate (ind_atm )
    if (associated(ind_soil)) deallocate (ind_soil)
    if (associated(ind_surf)) deallocate (ind_surf)
    call set_pointers (atmout)

    !----------------------------------------------------
    ! ensure consistency for interpolated Flake variables
    !----------------------------------------------------
    if (associated (atmout% t_mnw_lk) .and. &
        associated (atmout% t_wml_lk) .and. &
        associated (atmout% h_ml_lk ) .and. &
        associated (atmout% t_bot_lk) .and. &
        associated (atmout% c_t_lk  )       ) &
      call flake_init_atm (atmout, fr_lake(2))

  end subroutine hor_intp_atm
!------------------------------------------------------------------------------
  subroutine flake_init_atm (a, fr_thr)
  type (t_atm) ,intent(inout) :: a       ! atmospheric state
  real(wp)     ,intent(in)    :: fr_thr  ! threshold for Flake grid-points
  !-------------------------------------------------------------
  ! interface to flake_init
  ! for consistency check after interpolation of Flake variables
  !-------------------------------------------------------------
    type(t_grid) ,pointer :: g
    integer               :: nflkgb       ! number of Flake grid-points
    logical  ,allocatable :: m  (:,:,:,:) ! mask for valid Flake points
    real(wp) ,allocatable :: fr_lake  (:) ! invariant meta data
    real(wp) ,allocatable :: depth_lk (:) ! invariant meta data
    real(wp) ,allocatable :: fetch_lk (:) ! not used
    real(wp) ,allocatable :: dp_bs_lk (:) ! not used
    real(wp) ,allocatable :: t_bs_lk  (:) ! not used
    real(wp) ,allocatable :: gamso_lk (:) ! not used
    real(wp) ,allocatable :: t_snow   (:) ! set consistently with Flake
    real(wp) ,allocatable :: h_snow   (:) ! set consistently with Flake
    real(wp) ,allocatable :: t_ice    (:) ! set consistently with Flake
    real(wp) ,allocatable :: h_ice    (:) ! set consistently with Flake
    real(wp) ,allocatable :: t_mnw_lk (:) ! prognostic Flake variables
    real(wp) ,allocatable :: t_wml_lk (:) ! prognostic Flake variables
    real(wp) ,allocatable :: t_bot_lk (:) ! prognostic Flake variables
    real(wp) ,allocatable :: c_t_lk   (:) ! prognostic Flake variables
    real(wp) ,allocatable :: h_ml_lk  (:) ! prognostic Flake variables
    real(wp) ,allocatable :: t_b1_lk  (:) ! not used
    real(wp) ,allocatable :: h_b1_lk  (:) ! not used
    real(wp) ,allocatable :: t_scf_lk (:) ! not used

    logical  ,allocatable :: use_iceanalysis (:) ! dummy parameters,
    real(wp) :: fr_ice (0)                       ! currently not used

    g => a% grid
    allocate (m (g% shape(1),g% shape(2),1,g% shape(4)))
    m      = g% fr_lake >= fr_thr
    nflkgb = count (m)
    allocate (fr_lake         (nflkgb))
    allocate (depth_lk        (nflkgb))
    allocate (fetch_lk        (nflkgb))
    allocate (dp_bs_lk        (nflkgb))
    allocate (t_bs_lk         (nflkgb))
    allocate (gamso_lk        (nflkgb))
    allocate (t_snow          (nflkgb))
    allocate (h_snow          (nflkgb))
    allocate (t_ice           (nflkgb))
    allocate (h_ice           (nflkgb))
    allocate (t_mnw_lk        (nflkgb))
    allocate (t_wml_lk        (nflkgb))
    allocate (t_bot_lk        (nflkgb))
    allocate (c_t_lk          (nflkgb))
    allocate (h_ml_lk         (nflkgb))
    allocate (t_b1_lk         (nflkgb))
    allocate (h_b1_lk         (nflkgb))
    allocate (t_scf_lk        (nflkgb))
    allocate (use_iceanalysis (nflkgb))

    h_snow          =   0._wp
    t_snow          = 273._wp
    h_ice           =   0._wp
    t_ice           = 273._wp
    use_iceanalysis = .false.

    if (associated (a% t_snow)) t_snow = pack (a% t_snow   ,m)
    if (associated (a% h_snow)) h_snow = pack (a% h_snow   ,m)
    if (associated (a% t_ice )) t_ice  = pack (a% t_ice    ,m)
    if (associated (a% h_ice )) h_ice  = pack (a% h_ice    ,m)
    fr_lake                            = pack (g% fr_lake  ,m)
    depth_lk                           = pack (g% depth_lk ,m)
    t_mnw_lk                           = pack (a% t_mnw_lk ,m)
    t_wml_lk                           = pack (a% t_wml_lk ,m)
    t_bot_lk                           = pack (a% t_bot_lk ,m)
    c_t_lk                             = pack (a% c_t_lk   ,m)
    h_ml_lk                            = pack (a% h_ml_lk  ,m)

    call flake_init (nflkgb, use_iceanalysis,                &
                     fr_lake, depth_lk, fr_ice,              &
                     fetch_lk, dp_bs_lk, t_bs_lk, gamso_lk,  &
                     t_snow,   h_snow,                       &
                     t_ice, h_ice,                           &
                     t_mnw_lk, t_wml_lk, t_bot_lk,           &
                     c_t_lk, h_ml_lk,                        &
                     t_b1_lk, h_b1_lk,                       &
                     t_scf_lk                                )

    if (associated (a% t_snow)) a% t_snow = unpack (t_snow   ,m, a% t_snow  )
    if (associated (a% h_snow)) a% h_snow = unpack (h_snow   ,m, a% h_snow  )
    if (associated (a% t_ice )) a% t_ice  = unpack (t_ice    ,m, a% t_ice   )
    if (associated (a% h_ice )) a% h_ice  = unpack (h_ice    ,m, a% h_ice   )
    a% t_mnw_lk                           = unpack (t_mnw_lk ,m, a% t_mnw_lk)
    a% t_wml_lk                           = unpack (t_wml_lk ,m, a% t_wml_lk)
    a% t_bot_lk                           = unpack (t_bot_lk ,m, a% t_bot_lk)
    a% c_t_lk                             = unpack (c_t_lk   ,m, a% c_t_lk  )
    a% h_ml_lk                            = unpack (h_ml_lk  ,m, a% h_ml_lk )

  end subroutine flake_init_atm
!------------------------------------------------------------------------------
  subroutine vert_intp (atmin, atmout, lin, cubic, inplace)
  type (t_atm) ,intent(inout) TARGET :: atmin   ! input  atmospheric state
  type (t_atm) ,intent(inout) TARGET :: atmout  ! output atmospheric state
  logical      ,intent(in) ,optional :: lin     ! linear interpolation flag
  logical      ,intent(in) ,optional :: cubic   ! Cubic  interpolation flag
  logical      ,intent(in) ,optional :: inplace ! deallocate source
  !----------------------------------------------------------------
  ! vertically interpolate atmospheric state to new pressure levels
  !----------------------------------------------------------------
    real(wp) ,allocatable :: lnphi (:,:,:,:), lnpho (:,:,:,:) ! half-level pres.
    real(wp) ,allocatable :: lnpfi (:,:,:,:), lnpfo (:,:,:,:) ! full-level pres.
    real(wp) ,allocatable :: d2f   (:,:,:,:), d2h   (:,:,:,:)
    real(wp) ,allocatable :: dfh   (:,:,:),   df    (:,:,:)
    real(wp) ,allocatable :: a     (:,:,:),   b     (:,:,:)
    real(wp) ,allocatable :: temp  (:,:,:,:) ! to adjust reference pressure
    real(wp)     ,pointer :: pti   (:,:,:,:)
    real(wp)     ,pointer :: pto   (:,:,:,:)
    type(t_grid) ,pointer :: gi, go

    integer  :: l1,l2,l3,l4,u1,u2,u3hi,u3fi,u3ho,u3fo,u4
    integer  :: k, kei, l, luv, nzi, keo, nzo
    logical  :: li, linp, lcub
    logical  :: ltm        ! T is missing
    logical  :: lgeo       ! Current field to interpolate is geof or geoh
    logical  :: lintz      ! Interpolate over height (geopotential) not ln(p)
    logical  :: lghi, lgfi ! atmin%  geoh,geof allocated ?
    logical  :: lgho, lgfo ! atmout% geoh,geof allocated ?
    !---------------------------------------------------------------------
    ! Upper temperature bound for geopotential extrapolation below surface
    ! cf. IFS documentation (CY25R1), Part II, Paragraph 5.3.2
    !---------------------------------------------------------------------
    real(wp) ,parameter   :: Tx = 290.5_wp
    !--------------------
    ! derive array bounds
    !--------------------
    gi  => atmin % grid
    go  => atmout% grid
    l1   = atmin % lb(1); u1 = atmin% ub(1)
    l2   = atmin % lb(2); u2 = atmin% ub(2)
    l3   = atmin % lb(3)
    l4   = atmin % lb(4); u4 = atmin% ub(4)
    u3fi = atmin % ub(3)
    u3hi = atmin % ub(3) ! + 1
    u3fo = atmout% ub(3)
    u3ho = atmout% ub(3) ! + 1
    kei  = u3fi - l3 + 1
    keo  = u3fo - l3 + 1
    li   = .false.; if (present(lin))     li   = lin
    linp = .false.; if (present(inplace)) linp = inplace
    lcub = .false.; if (present(cubic))   lcub = cubic
    lintz= .false.
    !----------------------
    ! allocate local arrays
    !----------------------
    allocate   (lnphi (l1:u1,l2:u2,l3:u3hi,l4:u4))
    allocate   (lnpho (l1:u1,l2:u2,l3:u3ho,l4:u4))
    allocate   (lnpfi (l1:u1,l2:u2,l3:u3fi,l4:u4))
    allocate   (lnpfo (l1:u1,l2:u2,l3:u3fo,l4:u4))
    allocate   (d2f   (l1:u1,l2:u2,l3:u3fi,l4:u4))
    allocate   (d2h   (l1:u1,l2:u2,l3:u3hi,l4:u4))
    allocate   (dfh   (l1:u1,l2:u2,        l4:u4))
    allocate   (df    (l1:u1,l2:u2,        l4:u4))
    allocate   (a     (l1:u1,l2:u2,        l4:u4))
    allocate   (b     (l1:u1,l2:u2,        l4:u4))
    if (gi% vct == VCT_P_HYB .and. &
        go% vct == VCT_P_HYB       ) then
      allocate (temp  (l1:u1,l2:u2,      1,    1))
    else
      allocate (temp  (l1:u1,l2:u2,l3:u3fi,    1))
    endif
    !----------------------------------
    ! Surface gradient of geopotential:
    ! fallback value for t missing.
    ! See also explanation below.
    !----------------------------------
    dfh = -R * Tx
    ltm = .true.
    b   = lapse_cl                              ! Lapse rate gamma [K/m]
    a   = b * R / atmin% grid% g                ! Exponent a=R*gamma/g
    !--------------------------------
    ! loop over fields to interpolate
    !--------------------------------
    ! loop: 1: u,v on Arakawa-C-grid
    !       2: everything else
    !--------------------------------
    if (gi% arakawa    == 'C' .and. &
        atmin% gp_wind /= 'A' .and. &
       (gi% vct  /= VCT_P_ISO .or.  &
        go% vct  /= VCT_P_ISO     ) ) then
      luv = 1
    else
      luv = 2
    endif
    !------------------------------------------------------------
    ! derive vertical coordinates (ln p) for half and full levels
    ! (Arakawa A grid points)
    !------------------------------------------------------------
    if (gi% vct == VCT_P_HYB .or. &
        gi% vct == VCT_Z_HYB .or. &
        gi% vct == VCT_Z_GEN      ) call set_p(atmin)

    if (go% vct == VCT_P_HYB      ) call set_p(atmout)

    if (go% vct == VCT_Z_GEN .or. &
        go% vct == VCT_Z_HYB      ) lintz = .true.

    if (lintz) then
      lghi = associated (atmin % geoh)
      lgfi = associated (atmin % geof)
      lgho = associated (atmout% geoh)
      lgfo = associated (atmout% geof)
      if (.not. (lghi .and. lgfi)) call set_geo (atmin , geof= .not. lgfi)
      if (.not. (lgho .and. lgfo)) call set_geo (atmout, geof= .not. lgfo)
    end if

    do l=1,2
     if (l==1 .and. (atmout% grid% arakawa /= 'C' &
                     .or.  atmout% gp_wind == 'A')) cycle

     call init_spline1 (atmin,  lnphi, lnpfi, lgeo=lintz)
     call init_spline1 (atmout, lnpho, lnpfo, lgeo=lintz)
     !-----------------
     ! loop over fields
     !-----------------
     do k=1,size(atmout%m)
      !----------------------------------------------------
      ! keep vertical coordinate dependent fields untouched
      !----------------------------------------------------
      if (.not. atmin% m(k)%i% alloc) cycle
      if (atmout% m(k)%i% ref)        cycle
      if (linp) call allocate (atmout,atmin% m(k)%i%name)
      if (.not. atmout% m(k)%i% alloc) cycle
      !-------------------------------------------
      ! special handling for u,v on Arakawa-C-grid
      !-------------------------------------------
      select case (atmin% m(k)%i% name)
      case ('u')
        if (l/=luv) cycle
        if (l==1) then
          select case (gi% vct)
          case (VCT_P_ISO)
          case (VCT_P_HYB)
            temp(:,:,1:1,:) = atmin%psr
            call half_gp_w (atmin% psr, atmin% grid)
            call set_p(atmin)
            call init_spline1(atmin, lnphi, lnpfi)
            atmin%psr =  temp(:,:,1:1,:)
          case default
            if (lintz) then
              temp = atmin% geof
              call half_gp_w (atmin% geof, atmin% grid)
              call init_spline1 (atmin, lnphi, lnpfi, lgeo=lintz)
              atmin% geof =  temp
            else
              temp = atmin% pf
              call half_gp_w (atmin% pf, atmin% grid)
              call init_spline1 (atmin, lnphi, lnpfi)
              atmin% pf =  temp
            endif
          end select
          select case (go% vct)
          case (VCT_P_ISO)
          case (VCT_P_HYB)
            temp(:,:,1:1,:) = atmout%psr
            call half_gp_w (atmout% psr, atmout% grid)
            call set_p(atmout)
            call init_spline1(atmout, lnpho, lnpfo)
            atmout%psr = temp(:,:,1:1,:)
          case default
            if (lintz) then
              temp = atmout% geof
              call half_gp_w (atmout% geof, atmout% grid)
              call init_spline1 (atmout, lnpho, lnpfo, lgeo=lintz)
              atmout% geof = temp
            else
              temp = atmout% pf
              call half_gp_w (atmout% pf, atmout% grid)
              call init_spline1 (atmout, lnpho, lnpfo)
              atmout% pf = temp
            endif
          end select
        endif
      case ('v')
        if (l/=luv) cycle
        if (l==1) then
          select case (gi% vct)
          case (VCT_P_ISO)
          case (VCT_P_HYB)
            temp(:,:,1:1,:) = atmin%psr
            call half_gp_s (atmin% psr, atmin% grid)
            call set_p(atmin)
            call init_spline1(atmin, lnphi, lnpfi)
            atmin%psr = temp(:,:,1:1,:)
          case default
            if (lintz) then
              temp = atmin% geof
              call half_gp_s (atmin% geof, atmin% grid)
              call init_spline1 (atmin, lnphi, lnpfi, lgeo=lintz)
              atmin% geof =  temp
            else
              temp = atmin% pf
              call half_gp_s (atmin% pf, atmin% grid)
              call init_spline1 (atmin, lnphi, lnpfi)
              atmin% pf = temp
            endif
          end select
          select case (go% vct)
          case (VCT_P_ISO)
          case (VCT_P_HYB)
            temp(:,:,1:1,:) = atmout%psr
            call half_gp_s (atmout% psr, atmout% grid)
            call set_p(atmout)
            call init_spline1(atmout, lnpho, lnpfo)
            atmout%psr = temp(:,:,1:1,:)
          case default
            if (lintz) then
              temp = atmout% geof
              call half_gp_s (atmout% geof, atmout% grid)
              call init_spline1 (atmout, lnpho, lnpfo, lgeo=lintz)
              atmout% geof = temp
            else
              temp = atmout% pf
              call half_gp_s (atmout% pf, atmout% grid)
              call init_spline1 (atmout, lnpho, lnpfo)
              atmout%pf = temp
            endif
          end select
        endif
      case default
        if (l==1) cycle
      end select
      !---------------------------------------------
      ! deallocate output fields if no input present
      !---------------------------------------------
      if (.not. atmin % m(k)%i% alloc) then
        call deallocate (atmout% m(k))
        cycle
      endif
      nzi = atmin % m(k)%i% ub(3)-atmin % m(k)%i% lb(3)+1
      nzo = atmout% m(k)%i% ub(3)-atmout% m(k)%i% lb(3)+1

!!$      if (atmout%grid%arakawa == 'C') then
!!$        if (atmout%m(k)%i%name == 'u') then
!!$        else if (atmout%m(k)%i%name == 'v') then
!!$        endif
!!$      endif

      !----------------------------------------------------------
      ! Set vertical gradient at the surface.
      ! Surface derivatives are estimated from an atmosphere with
      ! uniform (climatological?) lapse rate gamma=-dT/dz with
      !   z(p) = z(p_0) + T_0/gamma * [1 - (p/p_0)^(R*gamma/g)],
      !   d(g*z)/dln(p)|_{p=p_0} = -R*T_0,
      !   dT/dln(p)|_{p=p_0} = -gamma*dz/dln(p) = R*gamma/g*T_0.
      ! (Climatologically, R*gamma/g ~ 0.19).
      !----------------------------------------------------------
      lgeo = atmin% m(k)%i% name(1:3) == 'geo'
      if (.not. lintz) then
       if (atmin% m(k)%i% name == 't') then
!        a   = atmin % m(k)% ptr(:,:,u3fi,:) &
!            + b * atmin% grid% geosp(atmin% lb(1):atmin% ub(1),     &
!                                     atmin% lb(2):atmin% ub(2),1,:) &
!                / atmin% grid% g
!        where (a > 290.5) &
!          b = (290.5 - a) * atmin% grid% g / atmin% grid% geosp(:,:,1,:)
!        a   = b * R / atmin% grid% g
        df  = lapse_cl * R / atmin% grid% g * atmin % m(k)% ptr(:,:,u3fi,:)
        dfh = -R * atmin % m(k)% ptr(:,:,u3fi,:)    ! Geop. surf.gradient, keep
        ltm = .false.
       else if (lgeo) then
        if (ltm) then
!          write (0,*) "vert_intp: WARNING: t not found, using fallback dfh"
        end if
        df  = dfh
       else
        df  = 0._wp
       endif
      else ! lintz==.true.
      !----------------------------------------------------------
      ! Vertical interpolation using height (-geopotential):
      ! use climatological vertical gradient for T extrapolation,
      !   dT/d(-g*z) = gamma/g,
      !   dln(p)/d(-g*z)|{p=p_0} = 1 / (R*T_0)
      !----------------------------------------------------------
       if (atmin% m(k)%i% name == 't') then
        df  = lapse_cl / atmin% grid% g
        dfh = (1/R) / atmin % m(k)% ptr(:,:,u3fi,:)
        ltm = .false.
       else if (atmin% m(k)%i% name == 'pf') then
        if (ltm) then
!        write (0,*) "vert_intp: WARNING: t not found, using fallback dfh"
         df = -1 / dfh
        else
         df = dfh
        end if
        !----------------------------------
        ! Keep pf, apply log to input field
        !----------------------------------
        temp = atmin % m(k)% ptr
        atmin % m(k)% ptr = log (temp)
       else if (atmin% m(k)%i% name == 'ph') then
        if (ltm) then
!        write (0,*) "vert_intp: WARNING: t not found, using fallback dfh"
         df = -1 / dfh
        else
         df = dfh
        end if
        !----------------------------------
        ! Keep ph, apply log to input field
        !----------------------------------
        temp = atmin % m(k)% ptr(:,:,2:,:)
        atmin % m(k)% ptr(:,:,2:,:) = log (temp)
       else if (lgeo) then
        cycle           ! Geopotential already present in output, skip
       else
        df  = 0._wp
       end if
      end if ! lintz
      !------------------------------
      ! just copy single level fields
      !------------------------------
      if (nzi==1) then
        atmout% m(k)% ptr = atmin% m(k)% ptr
      !--------------------------------
      ! ke vertical levels: interpolate
      !--------------------------------
      else if (nzi==kei) then
        if (nzo==keo) then
          call init_spline2 (atmin% m(k)% ptr, lnpfi, d2f, df)
          if (li) d2f = 0._wp
          call vintp_spline (lnpfi, atmin % m(k)% ptr, d2f, &
                             lnpfo, atmout% m(k)% ptr,      &
                             df, a, lgeo, lcub)
        else if (nzo==keo+1) then
          if (go% vct /= VCT_P_HYB) &
            call finish ('vert_intp','nzo==keo+1, COSMO')
          call init_spline2 (atmin% m(k)% ptr, lnpfi, d2f, df)
          if (li) d2f = 0._wp
          call vintp_spline (lnpfi, atmin % m(k)% ptr, d2f, &
                             lnpho, atmout% m(k)% ptr,      &
                             df, a, lgeo, lcub)
        else
          call deallocate (atmout% m(k))
        endif
      !----------------------------------
      ! ke+1 vertical levels: interpolate
      !----------------------------------
      else if (nzi==kei+1) then
!       if (gi% vct /= VCT_P_HYB .and. &
!           gi% vct /= VCT_Z_GEN     ) &
!         call finish ('vert_intp','nzi==kei+1, COSMO')
        pti => atmin % m(k)% ptr (:,:,2:,:)
        call init_spline2 (pti, lnphi, d2h, df)
        if (li) d2h = 0._wp
        if (nzo==keo) then
          call vintp_spline (lnphi,               pti, d2h, &
                             lnpfo, atmout% m(k)% ptr,      &
                             df, a, lgeo, lcub)
        else if (nzo==keo+1) then
          pto => atmout% m(k)% ptr (:,:,2:,:)
          call vintp_spline (lnphi,               pti, d2h, &
                             lnpho,               pto,      &
                             df, a, lgeo, lcub)
          atmout% m(k)% ptr (:,:,1,:) = atmin % m(k)% ptr (:,:,1,:)
        else
          call deallocate (atmout% m(k))
        endif
      else
      !-------------------------------------------------------------
      ! unexpected number of levels, copy if vertical extent matches
      !-------------------------------------------------------------
        if (nzi==nzo) then
          atmout% m(k)% ptr = atmin% m(k)% ptr
        else
          call deallocate (atmout% m(k))
        endif
      end if
      if (lintz) then
        if (atmin% m(k)%i% name == "pf") then
          !-----------------------------
          ! Undo logarithm applied to pf
          !-----------------------------
          if (.not. linp) atmin% m(k)% ptr = temp
          atmout% m(k)% ptr = exp (atmout% m(k)% ptr)
        else if (atmin% m(k)%i% name == "ph") then
          !-----------------------------
          ! Undo logarithm applied to ph
          !-----------------------------
          if (.not. linp) atmin% m(k)% ptr(:,:,2:,:) = temp
          atmout% m(k)% ptr(:,:,2:,:) = exp (atmout% m(k)% ptr(:,:,2:,:))
        end if
      end if
      if (linp) call deallocate (atmin, atmin% m(k)%i% name)
     end do
    end do

    !---------
    ! clean up
    !---------
    if (lintz) then
      if (.not. lghi) call deallocate (atmin,  'geoh')
      if (.not. lgfi) call deallocate (atmin,  'geof')
      if (.not. lgho) call deallocate (atmout, 'geoh')
      if (.not. lgfo) call deallocate (atmout, 'geof')
    endif
    call set_pointers (atmout)
    deallocate (lnphi, lnpho, lnpfi, lnpfo, d2f, d2h, dfh, df, a, b)

  end subroutine vert_intp
!------------------------------------------------------------------------------
  subroutine init_spline1 (atm, lnph, lnpf, lgeo)
  !-----------------------------
  ! derives argument grid (ln p)
  !-----------------------------
  type (t_atm)      ,intent(inout) :: atm
  real(wp)          ,intent(out)   :: lnph (:,:,:,:)
  real(wp)          ,intent(out)   :: lnpf (:,:,:,:)
  logical, optional ,intent(in)    :: lgeo      ! Use (-geopot.) not ln(p)

    type (t_grid) ,pointer :: g
    logical                :: lz

    lz = .false.; if (present (lgeo)) lz = lgeo

    g => atm% grid
    if (size (lnph,3) /= g% nz  ) call finish ('init_spline1','dim3/=nz')
    if (size (lnpf,3) /= g% nz  ) call finish ('init_spline1','dim3/=nz')

    if (lz) then
       if (associated (atm% geof)) then
          lnpf = - atm% geof
       else
          call finish ('init_spline1','geof not allocated')
       end if
       if (associated (atm% geoh)) then
          lnph = - atm% geoh(:,:,2:,:)
       else
          lnph = 0
          select case (atm% grid% vct)
          case (VCT_P_ISO)
             ! Probably ok
          case default
             call finish ('init_spline1','geoh not allocated')
          end select
       end if
       return
    end if

    select case (atm% grid% vct)
    case default
      call finish ('init_spline1','unknown leveltype')
    case (VCT_P_HYB)
      !------------------------------------
      ! GME/HRM hybrid pressure coordinates
      !------------------------------------
      if (.not. associated (atm% ph) .or. .not. associated (atm% pf)) &
        call set_p (atm)

      lnpf           = log (atm% pf)
      lnph(:,:,:,:)  = log (atm% ph (:,:,2:,:))
!     lnph(:,:,1 ,:) = log (atm% ph (:,:,2 ,:)/2._wp)

!!$    case (VCT_Z_HYB)
!!$      !---------------------------
!!$      ! COSMO hybrid z-coordinates
!!$      !---------------------------
!!$      if (.not. associated (atm% pf)) call set_p (atm)
!!$      lnpf = log (atm% pf)
!!$      lnph = 0._wp         ! ke+1 currently not supported

    case (VCT_Z_GEN, VCT_Z_HYB)
      !------------------------------------------------------------
      ! ICON height-based hybrid or generalized vertical coordinate
      ! COSMO hybrid z-coordinate (as used e.g. by MeteoSwiss)
      !------------------------------------------------------------
      if (.not. associated (atm% ph) .or. .not. associated (atm% pf)) &
        call set_p (atm)

      lnpf           = log (atm% pf)
      lnph(:,:,:,:)  = log (atm% ph (:,:,2:,:))

    case (VCT_P_ISO)
      !---------------------
      ! pressure coordinates
      !---------------------
      if (.not. associated (atm% pf)) call set_pf (atm)

      lnpf           = log (atm% pf)
      lnph           = lnpf

    end select

  end subroutine init_spline1
!==============================================================================
  subroutine average_pole_scalar (atm, s, nn)
  type (t_atm) ,intent(in)    :: atm              ! atmospheric state
  real(wp)     ,pointer       :: s (:,:,:,:)      ! scalar quantity
  integer      ,intent(in)    :: nn               ! neigbour
  !------------------------------------------------------
  ! average the scalar quantity S at the poles
  ! taking the 5 neighbour points on the icosahedral grid
  ! nn  = 1: next   neigbour
  ! nn >= 2: second neigbours
  !------------------------------------------------------
    integer  :: i,j, d, n
    real(wp) :: tn (size(s,3)), ts (size(s,3))
    !--------------------------------
    ! return if no action is required
    !--------------------------------
    if (nn <= 0)                                 return  ! no averaging
    if (atm% grid% gridtype /= DWD6_ICOSAHEDRON) return  ! no icosahedron
    !----------------------------
    ! determine indices of source
    !----------------------------
    select case (nn)
    case (1)
      i=1; j=1
    case (2:)
      i=1; j=2
    end select
    !-----------
    ! average NP
    !-----------
    tn = 0
    n  = 0
    do d=1,5
      if (atm% grid% marr(1,i,j,d) == dace% pe) then
        n = n + 1
        tn = tn + s(atm% grid% marr(2,i,j,d), &
                    atm% grid% marr(3,i,j,d), &
                    :                       , &
                    atm% grid% marr(4,i,j,d))
      endif
    end do
    n  = p_sum (n)
    tn = p_sum (tn) / 5
    if (n/=5) call finish ('average_pole_scalar','n/=5')
    if (all (atm% lb(1:2) == (/0,1/))) then
      do d=1,5
        s(0,1,:,d) = tn
      end do
    endif
    !-----------
    ! average SP
    !-----------
    ts = 0
    n  = 0
    do d=6,10
      if (atm% grid% marr(1,i,j,d) == dace% pe) then
        n = n + 1
        ts = ts + s(atm% grid% marr(2,i,j,d), &
                    atm% grid% marr(3,i,j,d), &
                    :                       , &
                    atm% grid% marr(4,i,j,d))
      endif
    end do
    n  = p_sum (n)
    ts = p_sum (ts) / 5
    if (n/=5) call finish ('average_pole_scalar','n/=5')
    if (all (atm% lb(1:2) == (/0,1/))) then
      do d=6,10
        s(0,1,:,d) = ts
      end do
    endif
    if (nn >= 2) then
      !---------------------------------
      ! set next neighbour at NP as well
      !---------------------------------
      i=1; j=1
      n = 0
      do d=1,5
        if (atm% grid% marr(1,i,j,d) == dace% pe) then
          n = n + 1
          s(atm% grid% marr(2,i,j,d), &
            atm% grid% marr(3,i,j,d), &
            :                       , &
            atm% grid% marr(4,i,j,d)) = tn
        endif
      end do
      n = p_sum (n)
      if (n/=5) call finish ('average_pole_scalar','n/=5')
      !---------------------------------
      ! set next neighbour at SP as well
      !---------------------------------
      i=1; j=1
      n = 0
      do d=6,10
        if (atm% grid% marr(1,i,j,d) == dace% pe) then
          n = n + 1
          s(atm% grid% marr(2,i,j,d), &
            atm% grid% marr(3,i,j,d), &
            :                       , &
            atm% grid% marr(4,i,j,d)) = ts
        endif
      end do
      n = p_sum (n)
      if (n/=5) call finish ('average_pole_scalar','n/=5')
    end if

  end subroutine average_pole_scalar
!------------------------------------------------------------------------------
  subroutine average_pole_wind (atm, n)
  type (t_atm) ,intent(inout) :: atm             ! atmospheric state
  integer      ,intent(in)    :: n               ! number of points
  !------------------------------------------------------
  ! average the wind components at the poles
  ! taking the 5 neighbour points on the icosahedral grid
  !------------------------------------------------------
!   real(wp) :: s, c           ! sin, cos of the rotation angle
!   real(wp) :: a              ! rotation angle
!   real(wp) :: u (atm% ub(3)) ! old wind component u
!   real(wp) :: v (atm% ub(3)) ! old wind component v
!   integer  :: l              ! diamond index
!   integer  :: i              ! horizontal index

    if (atm% grid% gridtype /= DWD6_ICOSAHEDRON) return
    if (n<=0)                                    return

!   if(dace% lpio) write(6,*) repeat('-',79)
    if(dace% lpio) write(6,'(/a,i2)') '  average_pole_wind: n =',n
    !--------------------------------------------
    ! rotate winds near the pole to global system
    !--------------------------------------------
    select case (n)
    case (1)
      call rotate_wind_global (atm,1,1)
    case (2:)
      call rotate_wind_global (atm,1,2)
    end select
    !----------------------
    ! average scalar fields
    !----------------------
    call average_pole_scalar (atm, atm% u, n)
    call average_pole_scalar (atm, atm% v, n)
    !------------------------------------------------
    ! rotate back winds near the pole to local system
    !------------------------------------------------
    call rotate_wind_local   (atm,1,1)
    select case (n)
    case (2:)
      call rotate_wind_local (atm,1,2)
    end select

  end subroutine average_pole_wind
!-------------------------------------------------------------------------------
  pure subroutine rotate_pole_global1 (atms)
    !------------------------------------------
    ! Rotate wind at the poles (vector version)
    ! (workaround for bug in gfortran, nfort).
    !------------------------------------------
    type(t_atm) ,intent(inout) :: atms (:)
    integer                    :: i        ! loop index
    do i = 1, size (atms)
      call rotate_pole_global (atms (i))
    end do
  end subroutine rotate_pole_global1
!------------------------------------------------------------------------------
  pure subroutine rotate_pole_global (atm)
  type (t_atm) ,intent(inout) :: atm
  !---------------------------------------------------------------------------
  ! Rotate the winds at the poles to the 'global' coordinate system:
  !   u is orientated to the east if the pole is approached on the 0-meridian.
  !   For the north pole this corresponds with the local coordinate system
  !   of diamond 1.
  ! The 'local'  coordinate system corresponds with the GME convention.
  ! The 'global' coordinate system corresponds with the 3D-Var convention.
  !---------------------------------------------------------------------------
    real(wp) :: s, c           ! sin, cos of the rotation angle
    real(wp) :: a              ! rotation angle
    real(wp) :: u (atm% ub(3)) ! old wind component u
    real(wp) :: v (atm% ub(3)) ! old wind component v
    integer  :: l              ! diamond index
    !---------------------------------
    ! check if this pe holds the poles
    !---------------------------------
    if (all (atm% lb(1:2) == (/0,1/))) then
      !------------
      ! north poles
      !------------
      do l=1,5
        a = 2._wp * pi * (l-1._wp) / 5
        c = cos(a)
        s = sin(a)
        u = atm% u (0,1,:,l)
        v = atm% v (0,1,:,l)
        atm% u(0,1,:,l) = c * u - s * v
        atm% v(0,1,:,l) = c * v + s * u
      end do
      !------------
      ! south poles
      !------------
      do l=6,10
        a = 2._wp * pi * (l-5.5_wp) / 5
        c = cos(a)
        s = sin(a)
        u = atm% u (0,1,:,l)
        v = atm% v (0,1,:,l)
        atm% u(0,1,:,l) = c * u + s * v
        atm% v(0,1,:,l) = c * v - s * u
      end do
    end if
  end subroutine rotate_pole_global
!-------------------------------------------------------------------------------
  pure subroutine rotate_pole_local1 (atms)
    !------------------------------------------
    ! Rotate wind at the poles (vector version)
    ! (workaround for bug in gfortran, nfort).
    !------------------------------------------
    type(t_atm) ,intent(inout) :: atms (:)
    integer                    :: i        ! loop index
    do i = 1, size (atms)
      call rotate_pole_local (atms (i))
    end do
  end subroutine rotate_pole_local1
!------------------------------------------------------------------------------
  pure subroutine rotate_pole_local (atm)
  type (t_atm) ,intent(inout) :: atm
  !---------------------------------------------------------------------------
  ! Rotate the winds at the poles to the 'local' coordinate system:
  !   u is orientated to the east if the pole is approached from the
  !   midpoint of the diamond. For the north pole this corresponds with
  !   the global coordinate system on diamond 1
  !---------------------------------------------------------------------------
    real(wp) :: s, c           ! sin, cos of the rotation angle
    real(wp) :: a              ! rotation angle
    real(wp) :: u (atm% ub(3)) ! old wind component u
    real(wp) :: v (atm% ub(3)) ! old wind component v
    integer  :: l              ! diamond index
    !---------------------------------
    ! check if this pe holds the poles
    !---------------------------------
    if (all (atm% lb(1:2) == (/0,1/))) then
      !------------
      ! north poles
      !------------
      do l=1,5
        a = 2._wp * pi * (l-1._wp) / 5
        c = cos(a)
        s = sin(a)
        u = atm% u (0,1,:,l)
        v = atm% v (0,1,:,l)
        atm% u(0,1,:,l) = c * u + s * v
        atm% v(0,1,:,l) = c * v - s * u
      end do
      !------------
      ! south poles
      !------------
      do l=6,10
        a = 2._wp * pi * (l-5.5_wp) / 5
        c = cos(a)
        s = sin(a)
        u = atm% u (0,1,:,l)
        v = atm% v (0,1,:,l)
        atm% u(0,1,:,l) = c * u - s * v
        atm% v(0,1,:,l) = c * v + s * u
      end do
    end if
  end subroutine rotate_pole_local
!------------------------------------------------------------------------------
  subroutine rotate_wind_global (atm, i, j)
  type (t_atm) ,intent(inout) :: atm
  integer      ,intent(in)    :: i, j
  !---------------------------------------------------------------------------
  ! Rotate the winds near the poles to the 'global' coordinate system:
  !   u is orientated to the east if the pole is approached on the 0-meridian.
  !   For the north pole this corresponds with the local coordinate system
  !   of diamond 1.
  ! The 'local'  coordinate system corresponds with the GME convention.
  ! The 'global' coordinate system corresponds with the 3D-Var convention.
  ! does not work for points with index (0,j)
  !---------------------------------------------------------------------------
    real(wp) :: s, c           ! sin, cos of the rotation angle
    real(wp) :: a              ! rotation angle
    real(wp) :: u (atm% ub(3)) ! old wind component u
    real(wp) :: v (atm% ub(3)) ! old wind component v
    integer  :: l              ! diamond index

    if (i==0) call finish ('rotate_wind_global','i==0')
    do l=1,10
      if  (atm%grid%marr(1,i,j,l) /= dace% pe) cycle
      a  = atm%grid%rlon(i,j,1,l)
      c = cos(a)
      s = sin(a)
      u = atm% u (i,j,:,l)
      v = atm% v (i,j,:,l)
      if (l<=5) then
        atm% u(i,j,:,l) = c * u - s * v
        atm% v(i,j,:,l) = c * v + s * u
      else
        atm% u(i,j,:,l) = c * u + s * v
        atm% v(i,j,:,l) = c * v - s * u
      endif
    end do
  end subroutine rotate_wind_global
!------------------------------------------------------------------------------
  subroutine rotate_wind_local (atm, i, j)
  type (t_atm) ,intent(inout) :: atm
  integer      ,intent(in)    :: i, j
  !---------------------------------------------------------------------------
  ! Rotate the winds near the poles to the 'local' coordinate system:
  !   u is orientated to the east if the pole is approached on the 0-meridian.
  !   For the north pole this corresponds with the local coordinate system
  !   of diamond 1.
  ! The 'local'  coordinate system corresponds with the GME convention.
  ! The 'global' coordinate system corresponds with the 3D-Var convention.
  ! does not work for points with index (0,j)
  !---------------------------------------------------------------------------
    real(wp) :: s, c           ! sin, cos of the rotation angle
    real(wp) :: a              ! rotation angle
    real(wp) :: u (atm% ub(3)) ! old wind component u
    real(wp) :: v (atm% ub(3)) ! old wind component v
    integer  :: l              ! diamond index

    if (i==0) call finish ('rotate_wind_local','i==0')
    do l=1,10
      if  (atm%grid%marr(1,i,j,l) /= dace% pe) cycle
      a  = atm%grid%rlon(i,j,1,l)
      c = cos(a)
      s = sin(a)
      u = atm% u (i,j,:,l)
      v = atm% v (i,j,:,l)
      if (l<=5) then
        atm% u(i,j,:,l) = c * u + s * v
        atm% v(i,j,:,l) = c * v - s * u
      else
        atm% u(i,j,:,l) = c * u - s * v
        atm% v(i,j,:,l) = c * v + s * u
      endif
    end do
  end subroutine rotate_wind_local
!------------------------------------------------------------------------------
  elemental subroutine scatter_wind_poles (atm)
  type (t_atm) ,intent(inout) :: atm
  !-----------------------------------------------------------------------
  ! Scatter (duplicate) the wind information of the poles to all diamonds.
  ! This makes sense only in the 'global' coordinate system.
  !-----------------------------------------------------------------------
    integer :: l, l1
    !---------------------------------
    ! check if this pe holds the poles
    !---------------------------------
    if (all (atm% lb(1:2) == (/0,1/))) then
      !------------
      ! north poles
      !------------
      l1 = atm% grid% marr (4,0,1,1)
      do l = 1, 5
        if (l == l1) cycle
        atm% u(0,1,:,l) = atm% u(0,1,:,l1)
        atm% v(0,1,:,l) = atm% v(0,1,:,l1)
      end do
      !------------
      ! south poles
      !------------
      l1 = atm% grid% marr (4,0,1,6)
      do l = 6, 10
        if (l == l1) cycle
        atm% u(0,1,:,l) = atm% u(0,1,:,l1)
        atm% v(0,1,:,l) = atm% v(0,1,:,l1)
      end do
    endif
  end subroutine scatter_wind_poles
!==============================================================================
  subroutine A_to_C_grid1 (atm, restore)
  type(t_atm) ,intent(inout)           :: atm (:)
  logical     ,intent(in)    ,optional :: restore
  !-------------------------------------------------------
  ! transform wind components ( u, v, u_10m, v_10m)
  ! from A to C grid (vector version)
  !-------------------------------------------------------
    integer :: i        ! loop index (states)
    do i = 1, size(atm)
      call A_to_C_grid (atm(i), restore)
    end do
  end subroutine A_to_C_grid1
!------------------------------------------------------------------------------
  subroutine A_to_C_grid (atm, restore)
  type(t_atm) ,intent(inout)           :: atm
  logical     ,intent(in)    ,optional :: restore
  !-------------------------------------------------------
  ! transform wind components ( u, v, u_10m, v_10m)
  ! from A to C grid.
  ! if 'restore' is true, the fields are restored from the
  ! saved values.
  !-------------------------------------------------------
    logical :: r
    r = .false.; if (present(restore)) r = restore
    if (atm% gp_wind == 'C') call finish('A_to_C_grid','state is on C grid')
    atm% gp_wind  = 'C'
    !---
    ! u
    !---
    if (associated (atm% u)) then
      if (r .and. associated (atm% u_c)) then
        atm% u = atm% u_c
      else
        call half_gp_w (atm% u, atm% grid)
      endif
    endif
    !---
    ! v
    !---
    if (associated (atm% v)) then
      if (r .and. associated (atm% v_c)) then
        atm% v = atm% v_c
      else
        call half_gp_s (atm% v, atm% grid)
      endif
    endif
    !------
    ! u_10m
    !------
    if (associated(atm% u_10m)) then
      if (r .and. associated (atm% u_10m_c)) then
        atm% u_10m = atm% u_10m_c
      else
        call half_gp_w (atm% u_10m, atm% grid)
      endif
    endif
    !------
    ! v_10m
    !------
    if (associated(atm% v_10m)) then
      if (r .and. associated (atm% v_10m_c)) then
        atm% v_10m = atm% v_10m_c
      else
        call half_gp_s (atm% v_10m, atm% grid)
      endif
    endif
  end subroutine A_to_C_grid
!-----------------------------------------------------------------------------
  subroutine C_to_A_grid1 (atm, restore)
  type(t_atm) ,intent(inout)           :: atm (:)
  logical     ,intent(in)    ,optional :: restore
  !-------------------------------------------------------
  ! transform wind components ( u, v, u_10m, v_10m)
  ! from C to A grid (vector version)
  !-------------------------------------------------------
    integer :: i        ! loop index (states)
    do i = 1, size(atm)
      call C_to_A_grid (atm(i), restore)
    end do
  end subroutine C_to_A_grid1
!------------------------------------------------------------------------------
  subroutine C_to_A_grid (atm, save)
  type(t_atm) ,intent(inout)        :: atm
  logical     ,intent(in) ,optional :: save
  !-------------------------------------------------------
  ! transform wind components ( u, v, u_10m, v_10m)
  ! from C to A grid.
  ! if 'save' is true, the fields are saved in the
  ! respective components ( u_c, v_c, u_10m_c, v_10m_c).
  !-------------------------------------------------------
    logical :: s
    s = .false.; if (present(save)) s = save
    if (atm% gp_wind == 'A') call finish('C_to_A_grid','state is on A grid')
    atm% gp_wind  = 'A'
    !---
    ! u
    !---
    if (associated (atm% u)) then
      if (s) then
        call allocate (atm, 'u_c')
        atm% u_c = atm% u
      endif
      call half_gp_e (atm% u, atm% grid)
    else
      call deallocate (atm, 'u_c')
    endif
    !---
    ! v
    !---
    if (associated (atm% v)) then
      if (s) then
        call allocate (atm, 'v_c')
        atm% v_c = atm% v
      endif
      call half_gp_n (atm% v, atm% grid)
    else
      call deallocate (atm, 'v_c')
    endif
    !------
    ! u_10m
    !------
    if (associated (atm% u_10m)) then
      if (s) then
        call allocate (atm, 'u_10m_c')
        atm% u_10m_c = atm% u_10m
      endif
      call half_gp_e (atm% u_10m, atm% grid)
    else
      call deallocate (atm, 'u_10m_c')
    endif
    !---
    ! v_10m
    !---
    if (associated (atm% v_10m)) then
      if (s) then
        call allocate (atm, 'v_10m_c')
        atm% v_10m_c = atm% v_10m
      endif
      call half_gp_n (atm% v_10m, atm% grid)
    else
      call deallocate (atm, 'v_10m_c')
    endif
  end subroutine C_to_A_grid
!------------------------------------------------------------------------------
  subroutine half_gp_e (x, g)
  real(wp)     ,intent(inout) :: x (:,:,:,:)  ! field
  type(t_grid) ,intent(in)    :: g            ! grid information
  !---------------------------------------------------------------
  ! interpolate a field to a location half a gridpoint to the East
  ! for regular field on a local domain
  !---------------------------------------------------------------
    real(wp) :: temp    (g% shape(1), g% shape(2))
    real(wp) :: slice_r (g% shape(2), size (x ,3))
    real(wp) :: slice_s (g% shape(2), size (x ,3))
    integer  :: i, j, k

    if (g% dc% nbaw /= -1) then
      slice_s = x(g% shape(1),:,:,1)
      call p_isend (slice_s, g% dc% nbaw, 1)
    endif
    if (g% dc% nbpe /= -1) then
      call p_recv  (slice_r, g% dc% nbpe, 1)
    else
      slice_r = x(1,:,:,1)
    endif
    call p_wait

    do k = 1, size (x ,3)
      do j = 1, g% shape(2)
        do i = 2, g% shape(1)
          temp(i,j) = (x(i-1,j,k,1) + x(i,j,k,1)) * 0.5_wp
        enddo
        temp (1,j) = (slice_r(j,k) + x (1,j,k,1)) * 0.5_wp
      enddo
      x (:,:,k,1) = temp (:,:)
    enddo

  end subroutine half_gp_e
!------------------------------------------------------------------------------
  subroutine half_gp_w (x, g)
  real(wp)     ,intent(inout) :: x (:,:,:,:)  ! field
  type(t_grid) ,intent(in)    :: g            ! grid information
  !---------------------------------------------------------------
  ! interpolate a field to a location half a gridpoint to the West
  ! for regular field on a local domain
  !---------------------------------------------------------------
    real(wp) :: temp    (g% shape(1), g% shape(2))
    real(wp) :: slice_r (g% shape(2), size (x ,3))
    real(wp) :: slice_s (g% shape(2), size (x ,3))
    integer  :: i, j, k

    if (g% dc% nbpe /= -1) then
      slice_s = x(1,:,:,1)
      call p_isend (slice_s, g% dc% nbpe, 1)
    endif
    if (g% dc% nbaw /= -1) then
      call p_recv  (slice_r, g% dc% nbaw, 1)
    else
      slice_r = x(g% shape(1),:,:,1)
    endif
    call p_wait

    do k = 1, size (x ,3)
      do j = 1, g% shape(2)
        do i = 1, g% shape(1)-1
          temp(i,j) = (x(i,j,k,1) + x(i+1,j,k,1)) * 0.5_wp
        enddo
        temp(g% shape(1),j) = (slice_r(j,k) + x(g% shape(1),j,k,1)) * 0.5_wp
      enddo
      x(:,:,k,1) = temp(:,:)
    enddo

  end subroutine half_gp_w
!------------------------------------------------------------------------------
  subroutine half_gp_n (x, g)
  real(wp)     ,intent(inout) :: x (:,:,:,:)  ! field
  type(t_grid) ,intent(in)    :: g            ! grid information
  !----------------------------------------------------------------
  ! interpolate a field to a location half a gridpoint to the North
  ! for regular field on a local domain
  !----------------------------------------------------------------
    real(wp) :: temp    (g% shape(1), g% shape(2))
    real(wp) :: slice_r (g% shape(1), size (x ,3))
    real(wp) :: slice_s (g% shape(1), size (x ,3))
    integer  :: i, j, k

     if (g% dc% nbae /= -1) then
       slice_s = x(:,g% shape(2),:,1)
       call p_isend (slice_s, g% dc% nbae, 1)
     endif
     if (g% dc% nbpw /= -1) then
       call p_recv  (slice_r, g% dc% nbpw, 1)
     else
      slice_r = x(:,1,:,1)
     endif
     call p_wait

    do k = 1, size (x ,3)
      do i = 1, g% shape(1)
        do j = 2, g% shape(2)
          temp(i,j) = (x(i,j-1,k,1) + x(i,j,k,1)) * 0.5_wp
        enddo
        temp(i,1) = (slice_r(i,k) + x(i,1,k,1)) * 0.5_wp
      enddo
      x(:,:,k,1) = temp(:,:)
    enddo

  end subroutine half_gp_n
!------------------------------------------------------------------------------
  subroutine half_gp_s (x, g)
  real(wp)     ,intent(inout) :: x (:,:,:,:)  ! field
  type(t_grid) ,intent(in)    :: g            ! grid information
  !----------------------------------------------------------------
  ! interpolate a field to a location half a gridpoint to the South
  ! for regular field on a local domain
  !----------------------------------------------------------------
    real(wp) :: temp    (g% shape(1), g% shape(2))
    real(wp) :: slice_r (g% shape(1), size (x ,3))
    real(wp) :: slice_s (g% shape(1), size (x ,3))
    integer  :: i, j, k

     if (g% dc% nbpw /= -1) then
       slice_s = x(:,1,:,1)
       call p_isend (slice_s, g% dc% nbpw, 1)
     endif
     if (g% dc% nbae /= -1) then
       call p_recv  (slice_r, g% dc% nbae, 1)
     else
      slice_r = x(:,g% shape(2),:,1)
     endif
     call p_wait

      do k = 1, size (x ,3)
        do i = 1, g% shape(1)
          do j = 1, g% shape(2)-1
            temp (i,j) = (x (i,j,k,1) + x (i,j+1,k,1)) * 0.5_wp
          enddo
          temp (i,g% shape(2)) = (slice_r(i,k) + x(i,g% shape(2),k,1)) * 0.5_wp
        enddo
        x (:,:,k,1) = temp(:,:)
      enddo

  end subroutine half_gp_s
!==============================================================================
  subroutine print_sel_par (dest, src1, par1, src2, par2, &
                                  src3, par3, src4, par4, dealloc, comment)
  !-----------------------------------------------------------
  ! to be called BEFORE 'select_params'
  ! in order to produce a printout for the selected parameters
  !-----------------------------------------------------------
  type (t_atm)     ,intent(in)           :: dest    ! state to prepare
  character(len=*) ,intent(in)           :: src1    ! source 1
  character(len=*) ,intent(in)           :: par1    ! list of params from src1
  character(len=*) ,intent(in) ,optional :: src2    ! source 2
  character(len=*) ,intent(in) ,optional :: par2    ! list of params from src2
  character(len=*) ,intent(in) ,optional :: src3    ! source 3
  character(len=*) ,intent(in) ,optional :: par3    ! list of params from src3
  character(len=*) ,intent(in) ,optional :: src4    ! source 4
  character(len=*) ,intent(in) ,optional :: par4    ! list of params from src4
  logical          ,intent(in) ,optional :: dealloc ! deallocate parameters
  character(len=*) ,intent(in) ,optional :: comment ! comment for printout

    integer ,parameter :: np = 4
    character(len=16)  :: pars (nm,np)
    integer            :: src  (nm)
    character(len=16)  :: nam
    character(len=16)  :: srcx (np)
    integer            :: n       (np)
    integer            :: i, m
    logical            :: de
    logical            :: p1, p2, p3, p4

    p1 = .true.
    p2 = present (src2)
    p3 = present (src3)
    p4 = present (src4)

    if (dace% lpio) then
      n   = 0
      src = 0
      de  = .true. ;if (present(dealloc)) de = dealloc

      write (6,'(a)') repeat('-',79)
      write (6,'()')
      if (present(comment)) then
        write (6,'(a,a,a,a)') '  Selecting parameters for ',&
          trim(dest% name)//' ('//trim(comment)//')'
      else
        write (6,'(a,a)') '  Selecting parameters for ',trim(dest% name)
      endif
      if (p1) write (6,'(a,a16,a,a)') '    from ',trim(src1),' : ',trim(par1)
      if (p2) write (6,'(a,a16,a,a)') '         ',trim(src2),' : ',trim(par2)
      if (p3) write (6,'(a,a16,a,a)') '         ',trim(src3),' : ',trim(par3)
      if (p4) write (6,'(a,a16,a,a)') '         ',trim(src4),' : ',trim(par4)
      write (6,'()')
      if (p1) call split (pars(:,1), par1, n(1))
      if (p2) call split (pars(:,2), par2, n(2))
      if (p3) call split (pars(:,3), par3, n(3))
      if (p4) call split (pars(:,4), par4, n(4))
      if (p1) srcx(1) = src1
      if (p2) srcx(2) = src2
      if (p3) srcx(3) = src3
      if (p4) srcx(4) = src4

      do i = 1, np
        if (n(i)==0) cycle
        do m = 1, size (dest% m)
          nam = dest% m(m)% i% name
          if (nam == 'u_c') nam = 'u'
          if (nam == 'v_c') nam = 'v'
          if(any(pars(:n(i),i)==nam)) src (m) = i
        end do
      end do

      do m = 1, size (dest% m)
        if (src (m) /= 0) then
          write (6,'(4x,a,a,a,"(",i1,")")') &
            dest% m(m)% i% name, 'taken from: ',srcx(src(m)), src(m)
        else if (dest% m(m)% i% alloc) then
          if (de) then
            write (6,'(4x,a,a)') dest% m(m)% i% name, 'DEALLOCATED'
          else
            write (6,'(4x,a,a)') dest% m(m)% i% name, 'KEPT'
          endif
        endif
      end do
      write (6,'()')

    endif

  end subroutine print_sel_par
!------------------------------------------------------------------------------
  elemental subroutine select_params (dest, src1, par1, src2, par2, &
                                            src3, par3, src4, par4, &
                                            smnd, padd, dealloc     )
  !-------------------------------------------------------------------------
  ! prepare atmospheric state (dest) by taking parameters (variables)
  ! from the sources (src1, src..)
  ! as specified by the parameter lists (par1, par..) .
  ! if 'dealloc' is .true. (the default)  all fields in 'dest' which
  ! are not present in any parameter list are deallocated .
  ! for soil variables MAGIC values (-huge(wp)) are not copied.
  !------------------------------------------------------------------------
  type (t_atm)     ,intent(inout)        :: dest    ! state to prepare
  type (t_atm)     ,intent(in)           :: src1    ! source 1
  character(len=*) ,intent(in)           :: par1    ! list of params from src1
  type (t_atm)     ,intent(in) ,optional :: src2    ! source 2
  character(len=*) ,intent(in) ,optional :: par2    ! list of params from src2
  type (t_atm)     ,intent(in) ,optional :: src3    ! source 3
  character(len=*) ,intent(in) ,optional :: par3    ! list of params from src3
  type (t_atm)     ,intent(in) ,optional :: src4    ! source 4
  character(len=*) ,intent(in) ,optional :: par4    ! list of params from src4
  type (t_atm)     ,intent(in) ,optional :: smnd    ! summand
  character(len=*) ,intent(in) ,optional :: padd    ! list of params to add
  logical          ,intent(in) ,optional :: dealloc ! deallocate parameters

    logical :: mask  (nm)
    logical :: de

    mask = .false.
    de   = .true. ;if (present(dealloc)) de = dealloc
    call                     select_params0 (dest, src1, par1, 0, mask)
    if (present (src2)) call select_params0 (dest, src2, par2, 0, mask)
    if (present (src3)) call select_params0 (dest, src3, par3, 0, mask)
    if (present (src4)) call select_params0 (dest, src4, par4, 0, mask)
    if (present (smnd)) call select_params0 (dest, smnd, padd, 1, mask)
    if (de)             call   keep_params  (dest,                mask)
    call set_pointers (dest)

  end subroutine select_params
!------------------------------------------------------------------------------
  pure subroutine select_params0 (dest, src, par, add, mask)
  type (t_atm)     ,intent(inout) :: dest    ! atmospheric state to prepare
  type (t_atm)     ,intent(in)    :: src     ! source
  character(len=*) ,intent(in)    :: par     ! list of parameters
  integer          ,intent(in)    :: add     ! 0:replace 1:add
  logical          ,intent(inout) :: mask(:) ! mark processed entries
  !-----------------------------------------------------------------------
  ! prepare atmospheric state by taking parameters (variables) from src
  ! as specified by the parameter list par .
  ! for soil variables MAGIC values (-huge(wp)) are not copied.
  !-----------------------------------------------------------------------
    integer           :: m     ! loop index (parameters)
    character(len=16) :: name  ! name of parameter processed
    character(len=16) :: pars(nm)
    integer           :: n

    if (par == '') return
    call split (pars, par, n)
    do m = 1, size (dest% m)
      name = dest% m(m)% i% name
      if (name == 'u_c') name = 'u'
      if (name == 'v_c') name = 'v'
      if(any(pars(:n)==name)) then
        if (src% m(m)% i% alloc) then
          mask (m) = .true.
          if (.not. dest% m(m)% i% alloc) then
            call allocate (dest, name)
            dest% m(m)% ptr = -huge(1._wp)   ! indication for not set so far
          endif
!!! TODO: unsure how to convert this to TR15581:
#ifdef TR15581
          if (.true.) then
!           call finish ("TR15581","not implemented")
#else
          if (.not. associated (dest% m(m)% ptr, src% m(m)% ptr)) then
#endif
            select case (src% m(m)% i% name)
            case ('t_so'  ,'w_so'  ,'tsurf' ,'t_ice','h_ice'   ,'fr_ice', &
                  'h_snow','t_snow','w_snow','w_i'  ,'freshsnw','rho_snow',&
                  'snowc')
              !-------------------------------------------
              ! surface fields: only copy non-MAGIC values
              !-------------------------------------------
              if (add == 0) then
!++++++++++++++++++++++++++++++++++++
! work around bug in NEC sxf90/rev441
!++++++++++++++++++++++++++++++++++++
!               where (src% m(m)% ptr /= -huge(1._wp)) &
!                     dest% m(m)% ptr = src% m(m)% ptr

                call conditional_assign (dest% m(m)% ptr, src% m(m)% ptr, &
                                    size( src% m(m)% ptr)                 )
              else if (add == 1) then
                where (src% m(m)% ptr /= -huge(1._wp)) &
                  dest% m(m)% ptr = dest% m(m)% ptr + src% m(m)% ptr
              endif
            case default
              !-----------------------
              ! other fields: copy all
              !-----------------------
              if (add == 0) then
                dest% m(m)% ptr = src% m(m)% ptr
              else if (add == 1) then
                dest% m(m)% ptr = dest% m(m)% ptr + src% m(m)% ptr
              endif
            end select
          endif
        endif
      endif
    end do
  contains
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!++++++++++++++++++++++++++++++++++++
! work around bug in NEC sxf90/rev441
!++++++++++++++++++++++++++++++++++++
  pure subroutine conditional_assign (d,s,n)
  integer  ,intent(in)    ::   n
  real(wp) ,intent(inout) :: d(n)
  real(wp) ,intent(in)    :: s(n)
    where (s/= -huge(1._wp)) d = s
  end subroutine conditional_assign
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine select_params0
!------------------------------------------------------------------------------
  subroutine derive_params1 (dest, par, dealloc, comment)
  !------------------------------------------------------
  ! transform fields (eg. q <-> rh etc.) (vector version)
  !------------------------------------------------------
  type (t_atm)     ,intent(inout)        :: dest(:) ! state to prepare
  character(len=*) ,intent(in)           :: par     ! list of parameters
  logical          ,intent(in)           :: dealloc ! deallocate unused params
  character(len=*) ,intent(in) ,optional :: comment ! for printout

    integer           :: i        ! loop index (states)
    integer           :: verbose

    do i = 1, size(dest)
      verbose           = 0
      if (i==1) verbose = 1
      call derive_params0 (dest(i), par, dealloc, verbose, comment)
    end do

  end subroutine derive_params1
!------------------------------------------------------------------------------
  subroutine derive_params0 (dest, par, dealloc, verbose, comment)
  !------------------------------------------------------
  ! transform fields (eg. q <-> rh etc.) (scalar version)
  !------------------------------------------------------
  type (t_atm)     ,intent(inout)        :: dest    ! state to prepare
  character(len=*) ,intent(in)           :: par     ! list of parameters
  logical          ,intent(in)           :: dealloc ! deallocate unused params
  integer          ,intent(in) ,optional :: verbose ! verbosity flag
  character(len=*) ,intent(in) ,optional :: comment ! for printout

    integer           :: m        ! loop index (parameters)
    character(len=16) :: name     ! name of parameter processed
    logical           :: mask(nm) ! variables to keep flag
    character(len=16) :: pars(nm) ! names of variables to derive
    integer           :: n        ! return code
    integer           :: verb     ! verbosity flag
    character(len=32) :: from     ! text to print
    integer           :: error    ! error occured for this variable

    verb = 0; if (present(verbose)) verb = verbose

    call split (pars, par, n)
    if(verb>0 .and. dace% lpio) then
      write (6,'(a)') repeat('-',79)
      write (6,'()')
      if (present(comment)) then
        write (6,'(a,a,a)') '  Deriving parameters for ',&
          trim(dest% name)//' ('//trim(comment)//')'
      else
        write (6,'(a,a)') '  Deriving parameters for ',trim(dest% name)
      endif
      write (6,'(a,l1)') '    dealloc = ',dealloc
      write (6,'()')
    endif
    mask  = .false.
    error = 0
    do m = 1, size (dest% m)
      name = dest% m(m)% i% name
      if(any(pars==name)) then
        mask (m) = dest% m(m)% i% alloc
        call allocate (dest, name)
        if (.not. mask (m)) then
          mask(m) = .true.
          from    = ''
          select case (name)
          case ('psr')
            dest% psr = dest% ps
          case ('pp')
            if (.not.(associated(dest% grid% p0).and. &
                      associated(dest%       pf)      )) &
              call finish('derive_params','pp: p0 or pf not present')
            dest% pp = dest% pf - dest% grid% p0
            from = 'derived from (pf-p0)'
          case ('ps')
            call set_ps (dest)
            from = 'derived from set_ps'
          case ('rh')
            if (.not.associated(dest% pf)) call set_p(dest)
!NEC$ inline
            dest% rh = rh_q (dest% q, dest% t, dest% pf)
            from = 'derived from rh_q(q)'
          case ('geoh')
            if (.not.associated(dest% ph)) call set_p(dest)
            call set_geo (dest)
            from = 'derived from set_geo (p)'
          case ('geof')
            if (.not.associated(dest% ph)) call set_p(dest)
            call set_geo (dest, geof=.true.)
            from = 'derived from set_geo (p)'
          case ('t')
            if (.not. associated (dest% ph)) call set_p      (dest)
            if (.not. associated (dest% tv)) call set_tv_geo (dest)
            call set_tq_tv (dest)
            from = 'derived from set_tq_tv(tv)'
          case ('tv')
            call set_tv (dest)
            from = 'derived from set_tv(t,q,p)'
          case ('den')
            call set_den (dest)
            from = 'derived from set_den(tv,p)'
          case ('q')
            if (.not.associated(dest% pf)) call set_p(dest)
            dest% q = q_rh (dest% rh, dest% t, dest% pf)
            from = 'derived from q_rh (rh, t, p)'
          case default
            error   =  m
            from = 'NOT IMPLEMENTED'
          end select
          if(verb>0 .and. dace% lpio) write (6,'(4x,a,a,a)') &
            dest% m(m)% i% name, from
        else
          if(verb>0 .and. dace% lpio) &
            write (6,'(4x,a,a)') dest% m(m)% i% name, 'KEPT (required)'
        endif
      else
        if(verb>0 .and. dace% lpio) then
          if (dest% m(m)% i% alloc) then
            if (dealloc) then
              write (6,'(4x,a,a)') dest% m(m)% i% name, 'DEALLOCATED'
            else
              write (6,'(4x,a,a)') dest% m(m)% i% name, 'KEPT'
            endif
          endif
        endif
      endif
    end do

    !---------------
    ! error handling
    !---------------
    if (error > 0)                                                &
      call finish ('derive_params',                               &
                   'cannot derive '//trim(dest% m(error)% i% name))

    !------------------------------------
    ! deallocate parameters not requested
    !------------------------------------
    if (dealloc) call keep_params (dest, mask)

  end subroutine derive_params0
!------------------------------------------------------------------------------
  pure subroutine keep_params (dest, mask)
  type (t_atm)     ,intent(inout) :: dest    ! atmospheric state to prepare
  logical          ,intent(in)    :: mask(:)
  !------------------------------------
  ! deallocate parameters not requested
  !------------------------------------
    integer           :: m     ! loop index (parameters)
    do m = 1, size (dest% m)
      if(.not.mask(m)) call deallocate (dest, dest%m(m)%i% name)
    end do
  end subroutine keep_params
!------------------------------------------------------------------------------
  pure subroutine keep_params_name (dest, par)
  type (t_atm)     ,intent(inout) :: dest    ! atmospheric state to prepare
  character(len=*) ,intent(in)    :: par
  !------------------------------------
  ! deallocate parameters not requested
  !------------------------------------
    integer           :: m        ! loop index (parameters)
    integer           :: n        ! return code
    character(len=16) :: pars(nm) ! names of variables to derive
    call split (pars, par, n)
    do m = 1, size (dest% m)
      if(all(pars(:n) /= dest% m(m)% i% name))  &
        call deallocate (dest, dest%m(m)%i% name)
    end do
  end subroutine keep_params_name
!=============================================================================
  subroutine hydro_balanc (fc, ana, bal_var)
  !----------------------------------------------------------
  ! hydrostatic balancing of pressure and temperature profile
  !----------------------------------------------------------
  type (t_atm) ,intent(in)    :: fc      ! reference forecast
  type (t_atm) ,intent(inout) :: ana     ! analysis to balance
  integer      ,intent(in)    :: bal_var ! 0: balance p; 1 balance t

    integer  :: i,j,n
    real(wp) :: zfbuyot       ! buyot term at level n
    real(wp) :: zfbuyob       ! buyot term at level n+1
    real(wp) :: qaninb        ! virt. temp. inc. at lev. n+1
    real(wp) :: qaninc        ! virt. temp. inc. at lev. n
!   real(wp) :: taninb        ! old temp. inc.
!   real(wp) :: taninc        ! new temp. inc
    integer  :: lb(4), ub(4)  ! Copy of grid local bounds
    logical  :: lqr, lqs      ! qr, qs present ?

    !------------------------------
    ! set enkf_an to ana increments
    !------------------------------
    !do k = 1, k_enkf
    !   enkf_an(k) = enkf_an(k) - enkf_fc(k)
    !end do
    !replaced by explicit loop:
    !later replace here (and in sat_adj, hydro_balanc)
    ! in i,j loop fc% grid by enkf_an/fc % grid

    if (.not. associated (ana% pf)) &
      call finish('hydro_balanc', 'ana% pf not allocated !')
    if (.not. associated ( fc% pf)) &
      call finish('hydro_balanc',  'fc% pf not allocated !')

    lqr = associated (fc% qr) .and. associated (ana% qr)
    lqs = associated (fc% qs) .and. associated (ana% qs)
    lb  = fc% grid% lb
    ub  = fc% grid% ub

!$omp parallel do private(i,j,n) schedule(static)
      do     n = 1, fc% grid% nz
        do   j = lb(2), ub(2)
!DIR$ IVDEP
          do i = lb(1), ub(1)
            ana% q  (i,j,n,1) = ana% q  (i,j,n,1) - fc% q  (i,j,n,1)
            ana% qcl(i,j,n,1) = ana% qcl(i,j,n,1) - fc% qcl(i,j,n,1)
            ana% pf (i,j,n,1) = ana% pf (i,j,n,1) - fc% pf (i,j,n,1)
            ana% t  (i,j,n,1) = ana% t  (i,j,n,1) - fc% t  (i,j,n,1)
          end do
        end do
      end do
!$omp end parallel do

    !------------------------------------------------
    ! hydrostatic balancing of pressure / temperature
    ! (up to now only for pressure implemented)
    !------------------------------------------------
      do     n = fc% grid% nz-1, 1, -1
!$omp parallel do private(i,j,zfbuyot,zfbuyob,qaninc,qaninb) schedule(static)
        do   j = lb(2), ub(2)
!DIR$ IVDEP
          do i = lb(1), ub(1)
            zfbuyot =  R * fc% grid% rho0(i,j,n  ,1) * fc% t(i,j,n  ,1)
            zfbuyob =  R * fc% grid% rho0(i,j,n+1,1) * fc% t(i,j,n+1,1)
            qaninc = rddrm1 * ana% q  (i,j,n  ,1) &
                            - ana% qcl(i,j,n  ,1)
            qaninb = rddrm1 * ana% q  (i,j,n+1,1) &
                            - ana% qcl(i,j,n+1,1)
            if (lqr) then
              qaninc = qaninc - (ana% qr (i,j,n  ,1) - fc% qr (i,j,n  ,1))
              qaninb = qaninb - (ana% qr (i,j,n+1,1) - fc% qr (i,j,n+1,1))
            endif
            if (lqs) then
              qaninc = qaninc - (ana% qs (i,j,n  ,1) - fc% qs (i,j,n  ,1))
              qaninb = qaninb - (ana% qs (i,j,n+1,1) - fc% qs (i,j,n+1,1))
            endif

            !select case(bal_var)
            !case(0)
            ana% pf(i,j,n,1) =                                   &
              1._wp / (2._wp + fc% grid% dp0(i,j,n+1,1)/zfbuyot) &
              *( ana% pf(i,j,n+1,1) *                            &
               (2._wp - fc% grid% dp0(i,j,n,1)/zfbuyob)          &
              +  ana% t(i,j,n,1)*fc% grid% dp0(i,j,n+1,1)        &
                 / fc% t(i,j,n,1)                                &
              +  ana% t(i,j,n+1,1)*fc% grid% dp0(i,j,n,1)        &
                 / fc% t(i,j,n+1,1)                              &
              +  qaninc*fc% grid% dp0(i,j,n+1,1)                 &
              +  qaninb*fc% grid% dp0(i,j,n,1)  )
            !case(1)
            !if (n==fc% grid% nz-1) taninc = ana% t(i,j,fc% grid% nz,1)
            !taninb = taninc
            !taninc = ana% t(i,j,n,1)
            !ana% t(i,j,n,1) =                                     &
            !     ( (2._wp + fc% grid% dp0(i,j,n+1,1)/zfbuyot) *   &
            !     ana% pf(i,j,n,1)                                 &
            !     - ana% pf(i,j,n+1,1) *                           &
            !     (2._wp - fc% grid% dp0(i,j,n,1)/zfbuyob)         &
            !     -  taninb*fc% grid% dp0(i,j,n,1)                 &
            !        / fc% t(i,j,n+1,1)                            &
            !     -  qaninc*fc% grid% dp0(i,j,n+1,1)               &
            !     -  qaninb*fc% grid% dp0(i,j,n,1) ) *             &
            !     fc% t(i,j,n,1)/fc% grid% dp0(i,j,n+1,1)
            !end select
          end do
        end do
!$omp end parallel do
      end do

    !--------------------------------------
    !add modified ana increments to enkf_fc
    !--------------------------------------
!   print *, 'add new incr to ana'
    !do k = 1, k_ens
    !   ana = ana + fc
    !end do
    !replace by explicit loop
!$omp parallel do private(i,j,n) schedule(static)
      do     n = 1, fc% grid% nz
        do   j = lb(2), ub(2)
!DIR$ IVDEP
          do i = lb(1), ub(1)
            ana% q  (i,j,n,1) = ana% q  (i,j,n,1) + fc% q  (i,j,n,1)
            ana% qcl(i,j,n,1) = ana% qcl(i,j,n,1) + fc% qcl(i,j,n,1)
            ana% pf (i,j,n,1) = ana% pf (i,j,n,1) + fc% pf (i,j,n,1)
            ana% t  (i,j,n,1) = ana% t  (i,j,n,1) + fc% t  (i,j,n,1)
          end do
        end do
      end do
!$omp end parallel do
  if (associated (ana% pp)) then
     ana % pp = ana % pf - ana % grid % p0
  end if

  end subroutine hydro_balanc
!==============================================================================
  subroutine keep_diamond (diamond, state, grid)
  integer     ,intent(in) :: diamond
  type(t_atm)  ,intent(inout) ,optional :: state
  type(t_grid) ,intent(inout) ,optional :: grid
  !---------------------------------------------------------------------------
  ! keep only the desired diamond
  ! for an ICON ensemble with gridpoints re-distributed to a GME partitioning.
  ! (for testing so that only a limeted amount of data is kept)
  !---------------------------------------------------------------------------

    integer               :: d_gme (0:10)
    integer               :: lb1, ub1
    integer               :: nd, m, i
    real(wp) ,ALLOCATABLE :: ptr (:,:,:,:)
    integer(i8)           :: siz

    if (diamond < 1) return

    if (present (state)) then

      !----------------------------------------------------
      ! calculate local bounds holding the required diamond
      !----------------------------------------------------
      d_gme = state% grid% d_gme
      if (d_gme(0) < 0) return
      nd  = d_gme (diamond)   - d_gme (diamond-1)
      lb1 = state% grid% lb(1) + d_gme (diamond-1)
      ub1 = lb1 + nd - 1

      !------------------------------------------------------------
      ! re-allocate fields of atmospheric state within local bounds
      !------------------------------------------------------------
      state% lb(1) = lb1
      state% ub(1) = ub1
      do m = 1, size(state% m)
        state% m(m)% i% lb(1) = lb1
        state% m(m)% i% ub(1) = ub1
        if (.not. state% m(m)% i% alloc) cycle
        allocate (ptr (state% m(m)% i% lb(1) : state% m(m)% i% ub(1), &
                       state% m(m)% i% lb(2) : state% m(m)% i% ub(2), &
                       state% m(m)% i% lb(3) : state% m(m)% i% ub(3), &
                       state% m(m)% i% lb(4) : state% m(m)% i% ub(4)) )
        ptr = state% m(m)% ptr (lb1:ub1,:,:,:)
        deallocate (state% m(m)% ptr)
        MOVE_ALLOC (ptr, state% m(m)% ptr)  ! state% m(m)% ptr => ptr
      end do
      call set_pointers (state)

      siz = state% size
      if (dace% lpio) write(6,*) '  keep_diamond:',diamond, siz, state% size

    endif

    if (present (grid)) then

      d_gme = grid% d_gme
      if (d_gme(0) < 0) return
      nd  = d_gme (diamond)   - d_gme (diamond-1)
      lb1 = grid% lb(1) + d_gme (diamond-1)
      ub1 = lb1 + nd - 1

      !-----------------------------------
      ! adjust grid% lb, ub (local bounds)
      !-----------------------------------
      grid% lb   (1)  = lb1
      grid% ub   (1)  = ub1
      grid% shape(1)  = nd
!     grid% dc% ilim1 = ip + 1
!     where (grid%m%i% lb(1) == lb1 .and. grid%m%i% ub(1) == ub1)
!       grid%m%i% lb(1) = grid% lb(1)
!       grid%m%i% ub(1) = grid% ub(1)
!     endwhere


      !------------------------------------------
      ! adjust grid% d_gme (local diamond offset)
      !------------------------------------------
      grid% d_gme (:diamond-1) = 0
      grid% d_gme ( diamond: ) = nd

      !----------------------------------------
      ! adjust grid% marr (permutation indices)
      !----------------------------------------
      do i=1, size (grid% marr, 1)
        if (grid% marr (4,i,1,1) /= diamond) then
          grid% marr ( 1 ,i,1,1) = -1
          grid% marr (2:4,i,1,1) =  0
        endif
      end do

!       pe = ix_pe(i)
!       d  = ix_d (i)
!       j  = id (d-1, pe  ) + 1
!       id      (d-1, pe  ) = j
!       j  = ip (     pe-1) + j
!       grid% marr (1,i,1,1) = pe
!       grid% marr (2,i,1,1) = j
!       grid% marr (3,i,1,1) = 1
!       grid% marr (4,i,1,1) = 1
!     end do
!
!     !--------------------------------------------------
!     ! adjust grid% rlon, rlat (coordinates already set)
!     !--------------------------------------------------
!     grid% rlat   (grid% marr (2,:,1,1),1,1,1) = grid% rlat   (:,1,1,1)
!     grid% rlon   (grid% marr (2,:,1,1),1,1,1) = grid% rlon   (:,1,1,1)
!     grid% xnglob (grid% marr (2,:,1,1),1,:,1) = grid% xnglob (:,1,:,1)

    endif

  end subroutine keep_diamond
!==============================================================================
end module mo_atm_state
