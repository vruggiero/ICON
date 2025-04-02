!
!+ Derived type to hold a model column (for inter-processor communication)
!
MODULE mo_t_col
!
! Description:
!   Definition of derived type to hold a model column. Used for
!   inter-processor communication.
!
! Current Maintainer: DWD, Harald Anlauf, Michael Bender
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
!  changes for COSMO vertical coordinate
! V1_5         2009/05/25 Harald Anlauf
!  gather_cols(extracted from mo_psas:get_cols):use send/recv or MPI_gather(v)
! V1_8         2009/12/09 Andreas Rhodin
!  Changes for COSMO: distribute pressure coordinates, dont infer from ak,bk
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_11        2010/06/16 Andreas Rhodin
!  COL_X: option to transfer water + ice load
! V1_13        2011/11/01 Andreas Rhodin
!  changed 3dvar revision numbers (CVS->SVN)
! V1_20        2012-06-18 Andreas Rhodin
!  allow interpolation of qcl,qci,qr,qs,qg
! V1_22        2013-02-13 Andreas Rhodin
!  handle t2m,rh2m,td2m,t_ice,h_ice,t_so
!  optimize MPI communication pattern to reduce latencies (H.Anlauf)
! V1_26        2013/06/27 Andreas Rhodin
!  option to interpolate inflation factor
! V1_28        2014/02/26 Andreas Rhodin
!  preparations for VarEnKF
! V1_31        2014-08-21 Andreas Rhodin
!  handle pp and w; check for presence of t,q,pf before calculating rhll
! V1_35        2014-11-07 Andreas Rhodin
!  atm2col: check for allocation of u,v before copying u10m,v10m
! V1_37        2014-12-23 Andreas Rhodin
!  changes for Variational Ensemble Kalman Filter (VarEnKF)
! V1_42        2015-06-08 Andreas Rhodin
!  implement temporal interpolation for COSMO MEC
! V1_45        2015-12-15 Harald Anlauf
!  use non-blocking communication; option to handle soil t_so, w_so
! V1_46        2016-02-05 Andreas Rhodin
!  base decisions on new flag 'vct', not 'ivctype'
! V1_47        2016-06-06 Andreas Rhodin
!  include cloud cover variables CLCT,CLCL,CLCM,CLCT,CLC; handle time range
! V1_48        2016-10-06 Harald Anlauf
!  Add roughness length z0 to t_spot, t_sl
! V1_51        2017-02-24 Valerio Cardinali
!  mo_t_col: add 'dpsdt' (for COMET)
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Andreas Rhodin  DWD  2002-2007
! Harald Anlauf   DWD  2007
! Michael Bender  DWD  2019
!==============================================================================
!
! Diagnostics for vectorization
!
#if defined (_FTRACE) && !defined (DISABLE_FTRACE_REGION)
#define FTRACE_BEGIN(text) CALL FTRACE_REGION_BEGIN (text)
#define FTRACE_END(text)   CALL FTRACE_REGION_END   (text)
#else
#define FTRACE_BEGIN(text)
#define FTRACE_END(text)
#endif
!==============================================================================
#include "tr15581.incf"
!==============================================================================
! Modules used
!=============
  use mo_kind,       only: wp,           &! working precision
                           i8             ! 64bit integer
  use mo_exception,  only: finish         ! abort in case of error
  use mo_physics,    only: pi,           &! pi
                           gacc,         &! gravity acceleration
                           p_ps,         &! derive p from ps
                           rhw_q          ! relative humidity over water
  use mo_t_obs,      only: t_mcol,       &! index array
                           t_mcols,      &! index array container
                           t_icol,       &!
                           t_coord,      &! coordinate data type
                           t_obs,        &! observation data type
                           t_spot,       &!
                           empty_spot,   &! empty spot
                           set_xuv        ! set up unit vectors on the sphere
  use mo_atm_state,  only: t_atm,              &! atmospheric state data type
                           VMAX_10M,TMIN_2M,   &! positions in t_atm
                           TMAX_2M,TOT_PREC,   &!
                           ASWDIR_S,ASWDIFD_S   !
  use mo_atm_grid,   only: t_grid,             &! model grid data type
                           same_horizontal_grid ! check for same grid
  use mo_physics,    only: tv_t_q         ! virtual temperature from t,q
  use mo_time,       only: t_time         ! time & date data type
  use mo_mpi_dace,   only: dace,         &! MPI group info
                           p_send,       &! generic send routine
                           p_recv,       &! generic receive routine
                           p_isend,      &! generic non-blocking send routine
                           p_irecv,      &! generic non-blocking receive routine
                           p_irequest,   &! index of next non-blocking request
                           p_mrequest,   &! maximum outstanding requests
                           p_waitall,    &! wait for MPI requests to complete
                           p_barrier      ! MPI barrier
  use mo_wmo_tables, only: WMO3_ISOBARIC  ! level type
  use mo_obs_sndrcv, only: p_mcol,       &! Pointers to model column descriptors
                           scatter_mcol, &! scatter model columns
                           alltoall_mcol  ! scatter model columns (alltoall)
  use mo_grid_intpol,only: idx_init,     &! interpolation indices and weights
                           t_h_idx        ! derived type for indices,weights

  use mo_fdbk_tables,only: VN_HEIGHT,    &! height               code
                           VN_HOSAG,     &! height above ground  code
                           VN_FLEV,      &! nominal flight level code
                           VN_P,         &! pressure
                           varno          ! 'varno' table
  use mo_t_table,    only: name_value     ! find name of table entry
  implicit none
!==============================================================================
! Public entities
!================
  private
  public :: t_sl              ! data type to hold single level data, part of t_col
  public :: t_col             ! data type to hold a column of the model
  public :: t_cols            ! data type to hold a bunch of columns
  public :: p_send            ! send columns to another processor
  public :: p_recv            ! receive columns from another processor
  public :: alloc_cols        ! allocate components of datatype t_cols
  public :: dealloc_cols      ! deallocate components of datatype t_cols
  public :: set_xuv           ! set up unit vectors on the sphere
  public :: assignment(=)     ! type(t_cols) = real(wp)
  public :: atm2col           ! gather columns from atmosphere
  public :: p_ps              ! derive p from ps
  public :: set_col_ptr
  public :: gather_cols       ! gather columns from other processors
  public :: igather_cols      ! gather_cols (using non-blocking comm.)
  public :: gather_cols_multi ! gather columns (all boxes, non-blocking comm.)
  public :: get_cols          ! redistribute model columns for obs processing
  public :: get_cols_adj      ! adjoint routine
  public :: set_cols_logp     ! store ln p instead of p
  public :: set_mcols_m       ! set t_mcols for grid to grid interpolation
  public :: mcol_idx          ! set t_mcols from t_h_idx derived type variable
  public :: nt                ! number of time ranges,
  public :: tr                ! time ranges in minutes
  public :: isend_cols        ! (non-blocking) send columns to another PE
  public :: irecv_cols        ! (non-blocking) receive columns from another PE
  public :: alltoall_cols     ! redistribute columns from all PEs
  public :: get_varlist       ! get a list of allocated identifiers
  public :: set_lev_insitu    ! set observation pressure levels from cols
  public :: u10_use_mlevel    ! use lowest model level as proxy for 10m wind
  !-------------------------
  ! identifier for variables
  !-------------------------
  public :: COL_T             ! temperature
  public :: COL_TV            ! virtual temperature
  public :: COL_Q             ! specific humidity
  public :: COL_RH            ! relative humidity
! public :: COL_U             ! u-wind component (deprecated)
! public :: COL_V             ! v-wind component (deprecated)
  public :: COL_UV            ! hor.wind u+v components
  public :: COL_P             ! pressure
  public :: COL_PH            ! pressure     (half levels)
  public :: COL_GEO           ! geopotential
  public :: COL_GEOH          ! geopotential (half levels)
  public :: COL_X             ! water load
  public :: COL_QCL           ! cloud water content
  public :: COL_QCI           ! cloud ice   content
  public :: COL_QR            ! rain
  public :: COL_QS            ! snow
  public :: COL_QG            ! graupel
  public :: COL_INFLAT        ! inflation factor
  public :: COL_TV2           ! tv from tv, not t,q
  public :: COL_PP            ! pressure deviation
  public :: COL_W             ! vertical velocity
  public :: COL_SOIL          ! soil temperature, humidity
  public :: COL_CLC           ! cloud cover
  public :: COL_RANGE         ! variables with time range indicator
  public :: COL_QVDIA         ! specific humidity, diagnostic
  public :: COL_QCDIA         ! spec. cloud water, diagnostic
  public :: COL_QIDIA         ! spec. cloud ice, diagnostic
  public :: COL_OZONE         ! ozone
  public :: COL_WSO           ! soil moisture only
  public :: COL_CO2           ! CO2
  public :: COL_REFF_QC       ! Effective radius of cloud water
  public :: COL_REFF_QI       ! Effective radius of cloud ice
  !-----------
  ! interfaces
  !-----------
  interface p_send
    module procedure send_cols
  end interface p_send

  interface p_recv
    module procedure recv_cols
  end interface p_recv

  interface dealloc_cols
    module procedure dealloc_col
    module procedure dealloc_cols
  end interface

  interface get_cols
    module procedure get_cols1
    module procedure get_cols
    module procedure get_cols_gg
  end interface get_cols

  interface get_cols_adj
    module procedure get_cols_adj
    module procedure get_cols_gg_adj
  end interface get_cols_adj

  interface set_cols_logp
    module procedure set_cols_logp
    module procedure set_cols_logp1
   end interface set_cols_logp

  interface set_xuv
    module procedure set_xuv_col
  end interface
!==============================================================================
! Data type definitions
!======================
  integer, parameter  :: nt     = 5                     ! number of time ranges
  integer, parameter  :: tr(nt) = [60,180,360,720,1440] ! ranges (minutes)
  !------------------
  ! single level data
  !------------------
  type t_sl
    !------------------
    ! single level data
    !------------------
    real(wp) :: ps     ! surface pressure
    real(wp) :: psr    ! reference surface pressure
    real(wp) :: dpsdt  ! surface pressure tendency
    real(wp) :: geosp  ! surface geopotential (orography)        [m2/s2]
    real(wp) :: geoid  ! geoid heights                           [m]
    real(wp) :: ts     ! surface temperature
!   real(wp) :: qs     ! surface humidity (currently unused)
    real(wp) :: lsm    ! land sea mask
    real(wp) :: ssd    ! SSO standard deviation                  [m]
    real(wp) :: z0     ! roughness length
    real(wp) :: tll    ! temperature at lowest model level       [K]
    real(wp) :: rhll   ! relative humidity at lowest model level [1]
    real(wp) :: ull    ! wind (u-comp.) at lowest model level    [m/s]
    real(wp) :: vll    ! wind (v-comp.) at lowest model level    [m/s]
    !--------------------------
    ! derived single level data
    !--------------------------
    real(wp) :: u10m
    real(wp) :: v10m
    real(wp) :: t2m
    real(wp) :: rh2m     ! 2m relative humidity                  [%]
    real(wp) :: td2m
    real(wp) :: dtdz     ! lapse rate (dT/dz) [K/m] for extrapolation
    real(wp) :: t2mland  ! t2m land tile average
    real(wp) :: rh2mland ! rh2m  "
    real(wp) :: td2mland ! td2m  "
    !--------------------------
    ! required for sst analysis
    !--------------------------
    real(wp) :: tsurf  ! temperature at surface                  [K]
    real(wp) :: fr_ice ! sea ice fraction
    real(wp) :: t_ice  ! ice temperature
    real(wp) :: h_ice  ! ice depth
    real(wp) :: t_so_0 ! surface temperature                     [K]
    !---------------------------
    ! required for snow analysis
    !---------------------------
    real(wp) :: h_snow   ! snow height
    real(wp) :: t_snow   ! snow temperature
    real(wp) :: w_snow   ! snow water equivalent
    real(wp) :: rho_snow ! snow density
    real(wp) :: freshsnw ! snow age
    real(wp) :: snowc    ! snow cover
    real(wp) :: w_i      ! interception reservoir water
    !------------------
    ! cloud diagnostics
    !------------------
    real(wp) :: clct     ! total  cloud cover
    real(wp) :: clcl     ! low    cloud cover
    real(wp) :: clcm     ! medium cloud cover
    real(wp) :: clch     ! high   cloud cover
    real(wp) :: ceiling  ! ceiling height
    real(wp) :: vis      ! hor. visibility
    !---------------------
    ! time range variables
    !---------------------
    real(wp) :: vmax_10m  (nt) ! max. 10m wind speed
    real(wp) :: tmin_2m   (nt) ! min.  2m temperature
    real(wp) :: tmax_2m   (nt) ! max.  2m temperature
    real(wp) :: tot_prec  (nt) ! total precipitation
    real(wp) :: aswdir_s  (nt) ! downw.direct SW rad
    real(wp) :: aswdifd_s (nt) ! down. diffusive r.
  end type t_sl
  !----------------------------------------------------------
  ! type t_col: metadata and single level data for one column
  !----------------------------------------------------------
  type t_col
    !-----------------------
    ! location and meta data
    !-----------------------
    type(t_coord)   :: c           ! coordinates:
    integer         :: i,j,l       ! horizontal/diamond index
    integer         :: pe          ! host processor element index
    integer         :: iom         ! index offset for multi level data
    integer         :: iot         ! index offset for derived data
    !------------------
    ! single level data
    !------------------
    type (t_sl)     :: s
    !-----------------
    ! multi level data
    !-----------------
    real(wp) _POINTER :: t    (:)   => NULL() ! temperature
    real(wp) _POINTER :: q    (:)   => NULL() ! specific humidity
    real(wp) _POINTER :: u    (:)   => NULL() ! wind component u
    real(wp) _POINTER :: v    (:)   => NULL() ! wind component v
!   real(wp) ,pointer :: tr   (:,:) => NULL() ! tracers
    !-------------------------
    ! derived multi level data
    !-------------------------
    real(wp) _POINTER :: p    (:)   => NULL() ! pressure
    real(wp) _POINTER :: ph   (:)   => NULL() ! pressure at half levels
    real(wp) _POINTER :: geo  (:)   => NULL() ! geopotential
    real(wp) _POINTER :: geoh (:)   => NULL() ! geopotential at half levels
    real(wp) _POINTER :: rh   (:)   => NULL() ! relative humidity
    real(wp) _POINTER :: tv   (:)   => NULL() ! virtual temperature
    real(wp) _POINTER :: x    (:)   => NULL() ! water load (cloud+ice+precip)
    real(wp) _POINTER :: qcl  (:)   => NULL() ! cloud water content
    real(wp) _POINTER :: qci  (:)   => NULL() ! cloud ice   content
    real(wp) _POINTER :: qr   (:)   => NULL() ! rain
    real(wp) _POINTER :: qs   (:)   => NULL() ! snow
    real(wp) _POINTER :: qg   (:)   => NULL() ! graupel
    real(wp) _POINTER :: infl (:)   => NULL() ! inflation factor
    real(wp) _POINTER :: irpc (:)   => NULL() ! IR PCA coefficients
    real(wp) _POINTER :: pp   (:)   => NULL() ! pressure deviation
    real(wp) _POINTER :: w    (:)   => NULL() ! vertical velocity
    real(wp) _POINTER :: clc  (:)   => NULL() ! cloud cover
    real(wp) _POINTER :: qv_dia(:)  => NULL() ! specific humidity, diag.
    real(wp) _POINTER :: qc_dia(:)  => NULL() ! spec. cloud water, diag.
    real(wp) _POINTER :: qi_dia(:)  => NULL() ! spec. cloud ice, diag.
    real(wp) _POINTER :: o3   (:)   => NULL() ! ozone
    real(wp) _POINTER :: co2  (:)   => NULL() ! carbon dioxide
    real(wp) _POINTER :: reff_qc(:) => NULL() ! effective radius cloud water
    real(wp) _POINTER :: reff_qi(:) => NULL() ! effective radius cloud ice
    !-----
    ! soil
    !-----
    real(wp) _POINTER :: t_so (:)   => NULL() ! soil temperature
    real(wp) _POINTER :: w_so (:)   => NULL() ! soil moisture
  end type t_col
!------------------------------------------------------------------------------
  !-----------------------------------------
  ! type t_cols: container for model columns
  !-----------------------------------------
  type t_cols
    integer              :: ncol        =  0      ! number of model columns
    integer              :: ke          =  0      ! number of levels
    integer              :: ns          =  0      ! number of soil layers
    integer              :: nt          =  nt     ! number of time ranges
    integer              :: tr (nt)     =  tr     ! range (hours)
    integer              :: nvar        =  0      ! number of variables
    integer(i8)          :: ids         =  0_i8   ! variables allocated
    integer              :: ntr         =  0      ! number of tracers
    integer              :: nm          =  0      ! multi level data size
    type (t_time)        :: time
    integer              :: levtyp      =  0      ! level type
    integer              :: vctype      =  0      ! vertical coordinate type
    real(wp)    _POINTER :: ak    (:)   => NULL() ! vertical coefficient table
    real(wp)    _POINTER :: bk    (:)   => NULL() ! vertical coefficient table
    type(t_col) ,pointer :: col   (:)   => NULL() ! columns
    real(wp)    _POINTER :: multi (:)   => NULL() ! multi level data
!   real(wp)    ,pointer :: tracs (:,:) => NULL() ! tracer
  end type t_cols
!==============================================================================
! Constants
!==========
  !-------------------------
  ! identifier for variables
  !-------------------------
  integer(i8), parameter :: COL_T      =         1_i8 ! temperature
  integer(i8), parameter :: COL_Q      =         2_i8 ! specific humidity
  integer(i8), parameter :: COL_RH     =         4_i8 ! relative humidity
! integer(i8), parameter :: COL_U      =         8_i8 ! u-wind component
! integer(i8), parameter :: COL_V      =        16_i8 ! v-wind component
  ! u,v shall always come in pairs.  Technically we currently must use 2 bits:
  integer(i8), parameter :: COL_UV     =   2_i8 ** 3 &! hor.wind (u+v components)
                                              + 2_i8 ** 4  ! hor.wind (u+v components)
  integer(i8), parameter :: COL_P      =        32_i8 ! pressure
  integer(i8), parameter :: COL_PH     =        64_i8 ! pressure     (half levels)
  integer(i8), parameter :: COL_GEO    =       128_i8 ! geopotential
  integer(i8), parameter :: COL_GEOH   =       256_i8 ! geopotential (half levels)
  integer(i8), parameter :: COL_TV     =       512_i8 ! virtual temperature
  integer(i8), parameter :: COL_X      =      1024_i8 ! water + ice load
  integer(i8), parameter :: COL_QCL    =      2048_i8 ! cloud water content
  integer(i8), parameter :: COL_QCI    =      4096_i8 ! cloud ice   content
  integer(i8), parameter :: COL_QR     =      8192_i8 ! rain
  integer(i8), parameter :: COL_QS     =     16384_i8 ! snow
  integer(i8), parameter :: COL_QG     =     32768_i8 ! graupel
  integer(i8), parameter :: COL_INFLAT =     65536_i8 ! inflation factor
  integer(i8), parameter :: COL_TV2    =    131072_i8 ! tv from tv, not from t, q
  integer(i8), parameter :: COL_IRPC   =    262144_i8 ! IR emissivity PCA coefficients
  integer(i8), parameter :: COL_PP     =    524288_i8 ! pressure deviation
  integer(i8), parameter :: COL_W      =   1048576_i8 ! vertical velocity
  integer(i8), parameter :: COL_SOIL   =   2097152_i8 ! soil moisture and temperature
  integer(i8), parameter :: COL_CLC    =   4194304_i8 ! cloud cover
  integer(i8), parameter :: COL_RANGE  =   8388608_i8 ! time range variables
  integer(i8), parameter :: COL_QVDIA  =  16777216_i8 ! specific humidity, diag.
  integer(i8), parameter :: COL_QCDIA  =  33554432_i8 ! spec. cloud water, diag.
  integer(i8), parameter :: COL_QIDIA  =  67108864_i8 ! spec. cloud ice, diag.
  integer(i8), parameter :: COL_OZONE  = 134217728_i8 ! ozone
  integer(i8), parameter :: COL_WSO    =     2_i8**28 ! soil moisture only
  integer(i8), parameter :: COL_CO2    =     2_i8**29 ! CO2
  integer(i8), parameter :: COL_REFF_QC=     2_i8**30 ! effective radius cloud water
  integer(i8), parameter :: COL_REFF_QI=     2_i8**31 ! effective radius cloud ice
  integer ,parameter :: N_COL      = 32      ! total number of identifiers
!==============================================================================
! Module variables
!=================
  !----------------------------------------------------------
  ! Use lowest model level as proxy for 10m wind?
  ! This is used for variational assimilation, until there is
  ! a forward model for deriving 10 minute averaged 10 wind.
  !----------------------------------------------------------
  logical :: u10_use_mlevel = .true.
  !-------------------------
  ! private module variables
  !-------------------------
  integer :: ncolsa  = 0      ! initialised variables of type t_cols
  integer :: ncola   = 0      ! allocated components t_cols% col
  integer :: nmulta  = 0      ! allocated components t_cols% multi
  integer :: ntraca  = 0      ! allocated components t_cols% tracs
  logical :: lprint = .false. ! flag for printing statistics
!==============================================================================
! Interfaces
!===========
  interface assignment (=)
    module procedure assign_cols_real
  end interface assignment (=)

  interface p_ps
    module procedure p_ps_cols
  end interface
!==============================================================================
! module Subroutines
!===================
contains
!------------------------------------------------------------------------------
  subroutine print_col_alloc (txt)
  character (len=*) ,intent(in) :: txt
    write(0,"(a,': pe=',i3,', ncols,ncol,nmult,ntrac=',4i10)") &
                      txt,dace% pe,ncolsa,ncola,nmulta,ntraca
  end subroutine print_col_alloc
!------------------------------------------------------------------------------
  subroutine send_cols (cols, dest)
  type (t_cols) ,intent(in) :: cols
  integer       ,intent(in) :: dest
    !------------------------------------------
    ! interfaces for derived type send routines
    !------------------------------------------
    interface
      subroutine p_send_derivedtype (buffer, count, dest, tag, comm)
        import :: t_cols
        type(t_cols)              :: buffer    ! variable to send
        integer, intent(in)       :: count     ! len(byte) of variable
        integer, intent(in)       :: dest      ! destination processor index
        integer, intent(in)       :: tag       ! tag
        integer, intent(in)       :: comm      ! communicator
      end subroutine p_send_derivedtype
      subroutine p_send_derivedtype2 (buffer, count, dest, tag, comm)
        import :: t_col
        type(t_col)                :: buffer(*) ! variable to send
        integer, intent(in)        :: count     ! len(byte) of variable
        integer, intent(in)        :: dest      ! destination processor index
        integer, intent(in)        :: tag       ! tag
        integer, intent(in)        :: comm      ! communicator
      end subroutine p_send_derivedtype2
    end interface
    !----------------
    ! local variables
    !----------------
    integer ,parameter :: tag = 1
    integer            :: size_cols
    integer            :: size_col
    !---------------
    ! send container
    !---------------
FTRACE_BEGIN("send_cols:container")
    size_cols = size(transfer(cols,['*']))
    call p_send_derivedtype (cols, size_cols, dest, tag, dace% comm)
FTRACE_END  ("send_cols:container")
    !----------------
    ! send components
    !----------------
    if (cols% ncol > 0) then
FTRACE_BEGIN("send_cols:components")
      size_col = size (transfer (cols% col(1),['*'])) * cols% ncol
      call p_send_derivedtype2 (cols% col, size_col, dest, tag, dace% comm)
      call p_send              (cols% multi,         dest, tag)
!     call p_send              (cols% tracs,         dest, tag)
FTRACE_END  ("send_cols:components")
    endif
  end subroutine send_cols
!------------------------------------------------------------------------------
  subroutine isend_cols (cols, dest, pass)
  type(t_cols) ,intent(in) :: cols
  integer      ,intent(in) :: dest      ! Destination pe
  integer      ,intent(in) :: pass      ! "Phase" of gather_cols
    !-------------------------------------------
    ! interfaces for derived type isend routines
    !-------------------------------------------
    interface
      subroutine p_isend_derivedtype (buffer, count, dest, tag, comm)
        import :: t_cols
        type(t_cols)               :: buffer    ! variable to send
        integer, intent(in)        :: count     ! len(byte) of variable
        integer, intent(in)        :: dest      ! destination processor index
        integer, intent(in)        :: tag       ! tag
        integer, intent(in)        :: comm      ! communicator
      end subroutine
      subroutine p_isend_derivedtype2 (buffer, count, dest, tag, comm)
        import :: t_col
        type(t_col)                :: buffer(*) ! variable to send
        integer, intent(in)        :: count     ! len(byte) of variable
        integer, intent(in)        :: dest      ! destination processor index
        integer, intent(in)        :: tag       ! tag
        integer, intent(in)        :: comm      ! communicator
      end subroutine
    end interface
    !----------------
    ! local variables
    !----------------
    integer ,parameter :: tag1 = 1, tag2 = 2, tag3 = 3
    integer            :: size_cols
    integer            :: size_col
    select case (pass)
    case (1)
       !---------------
       ! send container
       !---------------
       size_cols = size(transfer(cols,['*']))
       call p_isend_derivedtype (cols, size_cols, dest, tag1, dace% comm)
    case (2)
       !----------------
       ! send components
       !----------------
       if (cols% ncol > 0) then
          size_col = size (transfer (cols% col(1),['*'])) * cols% ncol
          call p_isend_derivedtype2 (cols% col, size_col, dest, tag2, dace% comm)
          call p_isend              (cols% multi,         dest, tag3)
!         call p_isend              (cols% tracs,         dest, tag)
       endif
    case default
       call finish ("isend_cols","invalid pass")
    end select
  end subroutine isend_cols
!------------------------------------------------------------------------------
  subroutine recv_cols (cols, src)
  type (t_cols) ,intent(inout) :: cols
  integer       ,intent(in)    :: src
    !----------------
    ! local variables
    !----------------
    integer ,parameter :: tag = 1
    integer            :: size_cols
    integer            :: size_col
    type (t_cols)      :: temp
    !---------------------------------------------
    ! include interfaces for external recv routine
    !---------------------------------------------
    interface
      subroutine p_recv_derivedtype (buffer, count, src, tag, comm)
        import :: t_cols
        type(t_cols)               :: buffer     ! variable to recv
        integer      ,intent(in)   :: count      ! len(byte) of variable
        integer      ,intent(in)   :: src        ! source processor el.
        integer      ,intent(in)   :: tag        ! tag
        integer      ,intent(in)   :: comm       ! communicator
      end subroutine p_recv_derivedtype
    end interface
    interface
      subroutine p_recv_derivedtype2 (buffer, count, src, tag, comm)
        import :: t_col
        type(t_col)                :: buffer (*) ! variable to recv
        integer      ,intent(in)   :: count      ! len(byte) of variable
        integer      ,intent(in)   :: src        ! source processor el.
        integer      ,intent(in)   :: tag        ! tag
        integer      ,intent(in)   :: comm       ! communicator
      end subroutine p_recv_derivedtype2
    end interface
    !------------------
    ! receive container
    !------------------
FTRACE_BEGIN("recv_cols:container")
    size_cols = size(transfer(cols,['*']))
    call p_recv_derivedtype (temp, size_cols, src, tag, dace% comm)
FTRACE_END  ("recv_cols:container")
    !--------------------
    ! allocate components
    !--------------------
    call alloc_cols (cols, tmp=temp, ncol=temp%ncol, ke=temp%ke, ns=temp%ns)
    !-------------------
    ! receive components
    !-------------------
    if (cols% ncol > 0) then
FTRACE_BEGIN("recv_cols:components")
      size_col = size (transfer (cols% col(1),['*'])) * cols% ncol
      call p_recv_derivedtype2 (cols% col, size_col, src, tag, dace% comm)
      call set_col_ptr (cols)
      call p_recv              (cols% multi,         src, tag)
!     call p_recv              (cols% tracs,         src, tag)
FTRACE_END  ("recv_cols:components")
    endif
  end subroutine recv_cols
!------------------------------------------------------------------------------
  subroutine irecv_cols (cols, src, pass, temp)
  type(t_cols) ,intent(inout) :: cols
  integer      ,intent(in)    :: src    ! Source pe
  integer      ,intent(in)    :: pass   ! "Phase" of gather_cols
  type(t_cols) ,intent(inout) :: temp   ! Temporary container between passes
    !----------------
    ! local variables
    !----------------
    integer ,parameter :: tag1 = 1, tag2 = 2, tag3 = 3
    integer            :: size_cols
    integer            :: size_col
    !----------------------------------------------
    ! include interfaces for external irecv routine
    !----------------------------------------------
    interface
      subroutine p_irecv_derivedtype (buffer, count, src, tag, comm)
        import :: t_cols
        type(t_cols)               :: buffer     ! variable to recv
        integer      ,intent(in)   :: count      ! len(byte) of variable
        integer      ,intent(in)   :: src        ! source processor el.
        integer      ,intent(in)   :: tag        ! tag
        integer      ,intent(in)   :: comm       ! communicator
      end subroutine p_irecv_derivedtype
    end interface
    interface
      subroutine p_irecv_derivedtype2 (buffer, count, src, tag, comm)
        import :: t_col
        type(t_col)                :: buffer (*) ! variable to recv
        integer      ,intent(in)   :: count      ! len(byte) of variable
        integer      ,intent(in)   :: src        ! source processor el.
        integer      ,intent(in)   :: tag        ! tag
        integer      ,intent(in)   :: comm       ! communicator
      end subroutine p_irecv_derivedtype2
    end interface

    select case (pass)
    case (1)
       !------------------
       ! receive container
       !------------------
       size_cols = size(transfer(cols,['*']))
       call p_irecv_derivedtype (temp, size_cols, src, tag1, dace% comm)
       !--------------------------------------------------------
       ! Make sure to mpi_wait on completion to get valid 'temp'
       !--------------------------------------------------------
    case (2)
       !--------------------
       ! allocate components
       !--------------------
       call alloc_cols (cols, tmp=temp, ncol=temp%ncol, ke=temp%ke)
       !-------------------
       ! receive components
       !-------------------
       if (cols% ncol > 0) then
          size_col = size (transfer (cols% col(1),['*'])) * cols% ncol
          call p_irecv_derivedtype2 (cols% col, size_col, src, tag2, dace% comm)
          call p_irecv              (cols% multi,         src, tag3)
!         call p_irecv              (cols% tracs,         src, tag)
!         call set_col_ptr (cols)   ! Defer until after completion
       endif
    case default
       call finish ("irecv_cols","invalid pass")
    end select
  end subroutine irecv_cols
!==============================================================================
! Gather columns.  The implementation using MPI_gather(v) might be faster,
! but tests with real data show that it is usually slightly slower than
! the straightforward implementation using explicit send/recv, since there
! are often only few partners exchanging larger amounts of data.
!------------------------------------------------------------------------------
#undef  USE_MPI_GATHER
#if defined (USE_MPI_GATHER)
  !-----------------------------------
  ! Implementation using MPI_gather(v)
  !-----------------------------------
  subroutine gather_cols (col, cols, root)
  use mo_mpi_dace, only: p_gather, p_gatherv
    type(t_cols) ,intent(in)  :: col
    type(t_cols) ,intent(out) :: cols(0:dace% npe-1)
    integer      ,intent(in)  :: root
    !----------------
    ! local variables
    !----------------
    integer                  :: n, nmulti
    integer                  :: j           ! Processor index
    integer                  :: rc  (0:dace% npe-1)
    integer                  :: ncol(0:dace% npe-1)
    type(t_cols)             :: temp(0:dace% npe-1)
    type(t_col), allocatable :: tcol(:)
    type(t_col),     pointer :: pcol(:)
    type(t_col),      target :: dummycol(0) ! Dummy array for bounds-checking
    real(wp),    allocatable :: multi(:)
    real(wp),        pointer :: pmulti(:)
    real(wp),         target :: dummyv(0)   ! Dummy array for bounds-checking
    !--------------------
    ! Gather container(s)
    !--------------------
FTRACE_BEGIN("gather_cols:p_gather_cols")
    call p_gather_cols ((/col/), temp, root=root)
FTRACE_END  ("gather_cols:p_gather_cols")
    !--------------------------------
    ! Allocate components on receiver
    !--------------------------------
    if (dace% pe == root) then
       do j = 0, dace% npe-1
          if (j == root) then
             cols(j) = col      ! Explicit copy on same PE
          else
             call alloc_cols (cols(j),                                     &
                  tmp=temp(j), ncol=temp(j)%ncol, ke=temp(j)%ke, ns=temp%ns)
          end if
       end do
       ncol = cols(:)% ncol
       allocate (tcol(sum (ncol)))
    else
       ncol = -1
       allocate (tcol(0))
    end if
    !---------------------------
    ! Gather relevant components
    !---------------------------
    if (col% ncol > 0) then
       pcol   => col% col
       pmulti => col% multi
       nmulti =  size (col% multi)
    else
       pcol   => dummycol
       pmulti => dummyv
       nmulti =  0
    end if
FTRACE_BEGIN("gather_cols:p_gather_col")
    call p_gather_col (pcol, tcol, recvcounts=ncol, root=root)
FTRACE_END  ("gather_cols:p_gather_col")
    if (dace% pe == root) then
       n = 0
       do j = 0, dace% npe-1
          if (ncol(j) > 0) then
             if (j /= root) then
                cols(j)%col = tcol(n+1:n+ncol(j))
                call set_col_ptr (cols(j))
             end if
          end if
          n = n + ncol(j)
       end do
    end if
    deallocate (tcol)
    call p_gather (nmulti, rc, root=root)
    if (dace% pe == root) then
       allocate (multi(sum (rc)))
    else
       allocate (multi(0))
    end if
FTRACE_BEGIN("gather_cols:p_gatherv")
    call p_gatherv (pmulti, multi, recvcounts=rc, root=root)
FTRACE_END  ("gather_cols:p_gatherv")
    if (dace% pe == root) then
       n = 0
       do j = 0, dace% npe-1
          if (rc(j) > 0) then
             if (j /= root) then
                cols(j)% multi = multi(n+1:n+rc(j))
             end if
          end if
          n = n + rc(j)
       end do
    end if
  end subroutine gather_cols
!------------------------------------------------------------------------------
#define DERIVED  type(t_cols)
#undef  MPI_TYPE
#define p_gather_DERIVED p_gather_cols
#include "p_gather_derived.incf"
#undef  DERIVED
#undef  MPI_TYPE
#undef  p_gather_DERIVED
!------------------------------------------------------------------------------
#define DERIVED  type(t_col)
#undef  MPI_TYPE
#define p_gather_DERIVED p_gather_col
#include "p_gather_derived.incf"
#undef  DERIVED
#undef  MPI_TYPE
#undef  p_gather_DERIVED
!------------------------------------------------------------------------------
#else
  !----------------------------------------
  ! Implementation using explicit send/recv
  !----------------------------------------
  subroutine gather_cols (col, cols, root)
    type(t_cols) ,intent(in)    :: col
    type(t_cols) ,intent(inout) :: cols(0:dace% npe-1)
    integer      ,intent(in)    :: root
    !----------------
    ! local variables
    !----------------
    integer :: j                ! loop index
#ifdef GATHER_COLS_ORIGINAL
    !-----------------------------
    ! loop over rank of partner PE
    !-----------------------------
    do j = 0,dace% npe-1
       !----------------------------------
       ! send if box is handled by partner
       !----------------------------------
       if (j == root) then
          if (j /= dace% pe) then
             call p_send (col, j)
          else if (j == dace% pe) then
             cols(j) = col
          end if
       end if
       !-------------------------------------
       ! receive if box is handled by this PE
       !-------------------------------------
       if (dace% pe == root) then
          if (j /= dace% pe) then
            call p_recv (cols(j), j)
          end if
       end if
    end do
#else /* simplified version of the above code */
    if (dace% pe == root) then
       cols(root) = col
       !-------------------------------------
       ! loop over rank of partner PE
       ! receive if box is handled by this PE
       !-------------------------------------
       do j = 0, dace% npe-1
          if (j /= root) then
             call p_recv (cols(j), j)
          end if
       end do ! j
    else
       !----------------------------------
       ! send if box is handled by partner
       !----------------------------------
       call p_send (col, root)
    end if
#endif
  end subroutine gather_cols
!------------------------------------------------------------------------------
#endif
!------------------------------------------------------------------------------
  !--------------------------------------------
  ! Implementation using non-blocking send/recv
  !--------------------------------------------
  subroutine igather_cols (col, cols, root)
    type(t_cols) ,intent(in)    :: col
    type(t_cols) ,intent(inout) :: cols(0:dace% npe-1)
    integer      ,intent(in)    :: root
    !----------------
    ! local variables
    !----------------
    integer      :: j                   ! loop index (pe)
    integer      :: pass                ! pass ("phase" of igather_cols)
    type(t_cols) :: temp(0:dace% npe-1)  ! Temporary container between passes

!write(0,*) dace% pe,root,0,"igather_cols::enter"
    if (dace% pe == root) then
       cols(root) = col
       !-------------------------------------
       ! receive if box is handled by this PE
       !-------------------------------------
       do pass = 1, 2
          do j = 0, dace% npe-1
             if (j /= root) then
                call irecv_cols (cols(j), j, pass=pass, temp=temp(j))
             end if
          end do ! j
!write(0,*) dace% pe,root,pass,"igather_cols::wait(irecv)"
          call p_waitall ()
       end do    ! pass
       do j = 0, dace% npe-1
          if (j /= root .and. cols(j)% ncol > 0) then
             call set_col_ptr (cols(j))
          end if
       end do ! j
    else
       !----------------------------------
       ! send if box is handled by partner
       !----------------------------------
       do pass = 1, 2
          call isend_cols (col, root, pass=pass)
!write(0,*) dace% pe,root,pass,"igather_cols::wait(isend)"
          call p_waitall ()
       end do    ! pass
    end if
!write(0,*) dace% pe,root,3,"igather_cols::leave"
!call p_barrier () ! Debug
  end subroutine igather_cols
!------------------------------------------------------------------------------
  subroutine gather_cols_multi (col, cols, root)
    type(t_cols) ,intent(in)    :: col (:)
    type(t_cols) ,intent(inout) :: cols(:,:)
    integer      ,intent(in)    :: root(:)
    !-----------------------------
    ! Gather columns for all boxes
    !-----------------------------
!!$    integer :: ib       ! "Box" index
    !---------------
    ! Sanity checks:
    !---------------
    if (size (col)    /= size (root) .or. &
        size (cols,2) /= size (root) .or. &
        size (cols,1) /= dace% npe        ) then
       write(0,*) "gather_cols_multi: size(col),shape(cols),size(root)=", &
            size (col), shape (cols), size (root)
       call finish ("gather_cols_multi","bad shape of args")
    end if

!!$    do ib = 1, size (root)
!!$!      call  gather_cols (col(ib), cols(:,ib), root(ib))
!!$       call igather_cols (col(ib), cols(:,ib), root(ib))
!!$    end do
    call igather_cols_multi (col, cols, root)
  end subroutine gather_cols_multi
!------------------------------------------------------------------------------
  !--------------------------------------------
  ! Implementation using non-blocking send/recv
  ! (use only by wrapper with argument checks.)
  !--------------------------------------------
  subroutine igather_cols_multi (col, cols, root)
    type(t_cols) ,intent(in)    :: col (   :)
    type(t_cols) ,intent(inout) :: cols(0:,:)
    integer      ,intent(in)    :: root(   :)
    !----------------
    ! local variables
    !----------------
    integer      :: ib                            ! box  index
    integer      :: pe                            ! root pe of box
    integer      :: j                             ! inner loop index (sender pe)
    integer      :: pass                          ! pass (c.f. igather_cols)
    type(t_cols) :: temp(0:dace% npe-1,size(col)) ! Temporary container
    !-----------------------------------------
    ! Limit number of outstanding MPI requests
    !-----------------------------------------
    integer, parameter :: maxwait = 32

    do pass = 1, 2
       do ib = 1, size (col)
          pe = root(ib)
          if (dace% pe == pe) then
             cols(pe,ib) = col(ib)
             !-------------------------------------
             ! receive if box is handled by this PE
             !-------------------------------------
             do j = 0, dace% npe-1
                if (j /= pe) then
                   call irecv_cols (cols(j,ib), j, pass=pass, temp=temp(j,ib))
                end if
             end do ! j
          else
             !----------------------------------
             ! send if box is handled by partner
             !----------------------------------
             call isend_cols (col(ib), pe, pass=pass)
          end if
!         if (mod (ib, maxwait) == 0) call p_waitall ()
          if (p_irequest > p_mrequest) call p_waitall ()
       end do ! ib
       call p_waitall ()
!      call p_barrier () ! Debug
    end do    ! pass

    do ib = 1, size (col)
       pe = root(ib)
       if (dace% pe == pe) then
          !--------------------------
          ! box is handled by this PE
          !--------------------------
          do j = 0, dace% npe-1
             if (j /= pe .and. cols(j,ib)% ncol > 0) then
                call set_col_ptr (cols(j,ib))
             end if
          end do ! j
       end if
    end do       ! ib

  end subroutine igather_cols_multi
!------------------------------------------------------------------------------
  subroutine scatter_cols (col, cols, root)
  !----------------------------------------
  ! send columns from all observation PEs
  ! to one of the model PEs
  ! Implementation using explicit send/recv
  ! Used in the adjoint code
  !----------------------------------------
    type(t_cols) ,intent(out) :: col
    type(t_cols) ,intent(in)  :: cols(0:dace% npe-1)
    integer      ,intent(in)  :: root
    !----------------
    ! local variables
    !----------------
    integer :: j                ! loop index
    !-----------------------------
    ! loop over rank of partner PE
    !-----------------------------
    do j = 0,dace% npe-1
       !----------------------------------
       ! send if box is handled by this PE
       !----------------------------------
       if (dace% pe == root) then
          if (j /= dace% pe) then
            call p_send (cols(j), j)
          end if
       end if
       !------------------------------------
       ! receve if box is handled by partner
       !------------------------------------
       if (j == root) then
          if (j /= dace% pe) then
             call p_recv (col, j)
          else if (j == dace% pe) then
             col = cols(j)
          end if
       end if
    end do
  end subroutine scatter_cols
!------------------------------------------------------------------------------
  subroutine alltoall_cols (in, out)
    !----------------------------------
    ! redistribute columns from all PEs
    !----------------------------------
    type(t_cols), intent(in)  :: in (0:) ! Send buffer (-> pe=0...)
    type(t_cols), intent(out) :: out(0:) ! Recv.buffer (<- pe=0...)
    !----------------
    ! local variables
    !----------------
    integer      :: pe                   ! loop index
    integer      :: pass                 ! pass ("phase" of alltoall_cols)
    type(t_cols) :: temp(0:dace% npe-1)  ! Temporary container between passes

    if (size (in)  /= dace% npe) call finish ("alltoall_cols","bad size(in)")
    if (size (out) /= dace% npe) call finish ("alltoall_cols","bad size(out)")
    !---------------------------
    ! Copy data local to this PE
    !---------------------------
    out(dace% pe) = in(dace% pe)
    !--------------------------------
    ! receive cols handled by this PE
    ! send cols handled by partner
    !--------------------------------
    do pass = 1, 2
       do pe = 0, dace% npe - 1
          if      (dace% pe <  pe) then
             call isend_cols (in(pe), dest=pe, pass=pass)
             call irecv_cols (out(pe), src=pe, pass=pass, temp=temp(pe))

          else if (dace% pe == pe) then ! Local data handled separately

          else if (dace% pe >  pe) then
             call irecv_cols (out(pe), src=pe, pass=pass, temp=temp(pe))
             call isend_cols (in(pe), dest=pe, pass=pass)
          endif
       end do ! pe
       call p_waitall ()
    end do    ! pass

    do pe = 0, dace% npe-1
       if (pe /= dace% pe .and. out(pe)% ncol > 0) then
          call set_col_ptr (out(pe))
       end if
    end do

  end subroutine alltoall_cols
!==============================================================================
  subroutine get_cols_gg (mc, atm, cols, iatm)
  type (t_mcols)    ,intent(in)           :: mc    ! model column meta data (in boxes)
  type (t_atm)      ,intent(in)           :: atm   ! atmospheric state (on model grid)
  type (t_cols)     ,intent(out)          :: cols  ! model columns, sorted in boxes
  integer(i8)       ,intent(in) ,optional :: iatm  ! fields required (forced)
  !---------------------------------------------------------------------------
  ! distribute model columns over PEs for interpolation-operator.
  !
  ! This routine is intended for GRID to GRID interpolation. In contrast
  ! to GRID to OBSERVATION interpolation the columns on the target PEs
  ! are not sorted in 'boxes'. In order to use the routine 'get_cols' which is
  ! tailored to the latter application we set up dummy 'boxes' (one for each
  ! PE) and call 'get_cols'.
  !---------------------------------------------------------------------------

    type (t_mcols) :: mc1   (0:dace% npe-1)
    type (t_cols)  :: cols1 (0:dace% npe-1)
    integer :: i
    mc1(dace% pe) = mc
    mc1% pe = (/(i,i=0,dace% npe-1)/)
    call get_cols (mc1, atm, cols1, iatm)
    cols = cols1(dace% pe)

  end subroutine get_cols_gg
!------------------------------------------------------------------------------
  subroutine get_cols1 (mc, atm, cols, iatm)
  type (t_mcols)   ,intent(in)           :: mc  (:) ! model column meta data (in boxes)
  type (t_atm)     ,intent(in)           :: atm (:) ! atmospheric state (on model grid)
  type (t_cols)    ,intent(out)          :: cols(:) ! model columns, sorted in boxes
  integer(i8)      ,intent(in) ,optional :: iatm    ! fields required (forced)
  !---------------------------------------------------------------------------
  ! distribute model columns over PEs for (observation-)interpolation-operator
  ! - version for single time slot -
  !---------------------------------------------------------------------------

    integer :: t ! time slot

    do t = 1, size (atm)
      call get_cols (mc, atm(t), cols, iatm, t)
    end do

  end subroutine get_cols1
!------------------------------------------------------------------------------
  subroutine get_cols (mc, atm, cols, iatm, tslot)
  type (t_mcols)             ,intent(in)    :: mc  (:) ! model column meta data (in boxes)
  type (t_atm)               ,intent(in)    :: atm     ! atmospheric state (on model grid)
  type (t_cols)              ,intent(inout) :: cols(:) ! model columns, sorted in boxes
  integer(i8)      ,optional ,intent(in)    :: iatm    ! fields required (forced)
  integer          ,optional ,intent(in)    :: tslot   ! time slot
  !---------------------------------------------------------------------------
  ! distribute model columns over PEs for (observation-)interpolation-operator
  ! - version for multiple time slots -
  !
  ! 1) The meta data (column indices required for interpolation in the
  !    'boxes' at the target PE) is gathered on the source PEs.
  ! 2) On the source PE the required model columns are stored in a array
  !    of derived type t_cols for redistribution.
  ! 3) The model columns are transmitted to the target PEs
  ! 4) The model columns received from the different PEs are sorted in
  !    desired order in the 'boxes'
  !
  ! The model columns on the target PEs are sorted in 'boxes' (corresponding
  ! with the array elements of the subroutine arguments 'mc' and 'cols') there
  ! a given 'box' is allocated only on one PE, as required by the PSAS
  ! algorithm. If present, the variable 'iatm' (bit flags) specifies the model
  ! variables (temperature, humidity,...) to be transmitted in all grid
  ! columns in addition to those specified individually within 'mc'.
  !---------------------------------------------------------------------------

    !----------------
    ! local variables
    !----------------
    integer               :: t          ! time slot
    integer               :: i, j, k, l ! indices
    integer               :: ntri       ! number of boxes
    integer               :: nmcol      ! number of model columns in a box
    integer(i8)           :: latm       ! bit flags
    integer(i8)           :: ids      (               size (mc)) ! Variables
    type (p_mcol)         :: mcols    (               size (mc))
    type (t_cols)         :: cols_send(               size (mc))
    type (t_cols)         :: cols_recv(0:dace% npe-1, size (mc))
    type (t_col) ,pointer :: cl, cr

    !--------------------------------------------------
    ! derive some quantities from subroutine parameters
    !--------------------------------------------------
    t    = 1     ;if (present (tslot)) t    = tslot
    latm = COL_P ;if (present (iatm )) latm = iatm
    ntri = size (mc)
    if (t==1) cols% ncol = 0
    !----------------------------------------------------------------
    ! 1) The meta data (column indices required for interpolation in
    !    the 'boxes' at the target PE) is gathered on the source PEs.
    !
    !    This is basically an alltoall communication pattern.
    !    Using alltoall communication directly may be more efficient.
    !----------------------------------------------------------------
    call alltoall_mcol (mc(:), mcols, t)
    !---------------------------------------------------------------
    ! 2) On the source PE the required model columns are stored into
    !    an array of derived type t_cols for redistribution.
    !---------------------------------------------------------------
    do i=1,ntri
      if (mc(i)% pe == dace% pe) then
         nmcol  = mc(i)% n
         ids(i) = latm
         do l=1,nmcol
            ids(i) = ior (ids(i), mc(i)% c(l)% iatm)
         end do
      end if
      call atm2col (cols_send(i), atm, mcols(i)%p, iatm=iatm)
      deallocate (mcols(i)%p)
    end do
    !----------------------------------------------------------------
    ! 3) The model columns are transmitted to the target PEs
    !
    !    This is basically an alltoall communication pattern.
    !    Using alltoall communication directly may be more efficient.
    !----------------------------------------------------------------
#ifdef GET_COLS_ORIGINAL
    do i=1,ntri
FTRACE_BEGIN("get_cols:gather_cols")
      call gather_cols (cols_send(i), cols_recv(:,i), root=mc(i)% pe)
!     call igather_cols (cols_send(i), cols_recv(:,i), root=mc(i)% pe)
FTRACE_END  ("get_cols:gather_cols")
    end do
#else
FTRACE_BEGIN("get_cols:gather_cols_multi")
    call gather_cols_multi (cols_send, cols_recv, root=mc(:)% pe)
FTRACE_END  ("get_cols:gather_cols_multi")
#endif
    !-----------------------------------------------------
    ! 4) The model columns received from the different PEs
    !    are sorted in desired order in the 'boxes'
    !-----------------------------------------------------
    do i=1,ntri
      !-------------------------------------
      ! process if box is handled by this PE
      !-------------------------------------
      if (mc(i)% pe == dace% pe) then
        nmcol = mc(i)% n
        !-----------------------------
        ! loop over rank of partner PE
        !-----------------------------
        do j=0,dace% npe-1
          if (cols(i)% ncol == 0 .and. cols_recv(j,i)% ncol> 0) then
            !-------------------------
            ! allocate result variable
            !-------------------------
            call alloc_cols (cols(i),              &
                             tmp = cols_recv(j,i), &
                             ncol= nmcol,          &
                             ids = ids(i)          )
          else
            cols(i)% time = cols_recv(j,i)% time
          endif
          l = 0
          do k=1,nmcol
            if (mc(i)% c(k)% ijdtp(5) == j .and. &
                mc(i)% c(k)% ijdtp(4) == t       ) then
              l = l + 1
              cl => cols(i)% col(k)        ! left  hand side
              cr => cols_recv(j,i)% col(l) ! right hand side
              cl% c     = cr% c
              cl% i     = cr% i
              cl% j     = cr% j
              cl% l     = cr% l
              cl% pe    = cr% pe
              cl% s     = cr% s
              call copy_col  (cl% t    ,cr% t   )
              call copy_col  (cl% q    ,cr% q   )
              call copy_col  (cl% u    ,cr% u   )
              call copy_col  (cl% v    ,cr% v   )
!             call copy_col2 (cl% tr   ,cr% tr  )
              call copy_col  (cl% p    ,cr% p   )
              call copy_col  (cl% ph   ,cr% ph  )
              call copy_col  (cl% geo  ,cr% geo )
              call copy_col  (cl% geoh ,cr% geoh)
              call copy_col  (cl% rh   ,cr% rh  )
              call copy_col  (cl% tv   ,cr% tv  )
              call copy_col  (cl% x    ,cr% x   )
              call copy_col  (cl% qcl  ,cr% qcl )
              call copy_col  (cl% qci  ,cr% qci )
              call copy_col  (cl% qr   ,cr% qr  )
              call copy_col  (cl% qs   ,cr% qs  )
              call copy_col  (cl% qg   ,cr% qg  )
              call copy_col  (cl% infl ,cr% infl)
              call copy_col  (cl% pp   ,cr% pp  )
              call copy_col  (cl% w    ,cr% w   )
              call copy_col  (cl% t_so ,cr% t_so)
              call copy_col  (cl% w_so ,cr% w_so)
              call copy_col  (cl% clc  ,cr% clc )
              call copy_col  (cl% qv_dia, cr% qv_dia )
              call copy_col  (cl% qc_dia, cr% qc_dia )
              call copy_col  (cl% qi_dia, cr% qi_dia )
              call copy_col  (cl% o3   ,cr% o3  )
              call copy_col  (cl% co2  ,cr% co2 )
              call copy_col  (cl% reff_qi,cr% reff_qi )
              call copy_col  (cl% reff_qc,cr% reff_qc )
            endif
          end do
          if (j /= dace% pe) call dealloc_cols (cols_recv(j,i))
        end do
      end if
      call dealloc_cols (cols_send(i))
    end do
  contains

    subroutine copy_col (y,x)
      real(wp) _POINTER :: y(:),x(:)
      if (associated(y)) then
        if(associated(x)) then
          y = x
        else
          nullify (y)
        endif
      endif
    end subroutine copy_col

!   subroutine copy_col2 (y,x)
!   real(wp), pointer :: y(:,:),x(:,:)
!     if (associated(y)) then
!       if(associated(x)) then
!         y = x
!       else
!         nullify (y)
!       endif
!     endif
!   end subroutine copy_col2

  end subroutine get_cols
!------------------------------------------------------------------------------
  subroutine get_cols_gg_adj (mc, atm, cols, iatm)
  type (t_mcols)    ,intent(in)           :: mc    ! model column meta data (in boxes)
  type (t_atm)      ,intent(inout)        :: atm   ! atmospheric state (on model grid)
  type (t_cols)     ,intent(in)           :: cols  ! model columns, sorted in boxes
  integer(i8)       ,intent(in) ,optional :: iatm  ! fields required (forced)

    type (t_mcols) :: mc1   (0:dace% npe-1)
    type (t_cols)  :: cols1 (0:dace% npe-1)
    integer :: i
    mc1% pe     = (/(i,i=0,dace% npe-1)/)
    mc1  (dace% pe) = mc
    cols1(dace% pe) = cols
    call get_cols_adj (mc1, atm, cols1, iatm)

  end subroutine get_cols_gg_adj
!------------------------------------------------------------------------------
  subroutine get_cols_adj (mc, atm, cols, iatm)
  type (t_mcols)   ,intent(in)           :: mc  (:) ! model column meta data (in boxes)
  type (t_atm)     ,intent(inout)        :: atm     ! atmospheric state (on model grid)
  type (t_cols)    ,intent(in)           :: cols(:) ! model columns, sorted in boxes
  integer(i8)      ,intent(in) ,optional :: iatm    ! fields required (forced)
  !-----------------------------------------
  ! distribute model columns over PEs
  ! for (observation-)interpolation-operator
  !-----------------------------------------

    !----------------
    ! local variables
    !----------------
    integer               :: i, j, k, l ! indices
    integer               :: ntri       ! number of boxes
    integer               :: nmcol      ! number of columns in a box
    integer(i8)           :: latm       ! bit flags
    integer(i8)           :: ids      (               size (mc)) ! Variables
    type (p_mcol)         :: mcols    (               size (mc))
    type (t_cols)         :: cols_send(               size (mc))
    type (t_cols)         :: cols_recv(0:dace% npe-1, size (mc))
    type (t_col) ,pointer :: cl, cr
    !----------------------------------------------
    ! send request for model columns to hosting PEs
    !----------------------------------------------
    ntri = size (mc)
    latm = COL_P; if (present (iatm)) latm = iatm
    !----------------------
    ! scatter the meta data
    !----------------------
    do i=1,ntri
FTRACE_BEGIN("get_cols_adj:scatter_mcol")
      nullify (mcols(i)%p)
      call scatter_mcol (mc(i), mcols(i)%p, 1)
FTRACE_END  ("get_cols_adj:scatter_mcol")
    end do
    !----------------
    ! prepare columns
    !----------------
    do i=1,ntri
      if (mc(i)% pe == dace% pe) then
         nmcol  = mc(i)% n
         ids(i) = latm
         do l=1,nmcol
            ids(i) = ior (ids(i), mc(i)% c(l)% iatm)
         end do
      end if
    end do

    !-------------------
    ! reorganise columns
    !-------------------
    do i=1,ntri
      !-------------------------------------
      ! process if box is handled by this PE
      !-------------------------------------
      if (mc(i)% pe == dace% pe) then
        nmcol = mc(i)% n
        !-----------------------------
        ! loop over rank of partner PE
        !-----------------------------
        do j=0,dace% npe-1
!         if (cols(i)% ncol == 0 .and. cols_recv(j,i)% ncol> 0) then
            !-------------------------
            ! allocate result variable
            !-------------------------
            call alloc_cols (cols_recv(j,i),       &
                             tmp = cols(i),        &
                             ncol= nmcol,          &
                             ids = ids(i)          )
!         else
!           cols(i)% time = cols_recv(j,i)% time
!         endif
          l = 0
          do k=1,nmcol
            if (mc(i)% c(k)% ijdtp(5) == j) then
              l = l + 1
              cr => cols(i)% col(k)        ! right hand side
              cl => cols_recv(j,i)% col(l) ! left  hand side
              cl% c     = cr% c
              cl% i     = cr% i
              cl% j     = cr% j
              cl% l     = cr% l
              cl% pe    = cr% pe
              cl% s     = cr% s
              call copy_col  (cl% u    ,cr% u   )
              call copy_col  (cl% v    ,cr% v   )
              call copy_col  (cl% geo  ,cr% geo )
              call copy_col  (cl% rh   ,cr% rh  )
              call copy_col  (cl% tv   ,cr% tv  )
            endif
          end do
        end do
      end if
    end do

    !-------------
    ! scatter data
    !-------------
    do i=1,ntri
FTRACE_BEGIN("get_cols_adj:scatter_cols")
      call scatter_cols (cols_send(i), cols_recv(:,i), root=mc(i)% pe)
FTRACE_END  ("get_cols_adj:scatter_cols")
        do j=0,dace% npe-1
          if (j /= dace% pe) call dealloc_cols (cols_recv(j,i))
        end do
    end do

    !-------------------------
    ! add to atmospheric state
    !-------------------------
    do i=1,ntri
      call atm2col_adj (cols_send(i), atm, mcols(i)%p, iatm=iatm)
      deallocate (mcols(i)%p)
      call dealloc_cols (cols_send(i))
    end do


  contains

    subroutine copy_col (y,x)
    real(wp), pointer :: y(:),x(:)
      if (associated(y)) then
        if(associated(x)) then
          y = x
        else
          nullify (y)
        endif
      endif
    end subroutine copy_col

  end subroutine get_cols_adj
!==============================================================================
  subroutine atm2col (cols, atm, idx, iatm)
  type (t_cols)     ,intent(out)           :: cols    ! columns to set
  type (t_atm)      ,intent(in)            :: atm     ! gridded atmospheric field
  type (t_mcol)     ,intent(in)            :: idx (:) ! index array
  integer(i8)       ,intent(in) ,optional  :: iatm
    !----------------
    ! local variables
    !----------------
    type (t_grid) ,pointer :: g                ! model grid
    type (t_col)  ,pointer :: c                ! column pointer
    integer                :: ncol             ! number of model columns to set
    integer                :: i1,i2,id,l       ! indices
    integer(i8)            :: ids              ! variables to allocate
    integer(i8)            :: latm
    real(wp)    ,parameter :: rdp = 180._wp/pi ! conversion factor
    g    => atm% grid
    ncol =  size (idx)
    latm = COL_P; if (present (iatm)) latm = iatm
    !+++++++++++++++++++++++++++++++
    ! allocate all columns requested
    !+++++++++++++++++++++++++++++++
    ids = latm
    do l=1,ncol
      ids  = ior (ids, idx(l)% iatm)
    end do
    !-----------------
    ! allocate columns
    !-----------------
    call alloc_cols (cols, ncol=ncol, ke=g% nz, ns=g% ns, ids=ids)
    !--------------
    ! copy metadata
    !--------------
    cols% time    = atm%   time
    cols% ke      = g%     nz
    cols% levtyp  = g%     levtyp
    cols% vctype  = g%     vct
    select case (g% levtyp)
    case (WMO3_ISOBARIC)
      cols% ak = g% akf (:g% nz+1)
      cols% bk = 0._wp
    case default
      cols% ak = g% ak (:g% nz+1)
      cols% bk = g% bk (:g% nz+1)
    end select
    !----------------------------
    ! check for allocation status
    !----------------------------
    call  check_alloc (COL_T       ,atm% t       ,'t   '  )
    call  check_alloc (COL_TV2     ,atm% tv      ,'tv  '  )
    call  check_alloc (COL_TV      ,atm% t       ,'t   '  )
    call  check_alloc (COL_TV      ,atm% q       ,'q   '  )
    call  check_alloc (COL_Q       ,atm% q       ,'q   '  )
    call  check_alloc (COL_RH      ,atm% rh      ,'rh  '  )
    call  check_alloc (COL_UV      ,atm% u       ,'u   '  )
    call  check_alloc (COL_UV      ,atm% v       ,'v   '  )
    call  check_alloc (COL_P       ,atm% pf      ,'p   '  )
    call  check_alloc (COL_PH      ,atm% ph      ,'ph  '  )
    call  check_alloc (COL_GEO     ,atm% geof    ,'geo '  )
    call  check_alloc (COL_GEOH    ,atm% geoh    ,'geoh'  )
    call  check_alloc (COL_QCL     ,atm% qcl     ,'qcl'   )
    call  check_alloc (COL_QCI     ,atm% qci     ,'qci'   )
    call  check_alloc (COL_QR      ,atm% qr      ,'qr'    )
    call  check_alloc (COL_QS      ,atm% qs      ,'qs'    )
    call  check_alloc (COL_QG      ,atm% qg      ,'qg'    )
    call  check_alloc (COL_INFLAT  ,atm% f_inflat,'inflat')
    call  check_alloc (COL_PP      ,atm% pp      ,'pp'    )
    call  check_alloc (COL_W       ,atm% w       ,'w'     )
    call  check_alloc (COL_SOIL    ,atm% t_so    ,'t_so'  )
    call  check_alloc (COL_SOIL    ,atm% w_so    ,'w_so'  )
    call  check_alloc (COL_WSO     ,atm% w_so    ,'w_so'  )
    call  check_alloc (COL_CLC     ,atm% clc     ,'clc'   )
    call  check_alloc (COL_QVDIA   ,atm% qv_dia  ,'qv_dia')
    call  check_alloc (COL_QCDIA   ,atm% qc_dia  ,'qc_dia')
    call  check_alloc (COL_QIDIA   ,atm% qi_dia  ,'qi_dia')
    call  check_alloc (COL_OZONE   ,atm% o3      ,'o3'    )
    call  check_alloc (COL_CO2     ,atm% co2     ,'co2'   )
    call  check_alloc (COL_REFF_QC ,atm% reff_qc ,'reff_qc')
    call  check_alloc (COL_REFF_QI ,atm% reff_qi ,'reff_qi')
    !----------
    ! copy data
    !----------
    do l =  1, ncol
      i1 =  idx(l)% ijdtp (1)
      i2 =  idx(l)% ijdtp (2)
      id =  idx(l)% ijdtp (3)
      c  => cols% col (l)
      c%c% dlat = rdp * g% rlat (i1, i2, 1, id)
      c%c% dlon = rdp * g% rlon (i1, i2, 1, id)
      c%   i    = i1
      c%   j    = i2
      c%   l    = id
      call set_xuv (c)
      c% s% ps       =   0._wp
      c% s% psr      =   0._wp
      c% s% dpsdt    =   0._wp
      c% s% geoid    =   0._wp
      c% s% geosp    =   0._wp
      c% s% ts       =   0._wp
!     c% s% qs       =   0._wp
      c% s% lsm      = -99._wp
      c% s% fr_ice   =   0._wp
      c% s% ssd      = -99._wp
      c% s% z0       = -99._wp
      c% s% tll      =   0._wp
      c% s% rhll     = -99._wp
      c% s% tsurf    =   0._wp
      c% s% t2m      = -99._wp
      c% s% rh2m     = -99._wp
      c% s% td2m     = -99._wp
      c% s% t2mland  = -99._wp
      c% s% rh2mland = -99._wp
      c% s% td2mland = -99._wp
      c% s% dtdz     = -99._wp
      c% s% t_ice    = -99._wp
      c% s% h_ice    = -99._wp
      c% s% t_so_0   = -99._wp
      c% s% h_snow   = -99._wp
      c% s% t_snow   = -99._wp
      c% s% w_snow   = -99._wp
      c% s% rho_snow = -99._wp
      c% s% freshsnw = -99._wp
      c% s% snowc    = -99._wp
      c% s% w_i      = -99._wp
      c% s% ull      = -99._wp
      c% s% vll      = -99._wp
      c% s% u10m     = -99._wp
      c% s% v10m     = -99._wp
      c% s% clct     = -99._wp
      c% s% clcl     = -99._wp
      c% s% clcm     = -99._wp
      c% s% clch     = -99._wp
      c% s% ceiling  = -99._wp
      c% s% vis      = -99._wp

      c% s% vmax_10m = -99._wp
      c% s% tmin_2m  = -99._wp
      c% s% tmax_2m  = -99._wp
      c% s% tot_prec = -99._wp
      c% s% aswdir_s = -99._wp
      c% s% aswdifd_s= -99._wp

      if (associated (atm%       u_10m)) &
        c% s% u10m  = atm%       u_10m (i1, i2, 1,     id)
      if (associated (atm%       v_10m)) &
        c% s% v10m  = atm%       v_10m (i1, i2, 1,     id)
      if (associated (atm%       u)) &
        c% s% ull   = atm%       u     (i1, i2, g% nz, id)
      if (associated (atm%       v)) &
        c% s% vll   = atm%       v     (i1, i2, g% nz, id)

      if ((u10_use_mlevel .and. c% s% ull  /= -99._wp) .or. &
                                c% s% u10m == -99._wp     ) &
        c% s% u10m  = c% s% ull
      if ((u10_use_mlevel .and. c% s% vll  /= -99._wp) .or. &
                                c% s% v10m == -99._wp     ) &
        c% s% v10m  = c% s% vll

      if (associated (atm%       ps)) &
        c% s% ps    = atm%       ps    (i1, i2, 1,     id)
      if (associated (atm%       psr)) &
        c% s% psr   = atm%       psr   (i1, i2, 1,     id)
      if (associated (atm%       dpsdt)) &
        c% s% dpsdt = atm%       dpsdt (i1, i2, 1,     id)
      if (associated (atm% grid% lsm)) &
        c% s% lsm   = atm% grid% lsm   (i1, i2, 1,     id)
      if (associated (atm% grid% sso_stdh)) &
        c% s% ssd   = atm% grid% sso_stdh (i1, i2, 1,  id)
      if (associated (atm%       fr_ice)) &
        c% s% fr_ice= atm%       fr_ice(i1, i2, 1,     id)
      if (associated (atm%       z0)) &
        c% s% z0    = atm%       z0    (i1, i2, 1,     id)
      if (associated (atm%       t2m)) &
        c% s% t2m   = atm%       t2m   (i1, i2, 1,     id)
      if (associated (atm%       rh2m)) &
        c% s% rh2m  = atm%       rh2m  (i1, i2, 1,     id)
      if (associated (atm%       td2m)) &
        c% s% td2m  = atm%       td2m  (i1, i2, 1,     id)

      if (associated (atm%       t2m_land )) &
        c% s% t2mland  = atm%    t2m_land (i1, i2, 1,  id)
      if (associated (atm%       rh2m_land)) &
        c% s% rh2mland = atm%    rh2m_land(i1, i2, 1,  id)
      if (associated (atm%       td2m_land)) &
        c% s% td2mland = atm%    td2m_land(i1, i2, 1,  id)

      if (associated (atm%       t)) then
        c% s% tll   = atm%       t     (i1, i2, g% nz, id)
        if(g% kmlbot /= 0 .and. g% kmltop /= 0 .and. &
           associated(atm%       pf)           .and. &
           associated(atm%       t)            .and. &
           associated(  g%       hhl)                ) then
        if (9999._wp /= g%       hhl   (i1, i2, g% kmltop,  id))&
          c%s%dtdz  =(atm%       t     (i1, i2, g% kmltop,  id) &
                    - atm%       t     (i1, i2, g% kmlbot,  id))&
                    /(  g%       hhl   (i1, i2, g% kmltop,  id) &
                       +g%       hhl   (i1, i2, g% kmltop+1,id) &
                       -g%       hhl   (i1, i2, g% kmlbot,  id) &
                       -g%       hhl   (i1, i2, g% kmlbot+1,id))&
                    * 2._wp
        endif
        if(associated(atm%       q ) .and. &
           associated(atm%       t ) .and. &
           associated(atm%       pf)     ) then
          ! dont apply to analysis increments:
          if (atm% t(i1, i2, g% nz, id) > 100._wp)         &
          c%s%rhll  = rhw_q(atm% q     (i1, i2, g% nz, id),&
                            atm% t     (i1, i2, g% nz, id),&
                            atm% pf    (i1, i2, g% nz, id))
       endif
      endif
      if (associated    (atm%    tsurf)) &
        c% s% tsurf    = atm%    tsurf    (i1, i2, 1,     id)
      if (associated    (atm%    t_ice)) &
        c% s% t_ice    = atm%    t_ice    (i1, i2, 1,     id)
      if (associated    (atm%    h_ice)) &
        c% s% h_ice    = atm%    h_ice    (i1, i2, 1,     id)
      if (associated    (atm%    h_snow)) &
        c% s% h_snow   = atm%    h_snow   (i1, i2, 1,     id)
      if (associated    (atm%    t_snow)) &
        c% s% t_snow   = atm%    t_snow   (i1, i2, 1,     id)
      if (associated    (atm%    w_snow)) &
        c% s% w_snow   = atm%    w_snow   (i1, i2, 1,     id)
      if (associated    (atm%    rho_snow)) &
        c% s% rho_snow = atm%    rho_snow (i1, i2, 1,     id)
      if (associated    (atm%    freshsnw)) &
        c% s% freshsnw = atm%    freshsnw (i1, i2, 1,     id)
      if (associated    (atm%    snowc)) &
        c% s% snowc    = atm%    snowc    (i1, i2, 1,     id)
      if (associated    (atm%    w_i)) &
        c% s% w_i      = atm%    w_i      (i1, i2, 1,     id)
      if (associated    (atm%    t_so )) &
        c% s% t_so_0   = atm%    t_so     (i1, i2, 0,     id)
      if (associated    (atm%    clct)) &
        c% s% clct     = atm%    clct     (i1, i2, 1,     id)
      if (associated    (atm%    clcl)) &
        c% s% clcl     = atm%    clcl     (i1, i2, 1,     id)
      if (associated    (atm%    clcm)) &
        c% s% clcm     = atm%    clcm     (i1, i2, 1,     id)
      if (associated    (atm%    clch)) &
        c% s% clch     = atm%    clch     (i1, i2, 1,     id)
      if (associated    (atm%    ceiling)) &
        c% s% ceiling  = atm%    ceiling  (i1, i2, 1,     id)
      if (associated    (atm%    vis)) &
        c% s% vis      = atm%    vis      (i1, i2, 1,     id)
      if (associated    (atm% grid% geosp)) &
        c% s% geosp =    atm% grid% geosp    (i1, i2, 1,     id)
      if (associated    (atm% grid% geoid)) &
        c% s% geoid =    atm% grid% geoid    (i1, i2, 1,     id)
      if (iand(ids,COL_T     )/=0_i8) c% t    = atm% t        (i1,i2, :,id)
      if (iand(ids,COL_TV2   )/=0_i8) c% tv   = atm% tv       (i1,i2, :,id)
      if (iand(ids,COL_TV    )/=0_i8) c% tv   = tv_t_q        (             &
                                                atm% t        (i1,i2, :,id),&
                                                atm% q        (i1,i2, :,id) )
      if (iand(ids,COL_Q     )/=0_i8) c% q    = atm% q        (i1,i2, :,id)
      if (iand(ids,COL_RH    )/=0_i8) c% rh   = atm% rh       (i1,i2, :,id)
      if (iand(ids,COL_UV    )/=0_i8) c% u    = atm% u        (i1,i2, :,id)
      if (iand(ids,COL_UV    )/=0_i8) c% v    = atm% v        (i1,i2, :,id)
      if (iand(ids,COL_P     )/=0_i8) c% p    = atm% pf       (i1,i2, :,id)
      if (iand(ids,COL_PH    )/=0_i8) c% ph   = atm% ph       (i1,i2, :,id)
      if (iand(ids,COL_GEO   )/=0_i8) c% geo  = atm% geof     (i1,i2, :,id)
      if (iand(ids,COL_GEOH  )/=0_i8) c% geoh = atm% geoh     (i1,i2, :,id)
      if (iand(ids,COL_QCL   )/=0_i8) c% qcl  = atm% qcl      (i1,i2, :,id)
      if (iand(ids,COL_QCI   )/=0_i8) c% qci  = atm% qci      (i1,i2, :,id)
      if (iand(ids,COL_QR    )/=0_i8) c% qr   = atm% qr       (i1,i2, :,id)
      if (iand(ids,COL_QS    )/=0_i8) c% qs   = atm% qs       (i1,i2, :,id)
      if (iand(ids,COL_QG    )/=0_i8) c% qg   = atm% qg       (i1,i2, :,id)
      if (iand(ids,COL_INFLAT)/=0_i8) c% infl = atm% f_inflat (i1,i2, :,id)
      if (iand(ids,COL_PP    )/=0_i8) c% pp   = atm% pp       (i1,i2, :,id)
      if (iand(ids,COL_W     )/=0_i8) c% w    = atm% w        (i1,i2, :,id)
      if (iand(ids,COL_SOIL  )/=0_i8) c% t_so = atm% t_so     (i1,i2,1:,id)
      if (iand(ids,COL_SOIL +    &
                   COL_WSO   )/=0_i8) c% w_so = atm% w_so     (i1,i2, :,id)
      if (iand(ids,COL_CLC   )/=0_i8) c% clc  = atm% clc      (i1,i2, :,id)
      if (iand(ids,COL_QVDIA )/=0_i8) c% qv_dia = atm% qv_dia (i1,i2, :,id)
      if (iand(ids,COL_QCDIA )/=0_i8) c% qc_dia = atm% qc_dia (i1,i2, :,id)
      if (iand(ids,COL_QIDIA )/=0_i8) c% qi_dia = atm% qi_dia (i1,i2, :,id)
      if (iand(ids,COL_REFF_QC)/=0_i8)c% reff_qc= atm% reff_qc(i1,i2, :,id)
      if (iand(ids,COL_REFF_QI)/=0_i8)c% reff_qi= atm% reff_qi(i1,i2, :,id)
      if (iand(ids,COL_X     )/=0_i8) then
                                 c% x    = 0._wp
        if(associated(atm% qcl)) c%x=c%x + atm% qcl  (i1,i2,:,id)
        if(associated(atm% qci)) c%x=c%x + atm% qci  (i1,i2,:,id)
        if(associated(atm% qr )) c%x=c%x + atm% qr   (i1,i2,:,id)
        if(associated(atm% qs )) c%x=c%x + atm% qs   (i1,i2,:,id)
        if(associated(atm% qg )) c%x=c%x + atm% qg   (i1,i2,:,id)
      endif
      if (iand(ids,COL_RANGE )/=0_i8) then
        if(associated(atm% vmax_10m )) call tr2col (c%s% vmax_10m  ,VMAX_10M )
        if(associated(atm% tmin_2m  )) call tr2col (c%s% tmin_2m   ,TMIN_2M  )
        if(associated(atm% tmax_2m  )) call tr2col (c%s% tmax_2m   ,TMAX_2M  )
        if(associated(atm% tot_prec )) call tr2col (c%s% tot_prec  ,TOT_PREC )
        if(associated(atm% aswdir_s )) call tr2col (c%s% aswdir_s  ,ASWDIR_S )
        if(associated(atm% aswdifd_s)) call tr2col (c%s% aswdifd_s ,ASWDIFD_S)
      endif
      if (iand(ids,COL_OZONE )/=0_i8) c% o3  = atm% o3      (i1,i2, :,id)
      if (iand(ids,COL_CO2   )/=0_i8) c% co2 = atm% co2     (i1,i2, :,id)
    end do
  contains
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    subroutine tr2col (c, k)
    !--------------------------------------------------------------------------
    ! assign variable valid over a time range
    ! constrain to values >= 0 (clip rounding errors, quantities should be >=0)
    !--------------------------------------------------------------------------
    real(wp) ,intent(inout) :: c (nt)  ! column variable
    integer  ,intent(in)    :: k       ! index in atmospheric state vector

      integer :: i, j
      do j = 1, atm% m(k)% i% ub(3)
        do i = 1, nt
          if (tr(i) == atm% m(k)% i% tr(j)) &
            c(i) = max (atm% m(k)% ptr (i1,i2,j,id), 0._wp)
        end do
      end do
    end subroutine tr2col
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    subroutine check_alloc (flag, ptr, name)
    !---------------------------------------------
    ! check allocation status of atmospheric field
    !---------------------------------------------
    integer(i8)      ,intent(in) :: flag
    real(wp)         ,pointer    :: ptr (:,:,:,:)
    character(len=*) ,intent(in) :: name

      if (iand(ids, flag)/=0_i8 .and. .not. (associated (ptr))) &
        call finish ('atm2col','not associated: '//name)
    end subroutine check_alloc
  end subroutine atm2col
!------------------------------------------------------------------------------
  subroutine atm2col_adj (cols, atm, idx, iatm)
  type (t_cols)     ,intent(in)           :: cols    ! columns to set
  type (t_atm)      ,intent(inout)        :: atm     ! gridded atmospheric field
  type (t_mcol)     ,intent(in)           :: idx (:) ! index array
  integer(i8)       ,intent(in) ,optional :: iatm
    !----------------
    ! local variables
    !----------------
!   type (t_grid) ,pointer :: g                ! model grid
    type (t_col)  ,pointer :: c                ! column pointer
    integer                :: ncol             ! number of model columns to set
    integer                :: i1,i2,id,l       ! indices
    integer(i8)            :: ids              ! variables to allocate
    integer(i8)            :: latm
!   real(wp)    ,parameter :: rdp = 180._wp/pi ! conversion factor
!   g    => atm% grid
    ncol =  size (idx)
    latm = COL_P; if (present (iatm)) latm = iatm
    !+++++++++++++++++++++++++++++++
    ! allocate all columns requested
    !+++++++++++++++++++++++++++++++
    ids = latm
    do l=1,ncol
      ids  = ior (ids, idx(l)% iatm)
    end do
    !----------------------------
    ! check for allocation status
    !----------------------------
    call  check_alloc (COL_TV2     ,atm% tv      ,'tv  '  )
    call  check_alloc (COL_RH      ,atm% rh      ,'rh  '  )
    call  check_alloc (COL_UV      ,atm% u       ,'u   '  )
    call  check_alloc (COL_UV      ,atm% v       ,'v   '  )
    call  check_alloc (COL_GEO     ,atm% geof    ,'geo '  )
    !----------
    ! copy data
    !----------
!print *,dace% pe,'### atm2col_adj ncol',ncol, size (cols% col)
    do l =  1, ncol
      i1 =  idx(l)% ijdtp (1)
      i2 =  idx(l)% ijdtp (2)
      id =  idx(l)% ijdtp (3)
      c  => cols% col (l)
      if (iand(ids,COL_TV2)/=0_i8) atm% tv   (i1,i2,:,id) =      &
                                   atm% tv   (i1,i2,:,id) + c% tv
      if (iand(ids,COL_RH )/=0_i8) atm% rh   (i1,i2,:,id) =      &
                                   atm% rh   (i1,i2,:,id) + c% rh
      if (iand(ids,COL_UV )/=0_i8) atm% u    (i1,i2,:,id) =      &
                                   atm% u    (i1,i2,:,id) + c% u
      if (iand(ids,COL_UV )/=0_i8) atm% v    (i1,i2,:,id) =      &
                                   atm% v    (i1,i2,:,id) + c% v
      if (iand(ids,COL_GEO)/=0_i8) atm% geof (i1,i2,:,id) =      &
                                   atm% geof (i1,i2,:,id) + c% geo
    end do
  contains
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    subroutine check_alloc (flag, ptr, name)
    integer(i8)      ,intent(in) :: flag
    real(wp)         ,pointer    :: ptr (:,:,:,:)
    character(len=*) ,intent(in) :: name
    !---------------------------------------------
    ! check allocation status of atmospheric field
    !---------------------------------------------
      if (iand(ids, flag)/=0 .and. .not. (associated (ptr))) &
        call finish ('atm2col_adj','not associated: '//name)
    end subroutine check_alloc
  end subroutine atm2col_adj
!==============================================================================
  subroutine alloc_cols (cols, ncol, ke, ns, ntr, ids, tmp)
  type(t_cols) ,intent(out)           :: cols ! columns to allocate
  integer      ,intent(in)  ,optional :: ncol ! number of model columns
  integer      ,intent(in)  ,optional :: ke   ! number of levels
  integer      ,intent(in)  ,optional :: ns   ! number of soil layers
  integer      ,intent(in)  ,optional :: ntr  ! number of tracers
  integer(i8)  ,intent(in)  ,optional :: ids  ! specify variables to allocate
  type(t_cols) ,intent(in)  ,optional :: tmp  ! template
    !----------------
    ! local variables
    !----------------
    integer              :: i         ! multi level data index
    integer              :: nm        ! size of multi level data / column
    integer              :: lke       ! number of model columns
    integer              :: lns       ! number of soil layers
!   integer              :: lntr      ! number of tracers
!   integer              :: lids      ! specify variables to allocate
    !--------------
    ! set metadata
    !--------------
    if (present (tmp)) then
      !------------------------
      ! use template if present
      !------------------------
      cols% ncol    = tmp% ncol
      cols% ke      = tmp% ke
      cols% ns      = tmp% ns
      cols% nvar    = tmp% nvar
      cols% ids     = tmp% ids
      cols% ntr     = tmp% ntr
      cols% nm      = tmp% nm
      cols% time    = tmp% time
      cols% levtyp  = tmp% levtyp
      cols% vctype  = tmp% vctype
    endif
    if (present (ncol)) cols% ncol = ncol
    if (present (ke))   cols% ke   = ke
    if (present (ns))   cols% ns   = ns
    if (present (ntr))  cols% ntr  = ntr
    lke  = cols% ke
    lns  = cols% ns
!   lntr = cols% ntr
    !------------------------------------------------
    ! determine total number of parameters per column
    !------------------------------------------------
    if (present (ids))  then
      cols% ids  = ids
      nm = 0
!     nm = 2 * (lke+1)
      if (btest(ids, 0)) nm = nm + lke     ! col_t
      if (btest(ids, 1)) nm = nm + lke     ! col_q
      if (btest(ids, 2)) nm = nm + lke     ! col_rh
      if (btest(ids, 3)) nm = nm + lke * 2 ! col_uv
      if (btest(ids, 5)) nm = nm + lke     ! col_p
      if (btest(ids, 6)) nm = nm + lke + 1 ! col_ph
      if (btest(ids, 7)) nm = nm + lke     ! col_geo
      if (btest(ids, 8)) nm = nm + lke + 1 ! col_geoh
      if (btest(ids, 9)) nm = nm + lke     ! col_tv
      if (btest(ids,10)) nm = nm + lke     ! col_x
      if (btest(ids,11)) nm = nm + lke     ! col_qcl
      if (btest(ids,12)) nm = nm + lke     ! col_qci
      if (btest(ids,13)) nm = nm + lke     ! col_qr
      if (btest(ids,14)) nm = nm + lke     ! col_qs
      if (btest(ids,15)) nm = nm + lke     ! col_qg
      if (btest(ids,16)) nm = nm + lke     ! col_inflat
      if (btest(ids,17)) nm = nm + lke     ! col_tv2
!     if (btest(ids,18)) nm = ?            ! col_irpc
      if (btest(ids,19)) nm = nm + lke     ! col_pp
      if (btest(ids,20)) nm = nm + lke + 1 ! col_w
      if (btest(ids,21)) nm = nm + lns     ! col_soil
      if (btest(ids,21) .or.              &! col_soil
          btest(ids,28)) nm = nm + lns     ! col_wso
      if (btest(ids,22)) nm = nm + lke     ! col_clc
      if (btest(ids,24)) nm = nm + lke     ! col_qvdia
      if (btest(ids,25)) nm = nm + lke     ! col_qcdia
      if (btest(ids,26)) nm = nm + lke     ! col_qidia
      if (btest(ids,27)) nm = nm + lke     ! col_ozone
      if (btest(ids,29)) nm = nm + lke     ! col_co2
      if (btest(ids,30)) nm = nm + lke     ! col_reff_qi
      if (btest(ids,31)) nm = nm + lke     ! col_reff_qc
      cols% nm   = nm
    endif
!   lids = cols% ids
    !----------------------------------
    ! allocate pointer array components
    !----------------------------------
    allocate (cols% col   (         cols% ncol            ))
    allocate (cols% multi (2*lke+2+ cols% ncol * cols% nm ))
!   allocate (cols% tracs (lke,     cols% ncol * cols% ntr))
    !-------------------
    ! associate pointers
    !-------------------
    call set_col_ptr (cols)
    !-------------------
    ! copy from template
    !-------------------
    if (present (tmp)) then
      !-------------------------------
      ! copy column specific meta data
      !-------------------------------
      if (.not. present (ncol)) then
        do i = 1, cols% ncol
          cols% col(i)% c    = tmp% col(i)% c
          cols% col(i)% i    = tmp% col(i)% i
          cols% col(i)% j    = tmp% col(i)% j
          cols% col(i)% l    = tmp% col(i)% l
          cols% col(i)% pe   = tmp% col(i)% pe
        end do
      endif
      !--------------------------------
      ! copy vertical coefficient table
      !--------------------------------
      if (.not. present (ke)) then
        if (associated (tmp% ak)) then
          cols% levtyp = tmp% levtyp
          cols% ak     = tmp% ak
          cols% bk     = tmp% bk
        endif
      endif
    endif
    !--------------------------------
    ! bookkeeping on allocated memory
    !--------------------------------
    ncolsa = ncolsa + 1
    ncola  = ncola  + cols% ncol
    nmulta = nmulta + cols% ncol * cols% nm
    ntraca = ntraca + cols% ntr
    if (lprint) call print_col_alloc ('  allocate cols')
  end subroutine alloc_cols
!------------------------------------------------------------------------------
  elemental subroutine set_col_ptr (cols)
  type (t_cols) ,intent(inout) :: cols

    integer               :: i, it, l, ntr, ke, ns
    integer(i8)           :: ids
    type (t_col) ,pointer :: c
    real(wp)     _POINTER :: p(:)
!   real(wp)     ,pointer :: t(:,:)

    i   = 0
    it  = 0
    ids = cols% ids
    ntr = cols% ntr
    ke  = cols% ke
    ns  = cols% ns
    p  => cols% multi
!   t  => cols% tracs
    cols% ak => p(i+1 : i+ke+1); i = i+ke+1
    cols% bk => p(i+1 : i+ke+1); i = i+ke+1
    do l=1,cols% ncol
      c => cols% col(l)
      c% iom = i
      c% iot = it
      if (btest(ids, 0)) then; c%t   =>p(i+1:i+ke  ); i=i+ke   ;endif
      if (btest(ids, 1)) then; c%q   =>p(i+1:i+ke  ); i=i+ke   ;endif
      if (btest(ids, 2)) then; c%rh  =>p(i+1:i+ke  ); i=i+ke   ;endif
      if (btest(ids, 3)) then; c%u   =>p(i+1:i+ke  ); i=i+ke
                               c%v   =>p(i+1:i+ke  ); i=i+ke   ;endif
      if (btest(ids, 5)) then; c%p   =>p(i+1:i+ke  ); i=i+ke   ;endif
      if (btest(ids, 6)) then; c%ph  =>p(i+1:i+ke+1); i=i+ke+1 ;endif
      if (btest(ids, 7)) then; c%geo =>p(i+1:i+ke  ); i=i+ke   ;endif
      if (btest(ids, 8)) then; c%geoh=>p(i+1:i+ke+1); i=i+ke+1 ;endif
      if (btest(ids, 9)) then; c%tv  =>p(i+1:i+ke  ); i=i+ke   ;endif
      if (btest(ids,10)) then; c%x   =>p(i+1:i+ke  ); i=i+ke   ;endif
      if (btest(ids,11)) then; c%qcl =>p(i+1:i+ke  ); i=i+ke   ;endif
      if (btest(ids,12)) then; c%qci =>p(i+1:i+ke  ); i=i+ke   ;endif
      if (btest(ids,13)) then; c%qr  =>p(i+1:i+ke  ); i=i+ke   ;endif
      if (btest(ids,14)) then; c%qs  =>p(i+1:i+ke  ); i=i+ke   ;endif
      if (btest(ids,15)) then; c%qg  =>p(i+1:i+ke  ); i=i+ke   ;endif
      if (btest(ids,16)) then; c%infl=>p(i+1:i+ke  ); i=i+ke   ;endif
      if (btest(ids,17)) then; c%tv  =>p(i+1:i+ke  ); i=i+ke   ;endif

      if (btest(ids,19)) then; c%pp  =>p(i+1:i+ke  ); i=i+ke   ;endif
      if (btest(ids,20)) then; c%w   =>p(i+1:i+ke+1); i=i+ke+1 ;endif
      if (btest(ids,21)) then; c%t_so=>p(i+1:i+ns  ); i=i+ns   ;endif
      if (btest(ids,21) .or. &
          btest(ids,28)) then; c%w_so=>p(i+1:i+ns  ); i=i+ns   ;endif
      if (btest(ids,22)) then; c%clc =>p(i+1:i+ke  ); i=i+ke   ;endif
      if (btest(ids,24)) then; c%qv_dia=>p(i+1:i+ke); i=i+ke   ;endif
      if (btest(ids,25)) then; c%qc_dia=>p(i+1:i+ke); i=i+ke   ;endif
      if (btest(ids,26)) then; c%qi_dia=>p(i+1:i+ke); i=i+ke   ;endif
      if (btest(ids,27)) then; c%o3  =>p(i+1:i+ke  ); i=i+ke   ;endif
      if (btest(ids,29)) then; c%co2 =>p(i+1:i+ke  ); i=i+ke   ;endif
      if (btest(ids,30)) then; c%reff_qc=>p(i+1:i+ke);i=i+ke   ;endif
      if (btest(ids,31)) then; c%reff_qi=>p(i+1:i+ke);i=i+ke   ;endif
!      if (ntr > 0) then
!        c%tr => t(:,it+1:it+ntr)
!      else  ! workaround for bug in IBM xlf 8.1.1.1 (-C option)
!        c%tr => t
!      endif ! workaround for bug in IBM xlf 8.1.1.1 (-C option)
      it=it+ntr
    end do
  end subroutine set_col_ptr
!------------------------------------------------------------------------------
  subroutine dealloc_col (cols)
  type(t_cols) ,intent(inout) :: cols ! columns to deallocate
    type(t_cols) ,save :: empty
    integer            :: stat
    !--------------------------------
    ! bookkeeping on allocated memory
    !--------------------------------
    ncolsa = ncolsa - 1
    ncola  = ncola  - cols% ncol
    nmulta = nmulta - cols% ncol * cols% nm
    ntraca = ntraca - cols% ntr
    if (lprint) call print_col_alloc ('deallocate cols')
    !-----------------------+++++++++++++++++++++++++++++++++++++++++++++
    ! deallocate components ! stat for NAG -mtrace option, compiler bug ?
    !-----------------------+++++++++++++++++++++++++++++++++++++++++++++
    if (associated (cols% col  )) deallocate (cols% col   ,stat=stat)
    if (associated (cols% multi)) deallocate (cols% multi ,stat=stat)
!   if (associated (cols% tracs)) deallocate (cols% tracs ,stat=stat)
    !---------------
    ! zero data type
    !---------------
    cols = empty
  end subroutine dealloc_col
!------------------------------------------------------------------------------
  subroutine dealloc_cols (cols)
  type(t_cols) ,intent(inout) :: cols(:) ! columns to deallocate
    integer :: i
    do i = 1, size (cols)
      call dealloc_col (cols(i))
    end do
  end subroutine dealloc_cols
!==============================================================================
  subroutine set_xuv_col (col)
  type (t_col) ,intent(inout) :: col
  !----------------------------------
  ! set up unit vectors on the sphere
  !----------------------------------
    call set_xuv (col% c)
  end subroutine set_xuv_col
!==============================================================================
  subroutine assign_cols_real (y,x)
  type (t_cols) ,intent(inout) :: y
  real (wp)     ,intent(in)    :: x
    !------------------
    ! single level data
    !------------------
    y% col% s% ps       = x
    y% col% s% geosp    = x
    y% col% s% geoid    = x
    y% col% s% ts       = x
!   y% col% s% qs       = x
    y% col% s% lsm      = x
    y% col% s% ssd      = x
    y% col% s% z0       = x
    y% col% s% tll      = x
    y% col% s% rhll     = x
    y% col% s% ull      = x
    y% col% s% vll      = x
    y% col% s% u10m     = x
    y% col% s% v10m     = x
    y% col% s% t2m      = x
    y% col% s% rh2m     = x
    y% col% s% td2m     = x
    y% col% s% t2mland  = x
    y% col% s% rh2mland = x
    y% col% s% td2mland = x
    y% col% s% dtdz     = x
    y% col% s% tsurf    = x
    y% col% s% fr_ice   = x
    y% col% s% t_ice    = x
    y% col% s% h_ice    = x
    y% col% s% t_so_0   = x
    y% col% s% h_snow   = x
    y% col% s% t_snow   = x
    y% col% s% w_snow   = x
    y% col% s% rho_snow = x
    y% col% s% freshsnw = x
    y% col% s% snowc    = x
    y% col% s% w_i      = x
    y% col% s% clct     = x
    y% col% s% clcl     = x
    y% col% s% clcm     = x
    y% col% s% clch     = x
    y% col% s% ceiling  = x
    y% col% s% vis      = x
!   y% col% s% vmax_10m = x
!   y% col% s% tmin_2m  = x
!   y% col% s% tmax_2m  = x
!   y% col% s% tot_prec = x
!   y% col% s% aswdir_s = x
!   y% col% s% aswdifd_s= x
    !-----------------
    ! multi level data
    !-----------------
    y% multi (2*y%ke+3:) = x
!   y% tracs (:,:)       = x
  end subroutine assign_cols_real
!==============================================================================
  subroutine p_ps_cols (cols)
  type (t_cols) ,intent(inout) :: cols

    real(wp) :: ph (cols% ke+1)
    real(wp) :: pf (cols% ke)

    integer :: i
    do i=1, cols% ncol
      if (associated (cols% col(i)% ph) .or. &
          associated (cols% col(i)% p )) then
        call p_ps (ph, pf, cols% ak, cols% bk, cols% col(i)% s% ps)
        if (associated (cols% col(i)% ph)) cols% col(i)% ph = ph
        if (associated (cols% col(i)% p )) cols% col(i)% p  = pf
      endif
    end do
  end subroutine p_ps_cols
!==============================================================================
  subroutine set_cols_logp (cols)
  !------------------------
  ! store ln p instead of p
  !------------------------
  type(t_cols), intent(inout) :: cols(:)
    integer :: i
    do i = 1, size (cols)
      call set_cols_logp1 (cols(i))
    end do
  end subroutine set_cols_logp
!------------------------------------------------------------------------------
  subroutine set_cols_logp1 (cols)
  !------------------------
  ! store ln p instead of p
  !------------------------
  type(t_cols), intent(inout) :: cols
    integer :: j
    do j = 1, cols% ncol
      if (associated (cols% col(j)% p)) &
        cols% col(j)% p  = log (cols% col(j)% p)
      if (associated (cols% col(j)% ph)) then
        if (cols% col(j)% ph(1) == 0._wp)                     &! GME/HRM
            cols% col(j)% ph(1) = cols% col(j)% ph(2)/10._wp
        cols% col(j)% ph = log (cols% col(j)% ph)
      endif
    end do
  end subroutine set_cols_logp1
!==============================================================================
  subroutine set_mcols_m (mc, ic, gy, gx, iatm)
  !-------------------------------------------
  ! set t_mcols for grid to grid interpolation
  !-------------------------------------------
  type (t_mcols)   ,intent(out) :: mc (dace% npe)
  type (t_icol)    ,pointer     :: ic (:,:,:)
  type (t_grid)    ,intent(in)  :: gy            ! destination grid
  type (t_grid)    ,intent(in)  :: gx            ! source      grid
  integer(i8)      ,intent(in)  :: iatm          ! field to gather flags

    !----------------
    ! local variables
    !----------------
    real(wp) ,parameter   :: rdp = 180._wp / pi ! degree - radiant factor
    integer               :: i,j,l              ! indices
    integer               :: pp,ii,jj,ll        ! indices
    integer               :: ly(4)              ! lower bounds destination grid
    integer               :: uy(4)              ! upper bounds destination grid
    integer               :: lx(4)              ! lower bounds source      grid
    integer               :: ux(4)              ! upper bounds source      grid
    integer               :: p1
    integer               :: n
    type(t_mcol) ,pointer :: c (:)
    logical               :: reordered          ! same grid but re-ordered points
    integer               :: pe_i1_i2_id(4)     ! re-ordered indices
    !---------------------------------------------
    ! determine columns required for interpolation
    !---------------------------------------------
    p1 =  dace% pe+1
    ly =  gy% lb
    uy =  gy% ub
    lx =  gx% lbg
    ux =  gx% ubg
    allocate (ic          (ly(1):uy(1),ly(2):uy(2),ly(4):uy(4)  ))
    allocate (mc(p1)% idx (lx(1):ux(1),lx(2):ux(2),lx(4):ux(4),1))
    allocate (mc(p1)% c (gx% nxny / dace% npe))
    mc(p1)% idx = 0
    mc    % n   = 0
    mc    % pe  = (/ (i,i=0,dace% npe-1) /)
    ic%c% dlon  = gy% rlon(ly(1):uy(1),ly(2):uy(2),1,ly(4):uy(4)) * rdp
    ic%c% dlat  = gy% rlat(ly(1):uy(1),ly(2):uy(2),1,ly(4):uy(4)) * rdp

    n = 0
    if (same_horizontal_grid (gx, gy, reordered=reordered)) then
      !------------------------------------
      ! for same grid the indices are known
      !------------------------------------
      if (.not. reordered) then
        do l = ly(4), uy(4)
          do j = ly(2), uy(2)
            do i = ly(1), uy(1)
              n            = n + 1
              if (size (mc(p1)% c) < n) then
                c => mc(p1)% c
                allocate (mc(p1)% c (int(1.5 * n)))
                mc(p1)% c (1:size(c)) = c
                deallocate (c)
              endif
              mc(p1)% n               = n
              mc(p1)% c(n)% icol      = n
              mc(p1)% c(n)% ijdtp(1)  = i
              mc(p1)% c(n)% ijdtp(2)  = j
              mc(p1)% c(n)% ijdtp(3)  = l
              mc(p1)% c(n)% ijdtp(4)  = 1
              mc(p1)% c(n)% ijdtp(5)  = dace% pe
              mc(p1)% c(n)% iatm      = iatm
              mc(p1)% c(n)% itrac     = 0

              ic(i,j,l)% h% imc(1,1)  = n
              ic(i,j,l)% h% w  (1)    = 1._wp
              ic(i,j,l)% h% ijdp(1:3) = mc(p1)% c(n)% ijdtp(1:3)
              ic(i,j,l)% h% ijdp(  4) = mc(p1)% c(n)% ijdtp(  5)

            end do
          end do
        end do
      else  ! re-ordered
        !-------------------------------------------------------------
        ! same grid, but grid-points of the source grid are re-ordered
        !-------------------------------------------------------------
        if (gy% d_gme(1) < 0 .or. gx%  d_gme(1) >= 0)          &
          call finish ('set_mcols_m','inconsistent re-ordering')
        !------------------------------------------------------------------
        ! same grid, but grid-points of the destination grid are re-ordered
        !------------------------------------------------------------------
        do l = gy% lbg(4), gy% ubg(4)
          do j = gy% lbg(2), gy% ubg(2)
            do i = gy% lbg(1), gy% ubg(1)
              if (gy% marr (1,i,j,l) == dace% pe) then
                pp          = gx% marr (1,i,j,l)
                pe_i1_i2_id = gy% marr (:,i,j,l)
                ii = pe_i1_i2_id(2)
                jj = pe_i1_i2_id(3)
                ll = pe_i1_i2_id(4)
                n            = n + 1
                if (size (mc(p1)% c) < n) then
                  c => mc(p1)% c
                  allocate (mc(p1)% c (int(1.5 * n)))
                  mc(p1)% c (1:size(c)) = c
                  deallocate (c)
                endif

                mc(p1)% n               = n
                mc(p1)% c(n)% icol      = n
                mc(p1)% c(n)% ijdtp     = [i,j,l,1,pp]
                mc(p1)% c(n)% iatm      = iatm
                mc(p1)% c(n)% itrac     = 0

                ic(ii,jj,ll)% h% imc(1,1)  = n
                ic(ii,jj,ll)% h% w  (1)    = 1._wp
                ic(ii,jj,ll)% h% ijdp(1:3) = mc(p1)% c(n)% ijdtp(1:3)
                ic(ii,jj,ll)% h% ijdp(  4) = mc(p1)% c(n)% ijdtp(  5)

              endif
            end do
          end do
        end do
      endif
    else
      !----------------------------------------------------
      ! for different grid search the neighbour grid-points
      !----------------------------------------------------
      do l = ly(4), uy(4)
        do j = ly(2), uy(2)
          do i = ly(1), uy(1)
            call idx_init ( &!
              ic(i,j,l)% c, &! <-  column descriptor
              ic(i,j,l)% h, &!  -> interpolation coefficients
              mc(p1),       &! <-> model column descriptors
              iatm,         &! <-  fields required
              0,            &! <-  tracers required
              gx,           &! <-  model grid
              1,            &! <-  time slot
              0._wp         )! <-  temporal weight
            if (ic(i,j,l)% h% imc(1,1) == 0) then
              write(6,*) dace% pe,'set_mcols_m: out of domain condition'
              write(6,*) dace% pe,'   gx% type =',gx% gridtype
              write(6,*) dace% pe,'   gy% type =',gy% gridtype
              write(6,*) dace% pe,'   gx% d_gme=',gx% d_gme(1)
              write(6,*) dace% pe,'   gy% d_gme=',gy% d_gme(1)
              write(6,*) dace% pe,'   i,j,l lb =',ly(1),ly(2),ly(4)
              write(6,*) dace% pe,'   i,j,l ub =',uy(1),uy(2),uy(4)
              write(6,*) dace% pe,'   i,j,l    =',i,j,l
              write(6,*) dace% pe,'   idx      =',ic(i,j,l)% h% imc
              write(6,*) dace% pe,'   lat, lon =',ic(i,j,l)% c% dlat, &
                                                  ic(i,j,l)% c% dlon
              write(6,*) dace% pe,'   iw12     =',ic(i,j,l)% h% iw12
              call finish('set_mcols_m','out of domain condition')
            endif
          end do
        end do
      end do
    endif
    deallocate (mc(p1)% idx)
  end subroutine set_mcols_m
!==============================================================================
  subroutine mcol_idx (mc, ic, idx, gx, iatm)
  !-----------------------------------------------
  ! set t_mcols from t_h_idx derived type variable
  !-----------------------------------------------
  type (t_mcols) ,intent(inout) :: mc  (dace% npe)! <-> list of colomns on PE
  type (t_icol)  ,pointer       :: ic  (:,:,:)    !  ->
  type (t_h_idx) ,pointer       :: idx (:,:,:)    !  <-
  type (t_grid)  ,intent(in)    :: gx             !  <- source      grid
  integer(i8)    ,intent(in)    :: iatm           !  <- parameter bit flag

    integer               :: i,j,l,k,n,m        ! indices
    integer               :: ly(4)              ! lower bounds destination grid
    integer               :: uy(4)              ! upper bounds destination grid
    integer               :: lx(4)              ! lower bounds source      grid
    integer               :: ux(4)              ! upper bounds source      grid
    integer               :: p1
    integer               :: ix
    type(t_mcol) ,pointer :: c, tmp(:)
    integer               :: ijdp (4)

    ly(1) = lbound (idx, 1)
    ly(2) = lbound (idx, 2)
    ly(4) = lbound (idx, 3)
    uy(1) = ubound (idx, 1)
    uy(2) = ubound (idx, 2)
    uy(4) = ubound (idx, 3)
    lx    =  gx% lbg
    ux    =  gx% ubg
    p1    =  dace% pe+1

    allocate   (ic      (ly(1):uy(1),ly(2):uy(2),ly(4):uy(4)  ))

    if (.not.associated (mc(p1)% idx)) then
      allocate (mc(p1)% idx (lx(1):ux(1),lx(2):ux(2),lx(4):ux(4),1))
      allocate (mc(p1)% c (size(ic)))
      mc(p1)% idx = 0
      mc    % n   = 0
      mc    % pe  = (/ (i,i=0,dace% npe-1) /)
    endif

    do l = ly(4), uy(4)
      do j = ly(2), uy(2)
        do i = ly(1), uy(1)
          ic (i,j,l)% h% w    = 0._wp
          ic (i,j,l)% h% imc  = 0
          ic (i,j,l)% h% ijdp = 0
          n  = idx(i,j,l)% n
          m  = sum (maxloc (idx(i,j,l)% w(:n)))
          if (n==0) cycle
          ic (i,j,l)% h% ijdp = idx(i,j,l)% ijdp (m,:)
          do k = 1, n
            ijdp = idx(i,j,l)% ijdp (k,:)
            ix = mc(p1)% idx (ijdp(1),ijdp(2),ijdp(3),1)
            if (ix == 0) then
              mc(p1)% n = mc(p1)% n + 1
              ix = mc(p1)% n
              mc(p1)% idx (ijdp(1),ijdp(2),ijdp(3),1) = ix
              if (size (mc(p1)% c) < ix) then
                tmp => mc(p1)% c
                allocate (mc(p1)% c (int(1.5 * ix)))
                mc(p1)% c (1:size(tmp)) = tmp
                deallocate (tmp)
              endif
              c => mc(p1)% c(ix)
              c% icol  = ix
              c% iatm  = iatm
              c% itrac = 0
              c% ijdtp = (/ijdp(:3),1,ijdp(4)/)
            else
              c => mc(p1)% c(ix)
              c% iatm  = ior (iatm,  c% iatm)
            endif
            ic (i,j,l)% h% w   (k)   = idx(i,j,l)% w (k)
            ic (i,j,l)% h% imc (k,1) = ix
          end do
        end do
      end do
    end do

  end subroutine mcol_idx
 !============================================================================
  subroutine get_varlist (ids, nid, idl)
    !---------------------------------------------------------
    ! get a list of identifiers for allocated multi level data
    !---------------------------------------------------------
    integer(i8),              intent(in)              :: ids    ! variables allocated (t_cols)
    integer,                  intent(out)             :: nid    ! number of allocated variables
    integer(i8), allocatable, intent(inout), optional :: idl(:) ! list of allocated variables

    integer(i8) :: idlist(N_COL)          ! local list of allocated variables
    integer     :: i

    nid = 0
    idlist = 0
    do i=0, N_COL-1
       if (btest(ids,i)) then
          nid = nid + 1
          idlist(nid) = ibset(idlist(nid),i)  ! should be equivalent to one of
       end if                                 ! the identifiers: COL_T, ...
    end do

    if (present(idl)) then
      if (allocated(idl)) deallocate(idl)
      allocate( idl(nid) )
      idl = idlist(1:nid)
    end if

  end subroutine get_varlist

!==============================================================================
! +++ The following routine belongs to a common module for low-level handling
!     of observations and their associated model columns (e.g. mo_obs_util) +++
!==============================================================================
  subroutine set_lev_insitu (spot, obs, cols)
    type(t_spot) ,intent(inout) :: spot     ! spot observations
    type(t_obs)  ,intent(inout) :: obs      ! observation data type
    type(t_cols) ,intent(in)    :: cols     ! model columns
    target :: cols
    !-----------------------------------------------------------
    ! interpolate pressure for in-situ-observations where height
    ! is the independent quantity (e.g. wind profiler); set plev
    !-----------------------------------------------------------
    real(wp)             :: olev            ! level of observation
    real(wp)             :: plev            ! pressure level
    real(wp)             :: zs              ! station height
    integer              :: ltyp            ! level type
    integer              :: i, j, k, n      ! indices
    type(t_col), pointer :: col
    logical              :: plev_only       ! all levels are plev?

    if (spot% o% i <  0      ) return
    if (spot% o% n == 0      ) return
    if (dace% pe   /= obs% pe) return

    plev_only = .true.
    do i = spot% o% i + 1, spot% o% i + spot% o% n
       select case (obs% body(i)% lev_typ)
       case (VN_P, -1)
       case default
          plev_only = .false.
          exit
       end select
    end do
    if (plev_only) return

    j = sum (maxloc (spot% col% h% w))
    k = spot% col% h% imc (j,1)

    col => cols% col(k)

    call check (col% geo, 'geo')
    call check (col% p,   'p')

    zs = spot% z
    if (zs == empty_spot% z) then
       ! Use model surface as proxy for (missing) station height
       zs = col% s% geosp / gacc
    end if

    n    = size (col% geo)
    olev = -huge(olev)
    plev = -1._wp

    do i = spot% o% i + 1, spot% o% i + spot% o% n
       if (olev /= obs% olev(i)) then
          plev = -1._wp
          olev = obs% olev(i)
          ltyp = obs% body(i)% lev_typ
          select case (ltyp)
          case (VN_HEIGHT, VN_FLEV)
             plev = pres_from_height (olev)
          case (VN_HOSAG)
             plev = pres_from_height (olev + zs)
          case (VN_P)
             plev = olev
          case default
             call finish ('set_lev_insitu',                     &
                          'lev_typ = '//name_value (varno, ltyp))
          end select

       endif

       select case (obs% body(i)% lev_typ)
       case (VN_HEIGHT, VN_FLEV, VN_HOSAG)
          obs% body(i)% plev = plev
       end select
    end do

    if (spot% z  /= empty_spot% z .and. &
        spot% ps == empty_spot% ps      ) spot% ps = pres_from_height (spot% z)

  contains

    subroutine check (ptr, name)
      !-------------------------------------------
      ! check allocation status of field in column
      !-------------------------------------------
      real(wp)         ,pointer    :: ptr(:)
      character(len=*) ,intent(in) :: name

      if (.not. (associated (ptr))) then
         write(0,*) "set_lev_insitu: statid, obstype, codetype: ", &
              spot% statid, spot% hd% obstype, spot% hd% codetype
         write(0,*) "set_lev_insitu: lev_typ =",             &
              obs% body(spot%o%i+1:spot%o%i+spot%o%n)% lev_typ
         call finish ('set_lev_insitu','not associated: '//trim (name))
      end if
    end subroutine check

    function pres_from_height (z) result (p)
      real(wp), intent(in) :: z
      real(wp)             :: p

      real(wp) :: gp    ! geopotential
      real(wp) :: w     ! interpolation weight
      integer  :: l     ! loop index

      gp = z * gacc
      if      (gp <= col% geo(n)) then
         p = exp (col% p(n))
      else if (gp >= col% geo(1)) then
         p = exp (col% p(1))
      else
         do l = 2, n
            if (gp > col% geo(l)) then
               w = (gp            - col% geo(l)) / &
                   (col% geo(l-1) - col% geo(l))
               p = exp (       w  * col% p(l-1) + &
                        (1._wp-w) * col% p(l  )   )
               exit
            endif
         end do
      end if
    end function pres_from_height

  end subroutine set_lev_insitu
!==============================================================================
end module mo_t_col
