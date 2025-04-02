!
!+ Radiance bias correction utilities
!
MODULE mo_radbiascor
!
! Description:
!   Radiance bias correction utilities
!
! Current Code Owner: DWD, Andreas Rhodin
!    phone: +49 69 8062 2722
!    fax:   +49 69 8062 3721
!    email: andreas.rhodin@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_8         2009/12/09 Mashrab Kuvatov
!  Radiance bias correction utilities
! V1_9         2010/04/20 Mashrab Kuvatov
!  new predictors, regularisation
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_11        2010/06/16 Andreas Rhodin
!  remove unused single-profile routines
! V1_13        2011/11/01 Mashrab Kuvatov
!  changes for IASI, HIRS, AMSUB
! V1_14        2011/11/08 Andreas Rhodin
!  add CDIR NOUNROLL directive (work around bug in sxf90/rev430)
! V1_15        2011/12/06 Harald Anlauf
!  workaround for sxf90 outerunroll bug; remove previous compiler directive
! V1_19        2012-04-16 Harald Anlauf
!  mo_radbiascor: read bias correction coefficients only on I/O pe
! V1_22        2013-02-13 Robin Faulwetter
!  Improved selection of FOVs to be bias corrected, and monitored
! V1_23        2013-03-26 Robin Faulwetter
!  Improve bias correction modules: t_decay = 0. : infinite memory
!                                   t_decay < 0. : static bias correction.
! V1_27        2013-11-08 Andreas Rhodin
!  implementation of Variational Bias Correction (VarBC)
! V1_28        2014/02/26 Andreas Rhodin
!  bias cor. coefficient calculation: saveguard for <=1 entries in statistics
!  set m_pred (max. number of predictors in BC file) to 16
! V1_29        2014/04/02 Robin Faulwetter
!  Fixed bug with two FOVs in the calculation of bias correction coefficients
! V1_31        2014-08-21 Robin Faulwetter
!  Unify mo_rad with COSMO. Improve mo_rttov_ifc. New write_rttov_prof
! V1_35        2014-11-07 Robin Faulwetter
!  New options for FOV bias correction for radiances
! V1_37        2014-12-23 Robin Faulwetter
!  Cleanup
! V1_47        2016-06-06 Harald Anlauf
!  force large namelist arrays on heap (reduces stack pressure with OpenMP)
! V1_48        2016-10-06 Robin Faulwetter
!  Allow different t_decay values for different channels
! V1_50        2017-01-09 Robin Faulwetter
!  Added task for modifying biascor files to rad_biascor.
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Mashrab Kuvatov   DWD            2009
! Marc Schwaerz     DWD/EUMETSAT,  2009
!==============================================================================

!-------------
! Modules used
!-------------
  use environment,  only: model_abort,    &! program abort in case of error
                          get_free_unit    ! get a free fortran unit
  use mo_kind,      only: wp               ! working precision kind parameter
  use mo_namelist,  only: position_nml,   &! routine to position nml group
                          POSITIONED       ! position_nml: OK return flag
  use mo_matrix,    only: inverse_rs       ! invert real symmetric matrix
  use mo_mpi_dace,  only: dace,           &! MPI group info
                          p_bcast          ! broadcast routine
  use mo_rad,       only: t_radv,         &! Type for radiance data
                          t_rad_set,      &
                          construct,      &
                          destruct,       &
                          assignment (=), &
                          m_instr,        &
                          m_chan,         &
                          chan_indx,      &
                          set_indx,       &
                          instr_chan,     &
                          rinvalid,       &
                          f_ins,          &
                          rad_set, n_set, &
                          USE_QCNOTUSE
  use mo_physics,   only: d2r,            &! factor degree -> radians
                          pi
  use mo_algorithms,only: fit_polynomial   ! fit a polynomial to given data
  use mo_t_obs,     only: usd
  use mo_dace_string,only:split2,         &!
                          tolower
  use mo_run_params,only: ana_time,       &!
                          fc_time
  use mo_time,      only: cyyyymmddhhmm,       &! convert time to string
                          time_cyyyymmddhhmm,  &! convert string to time
                          operator(-),         &! subtract times
                          days                  ! convert time to days
!--------------------------------
! modules to be used within COSMO
!--------------------------------

  use mo_physics,  only: tv_t_q,       &! calculate virtual temperature
                         gacc,         &! gravity acceleration
                         R_DRY => R,   &! Dry air gas constant  [J/(kg*K)]
                         EPS   => RdRd  ! R_WV/R_DRY   (R_WV: R_waterVapour)

  implicit none

!=======================================================================

  !================
  ! Public entities
  !================
  private

  public :: bccoef                 ! bias correction statistics & coefficients
  !--------------------
  ! constant parameters
  !--------------------
  public :: m_pred                 ! max.number of predictors used for bc.
  public :: P_OBS_BT               ! predictor: observed brightness temp.
  public :: P_BT                   ! predictor: model brightness temperature
  public :: P_THICK                ! predictor: model layer thickness
  public :: P_SPHER                ! predictor: spherical harmonic
  public :: P_PROD                 ! predictor: product of two predictors
  public :: P_SCANANG              ! predictor: scan angle
  public :: P_THICK_TR             ! predictor: model layer thickness weighted with
                                   !            transmission decrease within layer
  public :: P_THICK_TRF            ! predictor: model layer thickness weighted with
                                   !            top-of-layer transmission
  public :: P_TRANS_T              ! predictor: T integrated over transmission interval
  public :: P_TRANS_TS             ! predictor: T integrated over transmission interval
                                   !            but not below surface
  public :: P_TRANS_P              ! predictor: pressure of transmission value
  public :: P_TRANS_DP             ! predictor: pressure difference between transmission values
  public :: P_PLEV                 ! predictor: plevel of observation
  public :: P_INSTR_T              ! predictor: instrument temperature
  public :: P_ORB_PH               ! predictor: orbital phase angle
  public :: chan_dep_pred          ! array of channel dependent predictors
  !--------------------------
  ! parameters to be adjusted
  !--------------------------
  public :: min_entries            ! min. number of entries to keep
  public :: atm_thick_vers
  public :: offset_max_age         ! max. age allowed in offset calculations
  public :: offset_t_decay         ! t_decay applied to old offset, if offset can't be calculated
  !-------------------------
  ! derived type definitions
  !-------------------------
  public :: t_bcor_coef            ! bias correction coefficient data type
  public :: t_stats_data           ! bias correction statistics  data type
  !------------------------------------
  ! constructor and destructor routines
  !------------------------------------
  public :: construct              ! construct: t_stats_data
  public :: destruct               !  destruct: t_stats_data
  public :: assignment (=)         !    assign: t_rad_set
  public :: alloc_bcor_coef        ! allocate components in t_bcor_coef
  !-------------
  ! I/O routines
  !-------------
  public :: read_bcor_file         ! read  file with bc coeffs.+statistics
  public :: write_bcor_file        ! write file with bc coeffs.+statistics
  !-------------------------
  ! bias correction routines
  !-------------------------
  public :: calc_bcor_stats        ! bias corr. statistics (single channel)
  public :: calc_bcor_t_stats      ! bias corr. statistics (type t_rad)
  public :: merge_bcor_stats       ! merge statistics (online bias corr.)
  public :: scale_bcor_stats       ! rescale statistics
  public :: calc_bcor_coefs        ! bias corr. coefficients from statistics
  public :: calc_bcor_coefs_ac     ! coefficients from accumulated statistics
  public :: apply_bias_corr        ! apply bias correction
  !-------------------
  ! auxiliary routines
  !-------------------
  public :: chan_indx              ! get channel index in t_bcor_coef
  public :: set_indx               ! get bias correction set index
  public :: calc_predictors        ! calculate predictors given by t_bcor_coef
  public :: diagn_biasc            ! print bias coefficient diagnostic files
  public :: diagn_pred             ! determine most important predictors

!=======================================================================

!-----------------------------------------------------------------------
!--- error numbers of this module:

!> error code: everything was ok.
 integer, parameter :: &
                       NO_ERROR             = 0

!> error code: an allocation error occurred
 integer, parameter :: &
                       ERR_ALLOC            = 1

!> error code: missmatch in dimensions occurred
 integer, parameter :: &
                       ERR_DIM              = 2

!> error code: grid is not monotonous
 integer, parameter :: &
                       ERR_MONOTON          = 3

!-----------------------------------------------------------------------
!--- error check switch:

!> set to .true. causes extensive dimension and error checking:
logical, parameter :: checkInput = .true.

  !====================
  ! constant parameters
  !====================

  !--------------------------------------------------------
  ! max. size of derived type components within this module
  !--------------------------------------------------------
  integer, parameter :: m_pred  =  30 ! max. number of predictors
  integer, parameter :: m_pp    = m_pred*(m_pred+1)/2
  integer, parameter :: m_fov   = 100 ! max. number of fov
  !-----------
  ! predictors
  !-----------
  character(len=10), parameter :: P_THICK   = 'thickness' ! thickness
  character(len=10), parameter :: P_IWV     = 'iwv'       ! integrated water v.
  character(len=10), parameter :: P_BT      = 'bt'        ! brightness temper.
  character(len=10), parameter :: P_OBS_BT  = 'obs_bt'    ! observed br.temp.
  character(len=10), parameter :: P_TSKIN   = 'T_SurfSkin'! surface skin temp.
  character(len=10), parameter :: P_SCANANG = 'scan.ang'  ! scan angle
  character(len=10), parameter :: P_V10     = 'v10'       ! 10 m wind speed.
  character(len=10), parameter :: P_SPHER   = 'Y_lm'      ! spherical harmonics
  character(len=10), parameter :: P_PROD    = 'product'   ! predictor product
  character(len=10), parameter :: P_THICK_TR= 'thick_tr'  ! thickness weighted with transmission
  character(len=10), parameter :: P_THICK_TRF='thick_trf' ! thickness weighted with transmission
  character(len=10), parameter :: P_TRANS_T = 'trans_t'   ! temp. integrated over transmission
  character(len=10), parameter :: P_TRANS_TS= 'trans_ts'  ! temp. integrated over transmission
  character(len=10), parameter :: P_TRANS_P = 'trans_p'   ! pressure transmission value
  character(len=10), parameter :: P_TRANS_DP= 'trans_dp'  ! pressure between transmission value
  !TODO:
  character(len=10), parameter :: P_PLEV    = 'plev'      ! pressure transmission value
  character(len=10), parameter :: P_INSTR_T = 'instr_t'   ! instrument temperature
  character(len=10), parameter :: P_ORB_PH  = 'orbit_ph'  ! pressure between transmission value
  character(len=10), parameter :: chan_dep_pred(7) = (/P_THICK_TR,P_THICK_TRF,P_TRANS_T ,&
                                                       P_TRANS_TS,P_TRANS_P  ,P_TRANS_DP,&
                                                       P_PLEV     /)

  integer :: atm_thick_vers = 0
  ! 0: old version (constant Tv in extrapolation below ground)
  ! 1: emulates old version with new code
  ! 2: assumes the climatological T-gradient in extrap. below ground
  ! 3: Similar to 2, but the boundary layer is ignored. This makes sense, because the boundary layer
  !    temperatures are dominated by the daily-cycle and not by the synoptic situation in the the
  !    free atmosphere.
  !    The rationale behind this in connection to the Radiance-bias-correction should be explained:
  !    The original assumption behind the layer-depth-predictors is, that the main source of bias
  !    is the radiative transfer (RTTOV). With this assumption, option 3 does not make any sense.
  !    However, there is plenty of evidence, that the original explanation for the layer-thickness
  !    predictors is not appropriate. In my option (RF) they somehow represent the large-scale
  !    synoptic situation and correct mainly model biases. With this in mind, option 3 makes sense.

  !=====================
  ! parameters to adjust
  !=====================
  real(wp) :: min_entries = 100._wp ! min. number of entries in statistics
  real(wp) :: offset_max_age = 30._wp  ! maximum allowed age of statistics in offset calc. [days]
  real(wp) :: offset_t_decay = 30._wp  ! t_decay of offset [days], if offset can't be calculated

  !=========================
  ! derived type definitions
  !=========================


  !----------------------------------------------
  ! Derived type to hold data for statistics file
  !----------------------------------------------
  type t_stats_data
     ! FIXME:
     ! classify for:
     ! land / sea /coast
     ! clear / cloud (total overcast) / fractional cloud
     ! day (solar zenith < 83) / night (>90) / twilight
     ! longitude band
     ! latitude band
     logical           ::  alloc            = .false. ! flag for components allocated
     real(wp) ,pointer ::  n         (:)    => NULL() ! Number of entries / channel
     real(wp) ,pointer ::  btdscan   (:,:)  => NULL() ! Sum of bt difference
                                                      ! Size of FOV x channel
     real(wp) ,pointer ::  btdn      (:,:)  => NULL() ! Number of data
     real(wp) ,pointer ::  d2        (:,:)  => NULL() ! deviations squared
     real(wp) ,pointer ::  pred_i    (:,:,:)=> NULL() ! Sum of i-th predictor
                                                      ! Size: pred. x FOV x channel
     real(wp) ,pointer ::  btd_pred_i(:,:,:)=> NULL() ! Sum of bt difference * pred
                                                      ! Size: pred. x FOV x channel
     real(wp) ,pointer ::  pred_ij   (:,:,:)=> NULL() ! Sum of predictor products
                                                      ! example of 3 predictors
                                                      ! p11,p12,p13,p22,p23,p33
                                                      ! Size: pred. x FOV x channel
     integer           :: n_pred = -1                 ! number of predictors
  end type t_stats_data


  !-------------------------------------------------
  ! Parameters for types of fov bias correction
  !-------------------------------------------------
  integer, parameter :: fovbc_type_none    = 0  ! No FOV biascor
  integer, parameter :: fovbc_type_fov     = 1  ! independent FOV biascor for each individual FOV
  integer, parameter :: fovbc_type_poly    = 2  ! FOV biascor is a polynomial in the FOV number
  integer, parameter :: fovbc_const_none   = 0  ! Do not change mean FOV biascor
  integer, parameter :: fovbc_const_allobs = 1  ! FOV biascor = 0 averaged over all obs.
  integer, parameter :: fovbc_const_scanl  = 2  ! FOV biascor = 0 averaged over one scanline
  integer, parameter :: fovbc_const_nadir  = 3  ! FOV biascor = 0 at nadir
  !------------------------
  ! Parameters for mode_wgt
  !------------------------
  integer, parameter :: WGT_ERR            = 0
  integer, parameter :: WGT_LAT            = 1
  integer, parameter :: WGT_SZA            = 2
  !-------------------------------------------------
  ! Derived type to hold info on fov bias correction
  !-------------------------------------------------
  type t_fov_bc
     integer :: type       = fovbc_type_fov     ! type of FOV bcor (see parameters below)
     integer :: n          = 0                  ! degree of polynomial if type=fovbc_type_poly
     integer :: const_fov  = fovbc_const_allobs ! how to calulate constant in fov bcor (see parameters below)
     logical :: const_pred = .true.             ! use constant resulting from predictor "1"
  end type t_fov_bc
  type(t_fov_bc), save :: fov_bc_default

  !--------------------------------------------------
  ! Derived type to hold bias-correction coefficients
  !--------------------------------------------------
  type t_bcor_coef
     !----------
     ! meta data
     !----------
     integer                 :: exp      = -1    ! experiment number
     character(len=8)        :: model            ! model used
     type (t_rad_set)        :: i                ! set of instruments & channels used
     integer                 :: n_pred           ! number of predictors
     integer                 :: n_coeff          ! total number of coefficients
     character(len=12)       :: created          ! creation date 'yyyymmddhhmm'
     character(len=12)       :: modified         ! modification date
     character(len=12)       :: first_date       ! first training date
     character(len=12)       :: last_date        ! last  change date
     character(len=12)       :: last_update      ! last  training date
     real(wp)                :: entries          ! number of entries for training
     real(wp)                :: t_decay(m_chan)  ! background info decay time (days)
     integer                 :: mdlsfc_stat      ! surface mask used for statistics
     real(wp)                :: nu       = 0._wp ! regularisation parameter
     real(wp)                :: w_vbc    = 1._wp ! variational bias correction weight
     logical                 :: fov_bcorr=.true. ! bias correction for individual fov
     integer                 :: mode_wgt = 0     ! Mode for weighting of observations
                                                 ! (see WGT_* parameters)
                                                 !      0: no weighting (w=1.)
                                                 ! bit0=1: w=1./(e_fg^2+e_obs^2)
                                                 ! bit1=2: w=(1.+sin(lat)^2) i.e. counteract
                                                 !         the effect of thinning two orbits
                                                 !         (per 3hour ass. window) near the poles
                                                 ! bit2=4: w=1./cos(sza)?? i.e. resemble the
                                                 !         effect of thinning in case we have
                                                 !         unthinned data in biascor.
     integer                 :: constr_mode(m_chan)  = 0     ! mode for constrained bias correction,
                                                             ! 1: minimize <(d-b)^2 + alpha*(b-b0)^2>
                                                             ! 2: minimize <(d-b)^2> + alpha*(<b>-b0)^2
                                                             ! 3: minimize <(d-b)^2> under the constraint <b>=b0
                                                             ! <...> temporal mean, d=o-f, b=bcor
     real(wp)                :: constr_b0(m_chan)    = 0._wp ! "b0" for constrained bias correction,
                                                             ! "Target bias correction"
     real(wp)                :: constr_alpha(m_chan) = 0._wp ! "alpha" for constrained bias correction
                                                             ! "strength of constraint"
     character(len=300)      :: offset  (m_chan)  ! string describing offset
     real(wp)                :: offset_v(m_chan)  ! offset value (calculated in last analysis)
     character(len=10)       :: pred   (m_pred)  ! predictors
     integer                 :: pred_p1(m_pred)  ! predictor parameter 1
     integer                 :: pred_p2(m_pred)  ! predictor parameter 2
     logical                 :: chan_dep_pred  = .false.! whether the predictors depend on channels

     integer,        pointer :: pred_use (:,:) =>NULL() ! 1 for used predictors
     type(t_fov_bc), pointer :: fov_bc   (:)   =>NULL() !

     !-----------------------------
     ! bias correction coefficients
     !-----------------------------
     real(wp),       pointer :: coef_set1 (:,:) =>NULL() ! bias correction set 1
     real(wp),       pointer :: coef_set2 (:,:) =>NULL() ! bias correction set 2
     !----------------------------------------------------------
     ! bias coefficient background error covariance matrix
     ! (active predictors only, for Variational Bias Correction)
     !----------------------------------------------------------
     real(wp),       pointer :: n         (:)   =>NULL() ! number of entries
     real(wp),       pointer :: p_mean    (:,:) =>NULL() ! mean predictor value
     real(wp),       pointer :: p_scal    (:,:) =>NULL() ! normalisation
     real(wp),       pointer :: p_b       (:)   =>NULL() ! predictor B-matrix
     integer,        pointer :: i_p       (:,:) =>NULL() ! indices of active predict.
     integer,        pointer :: n_b       (:,:) =>NULL() ! :1) # of active predictors
                                                   ! :2) # squared (size of B)
                                                   ! :3) #**2 accumulated
     integer                 :: nb                       ! total size of B matrix
     integer                 :: mp              = 0      ! max # of active predict.
     !----------------
     ! bias statistics
     !----------------
     type(t_stats_data)      :: sti                      ! instantaneous statistics
     type(t_stats_data)      :: sta                      ! accumulated   statistics
!+   !--------------------
!+   ! 3D-Var internal use
!+   !--------------------
!+   integer                 :: n_ctr                    ! number of active coeffs.
!+   integer                 :: i_ctr                    ! controlvariable index
  end type t_bcor_coef

  type(t_bcor_coef) ,pointer :: bccoef  (:) => NULL()

  !-----------
  ! interfaces
  !-----------

  interface construct
    module procedure construct_stats_data
    module procedure construct_fov_bc
  end interface construct

  interface destruct
    module procedure destruct_stats_data
    module procedure destruct_bcor_coef
  end interface

  interface calc_predictors
    module procedure calc_predictors
    module procedure calc_predictors_rad
    module procedure calc_predictors_rads
  end interface

  interface apply_bias_corr
    module procedure apply_bias_corr
    module procedure apply_bias_corr_rads
  end interface apply_bias_corr

  interface calc_bcor_t_stats
    module procedure calc_bcor_t_stats
    module procedure calc_bcor_t_stat
  end interface calc_bcor_t_stats

  interface p_bcast
     module procedure p_bcast_fov_bc
  end interface

  interface operator(==)
    module procedure equal_fov_bc
  end interface


!==============================================================================
contains
!==============================================================================



  subroutine construct_stats_data (bcor_stats, n_pred, n_fov, n_chan)
  !====================================================
  ! Allocate memory for bias correction statistics data
  !====================================================
  type(t_stats_data) ,intent(inout) ::  bcor_stats ! Statistics data
  integer            ,intent(in)    ::  n_pred     ! Number of predictors
  integer            ,intent(in)    ::  n_fov      ! Number of FOVs
  integer            ,intent(in)    ::  n_chan     ! Number of channels

    bcor_stats% alloc      = .true.
    bcor_stats% n_pred     = n_pred
    allocate (bcor_stats% n          (                            n_chan))
    allocate (bcor_stats% btdscan    (                     n_fov, n_chan))
    allocate (bcor_stats% btdn       (                     n_fov, n_chan))
    allocate (bcor_stats% d2         (                     n_fov, n_chan))
    allocate (bcor_stats% pred_i     (n_pred,              n_fov, n_chan))
    allocate (bcor_stats% btd_pred_i (n_pred,              n_fov, n_chan))
    allocate (bcor_stats% pred_ij    (n_pred*(n_pred+1)/2, n_fov, n_chan))
    bcor_stats% n          = 0._wp
    bcor_stats% btdscan    = 0._wp
    bcor_stats% btdn       = 0._wp
    bcor_stats% d2         = 0._wp
    bcor_stats% pred_i     = 0._wp
    bcor_stats% btd_pred_i = 0._wp
    bcor_stats% pred_ij    = 0._wp
  end subroutine construct_stats_data

!------------------------------------------------------------------------------
!#define DEALLOC(var) if (associated(var)) deallocate(var, stat=stat)
#define DEALLOC(var) if (associated(var)) deallocate(var)

  elemental subroutine destruct_stats_data (bcor_stats)
  !======================================================
  ! Deallocate memory for bias correction statistics data
  !======================================================
  type(t_stats_data) ,intent(inout) ::  bcor_stats ! Statistics data

    bcor_stats% alloc   = .false.
    bcor_stats% n_pred  = -1
    DEALLOC(bcor_stats% n         )
    DEALLOC(bcor_stats% d2        )
    DEALLOC(bcor_stats% btdscan   )
    DEALLOC(bcor_stats% btdn      )
    DEALLOC(bcor_stats% pred_i    )
    DEALLOC(bcor_stats% btd_pred_i)
    DEALLOC(bcor_stats% pred_ij   )

  end subroutine destruct_stats_data

!------------------------------------------------------------------------------

  elemental subroutine destruct_bcor_coef (bcc)
  !======================================================
  ! Deallocate memory for bias correction statistics data
  !======================================================
    type(t_bcor_coef) ,intent(inout) ::  bcc ! Statistics data

    DEALLOC(bcc%pred_use)
    DEALLOC(bcc%fov_bc)
    DEALLOC(bcc%coef_set1)
    DEALLOC(bcc%coef_set2)
    DEALLOC(bcc%n)
    DEALLOC(bcc%p_mean)
    DEALLOC(bcc%p_scal)
    DEALLOC(bcc%p_b)
    DEALLOC(bcc%i_p)
    DEALLOC(bcc%n_b)
    call destruct_stats_data(bcc%sta)
    call destruct_stats_data(bcc%sti)

  end subroutine destruct_bcor_coef

!------------------------------------------------------------------------------

  elemental subroutine construct_fov_bc(x)
     type(t_fov_bc), intent(out) :: x
  end subroutine construct_fov_bc

!------------------------------------------------------------------------------

  subroutine alloc_bcor_coef (bcor)
  type (t_bcor_coef) ,intent(inout) :: bcor
  !----------------------------------------------------
  ! allocate components within derived type t_bcor_coef
  !----------------------------------------------------

    allocate (bcor% coef_set1 (bcor%    n_pred + 1 ,bcor% i% n_chan ))
    allocate (bcor% coef_set2 (bcor% i% n_chan     ,bcor% i% n_fov  ))

    bcor% coef_set1 = 0._wp
    bcor% coef_set2 = 0._wp

  end subroutine alloc_bcor_coef

!==============================================================================

  subroutine calc_bcor_coefs_ac (bcor_coef, mask_rs)
    type(t_bcor_coef) ,intent(inout)        :: bcor_coef(:) ! coefficients + statistics
    logical,           intent(in), optional :: mask_rs(:)
    !===================================================================
    ! Calculate bias correction coefficients from accumulated statistics
    !===================================================================
    integer :: i
    do i = 1, size (bcor_coef)
      if (present(mask_rs)) then
        if (.not.mask_rs(i)) CYCLE
      end if
      call calc_bcor_coefs (bcor_coef(i)% sta, bcor_coef(i))
    end do
  end subroutine calc_bcor_coefs_ac

!------------------------------------------------------------------------------

  subroutine calc_bcor_coefs (stats_data, bcor_coef)
  !=======================================
  ! Calculate bias correction coefficients
  !=======================================
  type(t_stats_data)    ,intent(in)    ::  stats_data ! Statistics data
  type(t_bcor_coef)     ,intent(inout) ::  bcor_coef  ! Bias cor. coefficients

    integer                 ::  n_chan         ! Number of channels
    integer                 ::  n_fov          ! Number of FOVs
    integer                 ::  n_pred         ! Number of predictors
    integer                 ::  ch_counter     ! Counter for channels
    integer                 ::  fov_counter    ! Counter for FOVs
    integer                 ::  pred_counter   ! Counter for predictors
    integer                 ::  pred_counter_i ! Counter for predictors
    integer                 ::  pred_counter_j ! Counter for predictors
    integer                 ::  i,j,n          ! matrix indices
    integer                 ::  stat
    real(wp)                ::  c, cn
    real(wp)  ,allocatable  ::  x(:)
    real(wp)  ,allocatable  ::  y(:)
    real(wp)  ,allocatable  ::  w(:)
    real(wp)  ,allocatable  ::  p(:)
    logical   ,allocatable  ::  l_fovbc_set(:)
    !----------------------------------------------
    ! Same as in t_stats_data, but sumed over FOVs
    !----------------------------------------------
    real(wp)  ,allocatable  ::  btdscan_ch    (:)   ! Tb difference
    real(wp)  ,allocatable  ::  btdscan_pi_ch (:,:) ! Tb diff. x predictor
    real(wp)  ,allocatable  ::  btdn_ch       (:)   ! number of entries
    real(wp)  ,allocatable  ::  pi_ch         (:,:) ! predictors
    real(wp)  ,allocatable  ::  pij_ch        (:,:) ! predictor products
    !----------------------------------------------
    real(wp)  ,allocatable  ::  sfov_pi       (:,:) ! rhs correction
    real(wp)  ,allocatable  ::  matrix_a      (:,:) ! matrix of pred. x pred.
    real(wp)  ,allocatable  ::  matrix_inv    (:,:) ! inverse matrix (B)
    real(wp)  ,allocatable  ::  vector_b      (:)   ! rhs
    real(wp)  ,allocatable  ::  tmp_b         (:)   ! lhs scaled coefficients
    real(wp)  ,allocatable  ::  scal          (:)   ! scaling factors
    !----------------------------------------------
    real(wp)                ::  offset              ! "offset" bias correction
    !-----------------------------------------
    ! if statistics file is empty, just bypass
    !-----------------------------------------
    bcor_coef% entries = 0._wp
    bcor_coef% mp      = 0
    if (.not. stats_data% alloc) return
    !------------------
    ! derive dimensions
    !-------------------
    n_chan = bcor_coef% i% n_chan
    n_fov  = bcor_coef% i% n_fov
    n_pred = stats_data%   n_pred
!   if (n_pred < 1) then
!      print *, 'Error: Wrong number of predictors.'
!      return
!   end if
    !--------------
    ! set meta data
    !--------------
    bcor_coef% entries = maxval (stats_data% n)
    !----------------------------------------------
    ! Allocate memory for coefficients and B matrix
    !----------------------------------------------
    if (associated (bcor_coef% coef_set1)) deallocate (bcor_coef% coef_set1)
    if (associated (bcor_coef% coef_set2)) deallocate (bcor_coef% coef_set2)
    if (associated (bcor_coef% p_mean   )) deallocate (bcor_coef% p_mean   )
    if (associated (bcor_coef% p_scal   )) deallocate (bcor_coef% p_scal   )
    if (associated (bcor_coef% p_b      )) deallocate (bcor_coef% p_b      )
    if (associated (bcor_coef% n_b      )) deallocate (bcor_coef% n_b      )
    if (associated (bcor_coef% i_p      )) deallocate (bcor_coef% i_p      )
    if (associated (bcor_coef% n        )) deallocate (bcor_coef% n        )

    allocate (bcor_coef% coef_set1 (n_pred + 1 ,n_chan))
    allocate (bcor_coef% coef_set2 (n_chan     ,n_fov ))
    allocate (bcor_coef% n_b       (n_chan     ,3     ))
    allocate (bcor_coef% n         (n_chan))

    bcor_coef% coef_set1 = 0._wp
    bcor_coef% coef_set2 = 0._wp
    bcor_coef% n         = 0._wp

    !----------------------------------------
    ! Allocate memory for auxiliary variables
    !----------------------------------------
    allocate (btdscan_ch    (n_chan))
    allocate (btdscan_pi_ch (n_pred, n_chan))
    allocate (btdn_ch       (n_chan))
    allocate (pi_ch         (n_pred, n_chan))
    allocate (pij_ch        (n_pred * (n_pred + 1)/2, n_chan))

    !----------------------------
    ! Sum up statistics over FOVs
    !----------------------------
    bcor_coef% n_b (1,3) = 0
    do i = 1, n_chan
       !----------
       ! mean bias
       !----------
       btdscan_ch(i)    =  sum(stats_data% btdscan(:, i))
       btdn_ch   (i)    =  sum(stats_data% btdn   (:, i))
       !----------------------------------
       ! predictors and predictor products
       !----------------------------------
       do pred_counter = 1, n_pred
          btdscan_pi_ch                (pred_counter,    i) =  &
            sum (stats_data% btd_pred_i(pred_counter, :, i))
          pi_ch                        (pred_counter,    i) =   &
            sum (stats_data% pred_i    (pred_counter, :, i))
       end do
       do pred_counter = 1, n_pred * (n_pred + 1)/2
          pij_ch                    (pred_counter,    i) =  &
            sum (stats_data% pred_ij(pred_counter, :, i))
       end do
       !-------------------------
       ! set indices for B-matrix
       !-------------------------
       bcor_coef%          n_b (i,1) = 0
       if (btdn_ch (i) > 0) &
         bcor_coef%        n_b (i,1) = count (bcor_coef% pred_use (:, i) == 1)
       bcor_coef%          n_b (i,2) = bcor_coef% n_b (i,1) ** 2
       if (i>1) bcor_coef% n_b (i,3) = bcor_coef% n_b (i-1,2) &
                                     + bcor_coef% n_b (i-1,3)
    end do

    !-----------------------------
    ! set up zero B matrix for VBC
    !-----------------------------
    bcor_coef% nb = sum         (bcor_coef% n_b (:,2))
    bcor_coef% mp = maxval      (bcor_coef% n_b (:,1))
    allocate (bcor_coef% p_b    (bcor_coef% nb        ))
    allocate (bcor_coef% i_p    (bcor_coef% mp, n_chan))
    allocate (bcor_coef% p_mean (bcor_coef% mp, n_chan))
    allocate (bcor_coef% p_scal (bcor_coef% mp, n_chan))
    bcor_coef% p_b = 0._wp
    bcor_coef% i_p = 0

    !-------------------------------------------------
    ! Calculate scan dependent bias and rhs correction
    !-------------------------------------------------
    allocate(sfov_pi(n_pred + 1, n_chan))
    sfov_pi              = 0._wp
    bcor_coef% coef_set2 = 0._wp
    if (bcor_coef% fov_bcorr) then
      !--------------------
      ! FOV bias correction
      !--------------------
      allocate(l_fovbc_set(n_fov))
      do ch_counter = 1, n_chan
        l_fovbc_set = .false.
        c = 0._wp
        select case(bcor_coef% fov_bc(ch_counter)% type)
        case(fovbc_type_none)
          CYCLE
        case(fovbc_type_fov )
          do fov_counter = 1, n_fov
            if (stats_data% btdn (fov_counter, ch_counter) > 0) then
              bcor_coef% coef_set2     (ch_counter, fov_counter) =   &
                   stats_data% btdscan  (fov_counter, ch_counter) /   &
                   stats_data% btdn     (fov_counter, ch_counter)
              l_fovbc_set(fov_counter) = .true.
            end if
          end do
        case(fovbc_type_poly)
          n = count(stats_data% btdn(:, ch_counter) > 0)
          if (n > 0 .and. bcor_coef% fov_bc(ch_counter)% n > 0) then
            allocate(x(n_fov),y(n),w(n),p(bcor_coef% fov_bc(ch_counter)% n+1))
            i = 0
            do j = 1, n_fov
              if (stats_data% btdn(j, ch_counter) > 0) then
                i = i + 1
                y(i) = stats_data% btdscan(j,ch_counter) / stats_data% btdn(j,ch_counter)
                x(i) = j-0.5*(n_fov+1)
                w(i) = stats_data% btdn   (j,ch_counter)
              end if
            end do
            call fit_polynomial(x(1:n), y, bcor_coef% fov_bc(ch_counter)% n, stat, &
                 wgts=w, p=p)
            x(:) = (/ (j-0.5*(n_fov+1),j=1,n_fov) /)
            bcor_coef% coef_set2(ch_counter, :) = p(1)
            do j = 1, bcor_coef% fov_bc(ch_counter)% n
              bcor_coef% coef_set2(ch_counter, :) = bcor_coef% coef_set2(ch_counter, :) + &
                   p(j+1) * x(:)**j
            end do
            c = p(1)
            deallocate(x,y,w,p)
            l_fovbc_set(:) = .true.
          else
            c = 0._wp
          end if
        end select

        !-------------------------------------
        ! constant part of FOV bias correction
        !-------------------------------------
        select case(bcor_coef% fov_bc(ch_counter)% const_fov)
        case(fovbc_const_none)
          ! Nothing to do
        case(fovbc_const_allobs)
          cn = 0._wp
          c  = 0._wp
          do fov_counter = 1, n_fov
            if (l_fovbc_set(fov_counter)) then
              c  = c  + stats_data% btdn (fov_counter, ch_counter) * &
                   bcor_coef% coef_set2(ch_counter, fov_counter)
              cn = cn + stats_data% btdn (fov_counter, ch_counter)
            end if
          end do
          if (cn > 0._wp) c = c / cn
        case(fovbc_const_scanl)
          n = 0
          c = 0._wp
          do fov_counter = 1, n_fov
            if (l_fovbc_set(fov_counter)) then
              c = c + bcor_coef% coef_set2(ch_counter, fov_counter)
              n = n + 1
            end if
          end do
          if (n > 0) c = c / (1._wp * n)
        case(fovbc_const_nadir)
          if (bcor_coef% fov_bc(ch_counter)% type /= fovbc_type_poly) then
            c = 0._wp
            if (mod(n_fov,2) == 1) then
              if (l_fovbc_set(n_fov/2+1)) then
                c = bcor_coef% coef_set2(ch_counter, n_fov/2 + 1)
              else
                j = 1
                do while (j <= n_fov/2)
                  if (l_fovbc_set(n_fov/2+1+j) .and. l_fovbc_set(n_fov/2+1-j)) then
                    c = 0.5_wp * (bcor_coef% coef_set2(ch_counter, n_fov/2+1+j) + &
                                  bcor_coef% coef_set2(ch_counter, n_fov/2+1-j))
                    exit
                  end if
                end do
              end if
            else
              j = 0
              do while (j < n_fov/2)
                if (l_fovbc_set(n_fov/2-j) .and. l_fovbc_set(n_fov/2+1+j)) then
                  c = 0.5_wp * (bcor_coef% coef_set2(ch_counter, n_fov/2-j) + &
                                bcor_coef% coef_set2(ch_counter, n_fov/2+1+j))
                  exit
                end if
              end do
            end if
          else
            ! use constant of polynomial (see above)
          end if
        end select
        where (l_fovbc_set(:)) &
             bcor_coef% coef_set2(ch_counter,:) = bcor_coef% coef_set2(ch_counter,:) - c
      end do
      deallocate(l_fovbc_set)

      !------------------------------------------------------------------
      ! rhs correction term accounting for scan dependent bias correction
      !------------------------------------------------------------------
      do ch_counter = 1, n_chan
        do fov_counter = 1, n_fov
          if (stats_data% btdn(fov_counter, ch_counter) > 0) then
            sfov_pi(1, ch_counter) = sfov_pi(1, ch_counter) +      &  ! predictor "1"
                 bcor_coef% coef_set2(ch_counter, fov_counter) *   &
                 stats_data% btdn    (fov_counter, ch_counter)
            do pred_counter = 1, n_pred
              sfov_pi     (pred_counter + 1, ch_counter) =         &
                   sfov_pi(pred_counter + 1, ch_counter) +         &
                   bcor_coef% coef_set2(ch_counter, fov_counter) * &
                   stats_data% pred_i(pred_counter, fov_counter, ch_counter)
            end do
          end if
        end do
      end do

    endif

    !----------------------------------------
    ! Calculate coefficients for each channel
    !----------------------------------------
    allocate (matrix_a   (n_pred + 1, n_pred + 1))
    allocate (matrix_inv (n_pred + 1, n_pred + 1))
    allocate (vector_b   (n_pred + 1            ))
    allocate (tmp_b      (n_pred + 1            ))
    allocate (scal       (n_pred                ))
    do ch_counter = 1, n_chan
       if (btdn_ch(ch_counter) > 2._wp) then
          ! A correlation matrix is inverted here. Since the correlation
          ! calculated from two timesteps is always one for any two variables
          ! the matrix is only invertible if the number of data is larger than 1.
          !--------------------------------------------------------
          ! a) build matrix and right hand side for predictors used
          !--------------------------------------------------------
          !---------------------------------------
          ! loop over rows, skip unused predictors
          !---------------------------------------
          pred_counter = 0
          i = 0
          do pred_counter_i = 1, n_pred
             if (bcor_coef% pred_use (pred_counter_i, ch_counter) < 1) then
               pred_counter = pred_counter + n_pred - pred_counter_i + 1
               cycle
             endif
             !------------------------------------------
             ! loop over columns, skip unused predictors
             !------------------------------------------
             i = i + 1
             j = i - 1
             bcor_coef% i_p (i, ch_counter) = pred_counter_i
             do pred_counter_j = pred_counter_i, n_pred
                pred_counter = pred_counter + 1
                if(bcor_coef% pred_use (pred_counter_j, ch_counter) < 1) cycle
                j = j + 1
                !---------------------------------
                ! set elements of symmetric matrix
                !---------------------------------
                matrix_a(i, j) = pij_ch(pred_counter, ch_counter)
                matrix_a(j, i) = matrix_a(i, j)
             end do
             !-------------------------------------
             ! row and column for the constant term
             !-------------------------------------
             matrix_a(i  , j+1) = pi_ch    (pred_counter_i, ch_counter)
             matrix_a(j+1, i  ) = matrix_a (i, j + 1)
             !----------------
             ! right hand side
             !----------------
             vector_b(i) = btdscan_pi_ch(pred_counter_i  ,ch_counter) - &
                           sfov_pi      (pred_counter_i+1,ch_counter)
          end do
          !------------------------------------
          ! remaining entries for constant term
          !------------------------------------
          n = i
          matrix_a (n+1, n+1) = btdn_ch   (ch_counter)
          vector_b (n+1)      = btdscan_ch(ch_counter) - sfov_pi(1, ch_counter)

          !--------------
          ! subtract mean
          !--------------
          do i = 1, n
            do j = 1, n
              matrix_a (i,j) = matrix_a (i  ,j  ) &
                             - matrix_a (i  ,n+1) &
                             * matrix_a (j  ,n+1) &
                             / matrix_a (n+1,n+1)
            end do
            vector_b (i) = vector_b(i) - matrix_a (i  ,n+1) * vector_b (n+1) &
                                       / matrix_a (n+1,n+1)
            bcor_coef% p_mean (i,ch_counter) = matrix_a (i  ,n+1) &
                                             / matrix_a (n+1,n+1)
          end do
          bcor_coef% n  (ch_counter) = btdn_ch (ch_counter)

          !------
          ! scale
          !------
          do i = 1, n
            if (matrix_a (i,i) > 0) then
              scal (i) = 1._wp / sqrt (matrix_a (i,i))
            else
              scal (i) = 0._wp                  ! saveguard for zero on diagonal
            endif
            bcor_coef% p_scal  (i,ch_counter) = scal(i)
          end do
          do i = 1, n
            do j = 1, n
              matrix_a (i,j) = matrix_a (i,j) * (scal(i) * scal(j))
            end do
            vector_b(i) = vector_b(i) * scal(i)
          end do

          !---------------
          ! constrained BC
          !---------------
          select case (bcor_coef%constr_mode (ch_counter))
          case(1)
            vector_b(n+1) = vector_b(n+1) + bcor_coef%constr_alpha(ch_counter) * &
                                            bcor_coef%constr_b0   (ch_counter) * &
                                            matrix_a (n+1, n+1)
            vector_b(1:n+1) = vector_b(1:n+1) / (1._wp + bcor_coef%constr_alpha(ch_counter))
          case(2)
            vector_b(n+1) = vector_b(n+1) + bcor_coef%constr_alpha(ch_counter) * &
                                            bcor_coef%constr_b0   (ch_counter) * &
                                            matrix_a (n+1, n+1)
            vector_b(n+1) = vector_b(n+1) / (1._wp + bcor_coef%constr_alpha(ch_counter))
          case(3)
            vector_b(n+1) = bcor_coef%constr_b0   (ch_counter) * matrix_a (n+1, n+1)
          end select

          !--------------------------
          ! regularisation (diagonal)
          !--------------------------
          if (bcor_coef% nu > 0._wp) then
            do i = 1, n
              matrix_a(i,i) = matrix_a(i,i) * (1._wp + bcor_coef% nu ** 2)
            end do
          endif

          !---------------------------------
          ! b) solve set of linear equations
          !---------------------------------
          if (any (scal (:n) == 0._wp)) then
            matrix_inv (1:n,1:n) = 0._wp
            tmp_b      (1:n+1)   = 0._wp
          else
            matrix_inv(1:n,1:n) = inverse_rs (matrix_a (1:n,1:n))
            tmp_b(1:n) =  matmul(matrix_inv(1:n,1:n), vector_b(1:n))

            !--------
            ! rescale
            !--------
            do i = 1, n
              tmp_b(i) = tmp_b(i) * scal(i)
            end do

            !-----------------
            ! correct for mean
            !-----------------
            tmp_b(n+1) = - sum (tmp_b(1:n) * matrix_a (1:n,n+1)) &
                                           / matrix_a (n+1,n+1)  &
                         + vector_b (n+1)  / matrix_a (n+1,n+1)
          endif

          !---------------------------------------------------
          ! c) redistribute coefficients used for this channel
          !---------------------------------------------------
          i = 0
          do pred_counter = 1, n_pred
             if (bcor_coef% pred_use (pred_counter, ch_counter) < 1) then
                bcor_coef% coef_set1 (pred_counter, ch_counter) = 0
             else
               i = i + 1
               bcor_coef% coef_set1 (pred_counter, ch_counter) =  tmp_b(i)
             endif
          end do
          bcor_coef% coef_set1 (n_pred+1, ch_counter) = tmp_b (i+1)

          !-----------------------------------------------------
          ! return B matrix (scaled and mean predictors removed)
          !-----------------------------------------------------
          bcor_coef% p_b ( bcor_coef% n_b(ch_counter,3) + 1 &
                         : bcor_coef% n_b(ch_counter,3) +   &
                           bcor_coef% n_b(ch_counter,2)     )&
                         = pack (matrix_inv (1:n,1:n), .true.)
       end if

       ! Offset
       if (trim(bcor_coef% offset(ch_counter)) /= '') then
         call calc_offset(bccoef(:), trim(bcor_coef%offset(ch_counter)), &
                          bcor_coef%offset_v(ch_counter), offset)
         bcor_coef%offset_v(ch_counter) = offset
         bcor_coef% coef_set1 (n_pred+1, ch_counter) = &
              bcor_coef% coef_set1 (n_pred+1, ch_counter) - offset
         if (dace%lpio) then
           write(*,'(2x,a,i3.3,a,i2.2,a,i5,a,F12.5)') &
                'Evaluated offset "'//trim(bcor_coef% offset(ch_counter))//'" &
                &(sat=',bcor_coef%i%satid,',grid=',bcor_coef%i%grid,',chan=',ch_counter,'): ',offset
         end if
       end if
    end do

  end subroutine calc_bcor_coefs

!==============================================================================

  subroutine write_bcor_file (bcor, unit, file, meta, templ)
  !======================================================================
  ! writes a bias correction file in new self-explanatory namelist-format
  !======================================================================
  type (t_bcor_coef)        ,intent(in) :: bcor (:) ! bias correction data
  integer                   ,intent(in) :: unit     ! I/O unit to use
  character(len=*),optional ,intent(in) :: file     ! file to open
  logical         ,optional ,intent(in) :: meta     ! write meta data only
  logical         ,optional ,intent(in) :: templ    ! write template

    target                       :: bcor
    integer                      :: iset
    integer                      :: i
    integer                      :: n_chan
    type (t_bcor_coef)  ,pointer :: bc
    type (t_stats_data) ,pointer :: st
    logical                      :: l_all, frame, zero
    character(len=*) ,parameter  :: form_b = '(5(1x,e20.12))' ! write coeffs.
    character(len=*) ,parameter  :: form_s = '(5(1x,e20.12))' ! write statist.
    character(len=20)            :: form_u = '(a,__i3,/(16x,__i3))'
    character(len=80)            :: form   = ''
    integer                      :: ios
    !---------------------------
    ! process optional arguments
    !---------------------------
    frame = .true.                ! write namelist header + trailer
    l_all = .true.                ! write all data
    zero  = .false.               ! write zero number of entries
    if (present (meta)) then      ! write meta data only (for job output)
      l_all   = .not.meta
      frame = .not.meta
    endif
    if (present (templ)) then     ! write namelist template (no entries)
      l_all   = .not.templ
      zero  =      templ
    endif
    if (present (file)) then
      open (unit, file=file, action='write',iostat=ios)
      if (ios/=0) call finish('write_bcor_file','cannot open '//trim(file))
    endif
    !-----------------------------------------------
    ! loop over sets of bias-correction coefficients
    !-----------------------------------------------
    do iset = 1, size(bcor)
      bc => bcor (iset)
      !-------------
      ! write header
      !-------------
      write (unit,'( )')
      if (frame) write (unit,'(a)') ' &BIASCOR_COEFF'
      !----------------
      ! write meta data
      !----------------
      write (unit,'(a, i8,a)')    '   exp         =' ,bc% exp,        '  ! experiment number'
      write (unit,'(a,  a,a)')    '   created     = "',trim(bc% created),'"'
      write (unit,'(a,  a,a)')    '   modified    = "',trim(bc% modified),'"'
      write (unit,'(a,  a,a)')    '   first_date  = "',trim(bc% first_date),'"'
      write (unit,'(a,  a,a)')    '   last_date   = "',trim(bc% last_date),'"'
      write (unit,'(a,  a,a)')    '   last_update = "',trim(bc% last_update),'"'
      if (zero) then
        write (unit,'(a,f9.0,a)') '   entries     =' ,0.,              ' ! number of entries used for training'
      else
        write (unit,'(a,f9.0,a)') '   entries     =' ,bc% entries,     ' ! number of entries used for training'
      endif
      if (all(bc% t_decay(2:bc% i% n_chan) == bc% t_decay(1))) then
        write (unit,'(a,f9.0,a)') '   t_decay     =' ,bc% t_decay(1),  ' ! background info decay time (days)'
      else
        write(form,'("(a,",I3,"(1x,f9.0),a)")') bc% i% n_chan
        write (unit,form)         '   t_decay     =' ,bc% t_decay(1:bc%i%n_chan),&
             ' ! background info decay time (days)'
      end if
      write (unit,'(a, i8,a)')    '   mdlsfc_stat =' ,bc% mdlsfc_stat,'  ! surface flag mask for statistics'
      write (unit,'(a,f9.2,a)')   '   nu          =' ,bc% nu,          ' ! regularisation parameter'
      write (unit,'(a,f9.2,a)')   '   w_vbc       =' ,bc% w_vbc,       ' ! VarBC weight'
      write (unit,'(a,7x,l1,a)')  '   fov_bcorr   =' ,bc% fov_bcorr,  '  ! individual FOV bias correction'
      write (unit,'(a, i8,a)')    '   mode_wgt    =' ,bc% mode_wgt ,  '  ! mode for weighting of obs'
      do i = 1, bc% i% n_chan
        if (bc% constr_mode(i) > 0) then
          write (unit,'(a,i5.5,a,1x,i1,a)') '   constr_mode(',i,') =' ,bc% constr_mode(i),&
               ' ! constrained BC mode'
        end if
      end do
      do i = 1, bc% i% n_chan
        if (bc% constr_mode(i) > 0) then
          write (unit,'(a,i5.5,a,1x,f7.4,a)') '   constr_b0(',i,') =' ,bc% constr_b0(i),&
               ' ! "b0" in constrained BC, "Target bcor"'
        end if
      end do
      do i = 1, bc% i% n_chan
        if (bc% constr_mode(i) > 0) then
          write (unit,'(a,i5.5,a,1x,E13.6,a)') '   constr_alpha(',i,') =' ,bc% constr_alpha(i),&
               ' ! "alpha" in constrained BC, "strength of constraint"'
        end if
      end do
      do i = 1, bc% i% n_chan
        if (trim(bc% offset(i)) /= '') then
          write (unit,'(a,i5.5,a,a,a)') '   offset(',i,') = "',trim(bc% offset(i)),&
               '" ! offset calculation rule'
          write (unit,'(a,i5.5,a,1x,E13.6,a)') '   offset_v(',i,') = ',bc% offset_v(i),&
               '  ! offset value (calculated in last analysis)'
        end if
      end do
      write (unit,'(a, i8,a)')    '   sat_id      =' ,bc% i% satid,   '  ! satellite ID'
      write (unit,'(a, i8,a)')    '   grid        =' ,bc% i% grid,    '  ! RTTOV-ID of target instrument'
      write (unit,'(a, i8,a)')    '   n_instr     =' ,bc% i% n_instr, '  ! number of instruments, channels'
      do i = 1, bc% i% n_instr
        write (unit,'(a,i1,a,i8)')'   instr(',i,')    =',bc% i% instr (i)
        write (unit,'(a,i1,a,i8)')'   n_chan(',i,')   =',bc% i% n_ch_i(i)
      end do
      write (unit,'(a, i8,a)')    '   n_coeff     =' ,bc% n_coeff,     '  ! coefficients'
      write (unit,'(a, i8,a)')    '   n_fov       =' ,bc% i% n_fov,    '  ! number of fov'
      write (unit,'(a, i8,a)')    '   n_pred      =' ,bc% n_pred,      '  ! number of predictors'
      do i = 1, bc% n_pred
        write (unit,'(a,i2.2,a,a,a)')'   pred(',i,')     = "',trim(bc% pred(i)),'"'
        write (unit,'(a,i2.2,a,i8)') '   pred_p1(',i,')  =' ,bc% pred_p1 (i)
        write (unit,'(a,i2.2,a,i8)') '   pred_p2(',i,')  =' ,bc% pred_p2 (i)
      end do
      if (bc% n_pred > 0) then
        write (form_u( 4: 5),'(i2)') bc% n_pred
        write (form_u(15:16),'(i2)') bc% n_pred
        write (unit,form_u)          '   pred_use    =', bc% pred_use(1:bc% n_pred,1:bc% i% n_chan)
      endif
      do i = 1, bc% i% n_chan
        if (.not. (bc% fov_bc(i) == fov_bc_default)) &
             write (unit,'(a,i4.4,a,i2,1x,i2,1x,i2,1x,l1)') '   fov_bc(',i,')=' ,bc% fov_bc(i)
      end do
      write (unit,'(a)')             '   chan        ='
      write (unit,'(10i12)') bc% i% chan
      if (l_all) then
        !-------------------
        ! write coefficients
        !-------------------
        write (unit,'(a)')          '   coef        ='
        n_chan = bc% i% n_chan
        do i = 1, n_chan
          write (unit,form_b) bc% coef_set1 (:,i), bc% coef_set2 (i,:)
        end do
      endif
      !-----------------------------
      ! write accumulated statistics
      !-----------------------------
      st => bc% sta
      if (st% alloc) then
        if (.not.zero) then
          write (unit,'(a)')          '                    !'
          write (unit,'(a)')          '                    ! statistics:'
          write (unit,'(a)')          '                    !'
          write (unit,'(a)')          '   n_c         ='
          write (unit,'(10f12.2)') st% n
        endif
        if (l_all) then
          write (unit,'(a)')          '   n_fc        ='
          do i = 1, n_chan
            write (unit,form_s) st% btdn (:,i)
          end do
          write (unit,'(a)')          '   d_fc        ='
          do i = 1, n_chan
            write (unit,form_s) st% btdscan (:,i)
          end do
          write (unit,'(a)')          '   d2_fc       ='
          do i = 1, n_chan
            write (unit,form_s) st% d2 (:,i)
          end do
          write (unit,'(a)')          '   p_pfc       ='
          do i = 1, n_chan
            write (unit,form_s) st% pred_i (:,:,i)
          end do
          write (unit,'(a)')          '   dp_pfc      ='
          do i = 1, n_chan
            write (unit,form_s) st% btd_pred_i (:,:,i)
          end do
          write (unit,'(a)')          '   pp_pfc      ='
          do i = 1, n_chan
            write (unit,form_s) st% pred_ij (:,:,i)
          end do
        endif
      endif
      !-------------------------------
      ! write instantaneous statistics
      !-------------------------------

      if (l_all) then
        st => bc% sti
        if (st% alloc) then

          write (unit,'(a)')          '   i_n_c         ='
          write (unit,'(10f12.2)') st% n
          write (unit,'(a)')          '   i_n_fc        ='
          do i = 1, n_chan
            write (unit,form_s) st% btdn (:,i)
          end do
          write (unit,'(a)')          '   i_d_fc        ='
          do i = 1, n_chan
            write (unit,form_s) st% btdscan (:,i)
          end do
          write (unit,'(a)')          '   i_d2_fc       ='
          do i = 1, n_chan
            write (unit,form_s) st% d2 (:,i)
          end do
          write (unit,'(a)')          '   i_p_pfc       ='
          do i = 1, n_chan
            write (unit,form_s) st% pred_i (:,:,i)
          end do
          write (unit,'(a)')          '   i_dp_pfc      ='
          do i = 1, n_chan
            write (unit,form_s) st% btd_pred_i (:,:,i)
          end do
          write (unit,'(a)')          '   i_pp_pfc      ='
          do i = 1, n_chan
            write (unit,form_s) st% pred_ij (:,:,i)
          end do
        endif
      endif
      !--------------
      ! write trailer
      !--------------
      if (frame) write (unit,'(a)') ' /'
      write          (unit,'( )')
    end do
    if (present (file)) close (unit)

  end subroutine write_bcor_file

!------------------------------------------------------------------------------

  subroutine read_bcor_file (bcor, unit, ierr, file)
  !------------------------------------------------------------------
  ! read bias correction file in new self-explanatory namelist-format
  !------------------------------------------------------------------
  type (t_bcor_coef)        ,pointer     :: bcor (:) ! bias correction data
  integer                   ,intent(in)  :: unit     ! I/O unit to use
  integer                   ,intent(out) :: ierr     ! error return argument
  character(len=*),optional ,intent(in)  :: file     ! file to open

    !------------------------------------------------------------
    ! namelist /BIASCOR_COEFF/ (mirrors derived type t_bcor_coef)
    !------------------------------------------------------------

    integer           :: sat_id          ! satellite ID
    integer           :: exp             ! experiment number
    integer           :: n_instr         ! number of instruments in set
    integer           :: instr (m_instr) ! instrument(s) (rttov ?)
    integer           :: grid            ! RTTOV-ID of grid instrument
    integer           :: n_chan(m_instr) ! number of channels per instrument
    integer           :: n_pred          ! number of predictors
    integer           :: n_fov           ! number of scan-angle depend.entries
    integer           :: n_coeff         ! total number of coefficients
    character(len=12) :: created         ! creation date 'yyyymmddhhmm'
    character(len=12) :: modified        ! modification date
    character(len=12) :: first_date      ! first training date
    character(len=12) :: last_date       ! last  change date
    character(len=12) :: last_update     ! last  training date
    real(wp)          :: entries         ! number of entries for training
    real(wp)          :: t_decay(m_chan) ! background info decay time (days)
    integer           :: mdlsfc_stat     ! surface flag mask for statistics
    real(wp)          :: nu              ! regularisation parameter
    real(wp)          :: w_vbc           ! variational bias correction weight
    logical           :: fov_bcorr       ! bias correction for individual fovs
    integer           :: mode_wgt        ! mode for weighting of obs
    integer           :: constr_mode(m_chan)  ! contrained bc mode
    real(wp)          :: constr_b0(m_chan)    ! target bc for constrained bcor
    real(wp)          :: constr_alpha(m_chan) ! strength of constrained bcor
    character(len=300):: offset(m_chan)       ! string describing offset
    real(wp)          :: offset_v(m_chan)     ! offset value

    character(len=10) :: pred   (m_pred) ! predictors
    integer           :: pred_p1(m_pred) ! predictor parameter 1
    integer           :: pred_p2(m_pred) ! predictor parameter 2
    integer           :: pred_use(m_pred * m_chan) ! 0:unused,+-1:stats,1:bcorr
    type(t_fov_bc)    :: fov_bc  (         m_chan) ! see declaration of t_fov_bc above

    integer           :: chan   (m_chan) ! list of channels
    real(wp)          :: coef   (m_chan*(m_pred+m_fov+1))  ! coefficients
    !-----------------------
    ! accumulated statistics
    !-----------------------
    real(wp)          :: n_c                 (m_chan) ! # of entries
    real(wp)          :: n_fc          (m_fov*m_chan) ! # of entries / fov
    real(wp)          :: d_fc          (m_fov*m_chan) ! sum of deviations
    real(wp)          :: d2_fc         (m_fov*m_chan) ! sum of deviations**2
    real(wp), allocatable :: p_pfc  (:)               ! deviations*predictors
    real(wp), allocatable :: dp_pfc (:)               ! predictor product sums
    real(wp), allocatable :: pp_pfc (:)               ! predictor product sums
    !-------------------------
    ! instantaneous statistics
    !-------------------------
    real(wp)          :: i_n_c                 (m_chan)
    real(wp)          :: i_n_fc          (m_fov*m_chan)
    real(wp)          :: i_d_fc          (m_fov*m_chan)
    real(wp)          :: i_d2_fc         (m_fov*m_chan)
    real(wp), allocatable :: i_p_pfc  (:)
    real(wp), allocatable :: i_dp_pfc (:)
    real(wp), allocatable :: i_pp_pfc (:)
    real(wp)          :: t_decay_def ! default background info decay time (days)

    namelist /BIASCOR_COEFF/                                                &
      sat_id, exp, n_instr, instr, grid, n_chan, n_pred, n_fov, n_coeff,    &
      created, modified, first_date, last_date, last_update, entries,       &
      pred, pred_p1, pred_p2, pred_use, fov_bc, chan, coef, t_decay,        &
      mdlsfc_stat, nu, w_vbc, fov_bcorr, mode_wgt,                          &
      constr_mode, constr_b0, constr_alpha, offset, offset_v,               &
      n_c, d2_fc, n_fc, d_fc, p_pfc, dp_pfc, pp_pfc,                        &
      i_n_c, i_d2_fc, i_n_fc, i_d_fc, i_p_pfc, i_dp_pfc, i_pp_pfc
    !----------------
    ! local variables
    !----------------
    integer                     :: iset
    integer                     :: i, j
    integer                     :: nchan, nfov, npred, npp
    integer                     :: ncof
    type (t_bcor_coef) ,pointer :: bc
    type (t_bcor_coef) ,pointer :: tmp (:)
    logical                     :: lexit

    ierr = 0
    !---------------------------
    ! process optional arguments
    !---------------------------
    if (dace% lpio .and. present (file))                              &
       open (unit, file=file, action='read', status='old', iostat=ierr)
    call p_bcast (ierr, dace% pio)
    if (ierr/=0) then
       if (dace% lpio) write(0,*) 'read_bcor_file: Cannot open file: ', file
       return
    endif
    allocate (p_pfc    (m_pred*m_fov*m_chan))
    allocate (dp_pfc   (m_pred*m_fov*m_chan))
    allocate (pp_pfc   (m_pp  *m_fov*m_chan))
    allocate (i_p_pfc  (m_pred*m_fov*m_chan))
    allocate (i_dp_pfc (m_pred*m_fov*m_chan))
    allocate (i_pp_pfc (m_pp  *m_fov*m_chan))
    !-----------------------------------------------
    ! loop over sets of bias-correction coefficients
    !-----------------------------------------------
    iset = 0; if (associated(bcor)) iset = size (bcor)
    do
      !-------------------
      ! set default values
      !-------------------
      sat_id      = -1
      exp         = -1
      n_instr     =  0
      instr       = -1
      grid        =  0
      n_chan      =  0
      n_pred      =  0
      n_fov       =  0
      n_coeff     =  0
      created     = '000000000000'
      modified    = '000000000000'
      first_date  = '000000000000'
      last_date   = '000000000000'
      last_update = '000000000000'
      entries     = -1._wp
      t_decay     = -huge(0._wp)
      mdlsfc_stat =  0
      nu          =  0._wp
      w_vbc       =  1._wp
      fov_bcorr   = .true.
      mode_wgt    =  0
      constr_mode =  0
      constr_b0   =  0._wp
      constr_alpha=  0._wp
      offset      = ''
      offset_v    =  0._wp
      pred        = ''
      pred_p1     =  0
      pred_p2     =  0
      pred_use    =  1
      call construct(fov_bc)
      chan        = -1
      coef        =  0
      n_c         = -1
      n_fc        = -1
      d_fc        = -1
      d2_fc       = -1
      p_pfc       = -1
      dp_pfc      = -1
      pp_pfc      = -1
      i_n_c       = -1
      i_n_fc      = -1
      i_d_fc      = -1
      i_d2_fc     = -1
      i_p_pfc     = -1
      i_dp_pfc    = -1
      i_pp_pfc    = -1
      !------------------------------
      ! read namelist /BIASCOR_COEFF/
      !------------------------------
      lexit = .false.
      if (dace% lpio) then
         call position_nml ('BIASCOR_COEFF', unit= unit,  &
                                           status= ierr,  &
                                          lrewind= iset==0)
         select case (ierr)
         case (POSITIONED)
#if defined(__ibm__)
            read (unit ,nml=BIASCOR_COEFF, iostat=ierr)
            if (ierr/=0) then
               write(0,*) 'read_bcor_file: ERROR in namelist /BIASCOR_COEFF/'
               lexit = .true.
            endif
#else
            read (unit ,nml=BIASCOR_COEFF)
#endif
         case default
            ierr  = 0
            lexit = .true.
         end select
      end if
      call p_bcast (ierr,  dace% pio)
      call p_bcast (lexit, dace% pio)
      if (lexit) exit
      !-------------------
      ! Broadcast namelist
      !-------------------
      call p_bcast (sat_id      ,dace% pio)
      call p_bcast (exp         ,dace% pio)
      call p_bcast (n_instr     ,dace% pio)
      call p_bcast (instr       ,dace% pio)
      call p_bcast (grid        ,dace% pio)
      call p_bcast (n_chan      ,dace% pio)
      call p_bcast (n_pred      ,dace% pio)
      call p_bcast (n_fov       ,dace% pio)
      call p_bcast (n_coeff     ,dace% pio)
      call p_bcast (created     ,dace% pio)
      call p_bcast (modified    ,dace% pio)
      call p_bcast (first_date  ,dace% pio)
      call p_bcast (last_date   ,dace% pio)
      call p_bcast (last_update ,dace% pio)
      call p_bcast (entries     ,dace% pio)
      call p_bcast (t_decay     ,dace% pio)
      call p_bcast (mdlsfc_stat ,dace% pio)
      call p_bcast (nu          ,dace% pio)
      call p_bcast (w_vbc       ,dace% pio)
      call p_bcast (fov_bcorr   ,dace% pio)
      call p_bcast (mode_wgt    ,dace% pio)
      call p_bcast (constr_mode ,dace% pio)
      call p_bcast (constr_b0   ,dace% pio)
      call p_bcast (constr_alpha,dace% pio)
      call p_bcast (offset      ,dace% pio)
      call p_bcast (offset_v    ,dace% pio)
      call p_bcast (pred        ,dace% pio)
      call p_bcast (pred_p1     ,dace% pio)
      call p_bcast (pred_p2     ,dace% pio)
      call p_bcast (pred_use    ,dace% pio)
      call p_bcast (fov_bc      ,dace% pio)
      call p_bcast (chan        ,dace% pio)
      call p_bcast (coef        ,dace% pio)
      call p_bcast (n_c         ,dace% pio)
      call p_bcast (n_fc        ,dace% pio)
      call p_bcast (d_fc        ,dace% pio)
      call p_bcast (d2_fc       ,dace% pio)
      call p_bcast (p_pfc       ,dace% pio)
      call p_bcast (dp_pfc      ,dace% pio)
      call p_bcast (pp_pfc      ,dace% pio)
      call p_bcast (i_n_c       ,dace% pio)
      call p_bcast (i_n_fc      ,dace% pio)
      call p_bcast (i_d_fc      ,dace% pio)
      call p_bcast (i_d2_fc     ,dace% pio)
      call p_bcast (i_p_pfc     ,dace% pio)
      call p_bcast (i_dp_pfc    ,dace% pio)
      call p_bcast (i_pp_pfc    ,dace% pio)
      !-----------------------------------
      ! For old biascor files without grid
      !-----------------------------------
      if ((sat_id == 784).and.(grid <= 0)) grid = 3

      nchan = sum (n_chan (1: n_instr))
      ncof  = n_pred + n_fov + 1
      if (n_fov == 1) fov_bcorr = .false.

      ! modify t_decay if necessary
      if (all(t_decay(:) == -huge(0._wp))) t_decay(1) = 0._wp
      do i = 1, nchan
        if (t_decay(i) /= -huge(0._wp)) then
          t_decay_def = t_decay(i)
        else
          t_decay(i)  = t_decay_def
        end if
      end do

      if (last_update == '000000000000') last_update = last_date

      !---------------------
      ! copy to derived type
      !---------------------
      iset = iset + 1
      tmp => bcor
      allocate (bcor (iset))
      if (associated (tmp)) then
        bcor (1:iset-1) = tmp
        deallocate (tmp)
      endif
      bc => bcor (iset)
      bc% i% satid              = sat_id
      bc% i% grid               = grid
      bc% i% n_chan             = nchan
      bc% i% n_instr            = n_instr
      bc% i% n_fov              = n_fov
      bc% i% instr  (1:n_instr) = instr  (1:n_instr)
      bc% i% n_ch_i (1:n_instr) = n_chan (1:n_instr)
      bc% i% o_ch_i (1)         = 0
      do i = 2, n_instr
        bc% i% o_ch_i (i) = bc% i% o_ch_i (i-1) + bc% i% n_ch_i (i-1)
      end do
      allocate (bc% i% chan (nchan))
      bc% i% chan               = chan (1:nchan)

      allocate (bc% coef_set1 (n_pred+1, nchan))
      allocate (bc% coef_set2 (nchan, n_fov)   )
      allocate (bc% pred_use  (n_pred, nchan)  )
      allocate (bc% fov_bc    (nchan)          )
      bc% exp                = exp
      bc% n_pred             = n_pred
      bc% n_coeff            = n_coeff
      bc% created            = created
      bc% modified           = modified
      bc% t_decay            = t_decay
      bc% mdlsfc_stat        = mdlsfc_stat
      bc% nu                 = nu
      bc% w_vbc              = w_vbc
      bc% fov_bcorr          = fov_bcorr
      bc% mode_wgt           = mode_wgt
      bc% constr_mode        = constr_mode
      bc% constr_b0          = constr_b0
      bc% constr_alpha       = constr_alpha
      do i = 1, size(offset)
        bc% offset(i)        = trim(offset(i))
        bc% offset_v(i)      = offset_v(i)
      end do
      bc% first_date         = first_date
      bc% last_date          = last_date
      bc% last_update        = last_update
      bc% entries            = entries
      bc% pred               = pred
      bc% pred_p1            = pred_p1
      bc% pred_p2            = pred_p2
      bc% pred_use           = reshape (pred_use(1:n_pred * nchan), &
                                                (/ n_pred,  nchan /))
      bc% fov_bc             = fov_bc(1:nchan)
      do i = 1, nchan
        j = i-1
        bc% coef_set1 (:,i) = coef (j*ncof +1       : j*ncof +1+n_pred)
        bc% coef_set2 (i,:) = coef (j*ncof +2+n_pred: i*ncof          )
      end do
      !----------------------------
      ! read accumulated statistics
      !----------------------------
      if (bc% sta% alloc) call destruct (bc% sta)
      if (any(n_c >= 0._wp)) then
        call construct (bc%    sta,    &
                        bc%    n_pred, &
                        bc% i% n_fov,  &
                        nchan          )
        nfov  = bc% i% n_fov
        npred = bc%    n_pred
        npp   = npred*(npred+1)/2

        bc% sta% n          = reshape (n_c   (1:size(bc% sta% n)),          shape(bc% sta% n))
        bc% sta% d2         = reshape (d2_fc (1:size(bc% sta% d2)),         shape(bc% sta% d2))
        bc% sta% btdn       = reshape (n_fc  (1:size(bc% sta% btdn)),       shape(bc% sta% btdn))
        bc% sta% btdscan    = reshape (d_fc  (1:size(bc% sta% btdscan)),    shape(bc% sta% btdscan))
        bc% sta% pred_i     = reshape (p_pfc (1:size(bc% sta% pred_i)),     shape(bc% sta% pred_i))
        bc% sta% btd_pred_i = reshape (dp_pfc(1:size(bc% sta% btd_pred_i)), shape(bc% sta% btd_pred_i))
        bc% sta% pred_ij    = reshape (pp_pfc(1:size(bc% sta% pred_ij)),    shape(bc% sta% pred_ij))
      endif
      !------------------------------
      ! read instantaneous statistics
      !------------------------------
      if (bc% sti% alloc) call destruct (bc% sti)
      if (any(i_n_c >= 0._wp)) then
        call construct (bc%    sti,    &
                        bc%    n_pred, &
                        bc% i% n_fov,  &
                        nchan          )
        nfov  = bc% i% n_fov
        npred = bc%    n_pred
        npp   = npred*(npred+1)/2

        bc% sti% n          = reshape (i_n_c   (1:size (bc% sti% n)),         &
                                shape                  (bc% sti% n))
        bc% sti% d2         = reshape (i_d2_fc (1:size (bc% sti% d2)),        &
                                shape                  (bc% sti% d2))
        bc% sti% btdn       = reshape (i_n_fc  (1:size (bc% sti% btdn)),      &
                                shape                  (bc% sti% btdn))
        bc% sti% btdscan    = reshape (i_d_fc  (1:size (bc% sti% btdscan)),   &
                                shape                  (bc% sti% btdscan))
        bc% sti% pred_i     = reshape (i_p_pfc (1:size (bc% sti% pred_i)),    &
                                shape                  (bc% sti% pred_i))
        bc% sti% btd_pred_i = reshape (i_dp_pfc(1:size (bc% sti% btd_pred_i)),&
                                shape                  (bc% sti% btd_pred_i))
        bc% sti% pred_ij    = reshape (i_pp_pfc(1:size (bc% sti% pred_ij)),   &
                                shape                  (bc% sti% pred_ij))
      endif
    end do
    !---------------------------------
    ! fallback: read 1dvar-file format
    !---------------------------------
    if (iset == 0 .and. ierr == 0) then
      call finish('read_bcor_file','error reading bcor_file')
    endif
    if (.not.associated( bcor)) allocate (bcor(0))
    !-----------
    ! close file
    !-----------
    if (dace% lpio .and. present (file)) close (unit)

    ! Testing
    ! do i = 1, bc%i%n_instr
    !   do ic = bc%i%o_ch_i(i)+1,bc%i%o_ch_i(i)+bc%i%n_ch_i(i)
    !     call calc_residual_bias(bc,ic,resid,n)
    !     write(*,*) 'residual',bc%i%instr(i),bc%i%chan(ic),resid,n
    !   end do
    ! end do


  end subroutine read_bcor_file
!==============================================================================

  subroutine calc_bcor_stats (stats_data, obs_fg, pred_i, pred_ij, pred_use, &
                              ifov, ichan, wgt                               )
  !-----------------------------------------
  ! Calculate statistics for bias correction
  !-----------------------------------------
  type(t_stats_data) ,intent(inout) :: stats_data   ! Bias cor. statistics data
  real(wp)           ,intent(in)    :: obs_fg       ! Observation - first guess
  real(wp)           ,intent(in)    :: pred_i   (:) ! Predictors
  real(wp)           ,intent(in)    :: pred_ij  (:) ! Predictor products
  integer            ,intent(in)    :: pred_use (:) ! Predictor uasge flag
  integer            ,intent(in)    :: ifov         ! FOV
  integer            ,intent(in)    :: ichan        ! Channel
  real(wp)           ,intent(in)    :: wgt          ! obs. weight

    !---------------------------
    ! check for valid predictors
    !---------------------------
    if (all (pred_i /= rinvalid .or. pred_use == 0)) then
      if (ifov > ubound(stats_data%d2,1) .or. ifov < lbound(stats_data%d2,1)) then
        write(0,*) 'ifov,bounds',ifov,lbound(stats_data%d2,1),ubound(stats_data%d2,1)
        call finish('calc_bcor_stats','ifov out of range')
      end if
      if (ichan> ubound(stats_data%d2,2) .or. ichan< lbound(stats_data%d2,2)) then
        write(0,*) 'ifov,bounds',ifov,lbound(stats_data%d2,2),ubound(stats_data%d2,2)
        call finish('calc_bcor_stats','ichan out of range')
      end if
      !------------------
      ! update statistics
      !------------------
      stats_data% d2            (ifov, ichan) = &
      stats_data% d2            (ifov, ichan) + wgt * obs_fg * obs_fg
      stats_data% btdscan       (ifov, ichan) = &
      stats_data% btdscan       (ifov, ichan) + wgt * obs_fg
      stats_data% btdn          (ifov, ichan) = &
      stats_data% btdn          (ifov, ichan) + wgt
      stats_data% btd_pred_i (:, ifov, ichan) = &
      stats_data% btd_pred_i (:, ifov, ichan) + wgt * obs_fg * pred_i(:)
      stats_data% pred_i     (:, ifov, ichan) = &
      stats_data% pred_i     (:, ifov, ichan) + wgt * pred_i(:)
      stats_data% pred_ij    (:, ifov, ichan) = &
      stats_data% pred_ij    (:, ifov, ichan) + wgt * pred_ij
    endif

  end subroutine calc_bcor_stats

!------------------------------------------------------------------------------

  subroutine calc_bcor_t_stats (stats, bcor, rad)
  !-----------------------------------------
  ! Calculate statistics for bias correction
  !-----------------------------------------
  type(t_stats_data) ,intent(inout) :: stats (:) ! Bias cor. statistics data
  type(t_bcor_coef)  ,intent(in)    :: bcor  (:) ! Bias correction data
  type(t_radv)       ,intent(inout) :: rad   (:)  ! model state + observation info

    integer :: i
    do i = 1, size (stats)
      call calc_bcor_t_stat (stats(i), bcor(i), rad(i))
    end do

  end subroutine calc_bcor_t_stats

!------------------------------------------------------------------------------

  subroutine calc_bcor_t_stat (stats, bcor, rad)
  !-----------------------------------------
  ! Calculate statistics for bias correction
  !-----------------------------------------
  type(t_stats_data) ,intent(inout) :: stats  ! Bias cor. statistics data
  type(t_bcor_coef)  ,intent(in)    :: bcor   ! Bias correction data
  type(t_radv)       ,intent(inout) :: rad    ! model state + observation info

    integer  :: ip
    integer  :: ic, ic_pr
    integer  :: ifov
    real(wp) :: pred_2 (bcor% n_pred * (bcor% n_pred+1) / 2)
    real(wp) :: wgt
    integer  :: i,j,k

    !----------------------------------
    ! loop over fovs in monitoring file
    !----------------------------------
    ifov = 1
    do ip = 1, rad% n_rec
      if (bcor% fov_bcorr) ifov = rad% fov (ip) + 1
      !---------------------------------
      ! calculate products of predictors
      !---------------------------------
      k = 0
      do i = 1, bcor% n_pred
        do j = i, bcor% n_pred
          k = k + 1
          if ( any(bcor% pred(i) == chan_dep_pred) .or. &
               any(bcor% pred(j) == chan_dep_pred)) cycle
          if (rad% pred (i, ip, 1) == rinvalid .or. &
              rad% pred (j, ip, 1) == rinvalid      ) then
            pred_2 (k) = rinvalid
          else
            pred_2 (k) = rad% pred (i, ip, 1) * rad% pred (j, ip, 1)
          endif
        end do
      end do
      !-------------------
      ! loop over channels
      !-------------------
      do ic = 1, rad% i% n_chan
        !-------------------------------------------
        ! update statistics for radiances,
        ! currently +++ hardcoded +++ :
        !   not rejected, no clouds, no land, no ice
        !-------------------------------------------
        if (      rad% not_rej (ic,ip)           .and. &
                  rad% cloudy  (ic,ip)     == 2  .and. &
            iand (rad% mdlsfc  (   ip), bcor% mdlsfc_stat) == 0 ) then
          wgt = 1._wp
          if (btest(bcor% mode_wgt, WGT_ERR)) then
            if (rad% e_fg (ic,ip) <= 0._wp) call finish('calc_bcor_t_stat','background error <= 0.')
            if (rad% e_obs(ic,ip) <= 0._wp) call finish('calc_bcor_t_stat','observation error <= 0.')
            wgt = wgt * (1._wp / (rad% e_fg(ic,ip)**2 + rad% e_obs(ic,ip)**2))
          end if
          if (btest(bcor% mode_wgt, WGT_LAT)) then
            wgt = wgt * (1._wp + sin(rad%dlat(ip)*d2r)**2)
          end if
          if (btest(bcor% mode_wgt, WGT_SZA)) then
            wgt = wgt * (1._wp / cos(min(abs(rad%stzen(ip)),80._wp)*d2r))
          end if

          if (bcor%chan_dep_pred) then
            ic_pr = ic
            k = 0
            do i = 1, bcor% n_pred
              do j = i, bcor% n_pred
                k = k + 1
                if (any(bcor% pred(i) == chan_dep_pred) .or.&
                    any(bcor% pred(j) == chan_dep_pred)) then
                  if (rad% pred (i, ip, ic_pr) == rinvalid .or. &
                      rad% pred (j, ip, ic_pr) == rinvalid      ) then
                    pred_2 (k) = rinvalid
                  else
                    pred_2 (k) = rad% pred (i, ip, ic_pr) * rad% pred (j, ip, ic_pr)
                  endif
                end if
              end do
            end do
          else
            ic_pr = 1
          end if

          call calc_bcor_stats (stats,                                   &
                                rad% bt_obs (ic,ip) - rad% bt_fg (ic,ip),&
                                rad% pred (:,ip,ic_pr), pred_2,          &
                                bcor% pred_use (:,ic),                   &
                                ifov,                                    &
                                ic,                                      &
                                wgt                                      )
        endif
      end do
    end do

  end subroutine calc_bcor_t_stat

!------------------------------------------------------------------------------

  subroutine calc_h_fg (rad)
  !------------------------------
  ! Calculate geopotential height
  !------------------------------
  type(t_radv) ,intent(inout) :: rad    ! model state + observation info

    real(wp) :: tv   (rad% n_lev, rad% n_rec)      ! virtual temperature
    real(wp) :: lnp  (rad% n_lev)                  ! ln(p)
    real(wp) :: dlnp (rad% n_lev-1)                ! ln(p) difference
    integer  :: k, i, i_p

    real(wp) ,parameter :: c = R_DRY / (2 * gacc)  ! constant factor

    !-------------------
    ! allocate component
    !-------------------
    if (rad% n_rec == 0) return
    if (.not.associated (rad%h_fg)) allocate (rad%h_fg (rad%n_lev, rad%n_rec))
    !-------------------
    ! calculate tv, dlnp
    !-------------------
    tv   =  tv_t_q (rad% t_fg, rad% q_fg)
    do i = 1, rad%n_rec
      i_p = min(i, ubound(rad%p,2))
      lnp  =  log    (rad% p(:,i_p))
      dlnp = (lnp (2:) - lnp (:rad% n_lev-1)) * c
      rad% h_fg (rad% n_lev, :) = 0._wp
      do k = rad% n_lev -1, 1, -1
        rad% h_fg (k,i) = rad% h_fg (k+1,i) + dlnp(k) * (tv(k,i) + tv(k+1,i))
      end do
    end do

    ! lnp  =  log    (rad% p(:,1))
    ! dlnp = (lnp (2:) - lnp (:rad% n_lev-1)) * c
    ! !------------
    ! ! integration
    ! !------------
    ! rad% h_fg (rad% n_lev, :) = 0._wp
    ! do k = rad% n_lev -1, 1, -1
    !   rad% h_fg (k,:) = rad% h_fg (k+1,:) + dlnp(k) * (tv(k,:) + tv(k+1,:))

    ! end do

  end subroutine calc_h_fg

!==============================================================================

  subroutine scale_bcor_stats (b, a_time, r_time, days)
  type (t_bcor_coef) ,intent(inout) :: b         ! pointer to bias corr.coeff.
  character(len=12)  ,intent(in)    :: a_time    ! Analysis time (hhhhmmddhhmm)
  character(len=12)  ,intent(in)    :: r_time    ! Run time      (hhhhmmddhhmm)
  real(wp)           ,intent(in)    :: days      ! days since last update
  !-------------------------------------
  ! Rescale bias correction statistics
  ! (adjust effective number of entries)
  !-------------------------------------

    real(wp)            :: fi        ! individual scaling factor for a channel
    integer             :: i         ! channel index
    real(wp) ,parameter :: m = 2._wp ! minimum number of entries to keep

    if (b% sta% alloc .and. days > 0._wp) then
      do i = 1, size (b% sta% n)
        if (b% t_decay(i) > 0._wp) then
          fi = exp ( - days / b% t_decay(i))
          if (b% sta% n(i) > 0._wp)  &
               fi = max (fi, min (min_entries/b% sta% n(i), 1._wp))
          b% sta% n              (i) = fi * b% sta% n              (i)
          b% sta% btdn         (:,i) = fi * b% sta% btdn         (:,i)
          b% sta% btdscan      (:,i) = fi * b% sta% btdscan      (:,i)
          b% sta% d2           (:,i) = fi * b% sta% d2           (:,i)
          b% sta% pred_i     (:,:,i) = fi * b% sta% pred_i     (:,:,i)
          b% sta% btd_pred_i (:,:,i) = fi * b% sta% btd_pred_i (:,:,i)
          b% sta% pred_ij    (:,:,i) = fi * b% sta% pred_ij    (:,:,i)
        end if
      end do
      !--------------------------
      ! adjust modification dates
      !--------------------------
      b% modified   = r_time
      b% last_date  = a_time
    endif

  end subroutine scale_bcor_stats

!==============================================================================

  subroutine merge_bcor_stats (b, st, a_time, r_time, days, force_add)
  !---------------------------------
  ! Merge bias correction statistics
  !---------------------------------
  type (t_bcor_coef)  ,pointer    :: b         ! pointer to bias corr.coeff.
  type (t_stats_data) ,pointer    :: st        ! pointer to statistics data
  character(len=12)   ,intent(in) :: a_time    ! Analysis time (hhhhmmddhhmm)
  character(len=12)   ,intent(in) :: r_time    ! Run time      (hhhhmmddhhmm)
  real(wp)            ,intent(in) :: days      ! days since last update
  logical   ,optional ,intent(in) :: force_add ! force adding of statistics

    logical                  :: l_force_add
    integer                  :: i

    if (present(force_add)) then
       l_force_add = force_add
    else
       l_force_add = .false.
    end if

    if (any (st% n > 0)) then
       !------------------------------
       ! update accumulated statistics
       !------------------------------
       if (.not. b% sta% alloc) then
          !-----------------------------
          ! old file is empty: just copy
          !-----------------------------
          call construct (b%    sta,    &
                          b%    n_pred, &
                          b% i% n_fov,  &
                          b% i% n_chan  )
          b% sta = st
          !------------------
          ! set creation date
          !------------------
          b% created    = r_time
          b% modified   = r_time
          b% first_date = a_time
          b% last_date  = a_time
          b% last_update= a_time
       else
          !-----------------
          ! rescale old file
          !-----------------
          call scale_bcor_stats (b, a_time, r_time, days)
          !----
          ! add
          !----
          do i = 1, b% i% n_chan
            if (b% t_decay(i) >= 0._wp .or. l_force_add) then
              b% sta% n              (i) = b% sta% n              (i) + st% n              (i)
              b% sta% btdn         (:,i) = b% sta% btdn         (:,i) + st% btdn         (:,i)
              b% sta% btdscan      (:,i) = b% sta% btdscan      (:,i) + st% btdscan      (:,i)
              b% sta% d2           (:,i) = b% sta% d2           (:,i) + st% d2           (:,i)
              b% sta% pred_i     (:,:,i) = b% sta% pred_i     (:,:,i) + st% pred_i     (:,:,i)
              b% sta% btd_pred_i (:,:,i) = b% sta% btd_pred_i (:,:,i) + st% btd_pred_i (:,:,i)
              b% sta% pred_ij    (:,:,i) = b% sta% pred_ij    (:,:,i) + st% pred_ij    (:,:,i)
            end if
          end do

          b% last_date  = a_time
          b% last_update= a_time
          !----------------------
          ! set modification date
          !----------------------
          b% modified   = r_time
       endif
       !--------------------------------
       ! update instantaneous statistics
       !--------------------------------
       if (.not. b% sti% alloc)      &
         call construct (b%    sti,    &
                         b%    n_pred, &
                         b% i% n_fov,  &
                         b% i% n_chan  )
       b% sti = st
    else
       !-----------------------------------------------------
       ! deallocate instantaneous statistics if no valid data
       !-----------------------------------------------------
       if (b% sti% alloc) call destruct (b% sti)
    endif
  end subroutine merge_bcor_stats
!==============================================================================

  subroutine apply_bias_corr_rads (rad, bcor, l_stat_equiv, mask_rs)
  use mo_mpi_dace, only: p_sum, p_or, p_max
  type(t_radv)      ,intent(inout)          :: rad (:)      ! model state + observation info
  type(t_bcor_coef) ,intent(in)             :: bcor(:)      ! bias cor. coefficients
  logical           ,intent(in)   ,optional :: l_stat_equiv ! selection of data to which the
                                                            ! bias correction is to be applied to:
                                                            ! T: equivalent to the data used for the
                                                            !    calculation of the bias correction
                                                            ! F: as many data as possible [default]
    logical,         intent(in), optional :: mask_rs(:)

    integer :: i,j,k
    integer :: ierr

    ! Check residual bcor
    real(wp), allocatable :: cmb(:)
    real(wp), allocatable :: wsum(:)
    integer,  allocatable :: n(:)
    real(wp) :: wgt
    integer  :: n_chan, satid, grid, chan
    logical  :: l_check



    do i = 1, size (rad)
      if (present(mask_rs)) then
        if (.not.mask_rs(i)) CYCLE
      end if

!       l_check = associated(rad(i)%e_fg)
!       l_check = p_or(l_check)
!      l_check = .true.
      l_check = .false.

      if (rad(i)% n_rec > 0) then

        call apply_bias_corr (bcor,               &
             bcor(i)% i% satid,  &
             bcor(i)% i% grid,   &
             rad (i)% fov,       &
             rad (i)% pred,      &
             rad (i)% bcor_,     &
             ierr,               &
             l_stat_equiv=l_stat_equiv,&
             ideb=rad(i)%ideb)

        if (ierr/=0) call finish('apply_bias_corr_rads','apply_bias_corr failed')

        !     !+++ where clause doesnt work on NEC SX9 (segmentation fault) +++
        !     where (rad(i)% bcor_ /= rinvalid .and. rad(i)% bt_obs > 0._wp)
        !       rad(i)% bt_bcor = rad(i)% bt_obs - rad(i)% bcor_
        !     elsewhere
        !       rad(i)% bt_bcor = rad(i)% bt_obs
        !     endwhere

        do j=1,rad(i)% n_rec
          do k=1,rad(i)% i% n_chan
            if (rad(i)% bcor_ (k,j) /= rinvalid .and. &
                 rad(i)% bt_obs(k,j)  > 0._wp          ) then
              rad(i)% bt_bcor(k,j) = rad(i)% bt_obs(k,j) - rad(i)% bcor_(k,j)
            else
              rad(i)% bt_bcor(k,j) = rad(i)% bt_obs(k,j)
            endif
          end do
        end do

      end if

      if (l_check) then
        n_chan = rad(i)% i% n_chan
        n_chan = p_max(n_chan)
        allocate(cmb(n_chan), n(n_chan), wsum(n_chan))
        cmb  = 0._wp
        wsum = 0._wp
        n    = 0

        do j=1,rad(i)% n_rec
          do k=1,rad(i)% i% n_chan
            if (rad(i)% bcor_ (k,j) /= rinvalid .and. &
                rad(i)% bt_obs(k,j)  > 0._wp    .and. &
                rad(i)% not_rej (k,j)           .and. &
                rad(i)% cloudy  (k,j) == 2      .and. &
                iand(rad(i)% mdlsfc  (  j), bcor(i)% mdlsfc_stat) == 0 ) then
              !write(2000+dace%pe,*) 'use',j,k,rad(i)% fov(j),rad(i)%obsnum(j),rad(i)% bt_fg (k,j),rad(i)% bt_obs (k,j)
              wgt = 1._wp
              if (btest(bcor(i)% mode_wgt, WGT_ERR)) then
                if (rad(i)% e_fg (k,j) <= 0._wp) call finish('calc_bcor_t_stat','background error <= 0.')
                if (rad(i)% e_obs(k,j) <= 0._wp) call finish('calc_bcor_t_stat','observation error <= 0.')
                wgt = wgt * (1._wp / (rad(i)% e_fg(k,j)**2 + rad(i)% e_obs(k,j)**2))
              end if
              if (btest(bcor(i)% mode_wgt, WGT_LAT)) then
                wgt = wgt * (1._wp + sin(rad(i)%dlat(j)*d2r)**2)
              end if
              if (btest(bcor(i)% mode_wgt, WGT_SZA)) then
                wgt = wgt * (1._wp / cos(min(abs(rad(i)%stzen(j)),80._wp)*d2r))
              end if
              cmb (k) = cmb (k) + wgt * (rad(i)% bt_bcor(k,j) - rad(i)% bt_fg(k,j))
              wsum(k) = wsum(k) + wgt
              n   (k) = n   (k) + 1
            end if
          end do
        end do

        satid = p_max(rad(i)%i%satid)
        grid  = p_max(rad(i)%i%grid )
        do k=1,n_chan
          cmb (k) = p_sum(cmb (k))
          wsum(k) = p_sum(wsum(k))
          n   (k) = p_sum(n   (k))
          chan = 0
          if (rad(i)%i%n_chan >= n_chan) chan = rad(i)%i%chan(k)
          chan  = p_max(chan)
          ! if (dace%lpio .and. n(k) > 0) then
          !   print*,'%bcor',i,satid,grid,k,chan,n(k),wsum(k),cmb(k)/(wsum(k))
          ! end if
        end do
        deallocate(cmb,n,wsum)
      end if

    end do

  end subroutine apply_bias_corr_rads

!------------------------------------------------------------------------------

  subroutine apply_bias_corr (bcor_coefs, cur_satid, cur_grid, fov, pred, &
       &bcor, ierr, l_stat_equiv, ideb)
    !--------------------------------------------------
    ! apply bias-correction
    ! 'bcorr' sign is 1dvar/ECMWF convention, not 3dvar
    !--------------------------------------------------
    type(t_bcor_coef),intent(in)             :: bcor_coefs(:) ! bias cor. coefficients
    integer          ,intent(in)             :: cur_satid     ! Given satellite ID
    integer          ,intent(in)             :: cur_grid      ! Given grid ID
    integer          ,intent(in)             :: fov (:)       ! FOV (starts with 0)
    real(wp)         ,intent(in)             :: pred(:,:,:)   ! predictors (n_pred, n_rec)
    real(wp)         ,intent(out)            :: bcor (:,:)    ! Bias correction (offsets)
    integer          ,intent(out)            :: ierr          ! Error flag
    logical          ,intent(in)   ,optional :: l_stat_equiv  ! selection of data to which the
                                                              ! bias correction is to be applied to:
                                                              ! T: equivalent to the data used for the
                                                              !    calculation of the bias correction
                                                              ! F: as many data as possible [default]
    integer,          intent(in),   optional :: ideb(:)

    integer                ::  nrec             ! number of records
    integer                ::  npred            ! Number of predictors
    integer                ::  cur_indx         ! Array index matching
                                                ! current satellite and grid
    integer                ::  cur_chan         ! Counter for channels
    integer                ::  cur_fov          ! Current field of view
    integer                ::  cur_rec          ! Counter for records
    integer                ::  cur_pred         ! Counter for predictors
    integer                ::  ich              ! index for chan.dep. predictors
    real(wp)               ::  cur_coef         ! Bias correction coefficient
    type(t_bcor_coef)      ::  sat_bcor_coef    ! Bias corr. coefs for current
                                                ! satellite
    logical ,allocatable   ::  valid_pred (:)   ! Flags for valid predictors
    logical                ::  l_equiv          ! see l_stat_equiv
    integer                ::  nchan            !+++ required for sxf90/rev430
    logical                ::  ldeb

    ierr  = 0
    bcor  = rinvalid
    if (present(l_stat_equiv)) then ; l_equiv=l_stat_equiv ; else ; l_equiv=.false. ; endif

    !----------------------------------
    ! Find record for current satellite
    !----------------------------------
    cur_indx = set_indx (bcor_coefs(:)% i, satid=cur_satid, grid=cur_grid)

    if (cur_indx < 1) then
      print *,                                                        &
        'apply_bias_corr: cannot find cur_indx for satid ',cur_satid,&
        ' and cur_grid ',cur_grid
          ierr = 1
          return
    endif

    sat_bcor_coef = bcor_coefs(cur_indx)
    nrec  = size (fov)
    npred = sat_bcor_coef% n_pred
    allocate(valid_pred(npred))
    do cur_rec = 1, nrec
       if (present(ideb)) then
          ldeb = any(ideb == cur_rec)
       else
          ldeb = .false.
       end if

       ! Check if all predictors are valid
       where (pred(:, cur_rec, 1) == rinvalid)
          valid_pred(:) = .false.
       elsewhere
          valid_pred(:) = .true.
       end where
       cur_fov = fov (cur_rec) + 1
       if (sat_bcor_coef% fov_bcorr .and. &
            (cur_fov < lbound(sat_bcor_coef% coef_set2,2) .or. &
             cur_fov > ubound(sat_bcor_coef% coef_set2,2))) then
         write(0,*) 'cur_fov=',cur_fov
         call finish('apply_bias_corr', 'Invalid FOV number')
       end if
       !-------------------
       ! Loop over channels
       !-------------------
       nchan = sat_bcor_coef% i% n_chan
       do cur_chan = 1, nchan
         if (sat_bcor_coef% fov_bcorr) then
           !---------------------------------------------------
           ! Scan angle dependent coefficient
           !---------------------------------------------------
           cur_coef = sat_bcor_coef% coef_set2 (cur_chan, cur_fov)
         else
           cur_coef = 0._wp
         end if

         if (ldeb) write(usd,*) 'debug_spot bc_calc',cur_rec,cur_chan,cur_coef

         if ((.not. l_equiv .and. &
              all(valid_pred(:) .or. sat_bcor_coef% pred_use(:,cur_chan)< 1)) .or.&
              (l_equiv .and. &
              all(valid_pred(:) .or. sat_bcor_coef% pred_use(:,cur_chan)==0))) then
!              !---------------------------------------------------
!              ! Scan angle dependent coefficient
!              !---------------------------------------------------
!              cur_coef = sat_bcor_coef% coef_set2(cur_chan, cur_fov)

           do cur_pred = 1, npred
             ich = 1
             if (sat_bcor_coef%chan_dep_pred) then
               ich = cur_chan
               valid_pred(cur_pred) = (pred(cur_pred, cur_rec, ich) /= rinvalid)
             end if
             if (valid_pred(cur_pred) .and. sat_bcor_coef% pred_use(cur_pred,cur_chan)==1) then
               cur_coef = cur_coef +                               &
                    sat_bcor_coef% coef_set1(cur_pred, cur_chan) * &
                    pred(cur_pred, cur_rec, ich)
               if (ldeb) write(usd,*) 'debug_spot bc_calc pred',cur_rec,cur_chan,cur_pred,cur_coef,&
                    sat_bcor_coef% coef_set1(cur_pred, cur_chan),pred(cur_pred, cur_rec, ich)
             end if
           end do
           if (sat_bcor_coef% fov_bc(cur_chan)% const_pred) &
                cur_coef = cur_coef + sat_bcor_coef% coef_set1(npred+1, cur_chan)
           if (ldeb) write(usd,*) 'debug_spot bc_calc pred_const',cur_rec,cur_chan,cur_coef,&
                sat_bcor_coef% coef_set1(npred+1, cur_chan)
         end if
         bcor(cur_chan, cur_rec) = cur_coef

       end do
     end do
   end subroutine apply_bias_corr

!=======================================================================

  subroutine calc_predictors_rads (rad, bcor, mask_rs)
    !----------------------------------------------
    ! calculate predictors from derived type t_radv
    !----------------------------------------------
    type(t_radv)      ,intent(inout)        :: rad (:) ! model state + observation info
    type(t_bcor_coef) ,intent(in)           :: bcor(:) ! Bias correction data
    logical,           intent(in), optional :: mask_rs(:)
    integer :: i
    do i = 1, size(rad)
      if (present(mask_rs)) then
        if (.not.mask_rs(i)) CYCLE
      end if
      call calc_predictors_rad (rad(i), bcor(i))
    end do
  end subroutine calc_predictors_rads

!-----------------------------------------------------------------------

  subroutine calc_predictors_rad (rad, bcor)
    !----------------------------------------------
    ! calculate predictors from derived type t_radv
    !----------------------------------------------
    type(t_radv)       ,intent(inout) :: rad     ! model state + observation info
    type(t_bcor_coef)  ,intent(in)    :: bcor    ! Bias correction data
    integer           :: nch
    integer, pointer  :: ib_ch(:)    ! bitmask for P_THICK_TR predictors
    logical           :: fg_present

    if (rad% n_rec == 0) return
    if (.not.associated (rad% pred)) then
      if (bcor%chan_dep_pred) then
        nch = rad%i%n_chan
      else
        nch = 1
      end if
      allocate (rad% pred (bcor% n_pred, rad% n_rec, nch))
    end if
    if (.not.associated(rad%tr))  allocate(rad%tr(0,0,0))
    if (associated(rad%i%bc%ib_ch)) then
      ib_ch => rad%i%bc%ib_ch
    else
      allocate(ib_ch(0))
    end if

    !-------------------------------------------
    ! check for presence of first guess profiles
    !-------------------------------------------
    fg_present = associated (rad% t_fg) .and. &
                 associated (rad% q_fg) .and. &
                 associated (rad% p   )

    if (fg_present) then

      !-----------------------------------
      ! profiles present, calculate height
      !-----------------------------------
      if (atm_thick_vers == 0) then
        if (rad%i%gopts%lev_mode > 0) then
          call finish('calc_predictors_rad','atm_thick_vers==0 only &
               &supported for calc. on RTTOV-levels (lev_mode==0)')
        end if
        call calc_h_fg(rad)
      end if
      call calc_predictors (rad% pred, bcor,                 &
                            p       = rad% p,                &
                            h       = rad% h_fg,             &
                            t       = rad% t_fg,             &
                            q       = rad% q_fg,             &
                            ps      = rad% ps_fg,            &
                            gp_surf = rad% gp_fg,            &
                            ts      = rad% ts_fg,            &
                            v10     = rad% v10_abs_fg,       &
                            bt      = rad% bt_fg,            &
                            obt     = rad% bt_obs,           &
                            dlon    = rad% dlon,             &
                            dlat    = rad% dlat,             &
                            scanang = rad% stzen,            & !!
                            tr      = rad% tr,               &
                            ib_ch   = ib_ch,                 &
                            plev    = rad% plev,             &
                            instr_t = rad% instr_temp(1,:),  &
                            orb_ph  = rad% orb_phase,        &
                            ideb    = rad% ideb              )
      !!!!!! rad%stzen is the satellite zenith angle, not the scan angle !!!!!!
    else
      !--------------------
      ! no profiles present
      !--------------------
      call calc_predictors (rad% pred, bcor, &
                             bt= rad% bt_fg, &
                            obt= rad% bt_obs )
    endif

    if (.not.associated(rad%i%bc%ib_ch)) deallocate(ib_ch)

  end subroutine calc_predictors_rad

!------------------------------------------------------------------------------

  subroutine calc_predictors (pred, biasc, p, h, t, q, ps, gp_surf, &
                              ts, v10, bt, obt, chan_id, dlon, dlat,&
                              scanang, tr, ib_ch, plev, instr_t, orb_ph, ideb)
    !---------------------------------------------------------
    ! Calculates predictors as specified in 'biasc' meta data.
    ! Vector version: 2nd index is number of profiles,
    !                 1st index is predictor number (pred)
    !                              level      (p, h, t, q)
    !                           or channel index (bt, obt)
    !---------------------------------------------------------
    real(wp)           ,intent(out) :: pred(:,:,:) ! predictors to calculate
    type(t_bcor_coef)  ,intent(in ) :: biasc      ! bias correction meta data
    real(wp) ,optional ,intent(in ) :: p    (:,:) ! pressure levels     (Pa)
    real(wp) ,optional ,intent(in ) :: h    (:,:) ! geopotential height (m)
    real(wp) ,optional ,intent(in ) :: t    (:,:) ! temperature profile (K)
    real(wp) ,optional ,intent(in ) :: q    (:,:) ! specif.hum. profile (kg/kg)
    real(wp) ,optional ,intent(in ) :: ps     (:) ! surface pressure
    real(wp) ,optional ,intent(in ) :: gp_surf(:) ! surface geopotential
    real(wp) ,optional ,intent(in ) :: ts     (:) ! surface temperature
    real(wp) ,optional ,intent(in ) :: v10    (:) ! 10 m wind speed
    real(wp) ,optional ,intent(in ) :: bt   (:,:) ! model brightness temperature
    real(wp) ,optional ,intent(in ) :: obt  (:,:) ! observed brightness temp.
    integer  ,optional ,intent(in ) :: chan_id(:) ! channel IDs
    real(wp) ,optional ,intent(in ) :: dlon   (:) ! longitude
    real(wp) ,optional ,intent(in ) :: dlat   (:) ! latitude
    real(wp) ,optional ,intent(in ) :: scanang(:) ! scan angle
    real(wp) ,optional ,intent(in ) :: tr (:,:,:) ! transm factors for P_THICK_TR
    integer  ,optional ,intent(in)  :: ib_ch  (:) ! bitmask for P_THICK_TR
    real(wp) ,optional ,intent(in ) :: plev (:,:) ! plevel of channel
    real(wp) ,optional ,intent(in ) :: instr_t(:) ! instr. temp.
    real(wp) ,optional ,intent(in ) :: orb_ph (:) ! orbit phase
    integer  ,optional ,intent(in ) :: ideb   (:) ! channel IDs

    real(wp), parameter   :: pi = dacos(-1._wp)
    integer               :: i                      ! predictor index variable
    integer               :: j                      ! channel   index variable
    integer               :: k                      ! level     index variable
    integer               :: i_p
    integer               :: i_pr_tr                ! index of P_THICK_TR pred.
    integer               :: ic                     ! channel index
    integer               :: npred                  ! number of predictors
    integer               :: nrec                   ! number of FOVs
    integer               :: nl                     ! number of levels
    integer               :: ndeb                   ! number of debug spots
    integer               :: ierr                   ! error status variable
    real(wp)              :: layer(2), l            ! layer boundary
    real(wp)              :: pl(100)                ! layer boundary
    real(wp), allocatable :: gp(:,:)   ! layer boundary
    integer               :: npl                     ! number of required levels
    integer               :: c_chan                 ! counter for channel
    integer               :: pred_chan              ! predictor channel index
    integer               :: lastlev (size(pred,2)) ! last valid level
    logical               :: pred_use               ! predictor is used ?
    logical               :: lchk

    character(len=100)    :: fname

    ierr  = NO_ERROR
    pred  = rinvalid
    npred = biasc% n_pred
    nrec  = size(pred,2)
    nl    = size(p,1)
    if (npred <= 0) return
    if (npred >  size(pred,1)) call finish ('calc_predictors',   &
                                            'n_pred > size(pred)')
    if (present(ideb)) then
      ndeb = count(ideb > 0)
    else
      ndeb = 0
    end if

    ! Check for monotonic pressure profiles
    lchk = .false.
    do i = 1, npred
      pred_use = any (biasc% pred_use (i,:) /= 0)
      if (pred_use .and. any(biasc%pred(i) == &
           (/P_THICK, P_THICK_TR, P_THICK_TRF, P_IWV/))) lchk = .true.
      if (lchk) exit
    end do
    if (lchk) then
      do i = 1, size(p,2)
        if (any(p(1:nl-1,i) >= p(2:nl,i))) then
          write(0,*) 'profile',i
          write(0,*) 'p',1,p(1,i)
          do j = 2, nl
            write(0,*) 'p',j,p(j,i),p(j,i) > p(j-1,i)
          end do
          call finish('calc_predictors', 'Non-monotonic pressure profile!')
        end if
      end do
    end if

    if (atm_thick_vers > 0) then
      ! Determine required layer heights
      npl = 0
      do i = 1, npred
        pred_use = any (biasc% pred_use (i,:) /= 0)
        if (pred_use .and. any(biasc% pred(i) == (/P_THICK, P_THICK_TR, P_THICK_TRF/))) then
          if (npl+2 > size(pl)) call finish('calc_predictors','too small array "pl"')
          if (.not. any(pl(1:npl) == biasc% pred_p1(i))) then
            npl = npl + 1
            pl(npl) = biasc% pred_p1(i)
          end if
          if (.not. any(pl(1:npl) == biasc% pred_p2(i))) then
            npl = npl + 1
            pl(npl) = biasc% pred_p2(i)
          end if
        end if
      end do
      if (npl > 0) then
        ! sort plevs
        do i = npl-1,1,-1
          do j = 1,i
            if (pl(j) > pl(j+1)) then
              l        = pl(j)
              pl(j)   = pl(j+1)
              pl(j+1) = l
            end if
          end do
        end do
        allocate(gp(npl,nrec))
        call calc_gp(p, t, q, ps, gp_surf, pl(1:npl), gp, ideb=ideb)
      end if
    end if

    i_pr_tr = 0
    do i = 1, npred
      pred_use = any (biasc% pred_use (i,:) /= 0)
      if (.not. pred_use) cycle
      select case (biasc% pred(i))
      case (P_THICK,P_THICK_TR,P_THICK_TRF)
        !----------
        ! thickness
        !----------
        layer = (/biasc% pred_p1(i), biasc% pred_p2(i)/)
        if (atm_thick_vers == 0) then
          call calc_atm_thick (layer, p, h, pred(i,:,1), ierr, ideb=ideb)
        else
          pred(i,:,1) = 0._wp
          do j = 1, npl
            if (biasc% pred_p2(i) == pl(j)) then
              pred(i,:,1) = pred(i,:,1) + gp(j,:)
            elseif (biasc% pred_p1(i) == pl(j)) then
              pred(i,:,1) = pred(i,:,1) - gp(j,:)
            end if
          end do
        end if
        if (any(biasc% pred(i) == (/P_THICK_TR,P_THICK_TRF/))) then
          i_pr_tr = i_pr_tr + 1
          pred(i,:,2:) = spread(pred(i,:,1),2,size(pred,3)-1)
          do ic = 1, size(ib_ch)
            if (btest(ib_ch(ic),i_pr_tr)) then
              pred(i,:,ic) = pred(i,:,ic) * tr(i_pr_tr,ic,:)
            end if
          end do
        end if
      case(P_TRANS_T,P_TRANS_TS,P_TRANS_P,P_TRANS_DP)
        i_pr_tr = i_pr_tr + 1
        do ic = 1, size(ib_ch)
          if (btest(ib_ch(ic),i_pr_tr)) then
            pred(i,:,ic) = tr(i_pr_tr,ic,:)
            if (ndeb > 0) write(usd,*) 'trans_t',ic,i_pr_tr,pred(i,ideb(1:ndeb),ic)
          end if
        end do
      case (P_PLEV)
        !---------------
        ! plevel of obs.
        !---------------
        do ic = 1, size(pred,3)
          pred (i,:,ic) = plev(ic,:)
          write(fname,'(A,"_",I5.5,"_",I5.5,"_",I2.2,"_",I3.3)') &
               trim(biasc% pred(i)),biasc% pred_p1(i),biasc% pred_p2(i),ic,dace%pe
          open(10,file=trim(fname))
          do j = 1, nrec
            write(10,*) pred(i,j,ic)
          end do
          close(10)
        end do
      case (P_IWV)
        !------------------------
        ! integrated water vapour
        !------------------------
        if (present (ps)) then
          do j = 1, nrec
            i_p = min(j, size(p,2))
            lastlev(j) = 1
            do k = 2, size (q,1)
              if (ps(j) > p (k,i_p)) lastlev(j) = k
            end do
          end do
        else
          lastlev = size (q,1)
        endif
        if (atm_thick_vers == 0) then
          call calc_col_cont_kgkg (t, q, h, p, q, lastlev, pred(i,:,1), ierr)
        else
          call calc_col_cont(p, q, ps, pred (i,:,1), ierr)
        end if
      case (P_TSKIN)
        !--------------------------
        ! model surface temperature
        !--------------------------
        pred (i,:,1) = ts
      case (P_SCANANG)
        !-----------
        ! scan angle
        !-----------
        pred(i,:,1) = scanang ** biasc% pred_p1(i) * &
                    dcos(biasc% pred_p2(i)*(dlat+90._wp)/360._wp * 2._wp*pi)
      case (P_SPHER)
        !--------------------
        ! spherical harmonics
        !--------------------
        pred (i,:,1) = y_lm (biasc% pred_p1(i), biasc% pred_p2(i), dlat, dlon)
      case (P_V10)
        !----------
        ! 10 m wind
        !----------
        if (present (v10)) pred (i,:,1) = v10
      case (P_INSTR_T)
        !-----------------
        ! instrument temp.
        !-----------------
        pred (i,:,1) = instr_t
      case (P_ORB_PH)
        !------------
        ! orbit phase
        !------------
        if (biasc% pred_p1(i) < 1) then
          pred (i,:,1) = orb_ph
        else
          pred (i,:,1) = cos(biasc% pred_p1(i) * orb_ph)
        end if
      case (P_PROD)
        !------------------------
        ! product of 2 predictors
        !------------------------
        if (biasc% pred_p1(i) >= i .or. biasc% pred_p2(i) >= i) &
          call finish('calc_predictors (product)','define arguments first')
        where (pred (biasc% pred_p1(i),:,1) /= rinvalid .and. &
               pred (biasc% pred_p2(i),:,1) /= rinvalid )     &
          pred (i,:,1) = pred (biasc% pred_p1(i),:,1) * pred (biasc% pred_p2(i),:,1)
      case (P_BT, P_OBS_BT)
        !-----------------------------
        ! brightness temperature
        !-----------------------------
        if (present(chan_id)) then
          !---------------------------------------------------
          ! Use channel IDs to find matching observation index
          ! (channel may be not present)
          !---------------------------------------------------
          j = -1
          pred_chan = biasc% pred_p1(i) * f_ins + biasc% pred_p2(i)
          do c_chan = 1, ubound(chan_id, 1)
            if(chan_id(c_chan) == pred_chan) then
              j = c_chan
              exit
            end if
          end do
        else
          !-------------------------------------------------------
          ! Use bias corr. data to find matching observation index
          !-------------------------------------------------------
          j = chan_indx (biasc% pred_p1(i), biasc% pred_p2(i), biasc% i)
          if (j < 1) then
            write(0,*)   'calc_predictors(2):  invalid instr/channel:',&
                         biasc% pred_p1(i), biasc% pred_p2(i)
            write(0,*)   '                     instruments:',&
                         biasc% i% instr(1:biasc%i%n_instr)
            write(0,*)   '                     channels:',&
                         biasc% i% chan(1:biasc%i%n_chan)
            call finish ('calc_predictors(2)','invalid instr/channel')
          endif
        end if
        !---------------------------------------------
        ! copy Tb for valid channel index and Tb value
        !---------------------------------------------
        if (j > 0) then
          if (biasc% pred(i) == P_BT)     pred (i,:,1) = bt  (j,:)
          if (biasc% pred(i) == P_OBS_BT) pred (i,:,1) = obt (j,:)
          where (pred (i,:,1) <= 0._wp .or. pred (i,:,1) >= 500._wp) &
                                          pred (i,:,1) = rinvalid
        endif
      case default
        write(0,*) 'predictor ',i
        write(0,*) 'dataset ',biasc% i% satid, biasc%i %grid
        call finish ('calc_predictors','invalid predictor "'//biasc% pred(i)//'"')
      end select
      if (biasc%chan_dep_pred .and. all(biasc% pred(i)/=chan_dep_pred)) then
        pred(i,:,2:) = spread(pred(i,:,1),2,size(pred,3)-1)
      end if

      if (ierr /= NO_ERROR) then
        write (0,*)  'calc_predictors:  error:', biasc% pred(i),', ierr=',ierr
        call finish ('calc_predictors','error:'//biasc% pred(i))
      endif
      if (present(ideb)) then
        if (any(ideb > 0)) then
          do j = 1, size(pred,2)
            if (any(ideb == j)) write(usd,*) 'debug_spot pred ',biasc% pred(i),&
                 biasc% pred_p1(i),biasc% pred_p2(i),pred(i,j,:)
          end do
        end if
      end if
    end do

  end subroutine calc_predictors

  recursive character(len=80) function pred_name(i, bc) result(name)
    integer,           intent(in) :: i
    type(t_bcor_coef), intent(in) :: bc
    character(len=20) :: str
    name = '***'
    if (i < 1) return
    if (bc% pred(i) == P_PROD) then
      name = trim(pred_name(bc%pred_p1(i),bc))//'*'//trim(pred_name(bc%pred_p2(i),bc))
    else
      name = bc% pred(i)
      select case(bc% pred(i))
      case(P_THICK)
        write(str,*) int(bc%pred_p1(i)) ; str = adjustl(str)
        name = trim(name)//'_'//trim(str)
        write(str,*) int(bc%pred_p2(i)) ; str = adjustl(str)
        name = trim(name)//'-'//trim(str)
      case(P_SCANANG)
        if (bc%pred_p1(i) /= 1) then
          write(str,*) int(bc%pred_p1(i)) ; str = adjustl(str)
          name = trim(name)//'^'//trim(str)
        end if
        if (bc%pred_p2(i) > 0) then
          write(str,*) int(bc%pred_p2(i)) ; str = adjustl(str)
          name = trim(name)//'*cos('//trim(str)//'*(lat+90)/360*2*pi)'
        end if
      case (P_BT, P_OBS_BT, P_SPHER)
        if (bc% pred(i) == P_SPHER) name = 'Y'
        write(str,*) int(bc%pred_p1(i)) ; str = adjustl(str)
        name = trim(name)//'('//trim(str)
        write(str,*) int(bc%pred_p2(i)) ; str = adjustl(str)
        name = trim(name)//','//trim(str)//')'
      end select
    end if

  end function pred_name

  !=======================================================================
  !> routine calculates the thickness of an atmospheric layer between
  !> two pressure levels for a set of profiles.\n

  !-----------------------------------------------------------------------
  ! - history:
  !-----------------------------------------------------------------------
  !> \version <b>1.0 | 2009-11-10 | M. Schwaerz</b>\n
  !> Initial release

  !-----------------------------------------------------------------------
  ! interface description:
  !-----------------------------------------------------------------------
  !> \param[in] layerBord  a vector conataining the top and the bottom
  !>                       pressure level between which\n the thickness of the
  !>                       atmospheric layer should be calculated;
  !>                       dim: (2)\n
  !>                       unit: [Pa]\n
  !> \param[in] presGrid   pressure grid\n
  !>                       dim: (nLevs)\n
  !>                       unit: [Pa]
  !> \param[in] hGrid      height grid\n
  !>                       dim: (nLevs)\n
  !>                       unit: [m]
  !> \param[out] atm_thick calculated atmospheric thickness\n
  !>                       unit: [m]\n
  !> \param[out] ierr      0 if everything was ok - else an error code > 0
  !> \param[in] lin_Ext    determines the kind of extrapolation outside the grid borders.
  !>                       Possible options:\n
  !>                       .true. : a linear extrapolation is performed\n
  !>                       .false.: a constant extrapolation is performed\n
  !>                       optional; default: .true. (i.e., linear extrapolation)
  !-----------------------------------------------------------------------
  subroutine calc_atm_thick (layerBord, presGrid, hGrid, atm_thick, ierr, lin_Ext, ideb)
    real (kind=wp), intent(in)           :: layerBord(2)
    real (kind=wp), intent(in)           :: presGrid(:,:)   ! dim (nLevs)
    real (kind=wp), intent(in)           :: hGrid   (:,:)   ! dim (nLevs)
    real (kind=wp), intent(out)          :: atm_thick(:)
    integer,        intent(out)          :: ierr
    logical,        intent(in), optional :: lin_Ext
    integer,        intent(in), optional :: ideb(:)
    ! local arguments:
    logical        :: loc_lin_Ext
    real (kind=wp) :: hBord  (2, size(atm_thick))
    real (kind=wp) :: h_noEx (2, size(atm_thick))
    integer        :: idx1   (2, size(atm_thick))
    integer        :: idx2   (2, size(atm_thick))
    real (kind=wp) :: int_Fac(2, size(atm_thick))
    integer        :: iRun1, iRun2, dimG, nProf, tmpIdx, i_p, n_p
    ierr   = NO_ERROR
    dimG          = size (hGrid,1)
    nProf         = size (hGrid,2)
    n_p           = size (presGrid,2)
    hBord   (:,:) = 0.0_wp
    h_noEx  (:,:) = 0.0_wp
    atm_thick (:) = 0.0_wp

    if (present(lin_Ext)) then
      loc_lin_Ext = lin_Ext
    else
      loc_lin_Ext = .true.
    end if
    !--- check input data:
    if ( checkInput ) then
      if ( ( size (hGrid,1)   /= dimG  ) .or. &
           ( size (hGrid,2)   /= nProf ) .or. &
           ( size (atm_thick) /= nProf ) ) then
        write(0,*) 'hGrid',shape(hGrid)
        write(0,*) 'presGrid',shape(presGrid)
        write(0,*) 'atm_thick',shape(atm_thick)
        ierr = ERR_DIM
        return
      endif
      do iRun1 = 1, nProf
        i_p = min(n_p, iRun1)
        if ( .not. ( all (presGrid (1:dimG-1,i_p  ) < presGrid (2:dimG,i_p  )) .or. &
                     all (presGrid (1:dimG-1,i_p  ) > presGrid (2:dimG,i_p  )) ) .or. &
             .not. ( all (hGrid    (1:dimG-1,iRun1) < hGrid    (2:dimG,iRun1)) .or. &
                     all (hGrid    (1:dimG-1,iRun1) > hGrid    (2:dimG,iRun1)) ) ) then
          write(0,*) 'ERR_MONOTON',iRun1,dimG
          do irun2 = 1, dimG
            write(0,*) 'h,p',irun2,presGrid(irun2,i_p),hgrid(irun2,irun1)
          end do
          ierr =  ERR_MONOTON
          return
        endif
        if (any(irun1 == ideb)) then
          do irun2 = 1, dimG
            write(usd,*) 'debug_spot bc_calc_pred atm_thick h,p',irun2,presGrid(irun2,i_p),hgrid(irun2,irun1)
          end do
        end if

      enddo
    endif
    if (layerBord (1) == layerBord (2)) then
      return
    endif
    !--- perform calc:
    !--- perform interpolation for both heights:
    do iRun1 = 1, nprof
      i_p = min(n_p, iRun1)
      do iRun2 = 1, 2
        h_noEx(irun2,irun1) = min(max(layerBord(irun2), min( presGrid (1,i_p), presGrid (dimG,i_p))), &
             max(presGrid(1,i_p), presGrid(dimG,i_p)))
      enddo
      if (presGrid (dimG,i_p) < layerBord(2)) atm_thick(irun1) = rinvalid
      if (presGrid (   1,i_p) > layerBord(1)) atm_thick(irun1) = rinvalid
    end do
    idx1(:,:) = 1
    idx2(:,:) = dimG


    do iRun2 = 1, nProf
      i_p = min(n_p, iRun2)
      if (atm_thick(iRun2) > rinvalid) then
        do iRun1 = 1, 2
          do while (idx2 (iRun1,iRun2) - idx1 (iRun1,iRun2) > 1)
            tmpIdx           = (idx2 (iRun1,iRun2) + idx1 (iRun1,iRun2)) / 2
            if (presGrid (tmpIdx,i_p) > h_noEx (iRun1,iRun2)) then
              idx2 (iRun1,iRun2)  = tmpIdx
            else
              idx1 (iRun1,iRun2)  = tmpIdx
            endif
          enddo
        enddo
      endif
    enddo

    if (loc_lin_ext) then
      do iRun2 = 1, nProf
        i_p = min(n_p, iRun2)
        if (atm_thick(iRun2) > rinvalid) then
          do iRun1 = 1, 2
            int_Fac(iRun1, iRun2)  = (layerBord (iRun1) - presGrid (idx1 (iRun1,iRun2), i_p)) / &
                                     (presGrid (idx2(iRun1,iRun2), i_p) - &
                                      presGrid (idx1(iRun1,iRun2), i_p))
            hBord  (iRun1, iRun2)  = (1.0_wp - int_Fac  (iRun1, iRun2)) * hGrid (idx1(iRun1,iRun2), iRun2) + &
                                      int_Fac  (iRun1, iRun2)  * hGrid (idx2(iRun1,iRun2), iRun2)
          enddo
        endif
      enddo
    else
      do iRun2 = 1, nProf
        if (atm_thick(iRun2) > rinvalid) then
          do iRun1 = 1, 2
            int_Fac(iRun1, iRun2)  = (h_noEx (iRun1,iRun2) - presGrid (idx1 (iRun1,iRun2), i_p)) / &
                                     (presGrid (idx2(iRun1,iRun2), i_p) - &
                                      presGrid (idx1(iRun1,iRun2), i_p))
            hBord  (iRun1, iRun2)  = (1.0_wp - int_Fac  (iRun1, iRun2)) * hGrid (idx1(iRun1,iRun2), iRun2) + &
                                      int_Fac  (iRun1, iRun2)  * hGrid (idx2(iRun1,iRun2), iRun2)
          enddo
        endif
      enddo
    endif
    do iRun2 = 1, nProf
      if (atm_thick (iRun2) > rinvalid) then
        atm_thick(iRun2) = abs (hBord(2, iRun2) - hBord(1, iRun2))
      endif
      if (any(irun2 == ideb)) write(usd,*) 'debug_spot bc_calc_pred hbord',atm_thick(irun2),hbord(:,irun2)
    enddo

  end subroutine calc_atm_thick


  subroutine calc_gp(p, t, q, psurf, gp_surf, pl, hgt, ideb)
    real(wp), intent(in)  :: p      (:,:)
    real(wp), intent(in)  :: t      (:,:)
    real(wp), intent(in)  :: q      (:,:)
    real(wp), intent(in)  :: psurf  (:)
    real(wp), intent(in)  :: gp_surf(:)
    real(wp), intent(in)  :: pl     (:)
    real(wp), intent(out) :: hgt    (:,:)
    integer,  intent(in), optional :: ideb(:)

    real(wp), parameter   :: c = R_DRY / gacc
    real(wp), parameter   :: tgrad = -0.0065_wp ! Clim.. temp. grad. [K/m]
    real(wp), allocatable :: tv     (:,:)
    real(wp), allocatable :: lnp    (:)
    real(wp), allocatable :: dlnp   (:)
    real(wp), allocatable :: h      (:)
    real(wp)              :: lnps, lnpl
    real(wp)              :: hsurf, w
    real(wp)              :: dlnp_, tsurf, qsurf, tvsurf, dz, t_, beta
    logical               :: ld, lfound
    integer               :: i,j,k
    integer               :: i_p,ksurf
    integer               :: nrec, nl, np
    ! atm_thick_vers ("remove" boundary layer)
    real(wp)              :: hbl, pbl, tbl, tvbl
    integer               :: kbl

    hgt = -1._wp

    nrec = size(t,2)
    nl   = size(p,1)
    np   = size(pl)

    allocate(tv(nl,nrec), lnp(nl), dlnp(nl-1), h(nl))
    tv = tv_t_q(t, q)

    do i = 1, nrec
      ld = .false.
      if (present(ideb)) ld = any(ideb == i)

      i_p = min(i, size(p,2))
      lnp  = log(p(:,i_p))
      lnps = log(psurf(i))
      dlnp = lnp(2:) - lnp(:nl-1)

      if (ld) write(usd,*) 'debug_spot calc_gp i',i,nrec,i_p

      ! Find surface
      ! ksurf is level below surface, ksurf-1 is lowest level above surface
      ksurf = nl + 1
      do k = nl, 1, -1
        if (p(k,i_p) >= psurf(i)) ksurf = k
      end do

      ! calc. surface values
      hsurf = gp_surf(i)/gacc
      if (ksurf > nl) then
        ! extrapolate to surface
        w = (lnps - lnp(nl)) / (lnp(nl) - lnp(nl-1))
        tsurf = t(nl,i) + w * (t(nl,i) - t(nl-1,i))
        qsurf = q(nl,i) + w * (q(nl,i) - q(nl-1,i))
      else
        ! interpolate to surface
        w = (lnp(ksurf) - lnps) / (lnp(ksurf) - lnp(ksurf-1))
        tsurf = t(ksurf,i) + w * (t(ksurf-1,i) - t(ksurf,i))
        qsurf = q(ksurf,i) + w * (q(ksurf-1,i) - q(ksurf,i))
      end if
      tvsurf = tv_t_q(tsurf, qsurf)
      if (ld) write(usd,*) 'debug_spot calc_gp surf',ksurf,psurf(i),hsurf,tsurf,qsurf,tvsurf

      ! Integrate hydrostatic equation
      h(ksurf-1) = hsurf + (lnps - lnp(ksurf-1)) * c * 0.5_wp * (tv(ksurf-1,i) + tvsurf)
      if (ld) write(usd,*) 'debug_spot calc_gp p,h',ksurf-1,p(ksurf-1,i_p), h(ksurf-1)
      do k = ksurf-2, 1, -1
        h(k) = h(k+1) + dlnp(k) * c * 0.5_wp * (tv(k,i) + tv(k+1,i))
        if (ld) write(usd,*) 'debug_spot calc_gp p,h',k,p(k,i_p), h(k), dlnp(k), (tv(k,i) + tv(k+1,i))*0.5_wp
      end do

      ! TODO
      if (atm_thick_vers == 3) then
        hbl = hsurf + 1000._wp ! RF: Assume "boundary layer" 1000m thick
                               ! RF: here the "boundary layer" is the layer that can be assumed to be
                               ! RF: free of daily cycle of temperature
        kbl = ksurf-1
        do k = ksurf-1, 1, -1
          if (h(k) >= hbl) then
            kbl = k
            exit
          end if
        end do
        if (kbl < nl) then
          w = (h(kbl) - hbl) / (h(kbl) - h(kbl+1))
          tvbl = (1._wp - w) * tv (kbl,i) + w * tv (kbl+1,i)
          tbl  = (1._wp - w) * t  (kbl,i) + w * t  (kbl+1,i)
          pbl  = (1._wp - w) * lnp(kbl  ) + w * lnp(kbl+1  )
          pbl  = exp(pbl)
        else ! should never happen
          tvbl = tv (kbl,i)
          tbl  = t  (kbl,i)
          pbl  = p  (kbl,i_p)
        end if
        if (ld) write(usd,*) 'debug_spot calc_gp bl_adapt',kbl,hbl,tvbl,tbl,pbl

        ! Integrate downwards with climatological temp. grad.
        beta = tvbl/tbl  ! Tv = beta * T
        do k = kbl + 1, nl
          t_ = tbl * (p(k,i_p)/pbl) ** (-c*tgrad*beta)
          dz = (t_ - tbl) / tgrad
          if (ld) write(usd,*) 'debug_spot calc_gp bl_adapt h',k,p(k,i_p),h(k),'->',hbl+dz,dz
          h(k) = hbl + dz
        end do
        ! For consistency we have to adapt hsurf,tsurf,tvsurf.
        ! The adapted values are purely hypothetical value.
        if (ksurf > nl) then
          ! extrapolate to surface
          w = (lnps - lnp(nl)) / (lnp(nl) - lnp(nl-1))
          hsurf = h(nl) + w * (h(nl) - h(nl-1))
        else
          ! interpolate to surface
          w = (lnp(ksurf) - lnps) / (lnp(ksurf) - lnp(ksurf-1))
          hsurf = h(ksurf) + w * (h(ksurf-1) - h(ksurf))
        end if
        tsurf  = tbl  + (hsurf - hbl) * tgrad
        tvsurf = tvbl + (hsurf - hbl) * tgrad
      end if

      ! Calculate height on required levels "pl"
      do j = 1, np
        lnpl = log(pl(j))
        if (ld) write(usd,*) 'debug_spot calc_gp pl',j,pl(j),p(1,i_p),p(nl,i_p)
        if  (pl(j) < p(1,i_p)) then
          ! upwards extrapolation: assume constant Tv
          dlnp_ = lnp(1) - lnpl
          hgt(j,i) = h(1) + dlnp_ * c * tv(1,i)
          if (ld) write(usd,*) 'debug_spot calc_gp hgt extrap_up',hgt(j,i),dlnp_,tv(1,i),c * tv(1,i)
        elseif (pl(j) > psurf(i)) then
          ! extrapolation below surface
          select case(atm_thick_vers)
          case(1)
            ! constant Tv
            dlnp_ = lnpl - lnps
            hgt(j,i) = hsurf - dlnp_ * c * tvsurf
            if (ld) write(usd,*) 'debug_spot calc_gp hgt extrap_down',hgt(j,i),hsurf,dlnp_,tvsurf,tvsurf*c
          case(2,3)
            ! climatological temp. gradient (0.65K/km = 0.00065K/m) and constant q
            ! analytic solution follows from integrating hydrostatic equ.
            ! dp = -rho*g*dz with rho = p /(R*Tv) and Tv = beta * T = beta * (tsurf  + tgrad * z)
            beta = tvsurf/tsurf ! Tv = beta * T
            t_ = tsurf * (pl(j)/psurf(i)) ** (-c*tgrad*beta)
            dz = (t_ - tsurf) / tgrad
            hgt(j,i) = hsurf + dz
            if (ld) write(usd,*) 'debug_spot calc_gp hgt extrap_down ana',hgt(j,i),hsurf,dz,t_
          end select
        elseif (pl(j) > p(ksurf-1,i_p)) then
          ! pl in lowest layer between surface and ksurf-1
          w = (lnpl - lnp(ksurf-1)) / (lnps - lnp(ksurf-1))
          hgt(j,i) = h(ksurf-1) + w * (hsurf - h(ksurf-1))
          if (ld) write(usd,*) 'debug_spot calc_gp hgt extrap_lowest_lev',hgt(j,i),w,hsurf,h(ksurf-1),&
               lnpl,lnp(nl),lnps
        else
          ! default
          lfound = .false.
          do k = 1, ksurf-2
            if (p(k,i_p) <= pl(j) .and. p(k+1,i_p) >= pl(j)) then
              w = (lnpl - lnp(k)) / (lnp(k+1) - lnp(k))
              hgt(j,i) = h(k) + w * (h(k+1) - h(k))
              if (ld) write(usd,*) 'debug_spot calc_gp hgt',hgt(j,i),w,h(k),h(k+1)
              lfound = .true.
              exit
            end if
          end do
          if (.not. lfound) then
            do k = 1, ksurf-1
              print*,'p',k,p(k,i_p)
            end do
            print*,'pl',pl(j)
            call finish('calc_gp','failed to find layer.')
          end if
        end if
      end do
    end do

  end subroutine calc_gp

  subroutine calc_col_cont(p, q, psurf, col_cont, ierr)
    real (kind=wp), intent(in)  :: p         (:,:)    ! (nLevs, nProf)
    real (kind=wp), intent(in)  :: q         (:,:)    ! (nLevs, nProf)
    real (kind=wp), intent(in)  :: psurf       (:)    ! (       nProf)
    real (kind=wp), intent(out) :: col_cont    (:)    ! (       nProf)
    integer,        intent(out) :: ierr

    real(wp) :: w, qsurf
    integer :: nrec, nl
    integer :: i, k, i_p, ksurf

    ! q = m_w/m ; m_w: mass of water, m: mass of moist air  (all masses per unit area)
    ! p = g * m -> dp = g * dm
    ! We have to compute the integral over (m_w/m) dm = q dp/g

    ierr  = NO_ERROR
    col_cont(:) = 0._wp
    nrec = size(q,2)
    nl   = size(p,1)
    do i = 1, nrec
      i_p = min(i, size(p,2))
      ! Find surface
      ! ksurf is level below surface, ksurf-1 is lowest level above surface
      ksurf = nl + 1
      do k = nl, 1, -1
        if (p(k,i_p) >= psurf(i)) ksurf = k
      end do
      if (ksurf > nl) then
        ! extrapolate to surface
        w = (log(psurf(i)) - log(p(nl,i_p))) / (log(p(nl,i_p)) - log(p(nl-1,i_p)))
        qsurf = q(nl,i) + w * (q(nl,i) - q(nl-1,i))
      else
        ! interpolate to surface
        w = (log(p(ksurf,i_p)) - log(psurf(i))) / (log(p(ksurf,i_p)) - log(p(ksurf-1,i_p)))
        qsurf = w * q(ksurf-1,i) + (1._wp - w) * q(ksurf,i)
      end if
      !col_cont(i) = 0.5_wp * (qsurf+q(ksurf-1,i)) * (psurf(i)-p(ksurf-1,i))
      col_cont(i) = 0._wp
      do k = ksurf-2, 1, -1
        col_cont(i) = col_cont(i) + 0.5_wp * (q(k,i)+q(k+1,i)) * (p(k+1,i_p)-p(k,i_p))
      end do
      col_cont(i) = col_cont(i)/gacc
    end do

  end subroutine calc_col_cont


  !=======================================================================
  !> routine calculates the column content of an atmospheric component
  !> for a set of profiles in kg*m$^{-2}$.\n

  !-----------------------------------------------------------------------
  ! - history:
  !-----------------------------------------------------------------------
  !> \version <b>1.0 | 2009-11-10 | M. Schwaerz</b>\n
  !> Initial release

  !-----------------------------------------------------------------------
  ! interface description:
  !-----------------------------------------------------------------------
  !> \param[in] temp         temperature profile\n
  !>                         dim: (nLevs,nProf)\n
  !>                         unit: [K]
  !> \param[in] humi         specific humidity profile\n
  !>                         dim: (nLevs,nProf)\n
  !>                         unit: [kg/kg]
  !> \param[in] presGrid     pressure grid\n
  !>                         dim: (nLevs,nProf)\n
  !>                         unit: [Pa]
  !> \param[in] hGrid        height grid\n
  !>                         dim: (nLevs,nProf)\n
  !>                         unit: [m]
  !> \param[in] comp         treated atmospheric trace gas component\n
  !>                         dim: (nLevs,nProf)\n
  !>                         unit: [kg/kg]
  !> \param[in] lastValidLev last valid (upper surface) levels of the treated profile\n
  !>                         dim: (nLevs,nProf)\n
  !>                         unit: [kg/kg]
  !> \param[out] col_cont    calculated column content of the atmospheric trace gas component\n
  !>                         dim: (nProf)\n
  !>                         unit: [m]\n
  !> \return 0 if everything was ok - else an error code > 0;
  !-----------------------------------------------------------------------
  subroutine calc_col_cont_kgkg (temp, humi, hGrid, presGrid, comp, lastValidLev, col_cont, ierr)
    real (kind=wp), intent(in)  :: temp      (:,:)    ! (nLevs, nProf)
    real (kind=wp), intent(in)  :: humi      (:,:)    ! (nLevs, nProf)
    real (kind=wp), intent(in)  :: hGrid     (:,:)    ! (nLevs, nProf)
    real (kind=wp), intent(in)  :: presGrid  (:,:)    ! (nLevs, nProf)
    real (kind=wp), intent(in)  :: comp      (:,:)    ! (nLevs, nProf)
    integer,        intent(in)  :: lastValidLev(:)    ! (       nProf)
    real (kind=wp), intent(out) :: col_cont    (:)    ! (       nProf)
    integer,        intent(out) :: ierr
    ! local arguments:
    real (kind=wp) :: comp_kgm3(size(comp,1)  , size(comp,2))
    real (kind=wp) :: layerCont(size(comp,1)-1, size(comp,2))
    integer :: iRun1, iRun2, dimG, nProf, n_p, i_p

    ierr  = NO_ERROR
    dimG  = size (comp,1)
    nProf = size (comp,2)
    n_p   = size (presGrid,2)
    col_cont  (:)   = 0.0_wp
    comp_kgm3 (:,:) = 0.0_wp
    layerCont (:,:) = 0.0_wp
    !--- check input data:
    if ( checkInput ) then
      if ( ( size(hGrid     ,1) /= dimG ) .or. &
           ( size(hGrid     ,2) /= nProf) .or. &
           ( size(temp      ,1) /= dimG ) .or. &
           ( size(temp      ,2) /= nProf) .or. &
           ( size(humi      ,1) /= dimG ) .or. &
           ( size(humi      ,2) /= nProf) .or. &
           ( size(col_cont    ) /= nProf) .or. &
           ( size(lastValidLev) /= nProf) .or. &
           ( size(presGrid  ,1) /= dimG ) ) then
        ierr = ERR_DIM
        return
      endif

      do iRun1 = 1, nProf
        i_p = min(n_p, iRun1)
        if (.not.(all(presGrid(1:dimG-1,i_p  ) < presGrid(2:dimG,i_p  )) .or. &
                  all(presGrid(1:dimG-1,i_p  ) > presGrid(2:dimG,i_p  )) ) .or. &
            .not.(all(hGrid   (1:dimG-1,iRun1) < hGrid   (2:dimG,iRun1)) .or. &
                  all(hGrid   (1:dimG-1,iRun1) > hGrid   (2:dimG,iRun1)) ) ) then
          ierr =  ERR_MONOTON
          return
        endif
      enddo
    endif

    !--- convert component into kg*m^-3:
    !--- (air density = p / (R_DRY * t_virt));  (t_virt = t * (1 + EPS * q));
    !--- with: t=temperature, q=specific humidity, p=pressure;t_virt= virtual temperature,
    !--- R_DRY=Dry air gas constant, R_WV=waterVapour gas constant, EPS=R_WV/R_DRY-1.0_wp;
    do iRun1 = 1, nProf
      i_p = min(n_p, iRun1)
      comp_kgm3 (1:lastValidLev(iRun1),iRun1) = comp (1:lastValidLev (iRun1),iRun1)     * &
                                                presGrid (1:lastValidLev (iRun1),i_p) / &
                                                (R_DRY * (temp (1:lastValidLev (iRun1),iRun1) * &
                                                (1.0_wp + EPS * humi (1:lastValidLev (iRun1),iRun1))))
    enddo

    !--- integrate over whole vertical range of each profile:
    do iRun2 = 1, nProf
      do iRun1 = 1, lastValidLev (iRun2) - 1
        layerCont (iRun1,iRun2) = abs      ( hGrid     (iRun1,iRun2) - hGrid     (iRun1+1,iRun2) ) * &
                                  0.5_wp * ( comp_kgm3 (iRun1,iRun2) + comp_kgm3 (iRun1+1,iRun2) )
      enddo
    enddo

!NEC$ nointerchange
    do iRun2 = 1, nProf
      col_cont (iRun2)    = sum (layerCont (1:lastValidLev(iRun2)-1 ,iRun2))
    enddo

  end subroutine calc_col_cont_kgkg

  subroutine calc_offset(bc, str, offset_old, offset)
    type (t_bcor_coef), intent(in), target :: bc(:)      ! coefficients
    character(len=*),   intent(in)         :: str        ! input string
    real(wp),           intent(in)         :: offset_old ! old offset value
    real(wp),           intent(out)        :: offset     ! offset
    ! Calculates an average over residual biases of selected channels/datasets.
    ! The input string "str" is a semicolon-separated list of instrument/channel/satellite
    ! specifications. Example:
    ! "instr=15 chan=3 ; instr=19 chan=22 ; instr=10 chan=11"
    ! Each instrument/channel/satellite specification (spec. in the following) is a blank separated
    ! list of assignments:
    ! instr=${instr} chan=${chan} [sat=${sat} nosat=${nosat} wgt=${wgt}]
    ! ${sat} and ${nosat} might be comma-separated lists of satellite IDs.
    ! For each spec. an age- and obs.number weighted average of the residual biases is calculated
    ! The different specs. are averaged either with age- and obs.number weights or with weights
    ! specified with "wgt=${wgt}". The wgt-averaging is used if one spec. contains a "wgt=${wgt}"
    ! entry. In this case for all specs. without "wgt=${wgt}" we assume "wgt=1.".
    ! If bias-correction is found that matches with any spec., we downweight the old offset
    ! "offset_old" (which is soted as offset_v in the bias-files).
    ! In the example above the resulting "offset" is an age- and obs.number weighted average of
    ! the residual biases for channel 3 of MHS(instr=15), channel 22 of ATMS(instr=19) and
    ! channel 11 of SSMIS(instr=10).
    ! If one of the specified channels has the "qcnotuse" bit set in the TOVS_OBS_CHAN namelist
    ! (for a particular satellite) it is not used in the offset calculations. Alternatively, bad
    ! channels (for a particular satellite) might be excluded with the nosat=badsat1,badsat2,...
    ! entry or by explicitely definining the good satellites with sat=sat1,sat2,...

    character(len=11),       parameter   :: proc = 'calc_offset'
    type(t_rad_set),         pointer     :: rs             ! bc(i)%i, dataset description
    character(len=len(str))              :: str_
    character(len=len(str)), allocatable :: s(:)
    character(len=len(str)), allocatable :: w(:)
    character(len=10),       allocatable :: ssat(:)
    integer,                 allocatable :: sat(:)
    real(wp)                             :: o(size(bc))
    real(wp)                             :: n(size(bc))
    real(wp)                             :: wgt(size(bc))
    real(wp)                             :: age(size(bc))
    real(wp)                             :: ws(3)
    real(wp)                             :: wgt_, w_age, dt_fc
    real(wp)                             :: age_, o_, n_
    integer                              :: nstr, nw, nsat
    integer                              :: nofs, no
    integer                              :: instr, chan
    integer                              :: ic, ic_       ! channel index
    integer                              :: ib, ib_       ! index in biasc(:)
    integer                              :: i,j,is,isat ,iw         ! loop indices
    integer                              :: stat
    logical                              :: l_tag(5), l_wgt, ldeb

    offset = 0._wp

    if (trim(str) == '') return

    read(str,*,iostat=stat) offset
    if (stat == 0) return

    ldeb = .false.

    ! prepare string array
    str_=tolower(str)
    call split2(str_, sep=';', n=nstr)
    allocate(s(nstr))
    call split2(str_, sep=';', array=s)

    dt_fc  = days(fc_time)
    offset = 0._wp
    nofs   = 0
    no     = 0
    l_wgt  = .false.
    do is = 1, nstr
      ! Interpret string
      call split2(s(is),sep=' =',n=nw)
      allocate(w(nw))
      call split2(s(is),sep=' =',array=w)
      l_tag =.false.
      instr = -1
      chan  = -1
      wgt_  = -1._wp
      nsat  = 0
      do iw = 1, nw
        select case (trim(w(iw)))
        case('instr','inst','instrument')
          l_tag =.false. ; l_tag(1)=.true.
        case('sat','satid')
          l_tag =.false. ; l_tag(2)=.true.
        case('nosat','nosatid','no_sat','no_satid')
          l_tag =.false. ; l_tag(3)=.true.
        case('chan','channel')
          l_tag =.false. ; l_tag(4)=.true.
        case('wgt','weight')
          l_tag =.false. ; l_tag(5)=.true.
        case default
          if (l_tag(1)) then
            read(w(iw),*,iostat=stat) instr
            if (stat/=0) call finish(proc,'invalid instr "'//trim(w(iw))//&
                 '" in offset-string "'//trim(s(is))//'"')
          end if
          if (l_tag(2) .or. l_tag(3)) then
            call split2(w(iw),sep=',',n=nsat)
            allocate(ssat(nsat),sat(nsat))
            call split2(w(iw),sep=',',array=ssat)
            do isat = 1, nsat
              read(ssat(isat),*,iostat=stat) sat(isat)
              if (stat/=0) call finish(proc,'invalid satellite "'//trim(ssat(isat))//&
                   '" in offset-string "'//trim(str)//'"')
              if (any(sat(1:isat-1) == abs(sat(isat)))) &
                   call finish(proc,'double satellite "'//trim(ssat(isat))//&
                   '" in offset-string "'//trim(str)//'"')
              if (l_tag(3)) sat(isat) = -abs(sat(isat))
            end do
          end if
          if (l_tag(4)) then
            read(w(iw),*,iostat=stat) chan
            if (stat/=0) call finish(proc,'invalid chan "'//trim(w(iw))//&
                 '" in offset-string "'//trim(s(is))//'"')
          end if
          if (l_tag(5)) then
            read(w(iw),*,iostat=stat) wgt_
            if (stat/=0) call finish(proc,'invalid weight "'//trim(w(iw))//&
                 '" in offset-string "'//trim(s(is))//'"')
            l_wgt = .true.
          end if
          l_tag = .false.
        end select
      end do
      deallocate(w)
      if (.not.allocated(sat)) allocate(sat(nsat))
      if (instr < 0) call finish(proc,'failed to determine instrument from string "'//trim(s(is))//'"')
      if (chan  < 0) call finish(proc,'failed to determine channel from string "'//trim(s(is))//'"')
      if (ldeb) print*,'offset sat',nsat,sat(1:nsat)
      if (ldeb) print*,'offset instr',instr,chan
      if (ldeb) print*,'offset wgt',wgt_

      ! Evaluate residual/offset
      ic_ = -1
      ic  = -1
      bc_loop: do ib = 1, size(bc)
        rs => bc(ib)%i
        if (any(sat(1:nsat) > 0)) then
          if (.not.any(sat(1:nsat) == rs%satid)) cycle bc_loop
          if (any(-sat(1:nsat) == rs%satid)) cycle bc_loop
        end if
        age_ = days (ana_time - time_cyyyymmddhhmm (bc(ib)% last_update))
        if (ldeb) print*,'offset',rs%satid,rs%grid
        if (ldeb) print*,'offset age',age_,offset_max_age
        if (age_ > offset_max_age) cycle bc_loop
        do i = 1, rs%n_instr
          if (rs%instr(i) == instr) then
            do j = rs%o_ch_i(i)+1,rs%o_ch_i(i)+rs%n_ch_i(i)
              if (chan == rs%chan(j)) then
                if (btest(get_flag(), USE_QCNOTUSE)) cycle bc_loop
                if (rs%instr(i) == rs%grid) then
                  call calc_residual_bias(bc(ib), j, o(nofs+1), n(nofs+1))
                  age(nofs+1) = age_
                  if (n(nofs+1) >= min_entries) then
                    nofs = nofs + 1
                    if (ldeb) print*,'offset bc sat',rs%satid,rs%grid
                    if (ldeb) print*,'offset res',nofs,o(nofs),n(nofs),age_,wgt_
                    ic  = j
                  end if
                else
                  ic_ = j
                  ib_ = ib
                end if
              end if
            end do
          end if
        end do
      end do bc_loop
      if (ic < 0) then
        if (ic_ > 0) then
          call calc_residual_bias(bc(ib_), ic_, o(nofs+1), n(nofs+1))
          age(nofs+1) = days (ana_time - time_cyyyymmddhhmm (bc(ib_)% last_update))
          if (n(nofs+1) >= min_entries) then
            nofs = nofs + 1
            if (ldeb) print*,'offset res_',nofs,o(nofs),n(nofs),age_,wgt_
          end if
          if (dace%lpio) write(0,*) 'WARNING: found only mapped instrument/channel &
               &corresponding to offset string "'//trim(s(is))//'"'
        else
          if (dace%lpio) write(0,*) 'WARNING: failed to find instrument/channel &
               &corresponding to offset string "'//trim(s(is))//'"'
        end if
      end if

      ! Average over data resulting from this part
      if (nofs > no) then
        o_   = 0._wp
        age_ = 0._wp
        n_   = 0._wp
        ws   = 0._wp
        do i = no+1, nofs
          w_age = 1._wp
          if (offset_max_age > 0._wp) &
               w_age = max(min(w_age - (age(i)-dt_fc)/offset_max_age,1._wp),0._wp)
          o_   = o_   + o(i)   * n(i) * w_age ; ws(1) = ws(1) + n(i) * w_age
          age_ = age_ + age(i) * n(i)         ; ws(2) = ws(2) + n(i)
          n_   = n_   + n(i)   * w_age        ; ws(3) = ws(3) + w_age
          if (ldeb) print*,'offset calc',i,o(i),age(i),n(i),w_age,n(i)*w_age
        end do
        if (ws(1) > 0._wp) o  (no+1) = o_   / ws(1)
        if (ws(2) > 0._wp) age(no+1) = age_ / ws(2)
        if (ws(3) > 0._wp) n  (no+1) = n_   !/ ws(3) ! We would like to weight with the total number
        wgt(no+1) = wgt_                             ! obs.number in the final average below
        if (ldeb) print*,'offset calc pre',no+1,o(no+1),age(no+1),n(no+1),wgt(no+1)
        no   = no+1
        nofs = no
      end if

      if (allocated(ssat)) deallocate(ssat)
      deallocate(sat)
    end do

    ! Calculate final average over all parts
    ws(1)   = 0._wp
    if (no > 0) then
      if (l_wgt) then
        where(wgt(1:no) < 0._wp) wgt(1:no) = 1._wp
      else
        wgt(1:no) = n(1:no)
      end if
      offset = 0._wp
      do i = 1, no
        w_age = 1._wp
        if (offset_max_age > 0._wp) &
             w_age = max(min(w_age - (age(i)-dt_fc)/offset_max_age,1._wp),0._wp)
        wgt_ = wgt(i) * w_age
        offset = offset + wgt_ * o(i)
        ws(1)  = ws(1)  + wgt_
        if (ldeb) print*,'offset calc final',i,wgt(i),w_age,wgt_,o(i),offset,ws(1),offset/ws(1)
      end do
    end if
    if (ws(1) > 0._wp) then
      if (ldeb) print*,'offset calc final',offset,ws(1),offset/ws(1)
      offset = offset / ws(1)
    else
      ! Downweight previous value
      ! (might happen if all instruments/channels used for offset calculation were not
      ! available for longer than offset_max_age
      offset = offset_old * exp(-dt_fc/offset_t_decay)
      if (ldeb) print*,'offset calc backup',offset_old,dt_fc,offset_t_decay,exp ( - age_ / offset_t_decay),offset
    end if

  contains

    function get_flag() result(flag)
      integer :: flag
      type(t_rad_set), pointer :: s            ! rad_set dataset description
      integer :: is,i,j
      flag = ibset(0,USE_QCNOTUSE)
      s => null()
      do is = 1, size(rad_set) ! size(rad_set) instead of n_set is crucial, since the
                               ! dataset might be inactive
        if (rad_set(is)%satid == rs%satid .and. rad_set(is)%grid == rs%grid) then
          s => rad_set(is)
          exit
        end if
      end do
      if (.not.associated(s)) return
      do i = 1, s%n_instr
        if (instr == s%instr(i)) then
          do j = s%o_ch_i(i)+1, s%o_ch_i(i)+s%n_ch_i(i)
            if (s%chan(j) == chan) then
              flag = s%flag(j)
              return
            end if
          end do
        end if
      end do
    end function get_flag

  end subroutine calc_offset

  subroutine calc_residual_bias(bc, ic, resid, n)
    type (t_bcor_coef), intent(in), target  :: bc    ! coefficients
    integer,            intent(in)          :: ic    ! channel index
    real(wp),           intent(out)         :: resid
    real(wp),           intent(out)         :: n

    character(len=18), parameter :: proc = 'calc_residual_bias'
    type(t_stats_data), pointer  :: st
    integer                      :: ifov, ipred
    real(wp)                     :: d_fov, b_fov

    resid = 0._wp
    n     = -1._wp
    st => bc%sta
    if (.not.st%alloc) return
    n = sum(st%btdn(:,ic))
    if (n > 0) then
      do ifov = 1, bc%i%n_fov
        if (st%btdn(ifov,ic) > 0._wp) then
          ! obs-fd
          d_fov = st%btdscan(ifov,ic)
          b_fov = 0._wp
          ! pred. biascorr.
          do ipred = 1, st%n_pred
            if (bc%pred_use(ipred, ic) == 1) then
              b_fov = b_fov + bc%coef_set1(ipred,ic) * st%pred_i(ipred,ifov,ic)
            end if
          end do
          if (bc%fov_bc(ic)%const_pred) b_fov = b_fov + bc%coef_set1(st%n_pred+1,ic) * st%btdn(ifov,ic)
          ! FOV biascorr.
          if (bc%fov_bcorr) b_fov = b_fov + bc%coef_set2(ic,ifov) * st%btdn(ifov,ic)
          resid = resid + (d_fov - b_fov)
        end if
      end do
      resid = resid / real(n,wp)
    end if

  end subroutine calc_residual_bias


  subroutine diagn_biasc (biasc, statistics, timeseries)
  type (t_bcor_coef)  ,intent(in) :: biasc (:)  ! coefficients
  character           ,intent(in) :: statistics ! statistics to use 'a' or 'i'
  logical             ,intent(in) :: timeseries ! create files for timeseries
  !---------------------------------------------
  ! print files for bias coefficient diagnostics
  ! to files to be interpreded by gnuplot
  !---------------------------------------------
    integer :: i
    do i=1,size(biasc)
      call diagn_biasc_1 (biasc(i), statistics, timeseries)
    end do
  end subroutine diagn_biasc

!------------------------------------------------------------------------------

  subroutine diagn_biasc_1 (biasc, statistics, timeseries)
  type (t_bcor_coef)  ,intent(in) :: biasc      ! coefficients
  character           ,intent(in) :: statistics ! statistics to use 'a' or 'i'
  logical             ,intent(in) :: timeseries ! create files for timeseries
  target                          :: biasc
  !---------------------------------------------
  ! print files for bias coefficient diagnostics
  ! to files to be interpreted by gnuplot
  !---------------------------------------------

    character(len=3)  :: csatid    ! satellite  id string
    character(len=3)  :: cinstr    ! instrument id string
    character(len=4)  :: cchan     ! channel    id string
    character(len=2)  :: cgrid     ! grid       id string
    character(len=64) :: file      ! output diagnostics file name
    integer           :: unit = -1 ! fortran unit number
    integer           :: i, j, k   ! indices
    integer           :: kchan     ! channel index variable
    integer           :: ichan     ! channel id
    integer           :: instr     ! instrument id
    integer           :: igrid     ! grid id
    real(wp)          :: f         ! factor (1 or 2)

    real(wp)          :: nn                        ! number of data points
    real(wp)          :: nn1                       !
    integer           :: np                        ! number of predictors
    real(wp)          :: gmcorr                    ! mean global correction
    real(wp)          :: corr   (biasc% i% n_fov)  ! fov debendent correction
    real(wp)          :: n      (biasc% i% n_fov)  ! # of data points   / fov
    real(wp)          :: d      (biasc% i% n_fov)  ! mean deviation     / fov
    real(wp)          :: rmse   (biasc% i% n_fov)  ! rmse of deviation  / fov
    real(wp)          :: gcorr  (biasc% i% n_fov)  ! global correction  / fov
    real(wp)          :: stdev  (biasc% i% n_fov)  ! stdev of deviation / fov
    real(wp)          :: crmse  (biasc% i% n_fov)  ! rmse after bias correction
    logical           :: msk    (biasc% i% n_fov)  ! mask for fov valid

    integer           :: ix     (biasc%    n_pred) ! matrix diagonal indices
    real(wp)          :: p_stdev(biasc%    n_pred,&! predictor stdev
                                 biasc% i% n_chan )!
    real(wp)          :: p_mean (biasc%    n_pred,&! predictor mean
                                 biasc% i% n_chan )!
    real(wp)          :: p_n    (biasc%    n_pred,&! number of entries
                                 biasc% i% n_chan )!

    real(wp)          :: g_n    (biasc% i% n_chan )! global number of data
    real(wp)          :: g_mean (biasc% i% n_chan )! global mean
    real(wp)          :: g_corr (biasc% i% n_chan )! global mean correction
    real(wp)          :: g_rmse (biasc% i% n_chan )! global rmse
    real(wp)          :: g_stdev(biasc% i% n_chan )! global stdev
    real(wp)          :: g_crmse(biasc% i% n_chan )! global corrected rmse

    real(wp)          :: beta_n (biasc%    n_pred,&! normalised coefficients
                                 biasc% i% n_chan )!
    logical           :: exist                     ! existence of file

    type(t_stats_data) ,pointer :: stat            ! statistics to use

    !----------------------------------------------
    ! chose accumulated or instantaneous statistics
    !----------------------------------------------
    select case (statistics)
    case ('a')
      stat => biasc% sta
    case ('i')
      stat => biasc% sti
    case default
      call finish('diagn_biasc',                                    &
                  'invalid value for statistics (a,i): '//statistics)
    end select

    !=======================
    ! scan angle diagnostics
    !=======================

    !-------------------
    ! loop over channels
    !-------------------
    do kchan = 1, biasc% i% n_chan               ! channel index
      ichan = biasc% i% chan (kchan)             ! channel id
      instr = instr_chan     (kchan, biasc% i)   ! instrument id
      igrid = biasc% i% grid

      !--------------------------------
      ! set some quantities in any case
      !--------------------------------
      nn      = 0._wp                      ! total number of data points
      nn1     = 1._wp                      !
      n       = 0._wp                      ! number of data points / fov
      d       = 0._wp                      ! sum of deviations     / fov
      rmse    = 0._wp                      ! rmse                  / fov
      gmcorr  = 0._wp                      ! correction by global coefficients
      np      = biasc% n_pred              ! number of predictors
      corr    = biasc% coef_set2 (kchan,:) ! fov debendent correction
      !-------------------------------------------
      ! set quantities if statistics are available
      !-------------------------------------------
      if (associated (stat% btdn)) then
        n    = stat% btdn (:,kchan)
        msk  = n > 0._wp
        nn   = sum (n, msk)
        nn1  = nn ;if (nn == 0._wp) nn1 = 1._wp
        where (.not. msk)
          !-----------------------------------------------------------------
          ! set invalid values to be printed as ***** for unused scan angles
          !-----------------------------------------------------------------
          n     = -huge (n)
          corr  = -huge (corr)
          d     = -huge (d)
          rmse  = -huge (rmse)
          stdev = -huge (stdev)
          gcorr = -huge (gcorr)
          crmse = -huge (crmse)
        elsewhere
          !----------------------------------
          ! valid values for used scan angles
          !----------------------------------
          d     =       stat% btdscan (:,kchan)  / n
          rmse  = sqrt (stat% d2      (:,kchan)  / n)
          stdev = sqrt (max(0._wp,rmse*rmse - d*d))
          gcorr = biasc% coef_set1 (np+1, kchan)
          crmse = stat% d2   (:,kchan)
          crmse = crmse + n * ( biasc% coef_set1 (np+1, kchan  ) &
                               +biasc% coef_set2 (      kchan,:) ) ** 2
          crmse = crmse - 2 * ( biasc% coef_set1 (np+1, kchan  ) &
                               +biasc% coef_set2 (      kchan,:))&
                            *   stat%  btdscan     (:  ,kchan  )
        endwhere
        !-----------------------------------------------------
        ! finish calculation of global correction for the mean
        !-----------------------------------------------------
        do i = 1, np
          where (msk)
            gcorr = gcorr + biasc% coef_set1 (i,  kchan) &
                          * stat%  pred_i    (i,:,kchan) / n
          endwhere
        end do
        !--------------------------------------
        ! finish global correction for the rmse
        !--------------------------------------
        if (nn > 0) then
          gmcorr = sum (biasc% coef_set1 (np+1, kchan) * n, msk)
          k = 1
          do i = 1, np
            ix (i) = k
            gmcorr = gmcorr +      biasc% coef_set1 (i,  kchan)     &
                            * sum (stat%  pred_i    (i,:,kchan), msk)
            where (n > 0)
              crmse = crmse - 2 *   stat% btd_pred_i   (i,:,kchan)  &
                                *   biasc% coef_set1   (i,  kchan)
              crmse = crmse + 2 * ( biasc% coef_set1 (np+1, kchan  ) &
                                   +biasc% coef_set2 (      kchan,:))&
                                *   biasc% coef_set1 (i   , kchan  ) &
                                *   stat%  pred_i    (i ,:, kchan  )
            endwhere
            do j = i, np
              f = 2._wp; if (i==j) f = 1._wp
              where (n > 0)
                crmse = crmse + f * stat%  pred_ij   (k,:,kchan) &
                                  * biasc% coef_set1 (i,  kchan) &
                                  * biasc% coef_set1 (j,  kchan)
              endwhere
              k = k + 1
            end do
          end do
          gmcorr = gmcorr / nn
           where (n > 0)
             crmse = sqrt ( MAX(crmse,0._wp) / n)
           endwhere
        endif
      endif

      !-----------------------------------------------
      ! write scan angle dependent diagnostics to file
      !-----------------------------------------------
      write (csatid,'(i3.3)') biasc% i% satid
      write (cchan,  '(i4.4)') ichan
      write (cinstr, '(i3.3)') instr
      write (cgrid,  '(i2.2)') igrid
      if (unit < 0) call get_free_unit (unit)
      file = 'biasc_scan_'//csatid//'_'//cinstr//'_'//cgrid//' '//cchan//'.gp'
      open (unit, file=file, action='write')
        write(unit,'(a)')   '#'
        write(unit,'(a)')   '# Scan Angle Dependent Bias Correction Statistics'
        write(unit,'(a)')   '#'
        write(unit,'(a,i5)')'# experiment =' ,biasc% exp
        write(unit,'(a,a)') '# file       = ',trim(file)
        write(unit,'(a,i5)')'# satid      =' ,biasc% i% satid
        write(unit,'(a,i5)')'# channel    =' ,ichan
        write(unit,'(a,i5)')'# instrument =' ,instr
        write(unit,'(a,a)') '# statistics = ',statistics
        write(unit,'(a)')   '#'
        write(unit,'(a)')   '# column 1   : field of view index'
        write(unit,'(a)')   '# column 2   : number of entries'
        write(unit,'(a)')   '# column 3   : fov         bias correction'
        write(unit,'(a)')   '# column 4   : global      bias correction'
        write(unit,'(a)')   '# column 5   : mean global bias correction'
        write(unit,'(a)')   '# column 6   : mean deviation'
        write(unit,'(a)')   '# column 7   : rmse'
        write(unit,'(a)')   '# column 8   : stdev'
        write(unit,'(a)')   '# column 9   : rmse after bias correction'
        write(unit,'(a)')   '#'
        do i = 1, biasc% i% n_fov
          write (unit,'(i4,1x,f10.2,7(1x,f8.3))') &
            i,                                    &
            n     (i),                            &
            corr  (i),                            &
            gcorr (i),                            &
            gmcorr,                               &
            d     (i),                            &
            rmse  (i),                            &
            stdev (i),                            &
            crmse (i)
        end do
      close (unit)

      !------------------------------------
      ! keep results for global diagnostics
      !------------------------------------
      g_n    (kchan) = nn
      if (nn > 0._wp) then
        g_mean (kchan) =       sum (d        * n, msk) / nn
        g_corr (kchan) = gmcorr
        g_rmse (kchan) = sqrt (sum (rmse **2 * n, msk) / nn)
        g_stdev(kchan) = sqrt (sum (stdev**2 * n, msk) / nn)
        g_crmse(kchan) = sqrt (sum (crmse**2 * n, msk) / nn)
      else
        g_mean (kchan) = -huge (g_mean )
        g_corr (kchan) = -huge (g_corr )
        g_rmse (kchan) = -huge (g_rmse )
        g_stdev(kchan) = -huge (g_stdev)
        g_crmse(kchan) = -huge (g_crmse)
      endif

      if (timeseries) then
        !-----------------------------------------------------------------
        ! append to time-series file: Scan Angle Dependent Bias Correction
        !-----------------------------------------------------------------
        file = 'biasc_scan_'//csatid//'_'//cinstr//'_'//cgrid//'_'//cchan//'_scancor.gp'
        inquire (file=file, exist=exist)
        open (unit, file=file, action='write', position='append')
        if (.not.exist) then
          write(unit,'(a)')   '#'
          write(unit,'(a)')   '# Scan Angle Dependent Bias Correction'
          write(unit,'(a)')   '#'
          write(unit,'(a,i5)')'# experiment     =' ,biasc% exp
          write(unit,'(a,a)') '# file           = ',trim(file)
          write(unit,'(a,i5)')'# satid          =' ,biasc% i% satid
          write(unit,'(a,i5)')'# channel        =' ,ichan
          write(unit,'(a,i5)')'# instrument     =' ,instr
          write(unit,'(a,i5)')'# number of FOVs =' ,biasc% i% n_fov
          write(unit,'(a,a)') '# statistics     = ',statistics
          write(unit,'(a)')   '#'
          write(unit,'(a)')   '# column 1   : time'
          write(unit,'(a)')   '# column n   : FOV bias correction'
          write(unit,'(a)')   '#'
        endif
        write(unit,'(a,100(1x,f8.3))') biasc% last_date, corr(:)
        close (unit)

        !----------------------------------------------------------------
        ! append to time-series file: Predictor Dependent Bias Correction
        !----------------------------------------------------------------
        file = 'biasc_scan_'//csatid//'_'//cinstr//'_'//cgrid//'_'//cchan//'_predcor.gp'
        inquire (file=file, exist=exist)
        open (unit, file=file, action='write', position='append')
        if (.not.exist) then
          write(unit,'(a)')   '#'
          write(unit,'(a)')   '# Predictor Dependent Bias Correction'
          write(unit,'(a)')   '#'
          write(unit,'(a,i5)')'# experiment     =' ,biasc% exp
          write(unit,'(a,a)') '# file           = ',trim(file)
          write(unit,'(a,i5)')'# satid          =' ,biasc% i% satid
          write(unit,'(a,i5)')'# channel        =' ,ichan
          write(unit,'(a,i5)')'# instrument     =' ,instr
          write(unit,'(a,i5)')'# number of FOVs =' ,biasc% i% n_fov
          write(unit,'(a,a)') '# statistics     = ',statistics
          write(unit,'(a)')   '#'
          write(unit,'(a)')   '# column 1   : time'
          write(unit,'(a)')   '# column n   : FOV bias correction'
          write(unit,'(a)')   '#'
        endif
        write(unit,'(a,100(1x,f8.3))') biasc% last_date, gcorr(:)
        close (unit)

        !------------------------------------------------------------
        ! append to time-series file: Mean Deviation after Correction
        !------------------------------------------------------------
        file = 'biasc_scan_'//csatid//'_'//cinstr//'_'//cgrid//'_'//cchan//'_resid.gp'
        inquire (file=file, exist=exist)
        open (unit, file=file, action='write', position='append')
        if (.not.exist) then
          write(unit,'(a)')   '#'
          write(unit,'(a)')   '# Mean Deviation after Correction'
          write(unit,'(a)')   '#'
          write(unit,'(a,i5)')'# experiment     =' ,biasc% exp
          write(unit,'(a,a)') '# file           = ',trim(file)
          write(unit,'(a,i5)')'# satid          =' ,biasc% i% satid
          write(unit,'(a,i5)')'# channel        =' ,ichan
          write(unit,'(a,i5)')'# instrument     =' ,instr
          write(unit,'(a,i5)')'# number of FOVs =' ,biasc% i% n_fov
          write(unit,'(a,a)') '# statistics     = ',statistics
          write(unit,'(a)')   '#'
          write(unit,'(a)')   '# column 1   : time'
          write(unit,'(a)')   '# column n   : FOV bias'
          write(unit,'(a)')   '#'
        endif
        write(unit,'(a,100(1x,f8.3))') biasc% last_date, d - corr - gcorr
        close (unit)

        !--------------------------------------------------
        ! append to time-series file: RMSE after Correction
        !--------------------------------------------------
        file = 'biasc_scan_'//csatid//'_'//cinstr//'_'//cgrid//'_'//cchan//'_rmse.gp'
        inquire (file=file, exist=exist)
        open (unit, file=file, action='write', position='append')
        if (.not.exist) then
          write(unit,'(a)')   '#'
          write(unit,'(a)')   '# RMSE after Correction'
          write(unit,'(a)')   '#'
          write(unit,'(a,i5)')'# experiment     =' ,biasc% exp
          write(unit,'(a,a)') '# file           = ',trim(file)
          write(unit,'(a,i5)')'# satid          =' ,biasc% i% satid
          write(unit,'(a,i5)')'# channel        =' ,ichan
          write(unit,'(a,i5)')'# instrument     =' ,instr
          write(unit,'(a,i5)')'# number of FOVs =' ,biasc% i% n_fov
          write(unit,'(a,a)') '# statistics     = ',statistics
          write(unit,'(a)')   '#'
          write(unit,'(a)')   '# column 1   : time'
          write(unit,'(a)')   '# column n   : FOV RMSE'
          write(unit,'(a)')   '#'
        endif
        write(unit,'(a,100(1x,f8.3))') biasc% last_date, crmse(:)
        close (unit)

        !-----------------------------------------------
        ! append to time-series file: Average Statistics
        !-----------------------------------------------
        file = 'biasc_scan_'//csatid//'_'//cinstr//'_'//cgrid//'_'//cchan//'_ave.gp'
        inquire (file=file, exist=exist)
        open (unit, file=file, action='write', position='append')
        if (.not.exist) then
          write(unit,'(a)')   '#'
          write(unit,'(a)')   '# Average Statistics'
          write(unit,'(a)')   '#'
          write(unit,'(a,i5)')'# experiment =' ,biasc% exp
          write(unit,'(a,a)') '# file       = ',trim(file)
          write(unit,'(a,i5)')'# satid      =' ,biasc% i% satid
          write(unit,'(a,i5)')'# channel    =' ,ichan
          write(unit,'(a,i5)')'# instrument =' ,instr
          write(unit,'(a,a)') '# statistics = ',statistics
          write(unit,'(a)')   '#'
          write(unit,'(a)')   '# column 1   : time'
          write(unit,'(a)')   '# column 2   : entries'
          write(unit,'(a)')   '# column 3   : obs-fg'
          write(unit,'(a)')   '# column 4   : sqrt((obs-fg)^2)'
          write(unit,'(a)')   '# column 5   : stdev'
          write(unit,'(a)')   '# column 6   : bcor'
          write(unit,'(a)')   '# column 7   : (obs-fg)-bcor'
          write(unit,'(a)')   '# column 8   : sqrt((obs-fg-bcor)^2)'
          write(unit,'(a)')   '#'
        endif
        write(unit,'(a,f12.2,6(1x,f8.3))') biasc% last_date, nn, &
          g_mean(kchan), g_rmse(kchan), g_stdev(kchan), gmcorr,  &
          g_mean(kchan)-gmcorr, g_crmse(kchan)
        close (unit)

      endif

    end do

    !===============================
    ! global coefficient diagnostics
    !===============================
    if (any (g_n > 0._wp)) then

      !-------------------------------
      ! calculate predictor bias, rmse
      !-------------------------------
      p_n  = spread (g_n,1,np)
      where (p_n > 0._wp .and. biasc% pred_use /= 0)
        p_mean  =       sum (stat% pred_i,         dim=2) / p_n
        p_stdev = sqrt (sum (stat% pred_ij(ix,:,:),dim=2) / p_n &
                        - p_mean ** 2                           )
      elsewhere
        p_mean  = -huge(p_mean)
        p_stdev = -huge(p_stdev)
      endwhere

      !----------------------
      ! print mean predictors
      !----------------------
      file = 'biasc_pred_mean_'//csatid//'_'//'.gp'
      open (unit, file=file, action='write')
        write(unit,'(a)')   '#'
        write(unit,'(a)')   '# predictor mean'
        write(unit,'(a)')   '#'
        write(unit,'(a,i5)')'# experiment =' ,biasc% exp
        write(unit,'(a,a)') '# file       = ',trim(file)
        write(unit,'(a,i5)')'# satid      =' ,biasc% i% satid
        write(unit,'(a,i5)')'# channels   =' ,biasc% i% n_chan
        write(unit,'(a,i5)')'# predictors =' ,np
        write(unit,'(a,a)') '# statistics = ',statistics
        write(unit,'(a)')   '#'
        write(unit,'(a)')   '# column  1  : index'
        write(unit,'(a)')   '# column  2  : instrument id'
        write(unit,'(a)')   '# column  3  : channel    id'
        write(unit,'(a)')   '# column  4  : number of entries'
        do i = 1, np
          write(unit,"(a,i2,a,i2,' : ',a,2i8)")                &
                            '# column ',i+4,'  : mean pred.',i,&
            biasc% pred(i), biasc% pred_p1(i), biasc% pred_p2(i)
        end do
        write(unit,'(a)')   '#'
        do kchan = 1, biasc% i% n_chan              ! channel index
          ichan = biasc% i% chan (kchan)            ! channel id
          instr = instr_chan     (kchan, biasc% i)  ! instrument id
          write (unit,'(3i5,f10.0,100(1x,e12.3))') &
            kchan, instr, ichan,                  & ! index, instrument, channel
            g_n    (  kchan),                     & ! number of entries
            p_mean (:,kchan)                        ! predictor mean value
        end do
      close (unit)

      !-------------------------
      ! print rmse of predictors
      !-------------------------
      file = 'biasc_pred_stdev_'//csatid//'.gp'
      open (unit, file=file, action='write')
        write(unit,'(a)')   '#'
        write(unit,'(a)')   '# predictor stdev'
        write(unit,'(a)')   '#'
        write(unit,'(a,i5)')'# experiment =' ,biasc% exp
        write(unit,'(a,a)') '# file       = ',trim(file)
        write(unit,'(a,i5)')'# satid      =' ,biasc% i% satid
        write(unit,'(a,i5)')'# channels   =' ,biasc% i% n_chan
        write(unit,'(a,i5)')'# predictors =' ,np
        write(unit,'(a,a)') '# statistics = ',statistics
        write(unit,'(a)')   '#'
        write(unit,'(a)')   '# column  1  : index'
        write(unit,'(a)')   '# column  2  : instrument id'
        write(unit,'(a)')   '# column  3  : channel    id'
        write(unit,'(a)')   '# column  4  : number of entries'
        do i = 1, np
          write(unit,"(a,i2,a,i2,' : ',a,2i8)")                &
                            '# column ',i+4,'  : stdev pred.',i,&
            biasc% pred(i), biasc% pred_p1(i), biasc% pred_p2(i)
        end do
        write(unit,'(a)')   '#'
        do kchan = 1, biasc% i% n_chan              ! channel index
          ichan = biasc% i% chan (kchan)            ! channel id
          instr = instr_chan     (kchan, biasc% i)  ! instrument id
          write (unit,'(3i5,f10.0,100(1x,f8.3))') &
            kchan, instr, ichan,                  & ! index, instrument, channel
            g_n     (  kchan),                    & ! number of entries
            p_stdev (:,kchan)                       ! predictor stdev
        end do
      close (unit)

      !-----------------------------------
      ! print bias correction coefficients
      !-----------------------------------

      where (p_n > 0._wp)
        beta_n = biasc% coef_set1 (1:np,:) * p_stdev
      elsewhere
        beta_n = -huge(beta_n)
      endwhere

      file = 'biasc_coef_p_stdev_'//csatid//'.gp'
      open (unit, file=file, action='write')
        write(unit,'(a)')   '#'
        write(unit,'(a)')   '# bias correction coefficients'
        write(unit,'(a)')   '# multiplied by predictor stdev'
        write(unit,'(a)')   '#'
        write(unit,'(a,i5)')'# experiment =' ,biasc% exp
        write(unit,'(a,a)') '# file       = ',trim(file)
        write(unit,'(a,i5)')'# satid      =' ,biasc% i% satid
        write(unit,'(a,i5)')'# channels   =' ,biasc% i% n_chan
        write(unit,'(a,i5)')'# predictors =' ,np
        write(unit,'(a,a)') '# statistics = ',statistics
        write(unit,'(a)')   '#'
        write(unit,'(a)')   '# column  1  : index'
        write(unit,'(a)')   '# column  2  : instrument id'
        write(unit,'(a)')   '# column  3  : channel    id'
        write(unit,'(a)')   '# column  4  : bt mean deviation'
        write(unit,'(a)')   '# column  5  : bt mean correction'
        write(unit,'(a)')   '# column  6  : bt rms deviation'
        write(unit,'(a)')   '# column  7  : bt stdev'
        write(unit,'(a)')   '# column  8  : corrected bt rmse'
        do i = 1, np
          write(unit,"(a,i2,a,i2,' : ',a,2i8)")                &
                            '# column ',i+8,'  : coefficient ',i,&
            biasc% pred(i), biasc% pred_p1(i), biasc% pred_p2(i)
        end do
        write(unit,'(a)')   '#'
        do k = 1, biasc% i% n_chan              ! channel index
          ichan = biasc% i% chan (k)            ! channel id
          instr = instr_chan     (k, biasc% i)  ! instrument id
          write (unit,'(3i5,5(1x,f7.3),100(1x,f8.3))') &
            k, instr, ichan,                  & ! index, instrument, channel
            g_mean (k), g_corr (k),           & ! mean deviation,correction
            g_rmse (k),g_stdev(k),g_crmse(k), & ! rms, stdev, corrected rms
            beta_n (:,k)                        ! coefficients
        end do
      close (unit)

      if (timeseries) then
        !--------------------------------------------------
        ! append to time-series file: Normalised Predictors
        !--------------------------------------------------
        do k = 1, biasc% i% n_chan              ! channel index
          ichan = biasc% i% chan (k)            ! channel id
          instr = instr_chan     (k, biasc% i)  ! instrument id
          igrid = biasc% i% grid
          write (cchan,  '(i4.4)') ichan
          write (cinstr, '(i3.3)') instr
          write (cgrid,  '(i2.2)') igrid

          file = 'biasc_scan_'//csatid//'_'//cinstr//'_'//cgrid//'_'//cchan//'_nbccof.gp'
          inquire (file=file, exist=exist)
          open (unit, file=file, action='write', position='append')
          if (.not.exist) then
            write(unit,'(a)')   '#'
            write(unit,'(a)')   '# Normalised Coefficients'
            write(unit,'(a)')   '#'
            write(unit,'(a,i5)')'# experiment           =' ,biasc% exp
            write(unit,'(a,a)') '# file                 = ',trim(file)
            write(unit,'(a,i5)')'# satid                =' ,biasc% i% satid
            write(unit,'(a,i5)')'# channel              =' ,ichan
            write(unit,'(a,i5)')'# instrument           =' ,instr
            write(unit,'(a,i5)')'# number of predictors =' ,biasc% n_pred
            write(unit,'(a,a)') '# statistics           = ',statistics
            write(unit,'(a)')   '#'
            write(unit,'(a)')   '# column 1   : time'
            write(unit,'(a)')   '# column 2   : obs-fg RMSE'
            write(unit,'(a)')   '# column n   : normalised coefficient'
            write(unit,'(a)')   '#'
            write(unit,'(a,100(1x,i8))')'# usage      :  ', biasc% pred_use(:,k)
            write(unit,'(a)')   '#'
          endif
          write(unit,'(a,100(1x,f8.4))') biasc% last_date,g_stdev(k),beta_n(:,k)
          close (unit)
        end do
      endif

    endif

  end subroutine diagn_biasc_1

!==============================================================================

  subroutine diagn_pred (biasc, use_all, keep_used, start_all, n_sel)
  type (t_bcor_coef) ,intent(in) :: biasc (:)  ! biascor meta data + statistics
  logical            ,intent(in) :: use_all    ! use all predictors
  logical            ,intent(in) :: keep_used  ! prevent active predictors from being deleted in the selection process
  logical            ,intent(in) :: start_all  ! start selection process with all pred.
  integer            ,intent(in) :: n_sel      ! brute force selection process: select n_select predictors.
  target                         :: biasc
  !-----------------------------------------
  ! Diagnose impact on individual predictors
  !-----------------------------------------
    integer                     :: is, ic, i, j, k, m, n, instr
    integer                     :: slot
    type (t_bcor_coef) ,pointer :: bc
    integer                     :: n_pred
    integer                     :: n_chan
    integer                     :: n_fov
    integer                     :: n_instr
    integer                     :: m_active
    integer                     :: m_store, ipr, ic0, n_chan_aux
    integer ,allocatable        :: n_active (:)
    integer ,allocatable        :: n_used   (:)
    integer ,allocatable        :: p_switch    (:)
    integer ,allocatable        :: p_ls     (:,:)
    integer ,allocatable        :: ind(:)
    integer ,allocatable        :: n_store(:)
    real(wp),allocatable        :: p_rmse   (:,:)
    integer ,allocatable        :: p_ls2    (:,:)
    real(wp),allocatable        :: p_rat2  (:,:)
    integer ,allocatable        :: pred_use (:,:)
    integer ,allocatable        :: pred   (:,:,:)
    integer ,allocatable        :: pred_ind(:)
    real(wp) ,allocatable       :: rmse   (:,:,:)
    real(wp) ,allocatable       :: tmp    (:)
    logical  ,allocatable       :: mask   (:)
    logical  ,allocatable       :: pset   (:,:)
    real(wp)                    :: rmse_aux
    logical                     :: l_found
    character(len=80)           :: pname
    character(len=64)           :: file                ! diagnostics file name
    character(len=3)            :: csatid              ! satellite id string
    character(len=2)            :: cgrid               ! grid      id string
    character(len=4)            :: cchan               ! channel   id string
    character(len=1)            :: cs                  ! method string
    character(len=1)            :: cp                  ! #predictors string
    integer                     :: unit = -1           ! fortran unit number

    integer :: icpr = -1

    type t_sum_sat
      integer              :: grid      = -1
      integer              :: n_chan1   = 0
      integer, allocatable :: chan      (:)
      integer, allocatable :: satellites(:)
      integer, allocatable :: pred_use  (:,:)
      integer, allocatable :: pred      (:,:,:)
      integer, allocatable :: n         (:,:)
      real(wp),allocatable :: rmse      (:,:)
      integer, allocatable :: n_entry   (:)
    end type t_sum_sat
    integer,         parameter :: m_grid = 5
    type(t_sum_sat), target    :: sum_arr(m_grid)
    type(t_sum_sat), pointer   :: s => null()

    !-----------------------------------
    ! loop over satellites / instruments
    !-----------------------------------
    do is=1,size(biasc)
      bc       => biasc(is)
      if (bc%entries < 2.) cycle
      n_pred   = bc%    n_pred
      n_chan   = bc% i% n_chan
      n_instr  = bc% i% n_instr
      n_fov    = bc% i% n_fov
      write(6,*)
      write(6,'(a,i3,a,i4,a,i4,a,10i4)') 'process set',is, &
           ': satid =', bc% i% satid,       &
           ', grid  =', bc% i% grid,         &
           ', instr =', bc% i% instr(:n_instr)
      write (csatid,'(i3.3)') bc% i% satid
      write (cgrid ,'(i2.2)') bc% i% grid

      !-----------------------------
      ! allocate arrays for this set
      !-----------------------------
      allocate (pred_use (n_pred,           n_chan))
      allocate (rmse     (n_pred+2, n_pred, n_chan))
      allocate (n_active (                  n_chan))
      allocate (n_used   (                  n_chan))
      allocate (p_switch (                  n_chan))
      allocate (tmp      (                  n_chan))
      allocate (p_rmse   (0:n_pred,         n_chan)); p_rmse = 0
      allocate (p_ls     (0:n_pred,         n_chan)); p_ls   = 0
      allocate (p_rat2   (0:n_pred,         n_chan)); p_rat2= 0
      allocate (p_ls2    (0:n_pred,         n_chan)); p_ls2  = 0
      allocate (mask     (n_pred                  ))

      !-----------------------------------------
      ! keep old 'pred_use' flags, set temporary
      !-----------------------------------------
      pred_use = bc% pred_use
      if (use_all) then
        where (bc% pred_use == -1) bc% pred_use = 1
      endif
      where   (bc% pred_use == -1) bc% pred_use = 0

      !--------------------------------------------------
      ! calculate bias correction for constant terms only
      !--------------------------------------------------
      write(6,*) 'calculate bias correction for constant terms only'
      rmse = 0._wp
      slot = 1
      p_rmse (0,:)      = scbc_stdev()
      rmse (slot, :, :) = spread (p_rmse (0,:),1,n_pred)

!     !------------------------------------------------------
!     ! calculate bias correction for empty set of predictors
!     !------------------------------------------------------
!     write(6,*) 'calculate bias correction for empty set of predictors'
!     where (bc% pred_use /= 0) bc% pred_use = -1
!     call calc_bcor_coefs (bc% sta, bc)
!     slot = 1
!     rmse (slot, :,:) = spread (scbc_rmse (),1,n_pred)

      !-----------------------------------------------------
      ! calculate bias correction for full set of predictors
      !-----------------------------------------------------
      write(6,*) 'calculate bias correction for full set of predictors'
      where (bc% pred_use /= 0) bc% pred_use = 1
      call calc_bcor_coefs (bc% sta, bc)
      slot = 2
      rmse (slot, :,:) = spread (scbc_rmse (),1,n_pred)
      n_active = count  (bc% pred_use == 1, dim=1)
      n_used   = n_active
      do ic = 1, n_chan
        if (rmse (slot,1,ic) == 0) n_used(ic) = 0
        p_rmse (n_used(ic),ic) = rmse (slot,1,ic)
      end do
      !----------------------------------------
      ! fill in rmse for all and none predictor
      !----------------------------------------

      p_rat2 = p_rmse

      if (n_sel <= 0) then
        if (start_all) then
          cs = 'a'
          !--------------------------------------
          ! loop over number of active predictors
          !--------------------------------------
          do
            n_active = count  (bc% pred_use == 1, dim=1)

            m_active = maxval (n_active)
            if (m_active < 1) exit
            slot = slot + 1

            !---------------------------------
            ! switch off one predictor in turn
            !---------------------------------
            write(6,*) 'switch off one predictor in turn out of',m_active
            p_switch    = 0
            do
              do ic = 1, n_chan
                if (p_switch (ic) > n_pred) cycle
                if (p_switch (ic) > 0     ) bc% pred_use (p_switch (ic), ic) = 1
                do i = p_switch (ic) + 1, n_pred
                  if (bc% pred_use (i,ic) == 1) then
                    bc% pred_use (i,ic) = -1
                    p_switch (ic) = i
                    exit
                  endif
                end do
                if (.not. any (bc% pred_use (:,ic) == -1)) then
                  p_switch(ic) = n_pred + 1
                  where (bc% pred_use (:,ic) == 1) bc% pred_use (:,ic) = -1
                endif
              end do
              if (all(p_switch > n_pred)) exit

              !--------------------------------------------------------
              ! calculate bias correction for reduced set of predictors
              !--------------------------------------------------------
              call calc_bcor_coefs (bc% sta, bc)
              tmp = scbc_rmse ()
              do ic = 1, n_chan
                if (bc% sta% n(ic) < 2) cycle
                if (p_switch(ic) <= n_pred) rmse (slot, p_switch(ic), ic) = tmp (ic)
              end do
            end do

            !------------------------------------
            ! determine least important predictor
            !------------------------------------
            do ic = 1, n_chan
              if (n_active(ic) > 0) then
                mask(:) = (bc% pred_use (:,ic) /= 0)
                if (keep_used) then
                  where(pred_use(:,ic) == 1) mask(:) = .false.
                  if (.not.any(mask)) mask(:) = (bc% pred_use (:,ic) /= 0)
                end if
                i = sum (minloc (rmse (slot,:,ic), mask = mask(:)))
                bc% pred_use (i,ic) = 0
                p_ls   (n_active(ic)  ,ic) = i
                p_rmse (n_active(ic)-1,ic) = rmse (slot,i,ic)
                if (bc% sta% n(ic) >= 2) then
                  mask(i) = .false.
                  if (any(mask(:))) i = sum(minloc (rmse (slot,:,ic), mask = mask(:)))
                  p_ls2  (n_active(ic)  ,ic) = i
                  p_rat2 (n_active(ic)  ,ic) = (rmse (slot,i,ic)          - p_rmse(n_active(ic),ic)) / &
                       (p_rmse(n_active(ic)-1,ic) - p_rmse(n_active(ic),ic)) - 1.
                end if
              endif
            end do
            where (bc% pred_use == -1) bc% pred_use = 1
          end do
        else

          cs = 'n'
          bc% pred_use = pred_use
          where(pred_use == 1) bc%pred_use = -1
          do
            n_active = count  (bc% pred_use == 1, dim=1)

            if (all(n_active >= count(pred_use /= 0, dim=1))) exit
            slot = slot + 1

            !------------------------
            ! switch on one predictor
            !------------------------
            write(6,*) 'switch on predictor ',n_active(1)+1
            p_switch    = 0
            do
              do ic = 1, n_chan
                if (p_switch(ic) > 0 .and. p_switch(ic) <= n_pred) bc% pred_use(p_switch(ic), ic) = -1
                if (.not. any(bc% pred_use(p_switch(ic)+1:n_pred, ic) == -1)) then
                  p_switch(ic) = n_pred + 1
                  cycle
                end if
                do j = p_switch(ic) + 1, n_pred
                  if (bc% pred_use(j, ic) == -1) then
                    p_switch(ic) = j
                    bc% pred_use(j, ic) = 1
                    exit
                  end if
                end do
              end do
              if (all(p_switch > n_pred)) exit

              !--------------------------------------------------------
              ! calculate bias correction for reduced set of predictors
              !--------------------------------------------------------
              call calc_bcor_coefs (bc% sta, bc)
              tmp = scbc_rmse ()
              do ic = 1, n_chan
                if (bc% sta% n(ic) < 2) cycle
                if (p_switch(ic) <= n_pred) rmse (slot, p_switch(ic), ic) = tmp (ic)
              end do
            end do

            !------------------------------------
            ! determine most important predictor
            !------------------------------------
            do ic = 1, n_chan
              if (n_active(ic) < count(pred_use(:,ic) /= 0)) then
                mask(:) = (bc% pred_use (:,ic) == -1)
                if (keep_used) then
                  where(pred_use(:,ic) /= 1) mask(:) = .false.
                  if (.not.any(mask)) mask(:) = (bc% pred_use (:,ic) == -1)
                end if
                i = sum (minloc (rmse (slot,:,ic), mask = mask(:)))
                bc% pred_use (i,ic) = 1
                p_ls   (n_active(ic)+1,ic) = i
                p_rmse (n_active(ic)+1,ic) = rmse (slot,i,ic)
                mask(i) = .false.
                if (bc% sta% n(ic) >= 2) then
                  if (any(mask(:))) i = sum(minloc (rmse (slot,:,ic), mask = mask(:)))
                  p_ls2  (n_active(ic)+1,ic) = i
                  p_rat2 (n_active(ic)+1,ic) = (p_rmse(n_active(ic),ic) - p_rmse(n_active(ic)+1,ic)) / &
                       (p_rmse(n_active(ic),ic) - rmse (slot,i,ic)) - 1.
                  if (p_rat2(n_active(ic)+1,ic)<0.) then
                    p_rat2(n_active(ic)+1,ic) = 0.
                  end if
                end if
              endif
              if (n_active(ic) == 0 .and. bc% sta% n(ic) > 2) then
                if (unit < 0) call get_free_unit (unit)
                write (cchan, '(i4.4)') ic
                n = count(bc% pred_use(:,ic) /= 0)
                allocate(ind(n))
                ind = pack((/ (j, j=1, n_pred) /), mask=bc% pred_use(:,ic) /= 0)
                do j = n-1, 1, -1
                  do k=1, j
                    if (rmse(slot,ind(k),ic) > rmse(slot,ind(k+1),ic)) then
                      i = ind(k)
                      ind(k) = ind(k+1)
                      ind(k+1) = i
                    end if
                  end do
                end do
                file = 'biasc_pred_rmse1_'//cs//'_'//csatid//'_'//cgrid//'_'//cchan//'.gp'
                open (unit, file=file, action='write')
                write(unit,'(a)')   '#'
                write(unit,'(a)')   '# rmse for predictor choices'
                write(unit,'(a)')   '#'
                write(unit,'(a,a)') '# file       = ',trim(file)
                write(unit,'(a,i5)')'# satid      =' ,bc% i% satid
                write(unit,'(a,i5)')'# grid       =' ,bc% i% grid
                write(unit,'(a,i5)')'# channel id :' ,ic
                write(unit,'(a)')   '#'
                write(unit,'(a)')   '# column  1  : number of predictor'
                write(unit,'(a)')   '# column  2  : rmse'
                write(unit,'(a)')   '# column  3  : rmse / rmse with no  predictors'
                write(unit,'(a)')   '# column  4  : predictor index'
                write(unit,'(a)')   '# column  5  : predictor name'
                write(unit,'(a)')   '#'
                do i = 1, n
                  j = ind(i)
                  write(unit,'(I3,2(1x,F13.8),1x,I3,1x,A50)') i, rmse(slot,j,ic), rmse(slot,j,ic)/rmse(1,j,ic),&
                       j, pred_name(j, bc)
                end do
                close(unit)
                deallocate(ind)
              end if
            end do
          end do
        end if

        !--------------------------------------
        ! restore old flags and bias correction
        !--------------------------------------
        write(6,*) 'restore old flags and bias correction'
        bc% pred_use = pred_use
        call calc_bcor_coefs (bc% sta, bc)

        !---------
        ! printout
        !---------
        write(6,*) 'printout'
        do ic = 1, n_chan
          if (bc% sta% n(ic) < 2) cycle
          if (n_used(ic) == 0) cycle
          write (cchan, '(i4.4)') ic

          if (unit < 0) call get_free_unit (unit)
          file = 'biasc_pred_rmse_'//cs//'_'//csatid//'_'//cgrid//'_'//cchan//'.gp'
          open (unit, file=file, action='write')
          write(unit,'(a)')   '#'
          write(unit,'(a)')   '# rmse for predictor choices'
          write(unit,'(a)')   '#'
          write(unit,'(a,a)') '# file       = ',trim(file)
          write(unit,'(a,i5)')'# satid      =' ,bc% i% satid
          write(unit,'(a,i5)')'# grid       =' ,bc% i% grid
          write(unit,'(a,i5)')'# channel id :' ,ic
          write(unit,'(a,i5)')'# predictors :' ,n_pred
          write(unit,'(a)')   '#'
          write(unit,'(a)')   '# column  1  : number of predictors used'
          write(unit,'(a)')   '# column  2  : rmse'
          write(unit,'(a)')   '# column  3  : rmse / rmse with all predictors'
          write(unit,'(a)')   '# column  4  : rmse / rmse with no  predictors'
          write(unit,'(a)')   '# column  5  : rmse / rmse without predictor'
          write(unit,'(a)')   '# column  6  : predictor index'
          write(unit,'(a)')   '# column  7  : predictor name'
          write(unit,'(a)')   '# column  8  : (rmse - rmse(1 pred. less))/(rmse2 - rmse(1 pred. less))'
          write(unit,'(a)')   '# column  9  : alternative predictor index (with rmse2)'
          write(unit,'(a)')   '# column 10  : alternative predictor name'
          write(unit,'(a)')   '#'
          write(unit,'(I3,4(1x,F13.8),1x,I3,1x,A50,1x,F13.8,1x,I3,1x,A50)') &
               0, p_rmse(0,ic),                      &
               p_rmse(0,ic) / p_rmse(n_used(ic),ic), &
               p_rmse(0,ic) / p_rmse(         0,ic), &
               p_rmse(0,ic) / p_rmse(         0,ic), &
               0, 'none', 1., 0, 'none'

          do i = 1, n_used(ic)
            j = p_ls(i,ic)
            if (p_ls2(i,ic) /= j) then
              pname = pred_name(p_ls2(i,ic), bc)
            else
              pname = 'none'
            end if
            if (ic==icpr) print*,i,n_used(ic),p_rmse(i,ic),p_rmse(n_used(ic),ic),p_rmse(i,ic) / p_rmse(n_used(ic),ic)
            write(unit,'(I3,4(1x,F13.8),1x,I3,1x,A50,1x,F13.8,1x,I3,1x,A50)') &
                 i, p_rmse(i,ic),                      &
                 p_rmse(i,ic) / p_rmse(n_used(ic),ic), &
                 p_rmse(i,ic) / p_rmse(         0,ic), &
                 p_rmse(i,ic) / p_rmse(       i-1,ic), &
                 j, trim(pred_name(j,bc)),             &
                 p_rat2(i,ic), p_ls2(i,ic),            &
                 trim(pname)
          end do
          close (unit)
        end do


      else
        !------------------------------
        ! brute force selection process
        !------------------------------
        bc% pred_use = pred_use
        where(pred_use == 1) bc%pred_use = -1

        n_active = count  (bc% pred_use /= 0, dim=1)
        m_active = maxval (n_active)
        if (m_active < 1) exit
        do icpr = 1, n_chan
          if (n_active(icpr) == m_active) exit
        end do

        allocate (pset   (m_active, n_chan))        ; pset    = .false.

        ! Count number of permutations
        do i = 1, n_sel
          pset(i,1) = .true.
        end do
        m = 0
        perm_loop0: do
          l_found = .false.
          do j = m_active, 2, -1
            if (.not.pset(j,1) .and. pset(j-1,1)) then
              pset(j,1) = .true.
              pset(j-1,1) = .false.
              pset(j+1:m_active,1) = .false.
              k = n_sel - count(pset(1:m_active,1))
              pset(j+1:j+k,1) = .true.
              l_found = .true.
              exit
            end if
          end do
          m = m + 1
          if (.not.l_found) exit perm_loop0
        end do perm_loop0

        ! Setup/find arrays for sums over different satellites
        s => null()
        do i = 1, size(sum_arr)
          if (sum_arr(i)% grid < 0) then
            sum_arr(i)% grid    = bc% i% grid
            sum_arr(i)% n_chan1 = bc% i% n_ch_i(1)
            s => sum_arr(i)
            allocate(s% chan(n_chan), s% satellites(size(biasc)),            &
                     s% pred_use(n_pred, n_chan), s% pred(m, 0:n_chan, n_sel), &
                     s% n(m, 0:n_chan), s% rmse(m, 0:n_chan), s% n_entry(n_chan))
            s% satellites = -1
            s% pred_use   = pred_use
            s% n_entry    = 0
            s% n          = 0
            s% rmse       = 0._wp
            instr = 1
            do j = 1, n_chan
              if (j > bc% i% o_ch_i(instr) + bc% i% n_ch_i(instr)) instr = instr + 1
              s% chan(j) = 100000 * bc% i% instr(instr) + bc% i% chan(j)
            end do
          end if
          if (bc% i% grid == sum_arr(i)% grid) then
            s => sum_arr(i)
            do j = 1, size(s% satellites)
              if (s% satellites(j) < 0) s%satellites(j) = bc% i% satid
              if (s% satellites(j) == bc% i% satid) exit
            end do
            exit
          end if
        end do


        write(6,*) 'calculate bias correction for all ',m,' permutations of ',n_sel,&
             ' out of ',m_active,' predictors'

        m_store = min(10000, m)

        deallocate (p_rmse  )
        allocate (pred   (m_store,  n_chan, n_sel)) ; pred    = 0
        allocate (p_rmse (m_store,  n_chan))        ; p_rmse  = huge(p_rmse(1,1))
        allocate (n_store(          n_chan))        ; n_store = 0
        allocate (pred_ind(n_sel))
        pset = .false.
        do i = 1, n_sel
          pset(i,:) = .true.
        end do

        n = m
        m = 0
        ipr=1
        perm_loop: do
          !------------------------------------------------
          ! calculate bias correction for set of predictors
          !------------------------------------------------
          do ic = 1, n_chan
            if (bc% sta% n(ic) < 2) cycle
            where(bc%pred_use(:,ic) == 1) bc%pred_use(:,ic) = -1
            j = 0
            do i = 1, n_pred
              if (bc%pred_use(i,ic) /= 0) then
                j = j + 1
                if (pset(j,ic)) bc%pred_use(i,ic) = 1
              end if
            end do
            if (count(bc%pred_use(:,ic) == 1) /= n_sel) stop 13
          end do

          call calc_bcor_coefs (bc% sta, bc)
          tmp = scbc_rmse ()

          !------------------------------------------------
           ! store permutation if rmse is small enough
          !------------------------------------------------
          do ic = 1, n_chan
            if (bc% sta% n(ic) < 2) cycle
            do i = 1, n_store(ic)
              if (tmp(ic) < p_rmse(i,ic)) exit
            end do
            if (i <= m_store) then
              if (i <= n_store(ic)) then
                p_rmse(i+1:m_store,ic) = p_rmse(i:m_store-1,ic)
                pred(i+1:m_store,ic,:) = pred(i:m_store-1,ic,:)
              end if
              p_rmse(i,ic) = tmp(ic)
              k = 0
              do j = 1, n_pred
                if (bc% pred_use(j,ic) == 1) then
                  k = k + 1
                  pred(i,ic,k) = j
                end if
              end do
              if (n_store(ic) < m_store) n_store(ic) = n_store(ic) + 1
            end if
          end do

          l_found = .false.
          do ic = 1, n_chan
            if (bc% sta% n(ic) < 2) cycle
            ! Find last true entry, that can be moved to the right
            do j = n_active(ic), 2, -1
              if (.not.pset(j,ic) .and. pset(j-1,ic)) then
                pset(j,ic) = .true.
                pset(j-1,ic) = .false.
                pset(j+1:n_active(ic),ic) = .false.
                k = n_sel - count(pset(1:n_active(ic),ic))
                pset(j+1:j+k,ic) = .true.
                l_found = .true.
                exit
              end if
            end do
          end do
          if (.not.l_found) exit perm_loop

          m = m + 1

          ! Store results for sum over different satellites
          if (associated(s)) then
            instr = 1
            do ic = 1, n_chan
              if (ic > bc% i% o_ch_i(instr) + bc% i% n_ch_i(instr)) instr = instr + 1
              k = 100000 * bc% i% instr(instr) + bc% i% chan(ic)
              do j = 1, size(s% chan)
                if (s% chan(j) == k) exit
              end do
              if (j <= size(s% chan)) then
                if (all(pred_use(:,ic) == s% pred_use(:,j))) then
                  if (bc% sta% n(ic) >= 2) then
                    i = 0
                    do k = 1, n_pred
                      if (bc% pred_use(k,ic) == 1) then
                        i = i + 1
                        pred_ind(i) = k
                      end if
                    end do
                    if (m <= s% n_entry(j)) then
                      if (.not.all(s% pred(m,j,:) == pred_ind(:))) call finish('diagn_pred','internal error')
                    else
                      s% pred(m,j,:) = pred_ind(:)
                    end if
                    if (tmp(ic) > 0._wp .and. rmse(1,1,ic) > 0._wp) then
                      tmp(ic) = min(tmp(ic), rmse(1,1,ic)) ! Only necessary for invalid values in tmp
                      s% rmse(m, j) = s% rmse(m, j) + tmp(ic)/rmse(1,1,ic) * bc% sta% n(ic)
                      s% n   (m, j) = s% n   (m, j) + bc% sta% n(ic)
                    else
                      s% rmse(m, j) = 9999._wp
                    end if
                    s% n_entry(j) = max(m, s% n_entry(j))
                  end if
                end if
              end if
            end do
          end if

          if (m*100/n >= ipr) then
            write(6,'(3x,I3,"%")') ipr
            ipr = m*100/n + 1
          end if
        end do perm_loop

        cs = 'p'
        write(cp,'(I1)') n_sel
        do ic = 1, n_chan
          if (bc% sta% n(ic) < 2) cycle
          if (n_used(ic) == 0) cycle
          write (cchan, '(i4.4)') ic

          if (unit < 0) call get_free_unit (unit)
          file = 'biasc_pred_rmse_'//cs//cp//'_'//csatid//'_'//cgrid//'_'//cchan//'.gp'
          open (unit, file=file, action='write')
          write(unit,'(a)')   '# rmse for best predictor selections'
          write(unit,'(a)')   '#'
          write(unit,'(a,a)') '# file       = ',trim(file)
          write(unit,'(a,i5)')'# satid      =' ,bc% i% satid
          write(unit,'(a,i5)')'# grid       =' ,bc% i% grid
          write(unit,'(a,i5)')'# channel id :' ,ic
          write(unit,'(a,i5)')'# predictors :' ,n_pred
          write(unit,'(a)')   '#'
          write(unit,'(a)')   '# column  1  : ranking'
          write(unit,'(a)')   '# column  2  : rmse'
          write(unit,'(a)')   '# column  3  : rmse / rmse with no  predictors'
          do i = 1, m_store
            write(unit,'(1x,I5,2(1x,F13.8),100(1x,I3,1x,A50))')  &
                 i,p_rmse(i,ic),p_rmse(i,ic)/rmse(1,1,ic),       &
                 (pred(i,ic,j), pred_name(pred(i,ic,j), bc), j=1,n_sel)
          end do
          close(unit)
        end do

        deallocate (pred, pset, n_store, pred_ind)
      end if

      !------------------------------
      ! deallocate arrays for this set
      !------------------------------
      deallocate (pred_use)
      deallocate (rmse    )
      deallocate (n_active)
      deallocate (n_used  )
      deallocate (p_switch)
      deallocate (tmp     )
      deallocate (p_rmse  )
      deallocate (p_ls    )
      deallocate (mask    )
      deallocate (p_rat2  )
      deallocate (p_ls2   )

    end do

    if (n_sel > 0) then
      allocate (pred_ind(n_sel))
      do i = 1, size(sum_arr)
        s => sum_arr(i)
        if (s% grid < 0) exit
        write(6,*) 'printout grid ',s% grid
        write (cgrid ,'(i2.2)') s% grid
        n_chan = size(s% chan)
        where(s% n > 0) s% rmse = s% rmse / s% n

        n_chan_aux = s% n_chan1
        l_found = .false.
        do j = 1, m
          do k = 1, n_sel
            if (.not.all(s% pred(j,1,k) == s% pred(j,2:n_chan_aux,k))) l_found = .true.
          end do
        end do
        if (.not.l_found) then
          s% rmse(:,0) = sum(s% rmse(:,1:n_chan_aux), dim=2)
          s% n   (:,0) = sum(s% n   (:,1:n_chan_aux), dim=2)
          s% pred(:,0,:) = s% pred(:,1,:)
          ic0 = 0
        else
          ic0 = 1
        end if

        do ic = ic0, n_chan
          if (all(s%n(:,ic) <= 0)) cycle
          if (ic == 0) then
            cchan='summ'
          else
            write (cchan, '(i4.4)') ic
          end if
          ! Sort
          write(6,*) '   sort channel ',ic
          do k = m-1, 1, -1
            do j = 1, k
              if (s% rmse(j, ic) > s% rmse(j+1,ic)) then
                rmse_aux          = s% rmse(j,ic)
                pred_ind          = s% pred(j,ic,:)
                s% rmse(j,ic)     = s% rmse(j+1,ic)
                s% pred(j,ic,:)   = s% pred(j+1,ic,:)
                s% rmse(j+1,ic)   = rmse_aux
                s% pred(j+1,ic,:) = pred_ind
              end if
            end do
          end do

          if (unit < 0) call get_free_unit (unit)
          file = 'biasc_pred_rmse_'//cs//cp//'_all_'//cgrid//'_'//cchan//'.gp'
          open (unit, file=file, action='write')
          write(unit,'(a)')   '# rmse for best predictor selections'
          write(unit,'(a)')   '#'
          write(unit,'(a,a)') '# file       = ',trim(file)
          write(unit,'(a,20(1x,I3.3))')'# satids     =' ,pack(s% satellites(:), mask=s% satellites(:)>0)
          write(unit,'(a,i5)')'# grid       =' ,s% grid
          if (ic >= 1) then
            write(unit,'(a,i5)')'# instrument :' ,s% chan(ic)/100000
            write(unit,'(a,i5)')'# channel id :' ,mod(s% chan(ic),100000)
          end if
          write(unit,'(a,i5)')'# predictors :' ,n_pred
          write(unit,'(a)')   '#'
          write(unit,'(a)')   '# column  1  : ranking'
          write(unit,'(a)')   '# column  2  : rmse / rmse with no  predictors'
          do j = 1, m
            write(unit,'(1x,I9,1x,F13.8,100(1x,I3,1x,A50))') &
                 j,s% rmse(j,ic),                            &
                 (s% pred(j,ic,k), pred_name(s% pred(j,ic,k), bc), k=1,n_sel)
          end do
          close(unit)

        end do
      end do
      deallocate (pred_ind)
    end if

  contains
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function scbc_rmse ()
    real(wp) :: scbc_rmse (n_chan)
    !-------------------------------------
    ! calculate rmse after bias correction
    !-------------------------------------

      real(wp) :: n     ! number of entries per Field of view
      real(wp) :: nn    ! total number of entries
      real(wp) :: s     ! rmse after bias correction
      real(wp) :: b     ! constant bias correction coefficient
      integer  :: ic    ! channel index
      integer  :: if    ! field of view index
      integer  :: i,j,k ! predictor index

      do ic = 1, n_chan
        s  = 0._wp
        nn = 0._wp
        if (associated (bc% sta% btdn)) then
          do if = 1, n_fov
            n = bc% sta% btdn (if,ic)
            b = bc% coef_set1 (n_pred+1, ic) + bc% coef_set2 (ic,if)
            if (n > 0._wp) then
              nn = nn + n
              s  = s + bc% sta% d2 (if,ic)                            &!
                 - 2 * bc% sta% btdscan (if,ic) * b                   &! d * b
                 + n * b * b                                           ! b * b
              k = 1
              do i = 1, n_pred
                if (bc% pred_use(i,ic) == 1) then
                  s = s - 2 * bc% sta% btd_pred_i (i,if,ic)           &! d * i
                            * bc% coef_set1       (i,   ic)           &!
                        + 2 * bc% sta% pred_i     (i,if,ic)           &! i * b
                            * bc% coef_set1       (i,   ic) * b       &!
                        +     bc% coef_set1       (i,   ic) ** 2      &!
                            * bc% sta% pred_ij    (k,if,ic)            ! i * i
                endif
                do j = i, n_pred
                  if (bc% pred_use(j,ic) == 1 .and. &
                      bc% pred_use(i,ic) == 1 .and. &
                                       i /= j       ) then
                    s = s + 2 * bc% coef_set1       (i,   ic)         &! i * j
                              * bc% coef_set1       (j,   ic)         &!
                              * bc% sta% pred_ij    (k,if,ic)
                  endif
                  k = k + 1
                end do
              end do
            endif
          end do
        endif
        if (s >= 0._wp) then
          if (nn > 0._wp) s = sqrt (s / nn)
          scbc_rmse (ic) = s
        else
          scbc_rmse (ic) = huge(0._wp)
        end if

      end do

    end function scbc_rmse
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function scbc_stdev ()
    real(wp) :: scbc_stdev (n_chan)
    !------------------------------------------
    ! calculate rmse with constant bias removed
    !------------------------------------------

      real(wp) :: n   ! number of entries per Field of view
      real(wp) :: nn  ! total number of entries
      real(wp) :: s   ! rmse with constant bias removed
      integer  :: ic  ! channel index
      integer  :: if  ! field of view index

      do ic = 1, n_chan
        s  = 0._wp
        nn = 0._wp
        if (associated (bc% sta% btdn)) then
          do if = 1, n_fov
            n = bc% sta% btdn (if,ic)
            if (n > 0._wp) then
              s  = s + bc% sta% d2 (if,ic) - bc% sta% btdscan (if,ic) ** 2 / n
              nn = nn + n
            endif
          end do
        endif
        if (nn > 0._wp) s = sqrt (s / nn)
        scbc_stdev (ic) = s
      end do

    end function scbc_stdev
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine diagn_pred

!==============================================================================

  function Y_lm (l, m, dlat, dlon)
  !-------------------------------------------------------
  ! calculation of spherical harmonics
  ! currently M=0 only
  ! no normalisation (not required if used as a predictor)
  !-------------------------------------------------------
  integer  ,intent(in) :: l
  integer  ,intent(in) :: m
  real(wp) ,intent(in) :: dlat (:)
  real(wp) ,intent(in) :: dlon (:)
  real(wp)             :: Y_lm (size(dlat))

   real(wp) ,parameter :: pi  = 3.141592653589793238_wp
   real(wp) ,parameter :: d2r = pi / 180._wp
   real(wp)            :: x (size(dlat))
   real(wp)            :: p (size(dlat), 0:l)
   integer             :: n

   if (m /= 0) call finish('Y_lm','m /= 0 currently not implemented')

   x = sin ( d2r * dlat )

   select case (l)
   case (0)
     Y_lm = 1._wp
   case (1)
     Y_lm = x
   case (2:)
     p (:,0) = 1
     p (:,1) = x
     do n = 1, l-1
       p (:,n+1) = ((2*n+1) * x * p (:,n) - n * p (:,n-1)) / (n+1)
     end do
     Y_lm = p (:,l)
   case default
     call finish ('Y_lm','l<1 not implemented')
   end select

  end function Y_lm
!==============================================================================
  logical function equal_fov_bc(a,b)
    type(t_fov_bc), intent(in) :: a, b
    equal_fov_bc = .true.
    if (a% type       /=     b% type      ) equal_fov_bc = .false.
    if (a% n          /=     b% n         ) equal_fov_bc = .false.
    if (a% const_fov  /=     b% const_fov ) equal_fov_bc = .false.
    if (a% const_pred .neqv. b% const_pred) equal_fov_bc = .false.
  end function equal_fov_bc
!==============================================================================
  subroutine finish (name, text)
  character(len=*) :: name
  character(len=*) :: text
    call model_abort (-1, 9999, text, name, 0)
  end subroutine finish
!==============================================================================

! specific MPI-bcast for derived type t_fov_bc
#undef  DERIVED
#define VECTOR
#define DERIVED type(t_fov_bc), dimension(:)
#define p_bcast_DERIVED p_bcast_fov_bc
#undef  MPI_TYPE
#include "p_bcast.incf"
!==============================================================================
end module mo_radbiascor
