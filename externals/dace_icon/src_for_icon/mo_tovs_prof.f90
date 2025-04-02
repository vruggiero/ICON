!
!+ provides TOVS routines with atmospheric profile etc.
!
MODULE mo_tovs_prof
!
! Description:
!   Provides TOVS routines with atmospheric profile, surface properties and other
!   (derived) quantities.
!
! Current Code Owner: DWD, Robin Faulwetter
!    phone: +49 69 8062 2746
!    fax:   +49 69 8062 3721
!    email: robin.faulwetter@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Robin Faulwetter  DWD  2019       original source
!==============================================================================

  !=============
  ! modules used
  !=============
  use iso_fortran_env,  only: stderr => error_unit               ! Standard error
  use mo_exception,     only: finish,                           &! abort routine
                              message                            ! write message to stderr
  use mo_kind,          only: wp, i1, sp                         ! kind parameters
  use mo_dace_string,   only: tolower
  use mo_run_params,    only: ana_time,                         &! analysis time
                              interp_strato
  use mo_cntrlvar,      only: tq_tvgh,                          &! t,spec.hum. from virt.t,gen.hum.
                              tq_tvgh_vec,                      &!
                              tq_tvgh_adj_vec,                  &!
                              rh_gh,                            &!
                              gh_min,                           &!
                              gh_max
  use mo_mpi_dace,      only: dace,                             &! MPI group info
                              p_gatherv,                        &! MPI gather routine
                              p_sum,                            &
                              p_max
  use mo_obs_set,       only: t_obs_set,                        &!
                              t_obs_block,                      &!
                              obs_block                          !
  use mo_instrid,       only: rttov_instr,                      &! RTTOV number from WMO instrument id
                              mw_instr                           ! MW instrument or not
  use mo_physics,       only: gacc,                             &! g
                              rh_q,                             &! rh <- t, q
                              Rdry => R,                        &! gas constant of dry air [J/(kg*K)]
                              gas_id_o3, gas_id_co2,            &!
                              tv_t_q,                           &! virtual temperat. <- temp., spec.hum., liq/ice
                              t_tv_q_adj,                       &! temp<-vir.temp.,spec.hum.;adjoint
                              t_tv_q,                           &! temp. <- virt.temp., spec.hum.
                              q_rh,                             &!
                              q_rh_adj
  use mo_time,          only: operator(-),                      &!
                              t_time,                           &!
                              localdate                          ! time as minutes
  use mo_t_obs,         only: t_obs,                            &!
                              t_spot,                           &!
                              po_context,                       &! process_obs context
                              po_ilns,                          &! process_obs context: number of call within linesearch
                              po_lmon,                          &! process_obs context: monitoring call?
                              ldeb,usd,dpref                     ! debug selected spot(s)
  use mo_sink,          only: bc_xd                              ! x from control variable
  use mo_t_col,         only: COL_T,                            &! T
                              COL_Q,                            &! q
                              COL_CLC,                          &! cloud cover
                              COL_QCDIA,                        &! spec. cloud water, diag.
                              COL_QIDIA,                        &! spec. cloud ice, diag.
                              COL_REFF_QC,                      &! eff. radius cloud water
                              COL_REFF_QI                        ! eff. radius cloud ice
  use mo_fdbk_tables,   only: OT_RAD,                           &! Radiances report type ID
                              TF_EMIS_SEA_MOD,                  &! Emiss. sea model
                              TF_REFL_BRDF,                     &! BRDF
                              TF_EMIS_FAILED,                   &! emiss. retrieval failed
                              TF_EMIS,                          &! all emiss. flags
                              SUR_LAND,                         &! surftype: land
                              SUR_SEA,                          &! surftype: sea
                              SUR_ICE,                          &! surftype: ice
                              SUR_NOICE,                        &! surftype: noice
                              SUR_SNOW,                         &! surftype: snow
                              SUR_NOSNOW                         ! surftype: nosnow
  use mo_t_use,         only: STAT_DEFAULT,                     &! invalid value flag
                              STAT_DISMISS,                     &! status flag: dismiss report
                              STAT_PASSIVE,                     &! status flag: passive report
                              STAT_REJECTED,                    &! status flag: rejected report
                              CHK_CLOUD                          ! cloud data flag
  use mo_hum_ana,       only: hum_ana_top                        ! log(hum_ana_top)
  use mo_rttov,         only: destruct,                         &! destruct radiances data type
                              construct,                        &! construct t_rttov_prof
                              call_rttov,                       &!
                              t_rttov_prof,                     &! rttov variable input data type
                              rttvi_called,                     &! whether call_rttvi was called
                              preslev,                          &! reference pressure levels (Pa)
                              qmin, qmax,                       &! min,max. q for valid RTTOVS calcs.
                              tmin, tmax,                       &! min,max. T for valid RTTOVS calcs.
                              rttov_bounds,                     &!
                              nsv,                              &! # of single level input variables
                              rt_hard_limits,                   &! Check RTTOV hard limits
                              q2rt_humi,                        &! convert humidity for RTTOV
                              rt_trg                             ! convert tracegases for RTTOV
  use mo_rad,           only: t_rad_set,                        &! specific. of instr.
                              t_radv,                           &! Complete satellite dataset
                              rad_set,                          &!
                              n_set,                            &!
                              print_rad_set,                    &!
                              rinvalid,                         &! invalid value
                              set_indx,                         &!
                              chan_indx,                        &!
                              m_chan,                           &! maximum number of channels
                              m_bd,                             &! maximum number of bands
                              n_styp,                           &! number of RTTOV surface types
                              destruct,                         &! destruct t_radv and t_rad_set
                              assignment(=),                    &! assignment of t_rad_set
                              amsua_chan_id,                    &!
                              amsub_chan_id,                    &!
                              lev2p,                            &! level to pressure
                              OPTV_CLD_FRC                       ! imager cloud fraction availability
  use mo_t_tovs,        only: t_tovs,                           &! observation operator specific type
                              t_tovs_instr,                     &!
                              tpp,                              &!
                              destruct,                         &! t_tovs  destructor routine
                              get_tovs_rs,                      &!
                              load,                             &! load  t_tovs
                              store,                            &! store t_tovs
                              mx_nav,                           &! max. size of av (2nd dim)
                              mx_nlev,                          &!
                              TTOVS_BASE_BIT,                   &!
                              TTOVS_CI_BIT,                     &!
                              TTOVS_AV_BIT,                     &!
                              TTOVS_L2C_BIT,                    &!
                              TTOVS_EMIS_BIT,                   &!
                              TTOVS_FLAG_BIT,                   &!
                              TTOVS_CLDLEV_BIT,                 &! 
                              TTOVS_BASE,                       &!
                              TTOVS_CI, TTOVS_AV,               &! Constants that determine, which parts of
                              TTOVS_CEMI,                       &! t_tovs are to be stored/loaded
                              TTOVS_EMIS, TTOVS_FLAG, TTOVS_SPEC
  use mo_ir_emis,       only: emis_pc,                          &! calculate emissivities
                              n_pc                               ! # of PC to use; -1: off
  use mo_emis,          only: set_emis,                         &!
                              print_emis_stat,                  &!
                              EMIS_CALC, EMIS_SUCCESS
  use mo_tskin,         only: set_tskin,                        &
                              print_tskin_stat,                 &
                              TSKIN_CALC
  use mo_rtifc,         only: rts_land, rts_sea, rts_ice         ! rttov surface types
  use mo_dec_matrix,    only: t_vector,                         &! vector data type
                              t_vector_segm,                    &! vector segment data type
                              t_matrix,                         &! matrix data type
                              t_matrix_block,                   &! matrix blocks
                              mp,                               &! kind parameter for matrix elements
                              construct,                        &! allocate data type components
                              destruct,                         &! deallocate data type components
                              deallocate,                       &! deallocate matrix block components
!                             gather,                           &! gather vector content on PE
!                               operator(*),                      &! matrix-vector multiply
!                               operator(+),                      &! add vectors
                              assignment(=),                    &! assign matrices or vectors
                              sub_block,                        &! extract from matrix block
                              FULL,                             &!
                              full_matrix                        ! store matrix in full   format
  use mo_matrix,        only: inverse_rs                         ! inverse of a real symmetric matrix
  use mo_fg_cov,        only: t_rowcol,                         &! data type to hold information
                              construct,                        &! construct type t_rowcol
                              destruct,                         &! destruct  type t_rowcol
                              get_corr,                         &! get correlations/covariances
                              init_spot_obs                      !
  use mo_fparse,        only: lv

  use mo_rtifc,         only: OUT_CSB,                          &! clear-sky BT
                              rts_land

  use mo_cloud_indices, only: mw_emiss, lwp_rain_thresh,        &! AMSUA-based lwp
                              lwp_QinZou2016                     ! MHS-based lwp

  implicit none

  private

  public :: fill_rad          ! fill radiances from t_obs_set into t_radv
  public :: prep_rttov_prof   ! prepare t_rttov_prof type swith profile info
  public :: destruct          ! t_jac, t_jac_arr
  public :: construct         ! t_jac
  public :: valid
  public :: t_jac
  public :: t_jac_arr
  public :: set_surface_type
  public :: surf_class_vers
  public :: check_recalc_stype
  public :: e_bg_ts_sea
  public :: e_bg_ts_land
  public :: e_bg_drts_land
  public :: e_bg_ts_ice
  public :: e_bg_t_top
  public :: e_bg_lnq_top
  public :: no_ps_dep
!  public :: no_t_dep
  public :: force_fastem
  public :: use_hum_top      ! deprecated
  public :: fg_prof_top
  public :: warn_gen_hum
  public :: p_top
  public :: min_dia_qc
  public :: min_dia_qi
  public :: scale_qi
  public :: snw_frc_mode
  public :: snw_frc_par
  ! Stuff for get_tovs_var
  public :: get_tovs_var
  public :: n_var
  public :: c_var
  public :: errmsg

  !----------------------------------------------------
  ! coefficients  for calculation of the Jacobian:
  ! d RTTOV_input_variables / d 3DVAR_control_variables
  !----------------------------------------------------
  type t_jac
    logical           :: ps_llev
    integer           :: il
    integer           :: iu
    real(wp)          :: wl
    real(wp)          :: wu
    real(wp)          :: d_lev
    real(wp)          :: d_tv
    real(wp)          :: d_rh
    real(wp), pointer :: dt_tv (:)  =>NULL() ! d T    / d Tv
    real(wp), pointer :: dt_rh (:)  =>NULL() ! d T    / d rh
    real(wp), pointer :: dq_tv (:)  =>NULL() ! d q    / d Tv
    real(wp), pointer :: dq_rh (:)  =>NULL() ! d q    / d rh
    real(wp), pointer :: de_dpc(:,:)=>NULL() ! d emis / d PC
    real(wp)          :: dt_tv_2m            ! d T_2m / d Tv_2m
    real(wp)          :: dt_rh_2m            ! d T_2m / d rh_2m
    real(wp)          :: dq_tv_2m            ! d q_2m / d Tv_2m
    real(wp)          :: dq_rh_2m            ! d q_2m / d rh_2m
    real(wp)          :: sigma_var_tskin     ! surface temperature error
    real(wp)          :: dclt_dum            ! d clt  / d dummy
    real(wp)          :: dclf_dum            ! d clf  / d dummy
    real(wp)          :: dsnf_dum            ! d snf  / d dummy
  end type t_jac

  type t_jac_arr
    type(t_jac), pointer  :: a(:) => NULL()
  end type t_jac_arr

  ! Options relevant in prep_rttov_prof
  logical, save :: force_fastem    = .true.   ! force calls to FASTEM. See note in prep_rttov_prof
  real(wp),save :: p_top           = -1._wp   ! model top level pressure
  integer, save :: surf_class_vers = 0
  logical, save :: check_recalc_stype = .false.
  real(wp),save :: e_bg_ts_sea     = 1._wp    ! surf.temp. bg.error sea
  real(wp),save :: e_bg_ts_land    = 3._wp    ! surf.temp. bg.error land
  real(wp),save :: e_bg_drts_land  = 1.5_wp    ! retrieved surf. temp. bg. error land
  real(wp),save :: e_bg_ts_ice     = 3._wp    ! surf.temp. bg.error ice
  logical, save :: warn_gen_hum    = .false.  ! Warning if gen. hum. bounds are exceeded
  logical       :: scale_qi        = .false.  ! scale cloud mass according to ratio of
                                              ! original and clipped diameters
  real(wp)      :: min_dia_qc      = 5.0_wp   ! minimum effective diameter of cloud droplets
  real(wp)      :: min_dia_qi      = 10.0_wp  ! minimum effective diameter of cloud ice
  ! stratospheric profiles
  real(wp),save :: e_bg_t_top      = 3.0_wp   ! bg. error top levels t (above level t_tovs%nld_t)
  real(wp),save :: e_bg_lnq_top    = 1.0_wp   ! bg. error top levels q (above level t_tovs%nld_h)
  integer, save :: use_hum_top     = 0        ! DEPRECATED, controls humidity profile above
                                              ! level_hum_dum/t_tovs%ld_h:
                                              !   1:fg-profile+e_bg_lnq_top=0., i.e. q is fixed
                                              !   2:fg-profile and q is a dummy variable
                                              ! We should better use fg_prof_top and e_bg_*_top
  integer, save :: fg_prof_top     = 0        ! 1:fg-profile
  !    temporary work arounds for discontinuities in RTTOV10
  logical, save :: no_ps_dep       = .false.  ! no dependence on surf.pres.
!  integer, save :: no_t_dep        = 0        ! no dependence on toplevel t
  integer       :: snw_frc_mode    = -1       ! Snow fraction handling
  real(wp)      :: snw_frc_par(5)  = 0._wp    ! Snow fraction parameterization parameters

  integer       :: npv                        ! no of profile inputs


  ! Stuff required for get_tovs_var
  integer, parameter :: n_var = 37
  character(len=lv)  :: c_var(n_var)
                        data c_var( 1) /'iwv'              /
                        data c_var( 2) /'rh_*'             /  ! _hgt [pressure level in hPa] must be added
                        data c_var( 3) /'q_*'              /  ! _hgt [pressure level in hPa] must be added
                        data c_var( 4) /'t_*'              /  ! _hgt [pressure level in hPa] must be added
                        data c_var( 5) /'sat_zen'          /
                        data c_var( 6) /'sat_azi'          /
                        data c_var( 7) /'sun_zen'          /
                        data c_var( 8) /'sun_azi'          /
                        data c_var( 9) /'rt_stype'         /
                        data c_var(10) /'mw_stype'         /
                        data c_var(11) /'ps'               /
                        data c_var(12) /'lat'              /
                        data c_var(13) /'lon'              /
                        data c_var(14) /'phase'            /
                        data c_var(15) /'fov'              /
                        data c_var(16) /'emis*'            /  ! _chan might be added
                        data c_var(17) /'l2c*'             /  ! _chan might be added
                        data c_var(18) /'p_l2c*'           /  ! _chan might be added
                        data c_var(19) /'plevel*'          /  ! _chan might be added
                        data c_var(20) /'z'                /
                        data c_var(21) /'mdlsfc'           /
                        data c_var(22) /'loctime'          /
                        data c_var(23) /'lwp'              /
                        data c_var(24) /'clw'              /
                        data c_var(25) /'lwp_qinzou2016'   /
                        data c_var(26) /'lwp_qinzou2016_s' /
                        data c_var(27) /'lwp_qin'          /
                        data c_var(28) /'lwp_qin_s'        /
                        data c_var(29) /'land_frac'        /
                        data c_var(30) /'lfr'              /
                        data c_var(31) /'tskin'            /
                        data c_var(32) /'ts'               / 
                        data c_var(33) /'abs_tobs_tana'    / 
                        data c_var(34) /'tobs_tana'        / 
                        data c_var(35) /'cldlev'           / 
                        data c_var(36) /'p_cldlev'         / 
                        

  character(len=300)    :: errmsg

  integer, parameter    :: mlev = 120

  type(t_tovs), target  :: tt_dum
  save                  :: tt_dum
  integer,      target  :: flag  (m_chan)
  real(tpp),    target  :: av    (mlev,10)
  real(sp),     target  :: l2c   (m_chan)
  real(wp),     target  :: emis  (m_chan)
  integer,      target  :: ci    (m_chan)
  real(sp),     target  :: cldlev(m_bd)

  interface destruct
    module procedure destruct_t_jac
    module procedure destruct_t_jac_arr
  end interface destruct

  interface construct
    module procedure construct_t_jac
  end interface construct

!==============================================================================
contains

#if defined (_FTRACE)
#define FTRACE_BEGIN(text) CALL FTRACE_REGION_BEGIN (text)
#define FTRACE_END(text)   CALL FTRACE_REGION_END   (text)
#else
#define FTRACE_BEGIN(text)
#define FTRACE_END(text)
#endif

  subroutine fill_rad(rad, set, obs, x, y, e_bg, lrttov, lprint, tj_set, min_stat, &
       ibox_s, ibox_e, pe_gather, lmon, lH, lB, lHBH, lmeta, HBHR, Roo,  zi, rhs, lpartial, lerr,&
       nadd, cp_opts, lbc, mask_rs)
  !-------------------------------------------
  ! fill type t_radv with data from type t_obs
  !-------------------------------------------
  ! TODO: more intuitive mechanism to select arrays that are to be filled.
  type (t_radv)    ,intent(out) ,target   :: rad(:)     ! radiance data
  type (t_rad_set) ,intent(in)  ,target   :: set(:)     ! meta data for 'rad'
  type (t_obs_set) ,intent(in)            :: obs        ! observation data (3dvar type)
  type (t_vector)  ,intent(in)  ,optional :: x          ! background, interpolation space
  type (t_vector)  ,intent(in)  ,optional :: y          ! background, observation   space
  type (t_vector)  ,intent(in)  ,optional :: e_bg       ! background error
  logical          ,intent(in)  ,optional :: lrttov
  logical          ,intent(in)  ,optional :: lprint
  type(t_jac_arr)  ,intent(out) ,optional :: tj_set(:)
  integer          ,intent(in)  ,optional :: min_stat   ! minimum status
  integer          ,intent(in)  ,optional :: ibox_s     ! start box (default: 1)
  integer          ,intent(in)  ,optional :: ibox_e     ! end   box (default: number of boxes)
  integer          ,intent(in)  ,optional :: pe_gather  ! pe on which all datasets shall be gathered
  logical          ,intent(in)  ,optional :: lmon
  logical          ,intent(in)  ,optional :: lH
  logical          ,intent(in)  ,optional :: lB
  logical          ,intent(in)  ,optional :: lHBH
  logical          ,intent(in)  ,optional :: lmeta
  type(t_matrix)   ,intent(in)  ,optional :: HBHR
  type(t_matrix)   ,intent(in)  ,optional :: Roo
  type(t_vector)   ,intent(in)  ,optional :: zi
  type(t_vector)   ,intent(in)  ,optional :: rhs
  logical          ,intent(in)  ,optional :: lpartial
  logical          ,intent(in)  ,optional :: lerr
  integer          ,intent(in)  ,optional :: nadd(:)    ! space for additional spots
  logical          ,intent(in)  ,optional :: cp_opts    ! copy i%Xopts from global rad_set to rad%i
  logical          ,intent(in)  ,optional :: lbc
  logical          ,intent(in)  ,optional :: mask_rs(:) ! mask rad_set, i.e. exclude some datasets

    type(t_spot),        pointer :: s
    type(t_radv),        pointer :: r => NULL ()
    type(t_rad_set),     pointer :: rs
    type(t_obs_block)            :: ob
    type(t_rttov_prof)           :: rtp
    type(t_tovs)                 :: tovs
    type(t_jac),         target  :: tj
    type(t_jac),         pointer :: tjp
    type(t_time)                 :: time
    integer                      :: ib
    integer                      :: is
    integer                      :: i
    integer                      :: j
    integer                      :: l
    integer                      :: m
    integer                      :: ii
    integer                      :: l2, i2, m2
    integer                      :: iset
    integer                      :: satid
    integer                      :: grid
    integer                      :: instr
    integer                      :: n_instr
    integer                      :: ichan
    integer                      :: pe_gather_loc
    integer                      :: i_rs
    integer                      :: nset (size(set))       ! # profiles in set
    integer                      :: nlev (size(set))       ! # levels in set
    integer                      :: nalloc                 ! # profiles to be allocated in set
    integer                      :: ib_s, ib_e             ! start/end box
    integer                      :: i_cloud(size(set))
    logical                      :: l_H, l_K, l_B, l_HBH
    logical                      :: l_rttov, l_print, l_mon, l_eb, l_meta, l_pc, l_pc_aux, &
                                    l_partial, l_err
    logical                      :: l_rt_lev(size(set))
    logical                      :: l_lamb(size(set))
    logical                      :: l_btcs(size(set))
    logical                      :: l_im_cfrac(size(set))
    integer                      :: n_pr_tr(size(set))
    logical                      :: l_x                    ! input profiles are present
    logical                      :: l_cp_opts
    logical                      :: l_bc
    integer                      :: ndeb
    integer                      :: max_nlev
    real(wp)                     :: fr, snf
    ! --- H ---
    real(wp),        allocatable :: H(:,:)                 ! temporary for H
    type(t_matrix_block)         :: Hb
    real(wp),        allocatable :: tmp(:)                 ! temporary
    real(wp),        allocatable :: jk(:,:,:)              ! jacobian cntrl->rttov
    real(wp),        allocatable :: kj(:,:,:)              ! inverse jk
    real(wp)                     :: edet                   ! 1 / determinant
    integer                      :: o0, o1, ox, on, i0, i1, in, ip
    integer                      :: ncv                    ! total no of inputs
    integer                      :: npc                    ! no of PC inputs (this spot)
    ! --- K ---
    type(t_matrix_block)         :: HBHRb                  !
    type(t_matrix_block)         :: Rb                     !
    type(t_rowcol)               :: rc                     ! vert. fg.correlation info
    type(t_vector)               :: v                      ! working variable
    type(t_rttov_prof)           :: rtp_k                  ! argument to rttov
    real(wp),        allocatable :: K (:,:)                ! (m_chan,2*nlev+2)
    real(wp),        allocatable :: HBHRi (:,:)            !
    real(mp),        allocatable :: Pb (:,:)               ! Pb in interpolation space
    ! --- B ---
    real(wp),        allocatable :: B(:,:)                 ! temporary for B
    type(t_matrix_block)         :: Bb

    FTRACE_BEGIN('fill_rad')

    FTRACE_BEGIN('fill_rad:setup')
    if (present(lrttov  )) then ; l_rttov   = lrttov     ; else ; l_rttov   = .false. ; endif
    if (present(lprint  )) then ; l_print   = lprint     ; else ; l_print   = .false. ; endif
    if (present(lmon    )) then ; l_mon     = lmon       ; else ; l_mon     = .false. ; endif
    if (present(lH      )) then ; l_H       = lH         ; else ; l_H       = .false. ; endif
    if (present(lB      )) then ; l_B       = lB         ; else ; l_B       = .false. ; endif
    if (present(lHBH    )) then ; l_HBH     = lHBH       ; else ; l_HBH     = .false. ; endif
    if (present(lmeta   )) then ; l_meta    = lmeta      ; else ; l_meta    = .true.  ; endif
    if (present(lpartial)) then ; l_partial = lpartial   ; else ; l_partial = .false. ; endif
    if (present(lerr    )) then ; l_err     = lerr       ; else ; l_err     = .false. ; endif
    if (present(lbc     )) then ; l_bc      = lbc        ; else ; l_bc      = .false. ; endif

    if (l_HBH .and. .not.present(Roo)) call finish('fill_rad','lHBH=T requires Roo argument')

    l_x  = present(x)
    l_eb = l_mon .and. present(e_bg)
    l_K  = l_H .and. present(HBHR) .and. present(zi) .and. present(rhs)
    if (present(ibox_s) .and. present(ibox_e)) then
      ib_s = ibox_s
      ib_e = ibox_e
    else
      ib_s = 1
      ib_e = size(obs% o)
    end if

    if (present(pe_gather)) then
      pe_gather_loc = pe_gather
    else
      pe_gather_loc = -99
    end if

    if (present(cp_opts)) then
      l_cp_opts = cp_opts
    else
      l_cp_opts = .false.
    end if


    npc  = max (n_pc, 0)     ! number of principal components
    l_pc = (npc > 0) .and. l_mon

    !-----------
    ! count data
    !-----------
    nset = 0
    nlev = 0
    do ib = ib_s, ib_e                        ! loop over observation 'boxes'
      if (obs% o(ib)% pe /= dace% pe) cycle   ! only on this PE
      do is=1,obs% o(ib)% n_spot              ! loop over reports/fovs
        s => obs% o(ib)%spot(is)              !
        if (s% hd% obstype /= OT_RAD) cycle   ! only treat radiances
        if (s% use% state  <= STAT_DISMISS) cycle
        if (present(min_stat)) then
          if (s% use% state  <= min_stat) cycle
        end if
        !------------------------------------------
        ! chose combination of satellite/instrument
        !------------------------------------------
        satid = s% hd% satid                  ! satellite id
        grid  = s% hd% grid_id                ! grid id
        iset = set_indx (set, satid=satid, grid=grid)
        if (iset < 1) then
          if (l_partial) then
            cycle
          else
            do j = 1, size(set)
              call print_rad_set(set(j), hint=dace% pe, unit=stderr)
            end do
            write(stderr,*) 'satid,grid',satid,grid
            call finish('fill_rad','no suitable dataset description (t_rad_set)')
          end if
        end if
        if (nset(iset) == 0) then
          call load(obs% o(ib), s, tovs, tovs_io=TTOVS_BASE)
          nlev(iset) = tovs%nlev
          call destruct(tovs)
        end if
        nset(iset) = nset(iset) + 1
      end do
    end do
    FTRACE_END('fill_rad:setup')

    ! Determine global rad_set and evaluate options in rad_set(:)%gopts and iopts.
    ! This is necessary, since the "set" supplied to fill_rad
    ! might be originating from e.g. the bias correction and might not contain the
    ! correct options (gopts, opts)
    do iset = 1, size(set)
      if (present(mask_rs)) then
        if (.not.mask_rs(iset)) CYCLE
      end if
      i_rs = set_indx(rad_set(1:n_set), satid=set(iset)%satid, grid=set(iset)%grid)
      if (i_rs < 1) then
        i = p_sum(nset(iset))
        if (.not.l_partial .and. i > 0) then
          do j = 1, n_set
            call print_rad_set(rad_set(j), hint=dace% pe, unit=stderr)
          end do
          write(stderr,*) 'satid,grid',set(iset)%satid,set(iset)%grid
          call finish('fill_rad','no suitable dataset description (t_rad_set)')
        end if
      end if
      r    => rad(iset)
      r% i =  set (iset)
      !-----------------
      ! evaluate options
      !-----------------
      i_cloud (iset)   = 0
      l_rt_lev(iset)   = .true.
      l_lamb  (iset)   = .false.
      l_im_cfrac(iset) = .false.
      n_pr_tr (iset) = 0
      r%n_trg = 0
      if (i_rs > 0) then
        rs => rad_set(i_rs)
        if (l_cp_opts) then
          r%i%gopts = rs%gopts
          r%i%iopts = rs%iopts
        end if
        if (any(rs%iopts(1:rs%n_instr)%cloud_mode > 0)) i_cloud(iset) = 1
        if (any(rs%iopts(1:rs%n_instr)%cloud_mode < 0)) i_cloud(iset) = -1
        if (i_cloud(iset) /= -1 .and. btest(rs%gopts%opt_vars, OPTV_CLD_FRC)) l_im_cfrac(iset) = .true.
        !   l_rt_lev(iset) = all(rs%iopts(1:rs%n_instr)%cloud_mode < 2)
        l_rt_lev(iset) = rs%gopts%lev_mode <= 0
        l_lamb  (iset) = any(rs%iopts(1:rs%n_instr)%do_lambertian)
        if (any(rs%iopts(1:rs%n_instr)%use_o3 > 0)) then
          r%n_trg = r%n_trg + 1
          r%i_o3  = r%n_trg
        end if
        if (any(rs%iopts(1:rs%n_instr)%use_co2 > 0)) then
          r%n_trg = r%n_trg + 1
          r%i_co2 = r%n_trg
        end if
        l_btcs (iset) = any(btest(rs%iopts(1:rs%n_instr)%rad_out, OUT_CSB))
        if (l_bc) n_pr_tr(iset) = rs%bc%n_tr
      end if
    end do

    !--------------------------------------------
    ! for non-empty set: allocate + preset arrays
    !--------------------------------------------
    FTRACE_BEGIN('fill_rad:init_sets_loop')
    max_nlev  = 0
    init_sets_loop: do iset = 1, size(set)
      if (present(mask_rs)) then
        if (.not.mask_rs(iset)) CYCLE
      end if
      r  => rad (iset)
      rs => set (iset)
      n_instr = rs%n_instr
      if (dace% lpio .and. l_print) write(6,'(3x,"fill_rad set ",I2," (sat. ",I3.3,", &
           &grid ",I3.3,"): ",I10," spots")') &
           &iset, rs% satid, rs% grid, nset(iset)
      r% n_rec = nset(iset)  ! number of records
      nalloc = r%n_rec
      if (present(nadd)) nalloc = nalloc + nadd(iset)
      if (present(pe_gather)) then
        i = p_sum(nset(iset))
        if (dace% pe == pe_gather) nalloc = i
        i = p_max(nlev(iset))
        if (dace% pe == pe_gather) nlev(iset) = i
      end if
      max_nlev = max(max_nlev,nlev(iset))
      l_pc_aux = l_pc .and. any(rs% instr_wmo(1:n_instr) == 221 &
                           .or. rs% instr_wmo(1:n_instr) == 620 )

      if (nalloc > 0 .or. dace% pe==pe_gather_loc) then
        !-----------------------
        ! set meta data in 'rad'
        !-----------------------
        ! allocate (r% i% chan (set(iset)% n_chan)) !????
        ! r% i     = set (iset)
        r% n_lev = nlev(iset)
        if (l_rttov) then
          r% p_unit   = 'hpa'
        else
          r% p_unit   = 'pa'
        end if
!        r% i% n_pc_emiss = npc
        if (l_pc_aux) then
          r% n_pc_emiss = npc
        else
          r% n_pc_emiss = 0
        end if
        r% ideb = -1
      end if

      if (nalloc > 0) then
        !---------
        ! allocate
        !---------
        if (l_meta) then
          if (present(pe_gather)) then
            allocate(r% pe       (                       nalloc))
            allocate(r% ind      (                       nalloc))
          end if
          allocate (r% obsnum    (                       nalloc))
          allocate (r% i_box     (                       nalloc))
          allocate (r% i_reprt   (                       nalloc))
          allocate (r% dlon      (                       nalloc))
          allocate (r% dlat      (                       nalloc))
          if (l_rt_lev(iset)) then
            allocate (r% p       (             r% n_lev,      1))
          else
            allocate (r% p       (             r% n_lev, nalloc))
          end if
          if (l_mon) then
            allocate (r% time(                           nalloc))
          end if
          allocate (r% nwc_flg   (                       nalloc))
          allocate (r% scanl     (                       nalloc))
          allocate (r% fov       (                       nalloc))
        end if
        if (l_x) then
          allocate (r% t_fg        (           r% n_lev, nalloc))
          allocate (r% q_fg        (           r% n_lev, nalloc))
          allocate (r% ps_fg       (                     nalloc))
          allocate (r% gp_fg       (                     nalloc))
          allocate (r% ts_fg       (                     nalloc))
          allocate (r% v10_abs_fg  (                     nalloc))
          allocate (r% snw_frc     (                     nalloc))
        endif
        if (l_rttov .or. l_pc) then
          allocate (r% emiss       (r%i% n_chan,           nalloc))
        end if
        if (l_pc_aux) then
          allocate(r% emis_pc    (r% n_pc_emiss,         nalloc))
        end if
        if (l_btcs(iset)) allocate(r%bt_fg_cs (r%i%n_chan, nalloc))
        if (.not.l_mon .or. l_H .or. l_pc .or. l_HBH) then
          allocate (r% valid     (r%i% n_chan,           nalloc))
        end if
        if (.not.l_mon .or. l_rttov) then
          allocate (r% stzen     (                       nalloc))
        end if
        if (.not.l_mon) then
          allocate (r% mdlsfc    (                       nalloc))
          allocate (r% op_na     (                       nalloc))
          allocate (r% i_body    (r%i% n_chan,           nalloc))
          allocate (r% bt_obs    (r%i% n_chan,           nalloc))
          allocate (r% bt_bcor   (r%i% n_chan,           nalloc))
          allocate (r% bcor_     (r%i% n_chan,           nalloc))
          allocate (r% bt_fg     (r%i% n_chan,           nalloc))
          allocate (r% not_rej   (r%i% n_chan,           nalloc))
          allocate (r% cloudy    (r%i% n_chan,           nalloc))
          allocate (r% state     (r%i% n_chan,           nalloc))
          allocate (r% flags     (r%i% n_chan,           nalloc))
          allocate (r% r_state   (                       nalloc))
        else
          allocate (r% tsm_fg      (                     nalloc))
        end if
        if (l_err) then
          allocate(r% e_fg       (r%i% n_chan,           nalloc))
          allocate(r% e_obs      (r%i% n_chan,           nalloc))
        end if
        if (present(tj_set)) allocate(tj_set(iset)% a(nalloc))
        if (l_rttov) then
          if (l_x) then
            allocate (r% t2m     (                       nalloc))
            allocate (r% q2m     (                       nalloc))
            allocate (r% u10_fg  (                       nalloc))
            allocate (r% v10_fg  (                       nalloc))
          endif
          allocate (r% stazi     (                       nalloc))
          allocate (r% sunzen    (                       nalloc))
          allocate (r% sunazi    (                       nalloc))
          allocate (r% shgt      (r% i% n_sens,          nalloc))
          allocate (r% stype     (r% i% n_sens,          nalloc))
          if (l_lamb(iset)) allocate (r% specularity(r%i% n_chan,nalloc))
        end if
        if (n_pr_tr(iset) > 0) allocate(r%tr(n_pr_tr(iset), r% i% n_chan, nalloc))
        if (l_bc) then
          allocate(r%plev        (r% i% n_chan,          nalloc))
          allocate(r%orb_phase   (                       nalloc))
          allocate(r%instr_temp  (1,                     nalloc))
        end if
        if (l_rttov.or.l_mon) then
          if(l_x .and. i_cloud(iset) /= 0) then
            if (mw_instr(rs%instr(1))) then
              allocate (r% clw      (r% n_lev,           nalloc))
            elseif (i_cloud(iset) > 0) then
              allocate (r% clwde    (r% n_lev-1,             nalloc))
              allocate (r% icede    (r% n_lev-1,             nalloc))
              allocate (r% cfrac    (r% n_lev-1,             nalloc))
              allocate (r% cld_frc  (r% n_lev-1,             nalloc))
              allocate (r% cld_fg   (6,          r% n_lev-1, nalloc))
            elseif (i_cloud(iset) < 0) then
              allocate (r% cld_top  (                        nalloc))
              allocate (r% cfraction(                        nalloc))
            end if
          endif
          if (l_x .and. l_im_cfrac(iset)) then
            allocate (r% cfraction(                      nalloc))
          end if
          if (r%n_trg > 0) then
            allocate (r% trg      (r%n_trg,    r% n_lev, nalloc))
          end if
        end if
        if (l_H) then
          allocate (r% H_t       (r%i% n_chan, r% n_lev, nalloc))
          allocate (r% H_q       (r%i% n_chan, r% n_lev, nalloc))
          allocate (r% H_ts      (r%i% n_chan,           nalloc))
          allocate (r% H_ps      (r%i% n_chan,           nalloc))
          if (l_pc_aux) then
            allocate (r% H_em_pc (r%i% n_chan, r% n_pc_emiss, nalloc))
          endif
          if (l_K) then
            allocate (r% K_t     (r%i% n_chan, r% n_lev, nalloc))
            allocate (r% K_q     (r%i% n_chan, r% n_lev, nalloc))
            allocate (r% K_ts    (r%i% n_chan,           nalloc))
            allocate (r% K_ps    (r%i% n_chan,           nalloc))
          end if
        end if
        if (l_B) then
          allocate (r% B_tt      (r%   n_lev,  r% n_lev, nalloc))
          allocate (r% B_qq      (r%   n_lev,  r% n_lev, nalloc))
          allocate (r% B_tq      (r%   n_lev,  r% n_lev, nalloc))
        endif
        if (l_HBH) then
          allocate (r% HBHR      (r%i% n_chan, r%i% n_chan, nalloc))
          allocate (r% R         (r%i% n_chan, r%i% n_chan, nalloc))
        endif
        if (l_eb) then
          allocate (r% t_eb      (             r% n_lev, nalloc))
          allocate (r% q_eb      (             r% n_lev, nalloc))
          allocate (r% ps_eb     (                       nalloc))
          allocate (r% ts_eb     (                       nalloc))
          allocate (r% cld_top_eb(                       nalloc))
          allocate (r% cld_frc_eb(             1,        nalloc))
          allocate (r% snw_frc_eb(                       nalloc))
          if (l_pc_aux) then
            allocate(r%emis_pc_eb(r% n_pc_emiss,         nalloc))
          end if
        end if

        !---------------------------
        ! preset with invalid values
        !---------------------------
        if (l_meta) then
          if (present(pe_gather)) then
            r% pe       = 0
            r% ind      = 0
          end if
          r% i_reprt    = 0
          r% i_box      = 0
          r% dlon       = rinvalid
          r% dlat       = rinvalid
          if (l_rt_lev(iset)) then
            r% p(:,1) = preslev    ! pressure levels (Pa) from rttvi
            if (l_rttov) r% p = r% p * 0.01_wp  !Pa -> hPa for RTTOV
          end if
          if (l_mon) then
            r% time     = 0
          end if
          r% nwc_flg  = 0
          r% scanl    = 0
          r% fov      = 0
        end if
        if (.not.l_mon .or. l_H .or. l_pc .or. l_HBH) then
          r% valid      = .false.
        end if
        if (l_x) then
          r% t_fg         = rinvalid
          r% q_fg         = rinvalid
          r% ps_fg        = rinvalid
          r% gp_fg        = rinvalid
          r% ts_fg        = rinvalid
          r% v10_abs_fg   = rinvalid
          r% snw_frc      = rinvalid
        endif
        if (l_rttov .or. l_pc) then
          r% emiss      = rinvalid
        end if
        if (l_pc_aux) then
          r% emis_pc    = rinvalid
        end if
        if (.not.l_mon .or. l_rttov) then
          r% stzen      = rinvalid
        end if
        if (l_btcs(iset)) r% bt_fg_cs = rinvalid
        if (.not.l_mon) then
          r% i_body     = 0
          r% bt_obs     = rinvalid
          r% bt_bcor    = rinvalid
          r% bcor_      = rinvalid
          r% bt_fg      = rinvalid
          r% not_rej    = .false.
          r% cloudy     = 0
          r% state      = STAT_DEFAULT
          r% flags      = 0
          r% r_state    = STAT_DEFAULT
          r% mdlsfc     = 0
          r% op_na      = 0
        else
          r% tsm_fg       = rinvalid
        end if
        if (l_err) then
          r% e_fg       = rinvalid
          r% e_obs      = rinvalid
        end if
        if (l_rttov) then
          if (l_x) then
            r% t2m      = rinvalid
            r% q2m      = rinvalid
            r% u10_fg   = rinvalid
            r% v10_fg   = rinvalid
          endif
          r% shgt       = rinvalid
          r% stype      = rinvalid
          r% stazi      = rinvalid
          r% sunzen     = rinvalid
          r% sunazi     = rinvalid
          if (l_lamb(iset)) r%specularity = rinvalid
        end if
        if (n_pr_tr(iset) > 0) r%tr        = rinvalid
        if (l_bc) then
          r%plev           = rinvalid
          r%orb_phase      = rinvalid
          r%instr_temp     = rinvalid
        end if
        if (l_rttov.or.l_mon) then
          if(l_x .and. i_cloud(iset) /= 0) then
            if (mw_instr(rs%instr(1))) then
              r% clw       = rinvalid
            elseif (i_cloud(iset) > 0) then
              r% clwde     = rinvalid
              r% icede     = rinvalid
              r% cfrac     = rinvalid
              r% cld_frc   = rinvalid
              r% cld_fg    = 0.0_wp
            elseif (i_cloud(iset) < 0) then
              r% cld_top   = rinvalid
              r% cfraction = rinvalid
            end if
          endif
          if (l_x .and. l_im_cfrac(iset)) then
            r% cfraction = rinvalid
          end if
          if (r%n_trg > 0) then
            r% trg       = rinvalid
          end if
        end if
        if (l_H) then
          r% H_t        = rinvalid
          r% H_q        = rinvalid
          r% H_ts       = rinvalid
          r% H_ps       = rinvalid
          if (l_pc_aux) then
            r% H_em_pc  = rinvalid
          endif
          if (l_K) then
            r% K_t      = rinvalid
            r% K_q      = rinvalid
            r% K_ts     = rinvalid
            r% K_ps     = rinvalid
          end if
        end if
        if (l_B) then
          r% B_tt       = rinvalid
          r% B_qq       = rinvalid
          r% B_tq       = rinvalid
        endif
        if (l_HBH) then
          r% HBHR       = rinvalid
          r% R          = rinvalid
        endif
        if (l_eb) then
          r% t_eb       = rinvalid
          r% q_eb       = rinvalid
          r% ps_eb      = rinvalid
          r% ts_eb      = rinvalid
          r% cld_top_eb = rinvalid
          r% cld_frc_eb = rinvalid
          r% snw_frc_eb = rinvalid
          if (l_pc_aux) then
            r% emis_pc_eb = rinvalid
          end if
        end if
      endif
    end do init_sets_loop
    FTRACE_END('fill_rad:init_sets_loop')

    FTRACE_BEGIN('fill_rad:alloc_aux')

    npv  = 2 * max_nlev               ! number of profile variables
    ncv  = npv + nsv + npc            ! number of input variables
    if (l_H) then
      allocate (H (m_chan, ncv), tmp (max (m_chan, ncv)))
      if (l_K) then
        allocate (HBHRi(m_chan,m_chan), Pb(ncv,ncv), K(m_chan, ncv))
        call construct(v, x% info)
      end if
    end if
    if (l_B) then
      allocate(B(ncv, ncv))
    endif

    if (l_H .or. l_B) then
      i = maxval(nlev(:))
      allocate(jk(i+1,2,2),kj(i+1,2,2))
    end if

    FTRACE_END('fill_rad:alloc_aux')

    FTRACE_BEGIN('fill_rad:box_loop')
    !-------------------------
    ! fill derived type t_radv
    !-------------------------
    nset = 0
!NEC$ nomove
    box_loop: do ib = ib_s, ib_e              ! loop over observation 'boxes'
      if (obs% o(ib)% pe /= dace% pe) cycle   ! only on this PE
      call obs_block (ob, obs, ib)
      if (l_K .or. l_HBH) then
        HBHRb = HBHR% b(ib,ib)
        call full_matrix(HBHRb)
      end if
      if (l_HBH) then
        Rb   = Roo% b(ib,ib)
        call full_matrix(Rb)
      end if

!NEC$ nomove
      spot_loop: do is=1,obs% o(ib)% n_spot   ! loop over reports/fovs
        s => obs% o(ib)%spot(is)              !
        if (s% hd% obstype /= OT_RAD) cycle   ! only treat radiances
        if (s% use% state  <= STAT_DISMISS) cycle
        if (present(min_stat)) then
          if (s% use% state  <= min_stat) cycle
        end if

        !------------------------------------------
        ! chose combination of satellite/instrument
        !------------------------------------------
        satid = s% hd% satid                  ! satellite id
        grid  = s% hd% grid_id                ! grid id
        iset  = set_indx (set, satid=satid, grid=grid)
        if (iset < 1 .and. l_partial) cycle spot_loop
        if (present(mask_rs)) then
          if (.not.mask_rs(iset)) cycle spot_loop
        end if
        nset (iset) = nset (iset) + 1
        j = nset (iset)
        !-------------
        ! spot indices
        !-------------
        i0 = s% i% i
        i1 = i0 + 1
        in = s% i% n
        ip = npv+nsv+1
        !------------
        ! fill header
        !------------
        r  => rad (iset)
        rs => set (iset)
        l_pc_aux = l_pc .and. (r% n_pc_emiss > 0)
        npv  = 2 * r%n_lev                ! number of profile variables
        ncv  = npv + nsv + npc            ! number of input variables

        if (l_rttov .or. l_mon .or..not.l_rt_lev(iset) .or. &
            l_meta  .or. n_pr_tr(iset) > 0 .or. l_bc) call load (ob% o, s, tovs)

        if (ldeb(s)) then
          ndeb = count(r%ideb > 0)
          if (ndeb < size(r%ideb)) r%ideb(ndeb+1) = j
        end if

        if (l_meta) then
          if (present(pe_gather)) then
            r% pe    (j) = dace%pe
            r% ind   (j) = j
          end if
          r% obsnum  (j) = s% hd% id
          r% i_box   (j) = ib
          r% i_reprt (j) = is
          r% dlon    (j) = s% col% c% dlon
          r% dlat    (j) = s% col% c% dlat
          if (l_mon) then
            time = s% hd% time - ana_time
            r% time(j) = real(time%days * 86400 + time%secs,wp) / 60._wp
          end if
          if (.not.l_rt_lev(iset)) then
            r% p(:,j) = tovs% av(:,tovs%i_p)
            if (l_rttov) r% p(:,j) = r% p(:,j) * 0.01_wp  !Pa -> hPa for RTTOV
          end if
          r% nwc_flg (j) = tovs% nwc_flg
          r% fov     (j) = s% phase
          r% scanl   (j) = tovs% scanl
        end if
        if (.not.l_mon .or. l_rttov) then
          r% stzen   (j) = s% stzen
        end if
        if (l_bc) then
          r%orb_phase(j)    = tovs% orb_phase
          r%instr_temp(1,j) = tovs% instr_temp
        end if

        if (.not.l_mon) then
          r% mdlsfc  (j) = s% mdlsfc
          if (s%o%n > 0) then
            r% op_na (j) = obs% o(ib)% body(s%o%i+1)% op_na
          else
            r% op_na (j) = 0
          end if
        else
          r% tsm_fg  (j) =   s%ts_bg
        end if
        !--------------------------------------------------------------
        ! transform variables: 3dvar interpolation space -> rttov input
        !--------------------------------------------------------------
        if (present(tj_set)) then
          tjp => tj_set(iset)% a(j)
        elseif (l_H .or. l_eb) then
          tjp => tj
        else
          tjp => null()
        end if

        if (l_x) then
          if (l_rttov) then

            if (present(tj_set) .or. l_H .or. l_eb .or. l_B) then
              call prep_rttov_prof(ob, s, x% s(ib), rtp=rtp, ttovs=tovs, tjac=tjp)
            else
              call prep_rttov_prof(ob, s, x% s(ib), rtp=rtp, ttovs=tovs)
            endif
          else
            if (present(tj_set) .or. l_H .or. l_eb .or. l_B) then
              call prep_rttov_prof(ob, s, x% s(ib), rtp=rtp, tjac=tjp)
            else
              call prep_rttov_prof(ob, s, x% s(ib), rtp=rtp)
            endif
          end if
        end if

        if (l_H .or. l_B) then
          !--------------------------------------------------------------------
          ! calculate transformation matrices for DACE <-> RTTOV representation
          ! DACE:  virtual temperature, generalised (relative) humidity
          ! RTTOV: temperature, specific humidity
          !--------------------------------------------------------------------
          on = s% o% n
          o0 = s% o% i
          o1 = o0 + 1
          ox = o0 + on
          jk (1:r%n_lev,1,1) = tjp% dt_tv(:)     ! dt_tv
          jk (1:r%n_lev,1,2) = tjp% dt_rh(:)     ! dt_rh
          jk (1:r%n_lev,2,1) = tjp% dq_tv(:)     ! dq_tv
          jk (1:r%n_lev,2,2) = tjp% dq_rh(:)     ! dq_rh
          jk (1+r%n_lev,1,1) = tjp% sigma_var_tskin
          jk (1+r%n_lev,1,2) = 0._wp
          jk (1+r%n_lev,2,1) = 0._wp
          jk (1+r%n_lev,2,2) = s% pz_bg
          kj  = 0._wp
          do i = 1, size (jk,1)
            edet = jk(i,1,1) * jk(i,2,2) - jk(i,2,1) * jk(i,1,2)
            if (edet /= 0._wp) then
              edet       = 1._wp/edet
              kj (i,1,1) =  edet * jk (i,2,2)
              kj (i,1,2) = -edet * jk (i,1,2)
              kj (i,2,1) = -edet * jk (i,2,1)
              kj (i,2,2) =  edet * jk (i,1,1)
            else
              if (jk (i,1,1) /= 0._wp) kj (i,1,1) = 1._wp / jk (i,1,1)
              if (jk (i,2,2) /= 0._wp) kj (i,2,2) = 1._wp / jk (i,2,2)
            endif
          end do
        endif

        if (l_H) then
          Hb = sub_block (ob%H, s%o%i, s%o%n, s%i%i, s%i%n, repr=FULL)
          H(1:on,1:in) = real(Hb%full(1:on,1:in), kind=wp)
          ! H (1:on   ,in+1:) = rinvalid  ! entries not used in this fof
          if (l_K) then
            call construct(rc, i0 + ncv, 1) !??????????
            call init_spot_obs (rc, s, ob% o)
            call get_corr(rc, rc, 1, 1, .true., .true., c=Pb)
            call destruct(rc)
            HBHRi (1:on,1:on) = inverse_rs(real(HBHRb% full(o1:ox,o1:ox),wp))
            do i=1,on
              HBHRi(i,1:on) = HBHRi(i,1:on) * rhs% s(ib)% x(o0+i)
            end do
            K (1:on,1:in) = matmul(HBHRi(1:on,1:on), H (1:on ,1:in))
            K (1:on,1:in) = matmul(K    (1:on,1:in), Pb(1:in,1:in))
          end if
          do i = 1, npv+2, 2
            tmp(1:on) = H(1:on,i)
            H(1:on,i  ) = kj(i/2+1,1,1) * tmp(1:on) + kj(i/2+1,2,1) * H(1:on,i+1)
            H(1:on,i+1) = kj(i/2+1,1,2) * tmp(1:on) + kj(i/2+1,2,2) * H(1:on,i+1)
          end do
        end if    ! l_H
        if (l_B) then
          !---------------------------------------------------------
          ! B in 'interpolation space' is written 'as used in RTTOV'
          ! i.e. in terms of t,q not tv,rh.
          !---------------------------------------------------------
          Bb = sub_block(ob%Bi, s%i%i, s%i%n, s%i%i, s%i%n, repr=FULL)
          B(1:in,1:in) = real(Bb%full(1:in,1:in), kind=wp)
          do i = 1, npv+2, 2
            tmp(1:in) = B(1:in,i)
            B(1:in,i  ) = jk(i/2+1,1,1) * tmp(1:in) + jk(i/2+1,1,2) * B(1:in,i+1)
            B(1:in,i+1) = jk(i/2+1,2,1) * tmp(1:in) + jk(i/2+1,2,2) * B(1:in,i+1)
          end do

          do i = 1, npv+2, 2
            tmp(1:in) = B(i,1:in)
            B(i  ,1:in) = jk(i/2+1,1,1) * tmp(1:in) + jk(i/2+1,1,2) * B(i+1,1:in)
            B(i+1,1:in) = jk(i/2+1,2,1) * tmp(1:in) + jk(i/2+1,2,2) * B(i+1,1:in)
          end do
        endif

        if (.not.l_mon .or. l_rttov .or. l_H .or. l_pc .or. l_err .or. l_HBH) then
          !----------
          ! fill body
          !----------
!NEC$ nomove
          do i = 1, s% o% n
            m = i + s% o% i
            ! Determine channel/instr index
            ! TODO: replace by tovs%ci?
            ichan =              nint(obs% o(ib)% olev (m))
            instr = rttov_instr (int (obs% o(ib)% body (m)% lev_sig), satid)
            l = chan_indx (instr, ichan, r% i, ii=ii)
            if (l <= 0) then
              write(0,*) i,m,r%i%satid, r%i%grid, l_mon, l_H, l_pc, dace% pe
              write(0,*) is,s%hd%id
              write(0,*) ichan, instr, int (obs% o(ib)% body (m)% lev_sig)
              write(0,*) l,instr,ichan
              do l = 1, r%i%n_instr
                write(0,*) 'instr',l,r%i%instr(l),instr
                if (r%i%instr(l)==instr) then
                  do m = r%i%o_ch_i(l)+1,r%i%o_ch_i(l)+r%i%n_ch_i(l)
                    write(0,*) 'chan',m,r%i%chan(m),ichan
                  end do
                end if
              end do
              write(0,*)
              call print_rad_set(r%i,hint=dace% pe,header='crash', unit=stderr)
              write(0,*)
              write(0,*) 'fill_rad: Failed to get channel index'
              write(0,*) '          instr, ichan =', instr, ichan
              call finish('fill_rad','Failed to get channel index')
            end if
            if (.not.l_mon .or. l_H .or. l_pc .or. l_HBH) then
              r% valid   (l,j) = .true.
            end if
            if (.not.l_mon) then
              r% bcor_   (l,j) = - obs% o(ib)% body(m)% bc
              r% bt_bcor (l,j) = obs% o(ib)% body(m)% o
              r% bt_obs  (l,j) = r% bt_bcor(l,j) + r% bcor_(l,j)
              r% not_rej (l,j) = valid (obs% o(ib)% body(m)% use% state)
              r% state   (l,j) = obs% o(ib)% body(m)% use% state
              r% flags   (l,j) = obs% o(ib)% body(m)% use% flags
              r% r_state (  j) = s% use% state
              r% i_body  (l,j) = m
              if (present(y)) r% bt_fg   (l,j) = y% s(ib)% x(m)
              if (btest  (obs% o(ib)% body(m)% use% flags, CHK_CLOUD)) then
                r% cloudy(l,j) = 1
              else
                r% cloudy(l,j) = 2
              endif
            end if
            if (l_err) then
              r% e_fg    (l,j) = obs% o(ib)% body(m)% eb
              r% e_obs   (l,j) = obs% o(ib)% body(m)% eo
            end if
            if (l_pc .or. l_rttov) then
              if (l_x) then
                ! rtp available
                r% emiss   (l,j) = rtp% emis (1,i)
              else
                r% emiss   (l,j) = tovs% emis(i)
              end if
            endif
            if (l_rttov .and. l_lamb(iset)) then
              if (rs%iopts(ii)%do_lambertian) then
                m = count(rs%iopts(1:ii)%do_lambertian)
                r% specularity(l,j) = max(min(real(tovs% spec(m),wp),1._wp),0._wp)
              end if
            end if
            if (n_pr_tr(iset) > 0) then
              r%tr(1:tovs%ntr,l,j) = tovs% tr(i,1:tovs%ntr)
              if (ldeb(s)) write(usd,*) dpref,' fill_rad tr',tovs% tr(i,1:tovs%ntr)
            end if
            if (l_bc) r%plev(l,j)  = obs% o(ib)% body(m)% plev
            if (l_btcs(iset) ) r% bt_fg_cs(l,j) = tovs% bt_cs(i)
            if (l_H) then
              r% H_t (l,:,j) = H(i, 1:npv:2)
              r% H_ts(l,  j) = H(i, npv+1)
              r% H_q (l,:,j) = H(i, 2:npv:2)
              r% H_ps(l,  j) = H(i, npv+2)
              if (l_pc_aux) then
                r% H_em_pc(l,:,j) = H(i, npv+6:)
              endif
              if (l_K) then
                v% s(ib)% x(i1:i0+in) = K(i,1:in) + x% s(ib)% x(i1:i0+in)
                call prep_rttov_prof(ob, s, v% s(ib), rtp=rtp_k)
                r% K_t (l,:,j) = rtp_k% av (:,1,1) - rtp% av (:,1,1)            ! T-FG
                r% K_q (l,:,j) = rtp_k% av (:,2,1) - rtp% av (:,2,1)            ! Q-FG
                r% K_ts(l,  j) = rtp_k% ssv(  1,1) - rtp% ssv(  1,1)            ! Tsurf
                r% K_ps(l,  j) =(rtp_k% sav(  3,1) - rtp% sav(  3,1)) * 100._wp ! Psurf (hPa)
                call destruct(rtp_k)
              end if
            end if
            if (l_HBH) then
              do i2 = 1, s% o% n
                m2 = i2 + s% o% i
                ! Determine channel/instr index
                ! TODO: replace by tovs%ci?
                ichan =              nint(obs% o(ib)% olev (m2))
                instr = rttov_instr (int (obs% o(ib)% body (m2)% lev_sig), satid)
                l2 = chan_indx (instr, ichan, r% i)
                !-----------------------------------------
                r% HBHR(l,l2,j) = real(HBHRb% full(m,m2),wp)
                r% R   (l,l2,j) = real(Rb  % full(m,m2),wp)
              enddo
            endif
          end do
          if (l_B) then
            r% B_tt(:,:,j) = B(1:npv:2, 1:npv:2)
            r% B_qq(:,:,j) = B(2:npv:2, 2:npv:2)
            r% B_tq(:,:,j) = B(1:npv:2, 2:npv:2)
          endif
        end if

        !--------------
        ! fill profiles
        !--------------
        if (l_x) then
          if (l_rttov .and. l_rt_lev(iset)) call rttov_bounds (rtp)

          r% t_fg  (:,j)    =      rtp% av  (:,1,1)            ! T-FG profile
          r% q_fg  (:,j)    =      rtp% av  (:,2,1)            ! Q-FG profile
          r% ts_fg   (j)    =      rtp% ssv (  1,1)            ! Tsurf
          r% ps_fg   (j)    =      rtp% sav (  3,1)            ! Psurf (hPa)
          r% gp_fg   (j)    =      s%gp_bg
          r% v10_abs_fg(j)  = sqrt(rtp% sav (  4,1) ** 2 + &
                                   rtp% sav (  5,1) ** 2   )   ! 10 m wind
          r% snw_frc (j)    =      rtp% snf
        endif

        if (l_rttov .and. l_x .and. i_cloud(iset) > 0) then
          ! Compute cloud parameters
          ! RTTOV cloud - cloud water/ice content on layers
          ! RTTOV cloud(1:6,nlayers): 6 different cloud types
          !                           1-5 five water cloud types
          !                           6   ice cloud
          ! Water cloud type 1 (2-5 not used)
          ! Scale cld_fg with cloud fraction within grid box (RTTOV12 and older)
          ! Cloud liquid /  ice water particle effective diameter in micro meters
          if (mw_instr(rs%instr(1))) then
            do m = 1, tovs%nav
              if (tovs%av_cont(m) == COL_QCDIA) then
                r%clw(:,j) = real(tovs%av(:,m),wp)
                exit
              end if
            end do
          else
            call cloud_features (            &
                 tovs,                       & ! <- tovs
                 r%p(:,min(j,ubound(r%p,2))),& ! <- p
                 r%t_fg(:,j),                & ! <- t
                 r%q_fg(:,j),                & ! <- q
                 scale_qi,                   & ! <- scaling of lwc/iwc if deff is clipped
                 r%cld_fg (1,:,j),           & ! <- cloud water content
                 r%cld_fg (6,:,j),           & ! <- cloud ice content
                 r % icede(:,j),             & ! <- effective radius of cloud water
                 r % clwde(:,j),             & ! <- effective radius of cloud ice
                 r % cfrac(:,j))
          end if
        end if

        if (l_pc_aux) then
          r% emis_pc(:in-ip+1,j) = x% s(ib)% x(i0+ip:i0+in)
        end if

        if (l_rttov) then
          if (l_x) then
            r% q_fg   (:,j) = q2rt_humi(r% q_fg(:,j))
            r% q2m      (j) = q2rt_humi(rtp% sav(2,1))
            r% t2m      (j) = rtp% sav(1,1)              ! 2m temp
            r% u10_fg   (j) = rtp% sav(4,1)              ! 10 m zonal wind
            r% v10_fg   (j) = rtp% sav(5,1)              ! 10 m meridional wind
          endif
          r% shgt     (:,j) = rtp% hsurf
          r% stype    (:,j) = tovs% rt_stype(1)
          r% stazi      (j) = tovs% boa
          r% sunzen     (j) = tovs% soza
          r% sunazi     (j) = tovs% soa
        else
          if (l_x) then
            r% ps_fg (j) = r% ps_fg (j) * 100._wp !hPa -> Pa
          endif
        end if
        if (l_rttov.or.l_mon) then
          if (l_x .and. i_cloud(iset) /= 0 .and..not.mw_instr(rs%instr(1))) then
            if (i_cloud(iset) > 0) then
              r% cld_frc(  1,j) = rtp% cv (2,1)              ! cloud fraction
            elseif (i_cloud(iset) < 0) then
              r% cld_top    (j) = rtp% cv (1,1)
              if (r% cld_top(j) > 0._wp) then
                r% cld_top(j) = r% cld_top(j) * 0.01_wp    ! cloud top height (hPa)
              else
                r% cld_top(j) = 1200._wp
              end if

              r% cfraction  (j) = rtp% cv (2,1) * 0.01_wp             ! simple cloud fraction
              if ( r% cfraction  (j) <= 0._wp .and. r% cld_top(j) < r% ps_fg (j)) then
                r% cfraction  (j) = 1._wp
              end if
              r% cfraction(j) = min(r%cfraction(j), 0.98_wp)
            end if
          endif
          if (l_x .and. l_im_cfrac(iset)) then
            r% cfraction  (j) = rtp% cv (2,1) * 0.01_wp             ! imager cloud fraction
          end if
          if (r%i_o3  > 0) r% trg(r%i_o3,1:r%n_lev,j) = &
               rt_trg(real(tovs% av(1:r%n_lev,tovs%i_o3),kind=wp), rtp% av(1:r%n_lev,2,1), gas_id_o3)
          if (r%i_co2 > 0) r% trg(r%i_co2,:,j) = &
               rt_trg(real(tovs% av(1:r%n_lev,tovs%i_co2),kind=wp), rtp% av(1:r%n_lev,2,1), gas_id_co2)
        end if

        if (l_eb) then
          r% t_eb        (:,j) = tjp% dt_tv(:)        * e_bg% s(ib)% x(i0+1:i0+npv:2)  ! T-EB profile
          r% q_eb        (:,j) = tjp% dq_rh(:)        * e_bg% s(ib)% x(i0+2:i0+npv:2)  ! Q-EB profile
          r% ts_eb         (j) = tjp% sigma_var_tskin * e_bg% s(ib)% x(i0+npv+1)       ! EB Tsurf
          r% ps_eb         (j) = s% pz_bg             * e_bg% s(ib)% x(i0+npv+2)       ! EB Psurf (hPa)
          r% cld_top_eb    (j) = tjp% dclt_dum        * e_bg% s(ib)% x(i0+npv+3)       ! EB cloud top height (hPa)
          r% cld_frc_eb  (1,j) = tjp% dclf_dum        * e_bg% s(ib)% x(i0+npv+4)       ! EB cloud fraction
          r% snw_frc_eb    (j) = tjp% dsnf_dum        * e_bg% s(ib)% x(i0+npv+5)       ! EB snow fraction
          if (l_pc_aux) then
            r% emis_pc_eb(:in-npv-5,j) =                e_bg% s(ib)% x(i0+npv+6:i0+in) ! EB emissivity PCs
          end if
        end if

        call destruct(rtp)
        call destruct(tovs)
        if (.not.present(tj_set) .and. (l_H .or. l_eb)) call destruct(tj)
        if (l_H) call destruct(Hb)
        if (l_B) call destruct(Bb)

      end do spot_loop

      if (l_K .or. l_HBH) call deallocate(HBHRb)
      if (l_HBH) then
         call deallocate(Rb)
      endif

    end do box_loop
    FTRACE_END('fill_rad:box_loop')

    call print_emis_stat
    call print_tskin_stat

    if (l_H .and. l_K .and. r%n_rec > 0) call destruct(v)

    if (present(pe_gather)) then
      do iset = 1, size(set)
        if (present(mask_rs)) then
          if (.not.mask_rs(iset)) CYCLE
        end if
        r  => rad (iset)
        call gather_radv(r, pe_gather, l_rt_lev(iset))
        if (dace% pe /= pe_gather) call destruct(r)
      end do
    end if
    FTRACE_END('fill_rad:finalize')


    FTRACE_END('fill_rad')
  end subroutine fill_rad

!------------------------------------------------------------------------------

  !--------------------------
  ! subroutine cloud_features
  !--------------------------
  !> Compute the cloud water/ice content and the effective diameter of droplets
  !>
  !> The cloud water content and the cloud ice content are estimated and scaled
  !> with the cloud fraction within each grid box.
  !>
  !> <b> call cloud_features (qc_dia, qi_dia, cfrac, p, t, q, lwc, ic, diameter) </b>
  !>
  !> @param[in]  qc_dia   Specific cloud liquid water content, diagnostic, in kg/kg
  !> @param[in]  qi_dia   Specific cloud ice content, diagnostic, in kg/kg
  !> @param[in]  cfrac    Cloud fraction per layer, 0 <= cfrac <= 1
  !> @param[in]  p        Pressure in hPa
  !> @param[in]  t        Temperature in K
  !> @param[in]  q        specific humidity in kg/kg
  !> @param[out] lwc      cloud water content (liquid water content) in g m**(-3)
  !> @param[out] ic       cloud ice content  in g m**(-3)
  !> @param[out] diameter effective cloud droplet diameter in micro meters
  !>                      (vertical profile), diameter >= 5 micro meter
  !>
  !>
  !> Density of liquid water rho_w:                                          \n
  !> Within the atmosphere rho_w does not change consideribly with pressure  \n
  !> or temperature and is assumed to be constant:                           \n
  !> rho_w = 1000 kg m**(-3)                                                 \n
  !>
  !> Liquid water content LWC                                                \n
  !> lwc = qc_dia * rho_l                                                    \n
  !> with                                                                    \n
  !> lwc    - liquid water content in kg m**(-3)                             \n
  !> qc_dia - specific cloud liquid water content in kg/kg                   \n
  !> rho_l  - density of air in kg m**(-3)                                   \n
  !>
  !> lwc [g m**(-3)] = 1000 * qc_dia [kg/kg] * rho_l [kg m**(-3)]
  !>
  !> Density of air rho_l                                                    \n
  !> rho_l = p / (Rdry * Tv)                                                 \n
  !> with                                                                    \n
  !> rho_l - density of air in kg m**(-3)                                    \n
  !> p     - pressure in Pa                                                  \n
  !> Rdry  - gas constant of dry air: Rdry [J/(kg*K)] = 287.05               \n
  !> Tv    - virtual temperature in K, provided by the function tv_t_q(t,q)  \n
  !>
  !> Pressure in hPa: rho_l = 100 * p / (Rdry * Tv)
  !>
  !>
  !------------------------------------------------------------------------
  subroutine cloud_features (tovs,p,t,q,scale_qi,lwc,iwc,deff_ice,deff_water,cfrac)
    type(t_tovs),intent(in)  :: tovs
    real(wp),    intent(in)  :: p(:)
    real(wp),    intent(in)  :: t(:)
    real(wp),    intent(in)  :: q(:)
    logical,     intent(in)  :: scale_qi
    real(wp),    intent(out) :: lwc(:)
    real(wp),    intent(out) :: iwc(:)
    real(wp),    intent(out) :: deff_ice  (:)
    real(wp),    intent(out) :: deff_water(:)
    real(wp),    intent(out) :: cfrac(:)

    real(wp) :: rho_l(size(p))
    integer  :: idx_qc
    integer  :: idx_qi
    integer  :: idx_rqc
    integer  :: idx_rqi
    integer  :: i, nlay
    real(wp),allocatable :: fac_qi(:)

    nlay = tovs%nlev - 1

    !Cloud variable indices
    idx_qc  = -1
    idx_qi  = -1
    idx_rqc = -1
    idx_rqi = -1
    do i=1, tovs%nav
      select case (tovs%av_cont(i))
      case (COL_CLC)  ! cloud fractional cover on layers = cloud cover
                      ! Model CLC: cloud cover in %, i.e. 0% <= CLC <= 100%
                      ! RTTOV cfrac: cloud cover between 0 and 1, i.e. 0 <= cfrac <= 1.0
                      ! => RTTOV cfrac = 0.01 * Model CLC
        !r% cfrac(:,j) = 0.01_wp * tovs%av(:,i)
        cfrac(1:nlay)  = 0.01_wp * tovs%av(1:nlay,i)
      case (COL_QCDIA)  ! Grid-scale + subgrid-scale cloud water
        idx_qc = i
      case (COL_QIDIA)  ! Grid-scale + subgrid-scale cloud ice
        idx_qi = i
      case (COL_REFF_QC)! Effective radius of cloud water
        idx_rqc = i
      case (COL_REFF_QI)! Effective radius of cloud ice
        idx_rqi = i
      end select
    end do

    !Effective diameter
    if (idx_rqc > 0 .and. idx_rqi > 0) then

      !Physical unit [reff -> deff, m -> micron]
      deff_water(1:nlay) = tovs%av(1:nlay,idx_rqc) * 2._wp * 1.0e+6_wp
      deff_ice  (1:nlay) = tovs%av(1:nlay,idx_rqi) * 2._wp * 1.0e+6_wp

      !Scale water mass according to clipping
      if ( scale_qi ) then
        allocate(fac_qi(nlay))
        fac_qi(:) = 1._wp
        where ( deff_ice < min_dia_qi .and. deff_ice > 0._wp)
          fac_qi = deff_ice/min_dia_qi
        end where
      end if

      !Clip effective diameter
      where ( deff_water < min_dia_qc .and. deff_water > 0._wp)
        deff_water = min_dia_qc
      end where
      where ( deff_ice < min_dia_qi .and. deff_ice > 0._wp)
        deff_ice = min_dia_qi
      end where

    else ! Use RTTOV parameterisation for diameter calculation
      deff_ice   = 0._wp
      deff_water = 0._wp
    end if

    !Liquid water content / ice water content
    if (idx_qc > 0 .and. idx_qi > 0) then
      ! qc_dia =  real(tovs%av(:,idx_qc),wp)
      ! qi_dia =  real(tovs%av(:,idx_qi),wp)

      ! Density of air
      ! => p in hPa
      ! => Rdry in J/(kg*K)
      ! => tv_t_q in K
      ! <= rho_l in kg m**(-3)
      rho_l = 100.0_wp * (p / (Rdry * tv_t_q(t,q)))

      ! Liquid water content of clouds
      ! => QC_DIA in kg/kg
      ! => rho_l in kg m**(-3)
      ! <= lwc in g m**(-3)
      lwc(1:nlay) = 1000.0_wp * real(tovs%av(1:nlay,idx_qc),wp) * rho_l(1:nlay)

#if (_RTTOV_VERSION <= 12)
      ! Scale lwc with cloud fraction within grid box
      where (cfrac(:) > 0.0_wp) lwc = lwc / cfrac
#endif
      ! Ice water content of clouds
      ! => QI_DIA in kg/kg
      ! => rho_l in kg m**(-3)
      ! <= ic in g m**(-3)
      iwc(1:nlay) = 1000.0_wp * real(tovs%av(1:nlay,idx_qi),wp) * rho_l(1:nlay)

#if (_RTTOV_VERSION <= 12)
      ! Scale iwc with cloud fraction within grid box
      where (cfrac > 0.0_wp) iwc = iwc / cfrac
#endif
      ! Scale iwc
      if ( scale_qi) then
        if (allocated(fac_qi)) then
          iwc = iwc/fac_qi
        end if
      end if
    end if
  end subroutine cloud_features

!------------------------------------------------------------------------------

  subroutine gather_radv (r, pe_gather, l_rt_lev)
    type(t_radv), intent(inout)         :: r
    integer,      intent(in)            :: pe_gather
    logical,      intent(in)            :: l_rt_lev

    integer :: n_send, n_rec_all

    n_rec_all = p_sum(r%n_rec)
    if (n_rec_all <= 0) return

    call gather_int    (r%pe        )
    call gather_int    (r%ind       )
    call gather_int    (r%i_box     )
    call gather_int    (r%j_box     )
    call gather_int    (r%i_reprt   )
    call gather_int_2d (r%i_body    )
    call gather_int    (r%obsnum    )
    call gather_int    (r%date      )
    call gather_int    (r%time      )
    call gather_int    (r%date_d    )
    call gather_int    (r%time_d    )
    call gather_int    (r%date_d    )
    call gather_int    (r%ntstep    )
    call gather_real   (r%dlat      )
    call gather_real   (r%dlon      )
    call gather_int    (r%fov       )
    call gather_int    (r%scanl     )
    call gather_real   (r%stzen     )
    call gather_real   (r%stazi     )
    call gather_real   (r%sunzen    )
    call gather_real   (r%sunazi    )
    call gather_real_2d(r%landfr    )
    call gather_int_2d (r%stype     )
    call gather_real_2d(r%shgt      )
    call gather_real_2d(r%specularity)
    call gather_int    (r%center    )
    call gather_int    (r%subcenter )
    ! call gather_real_2d(r%pred      )
    call gather_int    (r%mdlsfc    )
    call gather_int    (r%op_na     )
    call gather_log_2d (r%not_rej   )
    call gather_int_2d (r%cloudy    )
    call gather_int_2d (r%state     )
    call gather_int_2d (r%flags     )
    call gather_int    (r%r_state   )
    call gather_log_2d (r%valid     )
    call gather_real_2d(r%sinfl     )
    call gather_real_2d(r%bt_fg     )
    call gather_real_2d(r%bt_obs    )
    call gather_real_2d(r%bt_bcor   )
    call gather_real_2d(r%bcor_     )
    call gather_real_2d(r%emiss     )
    call gather_real_2d(r%bt_fg_cs  )
!    if (r% i% n_pc_emiss > 0) then
    if (n_pc > 0) then
      call gather_real_2d(r%emis_pc )
    end if
    if (l_rt_lev) then
      if (dace% pe == pe_gather .and. .not.associated(r% p)) then
        allocate(r% p(r%n_lev,1))
        r% p(:,1) = preslev
      end if
    else
      call gather_real_2d(r%p       )
    end if
    call gather_real_2d(r%h_fg      )
    call gather_real_2d(r%t_fg      )
    call gather_real_2d(r%q_fg      )
    call gather_real   (r%t2m       )
    call gather_real   (r%q2m       )
    call gather_real   (r%ps_fg     )
    call gather_real   (r%gp_fg     )
    call gather_real   (r%ts_fg     )
    call gather_real   (r%tsm_fg    )
    call gather_real   (r%u10_fg    )
    call gather_real   (r%v10_fg    )
    call gather_real   (r%v10_abs_fg)
    call gather_real   (r%snw_frc   )
    call gather_real   (r%cld_top   )
    call gather_real   (r%cfraction )
    call gather_real_2d(r%clwde     )
    call gather_real_2d(r%icede     )
    call gather_real_2d(r%cfrac     )
    call gather_real_3d(r%cld_fg    )
    call gather_real_2d(r%cld_frc   )
    call gather_real_3d(r%trg       )
    call gather_real_3d(r%H_t       )
    call gather_real_3d(r%H_q       )
    call gather_real_2d(r%H_ts      )
    call gather_real_2d(r%H_ps      )
    call gather_real_3d(r%K_t       )
    call gather_real_3d(r%K_q       )
    call gather_real_2d(r%K_ts      )
    call gather_real_2d(r%K_ps      )
    call gather_real_3d(r%B_tt      )
    call gather_real_3d(r%B_qq      )
    call gather_real_3d(r%B_tq      )
    call gather_real_3d(r%HBHR      )
    call gather_real_3d(r%R         )
    call gather_real_2d(r%t_eb      )
    call gather_real_2d(r%q_eb      )
    call gather_real   (r%ps_eb     )
    call gather_real   (r%ts_eb     )
    call gather_real   (r%cld_top_eb)
    call gather_real_2d(r%cld_frc_eb)
    call gather_real   (r%snw_frc_eb)
    if (n_pc > 0) then
      call gather_real_2d(r%emis_pc_eb)
      call gather_real_3d(r%H_em_pc   )
    endif

    if (dace% pe == pe_gather) r%n_rec = n_rec_all

  contains

    subroutine gather_int(dat)
      integer, pointer :: dat(:)
      integer, target  :: dum(0)
      integer, pointer :: send(:), recv(:)
      n_send = r%n_rec
      if (dace% pe == pe_gather) n_send = 0
      send => dum
      recv => dum
      if (associated(dat)) then
        if (n_send > 0)          send => dat(1:n_send)
        if (size(dat) > r%n_rec) recv => dat(r%n_rec+1:)
      else
        n_send =  0
      end if
      call p_gatherv(send, recv, pe_gather)
    end subroutine gather_int

    subroutine gather_int_2d(dat)
      integer, pointer :: dat(:,:)
      integer, target  :: dum(0,0)
      integer, pointer :: send(:,:), recv(:,:)
      n_send = r%n_rec
      if (dace% pe == pe_gather) n_send = 0
      send => dum
      recv => dum
      if (associated(dat)) then
        if (n_send > 0)            send => dat(:,1:n_send)
        if (size(dat,2) > r%n_rec) recv => dat(:,r%n_rec+1:)
      else
        n_send =  0
      end if
      call p_gatherv(send, recv, pe_gather)
    end subroutine gather_int_2d

    subroutine gather_real(dat)
      real(wp), pointer :: dat(:)
      real(wp), target  :: dum(0)
      real(wp), pointer :: send(:), recv(:)
      n_send = r%n_rec
      if (dace% pe == pe_gather) n_send = 0
      send => dum
      recv => dum
      if (associated(dat)) then
        if (n_send > 0)          send => dat(1:n_send)
        if (size(dat) > r%n_rec) recv => dat(r%n_rec+1:)
      else
        n_send =  0
      end if
      call p_gatherv(send, recv, pe_gather)
    end subroutine gather_real

    subroutine gather_real_2d(dat)
      real(wp), pointer :: dat(:,:)
      real(wp), target  :: dum(0,0)
      real(wp), pointer :: send(:,:), recv(:,:)
      n_send = r%n_rec
      if (dace% pe == pe_gather) n_send = 0
      send => dum
      recv => dum
      if (associated(dat)) then
        if (n_send > 0)            send => dat(:,1:n_send)
        if (size(dat,2) > r%n_rec) recv => dat(:,r%n_rec+1:)
      else
        n_send =  0
      end if
      call p_gatherv(send, recv, pe_gather)
    end subroutine gather_real_2d

    subroutine gather_real_3d(dat)
      real(wp), pointer :: dat(:,:,:)
      real(wp), target  :: dum(0,0,0)
      real(wp), pointer :: send(:,:,:), recv(:,:,:)
      n_send = r%n_rec
      if (dace% pe == pe_gather) n_send = 0
      send => dum
      recv => dum
      if (associated(dat)) then
        if (n_send > 0)            send => dat(:,:,1:n_send)
        if (size(dat,3) > r%n_rec) recv => dat(:,:,r%n_rec+1:)
      else
        n_send =  0
      end if
      call p_gatherv(send, recv, pe_gather)
    end subroutine gather_real_3d

    subroutine gather_log_2d(dat)
      logical, pointer :: dat(:,:)
      logical, target  :: dum(0,0)
      logical, pointer :: send(:,:), recv(:,:)
      n_send = r%n_rec
      if (dace% pe == pe_gather) n_send = 0
      send => dum
      recv => dum
      if (associated(dat)) then
        if (n_send > 0)            send => dat(:,1:n_send)
        if (size(dat,2) > r%n_rec) recv => dat(:,r%n_rec+1:)
      else
        n_send =  0
      end if
      call p_gatherv(send, recv, pe_gather)
    end subroutine gather_log_2d

  end subroutine gather_radv
!------------------------------------------------------------------------------
  subroutine prep_rttov_prof(obs, spot, xi, rtp, ierr, ttovs, tjac)
    type(t_obs_block),   intent(inout)                  :: obs  ! observation datatype
    type(t_spot),        intent(inout)                  :: spot ! SPOT observations
    type(t_vector_segm), intent(in)                     :: xi   ! interpolated values
    type(t_rttov_prof),  intent(out)                    :: rtp  ! arguments to rttov
    integer,             intent(out),  optional         :: ierr ! error return code
    type(t_tovs),        intent(out),  optional, target :: ttovs
    type(t_jac),         intent(out),  optional, target :: tjac

    type(t_jac),  target     :: tj_aux
    type(t_jac),  pointer    :: tj         => null()
    type(t_tovs), target     :: tovs_aux
    type(t_tovs), pointer    :: tovs       => null()
    logical                  :: l_debug    = .false.
    logical                  :: l_fill_jac = .false.
    integer                  :: ncv                       ! number of control variables
    integer                  :: npv                       ! number of profile variables
    integer                  :: npcl                      ! number of emissivity PCs
    integer                  :: nld_t                     ! max. T dummy level
    integer                  :: nld_h                     ! max. hum. dummy level
    integer                  :: k, k_lev, i1, i2, istat
    integer                  :: i_t, i_q
    integer                  :: tovs_io
    integer                  :: emis_flags
    integer,     allocatable :: i_fail(:)                 ! error return flag (tq_tvgh)
    real(wp)                 :: tmp(2), q_a, t_a, tv_a, p_a, rh_a, rh, drh_gh
    real(wp)                 :: tmin_, qmin_, tmax_, qmax_
    real(wp),    allocatable :: p(:)
    logical,     pointer     :: ps_llev         => NULL() ! ps is below lowest level -> true
    integer,     pointer     :: il              => NULL() ! indices for interpol. t2m, q2m
    integer,     pointer     :: iu              => NULL() ! indices for interpol. t2m, q2m
    real(wp),    pointer     :: wl              => NULL() ! weight for interpol. t2m, q2m
    real(wp),    pointer     :: wu              => NULL() ! weight for interpol. t2m, q2m
    real(wp),    pointer     :: d_lev           => NULL()
    real(wp),    pointer     :: d_tv            => NULL()
    real(wp),    pointer     :: d_rh            => NULL()
    real(wp),    pointer     :: dt_tv (:)       => NULL() ! d T    / d Tv
    real(wp),    pointer     :: dt_rh (:)       => NULL() ! d T    / d rh
    real(wp),    pointer     :: dq_tv (:)       => NULL() ! d q    / d Tv
    real(wp),    pointer     :: dq_rh (:)       => NULL() ! d q    / d rh
    real(wp),    pointer     :: dt_tv_2m        => NULL() ! d T_2m / d Tv_2m
    real(wp),    pointer     :: dt_rh_2m        => NULL() ! d T_2m / d rh_2m
    real(wp),    pointer     :: dq_tv_2m        => NULL() ! d q_2m / d Tv_2m
    real(wp),    pointer     :: dq_rh_2m        => NULL() ! d q_2m / d rh_2m
    real(wp),    pointer     :: de_dpc(:,:)     => NULL() ! d emis / d pc
    real(wp),    pointer     :: sigma_var_tskin => NULL() ! surface temperature error
    real(wp),    pointer     :: dclt_dum        => NULL() ! d clt / d dummy
    real(wp),    pointer     :: dclf_dum        => NULL() ! d clf / d dummy
    real(wp),    pointer     :: dsnf_dum        => NULL() ! d snf  / d dummy

    FTRACE_BEGIN('prep_rttov_prof')

    !-------------
    ! Preparations
    !-------------
    FTRACE_BEGIN('prep_rttov_prof:init')
    FTRACE_BEGIN('prep_rttov_prof:init1')
    l_debug = ldeb(spot)

    if (present(ierr)) ierr = 0

    ! rttov init must have been called before
    if (.not. rttvi_called) then
      call finish ('prep_rttov_prof','RTTOV must be initialized')
    end if

    ! TODO: here is optimization potential: we could use preallocated arrays to get
    !       results from t_tovs
    ! load specific parameters from t_tovs
    tovs_io = TTOVS_AV+TTOVS_CEMI+TTOVS_EMIS+TTOVS_CI+TTOVS_FLAG+TTOVS_SPEC+TTOVS_BASE
    if (present(ttovs)) then
      tovs => ttovs
      if (iand(tovs%init, tovs_io) /= tovs_io) then
        call destruct(tovs)
        call load (obs% o, spot, tovs)
      end if
    else
      tovs => tovs_aux
      call load (obs% o, spot, tovs, tovs_io=tovs_io)
    end if
    FTRACE_END('prep_rttov_prof:init1')

    FTRACE_BEGIN('prep_rttov_prof:init2')
    npv   = 2*tovs%nlev
    ncv   = spot% i% n
    npcl  = ncv - nsv - npv
    nld_t = tovs%nld_t
    nld_h = tovs%nld_h
    i_t   = tovs%i_t
    i_q   = tovs%i_q

    allocate(p(tovs%nlev))
    if (tovs%i_p > 0) then
      p(:) = tovs%av(:,tovs%i_p) !* 0.01_wp
    else
      p(:) = preslev(:)
    end if

    if (rt_hard_limits == 1) then
      qmin_ = 1.01_wp * qmin
      tmin_ = 1.01_wp * tmin
      qmax_ = 0.99_wp * qmax
      tmax_ = 0.99_wp * tmax
    else if (rt_hard_limits >= 2) then
      qmin_ = qmin
      tmin_ = tmin
      qmax_ = qmax
      tmax_ = tmax
    end if
    FTRACE_END('prep_rttov_prof:init2')

    FTRACE_END('prep_rttov_prof:init')

    FTRACE_BEGIN('prep_rttov_prof:prep_jac')
    ! prepare t_jac
    l_fill_jac = present(tjac)
    if (l_fill_jac) then
      tj => tjac
    else
      tj => tj_aux
    end if
    call construct(tj)
    ps_llev         => tj% ps_llev
    iu              => tj% iu
    il              => tj% il
    wu              => tj% wu
    wl              => tj% wl
    d_lev           => tj% d_lev
    sigma_var_tskin => tj% sigma_var_tskin
    if (l_fill_jac) then
      allocate ( tj% dt_tv  (tovs%nlev))               ; dt_tv => tj% dt_tv
      allocate ( tj% dt_rh  (tovs%nlev))               ; dt_rh => tj% dt_rh
      allocate ( tj% dq_tv  (tovs%nlev))               ; dq_tv => tj% dq_tv
      allocate ( tj% dq_rh  (tovs%nlev))               ; dq_rh => tj% dq_rh
      allocate ( tj% de_dpc (0:npcl, tovs% nchan)) ; de_dpc=> tj% de_dpc
      d_tv          => tj% d_tv
      d_rh          => tj% d_rh
      dt_tv_2m      => tj% dt_tv_2m
      dt_rh_2m      => tj% dt_rh_2m
      dq_tv_2m      => tj% dq_tv_2m
      dq_rh_2m      => tj% dq_rh_2m
      dclt_dum      => tj% dclt_dum
      dclf_dum      => tj% dclf_dum
      dsnf_dum      => tj% dsnf_dum
    end if
    FTRACE_END('prep_rttov_prof:prep_jac')

    FTRACE_BEGIN('prep_rttov_prof:fill1')
    ! Construct t_rttov_prof, which holds most results of this routine
    call construct (rtp, tovs% nchan, jakobi=0, n_av=2, nlev=tovs%nlev)

    !----------------------------------------------------------------
    ! fill atmospheric profile (tv, gh) and for dummy variables (T,q)
    !----------------------------------------------------------------
    rtp% av (:,1,1) = xi% x (spot%i%i+1:spot%i%i+npv:2) ! tv
    rtp% av (:,2,1) = xi% x (spot%i%i+2:spot%i%i+npv:2) ! gh
    if (warn_gen_hum) then
      do k = 1, tovs%nlev
        if (rtp% av(k,2,1) < gh_min) then
          write(6,'("WARNING: gen. hum. below lower bound in level ",I2,&
               &f6.3," spot ",I10," sat ",I3.3," grid ",I3.3,2x,f7.2,"N",f7.2,"E")') &
               k,rtp% av(k,2,1),spot%hd%id,spot%ident,rttov_instr(spot%sttyp,int(spot%hd%satid)),&
               spot%col%c%dlat,spot%col%c%dlon
        elseif (rtp% av(k,2,1) > gh_max) then
          write(6,'("WARNING: gen. hum. above upper bound in level ",I2,&
               &f6.3," spot ",I10," sat ",I3.3," grid ",I3.3,2x,f7.2,"N",f7.2,"E")') &
               k,rtp% av(k,2,1),spot%hd%id,spot%ident,rttov_instr(spot%sttyp,int(spot%hd%satid)),&
               spot%col%c%dlat,spot%col%c%dlon
        end if
      end do
    end if
    FTRACE_END('prep_rttov_prof:fill1')

    FTRACE_BEGIN('prep_rttov_prof:fill2')
    ! Do variable transformation for stratospheric dummy variables
    if (l_debug) then
      write(usd,*) dpref,'prep po',spot%hd%id,po_context,po_ilns,po_lmon
      write(usd,*) dpref,'prep nld_t,nld_h',nld_t,nld_h
      do k = 1, nld_h
        if (i_q > 0) then
          write(usd,*) dpref,'prep gh,q (0q)',k,rtp% av (k, 2,1),tovs% av (k,i_q)
        else
          write(usd,*) dpref,'prep gh,q (0q)',k,rtp% av (k, 2,1)
        end if
      end do
      do k = 1, nld_t
        if (i_t > 0) then
          write(usd,*) dpref,'prep gh,q (0t)',k,rtp% av (k, 1,1),tovs% av (k,i_t)
        else
          write(usd,*) dpref,'prep gh,q (0t)',k,rtp% av (k, 1,1)
        end if
      end do
    end if
    if (i_q > 0) then
      tovs% av(1:nld_h,i_q) = max(tovs% av(1:nld_h,i_q), real(qmin_,kind=tpp))
      tovs% av(1:nld_h,i_q) = min(tovs% av(1:nld_h,i_q), real(qmax_,kind=tpp))
      if (l_debug) then
        do k = 1, nld_h
          write(usd,*) dpref,'prep gh,q (0q)',k,rtp% av (k, 2,1),tovs% av (k,i_q)
        end do
      end if
    end if
    if (i_t > 0) then
      tovs% av(1:nld_t,i_t) = max(tovs% av(1:nld_t,i_t), real(tmin_,kind=tpp))
      tovs% av(1:nld_t,i_t) = min(tovs% av(1:nld_t,i_t), real(tmax_,kind=tpp))
    end if
    ! T dummy variable
    rtp% av (1:nld_t,1,1) = e_bg_t_top * rtp% av (1:nld_t,1,1)
    if (i_t > 0) then ! background profile available
      rtp% av (1:nld_t,1,1) = rtp% av (1:nld_t,1,1) + &
                                      tovs% av (1:nld_t,i_t)
    else
      rtp% av (1:nld_t,1,1) = rtp% av (1:nld_t,1,1) + &
                                      tmin_
    end if
    ! Q dummy variable
    rtp% av (1:nld_h,2,1) = e_bg_lnq_top * rtp% av (1:nld_h,2,1)
    if (i_q > 0) then ! background profile available
      rtp% av (1:nld_h,2,1) = rtp% av (1:nld_h,2,1) + &
                                       log(tovs% av (1:nld_h,i_q))
    else
      rtp% av (1:nld_h,2,1) = rtp% av (1:nld_h,2,1) + &
                                       log(qmin_)
    end if
    rtp% av (1:nld_h,2,1) = exp(rtp% av (1:nld_h,2,1))

    FTRACE_END('prep_rttov_prof:fill2')
    if (l_debug) then
      write(usd,*) dpref,'prep  level',nld_t,nld_h,k_lev
      do k = 1,nld_h
        write(usd,*) dpref,'prep  gh,q (1)',k,rtp% av (k,1:2,1)
      end do
    end if

    FTRACE_BEGIN('prep_rttov_prof:fill3')
    !---------------------------------
    ! Surface variables and parameters
    !---------------------------------
    if (no_ps_dep) then
      rtp% sav (3,1) = exp (obs% o% lev (spot% i% i+npv+2))
    else
      ! interpolate surface pressure
      !   -> with:
      !       obs% o% lev(spot% i% i+npv+2) ^= ln (ps_bg)
      !       exp(obs% o% lev(spot% i% i+npv+2)) ^= ps_bg, [obs% o% lev ^= ln(p)]
      !       spot% pz_bg ^= dp / dz
      !       xi% x(spot% i% i+npv+2) ^= geop(in ps at current iteration)
      !    => total ^= ps_bg + dp/dz * ( geop(iteration_i) - geop_bg ) = ps (iteration i)
      rtp% sav (3,1) = exp (obs% o% lev (spot% i% i+npv+2 )) +   &     !  ps (Pa)
                       spot% pz_bg * (xi% x(spot% i% i+npv+2) - spot% gp_bg/gacc)
      if (l_debug) write(usd,*) dpref,'prep ps',rtp% sav (3,1),exp (obs% o% lev (spot% i% i+npv+2 )),&
           obs% o% lev (spot% i% i+npv+2 ), spot% pz_bg * (xi% x(spot% i% i+npv+2) - spot% gp_bg/gacc), &
           spot% pz_bg, xi% x(spot% i% i+npv+2),spot% gp_bg/gacc
    endif
    ! prepare interpolation between levels to surface
    ps_llev  =  .true.
    do k = 1, size(p)
       if (p(k) .gt. rtp% sav(3,1)) then
          il = k
          iu = k-1
          ps_llev = .false.
          exit
       endif
    enddo
    if (.not. ps_llev) then
       d_lev = p(iu) - p(il)
       wl   = (p(iu) - rtp% sav  (3,1)) / d_lev
       wu   = (rtp% sav  (3,1) - p(il)) / d_lev
    else
       wl   = 1.0_wp
       wu   = 0.0_wp
       il       = size(p)
       iu       = size(p)
    endif
    rtp% sav  (1,1) = wl * rtp% av(il,1,1) + wu * rtp% av(iu,1,1)  ! tv2m
    rtp% sav  (2,1) = wl * rtp% av(il,2,1) + wu * rtp% av(iu,2,1)  ! rh2m
    if (l_fill_jac) then
      d_tv  =   rtp% av(iu,1,1) - rtp% av(il,1,1)
      d_rh  =   rtp% av(iu,2,1) - rtp% av(il,2,1)
    end if
    rtp% sav  (4,1) = spot% u10bg                          ! u10m
    rtp% sav  (5,1) = spot% v10bg                          ! v10m

    !----------------
    ! Skin parameters
    !----------------
!    if (no_ts_dep) rtp% ssv (2:,1) = 0._wp
    ! Variable transformation for skin temp. dummy variable
    ! initially Ts=ts_bg, at later iterations
    ! corrected by gradient (scaled with assumed variance sigma_var_tskin)
    FTRACE_END('prep_rttov_prof:fill3')

    FTRACE_BEGIN('prep_rttov_prof:bc_xd')
    !-----------------------------------------
    ! cloud top and fraction and snow fraction
    !-----------------------------------------
    ! TODO: what happens here if we did not initialize the sink variables
    if (l_fill_jac) then
      call bc_xd (rtp% cv(1,1), xi% x(spot%i%i+npv+3), obs%o% sink(spot%d%i+1), dx=dclt_dum)
      call bc_xd (rtp% cv(2,1), xi% x(spot%i%i+npv+4), obs%o% sink(spot%d%i+2), dx=dclf_dum)
      call bc_xd (rtp% snf    , xi% x(spot%i%i+npv+5), obs%o% sink(spot%d%i+3), dx=dsnf_dum)
    else
      call bc_xd (rtp% cv(1,1), xi% x(spot%i%i+npv+3), obs%o% sink(spot%d%i+1))
      call bc_xd (rtp% cv(2,1), xi% x(spot%i%i+npv+4), obs%o% sink(spot%d%i+2))
      call bc_xd (rtp% snf    , xi% x(spot%i%i+npv+5), obs%o% sink(spot%d%i+3))
    end if
    rtp% cv(1,1) = tovs% cloud_top
    rtp% cv(2,1) = tovs% cloud_imag
    if (l_debug) write(usd,*) dpref,'prep simpl_cloud',rtp% cv(1:2,1)

    FTRACE_END('prep_rttov_prof:bc_xd')
    ! TODO, ATTENTION: this is contradicting to the above use of snf as a sink variable.
    ! However, snf is only used in mo_ir_emis as a sink variable. Should be fixed whenever
    ! mo_ir_emis is used again.
    rtp%snf = min(max(real(tovs%snf,wp),0._wp),1._wp)
    select case(snw_frc_mode)
    case(1)
      if (all(snw_frc_par(1:4) > 0._wp)) then
        ! Taken from https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2012JD018178
        ! Niu and Yang 2007, J. Geophys. Res., 112, D21101
        ! (simplified: rho_snow = rho_new)
        rtp%snf = tanh(spot%hs_bg/ &
             (2.5_wp*snw_frc_par(1)*(snw_frc_par(2)/snw_frc_par(3))**snw_frc_par(4)))
        if (l_debug) write(usd,*) dpref,' prep set snf',spot%hs_bg,rtp%snf
        tovs%snf = rtp%snf
      else
        call finish('prep_rttov_prof', 'snw_frc_mode==1 requires snw_frc_par(1:2) > 0')
      end if
    end select

    call set_surface_type(spot, tovs)
    call store(obs% o, spot, tovs, tovs_io=TTOVS_BASE)
    sigma_var_tskin = tovs%e_bg_ts

    if (tovs% ts > 0) then ! tovs% ts stores retrieved tskin
      rtp% ssv(1,1) = tovs% ts
    else
      rtp% ssv(1,1) = spot% ts_bg
    end if
    if (l_debug) write(usd,*) dpref,'prep tskin',rtp% ssv(1,1),spot% ts_bg,tovs%ts
    rtp% ssv(1,1) = rtp% ssv(1,1) + sigma_var_tskin * xi% x(spot%i%i+npv+1) ! tsurf
    if (l_debug) write(usd,*) dpref,'prep tskin',rtp% ssv(1,1)

    !------------
    ! Coordinates
    !------------
    rtp% sv(1,1) = spot% col% c% dlat
    rtp% sv(2,1) = spot% col% c% dlon
    rtp% hsurf   = spot% gp_bg/(gacc*1.e3_wp)

    FTRACE_BEGIN('prep_rttov_prof:trans1')
    !-------------------------------------------------------------------------
    ! Transform variables to physical (RTTOV input) variables (Tv,gh) -> (T,q)
    !-------------------------------------------------------------------------
    ! Dummy variables
    if (l_debug) write(usd,*) dpref,'prep level',nld_t,nld_h
    allocate (i_fail(tovs%nlev))
    do k = 1, tovs%nlev
      if ((nld_t >= k) .and. (nld_h >= k)) then
        ! both t and q are dummy variables: nothing to convert, we have (T,q) already
        if (l_debug) write(usd,*) dpref,'prep gh,q t_q_dum',k
!        if (k <= no_t_dep) rtp% av (k,  1,1) = tovs% av (k,i_t)
        if (l_fill_jac) then
          dt_tv(k) = e_bg_t_top
          dt_rh(k) = 0.0_wp
          dq_tv(k) = 0.0_wp
          dq_rh(k) = rtp% av (k,2,1) * e_bg_lnq_top
        end if
      elseif ((nld_t >= k) .and.  &
              (nld_h < k)) then
        ! only t is dummy variable:
        ! call finish('prep_rttov_prof','T dummy without Q dummy not implemented')
        call rh_gh(rh, rtp% av (k,2,1))
        tmp(1) = q_rh(rh, rtp% av (k,1,1), p(k))
        if (l_fill_jac) then
          call rh_gh(rh, rtp% av (k,2,1), drh_gh)
          q_a  = 1._wp
          rh_a = 0._wp
          t_a  = 0._wp
          p_a  = 0._wp
          call q_rh_adj(q_a,rh_a,t_a,p_a,rh, rtp% av (k,1,1), p(k))
          dt_tv(k) = e_bg_t_top
          dt_rh(k) = 0.0_wp
          dq_tv(k) = t_a * e_bg_t_top
          dq_rh(k) = rh_a * drh_gh
        end if
        rtp% av (k,2,1) = tmp(1)
      elseif ((nld_t < k) .and. &
              (nld_h >= k)) then
        ! only q is (IFS) dummy variable:
        tmp(1) = t_tv_q (rtp% av (k,1,1), rtp% av (k,2,1))
        if (l_debug) write(usd,*) dpref,'prep gh,q q_dum',k,rtp% av (k,1:2,1),tmp(1),0.
        if (l_fill_jac) then
          t_a    = 1._wp
          tv_a   = 0._wp
          q_a    = 0._wp
          call t_tv_q_adj (t_a, tv_a, q_a, tv=rtp% av (k,1,1), q=rtp% av (k,2,1))
          dt_tv(k) = tv_a                                  ! part. deriv. dt/dtv
          dq_tv(k) = 0._wp                                 ! part. deriv. dq/dtv
          dt_rh(k) = q_a * rtp% av (k,2,1) * e_bg_lnq_top  ! dt/drh
          dq_rh(k) =       rtp% av (k,2,1) * e_bg_lnq_top  ! dq/drh
        end if
        rtp% av (k,1,1) = tmp(1)
      endif
    enddo
    FTRACE_END('prep_rttov_prof:trans1')

    FTRACE_BEGIN('prep_rttov_prof:trans2')
    ! Non-dummy atmospheric profile variables
    k_lev = max(nld_h,nld_t,0) + 1
    i1    = spot%i%i+2*(k_lev-1)
    i2    = spot%i%i+npv
    if (l_fill_jac) then
#ifdef __NEC__
      call tq_tvgh_adj_vec (tovs%nlev-k_lev+1,                            &
#else
      call tq_tvgh         (                                          &
#endif
                            rtp% av(k_lev:,1,1), rtp% av(k_lev:,2,1), & ! T , q
                            xi% x(i1+1:i2:2),    xi% x(i1+2:i2:2),    & ! Tv, gh
                            p(k_lev:), i_fail(k_lev:),          &
                            dt_tv(k_lev:), dt_rh(k_lev:), dq_tv(k_lev:), dq_rh(k_lev:))
    else
#ifdef __NEC__
      call tq_tvgh_vec (tovs%nlev-k_lev+1,                            &
#else
      call tq_tvgh     (                                          &
#endif
                        rtp% av(k_lev:,1,1), rtp% av(k_lev:,2,1), & ! T , q
                        xi% x(i1+1:i2:2),    xi% x(i1+2:i2:2),    & ! Tv, gh
                        p(k_lev:), i_fail(k_lev:))
    end if
    if(any(i_fail(k_lev:)<0)) then
      if (present (ierr)) then
        call message ('prep_rttov_prof','tq_tvgh (profile) did not converge')
        k = sum(minloc(i_fail(k_lev:)))+k_lev-1
        write(0,*)'t,q,tv,gh,p,ifail,n=',rtp% av (k,1,1),rtp% av (k,2,1),&
             xi% x (spot%i%i+2*k-1),xi% x (spot%i%i+2*k),p(k),i_fail(k),  &
             count(i_fail(k_lev:)<0)
        ierr = -1
        goto 998
      else
        write(0,*) 'spot',spot%hd%id
        write(0,*) 'i_fail',i_fail(k_lev:)
        write(0,*) 'tv',xi% x(i1+1:i2:2)
        write(0,*) 'gh',xi% x(i1+2:i2:2)
        call finish('prep_rttov_prof','tq_tvgh (profile) did not converge')
      endif
    endif
    FTRACE_END('prep_rttov_prof:trans2')

    if (l_debug) then
      do k = k_lev, tovs%nlev
        write(usd,*) dpref,'prep gh,q (2)',k,p(k),xi% x (spot%i%i+k*2-1:spot%i%i+k*2),rtp% av (k,1:2,1)
      end do
    end if


    FTRACE_BEGIN('prep_rttov_prof:trans4')
    ! Surface variables
    tmp(1:2) = rtp% sav (1:2,1)
    if (l_fill_jac) then
      call tq_tvgh(rtp% sav(1,1), rtp% sav(2,1), & ! Tv,gh
                   tmp(1),        tmp(2),        & ! T ,q
                   spot% ps_bg, i_fail(1),       &
                   dt_tv_2m, dt_rh_2m, dq_tv_2m, dq_rh_2m)
    else
      call tq_tvgh(rtp% sav(1,1), rtp% sav(2,1), & ! Tv,gh
                   tmp(1),        tmp(2),        & ! T ,q
                   spot% ps_bg, i_fail(1))
    end if
    if(i_fail(1)<0) then
      if (present (ierr)) then
        call message ('prep_rttov_prof','tq_tvgh (surface 2m) did not converge')
        write(0,*)'t,q,tv,gh,p,ifail=',rtp% sav (1:2,1),tmp(1:2),rtp% sav (3,1), i_fail(1)
        ierr = -1
        goto 998
      else
        call finish('prep_rttov_prof','tq_tvgh (surface 2m) did not converge')
      endif
    endif
    rtp% sav  (3,1) = rtp% sav  (3,1) / 100._wp    ! ps  (hPa)
    FTRACE_END('prep_rttov_prof:trans4')

    !-----------
    ! Emissivity
    !-----------
    FTRACE_BEGIN('prep_rttov_prof:emis')
    if (npcl > 0) then
      i2 = spot% i% i + ncv
      i1 = i2 - npcl + 1
      if (l_fill_jac) then
        call emis_pc (spot, tovs% cemi, obs% o, rtp% snf, ttovs%rt_stype(1), &
                      xi% x(i1:i2), rtp% emis(1,:), de_dpc)
      else
        call emis_pc (spot, tovs% cemi, obs% o, rtp% snf, ttovs%rt_stype(1), &
                      xi% x(i1:i2), rtp% emis(1,:))
      end if
    elseif (any(iand(tovs% flag(:),EMIS_CALC) == 0) .and. obs% o% pe == dace% pe) then
      call set_emis(spot, tovs, rtp, obs% o, istat)
      if ( istat /= 0 ) call finish('process_tovs','set_emis failed')
      if (l_fill_jac) de_dpc = 0._wp
    else
      ! Usually the emissivity from tovs%emis should be used, but RTTOV with FASTEM
      ! does not give the same results in this case, i.e. a call with emiss<=0.
      ! (call to FASTEM) does not give the same result as a call with the
      ! emissivity calculated by FASTEM.  (This is due to the diffuse reflection for
      ! polarized channels.) Thus, we have to force a call to FASTEM again and again
      ! for each RTTOV call.
      ! For mec we want to use the emissivity calculated within the analyses and do
      ! not want to use a (usually much worse) emissivity calculated from the forecasts.
      ! Thus, for mec we do not force FASTEM calls and take into account the small inaccuracy
      ! of the result.
      do k = 1, tovs%nchan
        emis_flags = iand(tovs%flag(k), EMIS_SUCCESS) ! Get all the emissivity calc. method flags
        if (l_debug) write(usd,*) dpref,'prep emis tf',k,rtp%emis(1,k),&
             tovs%flag(k),emis_flags,2**TF_EMIS_SEA_MOD,tovs% emis(k)
        if (emis_flags == 2**TF_EMIS_SEA_MOD) then
          ! Pure FASTEM - no coastal mixture
          if (force_fastem) then
            ! var3d
            rtp% emis(1,k) = -1.
          else
            ! mec: do not call FASTEM, use emissivity from analysis and take into account
            !      the little error due to the incorrect diffuse reflection
            rtp% emis(1,k) = abs(tovs%emis(k))
          end if
        else
          rtp% emis(1,k) = tovs% emis(k)
        end if
      end do
      if (l_fill_jac) de_dpc = 0._wp
    endif
    if (l_debug) then
      do k = 1, tovs%nchan
        write(usd,*) dpref,'prep emis',k,rtp%emis(1,k)
      end do
    end if
    FTRACE_END('prep_rttov_prof:emis')

    !-----------
    ! Tskin
    !-----------
    FTRACE_BEGIN('prep_rttov_prof:tskin')
    if (l_debug) write(usd,*) dpref,'prep ts',spot% ts_bg,tovs%ts,any(iand(tovs% flag(:), TSKIN_CALC) == 0),&
         tovs% flag(:), TSKIN_CALC
    if (any(iand(tovs% flag(:), TSKIN_CALC) == 0) .and. tovs%ts < 0._wp .and. obs% o% pe == dace% pe) then
      call set_tskin(spot, tovs, rtp, obs% o, istat)
      if (l_debug) write(usd,*) dpref, ' set_tskin done ',rtp%ssv(1,1)
      if ( istat /= 0 ) call finish('process_tovs','set_tskin failed')
    end if
    if (l_debug) write(usd,*) dpref,'prep ts',rtp%ssv(1,1),tovs%ts,spot%ts_bg
    FTRACE_END('prep_rttov_prof:tskin')
    FTRACE_END('prep_rttov_prof')

    ! Cleanup
998 if (.not.present(ttovs)) call destruct(tovs)

  end subroutine prep_rttov_prof


  subroutine set_surface_type(spot, tovs, l_check)
    type(t_spot), intent(inout)           :: spot
    type(t_tovs), intent(inout)           :: tovs
    logical,      intent(inout), optional :: l_check ! input: if .true. then compare new results
                                                     !        with old results in tovs. Do not
                                                     !        update tovs
                                                     ! output:.true. if new result are equal,
                                                     !        .false.
    character(len=16), parameter :: proc = 'set_surface_type'
    real(wp),          parameter :: eps_lf = 0.01_wp
    real(sp),          parameter :: eps    = 3 * epsilon(1._sp)
    integer,           parameter :: SUR_ALL = 2**SUR_LAND  + 2**SUR_SEA  + 2**SUR_ICE + &
                                              2**SUR_NOICE + 2**SUR_SNOW + 2**SUR_NOSNOW
    integer                      :: rts(n_styp)      ! resulting surface types
    real(sp)                     :: ws (n_styp)      ! resulting weights
    integer                      :: ns               ! resulting number of surface types
    real(sp)                     :: e_bg_ts          ! resulting e_bg_ts
    integer                      :: styp             ! auxiliary surface type
    integer                      :: stlsf            ! auxiliary stlsf
    integer                      :: i,j
    integer                      :: ii(1)
    real(wp)                     :: fr, fr_ice, snf_
    real(sp)                     :: land_frac, w
    logical                      :: mask(n_styp)
    logical                      :: lchk
    real(wp)                     :: e_bg_ts_land_sel

    ! TODO: inconsistent usage of tovs%land_frac vs. spot%sl_bg

    ! if (tovs%rt_stype >= 0) return  ! Do only for the first call (tovs%rt_stype == -1)
    ! This was uncommented in order to allow a new calculation of surface type in the global
    ! LETKF (Was important for backwards compatibility)

    if (ldeb(spot)) write(usd,*) dpref,'set_surface_type',surf_class_vers,spot%col%c%dlat,&
         spot%col%c%dlon,spot%stlsf,spot%sl_bg,spot%fi_bg,spot%ts_bg,tovs%mw_stype,tovs%snf

    if (present(l_check)) then
      lchk = l_check
    else
      lchk = .false.
    end if

    ns = 1
    ws(1) = 1._wp

    if (tovs% l_ts_dr) then
      e_bg_ts_land_sel = e_bg_drts_land
    else
      e_bg_ts_land_sel = e_bg_ts_land
    end if

    select case(surf_class_vers)
    case(0)
      if (btest(spot%stlsf, SUR_SEA) .and. .not.btest(spot%stlsf, SUR_LAND) .and. spot% sl_bg < 0.5 ) then
        ! sea
        e_bg_ts  = e_bg_ts_sea
        rts(1) = rts_sea
      else
        ! coeffs.land: open grass
        e_bg_ts  = e_bg_ts_land_sel
        rts(1) = rts_land
      endif
    case(1,2)
      ! 1. Decide whether land/sea/seaice/coast (1/2/3/4)
      styp=1  ! sea
      if (.not.btest(spot%stlsf,SUR_SEA)) then
        styp=2  ! land
      else
        if (btest(spot%stlsf,SUR_LAND) .or. spot%sl_bg > eps_lf) then
          styp=4  ! coast
        end if
        if ( (tovs%mw_stype >= 4 .and. tovs%mw_stype <= 9)                        .or. &
             ((spot%fi_bg > 0.01 .or. spot%ts_bg < 271.5)  .and. spot%sl_bg < 0.5) ) then
          styp=3   ! sea ice
        end if
      end if
      if (surf_class_vers == 2) then
        fr_ice = spot%fi_bg
        if (fr_ice < 0.1_wp .and. spot%ts_bg < 271.5) fr_ice = 1._wp
        fr = tovs%snf * tovs%land_frac + fr_ice * (1._wp - tovs%land_frac)
        if (fr > 0.8_wp) then !! questionable threshold !!
          styp = 3
        end if
      end if
      select case(styp)
      case(1) ! sea
        rts(1) = rts_sea
        e_bg_ts  = e_bg_ts_sea
      case(2) ! land
        rts(1) = rts_land
        e_bg_ts  = e_bg_ts_land_sel
      case(3) ! seaice/snow
        rts(1) = rts_ice
        e_bg_ts  = e_bg_ts_ice
      case(4) ! coast
        rts(1) = rts_land
        e_bg_ts  = sqrt(e_bg_ts_sea **2 * (1._wp - tovs%land_frac) + &
                             e_bg_ts_land_sel**2 * (tovs%land_frac))
      case default
        call finish('set_surface_type', 'invalid styp')
      end select
    case(3)
      if (tovs%land_frac <= 0.5) then
        rts(1) = rts_sea
        e_bg_ts  = e_bg_ts_sea
      else
        rts(1) = rts_land
        e_bg_ts  = e_bg_ts_land_sel
        if (spot%fi_bg > 0.5) then
          rts(1) = rts_ice
          e_bg_ts  = e_bg_ts_ice
        end if
      end if
    case(4)
      ! Allow for multiple surface types
      ns = 0
      land_frac = real(tovs%land_frac,wp)
      ! Take into account the model surface type, which might differ if the
      ! sea-land-mask(s) in icon and sat_pp differ.
      ! TODO: do something better
      if (spot%col%c%dlat <= -75._wp) then
        land_frac = max(land_frac, real(spot%sl_bg,sp))
      end if
      if (land_frac <= eps_lf) then
        if (spot%sl_bg >= eps_lf) then
          land_frac = 0.5_wp * (land_frac + spot%sl_bg)
        end if
      end if
      if (land_frac >= 1._wp-eps_lf) then
        if (spot%sl_bg <= 1._wp-eps_lf) then
          land_frac = 0.5_wp * (land_frac + spot%sl_bg)
        end if
      end if

      if (ldeb(spot)) write(usd,*) dpref, proc,' start',land_frac,spot%sl_bg,spot%fi_bg,spot%ts_bg,tovs%snf

      ! Basic classification
      if (land_frac >= 1._wp-eps_lf) then
        ! land
        ns = 1
        rts(ns) = rts_land
        ws (ns) = 1._wp
      elseif (land_frac <= eps_lf) then
        ! sea
        ns = 1
        rts(ns) = rts_sea
        ws (ns) = 1._wp
      else
        ! coast
        ns = 2
        rts(1:2) = (/rts_land,rts_sea/)
        ws (1:2) = (/land_frac,1._sp-land_frac/)
      end if
      do i = 1, ns
        if (ldeb(spot)) write(usd,*) dpref, proc,i,rts(i),ws(i)
      end do

      ! Sea ice fraction
      fr = 0._wp
      if (spot%sl_bg < 1._wp-eps_lf) then
        fr = spot%fi_bg
        !if (spot%ts_bg < 271.5) fr = 1._wp
      end if
      if (tovs%mw_stype >= 4 .and. tovs%mw_stype <= 9 .and. spot%ts_bg <= 275._wp) &
           fr = max(fr,0.2_wp)
           ! tovs%mw_stype might be ice also over clouds, therefore the safeguard spot%ts_bg < 275.
      fr = fr * (1._wp - land_frac)
      if (fr > eps_lf) then
        ! Add sea ice
        ns = ns+1
        if (ns > n_styp) call finish(proc,'ns too large (1)')
        rts(ns) = rts_ice
        ws(ns)   = fr
#if defined(__GFORTRAN__) && (__GNUC__ <= 8)    /* FINDLOC is Fortran 2008.  */
        do i = 1, ns
           if (rts(i) == rts_sea) exit
        end do
        if (i > ns) i = 0
#else
        ii = findloc(rts(1:ns), rts_sea) ; i = ii(1)
#endif
        if (i > 0) then
          ws(i) = ws(i) - fr
        else
          call finish(proc,'failed to find sea styp')
        end if
      end if

      ! ! Snow fraction
      ! snf_ = tovs%snf
      ! if (snf_ <= 0._wp .and. spot%fi_bg > 0.5_wp .and. &
      !     spot%sl_bg <= 0._wp .and. land_frac >= eps_lf) then
      !   ! We have no model land fraction and therefore no snow cover, but the spot
      !   ! (which is usually much bigger than the model resolution) sees land. If sea
      !   ! ice is dominating we can assume that the land is also covered by snow.
      !   ! This reduces the emissivity atlas value misses in polar regions.
      !   snf_ = 1._wp
      ! end if

      ! ! Replace snow
      ! fr = snf_ * land_frac
      ! if (fr >= eps_lf) then
      !   ii = findloc(rts(1:ns), rts_ice) ; i = ii(1)
      !   if (i <= 0) then
      !     ! Add ice
      !     ns = ns+1
      !     if (ns > n_styp) call finish(proc,'ns too large (2)')
      !     rts(ns) = rts_ice
      !     ws  (ns) = fr
      !   else
      !     ws  (i)  = ws  (i) + fr
      !   end if
      !   ii = findloc(rts(1:ns), rts_land) ; i = ii(1)
      !   if (i > 0) then
      !     ws  (i)  = ws  (i) - fr
      !   else
      !     call finish(proc,'failed to find land styp')
      !   end if
      ! end if

      ! Remove styp with too small weights
      mask(1:ns) = (ws(1:ns) >= eps_lf)
      i = count(mask(1:ns))
      if (i == 0) then
        call finish(proc,'number of surface types == 0')
      elseif (i < ns) then
        rts(1:i) = pack(rts(1:ns),mask=mask(1:ns))
        ws  (1:i) = pack(ws  (1:ns),mask=mask(1:ns))
        ns = i
      end if
      do i = 1, ns
        if (ldeb(spot)) write(usd,*) dpref, proc,i,rts(i),ws(i)
      end do

      ! Sort (decreasing ws)
      do i = ns-1,1,-1
        do j = 1, i
          if (ws(j+1) > ws(j)) then
            styp     = rts(j)   ; w       = ws(j)
            rts(j)   = rts(j+1) ; ws(j)   = ws(j+1)
            rts(j+1) = styp     ; ws(j+1) = w
          end if
        end do
      end do

      ! Re-normalize weights
      w = sum(ws(1:ns))
      if (ldeb(spot)) write(usd,*) dpref, 'w',w,ws(1:ns)
      if (w < 0.9_sp .or. w > 1._sp + eps) then
        write(0,*) 'spot%hd%id=',spot%hd%id,' w=',w,1._sp + tiny(w),tiny(w),epsilon(w),1._sp+epsilon(w)
        do i = 1, ns ; write(0,*) proc,i,rts(i),ws(i) ; end do
        call finish(proc,'invalid surface type weights')
      end if
      ws(1:ns) = ws(1:ns) / w

      ! Calc. e_bg_ts
      e_bg_ts = 0._wp
      do i = 1, ns
        select case(rts(i))
        case(rts_land)
          e_bg_ts = e_bg_ts + e_bg_ts_land_sel**2 * ws(i)
        case(rts_sea)
          e_bg_ts = e_bg_ts + e_bg_ts_sea **2 * ws(i)
        case(rts_ice)
          e_bg_ts = e_bg_ts + e_bg_ts_ice **2 * ws(i)
        end select
      end do
      e_bg_ts = sqrt(e_bg_ts)

      ns = ns
    case default
      call finish('set_surface_type', 'invalid surf_class_vers')
    end select

    if (ldeb(spot)) write(usd,*) dpref,'set_surface_type result',rts(1:ns), e_bg_ts

    if (lchk) then
      l_check = (ns == tovs%ns)
      if (l_check) l_check = all(rts(1:ns) == tovs%rt_stype(1:ns))
      ! if (l_check) l_check = all(ws (1:ns) == tovs%w_styp  (1:ns))
      if (.not.l_check .and. ldeb(spot)) then
        write(*,*) dpref,proc,' check_failed spot',spot%hd%id,spot%col%c%dlat,spot%col%c%dlon,tovs%mw_stype
        write(*,*) dpref,proc,' check_failed orig',tovs%ns,tovs%rt_stype(1:tovs%ns),tovs%w_styp(1:tovs%ns)
        write(*,*) dpref,proc,' check_failed new ',ns,rts(1:ns),ws(1:ns)
      end if
    else
      tovs%rt_stype(1:ns) = rts(1:ns)
      tovs%w_styp  (1:ns) = ws (1:ns)
      tovs%ns             = ns
      tovs%e_bg_ts        = e_bg_ts
      stlsf = spot%stlsf
      stlsf = stlsf - iand(stlsf, SUR_ALL)
      do i = 1, ns
        select case(rts(i))
        case(rts_land)
          stlsf = ibset(stlsf, SUR_LAND)
          if (tovs%snf < 1._sp) stlsf = ibset(stlsf, SUR_NOSNOW)
          if (tovs%snf > 0._sp) stlsf = ibset(stlsf, SUR_SNOW)
        case(rts_sea)
          stlsf = ibset(stlsf, SUR_SEA)
          stlsf = ibset(stlsf, SUR_NOICE)
        case(rts_ice)
          stlsf = ibset(stlsf, SUR_SEA)
          stlsf = ibset(stlsf, SUR_ICE)
        end select
      end do
      spot%stlsf = stlsf
      if (ldeb(spot)) write(usd,*) dpref,'set_surface_type stlsf',stlsf
    end if

  end subroutine set_surface_type


  function valid (state)
  !------------------------------------------------------------
  ! check of observation is valid (i.e passed all 3dvar-checks)
  !------------------------------------------------------------
  integer(i1) ,intent(in) :: state
  logical                 :: valid
    valid = (state >= STAT_PASSIVE .and. state /= STAT_REJECTED)
  end function valid


  elemental subroutine destruct_t_jac(tj)
    type(t_jac), intent(inout) :: tj
    if (associated(tj% dt_tv))  deallocate(tj% dt_tv)
    if (associated(tj% dt_rh))  deallocate(tj% dt_rh)
    if (associated(tj% dq_tv))  deallocate(tj% dq_tv)
    if (associated(tj% dq_rh))  deallocate(tj% dq_rh)
    if (associated(tj% de_dpc)) deallocate(tj% de_dpc)
  end subroutine destruct_t_jac

  elemental subroutine destruct_t_jac_arr(tj)
    type(t_jac_arr), intent(inout) :: tj
    call destruct(tj%a)
  end subroutine destruct_t_jac_arr

  subroutine construct_t_jac(tj)
    type(t_jac), intent(out) :: tj
  end subroutine construct_t_jac

  subroutine get_tovs_var(vars, values, status, spt, obs, fg, ttovs, chan, instr)
    character(len=*),   intent(in)                      :: vars(:)   ! variable names
    real(kind=wp),      intent(out)                     :: values(:) ! values of variables
    integer,            intent(out)                     :: status
    type(t_spot),       intent(inout)                   :: spt       ! intent(out) only for storing lwp_qin
    type(t_obs),        intent(inout), optional         :: obs       ! intent(out) only for storing lwp_qin
    type(t_vector_segm),intent(in),    optional         :: fg        ! background
    type(t_tovs),       intent(in),    optional, target :: ttovs
    integer,            intent(in),    optional         :: chan
    integer,            intent(in),    optional         :: instr

    character(len=*), parameter :: proc = 'get_tovs_var'
    type(t_rad_set), pointer :: rs
    type(t_tovs),    pointer :: tt_          => null()
    integer,         pointer :: flag_  (:)   => null()
    real(tpp),       pointer :: av_    (:,:) => null()
    real(sp),        pointer :: l2c_   (:)   => null()
    real(wp),        pointer :: emis_  (:)   => null()
    integer,         pointer :: ci_    (:)   => null()
    real(sp),        pointer :: cldlev_(:)   => null()
    character(len=lv)        :: vars_(size(vars))
    character(len=lv)        :: v
    type(t_time)             :: time_
    integer                  :: lhgt(size(vars))   ! whether a height/channel might/must be added
    integer                  :: ihgt(size(vars))   ! Added height/channel
    integer                  :: nvar, ivar
    integer                  :: i, j, stat, ic_, ltr, nc, ibd
    integer                  :: il, lev_ps, il1, il2
    integer                  :: tovs_req, tovs_io, tovs_arg
    logical                  :: l_ps, l_obs, l_rs, l_lfr, l_ts, l_fg
    logical                  :: l_sea
    logical                  :: ld
    real(wp)                 :: p(mx_nlev)
    real(wp)                 :: plev, ps, y, y1, y2, w1, w2
    real(wp)                 :: land_frac, ts
    real(wp)                 :: l2p(1)
    ! lwp calculation
    integer                  :: i_lwp(4),chan_lwp(4)
    integer                  :: id_amsua,ii,instr_lwp
    ! mhs index calculation
    integer                  :: i_mhs(5)
    integer                  :: id_amsub


#define SET_ERR(str,stat) call error(str,stat);return

    ld = ldeb(spt)

    tovs_req = 0
    l_ps    = .false.
    l_obs   = .false.
    l_rs    = .false.
    l_lfr   = .false.
    l_ts    = .false.
    l_fg    = .false.
    nvar = size(vars)
    vars_(:) = vars(:)
    if (size(values)<nvar)then;SET_ERR('values array too small',15);endif
    lhgt(:)  = 0
    ihgt(:)  = 0
    do ivar = 1, nvar
      if (len_trim(vars_(ivar)) > lv) then;SET_ERR('variable name "'//trim(vars_(ivar))//'" too long.', 1);endif
      v = trim(tolower(vars_(ivar)))
      if (ld) write(usd,*) dpref,proc,'var=',trim(v)
      if (v      == 'iwv') then
        tovs_req   = ibset(tovs_req, TTOVS_AV_BIT)
        l_ps       = .true.
      elseif (v(1:3) == 'rh_' .or. &
              v(1:2) == 'q_'  .or. &
              v(1:2) == 't_') then
        tovs_req   = ibset(tovs_req, TTOVS_AV_BIT)
        l_ps       = .true.
        lhgt(ivar) = 2  ! Height must be added
      else if (v == 'sat_zen' .or. &
               v == 'sat_azi' .or. &
               v == 'sun_zen' .or. &
               v == 'sun_azi' .or. &
               v == 'rt_stype'.or. &
               v == 'mw_stype') then
        tovs_req = ibset(tovs_req, TTOVS_BASE_BIT)
      else if (v == 'ps') then
        l_ps       = .true.
      else if (v(1:4) == 'emis') then
        tovs_req   = ibset(tovs_req, TTOVS_EMIS_BIT)
        tovs_req   = ibset(tovs_req, TTOVS_FLAG_BIT)
        lhgt(ivar) = 1  ! Channel might be added
      else if (v(1:3) == 'l2c') then
        tovs_req   = ibset(tovs_req, TTOVS_L2C_BIT)
        lhgt(ivar) = 1
      else if (v(1:5) == 'p_l2c') then
        tovs_req   = ibset(tovs_req, TTOVS_L2C_BIT)
        tovs_req   = ibset(tovs_req, TTOVS_AV_BIT)
        lhgt(ivar) = 1
      else if (v(1:6) == 'plevel') then
        l_obs      = .true.
        lhgt(ivar) = 1
      else if (v == 'cldlev') then
        tovs_req   = ibset(tovs_req, TTOVS_CLDLEV_BIT)
      else if (v == 'p_cldlev') then
        tovs_req   = ibset(tovs_req, TTOVS_CLDLEV_BIT)
        tovs_req   = ibset(tovs_req, TTOVS_AV_BIT)
      elseif (trim(v) == 'lwp' .or. trim(v) == 'clw') then
        tovs_req   = ibset(tovs_req, TTOVS_CI_BIT)
        l_obs      = .true.
        l_rs       = .true.
        l_lfr      = .true.
        l_ts       = .true.
      elseif (v(1:7) == 'lwp_qin') then !mhs-based lwp
        tovs_req   = ibset(tovs_req, TTOVS_CI_BIT)
        l_obs      = .true.
        l_rs       = .true.
        l_lfr      = .true.
        ! l_fg       = .true.
      elseif (trim(v) == 'land_frac' .or. trim(v) == 'lfr') then
        l_lfr      = .true.
      elseif (trim(v) == 'ts' .or. trim(v) == 'tskin') then
        l_ts       = .true.
      end if
      if (lhgt(ivar) > 0) then
        ! Determine height/channel number
        j = index(v, '_', back=.true.)
        if (j > 0) then
          read(v(j+1:),*,iostat=stat) ihgt(ivar)
          if (stat == 0) then
            v = v(1:j-1)
            if (lhgt(ivar) == 2) v = trim(v)//'_'
          else if (present(chan)) then
            ihgt(ivar) = chan
          else if (.not.present(chan)) then
            SET_ERR('Failed to determine height/channel from variable "'//trim(v)//'".',11)
          endif
        elseif (present(chan)) then
          ihgt(ivar) = chan
        else
          SET_ERR('No height/channel in "'//trim(v)//'". Add _X (X=height/channel) or chan argument',4)
        end if
        v = trim(v)//'*'
      end if
      if (.not.any(v == c_var(1:n_var))) then;SET_ERR('variable "'//trim(v)//'" not implemented.',2);endif
      vars_(ivar) = trim(v)
    end do

    if (l_fg) then
      if (.not.present(fg)) then;SET_ERR('Optional argument "fg" is required for this request.',13);endif
    end if
    if (any(lhgt(1:nvar) == 1)) then
      tovs_req = ibset(tovs_req, TTOVS_CI_BIT)
      l_rs     = .true.
    end if
    if (l_rs .or. l_lfr .or. l_ts) then
      tovs_req = ibset(tovs_req, TTOVS_BASE_BIT)
    end if
    l_obs = l_obs .or. l_ps
    if (btest(tovs_req, TTOVS_AV_BIT)) then
      if (mx_nlev > mlev) then; SET_ERR('Not enough levels, Recompile with increased mlev.',12);endif
    end if

    tovs_io  = tovs_req
    if (tovs_req > 0 .and. present(ttovs)) then
      ! Load only the quantities that are not already loaded in ttovs:
      tovs_io  = tovs_io  - iand(ttovs%init, tovs_req) ! Quantities that have to be loaded
      tovs_arg = tovs_req - tovs_io ! Quantities that will be used from the supplied ttovs
    else
      tovs_arg = 0
    end if
    if (ld) write(usd,*) dpref,proc,' tovs_req/io/arg',tovs_req,tovs_io,tovs_arg

    if (l_obs .or. tovs_io > 0) then
      if (.not.present(obs)) then;SET_ERR('Optional argument "obs" is required for this request.',3);endif
    end if

    if (tovs_io > 0) then
      if (mx_nav > size(av,2)) call finish('get_tovs_prof',&
           'array "av" to small. Recompile with increased size')
      call load(obs, spt, tovs=tt_dum, av=av, flag=flag, emis=emis, l2c=l2c, ci=ci, cldlev=cldlev, tovs_io=tovs_io)
      if (btest(tovs_io, TTOVS_BASE_BIT  )) tt_     => tt_dum
      if (btest(tovs_io, TTOVS_AV_BIT    )) av_     => av
      if (btest(tovs_io, TTOVS_FLAG_BIT  )) flag_   => flag
      if (btest(tovs_io, TTOVS_EMIS_BIT  )) emis_   => emis
      if (btest(tovs_io, TTOVS_L2C_BIT   )) l2c_    => l2c
      if (btest(tovs_io, TTOVS_CI_BIT    )) ci_     => ci
      if (btest(tovs_io, TTOVS_CLDLEV_BIT)) cldlev_ => cldlev
    endif
    if (tovs_arg > 0) then
      if (btest(tovs_arg, TTOVS_BASE_BIT  )) tt_     => ttovs
      if (btest(tovs_arg, TTOVS_AV_BIT    )) av_     => ttovs%av
      if (btest(tovs_arg, TTOVS_FLAG_BIT  )) flag_   => ttovs%flag
      if (btest(tovs_arg, TTOVS_EMIS_BIT  )) emis_   => ttovs%emis
      if (btest(tovs_arg, TTOVS_L2C_BIT   )) l2c_    => ttovs%l2c
      if (btest(tovs_arg, TTOVS_CI_BIT    )) ci_     => ttovs%ci
      if (btest(tovs_arg, TTOVS_CLDLEV_BIT)) cldlev_ => ttovs%cldlev
    end if

    if (l_ps) then
      npv  = 2*tt_%nlev                    ! number of profile variables
      ! TODO: Why not ps_bg?
      ps = exp (obs% lev (spt% i% i+npv+2 ))
    end if

    if (btest(tovs_req, TTOVS_AV_BIT)) then
      if (tt_%i_p > 0) then
        p(1:tt_%nlev) = tt_%av(1:tt_%nlev,tt_%i_p)
      else
        p(1:tt_%nlev) = preslev(1:tt_%nlev)
      end if
    end if

    if (l_rs) then
      call get_tovs_rs(tt_, rs=rs)
      if (.not.associated(rs)) then;SET_ERR('Failed to find rad_set entry.',8);endif
    end if

    if (l_lfr) then
      land_frac = 0._wp
      do j = 1, n_styp
        if (tt_%rt_stype(j) == rts_land) land_frac = tt_%w_styp(j)
      end do
    end if

    if (l_ts) then
      if (tt_%ts > 0) then
        ts = tt_%ts
      else
        ts = spt% ts_bg
      end if
    end if

    tovs_io = 0 ! here for storing

    ! Cycle through variables in order to calc. the values
    do ivar = 1, nvar
      v = trim(vars_(ivar))

      if (lhgt(ivar) == 2) then
        plev = ihgt(ivar) * 100._wp
        plev = min(plev, ps)
        ! Prepare interpolation to level
        il2 = -1
        do il = 1, tt_%nlev
          if (p(il) >= plev) then
            il2 = il
            EXIT
          end if
        end do
        if (il2 < 2) then;SET_ERR('Failed to find pressure level '//trim(v(j+1:)),6);endif
        il1 = il2 - 1
        w1 = (p(il2) - plev) / (p(il2) - p(il1))
        w2 = (plev - p(il1)) / (p(il2) - p(il1))
      elseif (lhgt(ivar) == 1) then
        ic_ = -1
        do i = 1, tt_%nchan
          if (rs%chan(tt_%ci(i)) == chan) then
            ic_ = i
            exit
          end if
        end do
        if (ic_ <= 0) then;SET_ERR('Failed to find channel in t_tovs',9);endif
      end if

      ! Calc. the variable values
      select case(trim(v))
      case('iwv')
        lev_ps = -1
        do i = 1, tt_%nlev
          if (p(i) >= ps) then
            lev_ps = i
            exit
          end if
        end do
        if (lev_ps < 2) then;SET_ERR('Failed to find surface level.',7);endif
        ! Calc. surface q:
        y1 = (av_(lev_ps,2) + (p(lev_ps-1)-ps) + &
             av_(lev_ps-1,2) * (ps-p(lev_ps))) / (p(lev_ps-1) - p(lev_ps-2))
        y = y1 * (p(lev_ps-1) - ps)
        ! Integrate
        do il = lev_ps-2,1,-1
          y = y + 0.5_wp * (av_(il,2) + av_(il+1,2)) * (p(il+1)-p(il))
        end do
        y = -100._wp / gacc
      case('rh_*')
        y1 = rh_q(real(av_(il1,2), kind=wp), real(av_(il1,1), kind=wp), p(il1))
        y2 = rh_q(real(av_(il2,2), kind=wp), real(av_(il2,1), kind=wp), p(il2))
        y = w1 * y1 + w2 * y2
      case('q_*')
        y = w1 * av_(il1,2) + w2 * av_(il1,2)
      case('t_*')
        y = w1 * av_(il1,1) + w2 * av_(il1,1)
      case('sat_zen')
        y = tt_% saza
      case('sat_azi')
        y = tt_% boa
      case('sun_zen')
        y = tt_% sasoa
      case('sun_azi')
        y = tt_% soa
      case('rt_stype')
        y = tt_% rt_stype(1)
      case('mw_stype')
        y = tt_% mw_stype
      case('ps')
        y = ps
      case('lat')
        y = spt% col% c% dlat
      case('lon')
        y = spt% col% c% dlon
      case('phase', 'fov')
        y = spt% phase
      case('emis*')
        if (iand(flag_(ic_),EMIS_CALC) == 0) then;SET_ERR('Emissivity not available.',10);endif
        y = emis_(ic_)
      case('l2c*')
        y = l2c_(ic_)
      case('p_l2c*')
        l2p = l2c_(ic_)
        call lev2p(p,l2p)
        y = l2p(1)
      case('cldlev','p_cldlev')
        if (.not.present(chan).or..not.present(instr)) then;SET_ERR('chan/instr argument missing',14);endif
        ibd = -1
        instr_loop: do i = 1, rs%n_instr
          if (rs%instr(i) == instr) then
            do j = rs%o_ch_i(i)+1,rs%o_ch_i(i)+rs%n_ch_i(i)
              if (rs%chan(j) == chan) then
                ibd = rs%band(j)
                exit instr_loop
              end if
            end do
          end if
        end do instr_loop
        if (ibd < 1) then
          write(0,*) 'instr=',instr,' chan=',chan
          SET_ERR('did not find instr/chan',15)
        endif
        if (.not.associated(cldlev_) .or. tt_%nband < ibd) then
          write(0,*) 'init=',associated(cldlev_),' nband=',tt_%nband,' ibd=',ibd
          SET_ERR('invalid band info',16)
        end if
        y = cldlev_(ibd)
        if (v(1:2) == 'p_') then
          l2p = y
          call lev2p(p,l2p)
          y = l2p(1)
        end if
      case('plevel*')
        y = obs% body(spt%o%i + ic_)% plev
      case('z')
        y = spt% gp_bg
      case('mdlsfc')
        y = spt% mdlsfc
      case('loctime')
        time_ =localdate(spt% actual_time, spt%col%c%dlon)
        y = time_%secs
      case('land_frac','lfr')
        y = land_frac
      case('ts','tskin')
        y = ts
      case('tobs_tana','abs_tobs_tana')
        time_ = spt% hd% time - ana_time
        y = real(time_%days * 86400 + time_%secs,wp) / 60._wp
        if (v(1:3) == 'abs') y = abs(y)
      case('lwp','clw')
        i_lwp = -1
        ii = 1
        do i = 1, tt_%nchan
          do while (tt_%ci(i) > rs%o_ch_i(ii) + rs%n_ch_i(ii) .and. ii < rs%n_instr)
            ii = ii + 1
          end do
          id_amsua = amsua_chan_id(rs%instr(ii), rs%chan(tt_%ci(i)))
          select case(id_amsua)
          case(1:3)
            i_lwp(id_amsua) = i
            instr_lwp = rs%instr(ii)
          case(15)
            i_lwp(4)        = i
          end select
        end do
        y = -1._wp
        if (all(i_lwp > 0)) then
          chan_lwp = rs%chan(tt_%ci(i_lwp))
          if (chan_lwp(4) == 1) chan_lwp(4) = 15 ! 89GHz: Fake MHS1 as AMSUA15
          call mw_emiss(instr_lwp,chan_lwp, real(obs%body(spt%o%i+i_lwp(:))%o,wp), spt%stzen, land_frac, ts, stat, lwp=y)
          if (stat /= 0) y = -1._wp
        end if
        if (y < 0._wp) y = lwp_rain_thresh
      case('lwp_qin','lwp_qin_s','lwp_qinzou2016','lwp_qinzou2016_s')
        i_mhs = -1
        ii = 1
        do i = 1, tt_%nchan
          do while (tt_%ci(i) > rs%o_ch_i(ii) + rs%n_ch_i(ii) .and. ii < rs%n_instr)
            ii = ii + 1
          end do
          id_amsub = amsub_chan_id(rs%instr(ii), rs%chan(tt_%ci(i)))
          if (id_amsub > 0 .and. id_amsub <=5) then
            if (i_mhs(id_amsub) < 0) i_mhs(id_amsub) = i
          end if
        end do
        y = -1._wp
        if (all(i_mhs(1:2) > 0)) then
          ltr = len_trim(v)
          if (v(ltr-1:ltr) == '_s') then
            l_sea = .true.
          else
            l_sea = (land_frac < 0.1)
          end if
          if (.not.present(fg)) then;SET_ERR('Optional argument "fg" is required for this request.',13);endif
          y = LWP_QinZou2016(real(obs%body(spt%o%i+i_mhs(1:2))%o,wp), fg%x(spt%o%i+i_mhs(1:2)), l_sea)
          tt_%lwp_qin = real(y, kind=sp)
          tovs_io = ibset(tovs_io, TTOVS_BASE_BIT)
        elseif (tt_%lwp_qin >= 0._sp) then
          y = real(tt_%lwp_qin, kind=wp)
        end if
        ! if (y < 0._wp) y = lwp_qin_rain_thresh
        if (y < 0._wp) y = 2._wp ! rainy scene
      case default
        SET_ERR('Unknown variable: "'//trim(v)//'"',11)
      end select
      values(ivar) = y
      if (ld) write(usd,*) dpref,proc,' var=',trim(v),' val=',values(ivar)
    end do

    if (tovs_io /= 0) call store(obs, spt, tt_, tovs_io=tovs_io)

    SET_ERR('',0)

  contains

    subroutine error(str, stat)
      character(len=*), intent(in) :: str
      integer,          intent(in) :: stat
      status=stat
      errmsg=str
    end subroutine error

  end subroutine get_tovs_var

end module mo_tovs_prof
