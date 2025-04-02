!
!+ Utilities for the PSAS system.
!
MODULE mo_psasutil
!
! Description:
!   Utilities for the PSAS system.
!   This module knows about:
!     - the atmospheric model data type (module mo_atm_state)
!     - block decomposed matrices       (module mo_dec_matrix)
!       and vectors
!     - the observation data type       (module t_obs)
!   It provides the following conversion routines:
!     - interpolate : interpolate atmospheric state to observation space
!     - get_obs     : get observed quantities from observation data type
!     - put_bg      : store background          in observation data type
!
! Current Code Owner: DWD, Andreas Rhodin
!    phone: +49 69 8062 2722
!    fax:   +49 69 8062 3721
!    email: andreas.rhodin@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_1         2008/11/05 Andreas Rhodin
!  First operational 3D-Var release
! V1_4         2009/03/26 Andreas Rhodin
!  changes for verification mode (TSK_YH)
! V1_8         2009/12/09 Andreas Rhodin
!  for COSMO get geop.height from full levels
!  changes for interpolation of IFS stratospheric background
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Andreas Rhodin
!  various technical changes
! V1_15        2011/12/06 Andreas Rhodin
!  option to interpolate stratospheric humidity in tv,q instead of tv,rh
! V1_19        2012-04-16 Andreas Rhodin
!  new subroutine get_plev: get pressure level from observation data type
! V1_22        2013-02-13 Andreas Rhodin
!  implementation of COSMO and GPSGB operators (no 3dvar, verification only)
! V1_27        2013-11-08 Andreas Rhodin
!  constant extrapolation of T above model top if IFS profile is not provided
! V1_29        2014/04/02 Andreas Rhodin
!  use PC fg coefficients for IR emissivity calculations
! V1_31        2014-08-21 Andreas Rhodin
!  consider snow fraction (derived from snow height) in IR emissivity model
! V1_35        2014-11-07 Andreas Rhodin
!  set hs_bg (snow height) when interpolating for COSMO
! V1_37        2014-12-23 Andreas Rhodin
!  implement RTTOV (Rochon) interpolation for model first guess (nwv_rad=4)
! V1_42        2015-06-08 Andreas Rhodin
!  implement temporal interpolation for COSMO MEC
! V1_44        2015-09-30 Harald Anlauf
!  Preparations for Jason-2
! V1_46        2016-02-05 Harald Anlauf
!  Fix DRIBU station id; extend SYNOP/DRIBU code for monitoring of T2m
! V1_47        2016-06-06 Harald Anlauf
!  STD: implementation of adjoint code
! V1_48        2016-10-06 Andreas Rhodin
!  changes for rh2m, t2m
!  Add roughness length z0 to t_spot, t_sl (H.Anlauf)
! V1_49        2016-10-25 Harald Anlauf
!  interpolate: derive z0_bg from sea points only
! V1_51        2017-02-24 Andreas Rhodin
!  use generalised humidity for GPSGB, GPSRO
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Andreas Rhodin  MPIfM/DWD  2001-2008  original code
!==============================================================================

#if defined(__ICON__) && !defined(__USE_RTTOV)
#undef _RTTOV_VERSION
#endif

  !=============
  ! modules used
  !=============
  !-------------------------
  ! general purpose routines
  !-------------------------
  use mo_kind,          only: wp, sp, i8        ! precision kinds
  use mo_exception,     only: finish            ! abort subroutine
  use mo_mpi_dace,      only: dace,            &! MPI group info
                              p_bcast,         &! generic MPI broadcast
                              p_sum,           &!
                              p_max,           &!
                              p_or
  use mo_run_params,    only: aux,             &! path for auxiliary files
                              input, output,   &!
                              path_file,       &! add path to file name
                              nproc1, nproc2,  &!
                              interp_strato,   &!
                              file_strato,     &!
                              grid_file_strato,&!
                              ana_time,        &!
                              fc_ref_time
  use mo_fortran_units, only: get_unit_number, &! reserve a unit number
                              return_unit_number! release the unit number
  use mo_physics,    only: gacc,        &! gravity acceleration
                           rearth,      &! earth radius
                           R,           &! gas constant
                           lapse_cl,    &! Climatological lapse rate [K/m]
                           ppmv_dry_q,  &!
                           q_ppmv_dry,  &!
                           rh_q,        &! rh -> q
                           t_tv_q,      &! tv,q -> t
                           tv_t_q,      &! t,q -> tv
                           tv_t_rh,     &! t,rh -> tv
                           rh_q,        &! rh -> q
                           tq_tvrh_vec, &!
                           tq_tvrh_vec2,&!
                           gas_id_o3,   &!
                           gas_id_co2,  &!
                           ppmv_dry_trg,&!
                           q_ppmv_dry_trg,&!
                           get_gas
  use mo_cntrlvar,   only: gh_rh,       &! generalized humidity from rel.hum.
                           tq_tvgh,     &! Tv,gh -> T,q
                           tq_tvgh_vec   ! Tv,gh -> T,q
  use mo_t_col,      only: t_cols,      &! atmospheric columns data type
                           t_col,       &! component of t_cols
                           get_cols,    &! collect columns
                           alloc_cols,  &! t_cols allocation routine
                           dealloc_cols,&! t_cols deallocation routine
                           COL_T,       &! temperature
                           COL_TV,      &! virtual temperature
                           COL_Q,       &! specific humidity
                           COL_RH,      &! relative humidity
                           COL_UV,      &! hor.wind (u,v)
                           COL_P,       &! pressure
                           COL_PH,      &! pressure     (half levels)
                           COL_GEO,     &! geopotential
                           COL_GEOH,    &! geopotential (half levels)
                           COL_X,       &! water load
                           COL_QCL,     &! cloud water content
                           COL_QCI,     &! cloud ice   content
                           COL_QR,      &! rain
                           COL_QS,      &! snow
                           COL_QG,      &! graupel
!                          COL_INFLAT,  &! inflation factor
!                          COL_TV2,     &! tv from tv, not t,q
!                          COL_PP,      &! pressure deviation
!                          COL_W,       &! vertical velocity
!                          COL_SOIL,    &! soil temperature, humidity
                           COL_CLC,     &! cloud cover
!                          COL_RANGE,   &! variables with time range indicator
                           COL_QVDIA,   &! specific humidity, diagnostic
                           COL_QCDIA,   &! spec. cloud water, diagnostic
                           COL_QIDIA,   &! spec. cloud ice, diagnostic
                           COL_OZONE,   &! O3 concentration
                           COL_CO2,     &! CO2 concentration
                           COL_REFF_QC, &! eff. radius cloud water
                           COL_REFF_QI, &! eff. radius cloud ice
                           u10_use_mlevel! lowest model level for 10m wind?
  !----------------------------
  ! acess observation data type
  !----------------------------
  use mo_t_obs,      only: t_obs,       &! observation data type
                           t_spot,      &! observation meta data type
                           t_mcols,     &!
                           t_hic,       &!
                           int_vh,      &! interpolation: 1.vert.,2.hor.
                           int_nn,      &! horiz. interp.: nearest neighbour
                           int_rad_hum, &! humidity variable for vertical interp. for radiances
                                         ! 0: RH, 1:QV, 2:RH below int_rad_hum_pblend + Q above
                                         ! 3: Q/T as in RTTOV - not necessary but useful for
                                         !    comparisons with e.g. radsim
                           int_rad_hum_pblend,&! Blending levels for radiance humidity interpolation
                           vint_lin_t,  &! linear vert.interp. for temperature
                           vint_lin_uv, &! linear vert.interp. for wind
                           vint_lin_z , &! linear vert.interp. for geopotential
                           vint_lin_tov,&! linear vert.interp. for RTTOV
                           nwv_rad,     &! radiance vert.intp. flag
                           vint_rttov,  &! RTTOV interpolation mode
                           mdlsfc_frac, &! model surface flags from fractions
                           invalid,     &! invalid value in t_spoz
                           TSK_YH,      &! run linear or forward operator
                           TSK_K,       &! evaluate linear operator
                           OBS_TV,      &! virtual temperature flag
                           OBS_T,       &!         temperature flag
                           OBS_RH,      &! relative humidity   flag
                           OBS_Q,       &!
                           OBS_U,       &! wind component u    flag
                           OBS_V,       &! wind component v    flag
                           OBS_FF,      &! wind speed          flag
                           OBS_H,       &! geopotential height flag
                           OBS_HS,      &! surface geop height flag
                           OBS_DUM,     &! dummy (sink) variable
                           TEMP,        &! TEMP  observation type identifier
                           TOVS,        &! ATOVS observation type identifier
                           SYNOP,       &! SYNOP observation type identifier
                           AMV,         &! AMV   observation type identifier
                           AIREP,       &! AIREP observation type identifier
                           GPSRO,       &! GPS RO (ray-tracer)    identifier
                           GPSGB,       &! GPS ground based
                           COSMO,       &! COSMO observation operators
                           SOIL,        &! ASCAT soil operator
                           SATEM,       &! SATEM  identifier
                           WLIDAR,      &! Wind Lidar identifier
                           SCATT,       &! SCATT observation type identifier
                           destruct,    &! t_obs destructor routine
                           release_mem, &! release memory
                           ldeb,        &! debug selected spot(s)
                           usd,dpref     ! ...
  use mo_t_use,      only: STAT_DISMISS  ! do not apply operator flag
  use mo_occ,        only: occ_col2xi    ! convert t,q ,gp  to  t,rh,gh
  use mo_std,        only: std_col2xi    ! convert t,q ,gp  to  t,rh,gh
  use mo_set_matrix, only: set_H,       &! extract linear obs.operator
!                          scatter_K,   &!
                           Pb_times_z    ! vector matrix product Pb * z
  use mo_synop,      only: dtdzp,       &! tempert.gradient for p extrapolation
                           zbs,         &! temperature extrapolation
                           version,     &! observation operator version
                           use_ps_model  ! invalid value: use model surf.pres.
  use mo_obs,        only: process_obs   ! generic observation handling routine
  use mo_obs_sndrcv, only: p_send,      &! generic MPI send
                           p_recv        ! generic MPI recv
  use mo_obs_set,    only: t_obs_set     !
  use mo_time,       only: t_time,      &! Data type to hold date and time
                           frac_year,   &!
                           operator(-), &!
                           operator(>), &!
                           operator(>=),&!
                           operator(<=),&!
                           operator(==),&!
                           cyyyymmddhhmm,&!
                           cyyyymmddhhmmss,&!
                           time_c,      &!
                           invalid_time,&!
                           days,        &!
                           date_tmpl
  use mo_t_bg_err_op,only: covm          !
  use mo_bg_err_2d,  only: apply_B_ii_2d
  use mo_rttov,      only: jplev,       &! number of rttov levels
                           lnp,         &! log pressure levels for rttov
                           preslev       ! pressure levels for rttov
  use mo_rtifc,      only: rtifc_coef_prop, &!
                           rt_gas_id_o3,&!
                           rt_gas_id_co2
  use mo_fdbk_tables,only: OT_SYNOP,    &! SYNOP observation type
                           OT_DRIBU,    &! BUOY  observation type
                           OT_RAD,      &! RADiances observ. type
                           VN_T2M,      &! 2m temperature code
                           VN_TD2M,     &! 2m dewpoint t. code
                           VN_RH,       &!    rh code
                           VN_RH2M,     &! 2m rh code
                           VN_HOSAG      ! height of sensor above surface code
  use mo_tovs_prof,  only: fg_prof_top   ! handling of TOVS profiles above level_*_dum
  use mo_tovs,       only: trg_file,    &! File with (current) tracegas info
                           trg_clim_file,&!
                           trg_clim_trend,&!
                           t_trg_clim_trend,& !
                           trg_hist_file, &!
                           trg_hist_inidate, &
                           glob_use_o3, &!
                           glob_use_co2,&!
                           use_reff,    &!
                           p_blend_ext, &!
                           nld_h_max,   &!
                           nld_t_max,   &!
                           iatm_tovs_all,&!
                           COL_CLD_TOVS
  use mo_t_tovs,     only: t_tovs,      &! type to store radiance specific stuff
                           store,       &! store t_tovs
                           load,        &! store t_tovs
                           destruct,    &! destruct t_tovs
                           add_av_cont, &! add entry to av_cont
                           get_tovs_rs, &! get rad_set corresponding to t_tovs
                           mx_nav,      &! highest value of t_tovs%n_av
                           mx_nlev,     &! highest value of t_tovs%nlev
                           tpp,         &! t_tovs%av precision
                           TTOVS_BASE,  &! I/O flag for t_tovs
                           TTOVS_AV,    &! I/O flag for t_tovs
                           TTOVS_AV_BIT,&! I/O flag for t_tovs
                           TTOVS_ALL     ! All I/O flags set for t_tovs
  use mo_rad,        only: rad_set,     &
                           n_set,       &
                           m_instr,     &
                           t_rad_set
  use mo_grib_handling,&
                     only: t_inventory,       &!
                           get_inventory,     &!
                           print_inventory     !
  use mo_grib,       only: read
  use mo_atm_state,  only: t_atm,             &! atmospheric state data type
                           construct,         &! t_grid initialisation routine
                           destruct,          &! t_grid cleanup routine
                           print,             &!
                           set_geo,           &!
                           set_p,             &!
                           set_ph,            &!
                           set_tv,            &!
                           allocate
  use mo_memory,     only: print
  use mo_atm_grid,   only: t_grid,            &! atmospheric grid data type
                           destruct
  use mo_grid_intpol,only: idx_init
  use mo_namelist,   only: position_nml,      &!
                           POSITIONED

  !------------------------
  ! access matrix data type
  !------------------------
  use mo_dec_matrix, only: t_vector,       &! vector data type
                           t_ivector,      &! integer vector data type
                           t_matrix,       &! matrix data type
                           construct,      &! allocate matrix block data type
                           destruct,       &! deallocate matrix block data type
                           gather,         &! gather vector content on PE
                           release_mem,    &! release memory of vector
                           sqrt,           &! square root of a matrix
                           operator (*),   &! matrix vector multiplication
                           operator (+),   &! matrix and vector operations
                           operator (-),   &! matrix and vector operations
                           assignment (=), &! assignment
                           random_gauss     ! normal distribution
  !------------------------------------
  ! Interface to interpolation routines
  !------------------------------------
  use mo_algorithms,   only: init_splinex  ! calculate coefficients (2nd deriv)
  use interpolation,   only: spline        ! perform spline interpolation
  use data_constants,  only: pi
  use mo_transform,    only: fft

  implicit none

  !-------------------------------------------------------------------
  ! col_field_pointer - pointer to one field in one column
  !-------------------------------------------------------------------
  ! fid  - field ID of atmospheric field, see COL_XX in mo_t_col.f90
  ! pt   - pointer to one column of an atmospheric field
  !-------------------------------------------------------------------
  type col_field_pointer
     integer(i8)        :: fid = -1_i8     ! field ID
     real (wp), pointer :: pt(:) => Null() ! pointer to field
  end type col_field_pointer

  !---------------------------------------------------------------------
  ! col_imc - fields in neighboring model columns: infos and pointer
  !---------------------------------------------------------------------
  ! ncol             - number of neighboring columns used for interpolation,
  !                    e.g. ncol=3 for ICON
  !                    maxcol = size(spt% col% h% imc,1) is the maximum
  !                    number of neighboring columns which is supported.
  ! nmax             - max. number of fields in any column
  ! nall             - number of fields in all columns
  ! nfield(1:maxcol) - number of fields per column
  ! ids(1:maxcol)    - bit field with info about fields in column
  ! idall            - bit field with fields available in all columns
  ! fidall(1:nall)   - field IDs for each bit set in "idall",
  !                    size(fidall)=nall = number of fields in all columns
  ! colp(1:maxcol,1:nmax) - pointer to fields in columns
  !                         colp(j,i): j - column index, j=1, ..., ncol
  !                                    i - field index, i=1, ..., nfield(j)
  ! idxall(1:ncol,1:nall) - index of given field in all columns
  !                         idxall(j,f): j - column index, as in "colp"
  !                                      f - field index, f=1, ..., nall
  !                         i=idxall(j,f)  - field index for "colp"
  !---------------------------------------------------------------------
  type col_imc
     integer                               :: ncol = 0
     integer                               :: nmax = 0
     integer                               :: nall = 0
     integer,                  allocatable :: nfield(:)
     integer(i8),              allocatable :: ids(:)
     integer(i8)                           :: idall = 0_i8
     integer(i8),              allocatable :: fidall(:)
     type (col_field_pointer), allocatable :: colp(:,:)
     integer,                  allocatable :: idxall(:,:)
  end type col_imc

  real(wp)              :: pblend(2)                    ! blending levels for int_rad_hum=2

!------------------------------------------------------------------------------
  !================
  ! public entities
  !================
  private
  public :: interpolate    ! interpolate atmospheric state to observation space
  public :: get_obs        ! get observed quantities from observation data type
  public :: get_plev       ! get pressure level from observation data type
  public :: get_obstype    ! get observation type    from observation data type
  public :: get_varno      ! get variable number     from observation data type
  public :: get_id         ! observation id from source/record
  public :: put_bg         ! store background          in observation data type
  public :: set_obs_random ! set random observations
  public :: check_ob_fg    ! check against background
  public :: write_intp     ! write information on interpolation space
  public :: set_new_x      ! derive x = Pb Ht z  .or.  x = H^-1 y
  public :: set_new_H      ! derive tl model
  public :: interpolate_trg! read and interp. atm. vars. from ext. sources
  public :: interpolate_strat   ! read and add/blend extra levels above model top
!==============================================================================
contains
!==============================================================================
  subroutine interpolate (col, obs, psi, bg)
  type (t_cols)   ,intent(in)              :: col(:) ! atmospheric columns
  type (t_obs_set),intent(inout)           :: obs    ! observations
  type (t_vector) ,intent(inout) ,optional :: psi    ! interpolated values
  logical                        ,optional :: bg     ! true for background

    !================
    ! local variables
    !================
    logical               :: lbg       ! true for background
    integer               :: l         ! vector segment index
    integer               :: j         ! spot           index
    integer               :: k         ! element        index
    integer               :: i         ! element        index
    integer               :: k2        ! rttov level index
    integer               :: ic        ! model column   index
    integer               :: ii        ! model column   index
    integer               :: it        ! model column   index (time)
    integer               :: nt        ! time slices (one or two)
    integer               :: n         ! number of points used for interpolation
    integer               :: ke        ! grid point     index
    integer               :: ml(1)     ! index of nearest neighbour
    integer               :: nl        ! number of profile levels (TOVS)
    real(wp)              :: z         ! last level (log p) processed
    real(wp) ,allocatable :: int  (:)  ! interpolated fields
    logical  ,allocatable :: lini (:)  ! spline coefficient initialisation flag
    integer  ,allocatable :: rti  (:)  ! RTTOV interpolation flag (see rti_* below)
    real(wp) ,allocatable :: rt_in (:,:)  ! RTTOV interpolated temperature
    real(wp) ,allocatable :: rt_out(:,:,:)! RTTOV interpolated relative humidity
    real(wp) ,allocatable :: dum_rt(:,:)  ! RTTOV humidity conversion
    real(wp) ,allocatable :: rt_wgt(:)    ! RTTOV humidity conversion
    real(wp) ,allocatable :: p(:)         ! Pressure profile
    type(t_cols)          :: d2        ! second derivative of spline
    type(t_col) ,pointer  :: c, d      ! temporary
    logical               :: lpe       ! flag: process segment on this pe
    type(t_spot) ,pointer :: spt       ! pointer to meta data element
    real(wp) ,allocatable :: cp  (:)   ! horizontally interpolated profiles: p
    real(wp) ,allocatable :: ctv (:)   ! horizontally interpolated profiles: t
    real(wp) ,allocatable :: crh (:)   ! horizontally interpolated profiles: rh
    real(wp) ,allocatable :: dtv (:)   ! horizont. interp. prof. 2nd deriv. t
    real(wp) ,allocatable :: drh (:)   ! horizont. interp. prof. 2nd deriv. rh
    logical               :: lin_t     ! local copy of vint_lin_t, vint_lin_tov
    integer               :: modtype   ! module type
    type(t_tovs)          :: ttovs     ! required to store background profile for radiances
    logical               :: l_tovs    ! store profile in t_tovs%av
    integer,  allocatable :: i_fail(:) ! exit status from tq_tvgh
    integer               :: i_rs      ! index in rad_set array
    logical               :: l_vint_rt ! perform vertival interpolation for radiance spot
    integer               :: i_p       ! index of p in idlist
    logical               :: l_cloud   ! cloud variables for TOVS
    logical               :: l_shift_cld ! shift cloud variable one "level" up. Useful for
                                         ! RTTOV, which expects clouds on layers (not levels)
    integer               :: nl_st     ! number of extrap. levels (above model top)
    logical               :: l_atm_lev ! whether level is in bg
    integer               :: i0,i1
    logical               :: l_debug
    !-------------------------------------
    ! parameters to interpolation routines
    !-------------------------------------
    real(wp)              :: tp            ! interpolated temperature
    real(wp)              :: tvp           ! interpolated virtual temp.
    real(wp)              :: qp            ! interpolated specific humidity
    real(wp)              :: rp            ! interpolated relative humidity
    real(wp)              :: up            ! interpolated wind component
    real(wp)              :: vp            ! interpolated wind component
    real(wp)              :: gp            ! interpolated geopotential
    real(wp)              :: t2p           ! interpolated 2m temperature
    real(wp)              :: t2l           ! interpolated 2m temperature  (land)
    real(wp)              :: r2l           ! interpolated 2m rel.humidity (land)
    real(wp)              :: t2c           ! height-adjusted 2m temperature
    real(wp)              :: r2c           ! height-adjusted 2m rel.humidity
    real(wp)              :: dt            ! temperature gradient (in ln p)
    real(wp)              :: gsp           ! interp. geopotential
    real(wp)              :: psp           ! interp. pressure
    real(wp)              :: slp           ! interp. sea land mask
    real(wp)              :: fip           ! interp. sea ice fraction
    real(wp)              :: z0p           ! interp. roughness length
    real(wp)              :: sdp           ! interp. SSO standard deviation
    real(wp)              :: tlp           ! interp. temp. at lowest level
    real(wp)              :: tsp           ! interp. temp. at surface
    real(wp)              :: rlp           ! interp. rel.hum. at lowest level
    real(wp)              :: dtz           ! dt/dz 850..700 hPa
    real(wp)              :: hsp           ! interp. snow height
    real(wp)              :: snf           ! interp. snow fraction
    real(wp)              :: twp           ! interp. sea/water temperature
    real(wp)              :: ulp           ! interp. u-wind at lowest level
    real(wp)              :: vlp           ! interp. v-wind at lowest level
    real(wp)              :: u10
    real(wp)              :: v10
    real(wp)              :: sws           ! sum of weights, sea points
    real(wp)              :: hs            ! sensor height above ground
    real(wp)              :: w1            ! interpolation weight (combined)
    real(wp)              :: w2            ! interpolation weight (2m/10m)
    real(wp)              :: wt            ! time interpolation coefficient
    real(wp)              :: int_          ! interpolated value
    ! interpolation of fields required for RTTOV
    integer               :: fi            ! index of field in column
    integer               :: nid           ! number of allocated vars in t_col
    integer(i8)           :: iatm          ! COL_* required for tovs
    logical               :: first_tovs = .true.
    logical               :: l_tt
    real(wp), allocatable :: tt_tv(:), tt_rh(:), tt_t(:), tt_q(:)
    real(wp)              :: rh_
    type (col_imc)        :: afp
    integer               :: fidx
    integer               :: rti_          ! rti and rti_ contain nl + RTI_*
    integer, parameter    :: RTI_CLOUD   = ibset(0,29)
    integer, parameter    :: RTI_VINT_RT = ibset(0,30)

    !======================
    ! executable statements
    !======================
    lbg = .false. ;if (present(bg)) lbg = bg

    ! Required for radiances:
    mx_nlev   = p_max(mx_nlev)
    nld_t_max = p_max(nld_t_max)
    nld_h_max = p_max(nld_h_max)
    if (fg_prof_top > 0 .and. (nld_t_max > 0 .or. nld_h_max > 0)) then
      allocate(tt_tv(mx_nlev),tt_rh(mx_nlev),tt_t(mx_nlev),tt_q(mx_nlev),i_fail(mx_nlev))
    end if

    !--------------------------
    ! loop over vector segments
    !--------------------------
    do l = 1, size(obs% o)
      !---------------------------------------
      ! handle only segments on this processor
      !---------------------------------------
      lpe = .false.
      if (present(psi)) lpe = lpe .or. (psi% s(l)% pe == dace% pe)
      if (.not. lpe) cycle
      ke = col(l)% ke
      !----------------------------------------------------
      ! check for correct size of interpolation space array
      !----------------------------------------------------
      if (present(psi)) then
        if (size(obs%o) /= psi% n_s) then
          write (6,*)  'interpolate:  nbox /= psi% n_s :',size(obs%o), psi% n_s
          call finish ('interpolate','nbox /= psi% n_s')
        endif
        if (obs%o(l)%n_int /= psi%s(l)% n) then
          write (6,*)  'interpolate:  n_int /= psi% n :',obs%o(l)%n_int, psi%s(l)% n
          call finish ('interpolate','n_int /= psi% n')
        endif
      endif
      !---------------------
      ! allocate temporaries
      !---------------------
      call alloc_cols (d2, tmp=col(l))
      allocate (lini (col(l)% ncol))
      lini = .false.
      !-----------------------------------------------------
      ! temporary arrays for radiance vertical interpolation
      !-----------------------------------------------------
      if (nwv_rad ==4) then
        if (int_rad_hum == 2) then
          allocate(rt_wgt(mx_nlev))
          if ( int_rad_hum_pblend(1) > int_rad_hum_pblend(2) .and. &
               int_rad_hum_pblend(2) > 0._wp) then
            pblend = log(int_rad_hum_pblend*100._wp)
          else
            call finish('interpolate','invalid values for int_rad_hum_pblend')
          endif
        end if
        if (int_rad_hum < 0 .or. int_rad_hum > 3) call finish('interpolate', &
             'invalid value for int_rad_hum (0<=int_rad_hum<=3).')
        if (int_rad_hum > 0) allocate(dum_rt(mx_nlev,2))
      endif
      first_tovs = .true.

      !-----------------------
      ! loop over observations
      !-----------------------
!NEC$ nomove
      do j = 1, obs% o(l)% n_spot
        spt => obs% o(l)% spot(j)
        l_debug = ldeb(spt)
        if (l_debug) write(usd,*) dpref//' int_spot'
        !-------------------
        ! only valid reports
        !-------------------
        if (spt% use% state <= STAT_DISMISS) then
          if (present (psi)) then
            psi% s(l)% x(spt% i% i + 1:spt% i% i + spt% i% n) = invalid
          endif
          cycle
        end if
        !-------------------
        ! time interpolation
        !-------------------
        nt = 2; if (spt% w_time==0) nt = 1
        !-----------------------------------------------
        !  allocate horizontally interpolated columns
        !-----------------------------------------------
        allocate (int (spt% i% n))
        if (.not. int_vh) then
          allocate (cp  (col(l)% ke))
          allocate (ctv (col(l)% ke))
          allocate (crh (col(l)% ke))
          allocate (dtv (col(l)% ke))
          allocate (drh (col(l)% ke))
        end if
        !============================================
        ! all observations: set entries in spot table
        !============================================
        gsp = 0._wp; psp = 0._wp; slp = 0._wp; tlp = 0._wp; tsp = 0._wp
        u10 = 0._wp; v10 = 0._wp; fip = 0._wp; hsp = 0._wp; snf = 0._wp
        dtz = 0._wp; z0p = 0._wp; sdp = 0._wp; twp = 0._wp; ulp = 0._wp
        vlp = 0._wp
        !--------------------------------------------------------------
        ! patch weights for horizontal nearest neigbour 'interpolation'
        !--------------------------------------------------------------
        if ( int_nn .and. spt% col% h% imc(1,1) /= 0) then
          ml = maxloc(spt% col% h% w)
          spt% col% h% w        = 0._wp
          spt% col% h% w(ml(1)) = 1._wp
          !RF: !!DEBUG!! For comparison with datool
          !spt% col% c% dlat = col(l)% col(spt% col% h% imc(ml(1),1))%c%dlat
          !spt% col% c% dlon = col(l)% col(spt% col% h% imc(ml(1),1))%c%dlon
        endif
        !---------------------------------------------------------------------------
        ! special handling for tovs: background for dummy variables stores in t_tovs
        !---------------------------------------------------------------------------
        l_tovs = lbg .and. spt%hd%modtype == TOVS ! .and. interp_strato <= 0)
        l_tt = l_tovs
        if (spt%hd%modtype == TOVS) then
          call load(obs%o(l), spt, ttovs, i_rs=i_rs)
          if (l_debug) write(usd,*) dpref,'int rs',rad_set(i_rs)%satid,rad_set(i_rs)%grid
          l_cloud     = any(rad_set(i_rs)%iopts(1:rad_set(i_rs)%n_instr)%cloud_mode >= 1)
          l_vint_rt   = rad_set(i_rs)%gopts%lev_mode <= 0
          l_shift_cld = any(rad_set(i_rs)%iopts(1:rad_set(i_rs)%n_instr)%cloud_mode == 3)
          if (l_debug) write(usd,*) dpref,'int l_tt',l_tt,ttovs%i_t,ttovs%i_q
          l_tt = l_tt .and. fg_prof_top > 0 .and. &
                            (ttovs%nld_t > ttovs%nl_st .or. ttovs%nld_h > ttovs%nl_st)
          if (l_debug) write(usd,*) dpref,'int l',l_cloud,l_vint_rt,l_tt
          nl_st = ttovs%nl_st
          nl    = ttovs%nlev - nl_st

          if (.not.l_vint_rt .and. nwv_rad < 4) call finish('interpolate',&
               'nwv_rad < 4 only supported with radiances on RTTOV-levels.')
          if (first_tovs) then
            first_tovs = .false.
            if (nwv_rad >= 4) then
              nid = popcnt(iatm_tovs_all)
              if (int_vh) then
                allocate (rt_out(nid,mx_nlev,0:col(l)% ncol))
                allocate (rti(col(l)% ncol))
                rti = 0
              else
                allocate (rt_out(nid,mx_nlev,0:0))
              endif
              if (.not.int_vh .or. int_rad_hum==3) allocate (rt_in(nid,ke))
              allocate(p(mx_nlev))
            else
              allocate(p(jplev))
              p(1:jplev) = preslev(1:jplev)
            end if
          end if
        end if

        !------------------------------
        ! interpolate background fields
        !------------------------------
        if (spt% col% h% imc(1,1) == 0) then
          !----------------------------
          ! no columns specified (MCOLS)
          !----------------------------
          if (lbg) then
            select case (spt% hd% modtype)
            case (COSMO, SOIL, GPSGB)
              !-------------------------------------------
              ! one column for COSMO observation operator
              ! no time interpolation so far for ITY_MCOLS
              !-------------------------------------------
              ii = spt% imcol(1)% imc(1)
              c => col(l)% col(ii)
              spt% ps_bg   = c% s% ps      ! surface pressure
              spt% gp_bg   = c% s% geosp   ! surface geopotential
              spt% sl_bg   = c% s% lsm     ! sea-land mask
              spt% fi_bg   = c% s% fr_ice  ! sea ice fraction
              spt% z0_bg   = c% s% z0      ! roughness length
              spt% ssd_bg  = c% s% ssd     ! SSO standard deviation
              spt% tl_bg   = c% s% tll     ! temp. at lowest level
              spt% ts_bg   = c% s% tsurf   ! temp. at surface
              spt% tw_bg   = c% s% t_so_0  ! sea/water temperature
              spt% dtdz_bg = c% s% dtdz    ! lapse rate 850..700 hPa
              spt% hs_bg   = c% s% h_snow  ! snow height
              spt% u10bg   = c% s% u10m    ! 10m wind
              spt% v10bg   = c% s% v10m    ! 10m wind
              spt% mdlsfc  = mdlsfc_frac (spt)
              spt% pz_bg   = gacc * c% s% ps / (R * c% s% tll)  ! d ps / d z
            case default
              !----------------------------------
              ! fill in reasonable values (GPSRO)
              !----------------------------------
              spt% sl_bg = 0.5  ! sea-land mask
              spt% fi_bg = 0.5  ! ice fraction
              spt% gp_bg = 0.0  ! surface geopotential
            end select
          endif
        else
          sws = 0._wp
          do it = 1, nt
          do ic = 1, size(spt% col% h% imc,1)
            ii = spt% col% h% imc(ic,it)
            if (ii==0) exit
            wt = spt% w_time
            if (it==1) wt = 1._wp - wt
            w1 = wt * spt% col% h% w(ic)    ! weight for this grid point/time
            c => col(l)% col(ii)
            psp = psp + w1 * c% s% ps
            gsp = gsp + w1 * c% s% geosp
            slp = slp + w1 * c% s% lsm
            fip = fip + w1 * c% s% fr_ice
            tlp = tlp + w1 * c% s% tll
            tsp = tsp + w1 * c% s% tsurf
            dtz = dtz + w1 * c% s% dtdz
            hsp = hsp + w1 * c% s% h_snow
            snf = snf + w1 * c% s% snowc
            ulp = ulp + w1 * c% s% ull
            vlp = vlp + w1 * c% s% vll
            u10 = u10 + w1 * c% s% u10m
            v10 = v10 + w1 * c% s% v10m
            sdp = sdp + w1 * c% s% ssd
            !------------------------------------------------------
            ! Interpolate roughness length over (mostly) sea points
            ! Interpolate sea temperature  over (mostly) sea points
            !------------------------------------------------------
            if (c% s% lsm < 0.05_wp) then
               sws = sws + w1
               z0p = z0p + w1 * c% s% z0
               twp = twp + w1 * c% s% t_so_0
            end if
          end do
          end do
          if (sws > 0._wp) then
             z0p = z0p / sws
             twp = twp / sws
          else
             z0p = 0._wp
             twp = 0._wp
          end if
          !------------------------
          ! fix for rounding errors
          !------------------------
          if(abs(slp-1._wp) < 1.e-5_wp) slp = 1._wp
          if(abs(slp      ) < 1.e-5_wp) slp = 0._wp
          if(abs(fip-1._wp) < 1.e-5_wp) fip = 1._wp
          if(abs(fip      ) < 1.e-5_wp) fip = 0._wp
          if(abs(snf      ) < 1.e-5_wp) snf = 0._wp
          snf = snf * 0.01_wp
          if (lbg) then
            spt% ps_bg   = psp  ! surface pressure
            spt% gp_bg   = gsp  ! surface geopotential
            spt% sl_bg   = slp  ! sea-land mask
            spt% fi_bg   = fip  ! sea ice fraction
            spt% z0_bg   = z0p  ! roughness length
            spt% ssd_bg  = sdp  ! SSO standard deviation
            spt% tl_bg   = tlp  ! temp. at lowest level
            spt% ts_bg   = tsp  ! temp. at surface
            spt% tw_bg   = twp  ! sea/water temperature
            spt% dtdz_bg = dtz  ! lapse rate 850..700 hPa
            spt% hs_bg   = hsp  ! snow height
            if (u10_use_mlevel) then
              spt% u10bg = ulp  ! 10m wind from proxy (lowest model level)
              spt% v10bg = vlp  ! 10m wind from proxy (lowest model level)
            else
              spt% u10bg = u10  ! 10m wind
              spt% v10bg = v10  ! 10m wind
            end if
            spt% mdlsfc  = mdlsfc_frac (spt)
            spt% pz_bg   = gacc * psp / (R * tlp) ! factor d ps / d z
!         else
!           psp = spt% ps_bg
!           gsp = spt% gp_bg
!           slp = spt% sl_bg
!           tlp = spt% tl_bg
!           tsp = spt% ts_bg
!           dtz = spt% dtdz_bg
!           hsp = spt% hs_bg
!           u10 = spt% u10bg
!           v10 = spt% v10bg
          endif
        endif
        !-----------------------------------
        ! observation operator specific part
        !-----------------------------------

        !-----------------------------------------------------
        ! special vertical interpolation for radiances
        ! 1: linear
        ! 2: spline
        ! 3: follow Rochon et al. (2007) QJRMS (for radiances)
        ! use modtype = -TOVS for choice 3
        !-----------------------------------------------------
        modtype = spt% hd% modtype
        if (modtype == TOVS) then
          if (nwv_rad == 4) modtype = -TOVS
        endif

        select case (modtype)
        !============================
        ! COSMO observation operators
        !============================
        case (COSMO, SOIL)
          !----------------------------------------------
          ! no interpolation, just copy model levels,
          ! currently not used; not integrated into 3dvar
          !----------------------------------------------
        !======================================
        ! GPSGB ground based GNSS (slant delay)
        !======================================
        case (GPSGB)
          !------------------------------------------
          ! no interpolation, just copy model levels,
          ! convert  t,q,gp  to  t,rh,gp
          !------------------------------------------
          n = spt% i% n
          i = 0
          do ic = 1, size(spt% imcol)
            ii = spt% imcol(ic)% imc(1)
            call std_col2xi (col(l)% col(ii)% t,    &! t        [K]
                             col(l)% col(ii)% q,    &! q        [kg/kg]
                             col(l)% col(ii)% geo,  &! geo      [m2/s2]
                         exp(col(l)% col(ii)% p),   &! p        [Pa]
                             int (i+1:i+3*ke)   ,   &! t,gh,geo [K][1][m2/s2]
                            .true.                  )! for gh
            i = i + 3*ke
          end do
        !======
        ! GPSRO
        !======
        case (GPSRO)
          !------------------------------------------
          ! no interpolation, just copy model levels,
          ! convert  t,q,gp  to  t,gh,z
          !------------------------------------------
          n = spt% i% n
          i = 0
          do ic = 1, size(spt% imcol)
            ii = spt% imcol(ic)% imc(1)
            call occ_col2xi (col(l)% col(ii)% t,        &! t       (K)
                             col(l)% col(ii)% q,        &! q       (kg/kg)
                             col(l)% col(ii)% s% geosp, &! geop    (m2/s2)
                         exp(col(l)% col(ii)% p),       &! p       (Pa)
                             int (i+1:i+2*ke+1) ,       &! t,gh,z  (K)(1)(gpm)
                            .true.                      )! for gh
            i = i + 2*ke+1
          end do
        !=======================
        ! TEMP, TOVS, AMV, AIREP
        !=======================
        case (TEMP, TOVS, AMV, AIREP, SATEM, WLIDAR)
          !-------------------------------
          ! set interpolation coefficients
          ! horizontal nn, not temporal
          !-------------------------------
          if ( int_nn ) then
            ml = maxloc(spt% col% h% w)
            spt% col% h% w        = 0._wp
            spt% col% h% w(ml(1)) = 1._wp
          endif

          lin_t = vint_lin_t
          if (spt% hd% modtype == TOVS) lin_t = vint_lin_tov

          !-------------------------------
          ! horizontal interpolation first
          !-------------------------------
          if ( .not. int_vh) then
            cp  = 0._wp
            ctv = 0._wp
            crh = 0._wp
            dtv = 0._wp
            drh = 0._wp
            !--------------------------------------------------------------
            !  loop over neighbouring grid points (horizontal and temporal)
            !--------------------------------------------------------------
            do it = 1, nt
            do ic = 1, size(spt% col% h% imc,1)
              ii = spt% col% h% imc(ic,it)
              if (ii==0) exit
              wt = spt% w_time
              if (it==1) wt = 1._wp - wt
              w1 = wt * spt% col% h% w(ic)  ! weight for this grid point/time
              c => col(l)% col(ii)
              d => d2    % col(ii)
              !----------------------------
              !  horizontal interpolation
              !----------------------------
              do i = 1, col(l)% ke
                cp(i)    = cp(i)  + w1 * c% p (i)
                if (associated(c% tv)) &
                  ctv(i) = ctv(i) + w1 * c% tv(i)
                if (associated(c% rh)) &
                  crh(i) = crh(i) + w1 * c% rh(i)
              enddo
            end do
            end do
            dt = lapse_cl * R/gacc * ctv(ke)
            call init_splinex (cp, ctv, dtv, ypn = dt   )
            call init_splinex (cp, crh, drh, ypn = 0._wp)
            drh = 0._wp                 !  linear interpolation rel. humidity
            if (lin_t) dtv = 0._wp      !  linear interpolation temperature
          endif

          !---------------------------
          ! set up spline coefficients
          !---------------------------
          do it = 1, nt
          do ic = 1, size(spt% col% h% imc,1)
            i = spt% col% h% imc(ic,it)
            if (i==0) exit
            if (.not. lini (i)) then
              lini (i) = .true.
              c => col(l)% col(i)
              d => d2    % col(i)
              if ( int_vh) then
                if (associated(c% tv   )) then
                  dt = lapse_cl * R/gacc * c% tv(ke)
                  call init_splinex (c% p, c% tv,  d% tv,  ypn = dt)
                  if (vint_lin_t) d% tv = 0._wp  ! linear interp. temperature
                endif
                if (associated(c% rh  )) then
!                 call init_splinex (c% p, c% rh,  d% rh, ypn = 0._wp)
                  d% rh=0 ! linear interpolation humidity
                end if
              end if
              if (associated(c% q   )) then
!               call init_splinex (c% p, c% q,   d% q,  ypn = 0._wp)
                d% q=0    ! linear interpolation humidity
              endif
              if (associated(c% u   )) then
                call init_splinex (c% p, c% u,   d% u,  ypn = 0._wp)
                if (vint_lin_uv) d% u = 0._wp  ! linear interp. wind
              endif
              if (associated(c% v   )) then
                call init_splinex (c% p, c% v,   d% v,  ypn = 0._wp)
                if (vint_lin_uv) d% v = 0._wp  ! linear interp. wind
              endif
              if (associated(c% geoh)) then
                !---------------------------------------------
                ! for GME/HRM get geop.height from half levels
                !---------------------------------------------
                if (associated(c% tv   )) then
                  call init_splinex (c% ph(2:),c% geoh(2:),d% geoh(2:), &
                                     ypn = -R * c% tv(ke))
                else if (associated(c% t   )) then
                  call init_splinex (c% ph(2:),c% geoh(2:),d% geoh(2:), &
                                     ypn = -R * c% t(ke))
                else
                  call init_splinex (c% ph(2:),c% geoh(2:),d% geoh(2:))
                endif
                d% geoh(1) = 0.
                if (vint_lin_z) d% geoh = 0._wp
              else if (associated(c% geo)) then
                !-------------------------------------------
                ! for COSMO get geop.height from half levels
                !-------------------------------------------
                if (associated(c% tv   )) then
                  call init_splinex (c% p ,c% geo ,d% geo, &
                                     ypn = -R * c% tv(ke))
                else if (associated(c% t   )) then
                  call init_splinex (c% p ,c% geo ,d% geo, &
                                     ypn = -R * c% t(ke))
                else
                  call init_splinex (c% p ,c% geo ,d% geo)
                endif
                if (vint_lin_z) d% geo = 0._wp
              endif
            endif
          end do
          end do
          !-------------------
          ! loop over elements
          !-------------------
          z = huge(z)
          do k = 1, spt% i% n
            i = spt% i% i + k
            k2 = (k+1)/2
            !----------------------------------------------
            ! call interpolation routine for each new level
            !----------------------------------------------
            if (obs% o(l)% lev(i) /= z) then
              z = obs% o(l)% lev(i)
            endif

            if (l_tt) l_atm_lev = k2>nl_st .and. k2<=nl+nl_st

            !----------------------------
            ! set interpolated quantities
            !----------------------------
            ! select case (obs% o(l)% t_int(i))
            ! case (OBS_TV)
            if (obs% o(l)% t_int(i) == OBS_TV .or. (l_tt.and.l_atm_lev)) then
              if ( .not. int_vh) then
              !------------------------------------------------------------------
              !  vertical interpolation (spline)
              !  of horizontally and temporally interpolated columns: temperature
              !------------------------------------------------------------------
                call spline  (      &
                     cp,            & ! <-- Argument grid
                     ctv,           & ! <-- Gridded function
                     dtv,           & ! <-- 2nd derivative of spline
                     z,             & ! <-- Interpolation point
                     tp)              ! --> Interpolated function value
                if (spt% hd% modtype == TOVS) then
                  !-----------------------------------------------------------
                  ! Radiances: constant extrapolation at top and below surface
                  !-----------------------------------------------------------
                  if     (z > cp (col(l)% ke)) then
                    tp = ctv (col(l)% ke)
                  elseif (z < cp (1))          then
                    tp = ctv (1)
                  endif
                endif
                int_ = tp
              else
              !------------------------------------------------
              ! vertical  ->  horizontal/temporal interpolation
              !------------------------------------------------
                int_ = 0._wp
                if (lin_t .and. .not. vint_lin_t) then
                  allocate (dtv(col(l)% ke))
                  dtv = 0._wp
                  do it = 1, nt
                  do ic = 1, size(spt% col% h% imc,1)
                    ii = spt% col% h% imc(ic,it)
                    if (ii==0) exit
                    wt = spt% w_time
                    if (it==1) wt = 1._wp - wt
                    c => col(l)% col(ii)
                    d => d2    % col(ii)
                    call spline  (c% p,  & ! <-- Argument grid
                                  c% tv, & ! <-- Gridded function
                                  dtv,   & ! <-- 2nd derivative of spline
                                  z,     & ! <-- Interpolation point
                                  tp)      ! --> Interpolated function value
                    if (spt% hd% modtype == TOVS) then
                      !--------------------------------------------------------
                      ! Radiances: constant extrapolation at top and below surf
                      !--------------------------------------------------------
                      if     (z > c% p (col(l)% ke)) then
                        tp = c% tv (col(l)% ke)
                      elseif (z < c% p (1))          then
                        tp = c% tv (1)
                      endif
                    endif
                    int_ = int_ + wt * spt% col% h% w(ic) * tp
                  end do
                  end do
                  deallocate (dtv)
                else
                  do it = 1, nt
                  do ic = 1, size(spt% col% h% imc,1)
                    ii = spt% col% h% imc(ic,it)
                    if (ii==0) exit
                    wt = spt% w_time
                    if (it==1) wt = 1._wp - wt
                    c => col(l)% col(ii)
                    d => d2    % col(ii)
                    call spline  (c% p,  & ! <-- Argument grid
                                  c% tv, & ! <-- Gridded function
                                  d% tv, & ! <-- 2nd derivative of spline
                                  z,     & ! <-- Interpolation point
                                  tp)      ! --> Interpolated function value
                    if (spt% hd% modtype == TOVS .or. vint_lin_t) then
                      !--------------------------------------------------------
                      ! Radiances: constant extrapolation at top and below surf
                      !--------------------------------------------------------
                      if     (z > c% p (col(l)% ke)) then
                        tp = c% tv (col(l)% ke)
                      elseif (z < c% p (1))          then
                        tp = c% tv (1)
                      endif
                    endif
                    int_ = int_ + wt * spt% col% h% w(ic) * tp
                  end do
                  end do
                endif
              end if

              if (obs% o(l)% t_int(i) == OBS_TV) int(k) = int_
              if (l_tt .and. l_atm_lev)          tt_tv(k2) = int_
            end if
            if (obs% o(l)% t_int(i) == OBS_RH .or. (l_tt.and.l_atm_lev)) then
              if ( .not. int_vh) then
                !---------------------------------------------------
                ! vertical interpolation (spline) of horizontally
                ! and temporally interpolated columns: rel. humidity
                !---------------------------------------------------
                call spline  (      &
                     cp,            &      ! <-- Argument grid
                     crh,           &      ! <-- Gridded function
                     drh,           &      ! <-- 2nd derivative of spline
                     z,             &      ! <-- Interpolation point
                     qp)                   ! --> Interpolated function value
                if (z < cp(1 )) qp=crh(1)  ! extrapolation lower boundary
                if (z > cp(ke)) qp=crh(ke) ! extrapolation upper boundary

                int_ = gh_rh (qp, .true., 200._wp, exp(z))
              else
              !----------------------------------------
              ! vertical  ->  horizontal interpolation
              !----------------------------------------

                int_ = 0._wp
                do it = 1, nt
                do ic = 1, size(spt% col% h% imc,1)
                  ii = spt% col% h% imc(ic,it)
                  if (ii==0) exit
                  wt = spt% w_time
                  if (it==1) wt = 1._wp - wt
                  c => col(l)% col(ii)
                  d => d2    % col(ii)
                  call spline  (c% p,  & ! <-- Argument grid
                                c% rh, & ! <-- Gridded function
                                d% rh, & ! <-- 2nd derivative of spline
                                z,     & ! <-- Interpolation point
                                qp)      ! --> Interpolated function value
                  if (z < c% p(1 )) qp=c% rh(1)  ! extrapolation
                  if (z > c% p(ke)) qp=c% rh(ke) !
                  int_ = int_ + wt * spt% col% h% w(ic) * qp
                end do
                end do
                int_ = gh_rh (int_, .true., 200._wp, exp(z))
              end if

              if (obs% o(l)% t_int(i) == OBS_RH) int(k)    = int_
              if (l_tt .and. l_atm_lev)          tt_rh(k2) = int_

            end if

            if (obs% o(l)% t_int(i) == OBS_Q) then

              if ( .not. int_vh) then
                call finish('interpolation',                         &
                            'OBS_Q .and..not. int_vh not implemented')
              else
              !----------------------------------------
              ! vertical  ->  horizontal interpolation
              !----------------------------------------

                int(k) = 0._wp
                do it = 1, nt
                do ic = 1, size(spt% col% h% imc,1)
                  ii = spt% col% h% imc(ic,it)
                  if (ii==0) exit
                  wt = spt% w_time
                  if (it==1) wt = 1._wp - wt
                  c => col(l)% col(ii)
                  d => d2    % col(ii)
                  call spline  (c% p,  & ! <-- Argument grid
                                c% q , & ! <-- Gridded function
                                d% q , & ! <-- 2nd derivative of spline
                                z,     & ! <-- Interpolation point
                                qp)      ! --> Interpolated function value
                  if(z<c% p(1 )) qp=c% q(1)  ! extrapolation
                  if(z>c% p(ke)) qp=c% q(ke) !
                  int(k) = int(k) + wt * spt% col% h% w(ic) * qp
                end do
                end do
              end if

            elseif (obs% o(l)% t_int(i) == OBS_U) then

              int(k) = 0._wp
              do it = 1, nt
              do ic = 1, size(spt% col% h% imc,1)
                ii = spt% col% h% imc(ic,it)
                if (ii==0) exit
                wt = spt% w_time
                if (it==1) wt = 1._wp - wt
                c => col(l)% col(ii)
                d => d2    % col(ii)
                call spline  (c% p, & ! <-- Argument grid
                              c% u, & ! <-- Gridded function
                              d% u, & ! <-- 2nd derivative of spline
                              z,    & ! <-- Interpolation point
                              up)     ! --> Interpolated function value
                int(k) = int(k) + wt * spt% col% h% w(ic) * up
              end do
              end do

            elseif (obs% o(l)% t_int(i) == OBS_V) then

              int(k) = 0._wp
              do it = 1, nt
              do ic = 1, size(spt% col% h% imc,1)
                ii = spt% col% h% imc(ic,it)
                if (ii==0) exit
                wt = spt% w_time
                if (it==1) wt = 1._wp - wt
                c => col(l)% col(ii)
                d => d2    % col(ii)
                call spline  (c% p, & ! <-- Argument grid
                              c% v, & ! <-- Gridded function
                              d% v, & ! <-- 2nd derivative of spline
                              z,    & ! <-- Interpolation point
                              vp)     ! --> Interpolated function value
                int(k) = int(k) + wt * spt% col% h% w(ic) * vp
              end do
              end do

            elseif (obs% o(l)% t_int(i) == OBS_H) then

              int(k) = 0._wp
              do it = 1, nt
              do ic = 1, size(spt% col% h% imc,1)
                ii = spt% col% h% imc(ic,it)
                if (ii==0) exit
                wt = spt% w_time
                if (it==1) wt = 1._wp - wt
                c => col(l)% col(ii)
                d => d2    % col(ii)
                if (associated(c% geoh)) then
                  call spline  (c% ph  (2:), & ! <-- Argument grid
                                c% geoh(2:), & ! <-- Gridded function
                                d% geoh(2:), & ! <-- 2nd derivative of spline
                                z,           & ! <-- Interpolation point
                                gp)            ! --> Interpolated function value
                else
                  call spline  (c% p  , & ! <-- Argument grid
                                c% geo, & ! <-- Gridded function
                                d% geo, & ! <-- 2nd derivative of spline
                                z,      & ! <-- Interpolation point
                                gp)       ! --> Interpolated function value
                endif
                int(k) = int(k) + wt * spt% col% h% w(ic) * gp
              end do
              end do
              int(k) = int(k) / gacc
            elseif (obs% o(l)% t_int(i) == OBS_HS) then
!             int(k) = (gsp + (log(psp) - obs% o(l)% lev(i)) * R * tlp) / gacc
              int(k) = gsp/gacc + (psp - exp(obs% o(l)% lev(i))) / spt%pz_bg
            elseif (obs% o(l)% t_int(i) == OBS_DUM) then
              if (associated (obs% o(l)% bgi)) then
                int(k) = obs% o(l)% bgi(i)
              else
                int(k) = 0._wp
              endif
            elseif (obs% o(l)% t_int(i) == OBS_TV .or. obs% o(l)% t_int(i) == OBS_RH) then
              ! Done above
            else
              write (0,*) 'interpolate(1): unknown interpolated quantity',&
                          obs% o(l)% t_int(i)
              write (0,*) dace% pe,'interpolate(1): modtype =',spt% hd% modtype
              write (0,*) 'interpolate(1): obstype =',spt% hd% obstype
              call finish ('interpolate(1)','unknown interpolated quantity')
            end if
          end do
        !=====================================================
        ! TOVS, vertical interpolation following Rochon et al.
        !=====================================================
        case (-TOVS)

          !---------------------------------------
          ! set interpolation coefficients
          ! horizontal, not temporal interpolation
          !---------------------------------------
          if ( int_nn ) then
            ml = maxloc(spt% col% h% w)
            spt% col% h% w        = 0._wp
            spt% col% h% w(ml(1)) = 1._wp
          endif

          ! Get neighboring columns and fields available in all columns
          call get_col_imc(spt, col(l), afp, nid, iatm_mask=iatm_tovs_all, int_rad_hum=int_rad_hum)
          if (l_debug) write(usd,*) dpref,'int fidall',afp%fidall(1:afp%nall)
          if (.not.l_vint_rt) then
            i_p = -1
            do i = 1, afp%nall
              if (afp%fidall(i) == COL_P) then
                i_p = i
                exit
              end if
            end do
          end if

          !-------------------------------
          ! horizontal interpolation first
          !-------------------------------
          if ( .not. int_vh) then
            cp  = 0._wp
            rt_in = 0._wp
!             ctv = 0._wp
!             crh = 0._wp
            !-------------------------------------
            !  loop over neighbouring grid points
            !-------------------------------------
            do it = 1, nt
            do ic = 1, size(spt% col% h% imc,1)
              ii = spt% col% h% imc(ic,it)
              if (ii==0) exit
              wt = spt% w_time
              if (it==1) wt = 1._wp - wt
              c => col(l)% col(ii)
              !--------------------------------------
              ! horizontal and temporal interpolation
              !--------------------------------------
              do fi = 1, afp%nall
                fidx = afp%idxall(ic,fi)
                do i=1, size( afp%colp(ic,fidx)%pt)
                  rt_in(fi,i) = rt_in(fi,i) + wt * spt% col% h% w(ic) * afp%colp(ic,fidx)%pt(i)
                end do
              end do
            end do
            end do
            !-----------------------------------
            ! vertical interpolation of T and RH
            !-----------------------------------
            if (l_vint_rt) then
              if (int_rad_hum == 3) call convert2rttov_units(afp, x=rt_in, mode=1)
              call interpolate_rttov(rt_out(:,1:nl,0), y_in=rt_in, afp=afp,&
                                           lnp_out=lnp(nl_st+1:nl_st+nl))
              if (int_rad_hum == 3) call convert2rttov_units(afp, x=rt_out(:,1:nl,0), mode=-1)
            else
              rt_out(:,1:nl,0) = rt_in
            end if
          else
            !---------------------------------------
            ! vertical interpolation at grid columns
            !---------------------------------------
            rti_ = nl
            if (l_cloud)   rti_ = rti_ + RTI_CLOUD
            if (l_vint_rt) rti_ = rti_ + RTI_VINT_RT
            rt_out(:,1:nl,0) = 0._wp
            do it = 1, nt
            do ic = 1, size(spt% col% h% imc,1)
              i = spt% col% h% imc(ic,it)
              if (i==0) exit
              wt = spt% w_time
              if (it==1) wt = 1._wp - wt
              if (l_debug) write(usd,*) dpref,'int rti',ic,i,rti_,rti(i),rti_ /= rti(i)
              if (rti_ /= rti(i)) then ! check whether we did this work already
                rti (i) = rti_
                c => col(l)% col(i)
                if (l_vint_rt) then
                  ! interpolate on RTTOV levels
                  if (int_rad_hum == 3) then
                    do fi = 1, afp%nall
                      fidx = afp%idxall(ic,fi)
                      k = size( afp%colp(ic,fidx)%pt)
                      rt_in(fi,1:k) = afp%colp(ic,fidx)%pt(1:k)
                    end do
                    call convert2rttov_units(afp, x=rt_in(1:afp%nall,:), mode=1)
                    call interpolate_rttov(rt_out(1:afp%nall,1:nl,i), y_in=rt_in(1:afp%nall,:),afp=afp,&
                         lnp_out=lnp(nl_st+1:nl_st+nl))
                    call convert2rttov_units(afp, x=rt_out(1:afp%nall,1:nl,i),mode=-1)
                  else
                    call interpolate_rttov(rt_out(:,1:nl,i), afp=afp, ic=ic,&
                                           lnp_out=lnp(nl_st+1:nl_st+nl))
                  end if
                else
                  ! copy model level
                  do fi = 1, afp%nall
                    fidx = afp%idxall(ic,fi)
                    rt_out(fi,1:nl,i) = afp%colp(ic,fidx)%pt(:)
                  end do
                end if
              end if
              !--------------------------------------
              ! horizontal and temporal interpolation
              !--------------------------------------
              rt_out(1:afp%nall,1:nl,0) = rt_out(1:afp%nall,1:nl,0) + &
                                          wt * spt% col% h% w(ic) * rt_out(1:afp%nall,1:nl,i)
              if (l_debug .and. spt%col%h%w(ic) > 0._wp) then
                do k = 1, nl
                  write(usd,*) dpref,'int out',dace%pe,ic,i,k,rt_out(:,k,0),spt% col% h% w(ic)
                end do
              end if
            end do
            end do
          endif
          if (l_vint_rt) then
            p(1:nl) = exp(lnp(nl_st+1:nl_st+nl))
          else
            p(1:nl) = exp(rt_out(i_p,1:nl,0))
          end if

          if (l_debug) then
            write(usd,*) dpref,'int i_p',i_p
            do k = 1, nl
              write(usd,*) dpref,'int p',k,p(k)
            end do
          end if

          select case(int_rad_hum)
          case(1,3)
            ! Convert (Tv,Q) to (Tv,RH)
            dum_rt(1:nl,1) = rt_out(2,1:nl,0)                       ! Q
            dum_rt(1:nl,2) = t_tv_q(rt_out (1,1:nl,0), dum_rt(1:nl,1)) ! T
            rt_out(2,1:nl,0) = rh_q  (dum_rt(1:nl,1), dum_rt(1:nl,2), p(1:nl))
          case(2)
            ! Convert (Tv,Q) to (Tv,RH)
            dum_rt(1:nl,1) = rt_out(3,1:nl,0) ! Q
            dum_rt(1:nl,2) = t_tv_q(rt_out (1,1:nl,0), dum_rt(1:nl,1)) ! T
            ! WARNING: Here index 3 (COL_Q) is filled with RH calc. from Q
            rt_out(3,1:nl,0) = rh_q  (dum_rt(1:nl,1), dum_rt(1:nl,2), p(1:nl))
            where(log(p(1:nl)) > pblend(1))
              rt_wgt(1:nl) = 1._wp
            elsewhere(log(p(1:nl)) > pblend(2))
              rt_wgt(1:nl) = cos((log(p(1:nl)) -pblend(1))/(pblend(2)-pblend(1))*pi*0.5_wp)**2
            elsewhere
              rt_wgt(1:nl) = 0._wp
            end where
            if (l_debug) then
              do k = 1, nl
                write(usd,*) dpref,'pblend',k,p(k),log(p(k)),log(p(k)) > pblend,rt_wgt(k),rt_out(2:3,k,0)
              end do
            end if
            ! Blend Q-interpolated profile with RH-interpolated profile
            rt_out(2,1:nl,0) = rt_out(2,1:nl,0) * rt_wgt(1:nl) + rt_out(3,1:nl,0) * (1._wp-rt_wgt(1:nl))
          end select

          ! Move p such that indices match with indices in t_tovs%av
          if (nl_st > 0) p(nl_st+1:nl_st+nl) = p(1:nl)

          !------------------------------
          ! loop over single observations
          !------------------------------
          z  = huge(z)
          do k = 1, spt% i% n
            i = spt% i% i + k
            k2 = (k+1)/2
            !----------------------------------------------
            ! variable transform for each new level
            !----------------------------------------------
            if (obs% o(l)% lev(i) /= z) then
              z = obs% o(l)% lev(i)
            endif
            l_atm_lev = k2>nl_st .and. k2<=nl+nl_st
            if (mod(k,2) == 1 .and. l_atm_lev) then
              rh_ = gh_rh (rt_out (2,k2-nl_st,0), .true., 200._wp, exp(z))
              if (l_tt) then
                tt_tv(k2) = rt_out  (1,k2-nl_st,0)
                tt_rh(k2) = rh_
              end if
            end if
            !----------------------------
            ! set interpolated quantities
            !----------------------------
            select case (obs% o(l)% t_int(i))
            case (OBS_TV)
              int(k) = rt_out  (1,k2-nl_st,0)
            case (OBS_RH)
              int(k) = rh_
            case (OBS_HS)
              int(k) = gsp/gacc + (psp - exp(obs% o(l)% lev(i))) / spt%pz_bg
            case (OBS_DUM)
              if (associated (obs% o(l)% bgi)) then
                int(k) = obs% o(l)% bgi(i)
              else
                int(k) = 0._wp
              endif
            case default
              write (0,*) 'interpolate(2): unknown interpolated quantity',obs% o(l)% t_int(i)
              write (0,*) dace% pe,'interpolate(2): modtype =',spt% hd% modtype
              write (0,*) 'interpolate(2): obstype =',spt% hd% obstype
              call finish ('interpolate(2)','unknown interpolated quantity')
            end select
          end do
          if (l_tovs) then
            i0 = nl_st + 1
            i1 = nl_st + nl  ! == ttovs%nlev
            do i = 1, ttovs%nav
              if (i == ttovs%i_t .or. i == ttovs%i_q) cycle  ! Done below
              do fi = 1, afp%nall
                if (ttovs%av_cont(i) == afp%fidall(fi)) then
                  if (l_debug) then
                    write(usd,*) dpref,'int fid',ttovs%av_cont(i),i,fi,spt%hd%id
                    do k = 1, nl
                      write(usd,*) dpref,'int store av',k,k+nl_st,rt_out(fi,k,0)
                    end do
                  end if
                  if (.not.l_vint_rt .and. ttovs%av_cont(i) == COL_P) then
                    ttovs%i_p = i
                    ttovs%av(i0:i1,i) = p(i0:i1)
                  else
                    ttovs%av(i0:i1,i) = rt_out(fi,1:nl,0)
                    if (iand(afp%fidall(fi), COL_CLD_TOVS) /= 0) then
                      if (l_shift_cld) ttovs%av(1:ttovs%nlev-1,i) = ttovs%av(2:ttovs%nlev,i)
                      ttovs%av(i0:i1,i) = max(ttovs%av(i0:i1,i), 0._tpp)
                    end if
                  end if
                end if
              end do
            end do

          end if

        case (SYNOP, SCATT)
          !--------------------------------------------
          ! SYNOP: for the time take lowest model level
          !--------------------------------------------
!          call w_intpol (atm% grid, spt% col% dlon, spt% col% dlat, &
!                         n, ix, iy, id, w)
!          if (n<=0) then
!            write(0,*) 'interpolate: w_intpol returned n<=0',&
!              spt%statid,spt%col%dlon,spt%col%dlat
!            spt% qcf = ior (spt% qcf, QC_INTPOL)
!          endif

          up = 0._wp; vp = 0._wp; t2p = 0._wp; rlp = 0._wp; tvp = 0._wp
          rp = 0._wp; t2l = 0._wp; r2l = 0._wp
          qp = -999._wp
          do it = 1, nt
          do ic = 1, size(spt% col% h% imc,1)
            ii = spt% col% h% imc(ic,it)
            if (ii==0) exit
            wt = spt% w_time
            if (it==1) wt = 1._wp - wt      ! first time weight
            w1 = wt * spt% col% h% w(ic)    ! weight for this grid point/time
            c => col(l)% col(ii)
            if (u10_use_mlevel) then
              up =  up  + w1 * c%    u(ke)
              vp =  vp  + w1 * c%    v(ke)
!             up =  up  + w1 * c% s% ull
!             vp =  vp  + w1 * c% s% vll
            else
              up =  up  + w1 * c% s% u10m
              vp =  vp  + w1 * c% s% v10m
            end if
            t2p  = t2p  + w1 * c% s% t2m
            rlp  = rlp  + w1 * c% s% rhll
            t2l  = t2l  + w1 * c% s% t2mland
            r2l  = r2l  + w1 * c% s% rh2mland
            !-------------------------------------------------------
            ! take rh2m first guess from 2m diagnostics or 10m level
            !-------------------------------------------------------
            if (c% s% rh2m < 0._wp) then
              !-------------------
              ! old: always 10m RH
              !-------------------
              rp =  rp  + w1 * c% s% rhll
            else
              select case (version)
              case (:0)
                !-------------------
                ! old: always 10m RH
                !-------------------
                rp =  rp  + w1 * c% s% rhll
              case (1,3:)
                if ((spt% hd% obstype  == OT_SYNOP .and. &
                     spt% hd% buf_type == 0      ) .or.  &
                     spt% hd% obstype  == OT_DRIBU       ) then
                  !-------------------
                  ! SYNOP land from 2m
                  !-------------------
                  rp =  rp  + w1 * c% s% rh2m
                else
                  !---------------
                  ! other from 10m
                  !---------------
                  rp =  rp  + w1 * c% s% rhll
                endif
              case (2)
                !-------------
                ! always 2m RH
                !-------------
                rp =  rp  + w1 * c% s% rh2m
              end select
            endif
          end do
          end do

          !-------------------------------------------------------------
          ! For SYNOP land, use land-tile average 2m values if available
          !-------------------------------------------------------------
          if (version >= 5 .and. t2l > 0._wp .and. r2l > 0._wp) then
             if (spt% hd% obstype  == OT_SYNOP .and. &
                 spt% hd% buf_type == 0              ) then
                t2p = t2l
                rp  = r2l
             end if
          end if

!         if (spt% qcf == 0) then
            !-------------------------------------------------------
            ! take model surface pressure if no pressure is reported
            !                          and surface pressure is valid
            !-------------------------------------------------------
            do k = 1, spt% o% n
              i = spt% o% i + k
              if (obs% o(l)% olev(i) == use_ps_model .and. &
                  spt%       ps_bg   >  0._wp        .and. &
                  spt%       ps_bg   <  150000._wp         ) then
                  obs% o(l)% olev(i) = spt% ps_bg
              endif
              if (obs% o(l)% body(i)% plev == use_ps_model .and. &
                  spt%       ps_bg   >  0._wp              .and. &
                  spt%       ps_bg   <  150000._wp               ) then
                  obs% o(l)% body(i)% plev = spt% ps_bg
              endif
            end do
            !----------------------------------------
            ! interpolation below lowest model level?
            !----------------------------------------
!           hs = -999._wp
            w2 =    0._wp
            if (version == 4 .or. version >= 6) then
               do k = 1, spt% o% n
                  i = spt% o% i + k
                  select case (obs% o(l)% varno (i))
                  case (VN_T2M, VN_TD2M, VN_RH, VN_RH2M)
                     if (obs% o(l)% body(i)% lev_typ == VN_HOSAG) then
                        hs = min (max (obs% o(l)% olev(i), 2._wp), 10._wp)
                        w2 = (hs - 2._wp) / (10._wp - 2._wp)
                        exit
                     end if
                  end select
               end do
            end if
            !--------------------
            ! interpolated values
            !--------------------
            do k = 1, spt% i% n
              i = spt% i% i + k
              select case (obs% o(l)% t_int(i))
              case (OBS_U )
                int(k) = up
                obs% o(l)% lev (i) = log (spt% ps_bg)
              case (OBS_V )
                int(k) = vp
                obs% o(l)% lev (i) = log (spt% ps_bg)
              case (OBS_FF)
                int(k) = sqrt (up**2 + vp**2)
                obs% o(l)% lev (i) = log (spt% ps_bg)
              case (OBS_H )
!               int(k) = (gsp + (log(psp) - obs% o(l)% lev(i)) * R * tvp) /gacc
!write(0,*) "OBS_H: gsp, dtdzp =", gsp, dtdzp
!write(0,*) "    psp, exp(lev) =", psp, exp(obs%o(l)% lev(i)), obs%o(l)% lev(i)
                if (dtdzp == 0._wp) then
                  int(k) = gsp/gacc + (psp - exp(obs% o(l)% lev(i))) / spt%pz_bg
                else
!write(0,*) "args:", exp(obs%o(l)%lev(i)), psp, tlp, gsp
                  int(k) = zbs (exp(obs%o(l)%lev(i)), psp, tlp, gsp) / gacc
!write(0,*) "zbs :", int(k), gacc
                endif
              case (OBS_RH)
                r2c = rp + w2*(rlp-rp)
                int(k) = gh_rh (r2c, fg=.true.)
                obs% o(l)% lev (i) = log (spt% ps_bg)
!write(0,*) "### hs,w,rlp,rp =", real ([hs,w2,rlp,rp])
              case (OBS_T,OBS_TV)
                !---------------------------------------
                ! T2M: apply simple height correction
                ! from model orography to station height
                !---------------------------------------
!write(0,*) "### hs,w,tlp,t2p=", real ([hs,w2,tlp,t2p])
                if (t2p > 0._wp) then
                  if (spt% dtdz_bg >  -98._wp .and. &
                      spt% dtdz_bg /= invalid       ) then
                    t2c = t2p - (gsp/gacc - spt% z) * spt% dtdz_bg
                  else
                    t2c = t2p + (gsp/gacc - spt% z) * dtdzp
                  endif
                else
                  t2c = t2p
                endif
                t2c = t2c + w2*(tlp-t2p)   ! correct below lowest level
                r2c = rp  + w2*(rlp-rp)    ! use linear interpolation 2m/10m

                if (obs% o(l)% t_int(i) == OBS_TV) then
                  if (rp <= 0._wp) then
                     write(0,*) 'interpolate(SYNOP,OBS_TV): rh=', rp
                     call finish ('interpolate','interpolated rh <= 0')
                  end if
                  tvp    = tv_t_rh (t2c, r2c, spt% ps_bg)
                  int(k) = tvp
!write(0,*) "OBS_TV: ", spt% statid, "t2p,t2c,r2c,tvp=", real([t2p,t2c,r2c,tvp])
                else
                  int(k) = t2c
!write(0,*) "OBS_T : ", spt% statid, "t2p,t2c,r2c   = ", real([t2p,t2c,r2c])
                end if
                obs% o(l)% lev (i) = log (spt% ps_bg)
              case default
                write(0,*) 'interpolate: t_int=', obs% o(l)% t_int(i)
                call finish ('interpolate','unknown interpolated quantity')
              end select
            end do
!         else
!           int = 0._wp
!         endif
        case default
          if (spt% i% n == 0) return
          write (0,*) 'interpolate : unknown observation operator', modtype
          call finish('interpolate','unknown observation operator')
        end select

        !------------------------------------
        ! set quantities in observation space
        !------------------------------------
        if (present (psi)) then
          psi% s(l)% x(spt% i% i +         1:  &
                       spt% i% i + spt% i% n ) &
                      = int (:)
        endif

        if (l_tovs) then
          if (l_tt) then
            ! We got Tv, gh but we want to store T,q in t_tovs.
            i0 = nl_st + 1
            i1 = nl_st + nl  ! == ttovs%nlev
#ifdef __NEC__
            call tq_tvgh_vec (nl, &
#else
            call tq_tvgh     (    &
#endif
                              tt_t(i0:i1), tt_q(i0:i1), tt_tv(i0:i1), tt_rh(i0:i1), &
                              p(i0:i1), i_fail(i0:i1))
            do i = i0,i1
              if (i_fail(i) < 0) write(0,*) 'tq_tvgh failed (interpolate):',spt%hd%id,i,p(i),&
                   tt_t(i), tt_q(i), tt_tv(i), tt_rh(i)
              if (l_debug) write(usd,*) dpref,'int ttovs%av(T,Q)',i,p(i),tt_t(i),tt_q(i),&
                   '<-Tv,GH:',tt_tv(i),tt_rh(i)
            end do
            if (ttovs%i_t > 0) ttovs%av(i0:i1,ttovs%i_t) = tt_t(i0:i1)
            if (ttovs%i_q > 0) ttovs%av(i0:i1,ttovs%i_q) = tt_q(i0:i1)
          end if
          if (snf >= 0._wp) ttovs%snf  = snf
          ttovs%init = ior(ttovs%init, TTOVS_AV)
          call store(obs%o(l), spt, ttovs)
        end if
        if (spt%hd%modtype == TOVS) call destruct(ttovs)


        deallocate ( int)
        !-----------------------------------------------
        !  deallocate horizontally interpolated columns
        !-----------------------------------------------
        if (.not. int_vh) then
          deallocate ( cp )
          deallocate ( ctv)
          deallocate ( crh)
          deallocate ( dtv)
          deallocate ( drh)
        end if
      end do

      !---------------------
      ! deallocate temporaries
      !---------------------
      deallocate (lini)
      call dealloc_cols (d2)
      if (nwv_rad ==4) then
        if (allocated(rt_out)) deallocate(rt_out)
        if (allocated(rti))    deallocate(rti)
        if (allocated(rt_in))  deallocate(rt_in)
        if (int_rad_hum > 0) then
          deallocate(dum_rt)
          if (int_rad_hum == 2) deallocate(rt_wgt)
        end if
      endif
      if (allocated(p)) deallocate(p)
    end do
  end subroutine interpolate
!---------------------------------------------------------------------------
  !--------------------------
  ! subroutine get_col_fields
  !--------------------------
  !> Get fields available in a given column
  !> The t_col structure doesn't provide any information about the fields
  !> available for a given column. There is a bit field "ids" in t_cols
  !> which lists all fields which may be found in any column, but in a
  !> certain column there are usually less fields allocated. Therefore, all
  !> fields defined in t_col are checked and pointers are set to the fields
  !> which are actually available.
  !>
  !> Fields are identified using the variables COL_XX defined in mo_t_col.f90,
  !> e.g. COL_T=1, COL_P=32, COL_QVDIA=16777216...
  !>
  !> \param [in]  c       pointer to column of type t_col
  !> \param [out] cfp     pointer array to fields in column
  !> \param [out] nfield  number of fields found in column, the first
  !>                      "nfield" elements of cfp(i)%pt point to the
  !>                      fields in the column.
  !> \param [in] int_rad_hum  set pointer in a specific order, required
  !>                          for satellite data
  !--------------------------------------------------------------------------
  subroutine get_col_fields(c, cfp, nfield, int_rad_hum, iatm)
    type(t_col) ,pointer              :: c
    type (col_field_pointer), pointer :: cfp(:)
    integer, intent(out)              :: nfield
    integer, intent(in)               :: int_rad_hum  ! see module header
    integer(i8), intent(in), optional :: iatm


    integer :: n, ncfp

    ! Re-initialize array of field IDs and nullify pointer
    cfp(:)%fid = -1_i8
    ncfp = size(cfp)
    do n=1, ncfp
       nullify(cfp(n)%pt)
    end do

    ! For RTTOV backward compatibility
    ! The order of some array elements is predefined
    n = 4  ! next free index for fields other than t, rh, q
    if (int_rad_hum /= 2) n = 3
    nfield = 0
    call check_field    (c%tv, COL_TV, nfield)  ! COL_TV => index 1
    if (int_rad_hum == 0 .or. int_rad_hum == 2) then
       nfield = 1
       call check_field (c%rh, COL_RH, nfield)  ! COL_RH => index 2
    end if
    if (int_rad_hum == 1 .or. int_rad_hum == 3) then
       nfield = 1
       call check_field (c%q,  COL_Q,  nfield)  ! COL_Q => index 2
    else if (int_rad_hum == 2) then
       nfield = 2
       call check_field (c%q,  COL_Q,  nfield)  ! COL_Q => index 3
    end if

    nfield = n - 1
    call check_field(c%t,      COL_T,      nfield)
    !call check_field(c%q,      COL_Q,      nfield)
    !call check_field(c%rh,     COL_RH,     nfield)
    call check_field(c%u ,     COL_UV,     nfield)
    call check_field(c%v ,     COL_UV,     nfield)
    call check_field(c%p ,     COL_P ,     nfield)
    call check_field(c%ph,     COL_PH,     nfield)
    call check_field(c%geo,    COL_GEO,    nfield)
    call check_field(c%geoh,   COL_GEOH,   nfield)
    !call check_field(c%tv,     COL_TV,     nfield)
    call check_field(c%x,      COL_X,      nfield)
    call check_field(c%qcl,    COL_QCL,    nfield)
    call check_field(c%qci,    COL_QCI,    nfield)
    call check_field(c%qr,     COL_QR,     nfield)
    call check_field(c%qs,     COL_QS,     nfield)
    call check_field(c%qg,     COL_QG,     nfield)
    !call check_field(c%infl,   COL_INFLAT, nfield)
    !call check_field(c%tv2,    COL_TV2,    nfield)
    !call check_field(c%irpc,   COL_IRPC,   nfield)
    !call check_field(c%pp,     COL_PP,     nfield)
    !call check_field(c%w,      COL_W,      nfield)
    call check_field(c%clc,    COL_CLC,    nfield)
    call check_field(c%qv_dia, COL_QVDIA,  nfield)
    call check_field(c%qc_dia, COL_QCDIA,  nfield)
    call check_field(c%qi_dia, COL_QIDIA,  nfield)
    call check_field(c%o3,     COL_OZONE,  nfield)
    call check_field(c%reff_qc,COL_REFF_QC,nfield)
    call check_field(c%reff_qi,COL_REFF_QI,nfield)
    !call check_field(c%, COL_, nfield)
    !call check_field(c%, COL_, nfield)
    !call check_field(c%, COL_, nfield)

  contains

    subroutine check_field(pt, id, nfield)
      real (wp),        pointer       :: pt(:)
      integer(i8),      intent(in)    :: id
      integer,          intent(inout) :: nfield

      if (present(iatm)) then
        if (iand(iatm, id) == 0) return
      end if
      if (associated(pt)) then
         nfield = nfield + 1
         if (nfield > ncfp) call finish('get_col_fields','nfield >ncfp')
         cfp(nfield)%fid = id
         cfp(nfield)%pt  => pt
      end if

    end subroutine check_field

  end subroutine get_col_fields
!---------------------------------------------------------------------------
  !-----------------------
  ! subroutine get_col_imc
  !-----------------------
  !> Get fields in neighboring columns
  !>
  !> Find all fields available in the neighboring columns, set poiters to
  !> the fields and get the subset of fields which is available in all
  !> columns. The latter is required for horizontal interpolation.
  !>
  !> Atmospheric fields can be specified in two ways:
  !> 1) Fields are identified using the variables COL_XX defined in mo_t_col.f90,
  !>    e.g. COL_T=1, COL_P=32, COL_QVDIA=16777216...
  !>    Fied IDS are stored in:
  !>    afp%colp(:,:)%fid, afp%fidall, iatm
  !> 2) Several fields can be specified in a bit field where each bit represents
  !>    a given field.
  !>    Bit fields are stored in:
  !>    afp%ids, afp%idall, atmid
  !>
  !> \param [in]    spt   spot
  !> \param [in]    cols  model columns container (t_cols)
  !> \param [inout] afp   structure with pointers to fields in
  !>                      neighboring columns
  !> \param [in]    nid   maximum number of fields per column, taken from
  !>                      t_cols/ids, required to allocate arrays
  !> \param [in]    iatm  bit field with required atmospheric fields, an error
  !>                      occours if one of these fields is missing, optional
  !> \param [in]    atmid array of field IDs with required atmospheric fields,
  !>                      an error occours if one of these fields is missing,
  !>                      optional
  !> \param [in] int_rad_hum  set pointer in a specific order, required
  !>                          for satellite data, optional
  !>
  !> The parameter "iatm" and "atmid" provide the same kind of information
  !> but in different formats. Only one of them should be used.
  !--------------------------------------------------------------------------
  subroutine get_col_imc(spt, cols, afp, nid, iatm, iatm_mask, atmid, int_rad_hum)
    type(t_spot),   pointer                   :: spt          ! pointer to meta data element
    type(t_cols),               intent(in)    :: cols         ! column data
    type (col_imc), target,     intent(inout) :: afp          ! structure with field info
    integer,                    intent(in)    :: nid          ! max. number of fields per column
    integer(i8),      optional, intent(in)    :: iatm         ! bit field with required atm. fields
    integer(i8),      optional, intent(in)    :: iatm_mask    ! get only fields in iatm_mask
    integer(i8),      optional, intent(in)    :: atmid(:)     ! required field IDs
    integer,          optional, intent(in)    :: int_rad_hum  ! specify order of fields

    type(t_col),             pointer      :: c                    ! pointer to column
    type(col_field_pointer), pointer      :: cfp(:)    ! pointer to fields in column nc
    integer(i8),             allocatable  :: fidall(:)
    integer                               :: ic, i, j, k, nc, it, nt
    integer(i8)                           :: iatm_
    logical                               :: TV_miss, RH_miss, Q_miss

    ! Check for time interpolation
    nt = 2; if (spt% w_time==0) nt = 1

    ! initialize arrays/pointers in afp
    if (allocated(afp%nfield)) deallocate(afp%nfield)
    if (allocated(afp%ids   )) deallocate(afp%ids   )
    if (allocated(afp%fidall)) deallocate(afp%fidall)
    if (allocated(afp%colp  )) deallocate(afp%colp  )
    if (allocated(afp%idxall)) deallocate(afp%idxall)
    afp%ncol  = 0
    afp%nmax  = 0
    afp%nall  = 0
    afp%idall = -1_i8   ! set all bits in integer

    afp%ncol = nt*size(spt% col% h% imc,1)
    allocate(afp%nfield(1:afp%ncol)      )
    allocate(afp%ids   (1:afp%ncol)      )
    allocate(afp%colp  (1:afp%ncol,1:nid))
    afp%nfield = -1
    afp%ids = 0_i8

    ! Set pointer to fields in all neighboring columns
    nc = 0

    do it = 1, nt
      do ic = 1, size(spt% col% h% imc,1)
        i = spt% col% h% imc(ic,it)
        if (i==0) exit
        c => cols%col(i)
        nc = nc + 1
        cfp => afp%colp(nc,:)
        call get_col_fields(c, cfp, afp%nfield(nc), int_rad_hum, iatm=iatm_mask)
        ! Create bit field with bits set for all fields in current column
        do j=1, afp%nfield(nc)
          afp%ids(nc) = ior(afp%ids(nc), afp%colp(nc,j)%fid)
        end do
        ! Create bit field with fields available in all columns
        afp%idall = iand(afp%idall, afp%ids(nc))
      end do
    end do
    afp%ncol = nc                  ! available number of neighboring columns
    afp%nmax = maxval(afp%nfield)  ! max. number of fields in any column
    afp%nall = popcnt(afp%idall)   ! number of fields available in all columns

    allocate(afp%fidall(1:afp%nall))
    allocate(afp%idxall(1:afp%ncol,1:afp%nall))
    afp%fidall = 0_i8
    afp%idxall = -1

    ! Convert bit field to single field IDs
    j = 1
    do i=0, digits(afp%idall)
       if (btest(afp%idall,i)) then
          afp%fidall(j) = ibset(afp%fidall(j),i)
          j = j + 1
       end if
    end do

    ! For RTTOV backward compatibility
    ! The order of some array elements is predefined
    if (present(int_rad_hum)) then

       allocate(fidall(1:afp%nall))
       fidall = afp%fidall

       ! Check if fields are available
       TV_miss = .false.
       RH_miss = .false.
       Q_miss  = .false.

       j = 4  ! next free index for fields other than t, rh, q
       if (int_rad_hum /= 2) j = 3
       afp%fidall(1) = COL_TV                           ! COL_TV => index 1
       TV_miss = .true.
       if (int_rad_hum == 0 .or. int_rad_hum == 2) then
          afp%fidall(2) = COL_RH                        ! COL_RH => index 2
          RH_miss = .true.
       end if
       if (int_rad_hum == 1 .or. int_rad_hum == 3) then
          afp%fidall(2) = COL_Q                         ! COL_Q => index 2
          Q_miss  = .true.
       else if (int_rad_hum == 2) then
          afp%fidall(3) = COL_Q                         ! COL_Q => index 3
          Q_miss  = .true.
       end if

       do i=1, afp%nall
          select case (fidall(i))
          case (COL_TV)
             TV_miss = .false.
          case (COL_RH)
             RH_miss = .false.
          case (COL_Q)
             Q_miss = .false.
          case default
             afp%fidall(j) = fidall(i)
             j = j + 1
          end select
       end do
       if (TV_miss) &
            call finish('interpolate','RTTOV interpolation, field TV missing')
       if (RH_miss) &
            call finish('interpolate','RTTOV interpolation, field RH missing')
       if (Q_miss) &
            call finish('interpolate','RTTOV interpolation, field Q missing')

       deallocate(fidall)

    end if    ! if (present(int_rad_hum)) then

    ! Check if all required fields are really available
    if (present(iatm) .or. present(atmid)) then
       if (present(iatm)) then
          ! Check bit field
          iatm_ = iatm
       else
          ! Check array of field IDs => convert array to bit field
          iatm_ = 0_i8
          do j=1, size(atmid)
             iatm_ = ior(iatm_, atmid(j))
          end do
       end if
       ! Check if required fields ar in "idall"
       if (iatm_ /= iand(iatm, afp%idall)) then
         write(0,*) 'iatm',iatm_,afp%idall,iand(iatm, afp%idall)
         iatm_ = iatm_ - iand(iatm, afp%idall)
         do j = 0, 63
           if (btest(iatm_, j)) write(0,*) 'Missing bit',j,2_i8**j
         end do
         call finish('interpolate','required fields missing')
       end if
    end if

    ! Get array indices for fields in "fidall"
    do j=1, afp%ncol
       do i=1, afp%nfield(j)
          do k=1, afp%nall
             if (iand(afp%colp(j,i)%fid,afp%fidall(k)) > 0_i8) then
               afp%idxall(j,k) = i
               if (.not.associated(afp%colp(j,i)%pt)) then
                 write(0,*) 'spot',spt%hd%id
                 write(0,*) 'fid',k,afp%fidall(k)
                 write(0,*) 'j,i,shape',j,i,shape(afp%colp)
                 call finish('get_col_imc','field not associated')
               end if
             end if
          end do
       end do
    end do

    ! Test: Check indices in "idxall"
    if (ldeb(spt)) then
      do i = 1, afp%nall
        write(usd,*) dpref,'fidall:',i,afp%fidall(i)
        write(usd,*) dpref,'idx: i, afp%idxall(:,i) = ', i, afp%idxall(:,i)
        do j=1, afp%ncol
          write(usd,*) dpref,'test idx: j, i, afp%idxall(j,i) = ', j, i, afp%idxall(j,i)
          if (afp%fidall(i) /= afp%colp(j,afp%idxall(j,i))%fid) then
            write(usd,*) dpref,'Wrong field ID: fidall(i), afp%colp%fid = ', &
                 afp%fidall(i), afp%colp(j,afp%idxall(j,i))%fid
          end if
        end do
      end do
    end if

  end subroutine get_col_imc

!------------------------------------------------------------------------------

  !-----------------------------
  ! subroutine interpolate_rttov
  !-----------------------------
  !> Vertical interpolation of atmospheric profiles on RTTOV levels
  !>
  !> Atmospheric profiles taken from the 3D fields are interpolated on
  !> RTTOV levels using the RTTOV routine "rttov_layeravg".
  !>
  !> The atmospheric profiles  can be provided in two ways:
  !> 1) Array of profiles, y_in
  !>    As the profiles in y_in are not specified, an additional pressure
  !>    profile needs to be given. This could be "lnp_in" or the pointer
  !>    array afp, where the pressure profile is selected automatically.
  !> 2) Array of pointers to several profiles, afp
  !>    If afp is provided, the current column index ic needs also to
  !>    be  given. If afp is available, there should always be a pointer
  !>    to the pressure profile and lnp_in is not required.
  !>
  !> \param [out] y_out   output array of interpolated profiles
  !> \param [in]  y_in    input array of atmospheric profiles, optional
  !> \param [in]  lnp_in    atmospheric profile of log(pressure), optional
  !> \param [in]  afp     structure with pointers to fields in
  !>                      neighboring columns, optional
  !> \param [in]  ic      current column index, optional
  !--------------------------------------------------------------------------
  subroutine interpolate_rttov (y_out, y_in, lnp_in, lnp_out, afp, ic, ldeb)

  real(wp),       intent(out)                  :: y_out  (:,:)
  real(wp),       intent(in), optional         :: y_in   (:,:) ! Array of input profiles
  real(wp),       intent(in), optional, target :: lnp_in (:)   ! input ln(p) profile
  real(wp),       intent(in), optional, target :: lnp_out(:)   ! output ln(p) profile
  type (col_imc), intent(in), optional         :: afp          ! Pointer array and column index
  integer,        intent(in), optional         :: ic           ! column index
  logical,        intent(in), optional         :: ldeb

#if (_RTTOV_VERSION <= 0)
      call finish ('interpolate_rttov','RTTOV not configured')
#else
#include "rttov_layeravg.interface"

    integer               :: nl_in, nl_out         ! in/out level numbers
    integer,  allocatable :: kstart(:)             ! start index ..
    integer,  allocatable :: kend  (:)             ! end   index ..
    real(wp), allocatable :: pz (:,:)              ! interpolation coefficients
    integer               :: k, j, n_int, fi, idx
    real(wp), pointer     :: cp(:)                 ! pointer to pressure field
    real(wp), pointer     :: lnp_(:)               ! pointer to pressure field
    logical               :: ld

    if (present(ldeb)) then
      ld = ldeb
    else
      ld = .false.
    end if

    ! Get logarithm of model pressure profile for vertical interpolation
    ! Warning: The COL_P profiles are already ln(p) profiles!
    !          There is no need to convert p to ln(p)!
    nullify(cp)
    if (present(y_in)) then
       ! Interpolate profiles in "y_in"
       if (present(lnp_in)) then
          cp => lnp_in
       else if (present(afp)) then
          do fi=1, afp%nall
             if (afp%fidall(fi) == COL_P) then
                ! set pointer to pressure field (full levels)
                cp => afp%colp(1,afp%idxall(1,fi))%pt
                exit
             end if
          end do
       end if
    else if(present(afp) .and. present(ic)) then
       ! Interpolate selected profiles in "afp"
       do fi=1, afp%nall
          if (afp%fidall(fi) == COL_P) then
          ! set pointer to pressure field (full levels)
             cp => afp%colp(ic,afp%idxall(ic,fi))%pt
             exit
          end if
       end do
    else
       ! Invalid input options
       call finish('interpolate_rttov','invalid combination of parameters')
    end if

    if (.not. associated(cp)) then

       call finish('interpolate_rttov','pressure profile not available')
    end if

    if (present(lnp_out)) then
      lnp_ => lnp_out
    else
      lnp_ => lnp
    end if

    nl_in  = size(cp)
    nl_out = size(lnp_)

    allocate(pz(nl_in,nl_out), kstart(nl_out), kend(nl_out))

    ! Get vertical interpolation coefficients from RTTOV
    ! lnp - ln(p) on RTTOV levels
    ! cp  - ln(p) on model levels
    call rttov_layeravg (lnp_,       &! <-  levels of output domain
                         cp,         &! <-  levels of input  domain
                         nl_out,     &! <-  output domain size
                         nl_in,      &! <-  input  domain size
                         pz(:,:),    &!  -> interpolation coefficients
                         kstart,     &!  -> start index
                         kend        &!  -> end   index
#if (_RTTOV_VERSION >= 12)
                         ,vint_rttov &! <-  RTTOV interpolation mode
#endif
                         )

    if (present(y_in)) then
       ! Interpolate profiles in "y_in"
       n_int = size(y_in,1)
       if (size(y_out,1) < n_int) then
          write(0,*) 'size(y_in,1)=',n_int,'size(y_out,1)=',size(y_out,1)
          call finish('interpolate_rttov','output array too small')
       end if
       y_out = 0._wp
       do k = 1, nl_out
          do j = kstart(k), kend(k)
             y_out(1:n_int,k) =  y_out(1:n_int,k) + pz(j,k) * y_in(1:n_int,j)
          end do
       end do
    else if(present(afp) .and. present(ic)) then
       ! Interpolate selected profiles in "afp"
       if (size(y_out,1) < afp%nall) then
         write(0,*) 'afp%nall=',afp%nall,'size(y_out,1)=',size(y_out,1)
         call finish('interpolate_rttov','output array too small')
       end if
       y_out = 0._wp
       do fi=1, afp%nall
          idx = afp%idxall(ic,fi)
          do k = 1, nl_out
             do j = kstart(k), kend(k)
                y_out(fi,k) =  y_out(fi,k) + pz(j,k) * afp%colp(ic,idx)%pt(j)
             end do
          end do
       end do
    end if

#endif
  end subroutine interpolate_rttov
!------------------------------------------------------------------------------
  !> Convert (Tv,Q) to RTTOV units (T,ppmv). Useful for comparison with other tools
  !> like RadSim, that call RTTOV with the internal RTTOV interpolation.
  subroutine convert2rttov_units(afp, ic, x, mode)
    type (col_imc), intent(inout)           :: afp         ! Pointer array and column index
    integer,        intent(in),    optional :: ic
    real(wp),       intent(inout), optional :: x(:,:)
    integer,        intent(in),    optional :: mode

    integer  :: lmode
    integer  :: fi, k, idx, nl
    integer  :: i_t, i_q
    real(wp) :: t, tv, q, ppmv

    if (present(mode)) then
      lmode = mode
    else
      lmode = 1
    end if
    lmode = sign(1, lmode)

    do fi=1, afp%nall
      select case(afp%fidall(fi))
      case(COL_TV)
        i_t = fi
      case(COL_Q)
        i_q = fi
      end select
    end do

    if (present(x)) then
      do k = 1, size(x,2)
        if (lmode == 1) then
          tv = x(i_t,k)
          q  = x(i_q,k)
          t = t_tv_q(tv, q)
          ppmv = ppmv_dry_q(q)
          x(i_t,k) = t
          x(i_q,k) = ppmv
        else
          t    = x(i_t,k)
          ppmv = x(i_q,k)
          q  = q_ppmv_dry(ppmv)
          tv = tv_t_q(t,q)
          x(i_t,k) = tv
          x(i_q,k) = q
        end if
      end do
    else if (present(ic)) then
      idx = afp%idxall(ic,i_t)
      nl = size(afp%colp(ic,idx)%pt)
      do k = 1, nl
        if (lmode == 1) then
          tv = afp%colp(ic,i_t)%pt(k)
          q  = afp%colp(ic,i_q)%pt(k)
          t = t_tv_q(tv, q)
          ppmv = ppmv_dry_q(q)
          afp%colp(ic,i_t)%pt(k) = t
          afp%colp(ic,i_q)%pt(k) = ppmv
        else
          t    = afp%colp(ic,i_t)%pt(k)
          ppmv = afp%colp(ic,i_q)%pt(k)
          q  = q_ppmv_dry(ppmv)
          tv = tv_t_q(t,q)
          afp%colp(ic,i_t)%pt(k) = tv
          afp%colp(ic,i_q)%pt(k) = q
        end if
      end do
    else
      call finish('convert2rttov_units','either x or ic missing')
    end if

  end subroutine convert2rttov_units
!------------------------------------------------------------------------------

  ! Interpolate tracegas fields to TOVS spots
  subroutine interpolate_trg (obs)
    type(t_obs_set)  ,intent(inout)  :: obs  ! observation data type

    ! NOTE: tracegas is abbreviated as "trg"
    ! We assume that all tracegase-gribs are assumed to contain "mass mixing ratio",
    ! i.e. kg/(kg dry air). Tracegas is stored in ppmv_dry

    character(len=15), parameter   :: proc = 'interpolate_trg'

    ! grid and fields
    type(t_inventory), pointer     :: invt(:)  ! GRIB file inventory
    type(t_grid),      pointer     :: grid     ! grid information
    type(t_atm)                    :: atm_trg  ! background

    ! spot handling
    type(t_spot),      pointer     :: spt => null()
    type(t_tovs),      target      :: ttovs
    integer                        :: tovs_io
    real(kind=tpp),    allocatable :: av   (:,:)           ! profiles in t_tovs
    integer                        :: nb, ib               ! box loops
    integer                        :: is, isp              ! spot loops
    ! spot counting
    integer                        :: nsp                  ! #spots on local PE
    integer                        :: nsp_all              ! #spots on all PEs
    integer                        :: n                    ! #spots auxiliary
    integer                        :: n_clim,n_curr,n_fixed! for final summary
    integer                        :: isp_box(size(obs%o)) ! offset of boxes in spt_mask
    integer                        :: nsp_box(size(obs%o)) ! number of spots per box
    logical                        :: l_data               ! Whether there is any data to be processed
    logical                        :: l_interp             ! Whether interpolations are to be done

    ! vertical interpolation
    real(kind=wp),     allocatable :: trg_grd(:,:)         ! trg profile on orig. grid
    real(kind=wp),     allocatable :: trg_vint(:,:)        ! trg profile on target grid
    real(kind=wp),     allocatable :: lnp_grd(:)           ! ln(p) orig. grid
    real(kind=wp),     allocatable :: lnp_vint(:)          ! ln(p) target grid
    integer                        :: nl                   ! number of levels target grid

    ! horizontal interpolation
    type(t_mcols),     allocatable :: mc(:)                ! column(s) from orig grid
    type(t_cols),      allocatable :: cbgb(:)              ! column(s) from orig grid
    type(t_col),       pointer     :: c                    ! column(s) from orig grid
    type(t_hic),       allocatable :: h(:)                 ! interpolation indices/weights
    integer(i8)                    :: trg_col              ! COL_* for current trg
    integer                        :: ih, ii, ic

    ! background profile
    real(wp),          allocatable :: trg_bkg(:)           ! background profile
    integer                        :: nl_rt                ! number of RTTOV levels

    ! Species, that might be processed:
    integer,           parameter   :: ntrg     = 2         ! number of trgs currently implemented
    integer,           parameter   :: I_O3     = 1         ! O3 index
    integer,           parameter   :: I_CO2    = 2         ! CO2 index
    integer,           parameter   :: I_XXX    = 3         ! future trg index
    character(len=20)              :: trg_name             ! name of current trg
    integer,           pointer     :: trg_ind              ! index of current trg in t_tovs%av ...
                                                           ! .... pointer to t_tovs%i_*

    ! handling of use_x options in TOVS_OBS_CHAN namelists
    type(t_rad_set),   pointer     :: rs  => null()        ! dataset options
    integer,           parameter   :: USE_CURR   = 0       ! use current profile
    integer,           parameter   :: USE_CLIM   = 1       ! use climatological profile
    integer,           parameter   :: USE_FIXED  = 2       ! use Fixed value profile
    integer                        :: gopts                ! "global" use_* option (ior over all ..
                                                           ! ... TOVS_OBS_CHAN namelists)

    ! Spot mask, that contains info about options for each individual spot
    integer,           allocatable :: spt_mask(:)          ! Bitmask for each spot
    integer,           parameter   :: O3_CURR   =  0       ! Current O3 profile
    integer,           parameter   :: O3_CLIM   =  1       ! Climatological O3 profile
    integer,           parameter   :: O3_FIXED  =  2       ! Fixed value O3 profile
    integer,           parameter   :: O3_DONE   =  3       ! O3 interpolatione done
    integer,           parameter   :: CO2_CURR  =  4       ! Current CO2 profile
    integer,           parameter   :: CO2_CLIM  =  5       ! Climatological CO2 profile
    integer,           parameter   :: CO2_FIXED =  6       ! Fixed value CO2 profile
    integer,           parameter   :: CO2_DONE  =  7       ! CO2 interpolatione done
    ! integer,         parameter   :: XXX_CURR  =  8       ! Current XXX profile
    ! integer,         parameter   :: XXX_CLIM  =  9       ! Climatological XXX profile
    ! integer,         parameter   :: XXX_FIXED = 10       ! Fixed value XXX profile
    ! integer,         parameter   :: XXX_DONE  = 11       ! XXX interpolatione done
    integer,           parameter   :: ALL_CURR  = 2**O3_CURR  + 2**CO2_CURR  ! + 2**XXX_CURR
    integer,           parameter   :: ALL_CLIM  = 2**O3_CLIM  + 2**CO2_CLIM  ! + 2**XXX_CLIM
    integer,           parameter   :: ALL_CUCL  = ALL_CURR + ALL_CLIM
    integer,           parameter   :: ALL_FIXED = 2**O3_FIXED + 2**CO2_FIXED ! + 2**XXX_FIXED
    integer                        :: trg_curr             ! *_CURR  for current trg
    integer                        :: trg_clim             ! *_CLIM  for current trg
    integer                        :: trg_fixed            ! *_FIXED for current trg
    integer                        :: trg_done             ! *_DONE  for current trg

    ! Temporal interpolation (for climatological profiles)
    real(kind=wp),     allocatable :: mm_hint(:,:,:)       ! monthly values horizontally interpolated
    real(kind=wp),     allocatable :: trg_aux(:,:,:)       ! horiz. interp. for surrounding months
    real(kind=wp),     allocatable :: trg_hint(:,:)        ! horiz. interp. for actual time
    real(kind=wp)                  :: t, fac
    integer                        :: i0, i1, im

    ! Climate trend adjustment ("ct")
    integer                        :: n_trg_ct             ! #trgs with ct info
    logical                        :: l_ct                 ! whether current trg has
    type(t_trg_clim_trend),pointer :: ct
    character(len=10), parameter   :: d_padding = '0701000000' ! padding for short dates in namelist
    real(kind=wp)                  :: avg                  ! average of climatology
    real(kind=wp)                  :: v_ct_spec            ! current (average) value according ...
    real(kind=wp)                  :: v_ct_ppmv            !  ... to climate trend
    integer                        :: gas_id               ! gas_id of current trg (DACE)
    integer                        :: rt_gas_id            ! gas_id of current trg (RTTOV)
    logical,            allocatable:: l_msg(:)             ! required for message printing
    integer                        :: satid

    ! trg history
    type t_trg_history
      type(t_time)                 :: last_update   = invalid_time  ! date of last run
      type(t_time)                 :: last_current  = invalid_time  ! date of last current field
      character(len=14)            :: gas_name      = ''            ! gas name
      real(wp)                     :: adapt_time    = 14._wp        ! adaptation time
      real(wp)                     :: max_age       =  1._wp        ! maximum age to be considered as "current"
      real(wp)                     :: wgt_clim      =  0._wp        ! weight of climate (vs. current) in last run
    end type t_trg_history
    type(t_trg_history), target    :: trg_hist(ntrg)       ! history
    type(t_trg_history), pointer   :: th                   ! history (current trg)
    logical,             target    :: l_hist(ntrg)         ! whether history is available
    logical,             pointer   :: lh                   ! whether history is available (current trg)
    type(t_time)                   :: hist_date            ! staring time for history
    real(wp)                       :: wgt                  ! weight of "climate" according to history
    real(wp)                       :: age                  ! Age of current field
    logical                        :: hist_upd, hist_upd_  ! whether history shall be updated

    ! general variables
    logical                        :: lread(ntrg)          ! whether trg fields were read
    logical                        :: ldum
    character(len=300)             :: msg
    integer                        :: i, j, l, instr, ierr

    !------------------
    ! Read history file
    !------------------
    l_hist = .false.
    hist_date = time_c(trg_hist_inidate)
    if (ana_time >= hist_date) then
      if (trg_hist_file /= '') then
        if (ana_time > hist_date) call read_hist_file
        do i = 1, ntrg
          call set_trg(i)
          if (.not.lh) then
            th%gas_name     = trg_name
            th%last_update  = time_c('1900010100')
            th%last_current = time_c('1900010100')
            lh = .true.
          end if
        end do
      else
        call finish(proc, 'ana_time >= trg_hist_inidate requires trg_hist_file to be set')
      end if
    end if

    nb = size (obs% o)

    !---------------------------------------------------
    ! Determine spots, that require external information
    !---------------------------------------------------
    ! Check namelist options
    gopts = convert_opts(use_o3=glob_use_o3, use_co2=glob_use_co2)
    hist_upd = (iand(gopts, ALL_CURR) > 0) .and. any(l_hist)
    l_data = .false.
    nl_rt  = 0
    do i = 1, n_set
      rs => rad_set(i)
      do instr = 1, rs%n_instr
        if (rs%iopts(instr)%use_o3  > 0 .or. rs%iopts(instr)%use_co2 > 0) then
          l_data = .true.
          nl_rt = max(nl_rt, rs%iopts(instr)%rt_nlevs)
        end if
      end do
    end do

    if (l_data .or. hist_upd) then
      ! Count OT_RAD spots
      nsp = 0
      do ib = 1, nb
        if (dace%pe /= obs%o(ib)%pe) cycle
        do is = 1, obs%o(ib)%n_spot
          spt => obs%o(ib)%spot(is)
          if (spt% hd% obstype /= OT_RAD) cycle
          nsp = nsp + 1
        end do
      end do
      nsp_all = p_sum(nsp)
      l_data = (nsp_all > 0)
    end if

    if (l_data .or. hist_upd) then
      ! Which spot requires which information?
      allocate(spt_mask(nsp))
      spt_mask = 0
      nsp_box = 0
      isp_box = 0
      isp = 0
      do ib = 1, nb
        if (dace%pe /= obs%o(ib)%pe) cycle
        isp_box(ib) = isp
        do is = 1, obs%o(ib)%n_spot
          spt => obs%o(ib)%spot(is)
          if (spt% hd% obstype /= OT_RAD) cycle
          isp = isp + 1
          call load(obs%o(ib), spt, tovs=ttovs, tovs_io=0, rs=rs)
          spt_mask(isp) = convert_opts(rs=rs)
          if (spt_mask(isp) > 0) nsp_box(ib) = nsp_box(ib) + 1
        end do
        nsp_box(ib) = isp - isp_box(ib)
      end do
      nsp = count(spt_mask > 0)
      nsp_all = p_sum(nsp)
      l_data = (nsp_all > 0)
      l_interp = any(iand(spt_mask,ALL_CUCL) > 0)
      l_interp = p_or(l_interp)
    end if

    if (l_data .or. hist_upd) then
      if (dace%lpio) then
        write(*,*)
        write(*,'(1x,A,1x,I9,a)') 'Interpolate tracegases for TOVS:',nsp_all,' spots'
      end if

      if (l_data) then
        !--------------------------------------
        ! Climate trend/Fixed value preparation
        !--------------------------------------
        n_trg_ct = count(trg_clim_trend(:)%gas_name /= '' .and. &
                         trg_clim_trend(:)%date     /= '' .and. &
                         trg_clim_trend(:)%val      >=  0._wp)
        trg_clim_trend(1:n_trg_ct) = pack(trg_clim_trend, &
                                          mask=trg_clim_trend(:)%gas_name /= '' .and. &
                                               trg_clim_trend(:)%date     /= '' .and. &
                                               trg_clim_trend(:)%val      >=  0._wp)
        do j = 1, n_trg_ct
          ct => trg_clim_trend(j)
          call get_gas(ct%gas_name, name=trg_name)
          if (trg_name == '') call finish(proc,'unknown gas "'//trim(ct%gas_name)//'"')
          ct%gas_name = trg_name
          l = len_trim(ct%date)
          if (l < 14 .and. l <= 4) then
            ct%date(l+1:14) = d_padding(l-3:10)
          else
            call finish(proc,'invalid date in climate trend directive for '//&
                 trim(trg_name)//': "'//trim(ct%date)//'"')
          end if
        end do
        !------------
        ! Allocations
        !------------
        allocate(av(mx_nlev,mx_nav+1))
        if (l_interp) then
          allocate(mc(nb), cbgb(nb))
          ! Count max. required size of h array
          n = 0
          do i = 1, ntrg
            call set_trg(i)
            nsp = count(btest(spt_mask, trg_curr))
            n = max(n, nsp)
            nsp = count(btest(spt_mask, trg_clim))
            n = max(n, nsp)
          end do
          allocate(h(n))
        end if
      end if

      !---------------
      ! Current fields
      !---------------
      if (iand(gopts, ALL_CURR) > 0) then
        if (trg_file == '') then
          write(*,*) '*** trg_file not given (file for current tracegas fields)'
          do i = 1, ntrg
            call set_trg(i)
            nsp = count(btest(spt_mask, trg_curr)) ; nsp_all = p_sum(nsp)
            if (lh) then
              call update_hist(th)
              if (th%wgt_clim < 1._wp) then
                write(msg,*) 'trg_file not defined but required for '//trim(trg_name)//&
                     ' if climate weight smaller than 1. wgt_clim=',th%wgt_clim
                if (nsp_all > 0) then
                  call finish(proc, trim(msg))
                else
                  write(*,*) '*** '//trim(msg)
                end if
              end if
            else if (nsp_all > 0) then
              call finish(proc,'trg_file required!')
            end if
          end do
        else
          lread = .false.
          trg_file = path_file(input, trg_file)
          nsp = count(iand(spt_mask, ALL_CURR)>0) ; nsp_all = p_sum(nsp)
          if (dace%lpio) then
            write(*,*)
            write(*,*) 'Read current tracegas fields from trg_file='//trim(trg_file)
            write(*,*) '   (required by ',nsp_all,' spots)'
          end if
          ! Read file
          invt => null()
          call get_inventory (invt, trg_file)
          if (size(invt) > 0) then
            ! read grid, setup_cols
            call print_inventory (invt,  first=.true.)
            if (dace%lpio) write(*,*)
            allocate(grid)
            call read(grid, trg_file, invt=invt, geosp=.false., lsm=.false., &
                      nproc1=nproc1, nproc2=nproc2, comm=dace% comm)
            call print(grid,  verbose=.true.)
            ! read trg fields
            do i = 1, ntrg
              call set_trg(i)
              nsp = count(btest(spt_mask, trg_curr)) ; nsp_all = p_sum(nsp)
              hist_upd_ = lh .and. btest(gopts, trg_curr)
              if (nsp_all > 0 .or. hist_upd_) then
                ! Read fields
                if (dace%lpio) then
                  write(*,*)
                  write(*,*) 'Read current '//trim(trg_name)//' fields'
                  write(*,*) '   (required by ',nsp_all,' spots)'
                end if
                call construct (atm_trg ,grid)
                call read(atm_trg, trg_file, invt, fields=trim(trg_name), ierr=ierr) !runtype='forecast|init_ana'
                if (ierr /= 0) call finish(proc, 'failed to read '//&
                     trim(trg_name)//' (current) from trg_file')
                lread(i) = .true.
                if (dace%lpio) write(*,*) 'time: ',cyyyymmddhhmmss(atm_trg%time)
                !call print(atm_trg, verbose=.true.)
                call print(atm_trg%m, verbose=.true.)
                age = days(ana_time - atm_trg%time)
                if (dace%lpio) write(*,'(1x,a,F7.2,a)') 'Found '//trim(trg_name)//&
                     ' fields with age=',age,' days'
              end if

              ! history file
              if (hist_upd_) call update_hist(th, atm_trg%time)

              ! Interpolation
              if (nsp_all > 0) then
                ! Columns for horizontal interpolation
                call setup_cols(trg_curr)
                call get_cols(mc, atm_trg, cbgb, iatm=trg_col)
              end if
              if (nsp > 0) then
                ! vertical interpolation and store in t_tovs
                allocate(trg_grd(1,grid%nz), trg_vint(1,mx_nlev), &
                         lnp_grd(  grid%nz), lnp_vint(  mx_nlev))
                lnp_grd(:) = log(grid%akf(1:grid%nz))
                isp = 0
                ih  = 0
                do ib = 1, nb
                  if (dace%pe /= obs%o(ib)%pe) cycle
                  do is = 1, obs%o(ib)%n_spot
                    spt => obs%o(ib)%spot(is)
                    if (spt% hd% obstype /= OT_RAD) cycle
                    isp = isp + 1
                    if (btest(spt_mask(isp),trg_curr)) then
                      ih = ih + 1
                      call prep_tovs
                      ! horizontal interpolation
                      trg_grd = 0._wp
                      do ic = 1, size(h(ih)%imc,1)
                        ii = h(ih)% imc(ic,1)
                        if (ii==0) exit
                        c => cbgb(ib)% col(ii)
                        trg_grd(1,:) = trg_grd(1,:) + h(ih)% w(ic) * get_col(c, i)
                      end do
                      ! Vertical interpolation
                      if (ttovs%i_p > 0) then
                        lnp_vint(1:nl) = log(real(av(1:nl,ttovs%i_p), wp))
                      else
                        lnp_vint(1:nl) = lnp(1:nl)
                      end if
                      call interpolate_rttov(trg_vint(1:1,1:nl), y_in=trg_grd(1:1,:), &
                                             lnp_out=lnp_vint(1:nl), lnp_in=lnp_grd)
                      trg_vint(1,1:nl) = ppmv_dry_trg(trg_vint(1,1:nl), 0._wp, gas_id)
                      call add_av(trg_vint(1,1:nl), 1._wp) ! Store ttovs%av
                      spt_mask(isp) = ibset(spt_mask(isp), trg_done)
                      ! Fallbacks required?
                      ldum = .true.
                      if (lh) ldum = th%wgt_clim <= 0._wp
                      if (ldum) then
                        if (btest(spt_mask(isp), trg_curr) .and. btest(spt_mask(isp), trg_clim )) &
                             spt_mask(isp) = spt_mask(isp) - 2**trg_clim
                        if (btest(spt_mask(isp), trg_curr) .and. btest(spt_mask(isp), trg_fixed)) &
                             spt_mask(isp) = spt_mask(isp) - 2**trg_fixed
                      end if
                    end if ! btest(spt_mask, trg_curr)
                  end do ! spots
                end do ! boxes
                deallocate(trg_grd, trg_vint, lnp_grd, lnp_vint)
                call destruct(mc)
                call dealloc_cols(cbgb)
              end if ! nsp > 0
              call destruct(atm_trg)
            end do ! ntrg
            call destruct(grid)
          else
            if (dace%lpio) write(*,*) '*** Emtpy file!'
          end if ! size(invt) > 0
          deallocate(invt)
          ! Update history, if no current fields were available
          do i = 1, ntrg
            call set_trg(i)
            nsp = count(btest(spt_mask, trg_curr)) ; nsp_all = p_sum(nsp)
            if (.not.lread(i) .and. lh .and. btest(gopts, trg_curr)) then
              call update_hist(th)
              if (th%wgt_clim < 1._wp) then
                write(msg,*) 'trg_file emtpy/unreadable but required for '//trim(trg_name)//&
                     ' if climate weight smaller than 1. wgt_clim=',th%wgt_clim,&
                     ' trg_file="'//trim(trg_file)//'"'
                if (nsp_all > 0) then
                  call finish(proc, trim(msg))
                else
                  write(*,*) '*** '//trim(msg)
                end if
              else if (nsp_all > 0) then
                call finish(proc,'trg_file required but emtpy/unreadable: &
                     &trg_file="'//trim(trg_file)//'"!')
              end if
            end if
          end do
        end if ! trg_file /= ''
      end if ! iand(gopts, ALL_CURR) > 0

      !-------------------
      ! Climatology fields
      !-------------------
      n = count(iand(spt_mask, ALL_CLIM ) > 0)
      n = p_sum(n)
      if (n > 0) then
        if (trg_clim_file == '') then
          call finish  (proc, 'trg_clim_file required but not defined!')
        else
          ! Climate annual cycle
          trg_clim_file = path_file(input, trg_clim_file)
          if (dace%lpio) then
            write(*,*)
            write(*,*) 'Read climatological tracegas fields from trg_clim_file='//trim(trg_clim_file)
            write(*,*) '   (required by ',n,' spots)'
          end if
          ! Read file
          invt => null()
          call get_inventory (invt, trg_clim_file)
          if (size(invt) > 0) then
            ! read grid
            call print_inventory (invt,  first=.true.)
            if (dace%lpio) write(*,*)
            allocate  (grid)
            call read (grid, trg_clim_file, invt=invt, geosp=.false., lsm=.false., &
                       nproc1=nproc1, nproc2=nproc2, comm=dace% comm)
            call print     (grid,  verbose=.true.)
            ! read trg fields
            do i = 1, ntrg
              call set_trg(i)
              nsp = count(btest(spt_mask,trg_clim))
              nsp_all = p_sum(nsp)
              if (nsp_all > 0) then
                if (dace%lpio) then
                  write(*,*)
                  write(*,*) 'Read climatological '//trim(trg_name)//' fields'
                  write(*,*) '   (required by ',nsp_all,' spots)'
                end if
                call get_ct
                if (l_ct) avg = 0._wp
                ! Read fields for each month and do horizontal interpolation
                if (nsp > 0) allocate(mm_hint(12,nsp,grid%nz)) ! monthy means horizontally interpolated
                call setup_cols(trg_clim)
                do im = 1, 12
                  call construct (atm_trg ,grid)
                  call read(atm_trg, trg_clim_file, invt, fields = trim(trg_name), runtype='average', &
                       month=im, ierr=ierr)
                  if (ierr /= 0) then
                    write(msg,'("failed to read ",A," month ",I2," from ",A," ierr=",I4)') &
                         trim(trg_name), im, trim(trg_clim_file), ierr
                    call finish(proc, trim(msg))
                  end if
                  if (dace%lpio) write(*,*) 'time: ',cyyyymmddhhmmss(atm_trg%time)
                  call print(atm_trg%m, verbose=.true.)
                  call get_cols(mc, atm_trg, cbgb, iatm=trg_col)
                  if (l_ct) avg = avg + avg_field() / 12._wp
                  call destruct(atm_trg)
                  if (nsp > 0) then
                    ! horizontal interpolation
                    isp = 0
                    ih  = 0
                    do ib = 1, nb
                      if (dace%pe /= obs%o(ib)%pe) cycle
                      do is = 1, obs%o(ib)%n_spot
                        spt => obs%o(ib)%spot(is)
                        if (spt% hd% obstype /= OT_RAD) cycle
                        isp = isp + 1
                        if (btest(spt_mask(isp),trg_clim)) then
                          ih = ih + 1
                          mm_hint(im,ih,:) = 0._wp
                          do ic = 1, size(h(ih)%imc,1)
                            ii = h(ih)% imc(ic,1)
                            if (ii==0) exit
                            c => cbgb(ib)% col(ii)
                            mm_hint(im,ih,:) = mm_hint(im,ih,:) + h(ih)% w(ic) * get_col(c, i)
                          end do
                        end if
                      end do
                    end do
                  end if
                  call dealloc_cols(cbgb)
                end do
                call destruct(mc)

                if (l_ct) then
                  ! Climate trend adjustment
                  avg = ppmv_dry_trg(avg,0._wp,gas_id)   ! average from climatology
                  fac = v_ct_ppmv/avg
                  if (dace%lpio) then
                    write(*,'(1x,a)') 'Climate trend adjustment for '//trim(trg_name)//':'
                    write(*,'(1x,2(a,e13.6,"ppmv"))') 'old mean: ',avg,',  new mean:',v_ct_ppmv
                  end if
                end if

                if (nsp > 0) then
                  ! Temporal interpolation to current analysis date
                  t = frac_year(ana_time)
                  i0 = floor(t*12+0.5)
                  if (i0 <= 0) i0 = 12
                  i1 = i0 + 1
                  if (i1 >= 13) i1 = 1
                  allocate(trg_aux(2,nsp, grid%nz))
                  trg_aux(1,:,:) = mm_hint(i0,:,:)
                  trg_aux(2,:,:) = mm_hint(i1,:,:)
                  ! FFT
                  call fft(mm_hint(1:12,:,:), isign=-1, norm=1)
                  ! backwards FFT to current time of year
                  t = t * 12 - 0.5_wp  ! t=0 corresponds to mid-January in fft -> subtract half a month
                  allocate(trg_hint(nsp,grid%nz))
                  trg_hint(:,:) = mm_hint(1,:,:) / sqrt(12._wp)
                  do j = 1, 5
                    trg_hint(:,:) = trg_hint(:,:) + mm_hint(j*2  ,:,:) * cos(2*pi*j*t/12._wp) / sqrt(6._wp)
                    trg_hint(:,:) = trg_hint(:,:) - mm_hint(j*2+1,:,:) * sin(2*pi*j*t/12._wp) / sqrt(6._wp)
                  end do
                  trg_hint(:,:) = trg_hint(:,:) + mm_hint(6*2  ,:,:) * cos(2*pi*j*t/12._wp) / sqrt(12._wp)
                  ! Avoid excessive overfitting, in particular negative values
                  trg_hint = min(trg_hint, 1.5_wp*maxval(trg_aux,1))
                  trg_hint = max(trg_hint, 0.5_wp*minval(trg_aux,1))
                  deallocate(mm_hint, trg_aux)

                  if (l_ct) trg_hint(:,:) = trg_hint(:,:) * fac

                  ! vertical interpolation and store in t_tovs
                  allocate(                    trg_vint(1,mx_nlev))
                  allocate(lnp_grd(  grid%nz), lnp_vint(  mx_nlev))
                  lnp_grd(:) = log(grid%akf(1:grid%nz))
                  isp = 0
                  ih  = 0
                  do ib = 1, nb
                    if (dace%pe /= obs%o(ib)%pe) cycle
                    do is = 1, obs%o(ib)%n_spot
                      spt => obs%o(ib)%spot(is)
                      if (spt% hd% obstype /= OT_RAD) cycle
                      isp = isp + 1
                      if (btest(spt_mask(isp),trg_clim)) then
                        ih = ih + 1
                        call prep_tovs
                        ! Vertical interpolation
                        if (ttovs%i_p > 0) then
                          lnp_vint(1:nl) = log(real(av(1:nl,ttovs%i_p), wp))
                        else
                          lnp_vint(1:nl) = lnp(1:nl)
                        end if
                        call interpolate_rttov(trg_vint(1:1,1:nl), y_in=trg_hint(ih:ih,:), &
                             lnp_out=lnp_vint(1:nl), lnp_in=lnp_grd)
                        trg_vint(1,1:nl) = ppmv_dry_trg(trg_vint(1,1:nl), 0._wp, gas_id)
                        if (lh .and. btest(spt_mask(isp), trg_curr)) then
                          wgt = th%wgt_clim
                        else
                          wgt = 1._wp
                        end if
                        call add_av(trg_vint(1,1:nl), wgt)
                        spt_mask(isp) = ibset(spt_mask(isp), trg_done)
                        ! Fallback not required
                        if (btest(spt_mask(isp), trg_curr) .and. btest(spt_mask(isp), trg_fixed)) &
                             spt_mask(isp) = spt_mask(isp) - 2**trg_fixed
                      end if ! btest(spt_mask, trg_clim)
                    end do ! spots
                  end do ! boxes
                  deallocate(trg_hint, trg_vint, lnp_grd, lnp_vint)
                end if ! nsp > 0
              end if ! nsp_all > 0
            end do ! ntrg
            call destruct(grid)
          else
            call finish(proc,'emtpy trg_clim_file: "'//trim(trg_clim_file)//'"')
          end if ! size(invt) > 0
          deallocate(invt)
        end if ! trg_clim_file
      end if ! n > 0

      !---------------------
      ! Fixed value profiles
      !---------------------
      n = count(iand(spt_mask, ALL_FIXED) > 0)
      n = p_sum(n)
      if (n > 0) then
        allocate(trg_bkg(nl_rt), trg_vint(1,mx_nlev), &
                 lnp_grd(nl_rt), lnp_vint(  mx_nlev))
        if (dace%lpio) then ; allocate(l_msg(n_set * m_instr)) ; l_msg = .false. ; endif

!NEC$ nomove
        do i = 1, ntrg
          call set_trg(i)
          nsp = count(btest(spt_mask,trg_fixed)) ; nsp_all = p_sum(nsp)
          if (nsp_all > 0) then
            call get_ct
            if (dace%lpio) then
              write(*,*)
              write(*,'(1x,a,F9.4,a,e13.6,a)') 'Set '//trim(trg_name)//&
                   ' to constant RTTOV default profile.'
              write(*,*) '   (for ',nsp_all,' spots)'
            end if
            isp = 0
!NEC$ nomove
            do ib = 1, nb
              if (dace%pe /= obs%o(ib)%pe) cycle
!NEC$ nomove
              do is = 1, obs%o(ib)%n_spot
                spt => obs%o(ib)%spot(is)
                if (spt% hd% obstype /= OT_RAD) cycle
                isp = isp + 1
                if (btest(spt_mask(isp),trg_fixed)) then
                  call prep_tovs
                  call get_tovs_rs(ttovs, rs=rs)
                  instr = -1
                  do j = 1, rs%n_instr
                    select case(gas_id)
                    case(gas_id_o3)
                      if (btest(rs%iopts(j)%use_o3, USE_FIXED)) then
                        instr = j
                        exit
                      end if
                    case(gas_id_co2)
                      if (btest(rs%iopts(j)%use_co2, USE_FIXED)) then
                        instr = j
                        exit
                      end if
                    end select
                  end do
                  call rtifc_coef_prop(rs%rttov_indx(instr), igas=rt_gas_id, gas_bkg=trg_bkg, &
                       nlevs=nl_rt, preslev=lnp_grd, satid=satid, instr=ii) ! misuse lnp_grd as p
                  lnp_grd = lnp_grd * 100._wp  ! hPa -> Pa
                  if (l_ct) then
                    ! Climate trend adjustment
                    avg = avg_vert(trg_bkg, lnp_grd)     ! Here lnp_grd is actually p - NOT log(p)
                    fac = v_ct_ppmv/avg
                    if (dace%lpio) then
                      if (.not.l_msg(rs%rttov_indx(instr))) then
                        write(*,'(1x,a," (sat. ",I3.3,", instr. ",I3.2,"): ")') &
                             'Climate trend adjustment for '//trim(trg_name),satid,ii
                        write(*,'(1x,2(a,e13.6,"ppmv"))') 'old mean: ',avg,',  new mean:',v_ct_ppmv
                        l_msg(rs%rttov_indx(instr)) = .true.
                      end if
                    end if
                    trg_bkg = trg_bkg * fac
                  end if

                  if (nl >= nl_rt-1 .and. nl <= nl_rt) then
                    ! Here the uppermost level is omitted (if we run without the uppermost level)
                    ! This might cause slight differences compared to RTTOV without tracegas
                    ! (furthermore the background profile is stored in single precision in DACE -
                    !  not in double precision)
                    trg_vint(1,1:nl) = trg_bkg(1+(nl_rt-nl):nl_rt)
                  else
                    ! Vertical interpolation
                    lnp_grd(1:nl_rt) = log(lnp_grd(1:nl_rt))
                    if (ttovs%i_p > 0) then
                      lnp_vint(1:nl) = log(real(av(1:nl,ttovs%i_p), wp))
                    else
                      lnp_vint(1:nl) = lnp(1:nl)
                    end if
                    call interpolate_rttov(trg_vint(1:1,1:nl), &
                                           y_in=reshape(trg_bkg(1:nl_rt),(/1,nl_rt/)), &
                                           lnp_out=lnp_vint(1:nl), lnp_in=lnp_grd(1:nl_rt))
                  end if
                  ! No units conversion necessary, since RTTOV default profile is in ppmv_dry

                  ! Store ttovs%av
                  if (lh .and. btest(spt_mask(isp), trg_curr)) then
                    wgt = th%wgt_clim
                  else
                    wgt = 1._wp
                  end if
                  call add_av(trg_vint(1,1:nl), wgt)
                  spt_mask(isp) = ibset(spt_mask(isp), trg_done)
                end if ! btest(spt_mask, trg_clim)
              end do ! spots
            end do ! boxes
          end if ! nsp_all > 0
        end do ! ntrg
        deallocate(trg_bkg, trg_vint, lnp_grd, lnp_vint)
        if (dace%lpio) deallocate(l_msg)
      end if ! n>0

      !-----------------------------------------------------------------
      ! Print summary and check, whether all spots got the required info
      !-----------------------------------------------------------------
      if (dace%lpio) write(*,'(1x,A)') 'Summary of spots with interpolated trace gas:'
      do i = 1, ntrg
        call set_trg(i)
        ii = 2**trg_curr  + 2**trg_done
        n_curr  = count(iand(spt_mask(:),ii) == ii) ; n_curr  = p_sum(n_curr )
        ii = 2**trg_clim  + 2**trg_done
        n_clim  = count(iand(spt_mask(:),ii) == ii) ; n_clim  = p_sum(n_clim )
        ii = 2**trg_fixed + 2**trg_done
        n_fixed = count(iand(spt_mask(:),ii) == ii) ; n_fixed = p_sum(n_fixed)
        if (dace%lpio) write(*,'(3x,A10,3(2x,I9))') trim(trg_name), n_curr, n_clim, n_fixed

        ii = 2**trg_curr + 2**trg_clim + 2**trg_fixed
        n = count(iand(spt_mask(:),ii) > 0 .and. .not.btest(spt_mask(:), trg_done))
        if (n > 0) then
          ! failed spots
          isp = 0
          do ib = 1, nb
            if (dace%pe /= obs%o(ib)%pe) cycle
            do is = 1, obs%o(ib)%n_spot
              spt => obs%o(ib)%spot(is)
              if (spt% hd% obstype /= OT_RAD) cycle
              isp = isp + 1
              if (iand(spt_mask(isp),ii) > 0 .and. .not.btest(spt_mask(isp), trg_done)) &
                   write(0,*) 'failed '//trim(trg_name),isp,spt%hd%id,spt_mask(isp),    &
                   btest(spt_mask(isp),trg_curr),btest(spt_mask(isp),trg_clim),         &
                   btest(spt_mask(isp),trg_fixed)
            end do
          end do
          write(msg,'(a,I9,a)') 'Failed to interpolate '//trim(trg_name)//' for ',n,' spots'
          call finish(proc, trim(msg))
        end if
      end do

      ! Write result in files that might be easily plotted
      ! do i = 1, ntrg
      !   call set_trg(i)
      !   write(msg,'("trg_",a,"_",I3.3,".dat")') trim(trg_name), dace%pe
      !   open(101,file=trim(msg))
      !   ii = 2**trg_curr + 2**trg_clim + 2**trg_fixed
      !   isp = 0
      !   do ib = 1, nb
      !     if (dace%pe /= obs%o(ib)%pe) cycle
      !     do is = 1, obs%o(ib)%n_spot
      !       spt => obs%o(ib)%spot(is)
      !       if (spt% hd% obstype /= OT_RAD) cycle
      !       isp = isp + 1
      !       if (iand(spt_mask(isp),ii) > 0) then
      !         call prep_tovs
      !         nl = ttovs%nlev
      !         write(101,'(100(1x,e13.6))') spt%col%c%dlat,spt%col%c%dlon,av(1:nl,trg_ind)
      !       end if
      !     end do
      !   end do
      !   close(101)
      ! end do

      if (l_data) then
        deallocate(av)
        if (l_interp) deallocate(mc,cbgb,h)
      end if

    end if ! l_data .or. l_hist_upd

    !-------
    ! Finish
    !-------
    if (any(l_hist)) call write_hist_file

  contains

    ! Convert use_* opts from namelist into internal bits
    function convert_opts(rs, use_o3, use_co2) result(bmask)
      type(t_rad_set), intent(in), optional :: rs
      integer,         intent(in), optional :: use_o3
      integer,         intent(in), optional :: use_co2
      integer :: bmask
      bmask = 0
      if (present(rs)) then
        if (any(btest(rs%iopts(1:rs%n_instr)%use_o3 , USE_CURR ))) bmask = ibset(bmask, O3_CURR  )
        if (any(btest(rs%iopts(1:rs%n_instr)%use_o3 , USE_CLIM ))) bmask = ibset(bmask, O3_CLIM  )
        if (any(btest(rs%iopts(1:rs%n_instr)%use_o3 , USE_FIXED))) bmask = ibset(bmask, O3_FIXED )
        if (any(btest(rs%iopts(1:rs%n_instr)%use_co2, USE_CURR ))) bmask = ibset(bmask, CO2_CURR )
        if (any(btest(rs%iopts(1:rs%n_instr)%use_co2, USE_CLIM ))) bmask = ibset(bmask, CO2_CLIM )
        if (any(btest(rs%iopts(1:rs%n_instr)%use_co2, USE_FIXED))) bmask = ibset(bmask, CO2_FIXED)
      end if
      if (present(use_o3)) then
        if (btest(use_o3 , USE_CURR )) bmask = ibset(bmask, O3_CURR  )
        if (btest(use_o3 , USE_CLIM )) bmask = ibset(bmask, O3_CLIM  )
        if (btest(use_o3 , USE_FIXED)) bmask = ibset(bmask, O3_FIXED )
      end if
      if (present(use_co2)) then
        if (btest(use_co2, USE_CURR )) bmask = ibset(bmask, CO2_CURR )
        if (btest(use_co2, USE_CLIM )) bmask = ibset(bmask, CO2_CLIM )
        if (btest(use_co2, USE_FIXED)) bmask = ibset(bmask, CO2_FIXED)
      end if
    end function convert_opts

    subroutine setup_cols(itest)
      integer, intent(in) :: itest
      integer :: ib, is
      integer :: nsp, n
      ih = 0
      do ib = 1, nb
        mc(ib)% pe = obs%o(ib)%pe
        if (dace%pe /= obs%o(ib)%pe) cycle
        isp = isp_box(ib)
        nsp = nsp_box(ib)
        n = count(iand(spt_mask(isp+1:isp+nsp), itest) > 0)
        allocate(mc(ib)% idx(grid%lbg(1):grid%ubg(1),grid%lbg(2):grid%ubg(2),grid%lbg(4):grid%ubg(4),1))
        allocate(mc(ib)% c(4*n))
        mc(ib)% idx = 0
        mc(ib)% n = 0
        do is = 1, obs%o(ib)%n_spot
          spt => obs%o(ib)%spot(is)
          if (spt% hd% obstype /= OT_RAD) cycle
          isp = isp + 1
          if (btest(spt_mask(isp),itest)) then
            ih = ih + 1
            call idx_init (   &
                 spt% col% c, &! <-  column descriptor
                 h(ih),       &!  -> interpolation coefficients
                 mc(ib),      &! <-> model column descriptors
                 0_i8,        &! <-  fields required
                 0,           &! <-  tracers required
                 grid,        &! <-  model grid
                 1,           &! <-  time slot
                 0._wp        )! <-  time interpolation weight
          end if
        end do
      end do
    end subroutine setup_cols

    ! Get column for tracegases
    function get_col(c, i) result(p)
      real(kind=wp), pointer :: p(:)
      type(t_col), intent(in) :: c
      integer,     intent(in) :: i
      nullify(p)
      select case(i)
      case(I_O3)
        p => c%o3
      case(I_CO2)
        p => c%co2
      end select
    end function get_col

    ! Load ttovs and prepare ttovs%av
    subroutine prep_tovs
      integer :: j
      tovs_io = 0
      call load(obs%o(ib), spt, tovs=ttovs, tovs_io=tovs_io, av=av)
      nl = ttovs%nlev
      call set_trg(i,ltt=.true.)
      do j = 1, ttovs%nav
        if (ttovs%av_cont(j) == trg_col) then
          trg_ind = j
          exit
        end if
      end do
      if (trg_ind <= 0) then
        write(0,*) 'spot%hd%id',spt%hd%id
        write(0,*) 'spt_mask',isp,spt_mask(isp)
        write(0,*) 't_tovs%av_cont',ttovs%av_cont(1:ttovs%nav)
        call finish(proc,'no place for '//trim(trg_name)//' in t_tovs%av')
      end if
    end subroutine prep_tovs

    ! Set values for current trg
    subroutine set_trg(i,ltt)
      integer, intent(in)           :: i
      logical, intent(in), optional :: ltt
      logical :: lt
      lt = .false.
      if (present(ltt)) lt = ltt
      if (lt) then
        select case(i)
        case(I_O3)
          trg_ind  => ttovs%i_o3
        case(I_CO2)
          trg_ind  => ttovs%i_co2
        end select
      else
        select case(i)
        case(I_O3)
          trg_name  = 'o3'
          trg_col   = COL_OZONE
          trg_curr  = O3_CURR
          trg_clim  = O3_CLIM
          trg_fixed = O3_FIXED
          trg_done  = O3_DONE
          gas_id    = gas_id_o3
          rt_gas_id = rt_gas_id_o3
          lh        => l_hist(i)
          th        => trg_hist(i)
        case(I_CO2)
          trg_name  = 'co2'
          trg_col   = COL_CO2
          trg_curr  = CO2_CURR
          trg_clim  = CO2_CLIM
          trg_fixed = CO2_FIXED
          trg_done  = CO2_DONE
          gas_id    = gas_id_co2
          rt_gas_id = rt_gas_id_co2
          lh        => l_hist(i)
          th        => trg_hist(i)
        end select
      end if
    end subroutine set_trg

    ! Get info on climate trend adjustment
    subroutine get_ct
      real(kind=wp) :: t_
      type(t_time)  :: tc, td
      integer       :: j
      l_ct = .false.
      do j = 1, n_trg_ct
        ct => trg_clim_trend(j)
        if (ct%gas_name == trg_name) then
          l_ct = .true.
          exit
        end if
      end do
      if (l_ct) then
        tc = time_c(ct%date)
        td = ana_time - tc
        t_ = (td%days + td%secs/86400._wp)/365.2425_wp !
        v_ct_ppmv = ct%val + t_ * ct%trend        ! "nominal" average (from namelist)
        v_ct_spec = q_ppmv_dry_trg(v_ct_ppmv, 0._wp, gas_id)
      end if
    end subroutine get_ct

    ! RF Comment: we need this in order to calc. a "representative" value for
    ! tracegases with climate trend adjustment. It is not fully clear, which
    ! heights/layers are best representing the publicly available numbers on
    ! average trace gas concentrations. We follow rttov_scale_ref_gas_prof.f90 here.
    function avg_vert(trg, p) result(avg)
      real(kind=wp)        :: avg
      real(wp), intent(in) :: trg(:)
      real(wp), intent(in) :: p(:)

      real(wp), parameter :: avg_top = 10000._wp  !(100hPa)
      real(wp), parameter :: avg_bot = 100000._wp !(1000hPa)
      real(wp) :: wsum, w, wh, v, v0
      integer  :: i

      avg  = 0._wp
      wsum = 0._wp
      do i = 1, size(p)-1
        if (p(i+1) > avg_top .and. p(i) < avg_top) then
          w = p(i+1) - avg_top
          wh = w / (p(i+1)-p(i))
          v0 = (trg(i) * wh + trg(i+1) * (1._wp-wh))
          v = 0.5_wp * (v0 + trg(i+1))
        elseif (p(i+1) > avg_top .and. p(i) >= avg_top .and. &
                p(i+1) <= avg_bot .and. p(i) < avg_bot) then
          w = (p(i+1)-p(i))
          v = 0.5_wp * (trg(i) + trg(i+1))
        elseif (p(i+1) > avg_bot .and. p(i) < avg_bot) then
          w = avg_bot - p(i)
          wh = w / (p(i+1)-p(i))
          v0 = (trg(i) * wh + trg(i+1) * (1._wp-wh))
          v = 0.5_wp * (v0 + trg(i))
        else
          w = 0._wp
          v = 0._wp
        end if
        avg = avg + w * v
        wsum = wsum + w
      end do
      if (wsum > 0._wp) avg = avg/wsum

    end function avg_vert


    ! Average current trg field
    ! WARNING: correct only for regular grids so far
    function avg_field() result(avg)
      real(kind=wp)          :: avg
      real(kind=wp)          :: avg_h(grid%nz)
      real(kind=wp)          :: wsum
      real(kind=wp), pointer :: field(:,:,:,:)
      real(kind=wp), allocatable :: wf(:,:,:)

      integer :: j,i1,i2,i3
      integer :: io1,io2,io3


      select case(trg_name)
      case('o3')
        field => atm_trg%o3
      case('co2')
        field => atm_trg%co2
      case default
        call finish('avg_field@interpolate_trg', trim(trg_name)//' not implemented.')
      end select
      ! horizontal
      io1 = lbound(field,1)-1
      io2 = lbound(field,2)-1
      io3 = lbound(field,4)-1
      allocate(wf(size(field,1),size(field,2),size(field,3)))
      wf = abs(cos(grid%rlat(:,:,1,:)))
      do j = 1, grid%nz
        avg_h(j) = 0._wp
        wsum     = 0._wp
        do i1 = 1, size(field,1)
          do i2 = 1, size(field,2)
            do i3 = 1, size(field,4)
              avg_h(j) = avg_h(j) + field(io1+i1,io2+i2,j,io3+i3)*wf(i1,i2,i3)
              wsum = wsum + wf(i1,i2,i3)
            end do
          end do
        end do
        avg_h(j) = p_sum(avg_h(j))
        wsum     = p_sum(wsum)
        avg_h(j) = avg_h(j)/wsum
      end do
      deallocate(wf)
      ! vertical
      avg = avg_vert(avg_h, grid%akf(1:grid%nz))
    end function avg_field

    ! Read trg history file
    subroutine read_hist_file
      character(len=300) :: fname
      logical            :: exist
      integer            :: iu, inml, itrg, i
      character(len=20)  :: name_aux
      ! Namelist
      character(len=14) :: last_update   = ''     ! date of last run
      character(len=14) :: last_current  = ''     ! date of last current field
      character(len=14) :: gas_name      = ''     ! gas name
      real(wp)          :: adapt_time    = -1._wp ! adaptation time
      real(wp)          :: max_age       = -1._wp ! maximum age to be considered as "current"
      real(wp)          :: wgt_clim      = -1._wp ! weight of climate (vs. current) in last run
      namelist /trg_history/ last_update, last_current, gas_name, &
                             adapt_time, max_age, wgt_clim
      if (dace%lpio) then
        fname = trim(trg_hist_file)
        call date_tmpl(fname, t=fc_ref_time)
        fname = path_file(input, fname)
        write(*,*) 'Read tracegas history file "'//trim(fname)//'"'
        inquire(file=trim(fname), exist=exist)
        inml = 0
        if (exist) then
          iu = get_unit_number()
          open(file=fname, unit=iu)
          loop_nml: do
            call position_nml ('TRG_HISTORY', unit=iu, lrewind=.false., status=ierr)
            if (ierr == POSITIONED) then
              read(iu, nml=trg_history)
              inml = inml + 1
              call get_gas(gas_name, name=name_aux)
              itrg = -1
              do i = 1, ntrg
                call set_trg(i)
                if (trim(name_aux) == trim(trg_name)) then
                  itrg = i
                  exit
                end if
              end do
              if (itrg > 0) then
                th => trg_hist(itrg)
                write(*,*) '  found history for '//trim(trg_name)
                th%gas_name     = trim(trg_name)
                th%last_update  = time_c(last_update)
                th%last_current = time_c(last_current)
                th%adapt_time   = max(adapt_time, tiny(adapt_time))
                th%max_age      = max_age
                th%wgt_clim     = wgt_clim
              else
                write(*,*) '  ignore namelist ',inml,' gas "'//trim(trg_name)//'" not implemented'
              end if
            else
              exit loop_nml
            end if
          end do loop_nml
          close(iu)
          call return_unit_number(iu)
        else
          call finish('read_hist_file@interpolate_trg', 'History file "'//&
               trim(fname)//' does not exist')
        end if
      end if
      do i = 1, ntrg
        call bcast_trg_history(trg_hist(i), dace%pio)
        l_hist(i) = (trg_hist(i)%gas_name /= '')
      end do
    end subroutine read_hist_file

    ! Update history
    subroutine update_hist(th, t_curr)
      type(t_trg_history), intent(inout) :: th
      type(t_time),        intent(in),   optional :: t_curr
      real(wp) :: dt, wgt_old, wgt_clim_min
      ! This should be done independently of the data amount (i.e. number of spots requiring trg info)
      ! Otherwise we could end up with different weights for different datasets, and this would mean
      ! different (inconsistent) trg fields for different data sets.
      wgt_old = th%wgt_clim
      dt  = days(ana_time - th%last_update)
      if (present(t_curr)) then
        age = days(ana_time - t_curr)
        if (age <=  th%max_age) then
          th%wgt_clim = max(th%wgt_clim - dt/th%adapt_time, 0._wp)
        else
          ! Increase
          th%wgt_clim = min(th%wgt_clim + dt/th%adapt_time, 1._wp)
          wgt_clim_min = min(0._wp + (age-th%max_age)/th%adapt_time, 1._wp)  ! Minimum climate weight for given age
          th%wgt_clim = max(th%wgt_clim, wgt_clim_min)
        end if
        th%last_current = t_curr
      else
        th%wgt_clim = min(th%wgt_clim + dt/th%adapt_time, 1._wp)
        age = days(ana_time - th%last_current)
        wgt_clim_min = min(0._wp + (age-th%max_age)/th%adapt_time, 1._wp)  ! Minimum climate weight for given age
        th%wgt_clim = max(th%wgt_clim, wgt_clim_min)
      end if
      if (dace%lpio) then
        if (wgt_old /= th%wgt_clim) then
          write(*,'(1x,"Adapt weights for ",a,": ",T27,F6.4," -> ",F6.4," (current)")') &
               trim(th%gas_name), 1._wp - wgt_old, 1._wp - th%wgt_clim
          write(*,'(1x,T27,F6.4," -> ",F6.4," (climate)")') wgt_old, th%wgt_clim
          write(*,'(1x,"corresponds to",T27,F9.4," days,  ",F9.4," hours")') &
               abs(wgt_old-th%wgt_clim) * th%adapt_time, abs(wgt_old-th%wgt_clim) * th%adapt_time * 24._wp
        else
          write(*,'(1x,"Current weight for ",a,":",T27,F6.4," (current), ",F6.4," (climate)")') &
               trim(th%gas_name), 1._wp - th%wgt_clim, th%wgt_clim
        end if
      end if
      th%last_update  = ana_time
    end subroutine update_hist

    ! Write trg history file
    subroutine write_hist_file
      character(len=300) :: fname
      integer :: i, iu
      if (dace%lpio) then
        fname = trim(trg_hist_file)
        call date_tmpl(fname, t=ana_time)
        fname = path_file(output, fname)
        write(*,*) 'Write tracegas history file "'//trim(fname)//'"'
        iu = get_unit_number()
        open(file=fname, unit=iu)
        do i = 1, ntrg
          call set_trg(i)
          if (lh) then
            write(iu,'(a)') ''
            write(iu,'(a)') '&trg_history'
            write(iu,'(2x,a,T24,"=",1x,"''",a,"''")') 'gas_name'    , trim(th%gas_name)
            write(iu,'(2x,a,T24,"=",1x,"''",a,"''")') 'last_update' , cyyyymmddhhmmss(th%last_update)
            write(iu,'(2x,a,T24,"=",1x,"''",a,"''")') 'last_current', cyyyymmddhhmmss(th%last_current)
            write(iu,'(2x,a,T24,"=",1x,F7.2)' ) 'max_age'     , th%max_age
            write(iu,'(2x,a,T24,"=",1x,F7.2)' ) 'adapt_time'  , th%adapt_time
            write(iu,'(2x,a,T24,"=",1x,e13.6)') 'wgt_clim'    , th%wgt_clim
            write(iu,'(a)') '/'
          end if
        end do
        close(iu)
        call return_unit_number(iu)
      end if
    end subroutine write_hist_file

    ! Add profile to t_tovs%av and store. Weighting of current and climate profile is done here
    subroutine add_av(trg, wgt)
      real(wp), intent(in) :: trg(:)
      real(wp), intent(in) :: wgt
      real(wp) :: wold
      integer  :: j
      wold = 1._wp - wgt
      if (wold > 0._wp) then
        av(1:nl,trg_ind) = real(wgt * trg(1:nl) + wold * av(1:nl,trg_ind), kind=sp)
      else
        av(1:nl,trg_ind) = real(trg(1:nl), kind=sp)
      end if
      if (ldeb(spt)) then
        do j = 1, nl
          write(usd,*) dpref,proc,trg_name,j,av(j,trg_ind)
        end do
      end if
      call store_av
    end subroutine add_av

    ! Store t_tovs%av
    subroutine store_av
      tovs_io = ior(tovs_io, TTOVS_BASE)
      call store(obs%o(ib), spt, tovs=ttovs, av=av(1:nl,1:ttovs%nav), tovs_io=tovs_io)
      call destruct(ttovs)
    end subroutine store_av

! specific MPI-bcast for derived type t_trg_history
#undef  VECTOR
#undef  DERIVED
#define DERIVED type(t_trg_history)
#define p_bcast_DERIVED bcast_trg_history
#undef  MPI_TYPE
#include "p_bcast.incf"
#undef  DERIVED
#undef  MPI_TYPE

  end subroutine interpolate_trg

!------------------------------------------------------------------------------

  subroutine interpolate_strat(obs)
    type(t_obs_set)  ,intent(inout)  :: obs  ! observation data type
    !-------------------------------------------------------------------------
    ! read extra levels above movel top and store in t_tovs%av
    !   interp_strato: how to fill the extra levels above model top
    !                  Bit 0 = 1 : get T,q from file_strato
    !                  Bit 1 = 2 : set T,q to extrap. constant
    !                  Bit 2 = 4 : set T,q to US standard atmosphere
    !                  If method 1 (bit 0, read from file file_strato) fails, the
    !                  other two method might serve as a fallback
    !-------------------------------------------------------------------------
    character(len=17), parameter    :: proc = 'interpolate_strat'
    ! General stuff
    type (t_spot),      pointer     :: spt      !
    logical                         :: ok
    integer                         :: nb       ! number of boxes
    integer                         :: ib, is   ! box, spot index
    integer                         :: i
    ! Spot counting
    integer                         :: nsp, nsp_all   ! Number of spots to be processed
    integer,            allocatable :: nsp_box(:)
    logical,            allocatable :: spt_mask(:)
    ! interp_strato bits
    integer,           parameter    :: ITS_FILE  = 0
    integer,           parameter    :: ITS_CONST = 1
    integer,           parameter    :: ITS_USS   = 2
    ! state of extra atmosphere
    type (t_grid),      pointer     :: grid     ! grid information
    type (t_atm)                    :: strat    ! background
    type (t_inventory), pointer     :: invt(:)  ! GRIB file inventory
    ! horizontal interpolation
    type(t_mcols),      allocatable :: mc(:)                ! column(s) from orig grid
    type(t_cols),       allocatable :: cbgb(:)              ! column(s) from orig grid
    type(t_col),        pointer     :: c                    ! column(s) from orig grid
    type(t_hic),        allocatable :: h(:)                 ! interpolation indices/weights
    integer                         :: ih                   ! index in h(:)
    integer                         :: ii, ic
    integer(i8)                     :: iatm
    ! t_tovs stuff
    type(t_tovs)                    :: ttovs     ! TOVS specific data type
    real(tpp),          allocatable :: av(:,:)
    type(t_rad_set),    pointer     :: rs
    integer                         :: it_t, it_q
    ! interpolation
    real(wp),           allocatable :: p_g(:), v_g(:,:) ! P and other vars on orig. grid
    real(wp),           allocatable :: p_i(:), v_i(:,:) ! P and other vars in interp. space
    integer,            parameter   :: i_tv = 1         ! Tv index in v_g/v_i
    integer,            parameter   :: i_rh = 2         ! RH index in v_g/v_i
    integer,            parameter   :: i_q  = 3         ! Q  index in v_g/v_i
    integer                         :: ng               ! number of levels in file
    integer                         :: nl               ! number levels in interp. space
    integer                         :: level_int        ! number of levels to filled
    logical                         :: l_vint_rt
    ! transformation
    integer,            allocatable :: i_fail(:)
    real(wp),           allocatable :: wgt(:)
    real(wp),           allocatable :: t(:)
    real(wp),           allocatable :: q(:)
    ! Blending
    integer                         :: nld_t, nld_h
    integer                         :: level_blend
    integer                         :: l0
    real(wp)                        :: w

    if (interp_strato == 0) return
    nld_t_max = p_max(nld_t_max)
    nld_h_max = p_max(nld_h_max)
    if (nld_t_max <= 0 .and. nld_h_max <= 0) return

    nb = size (obs% o)

    ! Count OT_RAD spots
    allocate(nsp_box(nb))
    nsp_box = 0
    do ib = 1, nb
      if (dace%pe /= obs%o(ib)%pe) cycle
      do is = 1, obs%o(ib)%n_spot
        spt => obs%o(ib)%spot(is)
        if (spt% hd% obstype /= OT_RAD) cycle
        nsp_box(ib) = nsp_box(ib) + 1
      end do
    end do
    nsp = sum(nsp_box(1:nb))
    nsp_all = p_sum(nsp)
    if (nsp_all <= 0) return
    ! Count OT_RAD spots, that require interpolation
    allocate(spt_mask(nsp))
    spt_mask = .false.
    nsp_box  = 0
    ih = 0
    do ib = 1, nb
      if (dace%pe /= obs%o(ib)%pe) cycle
      do is = 1, obs%o(ib)%n_spot
        spt => obs%o(ib)%spot(is)
        if (spt% hd% obstype /= OT_RAD) cycle
        ih = ih + 1
        call load(obs%o(ib), spt, tovs=ttovs, tovs_io=0)
        if (ttovs%nl_st > 0) then
          spt_mask(ih) = .true.
          nsp_box(ib) = nsp_box(ib) + 1
        end if
      end do
    end do
    nsp = sum(nsp_box(1:nb))

    ok = .false.
    ! Interpolate atmosphere from file
    if (btest(interp_strato, ITS_FILE)) then
      if (file_strato /= '') then
        if (dace%lpio) then
          write(*,*)
          write(*,*) 'Read file_strato='//trim(file_strato)
        end if
        nullify   (invt)
        call get_inventory (invt, file_strato)
        if (size(invt) > 0) then
          call print_inventory (invt,  first=.true.)
          if (dace%lpio) write(*,*)
          allocate(grid)
          call read(grid, file_strato, invt=invt, geosp=.false., lsm=.false., &
               nproc1=nproc1, nproc2=nproc2, comm=dace% comm, gridfile=grid_file_strato)
          call print(grid,  verbose=.true.)
          call construct(strat ,grid)
          call read     (strat, file_strato, invt, fields = 't q ps psr pf ph', optionals='pf ph')
          call print    (strat, verbose=.true.)
          if (.not.associated(strat%pf)) then
            call set_ph(strat)
            call set_p(strat)
          end if
          call set_geo  (strat, geof=.true.)
          call set_tv   (strat)
          call allocate (strat,'rh')
          strat% rh = rh_q (strat% q, strat% t, strat% pf)
          if (any(strat% rh < 0._wp)) call finish (proc,'rh < 0.')
          allocate(mc(nb), cbgb(nb), h(nsp))
          ih = 0
          do ib = 1, nb
            mc(ib)% pe = obs%o(ib)%pe
            if (dace%pe /= obs%o(ib)%pe) cycle
            allocate(mc(ib)% idx(grid%lbg(1):grid%ubg(1),grid%lbg(2):grid%ubg(2),grid%lbg(4):grid%ubg(4),1))
            allocate(mc(ib)% c(4*nsp_box(ib)))
            mc(ib)% idx = 0
            mc(ib)% n = 0
            do is = 1, obs%o(ib)%n_spot
              spt => obs%o(ib)%spot(is)
              if (spt% hd% obstype /= OT_RAD) cycle
              ih = ih + 1
              if (spt_mask(ih)) then
                call idx_init (   &
                     spt% col% c, &! <-  column descriptor
                     h(ih),       &!  -> interpolation coefficients
                     mc(ib),      &! <-> model column descriptors
                     0_i8,        &! <-  fields required
                     0,           &! <-  tracers required
                     grid,        &! <-  model grid
                     1,           &! <-  time slot
                     0._wp        )! <-  time interpolation weight
              end if
            end do
          end do
          iatm = COL_P + COL_TV + COL_RH + COL_Q
          call get_cols(mc, strat, cbgb, iatm=iatm)

          ng = grid%nz
          allocate(p_g(ng), v_g(3,ng), p_i(mx_nlev), v_i(3, mx_nlev), &
                   av(mx_nlev,mx_nav), i_fail(mx_nlev), wgt(mx_nlev), &
                   t(mx_nlev),q(mx_nlev))
          ih  = 0
          do ib = 1, nb
            if (dace%pe /= obs%o(ib)%pe) cycle
            do is = 1, obs%o(ib)%n_spot
              spt => obs%o(ib)%spot(is)
              if (spt% hd% obstype /= OT_RAD) cycle
              ih = ih + 1
              if (spt_mask(ih)) then
                call prep_spot
                ! horiz. interp.
                p_g = 0._wp
                v_g = 0._wp
                do ic = 1, size(h(ih)%imc,1)
                  ii = h(ih)% imc(ic,1)
                  if (ii==0) exit
                  c => cbgb(ib)% col(ii)
                  p_g(     :) = p_g(     :) + h(ih)% w(ic) * c%p (:)
                  v_g(i_tv,:) = v_g(i_tv,:) + h(ih)% w(ic) * c%tv(:)
                  v_g(i_rh,:) = v_g(i_rh,:) + h(ih)% w(ic) * c%rh(:)
                  v_g(i_q ,:) = v_g(i_q ,:) + h(ih)% w(ic) * c%q (:)
                end do
                ! vert. interp.
                call interpolate_rttov(v_i(:,1:nl), y_in=v_g(:,:), &
                                       lnp_out=log(p_i(1:nl)), lnp_in=log(p_g))
                ! humidity cases
                select case(int_rad_hum)
                case(0,2)
#ifdef __NEC__
                  call tq_tvrh_vec2 &
#else
                  call tq_tvrh_vec  &
#endif
                     (nl,t(1:nl),q(1:nl),v_i(i_tv,1:nl),v_i(i_rh,1:nl),p_i(1:nl), i_fail=i_fail(1:nl))
                  do i = 1, nl
                    if (i_fail(i) < 0) then
                      write(0,*) 'i,p,tv,rh',i,p_i(i),v_i(i_tv,i),v_i(i_rh,i),i_fail(i)
                      call finish(proc, 'tq_tvrh failed')
                    end if
                  end do
                  if (int_rad_hum == 2) then
                    where(log(p_i(1:nl)) > pblend(1))
                      wgt(1:nl) = 1._wp
                    elsewhere(log(p_i(1:nl)) > pblend(2))
                      wgt(1:nl) = cos((log(p_i(1:nl)) -pblend(1))/(pblend(2)-pblend(1))*pi*0.5_wp)**2
                    elsewhere
                      wgt(1:nl) = 0._wp
                    end where
                    v_i(i_rh,1:nl) = q(1:nl) ! i_rh is now Q calculated from RH
                    q(1:nl) = v_i(i_rh,1:nl) * wgt + v_i(i_q,1:nl) * (1._wp-wgt)
                    t(1:nl) = t_tv_q(v_i(i_tv,1:nl), q(1:nl))
                  end if
                case(1,3)
                  q(1:nl) = v_i(i_q,1:nl)
                  t(1:nl) = t_tv_q(v_i(i_tv,1:nl), q(1:nl))
                case default
                  call finish(proc, 'value for int_rad_hum not implemented so far')
                end select
                ! Add to/blend with interpolated model profile
                level_blend = nl+1
                do i = 1, nl
                  if (p_i(i) > p_blend_ext) then
                    level_blend = i
                    exit
                  end if
                end do
                if (ldeb(spt)) write(usd,*) dpref,proc,' levs',nl,level_int,level_blend
                do i = 1, min(level_int, level_blend-1)
                  av(i,it_t) = t(i)
                  av(i,it_q) = q(i)
                end do
                l0 = max(level_int+1,1)
                if (level_blend <= level_int) then
                  ! Blending
                  do i = level_blend, level_int
                    w = cos((log(p_i(i))      - log(p_i(l0))) / &
                            (log(p_blend_ext) - log(p_i(l0))) *pi*0.5_wp)**2
                    av(i,it_t) = w * av(l0,it_t) + (1._wp - w) * t(i)
                    av(i,it_q) = w * av(l0,it_q) + (1._wp - w) * q(i)
                  end do
                end if
                if (ldeb(spt)) then
                  do i = 1, level_int
                   write(usd,*) dpref,proc,' prof',i,p_i(i),av(i,it_t),av(i,it_q)
                  end do
                end if
                call store(obs%o(ib), spt, tovs=ttovs, tovs_io=0, av=av)
              end if
            end do
          end do

          call destruct(mc)
          call dealloc_cols(cbgb)
          call destruct(grid)
          call destruct(strat)
          ok = .true.
        else
          write(0,*) 'file_strat='//trim(file_strato)//' inventory empty !!!'
          write(6,*) 'file_strat='//trim(file_strato)//' inventory empty !!!'
        end if
      else
        call finish(proc, 'interp_strato/=0 but file_strato not set')
      end if
    endif

    if (.not.ok .and. btest(interp_strato, ITS_CONST)) then
      if (dace%lpio) then
        write(6,*) 'Extend TOVS-profiles with constant extrapolation'
      end if
      ih = 0
      do ib = 1, nb
        if (dace%pe /= obs%o(ib)%pe) cycle
        do is = 1, obs%o(ib)%n_spot
          spt => obs%o(ib)%spot(is)
          if (spt% hd% obstype /= OT_RAD) cycle
          ih = ih + 1
          if (spt_mask(ih)) then
            call prep_spot
            do i = 1, level_int
              av(i,it_t) = av(level_int+1,it_t)
              av(i,it_q) = av(level_int+1,it_q)
            end do
            call store(obs%o(ib), spt, tovs=ttovs, tovs_io=0, av=av)
          end if
        end do
      end do
      ok = .true.
    end if

    if (.not.ok .and. btest(interp_strato, ITS_USS)) then
      call finish(proc,'ITS_USS not implemented so far')
      ok = .true.
    end if

    if (.not.ok .and. interp_strato /= 0) &
         call finish(proc, 'FAILED! Either set interp_strato=0 or disable radiances.')

  contains

    subroutine prep_spot
      call load(obs%o(ib), spt, tovs=ttovs, tovs_io=0, av=av, rs=rs)
      it_t = ttovs%i_t
      if (it_t <= 0) call finish(proc, 'No T-entry in t_tovs%av')
      it_q = ttovs%i_q
      if (it_q <= 0) call finish(proc, 'No Q-entry in t_tovs%av')
      l_vint_rt = rs%gopts%lev_mode <= 0
      nl   = ttovs%nlev
      nld_t = ttovs%nld_t
      nld_h = ttovs%nld_h
      level_int = ttovs%nl_st
      if (l_vint_rt) then
        p_i(1:nl) = preslev(1:nl)
      else
        if (ttovs%i_p <= 0) call finish(proc, 'No P in t_tovs%av')
        p_i(1:nl) = real(av(1:nl, ttovs%i_p), wp)
      end if
    end subroutine prep_spot

  end subroutine interpolate_strat


!------------------------------------------------------------------------------
  subroutine get_obs (obs, phsp)
  !--------------------------------------------
  ! get observations from observation data type
  !--------------------------------------------
  type (t_obs)    ,intent(in)    :: obs (:)
  type (t_vector) ,intent(inout) :: phsp

    integer :: l      ! vector segment index
    !--------------------------
    ! loop over vector segments
    !--------------------------
    do l = 1, phsp% n_s
      !---------------------------------------
      ! handle only segments on this processor
      !---------------------------------------
      if (phsp% global .or. phsp% s(l)% pe == dace% pe) then
        if (.not. associated (phsp% s(l)% x)) &
          allocate (phsp%s(l)% x(phsp% s(l)% n))
        if (phsp% s(l)% n > 0) &
          phsp% s(l)% x(:) = obs(l)% body% o
      else
        if (associated(phsp% s(l)% x)) deallocate (phsp% s(l)% x)
      endif
    end do
  end subroutine get_obs
!------------------------------------------------------------------------------
  subroutine get_plev (obs, plev, fill)
  !----------------------------------------------
  ! get pressure level from observation data type
  !----------------------------------------------
  type (t_obs)    ,intent(in)           :: obs (:)
  type (t_vector) ,intent(inout)        :: plev
  real (wp)       ,intent(in) ,optional :: fill

    integer :: l      ! vector segment index
    !--------------------------
    ! loop over vector segments
    !--------------------------
    do l = 1, plev% n_s
      !---------------------------------------
      ! handle only segments on this processor
      !---------------------------------------
      if (plev% global .or. plev% s(l)% pe == dace% pe) then
        if (.not. associated (plev% s(l)% x)) &
          allocate (plev%s(l)% x(plev% s(l)% n))
        if (plev% s(l)% n > 0) then
          plev% s(l)% x(:) = obs(l)% body% plev
          if (present (fill)) then
            where (plev% s(l)% x(:) <= 0._wp)  plev% s(l)% x(:) = fill
          endif
        endif
      else
        if (associated(plev% s(l)% x)) deallocate (plev% s(l)% x)
      endif
    end do
  end subroutine get_plev
!------------------------------------------------------------------------------
  subroutine get_obstype (obs, phsp)
  !--------------------------------------------
  ! get observations from observation data type
  !--------------------------------------------
  type (t_obs)     ,intent(in)    :: obs (:)
  type (t_ivector) ,intent(inout) :: phsp

    integer :: l      ! vector segment index
    integer :: i      ! report index
    integer :: i1, in
    !--------------------------
    ! loop over vector segments
    !--------------------------
    do l = 1, phsp% n_s
      !---------------------------------------
      ! handle only segments on this processor
      !---------------------------------------
      if (phsp% global .or. phsp% s(l)% pe == dace% pe) then
        if (.not. associated (phsp% s(l)% x)) &
          allocate (phsp%s(l)% x(phsp% s(l)% n))
        if (phsp% s(l)% n > 0) then
          do i = 1, obs(l)% n_spot
            i1 = obs(l)% spot(i)% o% i + 1
            in = obs(l)% spot(i)% o% i + obs(l)% spot(i)% o% n
            phsp% s(l)% x(i1:in) = obs(l)% spot(i)% hd% obstype
          end do
        endif
      else
        if (associated(phsp% s(l)% x)) deallocate (phsp% s(l)% x)
      endif
    end do
  end subroutine get_obstype
!------------------------------------------------------------------------------
  subroutine get_id (obs, phsp)
  !--------------------------------------------
  ! get observations from observation data type
  !--------------------------------------------
  type (t_obs)     ,intent(in)    :: obs (:)
  type (t_ivector) ,intent(inout) :: phsp

    integer :: l      ! vector segment index
    integer :: i      ! report index
    integer :: i1, in
    !--------------------------
    ! loop over vector segments
    !--------------------------
    do l = 1, phsp% n_s
      !---------------------------------------
      ! handle only segments on this processor
      !---------------------------------------
      if (phsp% global .or. phsp% s(l)% pe == dace% pe) then
        if (.not. associated (phsp% s(l)% x)) &
          allocate (phsp%s(l)% x(phsp% s(l)% n))
        if (phsp% s(l)% n > 0) then
          do i = 1, obs(l)% n_spot
            i1 = obs(l)% spot(i)% o% i + 1
            in = obs(l)% spot(i)% o% i + obs(l)% spot(i)% o% n
            phsp% s(l)% x(i1:in) = obs(l)% spot(i)% hd% id
          end do
        endif
      else
        if (associated(phsp% s(l)% x)) deallocate (phsp% s(l)% x)
      endif
    end do
  end subroutine get_id
!------------------------------------------------------------------------------
  subroutine get_varno (obs, phsp)
  !--------------------------------------------
  ! get observations from observation data type
  !--------------------------------------------
  type (t_obs)     ,intent(in)    :: obs (:)
  type (t_ivector) ,intent(inout) :: phsp

    integer :: l      ! vector segment index
    !--------------------------
    ! loop over vector segments
    !--------------------------
    do l = 1, phsp% n_s
      !---------------------------------------
      ! handle only segments on this processor
      !---------------------------------------
      if (phsp% global .or. phsp% s(l)% pe == dace% pe) then
        if (.not. associated (phsp% s(l)% x)) &
          allocate (phsp%s(l)% x(phsp% s(l)% n))
        if (phsp% s(l)% n > 0) then
            phsp% s(l)% x(:) = obs(l)% varno
        endif
      else
        if (associated(phsp% s(l)% x)) deallocate (phsp% s(l)% x)
      endif
    end do
  end subroutine get_varno
!------------------------------------------------------------------------------
  subroutine put_bg (obs, bg, bge)
  !------------------------------------------
  ! store background in observation data type
  !------------------------------------------
  type (t_obs)    ,intent(inout) :: obs (:)  ! observation data type
  type (t_vector) ,intent(in)    :: bg       ! background value
  type (t_vector) ,intent(in)    :: bge      ! background error

    integer :: l      ! vector segment index
    !--------------------------
    ! loop over vector segments
    !--------------------------
    do l = 1, size(obs)
      !---------------------------------------
      ! handle only segments on this processor
      !---------------------------------------
      if (obs(l)% pe == dace% pe .and. obs(l)% n_obs > 0) then
        obs(l)% body% bg = bg % s(l)% x(:)
        obs(l)% body% eb = bge% s(l)% x(:)
      endif
    end do
  end subroutine put_bg
!------------------------------------------------------------------------------
  subroutine set_obs_random (obs, phsp, bg, R)
  type (t_obs)    ,intent(inout) :: obs (:) ! observation data type
  type (t_vector) ,intent(inout) :: phsp    ! observations
  type (t_vector) ,intent(in)    :: bg      ! background
  type (t_matrix) ,intent(in)    :: R       ! observation error matrix

    integer               :: l      ! vector segment index
    !--------------------------
    ! loop over vector segments
    !--------------------------
    do l = 1, phsp% n_s
      !---------------------------------------
      ! handle only segments on this processor
      !---------------------------------------
      if (phsp% s(l)% pe == dace% pe) then
        if (.not. associated (phsp% s(l)% x)) &
          allocate (phsp%s(l)% x(phsp% s(l)% n))
      else
        if (associated(phsp% s(l)% x)) deallocate (phsp% s(l)% x)
      endif
    end do
    call random_gauss(phsp)
    phsp = sqrt(R) * phsp + bg
    do l = 1, phsp% n_s
      if (phsp% s(l)% pe == dace% pe) obs(l)% body(:)% o = phsp% s(l)% x(:)
    end do
  end subroutine set_obs_random
!==============================================================================
  subroutine check_ob_fg (o, e_f, e_o, sig, idtwin)
  type (t_vector) ,intent(inout) :: o      ! observation - forecast
  type (t_vector) ,intent(in)    :: e_f    ! forecast    error (stddev.)
  type (t_vector) ,intent(in)    :: e_o    ! observation error (stddev.)
  real (wp)       ,intent(in)    :: sig    ! bound to reject observation
  integer         ,intent(in)    :: idtwin ! identical twin exp. flag

    integer               :: l, i ! vector segment index
    integer               :: nrej ! number of rejected observations
    real(wp) ,allocatable :: d (:)
    real(wp)              :: s2

    !-----------------------
    ! double twin experiment
    !-----------------------
    nrej = 0
    select case (idtwin)
    case (1)
      o = 1._wp
    case (2)
      call random_gauss (o)
    case default
      !-----------------
      ! background check
      !-----------------
      if (sig > 0._wp) then
        !--------------------------------------------
        ! background check, loop over vector segments
        !--------------------------------------------
        s2 = sig**2
        do l = 1, o% n_s
          !---------------------------------------
          ! handle only segments on this processor
          !---------------------------------------
          if (o% s(l)% pe == dace% pe) then
            allocate (d(o% s(l)% n))
            d = s2 * e_f% s(l)% x **2 + e_o% s(l)% x **2
            do i = 1, o% s(l)% n
              if (o% s(l)% x(i)**2 > d(i)) then
                o% s(l)% x(i) = 0._wp
                nrej = nrej + 1
              endif
            end do
            deallocate (d)
          endif
          !-----------------------
          ! broadcast to other PEs
          !-----------------------
          if (o% global) call p_bcast (o% s(l)% x, o% s(l)% pe)
        end do
      endif
    end select
!   write(6,*) 'check_ob_fg, pe',dace% pe,': rejected observations:',nrej
  end subroutine check_ob_fg
!=============================================================================
  subroutine write_intp (obs, name, bg, ana, e_bg)
  type(t_obs)      ,intent(in) :: obs(:) ! observations
  character(len=*) ,intent(in) :: name   ! file basename
  type(t_vector)   ,intent(in) :: bg     ! background
  type(t_vector)   ,intent(in) :: ana    ! analysis
  type(t_vector)   ,intent(in) :: e_bg   ! background error


    integer             :: ib, is, i, n, io
    character (len=128) :: file
    real(wp)            :: rlev
    character (len=12)  :: cdate
    real(wp)            :: dummy_ofc
!   real(wp)            :: dummy_aff
    real(wp)            :: dummy_z
    real(wp)            :: dummy_eo
    real(wp)            :: dummy_wvqc

    type(t_obs),pointer :: o
    type(t_obs),target  :: tmp
    target              :: obs

    if (dace% lpio) then
      io   = get_unit_number()
      file = trim(name)//'.info'
      file = path_file (aux,file)
      open  (io , file=file)
      n  = 0
      write(io,'(a)') &
'#      n  obs_nr     box    name  typ bft sbt obt       lat       lon&
&       lev              fc            o-fc            a-fc            af_f&
&               z            e_bg             e_o            w_qc  obs_id&
&         date  obs_id'

      dummy_ofc   = 0._wp
!     dummy_aff   = 0._wp
      dummy_z     = 0._wp
      dummy_eo    = 0.1_wp
      dummy_wvqc  = 1._wp
    endif

    do ib = 1,size(obs)
      !--------------------
      ! send / receive data
      !--------------------
      if (obs(ib)% bcast) then
        o => obs(ib)
      else
        if (dace% lpio) then
          if (obs(ib)% pe==dace% pio) then
            o => obs(ib)
          else
            call p_recv (tmp, obs(ib)% pe)
            o => tmp
          endif
        else
          if (obs(ib)% pe==dace% pe) call p_send (obs(ib), dace% pio)
        endif
      endif
      call gather(bg   ,dace% pio ,isgm=ib)
      call gather(ana  ,dace% pio ,isgm=ib)
      call gather(e_bg ,dace% pio ,isgm=ib)
      !------
      ! write
      !------
      if (dace% lpio) then
        do is=1,o% n_spot
!         if (o% spot(is)% hd% obstype /= RT_RAD) cycle
          cdate = cyyyymmddhhmm (o%spot(is)% actual_time)
          do i = o%spot(is)%i%i+1, &
                 o%spot(is)%i%i+o%spot(is)%i%n
            n = n + 1
            rlev = exp(o%lev(i))

!print *,'n          ',                  n
!print *,'id         ',o%spot(is)%       id
!print *,'ib         ',                  ib
!print *,'statid     ',o%spot(is)%       statid
!print *,'modtype    ',o%spot(is)%hd%    modtype
!print *,'buf_type   ',o%spot(is)%hd%    buf_type
!print *,'buf_subtype',o%spot(is)%hd%    buf_subtype
!print *,'t_int      ',o%                t_int (i)
!print *,'dlat       ',o%spot(is)%col%c% dlat
!print *,'dlon       ',o%spot(is)%col%c% dlon
!print *,'rlev       ',                  rlev
!print *,'bg         ',bg%   s(ib)% x(i)
!print *,'dummy_ofc  ',                  dummy_ofc
!print *,'ana-bg     ',ana%  s(ib)% x(i) -  bg%   s(ib)% x(i)
!print *,'ana-bg     ',ana%  s(ib)% x(i) -  bg%   s(ib)% x(i)
!print *,'dummy_z    ',                  dummy_z
!print *,'e_bg       ',e_bg% s(ib)% x(i)
!print *,'dummy_eo   ',dummy_eo
!print *,'dummy_wvqc ',                  dummy_wvqc
!print *,'cdate      ',                  cdate
!print *,'record     ',o%spot(is)%hd%    record

            write(io,'(3i8,1x,a,4i5,3f10.2,8f16.7,1x,a,1x,i8)')   &
                        n                                        ,&
                        o%spot(is)%       id                     ,&
                        ib                                       ,&
                        o%spot(is)%       statid                 ,&
                        o%spot(is)%hd%    modtype                ,&
                        o%spot(is)%hd%    buf_type               ,&
                        o%spot(is)%hd%    buf_subtype            ,&
                        o%                t_int (i)              ,&
                        o%spot(is)%col%c% dlat                   ,&
                        o%spot(is)%col%c% dlon                   ,&
                        rlev                                     ,&
                        bg%   s(ib)% x(i)                        ,&
                        dummy_ofc                                ,&
                        ana%  s(ib)% x(i) -  bg%   s(ib)% x(i)   ,&
                        ana%  s(ib)% x(i) -  bg%   s(ib)% x(i)   ,&
                        dummy_z                                  ,&
                        e_bg% s(ib)% x(i)                        ,&
                        dummy_eo                                 ,&
                        dummy_wvqc                               ,&
                        cdate                                    ,&
                        o%spot(is)%hd%    record
          end do
        end do
      end if
      if (dace% lpio) then
        if (o%pe /= dace% pe .and..not. o% bcast) then
          deallocate (tmp% spot)
          call destruct (tmp)
        endif
        call release_mem (bg)
        call release_mem (ana)
        call release_mem (e_bg)
      endif
    end do
    if (dace% lpio) then
      close (io)
      call return_unit_number(io)
    endif
  end subroutine write_intp
!==============================================================================
  subroutine set_new_x (new_H, x_xb, u_ub, Hi, Hinv, z, HPbHt, xi_xb, ui_ub, &
                        obs, time, pbaprx,e_f)
  integer         ,intent(in)           :: new_H ! flag
  type (t_vector) ,intent(inout)        :: x_xb  ! x - bg output (int. space)
  type (t_vector) ,intent(inout)        :: u_ub  ! y - bg   (observ. space)
  type (t_matrix) ,intent(in)           :: Hi    ! linear  observation operator
  type (t_matrix) ,intent(in)           :: Hinv  ! inverse observation operator
  type (t_vector) ,intent(in) ,optional :: z     ! solution vector (obsspace)
  type (t_matrix) ,intent(in) ,optional :: HPbHt ! backg. error covar. matrix
  type (t_vector) ,intent(in) ,optional :: xi_xb ! reference vector
  type (t_vector) ,intent(in) ,optional :: ui_ub ! reference vector
  type (t_obs_set),intent(inout),optional :: obs   ! observation data type
  type (t_time)   ,intent(in) ,optional :: time  ! analysis time
  integer         ,intent(in) ,optional :: pbaprx! Pb approximation flag
  type (t_vector) ,intent(in) ,optional :: e_f   ! forecast    error (stddev.)

  !===========================================================================
  ! The Routine derives a new x-bg (in interpolation vector space) from z
  ! (solution vector in dual space). y-bg = HPbHt z is updatated.
  !
  ! The two vector spaces X (interpolation space) and Y (observation
  ! space) are related by y-bg=H (x-bg) there H is the linearized
  ! observation operator. bg is the value of the background state in the
  ! respective vector space. x may be obtained from y either by inverting
  ! H (in the case that H is not invertable the generalised Inverse is
  ! used, leading to an approximation) or by calculating PbHt * z. The
  ! latter method may be expensive becouse the matrix Pb must be
  ! recalculated.
  !
  ! For the approximate solution reference values ui_ub, xi_xb (used for
  ! linearizing H) may be given to improve the conditioning.
  !
  ! The approximate solution is obtained if mod(new_H,10)==1
  ! The exact solution is optained       if mod(new_H,10)==2
  ! If B is available as operator (covm%valid>=5) the exact procedure is used
  !
  !   new_x=1: x = H^(-1) (H P_b H^t) z  (H^(-1): generalised inverse)
  !
  !         2: x =           P_b H^t  z  (set up new P_b)
  !
  !===========================================================================

    integer                  :: new_x ! flag for deriving x (the argument)
    type (t_vector)          :: xz    ! temporary
    type (t_vector)          :: dx    ! temporary
    type (t_vector)          :: dy    ! temporary
    type (t_vector)          :: du    ! temporary (dummy observations)
    type (t_vector) ,pointer :: Bz
    target                   :: u_ub

    !--------------------------------------------
    ! operator B matrix:
    ! always use exact method as Bii * z is cheap
    !--------------------------------------------
    if (present(z) .and. covm% valid >= 5) then
      allocate (Bz)
      call construct (Bz, u_ub% info)
      call construct (xz, x_xb% info)
      call construct (du, obs%  di)
      xz = z * Hi
      call apply_B_ii_2d (xz,xz,e_f,du)
      Bz = Hi * xz
      u_ub = Bz
      x_xb = xz
      call destruct (du)
      call destruct (xz)
      call destruct (Bz)
      deallocate    (Bz)
    else
      !---------------
      ! u_ub = HPbHt z
      !---------------
      if (present(z)) then
        allocate (Bz)
        call construct (Bz, u_ub% info)
        Bz = HPbHt * z
      else
        Bz => u_ub
      endif
      !-------------------------------------
      ! obtain vector in interpolation space
      !-------------------------------------
      new_x = mod (new_H , 10)
      select case (new_x)
      case (1)
        !---------------------
        ! approximate solution
        !---------------------
        if (present(ui_ub).and.present(xi_xb)) then
          call construct (dx, x_xb% info)
          call construct (dy, u_ub% info)
          dy = Bz - ui_ub
          dx = Hinv * dy
          x_xb  = dx + xi_xb
          call destruct (dx)
          call destruct (dy)
        else
          x_xb = Hinv * Bz
        endif
      case (2)
       if (.not.present(obs)   .or.&
           .not.present(time)  .or.&
           .not.present(pbaprx))   &
         call finish('set_new_x','obs,time,pbaprx not present for new_x==2')
       call construct (xz, x_xb% info)
       xz = z * Hi
       call Pb_times_z (x_xb, xz, obs% o, rearth, time, pbaprx)
       call destruct  (xz)
      case default
        write(0,*)  'set_new_x: new_x =',new_x
        call finish('set_new_x','invalid value of new_x')
      end select
      !---------------------
      ! store u_ub = HPbHt z
      !---------------------
      if (present(z)) then
        u_ub = Bz
        call destruct (Bz)
        deallocate    (Bz)
      endif
    endif
  end subroutine set_new_x
!------------------------------------------------------------------------------
  subroutine set_new_H (new_H, obs, xbg, xi_xb, ui_ub, Hi)
  integer         ,intent(in)    :: new_H ! flag
  type (t_obs_set),intent(inout) :: obs   ! observation data to update
  type (t_vector) ,intent(in)    :: xbg   ! background (in interpolation space)
  type (t_vector) ,intent(in)    :: xi_xb ! xi-bg      (reference value int.sp)
  type (t_vector) ,intent(inout) :: ui_ub ! yi-bg      (reference value obs.sp)
  type (t_matrix) ,intent(inout) :: Hi    ! linearized observation operator
  !===========================================================================
  ! Derive new Jakobi matrices for the nonlinear observation operators
  !   and store in 'obs'.
  !
  !   new_H = new_x + 10 * new_k
  !
  !   new_k=0: no update
  !         1: no update, derive      y=H(x)
  !         2: update          ui_ub, y
  !         3: update       H, ui_ub, y
  !===========================================================================

    integer                  :: new_K ! flag for deriving H (Jakobi matrix)
    type (t_vector)          :: x, y
!   target                   :: xi_xb

    !------------
    ! update H, y
    !------------
    new_K =      new_H / 10
    call construct (x, xi_xb% info)
    call construct (y, ui_ub% info)
    x     = xi_xb + xbg
    select case (new_k)
    case (0)
    case (1,2)
      !---------
      ! derive y
      !---------
      call process_obs (TSK_YH, obs, xi=x, y=y)
      if (new_k==2) then
        !------------
        ! set offsets
        !------------
        obs% l% x = x
        obs% l% y = y
      end if
    case (3)
      !------
      ! set H
      !------
      call process_obs (TSK_K, obs, xi=x)
!     call scatter_K   (obs)              ! not required for
!     call release_mem (obs% o)           ! operator B approach
      call set_H (Hi ,obs ,'f')
!     call apply_h (obs, x, y, 'o')
    case default
      call finish('set_new_H','invalid value for new_K')
    end select
    !-------------
    ! derive ui_ub
    !-------------
    select case (new_k)
    case (2:)
      ui_ub = Hi * xi_xb
    end select

    call destruct (x)
    call destruct (y)

  end subroutine set_new_H
!==============================================================================
end module mo_psasutil
