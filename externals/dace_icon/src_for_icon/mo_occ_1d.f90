!
!+ GNSS occultation 1D-observation-operator (Abel transform) routines.
!
MODULE mo_occ_1d
!
! Description:
!   GNSS occultation 1D-observation-operator (Abel transform) routines.
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
!  Changes for verification mode (TSK_SET_CHR)
! V1_9         2010/04/20 Andreas Rhodin
!  use same hydrostatic relationship as occ-library (consist nonlinear/adjoint)
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Harald Anlauf
!  improved diagnostics and checks on GPSRO fg and obs
! V1_15        2011/12/06 Harald Anlauf
!  option to remove observations in the assimilation step
!  modify status flag for 'borderline' rays
! V1_22        2013-02-13 Harald Anlauf
!  Add vertical coordinate type for ICON
!  Changes for GPSRO ray-tracer
!  Adjust PCD quality information to ROPP version 6
! V1_27        2013-11-08 Harald Anlauf
!  ICON horizontal interpolation, remove ak,bk, account for waterload
! V1_28        2014/02/26 Harald Anlauf
!  Bugfix (water load in adjoint for 1d operator)
! V1_42        2015-06-08 Harald Anlauf
!  horint_mode; adapt to MEC: temporal interpolation, feedback file I/O
! V1_43        2015-08-19 Andreas Rhodin
!  non-linear/adjoint code for compressibility correction
! V1_45        2015-12-15 Harald Anlauf
!  remove obs from tail of profile only when strict fg check enabled
! V1_50        2017-01-09 Andreas Rhodin
!  option for higher order interpolation in apply_L_m
! V1_51        2017-02-24 Andreas Rhodin
!  use generalised humidity for GPSRO, namelist flag for bug-compatibility
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Andreas Rhodin  DWD  2006-2007  original code
! Detlef Pingel   DWD  2006-2007  modifications
!==============================================================================
  !=============
  ! modules used
  !=============
  !------------------------
  ! general purpose modules
  !------------------------
  use mo_kind,      only: wp, i8          ! working precision
  use mo_exception, only: finish,        &! abort routine
                          message         ! print warning message
  use mo_mpi_dace,  only: dace
#ifdef ROPP
  use mo_mpi_dace,  only: p_send,        &! generic MPI send    routine
                          p_recv,        &! generic MPI receive routine
                          p_sum           ! generic MPI sum     routine
  use mo_version,   only: var3d_version, &! 3dvar version number
                          var3d_date      ! 3dvar modification date
  use mo_run_params,only: run_time,      &! processing time
                          ana_time,      &! analysis time
                          fc_hours,      &! forecast time (hours)
                          path_file,     &! concatenate path+basename
                          aux             ! path for diagnostic output
  use mo_time,      only: iyyyy,imm,idd, &! conversion: year, month, day
                          ihh, imi        !             hour, minute
#endif
  !----------------------
  ! observation data type
  !----------------------
  use mo_obs_set,    only: t_obs_set,    &! observation data type set, global
                           t_obs_block    ! obs data type
  use mo_t_obs,      only: t_obs,        &! observation data type
                           t_spot,       &! observation meta data type
                           t_imcol,      &! interpol. space model col. meta data
                           new_int,      &! reserve memory for interp.space
                           shrink_report,&! remove passive observations
                           GPSRO,        &! GNSS occultation type id
                           tsk_name,     &! name of task
                           TSK_INIT, TSK_READ, TSK_SET_CHR, TSK_SETUP_COLS,&
                           TSK_SETUP_FULL, TSK_SETUP_FUL0, TSK_R,          &
                           TSK_Y, TSK_K, TSK_SHRINK
! use mo_t_use,      only: STAT_OBS_ONLY  ! status flag: no model equivalent
#ifdef ROPP
  use mo_fdbk_tables,only: OT_GPSRO       ! report type value for GPSRO
#endif
  !---------------------------------------
  ! atmospheric state data type definition
  !---------------------------------------
  use mo_atm_state,  only: t_atm              ! atmospheric state data type
  use mo_atm_grid,   only: t_grid,           &! atmospheric grid  data type
                           VCT_P_ISO          ! isobaric coordinate
  use mo_t_col,      only: t_cols,           &! model columns data type
                           COL_T,            &! required variables: temperature
                           COL_Q,            &! specific humidity
                           COL_X,            &! water + ice load
#ifdef DEBUG_OCC
                           COL_GEOH,         &! geopotential (half levels)
#endif
                           COL_PH             ! pressure     (half levels)
  use mo_wmo_tables, only: WMO6_LATLON,      &!
                           WMO6_GAUSSIAN,    &!
                           DWD6_ICOSAHEDRON, &!
                           DWD6_ICON          !
  use mo_physics,    only: gacc,             &! gravity acceleration
!                          Rgas=>R,          &! gas constant
                           rh_q,             &! relative <- specific humid.
                           q_rh_adj           ! adjoint routine
  use mo_cntrlvar,   only: gh_rh,            &! generalised humidity from rh
                           rh_gh              ! relative humidity    from gh
  !-----------------
  ! matrix data type
  !-----------------
  use mo_dec_matrix, only: t_vector,         &! vector data type
                           t_vector_segm      ! segment data type
  !--------------------
  ! raytracing operator
  !--------------------
  use mo_occ,  only: process_occ_3d,  &! 3d-operator to ..
                     operator,        &! namelist parameter: 1=1d,3=3d,2=hybrid
!                    skip_invalid,    &! skip invalid rays in subsequent iter.
                     skip_lowest,     &! skip lowest ray
                     use_waterload,   &! Account for condensates (water/ice)
                     chk_1d_strict,   &! Strict checks on FG in 1d operator
                     z_oro_min,       &! Min. height above surface, 1d op. [km]
                     dz_duct,         &! Min. distance above ducting layer [km]
                     dxdz_min,        &! Lower bound on dx/dz, 1d operator
                     t_occ, t_ray,    &! GNNS occultation data type
                     load_occ,        &! load  t_occ from t_obs
                     store_occ,       &! store t_occ in   t_obs
                     invalid,         &! invalid value
                     use_gh,          &! use generalised humidity
                     jac_bugfix,      &! use proper state dependence of Jacobian
                     verbose,         &! 0:silent; >=3 verbose in 1d-operator
                     write_profile,   &! write N,t,q profile (1=1d, 3=3d)
                     horint_mode,     &! Horizontal interpolation mode
                     const_occ_data,  &! test: atmosphere to homogeneous state
                     restore_occ_data,&! test: restore atmosphere
                     refract_model     ! Refractivity model for forward operator
#ifdef ROPP
  use mo_occ,  only: set_occ,         &! set occ atmospheric data structures
                     destruct_occ,    &! deallocate occ atmospheric data
                     p_send, p_recv,  &! MPI send, recv type t_occ
                     institution,     &! NetCDF attribute
                     model,           &! NetCDF attribute
                     ztop_f,          &! upper bound feedback info
                     dz_f              ! increment   feedback info
#endif
  !---------------------------------------------------
  ! Interface to Michael Gorbunovs Raytracing operator
  !---------------------------------------------------
! use Earth,            only: geodetic       ! Geodetic coordinates data type
! use Coordinates,      only: cartesian      ! Cartesian coordinates data type
  use ICO_grid,         only: gf             !
  use ECHAM_fields,     only: ng1            !
!                             ECHAM_cleanup,&! deallocate module variables
!                             ECHAM_init     ! initialize global fields
  use mo_grid_intpol,   only: ECHAM_1d_idx_init, &
                              Grid_Indices,      &
                              add_index,         &
                              alloc_imcol         ! (re)allocate imcol
  use ECHAM_fields_adj, only: ECHAM_init_adj    ,&!
                              ECHAM_cleanup_adj ,&!
                              Interpolate_Refractivity_adj
  use ECHAM_1DVar,      only: ECHAM_Refraction_1D ! 1D-bending-angle-operator
  use Earth,            only: g_ave,             &! gravity acceleration
                              Rd                  ! gas constant
  use Occ_Meteoprofiles,only: Eps                 ! Rnu/Rd - 1
  use CIPM_2007,        only: Zeta_adj            ! Compressibility correction
#ifdef ROPP
  !------------------------------------
  ! use ROPP data type and I/O routines
  !------------------------------------
  use ropp_io_dwd,      only: ROprof,            &! RO profile data type
                              ropp_io_write,     &! write RO profile data type
                              ropp_io_init,      &! init  RO profile data type
                              ropp_io_free,      &! free  RO profile data type
                              ropp_io_occid,     &! derive occultation id
                              n_ropic_2,         &! no.levels for level 2  data
                              impact_heights,    &! ROPIC impact parameter set
                              assignment(=)       ! ROPP time = 3DVAR time
                                                  ! Product confidence data:
  use ropp_io,          only: PCD_nominal,       &!   nominal   / non-nominal
                              PCD_NRT,           &!   NRT       / off line
                              PCD_setting,       &!   Setting   / Rising
                              PCD_phase_nominal, &!   nominal   / non-nominal
                              PCD_bangle_nominal,&!   nominal   / non-nominal
                              PCD_refrac_nominal,&!   nominal   / non-nominal
                              PCD_met_nominal,   &!   nominal   / non-nominal
                              PCD_bg_nominal,    &!   nominal   / non-nominal
                              PCD_occultation,   &!   retreived / background
                              PCD_missing         !   valid     / invalid
#endif
  implicit none
!==============================================================================
  !================
  ! public entities
  !================
  private
  public :: process_occ      ! routine to process GNSS occultation data
  public :: write_gpsro_prof ! write profiles (t, q, n)

!==============================================================================
contains
!==============================================================================
  subroutine process_occ (task, spot, obs, atm, cols, xi, y, Jo, Jo_atm, state)
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

  !===========================================================
  ! generic routines for 1d and 3d radio occultation operator:
  !   1d operator is implemented here
  !   3d operator is called from mo_occ
  !===========================================================

    !----------------
    ! local variables
    !----------------
    type(t_occ)           :: o         ! occultation data type
    type(t_ray)  ,pointer :: rs(:)     ! rays data type
    type(t_ray)  ,pointer :: r         ! ray data type pointer
    integer               :: tsk       ! task, local copy
    logical               :: adjoint   ! flag for TSK_K set
    type(t_grid) ,pointer :: grid      ! model grid
    integer               :: i,j,k,l,m ! indices
    integer               :: n, nn     ! indices
    integer  ,allocatable :: stat(:)   ! error status
    real(wp) ,allocatable :: em  (:)   ! model bending angle
    real(wp) ,allocatable :: zp  (:)   ! geometric height for impact p.
!   real(wp)              :: Jo_o      ! d Jo / d y
    real(wp)              :: Jk        ! J contribution from level
    real(wp) ,allocatable :: rh  (:), q  (:), t  (:), p  (:) ! variables
    real(wp) ,allocatable :: rh_a(:), q_a(:), t_a(:), p_a(:) ! adjoints
    real(wp) ,allocatable :: gh_a(:), gh (:)                 ! gen.humidity
    real(wp) ,allocatable :: qx  (:)   ! water load (condensates)
    real(wp)              :: z_a       !
    integer               :: ke        ! number of model levels
!   integer               :: id        ! diamond index
    integer               :: ix        ! column index
    real(wp) ,allocatable :: Hnew(:,:) ! Jacobi matrix (part of)
    logical  ,pointer     :: msk(:)    ! true for observations kept
    logical               :: change    ! true for observations changed
    real(wp)              :: tvirt     ! virtual potential temperature
    integer               :: nlev, kp, ngp   ! dimensions
    integer               :: np
    integer  ,allocatable :: idx(:,:)      ! Subgrid indices [point, index]
    real(wp) ,allocatable :: w  (:)        ! horizontal interpolation weight
    real(wp) ,allocatable :: e_t (:,:,:)   ! d(EM(IP))/d(T(IGP,k))
    real(wp) ,allocatable :: e_q (:,:,:)   ! d(EM(IP))/d(Q(IGP,k))
    real(wp) ,allocatable :: e_p (:,:)     ! d(EM(IP))/d(Psur(IGP))
    real(wp) ,allocatable :: em_t(:,:,:)   ! d(EM(IP))/d(T(IGP,k))
    real(wp) ,allocatable :: em_q(:,:,:)   ! d(EM(IP))/d(Q(IGP,k))
    real(wp) ,allocatable :: em_p(:,:)     ! d(EM(IP))/d(Psur(IGP))
    integer               :: nimcol        ! elements of imcol in use
    integer(i8)           :: iatm          ! model columns required
    integer               :: natm          ! # of model columns required
    logical               :: strict        ! Strict checks in forward operator?
    logical               :: checked       ! Profile already fed thru FG check
    real(wp)              :: z_oro_tmp     ! Min. height above surface, 1d op.
    real(wp)              :: dxdz_tmp      ! Lower bound on dx/dz, 1d operator
    real(wp)              :: dz_duct_tmp   ! Min. distance above ducting layer
    logical               :: nneighb       ! nearest neighbour mode
    integer               :: order         ! interpolation order
    integer               :: ixmin, ixmax  ! min/max of indices to model cols.
    integer, allocatable  :: ir(:)         ! ray index list
    integer               :: nc            ! no. rays in chunk
    integer               :: nr            ! no. rays (check)
    real(wp),allocatable  :: em_  (:)      ! temporary model bending angle
    real(wp),allocatable  :: zp_  (:)      ! temporary geom. height for impact p.
    integer ,allocatable  :: stat_(:)      ! temporary error status
    logical ,allocatable  :: limcol(:)     ! Mask for tracking model columns
    integer               :: iii, iio      ! Auxiliary indices
    integer               :: nip           ! No. input parameter / model column
    type(t_imcol),pointer :: imcol_(:)     ! temporary imcol
    integer ,allocatable  :: jx    (:)     ! auxiliary map for column indices
    real(wp)              :: Zeta          ! Compressibility correction
    real(wp)              :: Zeta_p
    real(wp)              :: Zeta_t
    real(wp)              :: Zeta_q

    select case (operator)
    case (1)
      ! processed below, 1d operator (Abel transform)
    case (2)
      ! processed below, multiple 1d operator (ar perigee points)
    case (3)
      call process_occ_3d (task, spot, obs, atm, cols, xi, y, Jo, Jo_atm, state)
      return
    case default
      write(0,*)  'process_occ:  operator',operator,'not yet implemented!'
      call finish('process_occ','operator not yet implemented!')
    end select

    !--------------------------------------------------------------------------
    ! Main observation processing routine
    ! similar for all observation operators
    !-------------------------------------------------------------------------

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

    !==============================
    ! observation non_specific part
    !==============================
    tsk  =  task
    grid => atm% grid
    !-----------------------------------------------
    ! TSK_INIT : Initialization
    ! read the namelist, same for 1d and 3d operator
    !-----------------------------------------------
    if (iand (TSK_INIT,tsk) /= 0) then
      call process_occ_3d (TSK_INIT, spot, obs, atm, cols, xi, y, Jo, Jo_atm, state)
      tsk=tsk-TSK_INIT
    endif
    if (tsk==0) return
    !-----------------------------------------------
    ! TSK_READ: read data
    ! read NetCDF files, same for 1d and 3d operator
    !-----------------------------------------------
    if (iand (TSK_READ,tsk) /= 0) then
      call process_occ_3d (TSK_READ, spot, obs, atm, cols, xi, y, Jo, Jo_atm, state)
      tsk=tsk-TSK_READ
    endif
    if (tsk==0) return
    !---------------------------------------------
    ! TSK_SET_CHR: set observation characteristics
    !---------------------------------------------
    if (iand (TSK_SET_CHR,tsk) /= 0) then
      tsk=tsk-TSK_SET_CHR
    endif
    if (tsk==0) return

    if (verbose >= 2 .and. dace% lpio) print *, "process_occ: ", tsk_name (task)
    !==========================
    ! observation specific part
    !==========================
    if (.not.present(spot)) call finish('process_occ_1d','spot not present')
    if (.not.present(obs )) call finish('process_occ_1d','obs  not present')
    if (spot% hd% modtype /= GPSRO) return

    !--------------------------------------------------
    ! TSK_SETUP_FUL0: setup of PSAS interpolation space
    ! same for 1d and 3d operator
    !--------------------------------------------------
    if (iand (TSK_SETUP_FUL0,tsk) /= 0) then
      if (dace% pe==obs% o% pe) then
        call new_int (obs% o, spot, (2* atm% grid% nz + 1) * size(spot% imcol))
      endif
      tsk=tsk-TSK_SETUP_FUL0
    endif
    if (tsk==0) return
    !--------------------------------------------------
    ! TSK_SETUP_FULL: setup of PSAS interpolation space
    ! same for 1d and 3d operator:
    ! set obs% o% t_int       observed quantities
    !     obs% o% lev         levels
    !---------------------------------------------------
    if (iand (TSK_SETUP_FULL,tsk) /= 0) then
      if (dace% pe==obs% o% pe) then
        call process_occ_3d &
          (TSK_SETUP_FULL, spot, obs, atm, cols, xi, y, Jo, Jo_atm, state)
      endif
      tsk=tsk-TSK_SETUP_FULL
    endif
    if (tsk==0) return

    !===============================================
    ! tsk == TSK_R
    ! set up R (observation error covariance matrix)
    ! same for 1d and 3d operator:
    !     obs% R%             observational error
    !     obs% o% s_vqc       variational quality bounds
    !===============================================
    if (iand (TSK_R,tsk) /= 0) then
      if (dace% pe==obs% o% pe) then
        call process_occ_3d &
          (TSK_R, spot, obs, atm, cols, xi, y, Jo, Jo_atm, state)
      endif
      tsk = tsk - TSK_R
      if (tsk == 0) return
    endif

    !------------------------------------------
    ! TSK_Y: run nonlinear observation operator
    ! TSK_K: in addition provide Jacobi matrix
    !------------------------------------------
    if (iand   (TSK_K+TSK_Y, tsk) /= 0) then
      adjoint = iand (TSK_K, tsk) /= 0
      if (spot% pe_eval == dace% pe) then
        !----------------------------------
        ! allocate memory for Jacobi matrix
        !----------------------------------
        if (adjoint) then
          k = obs% H% ia (spot% i% i + 1)
          do j=1,spot% i% n
            obs% H% ia (spot% i% i +j) = k
          end do
          obs% H% ia (spot% i% i + spot% i% n + 1) = k
!         obs% xi% x (spot% i% i+1:spot% i% i + spot% i% n) = 0._wp
          obs% yi% x (spot% o% i+1:spot% o% i + spot% o% n) = invalid
        endif
        !------------------------------------------------
        ! copy atmospheric state to echam_fields module
        ! Initialization of echam_constituents_adj module
        !------------------------------------------------
!       call set_occ (gf, grid, cols, obs% o, xi)
        call echam_init_adj
        call const_occ_data (gf, spot)
        call load_occ  (obs% o, spot, o, rs)
        !====================
        ! run nonlinear model
        !====================
        !--------------------------------
        ! set dimensions, allocate arrays
        !--------------------------------
        nlev = gf% nz
        kp   = size(rs)
        ngp  = gf% ncol                 ! Number of actually used column(s)
!       select case (gf% grid_type)
!       case (WMO6_GAUSSIAN, WMO6_LATLON)
!         ngp = ng1**2
!       case (DWD6_ICOSAHEDRON, DWD6_ICON)
!         ngp = 3
!       end select
        nimcol = size(spot% imcol)
        allocate (em_t(kp, nimcol, nlev)); em_t = 0._wp
        allocate (em_q(kp, nimcol, nlev)); em_q = 0._wp
        allocate (em_p(kp, nimcol      )); em_p = 0._wp
        allocate (stat (kp)); stat  = HUGE (0)
        allocate (stat_(kp)); stat_ = HUGE (0)
        allocate (ir   (kp)); ir    = HUGE (0)
        allocate (em   (kp)); em    = HUGE (0._wp)
        allocate (em_  (kp)); em_   = HUGE (0._wp)
        allocate (zp   (kp)); zp    = HUGE (0._wp)
        allocate (zp_  (kp)); zp_   = HUGE (0._wp)
        allocate (idx (    ngp,    3   ))
        nneighb = (horint_mode == 0)
        if (nneighb .or. operator == 1) then
           allocate (e_t  (kp, ngp,    nlev))
           allocate (e_q  (kp, ngp,    nlev))
           allocate (e_p  (kp, ngp         ))
        else
           allocate (e_t  (1,  ngp,    nlev))
           allocate (e_q  (1,  ngp,    nlev))
           allocate (e_p  (1,  ngp         ))
        end if

        !-----------------------------------------------------------
        ! Strict validity checks for FG, diagnostics in minimization
        !-----------------------------------------------------------
        strict = chk_1d_strict .and. (.not. o% checked .or. &
                                      .not. o% critical     )

        if (chk_1d_strict .and. (.not. o% checked)) then
           z_oro_tmp = z_oro_min    ! Min. height above surface, 1d op. [km]
           dxdz_tmp  = dxdz_min     ! Lower bound on dx/dz, 1d operator
           dz_duct_tmp = dz_duct    ! Min. distance above ducting layer [km]
        else
           z_oro_tmp = 0._wp        ! Use loose checks (later iterations)
           dxdz_tmp  = 0._wp
           if (dz_duct >= 0._wp) then
              dz_duct_tmp =  0._wp
           else
              dz_duct_tmp = -1._wp
           end if
        end if

        !---------------------------
        ! run bending-angle operator
        !---------------------------
        select case (operator)
        case (1)
          !-----------------------------
          ! 1d operator (Abel transform)
          !-----------------------------
          if (kp /= o% nray) then
             write(0,*) "process_occ: ray count mismatch:", kp, o% nray
             call finish ("process_occ","ray count mismatch")
          end if
          nc = 0
          do k = 1, o% nray
!            if (obs% o% body(spot% o% i+k)% use% state <= STAT_OBS_ONLY) cycle
             nc     = nc + 1
             ir(nc) = k
          end do
          call ECHAM_Refraction_1D &
            (rs(ir(1:nc))% p, & ! <-- RO impact parameter [km]
             o % gp1d,        & ! <-- Occultation point (geodetic)
             o % rlc,         & ! <-- Curvature radius
             verbose,         & ! <-- Verbose level
             strict,          & ! <-- Strict checks on profile validity
             z_oro_tmp,       & ! <-- Min. height above surface [km]
             dxdz_tmp,        & ! <-- Lower bound on dx/dz
             dz_duct_tmp,     & ! <-- Min. dist. above ducting layer [km]
             adjoint,         & ! <-- evaluate adjoint
             em_,             & ! --> Model refraction angle
             zp_,             & ! --> geometric height for impact parameter
             idx,             & ! --> Subgrid indices [point, index]
             e_t,             & ! --> d(EM(IP))/d(T(IGP,k))
             e_q,             & ! --> d(EM(IP))/d(Q(IGP,k))
             e_p,             & ! --> d(EM(IP))/d(Psur(IGP))
             stat_            ) ! --> Error status
          do j = 1, nc
             i = ir(j)
             em  (i)     = em_  (j)
             zp  (i)     = zp_  (j)
             stat(i)     = stat_(j)
             em_t(i,:,:) = e_t  (j,:,:)
             em_q(i,:,:) = e_q  (j,:,:)
             em_p(i,:  ) = e_p  (j,:  )
          end do
        case (2)
          !-----------------------------------------
          ! multiple 1d operator (at perigee points)
          !-----------------------------------------
          if (nneighb) then
             !--------------------------------------------------------
             ! Optimized version with nearest-neighbor mode:
             ! for each model column, run operator on all nearby rays.
             !--------------------------------------------------------
             ixmin = minval (rs(:)% ix(1))      ! note: ix(2:)==0
             ixmax = maxval (rs(:)% ix(1))      ! for nneighb=T
             nr = 0
             do k = ixmin, ixmax
                nc = 0
                do i = 1, o% nray
                   if (rs(i)% ix(1) == k) then
                      nc     = nc + 1
                      ir(nc) = i
                   end if
                end do
                if (nc == 0) then
                   write(0,*) "process_occ: model column missing:", k, ixmin, ixmax
                   call message ("process_occ","model column missing")
                   cycle
                end if
                nr = nr + nc
                call ECHAM_Refraction_1D &
                  (rs(ir(1:nc))% p,    & ! <-- RO impact parameter [km]
                   rs(ir(1))   % geo,  & ! <-- Occultation point (geodetic)
                   o           % rlc,  & ! <-- Curvature radius
                   verbose,            & ! <-- Verbose level
                   strict,             & ! <-- Strict checks on profile validity
                   z_oro_tmp,          & ! <-- Min. height above surface [km]
                   dxdz_tmp,           & ! <-- Lower bound on dx/dz
                   dz_duct_tmp,        & ! <-- Min. dist. above ducting layer [km]
                   adjoint,            & ! <-- evaluate adjoint
                   em_(1:nc),          & ! --> Model refraction angle
                   zp_(1:nc),          & ! --> geometric height for impact parameter
                   idx,                & ! --> Subgrid indices [point, index]
                   e_t,                & ! --> d(EM(IP))/d(T(IGP,k))
                   e_q,                & ! --> d(EM(IP))/d(Q(IGP,k))
                   e_p,                & ! --> d(EM(IP))/d(Psur(IGP))
                   stat_(1:nc)         ) ! --> Error status
                do j = 1, nc
                   i = ir(j)
                   em  (i)                  = em_  (j)
                   zp  (i)                  = zp_  (j)
                   stat(i)                  = stat_(j)
                   em_t(i,rs(i)%ix(:ngp),:) = em_t(i,rs(i)%ix(:ngp),:) + e_t(j,:,:)
                   em_q(i,rs(i)%ix(:ngp),:) = em_q(i,rs(i)%ix(:ngp),:) + e_q(j,:,:)
                   em_p(i,rs(i)%ix(:ngp)  ) = em_p(i,rs(i)%ix(:ngp)  ) + e_p(j,:  )
                end do
             end do
             if (nr /= o% nray) then
                write(0,*) "process_occ: ray count mismatch:", nr, o% nray
                call finish ("process_occ","ray count mismatch")
             end if
          else
            !------------------------------------------
            ! Run forward model separately on each ray:
            !------------------------------------------
            do k = 1, o% nray
              call ECHAM_Refraction_1D &
                (rs(k:k)% p,     & ! <-- RO impact parameter [km]
                 rs(k)  % geo,   & ! <-- Occultation point (geodetic)
                 o      % rlc,   & ! <-- Curvature radius
                 verbose,        & ! <-- Verbose level
                 strict,         & ! <-- Strict checks on profile validity
                 z_oro_tmp,      & ! <-- Min. height above surface [km]
                 dxdz_tmp,       & ! <-- Lower bound on dx/dz
                 dz_duct_tmp,    & ! <-- Min. dist. above ducting layer [km]
                 adjoint,        & ! <-- evaluate adjoint
                 em(k:k),        & ! --> Model refraction angle
                 zp(k:k),        & ! --> geometric height for impact parameter
                 idx,            & ! --> Subgrid indices [point, index]
                 e_t,            & ! --> d(EM(IP))/d(T(IGP,k))
                 e_q,            & ! --> d(EM(IP))/d(Q(IGP,k))
                 e_p,            & ! --> d(EM(IP))/d(Psur(IGP))
                 stat(k:k)       ) ! --> Error status
              em_t(k,rs(k)%ix(:ngp),:) = em_t(k,rs(k)%ix(:ngp),:) + e_t(1,:,:)
              em_q(k,rs(k)%ix(:ngp),:) = em_q(k,rs(k)%ix(:ngp),:) + e_q(1,:,:)
              em_p(k,rs(k)%ix(:ngp)  ) = em_p(k,rs(k)%ix(:ngp)  ) + e_p(1,:  )
            end do
          end if
        case default
          write(0,*)  'process_occ(TSK_K+TSK_Y):  operator',&
                       operator,'not implemented!'
          call finish('process_occ(TSK_K+TSK_Y)','operator not implemented!')
        end select

        checked = o% checked
        !-------------------------------
        ! loop over observational levels
        !-------------------------------
        do k = 1, o% nray
          r => rs(k)
!         if (obs% o% body(spot% o% i+k)% use% state <= STAT_OBS_ONLY) then
!           if (present(y)) y% x (spot% o% i+k) = invalid
!           cycle
!         end if
          !-------------------------------------------------------
          ! remember validity of observation operator applications
          ! for later decisions in line search and minimisation
          !-------------------------------------------------------
          obs% o% body(spot% o% i+k)% op_na = stat(k)

          if(r% eps == invalid) then
          !-------------
          ! ray not used
          !-------------
            if (verbose > 1)                                                  &
            write(0,                                                          &
            '("process_occ_1d: pe=",i3,", ray=",i3,", z=",f7.3,", skipped")') &
            dace% pe,k,r% p-o% rlc
          else
          !---------
          ! ray used
          !---------
            if (verbose > 2 .or. (stat(k) == 1 .and. verbose > 0)) then
             if (stat(k)/=0) then
              write(0,                                                           &
               '("process_occ_1d: pe=",i3,", ray=",i3,", z=",f7.3,", stat=",i1)')&
               dace% pe,k,r% p-o% rlc,stat(k)
             else
              write(0,                                                           &
                '("process_occ_1d: pe=",i3,", ray=",i3,", z=",f7.3,", stat=",i1, &
                &", Em-o=",3f11.8)') dace% pe,k,r% p-o% rlc,stat(k),             &
                em(k) - r % eps, sqrt(r % var2), em(k)
             endif
            endif
            !---------------------------------------------------------
            ! mark invalid rays to be skipped in subsequent iterations
            !---------------------------------------------------------
            if(o% checked) then
              if(stat(k)/=0) then
                ! stat=1: bending angle could not be evaluated
                ! stat=2: b.a. unreliable (numerically unstable)
                ! stat=3: b.a. extrapolated/below minimum height
                if (verbose > 0) then
                  print *,'process_occ_1d: invalid ray for ', trim (spot% statid)
                  write(0,'(A,A,A,2F9.3,I4,3F11.6,2I3)') &
                       "process_occ_1d: invalid ray for ", spot% statid, ":", &
                       o%gp1d% phi, o%gp1d% lambda, k, r% p-o% rlc, r% eps,   &
                       em(k), stat(k), obs% o% body(spot% o% i+k)% use% state
                end if
!               call finish ('process_occ_1d','invalid ray')
                if (verbose > 1)                                      &
                  write(0,'("process_occ_1d: pe=",i3,", ray=",i3,", z=",f7.3, &
                           &", SKIP RAY !!!")') dace% pe,k,r% p-o% rlc
! Disabled because of destructive interference with FG checks...
!               if (skip_invalid) r% eps = invalid
                !---------------------------------------------------
                ! After FG check, diagnose bad rays once per profile
                !---------------------------------------------------
                if (checked)      o% critical = .true.
              endif
            else
              if(stat(k)/=0) then
                if (strict)      r% eps = invalid ! remove from tail of profile
              else
                if (skip_lowest) r% eps = invalid ! skip lowest ray
                o% checked = .true.
              endif
            endif

            !------------------------------------------
            ! return model bending angle, cost function
            !------------------------------------------
            if (present(y)) y% x (spot% o% i+k) = invalid
            if(r% eps /= invalid .and. stat(k)/=1) then
              Jk   = 0.5_wp * (em(k) - r % eps)**2 / r % var
!             Jo_o =          (em(k) - r % eps)    / r % var
              if (present(Jo))  Jo = Jo + Jk
              if (present(y))   y%  x (spot% o% i+k) = em(k)
              if (adjoint) obs% yi% x (spot% o% i+k) = em(k)
            endif
            !-------------------+++++++++++++
            ! calculate gradient (obsolete ?)
            !-------------------+++++++++++++
          end if
          if(r% eps /= invalid) rs(k)% geo% h = zp(k)
        end do

        !--------
        ! cleanup
        !--------
        call restore_occ_data (gf)
!       call destruct_occ
        call echam_cleanup_adj
!       call echam_cleanup

        deallocate (stat, stat_)
        deallocate (em)
        deallocate (zp)
        call store_occ (obs% o, spot, o, rs)
        deallocate (rs)

        !---------------
        ! setup K matrix
        !---------------
        if (adjoint) then
          ke = cols% ke
          allocate (gh_a(ke), rh_a(ke), q_a(ke), t_a(ke), p_a(ke))
          allocate (gh  (ke), rh  (ke), q  (ke), t  (ke), p  (ke), qx(ke))
          qx(:) = 0._wp
!         obs% xi% x (spot% i% i + 1 : spot% i% i + spot% i% n) = 0._wp
          k = obs% H% ia (spot% i% i + 1)
          do j=1,spot% i% n
            obs% H% ia (spot% i% i +j) = k
          end do
          obs% H% ia (spot% i% i + spot% i% n + 1) = k
          m = 0
          do l=1,size(spot% imcol)
            ix = spot% imcol(l)% imc(1)
!           i  = cols% col(ix)% i
!           j  = cols% col(ix)% j
!           id = cols% col(ix)% l
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
            if (use_waterload) qx = cols% col(ix)% x
!           z  =     cols% col(ix)% s% geosp / gacc
            rh = rh_q (q, t, p)
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
            !-------------------------------------------------
            ! use same hydrostatic relationship as occ-library
            ! for consistency of nonlinear/adjoint
            ! (including compressibility correction)
            !-------------------------------------------------
            if (refract_model == 1) then
               Zeta_p = 0._wp
            else
               call Zeta_adj (cols% col(ix)% s% ps, t(ke), q(ke), &
                              Zeta, Zeta_p, Zeta_t, Zeta_q        )
            end if
            tvirt = t(ke) * (1 + Eps * q(ke) - qx(ke))
!           z_a = cols% col(ix)% s% ps * gacc / (cols% col(ix)% t (ke) * Rgas)
!           z_a = cols% col(ix)% s% ps * g_ave / (Rd * tvirt)
            z_a = g_ave / Rd  * cols% col(ix)% s% ps /        &
                  (tvirt * (1 + cols% col(ix)% s% ps * Zeta_p))
            !---------------
            ! Jacobi matrix:
            !---------------
            allocate (Hnew(spot% o% n, nn       ))
            Hnew=0._wp
            do k=1,o% nray
              Hnew (k,1:2*ke:2) = em_t (k, l, :) +   &
                                  em_q (k, l, :) * t_a
!do i=1,ke
!print *,'d em/d t',k,l,i,em_t (k, l, i)
!end do
!print *
              Hnew (k,2:2*ke:2) = em_q (k, l, :) * rh_a
              Hnew (k,  nn    ) = em_p (k, l   ) * z_a
            end do
!Hnew=0._wp
!Hnew(:,2*ke-1) = 1
            k = obs% H% ia (spot% i% i + n + 1)
            do j=1,nn
              obs% H% ia (spot% i% i + n + j) = k
!write(6,'(a,2i3,(10f12.8))')'Hnew',l,j,Hnew (:,j)
!NEC$ ivdep
              do i=1,spot% o% n
                obs% H% packed (k) = Hnew (i,j)
                obs% H% ja     (k) = spot% o% i + i
                k = k + 1
              end do
            end do
            obs% H% ia (spot% i% i + n + nn + 1) = k
            deallocate (Hnew)
            !------------------
            ! final: set t,rh,z
            !------------------
!           obs% xi% x (spot% i% i+1:spot% i% i+2*ke:2) = t    ! temperature
!           obs% xi% x (spot% i% i+2:spot% i% i+2*ke:2) = rh   ! relative hum.
!           obs% xi% x (             spot% i% i+nn    ) = z    ! geopotential
          end do
          deallocate (gh_a, rh_a, q_a, t_a, p_a, gh, rh, q, t, p, qx)
        endif

        deallocate (idx)
        deallocate (em_t, em_q, em_p)
        deallocate (e_t, e_q, e_p)
        deallocate (em_, zp_)

      endif
      if (iand (TSK_K,tsk) /= 0) tsk=tsk-TSK_K
      if (iand (TSK_Y,tsk) /= 0) tsk=tsk-TSK_Y
    endif
    if (tsk==0) return

    !-------------------------------------------------
    ! TSK_SETUP_COLS:
    ! determine model columns required by the operator
    !-------------------------------------------------
    if (iand (TSK_SETUP_COLS,tsk) /= 0) then
      select case (grid% vct)
      case (VCT_P_ISO)
         call finish ("process_occ","isobaric model grid not supported")
      end select

      call load_occ (obs% o, spot, o, rs)
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
      if (associated (spot% imcol)) then
         deallocate (spot% imcol)
         nullify (spot% imcol)
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
      spot% n_spt = size (spot% imcol)
      spot% mke   = grid% nz
      select case (operator)
      case (1)
        ! processed above, 1d operator (Abel transform)
      case (2)
        !--------------------------------------
        ! true columns for multiple 1d-operator
        ! (need not contain mean tangent point)
        !--------------------------------------
        deallocate (spot% imcol)
        allocate   (spot% imcol(0))

        select case (grid% gridtype)
        case (WMO6_GAUSSIAN, WMO6_LATLON)
          ngp = ng1**2
        Case (DWD6_ICOSAHEDRON, DWD6_ICON)
          ngp = 3
        end select
        allocate (idx(ngp,4))
        allocate (w  (ngp  ))

        order  = 2; if (horint_mode == 0) order = 1
        nimcol = size (spot% imcol)
        do i = 1, size (rs)

          call Grid_Indices &
            (rs(i)%geo% lambda, & ! <-- geodetic longitude
             rs(i)%geo% phi,    & ! <-- geodetic latitude
             grid,              & ! <-- grid data type
             IDX,               & ! --> Grid point indices [Point, index]
             w,                 & ! --> Weight
             np,                & ! --> number of points returned
             order = order      ) ! <-- nearest neighbor mode

          call add_index (idx,         &! <-- Grid point indices [Point, index]
                          np,          &! <-- number of neigbour points
                          obs% o,      &! <-> observation data type
                          iatm,        &! <-- atmospheric parameters
                          natm,        &! <-- no. atmospheric parameters
                          spot% imcol, &! <-> model columns required
                          nimcol,      &! <-> elements of imcol in use
                          rs(i)% ix,   &! --> indices to imcol
                          grid,        &! <-- grid data type
                          spot% i_time,&! <-- time slot
                          spot% w_time )! <-- time interpolation weight
        end do
        call alloc_imcol (spot% imcol, nimcol)
        spot% n_spt = size (spot% imcol)
        deallocate (idx, w)
      case default
        call finish('process_occ','operator not yet implemented!')
      end select

      tsk=tsk-TSK_SETUP_COLS
      call store_occ (obs% o, spot, o, rs)
      deallocate (rs)
    endif
    if(tsk==0) return

    !------------------------------------------------------------------
    ! TSK_SHRINK:
    ! release unused observations in the report
    !------------------------------------------------------------------
    if (iand (TSK_SHRINK,tsk) /= 0) then
      call shrink_report (spot, obs%o, state, change, mask=msk)
      if (change) then
        call load_occ (obs% o, spot, o, rs)
        nip = 2 * spot% mke + 1                 ! no. input params./model column
        allocate (limcol(spot% i% n / nip))
        limcol(:) = .false.
        if (operator == 1) limcol(:) = .true.   ! Pure 1d op. retains all cols.
        j = 0
        do i=1,size(msk)
          if (msk(i)) then
            j = j + 1
            if (j < i) rs (j) = rs (i)
            !print *, "i=",i," rs(i)%ix =",rs(i)%ix
            do k = 1, size (rs(j)% ix)
               if (rs(j)% ix(k) == 0) exit
               limcol(rs(j)% ix(k)) = .true.    ! Remember required columns
            end do
          endif
        end do
        o%         nray = j
        spot% col% nlev = j
        n = count (limcol)
        if (n < size (limcol)) then
           !-------------------------------
           ! Shrink (pack) input parameters
           !-------------------------------
           allocate (jx(size (limcol)))
           jx = 0
           j  = 0
           do i = 1, size (limcol)
              if (limcol(i)) then
                 j     = j + 1
                 jx(i) = j      ! map from old to new imcols
                 if (j < i) then
                    iio = spot%i%i + nip*(j-1) ! to
                    iii = spot%i%i + nip*(i-1) ! from
                    obs% o% t_int(iio+1:iio+nip) = obs% o% t_int(iii+1:iii+nip)
                    obs% o% lev  (iio+1:iio+nip) = obs% o% lev  (iii+1:iii+nip)
                    obs% o% bgeri(iio+1:iio+nip) = obs% o% bgeri(iii+1:iii+nip)
                    if (associated (obs% o% bgi))                             &
                      obs% o% bgi(iio+1:iio+nip) = obs% o% bgi  (iii+1:iii+nip)
                 end if
              end if
           end do
           spot% i% n = j*nip   ! new no. of input parameters
           !-------------------
           ! Shrink spot%imcol?
           !-------------------
           if (associated (spot% imcol)) then
              imcol_ => spot% imcol
              allocate (spot% imcol(n))
              spot% imcol(:) = pack (imcol_, limcol)
              deallocate (imcol_)
              spot% n_spt = size (spot% imcol)
              !----------------------
              ! Update imcol indices:
              !----------------------
              do i = 1, o% nray
                 do k = 1, size (rs(i)% ix)
                    if (rs(i)% ix(k) == 0) exit
                    rs(i)% ix(k) = jx(rs(i)% ix(k))
                 end do
              end do
           end if
           deallocate (jx)
        end if
        call store_occ (obs% o, spot, o, rs(1:o% nray))
        deallocate (limcol)
        deallocate (msk)
        deallocate (rs)
      endif
      tsk = tsk - TSK_SHRINK
      if (tsk == 0) return
    endif

    !==========================
    ! abort if any task is left
    !==========================
!   if (iand (TSK_K,tsk) /= 0) tsk = tsk - TSK_K

    if (tsk /= 0) then
      write(0,*)  'process_occ:  unknown task =',tsk
      call finish('process_occ','unknown task')
    endif

  end subroutine process_occ
!==============================================================================
#if 0
  Subroutine Refraction_1D &
   (PRO,       & ! <-- RO impact parameter [km]
    GC,        & ! <-- Occultation point (geodetic)
    RC,        & ! <-- Curvature radius
    Vrb,       & ! <-- Verbose level
    strict,    & ! <-- Strict checks on profile validity
    EM,        & ! --> Model refraction angle
    Stat       ) ! --> Error status
  !--------------------------------
  ! wrapper for ECHAM_Refraction_1D
  !--------------------------------
  real(wp)       ,intent(in)  :: pro(:)      ! impact parameters [km]
  type(geodetic) ,intent(in)  :: gc          ! occultation point (geodetic)
  real(wp)       ,intent(in)  :: rc          ! Curvature radius
  integer        ,intent(in)  :: Vrb         ! Verbose level
  logical       , intent(in)  :: Strict      ! Strict checks on profile validity
  real(wp)       ,intent(out) :: em(:)       ! Model refraction angles
  integer        ,intent(out) :: stat(:)     ! Error status

    integer  ,allocatable :: idx(:,:)      ! Subgrid indices [point, index]
    real(wp) ,allocatable :: em_t(:,:,:)   ! d(EM(IP))/d(T(IGP,k))
    real(wp) ,allocatable :: em_q(:,:,:)   ! d(EM(IP))/d(Q(IGP,k))
    real(wp) ,allocatable :: em_p(:,:)     ! d(EM(IP))/d(Psur(IGP))

    integer :: nlev, kp, ngp
    !-----------
    ! dimensions
    !-----------
    nlev = gf% nz
    kp   = size(pro)
    select case (gf% grid_type)
    case (WMO6_GAUSSIAN, WMO6_LATLON)
      ngp = ng1**2
    case (DWD6_ICOSAHEDRON)
      ngp = 3
    end select
    !----------------------
    ! allocate local arrays
    !----------------------
    allocate (idx (    ngp, 3   ))
    allocate (em_t(kp, ngp, nlev))
    allocate (em_q(kp, ngp, nlev))
    allocate (em_p(kp, ngp      ))
    call ECHAM_Refraction_1D &
      (PRO,       & ! <-- RO impact parameter [km]
       GC,        & ! <-- Occultation point (geodetic)
       RC,        & ! <-- Curvature radius
       Vrb,       & ! <-- Verbose level
       Strict,    & ! <-- Strict checks on profile validity
       EM,        & ! --> Model refraction angle
       IDX,       & ! --> Subgrid indices [point, index]
       EM_T,      & ! --> d(EM(IP))/d(T(IGP,k))
       EM_Q,      & ! --> d(EM(IP))/d(Q(IGP,k))
       EM_P,      & ! --> d(EM(IP))/d(Psur(IGP))
       Stat)        ! --> Error status
    !------------------------
    ! deallocate local arrays
    !------------------------
    deallocate (idx)
    deallocate (em_t)
    deallocate (em_q)
    deallocate (em_p)
  end Subroutine Refraction_1D
#endif
!==============================================================================
  subroutine write_gpsro_prof (obs, grid, cols, bg_o, bg_o_e, step)
  type(t_obs_set) ,intent(in) :: obs     ! observations
  type(t_grid)    ,intent(in) :: grid    ! atmospheric grid description
  type(t_cols)    ,intent(in) :: cols(:) ! model columns
  type(t_vector)  ,intent(in) :: bg_o    ! background in observation space
  type(t_vector)  ,intent(in) :: bg_o_e  ! background error in obs.  space
  integer         ,intent(in) :: step    ! 1:monitoring, 2:analysis
  !---------------------------------------------
  ! write diagnostics for GPS radio occultations
  !---------------------------------------------

#ifndef ROPP
    !----------------------------------------------------
    ! no diagnostic output if ROPP library is not linked.
    ! give warning instead.
    !----------------------------------------------------
    if (iand (write_profile, step) == 0) return
    if (dace% lpio) &
    call message ('write_gpsro_prof',&
                  'libropp_io not linked, no GPSRO diagnostics written')
#else

    !================
    ! local variables
    !================
    logical                     :: append
    type(ROprof)                :: rp        ! radio occultation profile
    character(len=11)           :: filename  ! 'monGPSRP.nc' or 'monGPSRA.nc'
    !--------------------------------
    ! acess to 3dvar observation data
    !--------------------------------
    integer                 :: ib            ! box index
    integer                 :: is            ! occultation index
    integer                 :: i,j           ! level indices
    type(t_obs)  ,pointer   :: o             ! box pointer
    type(t_spot) ,pointer   :: s             ! spot pointer
    !-------------------
    ! indices and bounds
    !-------------------
    integer                 :: n1b           ! number of level 1 data
    integer                 :: ke            ! number of model levels
    integer                 :: iocc          ! occultation number in the box
    integer                 :: nocc          ! number of occultations per box
    integer                 :: mocc          ! total number of occultations
    !-----------------------------------------------
    ! parameters passed to the refractivity operator
    !-----------------------------------------------
    integer                 :: nz
    real(wp)                :: dz
    real(wp)   ,allocatable :: z(:)          ! geometric height (km)
    integer                 :: iz            ! level index
    real(wp)   ,allocatable :: ng  (:,:)     ! dN/d(alt,lat,lon)
    real(wp)                :: nh  (3,3)     ! hessian matrix of N
    integer    ,allocatable :: idx (:,:)     ! subgrid indices [point, index]
    real(wp)   ,allocatable :: np_t  (:,:)   ! d(np)/d(t(igp,k))
    real(wp)   ,allocatable :: np_q  (:,:)   ! d(np)/d(q(igp,k))
    real(wp)   ,allocatable :: np_p  (:)     ! d(np)/d(psur(igp))
    real(wp)   ,allocatable :: ng_t(:,:,:)   ! d(ng(i))/d(t(igp,k))
    real(wp)   ,allocatable :: ng_q(:,:,:)   ! d(ng(i))/d(q(IGP,k))
    real(wp)   ,allocatable :: ng_p(:,:)     ! d(ng(i))/d(psur(igp))
    real(wp)                :: zmin, zmax    ! bounds for valid refractivity
    integer                 :: ngp           ! subgrid size
    !----------------------------------------------------
    ! arrays for communication between processor elements
    !----------------------------------------------------
    type(t_ray)    ,pointer :: rs(:)
    type(t_occ)             :: occ
    type(t_occ),allocatable :: occs(:)
    real(wp)   ,allocatable :: lat (:,:)     ! latitude
    real(wp)   ,allocatable :: lon (:,:)     ! longitude
    real(wp)   ,allocatable :: np  (:,:)     ! refractivity
    real(wp)   ,allocatable :: t   (:,:)     ! temperature
    real(wp)   ,allocatable :: q   (:,:)     ! specific humidity
    real(wp)   ,allocatable :: gp  (:,:)     ! geopotential
    real(wp)   ,allocatable :: pf  (:,:)     ! Pressure
    real(wp)   ,allocatable :: ps  (:)       ! surface pressure
    real(wp)   ,allocatable :: hs  (:)       ! surface geopotential height
    real(wp)   ,allocatable :: gu  (:)       ! geoid undulation
    real(wp)   ,allocatable :: p   (:,:)     ! impact parameter
    real(wp)       ,pointer :: pr  (:)       ! ROPIC impact parameter set
    real(wp)   ,allocatable :: em  (:,:)     ! bending angle
    real(wp)   ,allocatable :: em_e(:,:)     ! bending angle error (stdev)

    !-----------------------------------
    ! check for diagnostic printout flag
    !-----------------------------------
    if (iand (write_profile, step) == 0) return

    !-----------------------
    ! count no. occultations
    !-----------------------
    mocc = 0
    do ib = 1, size(obs% o)
      o => obs% o(ib)
      if (o% pe == dace% pe) &
        mocc = mocc + count (o% spot(1:o% n_spot)% hd% obstype == OT_GPSRO)
    end do
    mocc = p_sum(mocc)
    if (mocc == 0) return

    !-------------------------------
    ! get ROPIC impact parameter set
    !-------------------------------
    call impact_heights (pr)
    n1b = size (pr)
    dz  = 0.2_wp
    nz  = ztop_f/dz + 1
    allocate (z(nz))
    allocate (ng(3,nz))
    !-------------------------------------------------------------------------
    ! derive profiles (diagnostics) of atmospheric parameters at perigee point
    !-------------------------------------------------------------------------
    append = .false.
    do ib = 1, size(obs% o)
      o => obs% o(ib)
      if (o% pe == dace% pe) then
        nocc = count (o% spot(1:o% n_spot)% hd% obstype == OT_GPSRO)
        if (nocc > 0) then

          !------------------------------------------------
          ! copy atmospheric state to echam_fields module
          !------------------------------------------------
          call set_occ (gf, grid, cols(ib))
          !-------------------------------------------------
          ! allocate arrays passed to interpolation routines
          !-------------------------------------------------
          ngp  = gf   % ngp
          ke   = grid % nz
          allocate (idx    (ngp, 3))
          allocate (np_t   (ngp,ke))
          allocate (np_q   (ngp,ke))
          allocate (np_p   (ngp)   )
          allocate (ng_t (3,ngp,ke))
          allocate (ng_q (3,ngp,ke))
          allocate (ng_p (3,ngp)   )

          !---------------------------------
          ! allocate arrays passed to I/O PE
          !---------------------------------
          allocate (np    (nz,nocc))
          allocate (t     (nz,nocc))
          allocate (q     (nz,nocc))
          allocate (gp    (nz,nocc))
          allocate (pf    (nz,nocc))
          allocate (ps       (nocc))
          allocate (hs       (nocc))
          allocate (gu       (nocc))
          allocate (occs     (nocc))
          allocate (p    (n1b,nocc))
          allocate (lat  (n1b,nocc))
          allocate (lon  (n1b,nocc))
          allocate (em   (n1b,nocc))
          allocate (em_e (n1b,nocc))

          !-----------------------------------
          ! loop over occultations in this box
          !-----------------------------------
          iocc = 0
          do is = 1, o% n_spot
            s => o% spot(is)
            if (s% hd% obstype /= OT_GPSRO) cycle
            iocc = iocc + 1

            !------------------------------------------------
            ! copy atmospheric state to echam_fields module
            ! Initialization of echam_constituents_adj module
            !------------------------------------------------
            call echam_init_adj
            call load_occ  (o, s, occ, rs)
            occs (iocc) = occ

            !---------------------------------------------------
            ! gather bending angles and impact height (level 1b)
            !---------------------------------------------------
            p    (:,iocc) = pr
            em   (:,iocc) = invalid
            em_e (:,iocc) = invalid
            lat  (:,iocc) = invalid
            lon  (:,iocc) = invalid
            do i = 1, s% o% n
              j  = s% o% i + i
              iz = sum(minloc(abs(pr-o% olev(j))))
              if (bg_o% s(ib)% x(j) > 0._wp) then
                p    (iz,iocc) = o%             olev(j)
                em   (iz,iocc) = bg_o%   s(ib)% x   (j)
                em_e (iz,iocc) = bg_o_e% s(ib)% x   (j)
                lat  (iz,iocc) = rs(i)% geo% Phi
                lon  (iz,iocc) = rs(i)% geo% Lambda
              endif
            end do

            !------------------------------------------------------
            ! calculate refractivity and its derivatives (level 2a)
            !------------------------------------------------------
            do iz = 1, nz
               z (iz) = (iz-1) * dz_f
              Call Interpolate_Refractivity_adj  &
                (s% col% c% dLon,       & ! <- Longiude of point [deg]
                 s% col% c% dLat,       & ! <- Latitude of point [deg]
                            Z(iz),      & ! <- Altitude of point [km]
                            Zmin,       & ! -> Minimum model Z for this lon/lat
                            Zmax,       & ! -> Maximum model Z for this lon/lat
                            NP(iz,iocc),& ! -> Interpolated N
                            NG,         & ! -> Interpolated dN/d(alt,lat,lon)
                            NH,         & ! -> Interpolated hessian matrix of N
                            IDX,        & ! -> Subgrid indices [point, index]
                            NP_T,       & ! -> d(NP)/d(T(IGP,k))
                            NP_Q,       & ! -> d(NP)/d(Q(IGP,k))
                            NP_P,       & ! -> d(NP)/d(Psur(IGP))
                            NG_T,       & ! -> d(NG(i))/d(T(IGP,k))
                            NG_Q,       & ! -> d(NG(i))/d(Q(IGP,k))
                            NG_P,       & ! -> d(NG(i))/d(Psur(IGP))
                            T (iz,iocc),& ! ~> interpolated temperature
                            Q (iz,iocc),& ! ~> interpolated humidity
                            Gp(iz,iocc),& ! ~> interpolated geopotential
                            Pf(iz,iocc),& ! ~> interpolated Pressure
                            Ps(iocc),   & ! ~> surface pressure
                            Hs(iocc),   & ! ~> surface geopotential height
                            Gu(iocc)    ) ! ~> geoid undulation
              if (Z(iz) < Zmin .or. Z(iz) > Zmax) NP(iz,iocc) = invalid
            end do

            !--------
            ! cleanup
            !--------
            call echam_cleanup_adj
            deallocate (rs)

            !------------------------
            ! interpolate temperature
            !------------------------

            !------------------------------
            ! interpolate specific humidity
            !------------------------------

            !----------------------------------------------------------
            ! calculate bg errors for refractivity,temperature,humidity
            !----------------------------------------------------------

          end do

          !--------
          ! cleanup
          !--------
          call destruct_occ
!         call echam_cleanup
          deallocate (idx, np_t, np_q, np_p, ng_t, ng_q, ng_p)
        endif

        !------------------------------
        ! send results to I/O processor
        !------------------------------
        if (.not.dace% lpio) then
          call p_send (nocc, dace% pio, 1)
          if (nocc > 0) then
            call p_send (occs, dace% pio   )
            call p_send (np  , dace% pio, 1)
            call p_send (t   , dace% pio, 1)
            call p_send (q   , dace% pio, 1)
            call p_send (gp  , dace% pio, 1)
            call p_send (pf  , dace% pio, 1)
            call p_send (ps  , dace% pio, 1)
            call p_send (hs  , dace% pio, 1)
            call p_send (gu  , dace% pio, 1)
            call p_send (p   , dace% pio, 1)
            call p_send (em  , dace% pio, 1)
            call p_send (em_e, dace% pio, 1)
            call p_send (lat , dace% pio, 1)
            call p_send (lon , dace% pio, 1)
            deallocate  (np,t,q,gp,pf,ps,hs,gu,occs,p,em,em_e,lat,lon)
          endif
        endif
      else

        !---------------------------------
        ! receive results at I/O processor
        !---------------------------------
        if (dace% lpio) then
          call p_recv (nocc, o% pe, 1)
          if (nocc > 0) then
            allocate (np   (nz ,nocc))
            allocate (t    (nz ,nocc))
            allocate (q    (nz ,nocc))
            allocate (gp   (nz ,nocc))
            allocate (pf   (nz ,nocc))
            allocate (ps   (    nocc))
            allocate (hs   (    nocc))
            allocate (gu   (    nocc))
            allocate (occs (    nocc))
            allocate (p    (n1b,nocc))
            allocate (em   (n1b,nocc))
            allocate (em_e (n1b,nocc))
            allocate (lat  (n1b,nocc))
            allocate (lon  (n1b,nocc))
            call p_recv (occs, o% pe)
            call p_recv (np,   o% pe, 1)
            call p_recv (t,    o% pe, 1)
            call p_recv (q,    o% pe, 1)
            call p_recv (gp,   o% pe, 1)
            call p_recv (pf,   o% pe, 1)
            call p_recv (ps,   o% pe, 1)
            call p_recv (hs,   o% pe, 1)
            call p_recv (gu,   o% pe, 1)
            call p_recv (p,    o% pe, 1)
            call p_recv (em,   o% pe, 1)
            call p_recv (em_e, o% pe, 1)
            call p_recv (lat,  o% pe, 1)
            call p_recv (lon,  o% pe, 1)
          endif
        endif
      end if

      !--------------------------------
      ! write NetCDF file (ROPP format)
      !--------------------------------
      if (dace% lpio) then

        do iocc = 1, nocc
          !--------------------------------------
          ! set ROPP I/O data type (ROPIC levels)
          !--------------------------------------
          call ropp_io_init            &
            (rp,                       &
             n_lev_1a = 0,             &
             n_lev_1b = n1b,           &
             n_lev_2a = n_ropic_2,     &
             n_lev_2b = n_ropic_2,     &
             n_lev_2c = 1,             &
             n_lev_2d = 0)

          !--------------------------
          ! deallocate DWD extensions
          !--------------------------
          deallocate (rp% Lev1b% r_gns)
          deallocate (rp% Lev1b% v_gns)
          deallocate (rp% Lev1b% r_leo)
          deallocate (rp% Lev1b% v_leo)

          !------------------------------------------
          ! set header (occultation independent part)
          !------------------------------------------
          rp% PCD = 0                                  ! 1: valid
          rp% PCD = ibset (rp% PCD, PCD_NRT)           ! 1: off line
!         rp% PCD = ibset (rp% PCD, PCD_setting)       ! 0: Setting
          rp% PCD = ibset (rp% PCD, PCD_phase_nominal) ! 1: non-nominal
!         rp% PCD = ibset (rp% PCD, PCD_missing)       ! 0: valid
          rp% processing_centre = institution
          if (fc_hours == 0 .or. step==2) then
            rp% refrac_method   = 'interpolation of '//trim(model)//' analysis'
            rp% meteo_method    = 'interpolation of '//trim(model)//' analysis'
          else
            rp% refrac_method   = 'interpolation of '//trim(model)//' forecast'
            rp% meteo_method    = 'interpolation of '//trim(model)//' forecast'
          endif
          if (step==2) then
            filename            = 'monGPSRA.nc'
          else
            filename            = 'monGPSRP.nc'
            rp% PCD  = ibset(rp% PCD,PCD_occultation) ! background
          endif
          rp% software_version  = '3DVAR '//trim(var3d_version()) &
                                          //' '//var3d_date()
          rp% DTpro             =  run_time
          select case (operator)
          case (1)
            rp% bangle_method   = 'Abel-Transform'
          case (2)
            rp% bangle_method   = 'Abel-Transform (multiple tangent points)'
          case (3)
            rp% bangle_method   = 'Ray-Tracing'
          case default
            call finish ('write_gpsro_prof','invalid operator value')
          end select
          rp% bg% source   = institution
          rp% bg% year     = iyyyy (ana_time)
          rp% bg% month    = imm   (ana_time)
          rp% bg% day      = idd   (ana_time)
          rp% bg% hour     = ihh   (ana_time)
          rp% bg% minute   = imi   (ana_time)
          rp% bg% fcperiod = fc_hours

          !-------------------------------------------
          ! prepare header (occultation specific part)
          !-------------------------------------------
          rp% leo_id             = occs(iocc)% leoid
          rp% gns_id             = occs(iocc)% gnsid
!         rp% stn_id             =
          rp% DTocc              = occs(iocc)% time
!         rp% overall_qual       =
          call ropp_io_occid (rp) ! set rp% occ_id
          rp% georef% r_coc      = occs(iocc)% xlc%x * 1000._wp
          rp% georef% roc        = occs(iocc)% rlc   * 1000._wp
          rp% georef% lat        = occs(iocc)% gp1d% phi
          rp% georef% lon        = occs(iocc)% gp1d% lambda
          rp% georef% undulation = gu  (iocc)

          !--------------------------------------
          ! fill in level1b data (bending angles)
          !--------------------------------------
!         rp% azimuth_tp
!         rp% bangle_qual
          rp% Lev1b% impact         = p    (:, iocc) + rp% georef% roc
          where (em(:, iocc) /= invalid)
            rp% Lev1b% bangle       = em   (:, iocc)
            rp% Lev1b% bangle_sigma = em_e (:, iocc)
            rp% Lev1b% lat_tp       = lat  (:, iocc)
            rp% Lev1b% lon_tp       = lon  (:, iocc)
          endwhere
          if (operator == 1) then
            rp% Lev1b% lat_tp       = occs(iocc)% gp1d% phi
            rp% Lev1b% lon_tp       = occs(iocc)% gp1d% lambda
          endif

          !------------------------------------
          ! fill in level2a data (refractivity)
          !------------------------------------
!         rp% Lev2a% geop_refrac
!         rp% Lev2a% refrac_sigma
!         rp% Lev2a% refrac_qual
          do iz = 1, size (rp% Lev2a% alt_refrac)
            rp% Lev2a% alt_refrac (iz) = (iz-1) * dz_f   * 1.e3_wp
            if (iz <= nz) then
              if (np (iz, iocc) /= invalid) then
                rp% Lev2a% refrac (iz) = np (iz, iocc) * 1.e6_wp
                !------------------------------------
                ! level2b data (atmospheric profiles)
                !------------------------------------
                rp% Lev2b% geop       (iz) = Gp (iz, iocc)
                rp% Lev2b% press      (iz) = Pf (iz, iocc) / 100._wp   ! hPa
                rp% Lev2b% temp       (iz) = T  (iz, iocc)             ! K
                rp% Lev2b% shum       (iz) = Q  (iz, iocc) * 1000._wp  ! g/kg
!               rp% Lev2b% geop_sigma (iz) =
!               rp% Lev2b% press_sigma(iz) =
!               rp% Lev2b% temp_sigma (iz) =
!               rp% Lev2b% shum_sigma (iz) =
!               rp% Lev2b% meteo_qual (iz) =
              endif
            endif
          end do

          !------------------------------------------
          ! fill in level2c data (surface parameters)
          !------------------------------------------
          rp% Lev2c% geop_sfc        = hs (iocc)
          rp% Lev2c% press_sfc       = ps (iocc) / 100._wp ! hPa
!         rp% Lev2c% press_sfc_sigma =
!         rp% Lev2c% press_sfc_qual  =



          !-----------
          ! write file
          !-----------
          call ropp_io_write (rp                      ,&
                        file= path_file(aux,filename) ,&
                        type= 'netcdf'                ,&
                      append= append                   )
          append = .true.
          call ropp_io_free (rp)
        end do

        !-----------------------
        ! deallocate temporaries
        !-----------------------
        if (nocc>0) deallocate (np,t,q,gp,pf,ps,hs,gu,occs,p,em,em_e,lat,lon)

      endif
    end do
#endif
  end subroutine write_gpsro_prof
!==============================================================================
end module mo_occ_1d
