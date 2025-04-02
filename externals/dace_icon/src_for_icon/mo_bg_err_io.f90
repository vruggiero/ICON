!
!+ Derived type definitions for interpolation operator meta data
!
MODULE mo_bg_err_io
!
! Description:
!   Derived type definitions for interpolation operator meta data
!   used in the background error covariance model.
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
!  bugfix for zero number of observations in a box
! V1_5         2009/05/25 Harald Anlauf
!  Cleanup of interface of init_spot
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Andreas Rhodin
!  changed 3dvar revision numbers (CVS->SVN)
! V1_22        2013-02-13 Andreas Rhodin
!  extend error message
! V1_28        2014/02/26 Andreas Rhodin
!  preparations for VarEnKF
! V1_31        2014-08-21 Andreas Rhodin
!  revised verical differentiation in rttov specific interpolation
! V1_37        2014-12-23 Andreas Rhodin
!  implement RTTOV (Rochon) interpolation for model first guess (nwv_rad=4)
!  changes for Variational Ensemble Kalman Filter (VarEnKF)
! V1_42        2015-06-08 Andreas Rhodin
!  implement temporal interpolation for COSMO MEC
! V1_50        2017-01-09 Andreas Rhodin
!  option for higher order interpolation in apply_L_m
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Andreas Rhodin  DWD  2007-2008 original code
!==============================================================================
  !=============
  ! Modules used
  !=============
  use mo_kind,        only: wp,           &! working precision kind parameter
                            i8             ! 8-byte integer
  use mo_mpi_dace,    only: dace,         &! MPI group info
                            p_bcast        ! generic MPI bcast routine
  use mo_exception,   only: finish         ! abort in case of error
  use mo_t_obs,       only: t_obs,        &! observations derived type
                            t_spot,       &! report meta data
                            t_vic,        &! vertical   interpolation coef.type
                            t_hic,        &! horizontal interpolation coef.type
                            t_mcols,      &! grid column descriptor set
                            t_mcol,       &! component of t_mcols
                            col2lev,      &! set component 'levs' of 't_obs'
                            nwv_rad,      &! radiance vert.intp. flag
                            OBS_H,        &! observation id (surf.press.)
                            ITY_ICOL,     &! report type specifier
                            ITY_MCOLS      ! report type specifier
  use mo_fdbk_tables, only: OT_RAD         ! RADIANCE observation type flag
  use mo_atm_grid,    only: t_grid,       &! model grid  derived type
                            construct,    &! constructor for type t_grid
                            destruct       ! destructor  for type t_grid
  use mo_atm_state,   only: t_atm          ! atmospheric state derived type
  use mo_t_col,       only: get_cols,     &! redistribute model columns
                            set_cols_logp,&! store ln p instead of p
                            dealloc_cols, &!
                            t_cols,       &! data type to hold columns
                            t_col          ! component of t_cols
  use mo_t_bg_err_op, only: covm,         &! covariance matrix meta data
                            vars,         &! mapping: index in covm% hc2
                            set_vertinpc, &! set vertical interpolation coefs.
                            compress_cov, &! store covm only once
                            uncompress_cov ! store covm on each PE
  use mo_bg_err_op,   only: horz_inp,     &! test: 1=next-neighbour,2=bi-linear
                            vert_inp,     &! test: 1=next-neighbour,2=interpol.
                            horz_nmc_inp, &! test: 1:next neighbour,2=baryz.
                            horz_grid_inp  ! test: 1:next neighbour,2=baryz.old,
                                           !       3=baryz_3pt-new, 5=baryz_6pt-new
  use mo_grid_intpol, only: idx_init       ! determine grid indices
  use mo_fg_cov,      only: t_rowcol,     &! OI bg-error meta data type
                            construct,    &!
                            destruct,     &!
                            init_spot,    &! set correlation data
                            init_spot_obs  ! set corr. data type from obs.
  use mo_physics,     only: x2d            ! m -> degree (earth sphere)
  use mo_wmo_tables,  only: WMO6_GAUSSIAN,&!
                            WMO6_LATLON,  &!
                            WMO6_ROTLL,   &!
                            WMO8_J_POSITIVE! Flag: Points scan in +j direction
  use mo_t_tovs,      only: t_tovs,       &!
                            load,         &!
                            TTOVS_BASE

  implicit none

  !================
  ! Public entities
  !================
  private
  public :: t_intop     ! derived type to hold interpolation operators
  public :: t_pe        !   component  of t_intop
  public :: t_box       !   component  of t_intop
  public :: t_vic       !   component  of t_box
  public :: intop_obs   ! interpolation operator (to obs.locations) meta data
  public :: construct   ! constructor for t_intop
  public :: destruct    ! destructor  for t_intop
  public :: set_vic_ps  ! set vert.interp. coefficients for surface pressure
  public :: set_vic_atm ! set vert.interp. coefficients for the atmosphere
  public :: IL_H, IL_RH, IL_U,   IL_V   ! positions in array xl(,:,)
  public :: IV_H, IV_RH, IV_PSI, IV_CHI ! positions in array xv(,:,)
  public :: IV_TV,IL_TV, IV_U,   IV_V   !   for ensemble B
  !=============================
  ! derived data type definition
  !=============================
  integer, parameter :: max_id = 4 ! h, psi chi rh u v

  type t_box
    integer                :: ib      ! box index
    integer                :: nc_l    ! number of columns in this box
    integer                :: n_vic   ! number of 'levels' in this box
    type (t_vic)  ,pointer :: vic (:) ! vertical   interpolation coefficients
    type (t_hic)  ,pointer :: hic (:) ! horizontal interpolation coefficients
    integer       ,pointer :: ip  (:) ! parameter indices
  end type t_box

  type t_pe
    integer          :: lb (max_id) ! bounds for layers allocated on the PE, ..
    integer          :: ub (max_id) ! .. for IV_H, IV_RH, IV_PSI, IV_CHI
    integer          :: lbl(max_id) ! local bounds for layers
    integer          :: ubl(max_id) ! .. for IV_H, IV_RH, IV_PSI, IV_CHI
    integer          :: nl (max_id) ! no. layers
    integer          :: nc          ! number of columns required by the PE
    integer ,pointer :: ij (:,:)    ! lat lon indices   required by the PE
  end type t_pe

  type t_intop
    !--------------------
    ! general information
    !--------------------
    integer :: hc2size        ! check allocation status
    integer :: nlev           ! number of levels (64)
    integer :: np_v           ! number of fields  on grid
    integer :: np_l           ! number of fields  at locations of observations
    !-----------------------------
    ! information specific to a PE
    !-----------------------------
    integer                :: nc_v  ! number of columns on grid
    integer                :: n_box ! number of boxes
    type(t_mcols)          :: mc    ! grid column descriptor set
    type(t_box)   ,pointer :: bx(:) ! vertical interpolation coefficients
    integer          :: lbl(max_id) ! local bounds for slices allocated
    integer          :: ubl(max_id) ! .. for IV_H, IV_RH, IV_PSI, IV_CHI
    integer          :: nl (max_id) ! no. layers
    !----------------------------------
    ! information gathered from all PEs
    !----------------------------------
    type(t_pe)    ,pointer :: pe(:) ! indices required by the PEs
  end type t_intop

  !=================
  ! Module variables
  !=================
  type (t_intop) ,save :: intop_obs ! interpolation operator meta data
  !-------------------------------------------------
  ! index variables in 'model' space
  ! (before vertical interpolation/differentiation )
  !-------------------------------------------------
  integer, parameter :: IL_H   = 1  ! height
  integer, parameter :: IL_RH  = 2  ! relative humidity
  integer, parameter :: IL_U   = 3  ! u-wind component
  integer, parameter :: IL_V   = 4  ! v-wind component
  !-----------------------------------------------------
  ! index variables in B-matrix representation for NMC-B
  ! (before horizontal interpolation/differentiation )
  !-----------------------------------------------------
  integer, parameter :: IV_H   = 1  ! height
  integer, parameter :: IV_RH  = 2  ! relative humidity
  integer, parameter :: IV_PSI = 3  ! stream function
  integer, parameter :: IV_CHI = 4  ! velocity potential
  !-------------------------------------------
  ! for VarEnKF B u,v are represented directly
  ! tv may be represented directly
  !-------------------------------------------
  integer, parameter :: IV_U   = 3  ! u-wind component
  integer, parameter :: IV_V   = 4  ! v-wind component
  integer, parameter :: IV_TV  = 5  ! virtual temperature
  integer, parameter :: IL_TV  = 5  ! virtual temperature

  integer, parameter :: NTS = 1     ! number of time slots
  integer, parameter ::  TS = 1     ! time slot
  real(wp),parameter ::  TW = 0._wp ! time weight
  !===========
  ! Interfaces
  !===========
  interface construct
    module procedure construct_intop_obs   ! to observation space
    module procedure construct_intop_grid  ! to model grid
  end interface construct

  interface destruct
    module procedure destruct_intop
  end interface destruct

contains
!==============================================================================

  subroutine destruct_intop (intop)
  type (t_intop) ,intent(inout) :: intop
  !=================================================
  ! deallocate meta data for interpolation operators
  ! to locations of observation.
  !=============================================
    integer :: ibl, pe
    if (intop% hc2size == 0) return
    do ibl = 1, size (intop% bx)
      deallocate (intop% bx(ibl)% vic)
      deallocate (intop% bx(ibl)% hic)
      deallocate (intop% bx(ibl)% ip )
    end do
    deallocate (intop% bx)
    deallocate (intop% mc% c)
    if (associated (intop% pe)) then
      do pe = 0, dace% npe-1
        deallocate (intop% pe(pe)% ij)
      end do
      deallocate (intop% pe)
    endif
  end subroutine destruct_intop

!==============================================================================

  subroutine construct_intop_obs (intop, obs, atm)
  type (t_intop) ,intent(out)          :: intop   ! interpolation operator
  type (t_obs)   ,intent(in)           :: obs (:) ! observation meta data
  type (t_atm)   ,intent(in) ,optional :: atm     ! atmospheric state
  !===========================================================================
  ! derive meta data for interpolation operators to locations of observations.
  !
  ! if 'atm' is present, background error covariances are given by the
  ! background ensemble on the respective model grid.
  !
  ! if 'atm' is not present, background error covariances are given by
  ! the NMC method on a gaussian or lat-lon grid. In this case derived wind
  ! covariances are derived from streamfunction and vorticity covariances
  ! by horizontal differentiation. Variances are accounted for in observation
  ! space.
  !===========================================================================
    target                :: obs
    integer               :: ibl, ib  ! box index: local, global
    integer               :: is, i    ! spot index
    integer               :: ic       ! column index
    type(t_obs)  ,pointer :: o        ! pointer to obs currently processed
    type(t_box)  ,pointer :: b        ! pointer to box currently processed
    type(t_spot) ,pointer :: s        ! pointer to report
    type(t_grid) ,target  :: grid     ! definition of grid (h.cov.matrices)
    type(t_grid) ,pointer :: g        ! pointer to grid or g_ens
    type(t_rowcol)        :: col      ! OI bg model meta data
    integer               :: k
    integer               :: ncol     ! no.columns per spot
    integer               :: nlev     ! no.levels  per 'icol'
    integer               :: olev     ! level offset
    integer               :: ii,ix,j  ! derive latitudinal index
    integer               :: l1, ln   ! level index bounds
    real(wp)              :: lh       ! horz.length scale (temporary)
    integer               :: iw12     ! flag for horizontal differentiation
    logical               :: lvarenkf ! EnKF B, not NMC
    type(t_cols)          :: cols     ! pressure level columns
    real(wp) ,allocatable :: logp(:)  ! model levels
    integer               :: nwv      ! vertical interpolation flag
    type(t_mcol) ,pointer :: c (:)    ! temporary
    integer               :: lorder   ! local copy of 'order'
    type(t_tovs)          :: ttovs

    lvarenkf = present (atm)
    if (lvarenkf) then
    !----------------------------------------------------------
    ! set up coefficients for ensemble background on model grid
    !----------------------------------------------------------
      g => atm% grid
      intop% hc2size = g% size
      intop% np_v    = 5
      intop% np_l    = 5
      if (intop% hc2size == 0) return
      iw12 = 0
      allocate (logp (g% nz))
    else
    !-------------------------------------------------------
    ! set up coefficients for NMC background on lat-lon grid
    !-------------------------------------------------------
      intop% hc2size = covm% hc2_totalsize
      intop% np_v    = 4
      intop% np_l    = 4
      if (intop% hc2size == 0) return
      call construct (grid, gridtype= covm% gridtype, &
                                  ny= covm% ny      , &
                            scanmode= WMO8_J_POSITIVE )
      g    => grid
      iw12 = 1
    endif

    select case (g% gridtype)
    case (WMO6_GAUSSIAN, WMO6_LATLON, WMO6_ROTLL)
       lorder =      horz_inp
    case default
       lorder = min (horz_inp, 2)   ! no higher-order interpolation available
    end select

    call uncompress_cov
    intop% n_box = count (obs% pe == dace% pe)
    allocate (intop% bx(intop% n_box))

    !==========================================
    ! set horizontal interpolation coefficients
    !==========================================
    allocate (intop% mc% c   (100))
    allocate (intop% mc% idx (g% lbg(1):g% ubg(1), &
                              g% lbg(2):g% ubg(2), &
                              g% lbg(4):g% ubg(4), &
                              NTS                 ))
    intop% mc% idx = 0
    intop% mc% n   = 0
    intop% mc% pe  = dace% pe
    !------------------
    ! loop over 'boxes'
    !------------------
    ibl = 0
    do ib = 1, size(obs)
      o => obs(ib)
      if (o% pe /= dace% pe) cycle
      call col2lev (o)                 ! set column meta data in obs
      ibl = ibl + 1
#if defined (__SX__)
#define BX intop%bx(ibl)
#else
#define BX b
      b => intop% bx(ibl)
#endif
      allocate (BX% hic (o% n_spt))
      i = 0
      do is = 1, o% n_spot
        s => o% spot(is)
        select case (s% int_type)
        case (ITY_ICOL)
        !--------------------------------------------------
        ! horizontal interpolation in between model columns
        !--------------------------------------------------
          i = i + 1
          BX% hic(i)% iw12 = iw12
          call idx_init (s% col% c, BX% hic(i), intop% mc, &
                         0_i8, 0, g, TS, TW, order=lorder  )
          !---------------------------------------------------
          ! set horizontal length scale for gradient operators
          !---------------------------------------------------
          if (iw12 == 1) then
            lh = 0._wp
            do ii = 1, size(BX% hic(i)% imc,1)
              ix = BX% hic(i)% imc(ii,1)
              if (ix == 0) exit
              j = intop% mc% c(ix)% ijdtp(2)
              lh = lh + BX% hic(i)% w(ii) * covm% L_h(j)
            end do
            BX% hic(i)% w1   = BX% hic(i)% w1 * (lh * x2d)
            BX% hic(i)% w2   = BX% hic(i)% w2 * (lh * x2d)
          endif
          !---------------
          ! adjust indices
          !---------------
          BX% hic(i)% l    = s% l
          !------------------------------------------
          ! nearest neighbour interpolation (testing)
          !------------------------------------------
          if (horz_inp==1)then
            BX% hic(i)% w (1)  = 1
            BX% hic(i)% w (2:) = 0
            BX% hic(i)% w1(1)  = 1
            BX% hic(i)% w1(2:) = 0
            BX% hic(i)% w2(1)  = 1
            BX% hic(i)% w2(2:) = 0
          endif
        case (ITY_MCOLS)
        !-------------------------------------------------
        ! horizontal to model column positions (for GPSRO)
        !-------------------------------------------------
          ncol = size(s% imcol)
          if (mod(s% l% n, ncol) /= 0) call finish('construct_intop_obs', &
                                                   'mod(s%l%n, ncol) /= 0')
          nlev = s% l% n / ncol
          olev = 0
          do ic = 1, ncol
            i = i + 1
            BX% hic(i)% iw12 = iw12
            call idx_init (s% imcol(ic)% c, BX% hic(i), intop% mc, &
                           0_i8, 0, g, TS, TW, order=lorder        )
            !---------------------------------------------------
            ! set horizontal length scale for gradient operators
            !---------------------------------------------------
            if (iw12 == 1) then
              lh = 0._wp
              do ii = 1, size(BX% hic(i)% imc,1)
                ix = BX% hic(i)% imc(ii,1)
                if (ix == 0) exit
                j = intop% mc% c(ix)% ijdtp(2)
                lh = lh + BX% hic(i)% w(ii) * covm% L_h(j)
              end do
              BX% hic(i)% w1   = BX% hic(i)% w1 * (lh * x2d)
              BX% hic(i)% w2   = BX% hic(i)% w2 * (lh * x2d)
            endif
            !---------------
            ! adjust indices
            !---------------
            BX% hic(i)% l% i = olev + s% l% i
            BX% hic(i)% l% n = nlev
            olev             = olev + nlev
            !------------------------------------------
            ! nearest neighbour interpolation (testing)
            !------------------------------------------
            if (horz_inp==1)then
              BX% hic(i)% w (1)  = 1
              BX% hic(i)% w (2:) = 0
              BX% hic(i)% w1(1)  = 1
              BX% hic(i)% w1(2:) = 0
              BX% hic(i)% w2(1)  = 1
              BX% hic(i)% w2(2:) = 0
            endif
          end do
        end select
      end do
      !------------------
      ! consistency check
      !------------------
      if (i /= o% n_spt) then
        write(0,*)  'construct_intop_obs:  i /= n_spt :', i, o% n_spt
        call finish('construct_intop_obs','i /= n_spt')
      endif
    end do

    !============================================================
    ! gather model levels (logp) for vertical interpolation setup
    !============================================================
    if (lvarenkf) then
      call get_cols (intop% mc, atm, cols)
      call set_cols_logp (cols)
    endif

    !========================================
    ! set vertical interpolation coefficients
    !========================================
    !------------------
    ! loop over 'boxes'
    !------------------
    ibl = 0
    do ib = 1, size (obs)
      o => obs(ib)
      if (o% pe /= dace% pe) cycle
      ibl = ibl + 1
#if defined (__SX__)
#define BX intop%bx(ibl)
#else
#define BX b
      b => intop% bx(ibl)
#endif
      BX% ib    = ib                   ! box index
      BX% nc_l  = o% n_spt             ! no.columns in the box
      !========================================
      ! set vertical interpolation coefficients
      !========================================
      BX% n_vic = o% n_lev             ! number of vertical intpol. coefs.
      o% n_lev = 0
      allocate (BX% vic (BX% n_vic))   ! vertical interpolation coefficients
      if (lvarenkf) then
        do i = 1, o% n_spt
          l1 = BX% hic(i)% l% i + 1
          ln = BX% hic(i)% l% i + BX% hic(i)% l% n

          !-----------------------------------------
          ! derive mean logp of neighbour gridpoints
          !-----------------------------------------
          logp = 0._wp
          do k = 1, size(BX% hic(i)% imc,1)
            if (BX% hic(i)% imc(k,1) == 0) exit
            logp = logp + BX% hic(i)% w(k) * cols% col (BX% hic(i)% imc(k,1))% p
          end do
          call set_vertinpc (BX% vic(l1:ln),    &!  -> coefficients to set
                             o% levs(l1:ln)% z, &! <-  ln(p)
                             logp,              &! <-  covariance matrix meta data
                             covm% nwv )         ! <-  covariance matrix meta data
        end do
      else
        do is = 1, o% n_spot
          l1  = o% spot(is)% l% i + 1
          ln  = o% spot(is)% l% i + o% spot(is)% l% n
          nwv = -1
          if (o% spot(is)% hd% obstype == OT_RAD) nwv = min (nwv_rad, 3)
          if (nwv < 0)                            nwv = covm% nwv
          if (nwv == 3) then
            call load  (o, o%spot(is), tovs=ttovs, tovs_io=TTOVS_BASE)
            ln = l1 + ttovs%nlev - 1
          end if
          call set_vertinpc          &!
                 (BX% vic(l1:ln),    &!  -> coefficients to set
                  o% levs(l1:ln)% z, &! <-  ln(p)
                  covm% logp,        &! <-  covariance matrix meta data
            nwv = nwv                )! <-  covariance matrix meta data
        end do
      endif
      !------------------------------------------
      ! nearest neighbour interpolation (testing)
      !------------------------------------------
      if (vert_inp==1)then
        do k = 1, BX% n_vic
          BX% vic(k)% wh(1) =1
          BX% vic(k)% wh(2:)=0
          BX% vic(k)% wt(1) =1
          BX% vic(k)% wt(2:)=0
        end do
      endif
      !-----------------------------------
      ! error (stdev) in observation space
      !-----------------------------------
      if (iw12 == 1) then
        call construct (col, o% n_int, o% n_spt)
        do is = 1, o% n_spot
          s => o% spot(is)
          call init_spot_obs (col, s, o)
        end do
        do k = 1, BX% n_vic
          BX% vic(k)% i    = o% levs(k)% i
          BX% vic(k)% ezn  = col% cols(o% levs(k)%i%i+1)% ezn
          BX% vic(k)% exzn = col% cols(o% levs(k)%i%i+1)% exzn
        end do
        call destruct (col)
      else
        do k = 1, BX% n_vic
          BX% vic(k)% i    = o% levs(k)% i
        end do
      endif
      !----------------------
      ! set parameter indices
      !----------------------
      allocate (BX% ip  (o% n_int))    ! parameter indices
      if (o% n_int > 0) BX% ip = o% t_int
    end do

    !=========
    ! clean up
    !=========
    intop% nc_v = intop% mc% n  ! number of columns used
    if (size (intop% mc% c) /= intop% mc% n) then
      allocate (c (intop% mc% n))
      c = intop% mc% c (:intop% mc% n)
      deallocate (intop% mc% c)
      intop% mc% c => c
    endif
    deallocate (intop% mc% idx)
    if (lvarenkf) then
      call dealloc_cols (cols)
    else
      call destruct (grid)
    endif
    call compress_cov

    !=======================
    ! set transposition info
    !=======================
    if (lvarenkf) then
      nullify (intop% pe)
    else
      call set_vertinp_pe (intop)
    endif

  end subroutine construct_intop_obs
!------------------------------------------------------------------------------
  subroutine set_vertinp_pe (intop)
  type (t_intop) ,intent(inout) :: intop

    integer               :: pe
    integer               :: ic
    integer               :: k
    !--------------------------------
    ! gather information from all PEs
    !--------------------------------
    allocate (intop% pe(0:dace% npe-1))
    do pe = 0, dace% npe-1
      if (pe == dace% pe) intop% pe(pe)% nc = intop% nc_v
      call p_bcast (intop% pe(pe)% nc, pe)
      allocate (intop% pe(pe)% ij (2,intop% pe(pe)% nc))
      if (pe == dace% pe) then
        do ic=1,intop% nc_v
          intop% pe(pe)% ij (1:2,ic) = intop% mc% c(ic)% ijdtp(1:2)
        end do
      endif
      call p_bcast (intop% pe(pe)% ij, pe)
    end do

    !---------------------------------------------------------------------
    ! set local and global bounds for slices in covm% hc2_part for each PE
    !---------------------------------------------------------------------
    intop% lbl = 1
    intop% ubl = 0
    do pe = 0, dace% npe-1
      intop% pe(pe)% lb =  1
      intop% pe(pe)% ub =  0
      intop% pe(pe)% lbl=  1
      intop% pe(pe)% ubl=  0
      intop% pe(pe)% nl =  0
      do k = 1, size(vars)
        if (any(covm% hc2_part% pe == pe &
           .and.covm% hc2_part% i  == k )) then
          intop% pe(pe)% lb (k) = minval (covm% hc2_part% k,       &
                                   mask = covm% hc2_part% pe == pe &
                                    .and. covm% hc2_part% i  == k  )
          intop% pe(pe)% ub (k) = maxval (covm% hc2_part% k,       &
                                   mask = covm% hc2_part% pe == pe &
                                    .and. covm% hc2_part% i  == k  )
          intop% pe(pe)% lbl(k) = minval (covm% hc2_part% l,       &
                                   mask = covm% hc2_part% pe == pe &
                                    .and. covm% hc2_part% i  == k  )
          intop% pe(pe)% ubl(k) = maxval (covm% hc2_part% l,       &
                                   mask = covm% hc2_part% pe == pe &
                                    .and. covm% hc2_part% i  == k  )
          intop% pe(pe)% nl(k) = intop%pe(pe)% ub(k) - intop%pe(pe)% lb(k) + 1
        endif
      end do
    end do
    intop% lbl = intop% pe(dace% pe)% lbl
    intop% ubl = intop% pe(dace% pe)% ubl
    intop% nl  = intop% pe(dace% pe)% nl

  end subroutine set_vertinp_pe
!==============================================================================
  subroutine construct_intop_grid (intop, cbg, atm, order)
  type(t_intop) ,intent(out)          :: intop   ! interpolation coefficients
  type(t_cols)  ,intent(in)           :: cbg(:)  ! atmospheric reference state
  type(t_atm)   ,intent(in) ,optional :: atm     ! atmospheric state
  integer       ,intent(in) ,optional :: order   ! hor.intp: 2=lin, 4=higher

  !=============================================
  ! derive meta data for interpolation operators
  ! to locations of analysis grid.
  !=============================================
    type(t_grid) ,target  :: grid     ! definition of grid (h.cov.matrices)
    type(t_grid) ,pointer :: g        ! pointer to NMC or EnKF grid
    type(t_cols) ,pointer :: bg
    target                :: cbg
    integer               :: ib       ! box index
    integer               :: ic       ! column index
    integer               :: i,ix,j   ! derive latitudinal index
    real(wp)              :: lh       ! horz.length scale (temporary)
    logical               :: lvarenkf ! EnKF B, not NMC
    integer               :: lorder   ! local copy of 'order'
    integer               :: iw12     ! flag for horizontal differentiation
    type(t_mcol) ,pointer :: c (:)    ! temporary
!   type(t_box)  ,pointer :: b        ! pointer to box currently processed

    lvarenkf = present (atm)
    if (lvarenkf) then
    !----------------------------------------------------------
    ! set up coefficients for ensemble background on model grid
    !----------------------------------------------------------
      g => atm% grid
      intop% hc2size = g% size
      intop% np_v    = 5
      intop% np_l    = 5
      if (intop% hc2size == 0) then
        if(dace% lpio) write(6,*) 'construct_intop_grid: returning (intop% hc2size == 0)'
        return
      endif
      iw12 = 0
    else
    !-------------------------------------------------------
    ! set up coefficients for NMC background on lat-lon grid
    !-------------------------------------------------------
      intop% hc2size = covm% hc2_totalsize
      intop% np_v    = 4
      intop% np_l    = 4
      if (covm% hc2_totalsize == 0) then
        if(dace% lpio) write(6,*) 'construct_intop_grid: returning (covm% hc2_totalsize == 0)'
        return
      endif
      call construct (grid, gridtype= covm% gridtype, &
                                  ny= covm% ny,       &
                            scanmode= WMO8_J_POSITIVE )
      g    => grid
      iw12 = 1
    endif

    if (horz_grid_inp==-99) then
      lorder = horz_inp;      if (present (order)) lorder = order
    else
      lorder = horz_grid_inp; if (present (order)) lorder = order
    end if


    select case (g% gridtype)
    case (WMO6_GAUSSIAN, WMO6_LATLON, WMO6_ROTLL)
      if (lvarenkf) then
        lorder = horz_inp; if (present (order)) lorder = order
      else
        if (horz_nmc_inp == -99) then
          lorder = horz_inp; if (present (order)) lorder = order
        else
          lorder = horz_nmc_inp
        end if
      end if
       ! OK
    case default
       if (lorder == 4) then
         lorder = min(lorder, 2)     ! no higher-order interpolation available
       else
         lorder = lorder
       end if
    end select

    !------------------------
    ! allocate auxiliary data
    !------------------------
    allocate (intop% mc% c   (100))
    allocate (intop% mc% idx (g% lbg(1):g% ubg(1), &
                              g% lbg(2):g% ubg(2), &
                              g% lbg(4):g% ubg(4), &
                              NTS                 ))
    !------------------
    ! loop over 'boxes'
    !------------------
    intop% mc% idx = 0
    intop% mc% n   = 0
    intop% n_box   = size(cbg)
    allocate (intop% bx(intop% n_box))
    do ib = 1, size(cbg)
      bg => cbg(ib)
!     b => intop% bx(ib)
      intop% bx(ib)% ib    = ib                   ! box index
      intop% bx(ib)% nc_l  = bg% ncol             ! no.columns in the box
      intop% bx(ib)% n_vic = bg% ke + 1
      allocate (intop% bx(ib)% hic (bg% ncol))
      allocate (intop% bx(ib)% vic (2 *bg% ke + 1))
      allocate (intop% bx(ib)% ip  (0))
      !------------------
      ! loop over columns
      !------------------
      do ic = 1, bg% ncol
        !------------------------------------------
        ! set horizontal interpolation coefficients
        !------------------------------------------
        intop% bx(ib)% hic(ic)% iw12 = iw12
        call idx_init (bg% col(ic)% c, intop% bx(ib)% hic(ic), intop% mc,&
                       0_i8, 0, g, TS, TW, order=lorder                  )
        if (horz_inp==1)then
          intop% bx(ib)% hic(ic)% w(1) =1
          intop% bx(ib)% hic(ic)% w(2:)=0
        endif
        !------------------------------------
        ! interpolate horizontal length scale
        !------------------------------------
        if (iw12 == 1) then
          lh = 0._wp
          do i = 1, size(intop% bx(ib)% hic(ic)% imc,1)
            ix = intop% bx(ib)% hic(ic)% imc(i,1)
            if (ix == 0) exit
            j = intop% mc% c(ix)% ijdtp(2)
            lh = lh + intop% bx(ib)% hic(ic)% w(i) * covm% L_h(j)
          end do
          intop% bx(ib)% hic(ic)% w1   = intop% bx(ib)% hic(ic)% w1 * (lh * x2d)
          intop% bx(ib)% hic(ic)% w2   = intop% bx(ib)% hic(ic)% w2 * (lh * x2d)
        endif
        !-----------------------------------------------------------
        ! set up of vertical interpolation coefficients is postponed
        ! (calculated on the fly)
        !-----------------------------------------------------------
        intop% bx(ib)% hic(ic)% l% i = 0
        intop% bx(ib)% hic(ic)% l% n = 1                              ! PS only
      end do
    end do

    intop% nc_v = intop% mc% n  ! number of columns used
    if (size (intop% mc% c) /= intop% mc% n) then
      allocate (c (intop% mc% n))
      c = intop% mc% c (:intop% mc% n)
      deallocate (intop% mc% c)
      intop% mc% c => c
    endif
    if (lvarenkf) then
      nullify (intop% pe)
    else
      call set_vertinp_pe (intop)
    endif

    !--------------------------
    ! deallocate auxiliary data
    !--------------------------
    if (.not. lvarenkf) then
      call destruct (grid)
    endif
    deallocate    (intop% mc% idx)

  end subroutine construct_intop_grid
!------------------------------------------------------------------------------
  subroutine set_vic_ps (b, ehs, c)
  type(t_box)  ,intent(inout) :: b       ! interpolation coefficients
  real(wp)     ,intent(out)   :: ehs     ! surface geop.  observation error
  type(t_col)  ,intent(in)    :: c       ! atmospheric reference state
  !-------------------------------------------------------------
  ! set vertical interpolation coefficients for surface pressure
  !-------------------------------------------------------------
    type(t_rowcol)        :: col      ! OI bg model meta data
    real(wp)              :: z  (1)   ! vertical coordinates (ln p)
    real(wp)              :: err(1)   ! surface geop.  observation error

    z        = log (c% s% ps)
    b% n_vic = 1                      ! number of vertical intpol. coefs.
    call set_vertinpc (b% vic,       &!  -> coefficients to set
                       z,            &! <-  ln(p)
                       covm% logp,   &! <-  covariance matrix meta data
                       covm% nwv )    ! <-  covariance matrix meta data

    call construct (col, 1, 1)
    call init_spot (col,          &!
                    1,            &! spot index within the box
                    c% c,         &! coordinates
                    (/OBS_H/),    &! observation type identifiers
                    z,            &! vertical coordinates (ln p)
                    err=err)       ! observation error
    call destruct  (col)
    ehs = err(1)

  end subroutine set_vic_ps
!------------------------------------------------------------------------------
  subroutine set_vic_atm (b, e, z, t_int, c, logp)
  type(t_box) ,intent(inout)        :: b         ! interpolation coefficients
  real(wp)    ,intent(out)          :: e    (:)  ! observation error
  real(wp)    ,intent(in)           :: z    (:)  ! log(p) (destination)
  integer     ,intent(in)           :: t_int(:)  ! observation type
  type(t_col) ,intent(in)           :: c         ! atmospheric reference state
  real(wp)    ,intent(in) ,optional :: logp(:)   ! z (VarenKF B)
  !-----------------------------------------------------------
  ! set vertical interpolation coefficients for the atmosphere
  !-----------------------------------------------------------
    integer              :: nvc      ! number of variables / model column
    integer              :: ke       ! number of levels    / model column
    type(t_rowcol)       :: col      ! OI bg model meta data
    integer              :: k        ! index
    real(wp),allocatable :: zz (:)   ! levels

    nvc = size (t_int)
    ke  = (nvc-1)/5
    b% n_vic = 2 * ke + 1
    allocate (zz (b% n_vic))
    zz(1)       = z(1)
    zz(2:ke+1 ) = z(2:4*ke+1:4)
    zz(  ke+2:) = z(  4*ke+2: )

    if (present (logp)) then
      call set_vertinpc (b% vic,     &!  -> coefficients to set
                         zz    ,     &! <-  ln(p) (destination)
                         logp,       &! <-  ln(p) (source)
                         covm% nwv )  ! <-  covariance matrix meta data
      e = 0._wp                       ! not used
      b% vic(1) % i% i = 0
      b% vic(1) % i% n = 1
      do k = 1, ke
        b% vic(k+1)% i% i = (k-1)*4+1
        b% vic(k+1)% i% n = 4
      end do
      do k = 1, ke
        b% vic(ke + 1 + k)% i% i = 4*ke + k
        b% vic(ke + 1 + k)% i% n = 1
      end do
    else
      call set_vertinpc (b% vic,     &!  -> coefficients to set
                         zz    ,     &! <-  ln(p)
                         covm% logp, &! <-  covariance matrix meta data
                         covm% nwv )  ! <-  covariance matrix meta data
      call construct (col, nvc, 1)
      call init_spot (col,          &!
                      1,            &! spot index within the box
                      c% c,         &! coordinates
                      t_int,        &! observation type identifiers
                      z,            &! vertical coordinates (ln p)
                      err=e)         ! observation error
      b% vic(1) % i% i = 0
      b% vic(1) % i% n = 1
      b% vic(1)% ezn  = col% cols(1)% ezn
      b% vic(1)% exzn = col% cols(1)% exzn
      do k = 1, ke
        b% vic(k+1)% i% i = (k-1)*4+1
        b% vic(k+1)% i% n = 4
        b% vic(k+1)% ezn  = col% cols(b% vic(k+1)%i%i+1)% ezn
        b% vic(k+1)% exzn = col% cols(b% vic(k+1)%i%i+1)% exzn
      end do
      do k = 1, ke
        b% vic(ke + 1 + k)% i% i = 4*ke + k
        b% vic(ke + 1 + k)% i% n = 1
        b% vic(ke + 1 + k)% ezn  = col% cols(b% vic(ke + 1 + k)%i%i+1)% ezn
        b% vic(ke + 1 + k)% exzn = col% cols(b% vic(ke + 1 + k)%i%i+1)% exzn
      end do
    endif
    call destruct  (col)
    deallocate     (zz)
  end subroutine set_vic_atm
!==============================================================================
end module mo_bg_err_io
