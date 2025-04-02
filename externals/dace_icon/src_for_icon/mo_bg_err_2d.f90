!
!+ Covariance matrices, 2d(hor,wavelet)+1d(vert), high level entry routines
!
MODULE mo_bg_err_2d
!
! Description:
!   Operator representation of the wavelet based background error
!   covariance model based on the separability assumption regarding
!   vertical and horizontal covariance matrices.
!
!   This module holds high level entry routines:
!
!       - apply_B_oo:  multiply a vector in observation space by B
!       - apply_B_io:  multiply a vector by B, result is in interpolation space
!       - apply_B_oi:  multiply a vector by B, result is in observation space
!       - apply_B_mi:  multiply a vector by B, result is on model grid
!       - test_B_oo :  test the operator vs. the explicit representation
!
!   Low level routines and derived type definitions are placed
!   in module mo_t_bg_err_op.
!   Medium level routines are placed in module mo_bg_err_op.
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
!  changes for COSMO vertical coordinate, split off apply_lhlv
! V1_8         2009/12/09 Harald Anlauf
!  optimize for NEC SX-9
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Andreas Rhodin
!  new routine apply_B_oi
! V1_26        2013/06/27 Andreas Rhodin
!  minor changes
! V1_27        2013-11-08 Andreas Rhodin
!  implement Variational Bias Correction (VarBC)
! V1_28        2014/02/26 Andreas Rhodin
!  preparations for VarEnKF
! V1_31        2014-08-21 Andreas Rhodin
!  temporal correlations in random NMC B representation
!  preparations for rttov specific vertical interpolation
! V1_35        2014-11-07 Andreas Rhodin
!  for 1dvar-mode testing: explicitly multiply with B matrix (code disabled)
! V1_37        2014-12-23 Andreas Rhodin
!  changes for Variational Ensemble Kalman Filter (VarEnKF)
! V1_42        2015-06-08 Harald Anlauf
!  simple OpenMP parallelization of some routines; time interpolation for MEC
! V1_45        2015-12-15 Harald Anlauf
!  apply_L_m: option 'verbose' for controlling timing output
! V1_46        2016-02-05 Andreas Rhodin
!  base decisions on new flag 'vct', not 'ivctype'
! V1_47        2016-06-06 Andreas Rhodin
!  in 1dvar mode use explicit matrix if present (only works for radiances)
! V1_48        2016-10-06 Andreas Rhodin
!  cosmetic change
! V1_50        2017-01-09 Andreas Rhodin
!  option for higher order interpolation in apply_L_m
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Andreas Rhodin  DWD  2007-2008
! Harald Anlauf   DWD  2007-2008
!------------------------------------------------------------------------------
!
! Diagnostics for vectorization
!
#if defined (_FTRACE) && 0
#define FTRACE_BEGIN(text) CALL FTRACE_REGION_BEGIN (text)
#define FTRACE_END(text)   CALL FTRACE_REGION_END   (text)
#else
#define FTRACE_BEGIN(text)
#define FTRACE_END(text)
#endif
!==============================================================================
  !=============
  ! modules used
  !=============
  !-----------------
  ! general routines
  !-----------------
  use mo_kind,       only: wp               ! precision kind parameter
  use mo_exception,  only: finish           ! abort on error condition
  use mo_mpi_dace,   only: dace,           &! MPI group info
                           p_barrier,      &! MPI barrier routine
                           p_max,          &! generic MPI max over PEs
                           p_min,          &! generic MPI min over PEs
                           p_sum,          &! generic MPI sum
                           p_alltoall       ! generic MPI alltoall
  use mo_cpu_time,   only: stop_time        ! determine cpu and wall time
  !-------------
  ! observations
  !-------------
  use mo_t_obs,      only: t_spot,         &! report     meta data type
                           obsq,           &! table of observed quantities
                           OBS_TV,         &! virtual temperature  id
                           OBS_RH,         &! relative humidity    id
                           OBS_H,          &! geopotential height  id
                           OBS_HS,         &! surface geopotential id
                           OBS_U,          &! u-wind               id
                           OBS_V,          &! v-wind               id
                           OBS_DUM          ! dummy type
  use mo_obs_set,    only: t_obs_set        ! observation data type
  use mo_fdbk_tables,only: OT_RAD           ! Radiances report type ID
  !------------------
  ! atmospheric state
  !------------------
  use mo_t_col,      only: t_cols,         &! data type to hold columns
!                          get_cols,       &! redistribute model columns
!                          get_cols_adj,   &! adjoint routine
                           alloc_cols,     &! allocate components of t_cols
                           dealloc_cols,   &! deallocate t_cols components
                           assignment(=),  &! assign t_cols
                           COL_T, COL_UV,  &! identifier for temp., hor.wind,
                           COL_RH,         &!   rel.humidity
                           COL_GEOH         !   geopot.height (half levels)
  use mo_physics,    only: gacc,           &! gravity acceleration
                           R                ! gas constant of dry air
  use mo_wmo_tables, only: WMO3_ISOBARIC    ! level type
  use mo_atm_state,  only: t_atm,          &! atmospheric state derived type
                           assignment(=),  &! t_atm = real
                           construct,      &! set up t_atm
                           destruct         ! deallocate t_atm components
  use mo_atm_grid,   only: VCT_P_HYB        ! GME/HRM Vertical coordinate
  !-----------------------------------
  ! background error covariance module
  !-----------------------------------
  use mo_tovs,       only: lBii             ! use explicit matrix in intp.space
  use mo_t_bg_err_op,only: covm,           &! NMC fitted covariances
                           compress_cov,   &! store covm only once
                           uncompress_cov, &! store covm on each PE
                           random_gauss     ! Gaussian random numbers
  use mo_bg_err_op,  only: apply_B_ii_1d,  &! B * vector in intp.space
                           apply_B_mi_1d,  &! B * vector -> model grid
                           repr_2dh,       &! horizontal operator repres. flag
                           test_bg,        &! compare with explicit represent.
                           test_bg_b,      &! bound for failure of test_bg
                           test_adj,       &! compare linear vs adjoint
                           l_L_h,          &! apply L matrix in the horizontal
                           l_W_h,          &! apply W matrix in the horizontal
                           ltrace,         &! trace flag
                           l1dvar,         &! in 1dvar mode
                           n_zonfilt,      &! number of latitudes to filter
                           zonfilt          ! number of coefficients of filter
  use mo_bg_err_io,  only: t_box,          &! derived type definition
                           t_intop,        &! interpolation operator data type
                           construct,      &! set up type t_intop
                           destruct,       &! deallocate  t_intop
                           set_vic_ps,     &! set vert.interp.coeffs. for ps
                           set_vic_atm,    &! set vert.interp.coeffs. for atm.
                     io => intop_obs,      &! interpolation operator meta data
                           IL_H, IL_RH, IL_U,   IL_V, &! positions in xl(,:,)
                           IV_TV,IL_TV, IV_U,   IV_V, &!   in xv(,:,) (Ens-B)
                           IV_H, IV_RH, IV_PSI, IV_CHI !   in xv(,:,) (NMC_B)
  use mo_bg_err_ens, only: apply_B_mi_ensb,&! apply Im B Io^t
                           apply_B_ii_ensb  ! apply Io B Io^t

  use mo_varenkf,    only: w_nmc_b          ! weight of NMC      B
  !----------------------------
  ! variational bias correction
  !----------------------------
  use mo_radbias_3dv,only: biascor_mode,  &! bias correction mode
                           BC_VARBC        ! VBC flag value
  use mo_vbc,        only: Hvbc_x,        &! bias correction forward model
                           z_Hvbc,        &! bias correction adjoint model
                           Bvbc_z          ! bias correction B application
  !--------------------------------------
  ! wavelet representation and operations
  !--------------------------------------
  use mo_cov_wavelet,only: t_cov_wavelet   ! Hor. wavelet cov. matrix
  use mo_1dmra,      only: wave_1d,       &! 1d wavelet transform
                           wave_1dt,      &! 1d wavelet transform (transp.)
                           TR_SYN,        &! flag: wavelet synthesis
                           TR_ADJ          !       adjoint synthesis
  use mo_transform,  only: fft             ! fast Fourier transform
  !--------------------------------------
  ! parallel matrix and vector operations
  !--------------------------------------
  use mo_allocate,   only: enter_function,&! to be called at start of routines
                           leave_function  ! to be called at end of routines
  use mo_dec_matrix, only: t_vector,      &! partitioned vector data type
                           t_matrix,      &! partitioned matrix data type
                           construct,     &! construct a vector
                           destruct,      &! destruct  a vector
                           delete_storage,&! deallocate pointer components
                           operator(+),   &! vector sum
                           operator(*),   &! matrix vector product
                           operator(/),   &! vector / vector
                           operator(-),   &! vector difference
                           operator(/=),  &! vector element compare
                           assignment(=), &! vector assignment
                           sub_block,     &! get sub-matrix
!                          FULL,          &! full representation flag value
                           size,          &! size of a vector
                           minval,maxval, &! minval, maxval of a vector
                           count,         &! count true elements of a vector
                           random_gauss,  &! set vector to random field
                           norm2,         &! 2-norm of a vector
                           sum             ! sum of vector components
  use mo_matrix,     only: operator(.x.)   ! matrix-vector multiplication
  use mo_random,     only: random_gauss    ! normal distribution -> test_B_adj
  implicit none

  !================
  ! public entities
  !================
  private
  public :: apply_B_oo         ! multiply a vector in obs.space by B
  public :: apply_B_io         ! multiply a vector by B, result in intp.space
  public :: apply_B_oi         ! multiply a vector by B, result in obs. space
  public :: apply_B_mi         ! multiply a vector by B, result on model grid
  public :: apply_B_ii         ! multiply a vector in intp.space by B
  public :: apply_B_ii_2d      ! multiply a vector in intp.space by B
  public :: apply_L_o          ! multiply a vector by L, result in obs.space
  public :: apply_L_o_t        ! multiply a vector in obs.space by L_t
  public :: apply_L_i          ! multiply a vector by L, result in interp.space
  public :: apply_L_i_t        ! multiply a vector in interp.space by L_t
  public :: apply_L_m          ! multiply a vector by L, result on model grid
  public :: apply_WhLh         ! multiply a vector in interp.space by Wh*Lh
  public :: test_B_oo          ! test the operator implementation of B
  public :: test_B_adj         ! adjoint tests

  !===========
  ! Interfaces
  !===========
  interface trace
    module procedure trace_r   ! print min,max,variance of real
    module procedure trace_r1  ! print min,max,variance of real (:)
    module procedure trace_r2  ! print min,max,variance of real (:,:)
    module procedure trace_r3  ! print min,max,variance of real (:,:,:)
    module procedure trace_v   ! print min,max,variance of type t_vector
  end interface trace

contains
!==============================================================================
! top level entries
!==============================================================================
  subroutine apply_B_oo (obs, x, y, e_fi, yi)
  type (t_obs_set) ,intent(inout) :: obs  ! observation meta data
  type (t_vector)  ,intent(in)    :: x    ! input
  type (t_vector)  ,intent(inout) :: y    ! output
  type (t_vector)  ,intent(in)    :: e_fi ! background error
  type (t_vector)  ,intent(inout)        &!
                ,optional, target :: yi   ! B*x in interpolation space
  !--------------------------------------------------------
  ! multiply a vector in observation space by the B matrix.
  ! result is a vector in observation space.
  ! this routine also accounts for the VarBC B term.
  !--------------------------------------------------------
    type (t_vector)          :: xi   ! intermediate storage
    type (t_vector) ,target  :: yiii
    type (t_vector) ,pointer :: yii

    call enter_function
    !------------------------------
    ! allocate intermediate storage
    !------------------------------
    call construct   (xi, obs% ii, 'x_intp')
    if (present(yi)) then
      yii => yi
    else
      yii => yiii
      call construct (yii, obs% ii, 'y_intp')
    endif
    !--------------------------------------------------
    ! xi = H_t * x : observation -> interpolation space
    !--------------------------------------------------
    xi = x * obs% l% H
    call z_Hvbc (x, obs% vbc% z, obs)
    !---------------------------------------------------
    ! yi = B * xi : interpolation -> interpolation space
    !---------------------------------------------------
    call apply_B_ii (obs, xi, yii, e_fi)
    call Bvbc_z (obs% vbc% cf, obs% vbc% z, obs)
    !------------------------------------------------
    ! y = H * yi : interpolation -> observation space
    !------------------------------------------------
    y = obs% l% H * yii
    if (iand (BC_VARBC, biascor_mode) /= 0) then
      call Hvbc_x (obs% vbc% bc, obs% vbc% cf, obs)
      y = y + obs% vbc% bc
    endif
    !---------
    ! clean up
    !---------
    if (.not. present(yi)) call destruct (yii)
    call destruct (xi)
    call delete_storage (x)
    call delete_storage (y)
    call delete_storage (e_fi)
    call leave_function
  end subroutine apply_B_oo
!-----------------------------------------------------------------------------
  subroutine apply_B_io (obs, x, y, e_fi)
  type (t_obs_set) ,intent(inout) :: obs  ! observation meta data
  type (t_vector)  ,intent(in)    :: x    ! input
  type (t_vector)  ,intent(inout) :: y    ! output
  type (t_vector)  ,intent(in)    :: e_fi ! background error
  !-------------------------------------------------------
  ! multiply a vector in observation space by the B matrix
  ! result is a vector in interpolation space:
  ! y = B * H^t x
  !-------------------------------------------------------

    type (t_vector) :: xi     ! intermediate storage in interpolation space
    call enter_function
    call construct (xi, obs% ii, 'x_intp')
    !--------------------------------------------------
    ! xi = H_t * x : observation -> interpolation space
    !--------------------------------------------------
    xi = x * obs% l% H
    !--------------------------------------------------
    ! y = B * xi : interpolation -> interpolation space
    !--------------------------------------------------
    call apply_B_ii (obs, xi, y, e_fi)
    !---------
    ! clean up
    !---------
    call destruct (xi)
    call delete_storage (x)
    call delete_storage (y)
    call delete_storage (e_fi)
    call leave_function

  end subroutine apply_B_io
!-----------------------------------------------------------------------------
  subroutine apply_B_oi (obs, x, y, e_fi)
  type (t_obs_set) ,intent(inout) :: obs  ! observation meta data
  type (t_vector)  ,intent(in)    :: x    ! input
  type (t_vector)  ,intent(inout) :: y    ! output
  type (t_vector)  ,intent(in)    :: e_fi ! background error
  !-------------------------------------------------------
  ! multiply a vector in interpolation space by the B matrix
  ! result is a vector in observation space
  !-------------------------------------------------------

    type (t_vector) :: xi     ! intermediate storage in interpolation space
    call enter_function
    call construct (xi, obs% ii, 'x_intp')
    !---------------------------------------------------
    ! xi = B * xi : interpolation -> interpolation space
    !---------------------------------------------------
    call apply_B_ii (obs, x, xi, e_fi)
    !------------------------------------------------
    ! y = H * xi : interpolation -> observation space
    !------------------------------------------------
    y = obs% l% H * xi
    !---------
    ! clean up
    !---------
    call destruct (xi)
    call delete_storage (x)
    call delete_storage (y)
    call delete_storage (e_fi)
    call leave_function

  end subroutine apply_B_oi
!-----------------------------------------------------------------------------
  subroutine apply_L_o (x, y, obs, e_fi, dummy)
  real(wp)         ,intent(in)    :: x(:,:,:) ! input
  type (t_vector)  ,intent(inout) :: y        ! output
  type (t_obs_set) ,intent(inout) :: obs      ! observation meta data
  type (t_vector)  ,intent(in)    :: e_fi     ! background error
  type (t_vector)  ,intent(in)    :: dummy    ! dummy observations
  optional                        :: dummy
  !----------------------------------------
  ! multiply a vector by the L matrix
  ! result is a vector in observation space
  !----------------------------------------
    type (t_vector) :: xi     ! intermediate storage in interpolation space
    call enter_function
    call construct (xi, obs% ii, 'x_intp')
    call apply_L_i (x, xi, e_fi, dummy)
    y = obs% l% H * xi
    call destruct (xi)
    call delete_storage (y)
    call delete_storage (e_fi)
    call leave_function
  end subroutine apply_L_o
!-----------------------------------------------------------------------------
  subroutine apply_L_o_t (x, y, obs, e_fi, dummy)
  type (t_vector)  ,intent(in)    :: x        ! input
  real(wp)         ,intent(out)   :: y(:,:,:) ! output
  type (t_obs_set) ,intent(inout) :: obs      ! observation meta data
  type (t_vector)  ,intent(in)    :: e_fi     ! background error
  type (t_vector)  ,intent(inout) :: dummy    ! dummy observations
  optional                        :: dummy
  !----------------------------------------------
  ! multiply a vector in observation space by L_t
  !----------------------------------------------
    type (t_vector) :: xi     ! intermediate storage in interpolation space
    call enter_function
    call construct (xi, obs% ii, 'x_intp')
    xi = x * obs% l% H
    call apply_L_i_t (xi, y, e_fi, dummy)
    call destruct (xi)
    call delete_storage (x)
    call delete_storage (e_fi)
    call leave_function
  end subroutine apply_L_o_t
!-----------------------------------------------------------------------------
  subroutine apply_B_ii_ex (obs, x, y, e_fi)
  type (t_obs_set) ,intent(inout) :: obs  ! observation meta data
  type (t_vector)  ,intent(in)    :: x    ! input
  type (t_vector)  ,intent(inout) :: y    ! output
  type (t_vector)  ,intent(in)    :: e_fi ! background error
  !++++++++++++++++++++++++++++++++++++++++++++++++
  ! use explicit matrix for radiances in 1dvar mode
  !++++++++++++++++++++++++++++++++++++++++++++++++
    integer                :: ib, is
    type (t_spot), pointer :: si
    do ib = 1,size(obs% o)
      if (obs% o(ib)% pe /= dace% pe) cycle
        do is = 1,obs% o(ib)% n_spot
        if (obs% o(ib)% spot(is)% hd% obstype /= OT_RAD) cycle
        si => obs% o(ib)% spot(is)
        y% s(ib)% x(si%i%i+1:si%i%i+si%i%n)                                   &
          = sub_block (obs% b% Bii% b(ib,ib), si%i%i, si%i%n, si%i%i, si%i%n) &
          * x% s(ib)% x(si%i%i+1:si%i%i+si%i%n)
      end do
    end do

  end  subroutine apply_B_ii_ex
!-----------------------------------------------------------------------------
  subroutine apply_B_ii (obs, x, y, e_fi)
  type (t_obs_set) ,intent(inout) :: obs  ! observation meta data
  type (t_vector)  ,intent(in)    :: x    ! input
  type (t_vector)  ,intent(inout) :: y    ! output
  type (t_vector)  ,intent(in)    :: e_fi ! background error
  !---------------------------------------------------------
  ! multiply  a vector in interpolation space by the B matrix
  ! result is a vector in interpolation space
  !---------------------------------------------------------
    type (t_vector) :: sink     ! intermediate storage (sink variables)
    !------------------------------------------------
    ! use explicit matrix for radiances in 1dvar mode
    !------------------------------------------------
    if (l1dvar .and. lBii) then
      call apply_B_ii_ex (obs, x, y, e_fi)
    !----------------------------
    ! use operator representation
    !----------------------------
    else
      select case (repr_2dh)
      case (0)
        call apply_B_ii_1d (obs, x, y, e_fi)
      case (1)
        call trace ('apply_B_ii_1d: x =',x, io% bx)
        call apply_B_ii_1d (obs, x, y, e_fi)
        call trace ('apply_B_ii_1d: y =',y, io% bx)
      case (2)
        call construct (sink, obs% di, 'sink')
        call apply_B_ii_2d   (x, y, e_fi, sink)
        call apply_B_ii_ensb (x, y)
        call destruct (sink)
      case default
        write(0,*)  'apply_B_ii:  repr_2dh =',repr_2dh
        call finish('apply_B_ii','invalid value for repr_2dh')
      end select
    endif
  end subroutine apply_B_ii
!-----------------------------------------------------------------------------
  subroutine apply_B_mi (a_m, cbg, obs, x, lnewpl, e_fi, pert_B)
  type(t_cols)   ,intent(out)         :: a_m(:) ! ana. increment (model space)
  type(t_cols)   ,intent(in)          :: cbg(:) ! reference state
  type(t_obs_set),intent(inout)       :: obs    ! observation meta data
  type(t_vector) ,intent(in)          :: x      ! input
  logical        ,intent(in)          :: lnewpl ! ana. increment on new p-levs
  type(t_vector) ,intent(in)          :: e_fi   ! background error
  real(wp)       ,intent(in),optional :: pert_B ! model perturbation factor
  !-------------------------------------------------------------
  ! multiply a vector in interpolation space by the NMC B matrix
  ! result is on model grid
  !-------------------------------------------------------------
    select case (repr_2dh)
    case (0)
      call apply_B_mi_1d (a_m, cbg, obs, x, lnewpl, e_fi)
    case (1)
      call trace ('apply_B_mi_1d: x =',x, io% bx)
      call apply_B_mi_1d (a_m, cbg, obs, x, lnewpl, e_fi)
    case (2)
      call apply_B_mi_2d   (a_m, cbg, x, lnewpl, e_fi, pert_B)
      call apply_B_mi_ensb (a_m, cbg, x, lnewpl, e_fi)
    case default
      write(0,*)  'apply_B_mi:  repr_2dh =',repr_2dh
      call finish('apply_B_mi','invalid value for repr_2dh')
    end select
  end subroutine apply_B_mi
!-----------------------------------------------------------------------------
  subroutine apply_B_mi_2d (a_m, cbg, x, lnewpl, e_fi, pert_B)
  type(t_cols)  ,intent(out)         :: a_m(:) ! ana. increment (model space)
  type(t_cols)  ,intent(in)          :: cbg(:) ! reference state
  type(t_vector),intent(in)          :: x      ! input
  logical       ,intent(in)          :: lnewpl ! ana. increment on new p-levs
  type(t_vector),intent(in)          :: e_fi   ! background error
  real(wp)      ,intent(in),optional :: pert_B ! model perturbation factor
  !---------------------------------------------------------
  ! multiply a vector in interpolation space by the B matrix
  ! result is on model grid
  !---------------------------------------------------------
    real(wp)              :: xd (covm% nx, covm% ny, covm% hc2_part_size)
    real(wp) ,allocatable :: xp (:,:,:)
    !----------------
    ! multiply by L^T
    !----------------
    call stop_time ('L_i_t')
    call apply_L_i_t (x, xd, e_fi)
    !------------------------
    ! scale with B_nmc weight
    !------------------------
    xd = xd * w_nmc_b
    !-------------------------------------
    ! optionally: add random perturbations
    !-------------------------------------
    if (present (pert_B)) then
      if (pert_B > 0._wp) then
        allocate (xp (covm% nx, covm% ny, covm% hc2_part_size))
        call random_gauss (xp, covm)
        xd = xd + pert_B * xp
        deallocate (xp)
      endif
    endif
    !----------------
    ! multiply by L^T
    !----------------
    call stop_time ('L_m')
    call apply_L_m   (xd, a_m, cbg, lnewpl)

  end subroutine apply_B_mi_2d
!-----------------------------------------------------------------------------
! Test the Covariance matrix implementation
!-----------------------------------------------------------------------------
  subroutine test_B_oo (obs, e_f, e_fi, B)
  type (t_obs_set) ,intent(inout) :: obs  ! observation meta data
  type (t_matrix)  ,intent(in)    :: B    ! explicit implementation of B
  type (t_vector)  ,intent(in)    :: e_f  ! background error
  type (t_vector)  ,intent(in)    :: e_fi ! background error
  !-----------------------------------------
  ! compare the operator implementation of B
  !    with the explicit implementation
  !-----------------------------------------
    type (t_vector) :: x                ! random vector
    type (t_vector) :: s                ! random vector, scaled, masked
    type (t_vector) :: yo               ! B * x  (operator representation)
    type (t_vector) :: ye               ! B * x  (explicit representation)
    real(wp)        :: diff, no, ne, nx ! vector norm, difference
    integer         :: n1               ! max (vector size, 1)

    integer                :: ib      ! block index
    integer                :: is      ! spot index
    integer                :: io      ! observation index
    integer                :: it, it2 ! test index
    logical                :: failed  ! test failed
    type (t_spot) ,pointer :: si      ! spot index

    real(wp) :: diffs (0:size(obsq))
    real(wp) :: mval  (0:size(obsq),0:size(obsq))

    if (.not. test_bg  ) return
    if (covm% valid < 2) return
    !----------------------------
    ! allocate, set random vector
    !----------------------------
    call construct (x,  obs% oi, 'x' )
    call construct (s,  obs% oi, 's' )
    call construct (yo, obs% oi, 'yo')
    call construct (ye, obs% oi, 'ye')
    call random_gauss (x)
    mval = 0._wp
    !----------------------------
    ! loop over observation types
    !----------------------------
    do it = size (obsq), 0, -1
      !--------------------------
      ! scale, mask random vector
      !--------------------------
      s = x / e_f
      if (it>0) then
        do ib = 1,size(obs% o)
          if (obs% o(ib)% pe /= dace% pe) cycle
          do is = 1,obs% o(ib)% n_spot
            si => obs% o(ib)% spot(is)
            do io = si% o% i + 1, si% o% i + si% o% n
              if (obs% o(ib)% varno(io) /= obsq(it)% key) &
                s% s(ib)% x(io) = 0._wp
            end do
          end do
        end do
      endif
      !------------------------------
      ! matrix multiplications, scale
      !------------------------------
      select case (count (s /= 0._wp))
      case (0)
        yo = 0._wp
      case default
        call apply_B_oo (obs, s, yo, e_fi)
      end select
      yo     = yo      / e_f
      ye     = (B * s) / e_f
      !------------
      ! differences
      !------------
      n1 = max (1,ye% n)
      if (it==0) then
        nx     = norm2 (x)  / n1
        no     = norm2 (yo) / n1
        ne     = norm2 (ye) / n1
      endif
      s      = ye - yo
      diff   = norm2 (s)  / n1
      diffs(it) = diff
      !------------------------------------------------
      ! matrix of differences according to lhs,rhs type
      !------------------------------------------------
      do it2 = 0, size (obsq)
        do ib = 1,size(obs% o)
          if (obs% o(ib)% pe /= dace% pe) cycle
          do is = 1,obs% o(ib)% n_spot
            si => obs% o(ib)% spot(is)
            do io = si% o% i + 1, si% o% i + si% o% n
              if (it2==0) then
                mval(it2,it) = max(mval(it2,it), abs(s% s(ib)% x(io)))
              elseif (obs% o(ib)% varno(io) == obsq(it2)% key) then
                mval(it2,it) = max(mval(it2,it), abs(s% s(ib)% x(io)))
              endif
            end do
          end do
        end do
      end do
    end do
    mval   = p_max (mval)
    failed = diffs(0) > test_bg_b
    !--------------
    ! print results
    !--------------
    if (dace% lpio) then
      write(6,'(a)') repeat ('-',79)
      write(6,'(a)')
      write(6,'(a)') 'test_B_oo: test the operator implementation of B'
      write(6,'(a)')
      write(6,'(a)') &
      '   ib   is    i  typ       x      y_explicit  y_operator  difference'
      !--------------
      ! example block
      !--------------
      do ib = 1,size(obs% o)
        if (obs% o(ib)% pe /= dace% pe) cycle
        do is = 1,obs% o(ib)% n_spot
          do io = obs% o(ib)% spot(is)% o% i + 1, &
                  obs% o(ib)% spot(is)% o% i + obs% o(ib)% spot(is)% o% n
            write(6,'(4i5,4f12.7)') ib,is,io,obs% o(ib)% varno(io),     &
                                    x% s(ib)% x(io),                    &
                                   ye% s(ib)% x(io),  yo% s(ib)% x(io), &
                                   ye% s(ib)% x(io) - yo% s(ib)% x(io)
          end do
          write(6,'()')
        end do
        if (.not.failed) exit
        write(6,'()')
      end do
      !------------------------------------------------
      ! matrix of differences according to lhs,rhs type
      !------------------------------------------------
      write (6,'(5x,a8,a40,20f12.7)') 'all', 'all', mval(0,:)
      do it = 1, size (obsq)
        write (6,'(i4,1x,a8,a40,20f12.7)') &
          obsq(it)% key, obsq(it)% name, obsq(it)% desc, mval(it,:)
      end do
      write (6,'(5x,8x,a40,20f12.7)') 'norm', diffs(:)
      write (6,'(5x,8x,40x,20a12)') 'all',obsq(:)% name
      !-------------
      ! final result
      !-------------
      write(6,'(a)')
      write(6,'(a,4e15.3)') '  norm, diff =',nx,ne,no,diff
      write(6,'(a)')
    endif
    !-------------------------
    ! abort if the test failed
    !-------------------------
    if (failed) then
      call p_barrier
      call finish ('test_B_oo','test failed')
    endif
    !---------
    ! clean up
    !---------
    call destruct (x)
    call destruct (s)
    call destruct (yo)
    call destruct (ye)
  end subroutine test_B_oo
!==============================================================================
! 2d horizontal + 1d vertical (separable) operator approach
!==========================================================
  subroutine apply_B_ii_2d (x, y, e_fi, sink)
  type (t_vector)  ,intent(in)    :: x     ! input
  type (t_vector)  ,intent(inout) :: y     ! output
  type (t_vector)  ,intent(in)    :: e_fi  ! background error
  type (t_vector)  ,intent(inout) :: sink  ! sink observations
  optional                        :: sink
  !---------------------------------------------------------
  ! multiply  a vector in interpolation space by the B matrix
  ! result is a vector in interpolation space
  !---------------------------------------------------------
    real(wp), allocatable :: xd (:,:,:) ! vector in wavelet representation
    real(wp), allocatable :: xh (:,:,:) ! temporary argument / result
    real(wp), allocatable :: xv (:,:,:) ! vertical profiles
    integer  :: nx, ny, nz, nh

    !---------------
    ! set dimensions
    !---------------
    nx = covm% nx             ! Horizontal dimensions of interpolation space
    ny = covm% ny
    nz = covm% nz
    nh = covm% hc2_part_size  ! Size of local partition (our "slice")

    !-------------------------------------------------
    ! adjoint interpolation
    ! (and vertical transformation + balance operator)
    !-------------------------------------------------
    allocate (xv (nz, io% np_v, io% nc_v))
    xv(:,:,:) = 0._wp
    call apply_LvWvKIhIv_t (xv, x, e_fi, sink)

    !----------------------
    ! adjoint transposition
    !----------------------
    allocate (xh (nx,     ny,       nh   ))
    call apply_T_t (xh, xv, io)

    !-------------
    ! apply WhLh_t
    !-------------
    allocate (xd(nx,ny,nh))
    call apply_WhLh_t (xd, xh)

    !-----------
    ! apply WhLh
    !-----------
    call apply_WhLh (xd, xh)

    !----------------------
    ! transposition
    !----------------------
    call apply_T (xh, xv, io)

    !------------------------
    ! scale with B_nmc weight
    !------------------------
    xv = xv * w_nmc_b

    !-------------------------------------------------
    ! interpolation
    ! (and vertical transformation + balance operator)
    !-------------------------------------------------
    call apply_LvWvKIhIv (xv, y, e_fi, sink)

  end subroutine apply_B_ii_2d
!------------------------------------------------------------------------------
  subroutine apply_L_m (xd, a_m, cbg, lnewpl, verbose, order)
  real(wp)     ,intent(in)  :: xd (:,:,:) ! input
  type(t_cols) ,intent(out) :: a_m(:)     ! analysis increment (model space)
  type(t_cols) ,intent(in)  :: cbg(:)     ! reference state
  logical      ,intent(in)  :: lnewpl     ! analysis increment on new p-levs
  logical      ,intent(in)  :: verbose    ! verbose timing information?
  integer      ,intent(in)  :: order      ! hor.intp.: 2=linear, 4=higer order
  optional                  :: verbose, order

  !---------------------------------------------------------
  ! multiply  a vector by the L matrix
  ! result is a vector in interpolation space
  !---------------------------------------------------------
    real(wp), allocatable :: xh (:,:,:) ! temporary argument / result
    real(wp), allocatable :: xv (:,:,:) ! vertical profiles
    type (t_intop)        :: io         ! interpolation coefficients
    integer               :: nx, ny, nz, nh
    logical               :: verb

    verb = .true.; if (present (verbose)) verb = verbose
    !---------------
    ! set dimensions
    !---------------
    nx = covm% nx             ! Horizontal dimensions of interpolation space
    ny = covm% ny
    nz = covm% nz
    nh = covm% hc2_part_size  ! Size of local partition (our "slice")

    if(ltrace) call trace ('apply_L_m: x =',xd)

    !-----------
    ! apply WhLh
    !-----------
    if (verb)  call stop_time ('WhLh')
    allocate (xh (nx, ny, nh))
    call apply_WhLh (xd, xh)

    if(ltrace) call trace ('WhLh x =',xh)

    !---------------------------------------------
    ! set up horizontal interpolation coefficients
    ! and transposition
    !---------------------------------------------
    if (verb)  call stop_time ('construct (io)')
    call construct (io, cbg, order=order)

    !----------------------
    ! transposition
    !----------------------
    allocate (xv (nz,io% np_v,io% nc_v))
    if (verb)  call stop_time ('T')
    call apply_T (xh, xv, io)

    if(ltrace) call trace ('TWhLh x =',xv)

    !-------------------------------------------------
    ! interpolation
    ! (and vertical transformation + balance operator)
    !-------------------------------------------------
    if (verb)  call stop_time ('LvWvKIhIv_m')
    call apply_LvWvKIhIv_m (xv, io, a_m, cbg, lnewpl, verbose=verb)

    !---------
    ! clean up
    !---------
    if (verb)  call stop_time ('destruct (io)')
    call destruct (io)

  end subroutine apply_L_m
!------------------------------------------------------------------------------
  subroutine apply_L_i (x, y, e_fi, sink)
  real(wp)         ,intent(in)    :: x(:,:,:) ! input
  type (t_vector)  ,intent(inout) :: y        ! output
  type (t_vector)  ,intent(in)    :: e_fi     ! background error
  type (t_vector)  ,intent(in)    :: sink     ! sink variables
  optional                        :: sink
  !---------------------------------------------------------
  ! multiply  a vector by the L matrix
  ! result is a vector in interpolation space
  !---------------------------------------------------------
    real(wp), allocatable :: xh (:,:,:) ! temporary argument / result
    real(wp), allocatable :: xv (:,:,:) ! vertical profiles
    integer  :: nx, ny, nz, nh

    !---------------
    ! set dimensions
    !---------------
    nx = covm% nx             ! Horizontal dimensions of interpolation space
    ny = covm% ny
    nz = covm% nz
    nh = covm% hc2_part_size  ! Size of local partition (our "slice")

    if(ltrace) call trace ('apply_L_i: x =',x)

    !-----------
    ! apply WhLh
    !-----------
    allocate (xh (nx,     ny,       nh   ))
    call apply_WhLh (x, xh)

    if(ltrace) call trace ('WhLh x =',xh)

    !----------------------
    ! transposition
    !----------------------
    allocate (xv (nz,io% np_v,io% nc_v))
    call apply_T (xh, xv, io)

    if(ltrace) call trace ('TWhLh x =',xv)

    !-------------------------------------------------
    ! interpolation
    ! (and vertical transformation + balance operator)
    !-------------------------------------------------
    call apply_LvWvKIhIv (xv, y, e_fi, sink)

    if(ltrace) call trace ('LvWvKIhIvTWhLh x =',y)

  end subroutine apply_L_i
!------------------------------------------------------------------------------
  subroutine apply_L_i_t (x, y, e_fi, sink)
  type (t_vector)  ,intent(in)    :: x        ! input
  real(wp)         ,intent(out)   :: y(:,:,:) ! output
  type (t_vector)  ,intent(in)    :: e_fi     ! background error
  type (t_vector)  ,intent(inout) :: sink     ! sink variables
  optional                        :: sink
  !---------------------------------------------------------
  ! multiply  a vector in interpolation space by the L_t
  !---------------------------------------------------------
    real(wp), allocatable :: xh (:,:,:) ! temporary argument / result
    real(wp), allocatable :: xv (:,:,:) ! vertical profiles
    integer  :: nx, ny, nz, nh

    !---------------
    ! set dimensions
    !---------------
    nx = covm% nx             ! Horizontal dimensions of interpolation space
    ny = covm% ny
    nz = covm% nz
    nh = covm% hc2_part_size  ! Size of local partition (our "slice")

    if(ltrace) call trace ('apply_L_i_t: x =',x)

    !-------------------------------------------------
    ! adjoint interpolation
    ! (and vertical transformation + balance operator)
    !-------------------------------------------------
    allocate (xv (nz,io% np_v,io% nc_v))
    xv(:,:,:) = 0._wp
    call apply_LvWvKIhIv_t (xv, x, e_fi, sink)

    if(ltrace) call trace ('LvWvKIhIv_t x =',xv)

    !----------------------
    ! adjoint transposition
    !----------------------
    allocate (xh (nx,     ny,       nh   ))
    call apply_T_t (xh, xv, io)

    if(ltrace) call trace ('T_LvWvKIhIv_t x =',xh)

    !-------------
    ! apply WhLh_t
    !-------------
    call apply_WhLh_t (y, xh)

    if(ltrace) call trace ('LhWhT_LvWvKIhIv_t x =',y)

  end subroutine apply_L_i_t
!==============================================================================
  subroutine test_B_adj (obs, e_fi)
  type (t_obs_set) ,intent(in) :: obs
  type (t_vector)  ,intent(in) :: e_fi ! background error

    real(wp), allocatable :: x (:,:,:)
    real(wp), allocatable :: y (:,:,:)
    real(wp), allocatable :: v (:,:,:)
    real(wp), allocatable :: w (:,:,:)
    type (t_vector) :: xi, wi
    integer  :: nx, ny, nz, nh, nc, np
    real(wp) :: t_v, t_t, t_h
    logical  :: ok_v, ok_h, ok_t

    if (repr_2dh   == 0) return
    if (test_adj   == 0) return
    if (obs% oi% n == 0) return  ! no observations
    test_adj = test_adj - 1

    !---------------
    ! set dimensions
    !---------------
    nx = covm% nx             ! number of zonal      grid-points
    ny = covm% ny             ! number of meridional grid-points
    nz = covm% nz             ! number of vertical   grid-points
    nh = covm% hc2_part_size  ! Size of local partition (our "slice")
    np = io%   np_v           ! number of parameters (h,rh,psi,chi)
    nc = io%   nc_v           ! number of grid columns

    !===============
    ! test LvWvKIhIv
    !===============
    call stop_time ('test_B_adj: random_gauss')
    call construct (xi, obs% ii)
    call construct (wi, obs% ii)
    allocate (y (nz, np, nc))
    allocate (v (nz, np, nc))
    call random_gauss (xi)
    call random_gauss (v)       ! uses old random_gauss
    call stop_time ('test_B_adj: LvWvKIhIv')
    call apply_LvWvKIhIv_t (y, xi, e_fi)
    call apply_LvWvKIhIv   (v, wi, e_fi)
    t_v  = abs((sum(xi*wi)-p_sum(sum(y*v)))/sum(xi*wi))
    ok_v = t_v < 2.e-9_wp ! previously 1.e-10, maybe needs rethinking
    deallocate    (y,v)
    call destruct (xi)
    call destruct (wi)

    !=======
    ! test T
    !=======
    call stop_time ('test_B_adj: random_gauss')
    allocate (x (nz, np, nc))
    allocate (w (nz, np, nc))
    allocate (y (nx, ny, nh))
    allocate (v (nx, ny, nh))
    call random_gauss (x)       ! uses old random_gauss
    call random_gauss (v, covm)
    call stop_time ('test_B_adj: T')
    call apply_T_t (y, x, io)
    call apply_T   (v, w, io)
    t_t   = abs((p_sum(sum(x*w))-p_sum(sum(y*v)))/p_sum(sum(x*w)))
    ok_t = t_t < 1.e-10_wp
    deallocate (x,y,v,w)

    !==========
    ! test WhLh
    !==========
    call stop_time ('test_B_adj: random_gauss')
    allocate (x(nx,ny,nh))
    allocate (y(nx,ny,nh))
    allocate (v(nx,ny,nh))
    allocate (w(nx,ny,nh))
    call random_gauss (x, covm)
    call random_gauss (v, covm)
    call stop_time ('test_B_adj: WhLh')
    call apply_WhLh_t (y, x)
    call apply_WhLh   (v, w)
    t_h  = abs((sum(x*w)-sum(y*v))/sum(x*w))
    t_h  = p_max(t_h)
    ok_h = t_h < 1.e-9_wp       ! previously 1.e-10, maybe needs rethinking
    deallocate (x,y,v,w)

    !---------
    ! printout
    !---------
    if (dace% lpio) then
!     write(6,'(a)') repeat('-',79)
      write(6,'( )')
      write(6,'(a)') '  test_B_adj: Adjoint tests for B matrix'
      write(6,'(a)') '              (2d+1d separable operator approach)'
      write(6,'( )')
      if (ok_v) then
        write(6,*) 'adjoint test on LvWvKIhIv OK    :',t_v
      else
        write(6,*) 'adjoint test on LvWvKIhIv FAILED:',t_v
      endif
      if (ok_t) then
        write(6,*) 'adjoint test on T         OK    :',t_t
      else
        write(6,*) 'adjoint test on T         FAILED:',t_t
      endif
      if (ok_h) then
        write(6,*) 'adjoint test on WhLh      OK    :',t_h
      else
        write(6,*) 'adjoint test on WhLh      FAILED:',t_h
      endif
      if (.not.(ok_v.and.ok_t.and.ok_h)) &
        call finish('test_B_adj','adjoint test FAILED')
    endif

  end subroutine test_B_adj
!==============================================================================
  subroutine apply_LvWvKIhIv_t (xv, xi, e_fi, sink)
  real(wp)       ,intent(out)   :: xv(:,:,:)! result   in 'vertical mode' space
  type(t_vector) ,intent(in)    :: xi       ! argument in 'interpolation' space
  type(t_vector) ,intent(in)    :: e_fi     ! background error
  type(t_vector) ,intent(inout) :: sink     ! sink variables (output)
  optional                      :: sink
  !----------------------------------------------------------
  ! apply adjoint interpolations (+ vertical transformations)
  !----------------------------------------------------------
    integer  :: icv                         ! column index (in xv)
    integer  :: j                           ! meridional index
#if defined (__NEC__)
    integer  :: ii, jj                      ! Aux. loop indices
    integer  :: nz, n1                      ! Temporary array dim.
                                            ! Temporary arrays
#undef PREFER_ALLOCATABLE
#ifdef PREFER_ALLOCATABLE  /* Prefer allocatable temporaries for sqcXY? */
    real(wp), allocatable, dimension(:,:)         :: sqcvh, sqcvq, sqcvs, sqcvv
#else
    real(wp), dimension((covm%nz/2)*2+1,covm% nz) :: sqcvh, sqcvq, sqcvs, sqcvv
#endif
    real(wp), dimension(covm% nz)         :: xw_h, xw_q, xw_s, xw_v
#if __NEC_VERSION__ < 30300 /* work around nfort issue with -fmatrix-multiply */
!NEC$ vreg(xw_h)
!NEC$ vreg(xw_q)
!NEC$ vreg(xw_s)
!NEC$ vreg(xw_v)
#endif
#endif

    if (ltrace) call trace ('IvIhKLvWv_t: x =',xi, io% bx)

    !-----------------------------
    ! apply adjoint interpolations
    !-----------------------------
FTRACE_BEGIN("apply_LvWvKIhIv_t:IhIv_t")
!CDIR IEXPAND(apply_IhIv_t)
    call apply_IhIv_t (xv, xi, e_fi, sink)
FTRACE_END  ("apply_LvWvKIhIv_t:IhIv_t")

    !----------------------------------------------------------
    ! prepare/check vertical covariances and psi-h correlations
    !----------------------------------------------------------
FTRACE_BEGIN("apply_LvWvKIhIv_t:uncompress_cov")
    call uncompress_cov
    if (.not.associated(covm% c_h_psi)) &
      call finish('apply_LvWvKIhIv','c_h_psi is not associated !')

    if(ltrace) then
      call trace ('IvIh_t: x =',xv)
      call trace ('h   =',xv(:,IV_H,:))
      call trace ('rh  =',xv(:,IV_RH,:))
      call trace ('psi =',xv(:,IV_PSI,:))
      call trace ('chi =',xv(:,IV_CHI,:))
      if (dace% lpio) write(6,'(a32,3f15.3)')  'c_h_psi =',&
                      minval(covm% c_h_psi),maxval(covm% c_h_psi)
    endif
FTRACE_END  ("apply_LvWvKIhIv_t:uncompress_cov")

    !==============================
    ! K_t: adjoint balance operator
    !==============================
#if !defined (__NEC__)
!$omp parallel private(icv,j)
#else
    nz = covm% nz
    n1 = (nz/2)*2 + 1
!$omp parallel private(icv,j,sqcvh,sqcvq,sqcvs,sqcvv,xw_h,xw_q,xw_s,xw_v)
#ifdef PREFER_ALLOCATABLE
    allocate (sqcvh(n1,nz))
    allocate (sqcvq(n1,nz))
    allocate (sqcvs(n1,nz))
    allocate (sqcvv(n1,nz))
#endif
FTRACE_BEGIN("apply_LvWvKIhIv_t:K")
#endif

!$omp do schedule(static)
    do icv = 1, io% nc_v
      j = io% mc% c(icv)% ijdtp(2)
!NEC$ shortloop
      xv (:, IV_H,   icv) = xv (:, IV_H,   icv) &
                          + xv (:, IV_PSI, icv) * covm% c_h_psi (j)
!NEC$ shortloop
      xv (:, IV_PSI, icv) = xv (:, IV_PSI, icv) * covm% c1_h_psi(j)
    end do
!$omp end do nowait
FTRACE_END  ("apply_LvWvKIhIv_t:K")

    !=========================================
    ! LvWv_t: adjoint vertical transformations
    !=========================================
FTRACE_BEGIN("apply_LvWvKIhIv_t:LvWv_t")
!$omp do schedule(static)
    do icv = 1, io% nc_v
      j = io% mc% c(icv)% ijdtp(2)
#if defined (__NEC__)
!CDIR BEGIN SHORTLOOP
      ! Copy arrays to reduce bank conflicts
      sqcvh(1:nz,1:nz) = covm% sqcvh (:,:,j)
      sqcvq(1:nz,1:nz) = covm% sqcvq (:,:,j)
      sqcvs(1:nz,1:nz) = covm% sqcvs (:,:,j)
      sqcvv(1:nz,1:nz) = covm% sqcvv (:,:,j)
      xw_h (1:nz) = 0._wp
      xw_q (1:nz) = 0._wp
      xw_s (1:nz) = 0._wp
      xw_v (1:nz) = 0._wp
!NEC$ shortloop
      do ii = 1, nz
!NEC$ shortloop
!NEC$ ivdep
        do jj = 1, nz
          xw_h (jj) = xw_h (jj) + xv (ii, IV_H,   icv) * sqcvh(ii,jj)
          xw_q (jj) = xw_q (jj) + xv (ii, IV_RH,  icv) * sqcvq(ii,jj)
          xw_s (jj) = xw_s (jj) + xv (ii, IV_PSI, icv) * sqcvs(ii,jj)
          xw_v (jj) = xw_v (jj) + xv (ii, IV_CHI, icv) * sqcvv(ii,jj)
        end do
      end do
      xv (:, IV_H,   icv) = xw_h (1:nz)
      xv (:, IV_RH,  icv) = xw_q (1:nz)
      xv (:, IV_PSI, icv) = xw_s (1:nz)
      xv (:, IV_CHI, icv) = xw_v (1:nz)
!CDIR END
#else
      xv (:, IV_H,   icv) = xv (:, IV_H,   icv) .x. covm% sqcvh (:,:,j)
      xv (:, IV_RH,  icv) = xv (:, IV_RH,  icv) .x. covm% sqcvq (:,:,j)
      xv (:, IV_PSI, icv) = xv (:, IV_PSI, icv) .x. covm% sqcvs (:,:,j)
      xv (:, IV_CHI, icv) = xv (:, IV_CHI, icv) .x. covm% sqcvv (:,:,j)
#endif
    end do
!$omp end do
FTRACE_END  ("apply_LvWvKIhIv_t:LvWv_t")
#if defined (__NEC__)
#ifdef PREFER_ALLOCATABLE
    deallocate (sqcvh)
    deallocate (sqcvq)
    deallocate (sqcvv)
    deallocate (sqcvs)
#endif
#endif

!$omp end parallel

    call compress_cov

    if(ltrace) then
      call trace ('apply_LvWvKIhIv_t: y =',xv)
      call trace ('h   =',xv(:,IV_H,:))
      call trace ('rh  =',xv(:,IV_RH,:))
      call trace ('psi =',xv(:,IV_PSI,:))
      call trace ('chi =',xv(:,IV_CHI,:))
    endif

  end subroutine apply_LvWvKIhIv_t
!------------------------------------------------------------------------------
  subroutine apply_IhIv_t (xv, xi, e_fi, sink)
  real(wp)       ,intent(out)   :: xv(:,:,:)! result   in 'vertical mode' space
  type(t_vector) ,intent(in)    :: xi       ! argument in 'interpolation' space
  type(t_vector) ,intent(in)    :: e_fi     ! background error
  type(t_vector) ,intent(inout) :: sink     ! sink variables (output)
  optional                      :: sink
  !-----------------------------
  ! apply adjoint interpolations
  !-----------------------------
    integer  :: nc_l                    ! no. columns  at observation pts.
    integer  :: ib                      ! global 'box' index
    integer  :: ibl                     ! local  'box' index
    integer  :: ic                      ! column index (in xl)
    integer  :: icv                     ! column index (in xv)
    integer  :: id                      ! sink index
    integer  :: il                      ! parameter index (in xl)
    integer  :: j                       ! meridional index
    integer  :: l                       ! level index
    integer  :: i
    integer  :: k
    real(wp) :: su, sv                  ! u, v error
    real(wp) :: xl (covm% nz, io% np_l) ! columns at observation pts
    type(t_box) ,pointer :: bx

    xv = 0._wp
    !----------------
    ! loop over boxes
    !----------------
    do ibl = 1, io% n_box
      bx   => io% bx(ibl)
      ib   = bx% ib
      nc_l = bx% nc_l
      id   = 0

      !=====================================
      ! Iv_t: adjoint vertical interpolation
      !=====================================

      !------------------
      ! loop over columns
      !------------------
      do ic = 1, nc_l
        xl = 0._wp

        !-----------------
        ! loop over levels
        !-----------------
        do l = bx% hic(ic)% l% i + 1, bx% hic(ic)% l% i + bx% hic(ic)% l% n

          !-------------------------------------
          ! loop over observations at this level
          !-------------------------------------
          do i = bx% vic(l)% i% i + 1, &
                 bx% vic(l)% i% i + bx% vic(l)% i% n
            select case (bx% ip(i))
            case (OBS_H, OBS_HS, OBS_TV)
              il = IL_H
            case (OBS_RH)
              il = IL_RH
            case (OBS_U)
              il = IL_U
            case (OBS_V)
              il = IL_V
            case (OBS_DUM)
            case default
              write(0,*)'apply_LvWvKIhIv_t: invalid interpolation type:',&
                        bx% ip(i)
              write(0,*)'apply_LvWvKIhIv_t: box,l,i=',ib,l,i
              call finish ('apply_LvWvKIhIv_t','invalid interpolation type')
            end select
            select case (bx% ip(i))
            case (OBS_TV)
              do k = 1, bx% vic(l)% g% n
                xl  (bx% vic(l)% g% i+k-1 ,il) = &
                  xl(bx% vic(l)% g% i+k-1 ,il) - &
                      xi% s(ib)% x   (i)       * &
                    e_fi% s(ib)% x   (i)       * &
                    (bx% vic(l)% wt(k)         * &
                     bx% vic(l)% exzn          + &
                     bx% vic(l)% wh(k)         * &
                     bx% vic(l)% ezn             )
              end do
            case (OBS_DUM)
              if (present(sink)) then
                id = id + 1
                sink% s(ib)% x(id) = xi% s(ib)% x (i) * e_fi% s(ib)% x (i)
              endif
            case default
!             mr(ib)% lq(ic) = .true.
              do k = 1, bx% vic(l)% g% n
                xl (bx% vic(l)% g% i+k-1 ,il) = &
                 xl(bx% vic(l)% g% i+k-1 ,il) + &
                    bx% vic(l)% wh(k)         * &
                     xi% s(ib)% x   (i)       * &
                   e_fi% s(ib)% x   (i)
              end do
            end select
          end do

        end do

        !=======================================
        ! Ih_t: adjoint horizontal interpolation
        !=======================================

        !--------------------------------
        ! NMC-B: derive u,v from chi,psi
        !--------------------------------
        if (covm% oi_comp == 2) then
          !--------------
          ! rescale winds
          !--------------
          su = 0._wp
          sv = 0._wp
          do k = 1, size(bx% hic(ic)% imc,1)
            icv = bx% hic(ic)% imc(k,1)
            if (icv==0) cycle
            j   = io% mc% c(icv)% ijdtp(2)
            su = su     + covm% sdev_u(j) * bx% hic(ic)% w (k)
            sv = sv     + covm% sdev_v(j) * bx% hic(ic)% w (k)
          end do
          xl (:, IL_U)  = xl (:, IL_U) * (1._wp / su)
          xl (:, IL_V)  = xl (:, IL_V) * (1._wp / sv)
        endif

        do k = 1, size(bx% hic(ic)% imc,1)
          icv = bx% hic(ic)% imc(k,1)
          if (icv==0) cycle
          j   = io% mc% c(icv)% ijdtp(2)

          xv   (:, IV_H,   icv) = &
            xv (:, IV_H,   icv) + &
            xl (:, IL_H)        * &
            bx% hic(ic)% w (k)

          xv   (:, IV_RH,  icv) = &
            xv (:, IV_RH,  icv) + &
            xl (:, IL_RH)       * &
            bx% hic(ic)% w (k)

          xv   (:, IV_CHI, icv) = &
            xv (:, IV_CHI, icv) + &
           (xl (:, IL_U)        * &
            bx% hic(ic)% w1(k)  + &
            xl (:, IL_V)        * &
            bx% hic(ic)% w2(k) )* &
            covm% sqnu(j)

          xv   (:, IV_PSI, icv) = &
            xv (:, IV_PSI, icv) - &
           (xl (:, IL_U)        * &
            bx% hic(ic)% w2(k)  - &
            xl (:, IL_V)        * &
           bx% hic(ic)% w1(k)) * &
             covm% sq1nu
        end do

        !-------------------------
        ! end of loop over columns
        !-------------------------
      end do

      !-----------------------
      ! end of loop over boxes
      !-----------------------
    end do
  end subroutine apply_IhIv_t
!------------------------------------------------------------------------------
  subroutine apply_LvWvKIhIv (xv, xi, e_fi, sink)
  real(wp)       ,intent(in)    :: xv(:,:,:)! argument in 'vertical mode' space
  type(t_vector) ,intent(inout) :: xi       ! result   in 'interpolation' space
  type (t_vector),intent(in)    :: e_fi     ! background error
  type(t_vector) ,intent(in)    :: sink     ! sink variables (input)
  optional                      :: sink
  !--------------------------------------------------
  ! apply interpolations (+ vertical transformations)
  !--------------------------------------------------
    integer  :: icv                         ! column index (in xv)
    integer  :: j                           ! meridional index
    real(wp) :: xw (size(xv,1),size(xv,2),size(xv,3)) ! working array
#if defined (__NEC__)
    integer  :: ii, jj
#endif

    if(ltrace) then
      call trace ('apply_LvWvKIhIv: x =',xv)
      call trace ('h   =',xv(:,IV_H,:))
      call trace ('rh  =',xv(:,IV_RH,:))
      call trace ('psi =',xv(:,IV_PSI,:))
      call trace ('chi =',xv(:,IV_CHI,:))
    endif

    !----------------------------------------------------------
    ! prepare/check vertical covariances and psi-h correlations
    !----------------------------------------------------------
    call uncompress_cov
    if (.not.associated(covm% c_h_psi)) &
      call finish('apply_LvWvKIhIv','c_h_psi is not associated !')

!$omp parallel do private(icv,j) schedule(static)
    do icv = 1, io% nc_v
      !===============================
      ! LvWv: vertical transformations
      !===============================
      j = io% mc% c(icv)% ijdtp(2)
#if defined (__NEC__)
!NEC$ shortloop
      xw (:, IV_H,   icv) = 0.
!NEC$ shortloop
      xw (:, IV_RH,  icv) = 0.
!NEC$ shortloop
      xw (:, IV_PSI, icv) = 0.
!NEC$ shortloop
      xw (:, IV_CHI, icv) = 0.
!NEC$ shortloop
      do ii = 1, covm% nz
!NEC$ shortloop
!NEC$ ivdep
        do jj = 1, covm% nz
          xw (jj, IV_H,   icv) = xw (jj, IV_H,   icv) &
                               + covm% sqcvh(jj,ii,j) * xv (ii, IV_H,   icv)
          xw (jj, IV_RH,  icv) = xw (jj, IV_RH,  icv) &
                               + covm% sqcvq(jj,ii,j) * xv (ii, IV_RH,  icv)
          xw (jj, IV_PSI, icv) = xw (jj, IV_PSI, icv) &
                               + covm% sqcvs(jj,ii,j) * xv (ii, IV_PSI, icv)
          xw (jj, IV_CHI, icv) = xw (jj, IV_CHI, icv) &
                               + covm% sqcvv(jj,ii,j) * xv (ii, IV_CHI, icv)
        end do
      end do
#else
      xw (:, IV_H,   icv) = covm% sqcvh (:,:,j) .x. xv (:, IV_H,   icv)
      xw (:, IV_RH,  icv) = covm% sqcvq (:,:,j) .x. xv (:, IV_RH,  icv)
      xw (:, IV_PSI, icv) = covm% sqcvs (:,:,j) .x. xv (:, IV_PSI, icv)
      xw (:, IV_CHI, icv) = covm% sqcvv (:,:,j) .x. xv (:, IV_CHI, icv)
#endif
      !====================
      ! K: balance operator
      !====================
      j = io% mc% c(icv)% ijdtp(2)
      xw (:, IV_PSI, icv) = covm% c_h_psi (j) * xw (:, IV_H,   icv) &
                          + covm% c1_h_psi(j) * xw (:, IV_PSI, icv)
    end do
!$omp end parallel do
    call compress_cov

    if(ltrace) then
      call trace ('KLvWv: x =',xw)
      call trace ('h   =',xw(:,IV_H,:))
      call trace ('rh  =',xw(:,IV_RH,:))
      call trace ('psi =',xw(:,IV_PSI,:))
      call trace ('chi =',xw(:,IV_CHI,:))
      if (dace% lpio) write(6,'(a32,3f15.3)')  'c_h_psi =',&
                      minval(covm% c_h_psi),maxval(covm% c_h_psi)
    endif

    !---------------------
    ! apply interpolations
    !---------------------
!CDIR IEXPAND(apply_IhIv)
    call apply_IhIv  (xw, xi, e_fi, sink)

    if (ltrace) call trace ('IvIhKLvWv x =',xi, io% bx)

  end subroutine apply_LvWvKIhIv
!------------------------------------------------------------------------------
  subroutine apply_IhIv (xv, xi, e_fi, sink)
  real(wp)       ,intent(in)    :: xv(:,:,:)! argument in 'vertical mode' space
  type(t_vector) ,intent(inout) :: xi       ! result   in 'interpolation' space
  type(t_vector) ,intent(in)    :: e_fi     ! background error
  type(t_vector) ,intent(in)    :: sink     ! sink variables (input)
  optional                      :: sink
  !---------------------
  ! apply interpolations
  !---------------------
    integer  :: nc_l                    ! no. columns  at observation pts.
    integer  :: ib                      ! global 'box' index
    integer  :: ibl                     ! local  'box' index
    integer  :: ic                      ! column index (in xl)
    integer  :: icv                     ! column index (in xv)
    integer  :: id                      ! sink index
    integer  :: il                      ! parameter index (in xl)
    integer  :: j                       ! meridional index
    integer  :: l                       ! level index
    integer  :: i
    integer  :: k
    real(wp) :: su, sv                  ! u, v error
    real(wp) :: xl (covm% nz, io% np_l) ! columns at observation pts
    type(t_box) ,pointer :: bx

    xi  = 0._wp
    !----------------
    ! loop over boxes
    !----------------
    do ibl = 1, io% n_box
      bx   => io% bx(ibl)
      ib   = bx% ib
      nc_l = bx% nc_l
      id   = 0
      !------------------
      ! loop over columns
      !------------------
      do ic = 1, nc_l
        xl = 0._wp
        !=============================
        ! Ih: horizontal interpolation
        !=============================

        !--------------------------------
        ! NMC-B: derive u,v from chi,psi
        !--------------------------------
        do k = 1, size(bx% hic(ic)% imc,1)
          icv = bx% hic(ic)% imc(k,1)
          if (icv==0) cycle
          j   = io% mc% c(icv)% ijdtp(2)

          xl (:, IL_H)  = xl (:, IL_H)       &
                        + xv (:, IV_H,  icv) &
                        * bx% hic(ic)% w (k)

          xl (:, IL_RH) = xl (:, IL_RH)      &
                        + xv (:, IV_RH, icv) &
                        * bx% hic(ic)% w (k)

          xl (:, IL_U)  = xl (:, IL_U)        &
                        + xv (:, IV_CHI, icv) &
                        * bx% hic(ic)% w1(k)  &
                        * covm% sqnu(j)       &
                        - xv (:, IV_PSI, icv) &
                        * bx% hic(ic)% w2(k)  &
                        * covm% sq1nu

          xl (:, IL_V)  = xl (:, IL_V)        &
                        + xv (:, IV_CHI, icv) &
                        * bx% hic(ic)% w2(k)  &
                        * covm% sqnu(j)       &
                        + xv (:, IV_PSI, icv) &
                        * bx% hic(ic)% w1(k)  &
                        * covm% sq1nu
        end do

        if (covm% oi_comp == 2) then
          !--------------
          ! rescale winds
          !--------------
          su = 0._wp
          sv = 0._wp
          do k = 1, size(bx% hic(ic)% imc,1)
            icv = bx% hic(ic)% imc(k,1)
            if (icv==0) cycle
            j   = io% mc% c(icv)% ijdtp(2)
            su = su     + covm% sdev_u(j) * bx% hic(ic)% w (k)
            sv = sv     + covm% sdev_v(j) * bx% hic(ic)% w (k)
          end do
          xl (:, IL_U)  = xl (:, IL_U) * (1._wp / su)
          xl (:, IL_V)  = xl (:, IL_V) * (1._wp / sv)
        endif

        !===========================
        ! Iv: vertical interpolation
        !===========================
        !-----------------
        ! loop over levels
        !-----------------
        do l = bx% hic(ic)% l% i + 1, bx% hic(ic)% l% i + bx% hic(ic)% l% n

          !-------------------------------------
          ! loop over observations at this level
          !-------------------------------------
          do i = bx% vic(l)% i% i + 1, &
                 bx% vic(l)% i% i + bx% vic(l)% i% n
            select case (bx% ip(i))
            case (OBS_H, OBS_HS, OBS_TV)
              il = IL_H
            case (OBS_RH)
              il = IL_RH
            case (OBS_U)
              il = IL_U
            case (OBS_V)
              il = IL_V
            case (OBS_DUM)
            case default
              write(0,*)'apply_LvWvKIhIv: invalid interpolation type:',&
                        bx% ip(i)
              write(0,*)'apply_LvWvKIhIv: box,l,i=',ib,l,i
              call finish ('apply_LvWvKIhIv','invalid interpolation type')
            end select
            select case (bx% ip(i))
            case (OBS_TV)
              do k = 1, bx% vic(l)% g% n
                xi% s(ib)% x(i) = xi% s(ib)% x(i) - &
                  e_fi% s(ib)% x   (i)            * &
                  xl(bx% vic(l)% g% i+k-1 ,il)    * &
                    (bx% vic(l)% wt(k)            * &
                     bx% vic(l)% exzn             + &
                     bx% vic(l)% wh(k)            * &
                     bx% vic(l)% ezn                )
              end do
            case (OBS_DUM)
              if (present(sink)) then
                id = id + 1
                xi% s(ib)% x (i) = e_fi% s(ib)% x (i) * sink% s(ib)% x(id)
              endif
            case default
!             ml% lq(levs(l)% is) = .true.
              do k = 1, bx% vic(l)% g% n
                xi% s(ib)% x(i) = xi% s(ib)% x(i) + &
                  xl(bx% vic(l)% g% i+k-1 ,il)    * &
                  bx% vic(l)% wh(k)               * &
                  e_fi% s(ib)% x(i)
              end do
            end select
          end do
        end do

        !-------------------------
        ! end of loop over columns
        !-------------------------
      end do

      !-----------------------
      ! end of loop over boxes
      !-----------------------
    end do
  end subroutine apply_IhIv
!------------------------------------------------------------------------------
  subroutine apply_LvWvKIhIv_m (xv, io, a_m, cbg, lnewpl, verbose)
  real(wp)     ,intent(in)    :: xv (:,:,:) ! input, vertical profiles
  type(t_intop),intent(inout) :: io         ! interpolation coefficients
  type(t_cols) ,intent(out)   :: a_m(:)     ! analysis increment (model space)
  type(t_cols) ,intent(in)    :: cbg(:)     ! reference state
  logical      ,intent(in)    :: lnewpl     ! analysis increment on new p-levs
  logical      ,intent(in)    :: verbose    ! verbose timing information?
  optional                    :: verbose
  !-------------------------------------------------
  ! interpolation
  ! (and vertical transformation + balance operator)
  ! result is on model grid
  !-------------------------------------------------
    real(wp)              :: xw (size(xv,1), &! working array
                                 size(xv,2), &!
                                 size(xv,3))  !
    !-------------------------------------------
    ! vertical transformation + balance operator
    !-------------------------------------------
    call apply_LvWvK_m (xv, xw, io, verbose)
    !---------------------------------------
    ! interpolation, result is on model grid
    !---------------------------------------
    call apply_IhIv_m  (xw, io, a_m, cbg, lnewpl, verbose)
  end subroutine apply_LvWvKIhIv_m
!------------------------------------------------------------------------------
  subroutine apply_LvWvK_m (xv, xw, io, verbose)
  real(wp)     ,intent(in)    :: xv (:,:,:) ! input, vertical profiles
  real(wp)     ,intent(out)   :: xw (:,:,:) !
  type(t_intop),intent(in)    :: io         ! interpolation coefficients
  logical      ,intent(in)    :: verbose    ! verbose timing information?
  optional                    :: verbose
  !-------------------------------------------------
  ! vertical transformation + balance operator
  ! for result on model grid
  !-------------------------------------------------
    integer            :: icv                 ! column index (in xv)
    integer            :: j                   ! meridional index
#if defined (__NEC__)
    integer :: ii, jj
#endif
    logical            :: verb

    verb = .true.; if (present (verbose)) verb = verbose

    if(ltrace) then
      call trace ('apply_LvWvKIhIv_m: x =',xv)
      call trace ('h   =',xv(:,IV_H,:))
      call trace ('rh  =',xv(:,IV_RH,:))
      call trace ('psi =',xv(:,IV_PSI,:))
      call trace ('chi =',xv(:,IV_CHI,:))
    endif

    !----------------------------------------------------------
    ! prepare/check vertical covariances and psi-h correlations
    !----------------------------------------------------------
    call uncompress_cov
    if (.not.associated(covm% c_h_psi)) &
      call finish('apply_LvWvKIhIv','c_h_psi is not associated !')

    if (verb)  call stop_time ('KLvWv')
!$omp parallel do private(icv,j) schedule(static)
    do icv = 1, io% nc_v
      !===============================
      ! LvWv: vertical transformations
      !===============================
      j = io% mc% c(icv)% ijdtp(2)
#if defined (__NEC__)
!NEC$ shortloop
      xw (:, IV_H,   icv) = 0.
!NEC$ shortloop
      xw (:, IV_RH,  icv) = 0.
!NEC$ shortloop
      xw (:, IV_PSI, icv) = 0.
!NEC$ shortloop
      xw (:, IV_CHI, icv) = 0.
!NEC$ shortloop
      do ii = 1, covm% nz
!NEC$ shortloop
!NEC$ ivdep
        do jj = 1, covm% nz
          xw (jj, IV_H,   icv) = xw (jj, IV_H,   icv) &
                               + covm% sqcvh(jj,ii,j) * xv (ii, IV_H,   icv)
          xw (jj, IV_RH,  icv) = xw (jj, IV_RH,  icv) &
                               + covm% sqcvq(jj,ii,j) * xv (ii, IV_RH,  icv)
          xw (jj, IV_PSI, icv) = xw (jj, IV_PSI, icv) &
                               + covm% sqcvs(jj,ii,j) * xv (ii, IV_PSI, icv)
          xw (jj, IV_CHI, icv) = xw (jj, IV_CHI, icv) &
                               + covm% sqcvv(jj,ii,j) * xv (ii, IV_CHI, icv)
        end do
      end do
#else
      xw (:, IV_H,   icv) = covm% sqcvh (:,:,j) .x. xv (:, IV_H,   icv)
      xw (:, IV_RH,  icv) = covm% sqcvq (:,:,j) .x. xv (:, IV_RH,  icv)
      xw (:, IV_PSI, icv) = covm% sqcvs (:,:,j) .x. xv (:, IV_PSI, icv)
      xw (:, IV_CHI, icv) = covm% sqcvv (:,:,j) .x. xv (:, IV_CHI, icv)
#endif

      !====================
      ! K: balance operator
      !====================
      j = io% mc% c(icv)% ijdtp(2)
      xw (:, IV_PSI, icv) = covm% c_h_psi (j) * xw (:, IV_H,   icv) &
                          + covm% c1_h_psi(j) * xw (:, IV_PSI, icv)
    end do
!$omp end parallel do

    if(ltrace) then
      call trace ('KLvWv_m: x =',xw)
      call trace ('h   =',xw(:,IV_H,:))
      call trace ('rh  =',xw(:,IV_RH,:))
      call trace ('psi =',xw(:,IV_PSI,:))
      call trace ('chi =',xw(:,IV_CHI,:))
      if (dace% lpio) write(6,'(a32,3f15.3)')  'c_h_psi =',&
                      minval(covm% c_h_psi),maxval(covm% c_h_psi)
    endif

  end subroutine apply_LvWvK_m
!------------------------------------------------------------------------------
  subroutine apply_IhIv_m (xw, io, a_m, cbg, lnewpl, verbose)
  real(wp)     ,intent(in)    :: xw (:,:,:) ! input, vertical profiles
  type(t_intop),intent(inout) :: io         ! interpolation coefficients
  type(t_cols) ,intent(out)   :: a_m(:)     ! analysis increment (model space)
  type(t_cols) ,intent(in)    :: cbg(:)     ! reference state
  logical      ,intent(in)    :: lnewpl     ! analysis increment on new p-levs
  logical      ,intent(in)    :: verbose    ! verbose timing information?
  optional                    :: verbose
  !-------------------------------------------------
  ! interpolation
  ! result is on model grid
  !-------------------------------------------------
    integer ,parameter :: nm  = 4             ! number of multi-level fields
    integer ,parameter :: ns  = 1             ! number of single-level fields
    integer            :: nvc = 0             ! number of variables/column
    integer            :: nvcf= 0             ! number of variables, full level
    integer            :: ke                  ! number of model levels
    integer            :: icv                 ! column index (in xv)
    integer            :: ic                  ! column index (in xl)
    integer            :: nc_l                ! no. columns at observation pts.
    integer            :: ib                  ! global 'box' index
    integer            :: j                   ! meridional index
    integer            :: k                   !
    integer            :: l                   !
    integer            :: i                   !
    integer            :: il                  !
    real(wp)           :: ehs                 ! surface geopotential obs.error
    real(wp)           :: ps                  ! surface pressure
    real(wp)           :: su, sv              ! u, v error
    real(wp)           :: xl (covm% nz,      &!
                              io% np_l)       ! columns at observation pts
    integer  ,allocatable :: t_int (:)        !
    real(wp) ,allocatable :: e     (:)        !
    real(wp) ,allocatable :: z     (:)
    real(wp) ,allocatable :: y     (:)
    real(wp) ,allocatable :: ph    (:)
    real(wp) ,allocatable :: lpf   (:)
    type (t_box) ,pointer :: bx               !
    logical               :: verb

    verb = .true.; if (present (verbose)) verb = verbose

    if (verb)  call stop_time ('IhIv')
    !----------------
    ! loop over boxes
    !----------------
    do ib  = 1, size (cbg)
      bx   => io% bx(ib)
      nc_l = bx% nc_l

      !--------------------------------------
      ! allocate result variable, set to zero
      !--------------------------------------
      call alloc_cols (a_m(ib), tmp=cbg(ib), ids=COL_T+COL_UV+COL_RH+COL_GEOH)
      a_m(ib) = 0._wp

      !---------------------
      ! allocate temporaries
      !---------------------
      if (nvc == 0) then
        ke  = a_m(ib)% ke
        nvcf = nm * ke + ns
        nvc  = nvcf + ke
        allocate (e     (nvc ))
        allocate (t_int (nvc ))
        allocate (z     (nvc ))
        allocate (y     (nvc ))
        allocate (ph    (ke+1))
        allocate (lpf   (ke  ))

        !-----------------
        ! set up meta data
        !-----------------
        t_int (1        ) = OBS_H
        t_int (2:nvcf:nm) = OBS_TV
        t_int (3:nvcf:nm) = OBS_RH
        t_int (4:nvcf:nm) = OBS_U
        t_int (5:nvcf:nm) = OBS_V
        t_int (1+nvcf:  ) = OBS_H
      endif

      !------------------
      ! loop over columns
      !------------------
      do ic = 1, nc_l
        xl = 0._wp
        su = 0._wp
        sv = 0._wp

        !-------------------------
        ! horizontal interpolation
        !-------------------------
        do k = 1, size(bx% hic(ic)% imc,1)
          icv = bx% hic(ic)% imc(k,1)
          if (icv==0) cycle
          j   = io% mc% c(icv)% ijdtp(2)

!NEC$ shortloop
          xl (:, IL_H)  = xl (:, IL_H)       &
                        + xw (:, IV_H,  icv) &
                        * bx% hic(ic)% w (k)

!NEC$ shortloop
          xl (:, IL_RH) = xl (:, IL_RH)      &
                        + xw (:, IV_RH, icv) &
                        * bx% hic(ic)% w (k)

!NEC$ shortloop
          xl (:, IL_U)  = xl (:, IL_U)        &
                        + xw (:, IV_CHI, icv) &
                        * bx% hic(ic)% w1(k)  &
                        * covm% sqnu(j)       &
                        - xw (:, IV_PSI, icv) &
                        * bx% hic(ic)% w2(k)  &
                        * covm% sq1nu

!NEC$ shortloop
          xl (:, IL_V)  = xl (:, IL_V)        &
                        + xw (:, IV_CHI, icv) &
                        * bx% hic(ic)% w2(k)  &
                        * covm% sqnu(j)       &
                        + xw (:, IV_PSI, icv) &
                        * bx% hic(ic)% w1(k)  &
                        * covm% sq1nu

          if (covm% oi_comp == 2) then
            su = su     + covm% sdev_u(j) * bx% hic(ic)% w (k)
            sv = sv     + covm% sdev_v(j) * bx% hic(ic)% w (k)
          endif
        end do
        if (covm% oi_comp == 2) then
          !--------------
          ! rescale winds
          !--------------
!NEC$ shortloop
          xl (:, IL_U)  = xl (:, IL_U) * (1._wp / su)
!NEC$ shortloop
          xl (:, IL_V)  = xl (:, IL_V) * (1._wp / sv)
        endif

        if (lnewpl) then
          !---------------------------------------------------------------
          ! set up vertical interpolation coefficient for surface pressure
          !---------------------------------------------------------------
          call set_vic_ps (bx, ehs, cbg(ib)%col(ic))

          !-----------------------------
          ! interpolate surface pressure
          !-----------------------------

          ps = 0
          do k = 1, bx% vic(1)% g% n
            ps = ps + xl(bx% vic(1)% g% i+k-1 ,IL_H)  * &
                         bx% vic(1)% wh(k)            * &
                         ehs
          end do
          a_m(ib)% col(ic)% s% ps = &
          ps * (gacc/R) * (cbg(ib)% col(ic)% s% ps / cbg(ib)% col(ic)% s% t2m)
          a_m(ib)% col(ic)% s% psr = cbg(ib)% col(ic)% s% ps &
                                   + a_m(ib)% col(ic)% s% ps
        else
          a_m(ib)% col(ic)% s% psr = cbg(ib)% col(ic)% s% ps
        endif

        !----------------------------------------------------------
        ! set up vertical interpolation coefficients for atmosphere
        !----------------------------------------------------------
        select case (a_m(ib)% levtyp)
        case (WMO3_ISOBARIC)
          !------------------------------------------------
          ! pressure levels: same for t,rh,u,v, geop.height
          !------------------------------------------------
!NEC$ ivdep
          do k = 1, ke
            lpf (k) = log (a_m(ib)% ak(k))
            z (nm*(k-1)+2:nm*k+1) = lpf(k)     ! t,rh,u,v
            z (nvcf+k)            = lpf(k)     ! geop.height
          end do
        case default
          !------------------------------------
          ! GME/HRM hybrid pressure coordinates
          !------------------------------------
          if (a_m(ib)% vctype == VCT_P_HYB) then
            !----------------------------------------
            ! hybrid levels: full levels for t,rh,u,v
            !----------------------------------------
            do k = 1, ke + 1
              ph (k) = a_m(ib)% ak(k) + a_m(ib)% bk(k) * a_m(ib)% col(ic)% s% psr
            end do
            do k = 1, ke
              lpf (k) = log(0.5_wp * (ph (k) + ph (k+1)))
              z (nm*(k-1)+2:nm*k+1) = lpf(k)
            end do
            !----------------------------
            ! half levels for geop.height
            !----------------------------
            do k = 1, ke
              z (nvcf+k) = log (ph (k+1))
            end do
          !---------------------------++++++++++++++++++
          ! COSMO hybrid z-coordinates (and ICON so far)
          !---------------------------++++++++++++++++++
          else
!NEC$ ivdep
            do k = 1, ke
              lpf (k) = log (cbg(ib)% col(ic)% p(k))
              z (nm*(k-1)+2:nm*k+1) = lpf(k)     ! t,rh,u,v
              z (nvcf+k)            = lpf(k)     ! geop.height
            end do
          endif
        end select
        z (1) = log (a_m(ib)% col(ic)% s% psr)
        call set_vic_atm (bx, e, z, t_int, a_m(ib)%col(ic))

        !------------------------------------------------
        ! interpolate atmospheric parameters (vertically)
        !------------------------------------------------

        y = 0._wp
        bx% hic(ic)% l% n = 2 * ke + 1  ! PS + 2 columns
        do l = bx% hic(ic)% l% i + 1, bx% hic(ic)% l% i + bx% hic(ic)% l% n
          do i = bx% vic(l)% i% i + 1, &
                 bx% vic(l)% i% i + bx% vic(l)% i% n
            select case (t_int(i))
            case (OBS_H, OBS_HS, OBS_TV)
              il = IL_H
            case (OBS_RH)
              il = IL_RH
            case (OBS_U)
              il = IL_U
            case (OBS_V)
              il = IL_V
            case default
              call finish ('apply_LvWvKIhIv_m','invalid interpolation type')
            end select
            select case (t_int(i))
            case (OBS_TV)
!NEC$ shortloop
              do k = 1, bx% vic(l)% g% n
                 y(i) = y(i) - e(i)                  * &
                        xl(bx% vic(l)% g% i+k-1, il) * &
                          (bx% vic(l)% wt(k)         * &
                           bx% vic(l)% exzn          + &
                           bx% vic(l)% wh(k)         * &
                           bx% vic(l)% ezn           )
              end do
            case default
!NEC$ shortloop
               do k = 1, bx% vic(l)% g% n
                 y(i) = y(i) + e (i)                 * &
                        xl(bx% vic(l)% g% i+k-1, il) * &
                        bx% vic(l)% wh(k)
               end do  ! k
             end select
          end do       ! i
        end do         ! l

        !-------------
        ! store result
        !-------------
        if (.not.lnewpl) &
          a_m(ib)% col(ic)% s% ps = y (1) * (gacc/R) * &
         (cbg(ib)% col(ic)% s% ps / cbg(ib)% col(ic)% s% t2m)
        a_m(ib)%   col(ic)%    t           = y (2:nvcf:nm)
        a_m(ib)%   col(ic)%    rh          = y (3:nvcf:nm)
        a_m(ib)%   col(ic)%    u           = y (4:nvcf:nm)
        a_m(ib)%   col(ic)%    v           = y (5:nvcf:nm)
        a_m(ib)%   col(ic)%    geoh(1)     = 0._wp
        a_m(ib)%   col(ic)%    geoh(2:ke+1)= y (1+nvcf:  ) * gacc

      end do           ! ic

      !-----------------------
      ! deallocate temporaries
      !-----------------------
      if (nvc /= 0) then
        nvc = 0
        deallocate (e    )
        deallocate (t_int)
        deallocate (z    )
        deallocate (y    )
        deallocate (ph   )
        deallocate (lpf  )
      endif

    end do             ! ib

    if (verb)  call stop_time ('compress_cov')
    call compress_cov

  end subroutine apply_IhIv_m
!------------------------------------------------------------------------------
  subroutine apply_T_t (xh, xv, ipc)
  real(wp)       ,intent(out) :: xh (:,:,:)     ! regular grid
  real(wp)       ,intent(in)  :: xv (:,:,:)     ! columns near observation pts.
  type(t_intop)  ,intent(in)  :: ipc            ! interpolation column meta data
  !--------------------------------------------
  ! transposition: horizontal slices <- columns
  !--------------------------------------------
    integer               :: sc (dace% npe)     ! send counts
    integer               :: rc (dace% npe)     ! receive counts
    real(wp) ,allocatable :: sendbuf (:)        ! send    buffer
    real(wp) ,allocatable :: recvbuf (:)        ! receive buffer
    integer               :: pe                 ! processor element index
    integer               :: ic                 ! grid column index (in xv)
    integer               :: i0                 ! offset in send/receive buffer
    integer               :: is                 ! Aux. offset in send buffer
    integer               :: i, j               ! horizontal grid indices
    integer               :: nl                 ! number of slices
    integer               :: ip                 ! number of parameter
    integer               :: lb, ub             ! level bounds
    integer               :: lbl, ubl           ! slice bounds
#if defined (__NEC__)
    real(wp) ,allocatable :: xh_(:,:,:)         ! Temporary copy of xh
    integer               :: m1,m2,m3,n1,n2     ! Dimensions of xh, xh_
    integer               :: off(0:dace% npe-1) ! Auxiliary offsets
#endif
    !---------------------------------
    ! count elements to send / receive
    !---------------------------------
    do pe = 0, dace% npe-1
      sc (pe+1) = sum(ipc% pe(      pe)% nl) * ipc% pe(dace% pe)% nc
      rc (pe+1) = sum(ipc% pe(dace% pe)% nl) * ipc% pe(      pe)% nc
    end do
    !-------------------------------
    ! allocate send / receive buffer
    !-------------------------------
    allocate (sendbuf (sum (sc)))
    allocate (recvbuf (sum (rc)))
    !-----------------
    ! fill send buffer
    !-----------------
FTRACE_BEGIN("apply_T_t:set_sendbuf")
    i0 = 0
    do pe= 0, dace% npe-1
      nl = sum(ipc% pe(pe)% nl (:))
!NEC$ ivdep
      do ip = 1, ipc% np_v
        if (ipc% pe(pe)% nl (ip) > 0) then
          lb = ipc% pe(pe)% lb (ip)
          ub = ipc% pe(pe)% ub (ip)
          lbl= ipc% pe(pe)% lbl(ip)
          ubl= ipc% pe(pe)% ubl(ip)
!NEC$ ivdep
          do ic = 1, ipc% nc_v
            is = i0 + (ic-1)*nl
!NEC$ ivdep
            sendbuf (is+lbl:is+ubl) = xv (lb:ub,ip,ic)
          end do
        endif
      end do
      i0 = i0 + ipc% nc_v * nl
    end do
FTRACE_END  ("apply_T_t:set_sendbuf")
    !---------------
    ! send / receive
    !---------------
FTRACE_BEGIN("apply_T_t:p_alltoall")
    call p_alltoall (sendbuf, recvbuf, sendcounts=sc, recvcounts=rc)
FTRACE_END  ("apply_T_t:p_alltoall")
    !-----------------------------
    ! get data from receive buffer
    !-----------------------------
    nl = sum(ipc% pe(dace% pe)% nl)
#if defined (__NEC__)
    !-------------------------------------------------------
    ! Compute offsets needed for optimization of inner loops
    !-------------------------------------------------------
    off(0) = 0
    do pe = 0, dace% npe-2
       off(pe+1) = off(pe) + ipc% pe(pe)% nc * nl
    end do
    !------------------------------------------------------
    ! Reduce bank conflicts on NEC SX using temporary array
    !------------------------------------------------------
    m1 = size (xh,1)
    m2 = size (xh,2)
    m3 = size (xh,3)
    n1 = (m1/2)*2 + 1
    n2 = (m2/2)*2 + 1
    allocate (xh_(n1,n2,m3))
    xh_ = 0._wp
FTRACE_BEGIN("apply_T_t:get_recvbuf")
!NEC$ novector
    do pe = 0, dace% npe-1
!NEC$ ivdep
      do ic = 1, ipc% pe(pe)% nc
        i = ipc% pe(pe)% ij (1,ic)
        j = ipc% pe(pe)% ij (2,ic)
        i0 = off(pe) + (ic-1)*nl
        xh_(i,j,:) = xh_(i,j,:) + recvbuf (i0+1 : i0+nl)
      end do
    end do
FTRACE_END  ("apply_T_t:get_recvbuf")
    xh = xh_(1:m1,1:m2,1:m3)
#else
    xh = 0._wp
    i0 = 0
    do pe= 0, dace% npe-1
      do ic = 1, ipc% pe(pe)% nc
        i = ipc% pe(pe)% ij (1,ic)
        j = ipc% pe(pe)% ij (2,ic)
        xh(i,j,:) = xh(i,j,:) + recvbuf (i0+1 : i0+nl)
        i0 = i0 + nl
      end do
    end do
#endif
  end subroutine apply_T_t
!------------------------------------------------------------------------------
  subroutine apply_T (xh, xv, ipc)
  real(wp)       ,intent(in)  :: xh(:,:,:)      ! regular grid
  real(wp)       ,intent(out) :: xv(:,:,:)      ! columns near observation pts.
  type(t_intop)  ,intent(in)  :: ipc            ! interpolation column meta data
  !--------------------------------------------
  ! transposition: horizontal slices -> columns
  !--------------------------------------------
    integer               :: sc (dace% npe)     ! send counts
    integer               :: rc (dace% npe)     ! receive counts
    real(wp) ,allocatable :: sendbuf (:)        ! send    buffer
    real(wp) ,allocatable :: recvbuf (:)        ! receive buffer
    integer               :: pe                 ! processor element index
    integer               :: ic                 ! grid column index (in xv)
    integer               :: i0                 ! offset in send/receive buffer
    integer               :: ir                 ! Aux. offset in receive buffer
    integer               :: i, j               ! horizontal grid indices
    integer               :: nl                 ! number of slices
    integer               :: ip                 ! number of parameter
    integer               :: lb, ub             ! level bounds
    integer               :: lbl, ubl           ! slice bounds
#if defined (__NEC__)
    real(wp) ,allocatable :: xh_(:,:,:)         ! Temporary copy of xh
    integer               :: m1,m2,m3,n1,n2     ! Dimensions of xh, xh_
    integer               :: off(0:dace% npe-1) ! Auxiliary offsets
#endif
    !---------------------------------
    ! count elements to send / receive
    !---------------------------------
    do pe = 0, dace% npe-1
      sc (pe+1) = sum(ipc% pe(dace% pe)% nl) * ipc% pe(      pe)% nc
      rc (pe+1) = sum(ipc% pe(      pe)% nl) * ipc% pe(dace% pe)% nc
    end do
    !-------------------------------
    ! allocate send / receive buffer
    !-------------------------------
    allocate (sendbuf (sum (sc)))
    allocate (recvbuf (sum (rc)))
    !-----------------
    ! fill send buffer
    !-----------------
    nl = sum(ipc% pe(dace% pe)% nl)
#if defined (__NEC__)
    !-------------------------------------------------------
    ! Compute offsets needed for optimization of inner loops
    !-------------------------------------------------------
    off(0) = 0
    do pe = 0, dace% npe-2
       off(pe+1) = off(pe) + ipc% pe(pe)% nc * nl
    end do
    !------------------------------------------------------
    ! Reduce bank conflicts on NEC SX using temporary array
    !------------------------------------------------------
    m1 = size (xh,1)
    m2 = size (xh,2)
    m3 = size (xh,3)
    n1 = (m1/2)*2 + 1
    n2 = (m2/2)*2 + 1
    allocate (xh_(n1,n2,m3))
    xh_(1:m1,1:m2,1:m3) = xh(1:m1,1:m2,1:m3)
FTRACE_BEGIN("apply_T:set_sendbuf")
!NEC$ novector
    do pe = 0, dace% npe-1
!NEC$ ivdep
      do ic = 1, ipc% pe(pe)% nc
        i = ipc% pe(pe)% ij (1,ic)
        j = ipc% pe(pe)% ij (2,ic)
        i0 = off(pe) + (ic-1)*nl
        sendbuf (i0+1 : i0+nl) = xh_(i,j,:)
      end do
    end do
FTRACE_END  ("apply_T:set_sendbuf")
#else
    i0 = 0
    do pe= 0, dace% npe-1
      do ic = 1, ipc% pe(pe)% nc
        i = ipc% pe(pe)% ij (1,ic)
        j = ipc% pe(pe)% ij (2,ic)
        sendbuf (i0+1 : i0+nl) = xh(i,j,:)
        i0 = i0 + nl
      end do
    end do
#endif
    !---------------
    ! send / receive
    !---------------
FTRACE_BEGIN("apply_T:p_alltoall")
    call p_alltoall (sendbuf, recvbuf, sendcounts=sc, recvcounts=rc)
FTRACE_END  ("apply_T:p_alltoall")
    !-----------------------------
    ! get data from receive buffer
    !-----------------------------
FTRACE_BEGIN("apply_T:get_recvbuf")
    i0 = 0
    do pe= 0, dace% npe-1
      nl = sum(ipc% pe(pe)% nl (:))
!NEC$ ivdep
      do ip = 1, ipc% np_v
        if (ipc% pe(pe)% nl (ip) > 0) then
          lb = ipc% pe(pe)% lb (ip)
          ub = ipc% pe(pe)% ub (ip)
          lbl= ipc% pe(pe)% lbl(ip)
          ubl= ipc% pe(pe)% ubl(ip)
!NEC$ ivdep
          do ic = 1, ipc% nc_v
            ir = i0 + (ic-1)*nl
            xv (lb:ub,ip,ic) = recvbuf (ir+lbl:ir+ubl)
          end do
        end if
      end do
      i0 = i0 + ipc% nc_v * nl
    end do
FTRACE_END  ("apply_T:get_recvbuf")
  end subroutine apply_T
!------------------------------------------------------------------------------
  subroutine apply_WhLh (xd, xt)
    real(wp) ,intent(in)   :: xd(:,:,:) ! argument in 'diagonal' representation
    real(wp) ,intent(inout):: xt(:,:,:) ! result
    !----------------
    ! Local variables
    !----------------
    integer                      :: m, l, id, idx, basis1, basis2
    real(wp)                     :: x(size (xd,1),size (xd,2))
    real(wp)                     :: y(size (xd,1)*size (xd,2))
    type(t_cov_wavelet), pointer :: W
    integer                      :: j     ! zonal index
    integer                      :: ny    ! number of latitudes
    integer                      :: nc1   ! no.filter coefficients + 1

    if (size (xd,3) /= covm% hc2_part_size) then
       write (0,*) "apply_WhLh: size (xd,3),hc2_part_size =", &
                    size (xd,3), covm% hc2_part_size
       call finish ("apply_WhLh","inconsistent partitioning")
    end if
    !-----------------------
    ! Loop over partitioning
    !-----------------------
    ny = covm% ny
    do m = 1, covm% hc2_totalsize
       if (dace% pe /= covm% hc2_part(m)% pe) cycle
       l   =  covm% hc2_part(m)% l      ! Local index
       id  =  covm% hc2_part(m)% id     ! Variable id
       idx =  covm% hc2_part(m)% i      ! Variable index
       W   => covm% hc2(idx)
       if (W% id /= id) then
          call finish ("apply_WhLh","id mismatch")
       end if
       if (.not. W% valid) then
          call finish ("apply_WhLh","invalid covariance matrix")
       end if
       !-------------------------
       ! Apply L in wavelet space
       !-------------------------
       y = reshape (xd(:,:,l), shape (y))
       if (l_L_h) y = W% b * y
       x(:,:) = reshape (y, shape (x))
       if (l_W_h) then
         !------------------
         ! Wavelet synthesis
         !------------------
         basis1 = W% meta% wv_basis(1)    ! Transformation in x-y plane
         basis2 = W% meta% wv_basis(2)
         call wave_1d  (x, TR_SYN, basis1)
         call wave_1dt (x, TR_SYN, basis2)
         !-------------------------
         ! Filter in Fourier domain
         !-------------------------
         if (n_zonfilt > 0) then
           call fft      (x (:,   :n_zonfilt   ), -1)
           call fft      (x (:, ny-n_zonfilt+1:), -1)
           do j = 1, n_zonfilt
             nc1 = zonfilt (j) * 2
             x (nc1:,    j  ) = 0
             x (nc1:, ny-j+1) = 0
           end do
           call fft      (x (:,   :n_zonfilt   ),  1)
           call fft      (x (:, ny-n_zonfilt+1:),  1)
         endif
       endif
       xt(:,:,l) = x(:,:)
       nullify (W)
    end do

  end subroutine apply_WhLh
!------------------------------------------------------------------------------
  subroutine apply_WhLh_t (xd, xt)
    real(wp) ,intent(out) :: xd (:,:,:) ! result in 'diagonal' representation
    real(wp) ,intent(in)  :: xt (:,:,:) ! argument
    !----------------
    ! Local variables
    !----------------
    integer                      :: m, l, id, idx, basis1, basis2
    real(wp)                     :: x(size (xt,1),size (xt,2))
    real(wp)                     :: y(size (xt,1)*size (xt,2))
    type(t_cov_wavelet), pointer :: W
    integer                      :: j     ! zonal index
    integer                      :: ny    ! number of latitudes
    integer                      :: nc1   ! no.filter coefficients + 1

    if (size (xt,3) /= covm% hc2_part_size) then
       write (0,*) size (xt,3), covm% hc2_part_size
       call finish ("apply_WhLh_t","inconsistent partitioning")
    end if
    !-----------------------
    ! Loop over partitioning
    !-----------------------
    ny = covm% ny
    do m = 1, covm% hc2_totalsize
       if (dace% pe /= covm% hc2_part(m)% pe) cycle
       l   =  covm% hc2_part(m)% l      ! Local index
       id  =  covm% hc2_part(m)% id     ! Variable id
       idx =  covm% hc2_part(m)% i      ! Variable index
       W   => covm% hc2(idx)
       if (W% id /= id) then
          call finish ("apply_WhLh_t","id mismatch")
       end if
       if (.not. W% valid) then
          call finish ("apply_WhLh_t","invalid covariance matrix")
       end if
       x(:,:) = xt(:,:,l)
       if (l_W_h) then
         !-------------------------
         ! Filter in Fourier domain
         !-------------------------
         if (n_zonfilt > 0) then
           call fft      (x (:,   :n_zonfilt   ), -1)
           call fft      (x (:, ny-n_zonfilt+1:), -1)
           do j = 1, n_zonfilt
             nc1 = zonfilt (j) * 2
             x (nc1:,    j  ) = 0
             x (nc1:, ny-j+1) = 0
           end do
           call fft      (x (:,   :n_zonfilt   ),  1)
           call fft      (x (:, ny-n_zonfilt+1:),  1)
         endif
         !--------------------------
         ! Adjoint wavelet synthesis
         !--------------------------
         basis1 = W% meta% wv_basis(1)    ! Transformation in x-y plane
         basis2 = W% meta% wv_basis(2)
         call wave_1d  (x, TR_ADJ, basis1)
         call wave_1dt (x, TR_ADJ, basis2)
       endif
       !---------------------------
       ! Apply L^T in wavelet space
       !---------------------------
       y = reshape (x, shape (y))
       if (l_L_h) y = y * W% b
       xd(:,:,l) = reshape (y, shape (x))
       nullify (W)
    end do

  end subroutine apply_WhLh_t
!==============================================================================
  subroutine trace_r3 (comment, x)
  character(len=*) ,intent(in) :: comment
  real(wp)         ,intent(in) :: x (:,:,:)
    real(wp) :: mi,ma,va
    mi = p_min (minval(x))
    ma = p_max (maxval(x))
    va = p_sum (sum (x*x)) / p_sum (size(x))
    if (dace% lpio) write(6,'(a32,3f15.3)')  comment,mi,ma,va
  end subroutine trace_r3
!------------------------------------------------------------------------------
  subroutine trace_r2 (comment, x)
  character(len=*) ,intent(in) :: comment
  real(wp)         ,intent(in) :: x (:,:)
    real(wp) :: mi,ma,va
    mi = p_min (minval(x))
    ma = p_max (maxval(x))
    va = p_sum (sum (x*x)) / p_sum (size(x))
    if (dace% lpio) write(6,'(a32,3f15.3)')  comment,mi,ma,va
  end subroutine trace_r2
!------------------------------------------------------------------------------
  subroutine trace_r1 (comment, x)
  character(len=*) ,intent(in) :: comment
  real(wp)         ,intent(in) :: x (:)
    real(wp) :: mi,ma,va
    mi = p_min (minval(x))
    ma = p_max (maxval(x))
    va = p_sum (sum (x*x)) / p_sum (size(x))
    if (dace% lpio) write(6,'(a32,3f15.3)')  comment,mi,ma,va
  end subroutine trace_r1
!------------------------------------------------------------------------------
  subroutine trace_r  (comment, x)
  character(len=*) ,intent(in) :: comment
  real(wp)         ,intent(in) :: x
    real(wp) :: mi,ma,va
    mi = p_min (x)
    ma = p_max (x)
    va = p_sum (x*x)
    if (dace% lpio) write(6,'(a32,3f15.3)')  comment,mi,ma,va
  end subroutine trace_r
!------------------------------------------------------------------------------
  subroutine trace_v (comment, x, bx)
  character(len=*) ,intent(in)           :: comment
  type(t_vector)   ,intent(in)           :: x
  type(t_box)      ,intent(in) ,optional :: bx(:)

    integer, parameter :: I_H  = 1
    integer, parameter :: I_T  = 2
    integer, parameter :: I_RH = 3
    integer, parameter :: I_U  = 4
    integer, parameter :: I_V  = 5
    integer, parameter :: I_X  = 5
    real(wp)           :: mi, ma,va
    real(wp)           :: mio(I_X),mao(I_X),vao(I_X)
    integer            :: no (I_X)
    integer            :: ibl, ib, i, j

    mi = minval(x)
    ma = maxval(x)
    va = sum (x*x) / size(x)
    if (dace% lpio) write(6,'(a32,3f15.3)')  comment,mi,ma,va

    if (present(bx)) then
      mio =  999._wp
      mao = -999._wp
      vao =    0._wp
      no  =    0
      do ibl = 1, size(bx)
        ib = bx(ibl)% ib

        do j = 1, x% s(ib)% n

          select case (bx(ibl)% ip(j))
          case (OBS_H, OBS_HS)
            i = I_H
          case (OBS_TV)
            i = I_T
          case (OBS_RH)
            i = I_RH
          case (OBS_U)
            i = I_U
          case (OBS_V)
            i = I_V
          case default
            cycle
          end select

          no (i) = no (i) + 1
          mio(i) = min (mio(i), x% s(ib)% x(j))
          mao(i) = max (mao(i), x% s(ib)% x(j))
          vao(i) = vao(i) +     x% s(ib)% x(j) **2

        end do
      end do

      no  = p_sum (no)
      mio = p_min (mio)
      mao = p_max (mao)
      vao = p_sum (vao)
      where (no>0) vao = vao / no

      if (dace% lpio) then
        write (6, '(a32,3f15.3,i10)') 'H' ,mio(1), mao(1), vao(1), no(1)
        write (6, '(a32,3f15.3,i10)') 'T' ,mio(2), mao(2), vao(2), no(2)
        write (6, '(a32,3f15.3,i10)') 'RH',mio(3), mao(3), vao(3), no(3)
        write (6, '(a32,3f15.3,i10)') 'U' ,mio(4), mao(4), vao(4), no(4)
        write (6, '(a32,3f15.3,i10)') 'V' ,mio(5), mao(5), vao(5), no(5)
      endif

    endif
  end subroutine trace_v
!==============================================================================
end module mo_bg_err_2d
