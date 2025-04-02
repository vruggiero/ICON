!
!+ Covariance matrices, 2d(hor,wavelet)+1d(vert), medium level routines
!
#if __NEC_VERSION__ / 100 == 305    /* Work around issues with nfort 3.5.x */
!NEC$ options "-fno-reorder-logical-expression"
#endif

MODULE mo_bg_err_op
!
! Description:
!   Representation of covariance matrices by a sequence of operators :
!     Vertical covariances fitted by the NMC method are represented by
!     matrices.
!     Horizontal covariances are represented by the explicit OI model.
!
!   This module holds medium level routines.
!   Low level routines and derived type definitions are placed
!   in module mo_t_bg_err_op .
!   High level entries and horizontal operator representations are placed
!   in module mo_bg_err_2d .
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
! V1_5         2009/05/25 Harald Anlauf
!  modify_vcov: use global lower bound on RH error
!  new namelist parameter p_const_v
! V1_7         2009/08/24 Harald Anlauf
!  print_bg_err_op_nml: print i_vloc, p_vloc
! V1_8         2009/12/09 Harald Anlauf
!  modify_vcov: implement recognition of file version for stdv_ratio
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Andreas Rhodin
!  option to adjust vertical covariance grid to model top
! V1_22        2013-02-13 Andreas Rhodin
!  1dvar mode implementation
! V1_28        2014/02/26 Andreas Rhodin
!  preparations for VarEnKF
! V1_31        2014-08-21 Andreas Rhodin
!  changes for 1dvar-mode and rttov specific vertical interpolation
! V1_35        2014-11-07 Andreas Rhodin
!  bug fix for 1dvar-mode & IR PC emissivity
! V1_37        2014-12-23 Harald Anlauf
!  implement stdv_ratio format version 2
! V1_42        2015-06-08 Harald Anlauf
!  modify_vcov: implement stdv_ratio format version 3 (interpolate in ln(p))
! V1_48        2016-10-06 Andreas Rhodin
!  handle z-based vertical coordinate in all post-multiplication branches
! V1_50        2017-01-09 Andreas Rhodin
!  set_bg_err_op: check if NMC/OI covariance model is aready set up
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Andreas Rhodin  DWD  2006-2008  original code
! Harald Anlauf   DWD  2007       horizontal wavelet approach
! Harald Anlauf   DWD  2008       optimizations for SX8
!=============================================================================
  !=============
  ! modules used
  !=============
  use mo_kind,       only: wp, i8          ! precision kind parameter
  use mo_exception,  only: finish,        &! abort on error condition
                           message         ! print warning message
  use mo_dec_matrix, only: t_vector,      &! partitioned vector data type
                           construct,     &! construct a vector
                           destruct,      &! destruct  a vector
                           delete_storage,&!
                           assignment(=), &! vector assignment
                           mp              ! kind parameter
  use mo_allocate,   only: enter_function,&! to be called at start of routines
                           leave_function  ! to be called at end of routines
  use mo_namelist,   only: position_nml,  &! routine to position nml group
                           nnml,          &! namelist fortran unit number
                           POSITIONED,    &! position_nml: OK  return flag
                           MISSING         !               namelist missing
  use mo_mpi_dace,   only: dace,          &! MPI group info
                           p_bcast,       &! broadcast routine
                           p_barrier       ! MPI barrier routine
  use mo_time,       only: t_time,        &! date and time data type
                           imm             ! get month
  use mo_dwd,        only: nfcoun,        &! climatol. 6h forecast error cf.
                           read_oi_coef    ! read OI bg-error model coeffs.
  use mo_fg_cov,     only: init_fg_cov,   &! initialize OI covariance model
                        read_nml_pscmodel,&! read OI covariance model namelist
                           t_rowcol,      &! OI bg-error meta data type
                           construct,     &! construct type t_rowcol
                           zero,          &! zero type t_rowcol for reuse
                           destruct,      &! destruct  type t_rowcol
                           init_spot,     &! precalculate meta data
                           get_corr,      &! get correlations
                           IFC_CFC,       &! climatological fc error flag
!                          IFC_AFC,       &! actual         fc error flag
                           L_h, L_q,      &! const. hor.lengthscale (namelist)
                           L_max,         &! max.l-scale (unit sphere)
                           x_cut,         &! cutoff radius
                           rnu,           &! e2_x / e2_v
                           rmu,           &! (constant) height wind correlation
                           lsvv,          &! largescale contrib. to uv-uv cor.
                           corr,          &! horizontal correlation function
                           dcorr,         &! hor.corr.function, derivatives
                           corrh,         &! hor.corr.function, rel.humidity
                           e_h,           &! (constant) height error
                           e_v,           &! (constant) wind   error
                           e_q,           &! (constant) relative humidity error
                           emin_q,        &! global lower bound on RH error
                           init_spot_obs   ! set corr. data type from obs.
  use mo_1dmra,      only: wave_1d,       &! wavelet transform
                           gp_mr_mat,     &! inverse matrix wavelet transform
                           TR_SYN,        &! Synthesis flag
                           TR_ADJ,        &! Adjoint flag
                           basis,         &! wavelet basis to use
                           WV_CUB_LIFT_INT ! wavelet basis used
  use mo_matrix,     only: sqrt_rs,       &! Squareroot of realsymmetric matrix
                           cholesky,      &! Cholesky factorization
                           scale,         &! scale a matrix
                           diag            ! diagonal of a matrix
  use mo_t_obs,      only: t_spot,        &! report meta data type
                           t_coord,       &! coordinates meta data type
                           t_hic1,        &! horizontal interpolation data type
                           t_vic,         &! vertical   interpolation data type
                           t_lev,         &! component of t_spot
                           set_xuv,       &! set up unit vectors on the sphere
                           col2lev,       &! set component 'levs' of 't_obs'
                           OBS_TV,        &! virtual temperature  id
                           OBS_RH,        &! relative humidity    id
                           OBS_H,         &! geopotential height  id
                           OBS_HS,        &! surface geopotential id
                           OBS_U,         &! u-wind               id
                           OBS_V,         &! v-wind               id
                           OBS_CHI,       &! velocity potential   id
                           OBS_DUM,       &! dummy type
                           ITY_ICOL,      &! interpolation type: column
                           ITY_MCOLS       !               model columns
  use mo_obs_set,    only: t_obs_set       ! observation data type
  use mo_t_bg_err_op,only: covm,           &! NMC fitted covariances
                           t_cov,          &! NMC covariances data type
                           construct,      &! construct type t_cov (covm)
                           cov_from_nmc,   &! read NMC fitted covariances
                           set_vertinpc,   &! set vertical interpolation coefs.
                           set_horinpc,    &! set horizont.interpolation coefs.
                           compress_cov,   &! store covm only once
                           uncompress_cov, &! store covm on each PE
                           read_horizontal_covariances  ! Hor. wavelet covs.
  use mo_physics,    only: R,              &! gas constant of dry air
                           gacc,           &! gravity acceleration
                           pi,             &! 3.1415..
                           d2r              ! factor: degree -> radians
  use mo_t_col,      only: t_cols,         &! data type to hold columns
                           alloc_cols,     &! allocate components of t_cols
                           assignment(=),  &! type(t_cols) = real(wp)
                           COL_T, COL_UV,  &! identifier for temp., hor.wind,
                           COL_RH           !   rel.humidity
  use mo_mem_usage,  only: trace_mem_usage,&! report memory usage
                           print_mem_usage  ! print memory usage
  use mo_grads,      only: t_ctl,          &! data type: .ctl file information
                           init_ctl,       &! preset .ctl file
                           write_var,      &! write variable to grads file
                           write_ctl,      &! write .ctl file
                           destruct         ! destruct t_ctl data type
  use mo_run_params, only: data,           &! data path
                           aux,            &! diagnostics output path
                           path_file        ! concatenate path/filename
  use mo_dace_string,only: split            ! Split text string into array
  use mo_atm_grid,   only: VCT_P_ISO,      &!       isobaric p coordinate
                           VCT_P_HYB,      &! GME/HRM hybrid p coordinate
                           VCT_Z_HYB,      &!         hybrid z coordinate
                           VCT_Z_GEN        !    generalised z coordinate
  use mo_fortran_units, only: get_unit_number, &!
                              return_unit_number!
  implicit none

  !================
  ! public entities
  !================
  private
  !----------
  ! operators
  !----------
  public :: apply_B_ii_1d     ! multiply  a vector in intp. space by B
  public :: apply_B_mi_1d     ! post multiplication
  !------------------------------------
  ! utility routines and data type defs
  !------------------------------------
  public :: set_bg_err_op     ! set the covariance model
  public :: t_mode            ! data type definition
  public :: vertloc           ! set vertical localisation matrix
  public :: average_merid     ! meridionally average vertical covariance matrix
  !------------------------
  ! namelist steering flags
  !------------------------
  public :: repr_2dh          ! horizontal operator representation flag
  public :: n_av_pole         ! average anal.incr. at poles (no.neighb.points)
  public :: zonfilt           ! number of Fourier coefs.in zonal filter
  public :: n_zonfilt         ! number of latitudes near poles to filter
  public :: l1dvar            ! run in 1dvar mode (no horizontal correlations)
  public :: horz_inp          ! horizontal interpolation: 2:bi-linear 4:higher-order
  public :: horz_grid_inp     ! horiz. grid interp.: 1:nn 2:bi-lin.
                              !     3: bilin. new. 5: bilin. vertex
  public :: horz_nmc_inp      ! horizontal interpolation: 2:bi-linear 4:higher-order
  !-----------
  ! test flags
  !-----------
  public :: test_bg           ! flag: compare with explicit representation
  public :: test_bg_b         ! bound for failure of test_bg
  public :: test_adj          ! flag: compare linear vs adjoint
  public :: ltrace            ! trace minval, maxval, variance
  !-------------------
  ! for debugging only
  !-------------------
  public :: vert_inp          ! test: 1=next-neighbour, 2=interpolate
  public :: l_L_h             ! apply L    matrix in the horizontal
  public :: l_W_h             ! apply W    matrix in the horizontal
  public :: unit_v            ! apply unit matrix in the vertical
  !-----------------------------------
  ! first guess error estimation flags
  !-----------------------------------
  public :: fger_iter         ! no. iterations for bg-error estimate
  public :: fger_est          ! estimate first guess error
  public :: fger_Li1          ! test apply_L_i with e_fi=1, abort
  public :: fger_seed         ! generate seed (for parallel run)
  public :: fger_write        ! write fg error to file
  public :: fger_ana          ! estimate for analysis step


  !=================
  ! module variables
  !=================
  logical           :: l1dvar    = .false.! run in '1dvar' mode

  !=============================
  ! namelist /BG_ERROR_OPERATOR/
  !=============================
  !------------------------------
  ! number of gridpoints, bounds
  !------------------------------
  integer           :: nz        =  64    ! number of vertical levels
  integer           :: ny        =   0    ! number of latitudinal  gridpoints
  integer           :: nx        =   0    ! number of longitudinal gridpoints
  real(wp)          :: pbot      = 1070.  ! bottom (surface) pressure (hPa)
! real(wp)          :: lats (2)  =(/-90._wp,90._wp/) ! latitude range
  real(wp)          :: rnd_time  = 0._wp  ! random B correlation time scale (h)
  !-----------------------------------
  ! specification of covariance model
  !-----------------------------------
  integer           :: valid     = 0      ! 1:content is valid 2:use cov.model
  integer           :: repr_2dh  = 0      ! hor. operator repr: 1=test, 2=use
  character(len=128):: file      = ''     ! input file name, NMC, vert.cov.
  character(len=128):: file_2dh  = ''     ! input file name, NMC, hor. cov.
  integer           :: lclim     = 2      ! 1:clim.err.  valid 2:use clim.err.
  real(wp)          :: scale_rh  = 1._wp  ! scaling factor rel.hum. fg.error
  integer           :: horz_inp  = 2      ! 2:bi-linear 4:higher-order
  integer           :: horz_grid_inp = -99 ! horiz. grid interp.: 1:nn 2|3|5 :bilin.
  integer           :: horz_nmc_inp  = -99 ! 2:bi-linear 4:higher-order
  !--------------------------------
  ! covariance matrix manipulations (partly not implemented)
  !--------------------------------
  logical           :: modify_cor= .true. ! modify correlations, else variances
  real(wp)          :: lat_homog = 99._wp ! latitudinal homogeneous covariances
  real(wp)          :: smooth_h  = 0._wp  ! smooth horizontally
  real(wp)          :: smooth_v  = 0._wp  ! smooth vertically
  character(len=128):: stdv_ratio= ''     ! file with stdev.ratios 3h/NMC time
  integer           :: i_vloc    = 0      ! vertical localisation flag
  real(wp)          :: p_vloc(2) =(/0,9999/)! bounds of localisation interval
  integer           :: i_av_vcov = 0      ! cov.matrix meridion. average flag
  real(wp)          :: p_const_q = 0._wp  ! constant RH covariances below (hPa)
  real(wp)          :: p_const_v = 0._wp  ! constant wind covar.    below (hPa)
  real(wp)          :: p_adj     = 0._wp  ! adapt levs.to model top above (hPa)
  !---------
  ! details
  !---------
  character         :: transform = 'w'    ! 'w'avelet, 'n'one
  character         :: sqr       = 's'    ! 's'ymmetric, 'c'holesky
  integer           :: nwv       = 4      ! no.points in vertical interpolation
  integer           :: base      = WV_CUB_LIFT_INT ! wavelet basis
  logical           :: h2psi     = .true. ! derive streamf.    cov. from height
  logical           :: h2t       =.false. ! derive temperature.cov. from height
  logical           :: sqrtcor   =.false. ! derive sqrt from corr.,not cov.
  logical           :: sparse    =.false. ! sparse repres. of vertical matrices
  integer           :: n_av_pole = -1     ! average anal.incr.at poles (no.pts)
  integer           :: n_zonfilt = 3      ! no of latitudes at poles to filter
  integer           :: zonfilt(16)=     & ! # Fourier coefs., zonal filter
                       (/2,2,3,0,0,0,0,0,0,0,0,0,0,0,0,0/) ! <0: auto
#if defined (__SX__)
  logical           :: csr_symm  = .true. ! exploit symmetry of sqrt(B)
  character(len=8)  :: csr_conv  = "JAD"  ! convert from CSR to representation
#else
  logical           :: csr_symm  =.false. ! do not exploit symmetry of sqrt(B)
  character(len=8)  :: csr_conv  = ""     ! convert from CSR to representation
#endif
  !-----------------------------
  ! First Guess error estimation
  !-----------------------------
  integer           :: fger_iter = 25     ! no.iterations for bg-error estimate
  logical           :: fger_est  =.false. ! estimate first guess error
  logical           :: fger_Li1  =.false. ! test apply_L_i with e_fi=1, abort
  logical           :: fger_seed = .true. ! generate seed (for parallel run)
  logical           :: fger_write=.false. ! write fg error to file
  integer           :: fger_ana  =  1     ! estimate for analysis step
                                          ! 0: none
  !-------                                ! 1: take from monitoring scan
  ! tests                                 ! 2: recalculate in analysis scan
  !-------
  logical           :: test_bg   =.true. ! compare operator impl. with explicit
  real(wp)          :: test_bg_b = 1.e-5 ! bound for failure of test_bg
  integer           :: test_adj  = 2     ! compare linear vs adjoint
  logical           :: ltrace    =.false.! trace minval,maxval,variance
  !----------
  ! debugging
  !----------
  integer           :: vert_inp  = 2     ! 1:next-neighbour, 2:interpolate
  logical           :: l_L_h     = .true.! apply L    matrix in the horizontal
  logical           :: l_W_h     = .true.! apply W    matrix in the horizontal
  logical           :: unit_v    =.false.! apply unit matrix in the vertical

  namelist /BG_ERROR_OPERATOR/ nx, ny, nz, pbot, file, transform, sqr, valid, &
                               base, nwv, h2psi, h2t, sqrtcor, sparse, lclim, &
                               scale_rh, test_bg, repr_2dh, file_2dh,         &
                               modify_cor, lat_homog, smooth_h, smooth_v,     &
                               test_adj, horz_inp, l_L_h, l_W_h, vert_inp,    &
                               test_bg_b, unit_v, ltrace,                     &
                               fger_iter, fger_est, fger_seed, fger_write,    &
                               fger_Li1, fger_ana, n_av_pole, zonfilt,        &
                               stdv_ratio, i_vloc, p_vloc, i_av_vcov,         &
                               p_const_q, p_const_v, p_adj, csr_symm,         &
                               csr_conv, rnd_time, horz_grid_inp, horz_nmc_inp
  !======================
  ! data type definitions
  !======================
  type t_mode
    real(wp) ,pointer :: h  (:,:) ! geopotential height
    real(wp) ,pointer :: u  (:,:) ! u-wind component
    real(wp) ,pointer :: v  (:,:) ! v-wind component
    real(wp) ,pointer :: ux (:,:) ! u-wind component
    real(wp) ,pointer :: vx (:,:) ! v-wind component
    real(wp) ,pointer :: q  (:,:) ! humidity
    logical  ,pointer :: lh   (:) ! used flags
    logical  ,pointer :: luv  (:) !
    logical  ,pointer :: lq   (:) !
    integer(i8)       :: size
  end type t_mode

  !===========
  ! interfaces
  !===========
  interface construct
    module procedure construct_modes
    module procedure construct_mode
  end interface construct

  interface destruct
    module procedure destruct_modes
    module procedure destruct_mode
  end interface destruct

contains
!==============================================================================
  subroutine construct_modes (modes, nlev, ncol, mask)
  type (t_mode) ,intent(out)          :: modes (:)
  integer       ,intent(in)           :: nlev
  integer       ,intent(in)           :: ncol  (:)
  logical       ,intent(in) ,optional :: mask  (:)
    integer :: i
    logical :: m
    do i = 1, size (modes)
      m = .true.; if(present(mask)) m = mask (i)
      call construct_mode (modes(i), nlev, ncol(i), m)
    end do
  end subroutine construct_modes
!------------------------------------------------------------------------------
  subroutine construct_mode (modes, nlev, ncol, mask)
  type (t_mode) ,intent(out) :: modes
  integer       ,intent(in)  :: nlev
  integer       ,intent(in)  :: ncol
  logical       ,intent(in)  :: mask
    if (mask) then
      allocate (modes% h  (nlev, ncol)); modes%  h = 0._wp
      allocate (modes% u  (nlev, ncol)); modes%  u = 0._wp
      allocate (modes% v  (nlev, ncol)); modes%  v = 0._wp
      allocate (modes% ux (nlev, ncol)); modes%  ux = 0._wp
      allocate (modes% vx (nlev, ncol)); modes%  vx = 0._wp
      allocate (modes% q  (nlev, ncol)); modes%  q = 0._wp
      allocate (modes% lh       (ncol)); modes% lh  = .false.
      allocate (modes% luv      (ncol)); modes% luv = .false.
      allocate (modes% lq       (ncol)); modes% lq  = .false.
      modes% size = size(modes% h  ) + &
                    size(modes% u  ) + &
                    size(modes% v  ) + &
                    size(modes% ux ) + &
                    size(modes% vx ) + &
                    size(modes% q  ) + &
                    size(modes% lh ) + &
                    size(modes% luv) + &
                    size(modes% lq)
      modes% size = modes% size * size(transfer(0._wp,(/' '/)))
    else
      nullify (modes% h ,modes% u  ,modes% v ,modes% ux, modes% vx, modes% q, &
               modes% lh,modes% luv,modes% lq)
      modes% size = 0
    endif
  end subroutine construct_mode
!-----------------------------------------------------------------------------
  subroutine destruct_modes (modes)
  type (t_mode) ,intent(inout) :: modes (:)
    integer :: i
    do i = 1, size (modes)
      call destruct_mode(modes(i))
    end do
  end subroutine destruct_modes
!-----------------------------------------------------------------------------
  subroutine destruct_mode  (modes)
  type (t_mode) ,intent(inout) :: modes
    if (associated (modes% h  )) deallocate (modes% h  )
    if (associated (modes% u  )) deallocate (modes% u  )
    if (associated (modes% v  )) deallocate (modes% v  )
    if (associated (modes% ux )) deallocate (modes% ux )
    if (associated (modes% vx )) deallocate (modes% vx )
    if (associated (modes% q  )) deallocate (modes% q  )
    if (associated (modes% lh )) deallocate (modes% lh )
    if (associated (modes% luv)) deallocate (modes% luv)
    if (associated (modes% lq )) deallocate (modes% lq )
    modes% size = 0
  end subroutine destruct_mode
!-----------------------------------------------------------------------------
  subroutine scatter_modes (modes, obs)
  type (t_mode)    ,intent(inout) :: modes(:) ! vertical modes  rhs
  type (t_obs_set) ,intent(in)    :: obs      ! observation meta data
  !--------------------------------------
  ! scatter modes over processor elements
  !--------------------------------------
    integer :: ib ! box  index
    !----------------
    ! loop over boxes
    !----------------
    do ib = 1,size(modes)
      call p_bcast (modes(ib)% h   ,obs% o(ib)% pe)
      call p_bcast (modes(ib)% u   ,obs% o(ib)% pe)
      call p_bcast (modes(ib)% v   ,obs% o(ib)% pe)
      call p_bcast (modes(ib)% ux  ,obs% o(ib)% pe)
      call p_bcast (modes(ib)% vx  ,obs% o(ib)% pe)
      call p_bcast (modes(ib)% q   ,obs% o(ib)% pe)
      call p_bcast (modes(ib)% lh  ,obs% o(ib)% pe)
      call p_bcast (modes(ib)% luv ,obs% o(ib)% pe)
      call p_bcast (modes(ib)% lq  ,obs% o(ib)% pe)
    end do
  end subroutine scatter_modes
!=============================================================================
! Namelist processing
!====================

  subroutine read_bg_err_op_nml
  !==================================
  ! read namelist /BG_ERROR_OPERATOR/
  !==================================
  integer :: i, ierr, j, fmax, nc1

    !==========================
    ! preset namelist variables
    !==========================
    nz        = 64       ! number of vertical levels
    ny        =   0      ! number of latitudinal  gridpoints
    nx        =   0      ! number of longitudinal gridpoints
    pbot      = 1070.    ! bottom (surface) pressure (hPa)
    rnd_time  = 0._wp    ! correlation time scale for random B (h)
    file      = ''       ! input file name, vertical   NMC derived covariances
    file_2dh  = ''       ! input file name, horizontal NMC derived covariances
    transform = 'w'      ! 'w'avelet, 'n'one
    sqr       = 's'      ! 's'ymmetric, 'c'holesky
    valid     = 0        ! 1:content is valid;  2:use this covariance model
    repr_2dh  = 0        ! horizontal operator representation: 1=test, 2=use
    lclim     = 2        ! 1:clim.err.  valid 2:use clim.err.
    nwv       = 4        ! no.points in vertical interpolation
    base      = WV_CUB_LIFT_INT ! wavelet basis
    h2psi     = .true.   ! derive streamf.    cov. from height
    h2t       =.false.   ! derive temperature.cov. from height
    sqrtcor   =.false.   ! derive sqrt from corr.,not cov.
    sparse    =.false.   ! sparse repres. of vertical matrices
    n_av_pole = -1       ! average anal.incr. at poles (no.pts)
    zonfilt   =         &!no.Fourier coefs.in zonal filter
                (/2,2,3,0,0,0,0,0,0,0,0,0,0,0,0,0/)
#if defined (__SX__)
    csr_symm  = .true.   ! exploit symmetry of sqrt(B)
    csr_conv  = "JAD"    ! convert from CSR to representation
#else
    csr_symm  =.false.   ! do not exploit symmetry of sqrt(B)
    csr_conv  = ""       ! convert from CSR to representation
#endif
    scale_rh  = 1._wp    ! scaling factor for rel.hum. fg.error
    test_bg   =.true.    ! compare operator impl. with explicit
    test_bg_b = 1.e-5    ! bound for failure of test_bg
    test_adj  = 2        ! compare linear vs adjoint
    horz_inp  = 2        ! 2:bi-linear 4:higher-order
    horz_grid_inp = -99  ! horiz. grid interp.: (1):nn (2|3|5):bi-lin.
    horz_nmc_inp  = -99  ! 2:bi-linear 4:higher-order
    vert_inp  = 2        ! 1:next-neighbour, 2:interpolate
    l_L_h     = .true.   ! apply L    matrix in the horizontal
    l_W_h     = .true.   ! apply W    matrix in the horizontal
    unit_v    =.false.   ! apply unit matrix in the vertical
    ltrace    =.false.   ! trace minval,maxval,variance
!   lats (2)  =(/-90._wp,90._wp/) ! latitude range
    modify_cor= .true.   ! modify correlations, else variances
    lat_homog = 99._wp   ! latitudinal homogeneous covariances
    smooth_h  = 0._wp    ! smooth horizontally
    smooth_v  = 0._wp    ! smooth vertically
    stdv_ratio= ''       ! file with stdev. ratios 3h/NMC time
    i_vloc    = 0        ! vertical localisation flag
    p_vloc    =(/0,9999/)! bounds of localisation interval
    i_av_vcov = 0        ! covariance matrix meridional average flag
    p_const_q = 0._wp    ! constant RH   covariances below (hPa)
    p_const_v = 0._wp    ! constant wind covariances below (hPa)
    p_adj     = 0._wp    ! adapt levels to model top above (hPa)
    !-----------------------------
    ! First Guess error estimation
    !-----------------------------
    fger_iter = 25     ! no. iterations for bg-error estimate
    fger_est  =.false. ! estimate first guess error
    fger_Li1  =.false. ! test apply_L_i with e_fi=1, abort
    fger_seed = .true. ! generate seed (for parallel run)
    fger_write=.false. ! write fg error to file
    fger_ana  =  1     ! estimate for analysis step

    !--------------
    ! read namelist
    !--------------
    if (dace% lpio) then
      write (6,'()')
      write (6,'(a)') repeat ('-',79)
      write (6,'(a)') '  reading namelist /BG_ERROR_OPERATOR/'
      write (6,'()')
      call position_nml ('BG_ERROR_OPERATOR', status=ierr)
      select case (ierr)
      case (POSITIONED)
        write (6,'(a)') '  namelist /BG_ERROR_OPERATOR/ found'
#if defined(__ibm__)
        read (nnml ,nml=BG_ERROR_OPERATOR, iostat=ierr)
        if (ierr/=0) call finish ('read_bg_err_op_nml',                  &
                                  'ERROR in namelist /BG_ERROR_OPERATOR/')
#else
        read (nnml ,nml=BG_ERROR_OPERATOR)
#endif
      case (MISSING)
        write (6,'(a)') '  namelist /BG_ERROR_OPERATOR/ not present, using defaults'
      case default
        call finish ('read_bg_err_op_nml',           &
                     'ERROR in reading namelist file')
      end select
    endif
    !-----------------------------
    ! broadcast namelist variables
    !-----------------------------
    call p_bcast (nx         ,dace% pio)
    call p_bcast (ny         ,dace% pio)
    call p_bcast (nz         ,dace% pio)
    call p_bcast (pbot       ,dace% pio)
    call p_bcast (rnd_time   ,dace% pio)
    call p_bcast (file       ,dace% pio)
    call p_bcast (file_2dh   ,dace% pio)
    call p_bcast (transform  ,dace% pio)
    call p_bcast (sqr        ,dace% pio)
    call p_bcast (valid      ,dace% pio)
    call p_bcast (repr_2dh   ,dace% pio)
    call p_bcast (lclim      ,dace% pio)
    call p_bcast (base       ,dace% pio)
    call p_bcast (nwv        ,dace% pio)
    call p_bcast (h2psi      ,dace% pio)
    call p_bcast (h2t        ,dace% pio)
    call p_bcast (sqrtcor    ,dace% pio)
    call p_bcast (sparse     ,dace% pio)
    call p_bcast (n_av_pole  ,dace% pio)
    call p_bcast (zonfilt    ,dace% pio)
    call p_bcast (csr_symm   ,dace% pio)
    call p_bcast (csr_conv   ,dace% pio)
    call p_bcast (fger_iter  ,dace% pio)
    call p_bcast (fger_est   ,dace% pio)
    call p_bcast (fger_Li1   ,dace% pio)
    call p_bcast (fger_seed  ,dace% pio)
    call p_bcast (fger_write ,dace% pio)
    call p_bcast (fger_ana   ,dace% pio)
    call p_bcast (scale_rh   ,dace% pio)
    call p_bcast (test_bg    ,dace% pio)
    call p_bcast (test_bg_b  ,dace% pio)
    call p_bcast (test_adj   ,dace% pio)
    call p_bcast (horz_inp   ,dace% pio)
    call p_bcast (horz_grid_inp, dace% pio)
    call p_bcast (horz_nmc_inp,  dace% pio)
    call p_bcast (vert_inp   ,dace% pio)
    call p_bcast (l_L_h      ,dace% pio)
    call p_bcast (l_W_h      ,dace% pio)
    call p_bcast (unit_v     ,dace% pio)
    call p_bcast (ltrace     ,dace% pio)
!   call p_bcast (lats       ,dace% pio)
    call p_bcast (modify_cor ,dace% pio)
    call p_bcast (lat_homog  ,dace% pio)
    call p_bcast (smooth_h   ,dace% pio)
    call p_bcast (smooth_v   ,dace% pio)
    call p_bcast (stdv_ratio ,dace% pio)
    call p_bcast (i_vloc     ,dace% pio)
    call p_bcast (p_vloc     ,dace% pio)
    call p_bcast (i_av_vcov  ,dace% pio)
    call p_bcast (p_const_q  ,dace% pio)
    call p_bcast (p_const_v  ,dace% pio)
    call p_bcast (p_adj      ,dace% pio)
    !--------------
    ! set path/file
    !--------------
    if (file /= 'OI') &
      file     = path_file (data, file)
    file_2dh   = path_file (data, file_2dh)
    stdv_ratio = path_file (data, stdv_ratio)
    !-------------------------------------------
    ! set points used for averaging at the poles
    !-------------------------------------------
    if (n_av_pole < 0) then
      if (valid >= 5) then
        n_av_pole = 2
      else
        n_av_pole = 0
      endif
    endif
    !-----------------------------
    ! zonal filtering at the poles
    !-----------------------------
    if(zonfilt(1) < 0) then
      !---------------------------------------------------------------------
      ! Number of coefficients not explicitly given for each latitude.
      ! Estimate of maximum frequency to keep at latitude j (for j=1..ny/2):
      ! fmax(j) = ny - 4*(j-ny/2)^2/ny
      ! fmax'(0) = 4 (indep. of ny), fmax(ny/2) ~ nx/2;
      ! and empirically fmax(1) ~ 2 (last lines closest to poles).
      !---------------------------------------------------------------------
      zonfilt = 0
      do j=1,size(zonfilt)
        fmax = 4*j - 2
        nc1  = 2 * (1+fmax)
!       if (nc1 > nx) exit
        zonfilt(j)= nc1/2
      end do
    endif
    !-----------------------------------------------
    ! count no. latitudes subject to zonal filtering
    !-----------------------------------------------
    n_zonfilt = 0
    do i = 1, size(zonfilt)
      if (zonfilt(i) <= 0) exit
      n_zonfilt = n_zonfilt + 1
    end do
    !---------
    ! printout
    !---------
    call print_bg_err_op_nml
  end subroutine read_bg_err_op_nml
!------------------------------------------------------------------------------
  subroutine print_bg_err_op_nml
    if (dace% lpio) then
      write(6,'()')
      write(6,'(a)')      '   number of gridpoints, bounds :'
      write(6,'()')
      write(6,'(a,i5)')   '    nx        = ',nx
      write(6,'(a,i5)')   '    ny        = ',ny
      write(6,'(a,i5)')   '    nz        = ',nz
      write(6,'(a,f6.0)') '    pbot      = ',pbot
      write(6,'(a,f6.0)') '    rnd_time  = ',rnd_time
      write(6,'()')
      write(6,'(a)')      '   specification of covariance model :'
      write(6,'()')
      write(6,'(a,i5)')   '    valid         = ',valid
      write(6,'(a,i5)')   '    repr_2dh      = ',repr_2dh
      write(6,'(a,a)')    '    file          = ',trim(file)
      write(6,'(a,a)')    '    file_2dh      = ',trim(file_2dh)
      write(6,'(a,i5)')   '    lclim         = ',lclim
      write(6,'(a,f6.2)') '    scale_rh      = ',scale_rh
      write(6,'(a,i5)')   '    horz_inp      = ',horz_inp
      write(6,'(a,i5)')   '    horz_grid_inp = ',horz_grid_inp
      write(6,'(a,i5)')   '    horz_nmc_inp  = ',horz_nmc_inp
      write(6,'()')
      write(6,'(a)')      '   details :'
      write(6,'()')
      write(6,'(a,a)')    '    transform = ',transform
      write(6,'(a,a)')    '    sqr       = ',sqr
      write(6,'(a,i5)')   '    base      = ',base
      write(6,'(a,i5)')   '    nwv       = ',nwv
      write(6,'(a,l1)')   '    h2psi     = ',h2psi
      write(6,'(a,l1)')   '    h2t       = ',h2t
      write(6,'(a,l1)')   '    sqrtcor   = ',sqrtcor
      write(6,'(a,l1)')   '    sparse    = ',sparse
      write(6,'(a,l1)')   '    fger_est  = ',fger_est
      write(6,'(a,i5)')   '    fger_iter = ',fger_iter
      write(6,'(a,l1)')   '    fger_Li1  = ',fger_Li1
      write(6,'(a,l1)')   '    fger_seed = ',fger_seed
      write(6,'(a,l1)')   '    fger_write= ',fger_write
      write(6,'(a,i5)')   '    fger_ana  = ',fger_ana
      write(6,'(a,i5)')   '    n_zonfilt = ',n_zonfilt
      write(6,'(a,16i5)') '    zonfilt   = ',zonfilt(1:n_zonfilt)
      write(6,'(a,l1)')   '    csr_symm  = ',csr_symm
      write(6,'(a,a)')    '    csr_conv  = ',csr_conv
      write(6,'(a,a)')    '    stdv_ratio= ',trim (stdv_ratio)
      write(6,'(a,f6.0)') '    p_const_q = ',p_const_q
      write(6,'(a,f6.0)') '    p_const_v = ',p_const_v
      write(6,'(a,f6.0)') '    p_adj     = ',p_adj
      write(6,'(a,i9)')   '    i_vloc    = ',i_vloc
      if (i_vloc > 0) &
      write(6,'(a,2f9.3)')'    p_vloc    = ',p_vloc
      write(6,'()')
      write(6,'(a)')      '   tests :'
      write(6,'()')
      write(6,'(a,l1)')   '    test_bg   = ',test_bg
      write(6,'(a,g6.1)') '    test_bg_b = ',test_bg_b
      write(6,'(a,l1)')   '    l_L_h     = ',l_L_h
      write(6,'(a,l1)')   '    l_W_h     = ',l_W_h
      write(6,'(a,l1)')   '    unit_v    = ',unit_v
      write(6,'(a,l1)')   '    ltrace    = ',ltrace
    endif
    !------------------
    ! consistency check
    !------------------
    valid = max (valid, 0)
    lclim = max (lclim, 0)
    lclim = min (lclim, valid)
    if (.not.h2psi .and. valid < 5) &
      call finish('print_bg_err_op_nml','.not.h2psi .and. valid < 5')
    select case (horz_grid_inp)
    case (1,2,3,4,5)
      select case (horz_nmc_inp)
      case (-99)
        call finish('print_bg_err_op_nml', &
                    'horz_grid_inp and horz_nmc_inp must be both set')
      end select
    case (-99)
    case default
      call finish('print_bg_err_op_nml', &
                  'horz_grid_inp must be one of 1,2,3,4,5')
    end select
    select case (horz_nmc_inp)
    case (2,4)
      select case (horz_grid_inp)
      case (-99)
        call finish('print_bg_err_op_nml', &
                    'horz_grid_inp and horz_nmc_inp must be both set')
      end select
      if (horz_nmc_inp /= horz_inp) then
        if (dace% lpio) write(6,*) &
          "Info: horz_nmc_inp /= horz_inp => horz_nmc_inp used only in postmult"
      end if
    case (-99)
    case default
      call finish('print_bg_err_op_nml','horz_nmc_inp must be one of 2,4')
    end select
  end subroutine print_bg_err_op_nml
!=============================================================================
! Set up vertical covariance matrices
!====================================
  subroutine set_bg_err_op (time, p_top)
  type (t_time) ,intent(in) :: time  ! time and date of analysis
  real(wp)      ,intent(in) :: p_top ! model top pressure (Pa)
  !=====================================
  !   read NMC derived covariances
  !   or derive from OI covariance model
  !=====================================
    integer                  :: k, j, l, m
    integer                  :: pe
    real(wp)    ,allocatable :: tmp (:,:)
    real(wp)    ,allocatable :: dh(:), ds(:), dv(:), dt(:), dq(:)
    type(t_vic) ,allocatable :: vic (:)
    type(t_ctl)              :: ctl

    if (dace% lpio) then
!     write(6,'()')
!     write(6,'(a)') repeat('-',79)
      write(6,'()')
      write(6,'(a)') '  set_bg_err_op'
    endif
    !-------------------------------------------
    ! check if covariance model is aready set up
    !-------------------------------------------
    if (covm% hc2_totalsize /= 0) then
      if (dace% lpio) then
        write(6,'(a)') '  returning, covariance model is aready set up.'
        write(6,'()')
      endif
      return
    endif
    !-------------------------------------
    ! get vertical covariances
    !   read NMC derived covariances
    !   or derive from OI covariance model
    !-------------------------------------
    if (file  == '') call read_bg_err_op_nml
    if (file  == '') return
    if (valid == 0)  return
    if (file  == 'OI') then
      call cov_from_oi      (time)
    else
      call read_nml_pscmodel   ! read old covariance model namelist
      call cov_from_nmc     (file, nx, ny, nz, pbot, h2psi, nwv, time, lclim, &
                             scale_rh)
      call adjust_top
      call read_OI_coef        ! read old covariance model coefficients
      call h_psi_cov_from_oi(time)
    endif
    !-------------------------------
    ! get horizontal covariances
    !   read NMC derived covariances
    !-------------------------------
    if (repr_2dh >= 1) then
       call read_horizontal_covariances (file_2dh,          &
                                         rnd_time,          &
                                         csr_symm=csr_symm, &
                                         csr_conv=csr_conv)
    endif
    !--------------------------------
    ! check for OI compatability mode
    !--------------------------------
    covm% oi_comp = 1
    if (associated (covm% mu)         .and. &
        associated (covm% sdev_psi)   .and. &
        associated (covm% sdev_psi_u) .and. &
        associated (covm% sdev_chi)   .and. &
        associated (covm% sdev_u)     .and. &
        associated (covm% sdev_v)           ) covm% oi_comp = 2
    if (dace% lpio) then
      write(6,*) 'oi_comp =',covm% oi_comp
      write(6,*)
    endif
    select case (covm% oi_comp)
    case default
      call finish ('set_bg_err_op','invalid value for oi_comp')
    case (1)
      allocate (covm% sdev_u (ny)); covm% sdev_u = 1._wp
      allocate (covm% sdev_v (ny)); covm% sdev_v = 1._wp
    case (2)
      covm% c_h_psi  = covm% sdev_psi   * covm% mu
      covm% c1_h_psi = covm% sdev_psi_u
      covm% L_h      = 1._wp
      covm% sqnu     = covm% sdev_chi
      covm% sq1nu    = 1._wp
      if (rmu >= -1._wp .and. rmu <= 1._wp) covm% c_h_psi = rmu
    end select
    !------------------------------------
    ! modify vertical covariance matrices
    !------------------------------------
    call plot_cov (ctl, covm, t=1)
    call modify_vcov   (covm)
    call plot_cov (ctl, covm, t=2)
    !--------------------------
    ! print modified parameters
    !--------------------------
    if (dace% lpio) then
      write(6,'()')
      write(6,'(a)') '  modified parameters:'
      write(6,'()')
    endif
    call print_bg_err_op_nml
    if (h2t) then
      allocate (vic (nz    ))
      allocate (tmp (nz, nz))
      call set_vertinpc (vic, covm% logp, covm% logp, covm% nwv)
    endif
!---------------------
! preliminary printout
!---------------------
if (dace% lpio) then
j = ny*1/6
write(6,*)
write(6,*)'set_bg_err_op: sqrt(diag) of covariance matrices at',covm%dlat (j)
write(6,*)
write(6,*)'    pressure level         height streamfunction     veloc.pot.    &
          &temperature       rel.hum.'
do k=1,nz
write(6,'(i4,6f15.3)') k, exp(covm%logp(k))/100, sqrt(covm% sqcvh(k,k,j)), &
sqrt(covm% sqcvs(k,k,j)), sqrt(covm% sqcvv(k,k,j)), &
sqrt(covm% sqcvt(k,k,j)), sqrt(covm% sqcvq(k,k,j))
end do
endif

    !==========================================
    ! final preparation of vertical covariances
    !==========================================
    if (sparse .and. .not. associated (covm% sqcv1h)) then
      allocate (covm% sqcv1h (nz, nz, ny))
      allocate (covm% sqcv1s (nz, nz, ny))
      allocate (covm% sqcv1q (nz, nz, ny))
      allocate (covm% sqcv1t (nz, nz, ny))
      allocate (covm% sqcv1v (nz, nz, ny))
    endif

    call p_barrier ! improve readability of printout in case of error
    do j=1,ny
      pe = mod (j, dace% npe)
      if  (pe == dace% pe) then
        !----------------
        ! constant errors
        !----------------
        if (e_h >= 0._wp) then
          call scale (covm% sqcvh(:,:,j), 1._wp / covm% eh   (:,j))
          covm% eh  (:,j) = e_h
          call scale (covm% sqcvh(:,:,j),         covm% eh   (:,j))
        endif
        if (e_v > 0._wp) then
          covm% ev  (:,j) = e_v
        endif
        if (e_q >= 0._wp) then
          call scale (covm% sqcvq(:,:,j), 1._wp / covm% erh  (:,j))
          covm% erh (:,j) = e_q
          call scale (covm% sqcvq(:,:,j),         covm% erh  (:,j))
        endif
        !-------------------------------------------
        ! derive temperature covariances from height
        !-------------------------------------------
        if (h2t) then
          tmp = 0._wp
          do l = 1, nz
            do m = 1, nz
              do k = 1, vic(l)% g% n
                tmp (l,m) = tmp (l,m) &
                          + vic(l)% wt(k) * covm% sqcvh(vic(l)% g% i+k-1,m,j)
              end do
            end do
          end do
          covm% sqcvt(:,:,j) = 0._wp
          do l = 1, nz
            do m = 1, nz
              do k = 1, vic(m)% g% n
                covm% sqcvt(l,m,j) = covm% sqcvt(l,m,j) &
                                   + vic(m)% wt(k) * tmp (l,vic(m)% g% i+k-1)
              end do
              covm% sqcvt(l,m,j) = covm% sqcvt(l,m,j) * (gacc / R) ** 2
              if (l==m) covm% et (l,j) = sqrt (covm% sqcvt(l,m,j))
            end do
          end do
        else
          if (all (covm% sqcvt == 0._wp)) call finish ('set_bg_err_op',  &
                                          'h2t==T but sqcvt is not valid')
        endif
        !-----------------
        ! scale the matrix
        !-----------------
        if (sqrtcor) then
          if (e_h /= 0._wp) call scale (covm% sqcvh(:,:,j), 1._wp / covm% eh   (:,j))
                            call scale (covm% sqcvs(:,:,j), 1._wp / covm% epsi (:,j))
                            call scale (covm% sqcvv(:,:,j), 1._wp / covm% exi  (:,j))
          if (e_q /= 0._wp) call scale (covm% sqcvq(:,:,j), 1._wp / covm% erh  (:,j))
          if (e_h /= 0._wp) call scale (covm% sqcvt(:,:,j), 1._wp / covm% et   (:,j))
          !-------------------------------------------------
          ! check for diagonal dominated covariance matrices
          !-------------------------------------------------
          if (e_h /= 0._wp) call checkdd (covm% sqcvh(:,:,j), j, 'sqcvh')
                            call checkdd (covm% sqcvs(:,:,j), j, 'sqcvs')
                            call checkdd (covm% sqcvv(:,:,j), j, 'sqcvv')
          if (e_q /= 0._wp) call checkdd (covm% sqcvq(:,:,j), j, 'sqcvq')
          if (e_h /= 0._wp) call checkdd (covm% sqcvt(:,:,j), j, 'sqcvt', softfail=.true.)
        endif
        !------------------------------
        ! set unit matrix (for testing)
        !------------------------------
        if (unit_v) then
          covm% sqcvh(:,:,j) = 0._wp
          covm% sqcvs(:,:,j) = 0._wp
          covm% sqcvv(:,:,j) = 0._wp
          covm% sqcvq(:,:,j) = 0._wp
          covm% sqcvt(:,:,j) = 0._wp
          do l = 1, nz
            covm% sqcvh(l,l,j) = 1._wp
            covm% sqcvs(l,l,j) = 1._wp
            covm% sqcvv(l,l,j) = 1._wp
            covm% sqcvq(l,l,j) = 1._wp
            covm% sqcvt(l,l,j) = 1._wp
          end do
        endif
        !----------
        ! transform
        !----------
        select case (transform)
        case ('n') ! none
        case ('w') ! wavelet transform
          basis = base
          call gp_mr_mat (covm% sqcvh(:,:,j))
          call gp_mr_mat (covm% sqcvs(:,:,j))
          call gp_mr_mat (covm% sqcvv(:,:,j))
          call gp_mr_mat (covm% sqcvq(:,:,j))
          call gp_mr_mat (covm% sqcvt(:,:,j))
        case default
         call finish('set_bg_err_op','invalid value for transform:'//transform)
        end select
        !------------
        ! derive sqrt
        !------------
        select case (sqr)
        case ('s') ! symmetric sqrt
          covm% sqcvh(:,:,j) = sqrt_rs (covm% sqcvh(:,:,j), min_ev = 0._wp)
          covm% sqcvs(:,:,j) = sqrt_rs (covm% sqcvs(:,:,j), min_ev = 0._wp)
          covm% sqcvv(:,:,j) = sqrt_rs (covm% sqcvv(:,:,j), min_ev = 0._wp)
          covm% sqcvq(:,:,j) = sqrt_rs (covm% sqcvq(:,:,j), min_ev = 0._wp)
          covm% sqcvt(:,:,j) = sqrt_rs (covm% sqcvt(:,:,j), min_ev = 0._wp)
        case ('c') ! Cholesky factorization
          call cholesky (covm% sqcvh(:,:,j))
          call cholesky (covm% sqcvs(:,:,j))
          call cholesky (covm% sqcvv(:,:,j))
          call cholesky (covm% sqcvq(:,:,j))
          call cholesky (covm% sqcvt(:,:,j))
        case default
          call finish('set_bg_err_op','invalid value for sqr:'//sqr)
        end select
        !----------------------------------------------------
        ! keep transformed matrices for sparse representation
        !----------------------------------------------------
        if (sparse) then
          covm% sqcv1h(:,:,j) = covm% sqcvh(:,:,j)
          covm% sqcv1s(:,:,j) = covm% sqcvs(:,:,j)
          covm% sqcv1q(:,:,j) = covm% sqcvq(:,:,j)
          covm% sqcv1t(:,:,j) = covm% sqcvt(:,:,j)
          covm% sqcv1v(:,:,j) = covm% sqcvv(:,:,j)
        endif
        !------------
        ! retransform
        !------------
        select case (transform)
        case ('n') ! none
        case ('w') ! wavelet transform
          call wave_1d (covm% sqcvh(:,:,j), TR_SYN, base)
          call wave_1d (covm% sqcvs(:,:,j), TR_SYN, base)
          call wave_1d (covm% sqcvv(:,:,j), TR_SYN, base)
          call wave_1d (covm% sqcvq(:,:,j), TR_SYN, base)
          call wave_1d (covm% sqcvt(:,:,j), TR_SYN, base)
        end select
        !----------
        ! normalize
        !----------
        if (.not. sqrtcor) then
          do k=1,nz
            if (e_h /= 0._wp) covm% sqcvh(k,:,j) = covm% sqcvh(k,:,j) * (1._wp / covm% eh  (k,j))
                              covm% sqcvs(k,:,j) = covm% sqcvs(k,:,j) * (1._wp / covm% epsi(k,j))
            if (e_q /= 0._wp) covm% sqcvq(k,:,j) = covm% sqcvq(k,:,j) * (1._wp / covm% erh (k,j))
            if (e_h /= 0._wp) covm% sqcvt(k,:,j) = covm% sqcvt(k,:,j) * (1._wp / covm% et  (k,j))
            covm% sqcvv(k,:,j) = covm% sqcvv(k,:,j) * (1._wp / covm% exi (k,j))
          end do
        end if
      endif
    end do
    do j=1,ny
      pe = mod(j,dace% npe)
      !------------------------------
      ! broadcast matrices to all PEs
      !------------------------------
      call p_bcast (covm% sqcvh(:,:,j), pe)
      call p_bcast (covm% sqcvs(:,:,j), pe)
      call p_bcast (covm% sqcvq(:,:,j), pe)
      call p_bcast (covm% sqcvt(:,:,j), pe)
      call p_bcast (covm% sqcvv(:,:,j), pe)
      call p_bcast (covm% eh     (:,j), pe)
      call p_bcast (covm% erh    (:,j), pe)
      call p_bcast (covm% et     (:,j), pe)
      if (sparse) then
        call p_bcast (covm% sqcv1h(:,:,j), pe)
        call p_bcast (covm% sqcv1s(:,:,j), pe)
        call p_bcast (covm% sqcv1q(:,:,j), pe)
        call p_bcast (covm% sqcv1t(:,:,j), pe)
        call p_bcast (covm% sqcv1v(:,:,j), pe)
      endif
    end do
    if (h2t) deallocate (vic, tmp)

!---------------------
! preliminary printout
!---------------------
if (dace% lpio) then
j = ny*1/6
allocate (tmp (nz,nz), dh(nz), ds(nz), dv(nz), dt(nz), dq(nz))
dh = diag (matmul (covm% sqcvh(:,:,j), transpose(covm% sqcvh(:,:,j))))
ds = diag (matmul (covm% sqcvs(:,:,j), transpose(covm% sqcvs(:,:,j))))
dv = diag (matmul (covm% sqcvv(:,:,j), transpose(covm% sqcvv(:,:,j))))
dt = diag (matmul (covm% sqcvt(:,:,j), transpose(covm% sqcvt(:,:,j))))
dq = diag (matmul (covm% sqcvq(:,:,j), transpose(covm% sqcvq(:,:,j))))

if (any (abs (dh - 1._wp) > 100._wp * EPSILON (1._wp)) .or. &
    any (abs (ds - 1._wp) > 100._wp * EPSILON (1._wp)) .or. &
    any (abs (dv - 1._wp) > 100._wp * EPSILON (1._wp)) .or. &
    any (abs (dt - 1._wp) > 100._wp * EPSILON (1._wp)) .or. &
    any (abs (dq - 1._wp) > 100._wp * EPSILON (1._wp))      ) then
write(6,*)
write(6,*)'set_bg_err_op: diagonals of correlation matrices at',covm%dlat (j)
write(6,*)
write(6,*)'    pressure level         height streamfunction     veloc.pot.    &
          &temperature       rel.hum.'
do k=1,nz
write(6,'(i4,6f15.3)')k,exp(covm%logp(k))/100,dh(k),ds(k),dv(k),dt(k),dq(k)
end do
end if

write(6,*)
write(6,*)'set_bg_err_op: Errors'
write(6,*)
write(6,*)'    pressure level              h             rh              t&
          &             xi              v'
do k=1,nz
write(6,'(i4,6f15.3)') k, exp(covm%logp(k))/100, covm% eh (k,j), &
covm% erh(k,j), covm% et (k,j), covm% exi(k,j), covm% ev(k,j)
end do
k=nz-nz/8
write(6,*)
write(6,*)'latitudinal dependence at k=',k, exp(covm%logp(k))/100
write(6,*)'set_bg_err_op: Errors:'
write(6,*)
write(6,*)'          latitude              h             rh              t&
          &             xi              v'
do j=1,ny
write(6,'(i4,6f15.3)') j,covm% dlat(j), covm% eh (k,j), covm% erh(k,j), &
covm% et (k,j), covm% exi(k,j), covm% ev(k,j)
end do

j=ny/6+1
write(6,*)
write(6,*)'analysis increment for near surface pressure observation at ',covm% dlat(j)
write(6,*)
write(6,*)'    p          cor h       cor t       cov h       cov t&
     &       var h       var t'
write(6,*)
dh=0; dh(nz-1)=1
dh=matmul (transpose(covm% sqcvh(:,:,j)), dh)
dt=matmul (          covm% sqcvt(:,:,j) , dh)
dh=matmul (          covm% sqcvh(:,:,j) , dh)
!dt=dt/dh(nz-1)
!dh=dh/dh(nz-1)
do k=1,nz
write(6,'(f9.3,6f12.5)') exp(covm%logp(k))/100, dh(k), dt(k), &
     dh(k)*covm% eh (k,j)*covm% eh (nz-1,j), &
     dt(k)*covm% et (k,j)*covm% eh (nz-1,j), &
     covm% eh (k,j), covm% et (k,j)
end do

deallocate (tmp, dh, ds, dv, dt, dq)
endif

!call p_barrier
!call finish('set_bg_err_op','test')

    !---------------------------------------
    ! finally switch on the covariance model
    !---------------------------------------
    covm% valid = valid
    if (dace% lpio) then
      write(6,'()')
      write(6,'(a)') '  set_bg_err_op finished.'
      write(6,'()')
    endif

    !------------------------------
    ! init old OI covariance model
    ! e.g. for forecast error model
    !------------------------------
    call init_fg_cov (time)
    call plot_cov (ctl, covm, t=3)
    call destruct (ctl)
    !--------------------------------
    ! store 3-d arrays on one PE only
    !--------------------------------
    call compress_cov

    contains
!..............................................................................
      subroutine adjust_top
      !---------------------------------------------
      ! adjust covariance matrix levels to model top
      !---------------------------------------------
      real(wp) :: l_adj   ! lower bound, logarithm
      real(wp) :: l_top   ! log (model top pressure)
      real(wp) :: padj    ! lower bound of adjustment zone (Pa)
      real(wp) :: a       ! coefficient
      real(wp) :: x       ! l_adj - log(p)
      real(wp) :: d       ! adjustment
      real(wp) :: p, l    ! values to adjust
      real(wp) :: s       ! scaling factor for temperature error
      integer  :: i       ! index

      if (p_adj <= 0._wp     ) return
      if (p_top <  covm% p(1)) then
        padj  = p_adj * 100._wp          ! hPa --> Pa
        l_top = log (p_top)
        l_adj = log (padj)
        if (dace% lpio) then
          write (6,'()')
          write (6,'(a)') repeat ('-',79)
          write (6,'()')
          write (6,'(a)') '  adjusting covariance matrix to model top'
          write (6,'()')
          write (6,'(a)') '     p (old,   new) log(p) (old,new)  scale'
        endif
        !----------------------
        ! adjust top level only
        !----------------------
        if (l_adj <= l_top) then
          if (dace% lpio) &
          write (6,'(i3,2f8.1,2f8.3)') 1,covm%p(1),p_top,covm%logp(1),l_top
          covm% p   (1) = p_top
          covm% logp(1) = l_top
        else
        !---------------------------
        ! adjust levels above 'padj'
        !---------------------------
          d = covm% logp(1) - l_top
          x = l_adj         - covm% logp(1)
          a = d / x**2
          do i = 1, covm% nz
            x = l_adj - covm% logp(i)
            if (x < 0) cycle
            d = a*x*x
            l = covm% logp(i)-d
            p = exp (l)
            s = 1._wp / (1._wp + 2._wp * a * x)
            if (dace% lpio) &
            write (6,'(i3,2f8.1,3f8.3)') i, covm% p(i), p, covm% logp(i), l, s
            covm% logp(i)   = l
            covm%    p(i)   = p
            covm%   et(i,:) = s * covm% et(i,:)
          end do
          if (dace% lpio) write (6,'()')
        endif
      endif

      end subroutine adjust_top
!..............................................................................
      subroutine checkdd (x, j, name, softfail)
      real(wp)         ,intent(in) :: x (:, :)
      integer          ,intent(in) :: j
      character(len=*) ,intent(in) :: name
      logical, optional,intent(in) :: softfail  ! Soft failure instead of error
        integer             :: k,l,m(1)
        real(wp)            :: d, d1, d2
        logical             :: error, soft
        real(wp) ,parameter :: eps = 1.e-5_wp
        !-------------------------------------------------
        ! check for diagonal dominated covariance matrices
        !-------------------------------------------------
        error = .false.
        soft = .false.; if (present (softfail)) soft = softfail
        do k=1,size(x,2)
          l = k
          d1 = maxval (x(:,k))
          d2 = maxval (x(k,:))
          d  = max (d1,d2) - x(k,k)
          if (d > eps) then
             error = .true.
             exit
          end if
        end do
        if (error) then
          if (d1 >= d2) then
             m = maxloc (x(:,k))
          else
             m = maxloc (x(k,:))
          end if
          write (6,*) repeat('=',79)
          write (6,'(4x,12i6)') (k,k=1,size(x,2))
          do k=1,size(x,2)
            write (6,'(i4,128f6.2)') k,x(:,k)
          end do
          write (6,*)
          write (6,*) 'checkdd: error in position j,k1,k2:',d,name,j,l,m
          write (6,*) repeat('=',79)
          if (soft) then
            call message('set_bg_err_op','error in checkdd: '//trim (name))
          else
            call finish('set_bg_err_op','error in checkdd')
          end if
        endif
      end subroutine checkdd
  end subroutine set_bg_err_op
!------------------------------------------------------------------------------
  subroutine modify_vcov (c)
  type (t_cov) ,intent(inout) :: c  ! covariance matrix data type
  !-------------------------------------------------------------------------
  ! modify vertical covariance matrices according to the namelist parameters
  !
  ! modify_cor: modify correlations, else variances
  ! lat_homog : latitudinal homogeneous covariances
  ! smooth_h  : smooth horizontally
  ! smooth_v  : smooth vertically
  !-------------------------------------------------------------------------
    integer               :: rversion = 0   ! Version of stdv_ratio
    real(wp) ,allocatable :: ratio  (:,:)
    real(wp) ,allocatable :: lev_rat(:)     ! Levels where 'ratio' is given
    integer               :: lev_num        ! Number of levels
    integer               :: ny, nz
    real(wp) ,pointer     :: L_vloc (:,:)
    integer  ,parameter   :: nlat = 361
    real(wp) ,parameter   :: lev_dflt(1:8) = (/999999._wp, &
         85000._wp,50000._wp,20000._wp,10000._wp,5000._wp,1000._wp, 0._wp/)
    !-------------
    ! any action ?
    !-------------
    if ( (lat_homog < -90._wp .or.  lat_homog > 90._wp) .and. &
          smooth_h  ==  0._wp .and. smooth_v ==  0._wp  .and. &
          stdv_ratio ==''                               .and. &
          i_vloc     == 0                               .and. &
          i_av_vcov  == 0                               .and. &
          p_const_v  == 0._wp                           .and. &
          p_const_q  == 0._wp                           .and. &
          emin_q     == 0._wp                                 ) return
    !----------------------
    ! vertical localisation
    !----------------------
    nullify (L_vloc)
    call vertloc (L_vloc, c% logp, i_vloc, p_vloc)
    !----------------------------------------------
    ! read stdev. ratios: forecast(3h) / NMC(48-24)
    !----------------------------------------------
    if (stdv_ratio /= '') then
      if (dace% lpio) then
        call read_ratio ()
      endif
      call p_bcast (lev_num, dace% pio)
      if (.not.dace% lpio) then
        allocate (lev_rat(     1:lev_num))
        allocate (ratio  (nlat,1:lev_num))
      end if
      call p_bcast (ratio,   dace% pio)
      call p_bcast (lev_rat, dace% pio)
    else
      !---------------------------
      ! Trivial fallback (ratio=1)
      !---------------------------
      lev_num    = size (lev_dflt)
      allocate (lev_rat(     1:lev_num))
      allocate (ratio  (nlat,1:lev_num))
      ratio      = 1._wp
      lev_rat(:) = lev_dflt(:)
    endif
    !--------------------------------------------------------------
    ! apply modification to each of the covariance matrices in turn
    !--------------------------------------------------------------
    ny = c% ny
    nz = c% nz
    call modify_cov (c ,c% sqcvh ,c% eh  ,c% eh_cl  ,0._wp    ,ratio) ! height
    call modify_cov (c ,c% sqcvs ,c% epsi,c% epsi_cl,0._wp    ,ratio) ! streamf.
    call modify_cov (c ,c% sqcvv ,c% exi ,c% exi_cl ,0._wp    ,ratio) ! vel.pot.
    call modify_cov (c ,c% sqcvq ,c% erh ,c% erh_cl ,p_const_q       &! humidity
                                                    ,errmin=emin_q)   !
    call modify_cov (c ,c% sqcvt ,c% et  ,c% et_cl  ,0._wp    ,ratio) ! temp.
    call rescale_error (          c% ev             ,p_const_v,ratio) ! wind
    !---------
    ! clean up
    !---------
    if (associated (L_vloc)) deallocate (L_vloc)
  contains
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    subroutine read_ratio ()
      integer            :: iu, i, j, k, ns, ios
      real(wp)           :: rlat
      character(len=255) :: line
      character(len=12)  :: str(25)
      real(wp)           :: lev(25)

      write(6,'()')
      write(6,'(a,a)') '  Reading ', trim (stdv_ratio)
      iu = get_unit_number()
      open (iu, file=stdv_ratio, action='read')
      lev_num = size (lev_dflt)           ! Default
      !---------------------------------
      ! Interprete or skip comment lines
      !---------------------------------
      line = ""
      do
         read (iu,'(A)') line
         if (line(1:1) /= "#") then
            backspace (iu)
            exit
         end if
         !---------------------------------------
         ! Check for format version specification
         !---------------------------------------
         i = index (line, 'version:')
         if (i > 0) then
            read (line(i+8:),*) rversion
            if (rversion < 0 .or. rversion > 3) then
               write(0,'(a)') trim (line)
               call finish ("read_ratio", "bad version number")
            end if
            cycle
         end if
         if (rversion == 2 .or. rversion == 3) then
            lev_num = 0
            lev = 0._wp
            i = index (line, 'Lat')
            if (i > 0) then
               str = ""
               call split (str, line(i+4:), ns)
               k = 1
               do j = 1, ns
                  read(str(j),*,iostat=ios) lev(k)
                  if (ios==0) k = k + 1
               end do
               lev_num = k + 1
               allocate (lev_rat(1:lev_num))
               lev_rat(1)   = lev_dflt(1)
               lev_rat(2:k) = lev(1:k-1) * 100._wp   ! hPa -> Pa
               lev_rat(k+1) = 0._wp
               if (rversion == 3) lev_rat(k+1) = 1.e-6_wp
            end if
         end if
      end do
      write(6,'(a,i0)') '  File format version: ', rversion

      select case (rversion)
      case (0,1)
         allocate (lev_rat   (1:lev_num))
         lev_rat(:) = lev_dflt(:)          ! Work in progress...
      case (2,3)
         if (lev_num <= 2) then
            write(6,'(a,i0)') '  Insufficient number of levels found: ', lev_num
            call finish ("read_ratio", "Insufficient number of levels")
         else
            write(6,'(a,99f9.3)') '  Levels provided [hPa]: ', &
                 lev_rat(2:lev_num-1)/100._wp
         end if
      end select

      allocate (ratio(nlat,1:lev_num))
      do i=1,size(ratio,1)
         read (iu,*) rlat, ratio(i,2:lev_num-1)
      end do
      select case (rversion)
      case (0)
         !-------------------------------------------
         ! Used fixed scaling factors in stratosphere
         !-------------------------------------------
         ratio(:,6) = 1.0_wp
         ratio(:,7) = 1.5_wp
      case (1)
         !--------------------------------------------------------------
         ! Prevent scaling factors only from becoming unreasonably large
         !--------------------------------------------------------------
         ratio(:,6) = min (ratio(:,6), 1.0_wp)
         ratio(:,7) = min (ratio(:,7), 1.5_wp)
      case (2,3)
         !-----------------
         ! No limits so far
         !-----------------
      case default
      end select
      ratio(:,1)       = ratio(:,2)
      ratio(:,lev_num) = ratio(:,lev_num-1)
      close (iu)
      call return_unit_number (iu)
    end subroutine read_ratio
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    subroutine modify_cov (c, cov, err, clerr, p_const, ratio, errmin)
    !------------------------------------------------
    ! apply modifications to single covariance matrix
    !------------------------------------------------
    type(t_cov) ,intent(in) :: c              ! covariance meta data
    real(wp)    ,pointer    :: cov    (:,:,:) ! covariance matrix
    real(wp)    ,pointer    :: err    (:,:)   ! error
    real(wp)    ,pointer    :: clerr  (:,:)   ! climatological error
    real(wp)                :: p_const        ! constant below p (hPa)
    real(wp)    ,optional   :: ratio  (:,:)   ! ratio stdev(3h)/stdev(NMCtime)
    real(wp)    ,optional   :: errmin         ! global lower bound on error
      integer  :: i, j,  k, jj(1)
      real(wp) :: w (size(cov,3),2)
      !----------------
      ! check arguments
      !----------------
      if (.not.associated(cov)) return
      if (.not.associated(err)) return

      !-----------------------------------------
      ! horizontally average covariance matrices
      !-----------------------------------------
      if (i_av_vcov /= 0) then
        w(:,1) = cos (c% dlat * d2r)
        w(:,2) = w(:,1) * (2 * i_av_vcov + 1)
        call average_merid (w(:,2), i_av_vcov)
        do i = 1, nz
          clerr (i, :) = clerr (i, :) ** 2
          call average_merid (clerr (i, :), i_av_vcov, w)
          clerr (i, :) = sqrt (clerr (i, :))
          do k = 1, nz
            call average_merid (cov (k, i, :), i_av_vcov, w)
            if (i==k) err(i,:) = sqrt (cov (k, i, :))
          end do
        end do
      end if

      !-----------------------------------------------------
      ! apply modifications to correlations, not covariances
      !-----------------------------------------------------
!     if (modify_cor) then
        do j = 1, ny
          do i = 1, nz
            do k = 1, nz
              cov (k, i, j) = cov (k, i, j) / (err(i,j) * err(k,j))
            end do
          end do
        end do
!     endif

      !-------------------------------------------------------------
      ! rescale variances by stdev. ratio: forecast(3h) / NMC(48-24)
      !-------------------------------------------------------------
      call rescale_error (err, p_const, ratio)

      !------------------------------------
      ! limit error by global minimum value
      !------------------------------------
      if (present(errmin)) err = max (err, errmin)

      !------------------------------------
      ! latitudinal homogeneous covariances
      !------------------------------------
      if (abs (lat_homog) <= 90._wp) then
        jj = minloc (abs (lat_homog - c% dlat))
        j  = jj (1)
        do i = 1, nz
          err (i,:)                         = err   (i,j)
          if(associated(clerr)) clerr (i,:) = clerr (i,j)
          do k = 1, nz
            cov (k, i, :) = cov (k, i, j)
          end do
        end do
      endif
      !-----------------------
      ! rescale to covariances
      !-----------------------
!     if (modify_cor) then
        do j = 1, ny
          do i = 1, nz
            do k = 1, nz
              cov (k, i, j) = cov (k, i, j) * (err(i,j) * err(k,j))
            end do
          end do
        end do
!     else
!      do j = 1, ny
!         do k = 1, nz
!           err(k,j) = sqrt (cov (k, k, j))
!         end do
!       end do
!     endif
      !----------------------
      ! vertical localisation
      !----------------------
      if (i_vloc /= 0) then
        do j = 1, ny
          cov (:,:,j) = cov (:,:,j) * L_vloc
        end do
      endif
    end subroutine modify_cov
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    subroutine rescale_error (err, p_const, ratio)
    real(wp)    ,pointer    :: err   (:,:)   ! error
    real(wp)                :: p_const       ! constant below p (hPa)
    real(wp)    ,optional   :: ratio (:,:)   ! ratio stdev(3h)/stdev(NMCtime)
      !-------------------------------------------------------------
      ! rescale variances by stdev. ratio: forecast(3h) / NMC(48-24)
      !-------------------------------------------------------------
      integer  :: ik(c% nz), ij(c% ny)
      real(wp) :: wk(c% nz), wj(c% ny)
      integer  :: j, k, kref
      real(wp) :: r

      !---------------------------------
      ! constant variances below p_const
      !---------------------------------
      if (p_const > 0._wp) then
        kref = 0
        do k = 1, nz
          if (covm% p(k) < p_const * 100._wp) then
            kref = k
          else
            if (kref == 0) then
              write (0,*)  'modify_cov:  invalid p_const =',p_const
              call finish ('modify_cov','invalid p_const')
            endif
            err(k,:) = err(kref,:)
          endif
        end do
      endif

      if (.not.present(ratio)) return

      !----------------------------------
      ! meridional index (-90..90 dy=0.5)
      !----------------------------------
      do j=1,ny
        ij(j) = int ((c% dlat(j) + 90._wp) * 2) + 1
        wj(j) = (c% dlat(j) + 90._wp) * 2 - int ((c% dlat(j) + 90._wp) * 2)
      end do
      !--------------------------------------
      ! vertical index (ca.1000 hPa .. 0 hPa)
      !--------------------------------------
      do k=1,nz
        ik(k) = sum (minloc (lev_rat - c% p(k), (lev_rat-c%p(k))>=0._wp))
        !--------------------------------------------------------------
        ! Interpolation linear in p (version<3) or in ln(p) (version=3)
        !--------------------------------------------------------------
        if (rversion < 3) then
          wk(k) = (lev_rat(ik(k))-c%p(k)) / (lev_rat(ik(k))-lev_rat(ik(k)+1))
        else
          wk(k) = log (lev_rat(ik(k)) / c% p(k))        &
                / log (lev_rat(ik(k)) / lev_rat(ik(k)+1))
        end if
      end do
      !--------------------------
      ! rescale error (variances)
      !--------------------------
      do j=1,ny
        do k=1,nz
          r = ratio (ij(j)   ,ik(k)  ) * (1._wp-wj(j)) * (1._wp-wk(k)) &
            + ratio (ij(j)+1 ,ik(k)  ) *        wj(j)  * (1._wp-wk(k)) &
            + ratio (ij(j)   ,ik(k)+1) * (1._wp-wj(j)) *        wk(k)  &
            + ratio (ij(j)+1 ,ik(k)+1) *        wj(j)  *        wk(k)
          err (k,j) = r * err (k,j)
        end do
      end do

    end subroutine rescale_error
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine modify_vcov
!------------------------------------------------------------------------------
  subroutine average_merid (cov, i_av_vcov, weight)
  real(wp) ,intent(inout)        :: cov (:)     ! vertical covariance matrices
  integer  ,intent(in)           :: i_av_vcov   !+- number of points to average
  real(wp) ,intent(in) ,optional :: weight(:,:) ! precalculated weights
  !------------------------------------------------------------
  ! meridionally average the vertical covariance matrices using
  ! + - i_av_vcov points (i.e. 2*i_av_vcov+1 points in total)
  !------------------------------------------------------------
    integer  :: j, ny, n
    real(wp) :: x (size(cov))

    if (i_av_vcov == 0) return
    ny = size(cov)
    n  = 2 * i_av_vcov + 1

    if (present(weight)) then
      x = cov * weight (:,1)
      do j = 1, i_av_vcov
        cov (j) = sum (x(1:n)) / weight (j,2)
      end do
      do j = i_av_vcov + 1, ny - i_av_vcov
        cov (j) = sum (x(j-i_av_vcov:j+i_av_vcov)) / weight (j,2)
      end do
      do j = ny - i_av_vcov + 1, ny
        cov (j) = sum (x(ny-n+1:ny)) / weight (j,2)
      end do
    else
      x = cov
      do j = 1, i_av_vcov
        cov (j) = sum (x(1:n)) / n
      end do
      do j = i_av_vcov + 1, ny - i_av_vcov
        cov (j) = sum (x(j-i_av_vcov:j+i_av_vcov)) / n
      end do
      do j = ny - i_av_vcov + 1, ny
        cov (j) = sum (x(ny-n+1:ny)) / n
      end do
    endif
  end subroutine average_merid
!------------------------------------------------------------------------------
  subroutine vertloc (L, logp, i_vloc, p_vloc)
  real(wp) ,pointer     :: L  (:,:)  ! localisation matrix
  real(wp) ,intent(in)  :: logp (:)  ! vertical coordinate
  integer  ,intent(in)  :: i_vloc    ! flag
  real(wp) ,intent(in)  :: p_vloc(2) ! bounds of localisation interval
  !----------------------------------------
  ! set up (vertical) localisation matrix L
  !----------------------------------------

    real(wp) :: logp_low
    real(wp) :: logp_high
    real(wp) :: w (size(logp),2)
    real(wp) :: z (size(logp))
    integer  :: j, k
    integer  :: n

    if (associated (L)) deallocate (L)
    select case (i_vloc)
    case (0)
      !-----------------------
      ! no localisation: L = 1
      !-----------------------
    case (1,2)
      !------------------------------------
      ! derive bounds for localisation area
      !------------------------------------
      logp_low  = minval(p_vloc) * 100._wp
      logp_high = maxval(p_vloc) * 100._wp
      if (logp_low <=0._wp) then
        logp_low = minval (logp)
      else
        logp_low = max (log(logp_low),  minval(logp))
      endif
      logp_high  = min (log(logp_high), maxval(logp))
      !---------------------------
      ! derive vertical coordinate
      !---------------------------
      z = 2._wp * (logp - logp_low) / (logp_high - logp_low) - 1._wp
      if (i_vloc == 1) &
        where (z > -1._wp .and. z < 1._wp) z = sin (pi/2._wp * z)
      !------------------------------------
      ! derive weights for ensemble 1 and 2
      !------------------------------------
      n = size(logp)
      do k = 1, n
        if (z(k) <= -1._wp) then
          w (k,1) = 0._wp; w (k,2) = 1._wp
        else if (z(k) >= 1._wp) then
          w (k,1) = 1._wp; w (k,2) = 0._wp
        else
          w(k,1) = cos (PI/4*(1-z(k)))
          w(k,2) = sin (PI/4*(1-z(k)))
        endif
!if(dace% lpio) then
!write(6,'(a,i3,4f10.3)') 'vertloc: k,p,z,w12 =',k,exp(logp(k)),z(k),w(k,:)
!end if
      end do
      allocate (L (n,n))
      do k = 1, n
        do j = 1, n
          L (j,k) = w(j,1)*w(k,1) + w(j,2)*w(k,2)
        end do
      end do
      ! Some simple sanity checks:
      if (maxval (L) > 1 + EPSILON (0._wp)) then
         write (*,*) "vertloc: maxval (L) =", maxval (L), "@", maxloc (L)
      end if
      if (minval (L) < 0) then
         write (*,*) "vertloc: minval (L) =", minval (L), "@", minloc (L)
      end if
    case default
      call finish ('vertloc','invalid value of i_vloc')
    end select
  end subroutine vertloc
!------------------------------------------------------------------------------
  subroutine plot_cov (ctl, c, t)
  type(t_ctl) ,intent(inout) :: ctl ! GRADS .ctl file info
  type(t_cov) ,intent(in)    :: c   ! covariance matrix info
  integer     ,intent(in)    :: t   ! time level
  !----------------------------------------------
  ! write GRADS file with variances actually used
  !----------------------------------------------
    integer  :: nz
    logical  :: zrev
    real(wp) :: zlev(c% nz)         ! Temporary: pressure levels (hPa)

    nz   = c% nz
    zrev = c% p(1) < c% p(nz)
    if (zrev) then
       zlev = c% p(nz:1:-1)/100
    else
       zlev = c% p         /100
    end if

    if (dace% lpio) then
      if (t==1)                                               &
        call init_ctl(ctl, path_file(aux,'^3dvar_variances'), &
                           title='variances used in 3dvar',   &
                           nx=1, zlev=zlev, ngl= c% ny,       &
                           comment =                          &
                        (/'time slices: 1: original stddevs', &
                          '             2: modified stddevs', &
                          '             3: final stddevs   '/))
      call write_var (ctl, transpose(c% eh)   ,'stdev_h'   ,t=t,zrev=zrev)
      call write_var (ctl, transpose(c% epsi) ,'stdev_psi' ,t=t,zrev=zrev)
      call write_var (ctl, transpose(c% exi)  ,'stdev_chi' ,t=t,zrev=zrev)
      call write_var (ctl, transpose(c% erh)  ,'stdev_rh'  ,t=t,zrev=zrev)
      call write_var (ctl, transpose(c% et)   ,'stdev_t'   ,t=t,zrev=zrev)
      call write_var (ctl, transpose(c% ev)   ,'stdev_v'   ,t=t,zrev=zrev)
      call write_ctl (ctl)
    endif

  end subroutine plot_cov
!------------------------------------------------------------------------------
  subroutine cov_from_oi (time)
  !-----------------------------------------------------------
  ! derive the operator representation of vertical covariances
  ! from the explicit OI formulation
  !-----------------------------------------------------------
  type (t_time) ,intent(in) :: time ! time and date of analysis
  !-------------------------------------------
  ! set the covariance matrix meta information
  ! from the OI model
  !-------------------------------------------
    type(t_rowcol)     :: col       ! OI bg model meta data
    type(t_coord)      :: coord     ! coordinate data meta type
    integer            :: j         ! latitude index
    integer            :: l         ! quantity index
    integer            :: it(nz)    ! observation type index
    real(mp)           :: c (nz,nz) ! covariance matrix
    real(wp)           :: e (nz)    ! error (stdev)
    integer, parameter :: its (5) = (/OBS_H, OBS_TV, OBS_RH, OBS_U, OBS_CHI/)
    real(wp)           :: rlsvv     ! renormalisation due to large scale c.

    if (dace% lpio) then
      write(6,'()')
      write(6,'(a)') '  cov_from_oi'
    endif
    if (ny<=0) ny = 61
    call construct (covm, nx, ny, nz, pbot * 100._wp, nwv, time, 0)

    !-------------------------------------
    ! take parameters mu, nu from OI model
    !-------------------------------------
    if (associated (covm% c_h_psi)) &
      call finish('cov_from_oi','c_h_psi already associated !')
    allocate       (covm% c_h_psi  (covm%ny))
    allocate       (covm% c1_h_psi (covm%ny))
    allocate       (covm% L_h      (covm%ny))
    allocate       (covm% sqnu     (covm%ny))

    !-------------------------------
    ! initialize OI covariance model
    ! set up OI bg meta data
    !-------------------------------
    call init_fg_cov (time)
    call construct (col, covm% nz, 1)

    !++++++++++++++++++++++++++++++++++++++++++++++
    ! reduced wind mass correlation, testing only !
    !++++++++++++++++++++++++++++++++++++++++++++++
    rlsvv  = 1._wp / (1._wp + lsvv)
    if (dace% lpio .and. rlsvv /= 1._wp) then
      write(6,'()')
      write(6,'(a)') repeat('!',79)
      write(6,'(a)') ' lsvv /= 0 in wavelet B approach'
      write(6,'(a)') repeat('!',79)
    endif

    !-------------------------------------------
    ! derive covariance matrix for each latitude
    !-------------------------------------------
    coord% dlon = 0._wp
    do l = 1, size(its)
      it = its(l)
      do j=1,ny
        c           = 0
        coord% dlat = covm% dlat(j)
        call set_xuv (coord)
        call zero (col)
        call init_spot (col,         &
                        1           ,&! spot index within the box
                        coord       ,&! coordinates
                        it(:)       ,&! observation type identifiers
                        covm%  logp ,&! vertical coordinates (ln p)
                    lcl=IFC_CFC     ,&! flag for climatological fc error
                    err=e)            ! error (stdev)
        call get_corr (&
                       col,    &! lhs description
                       col,    &! rhs description
                       1,      &! row spot index
                       1,      &! column spot index
                       .true., &! .T.:covariance, .F.:correlation
                       .true., &! .T. for symmetric matrix
                       c=c)     ! covariance matrix
        select case (its(l))
        case (OBS_H)
          covm% sqcvh (:,:,j) = c
          covm% sqcvs (:,:,j) = c
          covm% eh      (:,j) = e
          covm% epsi    (:,j) = e
          covm% c_h_psi   (j) = col% spots(1)% mu  * rlsvv
          covm% L_h       (j) = col% spots(1)% L_h
        case (OBS_TV)
          covm% sqcvt (:,:,j) = c
          covm% et      (:,j) = e
        case (OBS_RH)
          covm% sqcvq (:,:,j) = c
          covm% erh     (:,j) = e
        case (OBS_U)
          covm% ev      (:,j) = e
        case (OBS_CHI)
          covm% sqcvv (:,:,j) = c
          covm% exi     (:,j) = e
        end select
      end do
    end do

    if (rmu >= -1._wp .and. rmu <= 1._wp) covm% c_h_psi =       rmu
    covm% sqnu     = sqrt (rnu)
    covm% sq1nu    = sqrt (1._wp - covm% sqnu(1) ** 2)
    covm% c1_h_psi = sqrt (1._wp - covm% c_h_psi ** 2)

    !--------
    ! cleanup
    !--------
    call destruct (col)

  end subroutine cov_from_oi
!-----------------------------------------------------------------------------
  subroutine h_psi_cov_from_oi(time)
  type (t_time) ,intent(in) :: time ! time and date of analysis
  !------------------------------------------
  ! just take parameters mu, nu from OI model
  !------------------------------------------
    integer            :: j             ! latitude index
    real(wp)           :: dlat, wlat    ! interpolation coefficients
    integer            :: jlat, jlat1   ! .. for mu
    integer            :: fg_month      ! month
    real(wp)           :: rlsvv         ! renormalisation due to large scale c.

    !++++++++++++++++++++++++++++++++++++++++++++++
    ! reduced wind mass correlation, testing only !
    !++++++++++++++++++++++++++++++++++++++++++++++
    rlsvv  = 1._wp / (1._wp + lsvv)
    if (dace% lpio) then
      write(6,'()')
      write(6,'(a)') '  h_psi_cov_from_oi'
      if (rlsvv /= 1._wp) then
        write(6,'()')
        write(6,'(a)') repeat('!',79)
        write(6,'(a)') ' lsvv /= 0 in wavelet B approach'
        write(6,'(a)') repeat('!',79)
      endif
    endif

    if (associated (covm% c_h_psi)) &
      call finish('h_psi_cov_from_oi','c_h_psi already associated !')
    allocate       (covm% c_h_psi  (covm%ny))
    allocate       (covm% c1_h_psi (covm%ny))
    allocate       (covm% L_h      (covm%ny))
    allocate       (covm% sqnu     (covm%ny))
    !-------------------------------
    ! set up OI bg meta data
    !-------------------------------
    fg_month  = imm (time)

    !----------------------------
    ! derive mu for each latitude
    !----------------------------
    do j=1,ny
      !---------------------------------------------------------------
      ! indices and weights for latitude dependent correlations
      ! 6degree x 6degree lat-lon grid for analysis error is hardcoded
      !---------------------------------------------------------------
      dlat = covm% dlat(j)
      wlat  =  min (30._wp, max (1._wp, (dlat+87._wp)/6._wp+1))
      jlat  =  int(wlat)
      wlat  =  wlat - jlat
      jlat1 =  min (30, jlat+1)
      covm% c_h_psi(j)  =(       wlat  * nfcoun% mu (jlat1, fg_month) &
                        + (1._wp-wlat) * nfcoun% mu (jlat,  fg_month))&
                        * rlsvv
      covm% L_h    (j)  =        wlat  * nfcoun% b  (jlat1, fg_month) &
                        + (1._wp-wlat) * nfcoun% b  (jlat,  fg_month)
    end do

    if (L_h >   0._wp                   ) covm% L_h     = L_h
    if (rmu >= -1._wp .and. rmu <= 1._wp) covm% c_h_psi = rmu
    covm% sqnu     = sqrt (rnu)
    covm% sq1nu    = sqrt (1._wp - covm% sqnu(1) ** 2)
    covm% c1_h_psi = sqrt (1._wp - covm% c_h_psi ** 2)

  end subroutine h_psi_cov_from_oi
!=============================================================================
! Apply the Covariance matrix
!=============================================================================
  subroutine set_ic_obs (obs)
  type (t_obs_set) ,intent(inout) :: obs  ! observation meta data
  !--------------------------------------------------------
  ! set interpolation coefficients in observation data type
  !--------------------------------------------------------
    integer               :: ib, is, ic, ics
    integer               :: k
    type(t_rowcol)        :: col              ! OI bg model meta data
    type(t_spot) ,pointer :: si               ! coordinate meta data

    !----------------
    ! loop over boxes
    !----------------
    do ib = 1, size(obs% o)
      if (obs% o(ib)% pe /= dace% pe) cycle
      if (obs% o(ib)% n_lev /= 0) cycle
      call construct (col, obs% o(ib)% n_int, obs% o(ib)% n_spt)
      !------------------
      ! loop over reports
      !------------------
      ics = 1
      do is = 1, obs% o(ib)% n_spot
        si => obs% o(ib)% spot(is)
        call init_spot_obs (col, si, obs% o(ib))
        !-----------------------------------------------------
        ! set horizontal correlation coefs. in spot% hic
        !                                   or spot% imcol% hi
        !-----------------------------------------------------
        select case (si% int_type)
        case (ITY_ICOL)
          call set_horinpc (si% hic,          &
                            si% col% c% dlat, &
                            covm)
          si% emod = col% spots(ics)% emod
          ics = ics + 1
        case (ITY_MCOLS)
          do ic = 1, size(si% imcol)
            call set_horinpc (si% imcol(ic)%   hi,   &
                              si% imcol(ic)%c% dlat, &
                              covm)
            si% imcol(ic)% emod = col% spots(ics)% emod
            ics = ics + 1
          end do
        case default
          call finish('set_ic_obs','invalid int_type')
        end select
      end do
      !----------------------------------------------
      ! set vertival correlation coefs. in spot% levs
      !----------------------------------------------
      call col2lev (obs% o(ib))
      do k = 1, obs% o(ib)% n_lev
        call set_vertinpc (obs% o(ib)% levs(k)% vi, &
                           obs% o(ib)% levs(k)% z,  &
                           covm% logp, covm% nwv)
        obs%o(ib)%levs(k)%vi% ezn  = col% cols(obs% o(ib)% levs(k)%i%i+1)% ezn
        obs%o(ib)%levs(k)%vi% exzn = col% cols(obs% o(ib)% levs(k)%i%i+1)% exzn
      end do
      call destruct (col)
    end do

  end subroutine set_ic_obs
!-----------------------------------------------------------------------------
  subroutine apply_B_ii_1d (obs, x, y, e_fi)
  type (t_obs_set) ,intent(inout) :: obs  ! observation meta data
  type (t_vector)  ,intent(in)    :: x    ! input
  type (t_vector)  ,intent(inout) :: y    ! output
  type (t_vector)  ,intent(in)    :: e_fi ! background error
  !---------------------------------------------------------
  ! multiply  a vector in interpolation space by the B matrix
  ! result is a vector in interpolation space
  !---------------------------------------------------------
    integer                :: ib
    type (t_mode)          :: ml (size(obs%o)) ! vertical modes      lhs
    type (t_mode)          :: mr (size(obs%o)) ! vertical modes      rhs
    call enter_function
    call trace_mem_usage ('B_ii_1d:entry')
    call uncompress_cov
    call trace_mem_usage ('B_ii_1d:uncompress')
    !--------------------------------------------------------
    ! set interpolation coefficients in observation data type
    !--------------------------------------------------------
    call set_ic_obs (obs)
    !--------------------------------
    ! allocate temporary storage, rhs
    !--------------------------------
    call construct (mr, covm% nz, obs% o% n_spt)
    call print_mem_usage ('B_ii_1d:rhs','rhs', sum (mr% size))
    !--------------------------------------------------
    ! adjoint interpolation, multiply obs errors stdev.
    !--------------------------------------------------
    call vertinpa (x, mr, obs, e_fi)
    !--------------------------
    ! project on vertical modes
    !--------------------------
    call int2vmode (mr, obs)
    !------------
    ! scatter rhs
    !------------
    call scatter_modes (mr, obs)
    !--------------------------------
    ! allocate temporary storage, lhs
    !   take usage flags from rhs
    !--------------------------------
    call construct (ml, covm% nz, obs% o% n_spt, obs% o% pe == dace% pe)
    call print_mem_usage ('B_ii_1d:lhs','lhs', sum (ml% size))
    do ib = 1, size(obs% o)
      if (obs% o(ib)% pe /= dace% pe) cycle
      ml(ib)% lh  = mr(ib)% lh
      ml(ib)% luv = mr(ib)% luv
      ml(ib)% lq  = mr(ib)% lq
      if (l1dvar) then
        ml(ib)% h  = mr(ib)% h
        ml(ib)% u  = mr(ib)% u
        ml(ib)% v  = mr(ib)% v
        ml(ib)% q  = mr(ib)% q
        ml(ib)% ux = mr(ib)% ux
        ml(ib)% vx = mr(ib)% vx
      endif
    end do
    !-----------------------------------
    ! apply horizontal covariance matrix
    !-----------------------------------
    if (.not. l1dvar) call horicov_oo (mr, ml, obs)
    !----------------------------
    ! project from vertical modes
    !----------------------------
    call vmode2int_o (ml, obs)
    !--------------------------------------------------
    ! forward interpolation, multiply obs errors stdev.
    !--------------------------------------------------
    call vertinp_o (y, ml, obs, e_fi, x)
    !-----------------------------
    ! deallocate temporary storage
    !-----------------------------
    call trace_mem_usage ('B_ii_1d:total')
    call destruct (ml)
    call destruct (mr)
    call delete_storage (x)
    call delete_storage (y)
    call delete_storage (e_fi)
    call compress_cov
    call leave_function
  end subroutine apply_B_ii_1d
!-----------------------------------------------------------------------------
  subroutine apply_B_mi_1d (a_m, cbg, obs, x, lnewpl, e_fi)
  type(t_cols)     ,intent(out)   :: a_m(:) ! analysis increment (model space)
  type(t_cols)     ,intent(in)    :: cbg(:) ! reference state
  type (t_obs_set) ,intent(inout) :: obs    ! observation meta data
  type (t_vector)  ,intent(in)    :: x      ! input
  logical          ,intent(in)    :: lnewpl ! analysis increment on new p-levs
  type (t_vector)  ,intent(in)    :: e_fi   ! background error
  !---------------------------------------------------------
  ! multiply a vector in interpolation space by the B matrix
  ! result is on model grid
  !---------------------------------------------------------
    integer         ,parameter :: nm = 4  ! number of multilevel fields to ass.
    integer         ,parameter :: ns = 1  ! number of singlelevel fields
    integer                    :: ke      ! no. model levels
    integer                    :: ncol    ! no. columns
    integer                    :: nvc     ! no. variables/column
    integer                    :: ib      ! box index
    integer                    :: is      ! spot index
    integer                    :: k       ! vertical index
    integer                    :: i1, in  ! observation index bounds
    integer                    :: l1, ln  ! level index bounds
    type (t_mode)              :: ml                ! vertical modes      lhs
    type (t_mode)              :: mr (size(obs%o))  ! vertical modes      rhs
    type (t_hic1) ,allocatable :: hic(:)            ! horizontal interp. coefs
    type(t_rowcol)             :: col               ! OI bg model meta data
    type(t_rowcol)             :: ps                ! OI bg model meta data
    type(t_lev)   ,allocatable :: levs  (:)
    integer       ,allocatable :: t_int (:)
    real(wp)      ,allocatable :: e     (:)
    real(wp)      ,allocatable :: z     (:)
    real(wp)      ,allocatable :: lpf   (:), ph(:)
    real(wp)      ,allocatable :: y     (:)
    real(wp)      ,pointer     :: psx   (:)
    type(t_coord) ,pointer     :: c

    call enter_function
    call trace_mem_usage ('B_mi:entry')
    call uncompress_cov
    call trace_mem_usage ('B_mi:uncompress')

    if (dace% lpio) then
      write(6,'(a)') repeat('-',79)
      write(6,'()')
      write(6,'(a,l1)')'  apply_B_mi_1d newpl =',lnewpl
      write(6,'()')
    endif
    !--------------------------------------------------------
    ! set interpolation coefficients in observation data type
    !--------------------------------------------------------
    call set_ic_obs (obs)
    !---------------------------
    ! allocate temporary storage
    !---------------------------
    call construct (mr, covm% nz, obs% o% n_spot)
    call print_mem_usage ('B_mi:rhs','rhs', sum (mr% size))
    !--------------------------------------------------
    ! adjoint interpolation, multiply obs errors stdev.
    !--------------------------------------------------
    call vertinpa (x, mr, obs, e_fi)
    !--------------------------
    ! project on vertical modes
    !--------------------------
    call int2vmode (mr, obs)
    !------------
    ! scatter rhs
    !------------
    call scatter_modes (mr, obs)

    call construct (col, nm*cbg(1)%ke+ns, 1)
    call construct (ps,  1,               1)

    do ib = 1, size(cbg)
      !------------------------
      ! printout, communication
      !------------------------
      call p_barrier
      if (dace% lpio) write(6,*) &
        'apply_B_mi_1d: box, nbox, ncol',ib,size(cbg),cbg(ib)% ncol
      call p_barrier
      if(cbg(ib)% ncol==0) cycle
      ke   = cbg(ib)% ke
      ncol = cbg(ib)% ncol
      nvc  = nm * ke + ns
      !--------------------------------
      ! allocate temporary storage, rhs
      !--------------------------------
      call construct (ml, covm% nz, ncol, .true.)
      ml% lh  = .true.
      ml% luv = .true.
      ml% lq  = .true.
      !-----------------------------------
      ! apply horizontal covariance matrix
      !-----------------------------------
      call horicov (mr, ml, obs, cbg(ib)% col% c)
      !----------------------------
      ! project from vertical modes
      !----------------------------
      allocate (hic (ncol))
      call set_horinpc (hic, cbg(ib)% col% c% dlat, covm)
      call vmode2int (ml, hic)
      !=================================================
      ! forward interpolation, multiply obs error stdev.
      !=================================================
      !---------------
      ! first pass: ps
      !---------------
      if (lnewpl) then
        !---------------------
        ! allocate temporaries
        !---------------------
        allocate (t_int (ncol))
        allocate (e     (ncol))
        allocate (z     (1)   )
        allocate (levs  (ncol))
        allocate (y     (ncol))
        !-----------------
        ! set up meta data
        !-----------------
        t_int = OBS_H
        do is = 1, ncol
          c =>      cbg(ib)% col(is)% c
          z =  log (cbg(ib)% col(is)% s% ps)
          call set_vertinpc (levs(is)% vi, z(1), covm% logp, covm% nwv)
          call zero (ps)
          call init_spot (ps,        &
                          1         ,&! spot index within the box
                          c         ,&! coordinates
                          (/OBS_H/) ,&! observation type identifiers
                          z         ,&! vertical coordinates (ln p)
                      err=e(is:is))   ! error (stdev)
          levs(is)%     is   = is
          levs(is)% i % i    = is - 1
          levs(is)% i % n    = 1
          levs(is)%     z    = z(1)
          levs(is)% vi% ezn  = ps% cols(1)% ezn
          levs(is)% vi% exzn = ps% cols(1)% exzn
        end do
        !------------
        ! interpolate
        !------------
        call vertinp (y, ml, t_int, levs, e, y)
        !-----------------------
        ! deallocate temporaries
        !-----------------------
        deallocate (t_int)
        deallocate (e    )
        deallocate (z    )
        deallocate (levs )
      endif
      !------------------------------------------------------------
      ! allocate result variable, set to zero, patch pressure level
      !------------------------------------------------------------
      call alloc_cols (a_m(ib), tmp=cbg(ib), ids=COL_T+COL_UV+COL_RH)
      a_m(ib) = 0._wp
      if (lnewpl) then
        a_m(ib)% col% s% ps = &
          y * (gacc/R) * (cbg(ib)% col% s% ps / cbg(ib)% col% s% t2m)
        deallocate (y    )
        psx => cbg(ib)% col% s% ps
        psx = psx + a_m(ib)% col% s% ps
      endif
      if (dace% lpio) write(6,*) 'apply_B_mi_1d: pass 2'
      !----------------------
      ! second pass: t,rh,u,v
      ! allocate temporaries
      !----------------------
      allocate (t_int (ncol*nvc    ))
      allocate (e     (ncol*nvc    ))
      allocate (z     (nvc         ))
      allocate (levs  (ncol*(ke+ns)))
      allocate (y     (ncol*nvc    ))
      allocate (ph    (ke+1        ))
      allocate (lpf   (ke          ))
      !-----------------
      ! set up meta data
      !-----------------
      t_int (1       ) = OBS_H
      t_int (2:nvc:nm) = OBS_TV
      t_int (3:nvc:nm) = OBS_RH
      t_int (4:nvc:nm) = OBS_U
      t_int (5:nvc:nm) = OBS_V
      do is = 1, ncol
        c    =>      cbg(ib)% col(is)% c
        i1 = (is-1) *  nvc    + 1
        in =  is    *  nvc
        l1 = (is-1) * (ke+ns) + 1
        ln =  is    * (ke+ns)
        select case (cbg(ib)% vctype)
        case (VCT_P_ISO, VCT_P_HYB)
          do k = 1, ke + 1
            ph (k) = a_m(ib)% ak(k) + a_m(ib)% bk(k) * cbg(ib)% col(is)% s% ps
          end do
!NEC$ ivdep
          do k = 1, ke
            lpf (k) = log(0.5_wp * (ph (k) + ph (k+1)))
            z (nm*(k-1)+2:nm*k+1) = lpf(k)
          end do
        case (VCT_Z_HYB, VCT_Z_GEN)
!NEC$ ivdep
          do k = 1, ke
            lpf (k) = log (cbg(ib)% col(is)% p(k))
            z (nm*(k-1)+2:nm*k+1) = lpf(k)
          end do
        case default
          call finish ('apply_B_mi_1d','invalid vctype')
        end select
        z (1) = log (cbg(ib)% col(is)% s% ps)
        call set_vertinpc (levs(l1)% vi, z(1), covm% logp, covm% nwv)
        call zero (col)
        call init_spot (col,          &
                        1            ,&! spot index within the box
                        c            ,&! coordinates
                        t_int(1:nvc) ,&! observation type identifiers
                        z            ,&! vertical coordinates (ln p)
                    err=e(i1:in))      ! error (stdev)
        t_int(i1:in)        = t_int(1:nvc)
        levs(l1:ln)  %     is   = is
        levs(l1)     % i%  i    = i1 - 1
        levs(l1)     % i%  n    = 1
        levs(l1)     % vi% ezn  = col% cols(1)% ezn
        levs(l1)     % vi% exzn = col% cols(1)% exzn
        levs(l1)     %     z    = z(1)
        levs(l1+1:ln)%     z    = lpf
        do k=1,ke
          call set_vertinpc (levs(l1+k)% vi, lpf(k), covm% logp, covm% nwv)
          levs(l1+k) % i%  i    = i1 + (k-1)*nm
          levs(l1+k) % i%  n    = nm
          levs(l1+k) % vi% ezn  = col% cols(2+(k-1)*nm)% ezn
          levs(l1+k) % vi% exzn = col% cols(2+(k-1)*nm)% exzn
        end do
      end do
      !------------
      ! interpolate
      !------------
      call vertinp (y, ml, t_int, levs, e, y)
      !-------------
      ! store result
      !-------------
      do is = 1, ncol
        i1 = (is-1) *  nvc    + 1
        in =  is    *  nvc
        if (.not.lnewpl) &
          a_m(ib)% col(is)% s% ps = y (i1) * (gacc/R) * &
         (cbg(ib)% col(is)% s% ps / cbg(ib)% col(is)% s% t2m)
        a_m(ib)%   col(is)%    t  = y (i1+1:in:nm)
        a_m(ib)%   col(is)%    rh = y (i1+2:in:nm)
        a_m(ib)%   col(is)%    u  = y (i1+3:in:nm)
        a_m(ib)%   col(is)%    v  = y (i1+4:in:nm)
      end do
      !-----------------------
      ! deallocate temporaries
      !-----------------------
      deallocate (t_int )
      deallocate (e     )
      deallocate (z     )
      deallocate (levs  )
      deallocate (y     )
      deallocate (hic   )
      deallocate (ph,lpf)
      call destruct (ml)
    end do
    !-----------------------------
    ! deallocate temporary storage
    !-----------------------------
    call trace_mem_usage ('B_mi:total')
    call destruct (mr)
    call destruct (col)
    call destruct (ps)
    call delete_storage (x)
    call delete_storage (e_fi)
    call compress_cov
    call leave_function
  end subroutine apply_B_mi_1d
!-----------------------------------------------------------------------------
  subroutine horicov_oo (mr, ml, obs)
  type (t_mode)    ,intent(inout) :: mr (:) ! vertical modes  rhs
  type (t_mode)    ,intent(inout) :: ml (:) ! vertical modes  lhs
  type (t_obs_set) ,intent(in)    :: obs    ! observation meta data
  !-----------------------------------
  ! apply horizontal covariance matrix
  !-----------------------------------
!   target                 :: ml
    type(t_coord) ,pointer :: cl(:)
    integer                :: ibl
    integer                :: is, ic, ics
    type(t_spot)  ,pointer :: si
    !---------------------
    ! loop over boxes, lhs
    !---------------------
    do ibl = 1, size(ml)
      if (obs% o(ibl)% pe /= dace% pe) cycle
      !--------------------------
      ! prepare coordinates (lhs)
      !--------------------------
      if (obs% o(ibl)% n_spt == obs% o(ibl)% n_spot) then
        cl => obs% o(ibl)% spot(:)% col% c
      else
        allocate (cl(obs% o(ibl)% n_spt))
        ic = 0
        do is = 1, obs% o(ibl)% n_spot
          si  => obs% o(ibl)% spot(is)
          select case (si% int_type)
          case (ITY_ICOL)
            ic = ic + 1
            cl(ic) = si% col% c
          case (ITY_MCOLS)
            do ics = 1, size(si% imcol)
              ic = ic + 1
              cl(ic) = si% imcol(ics)% c
            end do
          case default
            call finish('horicov_oo','invalid int_type')
          end select
        end do
      endif
      !------------------------
      ! apply covariance matrix
      !------------------------
      call horicov (mr, ml(ibl), obs, cl)
      if (obs% o(ibl)% n_spt /= obs% o(ibl)% n_spot) deallocate (cl)
    end do
  end subroutine horicov_oo
!-----------------------------------------------------------------------------
  subroutine horicov (mr, ml, obs, col)
  type (t_mode)    ,intent(inout) :: mr (:) ! vertical modes  rhs
  type (t_mode)    ,intent(inout) :: ml     ! vertical modes  lhs
  type (t_obs_set) ,intent(in)    :: obs    ! observation meta data
  type(t_coord)    ,intent(in)    :: col(:) !
  !-----------------------------------
  ! apply horizontal covariance matrix
  !-----------------------------------
    target                 :: mr
    target                 :: col
    integer                :: ibr       ! box  indices  lhs,rhs
    integer                :: isl, isr  ! spot indices  lhs,rhs
    integer                :: is, ic, nc! report index, column index
    type(t_spot)  ,pointer :: si        ! reference to report meta data
    real(wp)               :: w, w1     ! interpolation weights
    integer                :: j, j1     ! latitude indices
    real(wp)  ,allocatable :: mul(:)
    real(wp)               :: mur
    type(t_coord) ,pointer :: cl, cr
    integer                :: month
    type (t_mode) ,pointer :: r
!   real(wp)               :: edem1
    real(wp)               :: r1nu      ! 1 - rnu
    real(wp)               :: sq1nu     ! sqrt(1-rnu)

    !----------------------------------
    ! horizontal correlation parameters
    !----------------------------------
    real(wp) :: chq     ! horizontal correlation for relative humidity
    real(wp) :: chh     ! horizontal correlation for height
    real(wp) :: clsvv   ! large scale contribution to wind-wind correlation
    real(wp) :: rlsvv   ! renormalisation due to large scale contribution
    real(wp) :: d1chh   ! first  derivative of chh
    real(wp) :: d1chx   ! first  derivative of chh, truncated for h-w corr.
    real(wp) :: d2chh   ! second derivative of chh
    real(wp) :: d1chho1 ! first  derivative of chh / 1.

    !----------------------------------------------------------------
    ! scalar product of longitudinal and transversal (l,t) directions
    ! with local (u,v) coordinate system
    !----------------------------------------------------------------
    real(wp) :: ab(3) ! vector between locations on the unit sphere
    real(wp) :: d     ! distance between locations on the unit sphere
    real(wp) :: abn   ! normalization factor
    real(wp) :: elui  ! longitidinal * u ,right hand side
    real(wp) :: elvi  ! longitidinal * v ,right hand side
    real(wp) :: etui  ! transversal  * u ,right hand side
    real(wp) :: etvi  ! transversal  * v ,right hand side
    real(wp) :: eluj  ! longitidinal * u ,left  hand side
    real(wp) :: elvj  ! longitidinal * v ,left  hand side
    real(wp) :: etuj  ! transversal  * u ,left  hand side
    real(wp) :: etvj  ! transversal  * v ,left  hand side

    real(wp) :: c_h_vt   ! correlation height      - transv.wind
    real(wp) :: c_vt_h   ! correlation transv.wind - height
    real(wp) :: ci_vt_vt ! correlation transv.wind - transv.wind
    real(wp) :: ci_vl_vl ! correlation longit.wind - longit.wind
    real(wp) :: xi_vt_vt ! correlation transv.wind - transv.wind
    real(wp) :: xi_vl_vl ! correlation longit.wind - longit.wind

    !------------------------------
    ! set and check some parameters
    !------------------------------
!   edem1 = 1._wp / (exp(1._wp)-1._wp)
    r1nu  = 1._wp - rnu
    sq1nu = sqrt(r1nu)
    month = imm (covm% time)
    if (L_h <= 0._wp) call finish('horicov','L_h should be a constant')
    if (L_q <= 0._wp) call finish('horicov','L_q should be a constant')
    !---------------------
    ! loop over boxes, lhs
    !---------------------
!   do ibl = 1, size(ml)
!     if (obs% o(ibl)% pe /= dace% pe) cycle
      !------------------------------------------
      ! calculate quantities depending on the lhs
      !------------------------------------------
      allocate (mul(size(col)))
      if (rmu >= -1._wp .and. rmu <= 1._wp) then
        mul (1:size(ml% lh)) = rmu
      else
        do isl = 1, size(ml% lh)
          cl => col(isl)
          w = min (30._wp, max (1._wp, (cl% dlat+87._wp)/6._wp+1))
          j = int(w)
          w = w - j
          w1 = 1._wp - w
          j1 = min (30, j+1)
          mul (isl) = w * nfcoun% mu (j1, month) + w1 * nfcoun% mu (j, month)
        end do
      end if
      !---------------------
      ! loop over boxes, rhs
      !---------------------
      do ibr = 1, size(mr)
        !---------------------
        ! loop over spots, rhs
        !---------------------
        is = 0
        ic = 0
        nc = 0
        do isr = 1, size(mr(ibr)% lh)
          ic = ic + 1
          if (ic > nc) then
            is = is + 1
            si => obs% o(ibr)% spot(is)
            ic = 1
            nc = 1
            if (associated (si% imcol)) nc = size (si% imcol)
          endif
          cr => si% col% c
          if (associated (si% imcol)) cr => si% imcol(ic)% c
          !------------------------------------------
          ! calculate quantities depending on the rhs
          !------------------------------------------
          if (mr(ibr)% luv(isr)) then
            if (rmu >= -1._wp .and. rmu <= 1._wp) then
              mur = rmu
            else
              w = min (30._wp, max (1._wp, (cr% dlat+87._wp)/6._wp+1))
              j = int(w)
              w = w - j
              w1 = 1._wp - w
              j1 = min (30, j+1)
              mur = w * nfcoun% mu(j1,month) + w1 * nfcoun% mu(j,month)
            endif
          endif
          !---------------------
          ! loop over spots, lhs
          !---------------------
          do isl = 1, size(col)
            cl => col(isl)
            !---------------
            ! check distance
            !---------------
            ab = cr% x - cl% x
            d  = sqrt(sum(ab**2))
            if (d >= x_cut*l_max) cycle
            !------------------------------
            ! derive horizontal covariances
            !------------------------------
            chh    = corr (d, L_h)
            call dcorr (d1chh, d1chho1, d2chh, d1chx)
            clsvv  = lsvv * chh
            rlsvv  = 1._wp / (1._wp + lsvv)
            if (chh == 0._wp .and. d1chh == 0._wp) cycle
            chq = corrh (d,L_q)
            !---------------------------------
            ! coefficients for wind directions
            !---------------------------------
            if (ml% luv(isl) .or. mr(ibr)% luv(isr)) then
              if (d==0._wp) ab = cl% du           ! arbitrary direction for d=0
              elui =   sum (ab * cl% du)          ! longit.vector in u,v system
              elvi =   sum (ab * cl% dv)          !
              abn  =   1._wp / sqrt (elui**2 + elvi**2) ! normalize
              elui =   abn * elui                 !
              elvi =   abn * elvi                 !
              eluj =   sum (ab * cr% du)          ! same at left hand side
              elvj =   sum (ab * cr% dv)          !  (changed sign)
              abn  =   1._wp / sqrt (eluj**2 + elvj**2)
              eluj = - (abn * eluj)               !
              elvj = - (abn * elvj)               !
              etui = - elvi                       ! rotate for tang.vector
              etvi =   elui                       !
              etuj = - elvj                       !
              etvj =   eluj
              if (mr(ibr)% luv(isr)) c_h_vt  = mur      * sq1nu * d1chx * rlsvv
              if (ml     % luv(isl)) c_vt_h  = mul(isl) * sq1nu * d1chx * rlsvv
              if (ml% luv(isl) .and. mr(ibr)% luv(isr)) then
                ci_vt_vt = (r1nu * d2chh   - clsvv) * rlsvv
                xi_vt_vt =  rnu  * d1chho1          * rlsvv
                ci_vl_vl = (r1nu * d1chho1 - clsvv) * rlsvv
                xi_vl_vl =  rnu  * d2chh            * rlsvv
              endif
            endif
            !---------
            ! multiply
            !---------
            r => mr(ibr)
            if (r% lh(isr)) then
              !------
              ! h - h
              !------
              if (ml% lh(isl)) then
                ml% h(:,isl) = ml% h(:,isl) + chh * r% h(:,isr)
              endif
              !--------
              ! u,v - h
              !--------
              if (ml% luv(isl)) then
                ml% u(:,isl) = ml% u(:,isl) - etui * c_vt_h * r% h(:,isr)
                ml% v(:,isl) = ml% v(:,isl) - etvi * c_vt_h * r% h(:,isr)
              endif
            endif
            if (r% lq(isr)) then
              !------
              ! q - q
              !------
              if (ml% lq(isl)) then
                ml% q(:,isl) = ml% q(:,isl) + chq * r% q(:,isr)
              endif
            endif
            if (r% luv(isr)) then
              !--------
              ! h - u,v
              !--------
              if (ml% lh(isl)) then
                ml% h(:,isl) = ml% h(:,isl)                     &
                  - c_h_vt * (etuj * r% u(:,isr) + etvj * r% v(:,isr))
              endif
              !----------
              ! u,v - u,v
              !----------
              if (ml% luv(isl)) then
                ml% u (:,isl) = ml% u (:,isl)   &
                  + (ci_vt_vt * etui * etuj + &
                     ci_vl_vl * elui * eluj ) * r% u (:,isr) &
                  + (ci_vt_vt * etui * etvj + &
                     ci_vl_vl * elui * elvj ) * r% v (:,isr)
                ml% v (:,isl) = ml% v (:,isl)   &
                  + (ci_vt_vt * etvi * etuj + &
                     ci_vl_vl * elvi * eluj ) * r% u (:,isr) &
                  + (ci_vt_vt * etvi * etvj + &
                     ci_vl_vl * elvi * elvj ) * r% v (:,isr)
                ml% ux(:,isl) = ml% ux (:,isl)   &
                  + (xi_vt_vt * etui * etuj + &
                     xi_vl_vl * elui * eluj ) * r% ux(:,isr) &
                  + (xi_vt_vt * etui * etvj + &
                     xi_vl_vl * elui * elvj ) * r% vx(:,isr)
                ml% vx(:,isl) = ml% vx (:,isl)   &
                  + (xi_vt_vt * etvi * etuj + &
                     xi_vl_vl * elvi * eluj ) * r% ux(:,isr) &
                  + (xi_vt_vt * etvi * etvj + &
                     xi_vl_vl * elvi * elvj ) * r% vx(:,isr)
              endif
            endif
          end do
        end do
      end do
      deallocate (mul)
!   end do
  end subroutine horicov
!-----------------------------------------------------------------------------
  subroutine int2vmode (mr, obs)
  type (t_mode)    ,intent(inout) :: mr (:) ! vertical modes       rhs
  type (t_obs_set) ,intent(in)    :: obs    ! observation meta data
  !----------------------------
  ! transform to vertical modes
  !----------------------------
    integer                :: ib                    ! box  index
    integer                :: is                    ! spot index
    integer                :: ic, ics               ! column index
    type (t_spot) ,pointer :: si                    ! spot index
    real(wp)               :: m(covm% nz, covm% nz) ! transformation matrix
    real(wp)               :: e(covm% nz)
    type (t_hic1) ,pointer :: hic                   ! hor. interpolation coefs
    !----------------
    ! loop over boxes
    !----------------
    do ib = 1,size(mr)
      if (obs% o(ib)% pe /= dace% pe) cycle
      !----------------
      ! loop over spots
      !----------------
      ic = 0
      do is = 1, obs% o(ib)% n_spot
        si => obs% o(ib)% spot(is)
        !------------------
        ! loop over columns
        !------------------
        ics = 0
        do
          ics = ics + 1
          select case (si% int_type)
          case (ITY_ICOL)
            if (ics > 1 ) exit
            hic => si% hic
          case (ITY_MCOLS)
            if (ics > size(si% imcol)) exit
            hic => si% imcol(ics)% hi
          case default
            call finish('int2vmode','invalid int_type')
          end select
          ic = ic + 1

          !------------------------------
          ! transform geopotential height
          !------------------------------
          if (mr(ib)% lh(ic)) then
            mr(ib)% lh  (ic) = .true.
            if (sparse) then
              if (.not. sqrtcor) then
                e = hic% w(1) * covm% eh(:, hic% ix  ) + &
                    hic% w(2) * covm% eh(:, hic% ix+1)
                mr(ib)% h (:,ic) = mr(ib)% h (:,ic) / e
              endif
              call wave_1d(mr(ib)% h (:,ic:ic),TR_ADJ,base)
              m = hic% w(1) * covm% sqcv1h (:,:, hic% ix  ) + &
                  hic% w(2) * covm% sqcv1h (:,:, hic% ix+1)
              mr(ib)% h (:,ic) = matmul (mr(ib)% h (:,ic), m)
            else
              m = hic% w(1) * covm% sqcvh (:,:, hic% ix  ) + &
                  hic% w(2) * covm% sqcvh (:,:, hic% ix+1)
              mr(ib)% h (:,ic) = matmul (mr(ib)% h (:,ic), m)
            endif
          endif
          !-------------------
          ! transform moisture
          !-------------------
          if (mr(ib)% lq(ic)) then
            mr(ib)% lq  (ic) = .true.
            if (sparse) then
              if (.not. sqrtcor) then
                e = hic% w(1) * covm% erh(:, hic% ix  ) + &
                    hic% w(2) * covm% erh(:, hic% ix+1)
                mr(ib)% q (:,ic) = mr(ib)% q (:,ic) / e
              endif
              call wave_1d( mr(ib)% q (:,ic:ic),TR_ADJ,base)
              m = hic% w(1) * covm% sqcv1q (:,:, hic% ix  ) + &
                  hic% w(2) * covm% sqcv1q (:,:, hic% ix+1)
              mr(ib)% q (:,ic) = matmul (mr(ib)% q (:,ic), m)
            else
              m = hic% w(1) * covm% sqcvq (:,:, hic% ix  ) + &
                  hic% w(2) * covm% sqcvq (:,:, hic% ix+1)
              mr(ib)% q (:,ic) = matmul (mr(ib)% q (:,ic), m)
            endif
          endif
          !--------------------------
          ! transform wind components
          !--------------------------
          if (mr(ib)% luv(ic)) then
            mr(ib)% luv (ic) = .true.
            if (sparse) then
              if (.not. sqrtcor) then
                e = hic% w(1) * covm% exi(:, hic% ix  ) + &
                    hic% w(2) * covm% exi(:, hic% ix+1)
                mr(ib)% ux (:,ic) = mr(ib)% u (:,ic) / e
                mr(ib)% vx (:,ic) = mr(ib)% v (:,ic) / e
                e = hic% w(1) * covm% epsi(:, hic% ix  ) + &
                    hic% w(2) * covm% epsi(:, hic% ix+1)
                mr(ib)% u (:,ic) = mr(ib)% u (:,ic) / e
                mr(ib)% v (:,ic) = mr(ib)% v (:,ic) / e
                call wave_1d( mr(ib)% ux (:,ic:ic),TR_ADJ,base)
                call wave_1d( mr(ib)% vx (:,ic:ic),TR_ADJ,base)
                call wave_1d( mr(ib)% u (:,ic:ic),TR_ADJ,base)
                call wave_1d( mr(ib)% v (:,ic:ic),TR_ADJ,base)
              else
                call wave_1d( mr(ib)% u (:,ic:ic),TR_ADJ,base)
                call wave_1d( mr(ib)% v (:,ic:ic),TR_ADJ,base)
                mr(ib)% ux (:,ic) = mr(ib)% u (:,ic)
                mr(ib)% vx (:,ic) = mr(ib)% v (:,ic)
              endif
              m = hic% w(1) * covm% sqcv1v (:,:, hic% ix  ) + &
                  hic% w(2) * covm% sqcv1v (:,:, hic% ix+1)
              mr(ib)% ux(:,ic) = matmul (mr(ib)% ux (:,ic), m)
              mr(ib)% vx(:,ic) = matmul (mr(ib)% vx (:,ic), m)
              m = hic% w(1) * covm% sqcv1s (:,:, hic% ix  ) + &
                  hic% w(2) * covm% sqcv1s (:,:, hic% ix+1)
              mr(ib)% u (:,ic) = matmul (mr(ib)% u (:,ic), m)
              mr(ib)% v (:,ic) = matmul (mr(ib)% v (:,ic), m)
            else
              m = hic% w(1) * covm% sqcvv (:,:, hic% ix  ) + &
                  hic% w(2) * covm% sqcvv (:,:, hic% ix+1)
              mr(ib)% ux(:,ic) = matmul (mr(ib)% u (:,ic), m)
              mr(ib)% vx(:,ic) = matmul (mr(ib)% v (:,ic), m)
              m = hic% w(1) * covm% sqcvs (:,:, hic% ix  ) + &
                  hic% w(2) * covm% sqcvs (:,:, hic% ix+1)
              mr(ib)% u (:,ic) = matmul (mr(ib)% u (:,ic), m)
              mr(ib)% v (:,ic) = matmul (mr(ib)% v (:,ic), m)
            endif
          endif
        end do
      end do
    end do
  end subroutine int2vmode
!-----------------------------------------------------------------------------
  subroutine vmode2int_o (ml, obs)
  type (t_mode)    ,intent(inout) :: ml (:) ! vertical modes       lhs
  type (t_obs_set) ,intent(in)    :: obs    ! observation meta data
  !------------------------------
  ! transform from vertical modes
  !------------------------------
    integer                :: ib                    ! box  index
    integer                :: is                    ! spot index
    integer                :: ic, ics               ! column index
    type (t_spot) ,pointer :: si                    ! spot index
    type (t_hic1) ,pointer :: hic (:)               ! hor. interpolation coefs
    !----------------
    ! loop over boxes
    !----------------
    do ib = 1,size(ml)
      if (obs% o(ib)% pe /= dace% pe) cycle
      !--------------------
      ! prepare coordinates
      ! call vmode2int
      !--------------------
      if (obs% o(ib)% n_spt == obs% o(ib)% n_spot) then
        call vmode2int (ml(ib), obs% o(ib)% spot(:)% hic)
      else
        allocate (hic(obs% o(ib)% n_spt))
        ic = 0
        do is = 1, obs% o(ib)% n_spot
          si  => obs% o(ib)% spot(is)
          select case (si% int_type)
          case (ITY_ICOL)
            ic = ic + 1
            hic(ic) = si% hic
          case (ITY_MCOLS)
            do ics = 1, size(si% imcol)
              ic = ic + 1
              hic(ic) = si% imcol(ics)% hi
            end do
          case default
            call finish('vmode2int_o','invalid int_type')
          end select
        end do
        call vmode2int (ml(ib), hic)
        deallocate (hic)
      endif
    end do
  end subroutine vmode2int_o
!-----------------------------------------------------------------------------
  subroutine vmode2int (ml, hic)
  type (t_mode)    ,intent(inout) :: ml     ! vertical modes       lhs
  type (t_hic1)    ,intent(in)    :: hic(:) ! horizontal interp. coefs

    integer  :: is                    ! spot index
    real(wp) :: m(covm% nz, covm% nz) ! transformation matrix
    real(wp) :: e(covm% nz)
    !----------------
    ! loop over spots
    !----------------
    do is = 1, size (hic)
      !------------------------------
      ! transform geopotential height
      !------------------------------
      if (ml% lh(is)) then
        ml% lh  (is) = .true.
        if (sparse) then
          m = hic(is)% w(1) * covm% sqcv1h (:,:, hic(is)% ix  ) + &
              hic(is)% w(2) * covm% sqcv1h (:,:, hic(is)% ix+1)
          ml% h (:,is) = matmul (m, ml% h (:,is))
          call wave_1d( ml% h (:,is:is),TR_SYN,base)
          if (.not. sqrtcor) then
            e = hic(is)% w(1) * covm% eh(:, hic(is)% ix  ) + &
                hic(is)% w(2) * covm% eh(:, hic(is)% ix+1)
            ml% h (:,is) = ml% h (:,is) / e
          endif
        else
          m = hic(is)% w(1) * covm% sqcvh (:,:, hic(is)% ix  ) + &
              hic(is)% w(2) * covm% sqcvh (:,:, hic(is)% ix+1)
          ml% h (:,is) = matmul (m, ml% h (:,is))
        endif
      endif
      !-------------------
      ! transform moisture
      !-------------------
      if (ml% lq(is)) then
        ml% lq  (is) = .true.
        if (sparse) then
          m = hic(is)% w(1) * covm% sqcv1q (:,:, hic(is)% ix  ) + &
              hic(is)% w(2) * covm% sqcv1q (:,:, hic(is)% ix+1)
          ml% q (:,is) = matmul (m, ml% q (:,is))
          call wave_1d( ml% q (:,is:is),TR_SYN,base)
          if (.not. sqrtcor) then
            e = hic(is)% w(1) * covm% erh(:, hic(is)% ix  ) + &
                hic(is)% w(2) * covm% erh(:, hic(is)% ix+1)
            ml% q (:,is) = ml% q (:,is) / e
          endif
        else
          m = hic(is)% w(1) * covm% sqcvq (:,:, hic(is)% ix  ) + &
              hic(is)% w(2) * covm% sqcvq (:,:, hic(is)% ix+1)
          ml% q (:,is) = matmul (m, ml% q (:,is))
        endif
      endif
      !--------------------------
      ! transform wind components
      !--------------------------
      if (ml% luv(is)) then
        ml% luv (is) = .true.
        if (sparse) then
          m = hic(is)% w(1) * covm% sqcv1v (:,:, hic(is)% ix  ) + &
              hic(is)% w(2) * covm% sqcv1v (:,:, hic(is)% ix+1)
          ml% ux(:,is) = matmul (m, ml% ux (:,is))
          ml% vx(:,is) = matmul (m, ml% vx (:,is))
          m = hic(is)% w(1) * covm% sqcv1s (:,:, hic(is)% ix  ) + &
              hic(is)% w(2) * covm% sqcv1s (:,:, hic(is)% ix+1)
          ml% u (:,is) = matmul (m, ml% u (:,is))
          ml% v (:,is) = matmul (m, ml% v (:,is))


          if (.not. sqrtcor) then
            call wave_1d( ml% ux (:,is:is),TR_SYN,base)
            call wave_1d( ml% vx (:,is:is),TR_SYN,base)
            call wave_1d( ml% u  (:,is:is),TR_SYN,base)
            call wave_1d( ml% v  (:,is:is),TR_SYN,base)
            e = hic(is)% w(1) * covm% epsi(:, hic(is)% ix  ) + &
                hic(is)% w(2) * covm% epsi(:, hic(is)% ix+1)
            ml% u (:,is) = ml% u (:,is) / e
            ml% v (:,is) = ml% v (:,is) / e
            e = hic(is)% w(1) * covm% exi(:, hic(is)% ix  ) + &
                hic(is)% w(2) * covm% exi(:, hic(is)% ix+1)
            ml% u (:,is) = ml% u (:,is) + ml% ux (:,is) / e
            ml% v (:,is) = ml% v (:,is) + ml% vx (:,is) / e
          else
            ml% u (:,is) = ml% u (:,is) + ml% ux (:,is)
            ml% v (:,is) = ml% v (:,is) + ml% vx (:,is)
            call wave_1d( ml% u (:,is:is),TR_SYN,base)
            call wave_1d( ml% v (:,is:is),TR_SYN,base)
          endif
        else
          m = hic(is)% w(1) * covm% sqcvv (:,:, hic(is)% ix  ) + &
              hic(is)% w(2) * covm% sqcvv (:,:, hic(is)% ix+1)
          ml% ux(:,is) = matmul (m, ml% ux (:,is))
          ml% vx(:,is) = matmul (m, ml% vx (:,is))
          m = hic(is)% w(1) * covm% sqcvs (:,:, hic(is)% ix  ) + &
              hic(is)% w(2) * covm% sqcvs (:,:, hic(is)% ix+1)
          ml% u (:,is) = ml% ux(:,is) + matmul (m, ml% u (:,is))
          ml% v (:,is) = ml% vx(:,is) + matmul (m, ml% v (:,is))
        endif
      endif
    end do
  end subroutine vmode2int
!-----------------------------------------------------------------------------
  subroutine vertinpa (x, mr, obs, e_fi)
  type (t_vector)  ,intent(in)    :: x      ! input
  type (t_mode)    ,intent(inout) :: mr (:) ! vertical modes      rhs
  type (t_obs_set) ,intent(in)    :: obs    ! observation meta data
  type (t_vector)  ,intent(in)    :: e_fi ! background error
  !-------------------------------
  ! adjoint vertical interpolation
  !-------------------------------
    integer :: ib                 ! box    index
    integer :: is                 ! spot   index
    integer :: ic                 ! column index
    integer :: l, ll              ! level  index
    integer :: i                  ! observation index
    integer :: k                  ! level  index
    integer :: nl                 ! no.levels/column
    integer :: ml                 ! no.levels/column
    type (t_spot) ,pointer :: si  ! spot pointer
    !----------------
    ! loop over boxes
    !----------------
    do ib = 1,size(mr)

mr(ib)%  h = 0._wp
mr(ib)%  u = 0._wp
mr(ib)%  v = 0._wp
mr(ib)%  ux = 0._wp
mr(ib)%  vx = 0._wp
mr(ib)%  q = 0._wp
mr(ib)% lh  = .false.
mr(ib)% luv = .false.
mr(ib)% lq  = .false.

      if (obs% o(ib)% pe /= dace% pe) cycle
      !----------------
      ! loop over spots
      !----------------
      ic = 0
      do is = 1, obs% o(ib)% n_spot
        si => obs% o(ib)% spot(is)
        nl = 0
        ll = 0
        ml = si% l% n
        if (associated(si% imcol)) ml = ml / size (si% imcol)
        !-----------------
        ! loop over levels
        !-----------------
        do l = si% l% i + 1, si% l% i + si% l% n
          ll = ll + 1
          if (ll>nl) then
            nl = nl + ml
            ic = ic + 1
          endif
          !-------------------------------------
          ! loop over observations at this level
          !-------------------------------------
          do i = obs% o(ib)% levs(l)% i% i + 1, &
                 obs% o(ib)% levs(l)% i% i + obs% o(ib)% levs(l)% i% n
            select case (obs% o(ib)% t_int(i))
            case (OBS_H, OBS_HS)
              mr(ib)% lh(ic) = .true.
              do k = 1, obs% o(ib)% levs(l)% vi% g% n
                mr(ib)  % h(obs% o(ib)% levs(l)% vi% g% i+k-1 ,ic) = &
                  mr(ib)% h(obs% o(ib)% levs(l)% vi% g% i+k-1 ,ic) + &
                            obs% o(ib)% levs(l)% vi% wh(k)         * &
                              x% s(ib)% x   (i)                    * &
                           e_fi% s(ib)% x   (i)
              end do
            case (OBS_TV)
              mr(ib)% lh(ic) = .true.
              do k = 1, obs% o(ib)% levs(l)% vi% g% n

!               mr(ib)  % h(obs% o(ib)% levs(l)% vi% ix+k-1 ,ic)          = &
!                 mr(ib)% h(obs% o(ib)% levs(l)% vi% ix+k-1 ,ic)          + &
!                           obs% o(ib)% levs(l)% vi% wt(k)                * &
!                             x% s(ib)% x   (i)                           * &
!                                  (gacc / R * si% emod) * (si% hic% w(1) * &
!                 covm% eh (obs% o(ib)% levs(l)% vi% ix+k-1,si% hic% ix)  + &
!                                                           si% hic% w(2) * &
!                 covm% eh (obs% o(ib)% levs(l)% vi% ix+k-1,si% hic% ix+1))

                mr(ib)  % h(obs% o(ib)% levs(l)% vi% g% i+k-1 ,ic) = &
                  mr(ib)% h(obs% o(ib)% levs(l)% vi% g% i+k-1 ,ic) - &
                              x% s(ib)% x   (i)                    * &
                           e_fi% s(ib)% x   (i)                    * &
                         (  obs% o(ib)% levs(l)% vi% wt(k)         * &
                            obs% o(ib)% levs(l)% vi% exzn          + &
                            obs% o(ib)% levs(l)% vi% wh(k)         * &
                            obs% o(ib)% levs(l)% vi% ezn             )

              end do
            case (OBS_RH)
              mr(ib)% lq(ic) = .true.
              do k = 1, obs% o(ib)% levs(l)% vi% g% n
                mr(ib)  % q(obs% o(ib)% levs(l)% vi% g% i+k-1 ,ic) = &
                  mr(ib)% q(obs% o(ib)% levs(l)% vi% g% i+k-1 ,ic) + &
                            obs% o(ib)% levs(l)% vi% wh(k)         * &
                              x% s(ib)% x   (i)                    * &
                           e_fi% s(ib)% x   (i)
              end do
            case (OBS_U)
              mr(ib)% luv(ic) = .true.
              do k = 1, obs% o(ib)% levs(l)% vi% g% n
                mr(ib)  % u(obs% o(ib)% levs(l)% vi% g% i+k-1 ,ic) = &
                  mr(ib)% u(obs% o(ib)% levs(l)% vi% g% i+k-1 ,ic) + &
                            obs% o(ib)% levs(l)% vi% wh(k)         * &
                              x% s(ib)% x   (i)                    * &
                           e_fi% s(ib)% x   (i)
              end do
            case (OBS_V)
              mr(ib)% luv(ic) = .true.
              do k = 1, obs% o(ib)% levs(l)% vi% g% n
                mr(ib)  % v(obs% o(ib)% levs(l)% vi% g% i+k-1 ,ic) = &
                  mr(ib)% v(obs% o(ib)% levs(l)% vi% g% i+k-1 ,ic) + &
                            obs% o(ib)% levs(l)% vi% wh(k)         * &
                              x% s(ib)% x   (i)                    * &
                           e_fi% s(ib)% x   (i)
              end do
            case (OBS_DUM)
            case default
              write(0,*)'vertinpa: invalid interpolation type:',&
                        obs% o(ib)% t_int(i)
              write(0,*)'vertinpa: box,spot,l,i=',ib,is,l,i
              call finish ('vertinpa','invalid interpolation type')
            end select
          end do
        end do
      end do
    end do
  end subroutine vertinpa
!-----------------------------------------------------------------------------
  subroutine vertinp_o (y, ml, obs, e_fi, x)
  type (t_vector)  ,intent(inout) :: y      ! output lhs
  type (t_mode)    ,intent(in)    :: ml (:) ! vertical modes      rhs
  type (t_obs_set) ,intent(in)    :: obs    ! observation meta data
  type (t_vector)  ,intent(in)    :: e_fi   ! background error
  type (t_vector)  ,intent(in)    :: x      ! input rhs (shortcut)
  !-----------------------
  ! vertical interpolation
  !-----------------------
    integer                :: ib  ! box  index
    y=0._wp
    !----------------
    ! loop over boxes
    !----------------
    do ib = 1,size(ml)
      if (obs% o(ib)% pe /= dace% pe) cycle
      call vertinp (y%    s(ib)% x,        ml(ib),          &
                    obs%  o(ib)% t_int, obs%o(ib)% levs(:), &
                    e_fi% s(ib)% x,       x%s(ib)% x        )
    end do
  end subroutine vertinp_o
!-----------------------------------------------------------------------------
  subroutine vertinp (y, ml, t_int, levs, e_fi, x)
  real(wp)       ,intent(inout) :: y    (:)  ! output lhs
  type (t_mode)  ,intent(in)    :: ml        ! vertical modes      rhs
  integer        ,intent(in)    :: t_int(:)  ! observation type
  type (t_lev)   ,intent(in)    :: levs (:)  ! level data type
  real(wp)       ,intent(in)    :: e_fi (:)  ! background error
  real(wp)       ,intent(in)    :: x    (:)  ! input rhs (shortcut)
  !-----------------------
  ! vertical interpolation
  !-----------------------
!   integer                :: is  ! spot index
    integer                :: l   ! level index
    integer                :: i   ! observation index
    integer                :: k   ! level index
!   type (t_spot) ,pointer :: si  ! spot pointer

!   y=0._wp
!   !----------------
!   ! loop over boxes
!   !----------------
!   do ib = 1,size(ml)
!     if (obs% o(ib)% pe /= dace% pe) cycle
      !----------------
      ! loop over spots
      !----------------
!     do is = 1, size(ix)
!       si => obs% spot(is)
        !-----------------
        ! loop over levels
        !-----------------
!       do l = si% l% i + 1, si% l% i + si% l% n
        do l = 1, size(levs)
          !-------------------------------------
          ! loop over observations at this level
          !-------------------------------------
          do i = levs(l)% i% i + 1, &
                 levs(l)% i% i + levs(l)% i% n
            y (i) = 0
            select case (t_int(i))
            case (OBS_H, OBS_HS)
!             ml(ib)% lh(is) = .true.
              do k = 1, levs(l)% vi% g% n

                y (i) = y (i)                             + &
                ml% h(levs(l)% vi% g% i+k-1 ,levs(l)% is) * &
                      levs(l)% vi% wh(k)                  * &
                      e_fi(i)

!               y (i) = y (i)              + &
!               ml% h(levs(l)% vi% ix+k-1 ,levs(l)% is) * &
!                         levs(l)% vi% wh(k)       * &
!                                             si% emod * (si% hic% w(1) * &
!               covm% eh (levs(l)% vi% ix+k-1,si% hic% ix)  + &
!                                                         si% hic% w(2) * &
!               covm% eh (levs(l)% vi% ix+k-1,si% hic% ix+1))

!write(0,*)dace% pe,'vertinp: ib,levs(l)% is,l,i,k,  e_fi,eh,eh,w1,w2,emod',&
!ib,levs(l)% is,l,i,k,&
!e_fi(i), &
!covm% eh (levs(l)% vi% ix+k-1,si% hic% ix),&
!covm% eh (levs(l)% vi% ix+k-1,si% hic% ix+1),&
!si% hic% w(1),si% hic% w(2),si% emod,&
!covm% eh (:,si% hic% ix)

!           c% eh = c% eh                                                     &
!             + spt%h%w(1) * c%vi% wh(k2) * covm% eh (c%vi% ix+k2-1,spt%h%ix )&
!             + spt%h%w(2) * c%vi% wh(k2) * covm% eh (c%vi% ix+k2-1,spt%h%ix+1)


              end do
            case (OBS_TV)
!             ml% lt(levs(l)% is) = .true.
              do k = 1, levs(l)% vi% g% n

!               y (i) = y (i)                       + &
!               ml% h(levs(l)% vi% ix+k-1 ,levs(l)% is)          * &
!                         levs(l)% vi% wt(k)                * &
!                                (gacc / R * si% emod) * (si% hic% w(1) * &
!               covm% eh (levs(l)% vi% ix+k-1,si% hic% ix)  + &
!                                                         si% hic% w(2) * &
!               covm% eh (levs(l)% vi% ix+k-1,si% hic% ix+1))

                y (i) = y (i)                                   - &
                               e_fi (i)                         * &
                      ml% h(levs(l)% vi% g% i+k-1 ,levs(l)% is) * &
                             (  levs(l)% vi% wt(k)              * &
                                levs(l)% vi% exzn               + &
                                levs(l)% vi% wh(k)              * &
                                levs(l)% vi% ezn                  )

              end do
            case (OBS_RH)
!             ml% lq(levs(l)% is) = .true.
              do k = 1, levs(l)% vi% g% n
                y (i) = y (i)                             + &
                ml% q(levs(l)% vi% g% i+k-1 ,levs(l)% is) * &
                      levs(l)% vi% wh(k)                  * &
                      e_fi(i)
              end do
            case (OBS_U)
!             ml% lu(levs(l)% is) = .true.
              do k = 1, levs(l)% vi% g% n
                y (i) = y (i)                             + &
                ml% u(levs(l)% vi% g% i+k-1 ,levs(l)% is) * &
                      levs(l)% vi% wh(k)                  * &
                      e_fi(i)
              end do
            case (OBS_V)
!             ml% lv(levs(l)% is) = .true.
              do k = 1, levs(l)% vi% g% n
                y (i) = y (i)                             + &
                ml% v(levs(l)% vi% g% i+k-1 ,levs(l)% is) * &
                      levs(l)% vi% wh(k)                  * &
                      e_fi(i)

              end do
            case (OBS_DUM)
              !---------------------------------
              ! shortcut for dummy sink variable
              !---------------------------------
              y (i) = x (i) * e_fi (i) ** 2
            case default
              write(0,*)'vertinp: invalid interpolation type:',t_int(i)
              call finish ('vertinp','invalid interpolation type')
            end select
          end do
        end do
!     end do
!   end do
  end subroutine vertinp
!=============================================================================
end module mo_bg_err_op
