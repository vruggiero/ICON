!
!+ Representation of covariance matrices by a sequence of operators.
!
MODULE mo_t_bg_err_op
!
! Description:
!   Representation of covariance matrices by a sequence of operators.
!   Covariances are represented by sparse matrices in wavelet transformed
!   space.  Covariances are fitted by the NMC method.
!
!   This module holds the low level routines and derived type definitions.
!   High level routines are put in module mo_bg_err_op and mo_bg_err_2d.
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
! V1_5         2009/05/25 Harald Anlauf
!  set_vertinpc: optimizations recommended by NEC
! V1_8         2009/12/09 Harald Anlauf
!  set_vertinpc: improve optimization for SX-9 and cleanup
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Harald Anlauf
!  cov_from_nmc: print description attribute of vertical covariances
! V1_19        2012-04-16 Harald Anlauf
!  compress_cov/uncompress_cov: bugfix for nprocs>ny
! V1_28        2014/02/26 Andreas Rhodin
!  preparations for VarEnKF
! V1_31        2014-08-21 Andreas Rhodin
!  modified rttov specific vertical interpolation (nwv_rad=3 in /observations/)
!  temporal correlations in random NMC B representation
! V1_45        2015-12-15 Harald Anlauf
!  Change compress_cov/uncompress_cov for higher efficiency
! V1_48        2016-10-06 Robin Faulwetter
!  Implemented RTTOV12
! V1_51        2017-02-24 Harald Anlauf
!  tweak for Cray using OpenMP; use the 'contiguous' attribute
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Andreas Rhodin  DWD  2006-2008  original version
! Harald Anlauf   DWD  2007       support for horizontal wavelet correlations
! Harald Anlauf   DWD  2008-2009  optimizations for SX8,SX9
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
!=============================================================================
#include "tr15581.incf"
  !=============
  ! modules used
  !=============
  use mo_kind,       only: wp, i8           ! working precision kind parameter
  use mo_system,     only: flush            !
  use mo_mpi_dace,   only: dace,           &! MPI group info
                           p_barrier,      &! MPI barrier
                           p_bcast,        &! generic MPI bcast routine
                           p_ibcast,       &! generic non-blocking MPI bcast
                           p_waitall,      &! wait for MPI requests to complete
                           p_real_wp        ! MPI type for real(wp)
  use mo_exception,  only: finish,         &! abort in case of error
                           message          ! write message to stderr
  use mo_time,       only: t_time,         &! date&time data type
                           imm              ! derive month
  use mo_grads,      only: t_ctl,          &! GRADS .ctl meta data type
                           read_ctl,       &! read GRADS meta data
                           read_var,       &! read field from GRADS
                           destruct         ! destruct GRADS meta data
  use mo_dace_string,only: split,          &! Split text string into array
                           char2,          &! 2 digit integer to character
                           toupper          ! Convert to uppercase
  use mo_physics,    only: gacc             ! gravity acceleration
  !-----------------------
  ! observation processing
  !-----------------------
  use mo_t_obs,      only: t_vic,          &! vert.int.coeff.data type
                           t_hic1,         &! hor. int.coeff.data type
                           mwv,            &! max.no.coeff for vert.interpol.
                           ndv_rad,        &! radiance vertical interpol. flag
                           vint_rttov_b,   &! interpolation mode for radiances
                           INT_H,          &! geopotential height id
                           INT_RH,         &! relative humidity   id
                           INT_CHI,        &! velocity potential  id
                           INT_PSI          ! stream   function   id
  !-------------------
  ! numerical routines
  !-------------------
  use mo_random,     only: random_gauss,        &! Gaussian random number dist.
                           random_state_t,      &! random generator state type
                           destruct              ! destruct random_state_t
  use mo_random_seed,only: construct_stream      ! New random generator
  use mo_cov_wavelet,only: t_cov_wavelet,       &! Hor. wavelet cov. matrix
                           t_cov_meta,          &! Covariance matrix metadata
                           construct,           &! component initialization
                           destruct,            &! component deallocation
                           p_bcast,             &! overloaded MPI bcast routine
                           p_send,              &! overloaded MPI send  routine
                           p_recv,              &! overloaded MPI recv  routine
                           read_cov_netcdf,     &! Load cov.matrix from NetCDF
                           read_cov_netcdf_meta,&! Read cov.matrix metadata
                           read_profile_netcdf, &! Read meridional profile
                           read_profile_list_netcdf
  use mo_1dmra,      only: wv_name               ! Wavelet mnemonics
  use mo_transform,  only: gauaw            ! calculate Gaussian latitudes
  use mo_dec_matrix, only: t_matrix_block, &! Matrix blocks
!                          srep,           &! Mnemonic for storage repr.
                           destruct,       &! component deallocation
                           CSR,            &! compressed sparse rows
                           CSC,            &! compressed sparse columns
                           JAD,            &! jagged diagonal representation
                           SYMMETRIC,      &! symmetric matrix
                           convert,        &! Convert between representations
                           set_perm         ! Set CSR row permutation
  use mo_letkf_util, only: t_rndm_tcor,    &! temporal correlation meta data
                           construct        ! set up t_rndm_tcor
  !---------------------
  ! netCDF f90 interface
  !---------------------
  use netcdf,        only: nf90_open,             &! open   NetCDF file
                           nf90_close,            &! close  NetCDF file
                           nf90_Inquire_Dimension,&!
                           nf90_inq_dimid,        &!
                           nf90_inq_varid,        &!
                           nf90_get_var,          &!
                           nf90_get_att,          &!
                           NF90_GLOBAL,           &!
                           NF90_NOERR,            &!
                           NF90_NOWRITE            ! NetCDF read flag
  implicit none
  !================
  ! public entities
  !================
  private
  public :: cov_from_nmc    ! read NMC fitted covariances from GRADS file
  public :: t_cov           ! derived data type to hold covariance model
  public :: t_hc2_part      ! partitioning of hor. covariance matrices
  public :: covm            ! module variable to hold covariance model
  public :: construct       ! construct variables of type t_cov
  public :: destruct        ! destruct  variables of type t_cov
  public :: clean_bg_err_op ! clean up module variables
  public :: set_vertinpc    ! set vertical interpolation coefficients
  public :: set_horinpc     ! set horizont.interpolation coefficients
  public :: cov_mem         ! report memory allocated
  public :: compress_cov    ! store covm only once
  public :: uncompress_cov  ! store covm on each PE
  public :: vars            ! mapping: index in covm% hc2 vs. INT_H .. INT_PSI
  public :: read_horizontal_covariances ! Read horizontal wavelet covariances
  public :: random_gauss    ! Gaussian random number distribution

  !============================================================================
  !===================
  ! module parameters
  !===================
  integer, parameter :: nvar       = 4
  integer, parameter :: vars(nvar) = (/ INT_H, INT_RH, INT_PSI, INT_CHI /)

  !===================
  ! derived data types
  !===================

  !--------------------------------------------------
  ! partitioning of 2d horizontal covariance matrices
  !--------------------------------------------------
  type t_hc2_part
     integer :: pe              ! processor id
     integer :: k               ! level index
     integer :: l               ! index in local array
     integer :: id              ! variable id: INT_H, _RH, _PSI, _CHI
     integer :: i               ! variable index (1..nvar)
  end type t_hc2_part

  !----------------------------
  ! covariance matrix meta data
  !----------------------------
  type t_cov
    !=========================================================
    ! Preliminary implementation of separable covariance model
    !=========================================================
    !-------------------
    ! general parameters
    !-------------------
    integer           :: valid = 0     ! 1:content is valid  2:use covar.model
    integer           :: lclim = 0     ! 1:clim.err.  valid  2:use clim.err.
    integer           :: nz            ! number of vertical levels
    integer           :: ny            ! number of latitudinal  gridpoints
    integer           :: nx            ! number of longitudinal gridpoints
    real(wp)          :: pbot  =-9._wp ! bottom (surface) pressure (Pa)
    type(t_time)      :: time          ! time of analysis
    real(wp) ,pointer :: p     (:)     ! levels    (Pa)
    real(wp) ,pointer :: logp  (:)     ! ln(p)
    real(wp) ,pointer :: dlat  (:)     ! latitudes (degree)
    !------------------------------
    ! vertical correlation matrices
    !------------------------------
    real(wp) _POINTER :: sqcvh (:,:,:) ! sqrt of vert. cov. height  (np,np,ny)
    real(wp) _POINTER :: sqcvs (:,:,:) ! sqrt of vert. cov. streamf (np,np,ny)
    real(wp) _POINTER :: sqcvv (:,:,:) ! sqrt of vert. cov. vel.pot (np,np,ny)
    real(wp) _POINTER :: sqcvq (:,:,:) ! sqrt of vert. cov. humidity(np,np,ny)
    real(wp) _POINTER :: sqcvt (:,:,:) ! sqrt of vert. cov. temper. (test)
    !----------------
    ! temporary stuff
    !----------------
    logical           :: cmpr  =.false.! compress 3-d arrays
    real(wp) _POINTER :: sqcv1h(:,:,:) ! sqrt of vert. cov. height  (np,np,ny)
    real(wp) _POINTER :: sqcv1s(:,:,:) ! sqrt of vert. cov. streamf (np,np,ny)
    real(wp) _POINTER :: sqcv1v(:,:,:) ! sqrt of vert. cov. vel.pot (np,np,ny)
    real(wp) _POINTER :: sqcv1q(:,:,:) ! sqrt of vert. cov. humidity(np,np,ny)
    real(wp) _POINTER :: sqcv1t(:,:,:) ! sqrt of vert. cov. temper. (test)
    !-----------------------------
    ! errors (standard deviations)
    !-----------------------------
    real(wp) ,pointer :: eh      (:,:) ! height                   error (stdev)
    real(wp) ,pointer :: et      (:,:) ! temperature              error (stdev)
    real(wp) ,pointer :: erh     (:,:) ! rel.humidity             error (stdev)
    real(wp) ,pointer :: ev      (:,:) ! wind                     error (stdev)
    real(wp) ,pointer :: epsi    (:,:) ! streamfunction           error (stdev)
    real(wp) ,pointer :: exi     (:,:) ! velocity.pot.            error (stdev)
    real(wp) ,pointer :: eh_cl   (:,:) ! climatol. height         error (stdev)
    real(wp) ,pointer :: et_cl   (:,:) ! climatol. temperature    error (stdev)
    real(wp) ,pointer :: erh_cl  (:,:) ! climatol. rel.humidity   error (stdev)
    real(wp) ,pointer :: ev_cl   (:,:) ! climatol. wind           error (stdev)
    real(wp) ,pointer :: epsi_cl (:,:) ! climatol. streamfunction error (stdev)
    real(wp) ,pointer :: exi_cl  (:,:) ! climatol. velocity.pot.  error (stdev)
    !-----------------------------------------------
    ! parameters used by the interpolation operators
    !-----------------------------------------------
    integer           :: nwv           ! no.points used for vert. interpolation
    integer           :: nwh           ! no.points used for horz. interpolation
    !--------------------------------------------------------------------
    ! Wavelet representation of the horizontal (auto-)covariance matrices
    !--------------------------------------------------------------------
    integer                      :: hc2_totalsize   ! Total size = nvar*levels
    integer                      :: hc2_part_size   ! Size of local arrays
    integer                      :: gridtype        ! WMO grid type
    type(t_hc2_part)    ,pointer :: hc2_part  (:)   ! Partitioning table
    type(t_cov_wavelet)          :: hc2(nvar)       ! Covariance matrices
    !---------------------------------------------------------------
    ! ensemble representation of B matrix with temporal correlations
    !---------------------------------------------------------------
    type(t_rndm_tcor)            :: ct              ! time correlation data
    type(random_state_t),pointer :: rndm_state(:,:) ! random no.generator state
    !-----------------------------------
    ! multivariate covariance parameters
    !-----------------------------------
    real(wp)         ,pointer :: c_h_psi  (:)   ! correlation geop.-streamf.
    real(wp)         ,pointer :: c1_h_psi (:)   ! 1-abs(c1_h_psi)
    real(wp)         ,pointer :: L_h      (:)   ! hor.length scale for geop.h.
    real(wp)         ,pointer :: sqnu     (:)   ! sqrt(nu)   wind-velocity.p
    real(wp)                  :: sq1nu   =-9._wp! sqrt(1-nu) wind-streamf.
    integer                   :: oi_comp = 0    ! OI compatibility flag
    !----------------------------------------------------------------------
    ! NMC-derived meridional profiles: standard deviations, coupling z-psi
    !----------------------------------------------------------------------
    real(wp), pointer :: mu        (:) => NULL () ! Geostrophic coupling
    real(wp), pointer :: sdev_h    (:) => NULL () ! Height
    real(wp), pointer :: sdev_rh   (:) => NULL () ! Rel.humidity
    real(wp), pointer :: sdev_chi  (:) => NULL () ! Velocity potential
    real(wp), pointer :: sdev_psi  (:) => NULL () ! Streamfunction (full)
    real(wp), pointer :: sdev_psi_u(:) => NULL () ! Streamfunction (unbalanced)
    real(wp), pointer :: sdev_u    (:) => NULL () ! u-wind
    real(wp), pointer :: sdev_v    (:) => NULL () ! v-wind (usually == sdev_u)
  end type t_cov

  !===========
  ! Interfaces
  !===========
  interface destruct
    module procedure destruct_t_cov
  end interface

  interface construct
    module procedure construct_t_cov
  end interface

  interface set_horinpc
    module procedure set_horinpc0
    module procedure set_horinpc1
  end interface set_horinpc

  interface set_vertinpc
    module procedure set_vertinpc0
    module procedure set_vertinpc1
  end interface set_vertinpc

  interface random_gauss
    module procedure random_gauss_covm
  end interface random_gauss

  !=================
  ! module variables
  !=================
  type (t_cov) ,target ,save :: covm ! covariance matrix meta information
contains
!==============================================================================
  subroutine clean_bg_err_op
  !--------------------------------------------------
  ! deallocate pointer components of module variables
  !--------------------------------------------------
    call destruct (covm)
  end subroutine clean_bg_err_op
!------------------------------------------------------------------------------
  subroutine destruct_t_cov (cov)
  !-----------------------------------
  ! destruct variables of type t_cov .
  ! deallocate pointer components
  !-----------------------------------
  type (t_cov) ,intent(inout) :: cov
    integer :: j
    if (cov% valid == 0) return
    if (associated (cov% sqcv1h )) deallocate (cov% sqcv1h )
    if (associated (cov% sqcv1s )) deallocate (cov% sqcv1s )
    if (associated (cov% sqcv1v )) deallocate (cov% sqcv1v )
    if (associated (cov% sqcv1q )) deallocate (cov% sqcv1q )
    if (associated (cov% sqcv1t )) deallocate (cov% sqcv1t )
    if (associated (cov% eh_cl  )) deallocate (cov% eh_cl  )
    if (associated (cov% et_cl  )) deallocate (cov% et_cl  )
    if (associated (cov% erh_cl )) deallocate (cov% erh_cl )
    if (associated (cov% ev_cl  )) deallocate (cov% ev_cl  )
    if (associated (cov% epsi_cl)) deallocate (cov% epsi_cl)
    if (associated (cov% exi_cl )) deallocate (cov% exi_cl )
    deallocate (cov% sqcvh)
    deallocate (cov% sqcvs)
    deallocate (cov% sqcvv)
    deallocate (cov% sqcvq)
    deallocate (cov% sqcvt)
    deallocate (cov% eh   )
    deallocate (cov% et   )
    deallocate (cov% erh  )
    deallocate (cov% ev   )
    deallocate (cov% epsi )
    deallocate (cov% exi  )
    deallocate (cov% p    )
    deallocate (cov% logp )
    deallocate (cov% dlat )
    cov% valid         =  0
    cov% lclim         =  0
    cov% hc2_totalsize =  0
    cov% hc2_part_size =  0
    cov% gridtype      = -1
    if (associated (cov% c_h_psi   )) deallocate (cov% c_h_psi)
    if (associated (cov% c1_h_psi  )) deallocate (cov% c1_h_psi)
    if (associated (cov% hc2_part  )) deallocate (cov% hc2_part)
    if (associated (cov% L_h       )) deallocate (cov% L_h)
    if (associated (cov% sqnu      )) deallocate (cov% sqnu)
    if (associated (cov% rndm_state)) then
      call destruct (cov% rndm_state)
      deallocate    (cov% rndm_state)
    endif
    do j=1, size (cov% hc2)
       call destruct (cov% hc2(j))
    end do
    if (associated (cov% mu))       deallocate (cov% mu)
    if (associated (cov% sdev_h    )) deallocate (cov% sdev_h)
    if (associated (cov% sdev_rh   )) deallocate (cov% sdev_rh)
    if (associated (cov% sdev_chi  )) deallocate (cov% sdev_chi)
    if (associated (cov% sdev_psi  )) deallocate (cov% sdev_psi)
    if (associated (cov% sdev_psi_u)) deallocate (cov% sdev_psi_u)
    if (associated (cov% sdev_u    )) deallocate (cov% sdev_u)
    if (associated (cov% sdev_v    )) deallocate (cov% sdev_v)
  end subroutine destruct_t_cov
!------------------------------------------------------------------------------
  subroutine construct_t_cov (cov, nx, ny, nz, pbot, nwv, time, lclim, plevels)
  type (t_cov) ,intent(out) :: cov        ! covariance data to allocate
  integer      ,intent(in)  :: nx, ny, nz ! no. longitudes,latitudes,levels
  real(wp)     ,intent(in)  :: pbot       ! pressure at bottom (Pa)
  integer      ,intent(in)  :: nwv        ! no.points in vertical interpolation
  type (t_time),intent(in)  :: time       ! date&time of analysis
  integer      ,intent(in)  :: lclim      ! climatological error flag
  real(wp)     ,intent(in), optional :: plevels(:) ! pressure (Pa), full levels
  !-------------------------------------
  ! construct covariance model data type
  !-------------------------------------
    real(wp) :: dz, dy, ptop, ltop, lbot
    integer  :: k, j
    !---------------
    ! set dimensions
    !---------------
    cov% nx   = nx
    cov% ny   = ny
    cov% nz   = nz
    cov% pbot = pbot
    cov% nwv  = nwv
    cov% nwh  = 2
    cov% time = time
    !----------------------------
    ! allocate pointer components
    !----------------------------
    nullify  (cov% sqcv1h)
    nullify  (cov% sqcv1s)
    nullify  (cov% sqcv1v)
    nullify  (cov% sqcv1q)
    nullify  (cov% sqcv1t)
    allocate (cov% sqcvh (nz, nz, ny))
    allocate (cov% sqcvs (nz, nz, ny))
    allocate (cov% sqcvv (nz, nz, ny))
    allocate (cov% sqcvq (nz, nz, ny))
    allocate (cov% sqcvt (nz, nz, ny))
    allocate (cov% eh        (nz, ny))
    allocate (cov% et        (nz, ny))
    allocate (cov% erh       (nz, ny))
    allocate (cov% ev        (nz, ny))
    allocate (cov% epsi      (nz, ny))
    allocate (cov% exi       (nz, ny))
    allocate (cov% p     (nz)        )
    allocate (cov% logp  (nz)        )
    allocate (cov% dlat          (ny))
    nullify  (cov% hc2_part)
    nullify  (cov% rndm_state)
    cov% lclim = lclim
    if (lclim == 0) then
      nullify  (cov% eh_cl  )
      nullify  (cov% et_cl  )
      nullify  (cov% erh_cl )
      nullify  (cov% ev_cl  )
      nullify  (cov% epsi_cl)
      nullify  (cov% exi_cl )
    else
      allocate (cov% eh_cl   (nz, ny))
      allocate (cov% et_cl   (nz, ny))
      allocate (cov% erh_cl  (nz, ny))
      allocate (cov% ev_cl   (nz, ny))
      allocate (cov% epsi_cl (nz, ny))
      allocate (cov% exi_cl  (nz, ny))
    endif
    if (present (plevels)) then
       !-------------------------------------------
       ! check consistency of given pressure levels
       !-------------------------------------------
       if (size (plevels) /= nz) then
          write (0,*) size (plevels), nz
          call finish ("construct_t_cov","size (plevels) /= nz")
       end if
       if (any (plevels (1:nz-1) >= plevels (2:nz))) then
           write (0,'(999f10.2)') plevels
          call finish ("construct_t_cov","plevels not monotonously increasing")
       end if
       if (abs (maxval (plevels) - pbot) / pbot > 1.e-6) then
          write (0,*) maxval (plevels), pbot
          call finish ("construct_t_cov","maxval (plevels) /= pbot")
       end if
       cov% p   (:) =      plevels(:)
       cov% logp(:) = log (plevels(:))
    else
       !-------------------------------------------
       ! fallback: levels are equidistant in log(p)
       !-------------------------------------------
       ptop = 1000._wp             ! top level pressure = 10 hPa
       ltop = log(ptop)
       lbot = log(pbot)
       dz   = (lbot-ltop) / (nz-1) ! log(p) increment
!NEC$ ivdep
       do k = 1, nz
          cov% logp(k) = ltop + (k-1) * dz
          cov% p   (k) = exp(cov% logp(k))
       end do
    end if
    !------------------------------------------------------------
    ! derive latitudes (regular for odd ny, Gaussian for even ny)
    !------------------------------------------------------------
    select case (mod(ny,2))
    case (1)
      dy = 180._wp / (ny-1)
      do j = 1, ny
        cov% dlat(j) = -90._wp + dy * (j-1)
      end do
    case (0)
      call gauaw (ga=cov% dlat)
!NEC$ ivdep
      cov% dlat = - 90._wp/asin(1._wp) * asin (cov% dlat)
    end select
    !------------------------------------------------------------
    ! We currently defer the initialization of the horizontal
    ! wavelet covariances until they are read from file
    !------------------------------------------------------------
    nullify (cov% c_h_psi)
    nullify (cov% c1_h_psi)
    nullify (cov% hc2_part)
    nullify (cov% rndm_state)
    nullify (cov% L_h)
    nullify (cov% sqnu)
    cov% hc2_totalsize = 0
    cov% hc2_part_size = 0
    do j=1, size (cov% hc2)
       call construct (cov% hc2(j))
    end do
    nullify (cov% mu)
    nullify (cov% sdev_h)
    nullify (cov% sdev_rh)
    nullify (cov% sdev_chi)
    nullify (cov% sdev_psi)
    nullify (cov% sdev_psi_u)
    nullify (cov% sdev_u)
    nullify (cov% sdev_v)
  end subroutine construct_t_cov
!------------------------------------------------------------------------------
  subroutine cov_mem (bytes)
  !------------------------
  ! report allocated memory
  !------------------------
  integer(i8), intent(out) :: bytes

    integer :: n = 0
!   if (n==0) n = size(transfer(1._wp,(/''/)))
    n = 8

    bytes = 0

    if (covm% valid <= 0) return
    if (associated (covm% sqcvh )) bytes =          5 * n * size (covm% sqcvh)
    if (associated (covm% sqcv1h)) bytes = bytes  + 5 * n * size (covm% sqcv1h)
    if (associated (covm% eh))     bytes = bytes  + 6 * n * size (covm% eh)
    if (associated (covm% eh_cl))  bytes = bytes  + 6 * n * size (covm% eh)
    bytes = bytes + n*(size (covm% p) + size (covm% logp) + size (covm% dlat))
  end subroutine cov_mem
!------------------------------------------------------------------------------
  subroutine compress_cov

    real(wp) _POINTER :: tmp (:,:,:)

    integer :: nz, ny
!   integer :: p1, np
!   integer :: nt                ! Compressed section size (this MPI rank)
    integer :: pm                ! no. processors keeping "compressed" data
    integer :: lb(0:dace% npe-1) ! Lower bounds
    integer :: ub(0:dace% npe-1) ! Upper bounds
#ifdef NO_OPT_MPI_ATOA
    return
#endif
    if (dace% npe   <= 1) return
    if (covm% valid == 0) return
    if (covm% cmpr) return
    covm% cmpr = .true.

    nz = covm% nz
    ny = covm% ny
!   p1 = dace% pe+1
!   np = dace% npe
!   nt = (ny-p1+np)/np  ! Be careful about p1 > ny!

    call get_compress_params (ny, pm, lb, ub)

    call compress (covm% sqcvh)
    call compress (covm% sqcvs)
    call compress (covm% sqcvv)
    call compress (covm% sqcvq)
    call compress (covm% sqcvt)
    call compress (covm% sqcv1h)
    call compress (covm% sqcv1s)
    call compress (covm% sqcv1v)
    call compress (covm% sqcv1q)
    call compress (covm% sqcv1t)

  contains

!!$    subroutine compress (x)
!!$      real(wp), pointer :: x(:,:,:)
!!$
!!$      if (associated (x)) then
!!$         allocate (tmp (nz, nz, nt))
!!$         tmp =  x (:, :, p1:ny:np)
!!$         deallocate (x)
!!$         x   => tmp
!!$      endif
!!$    end subroutine compress

    subroutine compress (x)
      real(wp) _POINTER :: x(:,:,:)

      if (associated (x)) then
         allocate (tmp (nz, nz, lb(dace% pe):ub(dace% pe)))
         tmp =       x ( :,  :, lb(dace% pe):ub(dace% pe))
         deallocate (x)
         x   => tmp
      endif
    end subroutine compress

  end subroutine compress_cov
!------------------------------------------------------------------------------
  subroutine get_compress_params (ny, pm, lb, ub)
    integer, intent(in)  :: ny        ! No. elements to "compress"
    integer, intent(out) :: pm        ! no. processors keeping "compressed" data
    integer, intent(out) :: lb(0:)    ! Lower bounds
    integer, intent(out) :: ub(0:)    ! Upper bounds

    integer  :: pe            ! Processor index
    real     :: inc           ! "Slices" per pe

    pm    = min (dace% npe, ny)
    inc   = real (ny) / dace% npe
    lb(0) = 1
    do pe = 0, dace% npe-2
       ub(pe)   = nint ((pe+1) * inc)
       lb(pe+1) = ub(pe) + 1
    end do
    ub(dace% npe-1) = ny
    if (any (ub - lb < -1) .or. (sum (ub - lb + 1) /= ny)) then
       if (dace% lpio) then
          write(0,*) "ny:", ny
          write(0,*) "lb:", lb
          write(0,*) "ub:", lb
       end if
       call p_barrier ()
       call finish ("get_compress_params","internal bounds checking error")
    end if
  end subroutine get_compress_params
!------------------------------------------------------------------------------
! Uncompression of covariance data can be done either using multiple (i)bcasts
! or using alltoallv.  On NEC SX, alltoallv was empirically much slower.
!#define PREFER_ALLTOALLV
!------------------------------------------------------------------------------
  subroutine uncompress_cov

    real(wp) _POINTER     :: tmp(:,:,:)
    real(wp) ,allocatable :: buf(:,:,:) ! Contiguous work area/bcast buffer
#ifdef _CRAYFTN
    real(wp) ,allocatable :: x_ (:,:,:) ! Temporary for crayftn, vectorization
#endif

    integer :: nz, ny, pm
    integer :: lb(0:dace% npe-1) ! Lower bounds
    integer :: ub(0:dace% npe-1) ! Upper bounds

#ifdef PREFER_ALLTOALLV
    integer :: rcount(0:dace% npe-1), scount(0:dace% npe-1) ! Counts
    integer :: rdispl(0:dace% npe-1), sdispl(0:dace% npe-1) ! Displacements
#endif

    if (dace% npe   <= 1) return
    if (covm% valid == 0) return
    if (.not. covm% cmpr) return
    covm% cmpr = .false.

    nz = covm% nz
    ny = covm% ny

    call get_compress_params (ny, pm, lb, ub)

    allocate (buf(nz,nz,ny))
#ifdef _CRAYFTN
    allocate (x_ (nz,nz,ny))
#endif

    call uncompress (covm% sqcvh)
    call uncompress (covm% sqcvs)
    call uncompress (covm% sqcvv)
    call uncompress (covm% sqcvq)
    call uncompress (covm% sqcvt)
    call uncompress (covm% sqcv1h)
    call uncompress (covm% sqcv1s)
    call uncompress (covm% sqcv1v)
    call uncompress (covm% sqcv1q)
    call uncompress (covm% sqcv1t)

    deallocate (buf)

  contains

    subroutine uncompress (x)
      real(wp) _POINTER :: x(:,:,:)

#ifdef PREFER_ALLTOALLV
      integer  :: nn, ierr
#else
      integer  :: pe
#endif

      if (.not. associated (x)) return

      tmp => x
      allocate (x (nz, nz, ny))

#ifdef PREFER_ALLTOALLV
      !----------------------------------------------
      ! Emulate "all-to-all-bcast" for data expansion
      !----------------------------------------------
      nn     = size (x, dim=1) * size (x, dim=2)
      sdispl = 0
      scount = (ub(dace% pe) - lb(dace% pe) + 1) * nn   ! "bcast"
      rcount = (ub(    :   ) - lb(    :   ) + 1) * nn
      rdispl = (lb(    :   )                - 1) * nn
      if (sum (rcount) /= size (x)) &
          call finish ("uncompress_cov","sum (rcount) /= size (x)")
#ifndef NOMPI
      call MPI_Alltoallv (tmp, scount, sdispl, p_real_wp, &
                          x,   rcount, rdispl, p_real_wp, &
                          dace% comm, ierr                )
#endif
      deallocate (tmp)

#else /* !PREFER_ALLTOALLV */
      !-------------------------
      ! Prepare broadcast buffer
      !-------------------------
      if (lb(dace% pe) <= ub(dace% pe)) then
#ifdef _CRAYFTN
!DIR$ IVDEP
#endif
         buf(:,:,lb(dace% pe):ub(dace% pe)) = tmp
      end if
      deallocate (tmp)
      !------------
      ! Expand data
      !------------
!     do pe = 0, pm-1
      do pe = 0, dace% npe-1
         if (lb(pe) > ub(pe)) cycle
!if (dace% lpio) write(0,*) dace% pe, "uncompress: pe,lb,ub=", pe,lb(pe),ub(pe)

#ifndef NO_MPI3        /* we do have MPI >=3.0 */
         call p_ibcast (buf(:,:,lb(pe):ub(pe)), pe)
#else                  /* fallback for MPI 2.x */
         call p_bcast  (buf(:,:,lb(pe):ub(pe)), pe)
#endif

      end do
      call p_waitall ()
      !---------------------------------------
      ! Extract/reorder data from bcast buffer
      !---------------------------------------
#ifdef _CRAYFTN     /* Argh! Crayftn can be sometimes rather dull... */
!$omp parallel workshare
!DIR$ IVDEP
      x_ = buf   ! <-- ok
!DIR$ IVDEP
      x  = x_    ! <-- ok (if x is contiguous)
!$omp end parallel workshare
#else
      x  = buf   ! <-- slow with crayftn 8.3.6/8.4.0
#endif

#endif /* !PREFER_ALLTOALLV */
    end subroutine uncompress

  end subroutine uncompress_cov
!------------------------------------------------------------------------------
  subroutine cov_from_nmc (file, nx, ny, nz, pbot, h2psi, nwv, time, lclim, &
                           scale_rh)
  character(len=*) ,intent(in)    :: file  ! file to read covariances from
  integer          ,intent(in)    :: nx    ! no. longitudes
  integer          ,intent(inout) :: ny    ! no. latitudes
  integer          ,intent(inout) :: nz    ! no. levels
  real(wp)         ,intent(inout) :: pbot  ! pressure at bottom (hPa)
  logical          ,intent(in)    :: h2psi ! derive streamf.cov. from height
  integer          ,intent(in)    :: nwv   ! no. points in vertical interpolat.
  type (t_time)    ,intent(in)    :: time  ! date&time of analysis
  integer          ,intent(in)    :: lclim ! climatological error use flag
  real(wp)         ,intent(in)    :: scale_rh ! scaling factor for RH fg.error
  !--------------------------------------------------------------------
  ! read the NMC derived vertical covariances from GRADS or NetCDF file
  !--------------------------------------------------------------------
    !----------------
    ! local variables
    !----------------
    type (t_ctl)          :: ctl           ! GRADS meta data
    character(len=256)    :: fname         ! full filename
    real(wp) ,allocatable :: tmp(:,:,:)    ! temporary (transposed array)
    integer               :: i, j, k, m    ! indices
    logical               :: zrev = .true. ! true for p(1) > p(nz)
    logical               :: has_pl        ! GrADS file provides p levels
    real(wp), allocatable :: plevels(:)    ! pressure levels in hPa
    character(len=256)    :: description
    !---------------------------
    ! local variables for NetCDF
    !---------------------------
    logical               :: lnetcdf       ! true for NetCDF file to read
    integer               :: status        ! status return variable
    integer               :: ncid          ! NetCDF file id
    integer               :: d_z, d_y, d_t ! NetCDF dimension id
    integer               :: varid         ! NetCDF variable id
    !------------
    ! read header
    !------------
    i = index (file, '.nc',.true.)
    lnetcdf = i > 0 .and. i == len_trim(file) - 2
    if (dace% lpio) then
      !---------
      ! printout
      !---------
      write(6,'()')
      write(6,'(a)') repeat('-',79)
      write(6,'()')
      write(6,'(a,a)') '  cov_from_nmc: ',trim(file)
      if (lnetcdf) then
        !----------------------
        ! read from NetCDF file
        !----------------------
        status = nf90_open (file, NF90_NOWRITE, ncid)
        if (status /= NF90_NOERR) call finish &
                                  ('cov_from_nmc','cannot open '//trim(file))
        status = nf90_get_att (ncid, NF90_GLOBAL, 'description', description)
        if (status == NF90_NOERR) then
           write(6,'(A,/,A)') "  Description:", trim (description)
        end if
        call chk (nf90_inq_dimid (ncid, 'latitude' ,d_y))
        call chk (nf90_inq_dimid (ncid, 'pressure' ,d_z))
        call chk (nf90_inq_dimid (ncid, 'month'    ,d_t))
        call chk (nf90_inquire_dimension (ncid, d_y, len=ny))
        call chk (nf90_inquire_dimension (ncid, d_z, len=nz))
        has_pl = .true.
        allocate (plevels (nz))
        call chk (nf90_inq_varid (ncid, 'pressure', varid))
        call chk (nf90_get_var   (ncid, varid, plevels   ))
        zrev   = plevels(1) > plevels(nz)
        if(zrev) call finish ('cov_from_nmc','zrev')
        pbot = maxval (plevels(:)) / 100._wp
      else
        !---------------------
        ! read from GRADS file
        !---------------------
        fname = trim(file)//'_geof.ctl'
        call read_ctl (ctl, fname)
        !------------------
        ! adjust parameters
        !------------------
        ny = ctl% ydefn
        nz = ctl% zdefn
        if (associated (ctl% zlevels)) then
          pbot = maxval (ctl% zlevels(:))
          has_pl = .true.
          allocate (plevels (nz))
          if (ctl% zlevels(1) > ctl% zlevels(nz)) zrev = .true.
          if (zrev) then
            !-----------------------------------------------------------------
            ! Reorder to match model levels: 1=top, nz=bottom; convert hPa->Pa
            !-----------------------------------------------------------------
            plevels(:) = ctl% zlevels(nz:1:-1) * 100
          else
            plevels(:) = ctl% zlevels(1:nz)    * 100
          end if
          write (*,*) "using pressure levels from ", trim (fname)
        else
          has_pl = .false.
          write (0,*) "cov_from_nmc: WARNING: pressure levels missing in GrADS!"
        end if
        if (all(ctl% var% name /= 'geof_varcl' .and. lclim > 0)) &
          call finish('cov_from_nmc',' lclim> 0 and climatology not in dataset')
        call destruct (ctl)
      endif
    endif
    call p_bcast (ny,     dace% pio)
    call p_bcast (nz,     dace% pio)
    call p_bcast (pbot,   dace% pio)
    call p_bcast (has_pl, dace% pio)
    if (has_pl) then
       if (.not.dace% lpio) allocate (plevels (nz))
       call p_bcast (plevels, dace% pio)
    end if
    !------------------------------------
    ! allocate vertical covariance arrays
    !------------------------------------
    if (has_pl) then
       call construct (covm, nx, ny, nz, pbot * 100._wp, nwv, time, lclim, &
            plevels)
    else
       call construct (covm, nx, ny, nz, pbot * 100._wp, nwv, time, lclim)
    end if
    !--------------------------
    ! read vertical covariances
    !--------------------------
    if (dace% lpio) then
      if (lnetcdf) then
        !-----------------
        ! read NetCDF file
        !-----------------
        m = imm (time)
        allocate (tmp         (nz, ny,1))
        !----------------------------------------------
        ! rel.humidity, check for reversed z-coordinate
        !----------------------------------------------
        call chk (nf90_inq_varid (ncid, 'var_rh', varid    ))
        call chk (nf90_get_var   (ncid, varid, covm% sqcvq, start=(/1,1,1,m/) ))
        covm% sqcvq = covm% sqcvq * scale_rh ** 2
        zrev = covm% sqcvq (2,2,ny/2) > covm% sqcvq (nz,nz,ny/2)
        if(zrev) call finish ('cov_from_nmc','zrev')
        forall (k=1:nz, j=1:ny) covm% erh (k,j) = sqrt(covm% sqcvq (k,k,j))
        if (covm% lclim > 0) then
          call chk (nf90_inq_varid (ncid, 'cler_rh', varid    ))
          call chk (nf90_get_var   (ncid, varid, covm% erh_cl, start=(/1,1,m/) ))
        endif
        !--------------------
        ! geopotential height
        !--------------------
        call chk (nf90_inq_varid (ncid, 'var_geof', varid    ))
        call chk (nf90_get_var   (ncid, varid, covm% sqcvh, start=(/1,1,1,m/) ))
        covm% sqcvh = covm% sqcvh / gacc**2     ! scale from (m**2/s**2)^2 to m^2
        forall (k=1:nz, j=1:ny) covm% eh (k,j) = sqrt(covm% sqcvh (k,k,j))
        if (covm% lclim > 0) then
          call chk (nf90_inq_varid (ncid, 'cler_geof', varid    ))
          call chk (nf90_get_var   (ncid, varid, covm% eh_cl, start=(/1,1,m/) ))
          covm% eh_cl = covm% eh_cl / gacc
        endif
        !--------------------------------------------
        ! streamfunction (copied from geopotential ?)
        !--------------------------------------------
        if (h2psi) then
          covm% sqcvs = covm% sqcvh
          covm% epsi  = covm% eh
          if (covm% lclim > 0) covm% epsi_cl = covm% eh_cl
        else
          call chk (nf90_inq_varid (ncid, 'var_psi', varid    ))
          call chk (nf90_get_var   (ncid, varid, covm% sqcvs, start=(/1,1,1,m/) ))
          forall (k=1:nz, j=1:ny) covm% epsi (k,j) = sqrt(covm% sqcvs (k,k,j))
          if (covm% lclim > 0) then
            call chk (nf90_inq_varid (ncid, 'cler_psi', varid    ))
           call chk (nf90_get_var   (ncid, varid, covm% epsi_cl, start=(/1,1,m/) ))
          end if
        endif
        !-------------------
        ! velocity potential
        !-------------------
        call chk (nf90_inq_varid (ncid, 'var_chi', varid    ))
        call chk (nf90_get_var   (ncid, varid, covm% sqcvv, start=(/1,1,1,m/) ))
        forall (k=1:nz, j=1:ny) covm% exi (k,j) = sqrt(covm% sqcvv (k,k,j))
        if (covm% lclim > 0) then
          call chk (nf90_inq_varid (ncid, 'cler_chi', varid    ))
         call chk (nf90_get_var   (ncid, varid, covm% exi_cl, start=(/1,1,m/) ))
        end if
        !---------------------------------
        ! temperature (for variances only)
        !---------------------------------
        call chk (nf90_inq_varid (ncid, 'err_t', varid    ))
        call chk (nf90_get_var   (ncid, varid, covm% et, start=(/1,1,m/) ))
        if (covm% lclim > 0) then
          call chk (nf90_inq_varid (ncid, 'cler_t', varid    ))
          call chk (nf90_get_var   (ncid, varid, covm% et_cl, start=(/1,1,m/) ))
        end if
        covm% sqcvt = 0._wp
        !-------------------------------------
        ! wind components (for variances only)
        !-------------------------------------
        call chk (nf90_inq_varid (ncid, 'err_u', varid    ))
        call chk (nf90_get_var   (ncid, varid, covm% ev, start=(/1,1,m/) ))
        call chk (nf90_inq_varid (ncid, 'err_v', varid    ))
        call chk (nf90_get_var   (ncid, varid, tmp, start=(/1,1,m/) ))
        covm% ev = sqrt (0.5_wp * (covm% ev **2 + tmp(:,:,1)**2))
        if (covm% lclim > 0) then
          call chk (nf90_inq_varid (ncid, 'cler_u', varid    ))
          call chk (nf90_get_var   (ncid, varid, covm% ev_cl, start=(/1,1,m/) ))
          call chk (nf90_inq_varid (ncid, 'cler_v', varid    ))
          call chk (nf90_get_var   (ncid, varid, tmp, start=(/1,1,m/) ))
          covm% ev_cl = sqrt (0.5_wp * (covm% ev_cl **2 + tmp(:,:,1)**2))
        endif
        call chk (nf90_close (ncid))

!call finish('cov_from_nmc','lnetcdf=T')

      else
        !----------------
        ! read GRADS file
        !----------------
        allocate (tmp         (nz, ny, nz))
        !----------------------------------------------
        ! rel.humidity, check for reversed z-coordinate
        !----------------------------------------------
        fname = trim(file)//'_rh.ctl'
        call read_ctl (ctl, fname)
        call read_var (ctl, tmp, 'rh_var')
        tmp = tmp * scale_rh ** 2
        if (.not. has_pl) then
           ! Fallback if pressure levels were not provided in GrADS
           zrev = tmp(2,ny/2,2) > tmp(nz,ny/2,nz)
        end if
        call transp (tmp, covm% sqcvq, zrev)
!NEC$ ivdep
        forall (k=1:nz, j=1:ny) covm% erh (k,j) = sqrt(covm% sqcvq (k,k,j))
        if (covm% lclim > 0) then
          call read_var (ctl, tmp, 'rh_varcl')
          forall (k=1:nz, j=1:ny)  covm% erh_cl (k,j) = sqrt (tmp (k,j,k))
          if (zrev) covm% erh_cl = covm% erh_cl (nz:1:-1,:)
        endif
        call destruct (ctl)
        !--------------------
        ! geopotential height
        !--------------------
        fname = trim(file)//'_geof.ctl'
        call read_ctl (ctl, fname)
        call read_var (ctl, tmp, 'geof_var')
        call transp (tmp, covm% sqcvh, zrev)
        covm% sqcvh = covm% sqcvh / gacc**2     ! scale from (m**2/s**2)^2 to m^2
        forall (k=1:nz, j=1:ny) covm% eh (k,j) = sqrt(covm% sqcvh (k,k,j))
        if (covm% lclim > 0) then
          call read_var (ctl, tmp, 'geof_varcl')
          forall (k=1:nz, j=1:ny)  covm% eh_cl (k,j) = sqrt (tmp (k,j,k)) / gacc
          if (zrev) covm%  eh_cl = covm% eh_cl (nz:1:-1,:)
        endif
        call destruct (ctl)
        !--------------------------------------------
        ! streamfunction (copied from geopotential ?)
        !--------------------------------------------
        if (h2psi) then
          covm% sqcvs = covm% sqcvh
          covm% epsi  = covm% eh
          if (covm% lclim > 0) covm% epsi_cl = covm% eh_cl
        else
          fname = trim(file)//'_psi.ctl'
          call read_ctl (ctl, fname)
          call read_var (ctl, tmp, 'psi_var')
          call transp (tmp, covm% sqcvs, zrev)
!NEC$ ivdep
          forall (k=1:nz, j=1:ny) covm% epsi (k,j) = sqrt(covm% sqcvs (k,k,j))
          if (covm% lclim > 0) then
            call read_var (ctl, tmp, 'psi_varcl')
            forall (k=1:nz, j=1:ny)   covm% epsi_cl (k,j) = sqrt (tmp (k,j,k))
            if (zrev) covm% epsi_cl = covm% epsi_cl (nz:1:-1,:)
          endif
          call destruct (ctl)
        endif
        !-------------------
        ! velocity potential
        !-------------------
        fname = trim(file)//'_chi.ctl'
        call read_ctl (ctl, fname)
        call read_var (ctl, tmp, 'chi_var')
        call transp (tmp, covm% sqcvv, zrev)
!NEC$ ivdep
        forall (k=1:nz, j=1:ny) covm% exi (k,j) = sqrt(covm% sqcvv (k,k,j))
        if (covm% lclim > 0) then
          call read_var (ctl, tmp, 'chi_varcl')
          forall (k=1:nz, j=1:ny)   covm% exi_cl (k,j) = sqrt (tmp (k,j,k))
          if (zrev) covm%  exi_cl = covm% exi_cl (nz:1:-1,:)
        endif
        call destruct (ctl)
        !---------------------------------
        ! temperature (for variances only)
        !---------------------------------
        fname = trim(file)//'_t.ctl'
        call read_ctl (ctl, fname)
        call read_var (ctl, tmp, 't_var')
        call transp (tmp, covm% sqcvt, zrev)
        forall (k=1:nz, j=1:ny) covm% et (k,j) = sqrt(covm% sqcvt (k,k,j))
        if (covm% lclim > 0) then
          call read_var (ctl, tmp, 't_varcl')
          forall (k=1:nz, j=1:ny)  covm% et_cl (k,j) = sqrt (tmp (k,j,k))
          if (zrev) covm%  et_cl = covm% et_cl (nz:1:-1,:)
        endif
        call destruct (ctl)
        !-------------------------------------
        ! wind components (for variances only)
        !-------------------------------------
        fname = trim(file)//'_u.ctl'
        call read_ctl (ctl, fname)
        call read_var (ctl, tmp, 'u_var')
        forall (k=1:nz, j=1:ny) covm% ev (k,j) = tmp (k,j,k)
        if (covm% lclim > 0) then
          call read_var (ctl, tmp, 'u_varcl')
          forall (k=1:nz, j=1:ny)   covm% ev_cl (k,j) = tmp (k,j,k)
        endif
        call destruct (ctl)
        fname = trim(file)//'_v.ctl'
        call read_ctl (ctl, fname)
        call read_var (ctl, tmp, 'v_var')
        forall (k=1:nz, j=1:ny) covm% ev (k,j) = covm% ev (k,j) + tmp (k,j,k)
!NEC$ ivdep
        covm% ev = sqrt(0.5_wp * covm% ev)
        if (zrev) covm% ev = covm% ev (nz:1:-1,:)
        if (covm% lclim > 0) then
          call read_var (ctl, tmp, 'v_varcl')
          forall (k=1:nz, j=1:ny) covm% ev_cl (k,j) = &
                                  covm% ev_cl (k,j) + tmp (k,j,k)
!NEC$ ivdep
          covm% ev_cl = sqrt(0.5_wp * covm% ev_cl)
          if (zrev) covm%   ev_cl = covm% ev_cl (nz:1:-1,:)
        endif
        call destruct (ctl)
      endif
      !---------
      ! clean up
      !---------
      deallocate    (tmp)
    endif
    !----------
    ! broadcast
    !----------
    call p_bcast   (covm% sqcvh,   dace% pio)
    call p_bcast   (covm% sqcvs,   dace% pio)
    call p_bcast   (covm% sqcvv,   dace% pio)
    call p_bcast   (covm% sqcvq,   dace% pio)
    call p_bcast   (covm% sqcvt,   dace% pio)
    call p_bcast   (covm% eh,      dace% pio)
    call p_bcast   (covm% et,      dace% pio)
    call p_bcast   (covm% erh,     dace% pio)
    call p_bcast   (covm% ev,      dace% pio)
    call p_bcast   (covm% epsi,    dace% pio)
    call p_bcast   (covm% exi,     dace% pio)
    if (covm% lclim > 0) then
      call p_bcast (covm% eh_cl,   dace% pio)
      call p_bcast (covm% et_cl,   dace% pio)
      call p_bcast (covm% erh_cl,  dace% pio)
      call p_bcast (covm% ev_cl,   dace% pio)
      call p_bcast (covm% epsi_cl, dace% pio)
      call p_bcast (covm% exi_cl,  dace% pio)
    endif

!..............................................................................
    contains
      subroutine transp (tmp, x, zrev)
      real(wp) ,intent(in)  :: tmp (nz, ny, nz)
      real(wp) ,intent(out) :: x   (nz, nz, ny)
      logical  ,intent(in)  :: zrev
      !-------------------------------------
      ! transpose array read from GRADS-file
      !-------------------------------------
        integer :: k,l
        if (zrev) then
          do k = 1, nz
          do l = 1, nz
            x (l,k,:) = tmp (nz+1-l,:,nz+1-k)
          end do
          end do
        else
          do k = 1, nz
          do l = 1, nz
            x (l,k,:) = tmp (l,:,k)
          end do
          end do
        end if
      end subroutine transp
!..............................................................................
    subroutine chk (status)
    integer ,intent(in) :: status
      if (status /= NF90_NOERR) &
      call finish('cov_from_nmc','NetCDF read')
    end subroutine chk
  end subroutine cov_from_nmc
!------------------------------------------------------------------------------
  subroutine set_horinpc1 (hic, dlat, covm)
  !------------------------------------------
  ! set horizontal interpolation coefficients
  !------------------------------------------
  type (t_hic1) ,intent(out) :: hic  (:)
  real(wp)      ,intent(in)  :: dlat (:)
  type (t_cov)  ,intent(in)  :: covm
    integer :: i
    do i=1,size(hic)
      call set_horinpc (hic(i), dlat(i), covm)
    end do
  end subroutine set_horinpc1
!------------------------------------------------------------------------------
  subroutine set_horinpc0 (hic, dlat, covm)
  !------------------------------------------
  ! set horizontal interpolation coefficients
  !------------------------------------------
  type (t_hic1) ,intent(out) :: hic
  real(wp)      ,intent(in)  :: dlat
  type (t_cov)  ,intent(in)  :: covm

    integer :: j

    if (covm% nwh /= 2) call finish ('set_horinpc','covm% nwh /= 2')
    hic% w = 0._wp
    if (dlat <= covm% dlat(1)) then
      hic% ix           = 1
      hic% w(1)         = 1._wp
    else if (dlat >= covm% dlat(covm% ny)) then
      hic% ix           = covm% ny - covm% nwh + 1
      hic% w(covm% nwh) = 1._wp
    else
      do j = covm% ny -1, 1, -1
        if (dlat >= covm% dlat(j)) then
          hic% ix       = j
          hic% w(2)     = (      dlat      - covm% dlat(j))&
                         / (covm% dlat(j+1) - covm% dlat(j))
          hic% w(1)     = 1._wp - hic% w (2)
          exit
        endif
      end do
    endif
    if (hic% ix == 0) call finish ('set_horinpc','dlat out of bounds')

  end subroutine set_horinpc0
!------------------------------------------------------------------------------
  subroutine set_vertinpc1 (vic, logp, logpg, nwv)
  !----------------------------------------
  ! set vertical interpolation coefficients
  !----------------------------------------
  type (t_vic) ,intent(out) :: vic   (:) ! interpolation coefficients to set
  real(wp)     ,intent(in)  :: logp  (:) ! vertical coordinates of observations
  real(wp)     ,intent(in)  :: logpg (:) ! vertical coordinates of grid
  integer      ,intent(in)  :: nwv       ! # of coefficients, 4:cubic interpol.
  !------------------------------------------------------
  ! define EXTRAPOLATE_TOP for extrapolation at model top
!#define EXTRAPOLATE_TOP
  !------------------------------------------------------
    integer :: k, n, nz, nv, nt

    n  = size (logp)
    nz = size (logpg)
    nv = size (vic)
    nt = 1; if (ndv_rad==1) nt = 2
    select case (nwv)
    case (4)
      !---------------------------------
      ! Cubic interpolation coefficients
      !---------------------------------
      call set_vertinpc_cub ()
    case (2)
      !---------------------------
      ! linear in lnp
      ! fallback to scalar version
      !---------------------------
      do k = 1, size (logp)
        call set_vertinpc (vic(k), logp(k), logpg, nwv)
      end do
    case (3)
      !-----------------------------------------------
      ! interpolation cosidering piecewise integration
      !-----------------------------------------------
      call set_vertinpc_rochon()
    case default
      call finish('set_vertinpc1','not implemented: nwv = '//char2(nwv))
    end select
  contains
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    subroutine set_vertinpc_cub ()
      !-----------------------------------------------------------------
      ! Vectorized version for setup of cubic interpolation coefficients
      !-----------------------------------------------------------------
      integer                     :: i, j, l, nbi
      integer,  dimension(n)      :: klo   !, khi
      real(wp), dimension(0:nz+1) :: lnp
      real(wp)                    :: dz, dz1, den, e2, f2, x(4), z(4)
      integer                     :: js, je     ! Indices for strip-mining
      integer  ,parameter         :: VL  = 256  ! Vector length
      integer  ,dimension(VL)     :: k_hi, k_lo ! Vector registers
!NEC$ vreg(k_lo)
!NEC$ vreg(k_hi)
      real(wp) ,parameter :: eps = 1.e-10_wp    ! tolerance for upper bound

      !------------------------------------------------------------------------
      ! Set up auxiliary vector for bisection search, taking care of boundaries
      ! and consistency with comparisons in subroutine set_vertinpc0.
      !------------------------------------------------------------------------
!cdirr on_adb(logpg)
!cdirr on_adb(klo)
!cdirr on_adb(lnp)
!cdirr on_adb(logp)
      lnp(0)      = -HUGE (0._wp)
      lnp(1)      =           logpg(1)      - eps
!NEC$ SHORTLOOP
      lnp(2:nz-2) =           logpg(2:nz-2)
      lnp(  nz-1) =  NEAREST (logpg(  nz-1), +1._wp)
      lnp(  nz  ) =  NEAREST (logpg(  nz  ), +1._wp)
      lnp(  nz+1) =  HUGE (0._wp)

FTRACE_BEGIN("set_vertinpc_cub:bisection")
      nbi = ceiling (log (nz + 1.0)/log (2.0) - 4*epsilon (0.0))
      !-------------------------
      ! Outer loop: strip-mining
      !-------------------------
      do js = 1, n, VL
         je = min (js+VL-1, n)
         l  = je - js + 1
         !--------------------------
         ! Bisection search on strip
         !--------------------------
!NEC$ SHORTLOOP
         k_lo(1:l) = 0
!NEC$ SHORTLOOP
         k_hi(1:l) = nz + 1
         do i = 1, nbi
!NEC$ SHORTLOOP
!NEC$ IVDEP
            do j = 1, l
               if (k_hi(j) - k_lo(j) > 1) then
                  k = (k_hi(j) + k_lo(j))/2
                  if (lnp(k) > logp(js-1+j)) then
                     k_hi(j) = k
                  else
                     k_lo(j) = k
                  end if
               end if
            end do
         end do
!NEC$ SHORTLOOP
         klo(js:js-1+l) = k_lo(1:l)
!NEC$ SHORTLOOP
!        khi(js:js-1+l) = k_hi(1:l)
      end do
FTRACE_END  ("set_vertinpc_cub:bisection")
!     if (any (khi-klo /= 1)) call finish ("set_vertinpc_cub","bad khi,klo!")

!NEC$ IVDEP
      do j = 1, n
!NEC$ UNROLL(mwv)
         vic(j)% wh   = 0._wp
!NEC$ UNROLL(mwv)
         vic(j)% wt   = 0._wp
         vic(j)% g% n = 4
         k = klo(j)
#ifndef EXTRAPOLATE_TOP
         if (k == 0) then
            vic(j)% g% i    = 1
            vic(j)% wh(1)   = 1._wp
         else if (k == nz  ) then
#else
         if      (k == nz  ) then
#endif
            vic(j)% g% i    = nz-3
            vic(j)% wh(4)   = 1._wp
         else if (k <= 1   ) then
            vic(j)% g% i    = 1
            dz1             = logp (j)  - logpg(1)
            dz              = logpg(2)  - logpg(1)
            vic(j)% wh(2)   = dz1   / dz
            vic(j)% wh(1)   = 1._wp - vic(j)% wh(2)
            vic(j)% wt(2)   = 1._wp / dz
            vic(j)% wt(1)   =       - vic(j)% wt(2)
         else if (k == nz-1) then
            vic(j)% g% i    = nz-3
            dz1             = logp (j)  - logpg(nz-1)
            dz              = logpg(nz) - logpg(nz-1)
            vic(j)% wh(4)   = dz1   / dz
            vic(j)% wh(3)   = 1._wp - vic(j)% wh(4)
            vic(j)% wt(4)   = 1._wp / dz
            vic(j)% wt(3)   =       - vic(j)% wt(4)
         else    ! 2 <= k <= nz-2
            vic(j)% g% i    = k - 1
            z(1:4)          = logpg(k-1:k+2)
            x(1:4)          = (logp(j) - z(1:4))
            dz              =  z(3) - z(2)

            den             = (z(3)-z(1)) * dz**2
            vic(j)% wh(1)   = - x(2) * x(3)**2          / den
            vic(j)% wt(1)   = - x(3) * (2*x(2) + x(3))  / den

            den             = (z(4)-z(2)) * dz**2
            vic(j)% wh(4)   =   x(3) * x(2)**2          / den
            vic(j)% wt(4)   =   x(2) * (2*x(3) + x(2))  / den

            den             = dz**3
            e2              =   x(3)**2 * (3*x(2)-x(3)) / den
            f2              =   6 * x(2) * x(3)         / den
            vic(j)% wh(2)   = - vic(j)% wh(4) +      e2
            vic(j)% wh(3)   = - vic(j)% wh(1) + (1 - e2)
            vic(j)% wt(2)   = - vic(j)% wt(4) +      f2
            vic(j)% wt(3)   = - vic(j)% wt(1) -      f2
         end if
      end do
    end subroutine set_vertinpc_cub
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    subroutine set_vertinpc_rochon
    !-------------------------------------------------
    ! Interpolation considering piecewise integration.
    ! Calls the respective RTTOV10 routine.
    ! Cf. Rochon et al. (2007) QJRMS 133 1547-1558
    ! and the RTTOV10 doku
    !-------------------------------------------------
#if (_RTTOV_VERSION <= 0)
      call finish ('set_vertinpc_rochon','RTTOV not configured')
#else
#include "rttov_layeravg.interface"
#include "rttov_layeravg_tl.interface"
      real(wp) :: pz (size(logpg), size(logp),nt) ! interpolation coefficients
      integer  :: kstart  (size(logp))            ! start index ..
      integer  :: kend    (size(logp))            ! end   index ..
      integer  :: kn      (size(logp))            ! .. no.of relevant coeffs.
      integer  :: kmax                            ! max no.of coefficients
      integer  :: k, kg, k3                       ! loop index
      real(wp) :: logp_tl (size(logp ))
      real(wp) :: logpg_tl(size(logpg))
      real(wp) :: d (size(logpg),-1:1)
      real(wp) :: z, z2, zz2, dd, edd
      !---------------------------------------
      ! call the RTTOV10 interpolation routine
      !---------------------------------------
      select case (nt)
      case (1)
        !--------------------------
        ! direct interpolation only
        !--------------------------
        call rttov_layeravg (logp         ,&! <-  levels of output domain
                             logpg        ,&! <-  levels of input  domain
                             size(logp)   ,&! <-  output domain size
                             size(logpg)  ,&! <-  input  domain size
                             pz(:,:,1)    ,&!  -> interpolation coefficients
                             kstart       ,&!  -> start index
                             kend          &!  -> end   index
#if (_RTTOV_VERSION >= 12)
                             ,vint_rttov_b &! <-  interpolation mode (see rttov users guide)
#endif
                             )
      case (2)
        !------------------------------------------
        ! plus vertical differentiation by rttov TL
        !------------------------------------------
        logp_tl  = 1._wp
        logpg_tl = 0._wp
        call rttov_layeravg_tl (logp         ,&! <-  levels of output domain
                                logp_tl      ,&! <-  tl input
                                logpg        ,&! <-  levels of input  domain
                                logpg_tl     ,&! <-  tl input
                                size(logp)   ,&! <-  output domain size
                                size(logpg)  ,&! <-  input  domain size
                                pz(:,:,1)    ,&!  -> interpolation coefficients
                                pz(:,:,2)    ,&!  -> tl output
                                kstart       ,&!  -> start index
                                kend          &!  -> end   index
#if (_RTTOV_VERSION >= 12)
                                ,vint_rttov_b &! <-  interpolation mode (see rttov users guide)
#endif
                                )
      end select

      !----------------------------------
      ! explicit vertical differentiation
      !----------------------------------
      if (nt==1) then
        select case (ndv_rad)
        case (0)
          d (:,-1) = 0
          d (:, 0) = 1
          d (:, 1) = 0
        case (1)
        case (2)
          kg = 1
          d (kg,  1) = 1._wp / (logpg(kg+1)-logpg(kg))
          d (kg,  0) = - d (kg,  1)
          d (kg, -1) = 0._wp
          do kg = 2, nz - 1
            d (kg,  1) = 1._wp / (logpg(kg+1)-logpg(kg-1))
            d (kg,  0) = 0._wp
            d (kg, -1) = - d (kg,  1)
          end do
          kg = nz
          d (kg,  1) = 0._wp
          d (kg,  0) = 1._wp / (logpg(kg)-logpg(kg-1))
          d (kg, -1) = - d (kg,  0)
        case default
          call finish('set_vertinpc_rochon',                       &
                      'not implemented: ndv_rad = '//char2(ndv_rad))
        case (3)
          kg = 1
          d (kg,  1) = 1._wp / (logpg(kg+1)-logpg(kg))
          d (kg,  0) = - d (kg,  1)
          d (kg, -1) = 0._wp
          kg = nz
          d (kg,  1) = 0._wp
          d (kg,  0) = 1._wp / (logpg(kg)-logpg(kg-1))
          d (kg, -1) = - d (kg,  0)
          do kg = 2, nz - 1
            dd  =  logpg(kg+1)-logpg(kg-1)
            z   = (logpg(kg)  -logpg(kg-1)) / dd
            z2  = z**2
            zz2 = z - z2
            edd = 1._wp / (zz2 * dd)
            d (kg, 1) = z2                * edd
            d (kg, 0) = (1._wp - 2._wp*z) * edd
            d (kg,-1) = - (d (kg,1) + d (kg,0))
          end do
        end select
      endif

      do k = 1, size(logp)
        !-------------------------
        ! ignore zero coefficients
        !-------------------------
        if (pz (kstart(k), k, 1 ) == 0._wp .and. &
            pz (kstart(k), k, nt) == 0._wp )  kstart(k) = kstart(k) + 1
        if (pz (kend  (k), k, 1 ) == 0._wp .and. &
            pz (kend  (k), k, nt) == 0._wp )  kend  (k) = kend  (k) - 1
        !------------------------------------------------------------
        ! add coefficient space for explicit vertical differentiation
        !------------------------------------------------------------
        if (nt==1) then
          if (kstart(k) > 1) then
            kstart(k) = kstart(k) - 1
            pz (kstart(k),k,1) = 0._wp
          endif
          if (kend  (k) < size(logpg)) then
            kend  (k) = kend  (k) + 1
            pz (kend  (k),k,1) = 0._wp
          endif
        endif
        kn(k) = kend(k) - kstart(k) + 1
      end do
      !-------------------------------------------------
      ! check for sufficient space to store coefficients
      !-------------------------------------------------
      kmax = maxval (kn)
      if (kmax > mwv) call finish ('set_vertinpc_rochon',           &
                                   'mwvf too small: < '//char2(kmax))
      !----------------------------------------
      ! copy coefficients to 3dvar derived type
      !----------------------------------------
      do k = 1, nv        ! loop over all RTTOV input variables
        if (k <= n) then  ! loop over profile variables only
          !---------------------------------------
          ! store coefficients from RTTOV routines
          !---------------------------------------
          vic(k)%   wh(:kn(k)) = pz (kstart(k):kend(k),k,1 )
          if (nt==2) then
            vic(k)% wt(:kn(k)) = pz (kstart(k):kend(k),k,nt)
          else
            !----------------------------------
            ! explicit vertical differentiation
            !----------------------------------
            vic(k)% wt(:kn(k)) = 0._wp
            do k3 = 1, kn(k)
              kg = k3 + kstart(k)-1
              if (pz (kg,k,1) == 0._wp) cycle
              if (d  (kg, -1) /= 0._wp) &
                vic(k)% wt(k3-1) = vic(k)% wt(k3-1) &
                                 + d  (kg, -1) * pz (kg,k,1)
              if (d  (kg,  0) /= 0._wp) &
                vic(k)% wt(k3  ) = vic(k)% wt(k3  ) &
                                 + d  (kg,  0) * pz (kg,k,1)
              if (d  (kg, +1) /= 0._wp) &
                vic(k)% wt(k3+1) = vic(k)% wt(k3+1) &
                                 + d  (kg, +1) * pz (kg,k,1)
            end do
          endif
          !---------------
          ! adjust indices
          !---------------
          vic(k)%   wh(kn(k)+1:) = huge(1._wp) ! dont use
          vic(k)%   wt(kn(k)+1:) = huge(1._wp) ! dont use
          vic(k)%   g% i         = kstart(k)
          vic(k)%   g% n         = kn    (k)
        else
          call set_vertinpc (vic(k), logp(n), logpg, 2)
        endif
      end do
#endif
    end subroutine set_vertinpc_rochon
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine set_vertinpc1
!------------------------------------------------------------------------------
! subroutine set_vertinpc1 (vic, logp, logpg, nwv)
! !----------------------------------------
! ! set vertical interpolation coefficients
! !----------------------------------------
! type (t_vic) ,intent(out) :: vic   (:) ! interpolation coefficients to set
! real(wp)     ,intent(in)  :: logp  (:) ! vertical coordinates of observations
! real(wp)     ,intent(in)  :: logpg (:) ! vertical coordinates of grig
! integer      ,intent(in)  :: nwv       ! # of coefficients, 4:cubic interpol.


  subroutine set_vertinpc0 (vic, logp, logpg, nwv)
  !----------------------------------------
  ! set vertical interpolation coefficients
  !----------------------------------------
  type (t_vic) ,intent(out) :: vic       ! interpolation coefficients to set
  real(wp)     ,intent(in)  :: logp      ! vertical coordinates of observations
  real(wp)     ,intent(in)  :: logpg (:) ! vertical coordinates of grig
  integer      ,intent(in)  :: nwv       ! # of coefficients, 4:cubic interpol.

    integer  :: k, nz
#define USE_BISECTION
#ifdef  USE_BISECTION
    integer  :: klo, khi                          ! Use bisection, not lin.search
#endif
#if 1
    real(wp) :: dz, dz1, den, e2, f2, x(4), z(4)  ! For generic pressure levels
#else
    real(wp) :: dz, dz1, dz2, dz3                 ! Equidistant in log(p)
#endif
    real(wp) ,parameter :: eps = 1.e-10_wp        ! tolerance for upper bound

    nz = size (logpg)
    vic% wh = 0._wp
    vic% wt = 0._wp
!   vic% ix = 0
    select case (nwv)
    !---------------------
    ! linear interpolation
    !---------------------
    case (2)
      vic% g% n = 2
#ifndef EXTRAPOLATE_TOP
     if (logp < logpg(1)) then
        vic% g% i    = 1
        vic% wh(1)   = 1._wp
      else if (logp > logpg(nz)) then
#else
      if      (logp > logpg(nz)) then
#endif
        vic% g% i    = nz - 1
        vic% wh(2)   = 1._wp
      else
        do k = nz-1, 1, -1
          if (logp >= logpg(k)) then
            exit
          endif
        end do
        k = max (1,k)  ! for extrapolation
        vic% g% i    = k
        dz1          = logp       - logpg(k)
        dz           = logpg(k+1) - logpg(k)
        vic% wh(2)   = dz1   / dz
        vic% wh(1)   = 1._wp - vic% wh(2)
        vic% wt(2)   = 1._wp / dz
        vic% wt(1)   =       - vic% wt(2)
      endif
!     if (vic% ix == 0) call finish ('set_vertinpc0','logp out of bounds')
    !--------------------
    ! cubic interpolation
    !--------------------
    case (4)
      vic% g% n = 4
#ifndef EXTRAPOLATE_TOP
      if (logp < logpg(1) - eps) then
        vic% g% i    = 1
        vic% wh(1)   = 1._wp
      else if (logp > logpg(nz)) then
#else
      if      (logp > logpg(nz)) then
#endif
        vic% g% i    = nz-3
        vic% wh(4)   = 1._wp
      else if (logp < logpg(2)) then
        !------------------------------------------------
        ! linear interpolation between topmost two levels
        !------------------------------------------------
        vic% g% i    = 1
        dz1          = logp     - logpg(1)
        dz           = logpg(2) - logpg(1)
        vic% wh(2)   = dz1   / dz
        vic% wh(1)   = 1._wp - vic% wh(2)
        vic% wt(2)   = 1._wp / dz
        vic% wt(1)   =       - vic% wt(2)
      else if (logp > logpg(nz-1)) then
        vic% g% i    = nz-3
        dz1          = logp      - logpg(nz-1)
        dz           = logpg(nz) - logpg(nz-1)
        vic% wh(4)   = dz1   / dz
        vic% wh(3)   = 1._wp - vic% wh(4)
        vic% wt(4)   = 1._wp / dz
        vic% wt(3)   =       - vic% wt(4)
      else
#ifdef  USE_BISECTION
        khi = nz-1
        klo = 2
        do while (khi-klo > 1)
           k = (khi+klo)/2
           if (logpg(k) > logp) then
              khi = k
           else
              klo = k
           endif
        end do
        k = klo
#else
        do k = nz-2, 2, -1
          if (logpg(k) <= logp) then
            exit
          endif
        end do
#endif
        vic% g% i    = k - 1
#if 1
        !------------------------------------------------------
        ! Determine coefficients of a 3rd order polynomial f(p)
        ! satisfying the following conditions:
        !   f(p2) = y2, f'(p2) = (y3-y1)/(x3-x1),
        !   f(p3) = y3, f'(p3) = (y4-y2)/(x4-x2),
        !   (where x1 < x2 < x3 < x4),
        ! such that the interpolating function and its first
        ! derivative (for the calculation of T) are continuous.
        !------------------------------------------------------
        z(1:4)       = logpg(k-1:k+2)
        x(1:4)       = (logp - z(1:4))
        dz           =  z(3) - z(2)

        den          = (z(3)-z(1)) * dz**2
        vic% wh(1)   = - x(2) * x(3)**2          / den
        vic% wt(1)   = - x(3) * (2*x(2) + x(3))  / den

        den          = (z(4)-z(2)) * dz**2
        vic% wh(4)   =   x(3) * x(2)**2          / den
        vic% wt(4)   =   x(2) * (2*x(3) + x(2))  / den

        den          = dz**3
        e2           =   x(3)**2 * (3*x(2)-x(3)) / den
        f2           =   6 * x(2) * x(3)         / den
        vic% wh(2)   = - vic% wh(4) +      e2
        vic% wh(3)   = - vic% wh(1) + (1 - e2)
        vic% wt(2)   = - vic% wt(4) +      f2
        vic% wt(3)   = - vic% wt(1) -      f2
#else
        !-------------------------------------------
        ! Interpolation formula for equidistant case
        !-------------------------------------------
        dz           =  logpg(k+1) - logpg(k)
        dz1          = (logp       - logpg(k)) / dz
        dz2          = dz1 * dz1
        dz3          = dz1 * dz2

        vic% wh(1)   =        -0.5_wp * dz1 + 1.0_wp * dz2 - 0.5_wp * dz3
        vic% wh(2)   =  1._wp               - 2.5_wp * dz2 + 1.5_wp * dz3
        vic% wh(3)   =         0.5_wp * dz1 + 2.0_wp * dz2 - 1.5_wp * dz3
        vic% wh(4)   =                      - 0.5_wp * dz2 + 0.5_wp * dz3

        vic% wt(1)   =        -0.5_wp       + 2.0_wp * dz1 - 1.5_wp * dz2
        vic% wt(2)   =                      - 5.0_wp * dz1 + 4.5_wp * dz2
        vic% wt(3)   =         0.5_wp       + 4.0_wp * dz1 - 4.5_wp * dz2
        vic% wt(4)   =                      - 1.0_wp * dz1 + 1.5_wp * dz2

        vic% wt      = vic% wt / dz
#endif
      endif
    case default
      call finish ('set_vertinpc','invalid nwv')
    end select

#if 0
    ! Debug, remove these checks later (TODO)
    dz1 = sum (vic% wh) - 1
    dz  = sum (vic% wt)
    if (abs (dz) > 1.e-10 .or. abs (dz1) > 1.e-10) then
       write (0,*) "logp, vic% g:", logp, vic% g
       write (0,*) "wh:", vic% wh
       write (0,*) "wt:", vic% wt
       call finish ("set_vertinpc0", "inconsistent interpolation coeffs.")
    end if
#endif

  end subroutine set_vertinpc0
!==============================================================================
  subroutine random_gauss_covm (harvest, cov)
  real(wp)             ,intent(out)   :: harvest (:,:,:)
  type(t_cov)          ,intent(inout) :: cov
    integer  :: l, j
    real(wp) :: h (size (harvest,1), size (harvest,2))
    do l = 1, size (harvest, 3)
      harvest (:,:,l) = 0._wp
      do j = 1, cov% ct% n
        call random_gauss (h, cov% rndm_state (l,j))
        harvest (:,:,l) = harvest (:,:,l) + cov% ct% w(j) * h
      end do
    end do
  end subroutine random_gauss_covm
!==============================================================================
  subroutine read_horizontal_covariances (filename, rnd_time, csr_symm, &
                                                              csr_conv  )
    !-----------------------------------------------------
    ! Read horizontal wavelet covariances from NetCDF file
    ! and store in the module variable covm.
    !-----------------------------------------------------
    character(len=*),  intent(in) :: filename
    real(wp),          intent(in) :: rnd_time ! random B correlation time scale
    logical, optional, intent(in) :: csr_symm ! Set symmetry flag for CSR
    character(len=*),  intent(in) :: csr_conv ! Convert sparse representation
    optional                      :: csr_conv
    !----------------
    ! Local variables
    !----------------
    integer,  parameter  :: MAXMAT = 16
    integer              :: j, n_fields, nf_file
    integer              :: i, id, k, l, m, pe, p_in, nz, n_data, i_h
    integer              :: lmax
    integer              :: map(nvar)
    logical              :: ok, perm, symm
    integer              :: u
    character(len=128)   :: field_list, description
    character(len=255)   :: profile_list
    character(len=16)    :: fields(MAXMAT)
    character(len=8)     :: conv
    type(t_cov_meta)     :: meta
    type(t_cov_wavelet),  pointer :: W => NULL ()
    type(t_matrix_block), pointer :: B => NULL ()
#if 0 /* new faster(?) version */
    type(t_cov_wavelet)  :: rcv_hc2     ! Dummy receive buffer for bcast
#endif
    !--------------------------------------------

    if (dace% lpio) then
       write (6,'()')
       write (6,'(a)') repeat('-',79)
       write (6,'()')
       write (6,'(a)') '  read_horizontal_covariances'
       write (6,'()')
    end if
    !--------------
    ! Sanity checks
    !--------------
    symm = .false.; if (present (csr_symm)) symm = csr_symm
    perm = .false.
    conv = "";      if (present (csr_conv)) conv = toupper (csr_conv)
    select case (conv)
    case ("CSR", "")
       ! OK
    case ("CSRPERM")
       perm = .true.
    case ("CSC")
       if (symm) then
          ! Not an error, just warn the user
          call message ("read_horizontal_covariances"," CSC and csr_symm")
       end if
    case ("JAD", "JAGGED", "JDS")
       if (.not. symm) then
          call message ("read_horizontal_covariances"," JAD but not csr_symm")
          symm = .true.
       end if
    case default
       call finish ("read_horizontal_covariances","csr_conv = "//conv)
    end select

    if (dace% lpio) then
       call read_cov_netcdf_meta (filename, meta, field_list, description)
       call read_profile_list_netcdf (filename, profile_list)
       write (6,*) "Reading NetCDF file   : ", trim (filename)
       write (6,*) "description           : ", trim (description)
       write (6,*) "gridinfo              : ", trim (meta% gridinfo)
       write (6,*) "gridtype (WMO table 6):" , meta% gridtype
       write (6,*) 'dimensions            :' , meta% dims(:)
       do k = 1, size (meta% dims)
          write (6,'(a,i2,a,i12,a,a)') '        wavelet basis', k, ":", &
               meta% wv_basis(k), " = ", trim (wv_name(meta% wv_basis(k)))
       end do
       write (6,*) "B representation type : ", meta% representation
       write (6,*) "field_list            : ", trim (field_list)
       if (profile_list == "") then
          write (6,*) "meridional profiles   : (none)"
       else
          write (6,*) "meridional profiles   : ", trim (profile_list)
       end if
       write (6,*)
       !-------------------------
       ! Validate global metadata
       !-------------------------
       if (meta% representation /= "L") then
          call finish ("read_horizontal_covariances",&
                       "cannot handle representation " // meta% representation)
       end if
       if (meta% dims(3) /= 1) then
          write (0,*) "Bad grid: dims(3) =", meta% dims(3)
          call finish ("read_horizontal_covariances","cannot continue")
       end if
       !-------------------------
       ! Overwrite other metadata
       !-------------------------
       covm% gridtype = meta% gridtype
       if (covm% nx /= meta% dims(1)) then
         if (covm% nx <= 0) then
           covm% nx = meta% dims(1)
         else
           write(0,*)'read_horizontal_covariances: nx /= meta% dims(1):',&
                     covm% nx,meta% dims(1)
           call finish('read_horizontal_covariances','nx /= meta% dims(1):')
         endif
       endif
       if (covm% ny /= meta% dims(2)) then
         if (covm% ny <= 0) then
           covm% ny = meta% dims(2)
         else
           write(0,*)'read_horizontal_covariances: ny /= meta% dims(2):',&
                     covm% ny,meta% dims(2)
           call finish('read_horizontal_covariances','ny /= meta% dims(2):')
         endif
       endif
    end if

    call p_bcast (covm% gridtype, dace% pio)
    call p_bcast (covm% nx,       dace% pio)
    call p_bcast (covm% ny,       dace% pio)
    call p_bcast (field_list,     dace% pio)
    call p_bcast (meta,           dace% pio)
    call p_bcast (profile_list,   dace% pio)

    call split (fields, field_list, n_fields)
    map = 0     ! 'map' holds an auxiliary index to the field to be read
    do j = 1, n_fields
       select case (fields(j))
       case ("geof") ; where (vars == INT_H)   map = j
       case ("rh")   ; where (vars == INT_RH)  map = j
       case ("chi")  ; where (vars == INT_CHI) map = j
       case ("psi_u"); where (vars == INT_PSI) map = j
       case default
          write (0, '(a)') "Not implemented: field: " // trim (fields(j))
          call finish ("read_horizontal_covariances","cannot continue")
       end select
    end do
    nf_file = count (map > 0)
    nz      = covm% nz
    !--------------------------------------------------------------------
    ! Correlations for geopotential and rel. humidity must be present.
    ! Correlations for streamfunction and velicity potential will be
    ! taken the same as for geopotential when not explicitly given.
    !--------------------------------------------------------------------
    ok = .true.
    if (nf_file < nvar) then
       ! Geopotential is mandatory
       if (.not. any (vars == INT_H .and. map > 0)) then
          ok = .false.
          if (dace% lpio) write (0,*) "geopotential is mandatory!"
       end if
       i_h = sum (map, vars == INT_H)
       ! Rel. humidity is mandatory
       if (.not. any (vars == INT_RH .and. map > 0)) then
          ok = .false.
          if (dace% lpio) write (0,*) "(relative) humidity is mandatory!"
       end if
       ! (Unbalanced) Streamfunction: fallback to geopotential
       if (.not. any (vars == INT_PSI .and. map > 0)) then
          if (dace% lpio) then
             write (6,*) "Streamfunction correlation: fallback to geopotential"
          end if
          where (vars == INT_PSI) map = i_h
       end if
       ! Velocity potential: fallback to geopotential
       if (.not. any (vars == INT_CHI .and. map > 0)) then
          if (dace% lpio) then
             write (6,*) "Velocity potential:         fallback to geopotential"
          end if
          where (vars == INT_CHI) map = i_h
       end if
    end if
    if (.not. ok) then
       call p_barrier ()
       call finish ("read_horizontal_covariances","missing correlations")
    end if
    if (dace% lpio) write (6,*)
    if (count (map > 0) < nvar) then
       if (dace% lpio) write (0,*) "nvar, map =", nvar, map(:)
       call p_barrier ()
       call finish ("read_horizontal_covariances","internal error")
    end if

    !------------------------------------------------------------------
    ! Set up data partitioning table, keeping adjacent model levels of
    ! the same variable together.
    ! The partitioning works best when nz*nvar is divisible by dace% npe
    !------------------------------------------------------------------
    n_data = nz * nvar
    covm% hc2_totalsize = n_data
    covm% hc2_part_size = (dace% pe + 1) * n_data / dace% npe &
                        -  dace% pe      * n_data / dace% npe

#if 0
    if (dace% lpio) then
       write (6,'(a,5i5)') &
            " Partitioning data: nz, nvar, n_data, nprocs, maxl =", &
            nz, nvar, n_data, dace% npe, covm% hc2_part_size
       write(6,*)
       write(6,'(a)') "     m    id   lvl  :    pe     l     i"
       write(6,'()')
    end if
#endif

    call construct (covm% ct, rnd_time, covm% time, verbose=.true.)
    allocate       (covm% rndm_state (covm% hc2_part_size,covm% ct% n))
    allocate       (covm% hc2_part   (n_data))
    m    = 0
    lmax = 0
    do i = 1, nvar
       do k = 1, nz
          m  = m + 1                        ! Total index: m=1..n
          pe = (m * dace% npe - 1)/ n_data  ! Map m=1..n onto pe=0..(dace% npe-1)
          l  =  m - pe * n_data / dace% npe ! Reduce index to local array index
          id =  vars(i)
          covm% hc2_part(m)% pe = pe
          covm% hc2_part(m)% k  = k
          covm% hc2_part(m)% l  = l
          covm% hc2_part(m)% id = id
          covm% hc2_part(m)% i  = i
          if (pe == dace% pe) lmax  = l

#if 0
          if (dace% lpio) write(6,'(3i6,"  :",3i6)') m,id,k,pe,l,i
#endif
          u = -1
          do j = 1, covm% ct% n
            if (dace% pe==pe) then
              call construct_stream (covm% rndm_state (l,j),       &
                                     hour= covm% ct% ht(j), uniq= u)
            else
              call construct_stream (hour= covm% ct% ht(j), uniq= u)
            endif
            u = 0
          end do
       end do
    end do
    if (dace% lpio) write (6,*)
    !----------------------------------
    ! Consistency check on partitioning
    !----------------------------------
    if (lmax /= covm% hc2_part_size) then
       write (0,*) lmax, " /=", covm% hc2_part_size
       call finish ("read_horizontal_covariances","hc2_part_size /= lmax")
    end if
    lmax = count(covm% hc2_part% pe == dace% pe)
    if (lmax /= covm% hc2_part_size) then
       write (0,*) lmax , " /=", covm% hc2_part_size
       call finish ("read_horizontal_covariances","hc2_part_size /= count")
    end if

    !--------------------------
    ! read and scatter matrices
    !--------------------------
    covm% hc2(:)% valid = .false.
    do i = 1, nvar
       id = vars(i)
       outer: do m = 1, n_data
          if (covm% hc2_part(m)% id /= id) cycle outer
          !------------------------------
          ! read data on chosen processor
          !------------------------------
          p_in = covm% hc2_part(m)% pe
          if (dace% pe == p_in) then
!            write (6,'("pe=",i3,"  id=",i3,a)') dace% pe, id, &
!                 " : reading field:          " // trim (fields(map(i)))
!            call flush (6)
             W => covm% hc2(i)
             call read_cov_netcdf (W, filename, field=trim (fields(map(i))))
             W% id = id
             B => W% b
             if (W% field /= fields(map(i))) then
                write (0,*) "???  Found: ", trim (W% field)
             end if
!            write (6,'("pe=",i3,"  id=",i3,a)') dace% pe, id,        &
!                 " : storage_representation: " // trim (srep(B% repr))
!            write (6,'("pe=",i3,"  id=",i3,a,i11,a,f8.2,a)') dace% pe,id,&
!                 " : nonzero elements      :" , B% nonzero,              &
!                 " =", B% nonzero/real (B% m), " coeffs./gridpoint"

             select case (conv)
             case ("", "CSR", "CSRPERM")
                if (conv /= "") then
!                  write (6,'("pe=",i3,"  id=",i3,a)') dace% pe, id, &
!                       " : converting to CSR"
                   call convert (B, repr=CSR)
                end if
                if (B% repr == CSR .and. perm) then
!                  write (6,'("pe=",i3,"  id=",i3,a)') dace% pe, id, &
!                       " : setting up row permutation for CSR"
                   call set_perm (B)
                end if
             case ("CSC")
!               write (6,'("pe=",i3,"  id=",i3,a)') dace% pe, id, &
!                    " : converting to CSC"
                call convert (B, repr=CSC)
             case ("JAD", "JAGGED", "JDS")
!               write (6,'("pe=",i3,"  id=",i3,a)') dace% pe, id, &
!                    " : converting to JAD"
                call convert (B, repr=JAD)
             end select

             if (symm) then
!               write (6,'("pe=",i3,"  id=",i3,a)') dace% pe, id, &
!                    " : assuming symmetry of sqrt(B)"
                B% qual = SYMMETRIC         ! Exploit symmetry of cov.matrix
             end if
             call flush (6)
             nullify (W, B)
          endif
#if 0 /* new faster(?) version */
          !--------------------------
          ! scatter using plain bcast
          !--------------------------
          if (any (      covm% hc2_part(:)% pe == dace% pe  &
                   .and. covm% hc2_part(:)% id == id       )) then
             !-----------------------------------------------------------
             ! This PE belongs to the group handling the current variable
             !-----------------------------------------------------------
             call p_bcast (covm% hc2(i), p_in)
          else
             !--------------------------
             ! Dummy receive and discard
             !--------------------------
             call p_bcast  (rcv_hc2, p_in)
             call destruct (rcv_hc2)
          end if
#else /* old slower(?) version */
          !------------------------------
          ! scatter (explicit send/recv)
          !         (bcast may be better)
          !------------------------------
          inner: do pe = 0, dace% npe-1
             if (pe == p_in) cycle inner
             if (any (      covm% hc2_part(:)% pe == pe  &
                      .and. covm% hc2_part(:)% id == id )) then
                !----------------------------------------------------------
                ! This var is also needed on *another* processor pe /= p_in
                !----------------------------------------------------------
                if (dace% pe == p_in) then
                   call p_send (covm% hc2(i), pe  )
                else if (dace% pe == pe) then
                   call p_recv (covm% hc2(i), p_in)
                endif
             endif
          end do inner
#endif
          !------------------------------------------------
          ! Move on to matrix for next var (outermost loop)
          !------------------------------------------------
          exit outer
       end do outer
    end do

#if 0
    !---------
    ! printout
    !---------
    call p_barrier
    if (dace% lpio) then
       write(6,*)
       write(6,*) "Verifying distribution of matrices onto processors:"
       write(6,*)
    end if
    do pe = 0, dace% npe-1
       call p_barrier
       if (pe /= dace% pe) cycle
       do i = 1, nvar
          if (.not. covm% hc2(i)% valid) cycle
          B => covm%hc2(i)% B
          write(6,*) "pe=", dace% pe, "id=", covm%hc2(i)% id, ":", B% nonzero, &
               real (minval (B% packed)), "..", real (maxval (B% packed))
       end do
    end do
#endif

    !----------------------------------------------------------
    ! Read meridional profiles and distribute to all processors
    !----------------------------------------------------------
    if (profile_list /= "") then
       call read_meridional_profiles (filename, profile_list)
    end if

    call p_barrier
    if (dace% lpio) then
       write (6,'()')
       write (6,*) "read_horizontal_covariances: done."
       write (6,'()')
       write (6,'(a)') repeat('-',79)
       write (6,'()')
    end if

  end subroutine read_horizontal_covariances
!==============================================================================
  subroutine read_meridional_profiles (filename, profile_list)
    !------------------------------------------------------
    ! Read meridional profiles from NetCDF covariances file
    ! and store in the module variable covm.
    !------------------------------------------------------
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: profile_list
    !----------------
    ! Local variables
    !----------------
    integer,  parameter :: MAXPROF = 32
    integer             :: ny, n_profiles, j, k, n
!   integer             :: pe
    character(len=16)   :: profiles(MAXPROF)
    character(len=128)  :: units, description
    !--------------------------------------------------------
    ! Normalization of horizontal wind components (variances)
    !--------------------------------------------------------
    type wind_var_t
       real(wp), pointer :: uu(:) => NULL (),    &! Sigma^2(uu)
                            vv(:) => NULL (),    &! Sigma^2(vv)
                            uv(:) => NULL ()      ! Sigma^2(uv)
    end type wind_var_t
    !-----------------------------------------------------
    ! Variances of balanced, unbalanced and divergent wind
    !-----------------------------------------------------
    type(wind_var_t) :: wind_var(3)

    if (profile_list == "") return

    if (dace% lpio) write (6,*)

    ! Decompose
    call split (profiles(:), profile_list, n_profiles)

    ny = covm% ny
    call get_profile ("mu",            covm% mu)
    call get_profile ("stddev_geof",   covm% sdev_h)
    call get_profile ("stddev_rh",     covm% sdev_rh)
    call get_profile ("stddev_chi",    covm% sdev_chi)
    call get_profile ("stddev_psi",    covm% sdev_psi)
    call get_profile ("stddev_psi_u",  covm% sdev_psi_u)
    call get_profile ("stddev_u",      covm% sdev_u)
    call get_profile ("stddev_v",      covm% sdev_v)
    !-------------------------------------------------------------------
    ! Variances of geostrophic, unbalanced and divergent wind components
    !-------------------------------------------------------------------
    call get_profile ("var_uu_g",     wind_var(1)% uu)
    call get_profile ("var_uu_u",     wind_var(2)% uu)
    call get_profile ("var_uu_div",   wind_var(3)% uu)
    call get_profile ("var_vv_g",     wind_var(1)% vv)
    call get_profile ("var_vv_u",     wind_var(2)% vv)
    call get_profile ("var_vv_div",   wind_var(3)% vv)
    call get_profile ("var_uv_g",     wind_var(1)% uv)
    call get_profile ("var_uv_u",     wind_var(2)% uv)
    call get_profile ("var_uv_div",   wind_var(3)% uv)

    if (associated (covm% mu)) then
      !----------------------------------------------
      ! (stddev_psi_u)^2 == (stddev_psi)^2 * (1-mu^2)
      !----------------------------------------------
      if (.not. associated (covm% sdev_psi_u) .and. associated (covm% sdev_psi)) then
         if (dace% lpio) then
            write (6,*) "stddev_psi_u not provided, deriving from stddev_psi"
         end if
         allocate (covm% sdev_psi_u(ny))
         covm% sdev_psi_u(:) = covm% sdev_psi(:) * sqrt (1 - covm% mu(:)**2)
      endif
      if (.not. associated (covm% sdev_psi) .and. associated (covm% sdev_psi_u)) then
         if (dace% lpio) then
            write (6,*) "stddev_psi not provided, deriving from stddev_psi_u"
         end if
         allocate (covm% sdev_psi(ny))
         covm% sdev_psi(:) = covm% sdev_psi_u(:) / sqrt (1 - covm% mu(:)**2)
      endif
    end if
    !-----------------------------------------------
    ! Check completeness of essential wind variances
    !-----------------------------------------------
    n = 0
    k = 0
    do j = 1, 3
       if (associated (wind_var(j)% uu)) n = n + 1
       if (associated (wind_var(j)% vv)) n = n + 1
       if (associated (wind_var(j)% uv)) k = k + 1
    end do
    !--------------------------------------------------------------
    ! (Re-)Derive isotropic expression for wind error normalization
    !--------------------------------------------------------------
    select case (n)
    case (0)
       if (associated (covm% sdev_u) .and. associated (covm% sdev_v)) then
          if (dace% lpio) then
             write (6,*)
             write (6,*) "read_meridional_profiles: ", &
                  "Warning: found only anisotropic NMC wind variances"
             write (6,*) "read_meridional_profiles: ", &
                  "Replacing NMC wind variances by isotropic mean"
          end if
          covm% sdev_u(:) = sqrt ((covm% sdev_u(:)**2 + covm% sdev_v(:)**2)/2)
          covm% sdev_v(:) = covm% sdev_u(:)
       else
          if (dace% lpio) then
             write (6,*)
             write (6,*) "read_meridional_profiles: ", &
                  "Warning: could not find NMC wind variances!"
          end if
       end if
    case (6)
       if (dace% lpio) then
          write (6,*)
          write (6,*) "read_meridional_profiles: ", &
               "Replacing NMC by recalculated wind variances"
       end if
       if (.not. associated (covm% sdev_u)) allocate (covm% sdev_u(ny))
       if (.not. associated (covm% sdev_v)) allocate (covm% sdev_v(ny))
       covm% sdev_u(:) = sqrt (( wind_var(1)% uu(:) + wind_var(1)% vv(:) &
                               + wind_var(2)% uu(:) + wind_var(2)% vv(:) &
                               + wind_var(3)% uu(:) + wind_var(3)% vv(:) )/2)
       covm% sdev_v(:) = covm% sdev_u(:)
       !----------------------------------------------------------
       ! Check presence of currently unused off-diagonal variances
       !----------------------------------------------------------
       select case (k)
       case (0)
          if (dace% lpio) then
             write (6,*) "read_meridional_profiles: ", &
                  "info: missing off-diagonal variances"
          end if
       case (3)
          ! OK
       case default
          if (dace% lpio) then
             write (6,*) "read_meridional_profiles: ", &
                  "Warning: unexpected count of off-diagonal variances", k
          end if
       end select
    case default
       write (0,*) "Count of wind diagonal variances (expected: 6):", n
       call finish ("read_meridional_profiles", "bad or corrupt input file?")
    end select
    !---------------------
    ! Clean up temporaries
    !---------------------
    do j = 1, 3
       if (associated (wind_var(j)% uu)) deallocate (wind_var(j)% uu)
       if (associated (wind_var(j)% vv)) deallocate (wind_var(j)% vv)
       if (associated (wind_var(j)% uv)) deallocate (wind_var(j)% uv)
       nullify (wind_var(j)% uu, wind_var(j)% vv, wind_var(j)% uv)
    end do

    if (dace% lpio) write (6,*)

#if 0
    ! Print some debugging info
    call p_barrier
    if (dace% lpio) then
       write(6,*)
       write(6,*) "Verifying distribution of profiles onto processors:"
       write(6,*)
    end if
    do pe = 0, dace% npe-1
       call p_barrier
       if (pe /= dace% pe) cycle
       write(6,*) "pe=", pe, "associated: ", associated (covm% mu), &
            associated (covm% sdev_h)  , associated (covm% sdev_rh), &
            associated (covm% sdev_chi), associated (covm% sdev_psi), &
            associated (covm% sdev_u),   associated (covm% sdev_v)
       if (associated (covm% mu)) write(6,*) "pe=", pe, &
            ": mu=", real (covm% mu(1)), "..", real (covm% mu(ny))
    end do
    call p_barrier
    if (dace% lpio) write(6,*)
#endif

    contains

      subroutine get_profile (varname, var)
        character(len=*), intent(in) :: varname
        real(wp),         pointer    :: var(:)

        if (any (profiles(:) == varname)) then
           allocate (var(ny))
           if (dace% lpio) then
              write (6,*) "pe=", dace% pe, ": reading ", varname
              call read_profile_netcdf (filename = filename,  &
                                        varname  = varname,   &
                                        var      = var(:),    &
                                        units    = units,     &
                                        desc     = description)
           end if
           call p_bcast (var, dace% pio)
        else
           if (dace% lpio) then
              write (6,*) varname, " not found in ", trim (filename)
           end if
           nullify (var)
        end if

      end subroutine get_profile
 end subroutine read_meridional_profiles
!==============================================================================
end module mo_t_bg_err_op
