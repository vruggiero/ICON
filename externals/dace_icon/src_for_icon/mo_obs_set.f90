!
!+ Defines set of data structures accompanying the observational data type
!
MODULE mo_obs_set
!
! Description:
!   Defines set of data structures accompanying the observational data
!   type T_OBS_SET. In addition to the basic observational data type T_OBS
!   observation error covariance matrix and linearised observation
!   operators are provided.
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
! V1_5         2009/05/25 Andreas Rhodin
!  subroutine destruct_obs_set: check association status of set% o
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Andreas Rhodin
!  technical changes for revised minimisation routine
! V1_15        2011/12/06 Andreas Rhodin
!  option to remove observations in the assimilation step
! V1_19        2012-04-16 Andreas Rhodin
!  define specific 'destruct' for arrays of derived type
! V1_23        2013-03-26 Andreas Rhodin
!  keep B in interpolation space for diagnostic printout in cofRTOVP.nc
! V1_27        2013-11-08 Andreas Rhodin
!  implement Variational Bias Correction (VarBC)
! V1_37        2014-12-23 Andreas Rhodin
!  t_obs_set: new components  Jo_spot, si
! V1_47        2016-06-06 Andreas Rhodin
!  option to write B in interpolation space to cofRTOVP.nc
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Andreas Rhodin  DWD  2005-2007
!========================================================================

!============
!Modules used
!============
use mo_kind,       only: wp               ! working precision kind parameter
use mo_exception,  only: finish           ! exit on error condition
use mo_dec_matrix, only: t_matrix,       &! block decomposed matrix
                         t_matrix_block, &! matrix blocks
                         t_vector,       &! segmented vector
                         t_vector_segm,  &! vector segment
                         t_ivector,      &! integer vector
                         t_dec_info,     &! decomposition info
                         t_info_block,   &! decomposition meta data
                         mp,             &! matrix precision
                         inverse,        &! invert matrix
                         operator(*),    &! multiply      matrices or vectors
                         operator(-),    &! difference of matrices or vectors
                         operator(+),    &! sum        of matrices or vectors
                         assignment(=),  &! assign vectors
                         destruct,       &! deallocate matrices, vectors
                         deallocate,     &! deallocate matrix block
                         allocate_block, &! allocate matrix block
                         FULL             ! matrix representation value
use mo_t_obs,      only: t_obs,          &! observation data type
                         t_spot,         &! observation data type
                         destruct,       &! deallocate t_obs
                         debug_spot,     &! debug selected spot(s)
                         ldeb_spot        ! whether spots are to be debugged
use mo_mpi_dace,   only: dace,           &! MPI group info
                         p_bcast          ! generic MPI broadcast
implicit none

!================
! Public entities
!================
private
public :: t_obs_set   ! observation data type set, global
public :: t_obs_block ! observation data type set, box level
public :: t_vqc       ! quadratic J_o data type
public :: t_ooj       ! observation operator Jakobian data type
public :: t_bg        ! background in interpolation and observation space
public :: t_vbc       ! variational bias correction data type
public :: t_vbc_meta  ! variational bias correction meta data
public :: obs_block   ! extraxt t_obs_block from t_obs_set
public :: apply_h     ! H*x, x*H, inv(H)*x, ..
public :: destruct    ! deallocate components
public :: debug_spot  ! debug selected spot(s)

!======================
! Data type definitions
!======================

  !-----------------------------------------------------------------
  ! quadratic approximation of VQC (Variational Quality Control) J_o
  !-----------------------------------------------------------------
  type t_vqc
    type(t_vector)            :: y         ! point of approximation
    real(wp)                  :: J         ! obs. cost function at yj
    type(t_vector)            :: Jo        ! obs. cost function
    type(t_vector)            :: Jo_spot   ! obs. cost function
    type(t_vector)            :: dJ        ! gradient of    Joj at yj
    type(t_matrix)            :: R         ! Hessian  of    Joj at yj
    type(t_vector)            :: w         ! VarQC weight       at yj
    type(t_vector)            :: active    ! mask for active variables (0/1)
  end type t_vqc

  !---------------------------------
  ! variational bias correction data
  !---------------------------------
  type t_vbc_meta
    !---------------------------------------------------
    ! meta data: one entry per satellite/instrument
    ! describes organisation of segments in t_vbc below:
    ! p (n_tot), n_tot = (n_fov + n_used) * n_chan
    !---------------------------------------------------
    integer                   :: codetype    = 0  ! code type (RAD / AIREP)
    character(len=8)          :: platform    = '' ! satellite / aircraft
    integer                   :: n_fov       = 0  ! number of FOVs used
    integer                   :: n_pred      = 0  ! number of predictors
    integer                   :: n_coeff     = 0  ! number of coefficients
    integer                   :: n_used      = 0  ! number of predictors used
    integer                   :: n_chan      = 0  ! number of channels used
    integer                   :: n_tot       = 0  ! (n_used + n_fov) * n_chan
    integer                   :: is          = 0  ! segment index in t_vbc
    integer                   :: o           = 0  ! offset in segment
  end type t_vbc_meta

  type t_vbc
    !----------------------------------------------------------------
    ! variational bias correction data
    ! one 'vector segment' per  satellite/instrument
    ! B is blockdiagonal
    ! Information on z is globally avaiylable (communicate via p_sum)
    ! Information on p is globally avaiylable
    ! B may be held locally
    !----------------------------------------------------------------
    type(t_vbc_meta) ,pointer :: meta(:) => Null() ! meta data
    type(t_ivector)           :: ic                ! channel index
    type(t_vector)            :: bc                ! bias correction
    type(t_vector)            :: z                 ! gradient in pred.space
    type(t_vector)            :: cf                ! coefficient update
    type(t_dec_info) ,pointer :: bi      => NULL() ! VBC space description
    real(wp)                  :: J       = 0._wp   ! VarBC cost function
  end type t_vbc

  !------------------------------
  ! observation operator Jakobian
  !------------------------------
  type t_ooj
    type(t_matrix)            :: H            ! linearized observation operator
    type(t_vector)            :: x            ! point of linearisation
    type(t_vector)            :: y            ! point of linearisation
  end type t_ooj

  !--------------------------------------------------
  ! background in interpolation and observation space
  !--------------------------------------------------
  type t_bg
    type(t_vector)            :: xb           ! background in interp. space
    type(t_vector)            :: yb           ! H  (xb)    in observ. space
    type(t_vector)            :: xl_xb        ! linearisation point - bg
    type(t_vector)            :: ul_ub        ! H * xl_xb
    type(t_vector)            :: zx           ! xl_xb = B * zx
    type(t_matrix)            :: Bii          ! B in int.space (optional)
    real(wp)                  :: J            ! bg. cost function at xl
  end type t_bg

  !------------------------
  ! set of observation data
  !------------------------
  type t_obs_set
    !---------------------------
    ! observations and meta data
    !---------------------------
    type(t_obs)      ,pointer :: o (:) => NULL() ! observation data type
    !----------------------
    ! description of spaces
    !----------------------
    type(t_dec_info) ,pointer :: si => NULL() ! report/fov    space description
    type(t_dec_info) ,pointer :: oi => NULL() ! observation   space description
    type(t_dec_info) ,pointer :: ii => NULL() ! interpolation space description
    type(t_dec_info) ,pointer :: di => NULL() ! dummy type    space description
    !---------------------------------
    ! linearised observation operators
    !---------------------------------
    type (t_ooj)              :: l
    !--------------------------------------------------
    ! background in interpolation and observation space
    !--------------------------------------------------
    type (t_bg)               :: b
    !-------------------------------------
    ! observation and observational errors
    !-------------------------------------
    type(t_vector)            :: obc          ! bias corrected observation
    type(t_matrix)            :: R            ! nominal observational error
    !-----------------------------------
    ! quadratic approcimation of VQC J_o
    !-----------------------------------
    type (t_vqc)              :: vc
    !----------------------------
    ! Variational bias correction
    !----------------------------
    type(t_vbc)               :: vbc
  end type t_obs_set

  !--------------------------
  ! block of observation data
  !--------------------------
  type t_obs_block
    type(t_obs)          ,pointer :: o  ! observation data type
    type(t_info_block)   ,pointer :: oi ! observation   space description
    type(t_info_block)   ,pointer :: ii ! interpolation space description
    type(t_info_block)   ,pointer :: di ! dummy type    space description
    type(t_matrix_block) ,pointer :: H  ! linearised observation operator
    type(t_vector_segm)  ,pointer :: xi ! linearisation point (argument)
    type(t_vector_segm)  ,pointer :: yi ! linearisation point (result)
    type(t_matrix_block) ,pointer :: R  ! observation error covariance matrix
    type(t_matrix_block) ,pointer :: Bi ! B matrix, interpolation space
  end type t_obs_block

  !----------
  ! interface
  !----------
  interface destruct
    module procedure destruct_obs_set
    module procedure destruct_obs_sets
    module procedure destruct_vc
    module procedure destruct_bg
    module procedure destruct_H
    module procedure destruct_vbc
  end interface destruct

  interface debug_spot
    module procedure debug_spot_obs_set
  end interface debug_spot

contains
!------------------------------------------------------------------------------
  subroutine obs_block (block, obs, i)
  type(t_obs_block) ,intent(out) :: block
  type(t_obs_set)   ,intent(in)  :: obs
  integer           ,intent(in)  :: i
  !-----------------------------------------------------
  ! extract     observation data type set on block level
  ! from global observation data type set
  !-----------------------------------------------------
    block% o => obs% o(i)
    nullify (block% R)
    nullify (block% H)
    nullify (block% oi)
    nullify (block% ii)
    nullify (block% di)
    nullify (block% Bi)
    if (associated (obs%    R%   b)) block% R  => obs%    R % b(i,i)
    if (associated (obs% l% H%   b)) block% H  => obs% l% H % b(i,i)
    if (associated (obs% l% x%   s)) block% xi => obs% l% x % s(i)
    if (associated (obs% l% y%   s)) block% yi => obs% l% y % s(i)
    if (associated (obs%    oi    )) block% oi => obs%   oi % b(i)
    if (associated (obs%    ii    )) block% ii => obs%   ii % b(i)
    if (associated (obs%    di    )) block% di => obs%   di % b(i)
    if (associated (obs% b% Bii% b)) block% Bi => obs% b%Bii% b(i,i)
  end subroutine obs_block
!------------------------------------------------------------------------------
  subroutine destruct_obs_sets (sets)
  type(t_obs_set) ,intent(inout) :: sets (:)
    integer :: i
    do i = 1, size (sets)
      call destruct (sets (i))
    end do
  end subroutine destruct_obs_sets
!------------------------------------------------------------------------------
  subroutine destruct_obs_set (set)
  type(t_obs_set) ,intent(inout) :: set
  !-----------------------
  ! deallocates components
  !-----------------------
    if (associated(set% o)) then
      call destruct (set% o)
      deallocate    (set% o)
    endif
    call destruct (set% R)
    call destruct (set% obc)
    call destruct (set% vc)
    call destruct (set% b)
    call destruct (set% l)
    if (associated (set% si)) then
      call destruct (set% si)
      deallocate (set% si)
    endif
    if (associated (set% oi)) then
      call destruct (set% oi)
      deallocate (set% oi)
    endif
    if (associated (set% ii)) then
      call destruct (set% ii)
      deallocate (set% ii)
    endif
    if (associated (set% di)) then
      call destruct (set% di)
      deallocate (set% di)
    endif
  end subroutine destruct_obs_set
!------------------------------------------------------------------------------
  subroutine destruct_vc (vc)
  type(t_vqc) ,intent(inout) :: vc
  !-----------------------
  ! deallocates components
  !-----------------------
    call destruct (vc% y)
    call destruct (vc% Jo)
    call destruct (vc% dJ)
    call destruct (vc% R)
    call destruct (vc% w)
    call destruct (vc% active)
  end subroutine destruct_vc
!------------------------------------------------------------------------------
  subroutine destruct_H (l)
  type(t_ooj) ,intent(inout) :: l
  !-----------------------
  ! deallocates components
  !-----------------------
    call destruct (l% H)
    call destruct (l% x)
    call destruct (l% y)
  end subroutine destruct_H
!------------------------------------------------------------------------------
  subroutine destruct_bg (bg)
  type(t_bg) ,intent(inout) :: bg
    call destruct (bg% xb)
    call destruct (bg% yb)
    call destruct (bg% xl_xb)
    call destruct (bg% ul_ub)
    call destruct (bg% zx)
    call destruct (bg% Bii)
  end subroutine destruct_bg
!------------------------------------------------------------------------------
  subroutine destruct_vbc (vbc)
  type (t_vbc) ,intent(inout) :: vbc
    if (associated   (vbc% meta)) deallocate (vbc% meta)
    call destruct    (vbc% ic)
    call destruct    (vbc% bc)
    call destruct    (vbc% z)
    call destruct    (vbc% cf)
    if (associated   (vbc% bi)) then
      call destruct  (vbc% bi)
      deallocate     (vbc% bi)
    endif
  end subroutine destruct_vbc
!------------------------------------------------------------------------------
  subroutine apply_H (obs, psi, pso, mode)
  type (t_obs_set),intent(in)    :: obs    ! observations
  type (t_vector) ,intent(inout) :: psi    ! interpolated values
  type (t_vector) ,intent(inout) :: pso    ! observables
  character       ,intent(in)    :: mode   ! 'o','l','a','i'
  !---------------------------------------------------------------
  ! apply the linearised observation operator
  !
  ! mode: 'o' account for background values used for linearisation
  !       'l' do not account for background values
  !       'a' adjoint mode
  !       'i' inverse mode
  !---------------------------------------------------------------
    integer                    :: l     ! observation box index
    integer                    :: j     ! observation index
    integer                    :: k,m   ! indices
    type (t_spot) ,pointer     :: spt   ! pointer to meta data element
    real(wp)      ,pointer     :: xi(:) ! interpolated values
    real(wp)      ,pointer     :: xo(:) ! observables
    type (t_matrix_block)      :: Hnew  ! Jakobi matrix
    type (t_matrix_block)      :: Hi    ! inverse H

    select case (mode)

    case ('o')
      pso = obs% l% H * (psi - obs% l% x) + obs% l% y

    case ('l')
      pso = obs% l% H * psi

    case ('a')
      psi = pso * obs% l% H

    case ('i')
      do l = 1, size(obs% o)
        if (psi% s(l)% pe /= dace% pe) cycle
        do j = 1, obs% o(l)% n_spot
          spt => obs% o(l)% spot(j)
          xi  => psi %s(l)% x (spt% i% i + 1: spt% i% i + spt% i% n )
          xo  => pso %s(l)% x (spt% o% i + 1: spt% o% i + spt% o% n )
          call allocate_block (Hnew,FULL,spt% o% n,spt% i% n)
          Hnew% full = 0._mp
          do k=1, spt% i% n
            do m=obs% l% H% b(l,l)% ia (spt% i% i + k),      &
                 obs% l% H% b(l,l)% ia (spt% i% i + k + 1) - 1
              Hnew% full(obs% l% H% b(l,l)% ja (m) - spt% o% i,k) &
                       = obs% l% H% b(l,l)% packed(m)
            end do
          end do
          Hi = inverse (Hnew)
          call deallocate (Hnew)
          xi = Hi * xo (:)
          call deallocate (Hi)
        end do
      end do

    case default
      call finish ('apply_H','invalid value for mode')
    end select

    do l = 1, size(obs% o)
      select case (mode)
      case ('a')
        if (psi% global) call p_bcast (psi% s(l)% x, psi% s(l)% pe)
      case default
        if (pso% global) call p_bcast (pso% s(l)% x, pso% s(l)% pe)
      end select
    end do

  end subroutine apply_H
!------------------------------------------------------------------------------

  subroutine debug_spot_obs_set(obs, flags, hint, msg)
    type(t_obs_set),  intent(in)           :: obs
    integer,          intent(in), optional :: flags(:)
    character(len=*), intent(in), optional :: hint
    character(len=*), intent(in), optional :: msg
    integer :: ib
    if (.not.ldeb_spot) return
    do ib = 1, size(obs%o)
      if (obs%o(ib)%pe == dace%pe) &
           call debug_spot(obs%o(ib), flags=flags, hint=hint, msg=msg)
    end do
  end subroutine debug_spot_obs_set

end module mo_obs_set
