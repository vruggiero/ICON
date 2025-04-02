!
!+ Interface to Discrete Fourier Transform library routines
!
MODULE mo_dft
!
! Description:
!   Module 'mo_dft' holds interface routines to FFTW library for
!   DCT, DST, and FFT of real data in one or two dimensions.
!
!   DCT/DST is in-place.
!
!   Conventions:
!     DCT/DST, basic transformation: type-I (forward == inverse)
!     DCT/DST, "the" forward transform: type-II
!     DCT/DST, "the" inverse transform: type-III
!
!   Normalization may be chosen as follows:
!     0: strictly orthonormal transformation (default)
!     1: Custom normalization for FFTW (1 for forward transform)
!
!  Definitions (excluding normalization factors; Fortran index conventions):
!
!     DCT-I:   Y(k) = 2*sum(j=2,n-1)*X(j)*cos((j-1)*(k-1)*pi/(n-1))+X(1)+(-1)**(k-1)*X(n)
!     DCT-II:  Y(k) = 2*sum(j=1,n)  *X(j)*cos((j-1/2)*(k-1)*pi/n)
!     DCT-III: Y(k) = 2*sum(j=2,n)  *X(j)*cos((j-1)*(k-1/2)*pi/n) + X(1)
!
!     DST-I:   Y(k) = 2*sum(j=1,n)  *X(j)*sin(j*k*pi/(n+1))
!     DST-II:  Y(k) = 2*sum(j=1,n)  *X(j)*sin((j-1/2)*k*pi/n)
!     DST-III: Y(k) = 2*sum(j=1,n-1)*X(j)*sin(j*(k-1/2)*pi/n) + (-1)**(k-1)*X(n)
!
! Current Code Owner: DWD, Harald Anlauf
!    phone: +49 69 8062 4941
!    fax:   +49 69 8062 3721
!    email: harald.anlauf@dwd.de
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Author:
! Harald Anlauf  DWD  2018
!----------------------------------------------------------------
  !-------------
  ! modules used
  !-------------
  use               mo_kind,       only: wp     ! working precision
  use               mo_exception,  only: finish ! exit on error condition
  use               data_constants,only: pi     ! 3.14159...
#ifdef HAVE_FFTW3
  use, intrinsic :: iso_c_binding               ! Fortran/C binding constants
#endif
  implicit none
  !----------------
  ! public entities
  !----------------
  private
  public :: dct_1d                  ! 1-d DCT transform
  public :: dct_1dt                 ! 1-d DCT transform, transposed version
  public :: dct_2d                  ! 2-d DCT transform
  public :: poisson_solve_dct       ! Poisson solve using DCT, 2-d field
  public :: setup_dct               ! Choose set of transformation pair
  public :: test_dct                ! Basic test of transformations
!------------------------------------------------------------------------------
  public :: NORM_ORTHO              ! Flag: Strictly orthonormal transforms
  public :: NORM_CUSTOM             ! Flag: custom normalization for FFTW
  public :: DFT_FORWARD             ! Flag: forward transform
  public :: DFT_INVERSE             ! Flag: inverse transform
!==============================================================================
#ifdef HAVE_FFTW3
#ifdef __NEC__
  include 'aslfftw3.f03'
#else
  include 'fftw3.f03'
#endif
#endif
  !-----------
  ! interfaces
  !-----------
  interface dct_1d
     module procedure dct_1d_1      ! 1-d DCT transform, 1-d field
     module procedure dct_1d        ! 1-d DCT transform, 2-d field
  end interface dct_1d
!------------------------------------------------------------------------------
  interface dct_1dt
     module procedure dct_1dt       ! 1-d DCT transform, 2-d field, transposed
  end interface dct_1dt
!------------------------------------------------------------------------------
  interface dct_2d
     module procedure dct_2d        ! 2-d DCT transform, 2-d field
     module procedure dct_2d_3      ! 2-d DCT transform, 3-d field
  end interface dct_2d
!------------------------------------------------------------------------------
  interface poisson_solve_dct
     module procedure poisson_solve_dct_2  ! Poisson solve using DCT, 2-d field
     module procedure poisson_solve_dct_3  ! Poisson solve using DCT, 3-d field
  end interface poisson_solve_dct
!------------------------------------------------------------------------------
  !==========
  ! Constants
  !==========
  integer, parameter :: NORM_ORTHO  = 0    ! Normalization flags for
  integer, parameter :: NORM_CUSTOM = 1    ! transformation routines
!------------------------------------------------------------------------------
#ifdef HAVE_FFTW3
  integer, parameter :: DFT_FORWARD = FFTW_FORWARD
  integer, parameter :: DFT_INVERSE = FFTW_BACKWARD
#else
  integer, parameter :: DFT_FORWARD = -1   ! Fallback definitions
  integer, parameter :: DFT_INVERSE = +1
#endif
!------------------------------------------------------------------------------
  !=================
  ! Module Variables
  !=================
! type(C_PTR), save :: plan
  !-------------------------------------------
  ! Pairing of forward/inverse transformations
  !-------------------------------------------
  integer :: trans_set = 2
#ifdef HAVE_FFTW3
  integer :: trans_fwd = FFTW_REDFT10
  integer :: trans_inv = FFTW_REDFT01
#endif
!==============================================================================
contains
!==============================================================================
  subroutine setup_dct (set)
    integer, intent(in) :: set
    !----------------------------------
    ! Choose set of transformation pair
    !   1 : fwd/inv: REDFT00 (DCT-I)
    !   2 : forward: REDFT10 (DCT-II)
    !       inverse: REDFT01 (DCT-III)
    !----------------------------------
#ifndef HAVE_FFTW3
    call finish ("setup_dct","FFTW not linked")
#else
    select case (set)
    case (1)
       trans_set = 1
       trans_fwd = FFTW_REDFT00
       trans_inv = FFTW_REDFT00
    case (2)
       trans_set = 2
       trans_fwd = FFTW_REDFT10
       trans_inv = FFTW_REDFT01
    case default
       call finish ("setup_dct","invalid set")
    end select
#endif
  end subroutine setup_dct
!==============================================================================
  subroutine test_dct (nx, ny)
    integer, intent(in) :: nx, ny
    !-------------------------------------------------
    ! Basic test of DCT routines on given regular grid
    !-------------------------------------------------
    integer               :: set, nrm, k
    real(wp), allocatable :: a(:,:), b(:,:)
    real(wp)              :: tmp
    integer               :: norm(2) = [ NORM_ORTHO, NORM_CUSTOM ]
    real(wp), parameter   :: tol = 1.e-15_wp
    real(wp)              :: t1, t2

    print *, "test_dct: nx,ny =", nx, ny
    allocate (a(nx,ny), b(nx,ny))
    call random_number (a)

    call cpu_time (t1)
    do set = 1, 2
       print *
       print *, "Testing set", set
       call setup_dct (set=set)
       do k = 1, size (norm)
          nrm = norm(k)
          print *
          print *, "Testing normalization type", nrm
          b = a
          print *, "Forward first, then inverse"
          call dct_2d (b, isign=DFT_FORWARD, norm=nrm)
          call dct_2d (b, isign=DFT_INVERSE, norm=nrm)
          tmp = sum ((a-b)**2) / sum (a**2)
          print *, "Residual =", tmp
          if (abs (tmp) > tol) stop "Failed!"
          b = a
          print *, "Inverse first, then forward"
          call dct_2d (b, isign=DFT_INVERSE, norm=nrm)
          call dct_2d (b, isign=DFT_FORWARD, norm=nrm)
          tmp = sum ((a-b)**2) / sum (a**2)
          print *, "Residual =", tmp
          if (abs (tmp) > tol) stop "Failed!"
       end do
    end do
    call cpu_time (t2)
    print *
    print *, "CPU time:", real (t2-t1), "s"
  end subroutine test_dct
!==============================================================================
  subroutine dct_1d_1 (x, isign, norm)
    real(wp) ,intent(inout)           :: x(:)
    integer  ,intent(in)              :: isign
    integer  ,intent(in)    ,optional :: norm
    !-------------------------------------------
    ! One-dimensional discrete cosine transform:
    ! 1-d input array.
    !-------------------------------------------
#ifndef HAVE_FFTW3
    call finish ("dct_1d_1","FFTW not linked")
#else
    integer        :: m, lno
    real(wp)       :: nrm, norm0, norm1
    real(C_DOUBLE) :: t1(size(x,1)), t2(size(x,1))
    type(C_PTR)    :: plan

    lno = 0; if (present (norm)) lno = norm
    if (size (x) == 0) return
    m = size (x,1)
    select case (trans_set)
    case (1)
       norm1 = 1._wp / (2*(m-1))
    case default
       norm1 = 1._wp / (2*m)
    end select
    norm0 = sqrt (norm1)
    select case (isign)
    case (DFT_FORWARD)  ! "The" DCT
       plan = fftw_plan_r2r_1d (m, t1, t2, trans_fwd, FFTW_ESTIMATE) ! DCT-II
       if (lno == 0) then; nrm = norm0; else; nrm = 1._wp; end if
    case (DFT_INVERSE)  ! "The" IDCT
       plan = fftw_plan_r2r_1d (m, t1, t2, trans_inv, FFTW_ESTIMATE) ! DCT-III
       if (lno == 0) then; nrm = norm0; else; nrm = norm1; end if
    case default
       return
    end select
    if (.not. c_associated (plan)) then
       write(0,*) "dct_1d_1: No plan for m, isign=", m,             &
            merge ("DFT_FORWARD","DFT_INVERSE", isign == DFT_FORWARD)
       call finish ("dct_1d_1","FFTW: no plan")
    end if

    t1   = x(:)
    call fftw_execute_r2r  (plan, t1, t2)
    x(:) = t2 * nrm
    call fftw_destroy_plan (plan)
#endif

  end subroutine dct_1d_1
!==============================================================================
  subroutine dct_1d (x, isign, norm)
    real(wp) ,intent(inout)           :: x(:,:)
    integer  ,intent(in)              :: isign
    integer  ,intent(in)    ,optional :: norm
    !-------------------------------------------
    ! One-dimensional discrete cosine transform:
    ! 2-d input array.
    !-------------------------------------------
#ifndef HAVE_FFTW3
    call finish ("dct_1d","FFTW not linked")
#else
    integer        :: i, m, lno
    real(wp)       :: nrm, norm0, norm1
    real(C_DOUBLE) :: t1(size(x,1)), t2(size(x,1))
    type(C_PTR)    :: plan

    lno = 0; if (present (norm)) lno = norm
    if (size (x) == 0) return
    m = size (x,1)
    select case (trans_set)
    case (1)
       norm1 = 1._wp / (2*(m-1))
    case default
       norm1 = 1._wp / (2*m)
    end select
    norm0 = sqrt (norm1)
    select case (isign)
    case (DFT_FORWARD)  ! "The" DCT
       plan = fftw_plan_r2r_1d (m, t1, t2, trans_fwd, FFTW_ESTIMATE) ! DCT-II
       if (lno == 0) then; nrm = norm0; else; nrm = 1._wp; end if
    case (DFT_INVERSE)  ! "The" IDCT
       plan = fftw_plan_r2r_1d (m, t1, t2, trans_inv, FFTW_ESTIMATE) ! DCT-III
       if (lno == 0) then; nrm = norm0; else; nrm = norm1; end if
    case default
       return
    end select
    if (.not. c_associated (plan)) then
       write(0,*) "dct_1d: No plan for m, isign=", m,               &
            merge ("DFT_FORWARD","DFT_INVERSE", isign == DFT_FORWARD)
       call finish ("dct_1d","FFTW: no plan")
    end if

    do i = 1, size (x,2)
       t1     = x(:,i)
       call fftw_execute_r2r (plan, t1, t2)
       x(:,i) = t2 * nrm
    end do
    call fftw_destroy_plan (plan)
#endif

  end subroutine dct_1d
!==============================================================================
  subroutine dct_1dt (x, isign, norm)
    real(wp) ,intent(inout)           :: x(:,:)
    integer  ,intent(in)              :: isign
    integer  ,intent(in)    ,optional :: norm
    !-------------------------------------------
    ! One-dimensional discrete cosine transform:
    ! 2-d input array.  Transposed version.
    !-------------------------------------------
#ifndef HAVE_FFTW3
    call finish ("dct_1dt","FFTW not linked")
#else
    integer        :: i, m, lno
    real(wp)       :: nrm, norm0, norm1
    real(C_DOUBLE) :: t1(size(x,2)), t2(size(x,2))
    type(C_PTR)    :: plan

    lno = 0; if (present (norm)) lno = norm
    if (size (x) == 0) return
    m = size (x,2)
    select case (trans_set)
    case (1)
       norm1 = 1._wp / (2*(m-1))
    case default
       norm1 = 1._wp / (2*m)
    end select
    norm0 = sqrt (norm1)
    select case (isign)
    case (DFT_FORWARD)  ! "The" DCT
       plan = fftw_plan_r2r_1d (m, t1, t2, trans_fwd, FFTW_ESTIMATE) ! DCT-II
       if (lno == 0) then; nrm = norm0; else; nrm = 1._wp; end if
    case (DFT_INVERSE)  ! "The" IDCT
       plan = fftw_plan_r2r_1d (m, t1, t2, trans_inv, FFTW_ESTIMATE) ! DCT-III
       if (lno == 0) then; nrm = norm0; else; nrm = norm1; end if
    case default
       return
    end select
    if (.not. c_associated (plan)) then
       write(0,*) "dct_1dt: No plan for m, isign=", m,              &
            merge ("DFT_FORWARD","DFT_INVERSE", isign == DFT_FORWARD)
       call finish ("dct_1dt","FFTW: no plan")
    end if

    do i = 1, size (x,1)
       t1     = x(i,:)
       call fftw_execute_r2r (plan, t1, t2)
       x(i,:) = t2 * nrm
    end do
    call fftw_destroy_plan (plan)
#endif

  end subroutine dct_1dt
!==============================================================================
  subroutine dct_2d (x, isign, norm)
    real(wp) ,intent(inout)           :: x(:,:)
    integer  ,intent(in)              :: isign
    integer  ,intent(in)    ,optional :: norm
    !-------------------------------------------
    ! Two-dimensional discrete cosine transform:
    ! 2-d input array.
    !-------------------------------------------
#ifndef HAVE_FFTW3
    call finish ("dct_2d","FFTW not linked")
#else
    integer        :: n1, n2, lno
    real(wp)       :: nrm, norm0, norm1
    real(C_DOUBLE) :: tmp(size(x,1),size(x,2))
    type(C_PTR)    :: plan

    lno = 0; if (present (norm)) lno = norm
    if (size (x) == 0) return
    n1 = size (x,1)
    n2 = size (x,2)
    select case (trans_set)
    case (1)
       norm1 = 1._wp / (4*(n1-1)*(n2-1))
    case default
       norm1 = 1._wp / (4*n1*n2)
    end select
    norm0 = sqrt (norm1)
    select case (isign)
    case (DFT_FORWARD)  ! "The" DCT
       plan = fftw_plan_r2r_2d (n2, n1, tmp, tmp, trans_fwd, trans_fwd, &
                                FFTW_ESTIMATE)
       if (lno == 0) then; nrm = norm0; else; nrm = 1._wp; end if
    case (DFT_INVERSE)  ! "The" IDCT
       plan = fftw_plan_r2r_2d (n2, n1, tmp, tmp, trans_inv, trans_inv, &
                                FFTW_ESTIMATE)
       if (lno == 0) then; nrm = norm0; else; nrm = norm1; end if
    case default
       return
    end select
    if (.not. c_associated (plan)) then
       write(0,*) "dct_2d: No plan for n1, n2, isign=", n1, n2,     &
            merge ("DFT_FORWARD","DFT_INVERSE", isign == DFT_FORWARD)
       call finish ("dct_2d","FFTW: no plan")
    end if

    tmp = x
    call fftw_execute_r2r  (plan, tmp, tmp)
    x   = tmp * nrm

    call fftw_destroy_plan (plan)
#endif

  end subroutine dct_2d
!==============================================================================
  subroutine dct_2d_3 (x, isign, norm)
    real(wp) ,intent(inout)           :: x(:,:,:)
    integer  ,intent(in)              :: isign
    integer  ,intent(in)    ,optional :: norm
    !-------------------------------------------
    ! Two-dimensional discrete cosine transform:
    ! 3-d input array.
    !-------------------------------------------
#ifndef HAVE_FFTW3
    call finish ("dct_2d_3","FFTW not linked")
#else
    integer        :: n1, n2, lno, k
    real(wp)       :: nrm, norm0, norm1
    real(C_DOUBLE) :: tmp(size(x,1),size(x,2))
    type(C_PTR)    :: plan

    lno = 0; if (present (norm)) lno = norm
    if (size (x) == 0) return
    n1 = size (x,1)
    n2 = size (x,2)
    select case (trans_set)
    case (1)
       norm1 = 1._wp / (4*(n1-1)*(n2-1))
    case default
       norm1 = 1._wp / (4*n1*n2)
    end select
    norm0 = sqrt (norm1)
    select case (isign)
    case (DFT_FORWARD)  ! "The" DCT
       plan = fftw_plan_r2r_2d (n2, n1, tmp, tmp, trans_fwd, trans_fwd, &
                                FFTW_ESTIMATE)
       if (lno == 0) then; nrm = norm0; else; nrm = 1._wp; end if
    case (DFT_INVERSE)  ! "The" IDCT
       plan = fftw_plan_r2r_2d (n2, n1, tmp, tmp, trans_inv, trans_inv, &
                                FFTW_ESTIMATE)
       if (lno == 0) then; nrm = norm0; else; nrm = norm1; end if
    case default
       return
    end select
    if (.not. c_associated (plan)) then
       write(0,*) "dct_2d_3: No plan for n1, n2, isign=", n1, n2,   &
            merge ("DFT_FORWARD","DFT_INVERSE", isign == DFT_FORWARD)
       call finish ("dct_2d_3","FFTW: no plan")
    end if

    do k = 1, size (x,3)
       tmp      = x(:,:,k)
       call fftw_execute_r2r  (plan, tmp, tmp)
       x(:,:,k) = tmp * nrm
    end do
    call fftw_destroy_plan (plan)
#endif

  end subroutine dct_2d_3
!==============================================================================
  subroutine poisson_solve_dct_2 (phi, rho, dx, dy)
    real(wp) ,intent(inout)           :: phi(:,:)   ! Function to solve
    real(wp) ,intent(in)              :: rho(:,:)   ! Right-hand side
    real(wp) ,intent(in)    ,optional :: dx         ! Grid spacing
    real(wp) ,intent(in)    ,optional :: dy         ! Grid spacing
    !------------------------------------------
    ! Solve two-dimensional Poisson's equation
    ! (\nabla^2 phi = rho) on rectangular grid.
    ! 2-d input array.
    !------------------------------------------
#ifndef HAVE_FFTW3
    call finish ("poisson_solve_dct_2","FFTW not linked")
#else
    integer  :: n1, n2, i, j, k
    real(wp) :: t1, t2, d1, d2
    real(wp) :: tmp(size(rho,1),size(rho,2))
    real(wp) :: ilp(size(rho,1),size(rho,2))  ! Inverse Laplacian

    if (size (phi) == 0) return
    if (any (shape (phi) /= shape (rho))) then
       print *, "shape (phi) =", shape (phi)
       print *, "shape (rho) =", shape (rho)
       call finish ("poisson_solve_dct_2","shape (psi) /= shape (rho)")
    end if

    d1 = 1._wp; if (present (dx)) d1 = dx
    d2 = d1   ; if (present (dy)) d2 = dy
    n1 = size (rho,1)
    n2 = size (rho,2)
    t1 = (pi/(n1*d1)) ** 2
    t2 = (pi/(n2*d2)) ** 2
    ilp(1,1) = 0._wp
    do j = 1, n2
       k = 1; if (j == 1) k = 2
       do i = k, n1
          ilp(i,j) = - 1._wp / (t1*(i-1)**2 + t2*(j-1)**2)
       end do
    end do

    tmp      = rho(:,:)
    call dct_2d (tmp, isign=DFT_FORWARD)
    tmp(1,1) = 0._wp
    tmp      = tmp * ilp
    call dct_2d (tmp, isign=DFT_INVERSE)
    phi(:,:) = tmp
#endif

  end subroutine poisson_solve_dct_2
!==============================================================================
  subroutine poisson_solve_dct_3 (phi, rho, dx, dy)
    real(wp) ,intent(inout)           :: phi(:,:,:) ! Function to solve
    real(wp) ,intent(in)              :: rho(:,:,:) ! Right-hand side
    real(wp) ,intent(in)    ,optional :: dx         ! Grid spacing
    real(wp) ,intent(in)    ,optional :: dy         ! Grid spacing
    !------------------------------------------
    ! Solve two-dimensional Poisson's equation
    ! (\nabla^2 phi = rho) on rectangular grid.
    ! 3-d input array.
    !------------------------------------------
#ifndef HAVE_FFTW3
    call finish ("poisson_solve_dct_3","FFTW not linked")
#else
    integer  :: n1, n2, i, j, k
    real(wp) :: t1, t2, d1, d2
    real(wp) :: tmp(size(rho,1),size(rho,2))
    real(wp) :: ilp(size(rho,1),size(rho,2))  ! Inverse Laplacian

    if (size (phi) == 0) return
    if (any (shape (phi) /= shape (rho))) then
       print *, "shape (phi) =", shape (phi)
       print *, "shape (rho) =", shape (rho)
       call finish ("poisson_solve_dct_3","shape (psi) /= shape (rho)")
    end if

    d1 = 1._wp; if (present (dx)) d1 = dx
    d2 = d1   ; if (present (dy)) d2 = dy
    n1 = size (rho,1)
    n2 = size (rho,2)
    t1 = (pi/(n1*d1)) ** 2
    t2 = (pi/(n2*d2)) ** 2
    ilp(1,1) = 0._wp
    do j = 1, n2
       k = 1; if (j == 1) k = 2
       do i = k, n1
          ilp(i,j) = - 1._wp / (t1*(i-1)**2 + t2*(j-1)**2)
       end do
    end do

    do k = 1, size (rho,3)
       tmp        = rho(:,:,k)
       call dct_2d (tmp, isign=DFT_FORWARD)
       tmp(1,1)   = 0._wp
       tmp        = tmp * ilp
       call dct_2d (tmp, isign=DFT_INVERSE)
       phi(:,:,k) = tmp
    end do
#endif

  end subroutine poisson_solve_dct_3
!==============================================================================
end MODULE mo_dft
