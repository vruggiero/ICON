!
!+ Donald E. Knuth's portable pseudo-random number generator
!
MODULE mo_random
!
! Description:
!   Donald E. Knuth's portable pseudo-random number generator,
!   from Seminumerical Algorithms by D E Knuth, 3rd edition (1997)
!   including the modifications made in the 9th printing (2002)
!   (Double precision version only).
!
!   The present implementation allows the creation of multiple
!   independent instances of the generator with differing seeds
!   for parallel harvesting and identified by a state variable.
!   The range of the (integer) seeds is: 0 <= seed <= (2**30-3).
!
!   For testing purposes there is also a static version of the
!   generator that generates the default sequence. (No seeding).
!
!   No facilities for saving states and restarting implemented yet.
!
!   This module is used commonly by the COSMO model and 3DVAR program packages !
!
! Current Code Owners:
!    For DWD 3DVAR:                        For COSMO:
!    DWD, Harald Anlauf                    DWD, Christoph Schraff
!    phone: +49 69 8062 4941               phone: +49 69 8062 2725
!    fax:   +49 69 8062 3721               fax:   +49 69 8062 3721
!    email: harald.anlauf@dwd.de           email: christoph.schraff@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_1         2008/11/05 Andreas Rhodin
!  First operational 3D-Var release
! V1_5         2009/05/25 Harald Anlauf
!  knuth_generate: minor optimizations, comment on non-vectorization
! V1_8         2009/12/09 Harald Anlauf
!  optimizations for NEC SX-9
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Harald Anlauf
!  Make FTRACE_REGIONs depending on macro DISABLE_FTRACE_REGION
! V1_26        2013/06/27 Andreas Rhodin
!  change "finish" to "model_abort" for compatibility with COSMO
! V1_45        2015-12-15 Harald Anlauf
!  fix for TR15581 support
! V1_51        2017-02-24 Andreas Rhodin
!  update shared files with COSMO 5.04d
!
! CAUTION: This module is used commonly by the 3DVAR and COSMO main programs.!!!
!!!        Therefore, anybody wanting to introduce a modification to this    !!!
!!!        module in the context of either of these programs must consult    !!!
!!!        the 'current code owner' of this module for the other program,    !!!
!!!        in order to allow for checking that the modification will comply  !!!
!!!        with both program packages. This must be done before the          !!!
!!!        modification is put into the Version Control System (VCS).        !!!
!
! Code Description:
! Language: Fortran 2003
! Software Standards:
!
! Author(s):
! Harald Anlauf (DWD) 2007-07
!
! Changes:
!==============================================================================
  !----------------------------------------------------------------------------
  ! Enhancements for vectorization
  !----------------------------------------------------------------------------
#if defined (__NEC__)
  ! Enable explicit strip-mining in subroutine knuth_generate
#define STRIP_MINING
#else
#undef  STRIP_MINING
#endif

#if defined (_FTRACE) && !defined (DISABLE_FTRACE_REGION)
#define FTRACE_BEGIN(text) CALL FTRACE_REGION_BEGIN (text)
#define FTRACE_END(text)   CALL FTRACE_REGION_END   (text)
  ! FTRACEing a user-defined region renders the function/subroutine impure
#define PURE
#else
#define FTRACE_BEGIN(text)
#define FTRACE_END(text)
#define PURE pure
#endif

#ifndef TR15581 /* should be available on all relevant platforms */
#define TR15581
#endif
  !============================================================================
  !-------------
  ! Modules used
  !-------------
  use kind_parameters,only: dp, i8         ! Double precision only (really!)
  use environment,    only: model_abort    ! abort in case of error
  !============================================================================
  implicit none
  private
  !============================================================================
  public :: random_number       ! Random number generator with state variable
  public :: random_state_t      ! Derived type: random generator state
  public :: construct           ! New random generator (state) from seed
  public :: destruct
  public :: assignment (=)      ! Clone a random generator state
  !----------------------------------------------------------------------------
  public :: test_random_knuth   ! Self-test of random number generator
  public :: test_sequences      ! Test independence of sequences
  public :: test_gaussian       ! Test implementations of normal distribution
  !----------------------------------------------------------------------------
  public :: random_static       ! Static random number generator
  !----------------------------------------------------------------------------
  public :: random_gauss        ! Generator for normal deviates
  public :: random_gauss_static ! Static version of normal deviate generator
  !----------------------------------------------------------------------------
  public :: DEFAULT_BUFFER_SIZE ! Default buffer size of generator
  public :: DEFAULT_SEED        ! Default seed (0)
  public :: MAX_SEED            ! Maximum seed (2^30-3)
  !----------------------------------------------------------------------------
  interface random_number
     module procedure random_knuth_scalar
     module procedure random_knuth_vector
     module procedure random_knuth_array2
     module procedure random_knuth_array3
  end interface

  interface random_static
     module procedure random_static_0d
     module procedure random_static_1d
     module procedure random_static_2d
     module procedure random_static_3d
  end interface

  interface random_gauss
     module procedure knuth_gauss_scalar
     module procedure knuth_gauss_array1
     module procedure knuth_gauss_array2
     module procedure knuth_gauss_array3
     module procedure random_gauss_static_0d
     module procedure random_gauss_static_1d
     module procedure random_gauss_static_2d
     module procedure random_gauss_static_3d
  end interface

  interface random_gauss_static
     module procedure random_gauss_static_0d
     module procedure random_gauss_static_1d
     module procedure random_gauss_static_2d
     module procedure random_gauss_static_3d
  end interface

  interface construct
     module procedure construct_knuth
  end interface

  interface destruct
     module procedure destruct_knuth
  end interface

  interface assignment (=)
     module procedure assign_random_state
  end interface
  !-------------------------------------
  ! General parameters of the generator:
  !-------------------------------------
  integer,  parameter :: DEFAULT_BUFFER_SIZE = 1009
  integer,  parameter :: DEFAULT_SEED = 0
  integer,  parameter :: MM           = 2**30
  integer,  parameter :: MAX_SEED     = MM - 3
  !--------------------------------------
  ! Internal parameters of the generator:
  !--------------------------------------
  integer,  parameter :: KK = 100
  integer,  parameter :: LL = 37
  integer,  parameter :: KKK = 2*KK-1
  real(dp), parameter :: ULP = 1/(2.0_dp**52)
  !------------------------------------
  ! The "raw" internal generator state:
  !------------------------------------
  type random_state_raw_t
     private
     real(dp)                 :: x(KK)
  end type random_state_raw_t
  !---------------------------------------
  ! The complete internal generator state:
  !---------------------------------------
  type random_state_t
     private
     integer                  :: last           ! Index of last random number
     integer                  :: buffer_end     ! Size of buffer
#ifdef TR15581
     real(dp), allocatable    :: buffer(:)              ! tr15581 & F2003
#else
     real(dp), pointer        :: buffer(:) => NULL ()   ! Fortran95
#endif
     type(random_state_raw_t) :: s
  end type random_state_t
  !
  ! (Note: using an allocatable buffer might result in faster code.)
  !
  !-----------------------------------------------------------------
  ! Static implementation for an almost replacement of RANDOM_NUMBER
  !-----------------------------------------------------------------
  logical                    :: s_init   = .false.
  type(random_state_t), save :: s_static
  !----------------------------------------------------------------------------
#if defined (__NEC__)
  integer :: gaussian_version = 2       ! 1=serial, 2=vectorized implementation
#else
  integer :: gaussian_version = 1       ! 1=serial, 2=vectorized implementation
#endif
  !----------------------------------------------------------------------------
contains
  !----------------------------------------------------------------------------
  PURE subroutine construct_knuth (state, seed, buffer_size)
    type(random_state_t), intent(inout) :: state
    integer,    optional, intent(in)    :: seed
    integer,    optional, intent(in)    :: buffer_size

    if (present (buffer_size)) then
       state% buffer_end = max (buffer_size, KKK)
    else
       state% buffer_end = DEFAULT_BUFFER_SIZE
    end if
#ifdef TR15581
    if (allocated  (state% buffer)) deallocate (state% buffer)
#else
    if (associated (state% buffer)) deallocate (state% buffer)
#endif
    allocate (state% buffer(state% buffer_end))
    state% last = state% buffer_end
    call construct_raw_state (state, seed)
  end subroutine construct_knuth
  !----------------------------------------------------------------------------
  PURE subroutine construct_raw_state (state, seed)
    type(random_state_t), intent(inout) :: state
    integer,    optional, intent(in)    :: seed

    call knuth_seed (state, seed)
  end subroutine construct_raw_state
  !----------------------------------------------------------------------------
  elemental subroutine destruct_knuth (state)
#ifdef TR15581
    type(random_state_t), intent(  out) :: state
#else
    type(random_state_t), intent(inout) :: state
    if(associated(state% buffer)) deallocate (state% buffer)
#endif
  end subroutine destruct_knuth
  !----------------------------------------------------------------------------
  PURE subroutine random_knuth_vector (harvest, state)
    real(dp),             intent(out)   :: harvest(:)
    type(random_state_t), intent(inout) :: state
    !----------------
    ! Local variables
    !----------------
    integer :: todo, done, chunk

    todo = size (harvest)
    done = 0
    do
       if (todo == 0) exit
       if (state% last >= state% buffer_end) then
          !---------------------
          ! Fill internal buffer
          !---------------------
          call knuth_generate (state)
          state% last = 0
       end if
       !---------------------------------------
       ! Fetch (next) chunk from current buffer
       !---------------------------------------
       chunk = min (todo, state% buffer_end - state% last)
       harvest(done+1:done+chunk) = &
            state% buffer(state% last+1:state% last+chunk)
       state% last = state% last + chunk
       done        = done + chunk
       todo        = todo - chunk
    end do
  end subroutine random_knuth_vector
  !----------------------------------------------------------------------------
  PURE subroutine knuth_generate (state, m)
    type(random_state_t), intent(inout) :: state
    integer, optional,    intent(in)    :: m
    !---------------------------------------------------------------
    ! Fortran version of "ranf_array"
    ! from Seminumerical Algorithms by D E Knuth, 3rd edition (1997)
    ! including the MODIFICATIONS made in the 9th printing (2002)
    ! ********* see the book for explanations and caveats! *********
    !---------------------------------------------------------------
    integer  :: j, n
    real(dp) :: y
#if defined (STRIP_MINING)
    integer  :: jj1, jj, mm
    integer, parameter :: LL3 = (KK-LL)/3
#endif

    n = state% buffer_end; if (present (m)) n = min (m, n)

    state% buffer(1:KK) = state% s% x(1:KK)

!!!FTRACE_BEGIN("knuth_generate:loop_1")
#if defined (STRIP_MINING)

#if __NEC_VERSION__ < 40000  /* Variant 1 */

    ! Explicit strip-mining for NEC SX to achieve partial vectorization
    ! (worked until nfort/3.5.1, started to fail with nfort/4.0.0)
    jj1 = max ((n-KK)/LL, 0)
    do jj = 1, jj1
       mm = KK+(jj-1)*LL
!NEC$ SHORTLOOP
!NEC$ ivdep
!NEC$ loop_count(LL)
       do j=mm+1, mm+LL
          y = state% buffer(j-KK) + state% buffer(j-LL)
          state% buffer(j) = y - int (y)
       end do
    end do
    mm = KK + jj1*LL
    ! The following loop is NOT vectorized with sxf90 revs.360,381,400
!NEC$ SHORTLOOP
!NEC$ IVDEP
!NEC$ loop_count(LL)
    do j=mm+1,n
       y = state% buffer(j-KK) + state% buffer(j-LL)
       state% buffer(j) = y - int (y)
    end do

#else  /* Variant 2 */

    ! Explicit strip-mining does NOT vectorize with sxf90 rev.400!?
    do mm = KK, n-1, LL
!NEC$ SHORTLOOP
!NEC$ IVDEP
       do j=mm+1, min (mm+LL, n)
          y = state% buffer(j-KK) + state% buffer(j-LL)
          state% buffer(j) = y - int (y)
       end do
    end do

#endif /* Variants */

#else
    ! Scalar machines: original code, no strip-mining
!DIR$ IVDEP SAFEVL=37
    do j=KK+1,n
       y = state% buffer(j-KK) + state% buffer(j-LL)
       state% buffer(j) = y - int (y)
    end do
#endif
!!!FTRACE_END  ("knuth_generate:loop_1")

!NEC$ SHORTLOOP
    do j=1,LL
       y = state% buffer(n+j-KK) + state% buffer(n+j-LL)
       state% s% x(j) = y - int (y)
    end do

!!!FTRACE_BEGIN("knuth_generate:loop_3")
#if defined (STRIP_MINING)
    ! Explicit strip-mining for NEC SX to achieve vectorization
    do mm = LL, KK-1, LL
!NEC$ SHORTLOOP
!NEC$ IVDEP
       do j=mm+1, min (mm+LL, KK)
          y = state% buffer(n+j-KK) + state% s% x(j-LL)
          state% s% x(j) = y - int (y)
       end do
    end do
#else
    ! Scalar machines: original code, no strip-mining
!DIR$ IVDEP SAFEVL=37
    do j=LL+1,KK
       y = state% buffer(n+j-KK) + state% s% x(j-LL)
       state% s% x(j) = y - int (y)
    end do
#endif
!!!FTRACE_END  ("knuth_generate:loop_3")
  end subroutine knuth_generate
  !----------------------------------------------------------------------------
  PURE subroutine random_knuth_scalar (harvest, state)
    real(dp),             intent(out)   :: harvest
    type(random_state_t), intent(inout) :: state

    real(dp) :: x(1)
    call random_knuth_vector (x, state)
    harvest = x(1)
  end subroutine random_knuth_scalar
  !----------------------------------------------------------------------------
  PURE subroutine random_knuth_array2 (harvest, state)
    real(dp),             intent(out)   :: harvest(:,:)
    type(random_state_t), intent(inout) :: state

    integer :: j
    do j = 1, size (harvest, dim=2)
       call random_knuth_vector (harvest(:,j), state)
    end do
  end subroutine random_knuth_array2
  !----------------------------------------------------------------------------
  PURE subroutine random_knuth_array3 (harvest, state)
    real(dp),             intent(out)   :: harvest(:,:,:)
    type(random_state_t), intent(inout) :: state

    integer :: j
    do j = 1, size (harvest, dim=3)
       call random_knuth_array2 (harvest(:,:,j), state)
    end do
  end subroutine random_knuth_array3
  !----------------------------------------------------------------------------
  PURE subroutine knuth_seed (state, seed)
    type(random_state_t), intent(inout) :: state
    integer, optional,    intent(in)    :: seed
    !----------------
    ! Local variables
    !----------------
    integer             :: j, s, t, seed_value
    integer, parameter  :: TT = 70
    real(dp)            :: ss, u(KKK), v

    if (present (seed)) then
       seed_value = seed
    else
       seed_value = DEFAULT_SEED
    end if
    if (seed_value < 0) then
       s = MM-1 - modulo (abs (seed_value)-1, MM)
    else
       s = modulo (seed_value, MM)
    end if
    ss = 2*ULP*(s+2)
    do j = 1, KK
       u(j) = ss
       ss = ss+ss
       if (ss >= 1d0) ss = ss-1+2*ULP
    end do
    u(2) = u(2) + ULP
    t = TT-1
    do
       do j = KK,2,-1
          u(2*j-1) = u(j)
          u(2*j-2) = 0
       end do
       do j = KKK,KK+1,-1
          v = u(j-(KK-LL)) + u(j)
          u(j-(KK-LL)) = v - int (v)
          v = u(j-KK) + u(j)
          u(j-KK) = v - int (v)
       end do
       if (modulo (s,2) == 1) then
          u(2:KK+1) = u(1:KK)
          u(1) = u(KK+1)
          v = u(LL+1) + u(KK+1)
          u(LL+1) = v - int (v)
       end if
       if (s /= 0) then
          s = s / 2
       else
          t = t - 1
       end if
       if (t <= 0) exit
    end do
    state% s% x(KK-LL+1:KK) = u(1:LL)
    state% s% x(1:KK-LL) = u(LL+1:KK)
    do j = 1,10
       call knuth_generate (state, KKK)
    end do
  end subroutine knuth_seed
  !============================================================================
  ! The code below contains static implementations of the generator that
  ! generate a single fixed sequence (no seeding implemented yet).
  !============================================================================
  subroutine random_static_0d (x)
    real(dp), intent(out) :: x
    if (.not. s_init) then
       call construct (s_static)
       s_init = .true.
    end if
    call random_number (x, s_static)
  end subroutine random_static_0d
  !----------------------------------------------------------------------------
  subroutine random_static_1d (x)
    real(dp), intent(out) :: x(:)
    if (.not. s_init) then
       call construct (s_static)
       s_init = .true.
    end if
    call random_number (x, s_static)
  end subroutine random_static_1d
  !----------------------------------------------------------------------------
  subroutine random_static_2d (x)
    real(dp), intent(out) :: x(:,:)
    if (.not. s_init) then
       call construct (s_static)
       s_init = .true.
    end if
    call random_number (x, s_static)
  end subroutine random_static_2d
  !----------------------------------------------------------------------------
  subroutine random_static_3d (x)
    real(dp), intent(out) :: x(:,:,:)
    if (.not. s_init) then
       call construct (s_static)
       s_init = .true.
    end if
    call random_number (x, s_static)
  end subroutine random_static_3d
  !============================================================================
  subroutine knuth_gauss_array1 (x, state)
    !-----------
    ! Parameters
    !-----------
    real(dp),             intent(out)   :: x(:)         ! Normal deviates
    type(random_state_t), intent(inout) :: state        ! Generator state

    select case (gaussian_version)
    case (1)
       call knuth_gauss_serial (x, state)
    case (2)
       call knuth_gauss_vector (x, state)
    case default
       write (0,*) "gaussian_version =", gaussian_version
       call model_abort (-1,-1,'invalid choice','knuth_gauss_array1')
       stop
    end select
  end subroutine knuth_gauss_array1
  !============================================================================
  PURE subroutine knuth_gauss_serial (x, state)
    !---------------------------------------------------------------------
    ! Gaussian random deviates using the Box-Muller method
    ! (See Knuth, Vol.2, 2nd ed., p.117, or Numerical Recipes, ch. 7.2)
    !
    ! Partially vectorized version by Harald Anlauf.
    !
    ! Beware: this implementation throws away superabundant random numbers,
    !         caching is not yet implemented.
    !-----------
    ! Parameters
    !-----------
    real(dp),             intent(out)   :: x(:)         ! Normal deviates
    type(random_state_t), intent(inout) :: state        ! Generator state
    !----------------
    ! Local variables
    !----------------
    integer,  parameter :: MAX_BUFFER_GAUSS = 512       ! Buffer size limit
    integer                             :: i, np, nbuf, done, todo
    real(dp)                            :: r, f
#ifndef __NEC__
    real(dp), allocatable, dimension(:) :: v1, v2
#else
    real(dp) :: v1(MAX_BUFFER_GAUSS*2), &! Fixed size buffers for
                v2(MAX_BUFFER_GAUSS)     ! efficient inlining on NEC
#endif

    todo = size (x)
    if (todo == 0) return
    !--------------------------------------------------------------
    ! Choice of size of the internal buffers: at each iteration
    ! we generate np pairs of uniform random numbers, with np
    ! determined by the algorithm's acceptance ratio (Pi/4 ~ 7/11).
    ! Tuning hint: adjust MAX_BUFFER_GAUSS to avoid cache thrashing.
    !--------------------------------------------------------------
    nbuf = min (((todo+1)/2)*14/11, MAX_BUFFER_GAUSS)
#ifndef __NEC__
    allocate (v1(nbuf), v2(nbuf))
#endif

    done = 0
!!!FTRACE_BEGIN("knuth_gauss_serial:outer")
!NEC$ novector
    outer: do
       np = min (((todo+1)/2)*14/11, nbuf)    ! Size of next chunk
#ifdef __NEC__
       call random_number (v1(1:2*np), state) ! Get np pairs of uniform deviates
       v1(1:np) = 2 * v1(   1:  np) - 1       ! Map to [-1:1,-1:1]
       v2(1:np) = 2 * v1(np+1:2*np) - 1
#else
       call random_number (v1(1:np), state)   ! Get np pairs of uniform deviates
       call random_number (v2(1:np), state)
       v1(1:np) = 2 * v1(1:np) - 1            ! Map to [-1:1,-1:1]
       v2(1:np) = 2 * v2(1:np) - 1
#endif
       do i = 1, np
          r = v1(i)**2 + v2(i)**2             ! Do they lie in the unit circle?
          if (r >= 1.0_dp) cycle              ! No, check next pair.
          f = sqrt (-2 * log(r) / r)
          done = done + 1
          x(done) = v1(i) * f
          todo = todo - 1
          if (todo == 0) exit outer
          done = done + 1
          x(done) = v2(i) * f
          todo = todo - 1
          if (todo == 0) exit outer
       end do
       if (todo <= 0) exit outer
    end do outer
!!!FTRACE_END  ("knuth_gauss_serial:outer")
  end subroutine knuth_gauss_serial
  !----------------------------------------------------------------------------
  PURE subroutine knuth_gauss_vector (x, state)
    !---------------------------------------------------------------------
    ! Gaussian random deviates using the Box-Muller polar method
    ! (See Knuth, Vol.2, 2nd ed., p.117, or Numerical Recipes, ch. 7.2)
    !
    ! Vectorized version borrowing from an F90 implementation by K.A.Hawick.
    ! Modified to produce the same sequence as the scalar implementation.
    ! Bug fixes and performance tweaks by Harald Anlauf.
    !
    ! Beware: this implementation throws away superabundant random numbers,
    !         caching is not yet implemented.
    !-----------
    ! Parameters
    !-----------
    real(dp),             intent(out)   :: x(:)         ! Normal deviates
    type(random_state_t), intent(inout) :: state        ! Generator state
    !----------------
    ! Local variables
    !----------------
    integer                             :: m, m1, m2, np, nbuf, done, todo
    integer,  parameter :: MAX_BUFFER_GAUSS = 512       ! Buffer size limit
#ifndef __NEC__
    real(dp), allocatable, dimension(:) :: v1, v2, r, s
    logical,  allocatable, dimension(:) :: mask
#else
    real(dp) :: v1(MAX_BUFFER_GAUSS*2), &! Fixed size buffers
                v2(MAX_BUFFER_GAUSS),   &! for efficient inlining
                r (MAX_BUFFER_GAUSS),   &! and reduced inlining overhead
                s (MAX_BUFFER_GAUSS)     ! on NEC SX Aurora
    integer  :: i
    real(dp) :: rtmp
#endif

    todo = size (x)
    if (todo == 0) return
    !--------------------------------------------------------------
    ! Choice of size of the internal buffers: at each iteration
    ! we generate np pairs of uniform random numbers, with np
    ! determined by the algorithm's acceptance ratio (Pi/4 ~ 7/11).
    ! Tuning hint: adjust MAX_BUFFER_GAUSS to avoid cache thrashing.
    !--------------------------------------------------------------
    nbuf = min (((todo+1)/2)*14/11, MAX_BUFFER_GAUSS)
#ifndef __NEC__
    allocate (v1(nbuf), v2(nbuf), r(nbuf), s(nbuf), mask(nbuf))
#endif

    done = 0
!!!FTRACE_BEGIN("knuth_gauss_vector:outer")
    outer: do while (todo > 0)
       np = min (((todo+1)/2)*14/11, nbuf)  ! Size of next chunk

#ifndef __NEC__
       call random_number (v1(1:np), state) ! Get np pairs of uniform deviates
       call random_number (v2(1:np), state)
       v1(1:np) = 2 * v1(1:np) - 1          ! Map to [-1:1,-1:1]
       v2(1:np) = 2 * v2(1:np) - 1
       r(1:np)  = min (v1(1:np)**2 + v2(1:np)**2, 1.0_dp)
       s(1:np)  = sqrt (-2 * log(r(1:np)) / r(1:np))
       v1(1:np) = v1(1:np) * s(1:np)
       v2(1:np) = v2(1:np) * s(1:np)
       mask(1:np) = r(1:np) < 1.0_dp        ! Do they lie in the unit circle?
       m = count (mask(1:np))               ! We'll obtain 2*m normal deviates
       if (m == 0) cycle                    ! "Bad luck", starting over...
       !-----------------------------------------------------------
       ! The buffers r,s receive the 2*m deviates.
       ! Interweave the buffers for the desired number of deviates.
       !-----------------------------------------------------------
       r(1:m)   = pack (v1(1:np), mask(1:np))
       s(1:m)   = pack (v2(1:np), mask(1:np))

#else  /* explicit code for what a good vectorizing compiler might achieve... */

       call random_number (v1(1:2*np), state) ! Get np pairs of uniform deviates
       m = 0
!NEC$ ivdep
!NEC$ nofma
       do i = 1, np
          v1(i) = 2 * v1(i   ) - 1            ! Map to [-1:1,-1:1]
          v2(i) = 2 * v1(i+np) - 1
          rtmp  = v1(i)**2 + v2(i)**2
          if (rtmp < 1.0_dp) then             ! Do they lie in the unit circle?
             s(i)  = sqrt (-2 * log(rtmp) / rtmp)
             v1(i) = v1(i) * s(i)
             v2(i) = v2(i) * s(i)
             m = m + 1                        ! We'll obtain 2*m normal deviates
             !-----------------------------------------------------------
             ! The buffers r,s receive the 2*m deviates.
             ! Interweave the buffers for the desired number of deviates.
             !-----------------------------------------------------------
             r(m) = v1(i)
             s(m) = v2(i)
          end if
       end do
       if (m == 0) cycle                    ! "Bad luck", starting over...
#endif
       m1 = min (m, (todo+1)/2)
       m2 = min (m, todo/2)
       x(done+1:done  +2*m1:2) = r(1:m1)
       x(done+2:done+1+2*m2:2) = s(1:m2)
       done = done + (m1+m2)
       todo = todo - (m1+m2)
    end do outer
!!!FTRACE_END  ("knuth_gauss_vector:outer")
  end subroutine knuth_gauss_vector
  !----------------------------------------------------------------------------
  subroutine knuth_gauss_scalar (x, state)
    real(dp),             intent(out)   :: x            ! Normal deviate
    type(random_state_t), intent(inout) :: state        ! Generator state

    real(dp) :: y(1)
    call random_gauss (y(:), state)
    x = y(1)
  end subroutine knuth_gauss_scalar
  !----------------------------------------------------------------------------
  subroutine knuth_gauss_array2 (x, state)
    real(dp),             intent(out)   :: x(:,:)       ! Normal deviates
    type(random_state_t), intent(inout) :: state        ! Generator state

    integer :: j
    do j = 1, size (x, dim=2)
       call random_gauss (x(:,j), state)
    end do
  end subroutine knuth_gauss_array2
  !----------------------------------------------------------------------------
  subroutine knuth_gauss_array3 (x, state)
    real(dp),             intent(out)   :: x(:,:,:)     ! Normal deviates
    type(random_state_t), intent(inout) :: state        ! Generator state

    integer :: j
    do j = 1, size (x, dim=3)
       call random_gauss (x(:,:,j), state)
    end do
  end subroutine knuth_gauss_array3
  !============================================================================
  subroutine random_gauss_static_0d (x)
    real(dp), intent(out) :: x
    if (.not. s_init) then
       call construct (s_static)
       s_init = .true.
    end if
    call random_gauss (x, s_static)
  end subroutine random_gauss_static_0d
  !----------------------------------------------------------------------------
  subroutine random_gauss_static_1d (x)
    real(dp), intent(out) :: x(:)
    if (.not. s_init) then
       call construct (s_static)
       s_init = .true.
    end if
    call random_gauss (x, s_static)
  end subroutine random_gauss_static_1d
  !----------------------------------------------------------------------------
  subroutine random_gauss_static_2d (x)
    real(dp), intent(out) :: x(:,:)
    if (.not. s_init) then
       call construct (s_static)
       s_init = .true.
    end if
    call random_gauss (x, s_static)
  end subroutine random_gauss_static_2d
  !----------------------------------------------------------------------------
  subroutine random_gauss_static_3d (x)
    real(dp), intent(out) :: x(:,:,:)
    if (.not. s_init) then
       call construct (s_static)
       s_init = .true.
    end if
    call random_gauss (x, s_static)
  end subroutine random_gauss_static_3d
  !============================================================================
  subroutine test_random_knuth ()
    !-------------------------------------------------
    ! Derived from official self-test of the generator
    !-------------------------------------------------
    !----------------
    ! Local variables
    !----------------
    integer             :: i, n0
    logical             :: success
    integer,  parameter :: TEST_SEED = 310952
    integer,  parameter :: DEFAULT_BUFFER_SIZE2 = DEFAULT_BUFFER_SIZE+1000
    integer,  parameter :: M = max (DEFAULT_BUFFER_SIZE,DEFAULT_BUFFER_SIZE2)
    real(dp)            :: a(M)
    real(dp), parameter :: R_2027082 = 1639783789831909_i8 * ULP
    real(dp), parameter :: A_2027082 = 0.36410514377569680455_dp
    real(dp), parameter :: EPS = epsilon (1.0_dp) / 2

    type(random_state_t) :: state, state2   ! Random generator state variables

    success = .true.
    n0 = 0

    if (abs (A_2027082 - R_2027082) > EPS) then
       write (0,'(2F23.20)') A_2027082, R_2027082
       call model_abort (-1,-1,'internal error','test_random_knuth')
    end if

    !------------------------------------------------------------
    ! Test 1: using default buffer size, output vector: same size
    !------------------------------------------------------------
    a = 0
    call construct (state, seed=TEST_SEED, buffer_size=DEFAULT_BUFFER_SIZE)
    do i = 1, (DEFAULT_BUFFER_SIZE2+1)
       call random_number (a(1:DEFAULT_BUFFER_SIZE),  state)
       n0 = n0 + count (a(1:DEFAULT_BUFFER_SIZE) == 0.0_dp)
    end do
    if (abs (a(1) - R_2027082) > EPS) then
       print *, "test_random_knuth: failure on test 1"
       print '(2F23.20)',a(1), R_2027082
       success = .false.
    end if
    call destruct  (state)
    !-------------------------------------------------------------
    ! Test 2: using enlarged buffer size, output vector: same size
    !-------------------------------------------------------------
    a = 0
    call construct (state, seed=TEST_SEED, buffer_size=DEFAULT_BUFFER_SIZE2)
    do i = 1, (DEFAULT_BUFFER_SIZE+1)
       call random_number (a(1:DEFAULT_BUFFER_SIZE2), state)
       n0 = n0 + count (a(1:DEFAULT_BUFFER_SIZE2) == 0.0_dp)
    end do
    if (abs (a(1) - R_2027082) > EPS) then
       print *, "test_random_knuth: failure on test 2"
       print '(2F23.20)',a(1), R_2027082
       success = .false.
    end if
    call destruct  (state)
    !----------------------------------------------------------------
    ! Test 3: using default buffer size, output vector: enlarged size
    !----------------------------------------------------------------
    a = 0
    call construct (state, seed=TEST_SEED, buffer_size=DEFAULT_BUFFER_SIZE)
    do i = 1, (DEFAULT_BUFFER_SIZE+1)
       call random_number (a(1:DEFAULT_BUFFER_SIZE2), state)
       n0 = n0 + count (a(1:DEFAULT_BUFFER_SIZE2) == 0.0_dp)
    end do
    if (abs (a(1) - R_2027082) > EPS) then
       print *, "test_random_knuth: failure on test 4"
       print '(2F23.20)',a(1), R_2027082
       success = .false.
    end if
    call destruct  (state)
    !-------------------------------------------------------
    ! Test 4: assignment (cloning) of random generator state
    !-------------------------------------------------------
    a = 0
    call construct (state, seed=TEST_SEED, buffer_size=DEFAULT_BUFFER_SIZE2)
    do i = 1, 100
       call random_number (a(1:DEFAULT_BUFFER_SIZE2), state)
       n0 = n0 + count (a(1:DEFAULT_BUFFER_SIZE2) == 0.0_dp)
    end do
    state2 = state
    call destruct  (state)
    do i = 101, (DEFAULT_BUFFER_SIZE+1)
       call random_number (a(1:DEFAULT_BUFFER_SIZE2), state2)
       n0 = n0 + count (a(1:DEFAULT_BUFFER_SIZE2) == 0.0_dp)
    end do
    if (abs (a(1) - R_2027082) > EPS) then
       print *, "test_random_knuth: failure on test 4"
       print '(2F23.20)',a(1), R_2027082
       success = .false.
    end if
    call destruct  (state2)

    if (n0 > 0) then
       print *, "test_random_knuth: number of unexpected zeros:", n0
       success = .false.
    end if
    if (success) then
       print *, "test_random_knuth: SUCCESS"
    else
       call model_abort (-1,-1,'FAIL','test_random_knuth')
    end if
  end subroutine test_random_knuth
  !============================================================================
  subroutine test_sequences ()
    !----------------------------------------------------
    ! Test independence of random sequences by setting up
    ! multiple independent generators...
    !----------------------------------------------------
    !----------------
    ! Local variables
    !----------------
    integer              :: i, j
    real(dp)             :: scal, sum1, sum2, sigma, tmp
    integer, parameter   :: N_gen   = 512       ! Number of generators
!   integer, parameter   :: vec_len = N_gen     ! Vector length for test
    integer, parameter   :: vec_len = N_gen*2   ! Vector length for test
    real(dp)             :: a(vec_len,N_gen)    ! Array of random sequences
    real(dp)             :: b(N_gen,N_gen)      ! "Covariance" matrix
    type(random_state_t) :: state(N_gen)        ! Array of generator states

    print *
    print '(1x,A,i0,A)', "test_sequences: creating ", N_gen, &
         " independent random generators"
    print *
    !-------------------------------------------------------------------
    ! Fill one row of array A with the first numbers from each generator
    !-------------------------------------------------------------------
    do i = 1, N_gen
!      call construct (state(i), seed=i)
       call construct (state(i), seed=i, buffer_size=max (N_gen,vec_len))
       call random_number (a(:,i), state(i))
    end do
    print *, "Overall check of first two moments (normalized):"
    sum1 = sum (a(:,:))    / (N_gen*vec_len) * 2
    sum2 = sum (a(:,:)**2) / (N_gen*vec_len) * 3
    print *, "sum1, sum2 =", sum1, sum2
    ! Expected error of the (scaled) mean:
    sigma = 1 / sqrt (12.0_dp*(N_gen*vec_len)) * 2
    print *, "Ratio of the error of the mean to the expected error:"
    print *, "|sum1-1|/sigma(sum1) =", abs (sum1-1) / sigma
    print *
    !-------------------
    ! Covariance matrix:
    !-------------------
    print *, "Calculating 'covariance matrix'"
    scal = sqrt (12.0_dp / vec_len)     ! Scaling factor for "unit variance"
    a = scal * (a - 0.5_dp)
    b = matmul (transpose (a), a)
    print *
    print *, "Upper left corner of covariance matrix:"
    print *
    j = 8
    do i = 1, j
       print '(9f8.4)', b(i,1:j)
    end do

    do i = 1, N_gen
       b(i,i) = b(i,i) - 1
    end do
    tmp = maxval (abs (b))
    print *, "maxval (abs (cov-1)):", tmp

    call destruct (state)
  end subroutine test_sequences
  !============================================================================
  subroutine test_gaussian ()
    !----------------
    ! Local variables
    !----------------
    integer              :: i, n0, k, saved_version
    logical              :: success
    integer,  parameter  :: TEST_SEED = 310952
    integer,  parameter  :: M = 512
    integer,  parameter  :: N = 512*16
    real(dp)             :: a(M)

    real(dp)             :: t0, t1, delta
    real(dp)             :: res(3,4,2)

    type(random_state_t) :: state ! Random generator state variable
    !------------------------------------------------------
    ! The reference output for the above set of parameters:
    !------------------------------------------------------
    real(dp), parameter  :: ref(3,4) = reshape ( (/ &
         -1.0329180943749412_dp,  0.4367814832920971_dp, &
         +1.4750718415862432_dp,  0.0385951392676669_dp, &
         -0.9493716438548104_dp,  0.4302813699551310_dp, &
         -1.0340399034907630_dp,  0.6976139403974778_dp, &
         +1.2362043267973737_dp,  1.8596221551872496_dp, &
         +0.0924647351019248_dp,  0.4898673204874738_dp  &
         /), (/3, 4/), order=(/1,2/) )
    real(dp), parameter  :: TOL = epsilon (1.0_dp)        ! Tolerance

    saved_version = gaussian_version
    success = .true.
    n0 = 0

    do k = 1, 2
       gaussian_version = k
       print *
       print *, "Testing gaussian_version =", gaussian_version
       a = 0
       call construct (state, seed=TEST_SEED, buffer_size=DEFAULT_BUFFER_SIZE)
       call cpu_time (t0)
       do i = 1, N
          call random_gauss (a(1:M), state)
          n0 = n0 + count (a(1:M) == 0.0_dp)
          if (i == 1) then
             res(1:3,1,k) = a(1:3)
             res(1:3,2,k) = a(M-2:M)
          else if (i == N) then
             res(1:3,3,k) = a(1:3)
             res(1:3,4,k) = a(M-2:M)
          end if
       end do
       call cpu_time (t1)
       print '("CPU time [s]:",f9.3)', t1-t0

       delta = maxval (abs (res(:,:,k)-ref))
       if (delta > TOL) then
          print *, "test_gaussian: test with reference run failed!"
          print *, "test_gaussian: max delta =", delta
          print *, "test_gaussian:        at :", &
               maxloc (transpose (abs (res(:,:,k)-ref)))
          success = .false.
          print *
          print *, "Result:"
          do i = 1, 4
             print '(3f20.16)', res(:,i,k)
          end do
          print *, "Reference:"
          do i = 1, 4
             print '(3f20.16)', ref(:,i)
          end do
          if (delta <= 2*EPSILON (1.0_dp)) then
             print *, "Note: very small differences can be due to FMA operations"
          end if
       end if
       call destruct  (state)

       if (n0 > 0) then
          print *, "test_gaussian: number of unexpected zeros:", n0
          success = .false.
       end if
    end do
    gaussian_version = saved_version

    if (success) then
       print *, "test_gaussian: SUCCESS"
    else
       call model_abort (-1,-1,'FAIL','test_gaussian')
    end if
  end subroutine test_gaussian
  !============================================================================
  subroutine assign_random_state (y, x)
    !-------------------------------
    ! Clone a random generator state
    !-------------------------------
    type(random_state_t), intent(inout) :: y
    type(random_state_t), intent(in)    :: x

#ifdef TR15581
    if (allocated  (y% buffer)) deallocate (y% buffer)
#else
    if (associated (y% buffer)) deallocate (y% buffer)
#endif
    allocate (y% buffer(x% buffer_end))

    y% last       = x% last
    y% buffer_end = x% buffer_end
    y% buffer(:)  = x% buffer(:)
    y% s          = x% s
  end subroutine assign_random_state
  !============================================================================
end module mo_random
