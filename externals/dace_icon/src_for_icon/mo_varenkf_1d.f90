!
!+ localisation operator in 1 dimension
!
module mo_varenkf_1d
!
! Description:
!
!  Localisation operator in 1 dimension. To be used for vertical localisation
!  by the VarEnKf scheme.
!
! Current Maintainer: DWD, Harald Anlauf
!    phone: +49 69 8062 4941
!    fax:   +49 69 8062 3721
!    email: harald.anlauf@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_27        2013-11-08 Andreas Rhodin
!  implementation of vertical Localisation operator for VarEnKF
! V1_28        2014/02/26 Andreas Rhodin
!  further implementation of VarEnKF
! V1_31        2014-08-21 Andreas Rhodin
!  generate 2d random fields with specified correlation length scale
! V1_37        2014-12-23 Andreas Rhodin
!  check for negative argument to gaspari_cohn
! V1_42        2015-06-08 Andreas Rhodin
!  optimisation for PSAS; vertical wavelet transform
! V1_45        2015-12-15 Harald Anlauf
!  mo_varenkf_1d: optimize replacing some pointer components by allocatables
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Andreas Rhodin  DWD  2013
!==============================================================================
  !=============
  ! modules used
  !=============
  use mo_kind,        only: wp           ! working precision kind parameter
  use mo_mpi_dace,    only: dace         ! MPI group info
  use mo_exception,   only: finish       ! abort in case of error
  use mo_dace_string, only: char1        ! convert integer to character
  use mo_atm_state,   only: t_atm        ! atmospheric state derived type
  use mo_letkf_util,  only: set_lv,     &! set levels for vertical localisation
                            gaspari_cohn ! Gaspary & Cohn function
  implicit none

!------------------------------------------------------------------------------
  !================
  ! public entities
  !================
  private
  public :: t_C1           ! 1D localisation operator parameters
  public :: setup_c1       ! set meta data for 1D localisation operator
  public :: construct      ! set meta data for 1D localisation operator
  public :: destruct       ! deallocate 1D localisation operator meta data
  public :: apply_c1t      ! apply adjoint sqrt of 1D localisation operator
  public :: apply_c1       ! apply sqrt of 1D localisation operator
  public :: apply_c1_c1t   ! apply 1D localisation operator
  public :: apply_c11      ! apply sqrt of 1D operator in 2 dimensions
  public :: apply_vloc     ! apply   vertical localisation operator (mode 3)
  public :: apply_vloc_t   ! adjoint vertical localisation operator (mode 3)
  public :: apply_vint     ! vertical interpolation
  public :: apply_vint_t   ! adjoint
  public :: apply_vrint    ! vertical (reverse) interpolation
  public :: test_c1_modes  ! test accuricy depending on # of modes
  public :: test_c1_bound  ! test accuricy depending on distance from bound

!------------------------------------------------------------------------------
  !===========
  ! Interfaces
  !===========
  interface construct
    module procedure construct_c1
  end interface construct

  interface destruct
    module procedure destruct_c1
  end interface destruct

!------------------------------------------------------------------------------
  !=========================
  ! Derived type definitions
  !=========================

  !---------------------------------
  ! localisation operator parameters
  !---------------------------------
  type t_C1
    !----------------------
    ! primary specification
    !----------------------
    integer               :: n          ! number of gridpoints
    real(wp)              :: l_loc      ! localisation length scale
    integer               :: mode       ! 1=explicit, 2=hierarchical grid
    logical               :: cyclic     ! cyclic grid
    logical               :: lsqrt      ! sqrt factorisation
    logical               :: sym        ! use symmetric square root
    !---------------------------
    ! parameters for sym = false
    !---------------------------
    integer               :: nm         ! number of 'modes'
    integer               :: mo         ! offset
    integer               :: mi         ! increment
    !------------------------
    ! parameters for mode = 1
    !------------------------
    integer               :: nsup       ! support (no points) of C
    real(wp) ,allocatable :: c(:)       ! coefficients
    real(wp) ,allocatable :: s(:)       ! scaling factor for normalisation
    !------------------------
    ! parameters for mode = 2
    !------------------------
    integer               :: nl         ! # of levels
    integer               :: ndc        ! # of application of D on coarse grid
    integer               :: ndf        ! # of application of D on fine grids
    real(wp)              :: dc         ! diffusion coefficent on coarse grid
    real(wp)              :: df         ! diffusion coefficent on fine grid
    !------------------------
    ! parameters for mode = 3
    !------------------------
    integer               :: nzr        ! number of height levels
    real(wp)              :: lv_surf    ! localisation length at surface
    real(wp)              :: lv_top     ! localisation length at top
    real(wp)              :: r_0        ! log(p) reference
    real(wp)              :: flr_d      ! localisation length increment
    real(wp)              :: fp         ! transformation parameter
    real(wp)              :: f0         ! transformation parameter
    real(wp)              :: b          ! transformation parameter
    real(wp)              :: d          ! transformation parameter
    real(wp) ,pointer     :: pr (:) => NULL() !   p-levels on reduced grid (Pa)
    real(wp) ,pointer     :: flr(:) => NULL() ! loc.length on reduced grid
    real(wp) ,pointer     :: fr (:) => NULL() ! levels in transformed coord.
    real(wp) ,allocatable :: cg (:,:,:,:,:)   ! coefficients on grid
    integer  ,allocatable :: l0   (:,:,:,:)   ! index to coarse gridpoints
    real(wp) ,allocatable :: cr (:,:,:,:,:)   ! reverse interpolation coef.
    real(wp) ,allocatable :: a    (:,:,:,:)   ! normalisation factors
    integer               :: lb   (4)         ! lbound of l0, cg
  end type t_C1

!==============================================================================
contains
!==============================================================================

  subroutine setup_c1 (c1, mean, lv_surf, lv_top, n_vert)
  type (t_c1)  ,intent(inout) :: c1       ! localisation parameters
  type (t_atm) ,intent(in)    :: mean     ! for pressure levels
  real(wp)     ,intent(in)    :: lv_top   ! localisation length at model top
  real(wp)     ,intent(in)    :: lv_surf  ! localisation length at bottom
  integer      ,intent(in)    :: n_vert   ! derive coefficients cr, ar
  optional                    :: lv_surf, lv_top, n_vert
  !------------------------------------------------------------
  ! Set up vertical localisation coefficients for mode == 3
  ! If lv_surf, lv_top are not given, vertical levels are kept,
  !   merely horizontal columns are recalculated
  ! If n_vert is given, lv_surf, lv_top are adapted so that
  ! the actual number of gridpoints 'n_vert' is used.
  !------------------------------------------------------------

    integer    ,parameter :: mc = 5          ! max.number of coefs. / gridpoint
    integer               :: lb(4), ub(4)    ! bounds of model grid
    integer               :: i,j,k,d,l,l0,l1 ! indices
    real(wp) ,allocatable :: f  (:)          ! transformed coordinate
    real(wp) ,allocatable :: df (:)          !
    real(wp)              :: a               ! normalisation factor
    real(wp) ,allocatable :: ar (:)          ! normalisation factor
    integer               :: nv              ! derive coefficients cr, ar
    real(wp) ,parameter   :: sq103 = sqrt (10._wp/3._wp) ! scaling gaspari-cohn

    nv = 0; if (present (n_vert)) nv = n_vert
    if (.not. (present(lv_surf) .and. present(lv_top))) then
      !------------------------------------------
      ! re-initialisation with new set of columns
      !------------------------------------------
      if (c1% mode /= 3) call finish('setup_c1','reinitialisation, mode /= 3')
      if (allocated (c1% a)) deallocate (c1% a,  c1% cr)
                             deallocate (c1% cg, c1% l0)
    else
      !------------------------------
      ! derive levels, linear in logp
      !------------------------------
      c1% lv_surf = lv_surf / sqrt(2._wp)
      c1% lv_top  = lv_top  / sqrt(2._wp)
!NEC$ nomove
      do
        c1% mode    =  3
        c1% nzr     = -1   ! special value to derive suitable spacing
        call set_lv (c1% nzr, c1% lv_surf, c1% lv_top,            &
                     c1% r_0, c1% flr_d,   c1% pr, c1% flr, mean, &
                     c1% fr,  c1% fp,      c1% f0, c1% b,   c1% d )
        if (nv == 0)       exit
        if (nv == c1% nzr) exit
        call destruct (c1)
        c1% lv_surf = c1% lv_surf * (c1% nzr / real(nv,wp))
        c1% lv_top  = c1% lv_top  * (c1% nzr / real(nv,wp))
      end do
    endif
    !--------------------
    ! derive coefficients
    !--------------------
    lb = mean% lb
    ub = mean% ub
    c1% lb = lb
    allocate  (c1% cg (mc, lb(1):ub(1), lb(2):ub(2), lb(3):ub(3), lb(4):ub(4)))
    allocate  (c1% l0     (lb(1):ub(1), lb(2):ub(2), lb(3):ub(3), lb(4):ub(4)))
    allocate  (f          (                          lb(3):ub(3)             ))
    if (nv>0) then
     allocate (df         (                          lb(3):ub(3)             ))
     allocate (c1% a      (lb(1):ub(1), lb(2):ub(2), lb(3):ub(3), lb(4):ub(4)))
     allocate (c1% cr (mc, lb(1):ub(1), lb(2):ub(2), lb(3):ub(3), lb(4):ub(4)))
     allocate (ar     (c1% nzr))
    endif
    do d = lb(4), ub(4)
    do j = lb(2), ub(2)
    do i = lb(1), ub(1)
    if (c1% lv_surf >= c1% lv_top) then
      do k = lb(3), ub(3)
        f(k) =  sq103 * sqrt(2._wp) * log(mean% pf (i,j,k,d))
      end do
    else
      do k = lb(3), ub(3)
        f(k) = c1% f0 + c1% fp * log (c1% r_0 - log(mean% pf (i,j,k,d)))
      end do
    end if
    if (nv>0) then
      do k = lb(3)+1, ub(3)-1
        df(k) = (f(k+1) - f(k-1)) / 2._wp
      end do
      df (lb(3)) = f(lb(3)+1) - f(lb(3)  )
      df (ub(3)) = f(ub(3)  ) - f(ub(3)-1)
      ar         = 0._wp
    endif
    do k = lb(3), ub(3)
      l0   = nint ((f(k) - c1% b) / c1% d) - 2
      c1% l0 (     i,j,k,d) = l0
      c1% cg (1:mc,i,j,k,d) = 0._wp
!NEC$ unroll_completely
      do l = 1, mc
        l1 = l + l0
        if (l1 < 1 .or. l1 > c1% nzr) cycle
        c1% cg (l,i,j,k,d) = gaspari_cohn (abs(f(k) - c1% fr(l1)), 1._wp)
        if (nv>0) then
          c1% cr (l,i,j,k,d) = c1% cg (l,i,j,k,d) * df (k)
          ar     (l1)        = ar     (l1)        + c1% cr (l,i,j,k,d)
        endif
      end do
      !----------
      ! Normalise
      !----------
      a = 1._wp / sqrt (sum (c1% cg (:,i,j,k,d)))
      c1% cg (:,i,j,k,d) = c1% cg (:,i,j,k,d) * a
      if (nv>0) c1% a  (i,j,k,d) = a
    end do
    if (nv>0) then
      where (ar /= 0._wp) ar = 1._wp / ar
      do k = lb(3), ub(3)
        l0 = c1% l0 (i,j,k,d)
        c1% cr (1:mc,i,j,k,d) = 0._wp
!NEC$ unroll_completely
        do l = 1, mc
          l1 = l + l0
          if (l1 < 1 .or. l1 > c1% nzr) cycle
          c1% cr (l,i,j,k,d) = c1% cr (l,i,j,k,d) * ar (l1)
        end do
      end do
    endif
    end do
    end do
    end do

  end subroutine setup_c1

!------------------------------------------------------------------------------

  subroutine construct_c1 (c1, n, l_loc, cyclic, lsqrt, stride, shift, modes)
  !--------------------------------------------------------
  ! Set up vertical localisation coefficients for mode == 1
  !--------------------------------------------------------
  type (t_c1) ,intent(out)          :: c1     ! localisation parameters
  integer     ,intent(in)           :: n      ! number of gridpoints
  real(wp)    ,intent(in)           :: l_loc  ! localisation length scale
  logical     ,intent(in) ,optional :: cyclic ! cyclic boundary conditions
  logical     ,intent(in) ,optional :: lsqrt  ! derive factorisation
  integer     ,intent(in) ,optional :: stride ! # modes / # gridpoints
  integer     ,intent(in) ,optional :: shift  ! shift modes vs. gridpoints
  integer     ,intent(in) ,optional :: modes  ! number of modes

    !----------------
    ! local variables
    !----------------
    real(wp) :: f_loc ! length scale parameter for Gaspari & Cohn function
    real(wp) :: f_sup ! maximum support        for Gaspari & Cohn function
    integer  :: i, j
    real(wp) ,allocatable :: x(:)

    !---------------------------
    ! set up mandatory arguments
    !---------------------------
    c1% n      = n
    c1% l_loc  = l_loc
    !---------------------------
    ! set up optional parameters
    !---------------------------
    c1% cyclic = .false. ;if (present (cyclic)) c1% cyclic = cyclic
    c1% lsqrt  = .true.  ;if (present (lsqrt )) c1% lsqrt  = lsqrt
    !----------------------------
    ! options not yet implemented
    !----------------------------
    c1% mode  = 1
    !-------------------------------
    ! set default derived parameters
    !-------------------------------
    c1% sym  = .true.
    c1% nm   = c1% n
    c1% mo   = 0
    c1% mi   = 1
    c1% nsup = 0
    !--------------------------
    ! set up derived parameters
    !--------------------------
    select case (c1% mode)
    case (1)
      !---------------------------------
      ! explicit Gaspary & Cohn function
      !---------------------------------
      f_loc    = sqrt (10._wp/3._wp) * c1% l_loc
      if (c1% lsqrt) f_loc = f_loc / sqrt(2._wp)
      f_sup    = 2._wp * f_loc
      c1% nsup = int (f_sup)
      allocate (c1% c (0 : c1% nsup))
      do i = 0, c1% nsup
        c1% c(i) = gaspari_cohn (real(i,wp), f_loc)
      end do
    case default
      call finish('construct_c1','mode = '//char1(c1%mode)//' not implemented')
    end select

    !------------------------------------------------
    ! options applicable for factorised matrices only
    !------------------------------------------------
    if (c1% lsqrt) then
      !--------------
      ! normalisation
      !--------------
      allocate (x (0:n-1))
      x = 0._wp; x(n/2) = 1._wp
      call apply_c1_c1t (c1, x)
      c1% c = c1% c / sqrt (x(n/2))
      !---------------
      ! default stride
      !---------------
      i = max (1, int(c1% l_loc))
      do
        c1% mi = i
        c1% nm = c1% n / i
        if (.not.c1% cyclic)          exit
        if (c1% mi * c1% nm == c1% n) exit
        i = i - 1
      end do
      !----------------
      ! explicit stride
      !----------------
      if (present (stride)) c1% mi = stride
      if (c1% cyclic) then
        c1% nm = c1% n / c1% mi
        if (c1% mi * c1% nm /= c1% n) then
          write(0,*)'invalid stride:', c1% mi,c1% nm,c1% mi * c1% nm,c1% n
          call finish ('construct_c1','invalid stride')
        endif
      endif
      !--------------------------------
      ! default offset, number of modes
      !--------------------------------
      if (c1% mi /= 1 .and. .not. c1% cyclic) then
        c1% nm = 1
        do
          i = c1% mi * c1% nm
          if (i>= c1% n) exit
          c1% nm = c1% nm + 1
        end do
        c1% nm = c1% nm + 2
        c1% mo = (c1% n - c1% mi * c1% nm) / 2
      endif
      !----------------
      ! explicit offset
      !----------------
      if (present (shift))  c1% mo = shift
      if (present (modes))  c1% nm = modes
      !----------------
      ! symetric sqrt ?
      !----------------
      if (c1% mi /= 1)      c1% sym  = .false.
      if (.not. c1% cyclic) c1% sym  = .false.
      !-------------------------
      ! individual normalisation
      !-------------------------
      if (.not. c1% sym) then
        allocate (c1% s (0:n-1)) ;c1% s = 1._wp
        if (c1% cyclic) then
          do i = 0, c1% mi -1
            x = 0._wp; x(i) = 1._wp
            call apply_c1_c1t (c1, x)
            c1% s(i) = 1._wp / sqrt (x(i))
            do j = i, n-1, c1% mi
              c1% s(j) = 1._wp / sqrt (x(i))
            end do
          end do
        else
          do i = 0, n-1
            x = 0._wp; x(i) = 1._wp
            call apply_c1_c1t (c1, x)
            if (x(i)>0._wp) c1% s(i) = 1._wp / sqrt (x(i))
          end do
        endif  ! not cyclic
      endif    ! sym
    endif      ! lsqrt
    !---------
    ! printout
    !---------
    if (dace% lpio) then
      write(6,'()')
      write(6,'(a)') '  1-D localisation operator parameters:'
      write(6,'()')
      write(6,'(a,i3,2x,a)') '    n      = ',c1% n     ,' ! number of gridpoints'
      write(6,'(a,f5.1 ,a)') '    l_loc  = ',c1% l_loc ,' ! localisation length scale'
      write(6,'(a,i3,2x,a)') '    mode   = ',c1% mode  ,' ! mode: 1=explicit, 2=hierarchical grid'
      write(6,'(a,l3,2x,a)') '    cyclic = ',c1% cyclic,' ! cyclic grid'
      write(6,'(a,l3,2x,a)') '    sqrt   = ',c1% lsqrt ,' ! use sqrt'
      write(6,'(a,l3,2x,a)') '    sym    = ',c1% sym   ,' ! use symmetric square root'
      write(6,'(a,i3,2x,a)') '    nm     = ',c1% nm    ," ! number of 'modes'"
      write(6,'(a,i3,2x,a)') '    mo     = ',c1% mo    ,' ! offset'
      write(6,'(a,i3,2x,a)') '    mi     = ',c1% mi    ,' ! increment'
      write(6,'(a,i3,2x,a)') '    nsup   = ',c1% nsup  ,' ! support (no points) of C'
      write(6,'()')
    endif
  end subroutine construct_c1

!------------------------------------------------------------------------------

  subroutine destruct_c1 (c1)
  !-----------------------------------------------
  ! deallocate 1D localisation operator parameters
  !-----------------------------------------------
  type (t_c1) ,intent(inout) :: c1
    if (allocated  (c1% c  )) deallocate (c1% c  )
    if (allocated  (c1% s  )) deallocate (c1% s  )
    if (associated (c1% pr )) deallocate (c1% pr )
    if (associated (c1% flr)) deallocate (c1% flr)
    if (associated (c1% fr )) deallocate (c1% fr )
    if (allocated  (c1% cg )) deallocate (c1% cg )
    if (allocated  (c1% l0 )) deallocate (c1% l0 )
    if (allocated  (c1% a  )) deallocate (c1% a  )
    if (allocated  (c1% cr )) deallocate (c1% cr )
    c1% mode = 0
  end subroutine destruct_c1

!------------------------------------------------------------------------------

  subroutine apply_c11 (cx, cy, lby, y, z)
  !------------------------------------------
  ! apply sqrt of 1D operator in 2 dimensions
  !------------------------------------------
  type (t_c1) ,intent(in)  :: cx         ! localisation parameters
  type (t_c1) ,intent(in)  :: cy         ! localisation parameters
  integer     ,intent(in)  :: lby(2)     ! y lower bound (offset)
  real(wp)    ,intent(out) :: y (:,:,:)  ! gridded field (MPI partitioned)
  real(wp)    ,intent(in)  :: z (:,:,:)  ! auxiliary coarse grid field

    real(wp) :: w (size(y,1),size(z,2)) !,size(z,3))
    integer  :: i,j,k

    do k = 1, size(z,3)
      do j = 1, size(z,2)
        call apply_c1 (cx, lby(1), w (:,j  ), z (:,j,k))
      end do
!   end do
!   do k = 1, size(z,3)
      do i = 1, size(y,1)
        call apply_c1 (cy, lby(2), y (i,:,k), w (i,:  ))
      end do
    end do

  end subroutine apply_c11

!------------------------------------------------------------------------------

  subroutine apply_c1_c1t (c1, x)
  type (t_c1) ,intent(in)    :: c1     ! localisation parameters
  real(wp)    ,intent(inout) :: x(0:)  ! argument / result
  !-----------------------------------------------
  ! apply the 1D localisation operator to a vector
  !-----------------------------------------------
    real(wp) ,allocatable :: y (:) ! temporary storage for result
    integer               :: i     ! result index
    integer               :: j     ! argument index
    integer               :: n     ! size  of argument and result

    n = size(x)
    !--------------------------------------------
    ! explicit formulation, no sqrt factorisation
    !--------------------------------------------
    if (.not. c1% lsqrt) then
      allocate (y(0:n-1))
      y = 0._wp
      if (c1% cyclic) then
        do i = 0, n-1
          y(i) = y(i) + c1% c(0) * x(i)
          do j = 1, c1% nsup
            y(i) = y(i) + c1% c(j) * (x(modulo(i+j,n)) + x(modulo(i-j,n)))
          end do
        end do
      else
        do i = 0, n-1
          y(i) = y(i) + c1% c(0) * x(i)
          do j = 1, min (c1% nsup, i)
            y(i) = y(i) + c1% c(j) *  x(i-j)
          end do
          do j = 1, min (c1% nsup, n-1-i)
            y(i) = y(i) + c1% c(j) *  x(i+j)
          end do
        end do
      end if
      x = y
    !------------------------
    ! factorisation C1 * C1^T
    !------------------------
    else
      allocate (y (c1% nm))
      call apply_c1t (c1,    y, x)
      call apply_c1  (c1, 0, x, y)
    end if

  end subroutine apply_c1_c1t

!------------------------------------------------------------------------------

  subroutine apply_c1t (c1, z, x)
  type (t_c1) ,intent(in)  :: c1      ! localisation parameters
  real(wp)    ,intent(out) :: z(0:)   ! result
  real(wp)    ,intent(in)  :: x(0:)   ! argument
  !------------------------------------------------
  ! apply the (transpose) of the sqrt factorisation
  !------------------------------------------------
    integer           :: n    ! size of argument x
    integer           :: m    ! size of result z
    integer           :: i    ! result index
    integer           :: k    ! corresponding argument index of the mode
    integer           :: j    ! offset to k
    real(wp) ,pointer :: y(:) ! temporary (for scaling)
    target            :: x

    if (.not. c1% lsqrt) call finish ('apply_c1t','lsqrt==F')

    n = size (x)
    m = size (z)

    !-------------------------------------------------
    ! individual scaling of gridpoints
    ! (required in # of modes /= number of gridpoints)
    !-------------------------------------------------
    y => x
    if (allocated (c1% s)) then
      allocate (y (0:n-1))
      y = x * c1% s
    endif

    !----------------------
    ! matrix multiplication
    !----------------------
    z = 0._wp
    if (c1% cyclic) then
      do i = 0, m-1
        k = c1% mo + c1% mi * i
        z(i) = z(i) + c1% c(0) * y(k)
        do j = 1, c1% nsup
          z(i) = z(i) + c1% c(j) * (y(modulo(k+j,n)) + y(modulo(k-j,n)))
        end do
      end do
    else
      do i = 0, m-1
        k = c1% mo + c1% mi * i
        do j = max (-c1% nsup,    -k) , &
               min ( c1% nsup, n-1-k)
          z(i) = z(i) + c1% c(abs(j)) *  y(k+j)
        end do
      end do
    end if

    if (allocated (c1% s)) deallocate (y)

  end subroutine apply_c1t

!------------------------------------------------------------------------------

  subroutine apply_c1 (c1, lby, y, z)

  type (t_c1) ,intent(in)  :: c1        ! localisation parameters
  integer                  :: lby       ! y lower bound
  real(wp)    ,intent(out) :: y(lby:)   ! result
  real(wp)    ,intent(in)  :: z(0:)     ! argument

    integer :: n   ! size of result   y
    integer :: m   ! size of argument z
    integer :: i   ! argument index
    integer :: k   ! corresponding result index of the mode
    integer :: j   ! offset to k
    integer :: l   ! result index
    integer :: uby ! y upper bound

    if (.not. c1% lsqrt) call finish ('apply_c1t','lsqrt==F')

    n   = size (y)
    m   = size (z)
    uby = lby + n - 1

    !----------------------
    ! matrix multiplication
    !----------------------
    y = 0._wp
    if (c1% cyclic) then
      do i = 0, m-1
        k = c1% mo + c1% mi * i
        y(k) = y(k) + c1% c(0) * z(i)
        do j = 1, c1% nsup
          l = modulo(k+j,n); y(l) = y(l) + c1% c(j) * z(i)
          l = modulo(k-j,n); y(l) = y(l) + c1% c(j) * z(i)
        end do
      end do
    else
      do i = 0, m-1
        k = c1% mo + c1% mi * i
        do j = max (-c1% nsup, lby-k) , &
               min ( c1% nsup, uby-k)
          y(k+j) = y(k+j) + c1% c(abs(j)) *  z(i)
        end do
      end do
    end if

    !-------------------------------------------------
    ! individual scaling of gridpoints
    ! (required in # of modes /= number of gridpoints)
    !-------------------------------------------------
    if (allocated (c1% s)) y = y * c1% s(lby:uby)

  end subroutine apply_c1

!------------------------------------------------------------------------------

  subroutine apply_vloc (c1, a, v)
  type (t_c1) ,intent(in)  :: c1
  real(wp)    ,intent(out) :: a (c1%lb(1):,c1%lb(2):,:,:)
  real(wp)    ,intent(in)  :: v (c1%lb(1):,c1%lb(2):,:,:)
#ifdef _CRAYFTN
  contiguous :: a, v
#endif

    integer  :: i,j,k,d,l,m,n,u
    real(wp) :: w

    n = size (c1% cg, 1)
    u = size (v, 3)
    a = 0._wp
    do d = lbound(a,4), ubound(a,4)
    do k = lbound(a,3), ubound(a,3)
    do j = lbound(a,2), ubound(a,2)
#ifdef __NEC__
      do l = 1, n
!NEC$ ivdep
    do i = lbound(a,1), ubound(a,1)
      m = c1% l0 (i,j,k,d)
        w = c1% cg (l,i,j,k,d)
        if ( (1 <= l+m) .and. (l+m <= u) ) &
                        a (i,j,k,d) = a (i,j,k,d) + v (i,j,l+m,d) * w
      end do
    end do
#else
    do i = lbound(a,1), ubound(a,1)
      m = c1% l0 (i,j,k,d)
      do l = 1, n
        w = c1% cg (l,i,j,k,d)
        if (w /= 0._wp) a (i,j,k,d) = a (i,j,k,d) + v (i,j,l+m,d) * w
      end do
    end do
#endif
    end do
    end do
    end do
  end subroutine apply_vloc

!------------------------------------------------------------------------------

  subroutine apply_vint (c1, a, v)
  type (t_c1) ,intent(in)  :: c1
  real(wp)    ,intent(out) :: a (c1%lb(1):,c1%lb(2):,:,:)
  real(wp)    ,intent(in)  :: v (c1%lb(1):,c1%lb(2):,:,:)

    integer  :: i,j,k,d,l,m,n
    real(wp) :: w

    n = size (c1% cg, 1)
    a = 0._wp
    do d = lbound(a,4), ubound(a,4)
    do k = lbound(a,3), ubound(a,3)
    do j = lbound(a,2), ubound(a,2)
    do i = lbound(a,1), ubound(a,1)
      m = c1% l0 (i,j,k,d)
      do l = 1, n
        w = c1% cg (l,i,j,k,d)
        if (w /= 0._wp) then
          w = w * c1% a (i,j,k,d)
          a (i,j,k,d) = a (i,j,k,d) + v (i,j,l+m,d) * w
        endif
      end do
    end do
    end do
    end do
    end do
  end subroutine apply_vint

!------------------------------------------------------------------------------

  subroutine apply_vint_t (c1, z, v)
  type (t_c1) ,intent(in)  :: c1
  real(wp)    ,intent(in)  :: z (c1%lb(1):,c1%lb(2):,:,:)
  real(wp)    ,intent(out) :: v (c1%lb(1):,c1%lb(2):,:,:)

    integer  :: i,j,k,d,l,m,n
    real(wp) :: w

    n = size (c1% cg, 1)
    v = 0._wp
    do d = lbound(z,4), ubound(z,4)
    do k = lbound(z,3), ubound(z,3)
    do j = lbound(z,2), ubound(z,2)
    do i = lbound(z,1), ubound(z,1)
      m = c1% l0 (i,j,k,d)
      do l = 1, n
        w = c1% cg (l,i,j,k,d)
        if (w /= 0._wp) then
          w = w * c1% a (i,j,k,d)
          v (i,j,l+m,d) = v (i,j,l+m,d) + z (i,j,k,d) * w
        endif
      end do
    end do
    end do
    end do
    end do
  end subroutine apply_vint_t

!------------------------------------------------------------------------------

  subroutine apply_vloc_t (c1, z, v)
  type (t_c1) ,intent(in)  :: c1
  real(wp)    ,intent(in)  :: z (c1%lb(1):,c1%lb(2):,:,:)
  real(wp)    ,intent(out) :: v (c1%lb(1):,c1%lb(2):,:,:)
#ifdef _CRAYFTN
  contiguous :: z, v
#endif

    integer  :: i,j,k,d,l,m,n,u
    real(wp) :: w

    n = size (c1% cg, 1)
    u = size (v, 3)
    v = 0._wp
    do d = lbound(z,4), ubound(z,4)
    do k = lbound(z,3), ubound(z,3)
    do j = lbound(z,2), ubound(z,2)
#ifdef __NEC__
      do l = 1, n
!NEC$ ivdep
    do i = lbound(z,1), ubound(z,1)
      m = c1% l0 (i,j,k,d)
        w = c1% cg (l,i,j,k,d)
        if ( (1 <= l+m) .and. (l+m <= u) ) &
                        v (i,j,l+m,d) = v (i,j,l+m,d) + z (i,j,k,d) * w
      end do
    end do
#else
    do i = lbound(z,1), ubound(z,1)
      m = c1% l0 (i,j,k,d)
      do l = 1, n
        w = c1% cg (l,i,j,k,d)
        if (w /= 0._wp) v (i,j,l+m,d) = v (i,j,l+m,d) + z (i,j,k,d) * w
      end do
    end do
#endif
    end do
    end do
    end do
  end subroutine apply_vloc_t

!------------------------------------------------------------------------------

  subroutine apply_vrint (c1, z, v)
  type (t_c1) ,intent(in)  :: c1
  real(wp)    ,intent(in)  :: z (c1%lb(1):,c1%lb(2):,:,:)
  real(wp)    ,intent(out) :: v (c1%lb(1):,c1%lb(2):,:,:)

    integer  :: i,j,k,d,l,m,n
    real(wp) :: w

    n = size (c1% cr, 1)
    v = 0._wp
    do d = lbound(z,4), ubound(z,4)
    do k = lbound(z,3), ubound(z,3)
    do j = lbound(z,2), ubound(z,2)
    do i = lbound(z,1), ubound(z,1)
      m = c1% l0 (i,j,k,d)
      do l = 1, n
        w = c1% cr (l,i,j,k,d)
        if (w /= 0._wp) v (i,j,l+m,d) = v (i,j,l+m,d) + z (i,j,k,d) * w
      end do
    end do
    end do
    end do
    end do
  end subroutine apply_vrint

!------------------------------------------------------------------------------

  subroutine test_c1_modes
  !-----------------------------------------------------------------
  ! test accuracy of SQRT factorisation compared to Gaspari & Cohn
  ! parameters tested: stride: 0.1 .. 2. * localisation length scale
  !                    shift : 0 ... stride
  !-----------------------------------------------------------------
    integer            :: stride
    integer            :: i, j, k, l
    integer ,parameter :: n = 2*3*2*5*7*2*3
    type (t_c1)        :: c1
    real(wp)           :: ref (0:n-1, 0:19)
    real(wp)           :: x   (0:n-1)
    real(wp)           :: norm
    real(wp)           :: dev
    real(wp)           :: asym

    integer     :: m (13) = &
                   (/1,2,3,4,5,6,7,8,9,10,12,15,20/)

    call construct (c1, n, 10._wp, cyclic=.true., lsqrt=.false.)
    do i = 0,19
      k = n/2 + i
      x = 0._wp; x(k) = 1._wp
      call apply_c1_c1t (c1, x)
      ref (:,i) = x
    end do
    norm = sum (x)
    call destruct (c1)

    do l = 1, 13
      stride =  m (l)
      call construct (c1, n, 10._wp, cyclic=.true., lsqrt=.true., &
                      stride=stride                               )
      do i = 0, stride - 1
      k = n/2 + i
        x = 0._wp; x(k) = 1._wp
        call apply_c1_c1t (c1, x)
        dev  = sum ( abs (x - ref (:,i))) / norm
        asym = 0._wp
        do j = 1, n/4
          asym = asym + abs (x (k-j) - x (k+j))
        end do
        asym = asym / norm
        print *,'stride =',stride,' offset =',i,' error =',dev,' asym =',asym
      end do
      print *
      call destruct (c1)
    end do

  end subroutine test_c1_modes

!------------------------------------------------------------------------------

  subroutine test_c1_bound

    integer ,parameter :: n = 100
    integer            :: i, j
    real(wp)           :: ref (0:n-1)
    real(wp)           :: x   (0:n-1)
    real(wp)           :: y   (0:n-1)
    type (t_c1)        :: c1
    real(wp)           :: norm
    real(wp)           :: dev
    real(wp)           :: asym

    call construct (c1, n, 10._wp, cyclic=.false., lsqrt=.false.)

    x   = 0._wp; x(n/2) = 1._wp
    ref = x
    call apply_c1_c1t (c1, ref)
    norm = sum (ref)

    do i = 0,n-1
      call construct (c1, n, 10._wp, cyclic=.false., lsqrt=.true., &
                      stride=10, shift = i                         )
      y = x
      call apply_c1_c1t (c1, y)
      dev  = sum ( abs (y - ref)) / norm
      asym = 0._wp
      do j = 1, n/2-1
        asym = asym + abs (y (n/2-j) - y (n/2+j))
      end do
      asym = asym / norm
      print *,'stride =',10,' offset =',i,' error =',dev,' asym =',asym
      call destruct (c1)
    end do

  end subroutine test_c1_bound

!==============================================================================
end module mo_varenkf_1d
