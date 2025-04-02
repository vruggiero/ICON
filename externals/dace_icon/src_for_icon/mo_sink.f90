!
!+ handling of (bound constrained) dummy sink variables
!
! $Id$
!
MODULE mo_sink
!
! Description:
!    Handling of (bound constrained) dummy sink variables.
!
! Current Code Owner: DWD, Andreas Rhodin
!    phone: +49 69 8062 2722
!    fax:   +49 69 8062 3721
!    email: andreas.rhodin@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_28        2014/02/26 Andreas Rhodin
!  New module for bound constrained dummy sink variables
! V1_35        2014-11-07 Andreas Rhodin
!  modify dummy sink variable error in case of logarithmic transform;
!  routines to account for Jakobian of sink variable transformation
! V1_51        2017-02-24 Andreas Rhodin
!  change meaning of sink variable mode (1/2 log/sqrt transform)
!  option to account for nonlinearity of generalised humidity in preconditioner
!
! Code Description:
! Language: Fortran 95.
! Software Standards:
!
!==============================================================================
  !=============
  ! modules used
  !=============
  use mo_kind,       only: wp               ! working precision kind parameter
  implicit none
  !================
  ! public entities
  !================
  private
  public :: t_sink     ! dummy sink variable meta data type
  public :: set_v_sink ! set derived quantities in t_sink
  !------------------------------------------------------
  ! Routines for bound constrained optimisation:
  ! Method: variable transformation
  !      x: bound constrained physical quantity
  !      v: transformed variable
  !      d: normalised transformed variable
  !     dx: dx  / dv  (or dx  / dd)   first derivative
  !    d2x: d2x / dv2 (or d2x / dd2) second derivative
  ! v is identical to x in the center of the domain
  ! but differs in the vicinity of the bounds.
  ! Additional transformations relate x to the normalised
  ! dummy sink variable d used in the psas minimisation.
  ! The relation between v and d is linear.
  !------------------------------------------------------
  public :: bc_xd  ! calculate x from (normally distributed) control variable
  public :: bc_xv  ! calculate x, dx  (optionally) from v
  public :: bc_vx  ! calculate v from x  (used for testing only)
  public :: bc_bv  ! derive B matrix accounting for change of Jakobian
  public :: bc_bd  ! derive B matrix accounting for change of Jakobian
  !========================
  ! derived type definition
  !========================
  type t_sink
    !------------------------------------------------------------
    ! specifications for (bound constrained) dummy sink variables
    !------------------------------------------------------------
    integer  :: ispot = -1            ! reference to report / FOV
    integer  :: iobs  = -1            ! reference to observation in report
    integer  :: mode  =  0            ! 0:no, 1:log transform, 2:quadratic
    real(wp) :: bg    =  0._wp        ! background value
    real(wp) :: err   =  1._wp        ! background error
    real(wp) :: min   = -huge(1._wp)  ! lower bound
    real(wp) :: max   =  huge(1._wp)  ! upper bound
    real(wp) :: lbnd  =  0._wp        ! nonlinear range
    real(wp) :: ubnd  =  0._wp        ! nonlinear range
    !-------------------
    ! derived quantities
    !-------------------
    real(wp) :: bg_v  =  0._wp        ! background value, transformed variable
    real(wp) :: err_v =  1._wp        ! background error, transformed variable
  end type t_sink

!==============================================================================
 contains
!==============================================================================
!
! Routines for bound constrained optimisation:
! Method: variable transformation
!      x: bound constrained parameter space
!      v: transformed parameter space
!     dx: dx / dv
!  bc_xv: calculate x, dx (optional) from v
!  bc_vx: calculate v from x
!  bc_xd: calculate x from (normally distributed) control variable
!  for x==xmin or x==xmax dx/dv vanishes
!  for xminr<x<xmaxr      dx/dv = 1
!  bound_constr = 0: no trensformation
!                 1: exponential/logarithmic transform
!                 2: quadratic/sqrt          transform
!==============================================================================

  subroutine bc_xv (x, v, xmin, xminr, xmaxr, xmax, bound_constr, dx, d2x, d2d)
  real(wp) ,intent(out)           :: x             ! physical quantity
  real(wp) ,intent(in)            :: v             ! transformed variable
  real(wp) ,intent(in)            :: xmin          ! lower bound for x
  real(wp) ,intent(in)            :: xminr         ! nonlinear range at lb
  real(wp) ,intent(in)            :: xmaxr         ! nonlinear range at ub
  real(wp) ,intent(in)            :: xmax          ! upper bound for x
  integer  ,intent(in)            :: bound_constr  ! no/sqrt/log transform
  real(wp) ,intent(out) ,optional :: dx            ! Jakobian
  real(wp) ,intent(out) ,optional :: d2x           ! second derivative
  real(wp) ,intent(out) ,optional :: d2d           ! second derivative / first
  !--------------------------------------------------------
  ! calculate physical quantity x from transformed variable
  !--------------------------------------------------------
    real(wp) :: r

    select case (bound_constr)
    !------------------
    ! no transformation
    !------------------
    case (0)
      x = v
      if (present(dx)) dx = 1._wp
    !------------------------------
    ! quadratic/sqrt transformation
    !------------------------------
    case (2)
      if      (v <= xmin-xminr) then
        x = xmin
        if (present(dx )) dx  = 0._wp
        if (present(d2x)) d2x = 0._wp
        if (present(d2d)) d2d = huge(0._wp)
      else if (v < xmin+xminr) then
        r = (v-xmin+xminr)/(2._wp*xminr)
        x = xmin + xminr * r**2
        if (present(dx )) dx  = r
        if (present(d2x)) d2x = 0.5_wp/xminr
        if (present(d2d)) d2d = 1._wp / (v-xmin+xminr)
      else if (v <= xmax-xmaxr) then
        x = v
        if (present(dx )) dx  = 1._wp
        if (present(d2x)) d2x = 0._wp
        if (present(d2d)) d2d = 0._wp
      else if (v < xmax+xmaxr) then
        r = (xmax+xmaxr-v)/(2._wp*xmaxr)
        x = xmax - xmaxr * r**2
        if (present(dx )) dx  = r
        if (present(d2x)) d2x = - 0.5_wp/xmaxr
        if (present(d2d)) d2d = - 1._wp / (xmax+xmaxr-v)
      else
        x = xmax
        if (present(dx )) dx  = 0._wp
        if (present(d2x)) d2x = 0._wp
        if (present(d2d)) d2d = - huge(0._wp)
      endif
    !---------------------------------------
    ! exponential/logarithmic transformation
    !---------------------------------------
    case (1)
      if      (v < xmin+xminr) then
        r = exp ((v-xmin) / xminr - 1._wp)
        x = xmin + xminr * r
        if (present(dx )) dx  = r
        if (present(d2x)) d2x = r / xminr
        if (present(d2d)) d2d = 1._wp / xminr
      else if (v <= xmax-xmaxr) then
        x = v
        if (present(dx )) dx  = 1._wp
        if (present(d2x)) d2x = 0._wp
        if (present(d2d)) d2d = 0._wp
      else
        r = exp ((xmax-v) / xmaxr - 1._wp)
        x = xmax - xmaxr * r
        if (present(dx )) dx  =   r
        if (present(d2x)) d2x = - r / xmaxr
        if (present(d2d)) d2d = - 1._wp / xmaxr
      endif
    end select

  end subroutine bc_xv
!------------------------------------------------------------------------------
  subroutine bc_vx (x, v, xmin, xminr, xmaxr, xmax, bound_constr)
  real(wp) ,intent(in)           :: x             ! physical quantity
  real(wp) ,intent(out)          :: v             ! transformed variable
  real(wp) ,intent(in)           :: xmin          ! lower bound
  real(wp) ,intent(in)           :: xminr         ! nonlinear range at lb
  real(wp) ,intent(in)           :: xmaxr         ! nonlinear range at ub
  real(wp) ,intent(in)           :: xmax          ! upper bound
  integer  ,intent(in)           :: bound_constr  ! kind of transformation

    real(wp) :: r
    real(wp) ,parameter :: epsr = 1.e-3_wp

    select case (bound_constr)
    case (0)
      v = x
    case (2)
      if      (x < xmin) then
        v = xmin - xminr
      else if (x < xmin+xminr) then
        r = (x-xmin)/xminr
        v = xmin - xminr + 2._wp * xminr * sqrt(r)
      else if (x < xmax-xmaxr) then
        v = x
      else if (x < xmax) then
        r = (xmax-x)/xmaxr
        v = xmax + xmaxr - 2._wp * xmaxr * sqrt(r)
      else
        v = xmax + xmaxr
      endif
    case (1)
      if      (x < xmin+xminr) then
        r = max (epsr, (x - xmin) / xminr)
        v = (log (r) + 1._wp) * xminr + xmin
      else if (x < xmax-xmaxr) then
        v = x
      else
        r = max (epsr, (xmax - x) / xmaxr)
        v = xmax - (log (r) + 1._wp) * xmaxr
      endif
    end select

  end subroutine bc_vx
!------------------------------------------------------------------------------
  subroutine bc_xd (x, d, s, dx, d2x, d2d)
  !-----------------------------------------------------------
  ! calculate x from (normally distributed) control variable d
  ! x: bound constrained physical quantity
  ! d: control variable used in the variational scheme with
  !    expectation value (first guess)       = 0.
  !    uncertainty       (first guess error) = 1.
  !-----------------------------------------------------------
  real(wp)      ,intent(out)           :: x   ! physical quantity
  real(wp)      ,intent(in)            :: d   ! normalised control variable
  type (t_sink) ,intent(in)            :: s   ! sink variable meta data
  real(wp)      ,intent(out) ,optional :: dx  ! Jakobian
  real(wp)      ,intent(out) ,optional :: d2x ! second derivativ
  real(wp)      ,intent(out) ,optional :: d2d ! second derivative / first

    real(wp) :: v

    v = s% bg_v + s% err_v * d
    call bc_xv (x, v, s% min, s% lbnd, s% ubnd, s% max, s% mode, &
                                          dx=dx, d2x=d2x, d2d=d2d)
    if (present (dx )) dx  = dx  * s% err_v
    if (present (d2x)) d2x = d2x * s% err_v ** 2
    if (present (d2d)) then
      if (d2d > -huge(0._wp) .and. d2d < huge(0._wp)) d2d = d2d * s% err_v
    endif

  end subroutine bc_xd
!------------------------------------------------------------------------------
  function bc_bd (d, s)
  !------------------------------------------------------------------
  ! Improved estimate of B in normalised control space (usually == 1)
  ! to be used in preconditioning,
  ! taking into account the Jakobian of the transformation x(v).
  !------------------------------------------------------------------
  real(wp)                   :: bc_bd  ! B estimate
  real(wp)      ,intent(in)  :: d      ! normalised control variable
  type (t_sink) ,intent(in)  :: s      ! sink variable meta data

    real(wp)            :: v              ! transformed variable
    real(wp)            :: x              ! physical quantity
    real(wp)            :: d2d            ! d2x/dd2 / dx/dd
    real(wp)            :: b              ! 1 / B
    real(wp) ,parameter :: edbmax = 1._wp ! bound for 1 / B

    v = s% bg_v + s% err_v * d

    if (v < s% min + s% lbnd .or. v > s% max - s% ubnd) then
      call bc_xd (x, d, s, d2d=d2d)
      if (d2d > -huge(0._wp) .and. d2d < huge(0._wp)) then
        b     = 1._wp - d * d2d
        bc_bd = 1._wp / max (b, edbmax)
      else
        bc_bd = 0._wp
      endif
    else
      bc_bd = 1._wp
    endif

  end function bc_bd
!------------------------------------------------------------------------------
  function bc_bv (v, vb, b, s)
  !------------------------------------------------------------------
  ! Improved estimate of B
  ! in case of nonlinear transformations of control variables.
  ! to be used in preconditioning,
  ! taking into account the Jakobian of the transformation x(v).
  !------------------------------------------------------------------
  real(wp)                   :: bc_bv  ! B estimate
  real(wp)      ,intent(in)  :: v      ! transformed control variable
  real(wp)      ,intent(in)  :: vb     ! control variable first guess
  real(wp)      ,intent(in)  :: b      ! control variable error
  type (t_sink) ,intent(in)  :: s      ! sink variable meta data

    real(wp)            :: x              ! physical quantity
    real(wp)            :: d2d            ! d2x/dd2 / dx/dd
    real(wp)            :: bi             ! 1 / B
    real(wp) ,parameter :: edbmax = 1._wp ! bound for 1 / B

    if (v < s% min + s% lbnd .or. v > s% max - s% ubnd) then
      call bc_xv (x, v, s% min, s% lbnd, s% ubnd, s% max, s% mode, d2d=d2d)
      if (d2d > -huge(0._wp) .and. d2d < huge(0._wp)) then
        bi    = 1._wp / b
        bi    = bi - (v-vb) * bi * d2d
        bc_bv = 1._wp / bi
      else
        bc_bv = 0._wp
      endif
    else
      bc_bv = b
    endif

  end function bc_bv
!==============================================================================
  subroutine set_v_sink (s)
  type(t_sink) ,intent(inout) :: s
  !---------------------------------
  ! set derived quantities in t_sink
  !---------------------------------
    real(wp) :: dx, dummy
    call bc_vx (s% bg, s% bg_v, s% min, s% lbnd, s% ubnd, s% max, s% mode)
    call bc_xv (dummy, s% bg_v, s% min, s% lbnd, s% ubnd, s% max, s% mode, dx)
    select case (s% mode)
!   case (2)
!     s% err_v = s% err / dx
    case default
      s% err_v = s% err
    end select
  end subroutine set_v_sink
!==============================================================================
end module mo_sink
