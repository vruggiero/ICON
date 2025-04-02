!
!+ generalised humidity transformations
!
MODULE mo_cntrlvar
!
! Description:
!   This module holds routines to transform from control variable
!   space (tv, virtual temperature and gh, generalized humidity) to
!   physical space (t, temperature and q, specific humidity)
!
!   Generalized humidity:
!   Above rh0 (default: 3%) gh is equal to rh (relative humidity over water).
!   version = 1 (exponential):
!     Below rh0: rh          = rh0/e * exp (gh / rh0)
!                d rh / d gh =   1/e * exp (gh / rh0)    1st derivative
!                            = 1/rh0 * rh
!                gh          = rh0   * ln  (rh * e/rh0)  inverse
!   version = 2 (quadratic):
!     Below rh0: rh          = rh0 * ((gh + rh0) / 2 rh0) ** 2
!                d rh / d gh =        (gh + rh0) / 2 rh0
!   Similar behaviour at saturation value.
!
! Current Maintainer: DWD, Robin Faulwetter
!    phone: +49 69 8062 2746
!    fax:   +49 69 8062 3721
!    email: robin.faulwetter@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_1         2008/11/05 Andreas Rhodin
!  First operational 3D-Var release
! V1_5         2009/05/25 Harald Anlauf
!  tq_tv_gh_vec: vectorized version of tq_tv_gh
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Andreas Rhodin
!  remove unused variables
! V1_26        2013/06/27 Andreas Rhodin
!  new version 2: quadratic transformation (not exponential)
! V1_51        2017-02-24 Andreas Rhodin
!  option to account for nonlinearity of generalised humidity in preconditioner
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Andreas Rhodin  DWD  2004 split from mo_physics, generalized humidity
! Detlef  Pingel  DWD  2004 new subroutine tq_tvrh
! Andreas Rhodin  DWD  2008 modification at saturation value
! Harald Anlauf   DWD  2008 preset saved variables, change printout
!=====================================================================

  !=============
  ! modules used
  !=============
  use mo_kind,     only: wp               ! working precision kind parameter
  use mo_exception,only: finish           ! abort on error condition
  use mo_namelist, only: position_nml,   &! routine to position nml group
                         nnml,           &! namelist fortran unit number
                         POSITIONED       ! position_nml: OK    return flag
  use mo_mpi_dace, only: dace,           &! MPI group info
                         p_bcast          ! broadcast routine
  use mo_sink,     only: t_sink           ! meta data derived type
  use mo_physics,  only: tq_tvrh,        &! conversion: t, q  <- tv, rh
#ifndef __NEC__
                         tq_tvrh_vec,    &!
#else
                         tq_tvrh_vec2,   &!
#endif
                         tq_tvrh_noxfail,&!
                         e                ! exp(1)
  implicit none
  !================
  ! public entities
  !================
  private
  !------------------------------
  ! conversion routines, adjoints
  !------------------------------
  public :: read_cntrlvar_nml ! read namelist /CNTRLVAR/
  public :: gh_rh             ! derive general. humidity from relative humidity
  public :: rh_gh             ! derive relative humidity from general. humidity
  public :: tq_tvgh           ! temp., spec.hum. <- virt.temp., gen.hum.
  public :: tq_tvgh_vec       !         same, but vectorization-adapted
  public :: tq_tvgh_adj_vec   !         same, but vectorization-adapted
  public :: trh_tvgh          ! temp.,  rel.hum. <- virt.temp., gen.hum.
  public :: rh0               ! <=0 indicates no generalized humidity
  public :: gh_meta           ! generalised humidity parameters
  public :: nl_prec           ! account for nonlinearity in preconditioner
  public :: gh_min            ! minimum value for generalized humidity
  public :: gh_max            ! maximum value for generalized humidity
  public :: disable_gh        ! disable generalized humidity transformation

#ifdef __NEC__
  ! Explicitly export to Work around inlining issue with NEC nfort:
  public :: a0,a1,edrh0,edb1,version,rhmax,edrh02,rh1,edb2,b1,b3
#endif


  integer,  parameter :: mlay             = 3          ! Max. number of rh0fg/rh1fg definitions
  real(wp), parameter :: r_init           = HUGE (999._wp)
  !-------------------
  ! namelist variables
  !-------------------
  real(wp)            :: rh0              = 0.03_wp    ! lower bound for linear range (rel.hum)
  real(wp)            :: rh0fg     (mlay) = 0.01_wp    ! lower bound for first guess rel.humidity
  real(wp)            :: rh0fg_ptop(mlay) = 0._wp      ! pressure until which rh0fg is applied
  real(wp)            :: q0fg             = 1.25e-6_wp ! lower bound for first guess spc.humidity
  real(wp)            :: rhmax            = r_init     ! upper bound for relative humidity
  real(wp)            :: rh1              = r_init     ! upper bound for linear range (rel.hum)
  real(wp)            :: rh1fg     (mlay) = r_init     ! upper bound for first guess rel.humidity
  real(wp)            :: rh1fg_ptop(mlay) = 0._wp      ! pressure until which rh0fg is applied
  integer             :: version          = 1          ! 1:exponential; 2:quadratic
  logical             :: warn             = .false.    ! print warnings if gen.humidity bounds are exceeded
  logical             :: nl_prec          = .false.    ! account for nonlinearity in preconditioner
  namelist /CNTRLVAR/ rh0, rh0fg, rh0fg_ptop, q0fg, rhmax, rh1, rh1fg, rh1fg_ptop, version, warn, nl_prec
  !-----------------
  ! module variables
  !-----------------
  integer             :: nlay0            =  0         ! Number of layers with rh0fg threshold
  integer             :: nlay1            =  0         ! Number of layers with rh1fg threshold
  real(wp)            :: a0               =  0._wp        ! rh0 / e
  real(wp)            :: edrh0            =  HUGE (0._wp) ! 1 / rh0
  real(wp)            :: rh02             =  0._wp        ! 2 * rh0
  real(wp)            :: edrh02           =  HUGE (0._wp) ! 1 / 2*rh0
  real(wp)            :: a1               =  0._wp        ! b1 / e
  real(wp)            :: b1               =  0._wp        ! rhmax - rh1
  real(wp)            :: b2               =  0._wp        ! 2 * b1
  real(wp)            :: b3               =  100._wp      ! rhmax + b1  TODO:999._wp ??
  real(wp)            :: edb1             =  HUGE (0._wp) ! 1 / b1
  real(wp)            :: edb2             =  HUGE (0._wp) ! 1 / b2
  real(wp)            :: gh_max           =  huge (0._wp)
  real(wp)            :: gh_min           = -huge (0._wp)
  type(t_sink),  save :: gh_meta          ! meta data
  !-----------
  ! interfaces
  !-----------
  interface gh_rh
    module procedure gh_rh
    module procedure gh_rh_tp
  end interface gh_rh

  interface rh_gh
    module procedure rh_gh_
    module procedure rh_gh_d
  end interface rh_gh

  interface tq_tvgh              !+++
    module procedure tq_tvgh     !+++ work around NEC SX bug:
    module procedure tq_tvgh_adj !+++ optional, intent(out) in elemental
  end interface tq_tvgh          !+++ subroutine

!  interface tq_tvgh_vec
!    module procedure tq_tvgh_vec
!    module procedure tq_tvgh_adj_vec
!  end interface tq_tvgh_vec

contains
!==============================================================================

  elemental function gh_rh (rh, fg) result (gh)
  !--------------------------------------------------------------------------
  ! derive generalized humidity from relative humidity,
  ! used after interpolation of first guess to the locations of observations.
  !--------------------------------------------------------------------------
  real(wp)             :: gh ! generalized humidity coordinate
  real(wp) ,intent(in) :: rh ! relative humidity
  logical  ,intent(in) :: fg ! first guess flag, account for rh0fg if .true.
  integer              :: lay

    gh = rh
    !---------------------
    ! restrict first guess
    !---------------------
    if (fg) then
      do lay = 1,nlay0
        gh = max(gh, rh0fg(lay))
      end do
      do lay = 1,nlay1
        gh = min(gh, rh1fg(lay))
      end do
    end if

    select case (version)
    case (1)
      if      (gh < rh0) then
        if (gh > 0._wp) then
          gh = rh0 * log (gh / a0)
        else
          gh = -huge(0._wp)
        end if
      else if (gh > rh1) then
        if (gh < rhmax) then
          gh = rhmax - gh
          gh = b1 * log (gh / a1)
          gh = rhmax - gh
        else
          gh = huge(0._wp)
        end if
      endif
    case (0,2)
      gh = max (min (gh, rhmax), 0._wp)
      if      (gh < rh0) then
        gh = rh02 * sqrt (gh*edrh0) - rh0
      else if (gh > rh1) then
        gh = b3 - b2 * sqrt ((rhmax-gh)*edb1)
      endif
    end select
  end function gh_rh
!------------------------------------------------------------------------------
  elemental function gh_rh_tp (rh, fg, t, p) result (gh)
  !--------------------------------------------------------------------------
  ! derive generalized humidity from relative humidity,
  ! used after interpolation of first guess to the locations of observations.
  !--------------------------------------------------------------------------
  real(wp)             :: gh ! generalized humidity coordinate
  real(wp) ,intent(in) :: rh ! relative humidity
  logical  ,intent(in) :: fg ! first guess flag, account for rh0fg if .true.
  real(wp) ,intent(in) :: t  ! temperature
  real(wp) ,intent(in) :: p  ! pressure
  integer              :: lay

    gh = rh
    !---------------------
    ! restrict first guess
    !---------------------
    if (fg) then
      do lay = 1,nlay0
        if (p > rh0fg_ptop(lay)) then
          gh = max(gh, rh0fg(lay))
          exit
        end if
      end do
      do lay = 1,nlay1
        if (p > rh1fg_ptop(lay)) then
          gh = min(gh, rh1fg(lay))
          exit
        end if
      end do
    end if

    select case (version)
    !------------------------
    ! exponential formulation
    !------------------------
    case (1)
      !---------
      ! gh <- rh
      !---------
      if (gh < rh0) then
        if (gh > 0._wp) then
          gh = rh0 * log (gh / a0)
        else
          gh = -huge(0._wp)
        end if
      else if (gh > rh1) then
        if (gh < rhmax) then
          gh = rhmax - gh
          gh = b1 * log (gh / a1)
          gh = rhmax - gh
        else
          gh = huge(0._wp)
        end if
      endif
    !----------------------
    ! quadratic formulation
    !----------------------
    case (0,2)
      gh = max (min (gh, rhmax), 0._wp)
      if      (gh < rh0) then
        gh = rh02 * sqrt (gh*edrh0) - rh0
      else if (gh > rh1) then
        gh = b3 - b2 * sqrt ((rhmax-gh)*edb1)
      endif
    end select
  end function gh_rh_tp

!------------------------------------------------------------------------------
  elemental subroutine rh_gh_ (rh, gh)
  real(wp) ,intent(out)           :: rh     ! relative humidity
  real(wp) ,intent(in)            :: gh     ! generalized humidity
  !---------------------------------------------------
  ! derive relative humidity from generalized humidity
  !---------------------------------------------------
    select case (version)
    case (1)
      if (gh < rh0) then
        rh = a0 * exp (gh * edrh0)
      else if (gh > rh1) then
        rh = rhmax - gh
        rh = a1 * exp (rh * edb1)
        rh = rhmax - rh
      else
        rh = gh
      endif
    case (0,2)
      if (gh < rh0) then
        rh =        rh0 * (edrh02 * max (gh + rh0, 0._wp)) ** 2
      else if (gh > rh1) then
        rh = rhmax - b1 * (edb2   * max (b3 - gh,  0._wp)) ** 2
      else
        rh = gh
      endif
    end select
  end subroutine rh_gh_
!------------------------------------------------------------------------------
  elemental subroutine rh_gh_d (rh, gh, drh_gh)
  real(wp) ,intent(out)           :: rh     ! relative humidity
  real(wp) ,intent(in)            :: gh     ! generalized humidity
  real(wp) ,intent(out)           :: drh_gh ! d rh / d gh
  !---------------------------------------------------
  ! derive relative humidity from generalized humidity
  ! optionally derive derivative
  !---------------------------------------------------
    select case (version)
    case (1)
      if (gh < rh0) then
        rh = a0 * exp (gh * edrh0)
        drh_gh = rh * edrh0
      else if (gh > rh1) then
        rh = rhmax - gh
        rh = a1 * exp (rh * edb1)
        drh_gh = rh * edb1
        rh = rhmax - rh
      else
        rh = gh
        drh_gh = 1._wp
      endif
    case (0,2)
      if (gh < rh0) then
        rh     = edrh02 * max (gh + rh0, 0._wp)
        drh_gh = rh
        rh     = rh0 * rh ** 2
      else if (gh > rh1) then
        rh     = edb2   * max (b3 - gh,  0._wp)
        drh_gh = rh
        rh     = rhmax - b1 * rh ** 2
      else
        rh = gh
        drh_gh = 1._wp
      endif
    end select
  end subroutine rh_gh_d
!------------------------------------------------------------------------------
  elemental subroutine trh_tvgh (t, rh, tv, gh, p, i_fail,  &
                                 dt_tv, dt_gh, drh_tv, drh_gh)
  !-----------------------------------------------------------------
  ! Transform virtual temperature (tv) and generalized humidity (gh)
  !                to temperature (t)  and specific humidity (q).
  ! Provide the Jacobi-matrix elements as well.
  !-----------------------------------------------------------------
  real(wp) ,intent (out)            :: t      ! temperature
  real(wp) ,intent (out)            :: rh     ! relative humidity
  real(wp) ,intent (in)             :: tv     ! virtual temperature
  real(wp) ,intent (in)             :: gh     ! generalized humidity
  real(wp) ,intent (in)             :: p      ! pressure
  integer  ,intent (out)            :: i_fail ! >0: number of iterations
                                              ! <0: no convergence
  real(wp) ,intent (out) ,optional  :: dt_tv  ! d t  / d tv
  real(wp) ,intent (out) ,optional  :: dt_gh  ! d t  / d gh
  real(wp) ,intent (out) ,optional  :: drh_tv ! d rh / d tv
  real(wp) ,intent (out) ,optional  :: drh_gh ! d rh / d gh
    !----------------
    ! local variables
    !----------------
    real(wp) :: q      ! specific
    real(wp) :: drhgh  ! d rh / d gh ,local variable
    real(wp) :: dt_rh  ! d t  / d gh
    !-------------------
    ! transform gh -> rh
    !-------------------
    call rh_gh (rh, gh, drhgh)
    !-----------------------
    ! transform tv,rh -> t,q
    !-----------------------
    call tq_tvrh (t, q, tv, rh, p, i_fail=i_fail, dt_tv=dt_tv, dt_rh=dt_rh)
    !-----------------------------
    ! adapt Jacobi-matrix elements
    !-----------------------------
    if (present (dt_gh )) dt_gh  = dt_rh  * drhgh ! d t  / d rh
    if (present (drh_tv)) drh_tv = 0._wp          ! d rh / d tv
    if (present (drh_gh)) drh_gh =          drhgh ! d rh / d gh
  end subroutine trh_tvgh
!------------------------------------------------------------------------------
  elemental subroutine tq_tvgh_adj (t, q, tv, gh, p, i_fail, &
                                    dt_tv, dt_gh, dq_tv, dq_gh, drh_gh)
  !-----------------------------------------------------------------
  ! Transform virtual temperature (tv) and generalized humidity (gh)
  !                to temperature (t)  and specific humidity (q).
  ! Provide the Jacobi-matrix elements as well.
  !-----------------------------------------------------------------
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! VERSION WITH ADJOINT
  ! NEC SX cannot handle optional intent(out) in elemental subroutine
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  real(wp) ,intent (out)            :: t      ! temperature
  real(wp) ,intent (out)            :: q      ! specific humidity
  real(wp) ,intent (in)             :: tv     ! virtual temperature
  real(wp) ,intent (in)             :: gh     ! generalized humidity
  real(wp) ,intent (in)             :: p      ! pressure
  integer  ,intent (out)            :: i_fail ! >0: number of iterations
                                              ! <0: no convergence
! real(wp) ,intent (out) ,optional  :: dt_tv  ! part. deriv. dt/dtv
! real(wp) ,intent (out) ,optional  :: dt_gh  ! part. deriv. dt/dgh
! real(wp) ,intent (out) ,optional  :: dq_tv  ! part. deriv. dq/dtv
! real(wp) ,intent (out) ,optional  :: dq_gh  ! part. deriv. dq/dgh
  real(wp) ,intent (out)            :: dt_tv  ! part. deriv. dt/dtv
  real(wp) ,intent (out)            :: dt_gh  ! part. deriv. dt/dgh
  real(wp) ,intent (out)            :: dq_tv  ! part. deriv. dq/dtv
  real(wp) ,intent (out)            :: dq_gh  ! part. deriv. dq/dgh
  real(wp) ,intent (out) ,optional  :: drh_gh
    !----------------
    ! local variables
    !----------------
    real(wp) :: rh    ! relative humidity
    real(wp) :: drhgh ! d rh / d gh
    !-------------------
    ! transform gh -> rh
    !-------------------
    call rh_gh (rh, gh, drhgh)
    !-----------------------
    ! transform tv,rh -> t,q
    !-----------------------
    call tq_tvrh (t, q, tv, rh, p, i_fail=i_fail, &
                  dt_tv=dt_tv, dt_rh=dt_gh, dq_tv=dq_tv, dq_rh=dq_gh)
    !-----------------------------
    ! adapt Jacobi-matrix elements
    !-----------------------------
!   if (present (dt_gh)) dt_gh = dt_gh * drhgh
!   if (present (dq_gh)) dq_gh = dq_gh * drhgh
    dt_gh = dt_gh * drhgh
    dq_gh = dq_gh * drhgh
    if (present(drh_gh)) drh_gh = drhgh
  end subroutine tq_tvgh_adj
!------------------------------------------------------------------------------
  elemental subroutine tq_tvgh (t, q, tv, gh, p, i_fail)
  !-----------------------------------------------------------------
  ! Transform virtual temperature (tv) and generalized humidity (gh)
  !                to temperature (t)  and specific humidity (q).
  ! Provide the Jacobi-matrix elements as well.
  !-----------------------------------------------------------------
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! VERSION WITHOUT ADJOINT
  ! NEC SX cannot handle optional intent(out) in elemental subroutine
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  real(wp) ,intent (out)            :: t      ! temperature
  real(wp) ,intent (out)            :: q      ! specific humidity
  real(wp) ,intent (in)             :: tv     ! virtual temperature
  real(wp) ,intent (in)             :: gh     ! generalized humidity
  real(wp) ,intent (in)             :: p      ! pressure
  integer  ,intent (out)            :: i_fail ! >0: number of iterations
                                              ! <0: no convergence
    !----------------
    ! local variables
    !----------------
    real(wp) :: rh    ! relative humidity
    !-------------------
    ! transform gh -> rh
    !-------------------
    call rh_gh (rh, gh)
    !-----------------------
    ! transform tv,rh -> t,q
    !-----------------------
    call tq_tvrh_noxfail (t, q, tv, rh, p, 0._wp,0._wp,i_fail=i_fail)
  end subroutine tq_tvgh
!------------------------------------------------------------------------------
  subroutine tq_tvgh_vec (n, t, q, tv, gh, p, i_fail)
  !-----------------------------------------------------------------
  ! Vectorization-adapted version of tq_tvgh
  !-----------------------------------------------------------------
  integer  ,intent(in)           :: n         ! Number of gridpoints
  real(wp) ,intent(out)          :: t     (n) ! temperature
  real(wp) ,intent(out)          :: q     (n) ! specific humidity
  real(wp) ,intent(in)           :: tv    (n) ! virtual temperature
  real(wp) ,intent(in)           :: gh    (n) ! generalized humidity
  real(wp) ,intent(in)           :: p     (n) ! pressure
  integer  ,intent(out)          :: i_fail(n) ! >0: number of iterations
                                              ! <0: no convergence
    !-------------
    ! local arrays
    !-------------
    real(wp) :: rh(n)    ! relative humidity
    !-------------------
    ! transform gh -> rh
    !-------------------
    call rh_gh (rh, gh)
    !-----------------------
    ! transform tv,rh -> t,q
    !-----------------------
#ifdef __NEC__
    call tq_tvrh_vec2 &
#else
    call tq_tvrh_vec  &
#endif
                     (n, t, q, tv, rh, p, i_fail=i_fail)
  end subroutine tq_tvgh_vec
!------------------------------------------------------------------------------
  subroutine tq_tvgh_adj_vec (n, t, q, tv, gh, p, i_fail, &
       dt_tv, dt_gh, dq_tv, dq_gh, drh_gh)
  !-----------------------------------------------------------------
  ! Vectorization-adapted version of tq_tvgh
  !-----------------------------------------------------------------
  integer  ,intent(in)             :: n         ! Number of gridpoints
  real(wp) ,intent(out)            :: t     (n) ! temperature
  real(wp) ,intent(out)            :: q     (n) ! specific humidity
  real(wp) ,intent(in)             :: tv    (n) ! virtual temperature
  real(wp) ,intent(in)             :: gh    (n) ! generalized humidity
  real(wp) ,intent(in)             :: p     (n) ! pressure
  integer  ,intent(out)            :: i_fail(n) ! >0: number of iterations
                                                ! <0: no convergence
  real(wp) ,intent (out)           :: dt_tv (n) ! part. deriv. dt/dtv
  real(wp) ,intent (out)           :: dt_gh (n) ! part. deriv. dt/dgh
  real(wp) ,intent (out)           :: dq_tv (n) ! part. deriv. dq/dtv
  real(wp) ,intent (out)           :: dq_gh (n) ! part. deriv. dq/dgh
  real(wp) ,intent (out) ,optional :: drh_gh(n) !-------------
    ! local arrays
    !-------------
    real(wp) :: rh(n)    ! relative humidity
    real(wp) :: drhgh(n) ! d rh / d gh
    !-------------------
    ! transform gh -> rh
    !-------------------
    call rh_gh (rh, gh, drhgh)
    !-----------------------
    ! transform tv,rh -> t,q
    !-----------------------
#ifdef __NEC__
    call tq_tvrh_vec2 &
#else
    call tq_tvrh_vec  &
#endif
                     (n, t, q, tv, rh, p, i_fail=i_fail, &
                      dt_tv=dt_tv, dt_rh=dt_gh, dq_tv=dq_tv, dq_rh=dq_gh)
    !-----------------------------
    ! adapt Jacobi-matrix elements
    !-----------------------------
!   if (present (dt_gh)) dt_gh = dt_gh * drhgh
!   if (present (dq_gh)) dq_gh = dq_gh * drhgh
    dt_gh = dt_gh * drhgh
    dq_gh = dq_gh * drhgh
    if (present(drh_gh)) drh_gh = drhgh

  end subroutine tq_tvgh_adj_vec
!------------------------------------------------------------------------------
  subroutine read_cntrlvar_nml
  !---------------------
  ! read namelist /CNTRLVAR/
  !---------------------
    logical :: mask(mlay)
    integer :: ierr
    !-------------
    ! set defaults
    !-------------
    rh0        = 0.03_wp    ! lower bound for linear range (rel.hum)
    rh0fg      = r_init
    rh0fg_ptop = 0._wp      ! lower bound for first guess rel.humidity
    q0fg       = 1.25e-6_wp ! lower bound for first guess spc.humidity
    rhmax      = r_init     ! upper bound for relative humidity
    rh1        = r_init     ! upper bound for linear range (rel.hum)
    rh1fg      = r_init     ! upper bound for first guess rel.humidity
    rh1fg_ptop = 0._wp      ! upper bound for first guess rel.humidity
    version    = 1          ! 1:exponential; 2:quadratic
    nl_prec    = .false.    ! account for nonlinearity in preconditioner
    !--------------
    ! read namelist
    !--------------
    if (dace% lpio) then
      write(6,'(a)') repeat('-',79)
      write(6,'( )')
      call position_nml ('CNTRLVAR', status=ierr)
      select case (ierr)
      case (POSITIONED)
        write(6,'(a)') '  reading namelist /CNTRLVAR/'
#if defined(__ibm__)
        read (nnml ,nml=CNTRLVAR, iostat=ierr)
        if (ierr/=0) call finish ('read_cntrlvar_nml',&
                                  'ERROR in namelist /CNTRLVAR/')
#else
        read (nnml ,nml=CNTRLVAR)
#endif
      case default
        write(6,'(a)') '  namelist /CNTRLVAR/ is not present'
        write(6,'(a)') '  using default parameters'
      end select

      select case (version)
      case default
        call finish('read_cntrlvar_nml','version not in 0,1,2')
      case (0)
        call disable_gh ()
      case (1,2)
        mask = (rh0fg(:) /= r_init)
        nlay0 = count(mask)
        rh0fg     (1:nlay0) = pack(rh0fg     (:), mask=mask)
        rh0fg_ptop(1:nlay0) = pack(rh0fg_ptop(:), mask=mask)
        if (version == 1 .and. nlay0 == 0) then
          ! Set default for backwards compatibility
          rh0fg     (1) = 0.01_wp
          rh0fg_ptop(1) = 0._wp
          nlay0         = 1
        end if
        where(rh0fg_ptop(1:nlay0) == r_init) rh0fg_ptop(1:nlay0) = 0._wp
        ! Avoid that levels are above these values due to rounding errors:
        rh0fg_ptop(1:nlay0) = rh0fg_ptop(1:nlay0) + 10 * spacing(rh0fg_ptop(1:nlay0))
        call sort(rh0fg_ptop(1:nlay0), rh0fg(1:nlay0))
        mask = (rh1fg(:) /= r_init)
        nlay1 = count(mask)
        rh1fg     (1:nlay1) = pack(rh1fg     (:), mask=mask)
        rh1fg_ptop(1:nlay1) = pack(rh1fg_ptop(:), mask=mask)
        where(rh1fg_ptop(1:nlay1) == r_init) rh1fg_ptop(1:nlay1) = 0._wp
        ! Avoid that levels are above these values due to rounding errors:
        rh1fg_ptop(1:nlay1) = rh1fg_ptop(1:nlay1) + 10 * spacing(rh1fg_ptop(1:nlay1))
        call sort(rh1fg_ptop(1:nlay0), rh1fg(1:nlay0))
      end select

      write(6,'( )')
      write(6,'(a,f10.2)')       '  rh0        = ',rh0
      write(6,'(a,3(es11.2))')   '  rh0fg      =' ,rh0fg     (1:nlay0)
      write(6,'(a,3(1x,f10.2))') '  rh0fg_ptop =' ,rh0fg_ptop(1:nlay0)
      write(6,'(a,es10.2)')      '  q0fg       = ',q0fg
      write(6,'(a,f10.2)')       '  rhmax      = ',rhmax
      write(6,'(a,f10.2)')       '  rh1        = ',rh1
      write(6,'(a,3(es11.2))')   '  rh1fg      =' ,rh1fg     (1:max(nlay1,1))
      write(6,'(a,3(1x,f10.2))') '  rh1fg_ptop =' ,rh1fg_ptop(1:nlay1)
      write(6,'(a,i7)')          '  version    =' ,version
      write(6,'(a,l1)')          '  nl_prec    = ',nl_prec
      write(6,'( )')
    endif
    !-----------------------------
    ! broadcast namelist variables
    !-----------------------------
    call p_bcast (rh0        ,dace% pio)
    call p_bcast (rh0fg      ,dace% pio)
    call p_bcast (rh0fg_ptop ,dace% pio)
    call p_bcast (nlay0      ,dace% pio)
    call p_bcast (q0fg       ,dace% pio)
    call p_bcast (rhmax      ,dace% pio)
    call p_bcast (rh1        ,dace% pio)
    call p_bcast (rh1fg      ,dace% pio)
    call p_bcast (rh1fg_ptop ,dace% pio)
    call p_bcast (nlay1      ,dace% pio)
    call p_bcast (version    ,dace% pio)
    call p_bcast (nl_prec    ,dace% pio)
    !------------------------
    ! derive other quantities
    !------------------------
    a0    = rh0 / e
    rh02  = rh0 * 2._wp
    if (rh0>0._wp) then
      edrh0  = 1._wp / rh0
      edrh02 = 1._wp / rh02
    endif
    if (rh1 < rhmax) then
      b1   = rhmax - rh1
      edb1 = 1._wp / b1
      a1   = b1 / e
      b2   = b1 * 2._wp
      b3   = rhmax + b1
      edb2 = 1._wp / b2
    endif

    gh_min = gh_rh(0._wp, .false.)
    if (rhmax < 2._wp) gh_max = gh_rh(rhmax, .false.)
    if (dace% lpio) then
      write(6,'(a,es10.2)') '  gh_min     = ',gh_min
      write(6,'(a,es10.2)') '  gh_max     = ',gh_max
      write(6,'()')
    end if

    !---------------------------------------------------------------
    ! store generalised humidity parameters in derived type variable
    !---------------------------------------------------------------
    gh_meta% min  = 0._wp
    gh_meta% max  = rhmax
    gh_meta% lbnd = rh0
    gh_meta% ubnd = rhmax - rh1
    gh_meta% mode = version

  end subroutine read_cntrlvar_nml

  subroutine sort(r1,r2)
    real(wp), intent(inout) :: r1(:)
    real(wp), intent(inout) :: r2(:)
    real(wp) :: r_
    integer  :: i, j
    do i = size(r1)-1, 1, -1
      do j = 1, i
        if (r1(j) < r1(j+1)) then
          r_      = r1(j)
          r1(j)   = r1(j+1)
          r1(j+1) = r_
          r_      = r2(j)
          r2(j)   = r2(j+1)
          r2(j+1) = r_
        end if
      end do
    end do
  end subroutine sort

!==============================================================================

  subroutine disable_gh ()
    !----------------------------------------------------
    ! Initialize module variables so that the generalized
    ! humidity transformation effectively gets disabled
    !----------------------------------------------------
    version    = 0
    nlay0      = 0
    nlay1      = 0
    rh0        = 0._wp        ! lower bound for linear range (rel.hum)
    rh0fg      = 0._wp
    a0         = 0._wp        ! rh0 / e
    edrh0      = HUGE (0._wp) ! 1 / rh0
    rh02       = 0._wp        ! 2 * rh0
    edrh02     = HUGE (0._wp) ! 1 / 2*rh0
    a1         = 0._wp        ! b1 / e
    b1         = 0._wp        ! rhmax - rh1
    b2         = 0._wp        ! 2 * b1
    b3         = 100._wp      ! rhmax + b1  TODO:999._wp ??
    edb1       = HUGE (0._wp) ! 1 / b1
    edb2       = HUGE (0._wp) ! 1 / b2
    rhmax      = r_init       ! upper bound for relative humidity
    rh1        = r_init       ! upper bound for linear range (rel.hum)
    rh1fg      = r_init       ! upper bound for first guess rel.humidity
    rh1fg_ptop = 0._wp        ! pressure until which rh0fg is applied
    gh_max     =  huge (0._wp)
    gh_min     = -huge (0._wp)

    gh_meta% min  = 0._wp
    gh_meta% max  = rhmax
    gh_meta% lbnd = rh0
    gh_meta% ubnd = rhmax - rh1
    gh_meta% mode = version
  end subroutine disable_gh

!==============================================================================
end module mo_cntrlvar
