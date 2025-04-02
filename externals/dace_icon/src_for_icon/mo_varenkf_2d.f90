!
!+ localisation operator in 2 dimensions
!
module mo_varenkf_2d
!
! Description:
!  Localisation operator in 2 dimensions. To be used for horizontal localisation
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
!  implementation of horizontal localisation operator for VarEnKF
! V1_28        2014/02/26 Andreas Rhodin
!  further implementation of VarEnKF
! V1_31        2014-08-21 Andreas Rhodin
!  generate 2d random fields with specified correlation length scale
! V1_42        2015-06-08 Andreas Rhodin
!  option: skip test on GME double points in hor.localisation
! V1_44        2015-09-30 Andreas Rhodin
!  fix ensemble disturbances (SST) for ICON
!  VARENKF: fix for large horizontal localisation radius
! V1_45        2015-12-15 Harald Anlauf
!  reorder some loops; minor optimizations
! V1_47        2016-06-06 Andreas Rhodin
!  optimise random pattern generator for ICON
! V1_50        2017-01-09 Harald Anlauf
!  apply_c2: improve cache use
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
  use mo_kind,        only: wp,         &! working precision kind parameter
                            sp           ! single  precision kind parameter
  use mo_mpi_dace,    only: dace,       &! MPI group info
                            p_max,      &! maximum over processors
                            p_sum        ! sum     over processors
  use mo_exception,   only: finish       ! abort in case of error
  use mo_dace_string, only: char3        ! convert integer to character string
  use mo_atm_grid,    only: t_grid,     &! grid meta data derived type
                            construct,  &! set up grid meta data
                            destruct,   &! deallocate grid meta data
                            print,      &! print grid information
                            fit_rotll    ! fit rotated lat-lon grid
  use mo_ico_grid,    only: factorize_nir! GME grid factorisation
  use mo_physics,     only: d2x,        &! conversion factors: degree -> meter
                            rearth       ! earth radius
  use mo_varenkf_1d,  only: t_C1,       &! 1D localisation operator parameters
                            construct,  &! set meta data for 1D localisation
                            destruct,   &! deallocate meta data
                            apply_c11    ! apply sqrt of 1D operator in 2 d
  use mo_wmo_tables,  only: WMO6_LATLON,&!         lat-lon grid type id
                            WMO6_ROTLL ,&! rotated lat-lon grid type id
                            DWD6_ICON  ,&!            ICON grid type id
                            DWD6_ICOSAHEDRON,&!        GME grid type id
                            DWD6_NONE    ! no grid, just a collection of columns
  use mo_letkf_util,  only: gaspari_cohn ! Gaspary & Cohn function
  implicit none

!------------------------------------------------------------------------------
  !================
  ! public entities
  !================
  private
  public :: t_icof         ! component of t_C2
  public :: t_C2           ! 2D localisation operator parameters
  public :: construct      ! set meta data for 2D localisation operator
  public :: destruct       ! deallocate 2D localisation operator meta data
  public :: apply_c2       ! apply sqrt of 2D localisation operator
  public :: apply_c2t      ! apply adjoint sqrt of 2D localisation operator
  public :: apply_c2_c2t   ! apply 2D localisation operator

!------------------------------------------------------------------------------
  !===========
  ! Interfaces
  !===========
  interface construct
    module procedure construct_c2
  end interface construct

  interface destruct
    module procedure destruct_c2
  end interface destruct

!------------------------------------------------------------------------------
  !=========================
  ! Derived type definitions
  !=========================

  !---------------------------
  ! interpolation coefficients
  !---------------------------
  type t_icof
    integer  :: i                           ! coarse grid x-index
    integer  :: j                           ! coarse grid y-index
    integer  :: d                           ! coarse grid diamond-index
    real(sp) :: w                           ! weight (coefficient)
  end type t_icof

  !------------------------------------
  ! 2D localisation operator parameters
  !------------------------------------
  type t_C2
    !----------------------
    ! primary specification
    !----------------------
    integer                   :: mode               ! 1=explicit, 2=coarse grid
    type(t_grid) ,pointer     :: grid => NULL()     ! grid type
    real(wp)                  :: r_loc              ! length scale (m)
    !------------------------
    ! parameters for mode = 1
    !------------------------
    type (t_C1)               :: cx                 ! x-direction  localisation
    type (t_C1)               :: cy                 ! y-direction  parameters
    !------------------------
    ! parameters for mode = 2
    !------------------------
    type(t_icof) ,allocatable :: c  (:,:,:,:)       ! coefficients (i,j,n,d)
    integer      ,allocatable :: nc (:,:  ,:)       ! number of coefficients
    integer                   :: mc                 ! max. number of coeff.
    !------------------------------
    ! parameters for mode = 1 and 2
    !------------------------------
    type(t_grid) ,pointer     :: gc => NULL()       ! coarse grid info
  end type t_C2

!==============================================================================
contains
!==============================================================================

  subroutine construct_c2 (c2, grid, r_loc, mode, verbose, refgrid)
  !-------------------------------------------
  ! set meta data for 2D localisation operator
  !-------------------------------------------
  type (t_c2)   ,intent(out),  target :: c2     ! localisation parameters
  type (t_grid) ,pointer              :: grid   ! grid meta data
  real(wp)      ,intent(in)           :: r_loc  ! localisation length scale (m)
  integer       ,intent(in) ,optional :: mode   ! 1:explicit 2,3:coarse grid
  integer       ,intent(in) ,optional :: verbose! 0:quiet 1:info 2:details
  type (t_grid) ,pointer    ,optional :: refgrid! ensemble grid meta data

    integer  :: ni, ni2, nir, ierror ! GME grid factorisation
    integer  :: i,j,d                ! grid indices
    integer  :: ic,jc,dc             ! coarse grid indices
    integer  :: mc                   ! max number of coefficients / gridpoint
    integer  :: nc                   ! number of coefficients / gridpoint
    real(wp) :: w                    ! weight (coefficient)
    real(wp) :: ws                   ! sum of weights
    real(wp) :: r                    ! normalised distance
    real(wp) :: s                    ! scale factor
    real(wp) :: s_                   ! 1 / scale factor
    real(wp) :: xd(3,10)             ! diamond central coordinate
    real(wp) :: rd                   ! diamond influence radius
    real(wp) :: xg(3)                ! gridpoint coordinates
    integer  :: lb(4), ub(4)         ! temporary loop bounds (crayftn+OpenMP)
    logical  :: lgcmarr              ! temporary
    integer  :: verb                 ! verbosity
    integer  :: nx, ny               ! grid size
    real(wp) :: di, dj               ! grid spacing
    real(wp) :: dlat,  dlon          ! S pole of axis of rotation
    real(wp) :: dlatr, dlonr         ! N pole of axis of rotation
    real(wp) :: dlat1, dlon1         ! latitude  bounds of rot.grid
    real(wp) :: dlatn, dlonn         ! longitude bounds of rot.grid
    type(t_icof) ,pointer :: c       ! coefficient pointer
    logical  :: global               ! localization on global grid?
#ifdef __NEC__
    integer  :: nc1
#endif
    !------------------------
    ! set mandatory arguments
    !------------------------
    c2% grid  => grid
    c2% r_loc =  r_loc
    !---------
    ! set mode
    !---------
    select case (grid% gridtype)
    case (WMO6_LATLON, WMO6_ROTLL)
      c2% mode = 1
      if (present (mode)) c2% mode = mode
    case default
      global = grid% global
      if (present (refgrid)) then
         if (associated (refgrid)) global = global .or. refgrid% global
      end if
      if (global) then
        c2% mode = 2
      else
        c2% mode = 3
        if (present (mode)) c2% mode = mode
        if (c2% mode <= 1) &
             call finish('construct_c2','conflicting mode '//char3(c2% mode))
      end if
    end select
    verb = 1; if (present (verbose)) verb = verbose
    !------------------
    ! set up parameters
    !------------------
    select case (c2% mode)
    case default
      call finish('construct_c2','invalid mode '//char3(c2% mode))
    !------------------------------------------------
    ! set up for mode == 1 (explicit, seperable grid)
    !------------------------------------------------
    case (1)
      !-----------------
      ! for lat-lon only
      !-----------------
      select case (grid% gridtype)
      case (WMO6_LATLON, WMO6_ROTLL)
      case default
        call finish('construct_c2',                                  &
                    'invalid mode 1 for grid '//char3(grid% gridtype))
      end select
      call construct ( c2% cx, grid% nx, r_loc/(grid% di*d2x))
      call construct ( c2% cy, grid% ny, r_loc/(grid% dj*d2x))
      allocate       ( c2% gc)
      call construct ( c2% gc,                            &
            gridtype = grid% gridtype,                    &
                  nx = c2% cx% nm,                        &
                  ny = c2% cy% nm,                        &
                 lo1 = grid% lo1 - c2% cx% mo * grid% di, &
                 la1 = grid% la1 - c2% cy% mo * grid% dj, &
                  di =             c2% cx% mi * grid% di, &
                  dj =             c2% cy% mi * grid% dj  )
    !---------------------------------------------------------
    ! set up for mode == 2 (non-symmetric sqrt on coarse grid)
    !---------------------------------------------------------
    case (2)
      !-----------------------------------
      ! estimate suitable GME grid spacing
      !-----------------------------------
      ni = ceiling (7054000._wp / r_loc)
      mc = 50             ! max. number of coarse points within loc.radius
      if (ni < 4) then
        ni =   4          ! constraint in construct_atm_grid ?
        mc = 250          ! total no.points for GME grid with ni=4
      endif
      do
        call factorize_nir (ni, ni2, nir, ierror)
        if (nir == 1 .or. nir == 3) exit
        ni = ni + 1
      end do
      allocate       (c2% gc)
      call construct (c2% gc, gridtype = DWD6_ICOSAHEDRON, ni = ni)
      !----------------------------------------
      ! allocate interpolation coefficient list
      !----------------------------------------
      allocate (c2% c (c2% grid% lb(1) : c2% grid% ub(1),&
                       c2% grid% lb(2) : c2% grid% ub(2),&
                       mc,                               &
                       c2% grid% lb(4) : c2% grid% ub(4)))
      allocate (c2%nc (c2% grid% lb(1) : c2% grid% ub(1),&
                       c2% grid% lb(2) : c2% grid% ub(2),&
                       c2% grid% lb(4) : c2% grid% ub(4)))
      !-----------------------
      ! calculate coefficients
      !-----------------------
      c2% c% w = 0._sp
      c2%   nc = 0
      c2%   mc = 0
      s        = rearth / (sqrt (10._wp/3._wp) * r_loc)
      s        = s * sqrt(2._wp)
      s_       = 2._wp / s
      lgcmarr  = associated (c2% gc% marr)
      !----------------------------------------------
      ! calculate diamond mean coordinates and radius
      !----------------------------------------------
      xd = 0._wp
      do dc = c2% gc% lbg(4), c2% gc% ubg(4)
        xd (:,dc) = ( c2%gc% xnglob (c2%gc% lbg(1),c2%gc% lbg(2),:,dc) &
                    + c2%gc% xnglob (c2%gc% ubg(1),c2%gc% lbg(2),:,dc) &
                    + c2%gc% xnglob (c2%gc% lbg(1),c2%gc% ubg(2),:,dc) &
                    + c2%gc% xnglob (c2%gc% ubg(1),c2%gc% ubg(2),:,dc) ) / 4._wp
      end do
      rd = max (sqrt (sum ((c2%gc% xnglob (c2%gc% lbg(1),c2%gc% lbg(2),:,1) - xd (:,1)) ** 2)),&
                sqrt (sum ((c2%gc% xnglob (c2%gc% ubg(1),c2%gc% lbg(2),:,1) - xd (:,1)) ** 2)) )
      rd = rd + s_
      !---------------------
      ! loop over gridpoints
      !---------------------
      lb = c2% grid% lb
      ub = c2% grid% ub
!$omp parallel do private(i,j,d,ic,jc,dc,nc,ws,r,w,c,xg) &
#ifdef __NEC__
!$omp             private(nc1)                           &
#endif
!$omp             collapse(3) schedule(dynamic)
      do d = lb(4), ub(4)
      do j = lb(2), ub(2)
      do i = lb(1), ub(1)
        !-----------------------
        ! skip GME double points
        !-----------------------
        if (c2% grid% gridtype == DWD6_ICOSAHEDRON) then
          if (c2% grid% marr (1,i,j,d) /= dace% pe) cycle
          if (c2% grid% marr (4,i,j,d) /= d       ) cycle
        endif
        !----------------------------
        ! loop over coarse gridpoints
        !----------------------------
        nc = 0
        ws = 0._wp
        xg = c2% grid% xnglob (i,j,1:3,d)
        do dc = c2% gc% lbg(4), c2% gc% ubg(4)
          r = sqrt ((xg(1) - xd(1,dc)) ** 2 + &
                    (xg(2) - xd(2,dc)) ** 2 + &
                    (xg(3) - xd(3,dc)) ** 2   )
          if (r > rd) cycle
          do jc = c2% gc% lbg(2), c2% gc% ubg(2)
#ifdef __NEC__
          nc1 = -9
#endif
!NEC$ ivdep
          do ic = c2% gc% lbg(1), c2% gc% ubg(1)
            !-----------------------
            ! skip GME double points
            !-----------------------
            if (lgcmarr) then
              if (c2% gc% marr (4,ic,jc,dc) /= dc) cycle
            endif
            !-----------------
            ! calculate weight
            !-----------------
            if (abs (xg(1) - c2% gc  % xnglob (ic,jc,1,dc)) > s_ ) cycle
            if (abs (xg(2) - c2% gc  % xnglob (ic,jc,2,dc)) > s_ ) cycle
            if (abs (xg(3) - c2% gc  % xnglob (ic,jc,3,dc)) > s_ ) cycle
            r = s * sqrt ((xg(1) - c2% gc  % xnglob (ic,jc,1,dc)) ** 2 + &
                          (xg(2) - c2% gc  % xnglob (ic,jc,2,dc)) ** 2 + &
                          (xg(3) - c2% gc  % xnglob (ic,jc,3,dc)) ** 2   )
            w = gaspari_cohn (r, 1._wp)
            !------------------
            ! store coefficient
            !------------------
            if (w == 0) cycle
            nc = nc + 1
#ifdef __NEC__
            if (nc1 == -9 .and. nc > mc) nc1 = nc
#else
            if (nc > mc) call finish                                            &
                         ('construct_c2','nc > mc : '//char3(nc)//' '//char3(mc))
#endif
            c => c2% c (i,j,nc,d)
            c% i = ic
            c% j = jc
            c% d = dc
            c% w = w
            ws   = ws + w * w
          end do ! ic
#ifdef __NEC__
          if (nc1 /= -9) call finish                                             &
                         ('construct_c2','nc > mc : '//char3(nc1)//' '//char3(mc))
#endif
          end do ! jc
        end do   ! dc
        if (nc > 0) then
          !-----------------------------------
          ! store # of coefficients, normalise
          !-----------------------------------
          c2% nc (i,j,d) = nc
          ws = 1._wp / sqrt (ws)
          c2% c (i,j,1:nc,d)% w = c2% c (i,j,1:nc,d)% w * ws
          c2% mc = max (c2% mc, nc)
        endif
      end do
      end do
      end do
!$omp end parallel do
      c2% mc = p_max (c2% mc)
    !---------------------------------------------------------
    ! set up for mode == 3 (non-symmetric sqrt on coarse grid)
    !---------------------------------------------------------
    case (3)
      !---------------------------------------------
      ! estimate using suitable rotated lat-lon grid
      !---------------------------------------------
      if (grid% gridtype == DWD6_NONE) then
         if (.not. present (refgrid)) &
              call finish ("construct_c2","gridtype==DWD6_NONE requires refgrid")
         if (.not. associated (refgrid)) &
              call finish ("construct_c2","gridtype==DWD6_NONE requires refgrid")
         call fit_rotll (refgrid,  dlonr, dlatr, dlon1, dlat1, dlonn, dlatn)
      else
         call fit_rotll (c2% grid, dlonr, dlatr, dlon1, dlat1, dlonn, dlatn)
      end if
      dlat = - dlatr
      dlon =   dlonr - 180._wp; if (dlon < -180._wp) dlon = dlonr + 180._wp
      !-------------------------------
      ! grid spacing, surplus boundary
      !-------------------------------
      di    = r_loc / d2x * sqrt (10._wp/3._wp)
      di    = di / sqrt (2._wp)
      dlon1 = dlon1 - di
      dlonn = dlonn + di
      dlat1 = dlat1 - di
      dlatn = dlatn + di
      dj    = di * 0.6_wp
      nx    = ceiling ( (dlonn - dlon1) / dj) + 1
      ny    = ceiling ( (dlatn - dlat1) / dj) + 1
      di    = (dlonn - dlon1) / (nx - 1)
      dj    = (dlatn - dlat1) / (ny - 1)
      allocate       (c2% gc)
      call construct (c2% gc,                &
                      gridtype = WMO6_ROTLL, &
                      nx       = nx,         &
                      ny       = ny,         &
                      di       = di,         &
                      dj       = dj,         &
                      lor      = dlon,       &
                      lar      = dlat,       &
                      lo1      = dlon1,      &
                      la1      = dlat1       )
      mc = 50   ! max. number of coarse points within loc.radius
      !----------------------------------------
      ! allocate interpolation coefficient list
      !----------------------------------------
      allocate (c2% c (c2% grid% lb(1) : c2% grid% ub(1),&
                       c2% grid% lb(2) : c2% grid% ub(2),&
                       mc,                               &
                       c2% grid% lb(4) : c2% grid% ub(4)))
      allocate (c2%nc (c2% grid% lb(1) : c2% grid% ub(1),&
                       c2% grid% lb(2) : c2% grid% ub(2),&
                       c2% grid% lb(4) : c2% grid% ub(4)))
      !-----------------------
      ! calculate coefficients
      !-----------------------
      c2% c% w = 0._sp
      c2%   nc = 0
      c2%   mc = 0
      s        = rearth / (sqrt (10._wp/3._wp) * r_loc)
      s        = s * sqrt(2._wp)
      s_       = 2._wp / s
      !---------------------
      ! loop over gridpoints
      !---------------------
      lb = c2% grid% lb
      ub = c2% grid% ub
!$omp parallel do private(i,j,d,ic,jc,dc,nc,ws,xg,r,w,c) &
#ifdef __NEC__
!$omp             private(nc1)                           &
#endif
!$omp             collapse(3) schedule(dynamic)
      do d = lb(4), ub(4)
      do j = lb(2), ub(2)
      do i = lb(1), ub(1)
         !----------------------------
         ! loop over coarse gridpoints
         !----------------------------
         nc = 0
         ws = 0._wp
         xg = c2% grid% xnglob (i,j,1:3,d)
         do dc = c2% gc% lbg(4), c2% gc% ubg(4)
         do jc = c2% gc% lbg(2), c2% gc% ubg(2)
#ifdef __NEC__
         nc1 = -9
#endif
!NEC$ ivdep
         do ic = c2% gc% lbg(1), c2% gc% ubg(1)
            !-----------------
            ! calculate weight
            !-----------------
            r = s * sqrt ((xg(1) - c2% gc  % xnglob (ic,jc,1,dc)) ** 2 + &
                          (xg(2) - c2% gc  % xnglob (ic,jc,2,dc)) ** 2 + &
                          (xg(3) - c2% gc  % xnglob (ic,jc,3,dc)) ** 2   )
            w = gaspari_cohn (r, 1._wp)
            !------------------
            ! store coefficient
            !------------------
            if (w == 0) cycle
            nc = nc + 1
#ifdef __NEC__
            if (nc1 == -9 .and. nc > mc) nc1 = nc
#else
            if (nc > mc) call finish                                            &
                         ('construct_c2','nc > mc : '//char3(nc)//' '//char3(mc))
#endif
            c => c2% c (i,j,nc,d)
            c% i = ic
            c% j = jc
            c% d = dc
            c% w = w
            ws   = ws + w * w
         end do ! ic
#ifdef __NEC__
          if (nc1 /= -9) call finish                                             &
                         ('construct_c2','nc > mc : '//char3(nc1)//' '//char3(mc))
#endif
         end do ! jc
         end do ! dc
         if (nc > 0) then
           !-----------------------------------
           ! store # of coefficients, normalise
           !-----------------------------------
           c2% nc (i,j,d) = nc
           ws = 1._wp / sqrt (ws)
           c2% c (i,j,1:nc,d)% w = c2% c (i,j,1:nc,d)% w * ws
           c2% mc = max (c2% mc, nc)
         endif
      end do
      end do
      end do
!$omp end parallel do
      c2% mc = p_max (c2% mc)
    end select
    !---------
    ! printout
    !---------
    if (dace% lpio .and. verb > 0) then
      write(6,'()')
      write(6,'(a)') '  2-D localisation operator parameters:'
      write(6,'()')
      write(6,'(a,i4,2x,a)') '    mode   = ',c2% mode  ,         &
                                           ' ! 1=explicit, 2,3=coarse grid'
      write(6,'(a,i4,2x,a)') '    grid   = ',c2% grid% gridtype
      write(6,'(a,f6.1 ,a)') '    r_loc  = ',c2% r_loc/1000._wp ,&
                                           ' ! localisation length scale (km)'
      select case (c2% mode)
      case (1)
      case (2,3)
      write(6,'(a,i4,2x,a)') '    gc     = ',c2% gc% gridtype,   &
                                           ' ! coarse grid'
      select case (c2% mode)
      case (2)
      write(6,'(a,i4,2x,a)') '    ni     = ',c2% gc% ni,         &
                                           ' ! cg resolution'
      case (3)
      write(6,'(a,i4,2x,a)') '    nx     = ',c2% gc% nx,         &
                                           ' ! cg nx gridpoints'
      write(6,'(a,i4,2x,a)') '    ny     = ',c2% gc% ny,         &
                                           ' ! cg ny gridpoints'
      end select
      write(6,'(a,i4,2x,a)') '    mc     = ',c2% mc,             &
                                           ' ! max.# of coefficients'
      end select
      write(6,'()')
    end if
    if (verb >= 2) call print (c2% gc)

  end subroutine construct_c2

!------------------------------------------------------------------------------

  subroutine destruct_c2 (c2)
  !-----------------------------------------------
  ! deallocate 2D localisation operator parameters
  !-----------------------------------------------
  type (t_c2) ,intent(inout) :: c2

    nullify         (c2% grid)
    call destruct   (c2% cx)
    call destruct   (c2% cy)
    if (associated  (c2% gc)) then
      call destruct (c2% gc)
      deallocate    (c2% gc)
      nullify       (c2% gc)
    endif
    if (allocated   (c2% c )) then
      deallocate    (c2% c )
      deallocate    (c2% nc)
    endif
    c2% mode = 0
  end subroutine destruct_c2

!------------------------------------------------------------------------------

  subroutine apply_c2 (c2, y, z)
  !---------------------------------------
  ! apply sqrt of 2D localisation operator
  !---------------------------------------
  type (t_c2) ,intent(in)  :: c2                   ! localisation parameters
  real(wp)    ,intent(out) :: y(c2% grid% lb(1):,& !
                                c2% grid% lb(2):,& !
                                               :,& !
                                c2% grid% lb(4):)  ! result
  real(wp)    ,intent(in)  :: z(c2% gc % lbg(1):,& !
                                c2% gc % lbg(2):,& !
                                               :,& !
                                c2% gc % lbg(4):)  ! argument
#ifdef _CRAYFTN
  contiguous :: y, z
#endif

    integer  :: i,j,d     ! grid indices
    integer  :: ic        ! weight index
    integer  :: offset(2) ! index offset for y
#ifndef __NEC__
    real(wp) :: tmp(size(y,3))
#else
    integer, parameter :: blksz = 256
    real(wp) :: tmp(blksz)
    integer  :: i__(blksz), j__(blksz)
    integer  :: ill, ilu, ilc
    integer  :: loopcount, ncmax, ij, nx, ny, ij_, ib, ist, ien
!$NEC vreg(tmp)
!$NEC vreg(i__)
!$NEC vreg(j__)
#endif

    y = 0._wp
    select case (c2% mode)
    case (1)
      offset = c2% grid% lb(1:2) - c2% grid% lbg(1:2)
      call apply_c11 (c2% cx, c2% cy, offset, y(:,:,:,1), z(:,:,:,1))
    case (2,3)
      !---------------------
      ! loop over gridpoints
      !---------------------
#ifndef __NEC__
      do d = c2% grid% lb(4), c2% grid% ub(4)
      do j = c2% grid% lb(2), c2% grid% ub(2)
      do i = c2% grid% lb(1), c2% grid% ub(1)
        if (c2% nc (i,j,d) == 0) cycle
        !-----------------------
        ! loop over coefficients
        !-----------------------
        tmp = 0._wp
        do ic = 1, c2% nc (i,j,d)
          tmp(:) =     tmp(:)                &
                     + c2% c (i,j,ic,d)% w * &
                    z (c2% c (i,j,ic,d)% i,  &
                       c2% c (i,j,ic,d)% j,:,&
                       c2% c (i,j,ic,d)% d   )
        end do
        y(i,j,:,d) = tmp(:)
      end do
      end do
      end do
#else
      ill = lbound(y,3)
      ilu = ubound(y,3)
      nx = c2% grid% ub(1) - c2% grid% lb(1) + 1
      ny = c2% grid% ub(2) - c2% grid% lb(2) + 1
      loopcount = nx*ny
      do d = c2% grid% lb(4), c2% grid% ub(4)
      do ij_ = 1, loopcount, blksz
        ist = ij_
        ien = min(ij_+blksz-1, loopcount)
        ncmax = -1
!$NEC shortloop
        do ij = ist, ien
          i = mod(ij-1,nx)+1; i = c2% grid% lb(1) + i -1
          j = (ij-1)/nx+1;    j = c2% grid% lb(2) + j -1
          ncmax = max(ncmax, c2% nc(i,j,d))
          ib = ij-ist+1
          i__(ib) = i
          j__(ib) = j
        end do

        do ilc = ill, ilu
!$NEC shortloop
          tmp = 0._wp
          do ic = 1, ncmax
!$NEC shortloop
!$NEC ivdep
            do ij = ist, ien
              ib = ij-ist+1
              i  = i__(ib)
              j  = j__(ib)
              !-----------------------
              ! loop over coefficients
              !-----------------------
              if(ic <= c2% nc(i,j,d)) then
                tmp(ib) = tmp(ib)                 &
                        + c2% c (i,j,ic,d)% w *   &
                       z (c2% c (i,j,ic,d)% i,    &
                          c2% c (i,j,ic,d)% j,ilc,&
                          c2% c (i,j,ic,d)% d   )
              end if
            end do
          end do

!$NEC shortloop
          do ij = ist, ien
            ib = ij-ist+1
            i  = i__(ib)
            j  = j__(ib)
            y(i,j,ilc,d) = tmp(ib)
          end do
        end do   ! ilc
      end do     ! ij_
      end do     ! d
#endif
    end select
  end subroutine apply_c2

!------------------------------------------------------------------------------

  subroutine apply_c2t (c2, z, x)
  !--------------------------------------------------
  ! apply adjoint of sqrt of 2D localisation operator
  !--------------------------------------------------
  type (t_c2) ,intent(in)  :: c2                   ! localisation parameters
  real(wp)    ,intent(out) :: z(c2% gc % lbg(1):,& !
                                c2% gc % lbg(2):,& !
                                               :,& !
                                c2% gc % lbg(4):)  ! argument
  real(wp)    ,intent(in)  :: x(c2% grid% lb(1):,& !
                                c2% grid% lb(2):,& !
                                               :,& !
                                c2% grid% lb(4):)  ! result
#ifdef _CRAYFTN
  contiguous :: x, z
#endif

    integer  :: i,j,d    ! grid indices
    integer  :: iz,jz,dz ! coarse grid indices
    integer  :: ic       ! weight index

    z = 0._wp
    select case (c2% mode)
    case (2,3)
      !---------------------
      ! loop over gridpoints
      !---------------------
      do d = c2% grid% lb(4), c2% grid% ub(4)
      do j = c2% grid% lb(2), c2% grid% ub(2)
      do i = c2% grid% lb(1), c2% grid% ub(1)
        if (c2% nc (i,j,d) == 0) cycle
        !-----------------------
        ! loop over coefficients
        !-----------------------
        do ic = 1, c2% nc (i,j,d)
          iz = c2% c (i,j,ic,d)% i
          jz = c2% c (i,j,ic,d)% j
          dz = c2% c (i,j,ic,d)% d
          z(iz,jz,:,dz) = z(iz,jz,:,dz) + c2% c (i,j,ic,d)% w * x(i,j,:,d)
        end do
      end do
      end do
      end do
    end select
    z = p_sum (z)
  end subroutine apply_c2t

!------------------------------------------------------------------------------

  subroutine apply_c2_c2t (c2, y, x)
  !-----------------------------------
  ! apply the 2D localisation operator
  !-----------------------------------
  type (t_c2) ,intent(in)  :: c2          ! localisation parameters
  real(wp)    ,intent(out) :: y(:,:,:,:)  ! result
  real(wp)    ,intent(in)  :: x(:,:,:,:)  ! argument

    real(wp), allocatable :: z (:,:,:,:)

    allocate (z (c2% gc% lbg(1) : c2% gc% ubg(1), &
                 c2% gc% lbg(2) : c2% gc% ubg(2), &
                 size(x, 3)                     , &
                 c2% gc% lbg(4) : c2% gc% ubg(4)) )

    call apply_c2t (c2, z, x)
    call apply_c2  (c2, y, z)
  end subroutine apply_c2_c2t

!==============================================================================
end module mo_varenkf_2d
