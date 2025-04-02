!
!+ Mathematical operations on gridded fields
!
MODULE mo_math_oper
!
! Description:
!   Routines for mathematical operations on gridded fields.
!   - Divergence of (horizontal) vector field
!   - Curl of (horizontal) vector field
!   - Gradient field (horizontal) of a scalar potential
!   - Vector field (horizontal) from a stream function
!
! Current Code Owner: DWD, Harald Anlauf
!    phone: +49 69 8062 4941
!    fax:   +49 69 8062 3721
!    email: harald.anlauf@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_48        2016-10-06 Harald Anlauf
!  new module: mathematical operations on gridded fields
!
! Code Description:
! Language: Fortran.
! Software Standards:
!
!------------------------------------------------------------------------------
  !-------------
  ! Modules used
  !-------------
  use mo_kind,         only: wp                 ! Working precision
  use mo_exception,    only: finish             ! abort on error condition
  !----------------
  ! Parallelization
  !----------------
  use mo_mpi_dace,     only: dace,             &! MPI communication info
                             p_alltoall,       &! Generic alltoall(v)
                             p_sum,            &! generic MPI sum
                             p_isend,          &! generic non-blocking send
                             p_recv,           &! generic receive
                             p_wait,           &! MPI_WAIT
                             p_barrier          ! MPI barrier
  !--------------
  ! Grid specific
  !--------------
  USE mo_wmo_tables, only:   WMO6_LATLON,      &!
                             WMO6_ROTLL,       &!
                             WMO6_GAUSSIAN,    &!
                             DWD6_ICOSAHEDRON, &! GME  triangular
                             DWD6_ICON          ! ICON triangular
  use mo_atm_grid,     only: t_grid,           &! Type definition for grid information
                             construct,        &! initialize t_grid
                             destruct,         &! free       t_grid
                             print              ! print      t_grid
  use mo_t_icon,       only: t_patch,          &! Derived type for ICON patch
                             t_grid_cells,     &! Derived type for grid cells
                             t_grid_edges,     &! Derived type for grid edges
                             t_grid_vertices,  &! Derived type for grid vertices
                             t_int_state,      &! Derived type: interpolation coeffs.
                             t_geographical_coordinates
! use mo_icon_grid,    only: t_grid_icon        ! Derived type for ICON grid metadata
  use mo_physics,      only: d2r,              &! factor degree->radians (pi/180)
                             Rearth             ! Earth radius

  implicit none
  !----------------
  ! Public entities
  !----------------
  private
  public :: hor_div             ! Derive divergence of horizontal vector field
  public :: hor_div_t           ! Dto., transposed version
  public :: hor_curl            ! Derive curl of horizontal vector field
  public :: hor_curl_t          ! Dto., transposed version
  public :: hor_grad            ! Derive horizontal gradient of scalar potential
  public :: hor_hrot            ! Derive vector field from scalar stream fct.
  public :: math_oper_init      ! Precompute operator coefficients
  public :: math_oper_cleanup   ! Clean up module
  public :: construct_icon_intp ! Precompute interpolation coefficients
  public :: destruct_icon_intp  ! Drop precomputed interpolation coefficients
  public :: set_boundlines      ! Set nboundlines for limited area grid
  !============================================================================
  !------------
  ! Parameters:
  !------------
  integer                    :: nboundlines = 3 ! nboundlines from COSMO
  !============================================================================
  !------------------------------------
  ! Derived types and module variables:
  !------------------------------------
  type t_stencil
     integer,    allocatable :: iidx(:,:,:)     !     line indices (m,i,j)
     integer,    allocatable :: iblk(:,:,:)     !    block indices (m,i,j)
     real(wp),   allocatable :: w_u (:,:,:)     ! Stencil: weights (m,i,j)
     real(wp),   allocatable :: w_v (:,:,:)     ! Stencil: weights (m,i,j)
     integer                 :: lb(4) = 0       ! Grid lower bounds (local)
     integer                 :: ub(4) = 0       ! Grid upper bounds (local)
  end type t_stencil

  type(t_stencil), save      :: s               ! Stencil for divergence
  type(t_stencil), save      :: r               ! Stencil for curl (rot)

  !---------------------
  ! Halo, administration
  !---------------------
  type t_halo
     integer,    allocatable :: iidx  (:)       ! (nhalo)
     integer,    allocatable :: iblk  (:)       ! (nhalo)
     real(wp),   allocatable :: u   (:,:)       ! (nlev,nhalo)
     real(wp),   allocatable :: v   (:,:)       ! (nlev,nhalo)
  end type t_halo

  type t_halo_adm
     integer                 :: nhalo = 0       ! Total halo points
     integer,    allocatable :: off   (:)       ! offset (0...(nprocs_1))
     integer,    allocatable :: cnt   (:)       ! count  (0...(nprocs_1))
     integer,    allocatable :: hidx(:,:)       ! (jc,jb) -> ihalo
     type(t_halo)            :: halo            !
  end type t_halo_adm

  type(t_halo_adm), save     :: d_send          ! Halo, points to send (div)
  type(t_halo_adm), save     :: d_recv          ! Halo, points to recv (div)

  type(t_halo_adm), save     :: r_send          ! Halo, points to send (rot)
  type(t_halo_adm), save     :: r_recv          ! Halo, points to recv (rot)

  !--------------------------------------
  ! Interpolation coefficients: ICON grid
  !--------------------------------------
  type(t_int_state), allocatable, save :: icon_intp

  !-----------------------
  ! Module is initialized?
  !-----------------------
  logical, save              :: init = .false.
  !============================================================================
! integer, save              :: dbg_level = 5
  integer, save              :: dbg_level = 0
  !============================================================================
CONTAINS
  !============================================================================
  subroutine math_oper_init (grid, only)
    !------------------------------------
    ! Initialize operators for given grid
    !------------------------------------
    type(t_grid), intent(in)           :: grid
    character(*), intent(in), optional :: only

    character(4) :: lonly

    lonly = ""; if (present (only)) lonly = only

    select case (grid% gridtype)
    case (DWD6_ICON)
       !--------------------------------
       ! Reordered ICON grid not handled
       !--------------------------------
       if (grid% d_gme(1) >= 0) &
            call finish ("math_oper_init", "grid% d_gme(1) >= 0")

       call construct_icon_intp (icon_intp, grid)
       select case (lonly)
       case default
          call init_stencil_icon_div  (s, grid, icon_intp)
          call init_stencil_icon_curl (r, grid, icon_intp)
       case ("div")
          call init_stencil_icon_div  (s, grid, icon_intp)
       case ("curl")
          call init_stencil_icon_curl (r, grid, icon_intp)
       end select
    case (WMO6_LATLON, WMO6_ROTLL)
       ! Nothing to do
    case default
       write (0,*) "math_oper_init: not implemented, gridtype =", grid% gridtype
       call finish ("math_oper_init", "gridtype not implemented!")
    end select
    init = .true.
  end subroutine math_oper_init
  !============================================================================
  subroutine math_oper_cleanup ()
    !----------------------
    ! Deallocate structures
    !----------------------
    call destruct_halo    (d_send)      ! Divergence
    call destruct_halo    (d_recv)
    call destruct_stencil (s)
    call destruct_halo    (r_send)      ! Curl
    call destruct_halo    (r_recv)
    call destruct_stencil (r)
    call destruct_icon_intp (icon_intp)
    init = .false.
  end subroutine math_oper_cleanup
  !============================================================================
  subroutine set_boundlines (nbl)
    integer, intent(in) :: nbl
    if (nbl < 0) then
       nboundlines = 3 ! default: nboundlines from COSMO
    else
       nboundlines = nbl
    end if
  end subroutine set_boundlines
  !============================================================================
  subroutine hor_div (u, v, div, grid)
    real(wp), pointer, intent(in)  :: u  (:,:,:,:)  ! Eastward  vector component
    real(wp), pointer, intent(in)  :: v  (:,:,:,:)  ! Northward vector component
    real(wp), pointer, intent(in)  :: div(:,:,:,:)  ! Divergence
    type(t_grid),      intent(in)  :: grid          ! Grid meta data
    !---------------------------------------------
    ! Calculate divergence of a given vector field
    ! with horizontal components u, v
    !---------------------------------------------
    logical :: linit          ! Coefficients already set up?

    linit = init
    if (.not. linit) call math_oper_init (grid)

    select case (grid% gridtype)
    case (DWD6_ICON)
       call div_icon   (u, v, div, grid)
    case (WMO6_LATLON, WMO6_ROTLL)
       call div_latlon (u, v, div, grid)
    case (WMO6_GAUSSIAN)
       call div_gauss  (u, v, div, grid)
    case default
       write (0,*) "hor_div: not implemented, gridtype =", grid% gridtype
       call finish ("hor_div", "gridtype not implemented!")
    end select

    if (.not. linit) call math_oper_cleanup ()
  end subroutine hor_div
  !============================================================================
  subroutine hor_div_t (u_t, v_t, div, grid)
    real(wp), pointer, intent(in)  :: u_t(:,:,:,:)  ! Eastward  vector component
    real(wp), pointer, intent(in)  :: v_t(:,:,:,:)  ! Northward vector component
    real(wp), pointer, intent(in)  :: div(:,:,:,:)  ! Divergence
    type(t_grid),      intent(in)  :: grid          ! Grid meta data
    !-----------------------------
    ! Transpose version of hor_div
    !-----------------------------
    logical :: linit          ! Coefficients already set up?

    linit = init
    if (.not. linit) call math_oper_init (grid)

    select case (grid% gridtype)
    case (DWD6_ICON)
       call div_icon_t   (u_t, v_t, div, grid)
    case (WMO6_LATLON, WMO6_ROTLL)
       call div_latlon_t (u_t, v_t, div, grid)
!!$    case (WMO6_GAUSSIAN)
!!$       call div_gauss_t  (u_t, v_t, div, grid)
    case default
       write (0,*) "hor_div_t: not implemented, gridtype =", grid% gridtype
       call finish ("hor_div_t", "gridtype not implemented!")
    end select

    if (.not. linit) call math_oper_cleanup ()
  end subroutine hor_div_t
  !============================================================================
  subroutine hor_curl (u, v, rot, grid)
    real(wp), pointer, intent(in)  :: u  (:,:,:,:)  ! Eastward  vector component
    real(wp), pointer, intent(in)  :: v  (:,:,:,:)  ! Northward vector component
    real(wp), pointer, intent(in)  :: rot(:,:,:,:)  ! Curl
    type(t_grid),      intent(in)  :: grid          ! Grid meta data
    !---------------------------------------
    ! Calculate curl of a given vector field
    ! with horizontal components u, v
    !---------------------------------------
    logical :: linit          ! Coefficients already set up?

    linit = init
    if (.not. linit) call math_oper_init (grid)

    select case (grid% gridtype)
    case (DWD6_ICON)
       call curl_icon   (u, v, rot, grid)
    case (WMO6_LATLON, WMO6_ROTLL)
       call curl_latlon (u, v, rot, grid)
!!$    case (WMO6_GAUSSIAN)
!!$       call curl_gauss  (u, v, rot, grid)
    case default
       write (0,*) "hor_curl: not implemented, gridtype =", grid% gridtype
       call finish ("hor_curl", "gridtype not implemented!")
    end select

    if (.not. linit) call math_oper_cleanup ()
  end subroutine hor_curl
  !============================================================================
  subroutine hor_curl_t (u_t, v_t, rot, grid)
    real(wp), pointer, intent(in)  :: u_t(:,:,:,:)  ! Eastward  vector component
    real(wp), pointer, intent(in)  :: v_t(:,:,:,:)  ! Northward vector component
    real(wp), pointer, intent(in)  :: rot(:,:,:,:)  ! Curl
    type(t_grid),      intent(in)  :: grid          ! Grid meta data
    !------------------------------
    ! Transpose version of hor_curl
    !------------------------------
    logical :: linit          ! Coefficients already set up?

    linit = init
    if (.not. linit) call math_oper_init (grid)

    select case (grid% gridtype)
    case (DWD6_ICON)
       call curl_icon_t   (u_t, v_t, rot, grid)
    case (WMO6_LATLON, WMO6_ROTLL)
       call curl_latlon_t (u_t, v_t, rot, grid)
!!$    case (WMO6_GAUSSIAN)
!!$       call curl_gauss_t  (u_t, v_t, rot, grid)
    case default
       write (0,*) "hor_curl_t: not implemented, gridtype =", grid% gridtype
       call finish ("hor_curl_t", "gridtype not implemented!")
    end select

    if (.not. linit) call math_oper_cleanup ()
  end subroutine hor_curl_t
  !============================================================================
  subroutine hor_grad (chi, u, v, grid)
    real(wp), pointer, intent(in)  :: chi(:,:,:,:)  ! Scalar potential
    real(wp), pointer, intent(in)  :: u  (:,:,:,:)  ! Eastward  vector component
    real(wp), pointer, intent(in)  :: v  (:,:,:,:)  ! Northward vector component
    type(t_grid),      intent(in)  :: grid          ! Grid meta data
    !----------------------------------------------------------
    ! Calculate horizontal gradient of a given scalar potential
    !----------------------------------------------------------
!   logical :: linit          ! Coefficients already set up?

!   linit = init
!   if (.not. linit) call math_oper_init (grid)

    select case (grid% gridtype)
!   case (DWD6_ICON)
!      call grad_icon   (chi, u, v, grid)
    case (WMO6_LATLON, WMO6_ROTLL)
       call grad_latlon (chi, u, v, grid)
!   case (WMO6_GAUSSIAN)
!      call grad_gauss  (chi, u, v, grid)
    case default
       write (0,*) "hor_grad: not implemented, gridtype =", grid% gridtype
       call finish ("hor_grad", "gridtype not implemented!")
    end select

!   if (.not. linit) call math_oper_cleanup ()
  end subroutine hor_grad
  !============================================================================
  subroutine hor_hrot (psi, u, v, grid)
    real(wp), pointer, intent(in)  :: psi(:,:,:,:)  ! Stream function
    real(wp), pointer, intent(in)  :: u  (:,:,:,:)  ! Eastward  vector component
    real(wp), pointer, intent(in)  :: v  (:,:,:,:)  ! Northward vector component
    type(t_grid),      intent(in)  :: grid          ! Grid meta data
    !---------------------------------------------------------
    ! Calculate vector field from given scalar stream function
    !---------------------------------------------------------
!   logical :: linit          ! Coefficients already set up?

!   linit = init
!   if (.not. linit) call math_oper_init (grid)

    select case (grid% gridtype)
!   case (DWD6_ICON)
!      call hrot_icon   (psi, u, v, grid)
    case (WMO6_LATLON, WMO6_ROTLL)
       call hrot_latlon (psi, u, v, grid)
!   case (WMO6_GAUSSIAN)
!      call hrot_gauss  (psi, u, v, grid)
    case default
       write (0,*) "hor_hrot: not implemented, gridtype =", grid% gridtype
       call finish ("hor_hrot", "gridtype not implemented!")
    end select

!   if (.not. linit) call math_oper_cleanup ()
  end subroutine hor_hrot
  !============================================================================
  subroutine div_gauss (u, v, div, grid)
    ! Implementation using spectral transform
    real(wp), pointer, intent(in)  :: u  (:,:,:,:)  ! Zonal      component
    real(wp), pointer, intent(in)  :: v  (:,:,:,:)  ! Meridional component
    real(wp), pointer, intent(in)  :: div(:,:,:,:)  ! Divergence
    type(t_grid),      intent(in)  :: grid

    write (0,*) "div_gauss: not yet implemented, gridtype =", grid% gridtype
    call finish ("div_gauss", "not yet implemented!")
  end subroutine div_gauss
  !============================================================================
  subroutine div_icon (u, v, div, grid)
    real(wp), pointer, intent(in)  :: u  (:,:,:,:)  ! Zonal      component
    real(wp), pointer, intent(in)  :: v  (:,:,:,:)  ! Meridional component
    real(wp), pointer, intent(in)  :: div(:,:,:,:)  ! Divergence
    type(t_grid),      intent(in)  :: grid

    integer                    :: jc, jb        ! Cell indices
    integer                    :: i, j, k       ! Loop indices
    integer                    :: ix            ! Halo index
    integer                    :: m             ! Stencil point index

    if (.not. associated (grid% icongrid)) &
         call finish ("div_icon", "icongrid not set up")
    if (.not. grid% global) &
         call finish ("div_icon", "non-global grid not yet implemented")
    if (.not. allocated (icon_intp)) &
         call finish ("div_icon", "icon_intp not set up")
    if (size (u,dim=3) /= size (div,dim=3) .or. &
        size (v,dim=3) /= size (div,dim=3)      ) then
       write(0,*) dace% pe, "div_icon: incompatible dimension 3: u,v,div=", &
            size (u,dim=3), size (v,dim=3), size (div,dim=3)
       call finish ("div_icon", "incompatible dimension 3")
    end if
    if (any (shape (u) /= shape (div)) .or. &
        any (shape (v) /= shape (div))      ) then
       write(0,*) dace% pe, "div_icon: incompatible spapes:"
       write(0,*) dace% pe, "shape (u)   =", shape (u)
       write(0,*) dace% pe, "shape (v)   =", shape (v)
       write(0,*) dace% pe, "shape (div) =", shape (div)
       call finish ("div_icon", "incompatible shapes")
    end if

    call halo_exchange_init (d_send, d_recv, u, v)
    !--------------
    ! Apply stencil
    !--------------
    div(:,:,:,:) = 0._wp
    do       j = lbound (div,2), ubound (div,2)
       do    i = lbound (div,1), ubound (div,1)
          do k = lbound (div,3), ubound (div,3)
             do m = 1, 4
                jc = s% iidx(m,i,j)
                jb = s% iblk(m,i,j)
                ix = d_recv% hidx(jc,jb)
                if (m == 1 .or. ix == 0) then
                   ! Non-halo point
                   div(i,j,k,1) = div(i,j,k,1)                 &
                                + s% w_u(m,i,j) * u(jc,jb,k,1) &
                                + s% w_v(m,i,j) * v(jc,jb,k,1)
                else
                   ! Halo point
                   div(i,j,k,1) = div(i,j,k,1)                          &
                                + s% w_u(m,i,j) * d_recv% halo% u(k,ix) &
                                + s% w_v(m,i,j) * d_recv% halo% v(k,ix)
                end if
             end do
          end do
       end do
    end do

    call halo_exchange_finish (d_send, d_recv, u, v)
  end subroutine div_icon
  !============================================================================
  subroutine halo_exchange_init (h_send, h_recv, u, v)
    !--------------------
    ! Scatter halo points
    !--------------------
    type(t_halo_adm),  intent(inout) :: h_send          ! Points to send
    type(t_halo_adm),  intent(inout) :: h_recv          ! Points to recv
    real(wp), pointer, intent(in)    :: u(:,:,:,:)      ! u component
    real(wp), pointer, intent(in)    :: v(:,:,:,:)      ! v component

    integer                :: k             ! Loop indices
    integer                :: jc, jb        ! Cell indices
    integer                :: pe            ! Processor index
    integer                :: nz            ! Number of levels
    integer                :: off           ! Offset

    real(wp),  allocatable :: sendbuf (:)   ! Buffers
    real(wp),  allocatable :: recvbuf (:)   ! for
    real(wp),  allocatable :: sendbuf2(:)   ! MPI_Alltoall
    real(wp),  allocatable :: recvbuf2(:)   !
    integer,   allocatable :: scounts (:)   ! send    counts
    integer,   allocatable :: rcounts (:)   ! receive counts
    integer                :: ns, nr        ! send, recv buffer sizes
    integer                :: is, ir        ! send, recv counter

    nz = size (u, dim=3)
    allocate (h_recv% halo% u(nz,h_recv% nhalo))
    allocate (h_recv% halo% v(nz,h_recv% nhalo))

    allocate (scounts(0:dace% npe-1))
    allocate (rcounts(0:dace% npe-1))
    scounts = h_send% cnt(0:dace% npe-1) * nz
    rcounts = h_recv% cnt(0:dace% npe-1) * nz
    !-----------------
    ! allocate buffers
    !-----------------
    ns = h_send% nhalo * nz
    nr = h_recv% nhalo * nz

    allocate (sendbuf (ns))
    allocate (recvbuf (nr))
    allocate (sendbuf2(ns))
    allocate (recvbuf2(nr))

    !-----------------
    ! Fill send buffer
    !-----------------
    k = 0
    do pe = 0, dace% npe-1
       off = h_send% off(pe)
       do is = 1, h_send% cnt(pe)
          jc = h_send% halo% iidx(off+is)
          jb = h_send% halo% iblk(off+is)
          sendbuf (k+1:k+nz) = u(jc,jb,:,1)
          sendbuf2(k+1:k+nz) = v(jc,jb,:,1)
          k = k + nz
       end do
    end do
    !--------------
    ! MPI alltoallv
    !--------------
    call p_alltoall (sendbuf,  recvbuf,  sendcounts=scounts, recvcounts=rcounts)
    call p_alltoall (sendbuf2, recvbuf2, sendcounts=scounts, recvcounts=rcounts)
    !-------------------------------
    ! fetch data from receive buffer
    !-------------------------------
    k = 0
    do pe = 0, dace% npe-1
       off = h_recv% off(pe)
       do ir = 1, h_recv% cnt(pe)
          jc = h_recv% halo% iidx(off+ir)
          jb = h_recv% halo% iblk(off+ir)
          h_recv% hidx   (jc,jb)    = off+ir
          h_recv% halo% u(:,off+ir) = recvbuf (k+1:k+nz)
          h_recv% halo% v(:,off+ir) = recvbuf2(k+1:k+nz)
          k = k + nz
       end do
    end do

    deallocate (sendbuf, sendbuf2, recvbuf, recvbuf2)

  end subroutine halo_exchange_init
  !----------------------------------------------------------------------------
  subroutine halo_exchange_finish (h_send, h_recv, u, v)
    type(t_halo_adm),  intent(inout) :: h_send          ! Points to send
    type(t_halo_adm),  intent(inout) :: h_recv          ! Points to recv
    real(wp), pointer, intent(in)    :: u(:,:,:,:)      ! u component
    real(wp), pointer, intent(in)    :: v(:,:,:,:)      ! v component

    deallocate (h_recv% halo% u)
    deallocate (h_recv% halo% v)
  end subroutine halo_exchange_finish
  !============================================================================
  subroutine halo_exchange_init_t (h_send, h_recv, u, v)
    !-----------------------------------------
    ! Scatter halo points (transposed version)
    !-----------------------------------------
    type(t_halo_adm),  intent(inout) :: h_send          ! Points to send
    type(t_halo_adm),  intent(inout) :: h_recv          ! Points to recv
    real(wp), pointer, intent(in)    :: u(:,:,:,:)      ! u component
    real(wp), pointer, intent(in)    :: v(:,:,:,:)      ! v component

    integer                :: nz            ! Number of levels
    integer                :: pe            ! Processor index
    integer                :: jc, jb        ! Cell indices
    integer                :: off           ! Offset
    integer                :: ir            ! recv counter

    do pe = 0, dace% npe-1
       off = h_recv% off(pe)
       do ir = 1, h_recv% cnt(pe)
          jc = h_recv% halo% iidx(off+ir)
          jb = h_recv% halo% iblk(off+ir)
          h_recv% hidx(jc,jb) = off+ir
       end do
    end do

    nz = size (u, dim=3)
    allocate (h_recv% halo% u(nz,h_recv% nhalo))
    allocate (h_recv% halo% v(nz,h_recv% nhalo))
    h_recv% halo% u = 0._wp
    h_recv% halo% v = 0._wp

  end subroutine halo_exchange_init_t
  !----------------------------------------------------------------------------
  subroutine halo_exchange_finish_t (h_send, h_recv, u, v)
    type(t_halo_adm),  intent(inout) :: h_send          ! Points to send
    type(t_halo_adm),  intent(inout) :: h_recv          ! Points to recv
    real(wp), pointer, intent(in)    :: u(:,:,:,:)      ! u component
    real(wp), pointer, intent(in)    :: v(:,:,:,:)      ! v component

    integer                :: k             ! Loop indices
    integer                :: jc, jb        ! Cell indices
    integer                :: pe            ! Processor index
    integer                :: nz            ! Number of levels
    integer                :: off           ! Offset

    real(wp),  allocatable :: sendbuf (:)   ! Buffers
    real(wp),  allocatable :: recvbuf (:)   ! for
    real(wp),  allocatable :: sendbuf2(:)   ! MPI_Alltoall
    real(wp),  allocatable :: recvbuf2(:)   !
    integer,   allocatable :: scounts (:)   ! send    counts
    integer,   allocatable :: rcounts (:)   ! receive counts
    integer                :: ns, nr        ! send, recv buffer sizes
    integer                :: is, ir        ! send, recv counter

    nz = size (u, dim=3)
    !-------------------
    ! Gather halo points
    ! (transpose of scatter)
    !-------------------
    ns = h_send% nhalo * nz
    nr = h_recv% nhalo * nz

    allocate (scounts(0:dace% npe-1))
    allocate (rcounts(0:dace% npe-1))
    scounts = h_send% cnt(0:dace% npe-1) * nz
    rcounts = h_recv% cnt(0:dace% npe-1) * nz
    !-----------------
    ! allocate buffers
    !-----------------
    allocate (sendbuf (ns))
    allocate (recvbuf (nr))
    allocate (sendbuf2(ns))
    allocate (recvbuf2(nr))
    !-------------------------------
    ! place data in "receive" buffer
    !-------------------------------
    k = 0
    do pe = 0, dace% npe-1
       off = h_recv% off(pe)
       do ir = 1, h_recv% cnt(pe)
          recvbuf (k+1:k+nz) = h_recv% halo% u(:,off+ir)
          recvbuf2(k+1:k+nz) = h_recv% halo% v(:,off+ir)
          k = k + nz
       end do
    end do
    !--------------
    ! MPI alltoallv
    !--------------
    call p_alltoall (recvbuf,  sendbuf,  sendcounts=rcounts, recvcounts=scounts)
    call p_alltoall (recvbuf2, sendbuf2, sendcounts=rcounts, recvcounts=scounts)
    !--------------------------
    ! Update from "send" buffer
    !--------------------------
    k = 0
    do pe = 0, dace% npe-1
       off = h_send% off(pe)
       do is = 1, h_send% cnt(pe)
          jc = h_send% halo% iidx(off+is)
          jb = h_send% halo% iblk(off+is)
          u(jc,jb,:,1) = u(jc,jb,:,1) + sendbuf (k+1:k+nz)
          v(jc,jb,:,1) = v(jc,jb,:,1) + sendbuf2(k+1:k+nz)
          k = k + nz
       end do
    end do

    deallocate (sendbuf, sendbuf2, recvbuf, recvbuf2)
    deallocate (h_recv% halo% u)
    deallocate (h_recv% halo% v)

  end subroutine halo_exchange_finish_t
  !============================================================================
  ! Transpose of div_icon
  !-------------------------
  subroutine div_icon_t (u, v, div, grid)
    real(wp), pointer, intent(in)  :: u  (:,:,:,:)  ! Zonal      component
    real(wp), pointer, intent(in)  :: v  (:,:,:,:)  ! Meridional component
    real(wp), pointer, intent(in)  :: div(:,:,:,:)  ! Divergence
    type(t_grid),      intent(in)  :: grid

    integer                        :: jc, jb        ! Cell indices
    integer                        :: i, j, k       ! Loop indices
    integer                        :: m             ! Stencil point index
    integer                        :: ix            ! Halo index

    if (.not. associated (grid% icongrid)) &
         call finish ("div_icon_t", "icongrid not set up")
    if (.not. grid% global) &
         call finish ("div_icon_t", "non-global grid not yet implemented")
    if (.not. allocated (icon_intp)) &
         call finish ("div_icon_t", "icon_intp not set up")
    if (size (u,dim=3) /= size (div,dim=3) .or. &
        size (v,dim=3) /= size (div,dim=3)      ) then
       write(0,*) dace% pe, "div_icon_t: incompatible dimension 3: u,v,div=", &
            size (u,dim=3), size (v,dim=3), size (div,dim=3)
       call finish ("div_icon_t", "incompatible dimension 3")
    end if
    if (any (shape (u) /= shape (div)) .or. &
        any (shape (v) /= shape (div))      ) then
       write(0,*) dace% pe, "div_icon_t: incompatible spapes:"
       write(0,*) dace% pe, "shape (u)   =", shape (u)
       write(0,*) dace% pe, "shape (v)   =", shape (v)
       write(0,*) dace% pe, "shape (div) =", shape (div)
       call finish ("div_icon_t", "incompatible shapes")
    end if

    call halo_exchange_init_t (d_send, d_recv, u, v)
    !-------------------------
    ! Apply transposed stencil
    !-------------------------
!   u(:,:,:,:) = 0._wp
!   v(:,:,:,:) = 0._wp
    do       j = lbound (div,2), ubound (div,2)
       do    i = lbound (div,1), ubound (div,1)
          do k = lbound (div,3), ubound (div,3)
             do m = 1, 4
                jc = s% iidx(m,i,j)
                jb = s% iblk(m,i,j)
                ix = d_recv% hidx(jc,jb)
                if (m == 1 .or. ix == 0) then
                   ! Update non-halo point
                   u(jc,jb,k,1) = u(jc,jb,k,1) + s% w_u(m,i,j) * div(i,j,k,1)
                   v(jc,jb,k,1) = v(jc,jb,k,1) + s% w_v(m,i,j) * div(i,j,k,1)
                else
                   ! Update halo point
                   d_recv% halo% u(k,ix) = d_recv% halo% u(k,ix)      &
                                         + s% w_u(m,i,j) * div(i,j,k,1)
                   d_recv% halo% v(k,ix) = d_recv% halo% v(k,ix)      &
                                         + s% w_v(m,i,j) * div(i,j,k,1)
                end if
             end do
          end do
       end do
    end do

    call halo_exchange_finish_t (d_send, d_recv, u, v)

  end subroutine div_icon_t
  !============================================================================
  subroutine init_stencil_icon_div (s, grid, intp)
    type(t_stencil),   intent(inout) :: s        ! Stencil for divergence
    type(t_grid),      intent(in)    :: grid
    type(t_int_state), intent(in)    :: intp

    type(t_patch),     pointer :: patch          ! global or top-level patch
    type(t_grid_edges),pointer :: edges
    integer                    :: ie, ib         ! Edge indices
    integer                    :: jc, jb, c      ! Cell indices
    integer                    :: i, j, e        ! Loop indices
    integer                    :: is, ir         ! Send, recv counter
    integer                    :: off
    integer                    :: pe, pe2        ! Processor indices
    integer                    :: m              ! Stencil point index
    integer                    :: lb (4), ub (4) ! Loop bounds (local)
    integer                    :: lbg(4), ubg(4) ! Loop bounds (global)
    integer,           pointer :: edge_idx(:,:,:)
    integer,           pointer :: edge_blk(:,:,:)
    integer,           pointer :: cell_idx(:,:,:)
    integer,           pointer :: cell_blk(:,:,:)
    real(wp)                   :: utmp, vtmp

    logical                    :: mine
    integer,       allocatable :: nsend    (:)   ! Estimate for send point count
    integer,       allocatable :: nrecv    (:)   ! Estimate for recv point count

    !-----------------------------------------------------------------
    ! Set up stencil:
    ! - 1st order stencil needs 4 points
    ! Topological information required:
    ! - cell centers to edges: patch% cells% edge_{idx,blk}(jl,jb,1:3)
    ! - edges to cell centers: patch% edges% cell_{idx,blk}(ie,ib,1:2)
    !-----------------------------------------------------------------
    ! Basics of computation:
    !
    ! 1) Derive normal components of vector field on edges:
    !        vn = c_lin_e(,1,) * (u(1)*n1(1) + v(1)*n2(1))
    !           + c_lin_e(,2,) * (u(2)*n1(2) + v(2)*n2(2))
    !    with (n1,n2) being the cartesian components of the edge normals,
    !        n1(1) = edges%primal_normal_cell(,1)%v1,
    !        n2(1) = edges%primal_normal_cell(,1)%v2,
    !    etc., and interpolation coefficients icon_intp%c_lin_e:
    !        c_lin_e(,1,) = edges%edge_cell_length(,2) / edges%dual_edge_length,
    !        c_lin_e(,2,) = 1 - c_lin_e(,1,)
    !
    ! 2) Apply divergence operator on normal components (theorem of Gauss):
    !        div = sum_{i=1..3} (vn(i) * w(i)),
    !    with w(i) = edges%primal_edge_length(i) * cells%edge_orientation(i)
    !              / cells%area
    !-----------------------------------------------------------------
    if (allocated (s% iidx)) then
       if (any (grid% lb /= s% lb) .or. any (grid% ub /= s% lb)) then
          write(0,*) "pe=", dace% pe, "new lb=", grid% lb, "ub=", grid% ub
          write(0,*) "pe=", dace% pe, "old lb=",    s% lb, "ub=",    s% ub
          call finish("init_stencil_icon_div","### Shape mismatch")
       end if
       return
    end if

    patch    => grid% icongrid% patch
    edge_idx => patch% cells% edge_idx
    edge_blk => patch% cells% edge_blk
    edges    => patch% edges
    cell_idx => patch% edges% cell_idx
    cell_blk => patch% edges% cell_blk

    lbg = grid% lbg
    ubg = grid% ubg
    lb  = grid% lb
    ub  = grid% ub

    allocate (s% iidx(4, lbg(1):ubg(1),lbg(2):ubg(2)))
    allocate (s% iblk(4, lbg(1):ubg(1),lbg(2):ubg(2)))
    allocate (s% w_u (4, lb (1):ub (1),lb (2):ub (2)))
    allocate (s% w_v (4, lb (1):ub (1),lb (2):ub (2)))
    s% iidx = 0
    s% iblk = 0
    s% w_u  = 0._wp
    s% w_v  = 0._wp

    allocate (nsend    (0:dace% npe-1))
    allocate (nrecv    (0:dace% npe-1))
    nsend     = 0
    nrecv     = 0
    !----------------
    ! Loop over cells
    !----------------
    do j = lbg(2), ubg(2)
    do i = lbg(1), ubg(1)
       pe             = grid% marr(1,i,j,1)
       mine           = (dace% pe == pe)
       s% iidx(1,i,j) = i
       s% iblk(1,i,j) = j
       m              = 1
       !----------------
       ! Loop over edges
       !----------------
       do e = 1, 3
          ie = edge_idx(i,j,e)
          ib = edge_blk(i,j,e)
          !----------------------------
          ! Loop over neighboring cells
          !----------------------------
          do c = 1, 2
             jc = cell_idx(ie,ib,c)
             jb = cell_blk(ie,ib,c)
             !-----------------------------------------------------------------
             ! Potential optimization: evaluate these only on owner pe of (i,j)
             !-----------------------------------------------------------------
             if (mine) then
                utmp = ( intp% geofac_div  (i, e,j ) * &
                         intp% c_lin_e     (ie,c,ib) ) &
                     * edges% primal_normal_cell(ie,ib,c)% v1
                vtmp = ( intp% geofac_div  (i, e,j ) * &
                         intp% c_lin_e     (ie,c,ib) ) &
                     * edges% primal_normal_cell(ie,ib,c)% v2
             else
                utmp = 0._wp
                vtmp = 0._wp
             end if
             if (jc == i .and. jb == j) then
                !-------------------------------------
                ! Originating cell same as destination
                !-------------------------------------
                if (mine) then
                   s% w_u (1,i,j) = s% w_u(1,i,j) + utmp
                   s% w_v (1,i,j) = s% w_v(1,i,j) + vtmp
                end if
             else
                !-----------------
                ! Neighboring cell
                !-----------------
                m              = m + 1
                s% iidx(m,i,j) = jc
                s% iblk(m,i,j) = jb
                if (mine) then
                   s% w_u (m,i,j) = utmp
                   s% w_v (m,i,j) = vtmp
                end if
                !------------------------------------------------
                ! Check if stencil crosses the processor boundary
                !------------------------------------------------
                pe2 = grid% marr(1,jc,jb,1)
                if      (      mine .and. pe2 /= dace% pe) then
                   nrecv(pe2) = nrecv(pe2) + 1          ! Request from pe2
                else if (.not. mine .and. pe2 == dace% pe) then
                   nsend(pe)  = nsend(pe)  + 1          ! Expected by pe
                end if
             end if
          end do
       end do
    end do
    end do

    if (dbg_level > 1) then
       do pe = 0, dace% npe-1
          call p_barrier()
          if (dace% pe == pe) print *, pe, "nsend    =", nsend
          if (dace% pe == pe) print *, pe, "nrecv    =", nrecv
       end do
       call p_barrier()
    end if
    if (nsend(dace% pe) /= 0 .or. nrecv(dace% pe) /= 0) then
       write(0,*) "### Error: nsend(self),nrecv(self) =",nsend(dace% pe),nrecv(dace% pe)
       call finish("init_stencil_icon_div","### Internal Error")
    end if

    !--------------------------------------------------------
    ! Determine grid points for halo exchange
    ! 1) Points to send
    ! 2) Points to receive
    ! The current implementation does not care for duplicates
    !--------------------------------------------------------
    d_send% nhalo = sum (nsend)
    allocate (d_send% halo% iidx(d_send% nhalo))
    allocate (d_send% halo% iblk(d_send% nhalo))
    allocate (d_send% off       (0:dace% npe-1))
    allocate (d_send% cnt       (0:dace% npe-1))
!   allocate (d_send% hidx      (lb(1):ub(1),lb(2):ub(2)))
    d_send% halo% iidx = 0
    d_send% halo% iblk = 0
    d_send% off        = 0
    d_send% cnt        = 0
!   d_send% hidx       = 0
    !-----------------------
    ! Loop over receiver pes
    !-----------------------
    do pe = 0, dace% npe-1
       if (pe > 0) d_send% off(pe) = d_send% off(pe-1) + d_send% cnt(pe-1)
       if (nsend(pe) == 0) cycle
       off = d_send% off(pe)
       is  = 0
       !--------------------
       ! Loop over all cells
       !--------------------
       do j    = lbg(2), ubg(2)
          do i = lbg(1), ubg(1)
             !-------------------------
             ! Cell handled by other pe
             !-------------------------
             if (pe == grid% marr(1,i,j,1)) then
                !-------------------------------
                ! Loop over other stencil points
                ! Check whether needed remotely.
                !-------------------------------
                do m = 2, 4
                   jc = s% iidx(m,i,j)
                   jb = s% iblk(m,i,j)
                   if (grid% marr(1,jc,jb,1) == dace% pe) then
                      is = is + 1
                      d_send% halo% iidx(off+is) = jc
                      d_send% halo% iblk(off+is) = jb
!                     d_send% hidx      (jc,jb)  = off+is
                   end if
                end do
             end if
          end do
       end do
       d_send% cnt (pe) = is
       if (dbg_level > 1) &
            write(0,*) "p_pe -> pe, is, nsend(pe)", dace% pe, pe, is, nsend(pe)
!      if (is > nsend(pe)) then
       if (is /= nsend(pe)) then
          write(0,*) "### Fatal Error: is > nsend(pe)", is, nsend(pe)
          call finish("init_stencil_icon_div","### Fatal Error: is > nsend")
       end if
    end do
    call p_barrier ()

    d_recv% nhalo = sum (nrecv)
    allocate (d_recv% halo% iidx(d_recv% nhalo))
    allocate (d_recv% halo% iblk(d_recv% nhalo))
    allocate (d_recv% off       (0:dace% npe-1))
    allocate (d_recv% cnt       (0:dace% npe-1))
    allocate (d_recv% hidx      (lbg(1):ubg(1),lbg(2):ubg(2)))
    d_recv% halo% iidx = 0
    d_recv% halo% iblk = 0
    d_recv% off        = 0
    d_recv% cnt        = 0
    d_recv% hidx       = 0
    !---------------------
    ! Loop over sender pes
    !---------------------
    do pe = 0, dace% npe-1
       if (pe > 0) d_recv% off(pe) = d_recv% off(pe-1) + d_recv% cnt(pe-1)
       if (nrecv(pe) == 0) cycle
       off = d_recv% off(pe)
       ir  = 0
       !---------------------
       ! Loop over "my" cells
       !---------------------
       do j    = lb(2), ub(2)
          do i = lb(1), ub(1)
             !------------------------
             ! Cell handled by this pe
             !------------------------
             if (dace% pe == grid% marr(1,i,j,1)) then
                !-------------------------------
                ! Loop over other stencil points
                !-------------------------------
                do m = 2, 4
                   jc = s% iidx(m,i,j)
                   jb = s% iblk(m,i,j)
                   if (grid% marr(1,jc,jb,1) == pe) then
                      ir = ir + 1
                      d_recv% halo% iidx(off+ir) = jc
                      d_recv% halo% iblk(off+ir) = jb
                      d_recv% hidx      (jc,jb)  = off+ir
                   end if
                end do
             end if
          end do
       end do
       d_recv% cnt (pe) = ir
       if (dbg_level > 1) &
            write(0,*) "p_pe <- pe, ir, nrecv(pe)", dace% pe, pe, ir, nrecv(pe)
!      if (ir > nrecv(pe)) then
       if (ir /= nrecv(pe)) then
          write(0,*) "### Fatal Error: ir > nrecv(pe)", ir, nrecv(pe)
          call finish("init_stencil_icon_div","### Fatal Error: ir > nrecv")
       end if
    end do

    ! TODO: Exchange halo point lists for cross-checks

  end subroutine init_stencil_icon_div
  !============================================================================
  subroutine construct_icon_intp (icon_intp, grid)
    !------------------------------------------------------
    ! Precompute interpolation coefficients cells <-> edges
    ! Needed for divergence and curl operator on ICON grid.
    !------------------------------------------------------
    type(t_int_state), allocatable :: icon_intp
    type(t_grid),       intent(in) :: grid

    integer                        :: jb, jc       ! Index variables (cells)
    integer                        :: je           ! Index variables (edges)
    integer                        :: jv           ! Index variables (vertices)
    integer                        :: ibe, ile     ! Edge indices
    integer                        :: nblks, nlen  ! Loop aux. variables
    integer                        :: nproma_e     ! Blocking sizes
    integer                        :: npromz       ! Size of last block
    type(t_patch),         pointer :: patch
    type(t_grid_cells),    pointer :: cells
    type(t_grid_edges),    pointer :: edges
    type(t_grid_vertices), pointer :: verts
    integer                        :: lb (4), ub (4) ! Loop bounds (cells)
    integer                        :: lbv(4), ubv(4) ! Loop bounds (vertices)

    if (allocated (icon_intp)) then
       return
    end if

    if (.not. allocated (icon_intp)) then
       allocate (icon_intp)
    end if

    patch => grid% icongrid% patch
    cells => patch% cells
    edges => patch% edges
    verts => patch% verts
    !-------------------------------------------
    ! Interpolation coefficients cells <-> edges
    !-------------------------------------------
    nproma_e = size (edges% edge_cell_length, 1)
    npromz   =       patch% npromz_e
    nblks    =       patch% nblks_e
    allocate (icon_intp% c_lin_e(nproma_e, 2, nblks))
    do jb = 1, nblks
       nlen = nproma_e; if (jb == nblks) nlen = npromz
       do je = 1, nlen
          icon_intp% c_lin_e(je,1,jb) = edges% edge_cell_length(je,jb,2) &
                                      / edges% dual_edge_length(je,jb)
          icon_intp% c_lin_e(je,2,jb) = 1._wp - icon_intp% c_lin_e(je,1,jb)
       end do
    end do
    !-------------------------------------
    ! Coefficients for divergence operator
    ! (keep coefficients for
    !-------------------------------------
    lb = grid% lb
    ub = grid% ub
    allocate (icon_intp% geofac_div(lb(1):ub(1), 3, lb(2):ub(2)))
    do    jb = lb(2),ub(2)
       do jc = lb(1),ub(1)
          do je = 1, 3
             ile = cells% edge_idx(jc,jb,je)
             ibe = cells% edge_blk(jc,jb,je)
             icon_intp% geofac_div(jc,je,jb) =         &
                    edges% primal_edge_length(ile,ibe) &
                  * cells% edge_orientation(jc,jb,je)  &
                  / cells% area            (jc,jb)
          end do
       end do
    end do
    !-------------------------------
    ! Coefficients for curl operator
    !-------------------------------
    lbv = grid% lb_d
    ubv = grid% ub_d
    allocate (icon_intp% geofac_rot(lbv(1):ubv(1), 6, lbv(2):ubv(2)))
    do    jb = lbv(2),ubv(2)
       do jv = lbv(1),ubv(1)
          do je = 1, 6
             if (je <= verts% num_edges(jv,jb)) then
                ile = verts% edge_idx(jv,jb,je)
                ibe = verts% edge_blk(jv,jb,je)
                icon_intp% geofac_rot(jv,je,jb) =        &
                       edges% dual_edge_length(ile,ibe)  &
                     * verts% edge_orientation(jv,jb,je) &
                     / verts% dual_area       (jv,jb)
             else
                ! Pentagon
                icon_intp% geofac_rot(jv,je,jb) = 0._wp
             end if
          end do
       end do
    end do
  end subroutine construct_icon_intp
  !============================================================================
  subroutine destruct_icon_intp (icon_intp)
    type(t_int_state), allocatable, intent(inout) :: icon_intp
    if (allocated (icon_intp)) deallocate (icon_intp)
  end subroutine destruct_icon_intp
  !============================================================================
  ! Stencil for discrete curl operator
  !-----------------------------------
  subroutine init_stencil_icon_curl (r, grid, intp)
    type(t_stencil),   intent(inout) :: r        ! Stencil for curl (rot)
    type(t_grid),      intent(in)    :: grid
    type(t_int_state), intent(in)    :: intp

    type(t_patch),     pointer :: patch          ! global or top-level patch
    type(t_grid_edges),pointer :: edges
    type(t_grid_vertices),pointer :: verts
    integer                    :: ie, ib         ! Edge indices
    integer                    :: jc, jb, c      ! Cell indices
    integer                    :: i, j, e        ! Loop indices
    integer                    :: is, ir         ! Send, recv counter
    integer                    :: off
    integer                    :: pe, pe2        ! Processor indices
    integer                    :: ne             ! Number of edges
    integer                    :: m              ! Stencil point index
    integer                    :: mm             ! Stencil point index
    integer                    :: lbg(4), ubg(4) ! Loop bounds (cells, global)
    integer                    :: lbv(4), ubv(4) ! Loop bounds (vertices)
    integer                    :: lb (4), ub (4) ! Loop bounds (local)
    integer,           pointer :: edge_idx(:,:,:)
    integer,           pointer :: edge_blk(:,:,:)
    integer,           pointer :: cell_idx(:,:,:)
    integer,           pointer :: cell_blk(:,:,:)
    real(wp)                   :: utmp, vtmp
!   integer                    :: j0, b0
    integer                    :: idx_tmp(12)
    integer                    :: blk_tmp(12)
    real(wp)                   :: w_u_tmp(12)
    real(wp)                   :: w_v_tmp(12)

    logical                    :: mine
    integer,       allocatable :: nsend    (:)   ! Estimate for send point count
    integer,       allocatable :: nrecv    (:)   ! Estimate for recv point count

    !-----------------------------------------------------------------
    ! Set up stencil:
    ! - Stencil needs 6 points (5 for hexagon cells)
    ! Topological information required:
    ! - vertices to edges:     patch% verts% edge_{idx,blk}(jl,jb,1:3)
    ! - edges to cell centers: patch% edges% cell_{idx,blk}(ie,ib,1:2)
    !-----------------------------------------------------------------
    ! Basics of computation:
    !
    ! 1) Derive normal components of vector field on edges:
    !        vn = c_lin_e(,1,) * (u(1)*n1(1) + v(1)*n2(1))
    !           + c_lin_e(,2,) * (u(2)*n1(2) + v(2)*n2(2))
    !    with (n1,n2) being the cartesian components of the edge normals,
    !        n1(1) = edges%primal_normal_cell(,1)%v1,
    !        n2(1) = edges%primal_normal_cell(,1)%v2,
    !    etc., and interpolation coefficients icon_intp%c_lin_e:
    !        c_lin_e(,1,) = edges%edge_cell_length(,2) / edges%dual_edge_length,
    !        c_lin_e(,2,) = 1 - c_lin_e(,1,)
    !
    !
    ! 2) Apply curl operator on normal components (theorem of Stokes):
    !        curl = sum_{i=1..6} (vn(i) * w(i)),
    !    with w(i) = geofac_rot(i)
    !              = edges%dual_edge_length(i) * verts%edge_orientation(i)
    !              / verts%dual_area
    !-----------------------------------------------------------------
    if (allocated (r% iidx)) then
       if (any (grid% lb /= r% lb) .or. any (grid% ub /= r% lb)) then
          write(0,*) "pe=", dace% pe, "new lb=", grid% lb, "ub=", grid% ub
          write(0,*) "pe=", dace% pe, "old lb=",    r% lb, "ub=",    r% ub
          call finish("init_stencil_icon_curl","### Shape mismatch")
       end if
       return
    end if

    patch    => grid% icongrid% patch
    verts    => patch% verts
    edge_idx => patch% verts% edge_idx
    edge_blk => patch% verts% edge_blk
    edges    => patch% edges
    cell_idx => patch% edges% cell_idx
    cell_blk => patch% edges% cell_blk

    lbv = grid% lbg_d
    ubv = grid% ubg_d
    lb  = grid% lb_d
    ub  = grid% ub_d

    allocate (r% iidx(6, lbv(1):ubv(1),lbv(2):ubv(2)))
    allocate (r% iblk(6, lbv(1):ubv(1),lbv(2):ubv(2)))
    allocate (r% w_u (6, lb (1):ub (1),lb (2):ub (2)))
    allocate (r% w_v (6, lb (1):ub (1),lb (2):ub (2)))
    r% iidx = 0
    r% iblk = 0
    r% w_u  = 0._wp
    r% w_v  = 0._wp

    allocate (nsend    (0:dace% npe-1))
    allocate (nrecv    (0:dace% npe-1))
    nsend     = 0
    nrecv     = 0
    !-------------------
    ! Loop over vertices
    !-------------------
    do j = lbv(2), ubv(2)
    do i = lbv(1), ubv(1)
       pe   = grid% pe_d(i,j,1)
       mine = (dace% pe == pe)
       idx_tmp = -1
       blk_tmp = -1
       !----------------
       ! Loop over edges
       !----------------
!      j0 = -1
!      b0 = -1
       m  = 0
       ne = verts% num_edges(i,j)
       do e = 1, ne
          ie = edge_idx(i,j,e)
          ib = edge_blk(i,j,e)
          !----------------------------
          ! Loop over neighboring cells
          !----------------------------
          do c = 1, 2
             m  = m + 1
             jc = cell_idx(ie,ib,c)
             jb = cell_blk(ie,ib,c)
             idx_tmp(m) = jc
             blk_tmp(m) = jb
             !-------------------------------------------
             ! Evaluate weights only on owner pe of (i,j)
             !-------------------------------------------
             if (mine) then
                utmp = ( intp% geofac_rot  (i, e,j ) * &
                         intp% c_lin_e     (ie,c,ib) ) &
                     * edges% primal_normal_cell(ie,ib,c)% v1
                vtmp = ( intp% geofac_rot  (i, e,j ) * &
                         intp% c_lin_e     (ie,c,ib) ) &
                     * edges% primal_normal_cell(ie,ib,c)% v2
                w_u_tmp(m) = utmp
                w_v_tmp(m) = vtmp
             else
                w_u_tmp(m) = 0._wp
                w_v_tmp(m) = 0._wp
             end if
          end do
       end do
       ! Reduce stencil (similar to insertion-sort)
       m  = 0           ! Index to update
       mm = 0           ! Last new index
       do c = 1, 2*ne
          do m = 1, mm
             if (idx_tmp(c) == r% iidx(m,i,j) .and. &
                 blk_tmp(c) == r% iblk(m,i,j)       ) exit
          end do
          if (m > mm) mm = m
          r% iidx(m,i,j) = idx_tmp(c)
          r% iblk(m,i,j) = blk_tmp(c)
          if (mine) then
           r% w_u(m,i,j) = r% w_u(m,i,j) + w_u_tmp(c)
           r% w_v(m,i,j) = r% w_v(m,i,j) + w_v_tmp(c)
          end if
       end do
!!$if (j == 1 .and. i < 10 .and. dace% lpio) then
!!$write(*,'(A,*(:,i8))') "i,idx_tmp=", i, idx_tmp
!!$end if
       if (mm /= ne) then
          if (dace% lpio) then
             write(0,*) "Internal error in stencil reduction:"
             write(0,*) "i,j   =", i, j
             write(0,*) "mm,ne =", mm, ne
             write(*,'(A,*(:,i8))') " idx_tmp=", idx_tmp
          end if
          call finish("init_stencil_icon_curl","### Internal Error")
       end if
       !------------------------------------------------
       ! Check if stencil crosses the processor boundary
       !------------------------------------------------
       do m = 1, ne
          jc  = r% iidx(m,i,j)
          jb  = r% iblk(m,i,j)
          pe2 = grid% marr(1,jc,jb,1)
          if      (      mine .and. pe2 /= dace% pe) then
             nrecv(pe2) = nrecv(pe2) + 1                ! Request from pe2
          else if (.not. mine .and. pe2 == dace% pe) then
             nsend(pe)  = nsend(pe)  + 1                ! Expected by pe
          end if
       end do
    end do
    end do

    if (dbg_level > 1) then
       do pe = 0, dace% npe-1
          call p_barrier()
          if (dace% pe == pe) print *, pe, "nsend    =", nsend
          if (dace% pe == pe) print *, pe, "nrecv    =", nrecv
       end do
       call p_barrier()
    end if
    if (nsend(dace% pe) /= 0 .or. nrecv(dace% pe) /= 0) then
       write(0,*) "### Error: nsend(self),nrecv(self) =",nsend(dace% pe),nrecv(dace% pe)
       call finish("init_stencil_icon_curl","### Internal Error")
    end if

    lbg = grid% lbg
    ubg = grid% ubg
    !--------------------------------------------------------
    ! Determine grid points for halo exchange
    ! 1) Points to send
    ! 2) Points to receive
    ! The current implementation does not care for duplicates
    !--------------------------------------------------------
    r_send% nhalo = sum (nsend)
    allocate (r_send% halo% iidx(r_send% nhalo))
    allocate (r_send% halo% iblk(r_send% nhalo))
    allocate (r_send% off       (0:dace% npe-1))
    allocate (r_send% cnt       (0:dace% npe-1))
!   allocate (r_send% hidx      (lbg(1):ubg(1),lbg(2):ubg(2)))
    r_send% halo% iidx = 0
    r_send% halo% iblk = 0
    r_send% off        = 0
    r_send% cnt        = 0
!   r_send% hidx       = 0
    !-----------------------
    ! Loop over receiver pes
    !-----------------------
    do pe = 0, dace% npe-1
       if (pe > 0) r_send% off(pe) = r_send% off(pe-1) + r_send% cnt(pe-1)
       if (nsend(pe) == 0) cycle
       off = r_send% off(pe)
       is  = 0
       !-----------------------
       ! Loop over all vertices
       !-----------------------
       do j    = lbv(2), ubv(2)
          do i = lbv(1), ubv(1)
             !---------------------------
             ! Vertex handled by other pe
             !---------------------------
             if (pe == grid% pe_d(i,j,1)) then
                !-------------------------------
                ! Loop over other stencil points
                ! Check whether needed remotely.
                !-------------------------------
                do m = 1, verts% num_edges(i,j)
                   jc = r% iidx(m,i,j)
                   jb = r% iblk(m,i,j)
                   if (grid% marr(1,jc,jb,1) == dace% pe) then
                      is = is + 1
                      r_send% halo% iidx(off+is) = jc
                      r_send% halo% iblk(off+is) = jb
!                     r_send% hidx      (jc,jb)  = off+is
                   end if
                end do
             end if
          end do
       end do
       r_send% cnt (pe) = is
       if (dbg_level > 1) &
            write(0,*) "p_pe -> pe, is, nsend(pe)", dace% pe, pe, is, nsend(pe)
!      if (is > nsend(pe)) then
       if (is /= nsend(pe)) then
          write(0,*) "### Fatal Error: is > nsend(pe)", is, nsend(pe)
          call finish("init_stencil_icon_curl","### Fatal Error: is > nsend")
       end if
    end do
    call p_barrier ()

    r_recv% nhalo = sum (nrecv)
    allocate (r_recv% halo% iidx(r_recv% nhalo))
    allocate (r_recv% halo% iblk(r_recv% nhalo))
    allocate (r_recv% off       (0:dace% npe-1))
    allocate (r_recv% cnt       (0:dace% npe-1))
    allocate (r_recv% hidx      (lbg(1):ubg(1),lbg(2):ubg(2)))
    r_recv% halo% iidx = 0
    r_recv% halo% iblk = 0
    r_recv% off        = 0
    r_recv% cnt        = 0
    r_recv% hidx       = 0
    !---------------------
    ! Loop over sender pes
    !---------------------
    do pe = 0, dace% npe-1
       if (pe > 0) r_recv% off(pe) = r_recv% off(pe-1) + r_recv% cnt(pe-1)
       if (nrecv(pe) == 0) cycle
       off = r_recv% off(pe)
       ir  = 0
       !------------------------
       ! Loop over "my" vertices
       !------------------------
       do j    = lb(2), ub(2)
          do i = lb(1), ub(1)
             !--------------------------
             ! Vertex handled by this pe
             !--------------------------
             if (dace% pe == grid% pe_d(i,j,1)) then
                !-------------------------------
                ! Loop over other stencil points
                !-------------------------------
                do m = 1, verts% num_edges(i,j)
                   jc = r% iidx(m,i,j)
                   jb = r% iblk(m,i,j)
                   if (grid% marr(1,jc,jb,1) == pe) then
                      ir = ir + 1
                      r_recv% halo% iidx(off+ir) = jc
                      r_recv% halo% iblk(off+ir) = jb
                      r_recv% hidx      (jc,jb)  = off+ir
                   end if
                end do
             end if
          end do
       end do
       r_recv% cnt (pe) = ir
       if (dbg_level > 1) &
            write(0,*) "p_pe <- pe, ir, nrecv(pe)", dace% pe, pe, ir, nrecv(pe)
!      if (ir > nrecv(pe)) then
       if (ir /= nrecv(pe)) then
          write(0,*) "### Fatal Error: ir > nrecv(pe)", ir, nrecv(pe)
          call finish("init_stencil_icon_curl","### Fatal Error: ir > nrecv")
       end if
    end do

    ! TODO: Exchange halo point lists for cross-checks

  end subroutine init_stencil_icon_curl
  !============================================================================
  subroutine destruct_halo (h)
    type(t_halo_adm), intent(out) :: h
  end subroutine destruct_halo
  !============================================================================
  subroutine destruct_stencil (s)
    type(t_stencil), intent(out) :: s
  end subroutine destruct_stencil
  !============================================================================
  subroutine curl_icon (u, v, rot, grid)
    real(wp), pointer, intent(in)  :: u  (:,:,:,:)  ! Zonal      component
    real(wp), pointer, intent(in)  :: v  (:,:,:,:)  ! Meridional component
    real(wp), pointer, intent(in)  :: rot(:,:,:,:)  ! Curl
    type(t_grid),      intent(in)  :: grid

    integer                        :: jc, jb        ! Cell indices
    integer                        :: i, j, k       ! Loop indices
    integer                        :: m             ! Stencil point index
    integer                        :: ix            ! Halo index
    type(t_grid_vertices), pointer :: verts

    if (.not. associated (grid% icongrid)) &
         call finish ("curl_icon", "icongrid not set up")
    if (.not. grid% global) &
         call finish ("curl_icon", "non-global grid not yet implemented")
    if (.not. allocated (icon_intp)) &
         call finish ("curl_icon", "icon_intp not set up")
    if (size (u,dim=3) /= size (rot,dim=3) .or. &
        size (v,dim=3) /= size (rot,dim=3)      ) then
       write(0,*) dace% pe, "curl_icon: incompatible dimension 3: u,v,rot=", &
            size (u,dim=3), size (v,dim=3), size (rot,dim=3)
       call finish ("curl_icon", "incompatible dimension 3")
    end if

    call halo_exchange_init (r_send, r_recv, u, v)

    verts => grid% icongrid% patch% verts
    !--------------
    ! Apply stencil
    !--------------
    rot(:,:,:,:) = 0._wp
    do       j = lbound (rot,2), ubound (rot,2)
       do    i = lbound (rot,1), ubound (rot,1)
          do k = lbound (rot,3), ubound (rot,3)
             do m = 1, verts% num_edges(i,j)
                jc = r% iidx(m,i,j)
                jb = r% iblk(m,i,j)
                ix = r_recv% hidx(jc,jb)
                if (ix == 0) then
                   ! Non-halo point (same pe)
                   rot(i,j,k,1) = rot(i,j,k,1)                 &
                                + r% w_u(m,i,j) * u(jc,jb,k,1) &
                                + r% w_v(m,i,j) * v(jc,jb,k,1)
                else
                   ! Halo point
                   rot(i,j,k,1) = rot(i,j,k,1)                          &
                                + r% w_u(m,i,j) * r_recv% halo% u(k,ix) &
                                + r% w_v(m,i,j) * r_recv% halo% v(k,ix)
                end if
             end do
          end do
       end do
    end do

    call halo_exchange_finish (r_send, r_recv, u, v)

  end subroutine curl_icon
  !============================================================================
  subroutine curl_icon_t (u, v, rot, grid)
    !-----------------------
    ! Transpose of curl_icon
    !-----------------------
    real(wp), pointer, intent(in)  :: u  (:,:,:,:)  ! Zonal      component
    real(wp), pointer, intent(in)  :: v  (:,:,:,:)  ! Meridional component
    real(wp), pointer, intent(in)  :: rot(:,:,:,:)  ! Curl
    type(t_grid),      intent(in)  :: grid

    integer                        :: jc, jb        ! Cell indices
    integer                        :: i, j, k       ! Loop indices
    integer                        :: m             ! Stencil point index
    integer                        :: ix            ! Halo index
    type(t_grid_vertices), pointer :: verts

    if (.not. associated (grid% icongrid)) &
         call finish ("curl_icon_t", "icongrid not set up")
    if (.not. grid% global) &
         call finish ("curl_icon_t", "non-global grid not yet implemented")
    if (.not. allocated (icon_intp)) &
         call finish ("curl_icon_t", "icon_intp not set up")
    if (size (u,dim=3) /= size (rot,dim=3) .or. &
        size (v,dim=3) /= size (rot,dim=3)      ) then
       write(0,*) dace% pe, "curl_icon_t: incompatible dimension 3: u,v,rot=", &
            size (u,dim=3), size (v,dim=3), size (rot,dim=3)
       call finish ("curl_icon_t", "incompatible dimension 3")
    end if

    call halo_exchange_init_t (r_send, r_recv, u, v)

    verts => grid% icongrid% patch% verts
    !-------------------------
    ! Apply transposed stencil
    !-------------------------
!   u(:,:,:,:) = 0._wp
!   v(:,:,:,:) = 0._wp
    do       j = lbound (rot,2), ubound (rot,2)
       do    i = lbound (rot,1), ubound (rot,1)
          do k = lbound (rot,3), ubound (rot,3)
             do m = 1, verts% num_edges(i,j)
                jc = r% iidx(m,i,j)
                jb = r% iblk(m,i,j)
                ix = r_recv% hidx(jc,jb)
                if (ix == 0) then
                   ! Update non-halo point (same pe)
                   u(jc,jb,k,1) = u(jc,jb,k,1)               &
                                + r% w_u(m,i,j) * rot(i,j,k,1)
                   v(jc,jb,k,1) = v(jc,jb,k,1)               &
                                + r% w_v(m,i,j) * rot(i,j,k,1)
                else
                   ! Update halo point
                   r_recv% halo% u(k,ix) = r_recv% halo% u(k,ix)      &
                                         + r% w_u(m,i,j) * rot(i,j,k,1)
                   r_recv% halo% v(k,ix) = r_recv% halo% v(k,ix)      &
                                         + r% w_v(m,i,j) * rot(i,j,k,1)
                end if
             end do
          end do
       end do
    end do

    call halo_exchange_finish_t (r_send, r_recv, u, v)

  end subroutine curl_icon_t
  !============================================================================
  subroutine div_latlon (u, v, div, grid)
    real(wp), pointer, intent(in)  :: u  (:,:,:,:)  ! -> Zonal      component
    real(wp), pointer, intent(in)  :: v  (:,:,:,:)  ! -> Meridional component
    real(wp), pointer, intent(in)  :: div(:,:,:,:)  ! <- Divergence
    type(t_grid),      intent(in)  :: grid

    integer                    :: i, j, k        ! Loop indices
    integer                    :: nbl            ! 0 or nboundlines
    integer                    :: l1, l2, l3     ! Effective local  lower bounds
    integer                    :: u1, u2, u3     ! Effective local  upper bounds
    integer                    :: l1g, l2g       ! Effective global lower bounds
    integer                    :: u1g, u2g       ! Effective global upper bounds
    integer                    :: lb (4), ub (4) ! Loop bounds (local)
    integer                    :: lbg(4), ubg(4) ! Loop bounds (global)
    integer                    :: pe_send        ! Sending pe
    integer                    :: pe_recv        ! Receiving pe
    real(wp)                   :: pre            ! Temporary
    real(wp)                   :: phi            ! Latitude (rad)
    real(wp)                   :: dlam           ! Delta lambda (rad)
    real(wp)                   :: dphi           ! Delta phi    (rad)
    real(wp), allocatable      :: w_u  (:,:)     ! Stencil weights
    real(wp), allocatable      :: w_v  (:,:)     ! Stencil weights
    real(wp), allocatable      :: uhalo(:,:)     ! Halo points (west)
    real(wp), allocatable      :: vhalo(:,:)     ! Halo points (south)
    real(wp), allocatable      :: usend(:,:)     ! Halo points to send
    real(wp), allocatable      :: vsend(:,:)     ! Halo points to send

    if (all (grid% gridtype /= [ WMO6_LATLON, WMO6_ROTLL ])) &
         call finish ("div_latlon", "gridtype must be (rotated) LATLON")
    if (grid% arakawa /= "C") &
         call finish ("div_latlon", "not implemented: Arakawa="//grid% arakawa)
    if (grid% global) &
         call finish ("div_latlon", "global latlon grid not yet implemented")

    if (size (u,dim=3) /= size (div,dim=3) .or. &
        size (v,dim=3) /= size (div,dim=3)      ) then
       write(0,*) dace% pe, "div_latlon: incompatible dimension 3: u,v,div=", &
            size (u,dim=3), size (v,dim=3), size (div,dim=3)
       call finish ("div_latlon", "incompatible dimension 3")
    end if
    if (any (shape (u) /= shape (div)) .or. &
        any (shape (v) /= shape (div))      ) then
       write(0,*) dace% pe, "div_latlon: incompatible shapes:"
       write(0,*) dace% pe, "shape (u)   =", shape (u)
       write(0,*) dace% pe, "shape (v)   =", shape (v)
       write(0,*) dace% pe, "shape (div) =", shape (div)
       call finish ("div_latlon", "incompatible shapes")
    end if

    l3  = lbound (div,3)
    u3  = ubound (div,3)
    lbg = grid% lbg
    ubg = grid% ubg
    lb  = grid% lb
    ub  = grid% ub
    !--------------------------------------------------
    ! Derive effective grid bounds, taking into account
    ! the number of grid lines which are not prognostic
    !--------------------------------------------------
    if (grid% global) then
       nbl = 0
    else
       nbl = max (nboundlines, 1)
    end if
    l1g = lbg(1) + nbl
    l2g = lbg(2) + nbl
    u1g = ubg(1) - nbl
    u2g = ubg(2) - nbl
    l1  = max (lb(1), l1g)
    l2  = max (lb(2), l2g)
    u1  = min (ub(1), u1g)
    u2  = min (ub(2), u2g)
    !--------------------------------------------------------------------
    ! Stencil weights for divergence on a spherical shell (radius: R)
    ! Div(u,v) = 1/(R*cos(phi)) * [d(u)/d(lambda) + d(v*cos(phi))/d(phi)]
    ! Valid for Arakawa-C grid only.
    !--------------------------------------------------------------------
    ! Stencil [convention: u(i,j) = u_n(i+1/2,j), v(i,j) = v_n(i,j+1/2) ]
    !
    !         +-----  v(i,j)  -----+
    !      u(i-1,j)  div(i,j)   u(i,j)
    !         +----- v(i,j-1) -----+
    !--------------------------------------------------------------------
    allocate (w_u(2,l2:u2))
    allocate (w_v(2,l2:u2))
    dlam = grid% di * d2r
    dphi = grid% dj * d2r
    do j = l2, u2
       phi      = grid% dlat(j) * d2r
       pre      = 1._wp / (Rearth * cos (phi))
       w_u(1,j) =   pre / dlam
       w_u(2,j) = - pre / dlam
       w_v(1,j) =   pre / dphi * cos (phi + 0.5_wp * dphi)
       w_v(2,j) = - pre / dphi * cos (phi - 0.5_wp * dphi)
    end do
    !--------------
    ! Halo exchange
    !--------------
    allocate (uhalo(l2:u2,l3:u3))
    allocate (usend(l2:u2,l3:u3))
    allocate (vhalo(l1:u1,l3:u3))
    allocate (vsend(l1:u1,l3:u3))

    if (l2 <= u2) then
       !-------------------------------------------------
       ! Send halo data @ i=u1 for right neighbor @ u1+1:
       !-------------------------------------------------
       if (u1 == ub(1) .and. u1+1 >= l1g .and. u1+1 <= u1g) then
          pe_recv    = grid% marr(1,u1+1,l2,1)
          usend(:,:) = u(u1,l2:u2,:,1)             ! Rightmost line @ i=u1
          call p_isend (usend, pe_recv, 1)
       end if
       !------------------------------------------
       ! Receive data @ i=l1-1 from left neighbor:
       !------------------------------------------
       if      (l1   == lb(1) .and. l1 <= u1) then
          pe_send    = grid% marr(1,l1-1,l2,1)
          call p_recv  (uhalo, pe_send, 1)
       else if (l1-1 <= ub(1) .and. l1 <= u1) then
          uhalo(:,:) = u(l1-1,l2:u2,:,1)           ! Halo elements on same pe
       else
          uhalo(:,:) = HUGE (0._wp)                ! No halo elements available
       end if
       call p_wait ()
    end if

    if (l1 <= u1) then
       !------------------------------------------------
       ! Send halo data @ j=u2 to upper neighbor @ u2+1:
       !------------------------------------------------
       if (u2 == ub(2) .and. u2+1 >= l2g .and. u2+1 <= u2g) then
          pe_recv    = grid% marr(1,l1,u2+1,1)
          vsend(:,:) = v(l1:u1,u2,:,1)             ! Topmost line @ j=u2
          call p_isend (vsend, pe_recv, 2)
       end if
       !-------------------------------------------
       ! Receive data @ j=l2-1 from lower neighbor:
       !-------------------------------------------
       if      (l2   == lb(2) .and. l2 <= u2) then
          pe_send    = grid% marr(1,l1,l2-1,1)
          call p_recv  (vhalo, pe_send, 2)
       else if (l2-1 <= ub(2) .and. l2 <= u2) then
          vhalo(:,:) = v(l1:u1,l2-1,:,1)           ! Halo elements on same pe
       else
          vhalo(:,:) = HUGE (0._wp)                ! No halo elements available
       end if
       call p_wait ()
    end if

    !--------------
    ! Apply stencil
    !--------------
    div(:,:,:,:) = 0._wp
    if (l1 <= u1 .and. l2 <= u2) then
       do    k = l3  , u3
!dir$ ivdep
          do j = l2  , u2
            div(l1,j,k,1) = div(l1,j,k,1) + w_u(1,j) * u    (l1 ,j,k,1) &! i=l1
                                          + w_u(2,j) * uhalo(    j,k  )
          end do ! j
          do j = l2  , u2
!dir$ ivdep
          do i = l1+1, u1
             div(i,j,k,1) = div(i ,j,k,1) + w_u(1,j) * u    (i  ,j,k,1) &! i>l1
                                          + w_u(2,j) * u    (i-1,j,k,1)
          end do ! i
          end do ! j
       end do    ! k

       do    k = l3  , u3
!dir$ ivdep
          do i = l1  , u1
            div(i,l2,k,1) = div(i,l2,k,1) + w_v(1,l2) * v    (i,l2 ,k,1) &! j=l2
                                          + w_v(2,l2) * vhalo(i    ,k  )
          end do ! i
          do j = l2+1, u2
!dir$ ivdep
          do i = l1  , u1
             div(i,j,k,1) = div(i,j ,k,1) + w_v(1,j)  * v    (i,j  ,k,1) &! j>l2
                                          + w_v(2,j)  * v    (i,j-1,k,1)
          end do ! i
          end do ! j
       end do    ! k
    end if

  end subroutine div_latlon
  !============================================================================
  subroutine div_latlon_t (u, v, div, grid)
    !------------------------
    ! Transpose of div_latlon
    !------------------------
    real(wp), pointer, intent(in)  :: u  (:,:,:,:)  ! Zonal      component
    real(wp), pointer, intent(in)  :: v  (:,:,:,:)  ! Meridional component
    real(wp), pointer, intent(in)  :: div(:,:,:,:)  ! Divergence
    type(t_grid),      intent(in)  :: grid

    integer                    :: i, j, k        ! Loop indices
    integer                    :: nbl            ! 0 or nboundlines
    integer                    :: l1, l2, l3     ! Effective local  lower bounds
    integer                    :: u1, u2, u3     ! Effective local  upper bounds
    integer                    :: l1g, l2g       ! Effective global lower bounds
    integer                    :: u1g, u2g       ! Effective global upper bounds
    integer                    :: lb (4), ub (4) ! Loop bounds (local)
    integer                    :: lbg(4), ubg(4) ! Loop bounds (global)
    integer                    :: pe_send        ! Sending pe
    integer                    :: pe_recv        ! Receiving pe
    real(wp)                   :: pre            ! Temporary
    real(wp)                   :: phi            ! Latitude (rad)
    real(wp)                   :: dlam           ! Delta lambda (rad)
    real(wp)                   :: dphi           ! Delta phi    (rad)
    real(wp), allocatable      :: w_u  (:,:)     ! Stencil weights
    real(wp), allocatable      :: w_v  (:,:)     ! Stencil weights
    real(wp), allocatable      :: uhalo(:,:)     ! Halo points (west)
    real(wp), allocatable      :: vhalo(:,:)     ! Halo points (south)
    real(wp), allocatable      :: urecv(:,:)     ! Halo points to recv
    real(wp), allocatable      :: vrecv(:,:)     ! Halo points to recv

    if (all (grid% gridtype /= [ WMO6_LATLON, WMO6_ROTLL ])) &
         call finish ("div_latlon_t", "gridtype must be (rotated) LATLON")
    if (grid% arakawa /= "C") &
         call finish ("div_latlon_t", "not implemented: Arakawa="//grid% arakawa)
    if (grid% global) &
         call finish ("div_latlon_t", "global latlon grid not yet implemented")

    if (size (u,dim=3) /= size (div,dim=3) .or. &
        size (v,dim=3) /= size (div,dim=3)      ) then
       write(0,*) dace% pe, "div_latlon_t: incompatible dimension 3: u,v,div=", &
            size (u,dim=3), size (v,dim=3), size (div,dim=3)
       call finish ("div_latlon_t", "incompatible dimension 3")
    end if
    if (any (shape (u) /= shape (div)) .or. &
        any (shape (v) /= shape (div))      ) then
       write(0,*) dace% pe, "div_latlon_t: incompatible shapes:"
       write(0,*) dace% pe, "shape (u)   =", shape (u)
       write(0,*) dace% pe, "shape (v)   =", shape (v)
       write(0,*) dace% pe, "shape (div) =", shape (div)
       call finish ("div_latlon_t", "incompatible shapes")
    end if

    l3  = lbound (div,3)
    u3  = ubound (div,3)
    lbg = grid% lbg
    ubg = grid% ubg
    lb  = grid% lb
    ub  = grid% ub
    !--------------------------------------------------
    ! Derive effective grid bounds, taking into account
    ! the number of grid lines which are not prognostic
    !--------------------------------------------------
    if (grid% global) then
       nbl = 0
    else
       nbl = max (nboundlines, 1)
    end if
    l1g = lbg(1) + nbl
    l2g = lbg(2) + nbl
    u1g = ubg(1) - nbl
    u2g = ubg(2) - nbl
    l1  = max (lb(1), l1g)
    l2  = max (lb(2), l2g)
    u1  = min (ub(1), u1g)
    u2  = min (ub(2), u2g)
    !--------------------------------------------------------------------
    ! Stencil weights for divergence on a spherical shell (radius: R)
    ! Div(u,v) = 1/(R*cos(phi)) * [d(u)/d(lambda) + d(v*cos(phi))/d(phi)]
    ! Valid for Arakawa-C grid only.
    !--------------------------------------------------------------------
    ! Stencil [convention: u(i,j) = u_n(i+1/2,j), v(i,j) = v_n(i,j+1/2) ]
    !
    !         +-----  v(i,j)  -----+
    !      u(i-1,j)  div(i,j)   u(i,j)
    !         +----- v(i,j-1) -----+
    !--------------------------------------------------------------------
    allocate (w_u(2,l2:u2))
    allocate (w_v(2,l2:u2))
    dlam = grid% di * d2r
    dphi = grid% dj * d2r
    do j = l2, u2
       phi      = grid% dlat(j) * d2r
       pre      = 1._wp / (Rearth * cos (phi))
       w_u(1,j) =   pre / dlam
       w_u(2,j) = - pre / dlam
       w_v(1,j) =   pre / dphi * cos (phi + 0.5_wp * dphi)
       w_v(2,j) = - pre / dphi * cos (phi - 0.5_wp * dphi)
    end do
    !--------------
    ! Halo exchange
    !--------------
    allocate (uhalo(l2:u2,l3:u3))
    allocate (urecv(l2:u2,l3:u3))
    allocate (vhalo(l1:u1,l3:u3))
    allocate (vrecv(l1:u1,l3:u3))

    !-------------------------
    ! Apply transposed stencil
    !-------------------------
    uhalo(:,:) = 0._wp
    vhalo(:,:) = 0._wp
!   u(:,:,:,:) = 0._wp
!   v(:,:,:,:) = 0._wp
    if (l1 <= u1 .and. l2 <= u2) then
       do    k = l3  , u3
!dir$ ivdep
          do j = l2  , u2
             u (l1 ,j,k,1) = u (l1 ,j,k,1) + w_u(1,j) * div(l1,j,k,1) ! i=l1
             uhalo( j,k  ) = uhalo( j,k  ) + w_u(2,j) * div(l1,j,k,1)
          end do ! j
!dir$ ivdep
          do j = l2  , u2
          do i = l1+1, u1
             u (i  ,j,k,1) = u (i  ,j,k,1) + w_u(1,j) * div(i ,j,k,1) ! i>l1
             u (i-1,j,k,1) = u (i-1,j,k,1) + w_u(2,j) * div(i ,j,k,1)
          end do ! i
          end do ! j
       end do    ! k

       do    k = l3  , u3
!dir$ ivdep
          do i = l1  , u1
             v (i, l2,k,1) = v (i, l2,k,1) + w_v(1,l2) * div(i,l2,k,1) ! j=l2
             vhalo(i ,k  ) = vhalo(i ,k  ) + w_v(2,l2) * div(i,l2,k,1)
          end do ! i
          do j = l2+1, u2
!dir$ ivdep
          do i = l1  , u1
             v (i,j  ,k,1) = v (i,j  ,k,1) + w_v(1,j)  * div(i,j ,k,1) ! j>l2
             v (i,j-1,k,1) = v (i,j-1,k,1) + w_v(2,j)  * div(i,j ,k,1)
          end do ! i
          end do ! j
       end do    ! k
    end if

    if (l2 <= u2) then
       !------------------------------------------
       ! Send halo data @ i=l1-1 to left neighbor:
       !------------------------------------------
       if      (l1   == lb(1) .and. l1 <= u1) then
          pe_recv    = grid% marr(1,l1-1,l2,1)
          call p_isend (uhalo, pe_recv, 3)
       else if (l1-1 <= ub(1) .and. l1 <= u1) then
          !--------------------------------
          ! Update halo elements on same pe
          !--------------------------------
          u(l1-1,l2:u2,:,1) = u(l1-1,l2:u2,:,1) + uhalo(:,:)
       else
          ! No halo elements available
       end if
       !----------------------------------------------------
       ! Receive data for i=u1 from right neighbor @ i=u1+1:
       !----------------------------------------------------
       if (u1 == ub(1) .and. u1+1 >= l1g .and. u1+1 <= u1g) then
          pe_send    = grid% marr(1,u1+1,l2,1)
          call p_recv  (urecv, pe_send, 3)
          u(u1  ,l2:u2,:,1) = u(u1  ,l2:u2,:,1) + urecv(:,:)
       end if
       call p_wait ()
    end if

    if (l1 <= u1) then
       !-------------------------------------------
       ! Send halo data @ j=l2-1 to lower neighbor:
       !-------------------------------------------
       if      (l2   == lb(2) .and. l2 <= u2) then
          pe_recv    = grid% marr(1,l1,l2-1,1)
          call p_isend (vhalo, pe_recv, 4)
       else if (l2-1 <= ub(2) .and. l2 <= u2) then
          !--------------------------------
          ! Update halo elements on same pe
          !--------------------------------
          v(l1:u1,l2-1,:,1) = v(l1:u1,l2-1,:,1) + vhalo(:,:)
       else
          ! No halo elements available
       end if
       !----------------------------------------------------
       ! Receive data for j=u2 from upper neighbor @ j=u2+1:
       !----------------------------------------------------
       if (u2 == ub(2) .and. u2+1 >= l2g .and. u2+1 <= u2g) then
          pe_send    = grid% marr(1,l1,u2+1,1)
          call p_recv  (vrecv, pe_send, 4)
          v(l1:u1,u2  ,:,1) = v(l1:u1,u2  ,:,1) + vrecv(:,:)
       end if
       call p_wait ()
    end if

  end subroutine div_latlon_t
  !============================================================================
  subroutine curl_latlon (u, v, rot, grid)
    real(wp), pointer, intent(in)  :: u  (:,:,:,:)  ! -> Zonal      component
    real(wp), pointer, intent(in)  :: v  (:,:,:,:)  ! -> Meridional component
    real(wp), pointer, intent(in)  :: rot(:,:,:,:)  ! <- Curl (z-component)
    type(t_grid),      intent(in)  :: grid

    integer                    :: i, j, k        ! Loop indices
    integer                    :: nbl            ! 0 or nboundlines
    integer                    :: l1, l2, l3     ! Effective local  lower bounds
    integer                    :: u1, u2, u3     ! Effective local  upper bounds
    integer                    :: l1g, l2g       ! Effective global lower bounds
    integer                    :: u1g, u2g       ! Effective global upper bounds
    integer                    :: lb (4), ub (4) ! Loop bounds (local)
    integer                    :: lbg(4), ubg(4) ! Loop bounds (global)
    integer                    :: pe_send        ! Sending pe
    integer                    :: pe_recv        ! Receiving pe
    real(wp)                   :: pre            ! Temporary
    real(wp)                   :: phi            ! Latitude (rad)
    real(wp)                   :: dlam           ! Delta lambda (rad)
    real(wp)                   :: dphi           ! Delta phi    (rad)
    real(wp), allocatable      :: w_u  (:,:)     ! Stencil weights
    real(wp), allocatable      :: w_v  (:,:)     ! Stencil weights
    real(wp), allocatable      :: uhalo(:,:)     ! Halo points (west)
    real(wp), allocatable      :: vhalo(:,:)     ! Halo points (south)
    real(wp), allocatable      :: usend(:,:)     ! Halo points to send
    real(wp), allocatable      :: vsend(:,:)     ! Halo points to send

    if (all (grid% gridtype /= [ WMO6_LATLON, WMO6_ROTLL ])) &
         call finish ("curl_latlon", "gridtype must be (rotated) LATLON")
    if (grid% arakawa /= "C") &
         call finish ("curl_latlon", "not implemented: Arakawa="//grid% arakawa)
    if (grid% global) &
         call finish ("curl_latlon", "global latlon grid not yet implemented")

    if (size (u,dim=3) /= size (rot,dim=3) .or. &
        size (v,dim=3) /= size (rot,dim=3)      ) then
       write(0,*) dace% pe, "curl_latlon: incompatible dimension 3: u,v,rot=", &
            size (u,dim=3), size (v,dim=3), size (rot,dim=3)
       call finish ("curl_latlon", "incompatible dimension 3")
    end if
    if (any (shape (u) /= shape (rot)) .or. &
        any (shape (v) /= shape (rot))      ) then
       write(0,*) dace% pe, "curl_latlon: incompatible shapes:"
       write(0,*) dace% pe, "shape (u)   =", shape (u)
       write(0,*) dace% pe, "shape (v)   =", shape (v)
       write(0,*) dace% pe, "shape (rot) =", shape (rot)
       call finish ("curl_latlon", "incompatible shapes")
    end if

    l3  = lbound (rot,3)
    u3  = ubound (rot,3)
    lbg = grid% lbg
    ubg = grid% ubg
    lb  = grid% lb
    ub  = grid% ub
    !--------------------------------------------------
    ! Derive effective grid bounds, taking into account
    ! the number of grid lines which are not prognostic
    !--------------------------------------------------
    if (grid% global) then
       nbl = 0
    else
       nbl = max (nboundlines, 1)
    end if
    l1g = lbg(1) + nbl
    l2g = lbg(2) + nbl
    u1g = ubg(1) - nbl
    u2g = ubg(2) - nbl
    l1  = max (lb(1), l1g)
    l2  = max (lb(2), l2g)
    u1  = min (ub(1), u1g)
    u2  = min (ub(2), u2g)
    !---------------------------------------------------------------------
    ! Stencil weights for curl on a spherical shell (radius: R)
    ! Curl(u,v) = 1/(R*cos(phi)) * [d(v)/d(lambda) - d(u*cos(phi))/d(phi)]
    ! Valid for Arakawa-C grid only.
    !---------------------------------------------------------------------
    ! Stencil [convention: u(i,j) = u_n(i+1/2,j), v(i,j) = v_n(i,j+1/2) ]
    !
    !         |             u(i,j+1)               |
    !         +-- v(i,j) -- rot(i,j) -- v(i+1,j) --+
    !         |              u(i,j)                |
    !---------------------------------------------------------------------
    allocate (w_u(2,l2:u2))
    allocate (w_v(2,l2:u2))
    dlam = grid% di * d2r
    dphi = grid% dj * d2r
    do j = l2, u2
       phi      = grid% dlat(j) * d2r + 0.5_wp * dphi           ! lat(j+1/2)
       pre      = 1._wp / (Rearth * cos (phi))
       w_v(1,j) =   pre / dlam
       w_v(2,j) = - pre / dlam
       w_u(1,j) = - pre / dphi * cos (phi + 0.5_wp * dphi)
       w_u(2,j) =   pre / dphi * cos (phi - 0.5_wp * dphi)
    end do
    !--------------
    ! Halo exchange
    !--------------
    allocate (uhalo(l1:u1,l3:u3))
    allocate (usend(l1:u1,l3:u3))
    allocate (vhalo(l2:u2,l3:u3))
    allocate (vsend(l2:u2,l3:u3))

    if (l2 <= u2) then
       !--------------------------------------------------
       ! Send halo data @ i=l1 for left neighbor @ i=l1-1:
       !--------------------------------------------------
       if (l1 == lb(1) .and. l1-1 >= l1g .and. l1-1 <= u1g) then
          pe_recv    = grid% marr(1,l1-1,l2,1)
          vsend(:,:) = v(l1,l2:u2,:,1)             ! Leftmost line @ i=l1
          call p_isend (vsend, pe_recv, 1)
       end if
       !-------------------------------------------
       ! Receive data @ i=u1+1 from right neighbor:
       !-------------------------------------------
       if      (u1   == ub(1) .and. l1 <= u1) then
          pe_send    = grid% marr(1,u1+1,l2,1)
          call p_recv  (vhalo, pe_send, 1)
       else if (u1+1 >= lb(1) .and. l1 <= u1) then
          vhalo(:,:) = v(u1+1,l2:u2,:,1)           ! Halo elements on same pe
       else
          vhalo(:,:) = HUGE (0._wp)                ! No halo elements available
       end if
       call p_wait ()
    end if

    if (l1 <= u1) then
       !--------------------------------------------------
       ! Send halo data @ j=l2 to lower neighbor @ j=l2-1:
       !--------------------------------------------------
       if (l2 == lb(2) .and. l2-1 >= l2g .and. l2-1 <= u2g) then
          pe_recv    = grid% marr(1,l1,l2-1,1)
          usend(:,:) = u(l1:u1,l2,:,1)             ! Lowermost line @ j=l2
          call p_isend (usend, pe_recv, 2)
       end if
       !-------------------------------------------
       ! Receive data @ j=u2+1 from upper neighbor:
       !-------------------------------------------
       if      (u2   == ub(2) .and. l2 <= u2) then
          pe_send    = grid% marr(1,l1,u2+1,1)
          call p_recv  (uhalo, pe_send, 2)
       else if (u2+1 >= lb(2) .and. l2 <= u2) then
          uhalo(:,:) = u(l1:u1,u2+1,:,1)           ! Halo elements on same pe
       else
          uhalo(:,:) = HUGE (0._wp)                ! No halo elements available
       end if
       call p_wait ()
    end if

    !--------------
    ! Apply stencil
    !--------------
    rot(:,:,:,:) = 0._wp
    if (l1 <= u1 .and. l2 <= u2) then
       do    k = l3  , u3
          do j = l2  , u2
!dir$ ivdep
          do i = l1  , u1-1
             rot(i,j,k,1) = rot(i ,j,k,1) + w_v(1,j) * v    (i+1,j,k,1) &! i<u1
                                          + w_v(2,j) * v    (i  ,j,k,1)
          end do ! i
          end do ! j
!dir$ ivdep
          do j = l2  , u2
            rot(u1,j,k,1) = rot(u1,j,k,1) + w_v(1,j) * vhalo(    j,k  ) &! i=u1
                                          + w_v(2,j) * v    (u1 ,j,k,1)
          end do ! j
       end do    ! k

       do    k = l3  , u3
          do j = l2  , u2-1
!dir$ ivdep
          do i = l1  , u1
             rot(i,j,k,1) = rot(i,j ,k,1) + w_u(1,j)  * u    (i,j+1,k,1) &! j<u2
                                          + w_u(2,j)  * u    (i,j  ,k,1)
          end do ! i
          end do ! j
!dir$ ivdep
          do i = l1  , u1
            rot(i,u2,k,1) = rot(i,u2,k,1) + w_u(1,u2) * uhalo(i    ,k  ) &! j=u2
                                          + w_u(2,u2) * u    (i,u2 ,k,1)
          end do ! i
       end do    ! k
    end if

  end subroutine curl_latlon
  !============================================================================
  subroutine curl_latlon_t (u, v, rot, grid)
    !-------------------------
    ! Transpose of curl_latlon
    !-------------------------
    real(wp), pointer, intent(in)  :: u  (:,:,:,:)  ! Zonal      component
    real(wp), pointer, intent(in)  :: v  (:,:,:,:)  ! Meridional component
    real(wp), pointer, intent(in)  :: rot(:,:,:,:)  ! Curl
    type(t_grid),      intent(in)  :: grid

    integer                    :: i, j, k        ! Loop indices
    integer                    :: nbl            ! 0 or nboundlines
    integer                    :: l1, l2, l3     ! Effective local  lower bounds
    integer                    :: u1, u2, u3     ! Effective local  upper bounds
    integer                    :: l1g, l2g       ! Effective global lower bounds
    integer                    :: u1g, u2g       ! Effective global upper bounds
    integer                    :: lb (4), ub (4) ! Loop bounds (local)
    integer                    :: lbg(4), ubg(4) ! Loop bounds (global)
    integer                    :: pe_send        ! Sending pe
    integer                    :: pe_recv        ! Receiving pe
    real(wp)                   :: pre            ! Temporary
    real(wp)                   :: phi            ! Latitude (rad)
    real(wp)                   :: dlam           ! Delta lambda (rad)
    real(wp)                   :: dphi           ! Delta phi    (rad)
    real(wp), allocatable      :: w_u  (:,:)     ! Stencil weights
    real(wp), allocatable      :: w_v  (:,:)     ! Stencil weights
    real(wp), allocatable      :: uhalo(:,:)     ! Halo points (west)
    real(wp), allocatable      :: vhalo(:,:)     ! Halo points (south)
    real(wp), allocatable      :: urecv(:,:)     ! Halo points to recv
    real(wp), allocatable      :: vrecv(:,:)     ! Halo points to recv

    if (all (grid% gridtype /= [ WMO6_LATLON, WMO6_ROTLL ])) &
         call finish ("curl_latlon_t", "gridtype must be (rotated) LATLON")
    if (grid% arakawa /= "C") &
         call finish ("curl_latlon_t", "not implemented: Arakawa="//grid% arakawa)
    if (grid% global) &
         call finish ("curl_latlon_t", "global latlon grid not yet implemented")

    if (size (u,dim=3) /= size (rot,dim=3) .or. &
        size (v,dim=3) /= size (rot,dim=3)      ) then
       write(0,*) dace% pe, "curl_latlon_t: incompatible dimension 3: u,v,rot=", &
            size (u,dim=3), size (v,dim=3), size (rot,dim=3)
       call finish ("curl_latlon_t", "incompatible dimension 3")
    end if
    if (any (shape (u) /= shape (rot)) .or. &
        any (shape (v) /= shape (rot))      ) then
       write(0,*) dace% pe, "curl_latlon_t: incompatible shapes:"
       write(0,*) dace% pe, "shape (u)   =", shape (u)
       write(0,*) dace% pe, "shape (v)   =", shape (v)
       write(0,*) dace% pe, "shape (rot) =", shape (rot)
       call finish ("curl_latlon_t", "incompatible shapes")
    end if

    l3  = lbound (rot,3)
    u3  = ubound (rot,3)
    lbg = grid% lbg
    ubg = grid% ubg
    lb  = grid% lb
    ub  = grid% ub
    !--------------------------------------------------
    ! Derive effective grid bounds, taking into account
    ! the number of grid lines which are not prognostic
    !--------------------------------------------------
    if (grid% global) then
       nbl = 0
    else
       nbl = max (nboundlines, 1)
    end if
    l1g = lbg(1) + nbl
    l2g = lbg(2) + nbl
    u1g = ubg(1) - nbl
    u2g = ubg(2) - nbl
    l1  = max (lb(1), l1g)
    l2  = max (lb(2), l2g)
    u1  = min (ub(1), u1g)
    u2  = min (ub(2), u2g)
    !---------------------------------------------------------------------
    ! Stencil weights for curl on a spherical shell (radius: R)
    ! Curl(u,v) = 1/(R*cos(phi)) * [d(v)/d(lambda) - d(u*cos(phi))/d(phi)]
    ! Valid for Arakawa-C grid only.
    !---------------------------------------------------------------------
    ! Stencil [convention: u(i,j) = u_n(i+1/2,j), v(i,j) = v_n(i,j+1/2) ]
    !
    !         |             u(i,j+1)               |
    !         +-- v(i,j) -- rot(i,j) -- v(i+1,j) --+
    !         |              u(i,j)                |
    !---------------------------------------------------------------------
    allocate (w_u(2,l2:u2))
    allocate (w_v(2,l2:u2))
    dlam = grid% di * d2r
    dphi = grid% dj * d2r
    do j = l2, u2
       phi      = grid% dlat(j) * d2r + 0.5_wp * dphi           ! lat(j+1/2)
       pre      = 1._wp / (Rearth * cos (phi))
       w_v(1,j) =   pre / dlam
       w_v(2,j) = - pre / dlam
       w_u(1,j) = - pre / dphi * cos (phi + 0.5_wp * dphi)
       w_u(2,j) =   pre / dphi * cos (phi - 0.5_wp * dphi)
    end do
    !--------------
    ! Halo exchange
    !--------------
    allocate (uhalo(l1:u1,l3:u3))
    allocate (urecv(l1:u1,l3:u3))
    allocate (vhalo(l2:u2,l3:u3))
    allocate (vrecv(l2:u2,l3:u3))

    !-------------------------
    ! Apply transposed stencil
    !-------------------------
    uhalo(:,:) = 0._wp
    vhalo(:,:) = 0._wp
!   u(:,:,:,:) = 0._wp
!   v(:,:,:,:) = 0._wp
    if (l1 <= u1 .and. l2 <= u2) then
       do    k = l3  , u3
!dir$ ivdep
          do j = l2  , u2
          do i = l1  , u1-1
             v (i+1,j,k,1) = v (i+1,j,k,1) + w_v(1,j) * rot(i ,j,k,1) ! i<u1
             v (i  ,j,k,1) = v (i  ,j,k,1) + w_v(2,j) * rot(i ,j,k,1)
          end do ! i
          end do ! j
!dir$ ivdep
          do j = l2  , u2
             vhalo( j,k  ) = vhalo( j,k  ) + w_v(1,j) * rot(u1,j,k,1) ! i=u1
             v (u1 ,j,k,1) = v (u1 ,j,k,1) + w_v(2,j) * rot(u1,j,k,1)
          end do ! j
       end do    ! k

       do    k = l3  , u3
          do j = l2  , u2-1
!dir$ ivdep
          do i = l1  , u1
             u (i,j+1,k,1) = u (i,j+1,k,1) + w_u(1,j)  * rot(i,j ,k,1) ! j<u2
             u (i,j  ,k,1) = u (i,j  ,k,1) + w_u(2,j)  * rot(i,j ,k,1)
          end do ! i
          end do ! j
!dir$ ivdep
          do i = l1  , u1
             uhalo(i ,k  ) = uhalo(i ,k  ) + w_u(1,u2) * rot(i,u2,k,1) ! j=u2
             u (i, u2,k,1) = u (i, u2,k,1) + w_u(2,u2) * rot(i,u2,k,1)
          end do ! i
       end do    ! k
    end if

    if (l2 <= u2) then
       !-------------------------------------------
       ! Send halo data @ i=u1+1 to right neighbor:
       !-------------------------------------------
       if      (u1   == ub(1) .and. l1 <= u1) then
          pe_recv    = grid% marr(1,u1+1,l2,1)
          call p_isend (vhalo, pe_recv, 3)
       else if (u1+1 >= lb(1) .and. l1 <= u1) then
          !--------------------------------
          ! Update halo elements on same pe
          !--------------------------------
          v(u1+1,l2:u2,:,1) = v(u1+1,l2:u2,:,1) + vhalo(:,:)
       else
          ! No halo elements available
       end if
       !---------------------------------------------------
       ! Receive data for i=l1 from left neighbor @ i=l1-1:
       !---------------------------------------------------
       if (l1 == lb(1) .and. l1-1 >= l1g .and. l1-1 <= u1g) then
          pe_send    = grid% marr(1,l1-1,l2,1)
          call p_recv  (vrecv, pe_send, 3)
          v(l1  ,l2:u2,:,1) = v(l1  ,l2:u2,:,1) + vrecv(:,:)
       end if
       call p_wait ()
    end if

    if (l1 <= u1) then
       !-------------------------------------------
       ! Send halo data @ j=u2+1 to upper neighbor:
       !-------------------------------------------
       if      (u2   == ub(2) .and. l2 <= u2) then
          pe_recv    = grid% marr(1,l1,u2+1,1)
          call p_isend (uhalo, pe_recv, 4)
       else if (u2+1 >= lb(2) .and. l2 <= u2) then
          !--------------------------------
          ! Update halo elements on same pe
          !--------------------------------
          u(l1:u1,u2+1,:,1) = u(l1:u1,u2+1,:,1) + uhalo(:,:)
       else
          ! No halo elements available
       end if
       !----------------------------------------------------
       ! Receive data for j=l2 from lower neighbor @ j=l2-1:
       !----------------------------------------------------
       if (l2 == lb(2) .and. l2-1 >= l2g .and. l2-1 <= u2g) then
          pe_send    = grid% marr(1,l1,l2-1,1)
          call p_recv  (urecv, pe_send, 4)
          u(l1:u1,l2  ,:,1) = u(l1:u1,l2  ,:,1) + urecv(:,:)
       end if
       call p_wait ()
    end if

  end subroutine curl_latlon_t
  !============================================================================
  subroutine grad_latlon (chi, u, v, grid)
    real(wp), pointer, intent(in)  :: chi(:,:,:,:)  ! -> Scalar potential
    real(wp), pointer, intent(in)  :: u  (:,:,:,:)  ! <- Zonal      component
    real(wp), pointer, intent(in)  :: v  (:,:,:,:)  ! <- Meridional component
    type(t_grid),      intent(in)  :: grid

    integer                    :: i, j, k        ! Loop indices
    integer                    :: nbl            ! 0 or nboundlines
    integer                    :: l1, l2, l3     ! Effective local  lower bounds
    integer                    :: u1, u2, u3     ! Effective local  upper bounds
    integer                    :: l1g, l2g       ! Effective global lower bounds
    integer                    :: u1g, u2g       ! Effective global upper bounds
    integer                    :: lb (4), ub (4) ! Loop bounds (local)
    integer                    :: lbg(4), ubg(4) ! Loop bounds (global)
    integer                    :: pe_send        ! Sending pe
    integer                    :: pe_recv        ! Receiving pe
    real(wp)                   :: pre            ! Temporary
    real(wp)                   :: phi            ! Latitude (rad)
    real(wp)                   :: dlam           ! Delta lambda (rad)
    real(wp)                   :: dphi           ! Delta phi    (rad)
    real(wp), allocatable      :: w_u  (:,:)     ! Stencil weights
    real(wp), allocatable      :: w_v  (:,:)     ! Stencil weights
    real(wp), allocatable      :: uhalo(:,:)     ! Halo points (east)
    real(wp), allocatable      :: vhalo(:,:)     ! Halo points (north)
    real(wp), allocatable      :: usend(:,:)     ! Halo points to send
    real(wp), allocatable      :: vsend(:,:)     ! Halo points to send

    if (all (grid% gridtype /= [ WMO6_LATLON, WMO6_ROTLL ])) &
         call finish ("grad_latlon", "gridtype must be (rotated) LATLON")
    if (grid% arakawa /= "C") &
         call finish ("grad_latlon", "not implemented: Arakawa="//grid% arakawa)
    if (grid% global) &
         call finish ("grad_latlon", "global latlon grid not yet implemented")

    if (size (u,dim=3) /= size (chi,dim=3) .or. &
        size (v,dim=3) /= size (chi,dim=3)      ) then
       write(0,*) dace% pe, "grad_latlon: incompatible dimension 3: chi,u,v=", &
            size (chi,dim=3), size (u,dim=3), size (v,dim=3)
       call finish ("grad_latlon", "incompatible dimension 3")
    end if
    if (any (shape (u) /= shape (chi)) .or. &
        any (shape (v) /= shape (chi))      ) then
       write(0,*) dace% pe, "grad_latlon: incompatible shapes:"
       write(0,*) dace% pe, "shape (chi) =", shape (chi)
       write(0,*) dace% pe, "shape (u)   =", shape (u)
       write(0,*) dace% pe, "shape (v)   =", shape (v)
       call finish ("grad_latlon", "incompatible shapes")
    end if

    l3  = lbound (chi,3)
    u3  = ubound (chi,3)
    lbg = grid% lbg
    ubg = grid% ubg
    lb  = grid% lb
    ub  = grid% ub
    !--------------------------------------------------
    ! Derive effective grid bounds, taking into account
    ! the number of grid lines which are not prognostic
    !--------------------------------------------------
    if (grid% global) then
       nbl = 0
    else
       nbl = max (nboundlines, 1)
    end if
    l1g = lbg(1) + nbl
    l2g = lbg(2) + nbl
    u1g = ubg(1) - nbl
    u2g = ubg(2) - nbl
    l1  = max (lb(1), l1g)
    l2  = max (lb(2), l2g)
    u1  = min (ub(1), u1g)
    u2  = min (ub(2), u2g)
    !---------------------------------------------------------------------
    ! Weights for horizontal gradient on a spherical shell (radius: R)
    ! Grad(chi) = [1/(R*cos(phi)) * d(chi)/d(lambda), 1/R * d(chi)/d(phi)]
    ! Valid for Arakawa-C grid only.
    !---------------------------------------------------------------------
    ! Stencil [convention: u(i,j) = u_n(i+1/2,j), v(i,j) = v_n(i,j+1/2) ]
    !          staggering: psi "lives" on mass points.
    !
    !       chi(i,j+1) --------------------+
    !         v(i,j)                       |
    !       chi(i,j) ----- u(i,j) --- chi(i+1,j)
    !---------------------------------------------------------------------
    allocate (w_u(2,l2:u2))
    allocate (w_v(2,l2:u2))
    dlam = grid% di * d2r
    dphi = grid% dj * d2r
    do j = l2, u2
       phi      = grid% dlat(j) * d2r
       pre      = 1._wp / Rearth
       w_u(1,j) =   pre / (dlam * cos (phi))
       w_u(2,j) = - w_u(1,j)
       w_v(1,j) =   pre /  dphi
       w_v(2,j) = - w_v(1,j)
    end do
    !--------------
    ! Halo exchange
    !--------------
    allocate (uhalo(l2:u2,l3:u3))
    allocate (usend(l2:u2,l3:u3))
    allocate (vhalo(l1:u1,l3:u3))
    allocate (vsend(l1:u1,l3:u3))

    if (l2 <= u2) then
       !------------------------------------------------
       ! Send halo data @ i=l1 for left neighbor @ l1-1:
       !------------------------------------------------
       if (l1 == lb(1) .and. l1-1 >= l1g .and. l1-1 <= u1g) then
          pe_recv    = grid% marr(1,l1-1,l2,1)
          usend(:,:) = chi(l1,l2:u2,:,1)           ! Leftmost line @ i=l1
          call p_isend (usend, pe_recv, 1)
       end if
       !-------------------------------------------
       ! Receive data @ i=u1+1 from right neighbor:
       !-------------------------------------------
       if      (u1   == ub(1) .and. l1 <= u1) then
          pe_send    = grid% marr(1,u1+1,l2,1)
          call p_recv  (uhalo, pe_send, 1)
       else if (u1+1 >= lb(1) .and. l1 <= u1) then
          uhalo(:,:) = chi(u1+1,l2:u2,:,1)         ! Halo elements on same pe
       else
          uhalo(:,:) = HUGE (0._wp)                ! No halo elements available
       end if
       call p_wait ()
    end if

    if (l1 <= u1) then
       !------------------------------------------------
       ! Send halo data @ j=l2 to lower neighbor @ l2-1:
       !------------------------------------------------
       if (l2 == lb(2) .and. l2-1 >= l2g .and. l2-1 <= u2g) then
          pe_recv    = grid% marr(1,l1,l2-1,1)
          vsend(:,:) = chi(l1:u1,l2,:,1)           ! Lowermost line @ j=l2
          call p_isend (vsend, pe_recv, 2)
       end if
       !-------------------------------------------
       ! Receive data @ j=u2+1 from upper neighbor:
       !-------------------------------------------
       if      (u2   == ub(2) .and. l2 <= u2) then
          pe_send    = grid% marr(1,l1,u2+1,1)
          call p_recv  (vhalo, pe_send, 2)
       else if (u2+1 >= lb(2) .and. l2 <= u2) then
          vhalo(:,:) = chi(l1:u1,u2+1,:,1)         ! Halo elements on same pe
       else
          vhalo(:,:) = HUGE (0._wp)                ! No halo elements available
       end if
       call p_wait ()
    end if

    !--------------
    ! Apply stencil
    !--------------
    u(:,:,:,:) = 0._wp
    v(:,:,:,:) = 0._wp
    if (l1 <= u1 .and. l2 <= u2) then
       do    k = l3  , u3
          do j = l2  , u2
!dir$ ivdep
          do i = l1, u1-1
             u(i ,j,k,1) = w_u(1,j) * chi(i+1,j,k,1) + w_u(2,j) * chi(i ,j,k,1)
          end do ! i
          end do ! j
!dir$ ivdep
          do j = l2  , u2
             u(u1,j,k,1) = w_u(1,j) * uhalo(  j,k  ) + w_u(2,j) * chi(u1,j,k,1)
          end do ! j
       end do    ! k

       do    k = l3  , u3
          do j = l2  , u2-1
!dir$ ivdep
          do i = l1  , u1
             v(i,j ,k,1) = w_v(1,j) * chi(i,j+1,k,1) + w_v(2,j) * chi(i,j ,k,1)
          end do ! i
          end do ! j
!dir$ ivdep
          do i = l1  , u1
             v(i,u2,k,1) = w_v(1,j) * vhalo(i  ,k  ) + w_v(2,j) * chi(i,u2,k,1)
          end do ! i
       end do    ! k
    end if

  end subroutine grad_latlon
  !============================================================================
  subroutine hrot_latlon (psi, u, v, grid)
    real(wp), pointer, intent(in)  :: psi(:,:,:,:)  ! -> Stream function
    real(wp), pointer, intent(in)  :: u  (:,:,:,:)  ! <- Zonal      component
    real(wp), pointer, intent(in)  :: v  (:,:,:,:)  ! <- Meridional component
    type(t_grid),      intent(in)  :: grid

    integer                    :: i, j, k        ! Loop indices
    integer                    :: nbl            ! 0 or nboundlines
    integer                    :: l1, l2, l3     ! Effective local  lower bounds
    integer                    :: u1, u2, u3     ! Effective local  upper bounds
    integer                    :: l1g, l2g       ! Effective global lower bounds
    integer                    :: u1g, u2g       ! Effective global upper bounds
    integer                    :: lb (4), ub (4) ! Loop bounds (local)
    integer                    :: lbg(4), ubg(4) ! Loop bounds (global)
    integer                    :: pe_send        ! Sending pe
    integer                    :: pe_recv        ! Receiving pe
    real(wp)                   :: pre            ! Temporary
    real(wp)                   :: phi            ! Latitude (rad)
    real(wp)                   :: dlam           ! Delta lambda (rad)
    real(wp)                   :: dphi           ! Delta phi    (rad)
    real(wp), allocatable      :: w_u  (:,:)     ! Stencil weights
    real(wp), allocatable      :: w_v  (:,:)     ! Stencil weights
    real(wp), allocatable      :: uhalo(:,:)     ! Halo points (west)
    real(wp), allocatable      :: vhalo(:,:)     ! Halo points (south)
    real(wp), allocatable      :: usend(:,:)     ! Halo points to send
    real(wp), allocatable      :: vsend(:,:)     ! Halo points to send

    if (all (grid% gridtype /= [ WMO6_LATLON, WMO6_ROTLL ])) &
         call finish ("hrot_latlon", "gridtype must be (rotated) LATLON")
    if (grid% arakawa /= "C") &
         call finish ("hrot_latlon", "not implemented: Arakawa="//grid% arakawa)
    if (grid% global) &
         call finish ("hrot_latlon", "global latlon grid not yet implemented")

    if (size (u,dim=3) /= size (psi,dim=3) .or. &
        size (v,dim=3) /= size (psi,dim=3)      ) then
       write(0,*) dace% pe, "hrot_latlon: incompatible dimension 3: psi,u,v=", &
            size (psi,dim=3), size (u,dim=3), size (v,dim=3)
       call finish ("hrot_latlon", "incompatible dimension 3")
    end if
    if (any (shape (u) /= shape (psi)) .or. &
        any (shape (v) /= shape (psi))      ) then
       write(0,*) dace% pe, "hrot_latlon: incompatible shapes:"
       write(0,*) dace% pe, "shape (psi) =", shape (psi)
       write(0,*) dace% pe, "shape (u)   =", shape (u)
       write(0,*) dace% pe, "shape (v)   =", shape (v)
       call finish ("hrot_latlon", "incompatible shapes")
    end if

    l3  = lbound (psi,3)
    u3  = ubound (psi,3)
    lbg = grid% lbg
    ubg = grid% ubg
    lb  = grid% lb
    ub  = grid% ub
    !--------------------------------------------------
    ! Derive effective grid bounds, taking into account
    ! the number of grid lines which are not prognostic
    !--------------------------------------------------
    if (grid% global) then
       nbl = 0
    else
       nbl = max (nboundlines, 1)
    end if
    l1g = lbg(1) + nbl
    l2g = lbg(2) + nbl
    u1g = ubg(1) - nbl
    u2g = ubg(2) - nbl
    l1  = max (lb(1), l1g)
    l2  = max (lb(2), l2g)
    u1  = min (ub(1), u1g)
    u2  = min (ub(2), u2g)
    !---------------------------------------------------------------------
    ! Weights for "horizontal rotation" on a spherical shell (radius: R)
    ! Rot(psi) = [-1/R * d(psi)/d(phi), 1/(R*cos(phi)) * d(psi)/d(lambda)]
    ! Valid for Arakawa-C grid only.  Sign conventions from meteorology.
    !---------------------------------------------------------------------
    ! Stencil [convention: u(i,j) = u_n(i+1/2,j), v(i,j) = v_n(i,j+1/2) ]
    !          staggering: psi(i,j) ~ Psi(i+1/2,j+1/2)
    !
    !   |     psi(i-1,j)   v(i,j)    psi(i,j)       |
    !   +---------------------+------- u(i,j) ------+
    !   |                     |      psi(i,j-1)     |
    !---------------------------------------------------------------------
    allocate (w_u(2,l2:u2))
    allocate (w_v(2,l2:u2))
    dlam = grid% di * d2r
    dphi = grid% dj * d2r
    do j = l2, u2
       phi      = grid% dlat(j) * d2r
       pre      = 1._wp / Rearth
       w_u(1,j) = - pre /  dphi
       w_u(2,j) = - w_u(1,j)
       w_v(1,j) =   pre / (dlam * cos (phi + 0.5_wp * dphi))
       w_v(2,j) = - w_v(1,j)
    end do
    !--------------
    ! Halo exchange
    !--------------
    allocate (uhalo(l2:u2,l3:u3))
    allocate (usend(l2:u2,l3:u3))
    allocate (vhalo(l1:u1,l3:u3))
    allocate (vsend(l1:u1,l3:u3))

    if (l2 <= u2) then
       !-------------------------------------------------
       ! Send halo data @ i=u1 for right neighbor @ u1+1:
       !-------------------------------------------------
       if (u1 == ub(1) .and. u1+1 >= l1g .and. u1+1 <= u1g) then
          pe_recv    = grid% marr(1,u1+1,l2,1)
          usend(:,:) = psi(u1,l2:u2,:,1)           ! Rightmost line @ i=u1
          call p_isend (usend, pe_recv, 1)
       end if
       !------------------------------------------
       ! Receive data @ i=l1-1 from left neighbor:
       !------------------------------------------
       if      (l1   == lb(1) .and. l1 <= u1) then
          pe_send    = grid% marr(1,l1-1,l2,1)
          call p_recv  (uhalo, pe_send, 1)
       else if (l1-1 <= ub(1) .and. l1 <= u1) then
          uhalo(:,:) = psi(l1-1,l2:u2,:,1)         ! Halo elements on same pe
       else
          uhalo(:,:) = HUGE (0._wp)                ! No halo elements available
       end if
       call p_wait ()
    end if

    if (l1 <= u1) then
       !------------------------------------------------
       ! Send halo data @ j=u2 to upper neighbor @ u2+1:
       !------------------------------------------------
       if (u2 == ub(2) .and. u2+1 >= l2g .and. u2+1 <= u2g) then
          pe_recv    = grid% marr(1,l1,u2+1,1)
          vsend(:,:) = psi(l1:u1,u2,:,1)           ! Topmost line @ j=u2
          call p_isend (vsend, pe_recv, 2)
       end if
       !-------------------------------------------
       ! Receive data @ j=l2-1 from lower neighbor:
       !-------------------------------------------
       if      (l2   == lb(2) .and. l2 <= u2) then
          pe_send    = grid% marr(1,l1,l2-1,1)
          call p_recv  (vhalo, pe_send, 2)
       else if (l2-1 <= ub(2) .and. l2 <= u2) then
          vhalo(:,:) = psi(l1:u1,l2-1,:,1)         ! Halo elements on same pe
       else
          vhalo(:,:) = HUGE (0._wp)                ! No halo elements available
       end if
       call p_wait ()
    end if

    !--------------
    ! Apply stencil
    !--------------
    u(:,:,:,:) = 0._wp
    v(:,:,:,:) = 0._wp
    if (l1 <= u1 .and. l2 <= u2) then
       do    k = l3  , u3
!dir$ ivdep
          do i = l1  , u1
             u(i,l2,k,1) = w_u(1,l2) * psi(i,l2,k,1) + w_u(2,l2) * vhalo( i,k  )
          end do ! i
          do j = l2+1, u2
!dir$ ivdep
          do i = l1  , u1
             u(i,j ,k,1) = w_u(1,j ) * psi(i,j ,k,1) + w_u(2,j) * psi(i,j-1,k,1)
          end do ! i
          end do ! j
       end do    ! k

       do    k = l3  , u3
!dir$ ivdep
          do j = l2  , u2
             v(l1,j,k,1) = w_v(1,j) * psi(l1,j,k,1) + w_v(2,j) * uhalo(  j,k  )
          end do ! j
          do j = l2  , u2
!dir$ ivdep
          do i = l1+1, u1
             v(i ,j,k,1) = w_v(1,j) * psi(i ,j,k,1) + w_v(2,j) * psi(i-1,j,k,1)
          end do ! i
          end do ! j
       end do    ! k
    end if

  end subroutine hrot_latlon
  !============================================================================
END MODULE mo_math_oper
