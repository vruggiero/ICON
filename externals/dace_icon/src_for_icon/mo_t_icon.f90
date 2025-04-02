!
!+ Derived type definitions and basic function for ICON patches
!
! $Id$
!
MODULE mo_t_icon
!
! Description:
!   Derived type definitions and basic function for ICON patches.
!   These are mostly taken from prep_icon (mo_utilities, mo_remap_shared).
!
! Current Code Owner: DWD, Harald Anlauf
!    phone: +49 69 8062 4941
!    fax:   +49 69 8062 3721
!    email: harald.anlauf@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_23        2013-03-26 Harald Anlauf
!  Initial version
! V1_26        2013/06/27 Harald Anlauf
!  Optimize loading ICON grid metadata for NEC SX-9
! V1_48        2016-10-06 Harald Anlauf
!  precomputation of coefficients cells<->edges; read data for curl operator
!
! Code Description:
! Language: Fortran 95.
! Software Standards:
!
!------------------------------------------------------------------------------
  !-------------
  ! Modules used
  !-------------
  use mo_kind,          only: wp, i2    ! working precision, short integer kind
  implicit none

  !----------------
  ! Public entities
  !----------------
  private
  public :: nproma                      ! Vector length for blocking
  public :: blk_no                      ! Block index from global index
  public :: idx_no                      ! Line  index within block
  public :: normalized_coord            ! Reduce lon/lat to [-180,180]/[-90,90]
  public :: gc2cc                       ! Geographical to cartesian coordinates
  public :: gvec2cvec                   ! Zonal/meridional vector to cartesian
  public :: cvec2gvec                   ! Cartesian vector to zonal/meridional
  public :: t_patch                     ! Derived type for ICON patch
  public :: t_grid_cells                ! Derived type for ICON grid cells
  public :: t_grid_vertices             ! Derived type for ICON grid vertices
  public :: t_grid_edges                ! Derived type for ICON grid edges
  public :: t_geographical_coordinates  ! Derived type (lon,lat)
  public :: t_cartesian_coordinates     ! Derived type (x(3))
  public :: t_tangent_vectors           ! Derived type (v1,v2)
  public :: t_int_state                 ! Derived type: interpolation coeffs.

  !------------------------
! integer :: nproma = 1         ! nproma blocking
  integer :: nproma = HUGE(0)   ! Default: no blocking
  !------------------------
!#ifdef __SX__
!  ! sxf90 rev.451:
!  ! Workaround for "Internal compiler error in optimization phase"
!#define ALLOCATABLE POINTER
!#endif
  !------------------------

  TYPE t_geographical_coordinates
    REAL(wp) :: lon                     ! -pi   .. pi   [rad]
    REAL(wp) :: lat                     ! -pi/2 .. pi/2 [rad]
  END TYPE t_geographical_coordinates

  TYPE t_cartesian_coordinates
    REAL(wp) :: x(3)
  END TYPE t_cartesian_coordinates

  TYPE t_tangent_vectors
    REAL(wp) :: v1
    REAL(wp) :: v2
  END TYPE t_tangent_vectors


  TYPE t_grid_cells

    ! line indices of edges of triangle:
    ! index1=1,nproma, index2=1,nblks_c, index3=1,cell_type
    INTEGER, ALLOCATABLE :: edge_idx(:,:,:)
    ! block indices of edges of triangle:
    ! index1=1,nproma, index2=1,nblks_c, index3=1,cell_type
    INTEGER, ALLOCATABLE :: edge_blk(:,:,:)

    ! line indices of verts of triangle:
    ! index1=1,nproma, index2=1,nblks_c, index3=1,cell_type
    INTEGER, ALLOCATABLE :: vertex_idx(:,:,:)
    ! block indices of verts of triangle:
    ! index1=1,nproma, index2=1,nblks_c, index3=1,cell_type
    INTEGER, ALLOCATABLE :: vertex_blk(:,:,:)

    ! longitude & latitude of centers of triangular cells
    ! index1=nproma, index2=1,nblks_c
    TYPE(t_geographical_coordinates), ALLOCATABLE :: center(:,:)

    ! area of cell
    ! index1=nproma, index2=1,nblks_c
    REAL(wp), ALLOCATABLE :: area(:,:)

    ! refinement control flag for cells
    ! index1=nproma, index2=1,nblks_c
    INTEGER, ALLOCATABLE :: c_ctrl(:,:)

    ! line/block indices of elements  next to each cell:
    ! index1=1,nproma, index2=1,nblks_c, index3=1,cell_type
    INTEGER, ALLOCATABLE :: neighbor_idx(:,:,:), neighbor_blk(:,:,:)

    ! global index of each cell: index=1,n_patch_cells
    INTEGER, ALLOCATABLE :: glb_index(:)

    ! dummy (unused in this code)
    INTEGER, ALLOCATABLE :: start_blk(:,:), end_blk(:,:)
    LOGICAL, ALLOCATABLE :: owner_mask(:,:)

    ! orientation (+/-1)
    ! index1=1,nproma, index2=1,nblks_c, index3=1,cell_type
    integer(i2), allocatable :: edge_orientation(:,:,:)
  END TYPE t_grid_cells


  TYPE t_grid_edges

    ! line indices of adjacent cells:
    ! index1=1,nproma, index2=1,nblks_e, index3=1,2
    INTEGER, ALLOCATABLE :: cell_idx(:,:,:)
    ! block indices of adjacent cells:
    ! index1=1,nproma, index2=1,nblks_e, index3=1,2
    INTEGER, ALLOCATABLE :: cell_blk(:,:,:)

    ! line indices of edge vertices:
    ! index1=1,nproma, index2=1,nblks_e, index3=1,2/4
    INTEGER, ALLOCATABLE :: vertex_idx(:,:,:)
    ! block indices of edge vertices:
    ! index1=1,nproma, index2=1,nblks_e, index3=1,2/4
    INTEGER, ALLOCATABLE :: vertex_blk(:,:,:)

    ! edges geometry (not necessarily initialized)

    ! longitude & latitude of edge midpoint
    ! index=1,nproma, index2=1,nblks_e
    TYPE(t_geographical_coordinates), ALLOCATABLE :: center(:,:)
    ! Cartesian normal to triangle edge
    ! index=1,nproma, index2=1,nblks_e
    TYPE(t_cartesian_coordinates),    ALLOCATABLE :: primal_cart_normal(:,:)
    ! normal to triangle edge
    ! index=1,nproma, index2=1,nblks_e
    TYPE(t_tangent_vectors),          ALLOCATABLE :: primal_normal(:,:)

    ! Edge length
    ! index1=nproma, index2=1,nblks_e
    REAL(wp),                         ALLOCATABLE :: primal_edge_length(:,:)

    ! Dual edge length
    ! index1=nproma, index2=1,nblks_e
    REAL(wp),                         ALLOCATABLE :: dual_edge_length(:,:)

    ! Edge to cell distance
    ! index1=nproma, index2=1,nblks_e, index3=1,2
    REAL(wp),                         ALLOCATABLE :: edge_cell_length(:,:,:)

    ! normal to triangle edge, projected to the location of the neighbor cells
    ! index=1,nproma, index2=1,nblks_e, index3=1,2
    TYPE(t_tangent_vectors),          ALLOCATABLE :: primal_normal_cell(:,:,:)
  END TYPE t_grid_edges


  ! grid_vertices class
  TYPE t_grid_vertices

    ! line indices of edges around a vertex:
    ! index1=1,nproma, index2=1,nblks_v, index3=1,6
    INTEGER, ALLOCATABLE :: edge_idx(:,:,:)
    ! block indices of edges around a vertex:
    ! index1=1,nproma, index2=1,nblks_v, index3=1,6
    INTEGER, ALLOCATABLE :: edge_blk(:,:,:)

    ! number of edges connected to vertex
    ! index1=1,nproma, index2=1,nblks_v
    INTEGER, ALLOCATABLE :: num_edges(:,:)

    ! longitude & latitude of vertex:
    ! index1=1,nproma, index2=1,nblks_v
    TYPE(t_geographical_coordinates), ALLOCATABLE :: vertex(:,:)

    ! line indices of neighbor vertices:
    ! index1=1,nproma, index2=1,nblks_v, index3=1,6
    INTEGER, ALLOCATABLE :: neighbor_idx(:,:,:)
    ! block indices of neighbor vertices:
    ! index1=1,nproma, index2=1,nblks_v, index3=1,6
    INTEGER, ALLOCATABLE :: neighbor_blk(:,:,:)

    ! line/block indices of cells around each vertex (index3=1,6):
    INTEGER, ALLOCATABLE :: cell_idx(:,:,:), cell_blk(:,:,:)

    ! area of hexagon/pentagon of which vertex is center:
    ! index1=1,nproma, index2=1,nblks_v
    REAL(wp), ALLOCATABLE :: dual_area(:,:)

    ! orientation (+/-1)
    ! index1=1,nproma, index2=1,nblks_v, index3=1,6
    integer(i2), ALLOCATABLE :: edge_orientation(:,:,:)
  END TYPE t_grid_vertices


  TYPE t_patch
    ! cell type (no. of cell edges)
    INTEGER :: cell_type

    ! total number of cells, edges and vertices
    INTEGER :: n_patch_cells
    INTEGER :: n_patch_edges
    INTEGER :: n_patch_verts

    ! ! values for the blocking
    !
    ! total...
    ! number of blocks and chunk length in last block
    ! ... for the cells
    INTEGER :: nblks_c
    INTEGER :: npromz_c
    ! ... for the edges
    INTEGER :: nblks_e
    INTEGER :: npromz_e
    ! ... for the vertices
    INTEGER :: nblks_v
    INTEGER :: npromz_v

    ! grid information on the patch
    TYPE(t_grid_cells)    ::  cells
    TYPE(t_grid_edges)    ::  edges
    TYPE(t_grid_vertices) ::  verts

    ! number of cells, edges and vertices in the global patch
    INTEGER :: n_patch_cells_g
    INTEGER :: n_patch_edges_g
    INTEGER :: n_patch_verts_g

    ! dummy (unused in this code)
    INTEGER :: n_childdom
  END TYPE t_patch


  !-------------------------------------
  ! Derived from ICON::mo_intp_data_strc
  !-------------------------------------
  TYPE t_int_state
    ! coefficient for interpolation from adjacent cells onto edge
    REAL(wp), ALLOCATABLE :: c_lin_e   (:,:,:)  ! (nproma,2,nblks_e)

    ! factor for divergence (nproma,cell_type,nblks_c)
    REAL(wp), ALLOCATABLE :: geofac_div(:,:,:)  ! (nproma,3,nblks_c)

    ! factor for curl (nproma,9-cell_type,nblks_v)
    REAL(wp), ALLOCATABLE :: geofac_rot(:,:,:)  ! (nproma,6,nblks_e)
  end type t_int_state

  !------------
  ! Parameters:
  !------------
  real(wp), parameter :: PI           = 3.1415926535897932384626433832795_wp
  real(wp), parameter :: TWOPI        = PI * 2
  real(wp), parameter :: PIHALF       = PI / 2

!------------------------------------------------------------------------------
CONTAINS
!------------------------------------------------------------------------------
  ELEMENTAL INTEGER FUNCTION blk_no(j)
    INTEGER, INTENT(in) :: j
    blk_no = MAX((ABS(j)-1)/nproma + 1, 1) ! i.e. also 1 for j=0, nproma=1
  END FUNCTION blk_no

  ELEMENTAL INTEGER FUNCTION idx_no(j)
    INTEGER, INTENT(in) :: j
    IF(j==0) THEN
      idx_no = 0
    ELSE
      idx_no = SIGN(MOD(ABS(j)-1,nproma)+1, j)
    ENDIF
  END FUNCTION idx_no

  ELEMENTAL INTEGER FUNCTION idx_1d(jl,jb)
    INTEGER, INTENT(in) :: jl, jb
    IF(jb<=0) THEN
      idx_1d = 0 ! This covers the special case nproma==1,jb=0,jl=1
      ! All other cases are invalid and get also a 0 returned
    ELSE
      idx_1d = SIGN((jb-1)*nproma + ABS(jl), jl)
    ENDIF
  END FUNCTION idx_1d

  SUBROUTINE get_indices_c (p_patch, i_blk, i_startblk, i_endblk, i_startidx, &
                            i_endidx, irl_start, opt_rl_end, opt_chdom)

    TYPE(t_patch), INTENT(IN) :: p_patch
    INTEGER, INTENT(IN) :: i_blk      ! Current block (variable jb in do loops)
    INTEGER, INTENT(IN) :: i_startblk ! Start block of do loop
    INTEGER, INTENT(IN) :: i_endblk   ! End block of do loop
    INTEGER, INTENT(IN) :: irl_start  ! refin_ctrl level where do loop starts
    INTEGER, INTENT(IN), OPTIONAL :: opt_rl_end ! refin_ctrl level where do loop ends
    INTEGER, INTENT(IN), OPTIONAL :: opt_chdom ! child domain position ID where
    ! negative refin_ctrl indices refer to
    INTEGER, INTENT(OUT) :: i_startidx, i_endidx ! Start and end indices (jc loop)

    i_startidx = 1
    i_endidx   = nproma
    IF (i_blk == i_endblk) i_endidx = p_patch%npromz_c
  END SUBROUTINE get_indices_c

  SUBROUTINE get_indices_e (p_patch, i_blk, i_startblk, i_endblk, i_startidx, &
                            i_endidx, irl_start, opt_rl_end, opt_chdom)

    TYPE(t_patch), INTENT(IN) :: p_patch
    INTEGER, INTENT(IN) :: i_blk      ! Current block (variable jb in do loops)
    INTEGER, INTENT(IN) :: i_startblk ! Start block of do loop
    INTEGER, INTENT(IN) :: i_endblk   ! End block of do loop
    INTEGER, INTENT(IN) :: irl_start  ! refin_ctrl level where do loop starts
    INTEGER, INTENT(IN) :: opt_rl_end ! refin_ctrl level where do loop ends
    INTEGER, INTENT(IN) :: opt_chdom ! child domain position ID where
    INTEGER, INTENT(OUT) :: i_startidx, i_endidx ! Start and end indices (je loop)
    i_startidx = 1
    i_endidx   = nproma
    IF (i_blk == i_endblk) i_endidx = p_patch%npromz_e
  END SUBROUTINE get_indices_e

  ! --------------------------------------------------------------------
  !> normalizes geographical coordinates (in radians)
  !  to the interval [-PI, PI]/[-PI/2,PI]
  !
  ELEMENTAL FUNCTION normalized_coord (p)
    TYPE (t_geographical_coordinates) :: normalized_coord
    TYPE (t_geographical_coordinates), INTENT(IN) :: p

    normalized_coord% lon = p% lon
    IF (p% lon < -PI) normalized_coord% lon = p% lon + TWOPI
    IF (p% lon >  PI) normalized_coord% lon = p% lon - TWOPI
    normalized_coord% lat = MIN (MAX (-PIHALF, p% lat), PIHALF)

! Old version of normalized_coord (for grid coordinates in degrees)
! normalizes geographical coordinates (in degrees)
! to the interval [-180., 180]/[-90.,90]
!   normalized_coord% lon = p% lon
!   IF (p% lon < -180._wp) normalized_coord% lon = p% lon + 360._wp
!   IF (p% lon >  180._wp) normalized_coord% lon = p% lon - 360._wp
!   normalized_coord% lat = MIN (MAX (-90._wp, p% lat), 90._wp)
  END FUNCTION normalized_coord

  !-------------------------------------------------------------------------
  !>
  !! Converts geographical to cartesian coordinates.
  !!
  !!
  !! @par Revision History
  !! Developed  by Luis Kornblueh  (2004).
  !!
  ELEMENTAL FUNCTION gc2cc (p_pos)  result(p_x)

    TYPE(t_geographical_coordinates), INTENT(in) :: p_pos     ! geo. coordinates

    TYPE(t_cartesian_coordinates)                :: p_x       ! Cart. coordinates

    REAL(wp)                                     :: z_cln, z_sln, z_clt, z_slt

    !-----------------------------------------------------------------------
    z_sln = SIN(p_pos%lon)
    z_cln = COS(p_pos%lon)
    z_slt = SIN(p_pos%lat)
    z_clt = COS(p_pos%lat)

    p_x%x(1) = z_cln*z_clt
    p_x%x(2) = z_sln*z_clt
    p_x%x(3) = z_slt

  END FUNCTION gc2cc
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Converts zonal @f$p\_gu@f$ and meridional vector components @f$p\_gv@f$ into cartesian.
  !!
  !! Converts zonal @f$p\_gu@f$ and meridional vector components @f$p\_gv@f$ into cartesian
  !! ones @f$(p\_cu, p\_cv, p\_cw)@f$ using the conversion
  !! @f{align*}{
  !! \begin{pmatrix} p\_cu \\ p\_cv \\ p\_cw \end{pmatrix} =
  !! \begin{pmatrix}
  !! -\sin p\_long & -\sin p\_lat \cdot \cos p\_long \\\
  !! \cos p\_long & -\sin p\_lat \cdot \sin p\_long \\\
  !! 0 & -\cos p\_lat
  !! \end{pmatrix} \cdot
  !! \begin{pmatrix} p\_gu \\ p\_gv \end{pmatrix}
  !! @f}
  !!
  !! @par Revision History
  !! Original version by Tobias Ruppert and Thomas Heinze, DWD (2006-11-14)
  !!
  ELEMENTAL SUBROUTINE gvec2cvec (p_gu, p_gv, p_long, p_lat, p_cu, p_cv, p_cw)
    !
    REAL(wp), INTENT(in)  :: p_gu, p_gv     ! zonal, meridional vec. component
    REAL(wp), INTENT(in)  :: p_long, p_lat  ! geo. coord. of data point
    REAL(wp), INTENT(out) :: p_cu, p_cv, p_cw            ! Cart. vector

    REAL(wp)              :: z_cln, z_sln, z_clt, z_slt  ! sin and cos of
    !-------------------------------------------------------------------------

    z_sln = SIN(p_long)
    z_cln = COS(p_long)
    z_slt = SIN(p_lat)
    z_clt = COS(p_lat)

    p_cu = z_sln * p_gu + z_slt * z_cln * p_gv
    p_cu = -1._wp * p_cu
    p_cv = z_cln * p_gu - z_slt * z_sln * p_gv
    p_cw = z_clt * p_gv

  END SUBROUTINE gvec2cvec
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Converts cartesian velocity vector @f$(p\_cu, p\_cv, p\_cw)@f$.
  !!
  !! Converts cartesian velocity vector @f$(p\_cu, p\_cv, p\_cw)@f$
  !! into zonal @f$p\_gu@f$ and meridional @f$g\_v@f$ velocity components
  !! using the conversion
  !! @f{align*}{
  !! \begin{pmatrix} p\_gu \\ p\_gv \end{pmatrix} =
  !! \begin{pmatrix}
  !! -\sin p\_long & \cos p\_long & 0 \\\
  !! -\sin p\_lat \cdot \cos p\_long & -\sin p\_lat \cdot \sin p\_long & \cos p\_lat
  !! \end{pmatrix} \cdot
  !! \begin{pmatrix} p\_cu \\ p\_cv \\ p\_cw \end{pmatrix}
  !! @f}
  !!
  !! @par Revision History
  !! Original version by Thomas Heinze, DWD (2006-11-16)
  !!
  ELEMENTAL SUBROUTINE cvec2gvec (p_cu, p_cv, p_cw, p_long, p_lat, p_gu, p_gv)
    !
    REAL(wp), INTENT(in)  :: p_cu, p_cv, p_cw  ! Cart. vector
    REAL(wp), INTENT(in)  :: p_long, p_lat     ! geo. coord. of data point
    REAL(wp), INTENT(out) :: p_gu, p_gv        ! zonal and meridional vec. comp.

    REAL(wp)              :: z_cln, z_clt, z_sln, z_slt  ! sin and cos of
    !-------------------------------------------------------------------------

    ! p_long and p_lat
    z_sln = SIN(p_long)
    z_cln = COS(p_long)
    z_slt = SIN(p_lat)
    z_clt = COS(p_lat)

    p_gu = z_cln * p_cv - z_sln * p_cu
    p_gv = z_cln * p_cu + z_sln * p_cv
    p_gv = z_slt * p_gv
    p_gv = z_clt * p_cw - p_gv

  END SUBROUTINE cvec2gvec
  !-------------------------------------------------------------------------

END MODULE mo_t_icon
