!
!+ Management of ICON grid structures, searching and interpolation.
!
MODULE mo_icon_grid
!
! Description:
!   Routines for allocation and reading ICON grid metadata,
!   searching within the ICON unstructured grid,
!   and calculation of interpolation coefficients.
!   The allocation and reading part is taken from
!   prep_icon(mo_remap_grid_icon).
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
!  Changes for GRIB2/GRIB_API/ICON
! V1_27        2013-11-08 Harald Anlauf
!  search_icon_global: warn for small negative weights due to rounding errors,
!                      but terminate when too large
! V1_42        2015-06-08 Harald Anlauf
!  ICON local patch
! V1_45        2015-12-15 Harald Anlauf
!  search_icon_global: allow for tolerance in vertex search
! V1_48        2016-10-06 Harald Anlauf
!  read/set up additional grid-related data (vertices/cells/edges)
! V1_50        2017-01-09 Andreas Rhodin
!  adapt COSMO-MEC to ICON-LAM:return the 6 surrounding grid-points
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
!------------------------------------------------------------------------------
  !-------------
  ! Modules used
  !-------------
  use mo_kind,         only: wp, i2
  use mo_exception,    ONLY: finish, message    ! abort on error condition
  USE mo_mpi_dace,     ONLY: dace,             &! DACE communicator info
              p_comm_work => d_comm,           &! Fixed communicator
                             p_bcast            ! Generic broadcast routine
  USE mo_t_icon,       ONLY: nproma,           &! Vector length for blocking
                             idx_no,           &! Block index from global index
                             blk_no,           &! Line  index within block
!                            gc2cc,            &! Geographical to cartesian coo.
                             gvec2cvec,        &! Zonal/meridional to cartesian
                             cvec2gvec,        &! Cartesian to zonal/meridional
                             t_patch,          &! Derived type for ICON patch
                             t_grid_cells,     &! Derived type for grid cells
                             t_grid_edges,     &! Derived type for grid edges
                             t_grid_vertices,  &! Derived type for grid vertices
                             normalized_coord, &! Normalize lon/lat
                             t_geographical_coordinates,&! Derived type (lon,lat)
                             t_cartesian_coordinates,   &! Derived type (x(3))
                             t_tangent_vectors           ! Derived type (v1,v2)
  use mo_dace_string,  only: decode_uuid        ! Decode uuid string
  use mo_algorithms,   only: index              ! Sort
  !---------------------
  ! netCDF f90 interface
  !---------------------
  use netcdf,      only: nf90_open,             &!
                         nf90_close,            &!
                         nf90_Inquire_Dimension,&!
                         nf90_inq_dimid,        &!
                         nf90_inq_varid,        &!
                         nf90_get_var,          &!
                         nf90_get_att,          &!
                         nf90_strerror,         &!
                         NF90_GLOBAL,           &! varid for global attributes
                         NF90_NOWRITE,          &! mode flag to open a dataset
                         NF90_NOERR              ! status return value: no error
  implicit none

  !----------------
  ! Public entities
  !----------------
  private
  public :: t_grid_icon         ! Derived type for ICON grid metadata
  public :: init_icongrid       ! Read metadata, set up interpolation
  public :: finish_icongrid     ! Clean up, release memory
  public :: search_icon_global  ! Determine interpolation points & weights
  public :: bench_search_icon   ! Simple benchmark code
  public :: compress_icon_grid  ! Compress redundant grid information
  public :: check_neighbours    ! check neighbour references for correct order
  public :: dist_to_bound       ! determine distance to boundary
  public :: set_search_grid     ! auxiliary grid for searching, interpolation
  !-------------------------------------
  ! Dummy publics for f95-only compilers
  !-------------------------------------
  public :: t_grid_reg
  public :: t_auxgrid
  public :: t_coord
  public :: t_uvec
  !============================================================================
  !------------
  ! Parameters:
  !------------
  real(wp), parameter :: PI           = 3.1415926535897932384626433832795_wp
  real(wp), parameter :: r2d          = 180/PI
  real(wp), parameter :: d2r          = PI/180
  real(wp), parameter :: earth_radius = 6.371229e6_wp   !! [m] average radius

! integer,  parameter :: p_comm_work  = dace% comm      ! Fixed communicator

  !============================================================================
  !------------------------------------
  ! Derived types and module variables:
  !------------------------------------

  !------------
  ! Unit-vector
  !------------
  type t_uvec
     real(wp) :: x(3)
  end type t_uvec

  !--------------------------
  ! (Local) Coordinate system
  !--------------------------
  type t_coord
     type(t_geographical_coordinates) :: coo
     real(wp) :: x(3)
     real(wp) :: u(3)
     real(wp) :: v(3)
  end type t_coord

  !-------------------------
  ! Simplified regular grids
  !-------------------------
  type t_grid_reg
    INTEGER           :: nx      = 0       ! number of grid points in x
    INTEGER           :: ny      = 0       ! number of grid points in y
    REAL(wp)          :: lon(2)  =   0._wp ! longitude bounds
    REAL(wp)          :: lat(2)  =   0._wp ! latitude  bounds
    REAL(wp)          :: lo1     =   0._wp ! longitude of first grid point
    REAL(wp)          :: la1     = -90._wp ! latitude  of first grid point
    REAL(wp)          :: di      =   0._wp ! longitudinal increment
    REAL(wp)          :: dj      =   0._wp ! latitudinal  increment
    REAL(wp)          :: dxi     =   0._wp ! longitudinal gridpoint density
    REAL(wp)          :: dyi     =   0._wp ! latitudinal  gridpoint density
    LOGICAL           :: global  = .true.  ! global grid
    LOGICAL           :: rot     = .false. ! rotated (currently not supported)
    LOGICAL           :: cyc_x   = .true.  ! cyclic boundary conditions in x
    LOGICAL           :: poly    = .true.  ! poles at y boundaries
  end type t_grid_reg

  !----------------------------
  ! Gridpoint of auxiliary grid
  !----------------------------
  type t_auxgrid
     real(wp) :: dlon    ! Longitude [degree]
     real(wp) :: dlat    ! Latitude  [degree]
     real(wp) :: dist    ! Distance
     real(wp) :: rvec(3)
     integer  :: idx(2)
  end type t_auxgrid

  !---------------------------------------------
  ! Derived type holding ICON grid information
  ! and auxiliary information for interpolation.
  !---------------------------------------------
  type t_grid_icon
     type(t_patch),    pointer :: patch => NULL() ! ICON grid structure
     type(t_grid_reg), pointer :: g               ! Regular mapping grid
     type(t_auxgrid),  pointer :: aux_grd  (:,:)  ! Auxiliary grid for search
     type(t_uvec),     pointer :: uvec_cell(:,:)  ! Cell   center coordinates
     type(t_uvec),     pointer :: uvec_vert(:,:)  ! Vertex center coordinates
     type(t_coord),    pointer :: vert     (:,:)  ! Vertices: local coordinates
     integer                   :: refcount   = 0  ! Reference counter (#links)
     integer                   :: grid_root  = -1 ! Primary subdivision
     integer                   :: grid_level = -1 ! Successive refinements
     integer                   :: grid_num   = -1 ! Number of grid used
     character(len=1)          :: uuid(16)        ! uuidOfHGrid
     logical                   :: global = .true. ! global grid?
  end type t_grid_icon

  !============================================================================
  integer                      :: dbg_level = 1
  CHARACTER(LEN=*), PARAMETER  :: modname   = 'mo_icon_grid'
  !============================================================================
CONTAINS
  !============================================================================
  subroutine init_icongrid (icongrid, gridfile, pio, verbose)
    type(t_grid_icon), intent(out) :: icongrid    ! ICON grid metadata
    character(len=*),  intent(in)  :: gridfile    ! NetCDF file: grid topology
    integer, optional, intent(in)  :: pio         ! Open file on I/O pe
    logical, optional, intent(in)  :: verbose     ! Enable debugging

    integer                :: lpio               ! Local copy of pio
    integer                :: ncid               ! NetCDF stream id
    integer                :: status             ! error status
    integer                :: ldbg               ! Debugging state
    integer                :: grid_root          ! Grid root parameter
    integer                :: grid_level         ! Grid subdivision level
    integer                :: grid_num           ! number_of_grid_used
    character(len=40)      :: uuid_str           ! uuidOfHGrid stored as string
    integer                :: euler              ! Euler characteristic
    logical                :: global             ! global grid?
    type(t_patch), pointer :: patch

    lpio = dace% pio; if (present (pio)) lpio = pio
    ncid = -1

    ldbg = dbg_level
    if (present (verbose)) then
       if (verbose) dbg_level = max (dbg_level, 5)
    end if

    if (dace% pe == lpio) then
       status = nf90_open (trim (gridfile), NF90_NOWRITE, ncid)
       if (status /= NF90_NOERR) write(0,*) "NetCDF error: ", nf90_strerror(status)
    end if
    call p_bcast (status, lpio)
    if (status /= NF90_NOERR) then
       call finish ("init_icongrid", "failure opening "// trim (gridfile))
    end if

    icongrid% uuid = achar (0)
    if (dace% pe == lpio) then
      status = nf90_get_att (ncid, NF90_GLOBAL, 'number_of_grid_used',grid_num)
      if (status /= NF90_NOERR) grid_num   = -1
      status = nf90_get_att (ncid, NF90_GLOBAL, 'grid_root',  grid_root)
      if (status /= NF90_NOERR) grid_root  = -1
      status = nf90_get_att (ncid, NF90_GLOBAL, 'grid_level', grid_level)
      if (status /= NF90_NOERR) grid_level = -1
      uuid_str = ""
      status = nf90_get_att (ncid, NF90_GLOBAL, 'uuidOfHGrid', uuid_str)
      if (status /= NF90_NOERR) uuid_str   = ""
      if (dbg_level > 0) then
        if (grid_num   /= -1) print *, "# number_of_grid_used = ", grid_num
        if (grid_root  /= -1) print *, "# grid_root           = ", grid_root
        if (grid_level /= -1) print *, "# grid_level          = ", grid_level
        if (uuid_str   /= "") print *, "# uuidOfHGrid         = ", &
                                       trim (uuid_str)
      end if
      icongrid% grid_num   = grid_num
      icongrid% grid_root  = grid_root
      icongrid% grid_level = grid_level
      call decode_uuid (trim (uuid_str), icongrid% uuid)
    end if
    call p_bcast (icongrid% grid_num   ,lpio)
    call p_bcast (icongrid% grid_root  ,lpio)
    call p_bcast (icongrid% grid_level ,lpio)
    call p_bcast (icongrid% uuid       ,lpio)

    allocate (icongrid% patch)
    call load_icon_grid (icongrid% patch, rank0=lpio, ncid=ncid)

    nullify  (icongrid% uvec_cell, icongrid% vert,  &
              icongrid% uvec_vert, icongrid% aux_grd)
    allocate (icongrid% g)

    if (dace% pe == lpio) then
       status = nf90_close (ncid)
       if (status /= NF90_NOERR) write(0,*) "NetCDF error: ", nf90_strerror(status)
    end if
    call p_bcast (status, lpio)
    if (status /= NF90_NOERR) then
       call finish ("init_icongrid", "failure closing "// trim (gridfile))
    end if

    patch => icongrid% patch
    euler =  patch% n_patch_verts_g - patch% n_patch_edges_g &
          +  patch% n_patch_cells_g
    select case (euler)
    case (2)
       global = .true.
    case (1)
       global = .false.
    case default
       call finish ("init_icongrid", "invalid grid file "// trim (gridfile))
    end select
    icongrid% global = global
    if (dace% pe == lpio) then
      if (dbg_level > 0) then
        print *, "# patch is global     = ", icongrid% global
      end if
    end if
    !--------------------------------------------------------------
    ! Set up coarser auxiliary grid for searching and interpolation
    !--------------------------------------------------------------
    call set_search_grid (icongrid, global)

    dbg_level = ldbg
  end subroutine init_icongrid
  !============================================================================
  subroutine finish_icongrid (icongrid)
    type(t_grid_icon), intent(inout) :: icongrid    ! ICON metadata

    if (associated (icongrid% aux_grd))   deallocate (icongrid% aux_grd)
    if (associated (icongrid% uvec_cell)) deallocate (icongrid% uvec_cell)
    if (associated (icongrid% uvec_vert)) deallocate (icongrid% uvec_vert)
    if (associated (icongrid% vert))      deallocate (icongrid% vert)
    if (associated (icongrid% g)) then
       call destruct_reg_grid (icongrid% g)
       deallocate  (icongrid% g)
    end if
    if (associated (icongrid% patch)) then
       call deallocate_icon_grid (icongrid% patch)
       deallocate  (icongrid% patch)
    end if
    nullify (icongrid% g, icongrid% aux_grd, icongrid% uvec_cell, &
             icongrid% uvec_vert, icongrid% vert, icongrid% patch)
  end subroutine finish_icongrid
  !============================================================================
  ! --------------------------------------------------------------------
  !> Load ICON grid to internal data structure
  !
  SUBROUTINE load_icon_grid (patch, rank0, ncid)
    TYPE (t_patch), INTENT(INOUT) :: patch
    INTEGER,        INTENT(IN)    :: rank0    !< MPI rank where data are read
    integer,        INTENT(IN)    :: ncid
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//'::load_icon_grid'
    INTEGER :: &
      &  i, j, idx, start_idx, end_idx, start_blk, end_blk, jb, jc,  &
      &  dimid, varID, n_patch_cells, n_patch_edges, n_patch_verts, ne
    REAL(wp), allocatable :: vlon(:), vlat(:), clon(:), clat(:),     &
                             area_of_c(:), area_of_v(:),             &
                             elon(:), elat(:), en_v1(:), en_v2(:),   &
                             len_e_p(:), len_e_d(:), dist_e_c(:,:)
    INTEGER,  allocatable :: refin_c_ctrl(:)
    INTEGER,  allocatable :: v_of_c(:,:), e_of_c(:,:), c_of_c(:,:),  &
      &                      v_of_e(:,:), c_of_e(:,:), c_of_v(:,:),  &
      &                      v_of_v(:,:), o_of_e(:,:), o_of_n(:,:),  &
      &                      e_of_v(:,:)
!   REAL(wp)              :: i_rr
    logical               :: dbg
    integer               :: euler              ! Euler characteristic
    logical               :: global             ! global grid?
    integer               :: status             ! error return parameter

    dbg = .false.
    IF (dace% pe == rank0) THEN
      IF (nf90_inq_dimid (ncid, 'cell',   dimid) /= NF90_NOERR) THEN
        CALL nf (nf90_inq_dimid (ncid, 'ncells', dimid), routine)
      END IF
      CALL nf (nf90_Inquire_Dimension (ncid, dimid, len=n_patch_cells), routine)
      IF (nf90_inq_dimid (ncid, 'edge',   dimid) /= NF90_NOERR) THEN
        CALL nf (nf90_inq_dimid(ncid, 'ncells2', dimid), routine)
      END IF
      CALL nf (nf90_Inquire_Dimension (ncid, dimid, len=n_patch_edges), routine)
      IF (nf90_inq_dimid (ncid, 'vertex', dimid) /= NF90_NOERR) THEN
        CALL nf (nf90_inq_dimid(ncid, 'ncells3', dimid), routine)
      END IF
      CALL nf (nf90_Inquire_Dimension (ncid, dimid, len=n_patch_verts), routine)
      IF (dbg_level >= 5) THEN
        dbg = .true.
        WRITE (0,*) "# n_patch_cells = ", n_patch_cells
        WRITE (0,*) "# n_patch_verts = ", n_patch_verts
        WRITE (0,*) "# n_patch_edges = ", n_patch_edges
      END IF
    END IF
    CALL p_bcast(n_patch_cells,rank0,p_comm_work)
    CALL p_bcast(n_patch_verts,rank0,p_comm_work)
    CALL p_bcast(n_patch_edges,rank0,p_comm_work)

    CALL allocate_icon_grid (patch, n_patch_cells, n_patch_edges, n_patch_verts)
    patch%n_patch_cells_g = n_patch_cells
    patch%n_patch_edges_g = n_patch_edges
    patch%n_patch_verts_g = n_patch_verts

    euler =  patch% n_patch_verts_g - patch% n_patch_edges_g &
          +  patch% n_patch_cells_g
    select case (euler)
    case (2)
       global = .true.
    case (1)
       global = .false.
    case default
       if (dace% pe == rank0) then
          WRITE (0,*) "# n_patch_cells = ", n_patch_cells
          WRITE (0,*) "# n_patch_edges = ", n_patch_edges
          WRITE (0,*) "# n_patch_verts = ", n_patch_verts
          WRITE (0,*) "# Euler number  = ", euler
       end if
       call finish ("load_icon_grid", "invalid grid file")
    end select

    IF (dbg) WRITE (0,*) "# load vertex coordinates"
    ALLOCATE(vlon(patch%n_patch_verts))
    ALLOCATE(vlat(patch%n_patch_verts))
    IF (dace% pe == rank0) THEN
      call chk ('vlon')
      CALL nf  (nf90_get_var (ncid, varid, vlon), routine)
      call chk ('vlat')
      CALL nf  (nf90_get_var (ncid, varid, vlat), routine)
    END IF
    CALL p_bcast(vlon,rank0,p_comm_work)
    CALL p_bcast(vlat,rank0,p_comm_work)
    DO i=1,patch%n_patch_verts
      jc = idx_no(i)
      jb = blk_no(i)
      patch%verts%vertex(jc,jb)%lon = vlon(i)   !*r2d
      patch%verts%vertex(jc,jb)%lat = vlat(i)   !*r2d
    END DO
    DEALLOCATE(vlon)
    DEALLOCATE(vlat)
    ! normalize coordinates
    start_blk = 1
    end_blk   = patch%nblks_v
!print *, "patch%npromz_v=", patch%npromz_v
    DO jb=start_blk,end_blk
      start_idx = 1
      end_idx   = nproma
      if (jb == end_blk) end_idx = patch%npromz_v
      DO jc=start_idx,end_idx
        patch%verts%vertex(jc,jb) = normalized_coord(patch%verts%vertex(jc,jb))
      END DO ! jc
    END DO !jb

    IF (dbg) WRITE (0,*) "# load cell-vertex indices"
    ALLOCATE(v_of_c(patch%n_patch_cells,patch%cell_type))
    IF (dace% pe == rank0) THEN
      call chk ('vertex_of_cell')
      CALL nf  (nf90_get_var (ncid, varid, v_of_c), routine)
    END IF
    CALL p_bcast(v_of_c,rank0,p_comm_work)
!NEC$ select_vector
!NEC$ ivdep
    DO i=1,patch%n_patch_cells
      jc = idx_no(i)
      jb = blk_no(i)
!NEC$ shortloop
      DO j=1,patch%cell_type
        idx = v_of_c(i,j)
        patch%cells%vertex_idx(jc,jb,j) = idx_no(idx)
        patch%cells%vertex_blk(jc,jb,j) = blk_no(idx)
      END DO
    END DO
    DEALLOCATE(v_of_c)

    IF (dbg) WRITE (0,*) "# load cell-edge indices"
    ALLOCATE(e_of_c(patch%n_patch_cells, patch%cell_type))
    IF (dace% pe == rank0) THEN
      call chk ('edge_of_cell')
      CALL nf  (nf90_get_var (ncid, varid, e_of_c), routine)
    END IF
    CALL p_bcast(e_of_c,rank0,p_comm_work)
!NEC$ select_vector
!NEC$ ivdep
    DO i=1,patch%n_patch_cells
      jc = idx_no(i)
      jb = blk_no(i)
!NEC$ shortloop
      DO j=1,patch%cell_type
        idx = e_of_c(i,j)
        patch%cells%edge_idx(jc,jb,j) = idx_no(idx)
        patch%cells%edge_blk(jc,jb,j) = blk_no(idx)
      END DO
    END DO
    DEALLOCATE(e_of_c)

    IF (dbg) WRITE (0,*) "# load vertex-edge indices"
    ALLOCATE(e_of_v(patch%n_patch_verts,6))
    IF (dace% pe == rank0) THEN
      call chk ('edges_of_vertex')
      CALL nf  (nf90_get_var (ncid, varid, e_of_v), routine)
    END IF
    CALL p_bcast(e_of_v,rank0,p_comm_work)
    patch%verts%edge_idx(:,:,6) = 0
!NEC$ ivdep
    DO i=1,patch%n_patch_verts
      jc = idx_no(i)
      jb = blk_no(i)
      ne = 0
!NEC$ unroll_completely
      DO j=1,6
        idx = e_of_v(i,j)
        ! take care of pentagon cells:
        if (idx > 0) then
          ne = ne + 1
          patch%verts%edge_idx(jc,jb,ne) = idx_no(idx)
          patch%verts%edge_blk(jc,jb,ne) = blk_no(idx)
        end if
      END DO
    END DO
    DEALLOCATE(e_of_v)

    IF (dbg) WRITE (0,*) "# load cell-neighbor indices"
    ALLOCATE(c_of_c(patch%n_patch_cells,patch%cell_type))
    IF (dace% pe == rank0) THEN
      call chk ('neighbor_cell_index')
      CALL nf  (nf90_get_var (ncid, varid, c_of_c), routine)
    END IF
    CALL p_bcast(c_of_c,rank0,p_comm_work)
!NEC$ select_vector
!NEC$ ivdep
    DO i=1,patch%n_patch_cells
      jc = idx_no(i)
      jb = blk_no(i)
!NEC$ shortloop
      DO j=1,patch%cell_type
        idx = c_of_c(i,j)
        patch%cells%neighbor_idx(jc,jb,j) = idx_no(idx)
        patch%cells%neighbor_blk(jc,jb,j) = blk_no(idx)
      END DO
    END DO
    DEALLOCATE(c_of_c)

    IF (dbg) WRITE (0,*) "# load vertex-neighbor indices"
    ALLOCATE(v_of_v(patch%n_patch_verts,6))
    IF (dace% pe == rank0) THEN
      call chk ('vertices_of_vertex')
      CALL nf  (nf90_get_var (ncid, varid, v_of_v), routine)
    END IF
    CALL p_bcast(v_of_v,rank0,p_comm_work)
    patch%verts%neighbor_idx(:,:,:) = -1
    patch%verts%neighbor_blk(:,:,:) = -1
#ifdef __NEC__
    patch%verts%num_edges   (:,:)   = -HUGE(0)
#endif
!NEC$ ivdep
    DO i=1,patch%n_patch_verts
      jc = idx_no(i)
      jb = blk_no(i)
      ne = 0
!NEC$ unroll_completely
      DO j=1,6
        idx = v_of_v(i,j)
        IF (idx > 0) THEN
           ne = ne + 1
           patch%verts%neighbor_idx(jc,jb,ne) = idx_no(idx)
           patch%verts%neighbor_blk(jc,jb,ne) = blk_no(idx)
        END IF
      END DO
      patch%verts%num_edges(jc,jb) = ne
#ifndef __NEC__   /* version for non-vector targets */
      IF (ne < 2 .and. dace% pe == rank0) THEN
         WRITE (0,*) "WARNING: vertex ",i,"has only", ne, "neighbors:"
         WRITE (0,*) v_of_v(i,:)
      END IF
#endif
    END DO
#ifdef __NEC__    /* version for vector targets */
    if (any (abs (patch%verts%num_edges) < 2)) then
       DO i=1,patch%n_patch_verts
          jc = idx_no(i)
          jb = blk_no(i)
          ne = patch%verts%num_edges(jc,jb)
          IF (ne < 2 .and. dace% pe == rank0) THEN
             WRITE (0,*) "WARNING: vertex ",i,"has only", ne, "neighbors:"
             WRITE (0,*) v_of_v(i,:)
          END IF
       END DO
    end if
#endif
    DEALLOCATE(v_of_v)

    IF (dbg) WRITE (0,*) "# load edge-vertex indices"
    ALLOCATE(v_of_e(patch%n_patch_edges,2))
    IF (dace% pe == rank0) THEN
      call chk ('edge_vertices')
      CALL nf  (nf90_get_var (ncid, varid, v_of_e), routine)
    END IF
    CALL p_bcast(v_of_e,rank0,p_comm_work )
!NEC$ ivdep
    DO i=1,patch%n_patch_edges
      jc = idx_no(i)
      jb = blk_no(i)
      DO j=1,2
        idx = v_of_e(i,j)
        patch%edges%vertex_idx(jc,jb,j) = idx_no(idx)
        patch%edges%vertex_blk(jc,jb,j) = blk_no(idx)
      END DO
    END DO
    DEALLOCATE(v_of_e)

    IF (dbg) WRITE (0,*) "# load edge-cell indices"
    ALLOCATE(c_of_e(patch%n_patch_edges,2))
    IF (dace% pe == rank0) THEN
      call chk ('adjacent_cell_of_edge')
      CALL nf  (nf90_get_var (ncid, varid, c_of_e), routine)
    END IF
    CALL p_bcast(c_of_e,rank0,p_comm_work )
!NEC$ ivdep
    DO i=1,patch%n_patch_edges
      jc = idx_no(i)
      jb = blk_no(i)
!NEC$ unroll(2)
      DO j=1,2
        idx = c_of_e(i,j)
        patch%edges%cell_idx(jc,jb,j) = idx_no(idx)
        patch%edges%cell_blk(jc,jb,j) = blk_no(idx)
      END DO
    END DO
    DEALLOCATE(c_of_e)

    IF (dbg) WRITE (0,*) "# load vertex-cell indices"
    ALLOCATE(c_of_v(patch%n_patch_verts,6))
    IF (dace% pe == rank0) THEN
      call chk ('cells_of_vertex')
      CALL nf  (nf90_get_var (ncid, varid, c_of_v), routine)
    END IF
    CALL p_bcast(c_of_v,rank0,p_comm_work)
    patch%verts%cell_idx(:,:,:) = -1
    patch%verts%cell_blk(:,:,:) = -1
!NEC$ ivdep
    DO i=1,patch%n_patch_verts
      jc = idx_no(i)
      jb = blk_no(i)
      ne = 0
!NEC$ unroll_completely
      DO j=1,6
        idx = c_of_v(i,j)
        ! take care of pentagon cells:
        IF (idx > 0) THEN
          ne = ne + 1
          patch%verts%cell_idx(jc,jb,ne) = idx_no(idx)
          patch%verts%cell_blk(jc,jb,ne) = blk_no(idx)
        END IF
      END DO
    END DO
    DEALLOCATE(c_of_v)

    IF (dbg) WRITE (0,*) "# load cell centers"
    ALLOCATE(clon(patch%n_patch_cells))
    ALLOCATE(clat(patch%n_patch_cells))
    IF (dace% pe == rank0) THEN
      call chk ('clon')
      CALL nf  (nf90_get_var (ncid, varid, clon), routine)
      call chk ('clat')
      CALL nf  (nf90_get_var (ncid, varid, clat), routine)
    END IF
    CALL p_bcast(clon,rank0,p_comm_work)
    CALL p_bcast(clat,rank0,p_comm_work)
    DO i=1,patch%n_patch_cells
      jc = idx_no(i)
      jb = blk_no(i)
      patch%cells%center(jc,jb)%lon = clon(i)   !*r2d
      patch%cells%center(jc,jb)%lat = clat(i)   !*r2d
    END DO
    DEALLOCATE(clon)
    DEALLOCATE(clat)
    ! normalize coordinates
    start_blk = 1
    end_blk   = patch%nblks_c
!print *, "patch%npromz_c=", patch%npromz_c
    DO jb=start_blk,end_blk
      start_idx = 1
      end_idx   = nproma
      if (jb == end_blk) end_idx = patch%npromz_c
      DO jc=start_idx,end_idx
        patch%cells%center(jc,jb) = normalized_coord(patch%cells%center(jc,jb))
      END DO ! jc
    END DO !jb

    IF (dbg) WRITE (0,*) "# load cell areas"
!   i_rr = 1 / earth_radius ** 2
    ALLOCATE(area_of_c(patch%n_patch_cells))
    IF (dace% pe == rank0) THEN
      call chk ('cell_area')
      CALL nf  (nf90_get_var (ncid, varid, area_of_c), routine)
    END IF
    CALL p_bcast(area_of_c,rank0,p_comm_work)
    DO i=1,patch%n_patch_cells
      jc = idx_no(i)
      jb = blk_no(i)
!     patch%cells%area(jc,jb) = area_of_c(i)*i_rr       ! [sr]
      patch%cells%area(jc,jb) = area_of_c(i)            ! [m^2]
    END DO
    DEALLOCATE(area_of_c)

    IF (dbg) WRITE (0,*) "# load edge-cell distances"
    ALLOCATE(dist_e_c(patch%n_patch_edges,2))
    IF (dace% pe == rank0) THEN
      call chk ('edge_cell_distance')
      CALL nf  (nf90_get_var (ncid, varid, dist_e_c), routine)
    END IF
    CALL p_bcast(dist_e_c,rank0,p_comm_work)
    DO i=1,patch%n_patch_edges
      jc = idx_no(i)
      jb = blk_no(i)
      patch%edges%edge_cell_length(jc,jb,1:2) = dist_e_c(i,1:2)
    END DO
    DEALLOCATE(dist_e_c)

    IF (dbg) WRITE (0,*) "# load primal and dual edge lengths"
    ALLOCATE(len_e_p(patch%n_patch_edges))
    ALLOCATE(len_e_d(patch%n_patch_edges))
    IF (dace% pe == rank0) THEN
      call chk ('edge_length')
      CALL nf  (nf90_get_var (ncid, varid, len_e_p), routine)
      call chk ('dual_edge_length')
      CALL nf  (nf90_get_var (ncid, varid, len_e_d), routine)
    END IF
    CALL p_bcast(len_e_p,rank0,p_comm_work)
    CALL p_bcast(len_e_d,rank0,p_comm_work)
!NEC$ ivdep
    DO i=1,patch%n_patch_edges
      jc = idx_no(i)
      jb = blk_no(i)
      patch%edges%primal_edge_length(jc,jb) = len_e_p(i)
      patch%edges%  dual_edge_length(jc,jb) = len_e_d(i)
    END DO
    DEALLOCATE(len_e_p)
    DEALLOCATE(len_e_d)

    IF (dbg) WRITE (0,*) "# load orientation of edges normals"
    ALLOCATE(o_of_n(patch%n_patch_cells,3))     ! Assume patch%cell_type==3
    IF (dace% pe == rank0) THEN
      call chk ('orientation_of_normal')
      CALL nf  (nf90_get_var (ncid, varid, o_of_n), routine)
    END IF
    CALL p_bcast(o_of_n,rank0,p_comm_work)
!NEC$ ivdep
    DO i=1,patch%n_patch_cells
      jc = idx_no(i)
      jb = blk_no(i)
      patch%cells%edge_orientation(jc,jb,1:3) = o_of_n(i,1:3)
    END DO
    DEALLOCATE(o_of_n)

    IF (dbg) WRITE (0,*) "# load edge orientation"
    ALLOCATE(o_of_e(patch%n_patch_verts,6))     ! Assume patch%cell_type==3
    IF (dace% pe == rank0) THEN
      call chk ('edge_orientation')
      CALL nf  (nf90_get_var (ncid, varid, o_of_e), routine)
    END IF
    CALL p_bcast(o_of_e,rank0,p_comm_work)
!NEC$ ivdep
    DO i=1,patch%n_patch_verts
      jc = idx_no(i)
      jb = blk_no(i)
      ! Take care of pentagon cells:
      ne = 0
!NEC$ unroll_completely
      DO j=1,6
        idx = o_of_e(i,j)
        IF (idx /= 0) THEN
          ne = ne + 1
          patch%verts%edge_orientation(jc,jb,ne) = idx
        END IF
      END DO
      if (ne == 5) patch%verts%edge_orientation(jc,jb,6) = 0
    END DO
    DEALLOCATE(o_of_e)

    IF (dbg) WRITE (0,*) "# load pentagon/hexagon areas"
    ALLOCATE(area_of_v(patch%n_patch_verts))
    IF (dace% pe == rank0) THEN
      call chk ('dual_area')
      CALL nf  (nf90_get_var (ncid, varid, area_of_v), routine)
    END IF
    CALL p_bcast(area_of_v,rank0,p_comm_work)
    DO i=1,patch%n_patch_verts
      jc = idx_no(i)
      jb = blk_no(i)
      patch%verts%dual_area(jc,jb) = area_of_v(i)       ! [m^2]
    END DO
    DEALLOCATE(area_of_v)

    !-------------------------------------------------------------
    ! Derive ancillary data required for computations on ICON grid
    !-------------------------------------------------------------
    CALL allocate_icon_grid_optionals (patch)

    IF (dbg) WRITE (0,*) "# load edges centers"
    ALLOCATE(elon(patch%n_patch_edges))
    ALLOCATE(elat(patch%n_patch_edges))
    IF (dace% pe == rank0) THEN
      call chk ('lon_edge_centre')
      CALL nf  (nf90_get_var (ncid, varid, elon), routine)
      call chk ('lat_edge_centre')
      CALL nf  (nf90_get_var (ncid, varid, elat), routine)
    END IF
    CALL p_bcast(elon,rank0,p_comm_work)
    CALL p_bcast(elat,rank0,p_comm_work)
    DO i=1,patch%n_patch_edges
      jc = idx_no(i)
      jb = blk_no(i)
      patch%edges%center(jc,jb)%lon = elon(i)   !*r2d
      patch%edges%center(jc,jb)%lat = elat(i)   !*r2d
    END DO
    DEALLOCATE(elon)
    DEALLOCATE(elat)

    IF (dbg) WRITE (0,*) "# load edges normals"
    ALLOCATE(en_v1(patch%n_patch_edges))
    ALLOCATE(en_v2(patch%n_patch_edges))
    IF (dace% pe == rank0) THEN
      call chk ('zonal_normal_primal_edge')
      CALL nf  (nf90_get_var (ncid, varid, en_v1), routine)
      call chk ('meridional_normal_primal_edge')
      CALL nf  (nf90_get_var (ncid, varid, en_v2), routine)
    END IF
    CALL p_bcast(en_v1,rank0,p_comm_work)
    CALL p_bcast(en_v2,rank0,p_comm_work)
    DO i=1,patch%n_patch_edges
      jc = idx_no(i)
      jb = blk_no(i)
      patch%edges%primal_normal(jc,jb)%v1 = en_v1(i)
      patch%edges%primal_normal(jc,jb)%v2 = en_v2(i)
    END DO
    DEALLOCATE(en_v1)
    DEALLOCATE(en_v2)

    if (.FALSE.)      then        ! not used, calculated by 'dist_to_bound'
      IF (dbg) WRITE (0,*) "# load refinement control flag for cells"
      IF (dace% pe == rank0) call chk ('refin_c_ctrl', status=status)
      CALL p_bcast (status, p_comm_work)
      IF (status == NF90_NOERR) then
        ALLOCATE (refin_c_ctrl (patch%n_patch_cells))
        IF (dace% pe == rank0) CALL nf (nf90_get_var (ncid, varid, refin_c_ctrl), routine)
        CALL p_bcast (refin_c_ctrl, rank0, p_comm_work)
        DO i=1,patch%n_patch_cells
          jc = idx_no(i)
          jb = blk_no(i)
          patch%cells%c_ctrl(jc,jb) = refin_c_ctrl(i)
        END DO
        DEALLOCATE(refin_c_ctrl)
      ELSE
        deallocate (patch%cells%c_ctrl)
      END IF
    else if (.not. global) then
      call dist_to_bound (patch, 50)
    else
      deallocate (patch%cells%c_ctrl)
    endif

    patch%edges%primal_normal_cell = t_tangent_vectors (HUGE(0._wp),HUGE(0._wp))
    ! Derivation of selected auxiliary data for edges
    IF (dbg) WRITE (0,*) "# derive auxiliary data for edges"
    if (global) then
       call complete_patchinfo (patch)
    else
       ! Not yet implemented (need refine_ctrl for grid boundaries)
    end if

    !-------------------------------------------------------------
    ! Deallocate temporaries which are usually not needed later on
    !-------------------------------------------------------------
    CALL deallocate_icon_grid_optionals (patch)

    !--------------------------------------------------------------------------
  contains
    !--------------------------------------------------------------------------
    subroutine chk (field, status)
      character(len=*)  ,intent(in)  :: field
      integer ,optional ,intent(out) :: status
      if (present (status)) then
        status = nf90_inq_varid (ncid, field, varid)
      else
        if (nf90_inq_varid (ncid, field, varid) /= NF90_NOERR) then
          CALL finish (routine, "Field <"//field//"> missing in grid file.")
        endif
      end if
    end subroutine chk
    !--------------------------------------------------------------------------
    SUBROUTINE nf (status, routine)
      INTEGER,          INTENT(in) :: status
      CHARACTER(len=*), INTENT(in) :: routine
      IF (status /= NF90_NOERR) THEN
         CALL finish (routine//' netCDF error', nf90_strerror(status))
      ENDIF
    END SUBROUTINE nf
  END SUBROUTINE load_icon_grid
  !============================================================================
  ! Allocate mandatory fields
  !--------------------------
  SUBROUTINE allocate_icon_grid (patch, n_cells, n_edges, n_verts)
    TYPE (t_patch), INTENT(INOUT) :: patch
    INTEGER,        INTENT(IN)    :: n_cells, n_edges, n_verts
    target                        :: patch
    ! local variables
!   CHARACTER(LEN=*), PARAMETER :: routine = modname//'::allocate_icon_grid'
    integer                     :: nproma_c, nproma_e, nproma_v
    type(t_patch),    pointer   :: p

    p => patch
    p%cell_type     = 3
    p%n_patch_cells = n_cells
    p%n_patch_edges = n_edges
    p%n_patch_verts = n_verts

    ! Handle case where nproma is ridiculously large (3D-Var)
    nproma_c   = min (nproma, n_cells)
    nproma_e   = min (nproma, n_edges)
    nproma_v   = min (nproma, n_verts)

    ! compute the no. of blocks:
    p%nblks_c  = blk_no(p%n_patch_cells)
    p%npromz_c = idx_no(p%n_patch_cells)
    p%nblks_e  = blk_no(p%n_patch_edges)
    p%npromz_e = idx_no(p%n_patch_edges)
    p%nblks_v  = blk_no(p%n_patch_verts)
    p%npromz_v = idx_no(p%n_patch_verts)

    ! create vertices and topology info:
    IF (dbg_level >=5 .and. dace% lpio) &
         WRITE (0,*) "# allocate data structures"
    ALLOCATE(p%verts%vertex            (nproma_v,p%nblks_v),             &! mand.
      &      p%verts%num_edges         (nproma_v,p%nblks_v),             &! mand.
      &      p%verts%neighbor_idx      (nproma_v,p%nblks_v, 6),          &! mand.
      &      p%verts%neighbor_blk      (nproma_v,p%nblks_v, 6),          &! mand.
      &      p%verts%cell_idx          (nproma_v,p%nblks_v, 6),          &! mand.
      &      p%verts%cell_blk          (nproma_v,p%nblks_v, 6),          &! mand.
      &      p%verts%edge_idx          (nproma_v,p%nblks_v, 6),          &
      &      p%verts%edge_blk          (nproma_v,p%nblks_v, 6),          &
      &      p%verts%dual_area         (nproma_v,p%nblks_v),             &
      &      p%verts%edge_orientation  (nproma_v,p%nblks_v, 6),          &
      &      p%cells%center            (nproma_c,p%nblks_c),             &! mand.
      &      p%cells%neighbor_idx      (nproma_c,p%nblks_c,p%cell_type), &! mand.
      &      p%cells%neighbor_blk      (nproma_c,p%nblks_c,p%cell_type), &! mand.
      &      p%cells%vertex_idx        (nproma_c,p%nblks_c,p%cell_type), &! needed
      &      p%cells%vertex_blk        (nproma_c,p%nblks_c,p%cell_type), &! needed
      &      p%cells%edge_idx          (nproma_c,p%nblks_c,p%cell_type), &
      &      p%cells%edge_blk          (nproma_c,p%nblks_c,p%cell_type), &
      &      p%cells%area              (nproma_c,p%nblks_c),             &
      &      p%cells%edge_orientation  (nproma_c,p%nblks_c,p%cell_type), &
!     &      p%edges%center            (nproma_e,p%nblks_e),             &! optional
      &      p%edges%cell_idx          (nproma_e,p%nblks_e, 2),          &
      &      p%edges%cell_blk          (nproma_e,p%nblks_e, 2),          &
      &      p%edges%vertex_idx        (nproma_e,p%nblks_e, 2),          &
      &      p%edges%vertex_blk        (nproma_e,p%nblks_e, 2),          &
!     &      p%cells%glb_index         (p%n_patch_cells),                &
!     &      p%edges%glb_index         (p%n_patch_edges),                &
      &      p%edges%  edge_cell_length(nproma_e,p%nblks_e, 2),          &
      &      p%edges%primal_edge_length(nproma_e,p%nblks_e),             &
      &      p%edges%  dual_edge_length(nproma_e,p%nblks_e),             &
!     &      p%edges%primal_cart_normal(nproma_e,p%nblks_e),             &! optional
!     &      p%edges%primal_normal     (nproma_e,p%nblks_e),             &! optional
      &      p%edges%primal_normal_cell(nproma_e,p%nblks_e, 2),          &
      &      p%cells%c_ctrl            (nproma_c,p%nblks_c)              )! may be missing

  END SUBROUTINE allocate_icon_grid
  !============================================================================
  ! Allocate optional fields
  !-------------------------
  SUBROUTINE allocate_icon_grid_optionals (patch)
    TYPE (t_patch), INTENT(INOUT) :: patch
    ! local variables
    integer                       :: nproma_e

    nproma_e  =  min (nproma, patch%n_patch_edges)

    allocate (patch%edges%center            (nproma_e,patch%nblks_e), &
         &    patch%edges%primal_cart_normal(nproma_e,patch%nblks_e), &
         &    patch%edges%primal_normal     (nproma_e,patch%nblks_e)  )
  end SUBROUTINE allocate_icon_grid_optionals
  !============================================================================
  ! Deallocate mandatory and optional fields
  !-----------------------------------------
  subroutine deallocate_icon_grid (patch)
    TYPE (t_patch), INTENT(INOUT) :: patch
    ! local variables
!   CHARACTER(LEN=*), PARAMETER :: routine = modname//'::deallocate_icon_grid'

    IF (dbg_level >=5 .and. dace% lpio) &
         WRITE (0,*) "# deallocate data structures"
    DEALLOCATE (patch%verts%vertex            , &! mandatory components
         &      patch%verts%num_edges         , &
         &      patch%verts%neighbor_idx      , &
         &      patch%verts%neighbor_blk      , &
         &      patch%verts%cell_idx          , &
         &      patch%verts%cell_blk          , &
         &      patch%cells%center            , &
         &      patch%cells%neighbor_idx      , &
         &      patch%cells%neighbor_blk        )
    if (allocated (patch%verts%edge_idx))       &! DACE stand-alone
    DEALLOCATE (patch%verts%edge_idx          , &
         &      patch%verts%edge_blk          , &
         &      patch%verts%dual_area         , &
         &      patch%verts%edge_orientation  , &
         &      patch%cells%vertex_idx        , &
         &      patch%cells%vertex_blk        , &
         &      patch%cells%edge_idx          , &
         &      patch%cells%edge_blk          , &
         &      patch%cells%area              , &
         &      patch%cells%edge_orientation  , &
!        &      patch%edges%center            , &! optional
         &      patch%edges%cell_idx          , &
         &      patch%edges%cell_blk          , &
         &      patch%edges%vertex_idx        , &
         &      patch%edges%vertex_blk        , &
         &      patch%edges%  edge_cell_length, &
         &      patch%edges%primal_edge_length, &
         &      patch%edges%  dual_edge_length, &
!        &      patch%edges%primal_cart_normal, &! optional
!        &      patch%edges%primal_normal     , &! optional
         &      patch%edges%primal_normal_cell  )
    if (allocated   (patch%cells%c_ctrl))       &! may be not present
         deallocate (patch%cells%c_ctrl)
    call deallocate_icon_grid_optionals (patch)
  end subroutine deallocate_icon_grid
  !============================================================================
  ! Deallocate optional fields
  !---------------------------
  subroutine deallocate_icon_grid_optionals (patch)
    TYPE (t_patch), INTENT(INOUT) :: patch

    if (allocated   (patch%edges%center)            ) &
         deallocate (patch%edges%center)
    if (allocated   (patch%edges%primal_cart_normal)) &
         deallocate (patch%edges%primal_cart_normal)
    if (allocated   (patch%edges%primal_normal)     ) &
         deallocate (patch%edges%primal_normal)
  end subroutine deallocate_icon_grid_optionals
  !============================================================================
  ! Compress redundant grid information using actual domain decomposition
  !----------------------------------------------------------------------
  subroutine compress_icon_grid (patch, lbc, ubc, lbv, ubv)
    type(t_patch), intent(inout) :: patch
    integer,       intent(in)    :: lbc(:)  ! Lower bounds for cells
    integer,       intent(in)    :: ubc(:)  ! Upper bounds for cells
    integer,       intent(in)    :: lbv(:)  ! Lower bounds for vertices
    integer,       intent(in)    :: ubv(:)  ! Upper bounds for vertices
    target                       :: patch

    type(t_grid_cells),    pointer :: cells
    type(t_grid_vertices), pointer :: verts
    integer                        :: lb2(2), lb3(3)
    integer                        :: ub2(2), ub3(3)
    real(wp),          allocatable :: tmp2d(:,:)
    integer(i2),       allocatable :: int3d(:,:,:)

    cells => patch% cells
    if (allocated (cells% area)) then
       lb2 = lbound (cells% area)
       ub2 = ubound (cells% area)
       if (any (lb2(1:2) > lbc(1:2)) .or. any (ub2(1:2) < ubc(1:2))) then
          write(0,*) dace% pe, "lb:", lb2(1:2), lbc(1:2)
          write(0,*) dace% pe, "ub:", ub2(1:2), ubc(1:2)
          call finish("compress_icon_grid", "bounds (cells %area)")
       end if
       if (any (lb2(1:2) < lbc(1:2)) .or. any (ub2(1:2) > ubc(1:2))) then
          allocate (tmp2d         (lbc(1):ubc(1),lbc(2):ubc(2)))
          tmp2d(:,:) = cells% area(lbc(1):ubc(1),lbc(2):ubc(2))
          call move_alloc (tmp2d, cells% area)
       end if
    end if

    if (allocated (cells% edge_orientation)) then
       lb3 = lbound (cells% edge_orientation)
       ub3 = ubound (cells% edge_orientation)
       if (any (lb3(1:2) > lbc(1:2)) .or. any (ub3(1:2) < ubc(1:2))) then
          write(0,*) dace% pe, "lb:", lb3(1:2), lbc(1:2)
          write(0,*) dace% pe, "ub:", ub3(1:2), ubc(1:2)
          call finish("compress_icon_grid", "bounds (cells %edge_orientation)")
       end if
       if (any (lb3(1:2) < lbc(1:2)) .or. any (ub3(1:2) > ubc(1:2))) then
          allocate (int3d(lbc(1):ubc(1),lbc(2):ubc(2),lb3(3):ub3(3)))
          int3d(:,:,:) = cells% edge_orientation(lbc(1):ubc(1),lbc(2):ubc(2),:)
          call move_alloc (int3d, cells% edge_orientation)
       end if
    end if

    verts => patch% verts
    if (allocated (verts% dual_area)) then
       lb2 = lbound (verts% dual_area)
       ub2 = ubound (verts% dual_area)
       if (any (lb2(1:2) < lbv(1:2)) .or. any (ub2(1:2) > ubv(1:2))) then
          allocate (tmp2d              (lbv(1):ubv(1),lbv(2):ubv(2)))
          tmp2d(:,:) = verts% dual_area(lbv(1):ubv(1),lbv(2):ubv(2))
          call move_alloc (tmp2d, verts% dual_area)
       end if
    end if

    if (allocated (verts% edge_orientation)) then
       lb3 = lbound (verts% edge_orientation)
       ub3 = ubound (verts% edge_orientation)
       if (any (lb3(1:2) < lbv(1:2)) .or. any (ub3(1:2) > ubv(1:2))) then
          allocate (int3d(lbv(1):ubv(1),lbv(2):ubv(2),lb3(3):ub3(3)))
          int3d(:,:,:) = verts% edge_orientation(lbv(1):ubv(1),lbv(2):ubv(2),:)
          call move_alloc (int3d, verts% edge_orientation)
       end if
    end if
  end subroutine compress_icon_grid
  !============================================================================
  !-------------------------------------------------------------------------
  !>
  !! Computes the local orientation of the edge primal normal and dual normal
  !! at the location of the cell centers and vertices.
  !!
  !! @par Revision History
  !!  developed by Guenther Zaengl, 2009-03-31
  !! Modification by Anurag Dipankar, MPIM, 2012-12-28
  !!
  SUBROUTINE complete_patchinfo (patch)
    TYPE(t_patch), INTENT(inout) :: patch
    !
    INTEGER  :: i, jb, je
    INTEGER  :: ilc1, ibc1, ilc2, ibc2
!   INTEGER  :: ilv1, ibv1, ilv2, ibv2  !,ilv3,ibv3,ilv4,ibv4
    REAL(wp) :: z_nu, z_nv, z_nx1(3)    !,z_nx2(3)
    REAL(wp) :: z_lon, z_lat, z_u, z_v  ! location and components of normal
    REAL(wp) :: z_norm                  ! norm of Cartesian normal
    TYPE(t_cartesian_coordinates) :: z_vec
#ifdef __NEC__
    integer  :: unresolved_3, unresolved_4
#endif

    CHARACTER(LEN=*), PARAMETER :: method_name = 'mo_intp_coeffs:complete_patchinfo'

    !-----------------------------------------------------------------------
    !
    ! First part derived from mo_grid_tools::calculate_patch_cartesian_positions
    !
!$omp parallel
!$omp do private(i,jb,je,z_lon,z_lat,z_u,z_v,z_norm,z_vec)
    DO i=1,patch% n_patch_edges
       je = idx_no(i)
       jb = blk_no(i)
       ! calculate edges positions

       ! location of edge midpoint
       z_lon = patch%edges%center(je,jb)%lon
       z_lat = patch%edges%center(je,jb)%lat

       ! zonal and meridional component of primal normal
       z_u = patch%edges%primal_normal(je,jb)%v1
       z_v = patch%edges%primal_normal(je,jb)%v2

       ! calculate Cartesian components of primal normal
       CALL gvec2cvec( z_u, z_v, z_lon, z_lat, z_vec%x(1), z_vec%x(2), z_vec%x(3) )

       ! compute unit normal to edge je
       z_norm = SQRT( DOT_PRODUCT(z_vec%x(1:3),z_vec%x(1:3)) )
       z_vec%x(1:3) = 1._wp / z_norm * z_vec%x(1:3)

       ! save the values in the according type structure of the patch
       patch%edges%primal_cart_normal(je,jb)%x(1:3) = z_vec%x(1:3)

    END DO
!$omp end do
!$omp end parallel

    !-----------------------------------------------------------------------
    !
    ! Second part derived from mo_intp_coeffs::complete_patchinfo
    !
#ifdef __NEC__
    unresolved_3 = 0
    unresolved_4 = 0
#endif
!$omp parallel
!$omp do private(i,jb,je,ilc1,ibc1,ilc2,ibc2, &
!$omp            z_nu,z_nv,z_lon,z_lat,z_nx1,z_norm)
!NEC$ ivdep
    DO i=1,patch% n_patch_edges
       je = idx_no(i)
       jb = blk_no(i)

       ! compute edge-vertex indices (and blocks) 3 and 4, which
       ! are the outer vertices of cells 1 and 2, respectively,
       ! and the inverse length bewtween vertices 3 and 4

       ilc1 = patch%edges%cell_idx(je,jb,1)
       ibc1 = patch%edges%cell_blk(je,jb,1)
       ilc2 = patch%edges%cell_idx(je,jb,2)
       ibc2 = patch%edges%cell_blk(je,jb,2)

!      ilv1 = patch%edges%vertex_idx(je,jb,1)
!      ibv1 = patch%edges%vertex_blk(je,jb,1)
!      ilv2 = patch%edges%vertex_idx(je,jb,2)
!      ibv2 = patch%edges%vertex_blk(je,jb,2)

       IF ((patch%cells%vertex_idx(ilc1,ibc1,1) /= &
            & patch%edges%vertex_idx(je,jb,1) .OR.  &
            & patch%cells%vertex_blk(ilc1,ibc1,1) /= &
            & patch%edges%vertex_blk(je,jb,1)) .AND.  &
            & (patch%cells%vertex_idx(ilc1,ibc1,1) /= &
            & patch%edges%vertex_idx(je,jb,2) .OR.  &
            & patch%cells%vertex_blk(ilc1,ibc1,1) /= &
            & patch%edges%vertex_blk(je,jb,2)) )        THEN

!!$          patch%edges%vertex_idx(je,jb,3) = patch%cells%vertex_idx(ilc1,ibc1,1)
!!$          patch%edges%vertex_blk(je,jb,3) = patch%cells%vertex_blk(ilc1,ibc1,1)

       ELSE IF ((patch%cells%vertex_idx(ilc1,ibc1,2) /= &
            & patch%edges%vertex_idx(je,jb,1) .OR.  &
            & patch%cells%vertex_blk(ilc1,ibc1,2) /= &
            & patch%edges%vertex_blk(je,jb,1)) .AND.  &
            & (patch%cells%vertex_idx(ilc1,ibc1,2) /= &
            & patch%edges%vertex_idx(je,jb,2) .OR.  &
            & patch%cells%vertex_blk(ilc1,ibc1,2) /= &
            & patch%edges%vertex_blk(je,jb,2)) )        THEN

!!$          patch%edges%vertex_idx(je,jb,3) = patch%cells%vertex_idx(ilc1,ibc1,2)
!!$          patch%edges%vertex_blk(je,jb,3) = patch%cells%vertex_blk(ilc1,ibc1,2)

       ELSE IF ((patch%cells%vertex_idx(ilc1,ibc1,3) /= &
            & patch%edges%vertex_idx(je,jb,1) .OR.  &
            & patch%cells%vertex_blk(ilc1,ibc1,3) /= &
            & patch%edges%vertex_blk(je,jb,1)) .AND.  &
            & (patch%cells%vertex_idx(ilc1,ibc1,3) /= &
            & patch%edges%vertex_idx(je,jb,2) .OR.  &
            & patch%cells%vertex_blk(ilc1,ibc1,3) /= &
            & patch%edges%vertex_blk(je,jb,2)) )        THEN

!!$          patch%edges%vertex_idx(je,jb,3) = patch%cells%vertex_idx(ilc1,ibc1,3)
!!$          patch%edges%vertex_blk(je,jb,3) = patch%cells%vertex_blk(ilc1,ibc1,3)

       ELSE
#ifdef __NEC__
!$omp atomic
          unresolved_3 = unresolved_3 + 1
#else
          CALL finish(method_name, "Unresolved edges%vertex(3)")
#endif
       ENDIF

       IF ((patch%cells%vertex_idx(ilc2,ibc2,1) /= &
            & patch%edges%vertex_idx(je,jb,1) .OR.  &
            & patch%cells%vertex_blk(ilc2,ibc2,1) /= &
            & patch%edges%vertex_blk(je,jb,1)) .AND.  &
            & (patch%cells%vertex_idx(ilc2,ibc2,1) /= &
            & patch%edges%vertex_idx(je,jb,2) .OR.  &
            & patch%cells%vertex_blk(ilc2,ibc2,1) /= &
            & patch%edges%vertex_blk(je,jb,2)) )        THEN

!!$          patch%edges%vertex_idx(je,jb,4) = patch%cells%vertex_idx(ilc2,ibc2,1)
!!$          patch%edges%vertex_blk(je,jb,4) = patch%cells%vertex_blk(ilc2,ibc2,1)

       ELSE IF ((patch%cells%vertex_idx(ilc2,ibc2,2) /= &
            & patch%edges%vertex_idx(je,jb,1) .OR.  &
            & patch%cells%vertex_blk(ilc2,ibc2,2) /= &
            & patch%edges%vertex_blk(je,jb,1)) .AND.  &
            & (patch%cells%vertex_idx(ilc2,ibc2,2) /= &
            & patch%edges%vertex_idx(je,jb,2) .OR.  &
            & patch%cells%vertex_blk(ilc2,ibc2,2) /= &
            & patch%edges%vertex_blk(je,jb,2)) )        THEN

!!$          patch%edges%vertex_idx(je,jb,4) = patch%cells%vertex_idx(ilc2,ibc2,2)
!!$          patch%edges%vertex_blk(je,jb,4) = patch%cells%vertex_blk(ilc2,ibc2,2)

       ELSE IF ((patch%cells%vertex_idx(ilc2,ibc2,3) /= &
            & patch%edges%vertex_idx(je,jb,1) .OR.  &
            & patch%cells%vertex_blk(ilc2,ibc2,3) /= &
            & patch%edges%vertex_blk(je,jb,1)) .AND.  &
            & (patch%cells%vertex_idx(ilc2,ibc2,3) /= &
            & patch%edges%vertex_idx(je,jb,2) .OR.  &
            & patch%cells%vertex_blk(ilc2,ibc2,3) /= &
            & patch%edges%vertex_blk(je,jb,2)) )        THEN

!!$          patch%edges%vertex_idx(je,jb,4) = patch%cells%vertex_idx(ilc2,ibc2,3)
!!$          patch%edges%vertex_blk(je,jb,4) = patch%cells%vertex_blk(ilc2,ibc2,3)

       ELSE
#ifdef __NEC__
!$omp atomic
          unresolved_4 = unresolved_4 + 1
#else
          CALL finish(method_name, "Unresolved edges%vertex(4)")
#endif
       ENDIF

!!$       ilv3 = patch%edges%vertex_idx(je,jb,3)
!!$       ibv3 = patch%edges%vertex_blk(je,jb,3)
!!$       ilv4 = patch%edges%vertex_idx(je,jb,4)
!!$       ibv4 = patch%edges%vertex_blk(je,jb,4)

       ! next step: compute projected orientation vectors for cells and vertices
       ! bordering to each edge (incl. vertices 3 and 4 introduced above)

       ! transform orientation vectors at local edge center to Cartesian space
       z_lon = patch%edges%center(je,jb)%lon
       z_lat = patch%edges%center(je,jb)%lat

       ! transform primal normal to cartesian vector z_nx1
       z_nx1(:)=patch%edges%primal_cart_normal(je,jb)%x(:)

!!$       ! transform dual normal to cartesian vector z_nx2
!!$       z_nx2(:)=patch%edges%dual_cart_normal(je,jb)%x(:)

       ! get location of cell 1

       z_lon = patch%cells%center(ilc1,ibc1)%lon
       z_lat = patch%cells%center(ilc1,ibc1)%lat

       ! compute local primal and dual normals at cell 1

       CALL cvec2gvec(z_nx1(1),z_nx1(2),z_nx1(3),z_lon,z_lat,z_nu,z_nv)
       z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

       patch%edges%primal_normal_cell(je,jb,1)%v1 = z_nu/z_norm
       patch%edges%primal_normal_cell(je,jb,1)%v2 = z_nv/z_norm

!!$       CALL cvec2gvec(z_nx2(1),z_nx2(2),z_nx2(3),z_lon,z_lat,z_nu,z_nv)
!!$       z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

!!$       patch%edges%dual_normal_cell(je,jb,1)%v1 = z_nu/z_norm
!!$       patch%edges%dual_normal_cell(je,jb,1)%v2 = z_nv/z_norm

       ! get location of cell 2

       z_lon = patch%cells%center(ilc2,ibc2)%lon
       z_lat = patch%cells%center(ilc2,ibc2)%lat

       ! compute local primal and dual normals at cell 2

       CALL cvec2gvec(z_nx1(1),z_nx1(2),z_nx1(3),z_lon,z_lat,z_nu,z_nv)
       z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

       patch%edges%primal_normal_cell(je,jb,2)%v1 = z_nu/z_norm
       patch%edges%primal_normal_cell(je,jb,2)%v2 = z_nv/z_norm

!!$       CALL cvec2gvec(z_nx2(1),z_nx2(2),z_nx2(3),z_lon,z_lat,z_nu,z_nv)
!!$       z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

!!$       patch%edges%dual_normal_cell(je,jb,2)%v1 = z_nu/z_norm
!!$       patch%edges%dual_normal_cell(je,jb,2)%v2 = z_nv/z_norm

       ! get location of vertex 1

!!$       z_lon = patch%verts%vertex(ilv1,ibv1)%lon
!!$       z_lat = patch%verts%vertex(ilv1,ibv1)%lat

       ! compute local primal and dual normals at vertex 1

!!$       CALL cvec2gvec(z_nx1(1),z_nx1(2),z_nx1(3),z_lon,z_lat,z_nu,z_nv)
!!$       z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)
!!$
!!$       patch%edges%primal_normal_vert(je,jb,1)%v1 = z_nu/z_norm
!!$       patch%edges%primal_normal_vert(je,jb,1)%v2 = z_nv/z_norm

!!$       CALL cvec2gvec(z_nx2(1),z_nx2(2),z_nx2(3),z_lon,z_lat,z_nu,z_nv)
!!$       z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

!!$       patch%edges%dual_normal_vert(je,jb,1)%v1 = z_nu/z_norm
!!$       patch%edges%dual_normal_vert(je,jb,1)%v2 = z_nv/z_norm

       ! get location of vertex 2

!!$       z_lon = patch%verts%vertex(ilv2,ibv2)%lon
!!$       z_lat = patch%verts%vertex(ilv2,ibv2)%lat

       ! compute local primal and dual normals at vertex 2

!!$       CALL cvec2gvec(z_nx1(1),z_nx1(2),z_nx1(3),z_lon,z_lat,z_nu,z_nv)
!!$       z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)
!!$
!!$       patch%edges%primal_normal_vert(je,jb,2)%v1 = z_nu/z_norm
!!$       patch%edges%primal_normal_vert(je,jb,2)%v2 = z_nv/z_norm

!!$       CALL cvec2gvec(z_nx2(1),z_nx2(2),z_nx2(3),z_lon,z_lat,z_nu,z_nv)
!!$       z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

!!$       patch%edges%dual_normal_vert(je,jb,2)%v1 = z_nu/z_norm
!!$       patch%edges%dual_normal_vert(je,jb,2)%v2 = z_nv/z_norm

       ! get location of vertex 3

!!$       z_lon = patch%verts%vertex(ilv3,ibv3)%lon
!!$       z_lat = patch%verts%vertex(ilv3,ibv3)%lat

       ! compute local primal and dual normals at vertex 3

!!$       CALL cvec2gvec(z_nx1(1),z_nx1(2),z_nx1(3),z_lon,z_lat,z_nu,z_nv)
!!$       z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)
!!$
!!$       patch%edges%primal_normal_vert(je,jb,3)%v1 = z_nu/z_norm
!!$       patch%edges%primal_normal_vert(je,jb,3)%v2 = z_nv/z_norm

!!$       CALL cvec2gvec(z_nx2(1),z_nx2(2),z_nx2(3),z_lon,z_lat,z_nu,z_nv)
!!$       z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

!!$       patch%edges%dual_normal_vert(je,jb,3)%v1 = z_nu/z_norm
!!$       patch%edges%dual_normal_vert(je,jb,3)%v2 = z_nv/z_norm

       ! get location of vertex 4

!!$       z_lon = patch%verts%vertex(ilv4,ibv4)%lon
!!$       z_lat = patch%verts%vertex(ilv4,ibv4)%lat

       ! compute local primal and dual normals at vertex 2

!!$       CALL cvec2gvec(z_nx1(1),z_nx1(2),z_nx1(3),z_lon,z_lat,z_nu,z_nv)
!!$       z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)
!!$
!!$       patch%edges%primal_normal_vert(je,jb,4)%v1 = z_nu/z_norm
!!$       patch%edges%primal_normal_vert(je,jb,4)%v2 = z_nv/z_norm

!!$       CALL cvec2gvec(z_nx2(1),z_nx2(2),z_nx2(3),z_lon,z_lat,z_nu,z_nv)
!!$       z_norm = SQRT(z_nu*z_nu+z_nv*z_nv)

!!$       patch%edges%dual_normal_vert(je,jb,4)%v1 = z_nu/z_norm
!!$       patch%edges%dual_normal_vert(je,jb,4)%v2 = z_nv/z_norm

    END DO
!$omp end do nowait

!$omp end parallel

#ifdef __NEC__
    if (unresolved_3 > 0) CALL finish(method_name, "Unresolved edges%vertex(3)")
    if (unresolved_4 > 0) CALL finish(method_name, "Unresolved edges%vertex(4)")
#endif

  END SUBROUTINE complete_patchinfo
  !----------------------------------------------------------------------------
  !============================================================================
  !------------------------------------------------------
  ! Set up auxiliary data for search within the ICON grid
  !------------------------------------------------------
  subroutine set_search_grid (icongrid, global)
    type(t_grid_icon), intent(inout) :: icongrid
    logical,           intent(in)    :: global    ! Global grid?

    integer                   :: ni                 ! grid resolution parameter
    integer                   :: nx, ny             ! auxiliary grid dimensions
    integer                   :: i, j, n, jb, jl
    integer                   :: nlen, nblks, npromz
    integer                   :: nproma_v, nproma_c
    integer                   :: ih, ix, iy, jx, jy, ixx
    integer                   :: nover, nset, nunset, nfar
    real(wp)                  :: dlat, dlon, glat, glon
    real(wp)                  :: dist
    real(wp)                  :: rvec(3), uvec(3)
    logical                   :: done
    type(t_patch),    pointer :: patch          ! global or top-level patch
    type(t_auxgrid),  pointer :: aux_grd(:,:)
    type(t_geographical_coordinates) :: coo
#ifdef __NEC__
    type(t_auxgrid), allocatable :: wrk_grd(:)
#endif
    integer, parameter        :: NI_MAX = 384   ! Limit size of auxiliary grid

    real(wp), parameter       :: INVAL = HUGE (0._wp)
    real                      :: t1, t2
!   logical                   :: debug = .true.
    logical                   :: debug = .false.

    patch => icongrid% patch
    !--------------------------
    ! Initialize auxiliary grid
    !--------------------------
    ni = nint (sqrt (patch% n_patch_cells_g / 20._wp))
    ni = min (ni, NI_MAX)
    nx = 3 * ni; if (mod(nx,2)==1) nx = nx + 1 ! even for checkerboard search
    ny = 2 * ni
    if (icongrid% g% nx == 0) then
       if (global) then
          call construct_reg_grid (icongrid% g, nx, ny)
       else
          call construct_reg_grid (icongrid% g, nx, ny, &
                                   lonw=minval (patch%verts%vertex%lon)*r2d, &
                                   lone=maxval (patch%verts%vertex%lon)*r2d, &
                                   lats=minval (patch%verts%vertex%lat)*r2d, &
                                   latn=maxval (patch%verts%vertex%lat)*r2d  )
       end if
    end if

    if (associated (icongrid% aux_grd)) deallocate (icongrid% aux_grd)
    allocate (icongrid% aux_grd(nx,ny))
    aux_grd => icongrid% aux_grd
    !===================================================================
    ! Initialize:
    !
    rvec(:) = INVAL
    do j = 1, ny
!NEC$ ivdep
       do i = 1, nx
          call gridpoint_coordinates (icongrid% g, i, j, dlon, dlat)
          aux_grd(i,j) = t_auxgrid (dlon, dlat, INVAL, rvec, (/ -1,-1 /) )
       end do
    end do
    !===================================================================
    if (debug) call cpu_time (t1)
    nproma_v = size (patch% verts% vertex, dim=1)
    npromz   =       patch% npromz_v
    nblks    =       patch% nblks_v
    allocate (icongrid% vert(nproma_v,nblks))
    do jb = 1, nblks
       nlen = nproma_v; if (jb == nblks) nlen = npromz
       do jl = 1, nlen
          !+++ sxf90 cannot vectorize this derived type assignment:
         !icongrid% vert(jl,jb)% coo      = patch% verts% vertex(jl,jb)
          !+++ but explicit assignment of components
          icongrid% vert(jl,jb)% coo% lon = patch% verts% vertex(jl,jb)% lon
          icongrid% vert(jl,jb)% coo% lat = patch% verts% vertex(jl,jb)% lat
          call set_xuv (icongrid% vert(jl,jb))
       end do
    end do
    if (debug) call cpu_time (t2)
    if (debug) print *, "Time for init of      vert [s]:", real (t2-t1)
    !===================================================================
    if (debug) call cpu_time (t1)
    nproma_c = size (patch% cells% center, dim=1)
    npromz   =       patch% npromz_c
    nblks    =       patch% nblks_c
    allocate (icongrid% uvec_cell(nproma_c,nblks))
    do jb = 1, nblks
       nlen = nproma_c; if (jb == nblks) nlen = npromz
       do jl = 1, nlen
          icongrid% uvec_cell(jl,jb)% x = &
               unit_vector (patch% cells% center(jl,jb))
       end do
    end do
    if (debug) call cpu_time (t2)
    if (debug) print *, "Time for init of uvec_cell [s]:", real (t2-t1)
    !===================================================================
    ! Phase 1: determine regular grid points which are closest to
    !          a vertex of the triangular grid.
    !
    ! Loop over vertices:
    call cpu_time (t1)
    nblks  = patch% nblks_v
    npromz = patch% npromz_v
    !
    nover  = 0
    nset   = 0
    nfar   = 0
#ifdef __NEC__
    allocate (wrk_grd(nproma_v))
#endif
    do jb = 1, nblks
       nlen = nproma_v; if (jb == nblks) nlen = npromz
#ifdef __NEC__  /* nfort may some day vectorize the following loop */
!NEC$ ivdep
       do jl = 1, nlen
          dlon = icongrid% vert(jl,jb)% coo% lon * r2d
          dlat = icongrid% vert(jl,jb)% coo% lat * r2d
          ! Map onto regular grid
          call nearest_gridpoint     (icongrid% g, dlon, dlat, i, j)
          ! Grid coordinates
          call gridpoint_coordinates (icongrid% g, i, j, glon, glat)
          dist = spdist (dlon, dlat, glon, glat)
          wrk_grd(jl) = t_auxgrid (dlon, dlat, dist, idx=[i,j], &
                                   rvec=icongrid% vert(jl,jb)% x)
       end do
       do jl = 1, nlen
          i    = wrk_grd(jl)% idx(1)
          j    = wrk_grd(jl)% idx(2)
          dist = wrk_grd(jl)% dist
          if (dist < aux_grd(i,j)% dist) then
             if (aux_grd(i,j)% dist == INVAL) then
                nset  = nset  + 1
             else
                nover = nover + 1
             end if
             aux_grd(i,j)         = wrk_grd(jl)
             aux_grd(i,j)% idx(1) = jl
             aux_grd(i,j)% idx(2) = jb
          else
             nfar = nfar + 1
          end if
       end do
#else
       do jl = 1, nlen
          dlon = icongrid% vert(jl,jb)% coo% lon * r2d
          dlat = icongrid% vert(jl,jb)% coo% lat * r2d
          ! Map onto regular grid
          call nearest_gridpoint     (icongrid% g, dlon, dlat, i, j)
          ! Grid coordinates
          call gridpoint_coordinates (icongrid% g, i, j, glon, glat)
          dist = spdist (dlon, dlat, glon, glat)
          if (dist < aux_grd(i,j)% dist) then
             if (aux_grd(i,j)% dist == INVAL) then
                nset  = nset  + 1
             else
                nover = nover + 1
             end if
             rvec(:) = icongrid% vert(jl,jb)% x
             aux_grd(i,j) = t_auxgrid (dlon, dlat, dist, rvec, (/ jl,jb /) )
          else
             nfar = nfar + 1
          end if
       end do
#endif
    end do

!!$    nunset = count (aux_grd(:,:)% dist == INVAL)
!!$    print *, "Total number of auxiliary grid points:", nx*ny
!!$    print *, "Number of unset auxiliary grid points:", nunset
!!$    print *, "Number of sets      :", nset
!!$    print *, "Number of overrides :", nover
!!$    print *, "Number of far points:", nfar
!!$    print *, "Totals (set)      :", nset, nx*ny-nunset
!!$    print *, "Totals (processed):", nset + nover + nfar, &
!!$         npromz + (nblks-1)*nproma
!!$
!!$    print *, "MV", maxval (count ((aux_grd(:,:)% dist == INVAL), dim=1))
!!$    print *, "ML", maxloc (count ((aux_grd(:,:)% dist == INVAL), dim=1))
!!$    print *, "mv", minval (count ((aux_grd(:,:)% dist == INVAL), dim=1))
!!$    print *, "ml", minloc (count ((aux_grd(:,:)% dist == INVAL), dim=1))
!!$    print *, count ((aux_grd(:,:30)% dist == INVAL), dim=1)

    if (debug) call cpu_time (t2)
    if (debug) print *, "Time for phase 1           [s]:", real (t2-t1)
    if (debug) call cpu_time (t1)

    ! Phase 2: recursively find vertices for remaining grid points
    !          where regular grid is apparently finer than triangular one.
    ih = 1
    n  = 1
    nset  = 0
    nover = 0
    do
!      print *, "n, ih =", n, ih
       !--------------------------------------------------------------
       ! "checkerboard" search (stride 2 in each direction) to prevent
       ! propagation of points in some preferred direction
       !--------------------------------------------------------------
       done = mod(n,4) == 1
       do jy = 1+mod(n,4)/2, ny,2
          if (all (aux_grd(:,jy)% dist /= INVAL)) cycle
          done = .false.
          do jx = 1+mod(n,2), nx, 2
             if (aux_grd(jx,jy)% dist /= INVAL) cycle
             ! Grid coordinates
             glon = aux_grd(jx,jy)% dlon
             glat = aux_grd(jx,jy)% dlat
             uvec = unitvector_d (glat, glon)
             ! Loop over halo of grid point (jx,jy) at "distance" ih
             do iy = jy-ih, jy+ih
                if (iy < 1 .or. iy > ny) cycle
                do ixx = jx-ih, jx+ih
                   if (abs (iy-jy)<ih .and. abs (ixx-jx)<ih) cycle
                   ix = ixx
                   if (global) then
                     if (ix <  1) ix = ix + nx
                     if (ix > nx) ix = ix - nx
                   else
                     if (ix <  1 .or. ix > nx)  cycle
                   endif
                   if (aux_grd(ix,iy)% dist == INVAL) cycle
                   ! Check this candidate
                   jl = aux_grd(ix,iy)% idx(1)
                   jb = aux_grd(ix,iy)% idx(2)
                   rvec = icongrid% vert(jl,jb)% x        !unit_vector (coo)
                   dist = gcdist (uvec, rvec)
                   if (dist < aux_grd(jx,jy)% dist) then
                      if (aux_grd(jx,jy)% dist == INVAL) then
                         nset  = nset  + 1
                      else
                         nover = nover + 1
                      end if
!                     coo  = patch   % verts% vertex(jl,jb)
                      coo  = icongrid% vert(jl,jb)% coo
                      dlat = coo% lat * r2d
                      dlon = coo% lon * r2d
                      if (dlon < 0) dlon = dlon + 360._wp
                      aux_grd(jx,jy) = &
                           t_auxgrid (dlon, dlat, dist, rvec, (/ jl,jb /) )
                   end if
                end do
             end do
          end do
       end do
       if (done) exit
       !ih = ih + 1
       n  = n  + 1
       if (n >= ny/2) then
          call finish ("set_search_grid","something went terribly wrong...")
       end if
    end do

!!$    print *, "mv", minval (count ((aux_grd(:,:)% dist == INVAL), dim=1))
!!$    print *, "ml", minloc (count ((aux_grd(:,:)% dist == INVAL), dim=1))
!!$    print *, "nset  =", nset
!!$    print *, "nover =", nover

    nunset = count (aux_grd(:,:)% dist == INVAL)
    if (nunset > 0) then
       write(0,*) "FATAL: n_unset =", nunset
       write(0,*) "maxloc (aux_grd(:,:)% dist)", maxloc (aux_grd(:,:)% dist)
       call finish ("set_search_grid","auxiliary grid not fully set up")
    end if
    if (debug) call cpu_time (t2)
    if (debug) print *, "Time for phase 2           [s]:", real (t2-t1)
  end subroutine set_search_grid
  !============================================================================
  pure function unit_vector (coo) result (x)
    real(wp)                                     :: x(3)
    type(t_geographical_coordinates), intent(in) :: coo ! Longitude, latitude
    !---------------------------------
    ! Return unit vector on the sphere
    !---------------------------------
    real(wp)                    :: rlon, rlat
    real(wp)                    :: sinlon, coslon, sinlat, coslat
    rlon = coo% lon !* d2r
    rlat = coo% lat !* d2r
    sinlon = sin (rlon)
    coslon = cos (rlon)
    sinlat = sin (rlat)
    coslat = cos (rlat)
    x(1)   = coslon * coslat    ! x
    x(2)   = sinlon * coslat    ! y
    x(3)   =          sinlat    ! z
  end function unit_vector
  !============================================================================
  pure function unitvector_d (dlat, dlon) result (x)
    real(wp), intent(in) :: dlat, dlon  ! Latitude, longitude [degree]
    real(wp)             :: x(3)
    !---------------------------------
    ! Return unit vector on the sphere
    !---------------------------------
    x = unit_vector (t_geographical_coordinates (dlon * d2r, dlat * d2r))
  end function unitvector_d
  !============================================================================
  !-------------------------------------
  ! Great circle distance on unit sphere
  !-------------------------------------
  pure function gcdist (uvec1, uvec2)
    real(wp)             :: gcdist
    real(wp), intent(in) :: uvec1(3), uvec2(3)
    real(wp) :: sprod
    sprod  = sum (uvec1(:) * uvec2 (:))
    gcdist = acos (sprod)
  end function gcdist
  !============================================================================
  elemental function spdist (lon1, lat1, lon2, lat2)
    real(wp)             :: spdist
    real(wp), intent(in) :: lon1, lat1  ! Latitude, longitude [degree]
    real(wp), intent(in) :: lon2, lat2  ! Latitude, longitude [degree]
    real(wp) :: sprod
    sprod  = sum (unitvector_d (lat1, lon1) * unitvector_d (lat2, lon2))
!   spdist = 1 - sprod
    spdist = acos (sprod)
  end function spdist
  !============================================================================
  elemental subroutine set_xuv (coord)
    type(t_coord), intent(inout)  :: coord
    !----------------------------------
    ! Set up unit vectors on the sphere
    !----------------------------------
    real(wp)                    :: rlon, rlat
    real(wp)                    :: sinlon, coslon, sinlat, coslat
    rlon = coord% coo% lon !* d2r
    rlat = coord% coo% lat !* d2r
    sinlon = sin (rlon)
    coslon = cos (rlon)
    sinlat = sin (rlat)
    coslat = cos (rlat)
    coord% x(1)   =  coslon * coslat ! x
    coord% x(2)   =  sinlon * coslat ! y
    coord% x(3)   =           sinlat ! z
    coord% u(1)   = -sinlon          ! u - vector
    coord% u(2)   =  coslon
    coord% u(3)   =  0._wp
    coord% v(1)   = -coslon * sinlat ! v - vector
    coord% v(2)   = -sinlon * sinlat
    coord% v(3)   =           coslat
  end subroutine set_xuv
  !============================================================================
  subroutine construct_reg_grid (g, nx, ny, lonw, lone, lats, latn)
    type(t_grid_reg),   intent(out) :: g
    integer,            intent(in)  :: nx, ny
    real(wp), optional, intent(in)  :: lonw, lone, lats, latn

    ! Global or regional grid, dedicated for triangle search:
    ! avoid poles and shift points relativ to regular grid by half grid spacing
    g% nx  = nx
    g% ny  = ny
    g% lon = (/  0._wp, 360._wp/) ! Longitude bounds
    g% lat = (/-90._wp,  90._wp/) ! Latitude  bounds
    if (present (lonw)) g% lon(1) = lonw
    if (present (lone)) g% lon(2) = lone
    if (present (lats)) g% lat(1) = lats
    if (present (latn)) g% lat(2) = latn
    if (abs (g% lon(2) - g% lon(1) - 360._wp) < 1.e-6_wp) then
       g% cyc_x  = .true.
    else
       g% cyc_x  = .false.
    end if
    if (g% lat(1) > -90._wp .or. g% lat(1) < 90._wp .or. .not. g% cyc_x) then
       g% global = .false.
    else
       g% global = .true.
    end if
    g% poly   = .false.
    g% di  = (g% lon(2) - g% lon(1)) / nx   ! Grid spacing
    g% dj  = (g% lat(2) - g% lat(1)) / ny
    g% lo1 =  g% lon(1) + g% di / 2         ! Reference longitude
    g% la1 =  g% lat(1) + g% dj / 2         ! Reference latitude
!   g% lo1 =  g% lon(1)                     ! Reference longitude
!   g% la1 =  g% lat(1)                     ! Reference latitude
    g% dxi = 1._wp / g% di                  ! Gridpoint density
    g% dyi = 1._wp / g% dj
!   print *, "construct_reg_grid:", nx, ny, &
!        g% lo1, g% la1 !, g% di, g% dj, g% dxi, g% dyi
  end subroutine construct_reg_grid
  !============================================================================
  subroutine destruct_reg_grid (g)
    type(t_grid_reg), intent(inout) :: g
    g% nx = 0
    g% ny = 0
  end subroutine destruct_reg_grid
  !============================================================================
  pure subroutine gridpoint_coordinates (g, i, j, dlon, dlat)
    type(t_grid_reg), intent(in)  :: g
    integer,          intent(in)  :: i, j
    real(wp),         intent(out) :: dlon, dlat
    dlon = g% lo1 + (i-1) * g% di
    dlat = g% la1 + (j-1) * g% dj
  end subroutine gridpoint_coordinates
  !============================================================================
  pure subroutine nearest_gridpoint (g, dlon, dlat, i, j)
    type(t_grid_reg), intent(in)  :: g
    real(wp),         intent(in)  :: dlon, dlat
    integer,          intent(out) :: i, j

    real(wp)  :: lon, lat
    lat = dlat
    lon = dlon
    if (lon > g% lon(2)) lon = lon - 360._wp
    if (lon < g% lon(1)) lon = lon + 360._wp
    i = nint ((lon - g% lo1) * g% dxi) + 1
    j = nint ((lat - g% la1) * g% dyi) + 1
    i = min (max (i, 1), g% nx)
    j = min (max (j, 1), g% ny)
    if (.not. g% global) then
       if (lon > g% lon(2) .or. lat < g% lat(1) .or. lat > g% lat(2)) then
          i = -1
          j = -1
       end if
    end if
  end subroutine nearest_gridpoint
  !============================================================================
  subroutine search_icon_global (icongrid, lon, lat, idx, blk, w, ngp, mode, &
                                 l_triangle_in)
  !----------------------------------------------------------------------
  ! mode = 3 (default) : return 3 neighbour grid-points with weights
  ! mode = 6           : return 6 surrounding grid-points with weights
  !----------------------------------------------------------------------
    type(t_grid_icon),  intent(in)  :: icongrid  ! ICON grid metadata
    real(wp),           intent(in)  :: lon(:)    ! Longitudes [degree]
    real(wp),           intent(in)  :: lat(:)    ! Latitudes  [degree]
    integer,            intent(out) :: idx(:,:)  ! Line  indices (mode,N)
    integer,            intent(out) :: blk(:,:)  ! Block indices (mode,N)
    real(wp),           intent(out) :: w  (:,:)  ! Weights       (mode,N)
    integer,            intent(out) :: ngp(:)    ! number of points returned
    integer,  optional, intent(in)  :: mode      ! mode (3 is default)
    logical,  optional, intent(in)  :: l_triangle_in ! use invariant triangles
                                                     ! (default: .false.)
    !----------------
    ! Local variables
    !----------------
    logical            :: debug
    integer            :: i, j, jj, n         ! Indices and counter
    integer            :: ip, np, ne          ! Number of points, edges
    integer            :: jl, jb              ! Vertex line & block index
    integer            :: kl, kb              ! Cell   line & block index
    real(wp)           :: dlon, dlat          ! Longitude, latitude [degree]
    real(wp)           :: wmin                ! Minimum weight
    real(wp)           :: point(3), uvec(3)   ! Unit vectors
    integer            :: nr                  ! local copy of 'mode'

    type(t_grid_cells),      pointer :: cells ! Cells of global patch
    type(t_grid_vertices),   pointer :: verts ! Vertices of global patch
    type(t_geographical_coordinates) :: coo   ! Coordinates (lon,lat)

    integer, parameter :: nv=6                ! no.vertices adjacent to vertex
    integer            :: jlc(nv), jbc(nv)    ! Cell   line & block indices
    integer            :: jlv(nv), jbv(nv)    ! Vertex line & block indices
    real(wp)           :: x1(3,nv), x2(3,nv)  ! Temporaries for ...
    real(wp)           :: x3(3,nv), z(3,nv)   !    weight calculations
    real(wp)           :: SW(3,nv)            ! Symplectic weights
    real(wp)           :: nvec(3,0:nv)        ! Unit vectors of vertices
    real(wp)           :: pdist(nv+1)         ! Distance point to vertices
    !-------------------------------------------------
    ! Triangle constellations for interpolation points
    !-------------------------------------------------
    integer, parameter :: ns=3                ! no.triangles considered
    integer, parameter :: djp(3) = (/  1,  2, -2 /)
    integer, parameter :: djm(3) = (/  2, -2, -1 /)
    integer, parameter,dimension(12) :: triangles = (/ 1, 3, 5, &
         1, 2, 3,&
         3, 4, 5,&
         5, 6, 1 /)
    integer            :: jp(ns), jm(ns)

#ifdef CHECK_CONSISTENCY
    integer            :: jle(nv), jbe(nv)    ! Edge line & block indices
    integer            :: jlc1(2), jbc1(2)    ! for connectivity/
    integer            :: jlc2(2), jbc2(2)    !     consistency checks
    integer            :: j1
    logical            :: ok1, ok2
    type(t_grid_edges), pointer :: edges      ! Edges of global patch
#endif
    !------------------------------------------------------
    ! Threshold for negative weights due to rounding errors
    !------------------------------------------------------
    real(wp),parameter :: maxnwgt = 1.e-9_wp  ! Max. negative weight
    real(wp),parameter :: tol = EPSILON(1._wp)! Vertex search tolerance
    integer            :: inext,ntriangles    ! ntriangles=to be tested nr.
                                              ! of triangles defined by
                                              ! cell centers of vertex
    logical            :: l_triangle

    cells => icongrid% patch% cells
    verts => icongrid% patch% verts
#ifdef CHECK_CONSISTENCY
    edges => icongrid% patch% edges
#endif

    nr         = 3; if (present(mode)) nr = mode
    if (.not.present(l_triangle_in)) then
      l_triangle = .false.
    else
      l_triangle = l_triangle_in
    end if
    if (size (idx,dim=1) < nr) then
       write(0,*) "mode, size(idx,1) =", nr, size (idx,dim=1)
       call finish ("search_icon_global","idx too small for mode!")
    end if
    idx(1:nr,:) = 0
    blk(1:nr,:) = 0
    w  (1:nr,:) = 0._wp
    ngp(     :) = 0
    np = size (lon)
    debug = .false.
!   debug = (np == 1)
l1: do ip = 1, np
       dlon = lon(ip)
       dlat = lat(ip)
       wmin = 0._wp
       !=======================================================================
       ! Represent point by unit vector in 3D space
       if (debug) then
         print *
         print *, "=============================================="
         print *, "New search at lon/lat =",ip,dlon,dlat
         print *
       endif
       coo   = t_geographical_coordinates (dlon * d2r, dlat * d2r)
       point = unit_vector (coo)
       do i = 1, nv
          z(:,i) = point
       end do
       !=======================================================================
       call nearest_gridpoint (icongrid% g, dlon, dlat, i, j)
       if (i < 0 .or. j < 0) then
          !------------------------
          ! Gridpoint out of domain
          !------------------------
          cycle l1
       end if
       kl = icongrid% aux_grd(i,j)% idx(1)
       kb = icongrid% aux_grd(i,j)% idx(2)

       ! Find closest vertex:
       n = 0
       do
          n         = n + 1
          ne        = verts% num_edges   (kl,kb)
          jlv(1:ne) = verts% neighbor_idx(kl,kb,1:ne)
          jbv(1:ne) = verts% neighbor_blk(kl,kb,1:ne)
          if (debug) then
             coo    = icongrid% vert(kl,kb)% coo
             print *, "Vertex at center:", coo% lon * r2d, coo% lat * r2d," , iteration", n
             print *, "Number of edges connected to vertex:", ne
             print *
             print *, "Neighbor vertices:"
          end if
          !--------------------------------------------------------------
          ! Compare scalar products of point with the current vertex and
          ! with its neighbors to locate the vertex closest to the point:
          !--------------------------------------------------------------
!         coo          = verts% vertex(kl,kb)
!         uvec         = unit_vector (coo)
          uvec         = icongrid% vert(kl,kb)% x
          pdist(ne+1)  = sum (point * uvec)
          do i = 1, ne
!            coo       = verts% vertex(jlv(i),jbv(i))
!            nvec(:,i) = unit_vector (coo)
             nvec(:,i) = icongrid% vert(jlv(i),jbv(i))% x
             pdist(i)  = sum (point * nvec(1:3,i))
             if (debug) then
                if (i == 1) print *, "point =", point
                print *, "nvec   ", i, ":", nvec(:,i)
                coo    = icongrid% vert(jlv(i),jbv(i))% coo
                print *, "Vertex ", i, ":", coo% lon*r2d, coo% lat*r2d, pdist(i), pdist(i)-pdist(ne+1)
             end if
          end do
          j = maxloc (pdist(:ne+1), dim=1)
          if (debug) print *, "pdist(ref)   =", pdist(ne+1)
          if (debug) print *, "maxloc (p*r) =", j

!         if (j > ne .or. pdist(j) == pdist(ne+1)) exit
          if (j > ne .or. abs (pdist(j) - pdist(ne+1)) <= tol) exit

          ! Continue with vertex closest to point
          kl = jlv(j)
          kb = jbv(j)
          if (n > 100) then
             write(6,*) dace% pe,"search_icon_global: lat, lon =", dlat, dlon
             write(6,*) dace% pe,"search_icon_global: pdist    =", pdist(:ne+1)
             write(6,*) dace% pe,"search_icon_global: ne,j,dist=", ne, j,     &
                        abs (pdist(j)-pdist(ne+1))
             !---------------------------------------------------------------
             ! Give up in pathological cases at grid boundary, else terminate
             !---------------------------------------------------------------
             if (ne < 5) exit
             call finish ("search_icon_global","search did not converge!")
          end if
          if (debug) print *, "Restarting from vertex", j
          ! cycle
       end do

       !---------------------------------------------------------------
       ! Vertex too close to boundary, has too few neighboring vertices
       !---------------------------------------------------------------
       if (ne < 5) cycle l1

       ! point lies in one of the (ne) triangles with a common vertex.
       ! Set up vertices of (ne) triangles adjacent to the central vertex:
       if (.not.l_triangle) then
       nvec(:,0) = nvec(:,ne)

!NEC$ shortloop
!NEC$ ivdep
       do i = 1, ne
!         z (:,i) = point
          x1(:,i) = uvec
          x2(:,i) = nvec(:,i-1)
          x3(:,i) = nvec(:,i)
       end do
       call Symplectic_Weights (z (:,:ne), x1(:,:ne), x2(:,:ne), x3(:,:ne), &
                                SW(:,:ne))
       j = maxloc (minval (SW(:,:ne),dim=1), dim=1)

       if (debug) then
          do i = 1, ne
             print *, "z (",i,"):",z (:,i)
             print *, "x1(",i,"):",x1(:,i)
             print *, "x2(",i,"):",x2(:,i)
             print *, "x3(",i,"):",x3(:,i)
          end do
          do i = 1, ne
             print *, "Weights", i, ":", SW(:,i)
          end do
          print *, "maxloc (minval (SW))=", j
       end if

       ! Save center of cell containing point

       jl = verts% cell_idx(kl,kb,j)
       jb = verts% cell_blk(kl,kb,j)
       if (debug) print *, "jl, jb =", jl, jb

#ifdef CHECK_CONSISTENCY
       ! Consistency checks of connectivity between vertices and edges
       j1 = j-1;  if (j1 == 0) j1 = ne
       jle(1:ne) = verts% edge_idx(kl,kb,1:ne)
       jbe(1:ne) = verts% edge_blk(kl,kb,1:ne)

       jlc1(1:2) = edges% cell_idx(jle(j1),jbe(j1),1:2)
       jbc1(1:2) = edges% cell_blk(jle(j1),jbe(j1),1:2)

       jlc2(1:2) = edges% cell_idx(jle(j), jbe(j), 1:2)
       jbc2(1:2) = edges% cell_blk(jle(j), jbe(j), 1:2)

       ok1 = .false.
       ok2 = .false.
       do i = 1, 2
          if (jl == jlc1(i) .and. jb == jbc1(i)) ok1 = .true.
          if (jl == jlc2(i) .and. jb == jbc2(i)) ok2 = .true.
       end do

       if (.not. (ok1 .and. ok2)) then
          print *, "Internal error during search:"
          print *, "jl, jb =", jl, jb
          do i = 1, 2
             print *, "i, jlc1, jbc1 =", i, jlc1(i), jbc1(i)
             print *, "i, jlc2, jbc2 =", i, jlc2(i), jbc2(i)
          end do
          call finish ("search_icon_global","Internal error during search")
       end if
#endif

       if (debug) then
          coo = cells% center(jl,jb)
          print *, "Center of Final Cell:", coo% lon * r2d, coo% lat * r2d
          print *, "------------------------------------------"
       end if
       end if
       !-----------------------------------------------------------
       ! Final step: search cell centers required for interpolation
       !-----------------------------------------------------------
       jlc(1:ne) = verts% cell_idx(kl,kb,1:ne)
       jbc(1:ne) = verts% cell_blk(kl,kb,1:ne)
       !---------------------------------------------------------
       ! Vertex too close to boundary, too few neighboring cells?
       !---------------------------------------------------------
       if (jlc(ne) < 0) cycle l1

!NEC$ shortloop
       do i = 1, ne
!         coo       = cells% center(jlc(i),jbc(i))
!         nvec(:,i) = unit_vector(coo)
          ! Use pre-calculated unit vectors
          nvec(:,i) = icongrid% uvec_cell(jlc(i),jbc(i))% x
       end do

       select case (nr)
       case default
         call finish ("search_icon_global","Invalid mode")
       case (6)
       !--------------------------------
       ! Return 6 surrounding gridpoints
       !--------------------------------
       if (l_triangle) then

!NEC$ shortloop
         do i = 1, ne
           z (:,1) = point
           x1(:,i) = uvec
           x2(:,i) = nvec(:,mod(i+ne-2,ne)+1)
           x3(:,i) = nvec(:,mod(i+ne-1,ne)+1)
         end do
         call Symplectic_Weights (z (:,:ne), x1(:,:ne), &
           x2(:,:ne), x3(:,:ne), SW(:,:ne))
         jj=1
!NEC$ shortloop
         do j=2,ne
           if ( minval(SW(:,j)) > minval(SW(:,jj)) ) then
              jj=j
           end if
         end do
         wmin = minval (SW(:,jj))

         if (SW(2,jj) > SW(3,jj)) then
            j=mod(jj+ne-2,ne)+1
            inext=jj
         else
            j=jj
            inext=mod(jj+ne-2,ne)+1
         end if

!NEC$ shortloop
           do i = 1, ne
             idx(i,ip) = jlc (mod (ne - 2 + i + j, ne) + 1 )
             blk(i,ip) = jbc (mod (ne - 2 + i + j, ne) + 1 )
             if ( mod(i+j+ne-2,ne)+1 == j) then
               w  (i,ip) = SW(1,jj)/real(ne,wp)+max(SW(3,jj),SW(2,jj))
             else if ( mod(i+j+ne-2,ne)+1 == inext) then
               w  (i,ip) = SW(1,jj)/real(ne,wp)+min(SW(2,jj),SW(3,jj))
             else
               w  (i,ip) = SW(1,jj)/real(ne,wp)
             end if
           end do

       else

         !----------------------------------------------------------
         ! Interpolation between central vertex and two cell centers
         !----------------------------------------------------------
#if 1
!NEC$ shortloop
         do i = 1, ne
            x1(:,i) = uvec
            x2(:,i) = nvec(:,i)
            x3(:,i) = nvec(:,mod (i,ne)+1)
         end do
         call Symplectic_Weights (z (:,:ne), x1(:,:ne), x2(:,:ne), x3(:,:ne), &
                                  SW(:,:ne))
         jj   = maxloc (minval (SW(:,:ne),dim=1), dim=1)
         wmin =         minval (SW(:,jj))
         w(     jj      ,ip) = SW(2,jj)
         w(mod (jj,ne)+1,ip) = SW(3,jj)
!NEC$ shortloop
         do i = 1, ne
            idx(i,ip) = jlc(i)
            blk(i,ip) = jbc(i)
            w  (i,ip) = w  (i,ip) + SW(1,jj) / ne
         end do

#else /* old version without weights */
!NEC$ shortloop
         do i = 1, ne
           idx(i,ip) = jlc (mod (j-2 + i, ne) + 1)
           blk(i,ip) = jbc (mod (j-2 + i, ne) + 1)
           w  (i,ip) = 0._wp
         end do
#endif

       end if ! l_triangle
       ngp  (ip) = ne

       case (3)
       !-------------------------------------------------
       ! Original version: Return 3 neighbour grid-points
       !                   and interpolation weights
       !-------------------------------------------------
       if (.not.l_triangle) then
       do i = 1, ns
          jp(i) = mod (j-1 + ne + djp(i), ne) + 1
          jm(i) = mod (j-1 + ne + djm(i), ne) + 1
          if (debug) print *, "jp,jm =", jp(i), jm(i)
       end do
!NEC$ ivdep
       do i = 1, ns
!         z (:,i) = point
          x1(:,i) = nvec(:,j)
          x2(:,i) = nvec(:,jp(i))
          x3(:,i) = nvec(:,jm(i))
       end do

       call Symplectic_Weights (z (:,:ns), x1(:,:ns), x2(:,:ns), x3(:,:ns), &
                                SW(:,:ns))
       jj   = maxloc (minval (SW(:,:ns),dim=1), dim=1)
       wmin = minval (SW(:,jj))

       !---------------------------------------------
       ! store cell indices and interpolation weights
       !---------------------------------------------
       idx(1,ip) = jl
       idx(2,ip) = jlc(jp(jj))
       idx(3,ip) = jlc(jm(jj))
       blk(1,ip) = jb
       blk(2,ip) = jbc(jp(jj))
       blk(3,ip) = jbc(jm(jj))
       else
         ntriangles=4; if (ne==5) ntriangles=3
!NEC$ shortloop
!NEC$ ivdep
         do i = 1, ntriangles
!           z (:,i) = point
            x1(:,i) = nvec(:,triangles((i-1)*3+1))
            x2(:,i) = nvec(:,triangles((i-1)*3+2))
            x3(:,i) = nvec(:,triangles((i-1)*3+3))
         end do
         call Symplectic_Weights (z (:,:ntriangles), x1(:,:ntriangles), &
                                  x2(:,:ntriangles), x3(:,:ntriangles), &
                                  SW(:,:ntriangles))
         jj=1
!NEC$ shortloop
         do j=2,ntriangles
           if (minval(SW(:,j)).gt.minval(SW(:,jj))) then
             jj=j
           end if
         end do

         idx(1,ip) = jlc(triangles((jj-1)*3+1))
         idx(2,ip) = jlc(triangles((jj-1)*3+2))
         idx(3,ip) = jlc(triangles((jj-1)*3+3))
         blk(1,ip) = jbc(triangles((jj-1)*3+1))
         blk(2,ip) = jbc(triangles((jj-1)*3+2))
         blk(3,ip) = jbc(triangles((jj-1)*3+3))
       end if
       w(1:3,ip) = SW (1:3, jj)
       wmin=minval (w(1:3,ip))
       ngp  (ip) = 3

       end select ! nr

!      if (wmin < 0) debug = .true.
       if (wmin < - maxnwgt) debug = .true.

       if (debug) then
          if (wmin < 0) print *, "#################################################"
          print*,l_triangle,nr,ne,ntriangles,ip,jj
          print*,w(:,ip)
          print*,minval(SW (:, jj)),minval(w(:,ip)),wmin
          select case (nr)
          case (6)
            do i = 1, ne
               print *, "Weights", i, ":", SW(:,i)
            end do
          case (3)
            if (l_triangle) then
              do i = 1, ntriangles
                 print *, "Weights", i, ":", SW(:,i)
              end do
            else
              do i = 1, ns
                 print *, "Weights", i, ":", SW(:,i)
              end do
            end if
          end select
          print *, "maxloc (minval (SW))=", jj
          print *
          print *, "Checkme:", j, jp(jj), jm(jj)
          print *, "Checkme (jp)", jp
          print *, "Checkme (jm)", jm
!         print *, "jlc:", jlc
          print *
          if (ip > 1) then
             print *, "ne =", ne
             print *, "Point: ", lon(ip), lat(ip)
             print *, "Vertex:", verts% vertex(kl,kb)% lon * r2d, &
                                 verts% vertex(kl,kb)% lat * r2d
          end if
          print *, "Resulting cell centers:"
          print *, jl,jb,cells% center(jl,jb)% lon * r2d, &
                         cells% center(jl,jb)% lat * r2d
          print *, jlc(jp(jj)),jbc(jp(jj)),                           &
                   cells% center(jlc(jp(jj)),jbc(jp(jj)))% lon * r2d, &
                   cells% center(jlc(jp(jj)),jbc(jp(jj)))% lat * r2d
          print *, jlc(jm(jj)),jbc(jm(jj)),                           &
                   cells% center(jlc(jm(jj)),jbc(jm(jj)))% lon * r2d, &
                   cells% center(jlc(jm(jj)),jbc(jm(jj)))% lat * r2d
          print *
          print *, "Weights:"
          print *, SW(:,jj)
          print *, "Scalar products:"
          print *, sum (point*nvec(:,j)), sum (point*nvec(:,jp(jj))), &
               sum (point*nvec(:,jm(jj)))
          print *
       end if

       if (wmin < - maxnwgt) then
!         print *, "Fatal: negative weight encountered!!!"
!         stop
          if (wmin < - 100*maxnwgt) then
             call finish  ("search_icon_global","negative weight encountered!")
          else
             call message ("search_icon_global","negative weight encountered!")
          end if
       end if
    end do l1

  end subroutine search_icon_global
  !----------------------------------------------------------------------------
  subroutine dist_to_bound (patch, n)
  type (t_patch) ,intent(inout) :: patch
  integer        ,intent(in)    :: n     ! max # of lines to check
  !----------------------------------------------------------------
  ! Determines the distance of grid-cells (in number of gridpoints)
  ! to the boundary of the domain. Same as read from 'refin_c_ctrl'
  ! but with more lines checked. Counting is as follows:
  !
  !  *===*===*===               ====== boundary
  !   \1/1\1/1\1                *      vertices
  !  --*---*---*--              1,2,3  cells
  !     \2/2\2/
  !    --*---*--
  !     /3\3/3\
  !
  !-----------------------------------------------------------------
    integer :: d       ! distance loop index
    integer :: jl      ! line & block index
    integer :: jb      ! block index
    integer :: il      ! line & block index
    integer :: ib      ! block index
    integer :: ne      ! number of neighbours
    integer :: j       ! loop over neighbours
    integer :: nproma_c

    if (patch% cell_type /= 3) call finish ("dist_to_bound", "cell_type /= 3")

    nproma_c = size (patch% cells% center, dim=1)
    if (.not. allocated (patch% cells% c_ctrl))              &
      allocate (patch%cells%c_ctrl (nproma_c, patch% nblks_c))

    patch% cells% c_ctrl = 0
    do d = 0, 2*n+1
      do jb = 1, size (patch% cells% neighbor_idx, 2)
      do jl = 1, size (patch% cells% neighbor_idx, 1)
        if (d == 0) then
          ne = count (patch% cells% neighbor_idx(jl,jb,1:3) > 0)
          if (ne <  patch% cell_type) patch% cells% c_ctrl (jl,jb) = 1
        else
          if (patch% cells% c_ctrl (jl,jb) == d) then
            do j = 1, 3  ! patch% cell_type
               il = patch% cells% neighbor_idx (jl,jb,j)
               if (il < 1) cycle
               ib = patch% cells% neighbor_blk (jl,jb,j)
               if (patch% cells% c_ctrl (il,ib) == 0) &
                   patch% cells% c_ctrl (il,ib) = d+1
            end do
          endif
        endif
      end do
      end do
    end do
    patch% cells% c_ctrl = (patch% cells% c_ctrl + 1) / 2

! do d = 0, n+1
! if(dace% pe==0) print *,'# d count',d,count(patch% cells% c_ctrl == d)
! end do

!   if (dace% lpio) then
!      do jb = 1, size (patch% cells% c_ctrl, 2)
!      do jl = 1, size (patch% cells% c_ctrl, 1)
!         d = patch% cells% c_ctrl(jl,jb)
!         if (d == 0) then
!            write(990,'(i8,i6,i6,2f9.3)') jl, jb, d,     &
!                 patch% cells% center(jl,jb)% lon * r2d, &
!                 patch% cells% center(jl,jb)% lat * r2d
!         else
!            write(991,'(i8,i6,i6,2f9.3)') jl, jb, d,     &
!                 patch% cells% center(jl,jb)% lon * r2d, &
!                 patch% cells% center(jl,jb)% lat * r2d
!         end if
!      end do
!      end do
!   end if

  end subroutine dist_to_bound
  !----------------------------------------------------------------------------
    subroutine check_neighbours (icongrid)
    type(t_grid_icon), intent(inout) :: icongrid ! ICON grid metadata
    !---------------------------------------------------
    ! loop over references in patch% verts% neighbor_id
    !                     and patch% verts% cell_idx
    ! and sort the references for a specific orientation
    !---------------------------------------------------

      integer :: jl ! line & block index
      integer :: jb ! block index
      integer :: ne ! number of neighbours

      do jb = 1, size (icongrid% patch% verts% cell_idx, 2)
      do jl = 1, size (icongrid% patch% verts% cell_idx, 1)
        ne  =          icongrid% patch% verts% num_edges (jl,jb)
        if (ne < 5 .or. ne > 6)                              cycle  ! exclude ..
        if (icongrid% patch% verts% cell_idx (jl,jb,ne) < 1) cycle  ! .. boundaries
        call check_neighbour (icongrid% patch% verts% neighbor_idx (jl, jb, 1:ne), &
                              icongrid% patch% verts% neighbor_blk (jl, jb, 1:ne), &
                              icongrid% patch% verts% cell_idx     (jl, jb, 1:ne), &
                              icongrid% patch% verts% cell_blk     (jl, jb, 1:ne), &
                              icongrid,                                       ne   )
      end do
      end do

    end subroutine check_neighbours
  !----------------------------------------------------------------------------
    subroutine check_neighbour (jlv, jbv, jlc, jbc, icongrid, ne)
    integer                        :: jlv (:)  ! Vertex line  indices
    integer                        :: jbv (:)  ! Vertex block indices
    integer                        :: jlc (:)  ! Cell   line  indices
    integer                        :: jbc (:)  ! Cell   block indices
    type(t_grid_icon), intent(in)  :: icongrid ! ICON grid metadata
    integer          , intent(in)  :: ne
    !---------------------------------------------------------------
    ! For local areas neighbour relationships may be stored in a
    ! different sequence (not consecutive) than for the global grid.
    ! This routine puts them into the expected order.
    !
    ! 1) sort list of neighbour vertices
    ! 2) ensure that list of neighbour cells is consistent
    !---------------------------------------------------------------
      integer  :: i, j
      real(wp) :: nvecv (3,ne)
      real(wp) :: nvec  (3)
      real(wp) :: nvecc (3,ne)
      real(wp) :: pdist   (ne)
      integer  :: ip1 (ne)
      integer  :: ip2 (ne)

      !-----------------------------
      ! store 3d-vectors of vertices
      !-----------------------------
!NEC$ shortloop
      do i = 1, ne
        nvecv(:,i) = icongrid% vert(jlv(i),jbv(i))% x
      end do
      !--------------------------------------------------------
      ! sort according to distance from first (arbitrary) point
      !--------------------------------------------------------
      pdist(1) = 0._wp
!NEC$ shortloop
      do i = 2, ne
        pdist(i) = sum ((nvecv(:,1) - nvecv(:,i)) **2)
      end do
      ip1 = index (pdist)
      do i = 4, 5
        pdist(i) = sum ((nvecv(:,ip1(2)) - nvecv(:,ip1(i))) **2)
      end do
      !------------------------------------------------
      ! ensure that second and third point are close by
      !------------------------------------------------
      ip2 (1:2) = ip1 (1:2)
      ip2 (ne)  = ip1 (3)
      select case (ne)
      case (6)
        ip2 (4) = ip1 (6)
        if (pdist(4) < pdist(5)) then
          ip2 (3) = ip1 (4)
          ip2 (5) = ip1 (5)
        else
          ip2 (3) = ip1 (5)
          ip2 (5) = ip1 (4)
        endif
      case (5)
        if (pdist(4) < pdist(5)) then
          ip2 (3) = ip1 (4)
          ip2 (4) = ip1 (5)
        else
          ip2 (3) = ip1 (5)
          ip2 (4) = ip1 (4)
        endif
      case default
        call finish('check_neighbour','ne /= 5 or 6')
      end select
      !----------------------------------------------------------
      ! now set references in derived type (subroutine arguments)
      !----------------------------------------------------------
!NEC$ shortloop
      jlv = jlv (ip2)
!NEC$ shortloop
      jbv = jbv (ip2)

      !----------------------------------------------------------
      ! finally ensure connectivity of cell centers with vertices
      !----------------------------------------------------------
!NEC$ shortloop
!NEC$ ivdep
      do j = 1, ne
        nvecv(:,j) = icongrid% vert      (jlv(j),jbv(j))% x
        nvecc(:,j) = icongrid% uvec_cell (jlc(j),jbc(j))% x
      end do

      ip2 = -huge (ip2)
      nvec = nvecv (:,ne)
!NEC$ shortloop
      do i = 1, ne
        nvec = (nvec + nvecv(:,i)) / 2._wp
!NEC$ shortloop
        do j = 1, ne
          pdist(j) = 1._wp - sum (nvec * nvecc(:,j))
        end do
!NEC$ shortloop
        ip2 (i) = minloc (pdist, dim=1)
        nvec = nvecv(:,i)
      end do
!NEC$ shortloop
      jlc = jlc (ip2)
!NEC$ shortloop
      jbc = jbc (ip2)

    end subroutine check_neighbour
  !============================================================================
  Subroutine Symplectic_Weights &
    (z,       & ! <-- Point in 3D space
     x1,      & ! <-- 1st vertices
     x2,      & ! <-- 2nd vertices
     x3,      & ! <-- 3rd vertices
     SW)        ! --> Array of symplectic weights
    !
    ! Calculation of symplectic weights for multiple (point/triangle) pairs.
    !----------------------------------------------------------
    ! Method:
    !            (z,x2,x3)     (x1,z,x3)     (x1,x2,z)
    ! SW(1:3) = ------------, ------------, ------------
    !                N             N             N
    ! N is normalization constant: Sum(SW(1:3)) = 1
    !----------------------------------------------------------
    ! Original code:
    ! (C) Copyright 2001, M. E. Gorbunov.
    !----------------------------------------------------------
    ! Input arguments:
    !
    Real(wp), Intent(In)  :: z(:,:)     ! Points
                                        ! [Component,point number]
    Real(wp), Intent(In)  :: x1(:,:)    ! 1st vertex direction
                                        ! [Component,point number]
    Real(wp), Intent(In)  :: x2(:,:)    ! 2nd vertex direction
    Real(wp), Intent(In)  :: x3(:,:)    ! 3rd vertex direction
    !
    ! Output arguments:
    !
    Real(wp), intent(out) :: SW(:,:)    ! Array of symplectic weights
                                        ! [Vertex,point number]
    !----------------------------------------------------------
    ! Local Scalars:
    !
    Integer  :: IP  ! Point index
    Real(wp) :: D   ! Symplectic denominator
    !
    ! Local Arrays:
    !
    real(wp) :: y1(3), y2(3), y3(3)     ! Temporaries
    !----------------------------------------------------------

    Do IP=1,Size(SW,2)

#if 0  /* original version */

       SW(1,IP) = (  z (1,IP)*x2(2,IP)*x3(3,IP)    &
                   - z (1,IP)*x2(3,IP)*x3(2,IP) +  &
                     z (2,IP)*x2(3,IP)*x3(1,IP)    &
                   - z (2,IP)*x2(1,IP)*x3(3,IP) +  &
                     z (3,IP)*x2(1,IP)*x3(2,IP)    &
                   - z (3,IP)*x2(2,IP)*x3(1,IP))

       SW(2,IP) = (  x1(1,IP)*z (2,IP)*x3(3,IP)    &
                   - x1(1,IP)*z (3,IP)*x3(2,IP) +  &
                     x1(2,IP)*z (3,IP)*x3(1,IP)    &
                   - x1(2,IP)*z (1,IP)*x3(3,IP) +  &
                     x1(3,IP)*z (1,IP)*x3(2,IP)    &
                   - x1(3,IP)*z (2,IP)*x3(1,IP))

       SW(3,IP) = (  x1(1,IP)*x2(2,IP)*z (3,IP)    &
                   - x1(1,IP)*x2(3,IP)*z (2,IP) +  &
                     x1(2,IP)*x2(3,IP)*z (1,IP)    &
                   - x1(2,IP)*x2(1,IP)*z (3,IP) +  &
                     x1(3,IP)*x2(1,IP)*z (2,IP)    &
                   - x1(3,IP)*x2(2,IP)*z (1,IP))

#else  /* alternative version for testing of sensitivity to rounding */

       !--------------------------------------------------------
       ! Rewrite determinants into numerically more stable form,
       ! using (z,x2,x3) == (z,x2-z,x3-z) etc.
       !--------------------------------------------------------
!      SW(1,IP) = (  x2(2,IP)*x3(3,IP)              &
!                  - x2(3,IP)*x3(2,IP) ) * z (1,IP) &
!               + (  x2(3,IP)*x3(1,IP)              &
!                  - x2(1,IP)*x3(3,IP) ) * z (2,IP) &
!               + (  x2(1,IP)*x3(2,IP)              &
!                  - x2(2,IP)*x3(1,IP) ) * z (3,IP)
!
!      SW(2,IP) = (  x1(1,IP)*x3(3,IP)              &
!                  - x1(3,IP)*x3(1,IP) ) * z (2,IP) &
!               + (  x1(2,IP)*x3(1,IP)              &
!                  - x1(1,IP)*x3(2,IP))  * z (3,IP) &
!               + (  x1(3,IP)*x3(2,IP)              &
!                  - x1(2,IP)*x3(3,IP) ) * z (1,IP)
!
!      SW(3,IP) = (  x1(1,IP)*x2(2,IP)              &
!                  - x1(2,IP)*x2(1,IP) ) * z (3,IP) &
!               + (  x1(2,IP)*x2(3,IP)              &
!                  - x1(3,IP)*x2(2,IP) ) * z (1,IP) &
!               + (  x1(3,IP)*x2(1,IP)              &
!                  - x1(1,IP)*x2(3,IP) ) * z (2,IP)

       y1(:) = x1(:,IP) - z(:,IP)
       y2(:) = x2(:,IP) - z(:,IP)
       y3(:) = x3(:,IP) - z(:,IP)
       SW(1,IP) = (  y2(2) * y3(3)              &
                   - y2(3) * y3(2) ) * z (1,IP) &
                + (  y2(3) * y3(1)              &
                   - y2(1) * y3(3) ) * z (2,IP) &
                + (  y2(1) * y3(2)              &
                   - y2(2) * y3(1) ) * z (3,IP)

       SW(2,IP) = (  y1(1) * y3(3)              &
                   - y1(3) * y3(1) ) * z (2,IP) &
                + (  y1(2) * y3(1)              &
                   - y1(1) * y3(2))  * z (3,IP) &
                + (  y1(3) * y3(2)              &
                   - y1(2) * y3(3) ) * z (1,IP)

       SW(3,IP) = (  y1(1) * y2(2)              &
                   - y1(2) * y2(1) ) * z (3,IP) &
                + (  y1(2) * y2(3)              &
                   - y1(3) * y2(2) ) * z (1,IP) &
                + (  y1(3) * y2(1)              &
                   - y1(1) * y2(3) ) * z (2,IP)

#endif

#if defined (__NEC__) || 1
       ! Explicit expansion of sum() for better vectorization
       D        = SW(1,IP) + SW(2,IP) + SW(3,IP)
#else
       D        = Sum (SW(1:3,IP), Dim=1)
#endif

!      SW(1,IP) = SW(1,IP)/D
!      SW(2,IP) = SW(2,IP)/D
!      SW(3,IP) = SW(3,IP)/D
       SW(1,IP) = SW(1,IP)*(1._wp/D)
       SW(2,IP) = SW(2,IP)*(1._wp/D)
       SW(3,IP) = SW(3,IP)*(1._wp/D)

    End Do

  End Subroutine Symplectic_Weights
  !============================================================================
  !----------------------------------------------------
  ! Little "benchmark" code for search within ICON grid
  !----------------------------------------------------
  subroutine bench_search_icon (icongrid)
    type(t_grid_icon), intent(in) :: icongrid
    !--
    integer, parameter :: N = 525000     ! Roughly T511
    integer   :: idx(3,N), blk(3,N), ngp(N)
    real(wp)  :: lon(N), lat(N), w(3,N)
    real(wp)  :: t1, t2
    integer   :: idx6(6,N), blk6(6,N)
    real(wp)  :: w6(6,N), sumw(N)
    integer   :: k

    call random_number (lon)
    call random_number (lat)
    ! Throw away first set...
    call random_number (lon)
    call random_number (lat)
    if (icongrid% global) then
       lon = 360._wp *  lon
       lat = 180._wp * (lat - 0.5_wp)
    else
       lon = icongrid% g% lo1 + (icongrid% g% di * icongrid% g% nx) *  lon
       lat = icongrid% g% la1 + (icongrid% g% dj * icongrid% g% ny) *  lat
    end if
    print *
    print *, "Testing random interpolation"
    print *, "range lon:", minval (lon), maxval (lon)
    print *, "range lat:", minval (lat), maxval (lat)
    print *, real (lon(1:6))
    print *, real (lon(N-5:N))
    print *, real (lat(1:6))
    print *, real (lat(N-5:N))
    print *
    print *, "Testing interpolation to triangle"
    call cpu_time (t1)
    call search_icon_global (icongrid, lon, lat, idx, blk, w, ngp)
    call cpu_time (t2)
    print *, real (w(1:3,1))
    print *, real (w(1:3,2))
    print *, real (w(1:3,N-1))
    print *, real (w(1:3,N))
    print *, "Time for search [s]:", real (t2-t1)
    print *, "Time for search [µs/point]:", real ((t2-t1)/N*1.e6_wp)
    print *
    sumw = sum (w, dim=1)
    where (ngp > 0) sumw = sumw - 1
    print *, "Checking sum of weights"
    print *, "min./max. of sum(w)-1 =", minval (sumw), maxval (sumw)
    print *
    print *, "Testing interpolation to hexagon/pentagon"
    call cpu_time (t1)
    call search_icon_global (icongrid, lon, lat, idx6, blk6, w6, ngp, mode=6)
    call cpu_time (t2)
    print *, real (w6(1:6,1))
    print *, real (w6(1:6,2))
    print *, real (w6(1:6,N-1))
    print *, real (w6(1:6,N))
    print *, "Time for search [s]:", real (t2-t1)
    print *, "Time for search [µs/point]:", real ((t2-t1)/N*1.e6_wp)
    print *
    sumw = sum (w6, dim=1)
    where (ngp > 0) sumw = sumw - 1
    print *, "Checking sum of weights"
    print *, "min./max. of sum(w)-1 =", minval (sumw), maxval (sumw)
    print *
    print *, "Checking plausibility of weights"
    do k = 1, 3
       call check_w6 (k)
    end do
    print *
    print *, "Done."
  contains
    subroutine check_w6 (k)
      integer, intent(in) :: k
      !---------------------------------------------------
      ! Check plausibility of weights for the hexagon case
      !---------------------------------------------------
      type(t_patch), pointer :: patch
      real(wp)               :: dlat(6), dlon(6), dist(6)
      integer                :: i, ld, lw
      if (all (w6(:,k) == 0._wp)) return
      print *
      write (*,'(a,2f12.6)') " Point:", lon(k), lat(k)
      patch => icongrid% patch
      do i = 1, ngp(k)
         dlon(i) = patch% cells% center(idx6(i,k),blk6(i,k))% lon * r2d
         dlat(i) = patch% cells% center(idx6(i,k),blk6(i,k))% lat * r2d
         dist(i) = spdist (lon(k), lat(k), dlon(i), dlat(i))      * r2d
         write(*,'(i7,4f12.6)') i, dlon(i), dlat(i), dist(i), w6(i,k)
      end do
      ld = minloc (dist(1:ngp(k)),1)
      lw = minloc (dist(1:ngp(k)),1)
      write(*,'(" minloc(dist) =",i2," maxloc(w) =",i2)') ld, lw
      if (ld /= lw) call finish ("bench_search_icon","checkme!")
    end subroutine check_w6
  end subroutine bench_search_icon
  !============================================================================
END MODULE mo_icon_grid
