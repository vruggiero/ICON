!
!+ Utility functions for atmospheric states
!
MODULE mo_atm_util
!
! Description:
!   Utility functions for atmospheric states (t_atm):
!   - conservative remapping
!   - interpolation
!   - derived type t_atm_set: set/collection of atmospheric states
!   - organize sets of states for parallel I/O
!
! Current Code Owner: DWD, Harald Anlauf
!    phone: +49 69 8062 4941
!    fax:   +49 69 8062 3721
!    email: harald.anlauf@dwd.de
!
! Code Description:
! Language: Fortran.
! Software Standards:
!
! Authors:
! Harald Anlauf      DWD  2019
! Thorsten Steinert  DWD  2020
!==============================================================================

  !=============
  ! Modules used
  !=============
  use mo_kind,       only: wp,                &! working precision kind
                           i8                  ! 64bit integer
  use mo_wmo_tables, only: DWD6_ICON           ! ICON unstructured grid id
  use mo_physics,    only: r2d,               &! conversion factors: radians -> degree
                           gacc,              &! gravity accel.
                           lapse_cl,          &! Clim. lapse rate
                           R,                 &! gas constant
                           q_rh,              &! spec.hum.from rh
                           rh_q                ! rel.hum. from  q
  use mo_physical_constants, only: p0sl_bg     ! [Pa]  sea level pressure
  use mo_exception,  only: finish              ! abort in case of error
  use mo_namelist,   only: position_nml,      &! routine to position nml group
                           nnml,              &! namelist fortran unit number
                           POSITIONED          ! position_nml: OK return flag
  use mo_mpi_dace,   only: dace,              &! MPI group info
                           p_bcast,           &! overloaded broadcast routine
                           p_min,             &! calculate minimum on all PEs
                           p_max,             &! calculate maximum on all PEs
                           p_sum               ! generic sum-function
  use mo_dace_string,only: split,             &! char string -> array
                           byte2hex,          &! char -> 'hexadecimal' CHAR(LEN=2)
                           decode_uuid,       &! Decode uuid (string -> char array)
                           eval_string         ! evaluate strings
  use mo_run_params, only: input,             &! input  directory path
                           path_file,         &! concatenate pathname/basename
                           method,            &! 'ENKF', 'LETKF', 'PSAS+LETKF'
                           ana_time            ! analysis time
  use mo_algorithms, only: init_splinex,      &! set coefficients (2nd deriv.) for spline interpol.
                           splint              ! Spline interpol.with optimization for subseq.calls
  use mo_atm_grid,   only: t_grid,            &! grid meta data type
                           VCT_P_HYB,         &! GME      vertical coordinate
                           VCT_P_ISO,         &! isobaric vertical coordinate
                           VCT_Z_HYB,         &! hybrid z coordinate           (COSMO)
                           VCT_Z_GEN,         &! generalised z coordinate (ICON,COSMO)
                           same_vertical_grid,&! compare vertical grids
                           same_horizontal_grid! compare horizontal grids
  use mo_atm_state,  only: t_atm,             &! atmospheric state derived type
                           construct,         &! initialize t_atm
                           destruct,          &! free memory of t_atm
                           print,             &! print atmospheric state
                           allocate,          &! allocate component of t_atm
                           deallocate,        &! deallocate component of t_atm
                           assignment (=),    &! t_atm = t_atm
                           set_geo,           &! set geopotential height
                           set_p,             &! set pressure
                           hor_intp_coef_atm, &!
                           T_SO                ! position of T_SO in t_atm
  use mo_t_col,      only: t_cols, t_col,     &!
                           t_sl,              &!
                           assignment(=),     &! type(t_cols) = real(wp)
                           alloc_cols,        &! allocate components of t_cols
                           dealloc_cols,      &! deallocate comps. of t_cols
                           get_cols,          &! redistribute model columns
                           mcol_idx,          &!
                           COL_T, COL_RH, COL_UV, COL_W,       &
                           COL_P, COL_Q, COL_QCL, COL_QCI,     &
                           COL_GEO, COL_GEOH, COL_PP,          &
                           COL_INFLAT, COL_SOIL, COL_WSO
  use mo_t_obs,      only: t_icol,            &! interpolation point descriptor
                           t_hic,             &! horizontal interpolation coefficients
                           t_mcols,           &! model column descriptor set
                           t_mcol,            &! model column descriptor
                           destruct,          &! destructor (t_mcols)
                           set_xuv,           &! set up unit vectors on the sphere
                           t_coord             ! unit vector data typ
  use mo_t_icon,     only: nproma,            &! Vector length for blocking
                           blk_no,            &! Block index from global index
                           idx_no              ! Line  index within block
  use mo_grid_intpol,only: t_h_idx             ! type for indices, weights
  use mo_soil,       only: n_sl,              &! number of soil layers
                           st,                &! soil thickness   (mm)
                     ca => cadp,              &! air dryness point   (lower bound)
                     cl => cpwp,              &! plant wilting point (lower bound)
                     cu => cfcap,             &! field capacity      (upper bound)
                     cp => cporv,             &! pore volume         (upper bound)
                           ST_ICE, ST_ROCK,   &!
                           ST_SEAWATER, ST_SEAICE
  !---------
  ! GRIB I/O
  !---------
  use mo_gribtables,   only: dwd_to_v3d,      &! convert to 3dvar convention
                             v3d_to_dwd        ! convert to DWD convention
  use mo_grib_handling,only: t_inventory,     &! Inventory Data Type
                             get_inventory,   &! read inventory table
                             print_inventory, &! print inventory
                             set_defaults      ! set center, subcenter, process
  use mo_grib,         only: read,            &! read atmospheric state or grid
                             write_grib        ! write atmospheric state
  !---------------------
  ! NetCDF f90 interface
  !---------------------
  use netcdf,        only: nf90_open,             &! open   NetCDF file
                           nf90_close,            &! close  NetCDF file
                           nf90_inquire_dimension,&!
                           nf90_get_att,          &! get    attribute
                           nf90_get_var,          &! get variable
                           nf90_inq_dimid,        &! inquire dimension id
                           nf90_inq_varid,        &! inquire variable  id
                           nf90_strerror,         &! error character string
                           NF90_GLOBAL,           &! global attribute flag
                           NF90_NOWRITE,          &! flag: readonly
                           NF90_NOERR              ! NetCDF code for no error

  implicit none

  !================
  ! Public entities
  !================
  private
  public :: init_remap            ! Initialize coefficients for remapping
  public :: remap_atm             ! Remap atmospheric state
  public :: cleanup_remap         ! Deallocate module variables
  public :: mix_q                 ! qv <- rh for lower atmosphere
  public :: interpolate_atm       ! interpolate atmospheric state
  !---------------------------------
  ! Handle sets of atmospheric state
  !---------------------------------
  public :: t_atm_set             ! type: set/collection of atmospheric states
  public :: construct             ! t_atm_set constructor routine
  public :: destruct              ! t_atm_set destructor routine
  public :: allocate              ! associate states (t_atm) with t_atm_set
  public :: deallocate            ! deallocate states in t_atm_set
  public :: write_grib            ! write t_atm_set in parallel
  !------------------
  ! external analysis
  !------------------
  public :: read_nml_extana       ! read namelist /EXTANA/
  public :: prep_extana           ! prepare external analysis
  public :: ext_snow              ! use external snow analysis
  public :: ext_sst               ! use external sst  analysis
  public :: ext_sma               ! use external sma  analysis
  public :: par_extana            ! parameters from external analysis
  public :: read_extana           ! parameters from external analysis
  public :: extana_file           ! external analysis filename
!------------------------------------------------------------------------------
  !================================
  ! Module variables and data types
  !================================
  integer                                :: src_size = -1        ! size of src grid
  integer                                :: dst_size = -1        ! size of dst grid
  integer,  dimension(:),    allocatable :: nlinks               ! number of links
  integer,  dimension(:,:),  allocatable :: ig_in                ! global src indices
  real(wp), dimension(:,:),  allocatable :: w_remap              ! remap coefficients
  character(len=256)                     :: remap_fname_mod = '' ! copy of remap_fname

  !-------------------------------
  ! namelist for external analysis
  !-------------------------------
  character(len=256) :: extana_file    = ''        ! external analysis
  logical            :: ext_snow       = .false.   ! external snow?
  logical            :: ext_sst        = .false.   ! external sst?
  logical            :: ext_sma        = .false.   ! external sma?
  character(len=128) :: par_extana     = ''        ! parameters from external analysis
  character(len=128) :: read_extana    = ''        ! parameters from external analysis
  character(len=64)  :: opt_extana     = ''        ! optional parameters/external analysis
  real(wp)           :: ext_sma_weight = 0.5_wp    ! relaxation factor for external SMA
  integer            :: ext_sma_top    = 1         ! top layer to apply external SMA
  integer            :: ext_sma_bot    = n_sl      ! bottom layer to apply external SMA

  namelist /EXTANA/ extana_file, ext_snow, ext_sst, ext_sma,  &
                    par_extana, opt_extana, ext_sma_weight,   &
                    ext_sma_top, ext_sma_bot

  !-------------------------------------
  ! Set/collection of atmospheric states
  !-------------------------------------
  type t_atm_set
     type(t_atm),        allocatable :: state(:) ! atmospheric states
     character(len=256), allocatable :: fname(:) ! associated file names
     integer,            allocatable :: gpid (:) ! generating process id
     integer,            allocatable :: ptype(:) ! type of gener.process
     integer                         :: size = 0 ! max. number of states
     integer                         :: used = 0 ! currently used states
  end type t_atm_set
!------------------------------------------------------------------------------
  !-----------
  ! interfaces
  !-----------
  interface mix_q
    module procedure mix_q
    module procedure mix_q1
  end interface mix_q

  interface construct
    module procedure construct_atm_set
  end interface construct

  interface destruct
    module procedure destruct_atm_set
  end interface destruct

  interface allocate
    module procedure allocate_patm_0
    module procedure allocate_patm_1
    module procedure allocate_patm_2
  end interface allocate

  interface deallocate
    module procedure deallocate_atm_set
  end interface deallocate

  interface write_grib
    module procedure write_atm_set
  end interface write_grib
!==============================================================================
contains
!==============================================================================
  subroutine init_remap (remap_file)
  !------------------------------------------------------
  ! Initialize coefficients for remapping.
  ! If remap_file is present, it is assumed to provide
  ! precomputed coefficients suitable for SCRIP remapping
  !------------------------------------------------------
  character(len=*),intent(in)              :: remap_file     ! Coefficient filename
    !----------------
    ! Local variables
    !----------------
    integer                                :: i,j            ! indices
    integer                                :: ncid, netcd_id ! for NetCDF
    integer                                :: num_links      ! remap dimension
    integer                                :: num_wgts       ! remap dimension
    integer                                :: dst, nlk       ! temporary
    real(wp)                               :: acc_lim        ! accuracy limit
    real(wp)                               :: renorm_w       ! renormalization
    logical                                :: is_there       ! aux logical variable
    integer, dimension(:),     allocatable :: src_address    ! address vector
    integer, dimension(:),     allocatable :: dst_address    ! address vector
    real(wp), dimension(:),    allocatable :: remap_vec      ! temporary
    real(wp), dimension(:,:),  allocatable :: remap_matrix   ! for read only
    !
    ! Check for remap file
    inquire(file=remap_file, exist=is_there)
    if( .NOT. is_there ) then
      call finish ('init_remap','remap file not found: '//trim(remap_file))
    endif
    !
    remap_fname_mod = remap_file
    !
    if(dace% lpio) then
      !
      ! Read mapping coefficients from file
      call check( nf90_open(remap_fname_mod, NF90_NOWRITE, ncid) )
      !
      call check( nf90_inq_dimid(ncid, "num_links", netcd_id) )
      call check( nf90_inquire_dimension(ncid, netcd_id, len = num_links) )
      !
      call check( nf90_inq_dimid(ncid, "num_wgts", netcd_id) )
      call check( nf90_inquire_dimension(ncid, netcd_id, len = num_wgts) )
      !
      call check( nf90_inq_dimid(ncid, "dst_grid_size", netcd_id) )
      call check( nf90_inquire_dimension(ncid, netcd_id, len = dst_size) )
      !
      call check( nf90_inq_dimid(ncid, "src_grid_size", netcd_id) )
      call check( nf90_inquire_dimension(ncid, netcd_id, len = src_size) )
      !
      allocate( src_address(num_links) )
      allocate( dst_address(num_links) )
      allocate( remap_matrix(num_wgts, num_links) )
      !
      call check( nf90_inq_varid(ncid, "src_address", netcd_id) )
      call check( nf90_get_var(ncid, netcd_id, src_address) )
      !
      call check( nf90_inq_varid(ncid, "dst_address", netcd_id) )
      call check( nf90_get_var(ncid, netcd_id, dst_address) )
      !
      call check( nf90_inq_varid(ncid, "remap_matrix", netcd_id) )
      call check( nf90_get_var(ncid, netcd_id, remap_matrix) )
      !
      call check( nf90_close(ncid) )
      ! prepare broadcast
      allocate( remap_vec(num_links) )
      remap_vec=reshape( remap_matrix, (/ num_links /) )
      deallocate( remap_matrix )
    endif
    !------------------------------------------------------------------------
    ! broadcast information to all pe's
    call p_bcast(num_links, dace% pio)
    call p_bcast(src_size, dace% pio)
    call p_bcast(dst_size, dace% pio)
    !
    if(.not. dace% lpio) then
      allocate( src_address(num_links) )
      allocate( dst_address(num_links) )
      allocate( remap_vec(num_links) )
    endif
    !
    call p_bcast(src_address, dace% pio)
    call p_bcast(dst_address, dace% pio)
    call p_bcast(remap_vec, dace% pio)
    !------------------------------------------------------------------------
    ! save in module variables
    allocate( nlinks  (dst_size) )
    allocate( w_remap (dst_size, 16) )
    allocate( ig_in   (dst_size, 16) )
    !
    ! adjust accuracy to keep number of remap-coeff. below 17
    acc_lim = 0.000325_wp
    do j = 1, 3
      nlinks=0
      ig_in=0
      w_remap=0._wp
      do i = 1, num_links
        dst=dst_address(i)
        nlk=nlinks(dst)+1
        !
        if( remap_vec(i) .lt. acc_lim ) cycle
        !
        nlinks(dst)=nlk
        if( nlk .gt. 16 ) exit
        ig_in(dst,nlk)=src_address(i)
        w_remap(dst,nlk)=remap_vec(i)
      end do
      !
      if( maxval(nlinks) .le. 16 ) exit
      if( j == 3 )  call finish('init_remap','nlk gt 16')
      acc_lim = acc_lim * 2._wp
    end do
    ! renormalization of remapping weights
    do i = 1, dst_size
      renorm_w=0._wp
      do j = 1, nlinks(i)
        renorm_w=renorm_w+w_remap(i,j)
      end do
      renorm_w=1._wp/renorm_w
      do j = 1, nlinks(i)
        w_remap(i,j)=w_remap(i,j)*renorm_w
      end do
    end do
    !
    deallocate( src_address )
    deallocate( dst_address )
    deallocate( remap_vec )

  contains

    !--------------------------------------------------------------------------
    subroutine check(status)
      ! check for errors in nf90-functions
      integer, intent ( in) :: status
      !
      if(status /= nf90_noerr) then
        call finish ("init_remap",trim(nf90_strerror(status)))
      end if
    end subroutine check
    !--------------------------------------------------------------------------

  end subroutine init_remap
!==============================================================================
  subroutine cleanup_remap ()
  !----------------------------
  ! Deallocate module variables
  !----------------------------
    deallocate ( nlinks )
    deallocate ( ig_in )
    deallocate ( w_remap )
  end subroutine cleanup_remap
!==============================================================================
  subroutine remap_atm (y, x, params)
  !------------------------------------------------
  ! Remap atmospheric state
  ! 1: horizontal conservative remapping
  !    for pressure only deviations from ref. atmo.
  ! 2: vertical interpolation of the difference
  !    between the remapped geopot and the
  !    true geopot of the target grid
  !------------------------------------------------
  type(t_atm)     ,intent(inout)         :: y        ! Destination
  type(t_atm)     ,intent(inout)         :: x        ! Source
  character(len=*),intent(in)            :: params   ! Parameters to interpolate

    !----------------
    ! Local variables
    !----------------
    integer               :: np                   ! no. of parameters to interpolate
    integer    ,parameter :: nm = 16              ! max. number of parameters
    character(len=16)     :: pars(nm)             ! names of variables to derive
    integer(i8)           :: iatm                 ! fields to gather
    integer(i8)           :: iatm_vert            ! fields for vertical interpolation
    integer               :: i,j,k,l              ! grid   indices
    integer               :: ke                   ! number of levels
    integer               :: m                    ! processor element
    integer               :: ic,n                 ! column indices
    integer               :: ly(4)                ! lower bounds destination grid
    integer               :: uy(4)                ! upper bounds destination grid
    logical               :: lghx, lgfx           ! geopotential height present in x
    logical               :: lghy, lgfy           ! geopotential height present in y
    type(t_mcols)         :: mc   (dace% npe)     ! model column descriptor set
    type(t_cols)          :: cols (dace% npe)     ! column sets
    real(wp)     ,pointer :: hhl(:)               ! hhl of single column
    type(t_icol) ,pointer :: icols(:,:,:)         ! interpolation point descriptors
    type(t_col)  ,pointer :: b, c                 ! temporary
    type(t_cols)          :: d                    ! temporary
    type(t_atm)           :: src_ref_atm          ! reference atmosphere for source
    type(t_atm)           :: dst_ref_atm          ! reference atmosphere for destination
    type(t_atm)           :: loc_src              ! local copy of source
    type(t_coord)         :: c_src, c_tgt         ! for wind recentering
    real(wp)              :: cuu,cuv,cvu,cvv      ! product of unit vectors
    real(wp) ,allocatable :: tinv(:)              ! temperature inversion
    real(wp) ,allocatable :: true_geof(:)         ! target geof for vert. intp.
    real(wp)              :: t_lapse              ! temperature lapse rate
    logical               :: lt, lrh, lu, lv, lw  ! logical for ...
    logical               :: lpf, lq, lqcl, lqci  !   ... variable interpolation

    ly =  y% lb
    uy =  y% ub
    ke =  uy(3)
    !----------------------------------------------------------------------
    ! compare grids of atmospheric fields with grids used in the remap file
    !----------------------------------------------------------------------
    call check_remap_grids(y% grid, x% grid)
    !-----------------------------------------
    ! set up vertical interpolation parameters
    !-----------------------------------------
    iatm  = COL_GEO
    lghx  = associated (x% geoh)
    lgfx  = associated (x% geof)
    lghy  = associated (y% geoh)
    lgfy  = associated (y% geof)
    if (.not. lgfx) call set_geo (x, geof=.true.)
    if (.not. lgfy) call set_geo (y, geof=.true.)

    !--------------------------------
    ! determine fields to interpolate
    ! convert to DWD conventions
    !--------------------------------
    call split (pars, params, np)
    call v3d_to_dwd (pars (1:np))

    lt    =  any (pars == 'T'      )
    lrh   =  any (pars == 'RELHUM' )
    lu    =  any (pars == 'U'      )
    lv    =  any (pars == 'V'      )
    lw    =  any (pars == 'W'      )
    lpf   =  any (pars == 'P'      )
    lq    =  any (pars == 'QV'     )
    lqcl  =  any (pars == 'QC'     )
    lqci  =  any (pars == 'QI'     )

    if( lu .neqv. lv ) call finish('remap_atm','only one wind component specified')

    !-----------------------------------------------------
    ! set key for atmospheric fields to get from other PEs
    !-----------------------------------------------------
    if (lt)    iatm = iatm + COL_T
    if (lrh)   iatm = iatm + COL_RH
    if (lu )   iatm = ior (iatm, COL_UV)
    if (lv )   iatm = ior (iatm, COL_UV)
    if (lw)    iatm = iatm + COL_W
    if (lpf)   iatm = iatm + COL_P
    if (lq)    iatm = iatm + COL_Q
    if (lqcl)  iatm = iatm + COL_QCL
    if (lqci)  iatm = iatm + COL_QCI

    !-------------------------------
    ! construct reference atmosphere
    !-------------------------------
    if ( lpf ) then
      if ( dace% lpio ) &
        write(6,'(a)') ' remap pressure deviations from ref. atmo.'
      call gen_ref_atm (dst_ref_atm, y% grid)
      call gen_ref_atm (src_ref_atm, x% grid)
    end if

    !-----------------------------------
    ! set key for vertical interpolation
    !-----------------------------------
    iatm_vert = iatm
    if (lpf) iatm_vert = iatm_vert - COL_P

    allocate( true_geof(ke) )

    !-----------------------------------------------
    ! set up model indices and weights for remapping
    !-----------------------------------------------
    nullify (icols)
    call set_mcols_m_remap (mc,      & !  -> list of columns on PEs
                            icols,   & !  -> interpolation point descriptors
                            y% grid, & ! <-  destination grid
                            x% grid, & ! <-  source grid
                            iatm     ) ! <-  parameter bit flag

    !--------------------------------------
    ! gather columns required for remapping
    !--------------------------------------
    call construct(loc_src, template=x)
    loc_src = x
    if (lpf) loc_src% pf = loc_src% pf - src_ref_atm% pf

    call get_cols (mc, loc_src, cols, iatm=iatm)

    call destruct (mc)
    call destruct (src_ref_atm)
    call destruct (loc_src)

    !-------------------------------------
    ! all cols(m) with m /= pe+1 are empty
    !-------------------------------------
    m  =  dace% pe+1

    !-------------------------
    ! horizontal interpolation
    !-------------------------
    call alloc_cols (d, tmp=cols(m), ncol=1)
    b => d% col(1)
    !----------------------------------
    ! loop over destination grid-points
    !----------------------------------
    do l = ly(4), uy(4)
    do j = ly(2), uy(2)
    do i = ly(1), uy(1)
      !-----------------------------------------
      ! intialize destination variables
      ! "loop" over variables, remap/interpolate
      !-----------------------------------------
      b% s% geosp = 0._wp
      b% geo      = 0._wp
      if (lt)   b% t   = 0._wp
      if (lpf)  b% p   = 0._wp
      if (lu)   b% u   = 0._wp
      if (lv)   b% v   = 0._wp
      if (lw)   b% w   = 0._wp
      if (lrh)  b% rh  = 0._wp
      if (lq)   b% q   = 0._wp
      if (lqcl) b% qcl = 0._wp
      if (lqci) b% qci = 0._wp

      do ic = 1, size (icols(i,j,l)% h% imc,1)
        n = icols(i,j,l)% h% imc(ic,1)
        if (n==0) exit

        c => cols(m)% col(n)

        b% s% geosp = b% s% geosp + icols(i,j,l)% h% w(ic) * c% s% geosp
        b% geo      = b% geo      + icols(i,j,l)% h% w(ic) * c% geo
        if (lt)   b% t   = b% t   + icols(i,j,l)% h% w(ic) * c% t
        if (lpf)  b% p   = b% p   + icols(i,j,l)% h% w(ic) * c% p
        if (lw)   b% w   = b% w   + icols(i,j,l)% h% w(ic) * c% w
        if (lrh)  b% rh  = b% rh  + icols(i,j,l)% h% w(ic) * c% rh
        if (lq)   b% q   = b% q   + icols(i,j,l)% h% w(ic) * c% q
        if (lqcl) b% qcl = b% qcl + icols(i,j,l)% h% w(ic) * c% qcl
        if (lqci) b% qci = b% qci + icols(i,j,l)% h% w(ic) * c% qci

        if ( lu .and. lv ) then

          c_src% dlon=c% c% dlon
          c_src% dlat=c% c% dlat
          call set_xuv(c_src)
          c_tgt% dlon=y% grid% rlon(i, j, 1, l) * r2d
          c_tgt% dlat=y% grid% rlat(i, j, 1, l) * r2d
          call set_xuv(c_tgt)
          !
          cuu=c_src% du(1) * c_tgt% du(1) + c_src% du(2) * c_tgt% du(2) &
                                               + c_src% du(3) * c_tgt% du(3)
          cvv=c_src% dv(1) * c_tgt% dv(1) + c_src% dv(2) * c_tgt% dv(2) &
                                               + c_src% dv(3) * c_tgt% dv(3)
          cuv=c_src% du(1) * c_tgt% dv(1) + c_src% du(2) * c_tgt% dv(2) &
                                               + c_src% du(3) * c_tgt% dv(3)
          cvu=c_src% dv(1) * c_tgt% du(1) + c_src% dv(2) * c_tgt% du(2) &
                                               + c_src% dv(3) * c_tgt% du(3)
          !
          b% u = b% u + icols(i,j,l)% h% w(ic) * ( c% u * cuu + c% v * cvu )
          b% v = b% v + icols(i,j,l)% h% w(ic) * ( c% u * cuv + c% v * cvv )

        end if
      end do

      !-----------------------------------------------------------
      ! vertical interpolation between remapped geof and true geof
      !-----------------------------------------------------------
      hhl => y% grid% hhl(i,j,:,l)
      do k = 1, ke
        true_geof(k) = gacc * 0.5_wp * ( hhl(k) + hhl(k+1) )
      end do

      call derive_tinv_sngl_col(b, tinv, t_lapse)
      b% t = b% t - tinv
      call vert_intp_sngl_col (b, true_geof, iatm_vert, d_t=t_lapse, geo=.true.)
      b% t = b% t + tinv
      deallocate(tinv)

      if (lt)   y%   t(i,j,:,l) = b% t
      if (lpf)  y%  pf(i,j,:,l) = b% p
      if (lu)   y%   u(i,j,:,l) = b% u
      if (lv)   y%   v(i,j,:,l) = b% v
      if (lw)   y%   w(i,j,:,l) = b% w
      if (lrh)  y%  rh(i,j,:,l) = b% rh
      if (lq)   y%   q(i,j,:,l) = b% q
      if (lqcl) y% qcl(i,j,:,l) = b% qcl
      if (lqci) y% qci(i,j,:,l) = b% qci

    end do ! i
    end do ! j
    end do ! l

    !-------------------------
    ! add reference atmosphere
    !-------------------------
    if (lpf) y% pf = y% pf + dst_ref_atm% pf

    !--------
    ! cleanup
    !--------
    call dealloc_cols (d)
    call dealloc_cols (cols)
    deallocate        (icols)
    deallocate        (true_geof)
    call destruct     (dst_ref_atm)

    if (.not. lghx) call deallocate (x,'geoh')
    if (.not. lgfx) call deallocate (x,'geof')
    if (.not. lghy) call deallocate (y,'geoh')
    if (.not. lgfy) call deallocate (y,'geof')

  end subroutine remap_atm
!==============================================================================
  subroutine check_remap_grids (grid_dst, grid_src)
  !------------------------------------------------------------------
  ! compare grids of atmospheric fields with grids used in remap file
  !------------------------------------------------------------------
  type(t_grid)     ,pointer              :: grid_dst           ! destination grid
  type(t_grid)     ,pointer              :: grid_src           ! source grid

    !----------------
    ! local variables
    !----------------
    integer                 :: ncid                ! for netcdf
    integer                 :: ierr                ! nf90 status
    character(len=40)       :: dst_grid_info       ! destination grid description
    character(len=40)       :: src_grid_info       ! source grid description

    !------------------
    ! uuid of ICON grid
    !------------------
    character(len=1)        :: uuid(16)            ! uuidOfHGrid

    !-------------------
    ! compare grid sizes
    !-------------------
    if(grid_dst% nxny /= dst_size) &
      call finish ('check_remap_grids','dst-grid: nxny mismatch')
    if(grid_src% nxny /= src_size) &
      call finish ('check_remap_grids','src-grid: nxny mismatch')

    !---------------------------------------
    ! read grid descriptions from remap file
    !---------------------------------------
    if ( dace% lpio ) then
      ierr = nf90_open(remap_fname_mod, NF90_NOWRITE, ncid)
      if(ierr == nf90_noerr) &
        ierr = nf90_get_att(ncid, NF90_GLOBAL, "dst_grid_info", dst_grid_info)
      if(ierr == nf90_noerr) &
        ierr = nf90_get_att(ncid, NF90_GLOBAL, "src_grid_info", src_grid_info)
      if(ierr == nf90_noerr) ierr = nf90_close(ncid)
      if(ierr /= nf90_noerr) &
        call finish ('check_remap_grids',trim(nf90_strerror(ierr)))
    end if
    call p_bcast(dst_grid_info, dace% pio)
    call p_bcast(src_grid_info, dace% pio)

    !-----------------------------------------------
    ! compare grid description with destination grid
    !-----------------------------------------------
    select case (grid_dst% gridtype)
    case (DWD6_ICON)
      ! ICON grid description is the uuid
      call decode_uuid (trim(dst_grid_info), uuid)
      if (any( grid_dst% icongrid% uuid /= uuid ) ) then
        write(6,'(a)') ' Expected uuid : '//byte2hex(grid_dst% icongrid% uuid)
        write(6,'(a)') ' Received uuid : '//byte2hex(uuid)
        call finish('check_remap_grids','dst-grid: uuid mismatch')
      end if
    case default
      call finish('check_remap_grids','unsupported destination grid')
    end select

    !------------------------------------------
    ! compare grid description with source grid
    !------------------------------------------
    select case (grid_src% gridtype)
    case (DWD6_ICON)
      ! ICON grid description is the uuid
      call decode_uuid (trim(src_grid_info), uuid)
      if (any( grid_src% icongrid% uuid /= uuid ) ) then
        write(6,'(a)') ' Expected uuid : '//byte2hex(grid_src% icongrid% uuid)
        write(6,'(a)') ' Received uuid : '//byte2hex(uuid)
        call finish('check_remap_grids','src-grid: uuid mismatch')
      end if
    case default
      call finish('check_remap_grids','unsupported source grid')
    end select

  end subroutine check_remap_grids
!==============================================================================
  subroutine set_mcols_m_remap (mc, ic, gy, gx, iatm)
  !---------------------------------------
  ! set t_mcols for grid to grid remapping
  !---------------------------------------
  type(t_mcols) ,intent(out) :: mc(dace% npe) ! set of model column descriptors
                                              ! (source grid)
  type(t_icol)  ,pointer     :: ic(:,:,:)     ! interpol. point (column) descriptor
                                              ! (destination grid)
  type(t_grid)  ,intent(in)  :: gy            ! destination grid
  type(t_grid)  ,intent(in)  :: gx            ! source      grid
  integer(i8)   ,intent(in)  :: iatm          ! field to gather flags

    !----------------
    ! local variables
    !----------------
    integer                               :: i,j,l ! indices
    integer                               :: ly(4) ! lower bounds destination grid
    integer                               :: uy(4) ! upper bounds destination grid
    integer                               :: lx(4) ! lower bounds source      grid
    integer                               :: ux(4) ! upper bounds source      grid
    integer                               :: p1
    integer                               :: ig    ! global index
    integer, dimension(:),    allocatable :: ig_a  ! global index array
    !---------------------------------------------
    ! determine columns required for interpolation
    !---------------------------------------------
    ly =  gy% lb
    uy =  gy% ub
    allocate (ic(ly(1):uy(1),ly(2):uy(2),ly(4):uy(4)))
    ! Coordinates of destination grids
    ic%c% dlon  = gy% rlon(ly(1):uy(1),ly(2):uy(2),1,ly(4):uy(4)) * r2d
    ic%c% dlat  = gy% rlat(ly(1):uy(1),ly(2):uy(2),1,ly(4):uy(4)) * r2d

    p1 =  dace% pe+1
    lx =  gx% lbg
    ux =  gx% ubg
    allocate (mc(p1)% idx (lx(1):ux(1),lx(2):ux(2),lx(4):ux(4),1))
    allocate (mc(p1)% c (gx% nxny / dace% npe))
    mc(p1)% idx = 0
    mc    % n   = 0
    mc    % pe  = (/ (i,i=0,dace% npe-1) /)

    ! Map (permutated) global index to global index used in remap file
    allocate( ig_a( gy% nxny ) )
    select case(gy% gridtype)
    case (DWD6_ICON)
      do i = 1, gy% nxny
        ig_a( gy% marr(2,i,1,1) + nproma * ( gy% marr(3,i,1,1) - 1 ) ) = i
      end do
    case default
      call finish('set_mcols_m_remap','unsupported grid')
    end select

    do l = ly(4), uy(4)
      do j = ly(2), uy(2)
        do i = ly(1), uy(1)

          ! Determine global index
          ig = ig_a ( i + nproma*(j-1) )

          call idx_init_remap ( &!
           ig,           &! <-  column descriptor (global index)
           ic(i,j,l)% h, &!  -> interpolation coefficients
           mc(p1),       &! <-> model column descriptors
           iatm,         &! <-  fields required
           0,            &! <-  tracers required
           gx            )! <-  model grid
        end do
      end do
    end do

    if (gy% global) then
      if ( any( ic% h% imc(1,1) == 0 ) ) &
        call finish('set_mcols_m_remap','out of domain condition')
    end if

    deallocate (ig_a)
    deallocate (mc(p1)% idx)

  end subroutine set_mcols_m_remap
!==============================================================================
  Subroutine idx_init_remap (ig_out, hic, mc, iatm, itrac, grid)
  !---------------------------------------------------------------------
  ! Determine model grid indices for a given location given by 'coord'.
  ! The horizontal interpolation coefficients 'hic' and the list
  ! of required model columns 'mc' and are updated accordingly.
  !---------------------------------------------------------------------
  integer       ,intent(in)    :: ig_out  ! global index (within output grid)
  type(t_hic)   ,intent(inout) :: hic     ! horizontal interpolation coefficients
  type(t_mcols) ,intent(inout) :: mc      ! model column descriptor
  integer(i8)   ,intent(in)    :: iatm    ! ids of atmospheric parameters required
  integer       ,intent(in)    :: itrac   ! tracers required
  type(t_grid)  ,intent(in)    :: grid    ! input model grid description

    !----------------
    ! local variables
    !----------------
    integer               :: idx  (16,4)  ! (point = 1...16 aktuell, ijdp)
    integer               :: i, ix, inn
    integer               :: ii, jj       ! line & block index in remap file
    integer               :: np
    type(t_mcol) ,pointer :: c, tmp(:)
    type(t_hic)  ,save    :: hic0

    hic = hic0

  !----------------------------
  ! determine neighbour columns
  !----------------------------

    np = nlinks(ig_out)
    do i = 1, np
      ii = idx_no (ig_in(ig_out,i))
      jj = blk_no (ig_in(ig_out,i))
      ! Map global index used in remap file to (permutated) global index
      idx(i,1) = grid% marr (2, ii, jj, 1)                    ! i
      idx(i,2) = grid% marr (3, ii, jj, 1)                    ! j
      idx(i,3) = grid% marr (4, ii, jj, 1)                    ! d = 1 for ICON
      idx(i,4) = grid% marr (1, ii, jj, 1)                    ! processor element
      hic%w(i) = w_remap (ig_out,i)
    end do

    !-------------------------------------------------
    ! copy indices and weights to subroutine arguments
    !-------------------------------------------------
    hic% imc (:,:)   = 0
    hic% w   (np+1:) = 0._wp

    ! Add columns required by remapping
    do i = 1, np
      ix = mc% idx (idx(i,1),idx(i,2),idx(i,3),1)
      if (ix == 0) then
        mc% n = mc% n + 1
        ix = mc% n
        mc% idx (idx(i,1),idx(i,2),idx(i,3),1) = ix
        if (size (mc% c) < ix) then
           tmp => mc% c
           allocate (mc% c (int(1.5 * ix)))
           mc% c (1:size(tmp)) = tmp
           deallocate (tmp)
        endif
        c => mc% c(ix)
        c% icol  = ix
        c% iatm  = iatm
        c% itrac = itrac
        c% ijdtp = (/idx(i,:3),1,idx(i,4)/)
      else
        c => mc% c(ix)
        c% iatm  = ior (iatm,  c% iatm)
        c% itrac = ior (itrac, c% itrac)
      endif
      hic% imc (i,1) = ix
    end do
    ! for now use first col to define position of the interpolation
    if (np > 0) then
      inn=1
      hic% ijdp(1:3) = mc% c (hic% imc (inn,1))% ijdtp(1:3)
      hic% ijdp(4)   = mc% c (hic% imc (inn,1))% ijdtp(  5)
    end if

  end Subroutine idx_init_remap
!==============================================================================
  subroutine gen_ref_atm (ref_atm, grid, lps)
  !------------------------------
  ! generate reference atmosphere
  !------------------------------
  type(t_atm)      ,intent(out)          :: ref_atm        ! reference atmosphere
  type(t_grid)     ,pointer              :: grid           ! grid information
  logical          ,intent(in)           :: lps            ! flag for PS
  optional                               :: lps

    !----------------
    ! local variables
    !----------------
    integer                 :: i,j,k,l                    ! indices
    integer                 :: lb(4)                      ! lower bounds grid
    integer                 :: ub(4)                      ! upper bounds grid
    real(wp)                :: t_ref                      ! ref. temperature [K]
    real(wp)                :: p_ref                      ! ref. pressure [Pa]
    real(wp)                :: z                          ! height [m]
    real(wp)                :: a, b                       ! temporary
    real(wp), PARAMETER     :: h_inv = 0.0001_wp          ! 1/ 10 km [1/m]
    real(wp), PARAMETER     :: t_00=213.15_wp             ! [K]
    real(wp), PARAMETER     :: delta_t=75._wp             ! [K]
    real(wp), PARAMETER     :: cst=-gacc/(h_inv*R*t_00)   ! -gacc*10km / (rd*t_00)
    logical                 :: lps_loc                    ! local copy of lps

    lps_loc = .false.
    if( present(lps) ) lps_loc = lps

    if(.not. associated( grid% hhl )) &
         call finish('gen_ref_atm','hhl is not associated')

    ! create new t_atm to store reference atmosphere
    call construct(ref_atm, grid, alloc = 't pf')
    if (lps_loc) call allocate (ref_atm, 'ps')

    lb=ref_atm% lb
    ub=ref_atm% ub

    do l = lb(4), ub(4)
      do k = lb(3), ub(3)
        do j = lb(2), ub(2)
          do i = lb(1), ub(1)

            z = 0.5_wp*(ref_atm% grid% hhl(i,j,k,l)+ref_atm% grid% hhl(i,j,k+1,l))

            a     = exp(-z*h_inv)
            t_ref = t_00 + delta_t * a

            b     = log( t_ref /( a * ( t_00 + delta_t) ) )
            p_ref = p0sl_bg * exp( cst * b )

            ref_atm% t (i,j,k,l) = t_ref
            ref_atm% pf(i,j,k,l) = p_ref

            if ( lps_loc .and. k == ub(3) ) then
              a     = exp( - ref_atm% grid% hhl(i,j,ub(3)+1,l) * h_inv)
              t_ref = t_00 + delta_t * a

              b     = log( t_ref /( a * ( t_00 + delta_t) ) )
              p_ref = p0sl_bg * exp( cst * b )

              ref_atm% ps(i,j,1,l) = p_ref
            end if

          end do
        end do
      end do
    end do

  end subroutine gen_ref_atm
!==============================================================================
  subroutine vert_intp_sngl_col (col, x_int, ids, d_t, d_p, geo)
  !---------------------------------------------
  ! vertical interpolation of t_col on log( pf )
  ! interpolation on geopotential for geo = true
  !---------------------------------------------
  type(t_col)          ,intent(inout)             :: col         ! column
  real(wp)             ,intent(in)                :: x_int(:)    ! new vertical coordinates
  integer(i8)          ,intent(in)                :: ids         ! fields to interpolate
  real(wp)             ,intent(in) ,optional      :: d_t         ! temperature surface derivative
  real(wp)             ,intent(in) ,optional      :: d_p         ! pressure surface derivative
  logical              ,intent(in) ,optional      :: geo         ! interpolate geopotential levels

    !----------------
    ! local variables
    !----------------
    integer                 :: ke                         ! number of levels
    integer                 :: k                          ! indices
    integer                 :: k0                         ! highest nonflat level
    logical                 :: lgeo                       ! local copy of geo
    real(wp)                :: dt                         ! derivative of surface temperature
    real(wp)                :: dp                         ! derivative of surface pressure
    real(wp)                :: T0                         ! surface temperature
    real(wp) ,allocatable   :: xi(:)                      ! abscissa
    real(wp) ,allocatable   :: yi(:)                      ! interpol. function value
    real(wp) ,allocatable   :: y2(:)                      ! second derivative of function value
    real(wp) ,allocatable   :: x (:)                      ! target abscissa
    real(wp) ,allocatable   :: y (:)                      ! interpolated fields
    real(wp), parameter     :: ACC = 1e-6_wp              !
    logical                 :: lt, lrh, lu, lv, lw, lpf   ! logical for ...
    logical                 :: lq, lqcl, lqci, lgp        !   ... variable interpolation
    !---------------------------------------------------------------------
    ! Upper temperature bound for geopotential extrapolation below surface
    ! cf. IFS documentation (CY25R1), Part II, Paragraph 5.3.2
    !---------------------------------------------------------------------
    real(wp) ,parameter   :: Tx = 290.5_wp

    lgeo = .false. ; if (present(geo)) lgeo = geo

    lt   = iand(ids,COL_T   )/=0_i8
    lrh  = iand(ids,COL_RH  )/=0_i8
    lu   = iand(ids,COL_UV  )/=0_i8
    lv   = iand(ids,COL_UV  )/=0_i8
    lw   = iand(ids,COL_W   )/=0_i8
    lpf  = iand(ids,COL_P   )/=0_i8
    lq   = iand(ids,COL_Q   )/=0_i8
    lqcl = iand(ids,COL_QCL )/=0_i8
    lqci = iand(ids,COL_QCI )/=0_i8
    lgp  = iand(ids,COL_GEO )/=0_i8

    !-------------------
    ! consistency checks
    !-------------------
    if (lgeo) then
      if(.not.associated(col% geo)) call finish('vert_intp_sngl_col','geo not associated')
      ke = size(col% geo)
    else
      if(.not.associated(col% p))   call finish('vert_intp_sngl_col','p not associated')
      ke = size(col% p)
    end if
    if ( size(x_int) /= ke )        call finish('vert_intp_sngl_col','ke mismatch')

    !--------------------
    ! surface derivatives
    !--------------------
    if (lt) then
      T0 = col% t(ke)
      if (present(d_t)) then
        dt = - d_t
      else
        if(lgeo) then
          dt = lapse_cl / gacc
        else
          dt = ( R * lapse_cl / gacc ) * T0
        end if
      end if
    else
      T0 = Tx
    end if
    if (present(d_p)) then
      dp = - d_p
    else
      if(lgeo) then
        dp = 1._wp / ( R * T0 )
      else
        dp = - R * T0
      end if
    end if

    allocate( xi(ke) )
    allocate( yi(ke) )
    allocate( y2(ke) )
    allocate( x (ke) )
    allocate( y (ke) )

    !-----------------
    ! set x-coordinate
    !-----------------
    if (lgeo) then
      xi = - col% geo
      x  = - x_int
    else
      xi = log(col% p)
      x  = log(x_int)
    end if

    !-----------------------------
    ! search highest nonflat level
    !-----------------------------
    k0 = ke + 1
    do k=1,ke
      if ( abs( xi(k) - x(k) ) > ACC ) then
        k0 = k
        exit
      end if
    end do

    !---------------------------------------------------------------------------
    ! interpolate pressure only if lgeo=true and geopotential only if lgeo=false
    !---------------------------------------------------------------------------
    if (lgeo) then
      if(lpf) then
        yi = log(col% p)
        call init_splinex(xi,yi,y2, ypn=dp)
        do k=k0,ke
          if( x(k) <= xi(ke) ) then
            call splint (xi, yi, y2, x(k), y(k), sequential = .true.)
          else
            y(k) = yi(ke) + (x(k) - xi(ke)) * dp
          end if
        end do
        col% p(k0:ke) = exp(y(k0:ke))
      end if
    else
      if(lgp) then
        yi = col% geo
        call init_splinex(xi,yi,y2, ypn=dp)
        do k=k0,ke
          if( x(k) <= xi(ke) ) then
            call splint (xi, yi, y2, x(k), y(k), sequential = .true.)
          else
            y(k) = yi(ke) + (x(k) - xi(ke)) * dp
          end if
        end do
        col% geo(k0:ke) = y(k0:ke)
      end if
    end if

    if(lt) then
      yi = col% t
      call init_splinex(xi,yi,y2, ypn=dt)
      do k=k0,ke
        if( x(k) <= xi(ke) ) then
          call splint (xi, yi, y2, x(k), y(k), sequential = .true.)
        else
          y(k) = yi(ke) + (x(k) - xi(ke)) * dt
        end if
      end do
      col% t(k0:ke) = y(k0:ke)
    end if

    if(lu) then
      yi = col% u
      call init_splinex(xi,yi,y2)
      do k=k0,ke
        if( x(k) <= xi(ke) ) then
          call splint (xi, yi, y2, x(k), y(k), sequential = .true.)
        else
          y(k) = yi(ke)
        end if
      end do
      col% u(k0:ke) = y(k0:ke)
    end if

    if(lv) then
      yi = col% v
      call init_splinex(xi,yi,y2)
      do k=k0,ke
        if( x(k) <= xi(ke) ) then
          call splint (xi, yi, y2, x(k), y(k), sequential = .true.)
        else
          y(k) = yi(ke)
        end if
      end do
      col% v(k0:ke) = y(k0:ke)
    end if

    if(lw) then
      yi = col% w
      call init_splinex(xi,yi,y2)
      do k=k0,ke
        if( x(k) <= xi(ke) ) then
          call splint (xi, yi, y2, x(k), y(k), sequential = .true.)
        else
          y(k) = yi(ke)
        end if
      end do
      col% w(k0:ke) = y(k0:ke)
    end if

    if(lq) then
      yi = col% q
      call init_splinex(xi,yi,y2)
      do k=k0,ke
        if( x(k) <= xi(ke) ) then
          call splint (xi, yi, y2, x(k), y(k), sequential = .true.)
        else
          y(k) = yi(ke)
        end if
      end do
      col% q(k0:ke) = y(k0:ke)
    end if

    if(lrh) then
      yi = col% rh
      call init_splinex(xi,yi,y2)
      do k=k0,ke
        if( x(k) <= xi(ke) ) then
          call splint (xi, yi, y2, x(k), y(k), sequential = .true.)
        else
          y(k) = yi(ke)
        end if
      end do
      col% rh(k0:ke) = y(k0:ke)
    end if

    if(lqcl) then
      yi = col% qcl
      call init_splinex(xi,yi,y2)
      do k=k0,ke
        if( x(k) <= xi(ke) ) then
          call splint (xi, yi, y2, x(k), y(k), sequential = .true.)
        else
          y(k) = yi(ke)
        end if
      end do
      col% qcl(k0:ke) = y(k0:ke)
    end if

    if(lqci) then
      yi = col% qci
      call init_splinex(xi,yi,y2)
      do k=k0,ke
        if( x(k) <= xi(ke) ) then
          call splint (xi, yi, y2, x(k), y(k), sequential = .true.)
        else
          y(k) = yi(ke)
        end if
      end do
      col% qci(k0:ke) = y(k0:ke)
    end if

    !-----------------------------
    ! save new vertical coordinate
    !-----------------------------
    if(lgeo) then
      col% geo = - x
    else
      col% p = exp(x)
    end if

    deallocate(xi)
    deallocate(yi)
    deallocate(y2)
    deallocate(x )
    deallocate(y )

  end subroutine vert_intp_sngl_col
!==============================================================================
  subroutine vert_intp_cols (cols, x_int, ids_col, d_t, d_p, geo)
  !----------------------------------------------
  ! currently not used !
  !----------------------------------------------
  ! vertical interpolation of t_cols on log( pf )
  ! interpolation on geopotential for geo = true
  !----------------------------------------------
  type(t_cols)         ,intent(inout)             :: cols        ! columns
  real(wp)             ,intent(in)                :: x_int(:,:)  ! new vertical coordinates
  integer(i8)          ,intent(in) ,optional      :: ids_col     ! specify variables to interpolate
  real(wp)             ,intent(in) ,optional      :: d_t(:)      ! temperature surface derivative
  real(wp)             ,intent(in) ,optional      :: d_p(:)      ! pressure surface derivative
  logical              ,intent(in) ,optional      :: geo         ! interpolate geopotential levels

    !----------------
    ! local variables
    !----------------
    integer                 :: ke                         ! number of levels
    integer                 :: ncol                       ! number of columns
    integer(i8)             :: ids                        ! fields to interpolate
    integer                 :: j,k                        ! indices
    integer                 :: k0                         ! highest nonflat level
    logical                 :: lgeo                       ! local copy of geo
    type (t_col) ,pointer   :: c                          ! column pointer
    real(wp)                :: dt(cols% ncol)             ! derivative of surface temperature
    real(wp)                :: dp(cols% ncol)             ! derivative of surface pressure
    real(wp)                :: T0(cols% ncol)             ! surface temperature
    real(wp)                :: d_surf                     ! surface derivative
    real(wp) ,allocatable   :: xi(:)                      ! abscissa
    real(wp) ,allocatable   :: yi(:)                      ! interpol. function value
    real(wp) ,allocatable   :: y2(:)                      ! second derivative of function value
    real(wp) ,allocatable   :: x (:)                      ! target abscissa
    real(wp) ,allocatable   :: y (:)                      ! interpolated fields
    real(wp), parameter     :: ACC = 1e-6_wp              !
    logical                 :: lt, lrh, lu, lv, lw, lpf   ! logical for ...
    logical                 :: lq, lqcl, lqci, lgp        !   ... variable interpolation
    !---------------------------------------------------------------------
    ! Upper temperature bound for geopotential extrapolation below surface
    ! cf. IFS documentation (CY25R1), Part II, Paragraph 5.3.2
    !---------------------------------------------------------------------
    real(wp) ,parameter   :: Tx = 290.5_wp

    ke   = cols% ke
    ncol = cols% ncol
    ids  = cols% ids ; if (present(ids_col)) ids  = ids_col
    lgeo = .false.   ; if (present(geo    )) lgeo = geo

    lt   = iand(ids,COL_T   )/=0_i8
    lrh  = iand(ids,COL_RH  )/=0_i8
    lu   = iand(ids,COL_UV  )/=0_i8
    lv   = iand(ids,COL_UV  )/=0_i8
    lw   = iand(ids,COL_W   )/=0_i8
    lpf  = iand(ids,COL_P   )/=0_i8
    lq   = iand(ids,COL_Q   )/=0_i8
    lqcl = iand(ids,COL_QCL )/=0_i8
    lqci = iand(ids,COL_QCI )/=0_i8
    lgp  = iand(ids,COL_GEO )/=0_i8

    !-------------------
    ! consistency checks
    !-------------------
    if ( size(x_int,1) /= ke )    call finish('vert_intp_cols','ke mismatch')
    if ( size(x_int,2) /= ncol )  call finish('vert_intp_cols','ncol mismatch')
    if( ncol == 0 ) return

    !--------------------
    ! surface derivatives
    !--------------------
    if (lt) then
      do j=1,ncol
        T0(j) = cols% col(j)% t(ke)
      end do
      !
      if (present(d_t)) then
        if( size(d_t) /= ncol ) call finish('vert_intp_cols','wrong dimension: d_t')
        dt = - d_t
      else
        if(lgeo) then
          dt = lapse_cl / gacc
        else
          dt = ( R * lapse_cl / gacc ) * T0
        end if
      end if
    else
      T0 = Tx
    end if
    !
    if (present(d_p)) then
      if( size(d_p) /= ncol ) call finish('vert_intp_cols','wrong dimension: d_p')
      dp = - d_p
    else
      if(lgeo) then
        dp = 1._wp / ( R * T0 )
      else
        dp = - R * T0
      end if
    end if

    allocate( xi(ke) )
    allocate( yi(ke) )
    allocate( y2(ke) )
    allocate( x (ke) )
    allocate( y (ke) )

    !-----------------------
    ! vertical interpolation
    !-----------------------
    do j=1,ncol
      c => cols% col(j)

      !-----------------
      ! set x-coordinate
      !-----------------
      if(lgeo) then
        if( .not. associated(c% geo) ) call finish('vert_intp_cols','geo not associated')
        xi = - c% geo
        x  = - x_int(:,j)
      else
        if( .not. associated(c% p) ) call finish('vert_intp_cols','p not associated')
        xi = log(c% p)
        x  = log(x_int(:,j))
      end if

      !-----------------------------
      ! search highest nonflat level
      !-----------------------------
      k0 = ke + 1
      do k=1,ke
        if ( abs( xi(k) - x(k) ) > ACC ) then
          k0 = k
          exit
        end if
      end do
      if ( k0 == ke + 1 ) cycle

      !---------------------------------------------------------------------------
      ! interpolate pressure only if lgeo=true and geopotential only if lgeo=false
      !---------------------------------------------------------------------------
      if (lgeo) then
        if(lpf) then
          yi = log(c% p)
          d_surf = dp(j)
          call init_splinex(xi,yi,y2, ypn=d_surf)
          do k=k0,ke
            if( x(k) <= xi(ke) ) then
              call splint (xi, yi, y2, x(k), y(k), sequential = .true.)
            else
              y(k) = yi(ke) + (x(k) - xi(ke)) * d_surf
            end if
          end do
          c% p(k0:ke) = exp(y(k0:ke))
        end if
      else
        if(lgp) then
          yi = c% geo
          d_surf = dp(j)
          call init_splinex(xi,yi,y2, ypn=d_surf)
          do k=k0,ke
            if( x(k) <= xi(ke) ) then
              call splint (xi, yi, y2, x(k), y(k), sequential = .true.)
            else
              y(k) = yi(ke) + (x(k) - xi(ke)) * d_surf
            end if
          end do
          c% geo(k0:ke) = y(k0:ke)
        end if
      end if

      if(lt) then
        yi = c% t
        d_surf = dt(j)
        call init_splinex(xi,yi,y2, ypn=d_surf)
        do k=k0,ke
          if( x(k) <= xi(ke) ) then
            call splint (xi, yi, y2, x(k), y(k), sequential = .true.)
          else
            y(k) = yi(ke) + (x(k) - xi(ke)) * d_surf
          end if
        end do
        c% t(k0:ke) = y(k0:ke)
      end if

      if(lu) then
        yi = c% u
        call init_splinex(xi,yi,y2)
        do k=k0,ke
          if( x(k) <= xi(ke) ) then
            call splint (xi, yi, y2, x(k), y(k), sequential = .true.)
          else
            y(k) = yi(ke)
          end if
        end do
        c% u(k0:ke) = y(k0:ke)
      end if

      if(lv) then
        yi = c% v
        call init_splinex(xi,yi,y2)
        do k=k0,ke
          if( x(k) <= xi(ke) ) then
            call splint (xi, yi, y2, x(k), y(k), sequential = .true.)
          else
            y(k) = yi(ke)
          end if
        end do
        c% v(k0:ke) = y(k0:ke)
      end if

      if(lw) then
        yi = c% w
        call init_splinex(xi,yi,y2)
        do k=k0,ke
          if( x(k) <= xi(ke) ) then
            call splint (xi, yi, y2, x(k), y(k), sequential = .true.)
          else
            y(k) = yi(ke)
          end if
        end do
        c% w(k0:ke) = y(k0:ke)
      end if

      if(lq) then
        yi = c% q
        call init_splinex(xi,yi,y2)
        do k=k0,ke
          if( x(k) <= xi(ke) ) then
            call splint (xi, yi, y2, x(k), y(k), sequential = .true.)
          else
            y(k) = yi(ke)
          end if
        end do
        c% q(k0:ke) = y(k0:ke)
      end if

      if(lrh) then
        yi = c% rh
        call init_splinex(xi,yi,y2)
        do k=k0,ke
          if( x(k) <= xi(ke) ) then
            call splint (xi, yi, y2, x(k), y(k), sequential = .true.)
          else
            y(k) = yi(ke)
          end if
        end do
        c% rh(k0:ke) = y(k0:ke)
      end if

      if(lqcl) then
        yi = c% qcl
        call init_splinex(xi,yi,y2)
        do k=k0,ke
          if( x(k) <= xi(ke) ) then
            call splint (xi, yi, y2, x(k), y(k), sequential = .true.)
          else
            y(k) = yi(ke)
          end if
        end do
        c% qcl(k0:ke) = y(k0:ke)
      end if

      if(lqci) then
        yi = c% qci
        call init_splinex(xi,yi,y2)
        do k=k0,ke
          if( x(k) <= xi(ke) ) then
            call splint (xi, yi, y2, x(k), y(k), sequential = .true.)
          else
            y(k) = yi(ke)
          end if
        end do
        c% qci(k0:ke) = y(k0:ke)
      end if

      !-----------------------------
      ! save new vertical coordinate
      !-----------------------------
      if(lgeo) then
        c% geo = - x
      else
        c% p = exp(x)
      end if
    end do

    deallocate(xi)
    deallocate(yi)
    deallocate(y2)
    deallocate(x )
    deallocate(y )

  end subroutine vert_intp_cols
!==============================================================================
  subroutine derive_tinv_sngl_col(col, tinv, t_lapse)
  !--------------------------------------------------------------
  ! Derive "near-surface temperature inversion" as departure from
  ! a constant lapse rate profile, whose structure shall be
  ! preserved! The lapse rate is derived from background T at
  ! 500m and 1000m above the surface, consistent with ICON's
  ! temperature_intp.
  ! Geopotential on full levels in col is required!
  !--------------------------------------------------------------
  type(t_col),           intent(in)  :: col           ! column
  real(wp), allocatable, intent(out) :: tinv(:)       ! temperature inversion
  real(wp)             , intent(out) :: t_lapse       ! lapse rate

    !----------------
    ! local variables
    !----------------
    integer                 :: ke                       ! number of levels
    integer                 :: i                        ! level below zpbl1
    integer                 :: k                        ! loop counter
    real(wp), pointer       :: t  (:)                   ! pointer to temperature
    real(wp), pointer       :: geo(:)                   ! pointer to geopotential
    real(wp)                :: gke                      ! surface geopotential
    real(wp)                :: t_zpbl1                  ! temperature at zpbl1
    real(wp)                :: t_zpbl2                  ! temperature at zpbl2
    real(wp)                :: zpbl1 =  500._wp * gacc  ! AGL heights used in
    real(wp)                :: zpbl2 = 1000._wp * gacc  ! lapse rate computations

    if(.not.associated(col% geo)) call finish('derive_tinv_sngl_col','geopotential missing')
    if(.not.associated(col% t))   call finish('derive_tinv_sngl_col','temperature missing')

    ke = size(col% t)

    allocate( tinv(ke) )
    tinv =  0._wp
    t    => col% t
    geo  => col% geo
    gke  =  col% s% geosp

    !------------------
    ! derive lapse rate
    !------------------
    i = 0
    do k=ke-1,1,-1
      if( geo(k) - gke >= zpbl1 ) then
        t_zpbl1 = t(k) + (zpbl1 - (geo(k)- gke)) * (t(k) - t(k+1)) / (geo(k) - geo(k+1))
        i = k + 1
        exit
      end if
    end do
    if ( i == 0 ) call finish ('derive_tinv_sngl_col','zpbl1 not found')
    t_zpbl2 = 0._wp
    do k=i-1,1,-1
      if( geo(k) - gke >= zpbl2 ) then
        t_zpbl2 = t(k) + (zpbl2 - (geo(k)- gke)) * (t(k) - t(k+1)) / (geo(k) - geo(k+1))
        exit
      end if
    end do
    if ( t_zpbl2 == 0._wp ) call finish ('derive_tinv_sngl_col','zpbl2 not found')
    t_lapse = (t_zpbl2 - t_zpbl1) / (zpbl2 - zpbl1)

    !---------------------------------------------------------------
    ! calculate deviation/departure from constant lapse rate profile
    !---------------------------------------------------------------
    do k=i,ke
      tinv(k) = t(k) - (t_zpbl1 + t_lapse * ( geo(k) - gke - zpbl1 ))
    end do

  end subroutine derive_tinv_sngl_col
!==============================================================================
  subroutine derive_tinv_cols(cols, tinv, t_lapse)
  !--------------------------------------------------------------
  ! currently not used !
  !--------------------------------------------------------------
  ! Derive "near-surface temperature inversion" as departure from
  ! a constant lapse rate profile, whose structure shall be
  ! preserved! The lapse rate is derived from background T at
  ! 500m and 1000m above the surface, consistent with ICON's
  ! temperature_intp.
  ! Geopotential on full levels in cols is required!
  !--------------------------------------------------------------
  type(t_cols),          intent(in)  :: cols          ! columns
  real(wp), allocatable, intent(out) :: tinv(:,:)     ! temperature inversion
  real(wp), allocatable, intent(out) :: t_lapse(:)    ! lapse rate

    !----------------
    ! local variables
    !----------------
    integer                 :: ke                       ! number of levels
    integer                 :: ncol                     ! number of columns
    integer                 :: i                        ! level below zpbl1
    integer                 :: j, k                     ! loop counter
    real(wp), pointer       :: t  (:)                   ! pointer to temperature
    real(wp), pointer       :: geo(:)                   ! pointer to geopotential
    real(wp)                :: gke                      ! surface geopotential
    real(wp)                :: t_zpbl1                  ! temperature at zpbl1
    real(wp)                :: t_zpbl2                  ! temperature at zpbl2
    real(wp)                :: zpbl1 =  500._wp * gacc  ! AGL heights used in
    real(wp)                :: zpbl2 = 1000._wp * gacc  ! lapse rate computations

    ke   = cols% ke
    ncol = cols% ncol

    if( ncol == 0 ) return
    if(.not.associated(cols% col(1)% geo)) call finish('derive_tinv_cols','geopotential missing')
    if(.not.associated(cols% col(1)% t))   call finish('derive_tinv_cols','temperature missing')

    allocate( t_lapse(ncol) )
    allocate( tinv(ke, ncol) )
    tinv = 0._wp

    do j=1,ncol
      t   => cols% col(j)% t
      geo => cols% col(j)% geo
      gke =  cols% col(j)% s% geosp

      !------------------
      ! derive lapse rate
      !------------------
      i = 0
      do k=ke-1,1,-1
        if( geo(k) - gke >= zpbl1 ) then
          t_zpbl1 = t(k) + (zpbl1 - (geo(k)- gke)) * (t(k) - t(k+1)) / (geo(k) - geo(k+1))
          i = k + 1
          exit
        end if
      end do
      if ( i == 0 ) call finish ('derive_tinv_cols','zpbl1 not found')
      t_zpbl2 = 0._wp
      do k=i-1,1,-1
        if( geo(k) - gke >= zpbl2 ) then
          t_zpbl2 = t(k) + (zpbl2 - (geo(k)- gke)) * (t(k) - t(k+1)) / (geo(k) - geo(k+1))
          exit
        end if
      end do
      if ( t_zpbl2 == 0._wp ) call finish ('derive_tinv_cols','zpbl2 not found')
      t_lapse(j) = (t_zpbl2 - t_zpbl1) / (zpbl2 - zpbl1)

      !---------------------------------------------------------------
      ! calculate deviation/departure from constant lapse rate profile
      !---------------------------------------------------------------
      do k=i,ke
        tinv(k,j) = t(k) - (t_zpbl1 + t_lapse(j) * ( geo(k) - gke - zpbl1 ))
      end do
    end do

    nullify(t  )
    nullify(geo)

  end subroutine derive_tinv_cols
!==============================================================================
  subroutine mix_q(a, h_mix_q)
  !-------------------------------------------------------
  ! derive humidity on lower levels from relative humidity
  ! and on higher levels from specific humidity
  !
  !  h_mix_q < 0: use only spec. hum.
  !  h_mix_q = 0: use only rel.  hum.
  !  h_mix_q > 0: mix hum. representations
  !-------------------------------------------------------
  type(t_atm)      ,intent(inout)                 :: a
  real(wp)         ,intent(in)                    :: h_mix_q

    !----------------
    ! local variables
    !----------------
    integer                 :: ke                       ! number of levels
    integer                 :: ktr                      ! transition level
    integer                 :: k                        ! loop counter
    real(wp)                :: hhl_top                  ! model top
    real(wp)                :: hhlmin, hhlmax           ! max/min hhl on each level

    !--------------------------------------
    ! check for presence of required fields
    !--------------------------------------
    if(.not.associated(a% q )) call finish('mix_q','q  is not associated')
    if(.not.associated(a% t )) call finish('mix_q','t  is not associated')
    if(.not.associated(a% rh)) call finish('mix_q','rh is not associated')
    if(.not.associated(a% pf)) call finish('mix_q','pf is not associated')

    ke = a% ub(3)

    if ( h_mix_q < 0._wp ) then
       a% rh = rh_q (a% q, a% t, a% pf)
      return
    end if

    !--------
    ! set ktr
    !--------
    ktr = 1
    if ( h_mix_q > 0._wp ) then
      select case (a% grid% vct)
      case (VCT_Z_HYB, VCT_Z_GEN)
        if( .not. associated(a% grid% hhl) ) call finish('mix_q','hhl not associated')

        hhl_top = maxval ( a% grid% hhl(:,:,1,:) )
        hhl_top = p_max (hhl_top)
        if ( h_mix_q > hhl_top ) call finish('mix_q','h_mix_q too large')

        do k=2,ke
          hhlmin = minval ( a% grid% hhl(:,:,k,:) )
          hhlmin = p_min (hhlmin)
          hhlmax = maxval ( a% grid% hhl(:,:,k,:) )
          hhlmax = p_max (hhlmax)

          if (hhlmin /= hhlmax) call finish('mix_q','levels are not flat at h_mix_q')

          if ( h_mix_q > hhlmax ) then
            ktr = k
            if ( dace% lpio ) write(6,'(a,i0,/,a,i0,a,i0)') &
              ' derive humidity from specific humidity on levels  1:', ktr-1, &
              ' and from relative humidity on levels ', ktr, ':', ke
            exit
          end if
        end do
      case default
        call finish('mix_q','unsupported vertical coordinate')
      end select
    end if

    a% q(:,:,ktr:ke,:) = q_rh (a% rh(:,:,ktr:ke,:), a% t(:,:,ktr:ke,:), a% pf(:,:,ktr:ke,:))
    if (ktr /= 1) &
      a% rh(:,:,1:ktr-1,:) = rh_q (a% q(:,:,1:ktr-1,:), a% t(:,:,1:ktr-1,:), a% pf(:,:,1:ktr-1,:))

  end subroutine mix_q
!==============================================================================
  subroutine mix_q1(states, h_mix_q)
  !-------------------------------------------------------
  ! derive humidity on lower levels from relative humidity
  ! and on higher levels from specific humidity
  !
  !  h_mix_q < 0: use only spec. hum.
  !  h_mix_q = 0: use only rel.  hum.
  !  h_mix_q > 0: mix hum. representations
  !-------------------------------------------------------
  type(t_atm)      ,intent(inout)                 :: states(:)
  real(wp)         ,intent(in)                    :: h_mix_q

    !----------------
    ! local variables
    !----------------
    integer                 :: i                        ! loop over ensemble members
    integer                 :: ke                       ! number of levels
    integer                 :: ktr                      ! transition level
    integer                 :: k                        ! loop counter
    real(wp)                :: hhl_top                  ! model top
    real(wp)                :: hhlmin, hhlmax           ! max/min hhl on each level

    if ( size(states) == 0 ) return
    !--------------------------------------
    ! check for presence of required fields
    !--------------------------------------
    if(.not.associated(states(1)% q )) call finish('mix_q1','q  is not associated')
    if(.not.associated(states(1)% t )) call finish('mix_q1','t  is not associated')
    if(.not.associated(states(1)% rh)) call finish('mix_q1','rh is not associated')
    if(.not.associated(states(1)% pf)) call finish('mix_q1','pf is not associated')

    ke = states(1)% ub(3)

    if ( h_mix_q < 0._wp ) then
      do i=1,size (states)
        states(i)% rh = rh_q ( states(i)% q, states(i)% t, states(i)% pf )
      end do
      return
    end if

    !--------
    ! set ktr
    !--------
    ktr = 1
    if ( h_mix_q > 0._wp ) then
      select case (states(1)% grid% vct)
      case (VCT_Z_HYB, VCT_Z_GEN)
        if( .not. associated(states(1)% grid% hhl) ) call finish('mix_q1','hhl not associated')

        hhl_top = maxval ( states(1)% grid% hhl(:,:,1,:) )
        hhl_top = p_max (hhl_top)
        if ( h_mix_q > hhl_top ) call finish('mix_q1','h_mix_q too large')

        do k=2,ke
          hhlmin = minval ( states(1)% grid% hhl(:,:,k,:) )
          hhlmin = p_min (hhlmin)
          hhlmax = maxval ( states(1)% grid% hhl(:,:,k,:) )
          hhlmax = p_max (hhlmax)

          if (hhlmin /= hhlmax) call finish('mix_q1','levels are not flat at h_mix_q')

          if ( h_mix_q > hhlmax ) then
            ktr = k
            if ( dace% lpio ) write(6,'(a,i0,/,a,i0,a,i0)') &
              ' derive humidity from specific humidity on levels  1:', ktr-1, &
              ' and from relative humidity on levels ', ktr, ':', ke
            exit
          end if
        end do
      case default
        call finish('mix_q1','unsupported vertical coordinate')
      end select
    end if

    do i=1,size (states)
      states(i)% q(:,:,ktr:ke,:) = q_rh ( states(i)% rh(:,:,ktr:ke,:), &
        states(i)% t(:,:,ktr:ke,:), states(i)% pf(:,:,ktr:ke,:) )
    end do

    if (ktr /= 1) then
      do i=1,size (states)
        states(i)% rh(:,:,1:ktr-1,:) = rh_q ( states(i)% q(:,:,1:ktr-1,:), &
          states(i)% t(:,:,1:ktr-1,:), states(i)% pf(:,:,1:ktr-1,:) )
      end do
    end if

  end subroutine mix_q1
!==============================================================================
  subroutine interpolate_atm (y, x, params, vertical)
  !------------------------------
  ! interpolate atmospheric state
  !------------------------------
  type (t_atm)    ,intent(inout)        :: y        ! destination
  type (t_atm)    ,intent(inout)        :: x        ! source
  character(len=*),intent(in)           :: params   ! parameters to interpolate
  character       ,intent(in) ,optional :: vertical ! y(es) n(o) p(res) g(eop)
    !----------------
    ! local variables
    !----------------
    integer               :: i,j,k,l,m,ic,n     ! indices
    integer               :: ly(4)              ! lower bounds destination grid
    integer               :: uy(4)              ! upper bounds destination grid
    type (t_mcols)        :: mc   (dace% npe)
    type (t_cols)         :: cols (dace% npe)
    type (t_icol),pointer :: icols     (:,:,:)
    type (t_icol),pointer :: ic_sea    (:,:,:)
    type (t_icol),pointer :: ic_land   (:,:,:)
    type (t_icol),pointer :: ic_soil   (:,:,:)
    type(t_h_idx),pointer :: idx       (:,:,:)
    type(t_h_idx),pointer :: idx_land  (:,:,:)
    type(t_h_idx),pointer :: idx_sea   (:,:,:)
    type(t_h_idx),pointer :: idx_soil  (:,:,:)
    type (t_cols)         :: d2       ! second derivative of spline
    logical  ,allocatable :: lini(:)  ! spline coefficient initialisation flag
    type (t_col) ,pointer :: c, d     ! temporary
    type (t_sl)  ,pointer :: s
    real(wp)              :: z
    real(wp)              :: a, b
    logical               :: lvert    ! flag for vertical interpolation
    logical               :: lvertp   ! vertical interpolation in p
    logical               :: lvertz   ! vertical interpolation in z
    integer               :: np       ! number of parameters to interpolate
    integer    ,parameter :: nm = 16  ! max. number of parameters
    character(len=16)     :: pars(nm) ! names of variables to derive
    integer(i8)           :: iatm     ! field to gather
    logical               :: lt, lrh, lu, lv, linfl
    logical               :: lpp, lq, lw, lqcl, lqci
    logical               :: ldpsdt, lpf
    logical               :: lwso
    type(t_grid) ,pointer :: gx, gy          ! pointers to grid info
    integer               :: n_land          !       land points
    integer               :: n_snow          !       points with snow
    integer               :: n_no_snow       !       points without snow
    logical               :: lsst            ! interpolate SST analysis results
    logical               :: lsnow           !            SNOW analysis results
    logical               :: lsma            !             SMA analysis results
    real(wp)     ,pointer :: h (:)           ! actual vertical coordinate
    logical               :: lghx, lgfx      ! geopotential height present in x
    logical               :: lghy, lgfy      ! geopotential height present in y
    logical               :: lphx, lpfx      ! pressure present in x
    logical               :: lphy, lpfy      ! pressure present in y
    integer               :: nzx             ! # of vertical gridpoints in x
    integer    ,parameter :: st_nosoil(4)   &! no soil moisture content
                             = [ST_ICE,ST_ROCK,ST_SEAWATER,ST_SEAICE]
    !--------------------------------------
    ! Surface derivatives for extrapolation
    ! Temperature: dT/d(-geo) = gamma/g,     (gamma = clim. lapse rate),
    !              dT/dlog(p) = dT/d(-geo) * d(-geo)/dlog(p)|_{p_0}
    !                         = gamma/g    * R*Tv_0
    !                        ~= gamma/g*R * T_0
    !--------------------------------------
    real(wp),   parameter :: dtz = lapse_cl / gacc
    real(wp),   parameter :: dtp = dtz * R
    real(wp), allocatable :: df1(:)          ! surface derivative of T
    real(wp), allocatable :: df2(:)          ! surface derivative of ln(p)
!   real(wp),   parameter :: Tx  = 290.5_wp  ! Fallback (virt.) temperature
    !-----------------------------------------
    ! set up vertical interpolation parameters
    !-----------------------------------------
    nzx    = x% grid% nz
    iatm   = 0_i8
    lvertp = .false.
    lvertz = .false.
    lvert  = .true.;
    if (present (vertical)) then
      lvert = .false.
      select case (vertical)
      case ('y')
        lvert  = .true.
      case ('n')
      case ('p')
        lvertp = .true.
      case ('g')
        lvertz = .true.
      end select
    endif
    if (lvert) then
      select case (y% grid% vct)
      case (VCT_P_ISO, VCT_P_HYB)
        lvertp = .true.
      case default
        lvertz = .true.
      end select
    endif
    if (lvertz) then
      iatm  = COL_GEO + COL_GEOH
      lvert = .true.
      lghx  = associated (x% geoh)
      lgfx  = associated (x% geof)
      lghy  = associated (y% geoh)
      lgfy  = associated (y% geof)
      if (.not. lgfx) call set_geo (x, geof=.true.)
      if (.not. lgfy) call set_geo (y, geof=.true.)
    endif
    if (lvertp) then
      iatm  = COL_P
      lvert = .true.
      lphx  = associated (x% ph)
      lpfx  = associated (x% pf)
      lphy  = associated (y% ph)
      lpfy  = associated (y% pf)
      if (.not. lpfx) call set_p (x)
      if (.not. lpfy) call set_p (y)
    endif
    !--------------------------------
    ! determine fields to interpolate
    !--------------------------------
    call split (pars, params, np)
    call v3d_to_dwd (pars (1:np))
    if (np <  0) call finish ('interpolate_atm','nm too small')
    if (np == 0) return
!   call message('interpolate_atm','====================================')
!   call message('interpolate_atm','check for lower boundary condition')
!   call message('interpolate_atm','check for consistency with vert_intp')
!   call message('interpolate_atm','====================================')
    !---------------------------------------------------------
    ! special keywords to handle SST/SNOW/SMA analysis results
    !---------------------------------------------------------
    lsst  = any (pars(1:np) == 'SST_ANA' )
    lsnow = any (pars(1:np) == 'SNOW_ANA')
    lsma  = any (pars(1:np) == 'SMA_ANA' )

    lt    =  any (pars == 'T'       )
    lrh   =  any (pars == 'RELHUM'  )
    lu    =  any (pars == 'U'       )
    lv    =  any (pars == 'V'       )
    lpf   =  any (pars == 'P'       )
    lpp   =  any (pars == 'PP'      )
    lw    =  any (pars == 'W'       )
    lq    =  any (pars == 'QV'      )
    lqcl  =  any (pars == 'QC'      )
    lqci  =  any (pars == 'QI'      )
    linfl =  any (pars == 'f_inflat')
    ldpsdt = any (pars == 'DPSDT'   )
    !--------------------------------------------------------
    ! We need to choose between external and LETKF-based SMA!
    ! If provided, an external analysis shall take precedence
    !--------------------------------------------------------
    lwso  =  any (pars == 'W_SO'    ) .and. .not. lsma

    !-----------------------------------------------------
    ! set key for atmospheric fields to get from other PEs
    !-----------------------------------------------------
    if (lt)    iatm = iatm + COL_T
    if (lrh)   iatm = iatm + COL_RH
    if (lu )   iatm = ior (iatm, COL_UV)
    if (lv )   iatm = ior (iatm, COL_UV)
    if (lpf)   iatm = ior (iatm, COL_P + COL_T)
    if (lpp)   iatm = iatm + COL_PP
    if (lw)    iatm = iatm + COL_W
    if (lq)    iatm = iatm + COL_Q
    if (lqcl)  iatm = iatm + COL_QCL
    if (lqci)  iatm = iatm + COL_QCI
    if (linfl) iatm = iatm + COL_INFLAT
    if (lwso)  iatm = iatm + COL_WSO
!   if (lwso)  call finish ("interpolate_atm","COL_WSO not fully implemented")

    !---------------------------------------------------
    ! set up model indices and weights for interpolation
    ! (obsolete, managed by 'hor_intp_coef_atm') below
    !---------------------------------------------------
!   call set_mcols_m (mc, icols, y% grid, x% grid, iatm)

    nullify (idx)
    nullify (idx_sea)
    nullify (idx_land)
    nullify (idx_soil)
    nullify (ic_sea)
    nullify (ic_land)
    nullify (ic_soil)

    !---------------------------------------------------
    ! set up model indices and weights for interpolation
    !---------------------------------------------------
    call hor_intp_coef_atm                                    &
      (x% grid, y% grid, lsma .or. lwso, lsnow, lsst, idx,    &
       ind_soil=idx_soil, ind_land=idx_land, ind_sea=idx_sea, &
       fr_lake=[.05_wp,.05_wp], fr_land=[.2_wp,.05_wp],       &
                            fr_land_sst=[.01_wp,.5_wp]        )

    call mcol_idx   (mc, icols,   idx,      x% grid, iatm)
    if (lsst) &
      call mcol_idx (mc, ic_sea,  idx_sea,  x% grid, 0_i8)
    if (lsnow) &
      call mcol_idx (mc, ic_land, idx_land, x% grid, 0_i8)
    if (lsma) &
      call mcol_idx (mc, ic_soil, idx_soil, x% grid, COL_SOIL)
    if (lwso .and. (.not. lsma)) &
      call mcol_idx (mc, ic_soil, idx_soil, x% grid, COL_WSO)

    !------------------------------------------
    ! gather columns required for interpolation
    !------------------------------------------
    call get_cols    (mc, x, cols, iatm=iatm)

    !------------
    ! interpolate
    !------------
    n_land          = 0
    n_snow          = 0
    n_no_snow       = 0

    ly =  y% lb
    uy =  y% ub
    gx => x% grid
    gy => y% grid

    do m=1,size(cols)
      if (cols(m)% ncol == 0) cycle
      !---------------------
      ! allocate temporaries
      !---------------------
      if (lvert) then
        call alloc_cols (d2, tmp=cols(m))
        allocate (lini (cols(m)% ncol))
        lini = .false.
        allocate (df1  (cols(m)% ncol))
        allocate (df2  (cols(m)% ncol))
        df1  = HUGE (0._wp)
        df2  = HUGE (0._wp)
      endif
      !-----------------------------
      ! check for presence of fields
      !-----------------------------
      if (lsst) then
        !--------------------------
        ! fields from SST analysis:
        !--------------------------
        call check_alloc (gx% lsm,     'SST_ANA: x% FR_LAND missing')
        call check_alloc (gy% lsm,     'SST_ANA: y% FR_LAND missing')
        call check_alloc ( x% t_so,    'SST_ANA: x% T_SO    missing')
        call check_alloc ( y% t_so,    'SST_ANA: y% T_SO    missing')
        call check_alloc ( x% tsurf,   'SST_ANA: x% T_G     missing')
        call check_alloc ( y% tsurf,   'SST_ANA: y% T_G     missing')
        call check_alloc ( x% t_ice,   'SST_ANA: x% T_ICE   missing')
        call check_alloc ( y% t_ice,   'SST_ANA: y% T_ICE   missing')
        call check_alloc ( x% h_ice,   'SST_ANA: x% H_ICE   missing')
        call check_alloc ( y% h_ice,   'SST_ANA: y% H_ICE   missing')
        call check_alloc ( x% fr_ice,  'SST_ANA: x% FR_ICE  missing')
        call check_alloc ( y% fr_ice,  'SST_ANA: y% FR_ICE  missing')
      endif
      if (lsnow) then
        !---------------------------
        ! fields from SNOW analysis:
        !---------------------------
        call check_alloc ( x% h_snow,   'SNOW_ANA: x% H_SNOW   missing')
        call check_alloc ( y% h_snow,   'SNOW_ANA: y% H_SNOW   missing')
        call check_alloc ( x% t_snow,   'SNOW_ANA: x% T_SNOW   missing')
        call check_alloc ( y% t_snow,   'SNOW_ANA: y% T_SNOW   missing')
        call check_alloc ( x% w_snow,   'SNOW_ANA: x% W_SNOW   missing')
        call check_alloc ( y% w_snow,   'SNOW_ANA: y% W_SNOW   missing')
        call check_alloc ( x% w_i,      'SNOW_ANA: x% W_I      missing')
        call check_alloc ( y% w_i,      'SNOW_ANA: y% W_I      missing')
        call check_alloc ( x% freshsnw, 'SNOW_ANA: x% FRESHSNW missing')
        call check_alloc ( y% freshsnw, 'SNOW_ANA: y% FRESHSNW missing')
        call check_alloc ( x% rho_snow, 'SNOW_ANA: x% RHO_SNOW missing')
        call check_alloc ( y% rho_snow, 'SNOW_ANA: y% RHO_SNOW missing')

!       call check_alloc ( x% t_so,     'SNOW_ANA: x% T_SO     missing')
!       call check_alloc ( y% t_so,     'SNOW_ANA: y% T_SO     missing')
!       call check_alloc ( x% tsurf,    'SNOW_ANA: x% T_G      missing')
!       call check_alloc ( y% tsurf,    'SNOW_ANA: y% T_G      missing')
      endif
      if (lsma) then
        !--------------------------
        ! fields from SMA analysis:
        !--------------------------
        call check_alloc ( x% t_so,    'SMA_ANA: x% T_SO    missing')
        call check_alloc ( y% t_so,    'SMA_ANA: y% T_SO    missing')
        call check_alloc ( x% w_so,    'SMA_ANA: x% W_SO    missing')
        call check_alloc ( y% w_so,    'SMA_ANA: y% W_SO    missing')
        call check_alloc (gx% soiltyp, 'SMA_ANA: x% soiltyp missing')
        call check_alloc (gy% soiltyp, 'SMA_ANA: y% soiltyp missing')
        call check_alloc (gx% lsm,     'SMA_ANA: x% FR_LAND missing')
        call check_alloc (gy% lsm,     'SMA_ANA: y% FR_LAND missing')
      endif
      !----------------------------------
      ! loop over destination grid-points
      !----------------------------------
      if (lsst .or. lsnow .or. lsma .or. lwso) then
        do l = ly(4), uy(4)
          do j = ly(2), uy(2)
src_grid:   do i = ly(1), uy(1)
              !-----------------------
              ! destination land point
              !-----------------------
              if (lsnow) then
                if (gy% lsm(i,j,1,l) >= 0.05_wp) then
                  n_land = n_land + 1
                  n = ic_land(i,j,l)% h% imc(1,1)
                  if (n==0) call finish('interpolate_atm','missing surf neighbour')
                  s => cols(m)% col(n)% s
                  !------------------
                  ! source land point
                  !------------------
!                 call intpl_snow (y, i,j,l, s% t_so_0, s% tsurf, &
                  call intpl_snow (y, i,j,l,                      &
                     s% h_snow,   s% t_snow, s% w_snow, s% w_i,   &
                     s% freshsnw, s% rho_snow                     )
                endif
              endif
              !----------------------
              ! destination sea point
              !----------------------
              if (lsst) then
                if (gy% lsm(i,j,1,l) < 0.5_wp) then
                  n = ic_sea(i,j,l)% h% imc(1,1)
                  if (n==0) call finish('interpolate_atm','missing sea neighbour')
                  s => cols(m)% col(n)% s
                  call intpl_sst (y, i, j, l,                         &
                    s% t_so_0, s% tsurf, s% t_ice, s% h_ice, s% fr_ice)
                endif
              endif
              !----------------------------
              ! destination land soil class
              !----------------------------
              if (lsma .or. lwso) then
                if (     gy% soiltyp (i,j,1,l) /= 9999._wp  ) then
                if (all (gy% soiltyp (i,j,1,l) /= st_nosoil)) then
                  n = ic_soil(i,j,l)% h% imc(1,1)
                  if (n==0) call finish('interpolate_atm','missing soil neighbour')
                  c => cols(m)% col(n)
                  if (lsma) then
                     call intpl_sm (y, i,j,l, c, x% grid)
                  elseif (lwso) then ! .not. lsma
                     call intpl_wso (y, i,j,l, c, x% grid)
                  endif
                endif
                endif
              endif
            end do src_grid
          end do
        end do
      endif

      !-------------------------
      ! loop over 'observations'
      !-------------------------
      do l = ly(4), uy(4)
        do j = ly(2), uy(2)
          do i = ly(1), uy(1)
            !-------------------------------
            ! set interpolation coefficients
            !-------------------------------
            if (lvert) then
              do ic = 1, size(icols(i,j,l)% h% imc,1)
                n = icols(i,j,l)% h% imc(ic,1)
                if (n==0) exit
                if (.not. lini (n)) then
                  lini (n) = .true.
                  c => cols(m)% col(n)
                  d => d2     % col(n)
                  !---------------------------------------------------------
                  ! Spline vertical coordinate must run from "low" to "high"
                  !---------------------------------------------------------
                  if (associated (c% p   )) c% p    = log (c% p)
                  if (associated (c% geo )) c% geo  =    - c% geo
                  if (associated (c% geoh)) c% geoh =    - c% geoh
                  if (lvertp) then
                    h => c% p
                    if (lt ) df1(n) = dtp   * c% t(nzx)  ! surface dT/dlog(p)
                    if (lpf) df2(n) = 1._wp              ! surface dlog(p)/dlog(p)
                  else
                    h => c% geo
                    if (lt ) df1(n) = dtz                ! surface dT/d(-geo)
                    if (lpf) df2(n) = (1/R) / c% t(nzx)  ! surface dlog(p)/d(-geo)
                  endif
                  if (lrh  ) call init_splinex (h,       c% rh,   d% rh  ,ypn=0._wp)
                  if (lu   ) call init_splinex (h,       c% u,    d% u   ,ypn=0._wp)
                  if (lv   ) call init_splinex (h,       c% v,    d% v   ,ypn=0._wp)
                  if (lq   ) call init_splinex (h,       c% q,    d% q   ,ypn=0._wp)
                  if (lqcl ) call init_splinex (h,       c% qcl,  d% qcl ,ypn=0._wp)
                  if (lqci ) call init_splinex (h,       c% qci,  d% qci ,ypn=0._wp)
                  if (linfl) call init_splinex (h,       c% infl, d% infl,ypn=0._wp)
                  if (lpp  ) call init_splinex (h,       c% pp,   d% pp  ,ypn=0._wp)
                  if (lt   ) call init_splinex (h,       c% t,    d% t   ,ypn=df1(n))
                  if (lpf  ) call init_splinex (h,       c% p,    d% p   ,ypn=df2(n))
                  if (lw   ) call init_splinex (c% geoh, c% w,    d% w   ,ypn=0._wp)
!                 if (associated(c% geoh)) then
!                   call init_splinex (c% ph(2:),c% geoh(2:),d% geoh(2:))
!                   d% geoh(1) = 0.
!                 endif
                endif
              end do
            endif
            !------------------------------
            ! loop over levels, interpolate
            !------------------------------
            do k = 1, y% grid% nz + 1
              !----------------
              ! w (half levels)
              !----------------
              if (lw) then
                if (lvertz) z = - y% geoh(i,j,k,l)
                b = 0._wp
                do ic = 1, size(icols(i,j,l)% h% imc,1)
                  n = icols(i,j,l)% h% imc(ic,1)
                  if (n==0) exit
                  c => cols(m)% col(n)
                  if (lvert) then
                    h => c% geoh
                    d => d2     % col(n)
                    call spline_ (h   , & ! <-- Argument grid
                                  c% w, & ! <-- Gridded function
                                  d% w, & ! <-- 2nd derivative of spline
                                  z,    & ! <-- Interpolation point
                                  a,    & ! --> Interpolated function value
                                  0._wp ) ! <-- surface derivative
                    b = b + icols(i,j,l)% h% w(ic) * a
                  else
                    b = b + icols(i,j,l)% h% w(ic) * c% w (k)
                  endif
                end do
                y% w(i,j,k,l) = b
              endif
              !------------
              ! full levels
              !------------
              if (k == y% grid% nz + 1) cycle
              if (lvertp) z = log (y% pf  (i,j,k,l))
              if (lvertz) z =    - y% geof(i,j,k,l)
              !------------
              ! temperature
              !------------
              if (lt) then
                b = 0._wp
                do ic = 1, size(icols(i,j,l)% h% imc,1)
                  n = icols(i,j,l)% h% imc(ic,1)
                  if (n==0) exit
                  c => cols(m)% col(n)
                  if (lvert) then
                    if (lvertp) then
                      h => c% p
                    else
                      h => c% geo
                    endif
                    d => d2   % col(n)
                    call spline_ (h,    & ! <-- Argument grid
                                  c% t, & ! <-- Gridded function
                                  d% t, & ! <-- 2nd derivative of spline
                                  z,    & ! <-- Interpolation point
                                  a,    & ! --> Interpolated function value
                                  df1(n)) ! <-- surface derivative
                    b = b + icols(i,j,l)% h% w(ic) * a
                  else
                    b = b + icols(i,j,l)% h% w(ic) * c% t (k)
                  endif
                end do
                y% t(i,j,k,l) = b
              endif
              !------------------
              ! specific humidity
              !------------------
              if (lq) then
                b = 0._wp
                do ic = 1, size(icols(i,j,l)% h% imc,1)
                  n = icols(i,j,l)% h% imc(ic,1)
                  if (n==0) exit
                  c => cols(m)% col(n)
                  if (lvert) then
                    if (lvertp) then
                      h => c% p
                    else
                      h => c% geo
                    endif
                    d => d2   % col(n)
                    call spline_ (h,    & ! <-- Argument grid
                                  c% q, & ! <-- Gridded function
                                  d% q, & ! <-- 2nd derivative of spline
                                  z,    & ! <-- Interpolation point
                                  a,    & ! --> Interpolated function value
                                  0._wp ) ! <-- surface derivative
                    b = b + icols(i,j,l)% h% w(ic) * a
                  else
                    b = b + icols(i,j,l)% h% w(ic) * c% q (k)
                  endif
                end do
                y% q(i,j,k,l) = b
              endif
              !------------
              ! cloud water
              !------------
              if (lqcl) then
                b = 0._wp
                do ic = 1, size(icols(i,j,l)% h% imc,1)
                  n = icols(i,j,l)% h% imc(ic,1)
                  if (n==0) exit
                  c => cols(m)% col(n)
                  if (lvert) then
                    if (lvertp) then
                      h => c% p
                    else
                      h => c% geo
                    endif
                    d => d2   % col(n)
                    call spline_ (h,      & ! <-- Argument grid
                                  c% qcl, & ! <-- Gridded function
                                  d% qcl, & ! <-- 2nd derivative of spline
                                  z,      & ! <-- Interpolation point
                                  a,      & ! --> Interpolated function value
                                  0._wp   ) ! <-- surface derivative
                    b = b + icols(i,j,l)% h% w(ic) * a
                  else
                    b = b + icols(i,j,l)% h% w(ic) * c% qcl (k)
                  endif
                end do
                y% qcl(i,j,k,l) = b
              endif
              !------------------
              ! specific humidity
              !------------------
              if (lqci) then
                b = 0._wp
                do ic = 1, size(icols(i,j,l)% h% imc,1)
                  n = icols(i,j,l)% h% imc(ic,1)
                  if (n==0) exit
                  c => cols(m)% col(n)
                  if (lvert) then
                    if (lvertp) then
                      h => c% p
                    else
                      h => c% geo
                    endif
                    d => d2   % col(n)
                    call spline_ (h,      & ! <-- Argument grid
                                  c% qci, & ! <-- Gridded function
                                  d% qci, & ! <-- 2nd derivative of spline
                                  z,      & ! <-- Interpolation point
                                  a,      & ! --> Interpolated function value
                                  0._wp   ) ! <-- surface derivative
                    b = b + icols(i,j,l)% h% w(ic) * a
                  else
                    b = b + icols(i,j,l)% h% w(ic) * c% qci (k)
                  endif
                end do
                y% qci(i,j,k,l) = b
              endif
              !------------------
              ! relative humidity
              !------------------
              if (lrh) then
                b = 0._wp
                do ic = 1, size(icols(i,j,l)% h% imc,1)
                  n = icols(i,j,l)% h% imc(ic,1)
                  if (n==0) exit
                  c => cols(m)% col(n)
                  if (lvert) then
                   if (lvertp) then
                     h => c% p
                   else
                     h => c% geo
                   endif
                   d => d2     % col(n)
                   call spline_ (h,     & ! <-- Argument grid
                                 c% rh, & ! <-- Gridded function
                                 d% rh, & ! <-- 2nd derivative of spline
                                 z,     & ! <-- Interpolation point
                                 a,     & ! --> Interpolated function value
                                 0._wp  ) ! <-- surface derivative
                    b = b + icols(i,j,l)% h% w(ic) * a
                  else
                    b = b + icols(i,j,l)% h% w(ic) * c% rh (k)
                  endif
                end do
                y% rh(i,j,k,l) = b
              endif
              !--
              ! u
              !--
              if (lu) then
                b = 0._wp
                do ic = 1, size(icols(i,j,l)% h% imc,1)
                  n = icols(i,j,l)% h% imc(ic,1)
                  if (n==0) exit
                  c => cols(m)% col(n)
                  if (lvert) then
                    if (lvertp) then
                      h => c% p
                    else
                      h => c% geo
                    endif
                    d => d2     % col(n)
                    call spline_ (h,    & ! <-- Argument grid
                                  c% u, & ! <-- Gridded function
                                  d% u, & ! <-- 2nd derivative of spline
                                  z,    & ! <-- Interpolation point
                                  a,    & ! --> Interpolated function value
                                  0._wp ) ! <-- surface derivative
                    b = b + icols(i,j,l)% h% w(ic) * a
                  else
                    b = b + icols(i,j,l)% h% w(ic) * c% u (k)
                  endif
                end do
                y% u(i,j,k,l) = b
              endif
              !--
              ! v
              !--
              if (lv) then
                b = 0._wp
                do ic = 1, size(icols(i,j,l)% h% imc,1)
                  n = icols(i,j,l)% h% imc(ic,1)
                  if (n==0) exit
                  c => cols(m)% col(n)
                  if (lvert) then
                    if (lvertp) then
                      h => c% p
                    else
                      h => c% geo
                    endif
                    d => d2     % col(n)
                    call spline_ (h   , & ! <-- Argument grid
                                  c% v, & ! <-- Gridded function
                                  d% v, & ! <-- 2nd derivative of spline
                                  z,    & ! <-- Interpolation point
                                  a,    & ! --> Interpolated function value
                                  0._wp ) ! <-- surface derivative
                    b = b + icols(i,j,l)% h% w(ic) * a
                  else
                    b = b + icols(i,j,l)% h% w(ic) * c% v (k)
                  endif
                end do
                y% v(i,j,k,l) = b
              endif
              !---
              ! pf (full-level pressure)
              !---
              if (lpf) then
                b = 0._wp
                do ic = 1, size(icols(i,j,l)% h% imc,1)
                  n = icols(i,j,l)% h% imc(ic,1)
                  if (n==0) exit
                  c => cols(m)% col(n)
                  if (lvert) then
                    if (lvertp) then
                      h => c% p
                    else
                      h => c% geo
                    endif
                    d => d2     % col(n)
                    call spline_ (h   , & ! <-- Argument grid
                                  c% p, & ! <-- Gridded function
                                  d% p, & ! <-- 2nd derivative of spline
                                  z,    & ! <-- Interpolation point
                                  a,    & ! --> Interpolated function value
                                  df2(n)) ! <-- surface derivative
                    b = b + icols(i,j,l)% h% w(ic) * a
                  else
                    b = b + icols(i,j,l)% h% w(ic) * c% p (k)
                  endif
                end do
                y% pf(i,j,k,l) = exp (b)  ! Invert log(p) transformation
              endif
              !---
              ! pp
              !---
              if (lpp) then
                b = 0._wp
                do ic = 1, size(icols(i,j,l)% h% imc,1)
                  n = icols(i,j,l)% h% imc(ic,1)
                  if (n==0) exit
                  c => cols(m)% col(n)
                  if (lvert) then
                    if (lvertp) then
                      h => c% p
                    else
                      h => c% geo
                    endif
                    d => d2     % col(n)
                    call spline_ (h    , & ! <-- Argument grid
                                  c% pp, & ! <-- Gridded function
                                  d% pp, & ! <-- 2nd derivative of spline
                                  z,     & ! <-- Interpolation point
                                  a,     & ! --> Interpolated function value
                                  0._wp  ) ! <-- surface derivative
                    b = b + icols(i,j,l)% h% w(ic) * a
                  else
                    b = b + icols(i,j,l)% h% w(ic) * c% pp (k)
                  endif
                end do
                y% pp(i,j,k,l) = b
              endif
              !-----------------
              ! inflation factor
              !-----------------
              if (linfl) then
                b = 0._wp
                do ic = 1, size(icols(i,j,l)% h% imc,1)
                  n = icols(i,j,l)% h% imc(ic,1)
                  if (n==0) exit
                  c => cols(m)% col(n)
                  if (lvert) then
                    if (lvertp) then
                      h => c% p
                    else
                      h => c% geo
                    endif
                   d => d2     % col(n)
                   call spline_ (h,       & ! <-- Argument grid
                                 c% infl, & ! <-- Gridded function
                                 d% infl, & ! <-- 2nd derivative of spline
                                 z,       & ! <-- Interpolation point
                                 a,       & ! --> Interpolated function value
                                 0._wp    ) ! <-- surface derivative
                   if (lvertz) then         ! no extrapolation
                     if (z > h (1)  ) a = c% infl (1)
                     if (z < h (nzx)) a = c% infl (nzx)
                   endif
                    b = b + icols(i,j,l)% h% w(ic) * a
                  else
                    b = b + icols(i,j,l)% h% w(ic) * c% infl (k)
                  endif
                end do
                y% f_inflat(i,j,k,l) = b
              endif
            end do ! k (levels)
            !---
            ! ps
            !---
            if (associated (y% ps)) then
              b = 0._wp
              do ic = 1, size(icols(i,j,l)% h% imc,1)
                n = icols(i,j,l)% h% imc(ic,1)
                if (n==0) exit
                b = b + icols(i,j,l)% h% w(ic) * cols(m)% col(n)% s% ps
              end do
              y% ps(i,j,1,l) = b
            endif
            !------
            ! dpsdt
            !------
            if (ldpsdt) then
              if (associated (y% dpsdt)) then
                b = 0._wp
                do ic = 1, size(icols(i,j,l)% h% imc,1)
                  n = icols(i,j,l)% h% imc(ic,1)
                  if (n==0) exit
                  b = b + icols(i,j,l)% h% w(ic) * cols(m)% col(n)% s% dpsdt
                end do
                y% dpsdt(i,j,1,l) = b
              endif
            endif
!           !----
!           ! psr
!           !----
!           if (associated (y% psr)) then
!             b = 0._wp
!             do ic = 1, size(icols(i,j,l)% h% imc,1)
!               n = icols(i,j,l)% h% imc(ic,1)
!               if (n==0) exit
!               b = b + icols(i,j,l)% h% w(ic) * cols(m)% col(n)% s% psr
!             end do
!             y% psr(i,j,1,l) = b
!           endif
          end do
        end do
      end do
      !---------------------
      ! deallocate temporaries
      !---------------------
      if (lvert) then
        deallocate (lini)
        call dealloc_cols (d2)
      endif
    end do

    !---------------------
    ! print out statistics
    !---------------------
    if (lsst .or. lsnow .or. lsma) then
      n_land          = p_sum (n_land)
      n_snow          = p_sum (n_snow)
      n_no_snow       = p_sum (n_no_snow)
      if (dace% lpio) then
        write (6,*)
        write (6,*)'  interpolate_atm:'
        write (6,*)
        write (6,*)'    params = ',trim(params)
        write (6,*)'    lsst   = ',lsst
        write (6,*)'    lsnow  = ',lsnow
        write (6,*)'    lsma   = ',lsma
        write (6,*)
        write (6,*) n_snow         ,'out of',n_land,'points with    snow.'
        write (6,*) n_no_snow      ,'out of',n_land,'points without snow.'
        write (6,*)
      endif
    endif
    !--------
    ! cleanup
    !--------
    deallocate        (icols)
    call dealloc_cols (cols)
    call destruct     (mc)
    deallocate        (idx)

    if(associated (idx_sea))  deallocate (idx_sea)
    if(associated (idx_land)) deallocate (idx_land)
    if(associated (idx_soil)) deallocate (idx_soil)
    if(associated (ic_sea))   deallocate (ic_sea)
    if(associated (ic_land))  deallocate (ic_land)
    if(associated (ic_soil))  deallocate (ic_soil)

    if (lvertz) then
      if (.not. lghx) call deallocate (x,'geoh')
      if (.not. lgfx) call deallocate (x,'geof')
      if (.not. lghy) call deallocate (y,'geoh')
      if (.not. lgfy) call deallocate (y,'geof')
    endif
    if (lvertp) then
      if (.not. lphx) call deallocate (x,'ph')
      if (.not. lpfx) call deallocate (x,'pf')
      if (.not. lphy) call deallocate (y,'ph')
      if (.not. lpfy) call deallocate (y,'pf')
    endif
  contains
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  subroutine spline_ (xa, ya, y2a, x, y, df)
    real(wp), intent(in)  :: xa(:), ya(:), y2a(:)
    real(wp), intent(in)  :: x             ! Abscissa
    real(wp), intent(out) :: y             ! Interpol. function value
    real(wp), intent(in)  :: df            ! Surface derivative
    !--------------------------------------------------------
    ! Wrapper to splint handling abscissa values out of range
    !--------------------------------------------------------
    integer :: nz
    nz = ubound (xa, 1)
    if (x < xa(nz)) then
       call splint (xa, ya, y2a, x, y)
    else
       y = ya(nz) + (x - xa(nz)) * df
    end if
  end subroutine spline_
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  subroutine intpl_snow (y, i,j,d,                                      &
                         h_snow, t_snow, w_snow, w_i, freshsnw, rho_snow)
  !--------------------------------------
  ! interpolate fields from SNOW analysis
  !--------------------------------------
  type (t_atm) ,intent(inout) :: y        ! destination state
  integer      ,intent(in)    :: i, j, d  ! destination grid indices
! real(wp)     ,intent(in)    :: t_so     ! soil surface temperature
! real(wp)     ,intent(in)    :: t_g      ! temperature of surface
  real(wp)     ,intent(in)    :: h_snow   ! snow height
  real(wp)     ,intent(in)    :: t_snow   ! snow temperature
  real(wp)     ,intent(in)    :: w_snow   ! snow water equivalent
  real(wp)     ,intent(in)    :: w_i      ! water in interception reservoir
  real(wp)     ,intent(in)    :: freshsnw ! snow age
  real(wp)     ,intent(in)    :: rho_snow ! snow density
    !----------------------------------------
    ! take parameters from closest land point
    !----------------------------------------
!print *,'# SNOW #',t_so, t_g, h_snow, t_snow, w_snow, w_i, freshsnw, rho_snow
    y% h_snow   (i,j,1,d) = h_snow
    y% t_snow   (i,j,1,d) = t_snow
    y% w_snow   (i,j,1,d) = w_snow
    y% freshsnw (i,j,1,d) = freshsnw
    y% rho_snow (i,j,1,d) = rho_snow
    if (h_snow > 0._wp) then
      !--------------------------------
      ! we have snow on the source grid
      !--------------------------------
      n_snow = n_snow + 1
      y% w_i   (i,j,1,d) = 0._wp ! ? water in interception reservoir
!     y% t_so  (i,j,0,d) = t_so  ! soil surface temperature
!     y% tsurf (i,j,1,d) = t_g   ! surface temperature to the atmosphere
    else
      !-----------------------------------
      ! we have no snow on the source grid
      !-----------------------------------
      n_no_snow = n_no_snow + 1
    endif
  end subroutine intpl_snow
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  subroutine intpl_sm (y, i,j,d, c, sg)
  !----------------------------
  ! interpolate soil quantities
  !----------------------------
  type (t_atm) ,intent(inout) :: y        ! destination state
  integer      ,intent(in)    :: i, j, d  ! destination grid indices
  type(t_col)  ,intent(in)    :: c        ! source
  type(t_grid) ,intent(in)    :: sg       ! source grid

    real(wp) :: w (sg% ns)
    integer  :: std, sts
    !----------------------------------------
    ! take parameters from closest land point
    !----------------------------------------
    std = int (y% grid% soiltyp (   i,   j, 1,    d))
    sts = int (sg%      soiltyp (c% i, c%j, 1, c% l))
    y%   t_so (i,j,1:,d) = c% t_so
    if (std == sts) then
      y% w_so (i,j, :,d) = c% w_so
    else
      select case (std)
      case (ST_ICE, ST_ROCK, ST_SEAWATER, ST_SEAICE)
      case default
        select case (sts)
        case (ST_ICE, ST_ROCK, ST_SEAWATER, ST_SEAICE)
        case default
          w = (c% w_so / st - cl(sts)) / (cu(sts) - cl(sts))
          w = min (1._wp, max (0._wp, w))
          y% w_so (i,j, :,d) = st * (w * (cu(std) - cl(std)) + cl(std))
        end select
      end select
    endif
  end subroutine intpl_sm
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  subroutine intpl_wso (y, i,j,d, c, sg)
    !----------------------------
    ! interpolate soil quantities
    !----------------------------
    type(t_atm)  ,intent(inout) :: y        ! destination state
    integer      ,intent(in)    :: i, j, d  ! destination grid indices
    type(t_col)  ,intent(in)    :: c        ! source
    type(t_grid) ,intent(in)    :: sg       ! source grid

    real(wp) :: w (sg% ns)
    integer  :: std, sts
    !----------------------------------------
    ! take parameters from closest land point
    !----------------------------------------
    std = int (y% grid% soiltyp (   i,   j, 1,    d))
    sts = int (sg%      soiltyp (c% i, c%j, 1, c% l))
    if (std == sts) then
      y% w_so (i,j, :,d) = c% w_so
    else
      select case (std)
      case (ST_ICE, ST_ROCK, ST_SEAWATER, ST_SEAICE)
      case default
        select case (sts)
        case (ST_ICE, ST_ROCK, ST_SEAWATER, ST_SEAICE)
        case default
          w = (c% w_so / st - ca(sts)) / (cp(sts) - ca(sts))
          y% w_so (i,j, :,d) = st * (w * (cp(std) - ca(std)) + ca(std))
        end select
      end select
    endif
  end subroutine intpl_wso
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  subroutine intpl_sst (y, i,j,d, t_so, t_g, t_ice, h_ice, fr_ice)
  !-------------------------------------
  ! interpolate fields from SST analysis
  !-------------------------------------
  type (t_atm) ,intent(inout) :: y       ! destination state
  integer      ,intent(in)    :: i, j, d ! destination grid indices
  real(wp)     ,intent(in)    :: t_so    ! sea surface temperature
  real(wp)     ,intent(in)    :: t_g     ! temperature of surface to atmosphere
  real(wp)     ,intent(in)    :: t_ice   ! ice temperature
  real(wp)     ,intent(in)    :: h_ice   ! ice height
  real(wp)     ,intent(in)    :: fr_ice  ! ice fraction
!print *,'# SST  #',t_so, t_g, t_ice, h_ice, fr_ice
    y% t_so   (i,j,0,d) = t_so
    y% tsurf  (i,j,1,d) = t_g
    y% t_ice  (i,j,1,d) = t_ice
    y% h_ice  (i,j,1,d) = h_ice
    y% fr_ice (i,j,1,d) = fr_ice
  end subroutine intpl_sst
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  subroutine intpl_sst_land (y, i,j,d, t_so, t_g)
  !---------------------------------------------------
  ! interpolate fields from SST analysis
  ! no land sea found in the vicinity of a sea point
  ! do something reasonable, +++ needs replacement +++
  !---------------------------------------------------
  type (t_atm) ,intent(inout) :: y       ! destination state
  integer      ,intent(in)    :: i, j, d ! destination grid indices
  real(wp)     ,intent(in)    :: t_so    ! sea surface temperature
  real(wp)     ,intent(in)    :: t_g     ! temperature of surface to atmosphere
    y% t_so   (i,j,0,d) = t_g
    y% tsurf  (i,j,1,d) = t_g
    y% t_ice  (i,j,1,d) = t_g
    y% h_ice  (i,j,1,d) = 0
    y% fr_ice (i,j,1,d) = 0
  end subroutine intpl_sst_land
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  subroutine check_alloc (field, text)
  !----------------------------------------------
  ! checks allocation status of atmospheric field
  ! aborts if not associated
  !----------------------------------------------
  real(wp)         ,pointer    :: field (:,:,:,:)
  character(len=*) ,intent(in) :: text
    if (.not.associated (field)) call finish ('interpolate_atm',text)
  end subroutine check_alloc
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  end subroutine interpolate_atm
!==============================================================================
  subroutine construct_atm_set (s, size)
    type(t_atm_set), intent(inout) :: s
    integer,         intent(in)    :: size
    allocate (s% state(size))
    allocate (s% fname(size)); s% fname = ""
    allocate (s% gpid (size)); s% gpid  = -1
    allocate (s% ptype(size)); s% ptype = -1
    s% size = size
    s% used = 0
  end subroutine construct_atm_set
  !----------------------------------------------------------------------------
  subroutine allocate_patm_0 (s, patm, fname, gpid, ptype)
    type(t_atm_set),  target,   intent(inout) :: s
    type(t_atm),      pointer,  intent(out)   :: patm
    character(len=*), optional, intent(in)    :: fname
    integer,          optional, intent(in)    :: gpid
    integer,          optional, intent(in)    :: ptype
    integer :: k
    k = s% used + 1
    if (k > s% size) call finish ("allocate_patm_0","atm_set too small")
    patm => s% state(k)
    if (present (fname)) then
       if (len_trim (fname) > len (s% fname)) then
          write(0,*) "fname too long: " // trim (fname)
          call finish ("allocate_patm_0","fname too long")
       end if
       s% fname(k) = fname
    end if
    if (present (gpid) ) s% gpid (k) = gpid
    if (present (ptype)) s% ptype(k) = ptype
    s% used = k
  end subroutine allocate_patm_0
  !----------------------------------------------------------------------------
  subroutine allocate_patm_1 (s, patm, n, fname, gpid, ptype)
    type(t_atm_set),  target,   intent(inout) :: s
    type(t_atm),      pointer,  intent(out)   :: patm(:)
    integer,                    intent(in)    :: n
    character(len=*), optional, intent(in)    :: fname(:)
    integer,          optional, intent(in)    :: gpid
    integer,          optional, intent(in)    :: ptype
    integer :: k, j
    j = s% used + 1
    k = s% used + n
    if (k > s% size) call finish ("allocate_patm_1","atm_set too small")
    patm(1:n) => s% state(j:k)
    if (present (fname)) then
       if (size (fname) /= n) call finish ("allocate_patm_1","bad size: fname")
       if (any (len_trim (fname) > len (s% fname))) then
          write(0,*) "fname too long: max.len=", maxval (len_trim (fname))
          call finish ("allocate_patm_1","fname too long")
       end if
       s% fname(j:k) = fname(:)
    end if
    if (present (gpid)) then
       s% gpid (j:k) = gpid
    end if
    if (present (ptype)) then
       s% ptype(j:k) = ptype
    end if
    s% used = k
  end subroutine allocate_patm_1
  !----------------------------------------------------------------------------
  subroutine allocate_patm_2 (s, patm, n1, n2, fname, gpid, ptype)
    type(t_atm_set),  target,   intent(inout) :: s
    type(t_atm),      pointer,  intent(out)   :: patm(:,:)
    integer,                    intent(in)    :: n1, n2
    character(len=*), optional, intent(in)    :: fname(:,:)
    integer,          optional, intent(in)    :: gpid
    integer,          optional, intent(in)    :: ptype
    integer :: k, j, n
    j = s% used + 1
    n = n1*n2
    k = s% used + n
    if (k > s% size) call finish ("allocate_patm_2","atm_set too small")
    patm(1:n1,1:n2) => s% state(j:k)
    if (present (fname)) then
       if (size (fname) /= n) call finish ("allocate_patm_2","bad size: fname")
       if (any (len_trim (fname) > len (s% fname))) then
          write(0,*) "fname too long: max.len=", maxval (len_trim (fname))
          call finish ("allocate_patm_2","fname too long")
       end if
       s% fname(j:k) = reshape (fname, [n])
    end if
    if (present (gpid)) then
       s% gpid (j:k) = gpid
    end if
    if (present (ptype)) then
       s% ptype(j:k) = ptype
    end if
    s% used = k
  end subroutine allocate_patm_2
  !----------------------------------------------------------------------------
  subroutine destruct_atm_set (s)
    type(t_atm_set), intent(inout) :: s
    if (s% used > 0) then
!!       call deallocate (s% state)
    end if
    deallocate (s% state)
    deallocate (s% fname)
    deallocate (s% gpid )
    deallocate (s% ptype)
    s% size = 0
    s% used = 0
  end subroutine destruct_atm_set
  !----------------------------------------------------------------------------
  subroutine deallocate_atm_set (s)
    type(t_atm_set), intent(inout) :: s
    if (s% used > 0) then
!!       call deallocate (s% state)
    end if
  end subroutine deallocate_atm_set
  !----------------------------------------------------------------------------
  subroutine write_atm_set (s, grid, parallel, ref, mask, check_grid, verbose)
    type(t_atm_set),        intent(inout) :: s
    character(*), optional, intent(in)    :: grid       ! grid-param. to write
    logical,      optional, intent(in)    :: parallel
    logical,      optional, intent(in)    :: ref
    logical,      optional, intent(in)    :: mask(:)    ! subset to write
    logical,      optional, intent(in)    :: check_grid ! check for ident.grids
    logical,      optional, intent(in)    :: verbose    ! verbose diagnostics

    integer              :: i, j, k, m
    logical              :: lpar
    logical              :: lgrd
    logical              :: reordered
    logical, allocatable :: lmask(:)

    lpar = .false.; if (present (parallel)  ) lpar = parallel
    lgrd = .true. ; if (present (check_grid)) lgrd = check_grid

    k = s% used
    if (k == 0) return

    allocate (lmask(k))
    if (present (mask)) then
       if (size (mask) < k) call finish ("write_atm_set", "invalid mask")
       lmask = mask(1:k)
       if (size (mask) > k) then
          if (any (mask(k+1:))) call finish ("write_atm_set", "invalid mask")
       end if
       m = count (lmask)
    else
       lmask = .true.
       m = k
    end if
    if (m == 0) return

    i = 0
    do j = 1, k
       if (.not. lmask(j)) cycle
       if (i == 0) i = j
       if (s% fname(j) == "") then
          call finish ("write_atm_set","file name(s) must be set")
       end if
       if (s% gpid (j) < 0) then
          call finish ("write_atm_set","gen. process identifier(s) must be set")
       end if
       if (s% ptype(j) < 0) then
          call finish ("write_atm_set","process type(s) must be set")
       end if
    end do

    !-------------------------
    ! Require consistent grids
    !-------------------------
    do j = i+1, k
       if (.not. lmask(j)) cycle
       if (.not. same_horizontal_grid (s% state(i)% grid, s% state(j)% grid, &
                                       reordered=reordered      )) then
          if (lgrd) then
             write(0,*) "inconsistent horizontal grids:", i, " /=", j
             call finish ("write_atm_set","inconsistent horizontal grids")
          end if
          if (dace% lpio) print *, "note: different horizontal grids:",i,"/=",j
       end if
       if (.not. same_vertical_grid (s% state(i)% grid, s% state(j)% grid)) then
          if (lgrd) then
             write(0,*) "inconsistent vertical grids:", i, " /=", j
             call finish ("write_atm_set","inconsistent vertical grids")
          end if
          if (dace% lpio) print *, "note: different vertical grids:",i,"/=",j
       end if
    end do

    if (lpar) then
       call write_grib (atm=s% state(i:k), file =s% fname(i:k),            &
                        pid=s% gpid (i:k), ptype=s% ptype(i:k), grid=grid, &
                        mask=lmask(i:k), mode="w", ref=ref, verbose=verbose)
    else
       ! Serialized version
       do j = 1, k
          if (.not. lmask(j)) cycle
          call set_defaults (process = s% gpid(j), process_type = s% ptype(j))
          call write_grib   (s% state(j), file=s% fname(j), &
                             mode="w", grid=grid, ref=ref)
       end do
    end if

  end subroutine write_atm_set
!==============================================================================
  subroutine read_nml_extana
  !=============================
  ! read namelist group /EXTANA/
  !=============================
    integer :: ierr
#if defined(__ibm__)
    integer :: ios
#endif

    !--------------
    ! read namelist
    !--------------
    if (dace% lpio) then
      call position_nml ('EXTANA', status=ierr)
      select case (ierr)
      case (POSITIONED)
#if defined(__ibm__)
        read (nnml ,nml=EXTANA, iostat=ios)
        if (ios/=0) call finish ('read_nml_extana','ERROR in namelist /EXTANA/')
#else
        read (nnml ,nml=EXTANA)
#endif
      case default
        write(6,'(/,a,/)') '  External analysis namelist /EXTANA/ not present'
      end select
    end if

    call p_bcast (ierr, dace% pio)
    if (ierr /= POSITIONED) return

    if (dace% lpio) then
      !--------------------------------------
      ! convert shortnames to DACE convention
      !--------------------------------------
      call dwd_to_v3d (par_extana)
      call dwd_to_v3d (opt_extana)

      if (ext_snow) call eval_string (par_extana,'+ SNOW_ANA')
      if (ext_sst)  call eval_string (par_extana,'+ SST_ANA')
      if (ext_sma)  call eval_string (par_extana,'+ SMA_ANA')

      if ( index(method,'LETKF') /= 0 ) then
        if (ext_snow) call eval_string (read_extana,'+ h_snow t_snow w_snow w_i freshsnw rho_snow')
        if (ext_sst)  call eval_string (read_extana,'+ t_so tsurf t_ice h_ice fr_ice')
      end if
      if (ext_sma) call eval_string (read_extana, '+ t_so w_so')

      ext_sma_weight = max( min(1._wp, ext_sma_weight), 0._wp)
      extana_file = path_file (input, extana_file)

      !---------
      ! printout
      !---------
      write(6,'(a)') repeat('_',79)
      write(6,'()')
      write(6,'(a)')      '  External analysis namelist /EXTANA/ :'
      write(6,'()')
      write(6,'(a,a)')    '    extana_file     = ',trim(extana_file)
      write(6,'(a,l1)')   '    ext_snow        = ',ext_snow
      write(6,'(a,l1)')   '    ext_sst         = ',ext_sst
      write(6,'(a,l1)')   '    ext_sma         = ',ext_sma
      write(6,'(a,a)')    '    par_extana      = ',trim(par_extana)
      write(6,'(a,a)')    '    read_extana     = ',trim(read_extana)
      write(6,'(a,a)')    '    opt_extana      = ',trim(opt_extana)
      write(6,'(a,f10.5)')'    ext_sma_weight  = ',ext_sma_weight
      write(6,'(a,i5)')   '    ext_sma_top     = ',ext_sma_top
      write(6,'(a,i5)')   '    ext_sma_bot     = ',ext_sma_bot
      write(6,'()')
    end if

    !-----------------------------
    ! broadcast namelist variables
    !-----------------------------
    call p_bcast (extana_file           ,dace% pio)
    call p_bcast (ext_snow              ,dace% pio)
    call p_bcast (ext_sst               ,dace% pio)
    call p_bcast (ext_sma               ,dace% pio)
    call p_bcast (par_extana            ,dace% pio)
    call p_bcast (read_extana           ,dace% pio)
    call p_bcast (opt_extana            ,dace% pio)
    call p_bcast (ext_sma_weight        ,dace% pio)
    call p_bcast (ext_sma_top           ,dace% pio)
    call p_bcast (ext_sma_bot           ,dace% pio)

    !------------------------
    ! check for allowed input
    !------------------------
    if( ext_sma_top <    1 ) call finish('read_nml_extana','ext_sma_top < 1')
    if( ext_sma_bot > n_sl ) call finish('read_nml_extana','ext_sma_bot > n_sl')
    if( ext_sma_top > ext_sma_bot ) &
        call finish('read_nml_extana','ext_sma_top > ext_sma_bot')

  end subroutine read_nml_extana
!==============================================================================
  subroutine prep_extana (extana, det_fc, extana_er, enkf_grid, inv)
  !---------------------------------------------------
  ! read and process (deterministic) external analysis
  !---------------------------------------------------
  type (t_atm)       ,pointer               :: extana        ! external analysis
  type (t_atm)       ,intent(in)            :: det_fc        ! deterministic forecast
  type (t_atm)       ,pointer    ,optional  :: extana_er     ! external analysis (ens.resolution)
  type (t_grid)      ,pointer    ,optional  :: enkf_grid     ! ensemble grid
  type (t_inventory) ,pointer    ,optional  :: inv(:)        ! GRIB file inventory

    !----------------
    ! Local variables
    !----------------
    type (t_inventory), pointer     :: invt(:)               ! GRIB file inventory
    real(wp)          ,allocatable  :: tmp_sst_extana(:,:,:) ! temporary storage

    nullify (extana)
    if ( present(extana_er) ) nullify (extana_er)
    if ( extana_file == "" ) return

    !-------------------------------------
    ! read external deterministic analysis
    !-------------------------------------

    nullify (invt)
    if (present(inv)) then
      if (.not.associated(inv)) call finish ('prep_extana','inventory not associated')
      if ( size(inv) == 0 )     call finish ('prep_extana','inventory is empty')
      invt => inv
    else
      call get_inventory (invt, extana_file)
    end if
    if (dace% lpio) then
      write(6,'()')
      call print_inventory (invt, first=.true.)
    endif

    allocate       (extana)
    call construct (extana, det_fc% grid)

    ! external sst contains only one level, expect only t_so(0)
    if (ext_sst) extana% m(T_SO)% i% ub(3) = extana% m(T_SO)% i% lb(3)

    !-----------------------
    ! read external analysis
    !-----------------------
    call read      (extana, extana_file, invt,     &
                    time      = ana_time,          &
                    runtype   = "analysis",        &
                    fields    = trim (read_extana),&
                    optionals = trim (opt_extana)  )
    call print     (extana, grid=.false., comment='external analysis')
    if (present(inv)) then
      nullify (invt)
    else
      if (associated(invt)) deallocate (invt)
    end if

    !----------------------
    ! restore shape of T_SO
    !----------------------
    if (ext_sst) then
      allocate( tmp_sst_extana (                                       &
                  extana% m(T_SO)% i% lb(1):extana% m(T_SO)% i% ub(1), &
                  extana% m(T_SO)% i% lb(2):extana% m(T_SO)% i% ub(2), &
                  extana% m(T_SO)% i% lb(4):extana% m(T_SO)% i% ub(4) ))
      tmp_sst_extana = extana% m(T_SO)% ptr(:,:,extana% m(T_SO)% i% lb(3),:)
      call deallocate (extana, 't_so')
      extana% m(T_SO)% i% ub(3) = extana% grid% ns
      call allocate (extana, 't_so')
      extana% m(T_SO)% ptr = -huge(1._wp)
      extana% m(T_SO)% ptr(:,:,extana% m(T_SO)% i% lb(3),:) = tmp_sst_extana
      deallocate( tmp_sst_extana )
    end if

    !------------
    ! adjust W_SO ( and maybe T_SO )
    !------------
    if (ext_sma) then
      if(.not.associated(det_fc% w_so)) call finish('prep_extana','W_SO missing in first guess')

      extana% w_so = ext_sma_weight * extana% w_so + (1._wp - ext_sma_weight) * det_fc% w_so
      if( ext_sma_top /= 1 ) &
        extana% w_so(:,:,1:ext_sma_top-1,:) = det_fc% w_so(:,:,1:ext_sma_top-1,:)
      if( ext_sma_bot /= n_sl ) &
        extana% w_so(:,:,ext_sma_bot+1:n_sl,:) = det_fc% w_so(:,:,ext_sma_bot+1:n_sl,:)
      !--------------------------------------------------------------
      ! take T_SO from first guess if extana contains no SST analysis
      !--------------------------------------------------------------
      if (ext_sst .and. associated(det_fc% t_so)) extana% t_so = det_fc% t_so
    end if

    if (ext_sst .or. ext_sma) call print (extana, grid=.false., comment='external analysis (adjusted)')

    !-----------------------------------------------------
    ! interpolate external analysis to ensemble resolution
    !-----------------------------------------------------
    if ( present(extana_er) ) then
      if ( .not. present(enkf_grid) ) call finish('prep_extana','missing enkf_grid')
      !---------------------------------------------------------
      ! Check for presence of mandatory fields for interpolation
      !---------------------------------------------------------
      if (ext_sst .and. .not. associated (extana% tsurf)) then
        if (associated (det_fc% tsurf)) then
          call allocate (extana, 'tsurf')
          extana% tsurf = det_fc% tsurf
        else
          call finish ('prep_extana','ext_sst, but tsurf not available')
        end if
      end if
      allocate             (extana_er)
      call construct       (extana_er, enkf_grid, alloc = read_extana,&
                                                    name='extana_er'  )
      extana_er = -huge(1._wp)
      call interpolate_atm (extana_er, extana, par_extana, vertical='n')
!     call print           (extana_er, grid=.false.,                        &
!                           comment='external analysis, ensemble resolution')
    end if

  end subroutine prep_extana
!==============================================================================
end module mo_atm_util
