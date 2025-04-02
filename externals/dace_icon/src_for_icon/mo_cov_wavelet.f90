!
!+ Utilities for reading and writing wavelet covariance matrices
!
MODULE mo_cov_wavelet
!
! Description:
!   Utilities for reading and writing wavelet covariance matrices
!   ("B" and their 'relatives', e.g. "L") as NetCDF files
!
! Current Code Owner: DWD, Harald Anlauf
!    phone: +49 69 8062 4941
!    fax:   +49 69 8062 3721
!    email: harald.anlauf@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_1         2008/11/05 Andreas Rhodin
!  First operational 3D-Var release
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Andreas Rhodin
!  cleanup (remove unused variables)
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Harald Anlauf  DWD  2007-2008
!-----------------------------------------------------------------------------

  !-------------
  ! Modules used
  !-------------
  use mo_kind,       only: wp
  use mo_constants,  only: PI
  use mo_exception,  only: finish
  use mo_ascii,      only: NUL              ! achar (00)
  use mo_dace_string,only: tolower,        &! Lowercase of string
                           split            ! Split text string into array

  use mo_mpi_dace,   only: p_bcast,        &! overloaded MPI broadcast routine
                           p_send,         &! overloaded MPI send routine
                           p_recv,         &! overloaded MPI receive routine
                           dace             ! MPI group info

  use mo_dec_matrix, only: t_matrix_block, &
                           allocate_block, &
                           construct,      &
                           destruct,       &
                           p_bcast,        &! overloaded MPI broadcast routine
                           p_send,         &! overloaded MPI send routine
                           p_recv,         &! overloaded MPI receive routine
                           CSR, CSC,       &! Compressed sparse row/column
                           srep,           &! Mnemonic for storage repr.
                           operator(*)      ! matrix vector multiply

  use netcdf,                              &! netcdf90
                     only: nf90_create,    &! create NetCDF file
                           nf90_open,      &! open   NetCDF file
                           nf90_close,     &! close  NetCDF file
                           nf90_get_att,   &! get    attribute
                           nf90_put_att,   &! define attribute
                           nf90_def_dim,   &! define dimension
                           nf90_inq_dimid, &! inquire dimension id
                           nf90_inq_varid, &! inquire variable  id
                           nf90_put_var,   &! write variable
                           nf90_enddef,    &! end of definition block
                           nf90_redef,     &! reenter define mode
                           nf90_abort,     &! abort   define mode
                           NF90_GLOBAL,    &! global attribute id
                           NF90_CLOBBER,   &! flag: overwrite
                           NF90_NOWRITE,   &! flag: readonly
                           NF90_WRITE,     &! flag: read/write(?)
                           NF90_NOERR,     &! status: no error
                           NF90_INT,       &! data type id
                           NF90_FLOAT       ! data type id

  use mo_t_netcdf,   only:                 &! NetCDF utilities
                                            ! Module variables:
                           ncid,           &! NetCDF file id
                           vname,          &! variable name or other comment
                           xtype,          &! type  of variable to be defined
                           dimid,          &! dimid of variable to be defined
                           stanc,          &! NetCDF start  parameter
                           counc,          &! NetCDF count  parameter
                           strnc,          &! NetCDF stride parameter
                           rname,          &! called routine+netcdf routine
!                          dimid2,         &! dimension of var(2) to be def.
!                          dimid3,         &! dimension of var(3) to be def.
!                                           ! Module procedures:
                           chk,            &! checks error status
                           def_var,        &! calls nf90_def_var
!                          def_var_2,      &! as def_var but with dimid2(2)
!                          def_var_3,      &! as def_var but with dimid3(3)
                           def_str,        &! define character string variable
                           get_var,        &! get integer   variable
                           get_attr,       &! get attribute
                           get_dim          ! get dimension from NetCDF file

  use mo_wmo_tables, only: WMO6_GAUSSIAN    ! Gaussian latitude/longitude grid
  use mo_legendre,   only: gauaw            ! Gaussian latitudes and weights

  use mo_1dmra,      only: wave_1d,        &! 1d wavelet transform
                           wave_1dt,       &! 1d wavelet transform (transp.)
                           TR_SYN,         &! flag: wavelet synthesis
                           TR_ADJ,         &!       adjoint synthesis
                           WV_DAUB_8,      &! wavelets for testing
                           WV_CDV_8,       &!
                           WV_NONE,        &!
                           wv_name          ! wavelet names
  !============================================================================
  implicit none
  private
  !----------------
  ! Public entities
  !----------------
  public :: t_cov_wavelet             ! Wavelet covariance matrix
  public :: t_cov_meta                ! Covariance matrix metadata
  public :: construct                 ! Initialization of covariance matrix
  public :: destruct                  ! Deallocation of covariance matrix
  public :: apply_wv_cov_2d           ! Apply covariance matrix to 2d-field
  public :: create_cov_netcdf         ! Create empty NetCDF file with header
  public :: append_cov_netcdf         ! Append covariance matrix to NetCDF file
  public :: read_cov_netcdf_meta      ! Read covariance metadata from NetCDF
  public :: read_cov_netcdf           ! Read covariance matrix from NetCDF file
  public :: find_cov_netcdf           ! Check presence of field in NetCDF file
  public :: test_create_cov_netcdf    ! Minimal test of implementation
  public :: append_profile_netcdf     ! Append (meridional) profile to NetCDF
  public :: read_profile_list_netcdf  ! Read list of profiles from NetCDF file
  public :: read_profile_netcdf       ! Read (meridional) profile from NetCDF
  public :: find_profile_netcdf       ! Check presence of profile in NetCDF
  public :: merge_netcdf              ! Merge selected covariances/profiles
  public :: isvalidname               ! Check variable name for validity
  !
  ! Overloaded versions of MPI communication routines
  ! for the derived types t_cov_wavelet and t_cov_meta
  !
  public :: p_send
  public :: p_recv
  public :: p_bcast
  !============================================================================
  ! Covariance matrix metadata (grid, wavelet basis and representation)
  type t_cov_meta
     character(len=64) :: gridinfo        ! Human-readable grid description
     integer           :: gridtype        ! WMO table 6 (Gaussian or lat/lon)
     integer           :: dims(3)         ! Grid dimensions (grid points)
     integer           :: wv_basis(3)     ! Wavelet basis (code)
     character         :: representation  ! "B" or "L" (=sqrt(B))
  end type t_cov_meta
  !============================================================================
  ! Wavelet covariance matrix
  type t_cov_wavelet
     type(t_matrix_block) :: b          ! Sparse representation
     type(t_cov_meta)     :: meta       ! Grid and wavelet basis information
     character(len=16)    :: field      ! Field name
     integer              :: id         ! Field identification (user def.)
     logical              :: valid      ! Information is valid
  end type t_cov_wavelet
  !============================================================================
  interface construct
     module procedure construct_cov_wavelet
  end interface

  interface destruct
     module procedure destruct_cov_wavelet
  end interface

  interface p_send
    module procedure send_cov_meta
    module procedure send_cov_wavelet
  end interface p_send

  interface p_recv
    module procedure recv_cov_meta
    module procedure recv_cov_wavelet
  end interface p_recv

  interface p_bcast
    module procedure bcast_cov_meta
    module procedure bcast_cov_wavelet
  end interface p_bcast
  !============================================================================
  ! Module local stuff
  integer, parameter :: MAXMAT = 16            ! Max. covariance matrices/file
  integer, parameter :: namlen = len (wv_name) ! String length of wavelet names
  integer, parameter :: MAXPROF = 32           ! Max. number of profiles/file
  ! MPI tags
  integer, parameter :: TAG_COV_WAVELET = 1
  integer, parameter :: TAG_COV_META    = 2
  !============================================================================
contains
  !============================================================================
  subroutine construct_cov_wavelet (W)
    type(t_cov_wavelet), intent(inout) :: W
    call construct (W% B)
    W% meta% representation = "?"
    W% field = ""
    W% id    = -1
    W% valid = .false.
  end subroutine construct_cov_wavelet
  !============================================================================
  subroutine destruct_cov_wavelet (W)
    type(t_cov_wavelet), intent(inout) :: W
    call destruct (W% B)
    W% meta% representation = "?"
    W% field = ""
    W% id    = -1
    W% valid = .false.
  end subroutine destruct_cov_wavelet
  !============================================================================
  subroutine apply_wv_cov_2d (W, x, ntimes)
    !---------------------------------------------------
    ! Apply wavelet covariance matrix to 2d grid, ntimes
    !---------------------------------------------------
    type(t_cov_wavelet), intent(in)    :: W
    real(wp)           , intent(inout) :: x(:,:)
    integer            , intent(in)    :: ntimes
    !----------------
    ! Local variables
    !----------------
    integer  :: i, basis1, basis2
    real(wp) :: y(size (x))
    !-------------------------------
    ! Determine transformation plane
    !-------------------------------
    if (W% meta% dims(3) == 1) then
       basis1 = W% meta% wv_basis(1)    ! Transformation in x-y plane
       basis2 = W% meta% wv_basis(2)
    else if (W% meta% dims(2) == 1) then
       basis1 = W% meta% wv_basis(1)    ! Transformation in x-z plane
       basis2 = W% meta% wv_basis(3)
    else if (W% meta% dims(1) == 1) then
       basis1 = W% meta% wv_basis(2)    ! Transformation in y-z plane
       basis2 = W% meta% wv_basis(3)
    else
       write (0,*) "Grid dimensions:", W% meta% dims(:)
       call finish ("apply_wv_cov_2d","cannot handle 3d-grid")
    end if
    !---------------------------
    ! Apply B/L in wavelet space
    !---------------------------
    call wave_1d  (x, TR_ADJ, basis1)
    call wave_1dt (x, TR_ADJ, basis2)
    y = reshape (x, shape (y))
    do i = 1, ntimes
       y = W% b * y
    end do
    x = reshape (y, shape (x))
    call wave_1d  (x, TR_SYN, basis1)
    call wave_1dt (x, TR_SYN, basis2)
  end subroutine apply_wv_cov_2d
  !============================================================================
  subroutine create_cov_netcdf (filename, meta, description)
    !-------------------------------------------------------
    ! Create NetCDF file to hold wavelet covariance matrices
    ! and initialize with specified metadata
    !-----------------------------------------------------
    character(len=*),     intent(in)           :: filename
    type(t_cov_meta),     intent(in)           :: meta
    character(len=*),     intent(in), optional :: description
    !----------------
    ! Local variables
    !----------------
    integer,  parameter   :: ndim = 3
    integer               :: k, m, n, dimid_m1, dimid_n1, dimid_ndim, &
                             dimid_slen
    integer               :: basis_id, descr_id
    integer               :: lon_id, lat_id, pres_id, &
                             dimid_lon, dimid_lat, dimid_pres
    real(wp), allocatable :: ga(:), gw(:)
    real,     allocatable :: lat(:), lon(:), pres(:)
    character(len=namlen) :: wv_descr(ndim)     ! Wavelet basis (name)
    !----------------------
    ! Set up misc. metadata
    !----------------------
    allocate (lon(meta% dims(1)))       ! Grid information
    allocate (lat(meta% dims(2)))
    allocate (pres(meta% dims(3)))
    pres(:) = 500
    lon(:)  =   0 + (/ (k-1, k=1,meta% dims(1)) /) * 360./ meta% dims(1)
    if (meta% gridtype == WMO6_GAUSSIAN) then
       ! Gaussian grid
       allocate (ga(meta% dims(2)), gw(meta% dims(2)))
       call gauaw (ga, gw, meta% dims(2))
       lat(:)  = (-180.0_wp/PI) * asin (ga)
    else
       ! Regular grid
       lat(:)  = -90 + (/ (k-1, k=1,meta% dims(2)) /) * 180./(meta% dims(2)-1)
    end if
    wv_descr(:) = wv_name(meta% wv_basis(:))    ! Keep wavelet names for checks
    !----------
    ! Open file
    !----------
    rname = 'create_cov_netcdf: nf90_create'
    call chk (nf90_create (filename, NF90_CLOBBER, ncid))
    !------------------
    ! Define attributes
    !------------------
    rname = 'create_cov_netcdf: nf90_put_att'
    call chk (nf90_put_att (ncid, NF90_GLOBAL, 'Conventions', 'CF-1.0'))
    call chk (nf90_put_att (ncid, NF90_GLOBAL, 'gridinfo', &
                            meta% gridinfo))
    call chk (nf90_put_att (ncid, NF90_GLOBAL, 'WMO_gridtype', &
                            meta% gridtype))
    call chk (nf90_put_att (ncid, NF90_GLOBAL, 'cov_representation', &
                            meta% representation))
    if (present (description)) then
       call chk (nf90_put_att (ncid, NF90_GLOBAL, 'description', &
                               description))
    end if
    !---------------------------------------
    ! Initialize list of covariance matrices
    !---------------------------------------
    call chk (nf90_put_att (ncid, NF90_GLOBAL, 'field_list', ""))
    !----------------------------
    ! Initialize list of profiles
    !----------------------------
    call chk (nf90_put_att (ncid, NF90_GLOBAL, 'profile_list', ""))
    !------------------
    ! Define dimensions
    !------------------
    rname = 'create_cov_netcdf: nf90_def_dim'
    call chk (nf90_def_dim (ncid, 'lon'  ,meta% dims(1),dimid_lon))
    call chk (nf90_def_dim (ncid, 'lat'  ,meta% dims(2),dimid_lat))
    call chk (nf90_def_dim (ncid, 'pres' ,meta% dims(3),dimid_pres))
    call chk (nf90_def_dim (ncid, 'ndim'  ,ndim  ,dimid_ndim))
    call chk (nf90_def_dim (ncid, 'slen'  ,namlen,dimid_slen))
    m = product (meta% dims(:))
    n = m
    call chk (nf90_def_dim (ncid, 'm1'    ,m+1   ,dimid_m1 ))
    call chk (nf90_def_dim (ncid, 'n1'    ,n+1   ,dimid_n1 ))
    !-----------------
    ! Define variables
    !-----------------
    rname = 'create_cov_netcdf: def_var'

    xtype = NF90_FLOAT
    dimid = dimid_lon
    call def_var('lon',  lon_id, 'degrees_east', 'longitude')
    dimid = dimid_lat
    call def_var('lat',  lat_id, 'degrees_north','latitude')
    dimid = dimid_pres
    call def_var('pres', pres_id,'hPa',          'pressure')

    xtype = NF90_INT
    dimid = dimid_ndim
    call def_var('wv_basis',basis_id,'',         'wavelet basis (code)')
    call def_str('wv_descr',descr_id,dimid_slen, 'wavelet basis (name)')
    !------------------------
    ! End of definition block
    !------------------------
    rname = 'create_cov_netcdf: nf90_enddef'
    call chk(nf90_enddef(ncid))
    !----------------
    ! Store variables
    !----------------
    rname = 'create_cov_netcdf: nf90_put_var'
    call chk (nf90_put_var (ncid, lon_id ,  lon ))
    call chk (nf90_put_var (ncid, lat_id ,  lat ))
    call chk (nf90_put_var (ncid, pres_id,  pres))
    call chk (nf90_put_var (ncid, basis_id, meta% wv_basis))
    call chk (nf90_put_var (ncid, descr_id, wv_descr))
    !-----------
    ! Close file
    !-----------
    call chk (nf90_close (ncid))
  end subroutine create_cov_netcdf
  !============================================================================
  subroutine append_cov_netcdf (W, filename, field)
    type(t_cov_wavelet),  intent(in), target   :: W
    character(len=*),     intent(in)           :: filename
    character(len=*),     intent(in)           :: field
    !----------------
    ! Local variables
    !----------------
    integer            :: n, repr
    integer            :: dimid_m1, dimid_n1, dimid_ns, a_id, ia_id, ja_id
    type(t_cov_meta)   :: tmpmeta
    character(len=8)   :: mnemonic           ! Sparse storage scheme
    character(len=255) :: field_list
    character(len=16)  :: fields(MAXMAT)
    !-------------------------------------
    type(t_matrix_block),  pointer :: B
    type(t_cov_meta),      pointer :: meta
    !-------------------------------------
    B    => W% b
    meta => W% meta
    !----------------------------------
    ! Preparations & consistency checks
    !----------------------------------
    if (.not. isvalidname (field)) then
       call finish ('append_cov_netcdf','bad field name: '//trim (field))
    end if
    repr     = B% repr  ! Array storage format
    mnemonic = srep(repr)
    select case (repr)
    case (CSR,CSC)      ! OK
    case default
       write (0,*) "B% repr, mnemonic =", repr, mnemonic
       call finish ('append_cov_netcdf','unsupported matrix representation')
    end select
    !----------
    ! Open file
    !----------
    rname = 'append_cov_netcdf: nf90_open'
    call chk (nf90_open (filename, NF90_WRITE, ncid))
    !-------------------------------------
    ! Get metadata of existing NetCDF file
    !-------------------------------------
    call read_cov_netcdf_common (filename, tmpmeta, field_list)
    call split (fields, field_list, n)
    if (any (fields(1:n) == field)) then
       call chk (nf90_abort (ncid))
       call finish ('append_cov_netcdf','cannot overwrite: '//trim (field))
    end if
    if ( meta% gridtype         /= tmpmeta% gridtype       .or. &
         meta% representation   /= tmpmeta% representation .or. &
         any (meta% dims(:)     /= tmpmeta% dims(:))       .or. &
         any (meta% wv_basis(:) /= tmpmeta% wv_basis(:)) ) then
       call finish ('append_cov_netcdf','incompatible metadata')
    end if
    !----------------------------
    ! Inquire existing dimensions
    !----------------------------
    rname = 'append_cov_netcdf: nf90_inq_dimid'
    call chk (nf90_inq_dimid  (ncid, 'm1', dimid_m1))
    call chk (nf90_inq_dimid  (ncid, 'n1', dimid_n1))
    !------------------------
    ! Reenter definition mode
    !------------------------
    rname = 'append_cov_netcdf: nf90_redef'
    call chk (nf90_redef (ncid))
    !----------------------
    ! Define new dimensions
    !----------------------
    rname = 'append_cov_netcdf: nf90_def_dim'
    call chk (nf90_def_dim (ncid, trim (field)//'_ns', B% nonzero, dimid_ns))
    !---------------------
    ! Define new variables
    !---------------------
    rname = 'append_cov_netcdf: def_var'
    select case (repr)
    case (CSR); dimid = dimid_m1
    case (CSC); dimid = dimid_n1
    case default
       call finish ('append_cov_netcdf','unsupported matrix representation')
    end select

    xtype = NF90_INT
    call def_var(trim (field)//'_ia'    ,ia_id,'',trim (field)//' ia indices')

    dimid = dimid_ns
    call def_var(trim (field)//'_ja'    ,ja_id,'',trim (field)//' ja indices')

    xtype = NF90_FLOAT
    call def_var(trim (field)//'_packed', a_id,'',trim(field)//' coefficients')
    !----------------------
    ! Add/modify attributes
    !----------------------
    rname = 'append_cov_netcdf: nf90_put_att'
    call chk (nf90_put_att (ncid, NF90_GLOBAL, 'field_list', &
                            adjustl (trim (field_list)//" "//trim (field))))
    call chk (nf90_put_att (ncid, NF90_GLOBAL, &
                            trim (field)//'_storage_representation', &
                            mnemonic))
    !------------------------
    ! End of definition block
    !------------------------
    rname = 'append_cov_netcdf: nf90_enddef'
    call chk (nf90_enddef (ncid))
    !--------------------
    ! Store new variables
    !--------------------
    rname = 'append_cov_netcdf: nf90_put_var'
    call chk (nf90_put_var (ncid, ia_id ,B% ia     ))
    call chk (nf90_put_var (ncid, ja_id ,B% ja     ))
    call chk (nf90_put_var (ncid,  a_id ,B% packed ))
    !-----------
    ! Close file
    !-----------
    call chk (nf90_close (ncid))
  end subroutine append_cov_netcdf
  !============================================================================
  subroutine read_cov_netcdf_meta (filename, meta, field_list, description)
    !------------------------------------------
    ! Read covariance metadata from NetCDF file
    !------------------------------------------
    character(len=*),     intent(in)            :: filename
    type(t_cov_meta),     intent(out), optional :: meta
    character(len=*),     intent(out), optional :: field_list
    character(len=*),     intent(out), optional :: description
    !----------
    ! Open file
    !----------
    vname = filename
    rname = 'read_cov_netcdf_meta: nf90_open'
    call chk (nf90_open (filename, NF90_NOWRITE, ncid))
    !-------------
    ! Get metadata
    !-------------
    call read_cov_netcdf_common (filename, meta, field_list, description)
    !-----------
    ! Close file
    !-----------
    call chk (nf90_close (ncid))
  end subroutine read_cov_netcdf_meta
  !============================================================================
  subroutine read_cov_netcdf (W, filename, field)
    !----------------------------------------------------
    ! Read covariance matrix for 'field' from NetCDF file
    !----------------------------------------------------
    type(t_cov_wavelet),  intent(out), target   :: W
    character(len=*),     intent(in)            :: filename
    character(len=*),     intent(in)            :: field
    !----------------
    ! Local variables
    !----------------
    integer            :: k, m, n, ns, repr
    character(len=8)   :: mnemonic
    character(len=255) :: field_list
    character(len=16)  :: fields(MAXMAT)
    !-------------------------------------
    type(t_matrix_block),  pointer :: B
    type(t_cov_meta),      pointer :: meta
    !-------------------------------------
    B    => W% b
    meta => W% meta
    W% field = trim (field)
    !----------
    ! Open file
    !----------
    rname = 'read_cov_netcdf: nf90_open'
    call chk (nf90_open (filename, NF90_NOWRITE, ncid))
    !-------------
    ! Get metadata
    !-------------
    call read_cov_netcdf_common (filename, meta, field_list)
    call split (fields, field_list, n)
    if (n == 0) then
       call finish ('read_cov_netcdf','empty NetCDF file')
    end if
    if (all (fields(1:n) /= field)) then
       call finish ('read_cov_netcdf','not in NetCDF file: '//trim (field))
    end if
    !------------------------
    ! Get specific attributes
    !------------------------
    rname = 'read_cov_netcdf: get_attr'
    call get_attr (mnemonic, '', trim (field)//'_storage_representation')
    ! Search in list of known mnemonics
    repr = HUGE (0)
    do k = lbound (srep,1), ubound (srep,1)
       if (mnemonic == srep(k)) then
          repr = k
          exit
       end if
    end do
    select case (repr)
    case (CSR,CSC)      ! OK
    case default
       write (0,*) "Found: mnemonic, repr = ", mnemonic, repr
       call finish ('read_cov_netcdf','unsupported matrix representation')
    end select
    !------------------------
    ! Get specific dimensions
    !------------------------
    rname = 'read_cov_netcdf: get_dim'
    call get_dim (m,  'm1'); m = m - 1
    call get_dim (n,  'n1'); n = n - 1
    call get_dim (ns, trim (field)//'_ns')
    !--------------------------------------
    ! Allocate arrays for data and metadata
    !--------------------------------------
    call allocate_block (B, repr, m=m, n=n, ns=ns)
    !----------
    ! Read data
    !----------
    stanc  = 1
    strnc  = 1
    counc  = 0
    rname = 'read_cov_netcdf: get_var'
    call get_var (B% ia    ,trim (field)//'_ia')
    call get_var (B% ja    ,trim (field)//'_ja')
    call get_var (B% packed,trim (field)//'_packed')
    W% valid = .true.
    W% id    = 0
    !-----------
    ! Close file
    !-----------
    call chk (nf90_close (ncid))
  end subroutine read_cov_netcdf
  !============================================================================
  subroutine read_cov_netcdf_common (filename, meta, field_list, description)
    !--------------------------------------------------
    ! Read covariance metadata from current NetCDF file
    !--------------------------------------------------
    character(len=*),     intent(in)            :: filename
    type(t_cov_meta),     intent(out), optional :: meta
    character(len=*),     intent(out), optional :: field_list
    character(len=*),     intent(out), optional :: description
    !----------------
    ! Local variables
    !----------------
    integer               :: ndim
    type(t_cov_meta)      :: tmpmeta
    integer,               allocatable :: dims(:), wv_basis(:)
    real,                  allocatable :: lon(:), lat(:), pres(:)
    character(len=namlen), allocatable :: units(:), wv_descr(:)
    !---------------
    ! Get attributes
    !---------------
    call get_attr (tmpmeta% gridinfo,      '', 'gridinfo')
    call get_attr (tmpmeta% gridtype,      '', 'WMO_gridtype')
    call get_attr (tmpmeta% representation,'', 'cov_representation')
    if (present (description)) then
       call get_attr (description,         '', 'description')
    end if
    if (present (field_list)) then
       call get_attr (field_list,          '', 'field_list')
    end if
    !---------------
    ! Get dimensions
    !---------------
    rname = 'read_cov_netcdf_common: get_dim'
    call get_dim (ndim, 'ndim')
    if (ndim /= 3) then
       write (0,*) trim (filename), " : ndim =", ndim
       call finish ('read_cov_netcdf_common','unsupported data format')
    end if
    allocate (dims(ndim))
    allocate (units(ndim))
    dims  = 0
    units = ""
    if (ndim >= 1) then
       call get_dim  (dims(1),  'lon')
       call get_attr (units(1), 'lon',  'units')
    end if
    if (ndim >= 2) then
       call get_dim  (dims(2),  'lat')
       call get_attr (units(2), 'lat',  'units')
    end if
    if (ndim >= 3) then
       call get_dim  (dims(3),  'pres')
       call get_attr (units(3), 'pres', 'units')
    end if
    !--------------------------------------
    ! Allocate arrays for data and metadata
    !--------------------------------------
    if (ndim >= 1) allocate (lon(dims(1)))
    if (ndim >= 2) allocate (lat(dims(2)))
    if (ndim >= 3) allocate (pres(dims(3)))
    allocate (wv_basis(ndim))
    allocate (wv_descr(ndim))
    wv_basis = 0
    wv_descr = ""
    !--------------
    ! Read metadata
    !--------------
    stanc  = 1
    strnc  = 1
    counc  = 0
    rname = 'read_cov_netcdf_common: get_var'
    if (allocated (lon))  call get_var (lon,  "lon")
    if (allocated (lat))  call get_var (lat,  "lat")
    if (allocated (pres)) call get_var (pres, "pres")
    call get_var (wv_basis, "wv_basis")
    call get_var (wv_descr, "wv_descr")
    !---------------
    ! Check metadata
    !---------------
    if (any (lon(2:) <= lon(:dims(1)-1))) then
       call finish ('read_cov_netcdf_common','longitudes not monotonous')
    end if
    if (.not. (all (lat(2:) >= lat(:dims(2)-1)) .or. &
               all (lat(2:) <= lat(:dims(2)-1)) )) then
       call finish ('read_cov_netcdf_common','latitudes not monotonous')
    end if
    if (.not. (all (pres(2:) >= pres(:dims(3)-1)) .or. &
               all (pres(2:) <= pres(:dims(3)-1)) )) then
       call finish ('read_cov_netcdf_common','pressure levels not monotonous')
    end if
    if (any (wv_descr(:) /= wv_name(wv_basis(:)))) then
       print *
       call finish ('read_cov_netcdf_common','wavelet codes/names mismatch')
    end if
    !-----------------
    ! Copy metadata
    !-----------------
    tmpmeta% dims(:)     = dims(:)
    tmpmeta% wv_basis(:) = wv_basis(:)
    if (present (meta)) meta = tmpmeta
#ifdef DEBUG
    !---------
    ! Printout
    !---------
!   print *,'dimensions     :', dims(:)
    print *,'longitudes     :', minval (lon), '...',maxval (lon), &
         trim (units(1))
    print *,'latitudes      :', minval (lat), '...',maxval (lat), &
         trim (units(2))
    print *,'pressure levels:', minval (pres),'...',maxval (pres), &
         trim (units(3))
#endif
  end subroutine read_cov_netcdf_common
  !============================================================================
  function find_cov_netcdf (filename, field) result (r)
    !-------------------------------------------
    ! Check for presence of field in NetCDF file
    !-------------------------------------------
    character(len=*), intent(in) :: filename, field
    logical                      :: r
    !----------------
    ! Local variables
    !----------------
    integer              :: n
    character(len=255)   :: field_list
    character(len=16)    :: fields(MAXMAT)
    !------------------------------------------
    ! Read covariance metadata from NetCDF file
    !------------------------------------------
    call read_cov_netcdf_meta (filename, field_list=field_list)
    call split (fields, field_list, n)
    r = (n > 0) .and. field /= "" .and. any (fields(1:n) == field)
  end function find_cov_netcdf
  !============================================================================
  subroutine test_create_cov_netcdf ()
    !------------------------------------------------------
    ! Simple test writing and reading back a nonsense array
    !------------------------------------------------------
    character(len=*), parameter :: filename = "test-dummy-matrix.nc"
    integer               :: j, k, m, n, ns, nlat
    type(t_cov_meta)      :: meta, meta1
    character(len=255)    :: field_list, description, desc
    character(len=255)    :: profile_list
    character(len=16)     :: fields(MAXMAT), units
    real(wp), allocatable :: profile(:)
    character(len=*), parameter :: profile_name = "testprofile"
    !-------------------------------
    type(t_cov_wavelet),  target  :: W, W2
!   type(t_matrix_block), pointer :: B

    meta% gridinfo = "horizontal-gaussian-latitudes"
    meta% dims      = (/ 16, 8, 1 /)
    meta% wv_basis  = (/ WV_DAUB_8, WV_CDV_8, WV_NONE /)
    meta% gridtype  = WMO6_GAUSSIAN
    meta% representation   = "B"        ! "B" or "L"

    description = "test case with invalid data, don't use"

    call create_cov_netcdf (filename, meta, description)

    ! Set up a small dummy array for testing
    W% meta = meta
    m = product (meta% dims)
    n = m
    ns = 10
    call allocate_block (W% B, repr=CSR, m=m, n=n, ns=ns)
    W% B% ia     = 1
    W% B% ja     = 2
    W% B% packed = 9
    W% valid = .true.
    W% id    = 0
    call append_cov_netcdf (W, filename, field="geof")
    call destruct (W% B)

    W2% meta = meta
    ns = 11
    call allocate_block (W2% B, repr=CSR, m=m, n=n, ns=ns)
    W2% B% ia     = 1
    W2% B% ja     = 3
    W2% B% packed = 8
    W2% valid = .true.
    W2% id    = 0
!   call append_cov_netcdf (W2, filename, field="geof") ! Should fail...
    call append_cov_netcdf (W2, filename, field="rh")
    call destruct (W2% B)

    ! Add a test profile
    nlat = meta% dims(2)
    allocate (profile (nlat))
    profile = (/ (0.1_wp * k, k=1,nlat) /)
    call append_profile_netcdf (profile, filename, profile_name, &
                                "nonsense/units", "Test profile")

!   print *, "Reading back:"
    call read_cov_netcdf_meta (filename, meta1, field_list, description)
!   print *,"description = ", trim (description)
!   print *,"field_list  = ", trim (field_list)
!   print *,"gridinfo              : ", trim (meta1% gridinfo)
!   print *,"gridtype (WMO table 6):" , meta1% gridtype
!   print *,'dimensions            :' , meta1% dims(:)
!   do k = 1, size (meta1% dims)
!      print '(a,i2,a,i12,a,a)',' wavelet basis', k, ":", &
!           meta1% wv_basis(k), " = ", trim (wv_name(meta1% wv_basis(k)))
!   end do
!   print *,"B representation type : ", meta1% representation

    call split (fields, field_list, k)
    do j = 1, k
       call read_cov_netcdf (W,  filename, field=trim (fields(j)))
       if (.not. W% valid) then
          print *, "Error: information is invalid for field ", trim (fields(j))
       end if
!      B     => W% b
!      print *,'Field: ', trim (W% field)
!      print *,'storage_representation: ', trim (srep(B%repr))
!      print *,'m              :' , B%m
!      print *,'n              :' , B%n
!      print *,'ns             :' , B%nonzero
!      print *,'ia             :' ,minval(B%ia),    '...',maxval(B%ia)
!      print *,'ja             :' ,minval(B%ja),    '...',maxval(B%ja)
!      print *,'b              :' ,minval(B%packed),'...',maxval(B%packed)
       call destruct (W% B)
    end do

    call read_profile_list_netcdf (filename, profile_list)
!   print *, "Profile list: ", trim (profile_list)
    if (.not. find_profile_netcdf (filename, profile_name)) then
       call finish ("test_create_cov_netcdf", &
                    "profile not found: " // profile_name)
    end if
    profile = 0
    call read_profile_netcdf (filename, profile_name, profile, units, desc)
!   print *, "Profile:"
!   print *, real (profile(:))
  end subroutine test_create_cov_netcdf
  !============================================================================
  subroutine append_profile_netcdf (var, filename, varname, units, desc, dim)
    real(wp),         intent(in)           :: var(:)    ! (Meridional) Profile
    character(len=*), intent(in)           :: filename  ! NetCDF filename
    character(len=*), intent(in)           :: varname   ! Name of profile
    character(len=*), intent(in)           :: units     ! SI units
    character(len=*), intent(in), optional :: desc      ! Profile description
    integer,          intent(in), optional :: dim       ! Dimension (reserved)
    !----------------
    ! Local variables
    !----------------
    integer            :: n, nlat, dimid_lat, var_id, rc
    character(len=255) :: profile_list
    character(len=16)  :: varnames(MAXPROF)
    character(len=64)  :: description
    !----------------------------------
    ! Preparations & consistency checks
    !----------------------------------
    if (.not. isvalidname (varname)) then
       call finish ('append_profile_netcdf',&
                    'bad variable name: '//trim (varname))
    end if
    if (present (desc)) then
       description = desc
    else
       description = "Profile of " // trim (varname)
    end if
    if (present (dim)) then
       call finish ("append_profile_netcdf", &
                    "optional argument dim not implemented yet!")
    end if
    !----------
    ! Open file
    !----------
    rname = 'append_profile_netcdf: nf90_open'
    call chk (nf90_open (filename, NF90_WRITE, ncid))
    !---------------------------------------
    ! Get old profile list (may be missing!)
    !---------------------------------------
    !call get_attr (profile_list,         '', 'profile_list')
    profile_list = ""
    rc = nf90_get_att (ncid, NF90_GLOBAL, 'profile_list', profile_list)
    if (rc /= NF90_NOERR) then
       write (6,*) "append_profile_netcdf: warning: missing profile_list!"
       profile_list = ""
    end if
    call split (varnames, profile_list, n)
    if (any (varnames(1:n) == varname)) then
       call chk (nf90_abort (ncid))
       call finish ('append_profile_netcdf','cannot modify: '//trim (varname))
    end if
    !----------------------
    ! Consistency of shapes
    !----------------------
    call get_dim (nlat, 'lat')
    if (size (var(:)) /= nlat) then
       call chk (nf90_abort (ncid))
       write (0,*) "nlat, size (var) =", nlat, size (var)
       call finish ('append_profile_netcdf','inconsistent shapes')
    end if
    !----------------------------
    ! Inquire existing dimensions
    !----------------------------
    rname = 'append_profile_netcdf: nf90_inq_dimid'
    call chk (nf90_inq_dimid  (ncid, 'lat', dimid_lat))
    !------------------------
    ! Reenter definition mode
    !------------------------
    rname = 'append_profile_netcdf: nf90_redef'
    call chk (nf90_redef (ncid))
    !---------------------
    ! Define new variables
    !---------------------
    rname = 'append_profile_netcdf: def_var'
    dimid = dimid_lat                   ! var(1:nlat)
    xtype = NF90_FLOAT
    call def_var (trim (varname), var_id, trim (units), trim (description))
    !----------------------
    ! Add/modify attributes
    !----------------------
    rname = 'append_profile_netcdf: nf90_put_att'
    call chk (nf90_put_att (ncid, NF90_GLOBAL, 'profile_list', &
                   adjustl (trim (profile_list)//" "//trim (varname))))
    !------------------------
    ! End of definition block
    !------------------------
    rname = 'append_profile_netcdf: nf90_enddef'
    call chk (nf90_enddef (ncid))
    !--------------------
    ! Store new variables
    !--------------------
    rname = 'append_profile_netcdf: nf90_put_var'
    call chk (nf90_put_var (ncid, var_id ,var(:)  ))
    !-----------
    ! Close file
    !-----------
    call chk (nf90_close (ncid))
  end subroutine append_profile_netcdf
  !----------------------------------------------------------------------------
  subroutine read_profile_list_netcdf (filename, profile_list)
    character(len=*), intent(in)  :: filename
    character(len=*), intent(out) :: profile_list
    !----------------
    ! Local variables
    !----------------
    integer :: i, rc
    !----------
    ! Open file
    !----------
    rname = 'read_profile_list_netcdf: nf90_open'
    call chk (nf90_open (filename, NF90_NOWRITE, ncid))
    !---------------
    ! Get attributes
    !---------------
    !call get_attr (profile_list, '', 'profile_list')
    profile_list = ""
    rc = nf90_get_att (ncid, NF90_GLOBAL, 'profile_list', profile_list)
    if (rc /= NF90_NOERR) then
       write (6,*) "read_profile_list_netcdf: warning: missing profile_list!"
       profile_list = ""
    end if
    !-------------------------------
    ! Remove spurious NUL-characters
    !-------------------------------
    do i=1, len (profile_list)
       if (profile_list(i:i) == NUL) profile_list(i:i) = ' '
!      if (llt (profile_list(i:i), " ")) profile_list(i:i) = " "
    end do
    !-----------
    ! Close file
    !-----------
    rname = 'read_profile_list_netcdf: nf90_close'
    call chk (nf90_close (ncid))
  end subroutine read_profile_list_netcdf
  !----------------------------------------------------------------------------
  subroutine read_profile_netcdf (filename, varname, var, units, desc)
    character(len=*), intent(in)            :: filename ! NetCDF filename
    character(len=*), intent(in)            :: varname  ! Name of profile
    real(wp),         intent(out)           :: var(:)   ! (Meridional) Profile
    character(len=*), intent(out), optional :: units    ! SI units
    character(len=*), intent(out), optional :: desc     ! Description
    !----------------
    ! Local variables
    !----------------
    integer            :: n, nlat, varid, ierr
    character(len=255) :: profile_list
    character(len=16)  :: varnames(MAXPROF)
    !----------
    ! Open file
    !----------
    rname = 'read_profile_netcdf: nf90_open'
    call chk (nf90_open (filename, NF90_NOWRITE, ncid))
    !---------------
    ! Get attributes
    !---------------
    profile_list = ""
    call get_attr (profile_list, '', 'profile_list')
    call split (varnames, profile_list, n)
    if (n == 0) then
       call finish ('read_profile_netcdf','no profiles found in NetCDF file')
    end if
    if (all (varnames(1:n) /= varname)) then
       call finish ('read_profile_netcdf', &
                    'profile not found: '//trim (varname))
    end if
    !------------------------
    ! Get specific dimensions
    !------------------------
    rname = 'read_profile_netcdf: get_dim'
    call get_dim (nlat, 'lat')
    if (size (var(:)) /= nlat) then
       call chk (nf90_abort (ncid))
       write (0,*) trim (varname), " : nlat, size (var) =", nlat, size (var)
       call finish ('read_profile_netcdf','inconsistent shapes')
    end if
    !----------
    ! Read data
    !----------
    var(:) = HUGE (0.0_wp)              ! Initialize to some "invalid" value
    stanc  = 1
    strnc  = 1
    counc  = 0
    rname = 'read_profile_netcdf: get_var'
    call get_var (var, trim (varname))
    !-----------
    ! Attributes
    !-----------
    if (present (units)) then
       units = ""
       vname = varname
       rname = 'read_profile_netcdf: nf90_inq_varid'
       call chk (nf90_inq_varid (ncid, trim (varname), varid))
       ierr = nf90_get_att (ncid, varid, "units", units)
       if (ierr /= 0) units = ""
    end if
    if (present (desc)) then
       desc = ""
       call get_attr (desc,  trim (varname), "long_name")
    end if
    !-----------
    ! Close file
    !-----------
    call chk (nf90_close (ncid))
  end subroutine read_profile_netcdf
  !----------------------------------------------------------------------------
  function find_profile_netcdf (filename, varname) result (r)
    !-------------------------------------------------
    ! Check for presence of a "profile" in NetCDF file
    !-------------------------------------------------
    character(len=*), intent(in)  :: filename
    character(len=*), intent(in)  :: varname
    logical                       :: r
    !----------------
    ! Local variables
    !----------------
    integer            :: n
    character(len=255) :: profile_list
    character(len=16)  :: varnames(MAXPROF)
    !---------------------------------------
    ! Read profile metadata from NetCDF file
    !---------------------------------------
    call read_profile_list_netcdf (filename, profile_list)
    call split (varnames, profile_list, n)
    r = (n > 0) .and. varname /= "" .and. any (varnames(1:n) == varname)
  end function find_profile_netcdf
  !============================================================================
  subroutine merge_netcdf (outfile, filelist, selection, title)
    !--------------------------------------------------------------
    ! Copy or merge NetCDF files containing covariances or profiles
    ! based on array "selection".  If selection(:) is empty or "*",
    ! all covariances and profiles are copied, excluding duplicates.
    ! A name in selection(:) prepended by a "-" will be removed.
    !--------------------------------------------------------------
    character(len=*), intent(in) :: outfile         ! Output NetCDF file
    character(len=*), intent(in) :: filelist(:)     ! Input  NetCDF files
    character(len=*), intent(in) :: selection(:)    ! Covariance/profile name
    character(len=*), intent(in), optional :: title ! Replacement title
    !----------------
    ! Local variables
    !----------------
    integer               :: nfiles, nselect, i, j, k, ny
    integer,  allocatable :: n_fields(:), n_profiles(:)
    real(wp), allocatable :: profile(:)
    logical               :: keep
    character(len=16)     :: current_field = ""
    character(len=16)     :: current_profile = ""
    character(len=16)     :: current_name = ""
    character(len=16)     :: units = ""
    character(len=128)    :: desc  = ""
    character(len=16),  allocatable :: fields(:,:)
    character(len=16),  allocatable :: profiles(:,:)
    character(len=255), allocatable :: field_list(:)
    character(len=255), allocatable :: description(:)
    character(len=255), allocatable :: profile_list(:)
    !--------------------------------------------------------------------------
    type(t_cov_meta),   allocatable :: meta(:)      ! Covariance metadata
    type(t_cov_wavelet)             :: W            ! Wavelet covariances
    !--------------------------------------------------------------------------
    ! Check presence of input files first
    nfiles  = size (filelist)
    nselect = size (selection)
    if (nfiles == 0) then
       call finish ("merge_netcdf","no input files")
    end if

    allocate (meta(nfiles))
    allocate (field_list(nfiles))
    allocate (profile_list(nfiles))
    allocate (description(nfiles))
    allocate (fields(MAXMAT,nfiles))
    allocate (n_fields(nfiles))
    allocate (n_profiles(nfiles))
    allocate (profiles(MAXPROF,nfiles))
    field_list   = ""
    profile_list = ""
    fields       = ""
    profiles     = ""
    n_fields     = 0
    n_profiles   = 0
    do i = 1, nfiles
       if (filelist(i) == "") cycle
       print *
       print *, "Reading metadata from file: ", trim (filelist(i))
       call read_cov_netcdf_meta     (filelist(i), meta(i), field_list(i), &
                                      description(i))
       call read_profile_list_netcdf (filelist(i), profile_list(i))

       call split (fields(:,i),   field_list(i),   n_fields(i))
       call split (profiles(:,i), profile_list(i), n_profiles(i))

       if (n_fields(i) == 0 .and. n_profiles(i) == 0) then
          print *, "*** File ", trim (filelist(i)), " is empty!?"
          cycle
       end if
       if (i <= nselect) then
          current_name = selection(i)
          keep = .true.
          if (current_name(1:1) == "-") then
             keep = .false.
             current_name = adjustl (current_name(2:))
          end if
          if (current_name /= "" .and. current_name /= "*") then
             if (any (fields(:,i) == current_name)) then
                if (keep) then
                   ! Retain only the selected field
                   fields(1,i)   = current_name
                   fields(2:,i)  = ""
                   n_fields(i)   = 1               ! Single field
                   n_profiles(i) = 0
!                  print *, "Found field ", trim (fields(1,i))
                else                               ! Don't keep
                   do j = 1, n_fields(i)
                      if (fields(j,i) == current_name) then
                         ! Found it, now remove it
                         fields(j:n_fields(i)-1,i) = fields(j+1:n_fields(i),i)
                         fields(  n_fields(i),  i) = ""
                         n_fields(i)               = n_fields(i) - 1
                         exit
                      end if
                   end do
                end if
             else if (any (profiles(:,i) == current_name)) then
                if (keep) then
                   ! Retain only the selected profile
                   profiles(1,i)  = current_name
                   profiles(2:,i) = ""
                   n_profiles(i)  = 1              ! Single profile
                   n_fields(i)    = 0
!                  print *, "Found profile ", trim (profiles(1,i))
                else                               ! Don't keep
                   do j = 1, n_profiles(i)
                      if (profiles(j,i) == current_name) then
                         ! Found it, now remove it
                         profiles(j:n_profiles(i)-1,i) = &
                              profiles(j+1:n_profiles(i),i)
                         profiles(  n_profiles(i),  i) = ""
                         n_profiles(i)                 = n_profiles(i) - 1
                         exit
                      end if
                   end do
                end if
             else
                write (0,*) "field_list:   '", trim (field_list(i)), "'"
                write (0,*) "profile_list: '", trim (profile_list(i)), "'"
                write (0,*) "selection:    '", trim (selection(i)), "'"
                call finish ("merge_netcdf","selection not in any list")
             end if
          else
             if (.not. keep) then
                n_fields(i)   = 0       ! Skip all fields and profiles
                n_profiles(i) = 0       ! contained in this file
                fields(:,i)   = ""
                profiles(:,i) = ""
             end if
          end if
       end if
    end do
    print *
    print *, "Checking compatibility of metadata"
    print *
    if (any (meta(2:)% gridtype /= meta(1)% gridtype)) then
       write (0,*) meta(:)% gridtype
       call finish ("merge_netcdf","incompatible gridtypes")
    end if
    if (any (meta(2:)% representation /= meta(1)% representation)) then
       write (0,*) meta(:)% representation
       call finish ("merge_netcdf","incompatible B representations")
    end if
    do k = 1, 3
       if (any (meta(2:)% dims(k) /= meta(1)% dims(k))) then
          write (0,*) "Dimension", k,":", meta(:)% dims(k)
          call finish ("merge_netcdf","incompatible grid size")
       end if
    end do
    do k = 1, 3
       if (any (meta(2:)% wv_basis(k) /= meta(1)% wv_basis(k))) then
          write (0,*) "Dimension", k,":", meta(:)% wv_basis(k)
          call finish ("merge_netcdf","incompatible wavelet bases")
       end if
    end do

    ny = meta(1)% dims(2)
    allocate (profile(ny))

    ! Create empty NetCDF file with metadata from first file in the input list
    print *, "Creating NetCDF file from template"
    print *
    if (present (title)) description(1) = title
    call create_cov_netcdf (outfile, meta(1), trim (description(1)))

    do i = 1, nfiles
       if (filelist(i) == "") cycle
       do j = 1, n_fields(i)
          current_field = fields(j,i)
          if (current_field == "") then
             call finish ("merge_netcdf", "internal error: current_field=''")
          end if
          if (any (fields(:,1:i-1) == current_field)) then
             print *, "*** Warning: ", trim (current_field), &
                  " already copied.  Skipping for file: ", trim (filelist(i))
             cycle
          end if
          print *, "Copying ", trim (current_field), " from ", &
               trim (filelist(i))
          call read_cov_netcdf   (W, filename=trim (filelist(i)), &
                                  field=trim (current_field))
          call append_cov_netcdf (W, filename=trim (outfile), &
                                  field=trim (current_field))
          call destruct (W)
       end do
       do j = 1, n_profiles(i)
          current_profile = profiles(j,i)
          if (current_profile == "") then
             call finish ("merge_netcdf", "internal error: current_profile=''")
          end if
          if (any (profiles(:,1:i-1) == current_profile)) then
             print *, "*** Warning: ", trim (current_profile), &
                  " already copied.  Skipping for file: ", trim (filelist(i))
             cycle
          end if
          print *, "Copying ", trim (current_profile), " from ", &
               trim (filelist(i))
          units = ""
          desc  = ""
          profile = 0
          call read_profile_netcdf (filename = trim (filelist(i)),     &
                                    varname  = trim (current_profile), &
                                    var      = profile(:),             &
                                    units    = units,                  &
                                    desc     = desc)
!         print *, "Units      : ", trim (units)
!         print *, "Description: ", trim (desc)
          call append_profile_netcdf (var      = profile,                &
                                      filename = trim (outfile),         &
                                      varname  = trim (current_profile), &
                                      units    = trim (units),           &
                                      desc     = trim (desc))
       end do
       print *
    end do
    deallocate (profile)
    deallocate (meta)
  end subroutine merge_netcdf
  !============================================================================
  function isvalidname (name) result (r)
    !---------------------------------
    ! Check variable name for validity
    !---------------------------------
    character(len=*), intent(in) :: name
    logical                      :: r
    !----------------
    ! Local variables
    !----------------
    integer               :: pos
    character(len (name)) :: string
    character(len=*), parameter :: LOWER   = "abcdefghijklmnopqrstuvwxyz"
    character(len=*), parameter :: DIGITS  = "0123456789"
    character(len=*), parameter :: CHARSET = LOWER // DIGITS // "_"

    string = tolower (name)
    !---------------------------------
    ! First character must be a letter
    !---------------------------------
    r = (len_trim (string) > 0)
    if (.not. r) return
    pos = verify (string(1:1),   set=LOWER)
    r = (pos == 0)
    if (.not. r) return
    !-----------------------------------------
    ! Other characters must be from chosen set
    !-----------------------------------------
    pos = verify (trim (string), set=CHARSET)
    r = (pos == 0)
  end function isvalidname
  !============================================================================
  ! MPI communication routines for wavelet covariance matrices and metadata
  ! (send receive bcast)
  !----------------------------------------------------------------------------
  subroutine send_cov_meta (meta, dest, tag, comm)
    type(t_cov_meta),  intent(in) :: meta
    integer,           intent(in) :: dest
    integer, optional, intent(in) :: tag
    integer, optional, intent(in) :: comm

    integer :: ltag, lcom
    integer :: count = 0

    if (count == 0) count = size (transfer (meta, (/' '/)))

    ltag = TAG_COV_META ;if (present (tag))  ltag = tag
    lcom = dace% comm   ;if (present (comm)) lcom = comm

    call p_send_derivedtype (meta, count, dest, ltag, lcom)

  end subroutine send_cov_meta
  !----------------------------------------------------------------------------
  subroutine send_cov_wavelet (cov, dest, tag, comm)
    type(t_cov_wavelet), intent(in) :: cov
    integer,             intent(in) :: dest
    integer, optional,   intent(in) :: tag
    integer, optional,   intent(in) :: comm

    integer :: ltag, lcom

    ltag = TAG_COV_WAVELET ;if (present (tag))  ltag = tag
    lcom = dace% comm      ;if (present (comm)) lcom = comm

    call p_send (cov% b,     dest, ltag, comm=lcom)
    call p_send (cov% meta,  dest, ltag, comm=lcom)
    call p_send (cov% field, dest, ltag, comm=lcom)
    call p_send (cov% id,    dest, ltag, comm=lcom)
    call p_send (cov% valid, dest, ltag, comm=lcom)

  end subroutine send_cov_wavelet
  !----------------------------------------------------------------------------
  subroutine recv_cov_meta (meta, source, tag, comm)
    type(t_cov_meta),  intent(inout) :: meta
    integer,           intent(in)    :: source
    integer, optional, intent(in)    :: tag
    integer, optional, intent(in)    :: comm

    integer :: ltag, lcom
    integer :: count = 0

    if (count == 0) count = size (transfer (meta, (/' '/)))

    ltag = TAG_COV_META ;if (present (tag))  ltag = tag
    lcom = dace% comm   ;if (present (comm)) lcom = comm

    call p_recv_derivedtype (meta, count, source, ltag, lcom)

  end subroutine recv_cov_meta
  !----------------------------------------------------------------------------
  subroutine recv_cov_wavelet (cov, source, tag, comm)
    type(t_cov_wavelet), intent(inout) :: cov
    integer,             intent(in)    :: source
    integer, optional,   intent(in)    :: tag
    integer, optional,   intent(in)    :: comm

    integer :: ltag, lcom

    ltag = TAG_COV_WAVELET ;if (present (tag))  ltag = tag
    lcom = dace% comm      ;if (present (comm)) lcom = comm

    call p_recv (cov% b,     source, ltag, comm=lcom)
    call p_recv (cov% meta,  source, ltag, comm=lcom)
    call p_recv (cov% field, source, ltag, comm=lcom)
    call p_recv (cov% id,    source, ltag, comm=lcom)
    call p_recv (cov% valid, source, ltag, comm=lcom)

  end subroutine recv_cov_wavelet
  !----------------------------------------------------------------------------
  subroutine bcast_cov_meta (meta, source, comm)
    type(t_cov_meta),  intent(inout) :: meta
    integer,           intent(in)    :: source
    integer, optional, intent(in)    :: comm

    integer :: lcom
    integer :: count = 0

    if (count == 0) count = size (transfer (meta, (/' '/)))

    lcom = dace% comm ;if (present (comm)) lcom = comm

    call p_bcast_derivedtype (meta, count, source, lcom)

  end subroutine bcast_cov_meta
  !----------------------------------------------------------------------------
  subroutine bcast_cov_wavelet (cov, source, comm)
    type(t_cov_wavelet), intent(inout) :: cov
    integer,             intent(in)    :: source
    integer, optional,   intent(in)    :: comm

    integer :: lcom

    lcom = dace% comm ;if (present (comm)) lcom = comm

    call p_bcast (cov% b,     source, comm=lcom)
    call p_bcast (cov% meta,  source, comm=lcom)
    call p_bcast (cov% field, source, comm=lcom)
    call p_bcast (cov% id,    source, comm=lcom)
    call p_bcast (cov% valid, source, comm=lcom)

  end subroutine bcast_cov_wavelet
  !============================================================================
end module mo_cov_wavelet
