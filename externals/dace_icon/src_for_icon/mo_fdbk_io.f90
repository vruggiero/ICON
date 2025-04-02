!
!+ routines and derived type for I/O of 3DVAR/COSMO feedback files
!
MODULE mo_fdbk_io
!
! Description:
!   Routines and derived type for I/O of 3DVAR/COSMO feedback files.
!   This module is jointly used by the COSMO model and 3DVAR program packages
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
! V1_13        2011/11/01 Andreas Rhodin
!  routines and derived type for feedback file I/O
! V1_15        2011/12/06 Andreas Rhodin
!  make error processing COSMO conform (use model_abort from module environment)
!  process (optional) body variable 'plevel'.
!  read_fdbk_head, read_fdbk_body : read complete set of variables.
!  write_fdbk_file                : new routine.
! V1_17        2011/12/21 Andreas Rhodin
!  default initialization of header and body tables with fillvalues
!  new optional argument 'close' to subroutine 'write_fdbk_file'
! V1_18        2012/01/06 Andreas Rhodin
!  Extended comments and documentation
! V1_19        2012-04-16 Andreas Rhodin
!  changes for Radiance fof files (LMU)
! V1_20        2012-06-18 Andreas Rhodin
!  account for possible fillvalues when writing 'e_o' to feedback file
! V1_22        2013-02-13 Andreas Rhodin
!  new routine 'scatter_fdbk_data' (preliminary, currently not portable).
!  read_fdbk_body: allow changes to read body after thinning of header.
! V1_26        2013/06/27 Andreas Rhodin
!  implement optional header entries: tracking, meas_type, rad_corr
!  ensure portability with COSMO
! V1_27        2013-11-08 Harald Anlauf
!  make compilable with MPI checksumming enabled (-DMPI_CHECKSUM, 3dvar only)
! V1_42        2015-06-08 Harald Anlauf
!  Read drifting lat/lon from feedback files
! V1_44        2015-09-30 Christoph Schraff
!  saveguard for no veri data present
! V1_45        2015-12-15 Andreas Rhodin
!  changed interface of generic subroutine get_var
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
! Authors:
! Andreas Rhodin  DWD  2010  original code
!------------------------------------------------------------------------------

!-------------
! Modules used
!-------------
  use kind_parameters, only: wp                 ! working precision kind parameter
  !----------------
  ! MPI environment
  !----------------
  use environment, only:               &!
                     model_abort,      &! abort in case of error
                     MPI_BYTE,         &! MPI specific parameters
                     MPI_INTEGER,      &!
                     MPI_SUCCESS        !

#ifdef __COSMO__
  use environment,   only:             &!
           d_comm => MPI_COMM_WORLD     ! default communicator
  use data_parallel, only:             &!
            d_npe => nproc              ! number of processor elements
#else
  use mo_mpi_dace, only: d_npe,        &! number of processor elements
                         d_comm         ! default communicator
#endif
  !-------------------
  ! NetCDF environment
  !-------------------
  use netcdf,  only:                   &! NetCDF f90 interface
                     nf90_get_var,     &! read NetCDF variable
                     nf90_get_att,     &! read NetCDF attributes
                     nf90_enddef,      &! close definition mode
                     nf90_strerror,    &! explanation of NetCDF error code
                     NF90_NOERR         ! NetCDF return code for no error
  !-------------------------------
  ! feedback file specific modules
  !-------------------------------
  use mo_fdbk, only: t_fdbk,           &! feedback file data type
                     t_fdbk_meta,      &! feedback file meta data
                     setup_fdbk,       &! set up feedback-file data structure
                     open_fdbk_read,   &! open feedback file for read access
                     read_meta,        &! read meta data
                     get_varid,        &! get NetCDF varid
                     get_veri_index,   &! return indices of verification runs
                     get_veri,         &! read verification entry
                     get_fillvalue,    &! get fillvalue
                     close_fdbk,       &! close feedback file
                     cleanup_fdbk,     &! deallocate components
!                    create_fdbk,      &! create new feedback file
                     add_history,      &! add history to feedback file
                     write_fdbk_var,   &! write variable to feedback file
                     write_global_attributes !
  use mo_t_netcdf_file,                &!
               only: create_netcdf_file !

implicit none

!----------------
! Public entities
!----------------
private
public :: t_fdbk_body        ! component of read_fdbk_data (body   table entry)
public :: t_fdbk_head        ! component of read_fdbk_data (header table entry)
public :: t_fdbk_data        ! derived type to hold feedback file content
public :: ifill, rfill       ! fillvalues used generally in this module
public :: read_fdbk_data     ! read feedback file, store in t_fdbk_data
public :: read_fdbk_head     ! read feedback file header
public :: read_fdbk_body     ! read feedback file body
public :: read_fdbk_veri     ! read verification data
public :: read_fdbk_var_head ! read variable from header
public :: read_fdbk_var_body ! read variable from body
public :: pack_fdbk          ! pack header and body according to masks
public :: destruct           ! deallocate components of t_fdbk_data
public :: associate_hb       ! associate related header and body entries
public :: get_var            ! directly read NetCDF variable
public :: write_fdbk_file    ! write complete feedback file
public :: scatter_fdbk_data  ! scatter derived type t_fdbk_data

!------------------
! Public interfaces
!------------------

interface destruct
  module procedure destruct_fdbk_data  ! deallocate components of t_fdbk_data
end interface destruct

interface get_var
  module procedure get_var_real
  module procedure get_var_real_2
  module procedure get_var_int
  module procedure get_var_char
end interface get_var

interface read_fdbk_var_head
  module procedure read_fdbk_ivar_head ! read integer variable from header
  module procedure read_fdbk_rvar_head ! read real    variable from header
end interface read_fdbk_var_head

interface read_fdbk_var_body
  module procedure read_fdbk_ivar_body ! read integer variable from body
  module procedure read_fdbk_rvar_body ! read real    variable from body
end interface read_fdbk_var_body

!--------------------------------
! Public derived type definitions
!--------------------------------

  integer  ,parameter :: ifill = -huge(0)
  real(wp) ,parameter :: rfill = -huge(0._wp)

  !-----------------------------
  ! feedback file header entries
  !-----------------------------
  type t_fdbk_head
    !------------------
    ! mandatory entries
    !------------------
    integer           :: i_body        = ifill
    integer           :: l_body        = ifill
    integer           :: n_level       = ifill ! number of levels in body
    integer           :: obstype       = ifill
    integer           :: codetype      = ifill
    integer           :: center        = ifill ! station processing center
    integer           :: sub_center    = ifill ! processing subcenter
    integer           :: data_category = ifill
    integer           :: sub_category  = ifill
    integer           :: ident         = ifill
    character(len=10) :: statid        = ''
    character(len=32) :: wsi           = ''
    real(wp)          :: lat           = rfill
    real(wp)          :: lon           = rfill
    integer           :: time          = ifill
    integer           :: time_dbase    = ifill ! data base time
    integer           :: time_nomi     = ifill ! nominal (synoptic) time
    real(wp)          :: sun_zenit     = rfill ! solar zenith angle
    integer           :: mdlsfc        = ifill
    integer           :: z_station     = ifill ! station height
    integer           :: z_modsurf     = ifill ! model surface height
    integer           :: sta_corr      = ifill ! correction indicator
    integer           :: instype       = ifill ! instrument type
    integer           :: r_state       = ifill
    integer           :: r_flags       = ifill
    integer           :: r_check       = ifill
    integer           :: index_x       = ifill
    integer           :: index_y       = ifill
    !----------------------
    ! some optional entries
    !----------------------
    integer           :: dbkz          = -1
    integer           :: record        = -1
    integer           :: subset        = -1
    integer           :: source        = -1
    integer           :: phase         = -1
    integer           :: obs_id        = -1
!   integer           :: flg_1dvar     = -1
    integer           :: flg_cld       = -1
    integer           :: index_d       = -1
    integer           :: surftype      = -1
    integer           :: soiltype      = -1
    integer           :: tracking      = -1
    integer           :: meas_type     = -1
    integer           :: rad_corr      = -1
    integer           :: scanline      = -1
    integer           :: fov           = -1
    real(wp)          :: sat_zenit     = -999._wp
    real(wp)          :: sat_azimuth   = -999._wp
    real(wp)          :: sun_azimuth   = -999._wp
    !-------------------------------------
    ! pointer to original position in file
    !-------------------------------------
    integer           :: pos           = 0        ! position in file
    !------------------------
    ! pointer to body entries
    !------------------------
    integer           :: ib            = 0        ! start index in body
    integer           :: nb            = 0        ! number of body entries
    integer           :: pe            = -1       ! processor element to use
  end type t_fdbk_head

  !---------------------------
  ! feedback file body entries
  !---------------------------
  type t_fdbk_body
    !-------------------------
    ! content of feedback file
    !-------------------------
    real(wp)          :: obs       = rfill
    real(wp)          :: bcor      = rfill
    real(wp)          :: e_o       = rfill
    real(wp)          :: level     = -99999._wp   ! fix missing value in fofs
    real(wp)          :: plevel    = -1._wp       ! optional
    real(wp)          :: lat       = rfill        ! optional (drifting) latitude
    real(wp)          :: lon       = rfill        ! optional (drifting) longitude
    real(wp)          :: plev_width= rfill        ! optional width of sensitivity fctn
    real(wp)          :: v_loc     = rfill        ! optional horizontal localisation
    real(wp)          :: h_loc     = rfill        ! optional vertical   localisation
    real(wp)          :: emiss     = rfill        ! optional emissivity
    integer           :: level_typ = ifill
    integer           :: level_sig = -99
    integer           :: varno     = ifill
    integer           :: state     = ifill
    integer           :: flags     = ifill
    integer           :: check     = ifill
    integer           :: qual      = -1
    !----------------
    ! pointer to file
    !----------------
    integer           :: pos       = 0 ! position in file
    !------------------------
    ! pointer to header entry
    !------------------------
    integer           :: ih        = 0 ! index of header entry
  end type t_fdbk_body

  !------------------------------------------
  ! container for feedback file data:
  !   mandatory header and body entries,
  !   verification meta data
  !------------------------------------------
  type t_fdbk_data
    type(t_fdbk)               :: f               ! feedback file meta data
    type(t_fdbk_head) ,pointer :: h (:) => NULL() ! feedback file header
    type(t_fdbk_body) ,pointer :: b (:) => NULL() ! feedback file body data
    real(wp)          ,pointer :: veri_data(:,:) => NULL() ! verification data
    type(t_fdbk_meta) ,pointer :: veri_meta  (:) => NULL() ! verific.meta data
  end type t_fdbk_data

!------------------------------------------------------------------------------
contains
!------------------------------------------------------------------------------
  subroutine read_fdbk_meta (data, file)
  type (t_fdbk_data) ,intent(inout) :: data
  character(len=*)   ,intent(in)    :: file
  !-----------------------------
  ! read feedback file mete data
  !-----------------------------

    !------------------------------
    ! clean up old data, close file
    !------------------------------
    call destruct (data)
    !------------------------------------
    ! set up feedback file data structure
    !------------------------------------
    call setup_fdbk (data% f% nc)
    !-------------------
    ! open feedback file
    !-------------------
    print *,'opening feedback file:',trim(file)
    call open_fdbk_read (data% f, file)
    !------------------------------------------------
    ! read global attributes + verification meta data
    !------------------------------------------------
    print *,'reading meta data'
    call read_meta (data% f, file)
    !-----------
    ! close file
    !-----------
    call close_fdbk (data% f)

  end subroutine read_fdbk_meta
!------------------------------------------------------------------------------
  subroutine read_fdbk_head (data, file)
  type (t_fdbk_data) ,intent(inout)        :: data
  character(len=*)   ,intent(in) ,optional :: file
  !----------------------------------
  ! read feedback file header entries
  !----------------------------------

    integer                     :: n_hdr
    integer                     :: i
    type (t_fdbk_head) ,pointer :: h (:)
    integer                     :: ierr

    !--------------------------
    ! read meta-data, open file
    !--------------------------
    if (present (file)) then
      call read_fdbk_meta (data, file)
    endif
    print *,'re-opening feedback file:',trim(data% f% nc% path)
    call open_fdbk_read (data% f, data% f% nc% path)

    !----------------
    ! allocate arrays
    !----------------
    if (associated (data% h)) deallocate (data% h)
    n_hdr  = data% f% n_hdr
    allocate (data% h (n_hdr))
    h => data% h

    !------------------------------
    ! read mandatory header entries
    !------------------------------
    print *,'reading header'
    call get_var_int  (data%f, 'i_body'       ,h% i_body        ,fill=ifill)
    call get_var_int  (data%f, 'l_body'       ,h% l_body        ,fill=ifill)
    call get_var_int  (data%f, 'obstype'      ,h% obstype                  )
    call get_var_int  (data%f, 'codetype'     ,h% codetype                 )
    call get_var_int  (data%f, 'data_category',h% data_category ,fill=ifill)
    call get_var_int  (data%f, 'sub_category' ,h% sub_category  ,fill=ifill)
    call get_var_int  (data%f, 'ident'        ,h% ident         ,fill=ifill)
    call get_var_int  (data%f, 'time'         ,h% time                     )
    call get_var_int  (data%f, 'time_nomi'    ,h% time_nomi     ,fill=ifill)
    call get_var_int  (data%f, 'time_dbase'   ,h% time_dbase    ,fill=ifill)
    call get_var_char (data%f, 'statid'       ,h% statid                   )
    call get_var_real (data%f, 'lat'          ,h% lat                      )
    call get_var_real (data%f, 'lon'          ,h% lon                      )
    call get_var_int  (data%f, 'mdlsfc'       ,h% mdlsfc        ,fill=ifill)
    call get_var_int  (data%f, 'r_state'      ,h% r_state                  )
    call get_var_int  (data%f, 'r_flags'      ,h% r_flags                  )
    call get_var_int  (data%f, 'r_check'      ,h% r_check                  )
    call get_var_int  (data%f, 'index_x'      ,h% index_x                  )
    call get_var_int  (data%f, 'index_y'      ,h% index_y                  )
    call get_var_int  (data%f, 'sub_center'   ,h% sub_center    ,fill=ifill)
    call get_var_int  (data%f, 'n_level'      ,h% n_level       ,fill=ifill)
    call get_var_int  (data%f, 'center'       ,h% center        ,fill=ifill)
    call get_var_int  (data%f, 'z_station'    ,h% z_station     ,fill=ifill)
    call get_var_int  (data%f, 'z_modsurf'    ,h% z_modsurf     ,fill=ifill)
    call get_var_int  (data%f, 'sta_corr'     ,h% sta_corr      ,fill=ifill)
    call get_var_int  (data%f, 'instype'      ,h% instype       ,fill=ifill)
    call get_var_real (data%f, 'sun_zenit'    ,h% sun_zenit     ,fill=rfill &
                                                              ,ierr=ierr )
    !---------------------------
    ! read some optional entries
    !---------------------------
    call get_var_int  (data%f, 'dbkz'         ,h% dbkz           ,fill=-1)
    call get_var_int  (data%f, 'subset'       ,h% subset         ,fill=-1)
    call get_var_int  (data%f, 'record'       ,h% record         ,fill=-1)
    call get_var_int  (data%f, 'source'       ,h% source         ,fill=-1)
    call get_var_int  (data%f, 'phase'        ,h% phase          ,fill=-1)
!   call get_var_int  (data%f, 'flg_1dvar'    ,h% flg_1dvar      ,fill=-1)
    call get_var_int  (data%f, 'flg_cld'      ,h% flg_cld        ,fill=-1)
    call get_var_int  (data%f, 'surftype'     ,h% surftype       ,fill=-1)
    call get_var_int  (data%f, 'soiltype'     ,h% soiltype       ,fill=-1)
    call get_var_int  (data%f, 'tracking'     ,h% tracking       ,fill=-1)
    call get_var_int  (data%f, 'meas_type'    ,h% meas_type      ,fill=-1)
    call get_var_int  (data%f, 'rad_corr'     ,h% rad_corr       ,fill=-1)
    call get_var_int  (data%f, 'scanline'     ,h% scanline       ,fill=-1)
    call get_var_int  (data%f, 'fov'          ,h% fov            ,fill=-1)
    call get_var_real (data%f, 'sat_zenit'    ,h% sat_zenit      ,fill=-999._wp)
    call get_var_real (data%f, 'sat_azimuth'  ,h% sat_azimuth    ,fill=-999._wp)
    call get_var_real (data%f, 'sun_azimuth'  ,h% sun_azimuth    ,fill=-999._wp)
    call get_var_int  (data%f, 'obs_id'       ,h% obs_id         ,fill=-1)
    call get_var_int  (data%f, 'index_d'      ,h% index_d        ,fill=-1)

    if (n_hdr > 0) then
      !-------------------------------------------------
      ! fix for old convention (i_body starting at zero)
      !-------------------------------------------------
      if (h(1)% i_body == 0) h% i_body = h% i_body + 1
      !--------------------------------------------
      ! consistency check on i_body, l_body, n_body
      !--------------------------------------------
      do i=1, n_hdr - 1
        if (h(i+1)% i_body - h(i)% i_body /= h(i)% l_body) then
          write (6,*)                                               &
          'inconsistent i, n_hdr, i_body, i_body+, diff, l_body :', &
           i, n_hdr, h(i)% i_body, h(i+1)% i_body,                  &
           h(i+1)% i_body - h(i)% i_body, h(i)% l_body
          call model_abort(-1,-1,'inconsistent i_body/l_body','read_fdbk_head')
        endif
      end do
      if (h(n_hdr)% i_body + h(n_hdr)% l_body - 1 /= data% f% n_body) then
        write (6,*)                                                  &
          'inconsistent n_hdr i_body l_body i+l-1 n_body:',          &
           n_hdr, h(n_hdr)% i_body, h(n_hdr)% l_body,                &
           h(n_hdr)% i_body + h(n_hdr)% l_body - 1, data% f% n_body
        call model_abort (-1,-1,'inconsistent i_body/l_body/n_body', &
                                'read_fdbk_head'                     )
      endif
    endif
    !------------
    ! set indices
    !------------
    h% pos = (/(i,i=1,n_hdr)/)
    h% nb  = 0
    h% ib  = 0
    !-----------
    ! close file
    !-----------
    call close_fdbk (data% f)

  end subroutine read_fdbk_head
!------------------------------------------------------------------------------
  subroutine read_fdbk_body (data, file)
  type (t_fdbk_data) ,intent(inout)        :: data
  character(len=*)   ,intent(in) ,optional :: file
  !--------------------------------
  ! read feedback file body entries
  !--------------------------------

    integer              :: n_body
    integer              :: i, ib
    logical ,allocatable :: mask (:)

    !--------------------------
    ! read meta-data, open file
    !--------------------------
    if (present (file)) then
      call read_fdbk_meta (data, file)
    endif
    print *,'re-opening feedback file:',trim(data% f% nc% path)
    call open_fdbk_read (data% f, data% f% nc% path)

    !----------------
    ! allocate arrays
    !----------------
    if (associated (data% b)) deallocate (data% b)
    n_body = sum   (data% h% l_body)
    allocate (mask (data% f% n_body))
    allocate (data% b (      n_body))

    !--------------------------
    ! set mask for partial read
    !--------------------------
    mask = .false.
    ib   = 0
    do i = 1, size(data% h)
      if (data%h(i)% l_body == 0) cycle
      mask (data%h(i)% i_body                        : &
            data%h(i)% i_body + data%h(i)% l_body - 1) = .true.
      data% h(i)%                         ib = ib+1
      data% h(i)%                         nb = data%h(i)% l_body
      data% b(ib+1:ib+data%h(i)% l_body)% ih = i
      ib = ib + data%h(i)% l_body
    end do

    !----------
    ! read body
    !----------
    print *,'reading body'
    call get_var_real (data%f, 'obs'       ,data% b% obs      ,mask=mask,fill=rfill     )
    call get_var_real (data%f, 'bcor'      ,data% b% bcor     ,mask=mask,fill=rfill     )
    call get_var_real (data%f, 'e_o'       ,data% b% e_o      ,mask=mask,fill=rfill     )
    call get_var_real (data%f, 'level'     ,data% b% level    ,mask=mask,fill=-99999._wp)
    call get_var_real (data%f, 'plevel'    ,data% b% plevel   ,mask=mask,fill=-1._wp    )
    call get_var_real (data%f, 'plev_width',data% b% plev_width,mask=mask,fill=-1._wp   )
    call get_var_real (data%f, 'dlat'      ,data% b% lat      ,mask=mask,fill=rfill     )
    call get_var_real (data%f, 'dlon'      ,data% b% lon      ,mask=mask,fill=rfill     )
    call get_var_int  (data%f, 'level_typ' ,data% b% level_typ,mask=mask)
    call get_var_int  (data%f, 'level_sig' ,data% b% level_sig,mask=mask,fill=-99       )
    call get_var_int  (data%f, 'varno'     ,data% b% varno    ,mask=mask)
    call get_var_int  (data%f, 'state'     ,data% b% state    ,mask=mask)
    call get_var_int  (data%f, 'flags'     ,data% b% flags    ,mask=mask)
    call get_var_int  (data%f, 'check'     ,data% b% check    ,mask=mask)
    call get_var_int  (data%f, 'qual'      ,data% b% qual     ,mask=mask,fill=-1        )
    data% b% pos = pack((/(i,i=1,data% f% n_body)/),mask)

    !-----------
    ! close file
    !-----------
    call close_fdbk (data% f)

  end subroutine read_fdbk_body
!------------------------------------------------------------------------------
  subroutine read_fdbk_veri (data, fill, ix, model, run_type, run_class,    &
                             initial_date, forecast_time, ens_member, exp_id)
  type (t_fdbk_data) ,intent(inout)        :: data
  real(wp)           ,intent(in) ,optional :: fill          ! new fillvalue
  integer            ,intent(in) ,optional :: ix (:)
  character(len=*)   ,intent(in) ,optional :: model
  integer            ,intent(in) ,optional :: run_type
  integer            ,intent(in) ,optional :: run_class
  character(len=*)   ,intent(in) ,optional :: initial_date
  integer            ,intent(in) ,optional :: forecast_time
  integer            ,intent(in) ,optional :: ens_member
  integer            ,intent(in) ,optional :: exp_id
  !------------------------------------------------------------
  ! read verification data.
  ! store in data% veri_data.
  ! select as specified by optional parameters
  ! or by explicit provision of verification indices 'ix'.
  ! filter as implied by filtered component index data% b% pos.
  !------------------------------------------------------------

    integer ,parameter    :: mv = 1200 ! max number of verification entries
    integer               :: nb        ! number of body entries used
    integer               :: nv        ! number of verification entries used
    integer               :: iv (mv)   ! indices of verification entries used
    integer               :: i         ! loop index
    real(wp)              :: fillvalue ! fillvalue in file
    real(wp) ,allocatable :: tmp (:)   ! temporary

    !--------------------------
    ! read meta-data, open file
    !--------------------------
    print *,'re-opening feedback file:',trim(data% f% nc% path)
    call open_fdbk_read (data% f, data% f% nc% path)

    !-------------------------------------------------------
    ! determine number of body and verification entries used
    !-------------------------------------------------------
    nb = data% f% n_body
    if (associated (data% b)) nb = size(data% b)
    if (present (ix)) then
      nv = size(ix)
      if (nv > mv)                                                            &
        call model_abort (-1,-1,"size of 'iv' is too small !","read_fdbk_veri")
      iv (1:nv) = ix
    else
      call get_veri_index (iv, nv, data% f, model, run_type, run_class,   &
                           initial_date, forecast_time, ens_member, exp_id)
      if (nv < 0) &
        call model_abort (-1,-1,"size of 'iv' is too small !","read_fdbk_veri")
    endif
    if (associated (data% veri_data)) then
      deallocate (data% veri_data)
      deallocate (data% veri_meta)
    endif
    allocate       (data% veri_data (nb, nv))
    allocate       (data% veri_meta (    nv))
    !-----------------------
    ! read verification data
    !-----------------------
    print *,'reading veri_data:',nv,'out of',data% f% n_veri
    print *,'  ',nv,'out of',data% f% n_veri,'verification entries'
    print *,'  ',nb,'out of',data% f% n_body,'body         entries'
    if (nv == data% f% n_veri .and. nb == data% f% n_body) then
      !---------------------------------------------
      ! no filtering: read verification data at once
      !---------------------------------------------
      call get_var_real_2 (data%f, 'veri_data' ,data% veri_data, fill=fill)
      data% veri_meta = data% f% veri
    else
      !-----------------------------------------------------------------
      ! filtering: 1) read data, store requested data in data% veri_data
      !-----------------------------------------------------------------
      allocate (tmp (data% f% n_body))
      do i = 1, nv
        tmp = REAL(get_veri (data% f, iv(i)), wp)
        if (nb /= data% f% n_body) then
          data% veri_data (:,i) = tmp (data% b% pos)
        else
          data% veri_data (:,i) = tmp (:)
        endif
        data% veri_meta (i) = data% f% veri (iv(i))
      end do
      deallocate (tmp)
      if (present(fill)) then
        fillvalue = REAL(get_fillvalue (data% f, 'veri_data'), wp)
        where (data% veri_data == fillvalue) data% veri_data = fill
      endif
    endif

    !-----------
    ! close file
    !-----------
    call close_fdbk (data% f)

  end subroutine read_fdbk_veri
!------------------------------------------------------------------------------
  subroutine associate_hb (data)
  type (t_fdbk_data) ,intent(inout) :: data
  !--------------------------------------------------------------
  ! associate corresponding header and body entries (set indices)
  ! remove header entries without body   entry
  ! remove body   entries without header entry
  !--------------------------------------------------------------
    integer :: ih   ! header index
    integer :: ib   ! body index
    integer :: nh   ! number of header       entries
    integer :: nb   ! number of body         entries
    integer :: nv   ! number of verification entries
    integer :: ip   ! body index in file
    integer :: ih2  ! header index after removal of unused entries
    integer :: ib2  ! body   index after removal of unused entries
    integer :: i, j ! loop indices

    type (t_fdbk_head) ,pointer :: h(:)   ! temporary header       table
    type (t_fdbk_body) ,pointer :: b(:)   ! temporary body         table
    real (wp)          ,pointer :: v(:,:) ! temporary verification data

    !-------------------------
    ! zero all index variables
    !-------------------------
    if (associated (data% h)) then
      nh  = size (data% h)
      data% h% ib = 0
      data% h% nb = 0
    end if
    if (associated (data% b)) data% b% ih = 0

    !------------------------------------
    ! if both header and body are present
    ! associate corresponding entries
    !------------------------------------
    if (associated (data% h) .and. associated (data% b)) then
      nb  = size (data% b)
      ib  = 1
      ih2 = 1
      ib2 = 1
lh:   do ih = 1, nh               ! loop over header entries
        do                        ! loop over corresponding body entries
          if (ib > nb) exit lh
          ip = data% b(ib)% pos
          if   (ifill == data% h(ih)% i_body .or. &
                ifill == data% h(ih)% l_body      ) then
            write (6,*)                                                      &
              'associate_hb: WARNING fillvalues in i_body, l_body; ih, nh =',&
              ih, nh
            exit lh
          endif
          if   (ip >= data% h(ih)% i_body       ) then
            if (ip >  data% h(ih)% i_body + data% h(ih)% l_body -1) exit
            !-----------------------------------------
            ! corresponding entries found, set indices
            !-----------------------------------------
            if (data% h(ih)% nb == 0) data% h(ih)% ib = ib2
            data% h(ih)% nb = data% h(ih)% nb + 1
            data% b(ib)% ih = ih2
            ib2             = ib2 + 1            ! body was valid
          endif
          ib                = ib + 1             ! end of body loop
        end do
        if (data% h(ih)% nb /= 0) then           ! header was valid
          ih2 = ih2 + 1
        else
          data% h(ih)% l_body = 0
        endif
      end do lh

      !--------------------------------------------
      ! remove unused body and verification entries
      !--------------------------------------------
      ib2 = count(data% b% ih /= 0)
      if (ib2 < nb) then
        nb = ib2
        !----------------------------
        ! remove verification entries
        !----------------------------
        if (associated (data% veri_data)) then
          nv  = size  (data% veri_data,2)
          allocate (v (nb, nv))
          j = 0
          do i = 1, size (data% b)
            if (data% b(i)% ih == 0) cycle
            j = j + 1
            v (j,:) = data% veri_data (i,:)
          end do
          deallocate (data% veri_data)
          data% veri_data => v
        endif
        !--------------------
        ! remove body entries
        !--------------------
        allocate (b (nb))
        b = pack (data% b, mask = data% b% ih /= 0)
        deallocate (data% b)
        data% b => b
      endif
    endif

    if (associated (data% h)) then
      !-----------------------------
      ! remove unused header entries
      !-----------------------------
      ih2 = count(data% h% l_body /= 0)
      if (ih2 < nh) then
        nh = ih2
        allocate (h (nh))
        h = pack (data% h, mask = data% h% l_body /= 0)
        deallocate (data% h)
        data% h => h
      endif
    endif

  end subroutine associate_hb
!------------------------------------------------------------------------------
  subroutine read_fdbk_data (data, file)
  type (t_fdbk_data) ,intent(inout) :: data
  character(len=*)   ,intent(in)    :: file
  !----------------------------------------
  ! read mandatory entries of feedback file
  !----------------------------------------

    call read_fdbk_head (data, file)

    call read_fdbk_body (data)

    call associate_hb   (data)

  end subroutine read_fdbk_data
!------------------------------------------------------------------------------
  subroutine destruct_fdbk_data (data)
  type (t_fdbk_data) ,intent(inout) :: data
  !----------------------------------------------------
  ! deallocate components of feedback file derived type
  !----------------------------------------------------

    if (associated (data% h))         deallocate (data% h)
    if (associated (data% b))         deallocate (data% b)
    if (associated (data% veri_data)) deallocate (data% veri_data)
    if (associated (data% veri_meta)) deallocate (data% veri_meta)
    call cleanup_fdbk (data% f)

  end subroutine destruct_fdbk_data
!------------------------------------------------------------------------------

  subroutine get_var_real (f, name, var, ierr, fill, mask)
  type(t_fdbk)       ,intent(in)  :: f       ! feedback file data type
  character(len=*)   ,intent(in)  :: name    ! variable name
  real(wp)           ,intent(out) :: var (:) ! variable returned
  integer  ,optional ,intent(out) :: ierr    ! error return parameter
  real(wp) ,optional ,intent(in)  :: fill    ! fillvalue to use
  logical  ,optional ,intent(in)  :: mask(:) ! mask to select elements
  !--------------------
  ! get NetCDF variable
  !--------------------

    integer               :: varid, status
    real(wp)              :: fillvalue
    real(wp) ,allocatable :: tmp (:)

    if (present(ierr)) ierr = 0
    varid = get_varid (f, name)
    if (varid < 1) then
      if (present(ierr)) then
        ierr = varid
        return
      else if (present(fill)) then
        var = fill
        return
      else
        call model_abort (-1,-1,'variable not in file: '//name,'get_var_real')
      endif
    endif

    if (present (mask)) then
      allocate (tmp (size(mask)))
      status = nf90_get_var (f% nc% ncid, varid, tmp)
      var = pack (tmp, mask)
    else
      status = nf90_get_var (f% nc% ncid, varid, var)
    endif
    if (status /= NF90_NOERR)                                         &
      call model_abort (-1,-1,                                        &
        'error reading : '//name//' - '//trim(nf90_strerror(status)), &
        'get_var_real'                                                )
    if (present(fill)) then
      status = nf90_get_att (f% nc% ncid, varid, '_FillValue', fillvalue)
      where (var == fillvalue) var = fill
    endif

  end subroutine get_var_real
!------------------------------------------------------------------------------

  subroutine get_var_real_2 (f, name, var, ierr, fill)
  type(t_fdbk)       ,intent(in)  :: f         ! feedback file data type
  character(len=*)   ,intent(in)  :: name      ! variable name
  real(wp)           ,intent(out) :: var (:,:) ! variable returned
  integer  ,optional ,intent(out) :: ierr      ! error return parameter
  real(wp) ,optional ,intent(in)  :: fill      ! fillvalue to use
  !--------------------
  ! get NetCDF variable
  !--------------------

    integer  :: varid, status
    real(wp) :: fillvalue

    if (present(ierr)) ierr = 0
    varid = get_varid (f, name)
    if (varid < 1) then
      if (present(ierr)) then
        ierr = varid
        return
      else if (present(fill)) then
        var = fill
        return
      else
        call model_abort (-1,-1,'variable not in file: '//name,&
                                'get_var_real_2'               )
      endif
    endif

    status = nf90_get_var (f% nc% ncid, varid, var)
    if (status /= NF90_NOERR)                                         &
      call model_abort (-1,-1,                                        &
        'error reading : '//name//' - '//trim(nf90_strerror(status)), &
        'get_var_real_2'                                              )
    if (present(fill)) then
      status = nf90_get_att (f% nc% ncid, varid, '_FillValue', fillvalue)
      where (var == fillvalue) var = fill
    endif

  end subroutine get_var_real_2
!------------------------------------------------------------------------------

  subroutine get_var_int (f, name, var, ierr, fill, mask)
  !--------------------
  ! get NetCDF variable
  !--------------------
  type(t_fdbk)       ,intent(in)  :: f       ! feedback file data type
  character(len=*)   ,intent(in)  :: name    ! variable name
  integer            ,intent(out) :: var (:) ! variable returned
  integer  ,optional ,intent(out) :: ierr    ! error return parameter
  integer  ,optional ,intent(in)  :: fill    ! fillvalue to use
  logical  ,optional ,intent(in)  :: mask(:) ! mask to select elements

    integer              :: varid, status, fillvalue
    integer ,allocatable :: tmp (:)

    if (present(ierr)) ierr = 0
    varid = get_varid (f, name)
    if (varid < 1) then
      if (present(ierr)) then
        ierr = varid
        return
      else if (present(fill)) then
        var = fill
        return
      else
        call model_abort (-1,-1,'variable not in file: '//name,'get_var_int')
      endif
    endif

    if (present (mask)) then
      allocate (tmp (size(mask)))
      status = nf90_get_var (f% nc% ncid, varid, tmp)
      var = pack (tmp, mask)
    else
      status = nf90_get_var (f% nc% ncid, varid, var)
    endif

    if (status /= NF90_NOERR)                                         &
      call model_abort (-1,-1,                                        &
        'error reading : '//name//' - '//trim(nf90_strerror(status)), &
        'get_var_int'                                                 )
    if (present(fill)) then
      status = nf90_get_att (f% nc% ncid, varid, '_FillValue', fillvalue)
      where (var == fillvalue) var = fill
    endif

  end subroutine get_var_int
!------------------------------------------------------------------------------
  subroutine get_var_char (f, name, var)
  type(t_fdbk)      ,intent(in)  :: f       ! feedback file data type
  character(len=*)  ,intent(in)  :: name    ! variable name
  character(len=*)  ,intent(out) :: var (:) ! variable returned
  !--------------------
  ! get NetCDF variable
  !--------------------

    integer                 :: varid, status
    character(len=len(var)) :: tmp (size(var))
    integer                 :: lc (2)
    integer                 :: i,j
    character  ,parameter   :: NUL = achar ( 00)

    varid = get_varid (f, name)
    if (varid < 1) &
      call model_abort (-1,-1,'variable not in file: '//name,'get_var_char')
    lc(1) = len (var)
    lc(2) = size(var)
    status = nf90_get_var (f% nc% ncid, varid, tmp, count=lc)
    if (status /= NF90_NOERR)                                         &
      call model_abort (-1,-1,                                        &
        'error reading : '//name//' - '//trim(nf90_strerror(status)), &
        'get_var_char'                                                )
    !--------------------------------------
    ! convert spurious NUL-character to ' '
    !--------------------------------------
    do j = 1, size (tmp)
      do i=1,len(tmp)
        if (tmp(j)(i:i) == NUL) tmp(j)(i:i) = ' '
      end do
    end do

    var = tmp

   end subroutine get_var_char
!------------------------------------------------------------------------------

  subroutine read_fdbk_ivar_head (fb, name, var, ierr, fill)
  !----------------------------------
  ! read integer variable from header
  !----------------------------------
  type(t_fdbk_data)  ,intent(in)  :: fb      ! feedback file data type
  character(len=*)   ,intent(in)  :: name    ! variable name
  integer            ,intent(out) :: var (:) ! variable returned
  integer  ,optional ,intent(out) :: ierr    ! error return parameter
  integer  ,optional ,intent(in)  :: fill    ! fillvalue to use

    integer :: tmp (fb% f% n_hdr)
    !--------------
    ! read variable
    !--------------
    call get_var_int (fb%f, name, tmp, ierr, fill)
    if (present (ierr)) then
      if (ierr /= 0) return
    endif
    !----------------------------------------------
    ! shrink to observations currently used in 'fb'
    !----------------------------------------------
    var = tmp (fb% h% pos)

  end subroutine read_fdbk_ivar_head

!------------------------------------------------------------------------------

  subroutine read_fdbk_rvar_head (fb, name, var, ierr, fill)
  !-------------------------------
  ! read real variable from header
  !-------------------------------
  type(t_fdbk_data)  ,intent(in)  :: fb      ! feedback file data type
  character(len=*)   ,intent(in)  :: name    ! variable name
  real(wp)           ,intent(out) :: var (:) ! variable returned
  integer  ,optional ,intent(out) :: ierr    ! error return parameter
  real(wp) ,optional ,intent(in)  :: fill    ! fillvalue to use

    real(wp) :: tmp (fb% f% n_hdr)
    !--------------
    ! read variable
    !--------------
    call get_var_real (fb%f, name, tmp, ierr, fill)
    if (present (ierr)) then
      if (ierr /= 0) return
    endif
    !----------------------------------------------
    ! shrink to observations currently used in 'fb'
    !----------------------------------------------
    var = tmp (fb% h% pos)

  end subroutine read_fdbk_rvar_head

!------------------------------------------------------------------------------

  subroutine read_fdbk_ivar_body (fb, name, var, ierr, fill)
  !----------------------------------
  ! read integer variable from header
  !----------------------------------
  type(t_fdbk_data)  ,intent(in)  :: fb      ! feedback file data type
  character(len=*)   ,intent(in)  :: name    ! variable name
  integer            ,intent(out) :: var (:) ! variable returned
  integer  ,optional ,intent(out) :: ierr    ! error return parameter
  integer  ,optional ,intent(in)  :: fill    ! fillvalue to use

    integer :: tmp (fb% f% n_body)
    !--------------
    ! read variable
    !--------------
    call get_var_int (fb%f, name, tmp, ierr, fill)
    if (present (ierr)) then
      if (ierr /= 0) return
    endif
    !----------------------------------------------
    ! shrink to observations currently used in 'fb'
    !----------------------------------------------
    var = tmp (fb% b% pos)

  end subroutine read_fdbk_ivar_body

!------------------------------------------------------------------------------

  subroutine read_fdbk_rvar_body (fb, name, var, ierr, fill)
  !-------------------------------
  ! read real variable from header
  !-------------------------------
  type(t_fdbk_data)  ,intent(in)  :: fb      ! feedback file data type
  character(len=*)   ,intent(in)  :: name    ! variable name
  real(wp)           ,intent(out) :: var (:) ! variable returned
  integer  ,optional ,intent(out) :: ierr    ! error return parameter
  real(wp) ,optional ,intent(in)  :: fill    ! fillvalue to use

    real(wp) :: tmp (fb% f% n_body)
    !--------------
    ! read variable
    !--------------
    call get_var_real (fb%f, name, tmp, ierr, fill)
    if (present (ierr)) then
      if (ierr /= 0) return
    endif
    !----------------------------------------------
    ! shrink to observations currently used in 'fb'
    !----------------------------------------------
    var = tmp (fb% b% pos)

  end subroutine read_fdbk_rvar_body

!------------------------------------------------------------------------------

  subroutine pack_fdbk (fb, mask_h, mask_b)
  !----------------------------------------
  ! pack header and body according to masks
  !----------------------------------------
  type(t_fdbk_data)  ,intent(inout) :: fb         ! feedback file data type
  logical, optional  ,intent(in)    :: mask_h (:) ! mask for header entries
  logical, optional  ,intent(in)    :: mask_b (:) ! mask for body   entries

    if (present (mask_h)) where (.not. mask_h) fb% h% l_body = 0
    if (present (mask_b)) where (.not. mask_b) fb% b% pos    = 0
    call associate_hb (fb)

  end subroutine pack_fdbk

!------------------------------------------------------------------------------

  subroutine write_fdbk_file (fb, file, comment, opt, close)
  !-----------------------------
  ! write complete feedback file
  !-----------------------------
  type(t_fdbk_data) ,intent(inout)        :: fb      ! feedback file data type
  character(len=*)  ,intent(in) ,optional :: file    ! file name
  character(len=*)  ,intent(in) ,optional :: comment ! history extension
  character(len=*)  ,intent(in) ,optional :: opt     ! options: model, obstypes
  logical           ,intent(in) ,optional :: close   ! if .false. dont close

    !----------------
    ! local variables
    !----------------
    integer      :: ir, iv
    logical      :: lclose
    lclose = .true.; if (present(close)) lclose = close

    !------------------------------------------------
    ! set up derived type, set file name, create file
    !------------------------------------------------
    fb% f% n_hdr    = size (fb% h)
    fb% f% n_body   = size (fb% b)
    if (present (file)) fb% f% nc% path = file
    where (fb%f%nc% dims% name=='d_hdr' ) fb%f%nc% dims% len = fb% f% n_hdr
    where (fb%f%nc% dims% name=='d_body') fb%f%nc% dims% len = fb% f% n_body
    call create_netcdf_file (fb% f%nc ,opt=opt)

    !------------------------------------
    ! write history and global attributes
    !------------------------------------
    if (present (comment)) then
      if (comment /= '') then
        call add_history (fb% f, '', '', comment, write = .false.)
      endif
    endif
    call write_global_attributes (fb% f)
    ir = nf90_enddef (fb% f% nc% ncid)
    if (ir /= 0) then
      call model_abort(-1,-1,                                            &
                       'ERROR in nf90_enddef: '//trim(nf90_strerror(ir)),&
                       'write_fdbk_file'                                 )
    endif
    fb% f% nc% defmode = .false.
    !--------------------
    ! loop over variables
    !--------------------
    do iv = 1, fb% f% nc% nvar
      if (.not. fb% f% nc% vars(iv)% opt_used) cycle
      select case (fb% f% nc% vars(iv)% name)
      !-----------------
      ! header variables
      !-----------------
      case ('i_body')
!       call write_fdbk_var (fb% f, iv, fb% h% i_body)
        call write_fdbk_var (fb% f, iv, fb% h% ib)
!
!print *,fb% h(1)% i_body, fb% h(1)% ib, fb% h(1)% pos
!print *,fb% h(1)% l_body, fb% h(1)% nb
!
      case ('l_body')
!       call write_fdbk_var (fb% f, iv, fb% h% l_body)
        call write_fdbk_var (fb% f, iv, fb% h% nb)
      case ('n_level')
        call write_fdbk_var (fb% f, iv, fb% h% n_level)
      case ('data_category')
        call write_fdbk_var (fb% f, iv, fb% h% data_category)
      case ('sub_category')
        call write_fdbk_var (fb% f, iv, fb% h% sub_category)
      case ('center')
        call write_fdbk_var (fb% f, iv, fb% h% center)
      case ('sub_center')
        call write_fdbk_var (fb% f, iv, fb% h% sub_center)
      case ('obstype')
        call write_fdbk_var (fb% f, iv, fb% h% obstype)
      case ('codetype')
        call write_fdbk_var (fb% f, iv, fb% h% codetype)
      case ('ident')
        call write_fdbk_var (fb% f, iv, fb% h% ident)
      case ('statid')
        call write_fdbk_var (fb% f, iv, fb% h% statid)
      case ('wsi')
        call write_fdbk_var (fb% f, iv, fb% h% wsi)
      case ('lat')
        call write_fdbk_var (fb% f, iv, fb% h% lat        ,fill=rfill)
      case ('lon')
        call write_fdbk_var (fb% f, iv, fb% h% lon        ,fill=rfill)
      case ('time')
        call write_fdbk_var (fb% f, iv, fb% h% time)
      case ('time_nomi')
        call write_fdbk_var (fb% f, iv, fb% h% time_nomi)
      case ('time_dbase')
        call write_fdbk_var (fb% f, iv, fb% h% time_dbase ,fill=ifill)
      case ('z_station')
        call write_fdbk_var (fb% f, iv, fb% h% z_station  ,fill=ifill)
      case ('z_modsurf')
        call write_fdbk_var (fb% f, iv, fb% h% z_modsurf  ,fill=ifill)
      case ('r_state')
        call write_fdbk_var (fb% f, iv, fb% h% r_state)
      case ('r_flags')
        call write_fdbk_var (fb% f, iv, fb% h% r_flags)
      case ('r_check')
        call write_fdbk_var (fb% f, iv, fb% h% r_check)
      case ('sta_corr')
        call write_fdbk_var (fb% f, iv, fb% h% sta_corr)
      case ('index_x')
        call write_fdbk_var (fb% f, iv, fb% h% index_x    ,fill=ifill)
      case ('index_y')
        call write_fdbk_var (fb% f, iv, fb% h% index_y    ,fill=ifill)
      case ('mdlsfc')
        call write_fdbk_var (fb% f, iv, fb% h% mdlsfc     ,fill=ifill)
      case ('instype')
        call write_fdbk_var (fb% f, iv, fb% h% instype    ,fill=ifill)
      case ('sun_zenit')
        call write_fdbk_var (fb% f, iv, fb% h% sun_zenit  ,fill=rfill)
      !-----------------
      ! optional entries
      !-----------------
      case ('dbkz')
        call write_fdbk_var (fb% f, iv, fb% h% dbkz      ,fill=-1)
      case ('subset')
        call write_fdbk_var (fb% f, iv, fb% h% subset    ,fill=-1)
      case ('record')
        call write_fdbk_var (fb% f, iv, fb% h% record    ,fill=-1)
      case ('source')
        call write_fdbk_var (fb% f, iv, fb% h% source    ,fill=-1)
      case ('index_d')
        call write_fdbk_var (fb% f, iv, fb% h% index_d   ,fill=-1)
      case ('obs_id')
        call write_fdbk_var (fb% f, iv, fb% h% obs_id    ,fill=-1)
!     case ('retrtype')
      case ('phase')
        call write_fdbk_var (fb% f, iv, fb% h% phase     ,fill=-1)
!     case ('flg_1dvar')
!       call write_fdbk_var (fb% f, iv, fb% h% flg_1dvar ,fill=-1)
      case ('flg_cld')
        call write_fdbk_var (fb% f, iv, fb% h% flg_cld   ,fill=-1)
      case ('surftype')
        call write_fdbk_var (fb% f, iv, fb% h% surftype  ,fill=-1)
      case ('soiltype')
        call write_fdbk_var (fb% f, iv, fb% h% soiltype  ,fill=-1)
      case ('tracking')
        call write_fdbk_var (fb% f, iv, fb% h% tracking , fill=-1)
      case ('meas_type')
        call write_fdbk_var (fb% f, iv, fb% h% meas_type, fill=-1)
      case ('rad_corr')
        call write_fdbk_var (fb% f, iv, fb% h% rad_corr , fill=-1)
      case ('sat_zenit')
        call write_fdbk_var (fb% f, iv, fb% h% sat_zenit ,fill=-999._wp)
      case ('sat_azimuth')
        call write_fdbk_var (fb% f, iv, fb% h% sat_azimuth ,fill=-999._wp)
      case ('sun_azimuth')
        call write_fdbk_var (fb% f, iv, fb% h% sun_azimuth ,fill=-999._wp)
      case ('scanline')
        call write_fdbk_var (fb% f, iv, fb% h% scanline    ,fill=-1)
      case ('fov')
        call write_fdbk_var (fb% f, iv, fb% h% fov         ,fill=-1)
      !-------------
      ! body entries
      !-------------
      case ('varno')
        call write_fdbk_var (fb% f, iv, fb% b% varno)
      case ('obs')
        call write_fdbk_var (fb% f, iv, fb% b% obs    ,fill=rfill)
      case ('bcor')
        call write_fdbk_var (fb% f, iv, fb% b% bcor)
      case ('level')
        call write_fdbk_var (fb% f, iv, fb% b% level ,fill=-99999._wp)
      case ('plevel')
        call write_fdbk_var (fb% f, iv, fb% b% plevel ,fill=-1._wp)
      case ('level_typ')
        call write_fdbk_var (fb% f, iv, fb% b% level_typ)
      case ('level_sig')
        call write_fdbk_var (fb% f, iv, fb% b% level_sig ,fill=-99)
      case ('state')
        call write_fdbk_var (fb% f, iv, fb% b% state)
      case ('flags')
        call write_fdbk_var (fb% f, iv, fb% b% flags)
      case ('check')
        call write_fdbk_var (fb% f, iv, fb% b% check)
      case ('e_o')
        call write_fdbk_var (fb% f, iv, fb% b% e_o  ,fill=rfill)
      case ('qual')
        call write_fdbk_var (fb% f, iv, fb% b% qual ,fill=ifill)
      case ('plev_width')
        call write_fdbk_var (fb% f, iv, fb% b% plev_width, fill=rfill)
      case ('v_loc')
        call write_fdbk_var (fb% f, iv, fb% b% v_loc,fill=rfill)
      case ('h_loc')
        call write_fdbk_var (fb% f, iv, fb% b% h_loc,fill=rfill)
      case ('emissiv')
        call write_fdbk_var (fb% f, iv, fb% b% emiss,fill=rfill)
      !------------------
      ! verification data
      !------------------
      case ('veri_data')
      case ('veri_model')
      case ('veri_run_type')
      case ('veri_run_class')
      case ('veri_initial_date')
      case ('veri_forecast_time')
      case ('veri_resolution')
      case ('veri_domain_size')
      case ('veri_description')
      case ('veri_ens_member')
      case ('veri_exp_id')
      case default
!       call model_abort (-1,-1,                                        &
!                         'not implemented: '//fb% f%nc%vars(iv)% name, &
!                         'write_fdbk_file'                             )
      end select
    end do

    !-------------------------------------
    ! close the file, deallocate meta data
    !-------------------------------------
    if (lclose) call close_fdbk   (fb% f)
!   call cleanup_fdbk (fb% f)

  end subroutine write_fdbk_file
!==============================================================================
  subroutine scatter_fdbk_data (data, src, com, npe, p_pe)
  !---------------------------------
  ! scatter derived type t_fdbk_data
  !---------------------------------
  type (t_fdbk_data) ,intent(inout) :: data
  integer            ,intent(in)    :: src
  integer            ,intent(in)    :: com
  integer            ,intent(in)    :: npe
  integer            ,intent(in)    :: p_pe

    integer           :: nh(0:npe-1), nb(0:npe-1)
    integer           :: ih(0:npe-1), ib(0:npe-1)
    type(t_fdbk_data) :: tmp
    integer           :: i, pe, nv
    integer           :: nbi, ib1, ibn, ob1, obn

    if (src==p_pe) then
      !-----------
      ! count data
      !-----------
      nh = 0
      nb = 0
      do i = 1, size(data% h)
        pe = data% h(i)% pe
        if (pe < 0) cycle
        nh(pe) = nh(pe) + 1
        nb(pe) = nb(pe) + data% h(i)% nb
      end do

      !-------------
      ! reorder data
      !-------------
      ib(0) = 0
      ih(0) = 0
      do pe = 1, npe-1
        ib(pe) = ib(pe-1) + nb(pe-1)
        ih(pe) = ih(pe-1) + nh(pe-1)
      end do
      nv = 0
      if (associated (data% veri_data)) nv = size (data% veri_data,2)
      allocate (tmp% h         (ih(npe-1)+nh(npe-1)   ))
      allocate (tmp% b         (ib(npe-1)+nb(npe-1)   ))
      allocate (tmp% veri_data (ib(npe-1)+nb(npe-1),nv))

      do i = 1, size(data% h)
        pe  = data% h(i)% pe
        if (pe < 0) cycle
        nbi = data% h(i)% nb
        ib1 = data% h(i)% ib
        ibn = ib1 + nbi -1
        ob1 = ib(pe) + 1
        obn = ib(pe) + nbi
        tmp% h           (ih(pe)+1   ) = data% h         (i          )
        tmp% b           (ob1:obn    ) = data% b         (ib1:ibn    )
        if (nv > 0)                                                  &
          tmp% veri_data (ob1:obn,:nv) = data% veri_data (ib1:ibn,:nv)
        ih(pe) = ih(pe) + 1
        ib(pe) = ib(pe) + nbi
      end do

      !---------------------------------
      ! deallocate data on src processor
      !---------------------------------
      deallocate (data% h)
      deallocate (data% b)
      if (associated (data% veri_data)) deallocate (data% veri_data)
    endif

    !----------------------
    ! broadcast array sizes
    !----------------------
    call bcast_int  (nv, src)
    call bcast_int1 (nh, src)
    call bcast_int1 (nb, src)

    !---------------------------------------------
    ! allocate destination derived type components
    !---------------------------------------------
    allocate     (data% h         (nh(p_pe)   ))
    allocate     (data% b         (nb(p_pe)   ))
    allocate     (data% veri_data (nb(p_pe),nv))

    !----------------------------
    ! allocate source tamporaries
    !----------------------------
    if (p_pe /= src) then
      allocate (data% veri_meta (  nv))
      allocate ( tmp% h         (0   ))
      allocate ( tmp% b         (0   ))
      allocate ( tmp% veri_data (0,nv))
    endif

    !--------------------------------------------
    ! broadcast / scatter derived type components
    !--------------------------------------------
    IF ( nv > 0 ) &
      call   bcast_fdbk_meta   (data% veri_meta,  src         )
    call   scatter_fdbk_head (tmp%  h, data% h, src, com, nh)
    call   scatter_fdbk_body (tmp%  b, data% b, src, com, nb)
    do i = 1, nv
      call scatter_fdbk_veri (tmp% veri_data(:,i),             &
                             data% veri_data(:,i), src, com, nb)
    end do

    !----------------
    ! adjust pointers
    !----------------
    ib1 = 1
    do i = 1, size(data% h)
      nbi                        = data% h(i)% nb
      data% h(i            )% ib = ib1
      data% b(ib1:ib1+nbi-1)% ih = i
      ib1                        = ib1 + nbi
    end do

    !-----------------------
    ! deallocate temporaries
    !-----------------------
    deallocate (tmp% h)
    deallocate (tmp% b)
    deallocate (tmp% veri_data)

  end subroutine scatter_fdbk_data
!------------------------------------------------------------------------------
! define basic bcast and scatter routines for derived type components
! by including the respective templates
!------------------------------------------------------------------------------
! MPI parameters are provided by this module, dont USE modules in included code
!
#undef   DEBUG
#ifndef  MPI_CHECKSUM   /* Checksumming of MPI communication (3dvar only) */
#define  MO_MPI_SOURCE
#endif
!------------------------------------------------------------------------------
#define  DERIVED  type(t_fdbk_head)
#undef   MPI_TYPE
#define  p_scatter_DERIVED scatter_fdbk_head
#include "p_scatter_derived.incf"
#undef   DERIVED
#undef   p_scatter_DERIVED
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#define  DERIVED  type(t_fdbk_body)
#undef   MPI_TYPE
#define  p_scatter_DERIVED scatter_fdbk_body
#include "p_scatter_derived.incf"
#undef   DERIVED
#undef   p_scatter_DERIVED
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#define  DERIVED  real(wp)
#undef   MPI_TYPE
#define  p_scatter_DERIVED scatter_fdbk_veri
#include "p_scatter_derived.incf"
#undef   DERIVED
#undef   p_scatter_DERIVED
!------------------------------------------------------------------------------
#define  VECTOR
#define  DERIVED type(t_fdbk_meta),dimension(:)
#define  p_bcast_DERIVED bcast_fdbk_meta
#undef   MPI_TYPE
#include "p_bcast.incf"
#undef   DERIVED
#undef   VECTOR
#undef   p_bcast_DERIVED
#undef   MPI_TYPE
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#define  DERIVED integer
#define  p_bcast_DERIVED bcast_int
#define  MPI_TYPE MPI_INTEGER
#include "p_bcast.incf"
#undef   p_bcast_DERIVED
#undef   MPI_TYPE
#undef   DERIVED
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#define  VECTOR
#define  DERIVED integer,dimension(:)
#define  p_bcast_DERIVED bcast_int1
#include "p_bcast.incf"
#undef   VECTOR
#undef   DERIVED
#undef   MPI_TYPE
#undef   p_bcast_DERIVED
!==============================================================================
end module mo_fdbk_io
