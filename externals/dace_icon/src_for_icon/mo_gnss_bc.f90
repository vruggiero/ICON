!
!+ Bias correction for ground based GNSS (STD / ZTD) data
!
MODULE mo_gnss_bc
!
! Description:
!   Bias correction for ground based GNSS (STD / ZTD) data
!
! Current Code Owner: DWD, Michael Bender
!    phone: +49 69 8062
!    fax:   +49 69 8062 3721
!    email: michael.bender@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_42        2015-06-08 Andreas Rhodin
!  Bias correction for ground based GNSS (STD / ZTD) data
! V1_43        2015-08-19 Michael Bender
!  modifications and error corrections in the STD operator.
! V1_45        2015-12-15 Andreas Rhodin
!  apply bias correction for GNSS ZTD/STD also in COSMO LETKF
! V1_46        2016-02-05 Michael Bender
!  Cleanup
! V1_48        2016-10-06 Robin Faulwetter
!  changed interface for multiple t_decay values
! V1_49        2016-10-25 Harald Anlauf
!  t_spot: generalize bc_airep -> bc_index
! V1_50        2017-01-09 Andreas Rhodin
!  rename modules STD_... to mo_std_...; adapt to 10 character station id
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
!-----------------------------------------------------------

!=============
! Modules used
!=============

  use mo_kind,       only: wp,              &! working precision kind parameter
                           i8                ! 8 byte integer kind parameter
  use mo_exception,  only: finish            ! abort in case of error
  use mo_time,       only: days,            &! derive days (real)
                           operator(-),     &! calculate time difference
                           mjd,             &! Modified Julian Date from time
                           cyyyymmddhhmm,   &! string from time
                           time_cyyyymmddhhmm! time from string
  use mo_obs_set,    only: t_obs_set         ! observation data derived type
  use mo_obs_tables, only: rept_use          ! use table entry
  use mo_dec_matrix, only: t_vector          ! vector data type
  use mo_mpi_dace,   only: dace,            &! MPI group info
                           p_sum,           &! sum over PEs
                           p_bcast,         &! generic MPI broadcast routine
                           p_gather          ! generic MPI gather routine
  use mo_run_params, only: ana_time,        &! analysis time
                           flag_biasc_gpsgb  ! GPSGB bias correction
  use mo_biasc_io,   only: t_bcor_head,     &! file header derived type
                           new_bc_head,     &! construct new file header
                           open_bc_read,    &! open  file for reading
                           open_bc_write,   &! open  file for writing
                           close_bc,        &! close file
                           bc_paths,        &! full pathnames
                           bc_obstyp,       &! observation type in files
                           nbcf,            &! number of files in list
                           verbose           ! Verbosity level of bias corr.
  use mo_fdbk_tables,only: OT_GPSGB,        &! GPSGB report type ID
                           VN_SPD,          &! observed quantity: slant path delay
                           VN_ZPD            ! observed quantity: zenith path delay

  use mo_t_use,      only: CHK_BIASCOR,     &! flag for no bias correction
                           decr_use,        &! decrease the state of a datum
                           STAT_PASSIVE,    &! observation status flags
                           STAT_ACTIVE_0I,  &!
                           STAT_ACCEPTED     !
  use mo_t_obs,      only: stat_hash         ! hash function for 'statid'
  use netcdf,        only:                  &! NetCDF f90 interface
                           nf90_def_dim,    &! define dimensions
                           nf90_def_var,    &! define variables
                           nf90_put_att,    &! define attribute
                           nf90_enddef,     &! end of definition mode
                           nf90_put_var,    &! write variable
                           nf90_strerror,   &! derive error character string
                           NF90_FLOAT,      &! float     type id
                           NF90_INT,        &! integer   type id
                           NF90_CHAR,       &! character type id
                           NF90_GLOBAL,     &! global attribute flag
                           NF90_NOERR,      &! NetCDF return code for no error
                           NF90_FILL_FLOAT   ! NetCDF fillvalue
  use mo_t_netcdf,   only: ncid,            &! NetCDF file id to use
                           stanc,           &! NetCDF start  parameter to use
                           counc,           &! NetCDF count  parameter to use
                           strnc,           &! NetCDF stride parameter to use
                           get_var,         &! get variable
                           get_dim           ! get dimension
  use slatec_module, only: sort              ! sort routine

  use mo_physics,    only: r2d,             &! factor radians -> degree
                           d2r               ! factor degree  -> radians

  use mo_std_coord,  only: gmf               ! Global Mapping Function GMF


  implicit none

!================
! Public entities
!================

  private
  !--------------
  ! derived types
  !--------------
  public :: t_gnss_bc        ! gnss bias correction data
  public :: t_bc             ! component of t_gnss_bc
  public :: t_plat_bc        ! component of t_gnss_bc
  !------------
  ! subroutines
  !------------
  public :: gnss_bc_init       ! initialize module: read namelist, biascor.coeff.
  public :: gnss_bc_bfg        ! apply bias correction, called before fg-check
  public :: gnss_bc_afg        ! update bias correction coefs., after first guess
  public :: gnss_bc_aan        ! update bias correction coefs., after analysis
  public :: read_gnss_bc_file  ! read gnss bias correction file
  !--------------------------------------------------
  ! namelist parameters to be set in namelist STD_OBS
  !--------------------------------------------------
  public :: biascor_mode     ! mode used for updating
  public :: t_decay          ! accumulation decay time (days)
  public :: n_required       ! number of entries required for correction
  public :: bc_fallback      ! fallback if biasc-file not present
  public :: BC_NOBC,BC_UP,BC_FG,BC_AN,BC_VARBC,BC_INVAL ! values for biascor_mode

!=========================
! Derived type definitions
!=========================

  !----------------------------------------
  ! station/center/processing specific data
  !----------------------------------------
  type t_bc
    integer  :: nb        = 0      ! number of entries (obs-bg)
    integer  :: na        = 0      ! number of entries (obs-ana)
    real(wp) :: ob        = 0._wp  ! obs-bg  mean deviation
    real(wp) :: oa        = 0._wp  ! obs-ana mean deviation
    real(wp) :: ob_ob     = 0._wp  ! obs-bg  variance
    real(wp) :: oa_oa     = 0._wp  ! obs-ana variance
    real(wp) :: nb_ac     = 0._wp  ! number of accumulated entries
    real(wp) :: na_ac     = 0._wp  ! number of accumulated entries
    real(wp) :: ob_ac     = 0._wp  ! obs-bg  mean     accumulated
    real(wp) :: oa_ac     = 0._wp  ! ana-bg  mean     accumulated
    real(wp) :: ob_ob_ac  = 0._wp  ! obs-bg  variance accumulated
    real(wp) :: oa_oa_ac  = 0._wp  ! ana-bg  variance accumulated
    real(wp) :: o_err     = 0._wp  ! nominal obs error
    real(wp) :: b_err     = 0._wp  ! nominal bg  error
    real(wp) :: o_err_ac  = 0._wp  ! nominal obs error
    real(wp) :: b_err_ac  = 0._wp  ! nominal bg  error
    real(wp) :: bc_b      = 0._wp  ! bias correction background value
    real(wp) :: bc_a      = 0._wp  ! bias correction analysis   value
!   real(wp) :: bc_b_err           ! bias correction background error
!   real(wp) :: bc_a_err           ! bias correction analysis   error
  end type t_bc

  !--------------------------
  ! observed quantity indices
  !--------------------------
  integer, parameter :: nob   = 1     ! number of observed quantities

  !-----------------------
  ! platform specific data
  !-----------------------
  type t_plat_bc
    character(len=10):: statid  = ''   ! gnss station id
    integer  (i8)    :: hash    = 0    ! integer representation of name
    integer          :: center  = 0    ! processing center
    integer          :: product = 0    ! product
  end type t_plat_bc

  !---------------------
  ! bias correction data
  !---------------------
  type t_gnss_bc
    type(t_bcor_head)        :: h                 ! file header data
    integer                  :: biascor_mode = 0  ! mode used for updating
    integer                  :: n            = 0  ! number of entries
    type (t_plat_bc),pointer :: plat (:)          ! platform specific metadata
    type (t_bc)     ,pointer :: data (:,:)      ! (nph,nob,:) data
  end type t_gnss_bc

  !-------------------------------------
  ! derived type used for cross-checking
  !-------------------------------------
  type t_tmp
    character(len=10):: statid ! station name
    integer  (i8)    :: hash   ! integer representation of name
    integer          :: center ! processing center
    integer          :: product! product
    integer          :: pe     ! processor       of report
    integer          :: i      ! old index
    integer          :: j      ! new index
  end type t_tmp

  !------------------------
  ! values for biascor_mode
  !------------------------
  integer ,parameter :: BC_INVAL = -9  ! not set so far
  integer ,parameter :: BC_NOBC  =  0  ! no bias correction
  integer ,parameter :: BC_UP    =  1  ! only update bias corrrection file
  integer ,parameter :: BC_FG    =  2  ! apply  bias corr. from first guess
  integer ,parameter :: BC_AN    =  3  ! apply  bias corr. from analysis
  integer ,parameter :: BC_VARBC =  4  ! variational bias correction

!=================
! Module variables
!=================

  type(t_gnss_bc),save :: gnss_bc              ! bias correction data
  type(t_gnss_bc),save :: empty_bc             ! empty bias correction
  !-----------------
  ! namelist entries
  !-----------------
  integer                  :: biascor_mode = BC_INVAL  ! mode used for updating
  real(wp)                 :: t_decay      =   -30._wp ! accumulation decaytime
  integer                  :: n_required   =    50     ! # of entries required
  logical                  :: bc_fallback  = .false.   ! biasc-file not present

contains
!==============================================================================

  subroutine gnss_bc_init
  !---------------------------------------
  ! Initialize this module:
  ! read bias-correction coefficient files
  !---------------------------------------

    integer :: i    ! loop index
    integer :: ierr ! error return value


    !--------------------------------------------------------------------
    ! set defaults for namelist parameters depending on 'ga3_biasc_gpsgb'
    ! if not set so far
    !--------------------------------------------------------------------
    if (biascor_mode == BC_INVAL) then
      select case (flag_biasc_gpsgb)
      case (-1)
        biascor_mode = BC_NOBC
        bc_fallback  = .false.
      case ( 0)
        biascor_mode = BC_FG
        bc_fallback  = .true.
      case ( 1)
        biascor_mode = BC_FG
        bc_fallback  = .false.
      end select
    endif

    if (biascor_mode /= 0) then
      !---------------------------------------
      ! read bias-correction coefficient files
      !---------------------------------------
      ierr = -1
      if (flag_biasc_gpsgb /= 0 .or. .not. bc_fallback) then
        do i = 1, nbcf
          if (bc_obstyp(i) == OT_GPSGB) then
            call read_gnss_bc_file (gnss_bc, bc_paths(i))
            ierr = 0
            exit
          endif
        end do
      end if
      !-----------------------------------------------
      ! create empty bias-correction coefficient files
      !-----------------------------------------------
      if (ierr /= 0) then
        if (.not.bc_fallback) &
          call finish ('gnss_bc_init','GPSGB: Coefficient file not present !')
        call new_bc_head (gnss_bc% h, OT_GPSGB, t_decay=(/abs(t_decay)/))
        gnss_bc% biascor_mode = biascor_mode
        gnss_bc% n            = 0
        allocate (gnss_bc% plat (0))
        allocate (gnss_bc% data (nob,0))
        if (dace% lpio) write(6,'(a,a)')'    created empty file ',&
                                             trim (gnss_bc%h% path)
      endif
      !-------------------------------
      ! prepare bias correction to use
      !-------------------------------
      select case (biascor_mode)
      case (BC_UP)
        gnss_bc% data% bc_b = 0
        gnss_bc% data% bc_a = 0
      case (BC_FG)
        where (gnss_bc% data% nb_ac >= n_required)
          gnss_bc% data% bc_b = gnss_bc% data% ob_ac
          gnss_bc% data% bc_a = gnss_bc% data% ob_ac
        elsewhere
          gnss_bc% data% bc_b = NF90_FILL_FLOAT
          gnss_bc% data% bc_a = NF90_FILL_FLOAT
        endwhere
      case (BC_AN)
        where (gnss_bc% data% na_ac >= n_required)
          gnss_bc% data% bc_b = gnss_bc% data% oa_ac
          gnss_bc% data% bc_a = gnss_bc% data% oa_ac
        elsewhere
          gnss_bc% data% bc_b = NF90_FILL_FLOAT
          gnss_bc% data% bc_a = NF90_FILL_FLOAT
        endwhere
      case (BC_VARBC)
        gnss_bc% data% bc_b = gnss_bc% data% bc_a
      end select
    endif

  end subroutine gnss_bc_init

!------------------------------------------------------------------------------

  subroutine read_gnss_bc_file (bc, file, inst)
  !-----------------------------------
  ! read gnss bias correction file
  !-----------------------------------
  type(t_gnss_bc) ,intent(out) :: bc    ! derived type variable to fill
  character(len=*)    ,intent(in)  :: file  ! name of file to read
  logical   ,optional ,intent(in)  :: inst  ! read instantanious data as well

    integer :: i
    logical :: linst

    linst = .false.; if (present(inst)) linst = inst

    if (dace% lpio) then
      !-------------------------------
      ! open NetCDF file, read header
      !-------------------------------
      bc% h% path = file
      call open_bc_read (bc% h)

      !----------------
      ! read dimensions
      !----------------
      ncid  = bc% h% ncid
      call get_dim (bc% n ,'platform')

      !--------------------
      ! allocate components
      !--------------------
      allocate (bc% plat (    bc% n))
      allocate (bc% data (nob,bc% n))

      !---------------
      ! read variables
      !---------------
      if (bc% n > 0) then
        stanc = 1
        counc = 0
        strnc = 1
        call get_var (bc% plat% statid   ,'statid'  )
        call get_var (bc% plat% center   ,'center'  )
        call get_var (bc% plat% product  ,'product' )
        call get_var (bc% data% nb_ac    ,'nb_ac'   )
        call get_var (bc% data% na_ac    ,'na_ac'   )
        call get_var (bc% data% ob_ac    ,'ob_ac'   )
        call get_var (bc% data% oa_ac    ,'oa_ac'   )
        call get_var (bc% data% ob_ob_ac ,'ob_ob_ac')
        call get_var (bc% data% oa_oa_ac ,'oa_oa_ac')
        call get_var (bc% data% o_err_ac ,'o_err_ac')
        call get_var (bc% data% b_err_ac ,'b_err_ac')
        call get_var (bc% data% bc_a     ,'bc_a'    )
        if (linst) then
          call get_var (bc% data% nb    ,'nb'   )
          call get_var (bc% data% na    ,'na'   )
          call get_var (bc% data% ob    ,'ob'   )
          call get_var (bc% data% oa    ,'oa'   )
          call get_var (bc% data% ob_ob ,'ob_ob')
          call get_var (bc% data% oa_oa ,'oa_oa')
          call get_var (bc% data% o_err ,'o_err')
          call get_var (bc% data% b_err ,'b_err')
        endif
!       bc% plat% hash = transfer  (bc% plat% statid, bc% plat% hash)
        bc% plat% hash = stat_hash (bc% plat% statid)
!       do i = 1, bc% n
!         bc% plat(i)% hash = transfer  (bc% plat(i)% statid, bc% plat(i)% hash)
!         bc% plat(i)% hash = stat_hash (bc% plat(i)% statid)
!       end do
      endif
      !-----------
      ! close file
      !-----------
      call close_bc (bc% h)
    endif

    !----------
    ! broadcast
    !----------
    call p_bcast (bc% n, dace% pio)
    if (.not.dace% lpio) then
      allocate (bc% plat (    bc% n))
      allocate (bc% data (nob,bc% n))
    endif
    call p_bcast_plat (gnss_bc% plat, dace% pio)
    call p_bcast_data (gnss_bc% data, dace% pio)

  end subroutine read_gnss_bc_file

!------------------------------------------------------------------------------

  subroutine destruct_bc_file (bc)
  type(t_gnss_bc) ,intent(inout) :: bc    ! derived type variable to fill

    !----------------------
    ! deallocate components
    !----------------------
    deallocate (bc% plat)
    deallocate (bc% data)

    bc = empty_bc
  end subroutine destruct_bc_file

!------------------------------------------------------------------------------

  subroutine write_gnss_bc_file (bc)
  !------------------------------------
  ! write gnss bias correction file
  !------------------------------------
  type(t_gnss_bc) ,intent(inout) :: bc

    integer :: status                            ! NetCDF return value
    integer :: dimid (2)                         ! NetCDF dimension id
    integer :: dimch                             ! NetCDF dimension id
    integer :: statid                            ! NetCDF variable ids
    integer :: center                            ! NetCDF variable ids
    integer :: product                           ! NetCDF variable ids
    integer :: nb, na                            ! .
    integer :: ob, oa                            ! .
    integer :: ob_ob, oa_oa                      ! .
    integer :: nb_ac, na_ac                      ! .
    integer :: ob_ac, oa_ac                      ! .
    integer :: ob_ob_ac, oa_oa_ac                ! .
    integer :: bc_b, b_err , b_err_ac ! bc_b_err ! .
    integer :: bc_a, o_err , o_err_ac ! bc_a_err ! .
    !++++++++++++++++++++++++++++++++++++++++++
    ! workaround for NEC SX9 nf90_put_var(text)
    !++++++++++++++++++++++++++++++++++++++++++
    character(len=10),allocatable :: statids(:)

    if (dace% lpio .and. bc%n > 0) then
      !-------------------------------
      ! open NetCDF file, write header
      !-------------------------------
      bc% h% t_decay(2:) = bc% h% t_decay(1) ! Write only one the first t_decay value
      call open_bc_write (bc% h, OT_GPSGB)

      !---------------
      ! set attributes
      !---------------
      status = nf90_put_att (bc%h% ncid, NF90_GLOBAL, 'biascor_mode', &
                                                       biascor_mode   )
      !------------------
      ! define dimensions
      !------------------
      call def_dim ('observation', nob,            dimid(1))
      call def_dim ('platform',    size(bc% plat), dimid(2))
      call def_dim ('chars',      10,              dimch   )

      !-----------------
      ! define variables
      !-----------------
      status = nf90_def_var (bc%h% ncid ,'statid' ,NF90_CHAR, &
                             (/dimch,dimid(2)/), statid)
      status = nf90_put_att (bc%h% ncid , statid, 'longname',&
                             'station id as character string')
      call def_var1('center'    ,center    ,'center')
      call def_var1('product'   ,product   ,'product')
      call def_var ('nb'        ,nb        ,'entries obs-bg')
      call def_var ('na'        ,na        ,'entries obs-ana')
      call def_var ('ob'        ,ob        ,'obs-bg  mean')
      call def_var ('oa'        ,oa        ,'obs-ana mean')
      call def_var ('ob_ob'     ,ob_ob     ,'obs-bg  variance')
      call def_var ('oa_oa'     ,oa_oa     ,'obs-ana variance')
      call def_var ('nb_ac'     ,nb_ac     ,'entries obs-bg accumulated')
      call def_var ('na_ac'     ,na_ac     ,'entries obs-ana accumulated')
      call def_var ('ob_ac'     ,ob_ac     ,'obs-bg  mean accumulated')
      call def_var ('oa_ac'     ,oa_ac     ,'obs-ana mean accumulated')
      call def_var ('ob_ob_ac'  ,ob_ob_ac  ,'obs-bg  variance accumulated')
      call def_var ('oa_oa_ac'  ,oa_oa_ac  ,'obs-ana variance accumulated')
      call def_var ('bc_b'      ,bc_b      ,'first guess bias correction')
      call def_var ('bc_a'      ,bc_a      ,'analysis bias correction')
      call def_var ('o_err'     ,o_err     ,'nominal obs error squared')
      call def_var ('b_err'     ,b_err     ,'nominal bg  error squared')
      call def_var ('o_err_ac'  ,o_err_ac  ,               &
                    'nominal obs error squared accumulated')
      call def_var ('b_err_ac'  ,b_err_ac  ,               &
                    'nominal bg  error squared accumulated')
!     call def_var ('bc_b'      ,bc_b      ,'bias correction background value')
!     call def_var ('bc_a'      ,bc_a      ,'bias correction analysis   value')
!     call def_var ('bc_b_err ' ,bc_b_err  ,                  &
!                   'bias correction background error squared')
!     call def_var ('bc_a_err ' ,bc_a_err  ,                 &
!                   'bias corrrection analysis error squared')
      status = nf90_enddef  (bc%h% ncid)
      !----------------
      ! write variables
      !----------------
      allocate (statids (bc% n))
      statids = bc% plat% statid
!+++  status  = nf90_put_var (bc%h% ncid, statid, bc% plat% statid)
      status  = nf90_put_var (bc%h% ncid, statid, statids)
      deallocate (statids)
      call put_var1 ('center'   ,center   ,bc% plat% center)
      call put_var1 ('product'  ,product  ,bc% plat% product)
      call put_vari ('nb'       ,nb       ,bc% data% nb)
      call put_vari ('na'       ,na       ,bc% data% na)
      call put_var  ('ob'       ,ob       ,bc% data% ob)
      call put_var  ('oa'       ,oa       ,bc% data% oa)
      call put_var  ('ob_ob'    ,ob_ob    ,bc% data% ob_ob)
      call put_var  ('oa_oa'    ,oa_oa    ,bc% data% oa_oa)
      call put_var  ('nb_ac'    ,nb_ac    ,bc% data% nb_ac)
      call put_var  ('na_ac'    ,na_ac    ,bc% data% na_ac)
      call put_var  ('ob_ac'    ,ob_ac    ,bc% data% ob_ac)
      call put_var  ('oa_ac'    ,oa_ac    ,bc% data% oa_ac)
      call put_var  ('ob_ob_ac' ,ob_ob_ac ,bc% data% ob_ob_ac)
      call put_var  ('oa_oa_ac' ,oa_oa_ac ,bc% data% oa_oa_ac)
      call put_var  ('bc_b'     ,bc_b     ,bc% data% bc_b)
      call put_var  ('bc_a'     ,bc_a     ,bc% data% bc_a)
      call put_var  ('o_err '   ,o_err    ,bc% data% o_err )
      call put_var  ('o_err_ac' ,o_err_ac ,bc% data% o_err_ac)
      call put_var  ('b_err '   ,b_err    ,bc% data% b_err )
      call put_var  ('b_err_ac' ,b_err_ac ,bc% data% b_err_ac)

      !-----------
      ! close file
      !-----------
      call close_bc (bc% h)
    endif

  contains
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    subroutine put_var1 (name, varid, values)
    character(len=*) ,intent(in) :: name
    integer          ,intent(in) :: varid
    integer          ,intent(in) :: values (:)

      status = nf90_put_var (bc%h% ncid ,varid, values)

      if (status /= NF90_NOERR) then
        call finish ('write_gnss_bc_file: put_var '//trim(name),&
                      trim(nf90_strerror(status))            )
      end if

    end subroutine put_var1
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    subroutine put_vari (name, varid, values)
    character(len=*) ,intent(in) :: name
    integer          ,intent(in) :: varid
    integer          ,intent(in) :: values (:,:)

      status = nf90_put_var (bc%h% ncid ,varid, values)

      if (status /= NF90_NOERR) then
        call finish ('write_gnss_bc_file: put_var '//trim(name),&
                      trim(nf90_strerror(status))            )
      end if

    end subroutine put_vari
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    subroutine put_var (name, varid, values)
    character(len=*) ,intent(in) :: name
    integer          ,intent(in) :: varid
    real(wp)         ,intent(in) :: values (:,:)

      status = nf90_put_var (bc%h% ncid ,varid, values)

      if (status /= NF90_NOERR) then
        call finish ('write_gnss_bc_file: put_var '//trim(name),&
                      trim(nf90_strerror(status))            )
      end if

    end subroutine put_var
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    subroutine def_dim (name, size, dimid)
    character(len=*) ,intent(in)  :: name
    integer          ,intent(in)  :: size
    integer          ,intent(out) :: dimid

      status = nf90_def_dim (bc%h% ncid ,name ,size ,dimid)

      if (status /= NF90_NOERR) then
        write(0,*)   'write_gnss_bc_file: def_dim '//trim(name)//' : size =',size
        call finish ('write_gnss_bc_file: def_dim '//trim(name),&
                      trim(nf90_strerror(status))            )
      end if

    end subroutine def_dim
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    subroutine def_var (name, varid, longname)
    character(len=*) ,intent(in)  :: name
    integer          ,intent(out) :: varid
    character(len=*) ,intent(in)  :: longname

      status = nf90_def_var (bc%h% ncid ,name ,NF90_FLOAT, dimid, varid)

      if (status /= NF90_NOERR) then
        call finish ('write_gnss_bc_file: def_var '//trim(name),&
                      trim(nf90_strerror(status))            )
      end if

      status = nf90_put_att (bc%h% ncid ,varid, 'longname', longname)

      if (status /= NF90_NOERR) then
        call finish ('write_gnss_bc_file: def_var '//trim(name),&
                      trim(nf90_strerror(status))            )
      end if
    end subroutine def_var
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    subroutine def_var1 (name, varid, longname)
    character(len=*) ,intent(in)  :: name
    integer          ,intent(out) :: varid
    character(len=*) ,intent(in)  :: longname

      status = nf90_def_var (bc%h% ncid ,name ,NF90_INT, dimid(2), varid)

      if (status /= NF90_NOERR) then
        call finish ('write_gnss_bc_file: def_var1 '//trim(name),&
                      trim(nf90_strerror(status))                )
      end if

      status = nf90_put_att (bc%h% ncid ,varid, 'longname', longname)

      if (status /= NF90_NOERR) then
        call finish ('write_gnss_bc_file: def_var1 '//trim(name),&
                      trim(nf90_strerror(status))                )
      end if
    end subroutine def_var1
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine write_gnss_bc_file

!==============================================================================

  subroutine gnss_bc_bfg (obs)
  !------------------------------------------------------------------------
  ! Radiance bias correction routine to be called before first guess check.
  ! Check for missing entries in bc file, extend file
  ! Apply bias correction to observations.
  !------------------------------------------------------------------------
  type (t_obs_set) ,intent(inout) :: obs  ! observation

    !----------------
    ! local variables
    !----------------
    integer ,parameter       :: mp = 10000   ! max. number of new planes
    type(t_tmp)              :: mis_pe (mp)  ! missing entries in bc file
    type(t_tmp)              :: mis_all(mp)  ! missing entries in bc file
    integer(i8)              :: hash         ! integer representation of name
    integer                  :: n            ! counter / PE
    integer                  :: m            ! counter all PEs
    integer                  :: ib           ! box index
    integer                  :: is           ! spot index
    integer                  :: i, j, k      ! loop indices
    type(t_plat_bc) ,pointer :: plat (:)     ! temporary
    type (t_bc)     ,pointer :: data (:,:)   ! temporary
    real(wp)                 :: rawobs       ! raw observation
    real(wp)                 :: bcor         ! bias correction
    integer                  :: ibc          ! bias correction index
    integer                  :: state        ! state to set if no BC available
    real(wp)                 :: gmfh         ! GMF: hydrostatic mapping function
    real(wp)                 :: gmfw         ! GMF: wet mapping function
    real(wp), parameter      :: pi05 = 1.57079632679489655D0 ! 0.5*Pi = 90 deg.
    !---------------------------------------------------------
    ! relate bias correction file entries to reports in 3dvar.
    ! check for missing entries in bc file.
    ! extend file if necessary.
    !---------------------------------------------------------

    if (biascor_mode /= 0) then
      !-----------------------------------------------
      ! relate bias correction file entries to reports
      ! check for missing entries on this PE
      !-----------------------------------------------
      n = 0
      do ib = 1, size (obs% o)
        if (obs% o(ib)% pe /= dace% pe) cycle
spot:   do is = 1, obs% o(ib)% n_spot
        if (obs% o(ib)% spot(is)% hd% obstype /= OT_GPSGB) cycle
          obs% o(ib)% spot(is)% bc_index = 0
          !do i = obs% o(ib)% spot(is)% o% i + 1,                       &
          !       obs% o(ib)% spot(is)% o% i + obs% o(ib)% spot(is)% o% n
!           !+++++++++++++++++++++++++++++++++++++++++++++
!           ! TODO: check for reasonable observed quantity
!           !+++++++++++++++++++++++++++++++++++++++++++++
!           select case (obs% o(ib)% varno(i))
!           case (VN_T, VN_RH)
!           case default
!             cycle
!           end select

            !-----------------------------------------------
            ! check for corresponding station/center/product
            !-----------------------------------------------
            ! hash               - station
            ! center = center_id - center
            ! stret              - product
            !-----------------------------------------------
            hash = stat_hash (obs% o(ib)% spot(is)% statid)
            do j = 1, size(gnss_bc% plat)
              if  (gnss_bc% plat(j)% hash    == hash        .and.       &
                   gnss_bc% plat(j)% center  ==                         &
                            obs% o(ib)% spot(is)% center_id .and.       &
                   gnss_bc% plat(j)% product ==                         &
                            obs% o(ib)% spot(is)% stret          ) then
                obs% o(ib)% spot(is)% bc_index = j
                cycle spot
              endif
            end do
            do j = 1, n
              if (mis_pe(j)% hash    == hash                 .and.       &
                  mis_pe(j)% center  ==                                  &
                             obs% o(ib)% spot(is)% center_id .and.       &
                  mis_pe(j)% product ==                                  &
                             obs% o(ib)% spot(is)% stret          ) then
                obs% o(ib)% spot(is)% bc_index = -j
                cycle spot
              endif
            end do
            !--------------------------------
            ! add to list of missing stations
            !--------------------------------
            n = n + 1
            if (n > mp) call finish('gnss_bc_bfg','n > mp')
            obs% o(ib)% spot(is)% bc_index = - n
            mis_pe (n)% statid  = obs% o(ib)% spot(is)% statid
            mis_pe (n)% hash    = hash
            mis_pe (n)% center  = obs% o(ib)% spot(is)% center_id
            mis_pe (n)% product = obs% o(ib)% spot(is)% stret
            mis_pe (n)% i       = n
            mis_pe (n)% j       = 0
            cycle spot
          !end do
        end do spot
      end do
      !---------------------------------------
      ! cross-check missing entries on all PEs
      !---------------------------------------
      m = p_sum (n)
      if (m > mp) call finish('gnss_bc_bfg','m > mp')
      if (m > 0) then
        call p_gather_tmp (mis_pe (1:n), mis_all(1:m), dace% pio)
        if (dace% lpio) then
          k             = 1
          i             = 1
          mis_all(i)% j = k
          mis_pe (k)    = mis_all(i)
outer:    do i = 2, m
            do j = 1, k
              !-----------------------------------------------
              ! Check for corresponding station/center/product
              !-----------------------------------------------
              if (mis_all(i)% hash    ==  mis_pe(j)% hash    .and.       &
                  mis_all(i)% center  ==  mis_pe(j)% center  .and.       &
                  mis_all(i)% product ==  mis_pe(j)% product      ) then
                mis_all(i)% j = j
                cycle outer
              endif
            end do
            k = k + 1
            mis_all(i)% j = k
            mis_pe (k)    = mis_all(i)
          end do outer
          !------------
          ! extend file
          !------------
          i = size (gnss_bc% plat)
          plat => gnss_bc% plat
          data => gnss_bc% data
          allocate (gnss_bc% plat (    i+k))
          allocate (gnss_bc% data (nob,i+k))
          gnss_bc% n = i+k
          gnss_bc% plat(  :i) = plat
          gnss_bc% data(:,:i) = data
          do j = 1, k
            gnss_bc% plat (i+j)% hash    = mis_pe (j)% hash
            gnss_bc% plat (i+j)% statid  = mis_pe (j)% statid
            gnss_bc% plat (i+j)% center  = mis_pe (j)% center
            gnss_bc% plat (i+j)% product = mis_pe (j)% product
          end do
          deallocate (plat, data)
          mis_all (1:m)% j = mis_all (1:m)% j + i
        endif
        !--------------------------------------
        ! redistribute information to other PEs
        !--------------------------------------
        call p_scatter_tmp (mis_all(1:m), mis_pe(1:n), dace% pio)
        if (dace% lpio) m = i+k
        call p_bcast (m, dace% pio)
        if (.not.dace% lpio) then
          deallocate (gnss_bc% plat)
          deallocate (gnss_bc% data)
          allocate   (gnss_bc% plat (m))
          allocate   (gnss_bc% data (nob,m))
        endif
        call p_bcast_plat (gnss_bc% plat, dace% pio)
        call p_bcast_data (gnss_bc% data, dace% pio)
        gnss_bc% n = m
        !-----------------------------------
        ! relate reports to new file entries
        !-----------------------------------
        do ib = 1, size (obs% o)
          if (obs% o(ib)% pe /= dace% pe) cycle
          do is = 1, obs% o(ib)% n_spot
            if (obs% o(ib)% spot(is)% hd% obstype /= OT_GPSGB) cycle
            if (obs% o(ib)% spot(is)% bc_index < 0)      &
                obs% o(ib)% spot(is)% bc_index = mis_pe( &
               -obs% o(ib)% spot(is)% bc_index)% j
          end do
        end do
      endif
    endif

    !-----------------------------------
    ! apply bias-correction within 3dvar
    !-----------------------------------
    if (biascor_mode >= BC_FG) then
      state = rept_use (OT_GPSGB)% use (CHK_BIASCOR)
      do ib = 1, size (obs% o)
        if (obs% o(ib)% pe /= dace% pe) cycle
        do is = 1, obs% o(ib)% n_spot
          if (obs% o(ib)% spot(is)% hd% obstype /= OT_GPSGB) cycle
          ibc = obs% o(ib)% spot(is)% bc_index
          if (ibc <= 0)                                      cycle
          do i = obs% o(ib)% spot(is)% o% i + 1,                       &
                 obs% o(ib)% spot(is)% o% i + obs% o(ib)% spot(is)% o% n
            !++++++++++++++++++++++++++++++++++++++++++++++
            ! TODO: check for reasonable quantity (ZTD/STD)
            !++++++++++++++++++++++++++++++++++++++++++++++
            j = 1
            rawobs = obs% o(ib)% body(i)% o  - &
                     obs% o(ib)% body(i)% bc
            bcor = gnss_bc% data(j,ibc)% bc_b
            if (bcor == NF90_FILL_FLOAT) then
              bcor = 0._wp
              call decr_use (obs% o(ib)% body(i)% use, &
                             check = CHK_BIASCOR,      &
                             state = state,            &
                             lflag = .true.            )
            endif

            ! The STD bias is almost independend from the elevation
            ! and can be applied without any mapping!
            obs% o(ib)% body(i)% bc =        - bcor
            obs% o(ib)% body(i)% o  = rawobs - bcor

!!$            write(*,*) 'GNSSBIAS ', obs% o(ib)% spot(is)% statid, &
!!$                 obs% o(ib)% varno(i),               &
!!$                 obs% o(ib)% olev(i),                &
!!$                 obs% o(ib)% spot(is)% hd % time,       &
!!$                 bcor, bcor/gmfh, gmfh, rawobs,  rawobs - bcor

          end do
        end do
      end do
    endif
  end subroutine gnss_bc_bfg

!==============================================================================

  subroutine gnss_bc_aan (obs, y)
  !-----------------------------------------------------------------------
  ! Bias correction routine to be called after analysis
  ! Update bias correction statistics.
  ! Write updated correction coefficient file.
  !-----------------------------------------------------------------------
  type (t_obs_set) ,intent(inout) :: obs  ! observation
  type (t_vector)  ,intent(in)    :: y    ! background, observation space

    integer             :: i, j      ! index
!   integer             :: n         ! number of gnsss in statistics
    integer             :: ib        ! box index
    integer             :: is        ! spot index
    integer             :: ibc       ! plane index
    real(wp)            :: an        ! analysis
    real(wp)            :: o         ! raw observation
    real(wp)            :: doa       ! obs - ana
!   real(wp)            :: dt        ! time since last update of statistics
!   real(wp)            :: f         ! weight for accumulated statistics
    integer             :: ier       ! error return flag
    integer,    pointer :: iperm (:) ! permutation index array for sorting
    type(t_bc), pointer :: pp  (:,:) ! pointer to statistics file entries
    type(t_bc), pointer :: p     (:) ! pointer to statistics file entries
    real(wp)            :: gmfh      ! GMF: hydrostatic mapping function
    real(wp)            :: gmfw      ! GMF: wet mapping function
    real(wp), parameter :: pi05 = 1.57079632679489655D0 ! 0.5*Pi = 90 deg.

    ! sorting by station/provider/product
    character (len=16), allocatable  :: sortarr(:) ! array to be sorted
    character (len=3)                :: prdstr     ! string with product ID
    character (len=3)                :: cntstr     ! string with center ID

    if (biascor_mode == 0) return
!   n = size(gnss_bc% plat)

    !-----------------
    ! set sums to zero
    !-----------------
    pp => gnss_bc% data(:,:)
    pp% na    = 0
    pp% oa    = 0._wp  ! obs-ana mean deviation
    pp% oa_oa = 0._wp  ! obs-ana variance

    !-----------------------
    ! Update bias statistics
    !-----------------------
    do ib = 1, size (obs% o)
      if (obs% o(ib)% pe /= dace% pe) cycle
      do is = 1, obs% o(ib)% n_spot
        if (obs% o(ib)% spot(is)% hd% obstype /= OT_GPSGB) cycle
        ibc = obs% o(ib)% spot(is)% bc_index
        if (ibc <= 0)                                      cycle
        do i = obs% o(ib)% spot(is)% o% i + 1,                       &
               obs% o(ib)% spot(is)% o% i + obs% o(ib)% spot(is)% o% n
          j = 1
          p => gnss_bc% data(:,ibc)
          select case (obs% o(ib)% body(i)% use% state)
          case (STAT_PASSIVE, STAT_ACTIVE_0I:STAT_ACCEPTED)
            an = y% s(ib)% x(i)
            o  = obs% o(ib)% body(i)% o - obs% o(ib)% body(i)% bc
            doa  = o  - an

            if (obs% o(ib)% varno(i) == VN_SPD) then
               ! map obsrvations to zenith using the Global Mapping Function GMF
               call gmf( mjd(obs% o(ib)% spot(is)% hd % time),       & ! MJD
                         obs% o(ib)% spot(is)% col % c % dlat * d2r, & ! lat
                         obs% o(ib)% spot(is)% col % c % dlon * d2r, & ! lon
                         obs% o(ib)% spot(is)% z,                    & ! height
                         pi05 - obs% o(ib)% olev(i) * d2r,           & ! Z
                         gmfh, gmfw)
               doa = doa / gmfh
            end if

            p(j)% na     = p(j)% na     + 1
            p(j)% oa     = p(j)% oa     + doa
            p(j)% oa_oa  = p(j)% oa_oa  + doa * doa
          end select
        end do
      end do
    end do

    !----------------
    ! sum up over PEs
    !----------------
    gnss_bc% data% na     = p_sum (gnss_bc% data% na)
    gnss_bc% data% oa     = p_sum (gnss_bc% data% oa)
    gnss_bc% data% oa_oa  = p_sum (gnss_bc% data% oa_oa)

    !-----------------------------------
    ! Update accumulated bias statistics
    !-----------------------------------
    pp => gnss_bc% data(:,:)
    !-------------------------------
    ! rescale accumulated statistics
    !-------------------------------
    pp% oa_ac      = pp% oa_ac    * pp% na_ac
    pp% oa_oa_ac   = pp% oa_oa_ac * pp% na_ac
    !-----------
    ! accumulate
    !-----------
    pp% na_ac      = pp% na_ac    + pp% na
    pp% oa_ac      = pp% oa_ac    + pp% oa
    pp% oa_oa_ac   = pp% oa_oa_ac + pp% oa_oa
    !--------
    ! rescale
    !--------
    where (pp% na > 0)
      pp% oa       = pp% oa       / pp% na
      pp% oa_oa    = pp% oa_oa    / pp% na
    endwhere
    where (pp% na_ac > 0)
      pp% oa_ac    = pp% oa_ac    / pp% na_ac
      pp% oa_oa_ac = pp% oa_oa_ac / pp% na_ac
    endwhere

    !-------------------------------
    ! sort by station/center/product
    !-------------------------------
    allocate (iperm  (gnss_bc% n))
    allocate (sortarr(gnss_bc% n))
    sortarr = gnss_bc% plat% statid
    do i=1, gnss_bc% n
       write(cntstr,'(i3.3)') gnss_bc% plat(i)% center
       write(prdstr,'(i3.3)') gnss_bc% plat(i)% product
       sortarr(i) = sortarr(i)(1:10) // cntstr // prdstr
    end do
    call sort (sortarr, iperm, 1, ier)
    gnss_bc% plat = gnss_bc% plat (  iperm)
    gnss_bc% data = gnss_bc% data (:,iperm)
    deallocate (iperm)
    deallocate (sortarr)

    !------------------------------------------
    ! Write updated correction coefficient file
    !------------------------------------------
    call write_gnss_bc_file (gnss_bc)

    !---------
    ! printout
    !---------
    if (dace% lpio) then
      write(6,'(a)') repeat('-',79)
      write(6,'( )')
      write(6,'(a)')      '  GNSS GPSGB STD/ZTD bias correction after analysis'
      write(6,'( )')
      write(6,'(a,a)')    '    last_date = ',gnss_bc% h% last_date
      write(6,'(a,a)')    '    ana_date  = ',cyyyymmddhhmm(ana_time)
      write(6,'(a,f11.2)')'    t_decay   =',    gnss_bc% h%    t_decay(1)
      !write(6,'(a,f13.4)')'    dt        =',                       dt
      !write(6,'(a,f13.4)')'    f         =',                       f
      write(6,'(a,i8)')   '    nb        =',sum(gnss_bc% data% nb)
      write(6,'(a,i8)')   '    na        =',sum(gnss_bc% data% na)
      write(6,'(a,f11.2)')'    nb_ac     =',sum(gnss_bc% data% nb_ac)
      write(6,'(a,f11.2)')'    na_ac     =',sum(gnss_bc% data% na_ac)
      write(6,'( )')
      if (verbose > 0) then
        !---------------------------------
        ! Verbosity of bias correction:
        ! 1 = print only non-empty entries
        ! 2 = print all
        !---------------------------------
        write(6,'( )')
      end if
    endif
    call destruct_bc_file (gnss_bc)
  end subroutine gnss_bc_aan

!==============================================================================

  subroutine gnss_bc_afg (obs, lwrite)
  !-----------------------------------------------------------------------
  ! Bias correction routine to be called after analysis
  ! Update bias correction statistics.
  ! Write updated correction coefficient file.
  !-----------------------------------------------------------------------
  type (t_obs_set) ,intent(inout) :: obs    ! observation
  logical          ,intent(in)    :: lwrite ! write bc file

    integer             :: i, j      ! index
!   integer             :: n         ! number of gnsss in statistics
    integer             :: ib        ! box index
    integer             :: is        ! spot index
    integer             :: ibc       ! plane index
    real(wp)            :: bg        ! background
    real(wp)            :: o         ! raw observation
    real(wp)            :: eo        ! observation error
    real(wp)            :: eb        ! background  error
    real(wp)            :: dob       ! obs - fg
    real(wp)            :: dt        ! time since last update of statistics
    real(wp)            :: f         ! weight for accumulated statistics
    integer             :: ier       ! error return flag
    integer,    pointer :: iperm (:) ! permutation index array for sorting
    type(t_bc), pointer :: pp  (:,:) ! pointer to statistics file entries
    type(t_bc), pointer :: p     (:) ! pointer to statistics file entries
    real(wp)            :: gmfh      ! GMF: hydrostatic mapping function
    real(wp)            :: gmfw      ! GMF: wet mapping function
    real(wp), parameter :: pi05 = 1.57079632679489655D0 ! 0.5*Pi = 90 deg.

    ! sorting by station/provider/product
    character (len=16), allocatable  :: sortarr(:) ! array to be sorted
    character (len=3)                :: prdstr     ! string with product ID
    character (len=3)                :: cntstr     ! string with center ID

    if (biascor_mode == 0) return
!   n = size(gnss_bc% plat)

    !-----------------
    ! set sums to zero
    !-----------------
    pp => gnss_bc% data(:,:)
    pp% nb    = 0
    pp% ob    = 0._wp  ! obs-fg mean deviation
    pp% ob_ob = 0._wp  ! obs-fg variance
    pp% o_err = 0._wp  !
    pp% b_err = 0._wp  !

    !-----------------------
    ! Update bias statistics
    !-----------------------
    do ib = 1, size (obs% o)
      if (obs% o(ib)% pe /= dace% pe) cycle
      do is = 1, obs% o(ib)% n_spot
        if (obs% o(ib)% spot(is)% hd% obstype /= OT_GPSGB) cycle
        ibc = obs% o(ib)% spot(is)% bc_index
        if (ibc <= 0)                                      cycle
        do i = obs% o(ib)% spot(is)% o% i + 1,                       &
               obs% o(ib)% spot(is)% o% i + obs% o(ib)% spot(is)% o% n
          j = 1
          p => gnss_bc% data(:,ibc)
          select case (obs% o(ib)% body(i)% use% state)
          case (STAT_PASSIVE, STAT_ACTIVE_0I:STAT_ACCEPTED)
            bg = obs% o(ib)% body(i)% bg
            o  = obs% o(ib)% body(i)% o - obs% o(ib)% body(i)% bc
            eo = obs% o(ib)% body(i)% eo       ! observational error
            eb = obs% o(ib)% body(i)% eb       ! background error
            dob  = o  - bg
            if (obs% o(ib)% varno(i) == VN_SPD) then
               ! map obsrvations to zenith using the Global Mapping Function GMF
               call gmf( mjd(obs% o(ib)% spot(is)% hd % time),       & ! MJD
                         obs% o(ib)% spot(is)% col % c % dlat * d2r, & ! lat
                         obs% o(ib)% spot(is)% col % c % dlon * d2r, & ! lon
                         obs% o(ib)% spot(is)% z,                    & ! height
                         pi05 - obs% o(ib)% olev(i) *d2r,           & ! Z
                         gmfh, gmfw)

!!$                write(*,*) 'GlobalMappingFkt : ',             &
!!$                     obs% o(ib)% spot(is)% hd % time,         &
!!$                     mjd(obs% o(ib)% spot(is)% hd % time),    &
!!$                     obs% o(ib)% spot(is)% col % c % dlon ,   &
!!$                     obs% o(ib)% spot(is)% col % c % dlat,    &
!!$                     obs% o(ib)% spot(is)% z,                 &
!!$                     obs% o(ib)% olev(i),              &
!!$                     pi05 - obs% o(ib)% olev(i) *d2r,  &
!!$                     gmfh, gmfw

               eo  = eo / gmfh
               eb  = eb / gmfh
               dob = dob / gmfh
            end if

            p(j)% nb     = p(j)% nb     + 1
            p(j)% ob     = p(j)% ob     + dob
            p(j)% ob_ob  = p(j)% ob_ob  + dob * dob
            p(j)% o_err  = p(j)% o_err  + eo
            p(j)% b_err  = p(j)% b_err  + eb
          end select
        end do
      end do
    end do

    !----------------
    ! sum up over PEs
    !----------------
    gnss_bc% data% nb     = p_sum (gnss_bc% data% nb)
    gnss_bc% data% ob     = p_sum (gnss_bc% data% ob)
    gnss_bc% data% ob_ob  = p_sum (gnss_bc% data% ob_ob)
    gnss_bc% data% o_err  = p_sum (gnss_bc% data% o_err )
    gnss_bc% data% b_err  = p_sum (gnss_bc% data% b_err )

    !-----------------------------------------
    ! weight factor for accumulated statistics
    !-----------------------------------------
    dt = days (ana_time - time_cyyyymmddhhmm (gnss_bc% h% last_date))
    if (gnss_bc% h% t_decay(1) > 0._wp .and. dt > 0._wp) then
      f = exp ( - dt / gnss_bc% h% t_decay(1))
    else
      f = 1._wp
    endif

    !-----------------------------------
    ! Update accumulated bias statistics
    !-----------------------------------
    pp => gnss_bc% data(:,:)
    !-------------------------------
    ! rescale accumulated statistics
    !-------------------------------
    pp% nb_ac      = pp% nb_ac    * f
    pp% ob_ac      = pp% ob_ac    * pp% nb_ac
    pp% ob_ob_ac   = pp% ob_ob_ac * pp% nb_ac

    pp% na_ac      = pp% na_ac    * f

    pp% o_err_ac   = pp% o_err_ac * pp% nb_ac
    pp% b_err_ac   = pp% b_err_ac * pp% nb_ac

    !-----------
    ! accumulate
    !-----------
    pp% nb_ac      = pp% nb_ac    + pp% nb
    pp% ob_ac      = pp% ob_ac    + pp% ob
    pp% ob_ob_ac   = pp% ob_ob_ac + pp% ob_ob
    pp% o_err_ac   = pp% o_err_ac + pp% o_err
    pp% b_err_ac   = pp% b_err_ac + pp% b_err

    !--------
    ! rescale
    !--------
    where (pp% nb > 0)
      pp% ob       = pp% ob       / pp% nb
      pp% ob_ob    = pp% ob_ob    / pp% nb
      pp% o_err    = pp% o_err    / pp% nb
      pp% b_err    = pp% b_err    / pp% nb
    endwhere
    where (pp% nb_ac > 0)
      pp% ob_ac    = pp% ob_ac    / pp% nb_ac
      pp% ob_ob_ac = pp% ob_ob_ac / pp% nb_ac
      pp% o_err_ac = pp% o_err_ac / pp% nb_ac
      pp% b_err_ac = pp% b_err_ac / pp% nb_ac
    endwhere

    if (lwrite) then
    !-------------------------------
    ! sort by station/center/product
    !-------------------------------
    allocate (iperm  (gnss_bc% n))
    allocate (sortarr(gnss_bc% n))
    sortarr = gnss_bc% plat% statid
    do i=1, gnss_bc% n
       write(cntstr,'(i3.3)') gnss_bc% plat(i)% center
       write(prdstr,'(i3.3)') gnss_bc% plat(i)% product
       sortarr(i) = sortarr(i)(1:10) // cntstr // prdstr
    end do
    call sort (sortarr, iperm, 1, ier)
    gnss_bc% plat = gnss_bc% plat (  iperm)
    gnss_bc% data = gnss_bc% data (:,iperm)
    deallocate (iperm)
    deallocate (sortarr)

    !------------------------------------------
    ! Write updated correction coefficient file
    !------------------------------------------
    call write_gnss_bc_file (gnss_bc)

    !---------
    ! printout
    !---------
    if (dace% lpio) then
      write(6,'(a)') repeat('-',79)
      write(6,'( )')
      write(6,'(a)')      '  GNSS STD/ZTD bias correction after first guess'
      write(6,'( )')
      write(6,'(a,a)')    '    last_date = ',gnss_bc% h% last_date
      write(6,'(a,a)')    '    ana_date  = ',cyyyymmddhhmm(ana_time)
      write(6,'(a,f11.2)')'    t_decay   =',    gnss_bc% h%    t_decay(1)
      write(6,'(a,f13.4)')'    dt        =',                       dt
      write(6,'(a,f13.4)')'    f         =',                       f
      write(6,'(a,i8)')   '    nb        =',sum(gnss_bc% data% nb)
      write(6,'(a,i8)')   '    na        =',sum(gnss_bc% data% na)
      write(6,'(a,f11.2)')'    nb_ac     =',sum(gnss_bc% data% nb_ac)
      write(6,'(a,f11.2)')'    na_ac     =',sum(gnss_bc% data% na_ac)
      write(6,'( )')
      if (verbose > 0) then
        !---------------------------------
        ! Verbosity of bias correction:
        ! 1 = print only non-empty entries
        ! 2 = print all
        !---------------------------------
        write(6,'( )')
      end if
    endif
    call destruct_bc_file (gnss_bc)
    endif
  end subroutine gnss_bc_afg

!==============================================================================
!----------------------------------------------------------------------------
! subroutine p_gather_tmp (sendbuffer,receivebuffer,root,[comm],[recvcounts])
!----------------------------------------------------------------------------
#define DERIVED type(t_tmp)
#undef  MPI_TYPE
#define p_gather_DERIVED p_gather_tmp
#include "p_gather_derived.incf"
#undef  DERIVED
#undef  p_gather_DERIVED
!==============================================================================
!--------------------------------------------------------------------
! subroutine p_scatter_tmp (sendbuf,recvbuf,root,[comm],[sendcounts])
!--------------------------------------------------------------------
#define DERIVED type(t_tmp)
#undef  MPI_TYPE
#define p_scatter_DERIVED p_scatter_tmp
#include "p_scatter_derived.incf"
#undef  DERIVED
!==============================================================================
!----------------------------------------------------
! subroutine p_bcast_plat (buffer, p_source, [comm])
!----------------------------------------------------
#define DERIVED type(t_plat_bc),dimension(:)
#define VECTOR
#undef  MPI_TYPE
#define p_bcast_DERIVED p_bcast_plat
#include "p_bcast.incf"
#undef  DERIVED
#undef  VECTOR
#undef  p_bcast_DERIVED
!----------------------------------------------------
! subroutine p_bcast_data (buffer, p_source, [comm])
!----------------------------------------------------
#define DERIVED type(t_bc),dimension(:,:)
#define VECTOR
#define RANK 2
#undef  MPI_TYPE
#define p_bcast_DERIVED p_bcast_data
#include "p_bcast.incf"
#undef  DERIVED
#undef  VECTOR
#undef  RANK
#undef  p_bcast_DERIVED
!==============================================================================
end module mo_gnss_bc
