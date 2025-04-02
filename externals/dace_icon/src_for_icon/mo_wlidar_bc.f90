!
!+ Bias correction for individual conventional platforms
!
MODULE mo_wlidar_bc
!
! Description:
!   Bias correction for Aeolus Wind Lidar HLOS winds
!
! Current Code Owner: DWD, Alexander Cress
!    phone: +49 69 8062 2716
!    fax:   +49 69 8062 3721
!    email: alexander.cress@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_49        2019-01-21 Alexander Cress
!  Online bias correction for WLIDAR
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
!-----------------------------------------------------------

!=============
! Modules used
!=============

  use mo_kind,       only: wp                ! working precision kind parameter
  use mo_exception,  only: finish            ! abort in case of error
  use mo_namelist,   only: position_nml,    &! position namelist
                           nnml,            &! namelist Fortran unit number
                           POSITIONED        ! ok    code from position_nml
  use mo_time,       only: days,            &! derive days (real)
                           operator(-),     &! calculate time difference
                           cyyyymmddhhmm,   &! string from time
                           time_cyyyymmddhhmm! time from string
  use mo_t_obs,      only: t_spot            ! observation data type
  use mo_obs_set,    only: t_obs_set         ! observation data derived type
  use mo_obs_tables, only: rept_use          ! use table entry
  use mo_dec_matrix, only: t_vector          ! vector data type
  use mo_mpi_dace,   only: dace,            &! MPI group info
                           p_sum,           &! sum over PEs
                           p_bcast,         &! generic MPI broadcast routine
                           p_gather          ! generic MPI gather routine
  use mo_run_params, only: ana_time,        &! analysis time
!                          run_type,        &! haupt=0, vor=1, ass=2
                           flag_biasc_wlidar ! wind lidar bias correction
  use mo_biasc_io,   only: t_bcor_head,     &! file header derived type
                           new_bc_head,     &! construct new file header
                           open_bc_read,    &! open  file for reading
                           open_bc_write,   &! open  file for writing
                           close_bc,        &! close file
                           bc_paths,        &! full pathnames
                           bc_obstyp,       &! observation type in files
                           nbcf,            &! number of files in list
                           verbose           ! Verbosity level of bias corr.
  use mo_fdbk_tables,only: VN_HLOS,         &! hlos observation flag
                           OT_WLIDAR         ! wind lidar report type ID
  use mo_t_use,      only: CHK_BIASCOR,     &! flag for no bias correction
                           decr_use,        &! decrease the state of a datum
                           STAT_PASSIVE,    &! observation status flags
                           STAT_ACTIVE_0I,  &!
                           STAT_ACCEPTED     !
  use mo_satid,      only: satid_t,         &! satellite identifier table
                           satname,         &! derive satellite name from satid
                           satid             ! derive satellite id from name
  use netcdf,        only:                  &! NetCDF f90 interface
                           nf90_def_dim,    &! define dimensions
                           nf90_def_var,    &! define variables
                           nf90_put_att,    &! define attribute
                           nf90_enddef,     &! end of definition mode
                           nf90_put_var,    &! write variable
                           nf90_strerror,   &! derive error character string
                           NF90_FLOAT,      &! float     type id
                           NF90_CHAR,       &! character type id
                           NF90_INT,        &! character type id
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
  implicit none

!================

!================
! Public entities
!================

  private
  !--------------
  ! derived types
  !--------------
  public :: t_wlidar_bc       ! hlos bias correction data
  public :: t_bc              ! components of t_wlidar_bc
  !------------
  ! subroutines
  !------------
  public :: wlidar_bc_init     ! initialize module: read namelist, biascor.coeff.
  public :: wlidar_bc_bfg      ! apply bias correction, called before fg-check
  public :: wlidar_bc_aan      ! update bias correction coefs., after analysis
  public :: read_bcor_file     ! read wind lidar bias correction file
  public :: read_wlidar_bc_nml ! read wind lidar bias correction namelist
  !--------------------------------------------------------
  ! namelist parameters to be set in namelist BIASCOR_SCATT
  !--------------------------------------------------------
  public :: biascor_mode      ! mode used for updating
  public :: t_decay           ! accumulation decay time (days)
  public :: n_required        ! number of entries required for correction
  public :: bc_fallback       ! fallback if biasc-file not present
  public :: BC_NOBC,BC_UP,BC_FG,BC_AN,BC_VARBC ! values for biascor_mode
  !-----------------------------------------------------
  ! Required in other modules/routines (e.g. monitoring)
  !-----------------------------------------------------
  public :: PH_MIE, PH_RAY
  public :: PH_ASC, PH_DES
  public :: bc_levs
  !--------------------------------
  ! Obsolete static bias correction
  !--------------------------------
  ! public :: wlidar_bc_init      ! Initialize wind lidar bias correction
  ! public :: wlidar_bc_bfg       ! Apply bias correction, called before FG check
  
  !=========================
  ! Derived type definitions
  !=========================
  
  !----------------------------
  ! wind lidar specific data
  !----------------------------
  type t_bc
    integer  :: n         = 0      ! number of entries
    real(wp) :: o_s       = 0._wp  ! obs  mean
    real(wp) :: b_s       = 0._wp  ! bg   mean
    real(wp) :: a_s       = 0._wp  ! an   mean
    real(wp) :: o_s2      = 0._wp  ! obs  square mean
    real(wp) :: b_s2      = 0._wp  ! bg   square mean
    real(wp) :: a_s2      = 0._wp  ! an   square mean
    real(wp) :: ob_s      = 0._wp  ! bg*obs mean
    real(wp) :: oa_s      = 0._wp  ! an*obs mean
    real(wp) :: ob        = 0._wp  ! obs-bg  mean deviation
    real(wp) :: oa        = 0._wp  ! obs-ana mean deviation
    real(wp) :: ab        = 0._wp  ! ana-obs mean deviation
    real(wp) :: ob_ob     = 0._wp  ! obs-bg  variance
    real(wp) :: oa_oa     = 0._wp  ! obs-ana variance
    real(wp) :: ab_ab     = 0._wp  ! ana-obs variance
    real(wp) :: ob_oa     = 0._wp  !
    real(wp) :: ob_ab     = 0._wp  !
    real(wp) :: oa_ab     = 0._wp  !
    real(wp) :: n_ac      = 0._wp  ! number of accumulated entries
    real(wp) :: o_s_ac    = 0._wp  ! obs  mean accumulated
    real(wp) :: b_s_ac    = 0._wp  ! bg   mean accumulated
    real(wp) :: a_s_ac    = 0._wp  ! an   mean accumulated
    real(wp) :: o_s2_ac   = 0._wp  ! obs  square mean accumulated
    real(wp) :: b_s2_ac   = 0._wp  ! bg   square mean accumulated
    real(wp) :: a_s2_ac   = 0._wp  ! an   square mean accumulated
    real(wp) :: ob_s_ac   = 0._wp  ! bg*obs mean accumulated
    real(wp) :: oa_s_ac   = 0._wp  ! an*obs mean accumulated
    real(wp) :: ob_ac     = 0._wp  ! obs-bg  mean     accumulated
    real(wp) :: oa_ac     = 0._wp  ! ana-bg  mean     accumulated
    real(wp) :: ab_ac     = 0._wp  ! ana-obs mean     accumulated
    real(wp) :: ob_ob_ac  = 0._wp  ! obs-bg  variance accumulated
    real(wp) :: oa_oa_ac  = 0._wp  ! ana-bg  variance accumulated
    real(wp) :: ab_ab_ac  = 0._wp  ! ana-obs variance accumulated
    real(wp) :: ob_oa_ac  = 0._wp  !
    real(wp) :: ob_ab_ac  = 0._wp  !
    real(wp) :: oa_ab_ac  = 0._wp  !
    real(wp) :: o_err     = 0._wp  ! nominal obs error
    real(wp) :: b_err     = 0._wp  ! nominal bg  error
    real(wp) :: o_err_ac  = 0._wp  ! nominal obs error
    real(wp) :: b_err_ac  = 0._wp  ! nominal bg  error
    real(wp) :: bc_b      = 0._wp  ! bias correction background value
    real(wp) :: bc_a      = 0._wp  ! bias correction analysis   value
  end type t_bc

  !---------------------------------
  ! Aeolus reciever channel indices
  !---------------------------------
  integer, parameter :: nph      = 2  ! number of aeolus receiver channels
  integer, parameter :: ntr      = 2  ! number of aeolus track directions (ascent or descent)
  integer, parameter :: nlv      = 5  ! number of aeolus bc levels
  integer ,parameter :: PH_MIE   = 1  ! Mie channel
  integer ,parameter :: PH_RAY   = 2  ! Rayleigh channel
  integer ,parameter :: PH_ASC   = 1  ! Ascending track
  integer ,parameter :: PH_DES   = 2  ! Descending track

  !-----------------------
  ! platform specific data
  !-----------------------
  type t_plat_bc
    integer          :: satid   = 0    ! wind lidar satellite id
    character(len=8) :: satids  =''
  end type t_plat_bc

  !-----------------------------------
  ! Aeolus Wind Lidar bias correction d ata
  !-----------------------------------
  type t_wlidar_bc
    type(t_bcor_head)        :: h                 ! file header data
    integer                  :: biascor_mode = 0  ! mode used for updating
    integer                  :: n            = 0  ! number of entries
    type (t_plat_bc),pointer :: plat (:)          ! platform specific metadata
    type (t_bc)     ,pointer :: data (:,:,:,:)    ! data
  end type t_wlidar_bc

  !-------------------------------------
  ! derived type used for cross-checking
  !-------------------------------------
  type t_tmp
    integer          :: satid =  0    ! wind_lidar satellite id
!   integer          :: pe    = -1    ! processor       of report
    integer          :: i     = -1    ! old index
    integer          :: j     = -1    ! new index
  end type t_tmp

  !------------------------
  ! values for biascor_mode
  !------------------------
  integer ,parameter :: BC_NOBC  = 0  ! no bias correction
  integer ,parameter :: BC_UP    = 1  ! only update bias corrrection file
  integer ,parameter :: BC_FG    = 2  ! apply  bias corr. from first guess
  integer ,parameter :: BC_AN    = 3  ! apply  bias corr. from analysis
  integer ,parameter :: BC_VARBC = 4  ! variational bias correction

  integer ,parameter :: BC_LEV   = 6  ! Number of BC Levels
  real(wp),parameter :: bc_levs(BC_LEV) = (/ 1050, 850, 500, 200, 70, 5 /) * 100._wp   ! hPa -> Pa

!------------------------------------------------------------------------------
  !----------------------------------------------
  ! type t_hlos_bc: wind liddar bias correction
  !----------------------------------------------
  type t_hlos_bc
     integer  :: satid =  0               ! Satellite id
     integer  :: instr = -1               ! Instrument id (currently unused)
     real(wp) :: scale =  1._wp           ! Scale factor of observation
     real(wp) :: const =  0._wp           ! Constant offset
  end type t_hlos_bc
  integer               :: nhlos = 0      ! Actual number of bias corrections
  integer,    parameter :: mhlos = 10     ! Max. number of bias corrections
  type(t_hlos_bc), save :: hlos_bc(mhlos)
!------------------------------------------------------------------------------

!=================
! Module variables
!=================

  type(t_wlidar_bc),save :: wlidar_bc              ! bias correction data
  !-----------------
  ! namelist entries
  !-----------------
  integer                  :: biascor_mode = BC_NOBC   ! mode used for updating
  real(wp)                 :: t_decay      =   -30._wp ! accumulation decaytime
  integer                  :: n_required   =  5000     ! # of entries required
  logical                  :: bc_fallback  = .false.   ! biasc-file not present

contains
!==============================================================================

  subroutine wlidar_bc_init
  !---------------------------------------
  ! Initialize this module:
  ! read bias-correction coefficient files
  !---------------------------------------

    integer :: i    ! loop index
    integer :: ierr ! error return value

    if (biascor_mode /= 0) then
      !---------------------------------------
      ! read bias-correction coefficient files
      !---------------------------------------
      ierr = -1
      if (flag_biasc_wlidar /= 0 .or. .not. bc_fallback) then
        do i = 1, nbcf
          if (bc_obstyp(i) == OT_WLIDAR) then
            call read_bcor_file (wlidar_bc, bc_paths(i))
            ierr = 0
            exit
          endif
        end do
      end if
      !-----------------------------------------------
      ! create empty bias-correction coefficient files
      !-----------------------------------------------
      if (ierr /= 0) then
!print*,'Bin in Create bias file: ', bc_fallback, biascor_mode
        if (.not.bc_fallback) &
          call finish ('wlidar_bc_init','coefficient file not present !')
        call new_bc_head (wlidar_bc% h, OT_WLIDAR, t_decay=[abs(t_decay)])
        wlidar_bc% biascor_mode = biascor_mode
        wlidar_bc% n            = 0
        allocate (wlidar_bc% plat (0))
        allocate (wlidar_bc% data (nph,ntr,nlv,0))
        if (dace% lpio) write(6,'(a,a,i2)')'    created empty file ', &
                        trim (wlidar_bc%h% path), wlidar_bc% biascor_mode
      endif
      !-------------------------------
      ! prepare bias correction to use
      !-------------------------------
      select case (biascor_mode)
      case (BC_UP)
        wlidar_bc% data% bc_b = 0
        wlidar_bc% data% bc_a = 0
      case (BC_FG)
        where (wlidar_bc% data% n_ac >= n_required)
          wlidar_bc% data% bc_b = wlidar_bc% data% ob_ac
          wlidar_bc% data% bc_a = wlidar_bc% data% ob_ac
        elsewhere
          wlidar_bc% data% bc_b = NF90_FILL_FLOAT
          wlidar_bc% data% bc_a = NF90_FILL_FLOAT
        endwhere
      case (BC_AN)
        where (wlidar_bc% data% n_ac >= n_required)
          wlidar_bc% data% bc_b = wlidar_bc% data% oa_ac
          wlidar_bc% data% bc_a = wlidar_bc% data% oa_ac
        elsewhere
          wlidar_bc% data% bc_b = NF90_FILL_FLOAT
          wlidar_bc% data% bc_a = NF90_FILL_FLOAT
        endwhere
      case (BC_VARBC)
        wlidar_bc% data% bc_b = wlidar_bc% data% bc_a
      end select
    endif
!   write(6,'(a,i2,i6,i6,i6)')'    wlidar_bc_init:', &
!        biascor_mode, size(wlidar_bc% data% n_ac), &
!        n_required,size(wlidar_bc% data% bc_b)

  end subroutine wlidar_bc_init

  !------------------------------------------------------------------------------

  subroutine read_bcor_file (bc, file, inst)
  !----------------------------------------
  ! read wind lidar bias correction file
  !----------------------------------------
  use mo_exception,   only: finish           ! abort routine

  use mo_head_netcdf, only: ncid             ! NetCDF file id

  !---------------------
  ! netCDF f90 interface
  !---------------------
  use netcdf,       only: nf90_Inquire,          &!
                          nf90_Inquire_Dimension,&!
                          nf90_Inquire_Variable, &!
                          nf90_inq_varid,        &!
                          nf90_get_var,          &!
                          nf90_strerror,         &!
                          NF90_MAX_NAME,         &!
                          NF90_NOERR              !

  type(t_wlidar_bc),intent(out) :: bc    ! derived type variable to fill
  character(len=*) ,intent(in)  :: file  ! name of file to read
  logical,optional ,intent(in)  :: inst  ! read instantaneous data as well

    integer  :: i, isat, j, k, l
    logical  :: linst
    real(wp) :: s_o, s_u, c_o, scale, intercept

  !---------------------------------------
  ! NetCDF variable IDs for body   section
  !---------------------------------------
  integer              :: varid_satid      ! Satelite id
  integer              :: varid_n_ac       ! Number of data (accumulated)
  integer              :: varid_o_s_ac     ! Number of data
  integer              :: varid_b_s_ac     ! Number of data
  integer              :: varid_a_s_ac     ! Number of data
  integer              :: varid_o_s2_ac    ! Number of data
  integer              :: varid_b_s2_ac    ! Number of data
  integer              :: varid_a_s2_ac    ! Number of data
  integer              :: varid_ob_s_ac    ! Number of data
  integer              :: varid_oa_s_ac    ! Number of data
  integer              :: varid_ob_ac    ! Number of data
  integer              :: varid_oa_ac    ! Number of data
  integer              :: varid_ab_ac    ! Number of data
  integer              :: varid_ob_ob_ac ! Number of data
  integer              :: varid_oa_oa_ac ! Number of data
  integer              :: varid_ob_oa_ac ! Number of data
  integer              :: varid_ab_ab_ac ! Number of data
  integer              :: varid_ob_ab_ac ! Number of data
  integer              :: varid_oa_ab_ac ! Number of data
  integer              :: varid_o_err_ac ! Number of data
  integer              :: varid_b_err_ac ! Number of data

  integer              :: varid_n       ! Number of data
  integer              :: varid_o_s     ! Number of data
  integer              :: varid_b_s     ! Number of data
  integer              :: varid_a_s     ! Number of data
  integer              :: varid_o_s2    ! Number of data
  integer              :: varid_b_s2    ! Number of data
  integer              :: varid_a_s2    ! Number of data
  integer              :: varid_ob_s    ! Number of data
  integer              :: varid_oa_s    ! Number of data
  integer              :: varid_ob    ! Number of data
  integer              :: varid_oa    ! Number of data
  integer              :: varid_ab    ! Number of data
  integer              :: varid_ob_ob ! Number of data
  integer              :: varid_oa_oa ! Number of data
  integer              :: varid_ob_oa ! Number of data
  integer              :: varid_ab_ab ! Number of data
  integer              :: varid_oa_ab ! Number of data
  integer              :: varid_ob_ab ! Number of data
  integer              :: varid_bc_b  ! Number of data
  integer              :: varid_bc_a  ! Number of data
  integer              :: varid_o_err ! Number of data
  integer              :: varid_b_err ! Number of data

  integer              :: status         ! NetCDF status variable


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
      allocate (bc% plat (bc% n))
      allocate (bc% data (nph,ntr,nlv,bc% n))

      !---------------
      ! read variables
      !---------------
      if (bc% n > 0) then
        stanc = 1
        counc = 0
        strnc = 1

        status = nf90_inq_varid (ncid, 'satid' , varid_satid )
        if (status == nf90_noerr) then
           status = nf90_get_var (ncid, varid_satid, bc% plat% satid, start=(/1/), count=(/bc% n/) )
           if (status /= nf90_noerr ) then
              call finish('get_var_satid cannot read  ')
           endif
        else
            call finish('get_var_satid error in status')
        endif

!       call get_var (bc% plat% satid    ,'satid'   )

        status = nf90_inq_varid (ncid, 'n_ac' , varid_n_ac )
        if (status == nf90_noerr) then
           status = nf90_get_var (ncid, varid_n_ac, bc% data% n_ac, start=(/1,1,1,1/), count=(/nph,ntr,nlv,bc% n/) )
           if (status /= nf90_noerr ) then
              call finish('get_var_n_ac cannot read  ')
           endif
        else
            call finish('get_var_n_ac error in status')
        endif

!       call get_var (bc% data% n_ac     ,'n_ac'    )

        status = nf90_inq_varid (ncid, 'o_s_ac' , varid_o_s_ac )
        if (status == nf90_noerr) then
           status = nf90_get_var (ncid, varid_o_s_ac, bc% data% o_s_ac, start=(/1,1,1,1/), count=(/nph,ntr,nlv,bc% n/) )
           if (status /= nf90_noerr ) then
              call finish('get_var_o_s_ac cannot read  ')
           endif
        else
            call finish('get_var_o_s_ac error in status')
        endif

!       call get_var (bc% data% o_s_ac   ,'o_s_ac'  )

        status = nf90_inq_varid (ncid, 'b_s_ac' , varid_b_s_ac )
        if (status == nf90_noerr) then
           status = nf90_get_var (ncid, varid_b_s_ac, bc% data% b_s_ac, start=(/1,1,1,1/), count=(/nph,ntr,nlv,bc% n/) )
           if (status /= nf90_noerr ) then
              call finish('get_var_b_s_ac cannot read  ')
           endif
        else
            call finish('get_var_b_s_ac error in status')
        endif

!       call get_var (bc% data% b_s_ac   ,'b_s_ac'  )

        status = nf90_inq_varid (ncid, 'a_s_ac' , varid_a_s_ac )
        if (status == nf90_noerr) then
           status = nf90_get_var (ncid, varid_a_s_ac, bc% data% a_s_ac, start=(/1,1,1,1/), count=(/nph,ntr,nlv,bc% n/) )
           if (status /= nf90_noerr ) then
              call finish('get_var_a_s_ac cannot read  ')
           endif
        else
            call finish('get_var_a_s_ac error in status')
        endif

!       call get_var (bc% data% a_s_ac   ,'a_s_ac'  )

        status = nf90_inq_varid (ncid, 'o_s2_ac' , varid_o_s2_ac )
        if (status == nf90_noerr) then
           status = nf90_get_var (ncid, varid_o_s2_ac, bc% data% o_s2_ac, start=(/1,1,1,1/), count=(/nph,ntr,nlv,bc% n/) )
           if (status /= nf90_noerr ) then
             call finish('get_var_o_s2_ac cannot read  ')
           endif
        else
            call finish('get_var_o_s2_ac error in status')
        endif

!       call get_var (bc% data% o_s2_ac  ,'o_s2_ac' )

        status = nf90_inq_varid (ncid, 'b_s2_ac' , varid_b_s2_ac )
        if (status == nf90_noerr) then
           status = nf90_get_var (ncid, varid_b_s2_ac, bc% data% b_s2_ac, start=(/1,1,1,1/), count=(/nph,ntr,nlv,bc% n/) )
           if (status /= nf90_noerr ) then
              call finish('get_var_b_s2_ac cannot read  ')
           endif
        else
            call finish('get_var_b_s2_ac error in status')
        endif

!       call get_var (bc% data% b_s2_ac  ,'b_s2_ac' )

        status = nf90_inq_varid (ncid, 'a_s2_ac' , varid_a_s2_ac )
        if (status == nf90_noerr) then
           status = nf90_get_var (ncid, varid_a_s2_ac, bc% data% a_s2_ac, start=(/1,1,1,1/), count=(/nph,ntr,nlv,bc% n/) )
           if (status /= nf90_noerr ) then
              call finish('get_var_a_s2_ac cannot read  ')
           endif
        else
            call finish('get_var_a_s2_ac error in status')
        endif

!       call get_var (bc% data% a_s2_ac  ,'a_s2_ac' )

        status = nf90_inq_varid (ncid, 'ob_s_ac' , varid_ob_s_ac )
        if (status == nf90_noerr) then
           status = nf90_get_var (ncid, varid_ob_s_ac, bc% data% ob_s_ac, start=(/1,1,1,1/), count=(/nph,ntr,nlv,bc% n/) )
           if (status /= nf90_noerr ) then
              call finish('get_var_ob_s_ac cannot read  ')
           endif
        else
            call finish('get_var_ob_s_ac error in status')
        endif

!       call get_var (bc% data% ob_s_ac  ,'ob_s_ac' )

        status = nf90_inq_varid (ncid, 'oa_s_ac' , varid_oa_s_ac )
        if (status == nf90_noerr) then
           status = nf90_get_var (ncid, varid_oa_s_ac, bc% data% oa_s_ac, start=(/1,1,1,1/), count=(/nph,ntr,nlv,bc% n/) )
           if (status /= nf90_noerr ) then
              call finish('get_var_oa_s_ac cannot read  ')
           endif
        else
            call finish('get_var_oa_s_ac error in status')
        endif

!       call get_var (bc% data% oa_s_ac  ,'oa_s_ac' )

        status = nf90_inq_varid (ncid, 'ob_ac' , varid_ob_ac )
        if (status == nf90_noerr) then
           status = nf90_get_var (ncid, varid_ob_ac, bc% data% ob_ac, start=(/1,1,1,1/), count=(/nph,ntr,nlv,bc% n/) )
           if (status /= nf90_noerr ) then
              call finish('get_var_ob_ac cannot read  ')
           endif
        else
            call finish('get_var_ob_ac','error in status')
        endif

!       call get_var (bc% data% ob_ac    ,'ob_ac'   )

        status = nf90_inq_varid (ncid, 'oa_ac' , varid_oa_ac )
        if (status == nf90_noerr) then
           status = nf90_get_var (ncid, varid_oa_ac, bc% data% oa_ac, start=(/1,1,1,1/), count=(/nph,ntr,nlv,bc% n/) )
           if (status /= nf90_noerr ) then
              call finish('get_var_oa_ac cannot read  ')
           endif
        else
            call finish('get_var_oa_ac error in status')
        endif

!       call get_var (bc% data% oa_ac    ,'oa_ac'   )

        status = nf90_inq_varid (ncid, 'ab_ac' , varid_ab_ac )
        if (status == nf90_noerr) then
           status = nf90_get_var (ncid, varid_ab_ac, bc% data% ab_ac, start=(/1,1,1,1/), count=(/nph,ntr,nlv,bc% n/) )
           if (status /= nf90_noerr ) then
              call finish('get_var_ab_ac cannot read  ')
           endif
        else
            call finish('get_var_ab_ac error in status')
        endif

!       call get_var (bc% data% ab_ac    ,'ab_ac'   )

        status = nf90_inq_varid (ncid, 'ob_ob_ac' , varid_ob_ob_ac )
        if (status == nf90_noerr) then
           status = nf90_get_var (ncid, varid_ob_ob_ac, bc% data% ob_ob_ac, start=(/1,1,1,1/), count=(/nph,ntr,nlv,bc% n/) )
           if (status /= nf90_noerr ) then
              call finish('get_var_ob_ob_ac cannot read  ')
           endif
        else
            call finish('get_var_ob_ob_ac error in status')
        endif

!       call get_var (bc% data% ob_ob_ac ,'ob_ob_ac')

        status = nf90_inq_varid (ncid, 'oa_oa_ac' , varid_oa_oa_ac )
        if (status == nf90_noerr) then
           status = nf90_get_var (ncid, varid_oa_oa_ac, bc% data% oa_oa_ac, start=(/1,1,1,1/), count=(/nph,ntr,nlv,bc% n/) )
           if (status /= nf90_noerr ) then
              call finish('get_var_oa_oa_ac cannot read  ')
           endif
        else
            call finish('get_var_oa_oa_ac error in status')
        endif

!       call get_var (bc% data% oa_oa_ac ,'oa_oa_ac')

        status = nf90_inq_varid (ncid, 'ab_ab_ac' , varid_ab_ab_ac )
        if (status == nf90_noerr) then
           status = nf90_get_var (ncid, varid_ab_ab_ac, bc% data% ab_ab_ac, start=(/1,1,1,1/), count=(/nph,ntr,nlv,bc% n/) )
           if (status /= nf90_noerr ) then
              call finish('get_var_ab_ab_ac cannot read  ')
           endif
        else
            call finish('get_var_ab_ab_ac error in status')
        endif

!       call get_var (bc% data% ab_ab_ac ,'ab_ab_ac')

        status = nf90_inq_varid (ncid, 'ob_oa_ac' , varid_ob_oa_ac )
        if (status == nf90_noerr) then
           status = nf90_get_var (ncid, varid_ob_oa_ac, bc% data% ob_oa_ac, start=(/1,1,1,1/), count=(/nph,ntr,nlv,bc% n/) )
           if (status /= nf90_noerr ) then
              call finish('get_var_ob_oa_ac cannot read  ')
           endif
        else
            call finish('get_var_ob_oa_ac error in status')
        endif

!       call get_var (bc% data% ob_oa_ac ,'ob_oa_ac')

        status = nf90_inq_varid (ncid, 'ob_ab_ac' , varid_ob_ab_ac )
        if (status == nf90_noerr) then
           status = nf90_get_var (ncid, varid_ob_ab_ac, bc% data% ob_ab_ac, start=(/1,1,1,1/), count=(/nph,ntr,nlv,bc% n/) )
           if (status /= nf90_noerr ) then
              call finish('get_var_ob_ab_ac cannot read  ')
           endif
        else
            call finish('get_var_ob_ab_ac error in status')
        endif

!       call get_var (bc% data% ob_ab_ac ,'ob_ab_ac')

        status = nf90_inq_varid (ncid, 'oa_ab_ac' , varid_oa_ab_ac )
        if (status == nf90_noerr) then
           status = nf90_get_var (ncid, varid_oa_ab_ac, bc% data% oa_ab_ac, start=(/1,1,1,1/), count=(/nph,ntr,nlv,bc% n/) )
           if (status /= nf90_noerr ) then
              call finish('get_var_oa_ab_ac cannot read  ')
           endif
        else
            call finish('get_var_oa_ab_ac error in status')
        endif

!       call get_var (bc% data% oa_ab_ac ,'oa_ab_ac')

        status = nf90_inq_varid (ncid, 'o_err_ac' , varid_o_err_ac )
        if (status == nf90_noerr) then
           status = nf90_get_var (ncid, varid_o_err_ac, bc% data% o_err_ac, start=(/1,1,1,1/), count=(/nph,ntr,nlv,bc% n/) )
           if (status /= nf90_noerr ) then
              call finish('get_var_o_err_ac cannot read  ')
           endif
        else
            call finish('get_var_o_err_ac error in status')
        endif

!       call get_var (bc% data% o_err_ac ,'o_err_ac')

        status = nf90_inq_varid (ncid, 'b_err_ac' , varid_b_err_ac )
        if (status == nf90_noerr) then
           status = nf90_get_var (ncid, varid_b_err_ac, bc% data% b_err_ac, start=(/1,1,1,1/), count=(/nph,ntr,nlv,bc% n/) )
           if (status /= nf90_noerr ) then
              call finish('get_var_b_err_ac cannot read  ')
           endif
        else
            call finish('get_var_b_err_ac error in status')
        endif

!       call get_var (bc% data% b_err_ac ,'b_err_ac')

        status = nf90_inq_varid (ncid, 'bc_b' , varid_bc_b )
        if (status == nf90_noerr) then
           status = nf90_get_var (ncid, varid_bc_b, bc% data% bc_b, start=(/1,1,1,1/), count=(/nph,ntr,nlv,bc% n/) )
           if (status /= nf90_noerr ) then
              call finish('get_var_bc_b cannot read  ')
           endif
        else
            call finish('get_var_bc_b error in status')
        endif

!       call get_var (bc% data% bc_b     ,'bc_b'    )

        status = nf90_inq_varid (ncid, 'bc_a' , varid_bc_a )
        if (status == nf90_noerr) then
           status = nf90_get_var (ncid, varid_bc_a, bc% data% bc_a, start=(/1,1,1,1/), count=(/nph,ntr,nlv,bc% n/) )
           if (status /= nf90_noerr ) then
              call finish('get_var_bc_a cannot read  ')
           endif
        else
            call finish('get_var_bc_a error in status')
        endif

!       call get_var (bc% data% bc_a     ,'bc_a'    )

        if (linst) then
        status = nf90_inq_varid (ncid, 'n' , varid_n )
        if (status == nf90_noerr) then
           status = nf90_get_var (ncid, varid_n, bc% data% n, start=(/1,1,1,1/), count=(/nph,ntr,nlv,bc% n/) )
           if (status /= nf90_noerr ) then
              call finish('get_var_n cannot read  ')
           endif
        else
            call finish('get_var_n',' error in status')
        endif

!         call get_var (bc% data% n     ,'n'    )

        status = nf90_inq_varid (ncid, 'o_s' , varid_o_s )
        if (status == nf90_noerr) then
           status = nf90_get_var (ncid, varid_o_s, bc% data% o_s, start=(/1,1,1,1/), count=(/nph,ntr,nlv,bc% n/) )
           if (status /= nf90_noerr ) then
             call finish('get_var_o_s','cannot read')
           endif
        else
           call finish('get_var_o_s','error in status')
        endif

!         call get_var (bc% data% o_s   ,'o_s'  )

        status = nf90_inq_varid (ncid, 'b_s' , varid_b_s )
        if (status == nf90_noerr) then
           status = nf90_get_var (ncid, varid_b_s, bc% data% b_s, start=(/1,1,1,1/), count=(/nph,ntr,nlv,bc% n/) )
           if (status /= nf90_noerr ) then
             call finish('get_var_b_s','cannot read  ')
           endif
        else
           call finish('get_var_b_s','error in status')
        endif

!         call get_var (bc% data% b_s   ,'b_s'  )

        status = nf90_inq_varid (ncid, 'a_s' , varid_a_s )
        if (status == nf90_noerr) then
           status = nf90_get_var (ncid, varid_a_s, bc% data% a_s, start=(/1,1,1,1/), count=(/nph,ntr,nlv,bc% n/) )
           if (status /= nf90_noerr ) then
             call finish('get_var_a_s','cannot read  ')
           endif
        else
           call finish('get_var_a_s','error in status')
        endif

!         call get_var (bc% data% a_s   ,'a_s'  )

        status = nf90_inq_varid (ncid, 'o_s2' , varid_o_s2 )
        if (status == nf90_noerr) then
           status = nf90_get_var (ncid, varid_o_s2, bc% data% o_s2, start=(/1,1,1,1/), count=(/nph,ntr,nlv,bc% n/) )
           if (status /= nf90_noerr ) then
             call finish('get_var_o_s2','cannot read  ')
           endif
        else
           call finish('get_var_o_s2','error in status')
        endif

!         call get_var (bc% data% o_s2  ,'o_s2' )

        status = nf90_inq_varid (ncid, 'b_s2' , varid_b_s2 )
        if (status == nf90_noerr) then
           status = nf90_get_var (ncid, varid_b_s2, bc% data% b_s2, start=(/1,1,1,1/), count=(/nph,ntr,nlv,bc% n/) )
           if (status /= nf90_noerr ) then
              call finish('get_var_b_s2','cannot read  ')
           endif
        else
            call finish('get_var_b_s2','error in status')
        endif

!         call get_var (bc% data% b_s2  ,'b_s2' )

        status = nf90_inq_varid (ncid, 'a_s2' , varid_a_s2 )
        if (status == nf90_noerr) then
           status = nf90_get_var (ncid, varid_a_s2, bc% data% a_s2, start=(/1,1,1,1/), count=(/nph,ntr,nlv,bc% n/) )
           if (status /= nf90_noerr ) then
               call finish('get_var_a_s2','cannot read  ')
           endif
        else
            call finish('get_var_a_s2','error in status')
        endif

!         call get_var (bc% data% a_s2  ,'a_s2' )

        status = nf90_inq_varid (ncid, 'ob_s' , varid_ob_s )
        if (status == nf90_noerr) then
           status = nf90_get_var (ncid, varid_ob_s, bc% data% ob_s, start=(/1,1,1,1/), count=(/nph,ntr,nlv,bc% n/) )
           if (status /= nf90_noerr ) then
             call finish('get_var_ob_s','cannot read  ')
           endif
        else
            call finish('get_var_ob_s','error in status')
        endif

!         call get_var (bc% data% ob_s  ,'ob_s' )

        status = nf90_inq_varid (ncid, 'oa_s' , varid_oa_s )
        if (status == nf90_noerr) then
           status = nf90_get_var (ncid, varid_oa_s, bc% data% oa_s, start=(/1,1,1,1/), count=(/nph,ntr,nlv,bc% n/) )
           if (status /= nf90_noerr ) then
             call finish('get_var_oa_s','cannot read  ')
           endif
        else
            call finish('get_var_oa_s','error in status')
        endif

!         call get_var (bc% data% oa_s  ,'oa_s' )

        status = nf90_inq_varid (ncid, 'ob' , varid_ob )
        if (status == nf90_noerr) then
           status = nf90_get_var (ncid, varid_ob, bc% data% ob, start=(/1,1,1,1/), count=(/nph,ntr,nlv,bc% n/) )
           if (status /= nf90_noerr ) then
               call finish('get_var_ob','cannot read  ')
           endif
        else
            call finish('get_var_ob','error in status')
        endif

!         call get_var (bc% data% ob    ,'ob'   )

        status = nf90_inq_varid (ncid, 'oa' , varid_oa )
        if (status == nf90_noerr) then
           status = nf90_get_var (ncid, varid_oa, bc% data% oa, start=(/1,1,1,1/), count=(/nph,ntr,nlv,bc% n/) )
           if (status /= nf90_noerr ) then
              call finish('get_var_oa','cannot read  ')
           endif
        else
            call finish('get_var_oa','error in status')
        endif

!         call get_var (bc% data% oa    ,'oa'   )

        status = nf90_inq_varid (ncid, 'ab' , varid_ab )
        if (status == nf90_noerr) then
           status = nf90_get_var (ncid, varid_ab, bc% data% ab, start=(/1,1,1,1/), count=(/nph,ntr,nlv,bc% n/) )
           if (status /= nf90_noerr ) then
              call finish('get_var_ab','cannot read  ')
           endif
        else
            call finish('get_var_ab','error in status')
        endif

!         call get_var (bc% data% ab    ,'ab'   )

        status = nf90_inq_varid (ncid, 'ob_ob' , varid_ob_ob )
        if (status == nf90_noerr) then
           status = nf90_get_var (ncid, varid_ob_ob, bc% data% ob_ob, start=(/1,1,1,1/), count=(/nph,ntr,nlv,bc% n/) )
           if (status /= nf90_noerr ) then
              call finish('get_var_ob_ob','cannot read  ')
           endif
        else
            call finish('get_var_ob_ob','error in status')
        endif

!         call get_var (bc% data% ob_ob ,'ob_ob')

        status = nf90_inq_varid (ncid, 'oa_oa' , varid_oa_oa )
        if (status == nf90_noerr) then
           status = nf90_get_var (ncid, varid_oa_oa, bc% data% oa_oa, start=(/1,1,1,1/), count=(/nph,ntr,nlv,bc% n/) )
           if (status /= nf90_noerr ) then
              call finish('get_var_oa_oa','cannot read  ')
           endif
        else
            call finish('get_var_oa_oa','error in status')
        endif

!         call get_var (bc% data% oa_oa ,'oa_oa')

        status = nf90_inq_varid (ncid, 'ab_ab' , varid_ab_ab )
        if (status == nf90_noerr) then
           status = nf90_get_var (ncid, varid_ab_ab, bc% data% ab_ab, start=(/1,1,1,1/), count=(/nph,ntr,nlv,bc% n/) )
           if (status /= nf90_noerr ) then
              call finish('get_var_ab_ab','cannot read  ')
           endif
        else
            call finish('get_var_ab_ab','error in status')
        endif

!         call get_var (bc% data% ab_ab ,'ab_ab')

        status = nf90_inq_varid (ncid, 'ob_oa' , varid_ob_oa )
        if (status == nf90_noerr) then
           status = nf90_get_var (ncid, varid_ob_oa, bc% data% ob_oa, start=(/1,1,1,1/), count=(/nph,ntr,nlv,bc% n/) )
           if (status /= nf90_noerr ) then
              call finish('get_var_ob_oa','cannot read  ')
           endif
        else
            call finish('get_var_ob_oa','error in status')
        endif

!         call get_var (bc% data% ob_oa ,'ob_oa')

        status = nf90_inq_varid (ncid, 'ob_ab' , varid_ob_ab )
        if (status == nf90_noerr) then
           status = nf90_get_var (ncid, varid_ob_ab, bc% data% ob_ab, start=(/1,1,1,1/), count=(/nph,ntr,nlv,bc% n/) )
           if (status /= nf90_noerr ) then
              call finish('get_var_ob_ab','cannot read  ')
           endif
        else
            call finish('get_var_ob_ab','error in status')
        endif

!         call get_var (bc% data% ob_ab ,'ob_ab')

        status = nf90_inq_varid (ncid, 'oa_ab' , varid_oa_ab )
        if (status == nf90_noerr) then
           status = nf90_get_var (ncid, varid_oa_ab, bc% data% oa_ab, start=(/1,1,1,1/), count=(/nph,ntr,nlv,bc% n/) )
           if (status /= nf90_noerr ) then
              call finish('get_var_oa_ab','cannot read  ')
           endif
        else
            call finish('get_var_oa_ab','error in status')
        endif

!         call get_var (bc% data% oa_ab ,'oa_ab')

        status = nf90_inq_varid (ncid, 'o_err' , varid_o_err )
        if (status == nf90_noerr) then
           status = nf90_get_var (ncid, varid_o_err, bc% data% o_err, start=(/1,1,1,1/), count=(/nph,ntr,nlv,bc% n/) )
           if (status /= nf90_noerr ) then
              call finish('get_var_o_err','cannot read  ')
           endif
        else
            call finish('get_var_o_err','error in status')
        endif

!         call get_var (bc% data% o_err ,'o_err')

        status = nf90_inq_varid (ncid, 'b_err' , varid_b_err )
        if (status == nf90_noerr) then
           status = nf90_get_var (ncid, varid_b_err, bc% data% b_err, start=(/1,1,1,1/), count=(/nph,ntr,nlv,bc% n/) )
           if (status /= nf90_noerr ) then
              call finish('get_var_b_err','cannot read  ')
           endif
        else
            call finish('get_var_b_err','error in status')
        endif

!         call get_var (bc% data% b_err ,'b_err')
        endif
        do i = 1, bc% n
          do isat = 1, size(satid_t)
            if (satid_t(isat)% code == bc% plat(i)% satid) then
              bc% plat(i)% satids = satid_t(isat)% mnem
            endif
          enddo
        enddo
      endif ! linst
      !-----------
      ! close file
      !-----------
      call close_bc (bc% h)
    endif ! bc%n > 0

    !----------
    ! broadcast
    !----------
    call p_bcast (bc% n, dace% pio)
    if (.not.dace% lpio) then
      allocate (bc% plat (bc% n))
      allocate (bc% data (nph,ntr,nlv,bc% n))
    endif
    call p_bcast_plat (wlidar_bc% plat, dace% pio)
    call p_bcast_data (wlidar_bc% data, dace% pio)

    if (dace% lpio .and. wlidar_bc% n > 0) then
       write(6,'(a)') "  HLOS Wind Lidarbias correction"
       write(6,'()')
       write(6,'(a)') &
             '   satid      n_ac   scale*obs + offset'
sat:   do i = 1, wlidar_bc% n
       do k = 1, nlv
       do l = 1, ntr
       do j = 1, nph
          scale     = 1
          intercept = 0
          if (wlidar_bc% data(j,l,k,i)% n_ac > 0._wp) then
             s_o = (wlidar_bc% data(j,l,k,i)% ob_s_ac                             - &
                    wlidar_bc% data(j,l,k,i)% o_s_ac * wlidar_bc% data(j,l,k,i)% b_s_ac    )
             s_u = (wlidar_bc% data(j,l,k,i)% o_s2_ac                             - &
                    wlidar_bc% data(j,l,k,i)% o_s_ac**2                             )
             c_o = (wlidar_bc% data(j,l,k,i)% b_s_ac * wlidar_bc% data(j,l,k,i)% o_s2_ac - &
                    wlidar_bc% data(j,l,k,i)% o_s_ac * wlidar_bc% data(j,l,k,i)% ob_s_ac   )
             if (s_u > 0._wp) then
                scale     = s_o / s_u
                intercept = c_o / s_u
             end if
             write(6,'(i8,f11.0,2x,2f9.3)') &
                  wlidar_bc% plat(i)% satid, &
                  wlidar_bc% data(j,l,k,i)% n_ac,  &
                  scale, intercept
          end if
       end do
       end do
       end do
       end do sat
       write(6,'()')
    end if

  end subroutine read_bcor_file

!------------------------------------------------------------------------------

  subroutine destruct_bc_file (bc)
    type(t_wlidar_bc) ,intent(inout) :: bc    ! derived type variable to fill

    type(t_wlidar_bc) :: empty_bc             ! empty bias correction
    !----------------------
    ! deallocate components
    !----------------------
    deallocate (bc% plat)
    deallocate (bc% data)

    bc = empty_bc
  end subroutine destruct_bc_file

!------------------------------------------------------------------------------

  subroutine write_bcor_file (bc)
  !-----------------------------------------
  ! write HLOS Win Lidar bias correction file
  !-----------------------------------------
  type(t_wlidar_bc) ,intent(inout) :: bc

    integer :: status                            ! NetCDF return value
    integer :: dimid (4)                         ! NetCDF dimension id
    integer :: satid                             ! NetCDF variable ids
    integer :: n                                 ! .
    integer :: ob, oa, ab                        ! .
    integer :: o_s, b_s, a_s                     ! .
    integer :: o_s2, b_s2, a_s2                  ! .
    integer :: ob_s, oa_s                        ! .
    integer :: ob_ob, oa_oa, ab_ab               ! .
    integer :: ob_oa, ob_ab, oa_ab               ! .
    integer :: n_ac                              ! .
    integer :: o_s_ac, b_s_ac, a_s_ac            ! .
    integer :: o_s2_ac, b_s2_ac, a_s2_ac         ! .
    integer :: ob_s_ac, oa_s_ac                  ! .
    integer :: ob_ac, oa_ac, ab_ac               ! .
    integer :: ob_ob_ac, oa_oa_ac, ab_ab_ac      ! .
    integer :: ob_oa_ac, ob_ab_ac, oa_ab_ac      ! .
    integer :: bc_b, b_err , b_err_ac ! bc_b_err ! .
    integer :: bc_a, o_err , o_err_ac ! bc_a_err ! .

    if (dace% lpio .and. bc%n > 0) then
      !-------------------------------
      ! open NetCDF file, write header
      !-------------------------------
      call open_bc_write (bc% h, OT_WLIDAR)

      !---------------
      ! set attributes
      !---------------
      status = nf90_put_att (bc%h% ncid, NF90_GLOBAL, 'biascor_mode', &
                                                       biascor_mode   )
      !------------------
      ! define dimensions
      !------------------
!     call def_dim ('observation', nob,            dimid(1))
      call def_dim ('rtrtype',     nph           , dimid(1))
      call def_dim ('track',       ntr           , dimid(2))
      call def_dim ('bclevel',     nlv           , dimid(3))
      call def_dim ('platform',    size(bc% plat), dimid(4))

      !-----------------
      ! define variables
      !-----------------
      status = nf90_def_var (bc%h% ncid ,'satid' ,NF90_INT, &
                             dimid(4), satid)
      status = nf90_put_att (bc%h% ncid , satid, 'longname',&
                             'satellite id ')

      call def_var ('n'         ,n         ,'entries')
      call def_var ('o_s'       ,o_s       ,'obs mean')
      call def_var ('b_s'       ,b_s       ,'bg  mean')
      call def_var ('a_s'       ,a_s       ,'ana mean')
      call def_var ('o_s2'      ,o_s2      ,'obs square mean')
      call def_var ('b_s2'      ,b_s2      ,'bg  square mean')
      call def_var ('a_s2'      ,a_s2      ,'ana square mean')
      call def_var ('ob_s'      ,ob_s      ,'obs x bg mean')
      call def_var ('oa_s'      ,oa_s      ,'obs x ana mean')
      call def_var ('ob'        ,ob        ,'obs-bg  mean')
      call def_var ('oa'        ,oa        ,'obs-ana mean')
      call def_var ('ab'        ,ab        ,'ana-bg  mean')
      call def_var ('ob_ob'     ,ob_ob     ,'obs-bg  variance')
      call def_var ('oa_oa'     ,oa_oa     ,'obs-ana variance')
      call def_var ('ab_ab'     ,ab_ab     ,'ana-bg  variance')
      call def_var ('ob_oa'     ,ob_oa     ,'obs-bg  obs-ana covariance')
      call def_var ('ob_ab'     ,ob_ab     ,'obs-bg  ana-bg  covariance')
      call def_var ('oa_ab'     ,oa_ab     ,'obs-ana ana-bg  covariance')
      call def_var ('n_ac'      ,n_ac      ,'entries accumulated')
      call def_var ('o_s_ac'    ,o_s_ac    ,'obs mean accumulated')
      call def_var ('b_s_ac'    ,b_s_ac    ,'bg  mean accumulated')
      call def_var ('a_s_ac'    ,a_s_ac    ,'ana mean accumulated')
      call def_var ('o_s2_ac'   ,o_s2_ac   ,'obs square mean accumulated')
      call def_var ('b_s2_ac'   ,b_s2_ac   ,'bg  square mean accumulated')
      call def_var ('a_s2_ac'   ,a_s2_ac   ,'ana square mean accumulated')
      call def_var ('ob_s_ac'   ,ob_s_ac   ,'obs x bg mean accumulated')
      call def_var ('oa_s_ac'   ,oa_s_ac   ,'obs x ana mean accumulated')
      call def_var ('ob_ac'     ,ob_ac     ,'obs-bg  mean accumulated')
      call def_var ('oa_ac'     ,oa_ac     ,'obs-ana mean accumulated')
      call def_var ('ab_ac'     ,ab_ac     ,'ana-bg  mean accumulated')
      call def_var ('ob_ob_ac'  ,ob_ob_ac  ,'obs-bg  variance accumulated')
      call def_var ('oa_oa_ac'  ,oa_oa_ac  ,'obs-ana variance accumulated')
      call def_var ('ab_ab_ac'  ,ab_ab_ac  ,'ana-bg  variance accumulated')
      call def_var ('ob_oa_ac'  ,ob_oa_ac  ,'obs-bg  obs-ana covariance acc.')
      call def_var ('ob_ab_ac'  ,ob_ab_ac  ,'obs-bg  ana-bg  covariance acc.')
      call def_var ('oa_ab_ac'  ,oa_ab_ac  ,'obs-ana ana-bg  covariance acc.')
      call def_var ('bc_b'      ,bc_b      ,'first guess bias correction')
      call def_var ('bc_a'      ,bc_a      ,'analysis bias correction')
      call def_var ('o_err'     ,o_err     ,'nominal obs error squared')
      call def_var ('b_err'     ,b_err     ,'nominal bg  error squared')
      call def_var ('o_err_ac'  ,o_err_ac  ,               &
                    'nominal obs error squared accumulated')
      call def_var ('b_err_ac'  ,b_err_ac  ,               &
                    'nominal bg  error squared accumulated')
      status = nf90_enddef  (bc%h% ncid)
      !----------------
      ! write variables
      !----------------
      status = nf90_put_var (bc%h% ncid, satid, bc% plat% satid)
      call put_vari ('n'        ,n        ,bc% data% n)
      call put_var  ('o_s'      ,o_s      ,bc% data% o_s)
      call put_var  ('b_s'      ,b_s      ,bc% data% b_s)
      call put_var  ('a_s'      ,a_s      ,bc% data% a_s)
      call put_var  ('o_s2'     ,o_s2     ,bc% data% o_s2)
      call put_var  ('b_s2'     ,b_s2     ,bc% data% b_s2)
      call put_var  ('a_s2'     ,a_s2     ,bc% data% a_s2)
      call put_var  ('ob_s'     ,ob_s     ,bc% data% ob_s)
      call put_var  ('oa_s'     ,oa_s     ,bc% data% oa_s)
      call put_var  ('ob'       ,ob       ,bc% data% ob)
      call put_var  ('oa'       ,oa       ,bc% data% oa)
      call put_var  ('ab'       ,ab       ,bc% data% ab)
      call put_var  ('ob_ob'    ,ob_ob    ,bc% data% ob_ob)
      call put_var  ('oa_oa'    ,oa_oa    ,bc% data% oa_oa)
      call put_var  ('ab_ab'    ,ab_ab    ,bc% data% ab_ab)
      call put_var  ('ob_oa'    ,ob_oa    ,bc% data% ob_oa)
      call put_var  ('ob_ab'    ,ob_ab    ,bc% data% ob_ab)
      call put_var  ('oa_ab'    ,oa_ab    ,bc% data% oa_ab)
      call put_var  ('n_ac'     ,n_ac     ,bc% data% n_ac)
      call put_var  ('o_s_ac'   ,o_s_ac   ,bc% data% o_s_ac)
      call put_var  ('b_s_ac'   ,b_s_ac   ,bc% data% b_s_ac)
      call put_var  ('a_s_ac'   ,a_s_ac   ,bc% data% a_s_ac)
      call put_var  ('o_s2_ac'  ,o_s2_ac  ,bc% data% o_s2_ac)
      call put_var  ('b_s2_ac'  ,b_s2_ac  ,bc% data% b_s2_ac)
      call put_var  ('a_s2_ac'  ,a_s2_ac  ,bc% data% a_s2_ac)
      call put_var  ('ob_s_ac'  ,ob_s_ac  ,bc% data% ob_s_ac)
      call put_var  ('oa_s_ac'  ,oa_s_ac  ,bc% data% oa_s_ac)
      call put_var  ('ob_ac'    ,ob_ac    ,bc% data% ob_ac)
      call put_var  ('oa_ac'    ,oa_ac    ,bc% data% oa_ac)
      call put_var  ('ab_ac'    ,ab_ac    ,bc% data% ab_ac)
      call put_var  ('ob_ob_ac' ,ob_ob_ac ,bc% data% ob_ob_ac)
      call put_var  ('oa_oa_ac' ,oa_oa_ac ,bc% data% oa_oa_ac)
      call put_var  ('ab_ab_ac' ,ab_ab_ac ,bc% data% ab_ab_ac)
      call put_var  ('ob_oa_ac' ,ob_oa_ac ,bc% data% ob_oa_ac)
      call put_var  ('ob_ab_ac' ,ob_ab_ac ,bc% data% ob_ab_ac)
      call put_var  ('oa_ab_ac' ,oa_ab_ac ,bc% data% oa_ab_ac)
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
    subroutine put_vari (name, varid, values)
    character(len=*) ,intent(in) :: name
    integer          ,intent(in) :: varid
    integer          ,intent(in) :: values (:,:,:,:)

      status = nf90_put_var (bc%h% ncid ,varid, values)

      if (status /= NF90_NOERR) then
        call finish ('write_bcor_file: put_var '//trim(name),&
                      trim(nf90_strerror(status))            )
      end if

    end subroutine put_vari
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    subroutine put_var (name, varid, values)
    character(len=*) ,intent(in) :: name
    integer          ,intent(in) :: varid
    real(wp)         ,intent(in) :: values (:,:,:,:)

      status = nf90_put_var (bc%h% ncid ,varid, values)

      if (status /= NF90_NOERR) then
        call finish ('write_bcor_file: put_var '//trim(name),&
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
        write(0,*)   'write_bcor_file: def_dim '//trim(name)//' : size =',size
        call finish ('write_bcor_file: def_dim '//trim(name),&
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
        call finish ('write_bcor_file: def_var '//trim(name),&
                      trim(nf90_strerror(status))            )
      end if

      status = nf90_put_att (bc%h% ncid ,varid, 'longname', longname)

      if (status /= NF90_NOERR) then
        call finish ('write_bcor_file: def_var '//trim(name),&
                      trim(nf90_strerror(status))            )
      end if
    end subroutine def_var
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine write_bcor_file

  !==============================================================================

  subroutine wlidar_bc_bfg (obs)
  !------------------------------------------------------------------------
  ! Radiance bias correction routine to be called before first guess check.
  ! Check for missing entries in bc file, extend file
  ! Apply bias correction to observations.
  !------------------------------------------------------------------------
  type (t_obs_set) ,intent(inout) :: obs  ! observation

    !----------------
    ! local variables
    !----------------
    integer ,parameter       :: mp = 10000      ! max. number of new wind lidars
    type(t_tmp)              :: mis_pe (mp)  ! missing entries in bc file
    type(t_tmp)              :: mis_all(mp)  ! missing entries in bc file
    integer                  :: satid        ! wind lidar satellite id
    integer                  :: n            ! counter / PE
    integer                  :: m            ! counter all PEs
    integer                  :: ib           ! box index
    integer                  :: is           ! spot index
    integer                  :: i, j, k, l   ! loop indices
    integer                  :: ip           ! Lidar Reciever index
    integer                  :: itr          ! Lidar Reciever track direction (ascent/descent)
    integer                  :: ilev_bc   ! Bias correction level
    integer                  :: ilev      ! loop index
    type(t_plat_bc) ,pointer :: plat (:)     ! temporary
    type(t_bc)      ,pointer :: data (:,:,:,:)     ! temporary
    integer                  :: ibc          ! bias correction index
    integer                  :: state        ! state to set if no BC available
    real(wp)                 :: raw_hlos     ! raw observation
    real(wp)                 :: hlat         ! latitude of raw observation
    real(wp)                 :: cor_hlos     ! corrected observation
    real(wp)                 :: c_o, s_o, s_u
    real(wp) :: scale     =  1._wp           ! Scale factor of observation
    real(wp) :: intercept =  0._wp           ! Constant offset

!     bc_levs(1:BC_LEV) = [ 1050, 850, 500, 200, 70, 5 ]
!     bc_levs = bc_levs * 100._wp   ! hPa -> Pa


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
        if (obs% o(ib)% spot(is)% hd% obstype /= OT_WLIDAR) cycle
          obs% o(ib)% spot(is)% bc_index = 0
          do i = obs% o(ib)% spot(is)% o% i + 1,                       &
                 obs% o(ib)% spot(is)% o% i + obs% o(ib)% spot(is)% o% n
             select case (obs% o(ib)% varno(i))
            case (VN_HLOS)
            case default
              cycle
            end select
            satid = obs% o(ib)% spot(is)% ident
!           write(*,*) 'VARNO/SATID WLIDAR: ', VN_HLOS, satid
            do j = 1, size(wlidar_bc% plat)
              if  (satid == wlidar_bc% plat(j)% satid) then
                obs% o(ib)% spot(is)% bc_index = j
                cycle spot
              endif
            end do
            do j = 1, n
              if (mis_pe (j)% satid == satid) then
                obs% o(ib)% spot(is)% bc_index = -j
                cycle spot
              endif
            end do
            n = n + 1
            if (n > mp) call finish('wlidar_bc_bfg','n > mp')
            obs% o(ib)% spot(is)% bc_index = - n
            mis_pe (n)% satid  = satid
            mis_pe (n)% i      = n
            mis_pe (n)% j      = 0
            cycle spot
          end do
        end do spot
      end do
      !---------------------------------------
      ! cross-check missing entries on all PEs
      !---------------------------------------
      m = p_sum (n)
      if (dace% lpio) write(6,'(a,i2,i6,i6)')'    print satid1 ', &
                                                    n,m,mp
      if (m > mp) call finish('wlidar_bc_bfg','m > mp')
      if (m > 0) then
        call p_gather_tmp (mis_pe (1:n), mis_all(1:m), dace% pio)
        if (dace% lpio) then
          k             = 1
          i             = 1
          mis_all(i)% j = k
          mis_pe (k)    = mis_all(i)
outer:    do i = 2, m
            do j = 1, k
              if (mis_all(i)% satid == mis_pe(j)% satid) then
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
          i = size (wlidar_bc% plat)
          plat => wlidar_bc% plat
          data => wlidar_bc% data
          allocate (wlidar_bc% plat (i+k))
          allocate (wlidar_bc% data (nph,ntr,nlv,i+k))
          wlidar_bc% n = i+k
          wlidar_bc% plat(:i) = plat
          wlidar_bc% data(:,:,:,:i) = data
          do j = 1, k
            wlidar_bc% plat (i+j)% satid = mis_pe (j)% satid
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
          deallocate (wlidar_bc% plat)
          deallocate (wlidar_bc% data)
          allocate   (wlidar_bc% plat (m))
          allocate   (wlidar_bc% data (nph,ntr,nlv,m))
        endif
        call p_bcast_plat (wlidar_bc% plat, dace% pio)
        call p_bcast_data (wlidar_bc% data, dace% pio)
        wlidar_bc% n = m
        !-----------------------------------
        ! relate reports to new file entries
        !-----------------------------------
        do ib = 1, size (obs% o)
          if (obs% o(ib)% pe /= dace% pe) cycle
          do is = 1, obs% o(ib)% n_spot
            if (obs% o(ib)% spot(is)% hd% obstype /= OT_WLIDAR) cycle
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
!   write(6,'(a,i2,i2)')'    Apply bc_bfg', biascor_mode, BC_FG
    if (biascor_mode >= BC_FG) then
      state = rept_use (OT_WLIDAR)% use (CHK_BIASCOR)
      do ib = 1, size (obs% o)
        if (obs% o(ib)% pe /= dace% pe) cycle
        do is = 1, obs% o(ib)% n_spot
          if (obs% o(ib)% spot(is)% hd% obstype /= OT_WLIDAR) cycle
          ibc = obs% o(ib)% spot(is)% bc_index
          if (ibc <= 0)                                      cycle
          select case (obs% o(ib)% spot(is)% stret)
          case (0)     ! MIE Receiver
            ip = PH_MIE
          case (1)     ! Rayleigh Receiver
            ip = PH_RAY
          case default ! Unsteady, reserved
            ip = 0
          end select

        do ilev = 1, size(bc_levs) -1
           if (obs% o(ib)% spot(is)% ps <= bc_levs(ilev) .and. obs% o(ib)% spot(is)% ps > bc_levs(ilev+1)) then
             ilev_bc = ilev
           endif
        enddo

!         raw_hlos = -1._wp
          do i = obs% o(ib)% spot(is)% o% i + 1,                       &
                 obs% o(ib)% spot(is)% o% i + obs% o(ib)% spot(is)% o% n
            select case (obs% o(ib)% varno(i))
            case (VN_HLOS)

              if (obs% o(ib)% body(i)% obs_par(1) < 180.0) then
                itr = PH_DES
              else
                itr = PH_ASC
              endif

            if (ip > 0) then
              k = ip
!             if (raw_hlos < 0._wp) raw_hlos = obs% o(ib)% body(i)% o  - &
!                                          obs% o(ib)% body(i)% bc
                raw_hlos = obs% o(ib)% body(i)% o  - &
                           obs% o(ib)% body(i)% bc
                hlat = obs% o(ib)% spot(is)% col% c% dlat
              !--------------------------------
              ! Apply transformation to raw obs
              !--------------------------------
              if (wlidar_bc% data(k,itr,ilev_bc,ibc)% n_ac <  n_required .or. &
                  wlidar_bc% data(k,itr,ilev_bc,ibc)% n_ac <= 0._wp           ) then

                 cor_hlos = raw_hlos
!      write(6,'(a,i6,2X,f8.3,f8.3,f8.3,2X,f12.3,2X,i3)')' Bin in wlidar_bc_bfg HLOS ', &
!            n_required, wlidar_bc% data(k,ilev_bc,ibc)% n_ac,raw_hlos, cor_hlos, obs% o(ib)% spot(is)% ps, ilev_bc
              else
!                s_o = (wlidar_bc% data(k,ibc)% ob_s_ac                               - &
!                       wlidar_bc% data(k,ibc)% o_s_ac * wlidar_bc% data(k,ibc)% b_s_ac    )
!                s_u = (wlidar_bc% data(k,ibc)% o_s2_ac                               - &
!                       wlidar_bc% data(k,ibc)% o_s_ac**2                               )
!                c_o = (wlidar_bc% data(k,ibc)% b_s_ac * wlidar_bc% data(k,ibc)% o_s2_ac - &
!                       wlidar_bc% data(k,ibc)% o_s_ac * wlidar_bc% data(k,ibc)% ob_s_ac   )
!
! New fo WLIDAR correction versus latitude
!
!      write(6,'(a,i6,2X,f10.3)')' Bin in HLOS ', n_required, wlidar_bc% data(k,ibc)% n_ac

                 s_o = (wlidar_bc% data(k,itr,ilev_bc,ibc)% ob_s_ac                               - &
                        wlidar_bc% data(k,itr,ilev_bc,ibc)% o_s_ac * wlidar_bc% data(k,itr,ilev_bc,ibc)% ob_ac    )
                 s_u = (wlidar_bc% data(k,itr,ilev_bc,ibc)% o_s2_ac                               - &
                        wlidar_bc% data(k,itr,ilev_bc,ibc)% o_s_ac**2                               )
                 c_o = (wlidar_bc% data(k,itr,ilev_bc,ibc)% ob_ac * wlidar_bc% data(k,itr,ilev_bc,ibc)% o_s2_ac - &
                        wlidar_bc% data(k,itr,ilev_bc,ibc)% o_s_ac * wlidar_bc% data(k,itr,ilev_bc,ibc)% ob_s_ac   )

!                if (s_u > 0._wp) then
                    scale     = s_o / s_u
                    intercept = c_o / s_u
!                else
!                   scale     = 1
!                   intercept = 0
!                endif
!                cor_hlos = scale * raw_hlos + intercept
                 if (abs(wlidar_bc% data(k,itr,ilev_bc,ibc)% ob_ac) > 0.3_wp) then
                 cor_hlos = raw_hlos - (scale * hlat + intercept)
                 else
                 cor_hlos = raw_hlos
                 endif

!                write(6,'(a,i3,f8.3,f8.3,f15.11,f8.3,f8.3,f10.3,2X,f12.3,2X,i3)')'    print wlidar TRANS ', &
!                      k,raw_hlos, scale, intercept, hlat, cor_hlos, wlidar_bc% data(k,ilev_bc,ibc)% ob_ac, obs% o(ib)% spot(is)% ps, ilev_bc
              endif
              obs% o(ib)% body(i)% bc = cor_hlos - raw_hlos
              obs% o(ib)% body(i)% o  = cor_hlos

!             if (dace% lpio) write(6,'(a,f8.3,f8.3,f15.11,f8.3,f8.3)')'    print wlidar TRANS ', &
!                      raw_hlos, scale, intercept, hlat, cor_hlos
!             if (dace% lpio) write(6,'(a,f8.3,f8.3,f8.3,f12.3)')'    print wlidar TRANS1 ', &
!                      wlidar_bc% data(j,ibc)% ob_s_ac, wlidar_bc% data(j,ibc)% o_s_ac, &
!                      wlidar_bc% data(j,ibc)% ob_ac, wlidar_bc% data(j,ibc)% o_s2_ac
            else
              call decr_use (obs% o(ib)% body(i)% use, &
                             check = CHK_BIASCOR,      &
                             state = state,            &
                             lflag = .true.            )
            endif
              cycle
            case default
            end select
            exit
          end do  ! i
        end do    ! is
      end do      ! ib
    endif         ! biascor_mode
  end subroutine wlidar_bc_bfg

!==============================================================================!

  subroutine wlidar_bc_aan (obs, y)
  !-----------------------------------------------------------------------
  ! Bias correction routine to be called after analysis
  ! Update bias correction statistics.
  ! Write updated correction coefficient file.
  !-----------------------------------------------------------------------
  type (t_obs_set) ,intent(inout) :: obs  ! observation
  type (t_vector)  ,intent(in)    :: y    ! background, observation space

    integer             :: i,ip,itr  ! index
!   integer             :: n         ! number of satellites in statistics
    integer             :: ib        ! box index
    integer             :: is        ! spot index
    integer             :: ibc       ! plane index
    real(wp)            :: bg        ! background
    real(wp)            :: an        ! analysis
    real(wp)            :: o         ! raw observation
    real(wp)            :: hlat      ! latitude of raw observation
    real(wp)            :: eo        ! observation error
    real(wp)            :: eb        ! background  error
    real(wp)            :: dob       ! obs - bg
    real(wp)            :: doa       ! obs - ana
    real(wp)            :: dab       ! ana - bg
    real(wp)            :: dt        ! time since last update of statistics
    real(wp)            :: f         ! weight for accumulated statistics
    integer             :: ier       ! error return flag
    integer             :: ilev_bc   ! Bias correction level
    integer             :: ilev      ! loop index
    integer,    pointer :: iperm (:) ! permutation index array for sorting
    type(t_bc), pointer :: pp    (:,:,:,:) ! pointer to statistics file entries
    type(t_bc), pointer :: p         ! pointer to statistics file entries
    logical             :: first     ! first obs. (hlos) of this spot?

    if (biascor_mode == 0) return

!    write(6,'(a)')'    Bin in wlidar_bc_aan'

    !-----------------
    ! set sums to zero
    !-----------------
    pp => wlidar_bc% data(:,:,:,:)
    pp% n     = 0
    pp% o_s   = 0._wp  ! observation mean
    pp% b_s   = 0._wp  ! background  mean
    pp% a_s   = 0._wp  ! analyses    mean
    pp% o_s2  = 0._wp  ! observation square mean
    pp% b_s2  = 0._wp  ! background  square mean
    pp% a_s2  = 0._wp  ! analyses    square mean
    pp% ob_s  = 0._wp  ! observation x background mean
    pp% oa_s  = 0._wp  ! observation x analyses mean
    pp% ob    = 0._wp  ! obs-bg  mean deviation
    pp% oa    = 0._wp  ! obs-ana mean deviation
    pp% ab    = 0._wp  ! ana-obs mean deviation
    pp% ob_ob = 0._wp  ! obs-bg  variance
    pp% oa_oa = 0._wp  ! obs-ana variance
    pp% ab_ab = 0._wp  ! ana-obs variance
    pp% ob_oa = 0._wp  !
    pp% ob_ab = 0._wp  !
    pp% oa_ab = 0._wp  !
    pp% o_err = 0._wp  !
    pp% b_err = 0._wp  !

    !-----------------------
    ! Update bias statistics
    !-----------------------
    do ib = 1, size (obs% o)
      if (obs% o(ib)% pe /= dace% pe) cycle
      do is = 1, obs% o(ib)% n_spot
        if (obs% o(ib)% spot(is)% hd% obstype /= OT_WLIDAR) cycle
        ibc = obs% o(ib)% spot(is)% bc_index
        if (ibc <= 0)                                       cycle
        first=.true.

        do ilev = 1, size(bc_levs) -1
           if (obs% o(ib)% spot(is)% ps <= bc_levs(ilev) .and. obs% o(ib)% spot(is)% ps > bc_levs(ilev+1)) then
             ilev_bc = ilev
           endif
        enddo

!          write(6,'(a,2x,f12.3,2X,i3)')'    Bin in wlidar_bc_aan HLOS', obs% o(ib)% spot(is)% ps, ilev_bc

        select case (obs% o(ib)% spot(is)% stret)
        case (0)     ! MIE Reciever
          ip = PH_MIE
        case (1)     ! Rayleigh Reciever
          ip = PH_RAY
        case default ! Unsteady, reserved
          ip = 0
        end select
        if (ip > 0) then
        do i = obs% o(ib)% spot(is)% o% i + 1,                       &
               obs% o(ib)% spot(is)% o% i + obs% o(ib)% spot(is)% o% n
           select case (obs% o(ib)% varno(i))
           case (VN_HLOS)
!           write(6,'(a,i2)')'    Bin in wlidar_bc_aan HLOS', ip
            if (obs% o(ib)% body(i)% obs_par(1) < 180.0) then
                itr = PH_DES
              else
                itr = PH_ASC
              endif
            select case (obs% o(ib)% body(i)% use% state)
            case (STAT_PASSIVE, STAT_ACTIVE_0I:STAT_ACCEPTED)
              bg = obs% o(ib)% body(i)% bg
              an = y% s(ib)% x(i)
              o  = obs% o(ib)% body(i)% o - obs% o(ib)% body(i)% bc
              eo = obs% o(ib)% body(i)% eo
              eb = obs% o(ib)% body(i)% eb
            hlat = obs% o(ib)% spot(is)% col% c% dlat
!write(6,'(a,i2,f8.3,f8.3,f8.3)')'    Bin in wlidar_bc_aan FF neu', j, &
!                                bg, an, o
            case default
              cycle
            end select
          case default
             cycle
         end select
           if (first) then
              p    => wlidar_bc% data(ip,itr,ilev_bc,ibc)
              dob  = o  - bg
              doa  = o  - an
              dab  = an - bg
              p% n      = p% n      + 1
!             p% o_s    = p% o_s    + o
              p% o_s    = p% o_s    + hlat
              p% b_s    = p% b_s    + bg
              p% a_s    = p% a_s    + an
!             p% o_s2   = p% o_s2   + (o * o)
              p% o_s2   = p% o_s2   + (hlat * hlat)
              p% b_s2   = p% b_s2   + (bg * bg)
              p% a_s2   = p% a_s2   + (an * an)
!             p% ob_s   = p% ob_s   + (o * bg)
!             p% ob_s   = p% ob_s   + (hlat * bg)
              p% ob_s   = p% ob_s   + (hlat * dob)
!             p% oa_s   = p% oa_s   + (o * an)
              p% oa_s   = p% oa_s   + (hlat * an)
              p% ob     = p% ob     + dob
              p% oa     = p% oa     + doa
              p% ab     = p% ab     + dab
              p% ob_ob  = p% ob_ob  + dob * dob
              p% oa_oa  = p% oa_oa  + doa * doa
              p% ab_ab  = p% ab_ab  + dab * dab
              p% ob_oa  = p% ob_oa  + dob * doa
              p% ob_ab  = p% ob_ab  + dob * dab
              p% oa_ab  = p% oa_ab  + doa * dab
              p% o_err  = p% o_err  + eo
              p% b_err  = p% b_err  + eb
              first = .false.
           end if
        end do ! i
       endif
      end do   ! is
    end do     ! ib

    !----------------
    ! sum up over PEs
    !----------------
    wlidar_bc% data% n      = p_sum (wlidar_bc% data% n)
    wlidar_bc% data% o_s    = p_sum (wlidar_bc% data% o_s)
    wlidar_bc% data% b_s    = p_sum (wlidar_bc% data% b_s)
    wlidar_bc% data% a_s    = p_sum (wlidar_bc% data% a_s)
    wlidar_bc% data% o_s2   = p_sum (wlidar_bc% data% o_s2)
    wlidar_bc% data% b_s2   = p_sum (wlidar_bc% data% b_s2)
    wlidar_bc% data% a_s2   = p_sum (wlidar_bc% data% a_s2)
    wlidar_bc% data% ob_s   = p_sum (wlidar_bc% data% ob_s)
    wlidar_bc% data% oa_s   = p_sum (wlidar_bc% data% oa_s)
    wlidar_bc% data% ob     = p_sum (wlidar_bc% data% ob)
    wlidar_bc% data% oa     = p_sum (wlidar_bc% data% oa)
    wlidar_bc% data% ab     = p_sum (wlidar_bc% data% ab)
    wlidar_bc% data% ob_ob  = p_sum (wlidar_bc% data% ob_ob)
    wlidar_bc% data% oa_oa  = p_sum (wlidar_bc% data% oa_oa)
    wlidar_bc% data% ab_ab  = p_sum (wlidar_bc% data% ab_ab)
    wlidar_bc% data% ob_oa  = p_sum (wlidar_bc% data% ob_oa)
    wlidar_bc% data% ob_ab  = p_sum (wlidar_bc% data% ob_ab)
    wlidar_bc% data% oa_ab  = p_sum (wlidar_bc% data% oa_ab)
    wlidar_bc% data% o_err  = p_sum (wlidar_bc% data% o_err )
    wlidar_bc% data% b_err  = p_sum (wlidar_bc% data% b_err )

    !-----------------------------------------
    ! weight factor for accumulated statistics
    !-----------------------------------------
    dt = days (ana_time - time_cyyyymmddhhmm (wlidar_bc% h% last_date))
    if (wlidar_bc% h% t_decay(1) > 0._wp .and. dt > 0._wp) then
      f = exp ( - dt / wlidar_bc% h% t_decay(1))
    else
      f = 1._wp
    endif

    !-----------------------------------
    ! Update accumulated bias statistics
    !-----------------------------------
    pp => wlidar_bc% data(:,:,:,:)
    !-------------------------------
    ! rescale accumulated statistics
    !-------------------------------
    pp% n_ac       = pp% n_ac     * f
    pp% o_s_ac     = pp% o_s_ac   * pp% n_ac
    pp% b_s_ac     = pp% b_s_ac   * pp% n_ac
    pp% a_s_ac     = pp% a_s_ac   * pp% n_ac
    pp% o_s2_ac    = pp% o_s2_ac  * pp% n_ac
    pp% b_s2_ac    = pp% b_s2_ac  * pp% n_ac
    pp% a_s2_ac    = pp% a_s2_ac  * pp% n_ac
    pp% ob_s_ac    = pp% ob_s_ac  * pp% n_ac
    pp% oa_s_ac    = pp% oa_s_ac  * pp% n_ac
    pp% ob_ac      = pp% ob_ac    * pp% n_ac
    pp% oa_ac      = pp% oa_ac    * pp% n_ac
    pp% ab_ac      = pp% ab_ac    * pp% n_ac
    pp% ob_ob_ac   = pp% ob_ob_ac * pp% n_ac
    pp% oa_oa_ac   = pp% oa_oa_ac * pp% n_ac
    pp% ab_ab_ac   = pp% ab_ab_ac * pp% n_ac
    pp% ob_oa_ac   = pp% ob_oa_ac * pp% n_ac
    pp% ob_ab_ac   = pp% ob_ab_ac * pp% n_ac
    pp% oa_ab_ac   = pp% oa_ab_ac * pp% n_ac
    pp% o_err_ac   = pp% o_err_ac * pp% n_ac
    pp% b_err_ac   = pp% b_err_ac * pp% n_ac
    !-----------
    ! accumulate
    !-----------
    pp% n_ac       = pp% n_ac     + pp% n
    pp% o_s_ac     = pp% o_s_ac   + pp% o_s
    pp% b_s_ac     = pp% b_s_ac   + pp% b_s
    pp% a_s_ac     = pp% a_s_ac   + pp% a_s
    pp% o_s2_ac    = pp% o_s2_ac  + pp% o_s2
    pp% b_s2_ac    = pp% b_s2_ac  + pp% b_s2
    pp% a_s2_ac    = pp% a_s2_ac  + pp% a_s2
    pp% ob_s_ac    = pp% ob_s_ac  + pp% ob_s
    pp% oa_s_ac    = pp% oa_s_ac  + pp% oa_s
    pp% ob_ac      = pp% ob_ac    + pp% ob
    pp% oa_ac      = pp% oa_ac    + pp% oa
    pp% ab_ac      = pp% ab_ac    + pp% ab
    pp% ob_ob_ac   = pp% ob_ob_ac + pp% ob_ob
    pp% oa_oa_ac   = pp% oa_oa_ac + pp% oa_oa
    pp% ab_ab_ac   = pp% ab_ab_ac + pp% ab_ab
    pp% ob_oa_ac   = pp% ob_oa_ac + pp% ob_oa
    pp% ob_ab_ac   = pp% ob_ab_ac + pp% ob_ab
    pp% oa_ab_ac   = pp% oa_ab_ac + pp% oa_ab
    pp% o_err_ac   = pp% o_err_ac + pp% o_err
    pp% b_err_ac   = pp% b_err_ac + pp% b_err

    !--------
    ! rescale
    !--------
    where (pp% n > 0)
      pp% o_s      = pp% o_s      / pp% n
      pp% b_s      = pp% b_s      / pp% n
      pp% a_s      = pp% a_s      / pp% n
      pp% o_s2     = pp% o_s2     / pp% n
      pp% b_s2     = pp% b_s2     / pp% n
      pp% a_s2     = pp% a_s2     / pp% n
      pp% ob_s     = pp% ob_s     / pp% n
      pp% oa_s     = pp% oa_s     / pp% n
      pp% ob       = pp% ob       / pp% n
      pp% oa       = pp% oa       / pp% n
      pp% ab       = pp% ab       / pp% n
      pp% ob_ob    = pp% ob_ob    / pp% n
      pp% oa_oa    = pp% oa_oa    / pp% n
      pp% ab_ab    = pp% ab_ab    / pp% n
      pp% ob_oa    = pp% ob_oa    / pp% n
      pp% ob_ab    = pp% ob_ab    / pp% n
      pp% oa_ab    = pp% oa_ab    / pp% n
      pp% o_err    = pp% o_err    / pp% n
      pp% b_err    = pp% b_err    / pp% n
    endwhere

    where (pp% n_ac > 0)
      pp% o_s_ac   = pp% o_s_ac   / pp% n_ac
      pp% b_s_ac   = pp% b_s_ac   / pp% n_ac
      pp% a_s_ac   = pp% a_s_ac   / pp% n_ac
      pp% o_s2_ac  = pp% o_s2_ac  / pp% n_ac
      pp% b_s2_ac  = pp% b_s2_ac  / pp% n_ac
      pp% a_s2_ac  = pp% a_s2_ac  / pp% n_ac
      pp% ob_s_ac  = pp% ob_s_ac  / pp% n_ac
      pp% oa_s_ac  = pp% oa_s_ac  / pp% n_ac
      pp% ob_ac    = pp% ob_ac    / pp% n_ac
      pp% oa_ac    = pp% oa_ac    / pp% n_ac
      pp% ab_ac    = pp% ab_ac    / pp% n_ac
      pp% ob_ob_ac = pp% ob_ob_ac / pp% n_ac
      pp% oa_oa_ac = pp% oa_oa_ac / pp% n_ac
      pp% ab_ab_ac = pp% ab_ab_ac / pp% n_ac
      pp% ob_oa_ac = pp% ob_oa_ac / pp% n_ac
      pp% ob_ab_ac = pp% ob_ab_ac / pp% n_ac
      pp% oa_ab_ac = pp% oa_ab_ac / pp% n_ac
      pp% o_err_ac = pp% o_err_ac / pp% n_ac
      pp% b_err_ac = pp% b_err_ac / pp% n_ac
    endwhere

    !-----
    ! sort
    !-----
    allocate (iperm (wlidar_bc% n))
    call sort (wlidar_bc% plat% satid, iperm, 1, ier)
    wlidar_bc% plat = wlidar_bc% plat (  iperm)
    wlidar_bc% data = wlidar_bc% data (:,:,:,iperm)
    deallocate (iperm)

    !---------
    ! printout
    !---------
    if (dace% lpio) then
      write(6,'(a)') repeat('-',79)
      write(6,'()')
      write(6,'(a)')      '  Aeolus wind lidar hlos wind bias correction'
      write(6,'()')
      write(6,'(a,a)')    '    last_date = ',wlidar_bc% h% last_date
      write(6,'(a,a)')    '    ana_date  = ',cyyyymmddhhmm(ana_time)
      write(6,'(a,f11.2)')'    t_decay   =',    wlidar_bc% h%    t_decay(1)
      write(6,'(a,f13.4)')'    dt        =',                    dt
      write(6,'(a,f13.4)')'    f         =',                    f
      write(6,'(a,i8)')   '    n         =',sum(wlidar_bc% data% n)
      write(6,'(a,f11.2)')'    n_ac      =',sum(wlidar_bc% data% n_ac)
      write(6,'()')
      if (verbose > 0) then
        !---------------------------------
        ! Verbosity of bias correction:
        ! 1 = print only non-empty entries
        ! 2 = print all
        !---------------------------------
        write(6,'(a)') &
             ' satid      n_ac   BIASCOR     RMSE'
        write(6,'()')
        do i = 1, wlidar_bc% n
          if (verbose > 1 .or. sum (wlidar_bc% data(:,:,:,i)% n_ac) > 0._wp) &
          write(6,'(i6,f11.0,2f9.2)')                         &
                                wlidar_bc% plat(i)% satid,       &
                                wlidar_bc% data(:,:,:,i)% n_ac,      &
                                wlidar_bc% data(:,:,:,i)% bc_a,      &
              sqrt (max (0._wp, wlidar_bc% data(:,:,:,i)% ob_ob_ac   &
                              - wlidar_bc% data(:,:,:,i)% ob_ac ** 2))
        end do
        write(6,'()')
      end if
    endif

    !------------------------------------------
    ! Write updated correction coefficient file
    !------------------------------------------
    call write_bcor_file  (wlidar_bc)
    call destruct_bc_file (wlidar_bc)
  end subroutine wlidar_bc_aan

!==============================================================================
  subroutine read_wlidar_bc_nml
  !------------------------------
  ! read namelist /BIASCOR_WLIDAR/
  !------------------------------
    integer :: ierr

    !-----------------------
    ! Namelist BIASCOR_WLIDAR
    !-----------------------
    namelist /BIASCOR_WLIDAR/ biascor_mode, n_required, t_decay, bc_fallback

    !--------------------------------------------
    ! set defaults depending on 'ga3_biasc_wlidar'
    !--------------------------------------------
    select case (flag_biasc_wlidar)
    case (-1)
      biascor_mode = BC_NOBC
      bc_fallback  = .false.
    case ( 0)
      biascor_mode = - BC_FG
      bc_fallback  = .true.
    case ( 1)
      biascor_mode = - BC_FG
      bc_fallback  = .false.
    end select

    !--------------
    ! read namelist
    !--------------
    if (dace% lpio) then
      call position_nml ('BIASCOR_WLIDAR', status=ierr)
      select case (ierr)
      case (POSITIONED)
#if defined(__ibm__)
        read (nnml ,nml=BIASCOR_WLIDAR, iostat=ierr)
        if (ierr/=0) call finish ('read_wlidar_bc_nml',               &
                                  'ERROR in namelist /BIASCOR_WLIDAR/')
#else
        read (nnml ,nml=BIASCOR_WLIDAR)
#endif
      end select
      !-----------------------------------------------------
      ! adjust 'biascor_mode' depending on 'ga3_biasc_wlidar'
      !-----------------------------------------------------
      if (biascor_mode < 0) then
        select case (flag_biasc_wlidar)
        case (-1)
          biascor_mode = BC_NOBC
        case (0:1)
          biascor_mode = - biascor_mode
        end select
      endif

    endif

    call p_bcast (biascor_mode,  dace% pio)
    call p_bcast (n_required,    dace% pio)
    call p_bcast (t_decay,       dace% pio)
    call p_bcast (bc_fallback,   dace% pio)

    !-------------------------------
    ! Print namelist /BIASCOR_SCATT/
    !-------------------------------
    if (dace% lpio) then
      write(6,'(a)') repeat('-',79)
      write(6,'()')
      write(6,'(a)')      ' Namelist /BIASCOR_WLIDAR/:'
      write(6,'()')
      write(6,'(a,i6)'  ) ' biascor_mode      = ', biascor_mode
      write(6,'(a,l6)'  ) ' bc_fallback       = ', bc_fallback
      write(6,'(a,i6)'  ) ' n_required        = ', n_required
      write(6,'(a,f6.0)') ' t_decay           = ', t_decay
      write(6,'()')
      if (flag_biasc_wlidar < 0 .and.  biascor_mode > 0) then
        write(0,*) '*** Warning: biascor_mode > 0 but flag_biasc_wlidar < 0 ***'
        write(6,*) '*** Warning: biascor_mode > 0 but flag_biasc_wlidar < 0 ***'
        write(6,'()')
      end if
    endif
  end subroutine read_wlidar_bc_nml
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
#define DERIVED type(t_bc),dimension(:,:,:,:)
#define VECTOR
#define RANK 4
#undef  MPI_TYPE
#define p_bcast_DERIVED p_bcast_data
#include "p_bcast.incf"
#undef  DERIVED
#undef  VECTOR
#undef  RANK
#undef  p_bcast_DERIVED
!==============================================================================
end module mo_wlidar_bc
