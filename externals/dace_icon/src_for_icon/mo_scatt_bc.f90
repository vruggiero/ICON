!
!+ Bias correction for individual conventional platforms
!
MODULE mo_scatt_bc
!
! Description:
!   Bias correction for individual conventional platforms.
!   Currently: Scatterometer and altimeter wind speed.
!
! Current Maintainer: DWD, Alexander Cress, Harald Anlauf
!    phone: +49 69 8062 2716
!    fax:   +49 69 8062 3721
!    email: alexander.cress@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_49        2016-10-25 Alexander Cress
!  Online bias correction for SCATT
! V1_50        2017-01-09 Harald Anlauf
!  Rename namelist BIASCOR_SCA to BIASCOR_SCATT, bugfixes, cleanup
! V1_51        2017-02-24 Harald Anlauf
!  read_scatt_bc_nml: warn if biascor_mode > 0 but flag_biasc_scatt < 0
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
  use mo_obs_set,    only: t_obs_set         ! observation data derived type
  use mo_obs_tables, only: rept_use          ! use table entry
  use mo_dec_matrix, only: t_vector          ! vector data type
  use mo_mpi_dace,   only: dace,            &! MPI group info
                           p_sum,           &! sum over PEs
                           p_bcast,         &! generic MPI broadcast routine
                           p_gather          ! generic MPI gather routine
  use mo_run_params, only: ana_time,        &! analysis time
!                          run_type,        &! haupt=0, vor=1, ass=2
                           flag_biasc_scatt  ! scatterometer bias correction
  use mo_biasc_io,   only: t_bcor_head,     &! file header derived type
                           new_bc_head,     &! construct new file header
                           open_bc_read,    &! open  file for reading
                           open_bc_write,   &! open  file for writing
                           close_bc,        &! close file
                           bc_paths,        &! full pathnames
                           bc_obstyp,       &! observation type in files
                           nbcf,            &! number of files in list
                           verbose           ! Verbosity level of bias corr.
  use mo_fdbk_tables,only: VN_FF,           &! wind speed observation flag
                           VN_U,            &! u-component observation flag
                           VN_V,            &! v-component observation flag
                           OT_SCATT          ! Scatterometer report type ID
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
                           NF90_INT,        &! integer   type id
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
! Public entities
!================

  private
  !--------------
  ! derived types
  !--------------
  public :: t_scatt_bc        ! scatterometer bias correction data
  public :: t_bc              ! components of t_scatt_bc
  !------------
  ! subroutines
  !------------
  public :: scatt_bc_init     ! initialize module: read namelist, biascor.coeff.
  public :: scatt_bc_bfg      ! apply bias correction, called before fg-check
  public :: scatt_bc_aan      ! update bias correction coefs., after analysis
  public :: read_bcor_file    ! read scatterometer bias correction file
  public :: read_scatt_bc_nml ! read scatterometer bias correction namelist
  !--------------------------------------------------------
  ! namelist parameters to be set in namelist BIASCOR_SCATT
  !--------------------------------------------------------
  public :: biascor_mode      ! mode used for updating
  public :: t_decay           ! accumulation decay time (days)
  public :: n_required        ! number of entries required for correction
  public :: bc_fallback       ! fallback if biasc-file not present
  public :: BC_NOBC,BC_UP,BC_FG,BC_AN,BC_VARBC ! values for biascor_mode
  !--------------------------------
  ! Obsolete static bias correction
  !--------------------------------
! public :: scat_bc_init      ! Initialize scatterometer bias correction
! public :: scat_bc_bfg       ! Apply bias correction, called before FG check

!=========================
! Derived type definitions
!=========================

  !----------------------------
  ! scatterometer specific data
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
!   real(wp) :: bc_b_err           ! bias correction background error
!   real(wp) :: bc_a_err           ! bias correction analysis   error
  end type t_bc

  !--------------------------
  ! observed quantity indices
  !--------------------------
! integer, parameter :: nob   = 2     ! number of observed quantities
! integer ,parameter :: OB_FF = 1     ! wind speed component
! integer ,parameter :: OB_U  = 2     ! u-wind component

  !-----------------------
  ! platform specific data
  !-----------------------
  type t_plat_bc
    integer          :: satid   = 0    ! scatterometer satellite id
    character(len=8) :: satids  =''
  end type t_plat_bc

  !-----------------------------------
  ! Scatterometer bias correction data
  !-----------------------------------
  type t_scatt_bc
    type(t_bcor_head)        :: h                 ! file header data
    integer                  :: biascor_mode = 0  ! mode used for updating
    integer                  :: n            = 0  ! number of entries
    type (t_plat_bc),pointer :: plat (:)          ! platform specific metadata
    type (t_bc)     ,pointer :: data (:)          ! data
  end type t_scatt_bc

  !-------------------------------------
  ! derived type used for cross-checking
  !-------------------------------------
  type t_tmp
    integer          :: satid =  0    ! scatterometer satellite id
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

!------------------------------------------------------------------------------
  !--------------------------------
  ! Obsolete static bias correction
  !--------------------------------
  !----------------------------------------------
  ! type t_scat_bc: scatterometer bias correction
  !----------------------------------------------
  type t_scat_bc
     integer  :: satid =  0               ! Satellite id
     integer  :: instr = -1               ! Instrument id (currently unused)
     real(wp) :: scale =  1._wp           ! Scale factor of observation
     real(wp) :: const =  0._wp           ! Constant offset
  end type t_scat_bc
  integer               :: nscat = 0      ! Actual number of bias corrections
  integer,    parameter :: mscat = 10     ! Max. number of bias corrections
  type(t_scat_bc), save :: scat_bc(mscat)
!------------------------------------------------------------------------------

!=================
! Module variables
!=================

  type(t_scatt_bc),save :: scatt_bc              ! bias correction data
  !-----------------
  ! namelist entries
  !-----------------
  integer                  :: biascor_mode = BC_NOBC   ! mode used for updating
  real(wp)                 :: t_decay      =   -30._wp ! accumulation decaytime
  integer                  :: n_required   =  5000     ! # of entries required
  logical                  :: bc_fallback  = .false.   ! biasc-file not present

contains
!==============================================================================

  subroutine scatt_bc_init
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
      if (flag_biasc_scatt /= 0 .or. .not. bc_fallback) then
        do i = 1, nbcf
          if (bc_obstyp(i) == OT_SCATT) then
            call read_bcor_file (scatt_bc, bc_paths(i))
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
          call finish ('scatt_bc_init','coefficient file not present !')
        call new_bc_head (scatt_bc% h, OT_SCATT, t_decay=[abs(t_decay)])
        scatt_bc% biascor_mode = biascor_mode
        scatt_bc% n            = 0
        allocate (scatt_bc% plat (0))
        allocate (scatt_bc% data (0))
        if (dace% lpio) write(6,'(a,a,i2)')'    created empty file ', &
                        trim (scatt_bc%h% path), scatt_bc% biascor_mode
      endif
      !-------------------------------
      ! prepare bias correction to use
      !-------------------------------
      select case (biascor_mode)
      case (BC_UP)
        scatt_bc% data% bc_b = 0
        scatt_bc% data% bc_a = 0
      case (BC_FG)
!write(6,'(a,i2,i6,i6)')' Bin in scat init:', BC_FG, n_required
        where (scatt_bc% data% n_ac >= n_required)
          scatt_bc% data% bc_b = scatt_bc% data% ob_ac
          scatt_bc% data% bc_a = scatt_bc% data% ob_ac
        elsewhere
          scatt_bc% data% bc_b = NF90_FILL_FLOAT
          scatt_bc% data% bc_a = NF90_FILL_FLOAT
        endwhere
      case (BC_AN)
        where (scatt_bc% data% n_ac >= n_required)
          scatt_bc% data% bc_b = scatt_bc% data% oa_ac
          scatt_bc% data% bc_a = scatt_bc% data% oa_ac
        elsewhere
          scatt_bc% data% bc_b = NF90_FILL_FLOAT
          scatt_bc% data% bc_a = NF90_FILL_FLOAT
        endwhere
      case (BC_VARBC)
        scatt_bc% data% bc_b = scatt_bc% data% bc_a
      end select
    endif
!   write(6,'(a,i2,i6,i6,i6)')'    scatt_bc_init:', &
!        biascor_mode, size(scatt_bc% data% n_ac), &
!        n_required,size(scatt_bc% data% bc_b)

  end subroutine scatt_bc_init

  !------------------------------------------------------------------------------

  subroutine read_bcor_file (bc, file, inst)
  !----------------------------------------
  ! read scatterometer bias correction file
  !----------------------------------------
  type(t_scatt_bc) ,intent(out) :: bc    ! derived type variable to fill
  character(len=*) ,intent(in)  :: file  ! name of file to read
  logical,optional ,intent(in)  :: inst  ! read instantaneous data as well

    integer  :: i, isat
    logical  :: linst
    real(wp) :: s_o, s_u, c_o, scale, intercept

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
      allocate (bc% data (bc% n))

      !---------------
      ! read variables
      !---------------
      if (bc% n > 0) then
        stanc = 1
        counc = 0
        strnc = 1
        call get_var (bc% plat% satid    ,'satid'   )
        call get_var (bc% data% n_ac     ,'n_ac'    )
        call get_var (bc% data% o_s_ac   ,'o_s_ac'  )
        call get_var (bc% data% b_s_ac   ,'b_s_ac'  )
        call get_var (bc% data% a_s_ac   ,'a_s_ac'  )
        call get_var (bc% data% o_s2_ac  ,'o_s2_ac' )
        call get_var (bc% data% b_s2_ac  ,'b_s2_ac' )
        call get_var (bc% data% a_s2_ac  ,'a_s2_ac' )
        call get_var (bc% data% ob_s_ac  ,'ob_s_ac' )
        call get_var (bc% data% oa_s_ac  ,'oa_s_ac' )
        call get_var (bc% data% ob_ac    ,'ob_ac'   )
        call get_var (bc% data% oa_ac    ,'oa_ac'   )
        call get_var (bc% data% ab_ac    ,'ab_ac'   )
        call get_var (bc% data% ob_ob_ac ,'ob_ob_ac')
        call get_var (bc% data% oa_oa_ac ,'oa_oa_ac')
        call get_var (bc% data% ab_ab_ac ,'ab_ab_ac')
        call get_var (bc% data% ob_oa_ac ,'ob_oa_ac')
        call get_var (bc% data% ob_ab_ac ,'ob_ab_ac')
        call get_var (bc% data% oa_ab_ac ,'oa_ab_ac')
        call get_var (bc% data% o_err_ac ,'o_err_ac')
        call get_var (bc% data% b_err_ac ,'b_err_ac')
        call get_var (bc% data% bc_b     ,'bc_b'    )
        call get_var (bc% data% bc_a     ,'bc_a'    )
        if (linst) then
          call get_var (bc% data% n     ,'n'    )
          call get_var (bc% data% o_s   ,'o_s'  )
          call get_var (bc% data% b_s   ,'b_s'  )
          call get_var (bc% data% a_s   ,'a_s'  )
          call get_var (bc% data% o_s2  ,'o_s2' )
          call get_var (bc% data% b_s2  ,'b_s2' )
          call get_var (bc% data% a_s2  ,'a_s2' )
          call get_var (bc% data% ob_s  ,'ob_s' )
          call get_var (bc% data% oa_s  ,'oa_s' )
          call get_var (bc% data% ob    ,'ob'   )
          call get_var (bc% data% oa    ,'oa'   )
          call get_var (bc% data% ab    ,'ab'   )
          call get_var (bc% data% ob_ob ,'ob_ob')
          call get_var (bc% data% oa_oa ,'oa_oa')
          call get_var (bc% data% ab_ab ,'ab_ab')
          call get_var (bc% data% ob_oa ,'ob_oa')
          call get_var (bc% data% ob_ab ,'ob_ab')
          call get_var (bc% data% oa_ab ,'oa_ab')
          call get_var (bc% data% o_err ,'o_err')
          call get_var (bc% data% b_err ,'b_err')
        endif
        do i = 1, bc% n
          do isat = 1, size(satid_t)
            if (satid_t(isat)% code == bc% plat(i)% satid) then
              bc% plat(i)% satids = satid_t(isat)% mnem
            endif
          enddo
        enddo
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
      allocate (bc% plat (bc% n))
      allocate (bc% data (bc% n))
    endif
    call p_bcast_plat (scatt_bc% plat, dace% pio)
    call p_bcast_data (scatt_bc% data, dace% pio)

    if (dace% lpio .and. scatt_bc% n > 0) then
       write(6,'(a)') "  Scatterometer bias correction"
       write(6,'()')
       write(6,'(a)') &
             '   satid      n_ac   scale*obs + offset'
sat:   do i = 1, scatt_bc% n
          scale     = 1
          intercept = 0
          if (scatt_bc% data(i)% n_ac > 0._wp) then
             s_o = (scatt_bc% data(i)% ob_s_ac                             - &
                    scatt_bc% data(i)% o_s_ac * scatt_bc% data(i)% b_s_ac    )
             s_u = (scatt_bc% data(i)% o_s2_ac                             - &
                    scatt_bc% data(i)% o_s_ac**2                             )
             c_o = (scatt_bc% data(i)% b_s_ac * scatt_bc% data(i)% o_s2_ac - &
                    scatt_bc% data(i)% o_s_ac * scatt_bc% data(i)% ob_s_ac   )
             if (s_u > 0._wp) then
                scale     = s_o / s_u
                intercept = c_o / s_u
             end if
             write(6,'(i8,f11.0,2x,2f9.3)') &
                  scatt_bc% plat(i)% satid, &
                  scatt_bc% data(i)% n_ac,  &
                  scale, intercept
          end if
       end do sat
       write(6,'()')
    end if

  end subroutine read_bcor_file

!------------------------------------------------------------------------------

  subroutine destruct_bc_file (bc)
    type(t_scatt_bc) ,intent(inout) :: bc    ! derived type variable to fill

    type(t_scatt_bc) :: empty_bc             ! empty bias correction
    !----------------------
    ! deallocate components
    !----------------------
    deallocate (bc% plat)
    deallocate (bc% data)

    bc = empty_bc
  end subroutine destruct_bc_file

!------------------------------------------------------------------------------

  subroutine remove_inactive (bc)
    type(t_scatt_bc) ,intent(inout) :: bc
    !----------------------------------------------------
    ! Remove inactive stations when accumulated effective
    ! number of observations falls below threshold
    !----------------------------------------------------
    integer             :: idx(bc% n)       ! Index list
    integer             :: i
    integer             :: n
    type(t_scatt_bc)    :: tmp
    real(wp), parameter :: thresh = 0.5_wp

    if (bc% n <= 0) return

    idx = 0
    n   = 0
    do i = 1, bc% n
       if (bc% data(i)% n_ac > thresh) then
          n      = n + 1
          idx(n) = i
       end if
    end do
    if (n == bc% n) return

    if (dace% lpio) then
       write(6,'(/,A,i0)') '    inactive entries removed: ', bc% n - n
    end if

    allocate (tmp% plat(n))
    allocate (tmp% data(n))
    tmp% plat = bc% plat(idx(1:n))
    tmp% data = bc% data(idx(1:n))
    deallocate (bc% plat, bc% data)
    bc% n     =  n
    bc% plat  => tmp% plat
    bc% data  => tmp% data

  end subroutine remove_inactive

!------------------------------------------------------------------------------

  subroutine write_bcor_file (bc)
  !-----------------------------------------
  ! write scatterometer bias correction file
  !-----------------------------------------
  type(t_scatt_bc) ,intent(inout) :: bc

    integer :: status                            ! NetCDF return value
    integer :: dimid (1)                         ! NetCDF dimension id
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
      call open_bc_write (bc% h, OT_SCATT)

      !---------------
      ! set attributes
      !---------------
      status = nf90_put_att (bc%h% ncid, NF90_GLOBAL, 'biascor_mode', &
                                                       biascor_mode   )
      !------------------
      ! define dimensions
      !------------------
!     call def_dim ('observation', nob,            dimid(1))
      call def_dim ('platform',    size(bc% plat), dimid(1))

      !-----------------
      ! define variables
      !-----------------
      status = nf90_def_var (bc%h% ncid ,'satid' ,NF90_INT, &
                             dimid(1), satid)
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
!     call def_var ('bc_b_err ' ,bc_b_err  ,                  &
!                   'bias correction background error squared')
!     call def_var ('bc_a_err ' ,bc_a_err  ,                 &
!                   'bias corrrection analysis error squared')
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
    integer          ,intent(in) :: values (:)

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
    real(wp)         ,intent(in) :: values (:)

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

  subroutine scatt_bc_bfg (obs)
  !------------------------------------------------------------------------
  ! Radiance bias correction routine to be called before first guess check.
  ! Check for missing entries in bc file, extend file
  ! Apply bias correction to observations.
  !------------------------------------------------------------------------
  type (t_obs_set) ,intent(inout) :: obs  ! observation

    !----------------
    ! local variables
    !----------------
    integer ,parameter       :: mp = 10000   ! max. number of new scatterometers
    type(t_tmp)              :: mis_pe (mp)  ! missing entries in bc file
    type(t_tmp)              :: mis_all(mp)  ! missing entries in bc file
    integer                  :: satid        ! scatterometer satellite id
    integer                  :: n            ! counter / PE
    integer                  :: m            ! counter all PEs
    integer                  :: ib           ! box index
    integer                  :: is           ! spot index
    integer                  :: i, j, k      ! loop indices
    type(t_plat_bc) ,pointer :: plat (:)     ! temporary
    type(t_bc)      ,pointer :: data (:)     ! temporary
    integer                  :: ibc          ! bias correction index
    integer                  :: state        ! state to set if no BC available
    real(wp)                 :: raw_ff       ! raw observation
    real(wp)                 :: raw_u, raw_v ! raw observation
    real(wp)                 :: cor_ff       ! corrected observation
    real(wp)                 :: cor_u, cor_v ! corrected observation
    real(wp)                 :: c_o, s_o, s_u
    real(wp) :: scale     =  1._wp           ! Scale factor of observation
    real(wp) :: intercept =  0._wp           ! Constant offset

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
        if (obs% o(ib)% spot(is)% hd% obstype /= OT_SCATT) cycle
          obs% o(ib)% spot(is)% bc_index = 0
          do i = obs% o(ib)% spot(is)% o% i + 1,                       &
                 obs% o(ib)% spot(is)% o% i + obs% o(ib)% spot(is)% o% n
             select case (obs% o(ib)% varno(i))
            case (VN_FF, VN_U)
            case default
              cycle
            end select
            satid = obs% o(ib)% spot(is)% ident
            do j = 1, size(scatt_bc% plat)
              if  (satid == scatt_bc% plat(j)% satid) then
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
            if (n > mp) call finish('scatt_bc_bfg','n > mp')
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
!     if (dace% lpio) write(6,'(a,i2,i6,i6)')'    print satid1 ', &
!                                                   n,m,mp
      if (m > mp) call finish('scatt_bc_bfg','m > mp')
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
          i = size (scatt_bc% plat)
          plat => scatt_bc% plat
          data => scatt_bc% data
          allocate (scatt_bc% plat (i+k))
          allocate (scatt_bc% data (i+k))
          scatt_bc% n = i+k
          scatt_bc% plat(:i) = plat
          scatt_bc% data(:i) = data
          do j = 1, k
            scatt_bc% plat (i+j)% satid = mis_pe (j)% satid
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
          deallocate (scatt_bc% plat)
          deallocate (scatt_bc% data)
          allocate   (scatt_bc% plat (m))
          allocate   (scatt_bc% data (m))
        endif
        call p_bcast_plat (scatt_bc% plat, dace% pio)
        call p_bcast_data (scatt_bc% data, dace% pio)
        scatt_bc% n = m
        !-----------------------------------
        ! relate reports to new file entries
        !-----------------------------------
        do ib = 1, size (obs% o)
          if (obs% o(ib)% pe /= dace% pe) cycle
          do is = 1, obs% o(ib)% n_spot
            if (obs% o(ib)% spot(is)% hd% obstype /= OT_SCATT) cycle
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
!if (dace% lpio) write(6,'(a)')'    Bin in biascor within 3dvar'
      state = rept_use (OT_SCATT)% use (CHK_BIASCOR)
      do ib = 1, size (obs% o)
        if (obs% o(ib)% pe /= dace% pe) cycle
        do is = 1, obs% o(ib)% n_spot
          if (obs% o(ib)% spot(is)% hd% obstype /= OT_SCATT) cycle
          ibc = obs% o(ib)% spot(is)% bc_index
          if (ibc <= 0)                                      cycle
          raw_ff = -1._wp
          do i = obs% o(ib)% spot(is)% o% i + 1,                       &
                 obs% o(ib)% spot(is)% o% i + obs% o(ib)% spot(is)% o% n
            select case (obs% o(ib)% varno(i))
            case (VN_FF)
!             j = OB_FF
              if (raw_ff < 0._wp) raw_ff = obs% o(ib)% body(i)% o  - &
                                           obs% o(ib)% body(i)% bc
              !--------------------------------
              ! Apply transformation to raw obs
              !--------------------------------
              if (scatt_bc% data(ibc)% n_ac <  n_required .or. &
                  scatt_bc% data(ibc)% n_ac <= 0._wp           ) then
!write(6,'(a)')'    Bin in FF '
                 cor_ff = raw_ff
              else
                 s_o = (scatt_bc% data(ibc)% ob_s_ac                               - &
                        scatt_bc% data(ibc)% o_s_ac * scatt_bc% data(ibc)% b_s_ac    )
                 s_u = (scatt_bc% data(ibc)% o_s2_ac                               - &
                        scatt_bc% data(ibc)% o_s_ac**2                               )
                 c_o = (scatt_bc% data(ibc)% b_s_ac * scatt_bc% data(ibc)% o_s2_ac - &
                        scatt_bc% data(ibc)% o_s_ac * scatt_bc% data(ibc)% ob_s_ac   )
                 if (s_u > 0._wp) then
                    scale     = s_o / s_u
                    intercept = c_o / s_u
                 else
                    scale     = 1
                    intercept = 0
                 end if
                 cor_ff = scale * raw_ff + intercept
              endif

              if (cor_ff <= 0._wp) then
                   cor_ff = max (raw_ff, 1.e-6_wp)
                   call decr_use (obs% o(ib)% body(i)% use, &
                                  check = CHK_BIASCOR,      &
                                  state = state,            &
                                  lflag = .true.            )
              end if
              obs% o(ib)% body(i)% bc = cor_ff - raw_ff
              obs% o(ib)% body(i)% o  = cor_ff
            case (VN_U)
!             j = OB_U
              raw_u  = obs% o(ib)% body(i  )% o  - &
                       obs% o(ib)% body(i  )% bc
              raw_v  = obs% o(ib)% body(i+1)% o  - &
                       obs% o(ib)% body(i+1)% bc
              raw_ff = sqrt (raw_u**2 + raw_v**2)
              !--------------------------------
              ! Apply transformation to raw obs
              !--------------------------------
              if (scatt_bc% data(ibc)% n_ac <  n_required .or. &
                  scatt_bc% data(ibc)% n_ac <= 0._wp           ) then
!write(6,'(a)')'    Bin in U/V '
                 cor_ff = raw_ff
              else
                 s_o = (scatt_bc% data(ibc)% ob_s_ac                               - &
                        scatt_bc% data(ibc)% o_s_ac * scatt_bc% data(ibc)% b_s_ac    )
                 s_u = (scatt_bc% data(ibc)% o_s2_ac                               - &
                        scatt_bc% data(ibc)% o_s_ac**2                               )
                 c_o = (scatt_bc% data(ibc)% b_s_ac * scatt_bc% data(ibc)% o_s2_ac - &
                        scatt_bc% data(ibc)% o_s_ac * scatt_bc% data(ibc)% ob_s_ac   )
                 if (s_u > 0._wp) then
                    scale     = s_o / s_u
                    intercept = c_o / s_u
                 else
                    scale     = 1
                    intercept = 0
                 end if
                 cor_ff = scale * raw_ff + intercept
              endif

              if (cor_ff <= 0._wp) then
                   cor_ff = max (raw_ff, 1.e-6_wp)
                   call decr_use (obs% o(ib)% body(i)% use, &
                                  check = CHK_BIASCOR,      &
                                  state = state,            &
                                  lflag = .true.            )
              end if
!      if (dace% lpio) write(6,'(a,f8.3,f8.3,f15.11,f8.3,f10.3,f10.3,f10.3)')'    print scatt TRANS ', &
!                      raw_ff, scale, intercept, cor_ff,s_o,s_u,c_o
!      if (dace% lpio) write(6,'(a,f8.3,f8.3,f8.3,f8.3)')'    print scatt TRANS1 ', &
!                      scatt_bc% data(j,ibc)% b_s_ac, scatt_bc% data(j,ibc)% o_s2_ac, &
!                      scatt_bc% data(j,ibc)% o_s_ac, scatt_bc% data(j,ibc)% ob_s_ac
          !-----------------------------------------------------------------------
          !--- End transformation
          !----------------------------------------------------------------------

              cor_u  = raw_u * (cor_ff / raw_ff)
              cor_v  = raw_v * (cor_ff / raw_ff)
!      if (dace% lpio) write(6,'(a,f8.3,f8.3,f8.3,f8.3,f8.3,f8.3,f8.3)')'    print scatt BFG1 ', &
!                      raw_u, raw_v, raw_ff, bcor, cor_ff, cor_u, cor_v
              obs% o(ib)% body(i  )% o  = cor_u
              obs% o(ib)% body(i  )% bc = cor_u - raw_u
              obs% o(ib)% body(i+1)% o  = cor_v
              obs% o(ib)% body(i+1)% bc = cor_v - raw_v
!!$      if (dace% lpio) write(6,'(a,f8.3,f8.3,f8.3,f8.3)')'    print scatt BFG1 ', &
!!$                      obs% o(ib)% body(i  )% o, obs% o(ib)% body(i  )% bc, &
!!$                      obs% o(ib)% body(i+1)% o, obs% o(ib)% body(i+1)% bc
              cycle
            case (VN_V) ! Handled by case (VN_U)
              cycle
            case default
            end select
            exit
          end do ! i
        end do   ! is
      end do     ! ib
    endif        ! biascor_mode
  end subroutine scatt_bc_bfg

!==============================================================================!

  subroutine scatt_bc_aan (obs, y)
  !-----------------------------------------------------------------------
  ! Bias correction routine to be called after analysis
  ! Update bias correction statistics.
  ! Write updated correction coefficient file.
  !-----------------------------------------------------------------------
  type (t_obs_set) ,intent(inout) :: obs  ! observation
  type (t_vector)  ,intent(in)    :: y    ! background, observation space

    integer             :: i         ! index
!   integer             :: n         ! number of satellites in statistics
    integer             :: ib        ! box index
    integer             :: is        ! spot index
    integer             :: ibc       ! plane index
    real(wp)            :: bg        ! background
    real(wp)            :: an        ! analysis
    real(wp)            :: o         ! raw observation
    real(wp)            :: eo        ! observation error
    real(wp)            :: eb        ! background  error
    real(wp)            :: dob       ! obs - bg
    real(wp)            :: doa       ! obs - ana
    real(wp)            :: dab       ! ana - bg
    real(wp)            :: dt        ! time since last update of statistics
    real(wp)            :: f         ! weight for accumulated statistics
    integer             :: ier       ! error return flag
    integer,    pointer :: iperm (:) ! permutation index array for sorting
    type(t_bc), pointer :: pp    (:) ! pointer to statistics file entries
    type(t_bc), pointer :: p         ! pointer to statistics file entries
    logical             :: first     ! first obs. (u/v or ff) of this spot?

    if (biascor_mode == 0) return

!   if (run_type < 2)      return    ! update file only in analysis cycle

!   n = size(scatt_bc% plat)

!write(6,'(a)')'    Bin in scatt_bc_aan', n

    !-----------------
    ! set sums to zero
    !-----------------
    pp => scatt_bc% data(:)
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
        if (obs% o(ib)% spot(is)% hd% obstype /= OT_SCATT) cycle
        ibc = obs% o(ib)% spot(is)% bc_index
        if (ibc <= 0)                                      cycle
        first = .true.
        do i = obs% o(ib)% spot(is)% o% i + 1,                       &
               obs% o(ib)% spot(is)% o% i + obs% o(ib)% spot(is)% o% n
           select case (obs% o(ib)% varno(i))
           case (VN_FF)
!             j = OB_FF
!write(6,'(a,i2)')'    Bin in scatt_bc_aan FF', j
            select case (obs% o(ib)% body(i)% use% state)
            case (STAT_PASSIVE, STAT_ACTIVE_0I:STAT_ACCEPTED)
              bg = obs% o(ib)% body(i)% bg
              an = y% s(ib)% x(i)
              o  = obs% o(ib)% body(i)% o - obs% o(ib)% body(i)% bc
              eo = obs% o(ib)% body(i)% eo
              eb = obs% o(ib)% body(i)% eb
!write(6,'(a,i2,f8.3,f8.3,f8.3)')'    Bin in scatt_bc_aan FF neu', j, &
!                                bg, an, o
            case default
              cycle
            end select
           case (VN_U)
!             j = OB_U
!write(6,'(a,i2)')'    Bin in scatt_bc_aan U/V', j
            select case (obs% o(ib)% body(i)% use% state)
            case (STAT_PASSIVE, STAT_ACTIVE_0I:STAT_ACCEPTED)
              bg = sqrt(obs% o(ib)% body(i)% bg**2 + obs% o(ib)% body(i+1)% bg**2)
              an = sqrt(y% s(ib)% x(i)**2 + y% s(ib)% x(i+1)**2)
              o  = sqrt((obs% o(ib)% body(i)% o - obs% o(ib)% body(i)% bc)**2 + &
                        (obs% o(ib)% body(i+1)% o - obs% o(ib)% body(i+1)% bc)**2)
!write(6,'(a,i2,f8.3,f8.3,f8.3)')'    print scatt AANFG ', &
!                                j, bg, an, o
              eo = obs% o(ib)% body(i)% eo
              eb = obs% o(ib)% body(i)% eb
            case default
              cycle
            end select
           case default
              cycle
           end select
           if (first) then
              p    => scatt_bc% data(ibc)
              dob  = o  - bg
              doa  = o  - an
              dab  = an - bg
              p% n      = p% n      + 1
              p% o_s    = p% o_s    + o
              p% b_s    = p% b_s    + bg
              p% a_s    = p% a_s    + an
              p% o_s2   = p% o_s2   + (o * o)
              p% b_s2   = p% b_s2   + (bg * bg)
              p% a_s2   = p% a_s2   + (an * an)
              p% ob_s   = p% ob_s   + (o * bg)
              p% oa_s   = p% oa_s   + (o * an)
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
      end do   ! is
    end do     ! ib

    !----------------
    ! sum up over PEs
    !----------------
    scatt_bc% data% n      = p_sum (scatt_bc% data% n)
    scatt_bc% data% o_s    = p_sum (scatt_bc% data% o_s)
    scatt_bc% data% b_s    = p_sum (scatt_bc% data% b_s)
    scatt_bc% data% a_s    = p_sum (scatt_bc% data% a_s)
    scatt_bc% data% o_s2   = p_sum (scatt_bc% data% o_s2)
    scatt_bc% data% b_s2   = p_sum (scatt_bc% data% b_s2)
    scatt_bc% data% a_s2   = p_sum (scatt_bc% data% a_s2)
    scatt_bc% data% ob_s   = p_sum (scatt_bc% data% ob_s)
    scatt_bc% data% oa_s   = p_sum (scatt_bc% data% oa_s)
    scatt_bc% data% ob     = p_sum (scatt_bc% data% ob)
    scatt_bc% data% oa     = p_sum (scatt_bc% data% oa)
    scatt_bc% data% ab     = p_sum (scatt_bc% data% ab)
    scatt_bc% data% ob_ob  = p_sum (scatt_bc% data% ob_ob)
    scatt_bc% data% oa_oa  = p_sum (scatt_bc% data% oa_oa)
    scatt_bc% data% ab_ab  = p_sum (scatt_bc% data% ab_ab)
    scatt_bc% data% ob_oa  = p_sum (scatt_bc% data% ob_oa)
    scatt_bc% data% ob_ab  = p_sum (scatt_bc% data% ob_ab)
    scatt_bc% data% oa_ab  = p_sum (scatt_bc% data% oa_ab)
    scatt_bc% data% o_err  = p_sum (scatt_bc% data% o_err )
    scatt_bc% data% b_err  = p_sum (scatt_bc% data% b_err )

    !-----------------------------------------
    ! weight factor for accumulated statistics
    !-----------------------------------------
    dt = days (ana_time - time_cyyyymmddhhmm (scatt_bc% h% last_date))
    if (scatt_bc% h% t_decay(1) > 0._wp .and. dt > 0._wp) then
      f = exp ( - dt / scatt_bc% h% t_decay(1))
    else
      f = 1._wp
    endif

    !-----------------------------------
    ! Update accumulated bias statistics
    !-----------------------------------
    pp => scatt_bc% data(:)
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

    where (pp% n_ac > 1.e-37_wp)
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
    allocate (iperm (scatt_bc% n))
    call sort (scatt_bc% plat% satid, iperm, 1, ier)
    scatt_bc% plat = scatt_bc% plat (iperm)
    scatt_bc% data = scatt_bc% data (iperm)
    deallocate (iperm)

    if (dace% lpio) then
      write(6,'(a)') repeat('-',79)
      write(6,'()')
      write(6,'(a)')      '  Scatterometer wind speed bias correction'
    end if

    call remove_inactive (scatt_bc)
    !---------
    ! printout
    !---------
    if (dace% lpio) then
      write(6,'()')
      write(6,'(a,a)')    '    last_date = ',scatt_bc% h% last_date
      write(6,'(a,a)')    '    ana_date  = ',cyyyymmddhhmm(ana_time)
      write(6,'(a,f11.2)')'    t_decay   =',    scatt_bc% h%    t_decay(1)
      write(6,'(a,f13.4)')'    dt        =',                    dt
      write(6,'(a,f13.4)')'    f         =',                    f
      write(6,'(a,i8)')   '    n         =',sum(scatt_bc% data% n)
      write(6,'(a,f11.2)')'    n_ac      =',sum(scatt_bc% data% n_ac)
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
        do i = 1, scatt_bc% n
          if (verbose > 1 .or. (scatt_bc% data(i)% n_ac) > 0._wp) &
          write(6,'(i6,f11.0,2f9.2)')                         &
                                scatt_bc% plat(i)% satid,     &
                                scatt_bc% data(i)% n_ac,      &
                                scatt_bc% data(i)% bc_a,      &
              sqrt (max (0._wp, scatt_bc% data(i)% ob_ob_ac   &
                              - scatt_bc% data(i)% ob_ac ** 2))
        end do
        write(6,'()')
      end if
    endif

    !------------------------------------------
    ! Write updated correction coefficient file
    !------------------------------------------
    call write_bcor_file  (scatt_bc)
    call destruct_bc_file (scatt_bc)
  end subroutine scatt_bc_aan

!==============================================================================
  subroutine read_scatt_bc_nml
  !------------------------------
  ! read namelist /BIASCOR_SCATT/
  !------------------------------
    integer :: ierr

    !-----------------------
    ! Namelist BIASCOR_SCATT
    !-----------------------
    namelist /BIASCOR_SCATT/ biascor_mode, n_required, t_decay, bc_fallback

    !--------------------------------------------
    ! set defaults depending on 'ga3_biasc_scatt'
    !--------------------------------------------
    select case (flag_biasc_scatt)
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
      call position_nml ('BIASCOR_SCATT', status=ierr)
      select case (ierr)
      case (POSITIONED)
#if defined(__ibm__)
        read (nnml ,nml=BIASCOR_SCATT, iostat=ierr)
        if (ierr/=0) call finish ('read_scatt_bc_nml',               &
                                  'ERROR in namelist /BIASCOR_SCATT/')
#else
        read (nnml ,nml=BIASCOR_SCATT)
#endif
      end select
      !-----------------------------------------------------
      ! adjust 'biascor_mode' depending on 'ga3_biasc_scatt'
      !-----------------------------------------------------
      if (biascor_mode < 0) then
        select case (flag_biasc_scatt)
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
      write(6,'(a)')      ' Namelist /BIASCOR_SCATT/:'
      write(6,'()')
      write(6,'(a,i6)'  ) ' biascor_mode      = ', biascor_mode
      write(6,'(a,l6)'  ) ' bc_fallback       = ', bc_fallback
      write(6,'(a,i6)'  ) ' n_required        = ', n_required
      write(6,'(a,f6.0)') ' t_decay           = ', t_decay
      write(6,'()')
      if (flag_biasc_scatt < 0 .and.  biascor_mode > 0) then
        write(0,*) '*** Warning: biascor_mode > 0 but flag_biasc_scatt < 0 ***'
        write(6,*) '*** Warning: biascor_mode > 0 but flag_biasc_scatt < 0 ***'
        write(6,'()')
      end if
    endif
  end subroutine read_scatt_bc_nml
!==============================================================================
  !--------------------------------
  ! Obsolete static bias correction
  !--------------------------------
  subroutine scat_bc_init ()
    !----------------------------------------------------
    ! Read bias-correction coefficient files, definitions
    !----------------------------------------------------
    character(len=8)      :: satellite
    integer               :: sat_id       ! satellite id
    integer               :: instr        ! instrument id
    real(wp)              :: scale        ! scaling factor
    real(wp)              :: const        ! constant offset

    integer               :: ierr
    logical               :: first
#if defined(__ibm__)
    integer               :: ios
#endif

    namelist /SCATTBIASCOR/ satellite, instr, scale, const

    if (dace% lpio) then
      write(6,'(a)') repeat('-',79)
      write(6,'()')
      write(6,'(a)') '  Namelist /SCATTBIASCOR/'
      write(6,'()')
    endif
    first = .true.
    nscat = 0

    do
       !-------------
       ! set defaults
       !-------------
       satellite = ""
       instr     = -1
       scale     =  1._wp
       const     =  0._wp
       !---------------------------------
       ! read namelist, consistency check
       !---------------------------------
       if (dace% lpio) then
          call position_nml ('SCATTBIASCOR' ,lrewind=first ,status=ierr)
          first = .false.
          select case (ierr)
          case (POSITIONED)
#if defined(__ibm__)
             read (nnml ,nml=SCATTBIASCOR, iostat=ios)
             if(ios/=0)call finish('scat_bc_init',                   &
                                   'ERROR in namelist /SCATTBIASCOR/')
#else
             read (nnml ,nml=SCATTBIASCOR)
#endif
          end select
       endif
       !-------------------------------------------
       ! exit if no further namelist group is found
       !-------------------------------------------
       call p_bcast (ierr, dace% pio)
       if (ierr /= POSITIONED) exit
       !-----------------------
       ! broadcast to other PEs
       !-----------------------
       call p_bcast (satellite, dace% pio)
       call p_bcast (instr,     dace% pio)
       call p_bcast (scale,     dace% pio)
       call p_bcast (const,     dace% pio)
       nscat = nscat + 1
       if (nscat > size (scat_bc))                                    &
         call finish ("scat_bc_init",                                 &
                      "too many occurences of namelist /SCATTBIASCOR/")
       !-------------------
       ! Consistency checks
       !-------------------
       if (scale  <= 0._wp) call finish ("scat_bc_init","scale <= 0")
       sat_id = satid (satellite)
       if (sat_id == 0    ) &
         call finish ('scat_bc_init','invalid satellite name: '//trim (satellite))

       scat_bc(nscat) = t_scat_bc (sat_id, instr, scale=scale, const=const)

       if (dace% lpio) then
         write(6,'(A,1x,A)') '  Satellite:', trim (satname (scat_bc(nscat)% satid))
         write(6,'(A,I10)')  '     satid =', scat_bc(nscat)% satid
         write(6,'(A,I10)')  '     instr =', scat_bc(nscat)% instr
         write(6,'(A,F10.4)')'     scale =', scale
         write(6,'(A,F10.4)')'     const =', const
         write(6,'()')
      end if
    end do

  end subroutine scat_bc_init
!------------------------------------------------------------------------------
  subroutine scat_bc_bfg (obs)
    type(t_obs_set) ,intent(inout) :: obs  ! observation
    !---------------------------------------------------------------
    ! Bias correction routine to be called before first guess check.
    ! Apply bias correction to observations.
    !---------------------------------------------------------------
    integer                  :: ib           ! box index
    integer                  :: is           ! spot index
    integer                  :: i, j         ! loop indices
    integer                  :: satid        ! satellite id
    integer                  :: state        ! state to set if no BC available
    real(wp)                 :: raw_ff       ! raw observation
    real(wp)                 :: raw_u, raw_v ! raw observation
    real(wp)                 :: cor_ff       ! corrected observation
    real(wp)                 :: cor_u, cor_v ! corrected observation

    if (nscat == 0) return

    state = rept_use (OT_SCATT)% use (CHK_BIASCOR)

    do ib = 1, size (obs% o)
       if (obs% o(ib)% pe /= dace% pe) cycle
spots: do is = 1, obs% o(ib)% n_spot
          if (obs% o(ib)% spot(is)% hd% obstype /= OT_SCATT) cycle
          do i = obs% o(ib)% spot(is)% o% i + 1,                       &
                 obs% o(ib)% spot(is)% o% i + obs% o(ib)% spot(is)% o% n
             satid = obs% o(ib)% spot(is)% ident
             do j = 1, nscat
                if (satid == scat_bc(j)% satid) exit
                if (j     == nscat            ) cycle spots
             end do
             select case (obs% o(ib)% varno(i))
             !------------------------
             ! Wind speed measurements
             !------------------------
             case (VN_FF)
                raw_ff = obs% o(ib)% body(i)% o  - &
                         obs% o(ib)% body(i)% bc
                !--------------------------------
                ! Apply transformation to raw obs
                !--------------------------------
                cor_ff = scat_bc(j)% scale * raw_ff + scat_bc(j)% const
                if (cor_ff < 0._wp) then
                   cor_ff = raw_ff
                   call decr_use (obs% o(ib)% body(i)% use, &
                                  check = CHK_BIASCOR,      &
                                  state = state,            &
                                  lflag = .true.            )
                end if
                obs% o(ib)% body(i)% o  = cor_ff
                obs% o(ib)% body(i)% bc = cor_ff - raw_ff
             !-------------------------
             ! Wind vector measurements
             !-------------------------
             case (VN_U)
                raw_u  = obs% o(ib)% body(i  )% o  - &
                         obs% o(ib)% body(i  )% bc
                raw_v  = obs% o(ib)% body(i+1)% o  - &
                         obs% o(ib)% body(i+1)% bc
                raw_ff = sqrt (raw_u**2 + raw_v**2)
                !--------------------------------
                ! Apply transformation to raw obs
                !--------------------------------
                cor_ff = scat_bc(j)% scale * raw_ff + scat_bc(j)% const
                if (cor_ff <= 0._wp) then
                   cor_ff = max (raw_ff, 1.e-6_wp)
                   call decr_use (obs% o(ib)% body(i)% use, &
                                  check = CHK_BIASCOR,      &
                                  state = state,            &
                                  lflag = .true.            )
                end if
                cor_u  = raw_u * (cor_ff / raw_ff)
                cor_v  = raw_v * (cor_ff / raw_ff)
                obs% o(ib)% body(i  )% o  = cor_u
                obs% o(ib)% body(i  )% bc = cor_u - raw_u
                obs% o(ib)% body(i+1)% o  = cor_v
                obs% o(ib)% body(i+1)% bc = cor_v - raw_v
             case default
                cycle
             end select
          end do
       end do spots
    end do
  end subroutine scat_bc_bfg
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
#define DERIVED type(t_bc),dimension(:)
#define VECTOR
#define RANK 1
#undef  MPI_TYPE
#define p_bcast_DERIVED p_bcast_data
#include "p_bcast.incf"
#undef  DERIVED
#undef  VECTOR
#undef  RANK
#undef  p_bcast_DERIVED
!==============================================================================
end module mo_scatt_bc
