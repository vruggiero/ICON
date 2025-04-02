!
!+ handle observations with time range indicator for MEC
!
module mo_obs_trange
!
! Description:
!   handle observations with time range indicator
!   i.e. averaged/accumulated/aggregated fields
!   for verification with MEC.
!
! Current Code Owner: DWD, Andreas Rhodin
!    phone: +49 69 8062 2722
!    fax:   +49 69 8062 3721
!    email: andreas.rhodin@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_47        2016-06-06 Andreas Rhodin
!  new module to handle observations with time range indicator in MEC
! V1_48        2016-10-06 Andreas Rhodin
!  bugfix for radiance averaging; check start time for time range variables
! V1_50        2017-01-09 Andreas Rhodin
!  improve printout
! V1_51        2017-02-24 Andreas Rhodin
!  RR verification: check for corresponding expid,runclass
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
!==============================================================================

  !=============
  ! Modules used
  !=============
  !------------------------
  ! general purpose modules
  !------------------------
  use mo_kind,        only: wp                 ! kind parameters
  use mo_time,        only: t_time,           &! date+time data type
                            init_time,        &! set t_time from hh, ...
                            operator(==),     &! compare times
                            operator(/=),     &! compare times
                            operator(>),      &! compare times
                            operator(-),      &! take time difference
                            maxval,           &! maximum time
                            seconds,          &! convert to seconds
                            ihhhmm,           &! convert to hours/minutes
                            cyyyymmddhh        ! convert to character string
  use mo_dace_string, only: char3              ! Conversion: INTEGER -> CHAR (LEN=3)
  use mo_mpi_dace,    only: dace,             &! MPI group info
                            self,             &! MPI group info (COMM_SELF)
                            p_sum              ! generic MPI sum
  use mo_exception,   only: finish             ! abort routine
  use mo_run_params,  only: p_readgrib,       &! PE used to read GRIB
                            path_file          ! derive full file name
  use mo_physics,     only: gacc               ! gravity acceleration
  use mo_algorithms,  only: index              ! sorting routine
  !-------------------------------
  ! observation related data types
  !-------------------------------
  use mo_t_obs,       only: t_obs,            &! observation data type
                            t_spot             ! report data type
  use mo_obs_set,     only: t_obs_block        ! observation data type
  use mo_t_datum,     only: rvind              ! invalid value indicator
  use mo_t_col,       only: nt,               &! number of time ranges,
                            tr                 ! time ranges in minutes
  !-------------------------------------
  ! atmospheric state related data types
  !-------------------------------------
  use mo_atm_grid,    only: t_grid             ! grid type
  use mo_atm_state,   only: t_atm,            &! atmospheric state type
                            allocate           ! allocate field
  use mo_atm_transp,  only: scatter_multi      ! scatter model levels
  use mo_memory,      only: t_m                ! type to hold 3d fields
  !-------------------
  ! GRIB file handling
  !-------------------
  use mo_emos_grib1,  only: t_grib1            ! GRIB record data type
  use mo_grib_handling,only:open_gribfile,    &! GRIB open wrapper
                            close_gribfile,   &!
                            t_inventory,      &! Inventory Data Type
                            get_inventory,    &! read inventory table
                            print_inventory,  &! print inventory
                            p_bcast            ! broadcast inventory via MPI
  use mo_grib_invt,   only: operator(==)       ! compare   inventory entries
  use mo_grib,        only: read_single,      &! read a single-level field
                            set_iope_ens,     &! set I/O PE pattern for ensemble
                            check_iope_ens,   &! Compare I/O PE patterns
!                           SCATTER,          &! flag to scatter fields
                            t_ctr_slot         ! slot container
  use mo_wmo_tables,  only: WMO5_RANGE, WMO5_AVERAGE, WMO5_ACCU, WMO5_DIFF, &
                            WMO8_J_POSITIVE
  !-----------------------
  ! Feedback file handling
  !-----------------------
  use mo_fdbk_tables, only: OT_SYNOP         ,&! SYNOP observation type id
                            VN_TRTR,          &! time period of information
                            VN_PTEND,         &! pressure tendency
                            VN_RR,            &! precipitation amount
                            VN_GUST,          &! wind gust
                            VN_RAD_GL,        &! global solar radiation
                            VN_RAD_DF,        &! diffuse solar radiation
                            VN_RAD_LW,        &! long-wave (downward) radiation
                            VN_TMAX,          &! maximum temperature
                            VN_TMIN,          &! minimum temperature
                            varno              ! varno table
  use mo_t_table,     only: name_value,       &! get name from table entry
                            text_value         ! get description from table

  implicit none

  !================
  ! Public Entities
  !================
  private
  public :: t_trange     ! table of variables with time range indicator
  public :: obs_trange   ! identify observations with time range indicator
  public :: grib_trange  ! derive required time ranges, store in atm.state
  public :: apply_trange ! apply observation operator

  interface grib_trange
    module procedure grib_trange_det1 ! deterministic/scalar version  (1 slot )
    module procedure grib_trange_det  ! deterministic/scalar version  (n slots)
    module procedure grib_trange_ens1 ! ensemble version         (1 time slot )
    module procedure grib_trange_ens  ! ensemble version         (n time slots)
  end interface

!------------------------------------------------------------------------------
  !----------------------------------------------------
  ! table entry for variables with time range indicator
  !----------------------------------------------------
  integer ,parameter :: MF = 256       ! max. file name length
  integer ,parameter :: MS =   8       ! max. number of time slots

  type t_trange
    integer           :: i             ! sort index for nice printout
    !--------------------------------------------
    ! specifications determined from observations
    !--------------------------------------------
    integer           :: varno
    type (t_time)     :: ver_time      ! verification time
    integer           :: rng_h         ! requested time range (h)
    type (t_time)     :: rng_time      ! requested time range
    type (t_time)     :: rng_start     ! start of period
    integer           :: n             ! number of observation entries
    !----------------------------------------
    ! specifications required from GRIB input
    !----------------------------------------
    character(len=24) :: shortnames(2) ! name of parameter required
    character(len=16) :: name          ! name (DACE convention)
    integer           :: erange        ! extended time range indicator
    !-----------------
    ! calculated field
    !-----------------
    logical           :: found         ! required fields are present
    real(wp) ,allocatable :: ptr (:,:,:)
  end type t_trange

!------------------------------------------------------------------------------
  !------------------------------------------------------------
  ! table entry to relate observation varnos to GRIB shortnames
  !------------------------------------------------------------
  type t_obs_grib
    integer           :: varno          ! observation variable number
    integer           :: erange         ! extended range indicator
    character(len=24) :: shortnames(2)  ! related GRIB shortname(s) (DWD or local)
    character(len=16) :: name           ! name (DACE internal convention)
  end type t_obs_grib

  !-------------------------------------------
  ! extended range indicator (DACE convention)
  !-------------------------------------------
  integer, parameter :: R_UNKNOWN = 0
  integer, parameter :: R_AVER    = 1
  integer, parameter :: R_ACCUM   = 2
  integer, parameter :: R_MIN     = 3
  integer, parameter :: R_MAX     = 4
  integer, parameter :: R_DIFF    = 5

  !------------------------------------------------------
  ! table to relate observation varnos to GRIB shortnames
  !------------------------------------------------------
  type(t_obs_grib) ,parameter :: obs_grib(8) = &
      [t_obs_grib (VN_PTEND  ,R_DIFF ,['PS       ','         '],'ps       '),    &
       t_obs_grib (VN_RR     ,R_ACCUM,['TOT_PREC ','         '],'tot_prec '),    &
       t_obs_grib (VN_GUST   ,R_MAX  ,['VMAX_10M ','         '],'vmax_10m '),    &
       t_obs_grib (VN_RAD_GL ,R_AVER ,['ASWDIR_S ','ASWDIFD_S'],'aswdir_s '),    &
       t_obs_grib (VN_RAD_DF ,R_AVER ,['ASWDIFD_S','         '],'aswdifd_s'),    &
       t_obs_grib (VN_RAD_LW ,R_AVER ,['ALWD_S   ','         '],'alwd_s   '),    &
       t_obs_grib (VN_TMAX   ,R_MAX  ,['TMAX_2M  ','         '],'tmax_2m  '),    &
       t_obs_grib (VN_TMIN   ,R_MIN  ,['TMIN_2M  ','         '],'tmin_2m  ')     ]

!==============================================================================
contains
!==============================================================================

  subroutine obs_trange (obs, table)
  type(t_obs)    ,intent(in) :: obs   (:)  ! list of observations
  type(t_trange) ,pointer    :: table (:)  ! list of fields required
  !----------------------------------------------------------
  ! analyse the set of observations,
  ! identify observations with time range indicator,
  ! relate the observations to the corresponding GRIB entries
  !----------------------------------------------------------

    integer                     :: io
    integer                     :: is
    integer                     :: ib
    integer                     :: n, i, j, m
    type(t_trange)              :: tr  (200)
    type(t_trange) ,allocatable :: tr1 (:)
    real(wp)                    :: range
    integer                     :: irange
    integer        ,allocatable :: k (:), l(:)
    type(t_time)                :: tmax

    !-----------------------
    ! loop over observations
    !-----------------------
    n = 0
    do io = 1, size(obs)
      if (obs(io)% pe /= dace% pe) cycle
      do is = 1, obs(io)% n_spot
        if (obs(io)% spot(is)% hd% obstype /= OT_SYNOP) cycle
l1:     do ib = obs(io)% spot(is)% o% i + 1,                    &
                obs(io)% spot(is)% o% i + obs(io)% spot(is)% o% n
          if (obs(io)% body (ib)% lev_typ /= VN_TRTR) cycle
          !---------------------------
          ! check for reasonable entry
          !---------------------------
          range = obs(io)% olev (ib)
          irange = int(range)
          if (range - irange /= 0._wp) then
            write(6,*) &
              'WARNING: obs_trange: time range observation skipped',range
            cycle l1
          endif
          !----------------------------------
          ! check if entry is already present
          !----------------------------------
          do i = 1, n
            if (tr(i)% varno    == obs(io)% varno(ib)          .and. &
                tr(i)% ver_time == obs(io)% spot(is)% hd% time .and. &
                tr(i)% rng_h    == irange                          ) then
              tr(i)% n = tr(i)% n + 1
              cycle l1
            endif
          end do
          !-----------------------
          ! insert new table entry
          !-----------------------
          n = n + 1
          if (n > size(tr)) call finish('obs_trange','n > size(tr)')
            tr(n)% varno    = obs(io)% varno(ib)
            tr(n)% ver_time = obs(io)% spot(is)% hd% time
            tr(n)% rng_h    = irange
            tr(n)% n        = 1
        end do l1
      end do
    end do
    !------------------
    ! gather on each PE
    !------------------
    m = p_sum (n)
    allocate (tr1 (2*m))                    ! factor 2 for extension of table
    call allgatherv_tr (tr(1:n), tr1(1:m))

    !-------------------------
    ! keep only unique entries
    !-------------------------
    n = 0
l2: do j = 1, m
      do i = 1, n
        if (tr1(j)% varno    == tr1(i)% varno    .and. &
            tr1(j)% ver_time == tr1(i)% ver_time .and. &
            tr1(j)% rng_h    == tr1(i)% rng_h          ) then
          tr1(i)% n = tr1(i)% n + tr1(j)% n
          cycle l2
        endif
      end do
      n = n + 1
      tr1(n) = tr1(j)
    end do l2

    !--------------------------------------
    ! calculate full time range information
    !--------------------------------------
    do i = 1, n
      call init_time (tr1(i)% rng_time, hh=tr1(i)% rng_h)
      tr1(i)% rng_start = tr1(i)% ver_time - tr1(i)% rng_time
    end do

    !-----------------------------------------
    ! relate observations to GRIB file entries
    !-----------------------------------------
    do j = 1, n
      tr1(j)% erange     = 0
      tr1(j)% shortnames = ''
      tr1(j)% name       = ''
      do i = 1, size (obs_grib)
        if (obs_grib(i)% varno == tr1(j)% varno) then
          tr1(j)% erange     = obs_grib(i)% erange
          tr1(j)% shortnames = obs_grib(i)% shortnames
          tr1(j)% name       = obs_grib(i)% name
          exit
        endif
      end do
    end do

    !-----------------------------------------------------------
    ! Cross-check entries that depend on multiple grib variables
    !-----------------------------------------------------------
    m = n
    do j = 1, m
       if (tr1(j)% shortnames(2) == "") cycle
       !------------------------------------------------------------
       ! Check for matching primary entry (shortname,range,ver_time)
       !------------------------------------------------------------
       if (any (tr1(1:n)% shortnames(1) == tr1(j)% shortnames(2) .and. &
                tr1(1:n)% rng_h         == tr1(j)% rng_h         .and. &
                tr1(1:n)% ver_time      == tr1(j)% ver_time           )) cycle
       !-----------------------------------
       ! Append entry for required variable
       !-----------------------------------
       do i = 1, size (obs_grib)
          if (obs_grib(i)% shortnames(1) == tr1(j)% shortnames(2)) then
             n = n + 1
             if (n > size(tr1)) call finish('obs_trange','n > size(tr1)')
             tr1(n)             = tr1(j)
             tr1(n)% varno      = obs_grib(i)% varno
             tr1(n)% erange     = obs_grib(i)% erange
             tr1(n)% shortnames = obs_grib(i)% shortnames
             tr1(n)% name       = obs_grib(i)% name
             tr1(n)% n          = 0                       ! artifical (no obs)
             exit
          end if
       end do
       if (i > size (obs_grib)) then
          call finish ('obs_trange', 'missing primary entry for shortname: ' &
                                     // trim (tr1(j)% shortnames(2)))
       end if
    end do

    !---------------------------
    ! sort (for proper printout)
    !---------------------------
    allocate (k (n), l(n))
    tmax        = maxval (tr1(1:n)% ver_time)
    k           = index  (tr1(1:n)% rng_h)
    l           = index  (tr1( k )% varno)
    k           = k  (l)
    l           = index (-seconds (tmax - tr1(k)% ver_time))
    k           = k  (l)
    tr1(1:n)% i = k  ! keep index for later use
    !---------
    ! printout
    !---------
    if (dace% lpio) then
      write(6,*)
      write(6,*) '  Observations with time range indicator:'
      write(6,*)
      write(6,*) '  no.  varno   range  date      entries  shortnames'
      write(6,*)
      do j = 1, n
        i = k(j)
        if (j>1) then
          if (cyyyymmddhh(tr1(i)% ver_time)&
           /= cyyyymmddhh(tr1(m)% ver_time)) write(6,*)
        endif
        m = i
        write (6,'(i6,2x,a8,i5,2x,a,i7,2x,2a12,a)') j, &
               name_value(varno, tr1(i)% varno),       &
                                 tr1(i)% rng_h   ,     &
                     cyyyymmddhh(tr1(i)% ver_time),    &
                                 tr1(i)% n,            &
                                 tr1(i)% shortnames,   &
          trim(text_value(varno, tr1(i)% varno))
      end do
      write(6,*)
    endif

    !----------------
    ! return argument
    !----------------
    allocate (table (n))
#ifdef __SX__
    do i = 1, n
    table(i) = tr1(i)
    enddo
#else
    table = tr1(1:n)
#endif

  end subroutine obs_trange

!==============================================================================

  subroutine grib_trange_det1 (fname, table, state, ref_time, &
                               expid, suffix, verbose, ctr    )
  character(len=*) ,intent(in)           :: fname     ! grib file name
  type(t_trange)   ,intent(inout)        :: table (:) ! table
  type(t_atm)      ,intent(inout)        :: state (:) ! atmospheric state
  type(t_time)     ,intent(in)           :: ref_time  ! fc start time
  integer          ,intent(in) ,optional :: expid     ! expid to read only
  character(len=*) ,intent(in) ,optional :: suffix    ! file name suffix
  logical          ,intent(in) ,optional :: verbose   ! steering of printout
  type(t_ctr_slot) ,intent(in) ,optional :: ctr       ! inventory container
  target :: state
  !----------------------------------------------------------------------
  ! get averaged/accumulated/aggregated fields available for verification
  ! for deterministic run and 1 input-file
  !----------------------------------------------------------------------
    type(t_atm) ,pointer           :: st (:,:)    ! atmospheric state
    character(len=MF)              :: fnames(1)   ! grib file name

    fnames                  =  fname
    st (1:1,1:size (state)) => state
    call grib_trange_ens (fnames, table, st, ref_time, expid=expid,   &
                          suffix=suffix, verbose=verbose, ctr=(/ctr/) )

  end subroutine grib_trange_det1

!==============================================================================

  subroutine grib_trange_det (fnames, table, state, ref_time, &
                              expid, suffix, verbose, ctr     )
  character(len=*) ,intent(in)           :: fnames(:) ! grib file name
  type(t_trange)   ,intent(inout)        :: table (:) ! table
  type(t_atm)      ,intent(inout)        :: state (:) ! atmospheric state
  type(t_time)     ,intent(in)           :: ref_time  ! fc start time
  integer          ,intent(in) ,optional :: expid     ! expid to read only
  character(len=*) ,intent(in) ,optional :: suffix    ! file name suffix
  logical          ,intent(in) ,optional :: verbose   ! steering of printout
  type(t_ctr_slot) ,intent(in) ,optional :: ctr(:)    ! inventory container
  target :: state
  !----------------------------------------------------------------------
  ! get averaged/accumulated/aggregated fields available for verification
  ! for deterministic run and multiple input-files (time slots)
  !----------------------------------------------------------------------
    type(t_atm) ,pointer           :: st (:,:)   ! atmospheric state

    st (1:1,1:size (state)) => state
    call grib_trange_ens (fnames, table, st, ref_time, expid=expid, &
                          suffix=suffix, verbose=verbose, ctr=ctr   )

  end subroutine grib_trange_det

!==============================================================================

  subroutine grib_trange_ens1 (fname, table, state, ref_time, members, &
                               expid, suffix, verbose, ctr             )
  !----------------------------------------------------------------------
  ! get averaged/accumulated/aggregated fields available for verification
  ! parallelized ensemble version for 1 input-file / member
  !----------------------------------------------------------------------
  character(len=*) ,intent(in)           :: fname        ! grib file names
  type(t_trange)   ,intent(inout)        :: table  (:)   ! table
  type(t_atm)      ,intent(inout)        :: state  (:,:) ! atmospheric state
  type(t_time)     ,intent(in)           :: ref_time     ! fc start time
  integer          ,intent(in) ,optional :: members(:)   ! renumbering
  integer          ,intent(in) ,optional :: expid        ! expid to read only
  character(len=*) ,intent(in) ,optional :: suffix       ! file name suffix
  logical          ,intent(in) ,optional :: verbose      ! steering of printout
  type(t_ctr_slot) ,intent(in) ,optional :: ctr          ! inventory container

    character(len=MF)              :: fnames(1)   ! grib file name

    fnames                  =  fname
    call grib_trange_ens (fnames, table, state, ref_time, members, &
                          expid, suffix, verbose, ctr=(/ctr/)      )

  end subroutine grib_trange_ens1

!==============================================================================
  subroutine grib_trange_ens (fnames, table, state, ref_time, members, &
                              expid, suffix, verbose, ctr              )
  character(len=*) ,intent(in)           :: fnames (:)   ! grib file names
  type(t_trange)   ,intent(inout)        :: table  (:)   ! table
  type(t_atm)      ,intent(inout)        :: state  (:,:) ! atmospheric state
  type(t_time)     ,intent(in)           :: ref_time     ! fc start time
  integer          ,intent(in) ,optional :: members(:)   ! renumbering
  integer          ,intent(in) ,optional :: expid        ! expid to read only
  character(len=*) ,intent(in) ,optional :: suffix       ! file name suffix
  logical          ,intent(in) ,optional :: verbose      ! steering of printout
  type(t_ctr_slot) ,intent(in) ,optional :: ctr(:)       ! inventory container
  target :: state
  !----------------------------------------------------------------------
  ! get averaged/accumulated/aggregated fields available for verification
  ! parallelized ensemble version for multiple file (time slots) / member
  !----------------------------------------------------------------------
    integer                    :: ns            ! number of files (time slots)
    integer                    :: j, l          ! indices
    integer                    :: i_s, i_e      ! counter
    integer                    :: n             ! inventory size
    integer                    :: k             ! ensemble index
    integer                    :: m             ! ensemble member
    integer                    :: nens          ! ensemble size
    integer                    :: i             ! start index of chunk to read
    integer                    :: nm            ! size        of chunk to read
    integer                    :: ie            ! ensemble member index in chunk
    integer                    :: ke            ! ensemble member index
    integer                    :: nchunk        ! number of chunks to read
    integer                    :: mcsize        ! maximum chunk size
    integer       ,allocatable :: pio(:)        ! PEs to use for reading
    integer                    :: idx           ! index of field in in t_atm
    integer                    :: ix(13)        ! indices
    real(wp)                   :: sj, sk        ! time periods (seconds)
    type(t_time)               :: last_time(13) ! times to concatenate ranges
    type(t_time)               :: next_time     ! time  to concatenate ranges
    logical       ,allocatable :: msk  (:)      ! mask
    logical                    :: verb          ! local copy of 'verbose'
    integer                    :: exp           ! expid to read only
    integer                    :: itr, ntr      ! time range count
    integer                    :: lb (4)        ! lower bounds of grib fields
    integer                    :: ub (4)        ! upper bounds of grib fields
    integer                    :: lbg(4)        ! lower bounds of grib fields
    integer                    :: ubg(4)        ! upper bounds of grib fields
    integer                    :: rng (size (state(1,1)% m(1)% i% tr))
    type(t_grib1)              :: grib          ! GRIB record
    type(t_inventory) ,pointer :: invt (:)      ! GRIB file inventory
    type(t_grid)      ,pointer :: grid          ! grid description
    real(wp)          ,pointer :: ptr(:,:,:,:)  ! temporary storage
    type(t_m)     ,allocatable :: t    (:)      ! temporary storage
    !-----------------------------------------------------------
    ! Container for handling ensemble file names and inventories
    !-----------------------------------------------------------
    type t_gribfile
       character(len=MF)          :: name (MS) ! file names / member
       type(t_inventory), pointer :: invt (:) => NULL ()
    end type t_gribfile
    type(t_gribfile), allocatable :: files(:)

    verb = .true. ;if (present(verbose)) verb = verbose
    exp  = -1     ;if (present(expid  )) exp  = expid
    ns   = size (fnames)
    nens = size (state, dim=1)

    if (ns > MS) call finish('grib_trange_ens','increase MS (number of time slots)')
    allocate (files(nens))

    do k = 1, nens
       m = k; if (present(members)) m = members(k)
       if (nens > 1) then
         do l = 1, ns
           files(k)% name(l) = path_file ('',fnames(l), iens=m, suffix=suffix)
         end do
       else
         do l = 1, ns
           files(k)% name(l) = path_file ('',fnames(l))
         end do
       end if
    end do

    allocate (pio(nens))
    call set_iope_ens (pio, nchunk, mcsize, stride=-1)  ! Optimized I/O stride

    if (present(ctr)) then
      !-------------------
      ! concat inventories
      !-------------------
      if (nens > 1) then
        do l = 1, ns
          call check_iope_ens( ctr(ns)% c(:)% pe, pio,         &
                              'grib_trange_ens ns='//char3(ns) )
        end do
      end if
      if ( size(ctr) /= ns ) &
        call finish('grib_trange_ens','number of slots do not match')
      allocate (msk(ns))
      do k = 1, nens
        if (dace% pe /= pio(k)) cycle
        do l = 1, ns
          if ( trim(files(k)% name(l)) /= trim(ctr(l)% c(k)% name) ) &
            call finish('grib_trange_ens','name mismatch: '// &
                       & trim(files(k)%name(l))//' /= '//trim(ctr(l)%c(k)%name))
          if (.not.associated(ctr(l)% c(k)% invt)) &
            call finish('grib_trange_ens','inventory not associated')
          if (size(ctr(l)% c(k)% invt) == 0) &
            call finish('grib_trange_ens','empty inventory')

          if (l == 1) then
            msk(l) = .false.
            n = ctr(l)% c(k)% size
          else
            if (ctr(l)% c(k)% name == ctr(l-1)% c(k)% name) then
              msk(l) = .true.
            else
              msk(l) = .false.
              n = n + ctr(l)% c(k)% size
            end if
          end if
        end do

        allocate (files(k)% invt(n))
        i_s = 0
        i_e = 0
        do l = 1, ns
          if (msk(l)) cycle
          i_s = i_e + 1
          i_e = i_e + ctr(l)% c(k)% size
          files(k)% invt(i_s:i_e) = ctr(l)% c(k)% invt
        end do

       ! if ( nens > 1 ) then
       !   m = k; if (present(members)) m = members(k)
       !   if (.not.all( files(k)% invt(:)% en% no == m )) &
       !     call finish('grib_trange_ens','wrong member in inventory')
       ! end if

      end do
      deallocate (msk)
    else
      !-----------------------------------------------------------------
      ! distribute read of inventories, 1 member per chosen PE at a time
      !-----------------------------------------------------------------
      do i = 1, nens, mcsize              ! loop over bunches to read
        nm = min (nens-i+1, mcsize)       ! number of members to read
        ie = -1
        ke = -1
        if (any (pio(i:i+nm-1) == dace% pe)) then
          ie = minloc (abs (pio(i:i+nm-1) - dace% pe),1) ! index within chunk
          ke = i - 1 + ie                                ! member to read on this PE
          !print *, "### grib_trange_ens: p_pe,ie,ke =", dace% pe, ie, ke
        endif
!       if (dace% lpio) then
!         write(6,*)
!         write(6,*) ' reading inventory of members',i,' to',i+nm-1
!         write(6,*)
!       endif
        if (ie > 0) then
          do l = 1, ns
            if (l > 1) then
              if (files(ke)% name(l) == files(ke)% name(l-1)) cycle
            endif
            call get_inventory (files(ke)% invt, files(ke)% name(l),               &
                               pio=pio(ke), comm=self% comm, ifile=l, append=.true.)
          end do
        end if
      end do
    end if
    !---------------------------------
    ! broadcast inventories to all PEs
    !---------------------------------
    do k = 1, nens
      call p_bcast (files(k)% invt, pio(k))
    end do
    !-----------------------------------------
    ! Select entries with time range indicator
    !-----------------------------------------
    do k = 1, nens
      invt => files(k)% invt
      allocate (msk (size (invt)))
      msk = .false.
      do i = 1, size (invt)
        if (invt(i)% ti% ref_time /= ref_time) cycle
        do j = 1, size(table)
          if (exp > 0 .and. exp /= invt(i)% pa% expid) cycle
          if (invt(i)% pa% iname     == table(j)% name         ) msk(i) = .true.
          if (invt(i)% pa% shortname == table(j)% shortnames(1)) msk(i) = .true.
          if (invt(i)% pa% shortname == table(j)% shortnames(2)) msk(i) = .true.
        end do
      end do
      n = count (msk)
      allocate  (files(k)% invt (n))
      if (n > 0) files(k)% invt = pack (invt, msk)
      deallocate (invt)
      deallocate (msk)
    end do

    invt => files(1)% invt
    n = size (invt)
    if (verb .and. dace% lpio) then
      write (6,*)
      write (6,*) '  GRIB records with time ranges considered: ',&
           trim(files(1)% name(1)),' start=',cyyyymmddhh(ref_time),'  count=',n
      do l = 2, ns
        write (6,*) '                                            ',&
           trim(files(1)% name(l)),' start=',cyyyymmddhh(ref_time),'  count=',n
      end do
      write (6,*)
    endif
    if (n == 0) goto 99
    if (verb) call print_inventory (invt)
    if (verb) call print_inventory (invt, liname=.true.)
    if (verb .and. dace% lpio) write (6,*)
    !-------------------------------------------------------
    ! Consistency check of entries with time range indicator
    !-------------------------------------------------------
    allocate (msk(n))
    do k = 2, nens
       if (size (files(k)% invt) /= n) then
          if (verb .and. dace% lpio) then
            write(6,*) 'Inventory size MISMATCH: ', &
            trim(files(k)% name(1)),'  count=', size (files(k)% invt)
            do l = 2, ns
              write(6,*) '                         ', &
              trim(files(k)% name(l)),'  count=', size (files(k)% invt)
            end do
          end if
          goto 999
       end if
       call same_order (files(k)% invt, invt, i)
       msk(:) = ( files(k)% invt% pa == invt% pa & ! parameter
            .and. files(k)% invt% ti == invt% ti & ! time
            .and. files(k)% invt% lv == invt% lv ) ! level
       if (count (msk) == n .and. i == 0) then
          cycle
       else
          if (verb .and. dace% lpio) then
            write(6,*) 'Inventory list MISMATCH: ', &
            trim(files(k)% name(1)),'  count (msk)=', count (msk)
            do l = 2, ns
              write(6,*) '                         ', &
              trim(files(k)% name(1)),'  count (msk)=', count (msk)
            end do
          end if
       end if
999    continue
       if (verb .and. dace% lpio) write (6,*)
       if (verb) call print_inventory (files(k)% invt)
       if (verb) call print_inventory (files(k)% invt, liname=.true.)
       goto 99
    end do
    deallocate (msk)

    !------------------
    ! read GRIB entries
    !------------------
    grid => state(1,1)% grid
    invt => files(1)% invt
    lb = grid% lb
    ub = grid% ub
    if (any (pio == dace% pe)) then
       lbg = grid% lbg
       ubg = grid% ubg
    else
       lbg = 1
       ubg = 0
    end if
    n = size (invt)
    allocate (t(nens))
    do k = 1, nens
       allocate (t(k)% ptr(lbg(1):ubg(1), lbg(2):ubg(2), n, lbg(4):ubg(4)))
    end do

    do i = 1, nens, mcsize              ! loop over bunches to read
      nm = min (nens-i+1, mcsize)       ! number of members to read
      ie = -1
      ke = -1
      if (any (pio(i:i+nm-1) == dace% pe)) then
        ie = minloc (abs (pio(i:i+nm-1) - dace% pe),1) ! index within chunk
        ke = i - 1 + ie                                ! member to read on this PE
      endif
      if (dace% lpio .and. nens > 1) then
        write(6,*)
        write(6,*) ' reading time range fields, members',i,' to',i+nm-1
      endif
      if (ke > 0) then
        l = -1
        do j = 1, size (invt)
          if (l /= files(ke)% invt(j)% ifile) then
            if (l > 0) call close_gribfile (grib)
            l = files(ke)% invt(j)% ifile
            call open_gribfile (grib, files(ke)% name(l), 'r')
          endif
          call read_single (t(ke)% ptr(:,:,j:j,:), files(ke)% invt, &
                            grid, grib, pos=j, pio=dace% pe,        &
                            scanmode=WMO8_J_POSITIVE, pmode=0       )
        end do
        call close_gribfile (grib)
      end if
    end do
    !---------------
    ! Scatter fields
    !---------------
    if (dace% lpio) then
      write(6,*) ' scatter fields'
    endif
    n = size (invt)
    do k = 1, nens
       ptr => t(k)% ptr
       allocate (t(k)% ptr(lb(1):ub(1), lb(2):ub(2), n, lb(4):ub(4)))
       call scatter_multi (ptr, t(k)% ptr, grid%dc, pio(k))
       deallocate (ptr)
    end do

    if (dace% lpio) then
      write(6,*) ' calculate required fields'
    endif
    n = size (invt)
    do k = 1, nens
     ptr => t(k)% ptr
     !------------------------------
     ! calculate the required fields
     !------------------------------
     do i = 1, size (table)
      table(i)% found = .false.
      !---------------------------------------
      ! take over without further calculations
      !---------------------------------------
      do j = 1, n
        if ((invt(j)% pa% iname     == table(i)% name  .or.   &
             invt(j)% pa% shortname == table(i)% shortnames(1)) .and. &
             invt(j)% ti% ver_time  == table(i)% ver_time       .and. &
             invt(j)% ti% rng_time  == table(i)% rng_time             ) then
!         if (dace% pe == p_readgrib) &
!           write(6,*) i, j, 'found: ', invt(j)% pa% shortname
          allocate (table(i)% ptr (lb(1):ub(1), lb(2):ub(2), lb(4):ub(4)))
          table(i)% ptr   = ptr (:,:,j,:)
          table(i)% found = .true.
        endif
      end do
      if (table(i)% found) cycle
      !----------------------------
      ! take difference of 2 fields
      !----------------------------
diff: do j = 1, n
        if ((invt(j)% pa% iname     == table(i)% name  .or.   &
             invt(j)% pa% shortname == table(i)% shortnames(1)) .and. &
             invt(j)% ti% ver_time  == table(i)% ver_time             ) then
          do l = 1, n
            if ((invt(l)% pa% iname     == table(i)% name  .or.   &
                 invt(l)% pa% shortname == table(i)% shortnames(1)) .and. &
                 invt(l)% ti% ver_time  == table(i)% rng_start      .and. &
                 invt(l)% ti% range     == invt (j)% ti% range      .and. &
                 invt(l)% pa% expid     == invt (j)% pa% expid      .and. &
                 invt(l)% pa% runclass  == invt (j)% pa% runclass         ) then
              if (invt(j)% ti% ver_time - invt(j)% ti% rng_time /= &
                  invt(l)% ti% ver_time - invt(l)% ti% rng_time    ) cycle
              select case (invt(l)% ti% range)
              case (WMO5_ACCU, WMO5_DIFF)
!               if (dace% pe == p_readgrib) &
!                 write(6,*) i, j, 'diff : ', invt(j)% pa% shortname
                allocate (table(i)% ptr (lb(1):ub(1), lb(2):ub(2), lb(4):ub(4)))
                table(i)% ptr   = ptr (:,:,j,:) - ptr (:,:,l,:)
                table(i)% found = .true.
                exit diff  ! for multiple matching combinations
              case (WMO5_AVERAGE)
!               if (dace% pe == p_readgrib) &
!                 write(6,*) i, j, l, 'aver : ', invt(j)% pa% shortname
                allocate (table(i)% ptr (lb(1):ub(1), lb(2):ub(2), lb(4):ub(4)))
                sj              = seconds (invt(j)% ti% rng_time)
                sk              = seconds (invt(l)% ti% rng_time)
                table(i)% ptr   = (sj * ptr (:,:,j,:) - sk * ptr (:,:,l,:)) / (sj - sk)
                table(i)% found = .true.
                exit diff  ! for multiple matching combinations
              case (WMO5_RANGE)
                if (dace% pe == p_readgrib) &
                     write(6,*) l, trim (invt(l)% pa% shortname),           &
                     ' : cannot derive range :',                            &
                     int (ihhhmm (invt(l)% ti% ver_time - ref_time)), "..", &
                     int (ihhhmm (invt(j)% ti% ver_time - ref_time))
              case default
                if (dace% pe == p_readgrib) &
                     write(6,*) l, trim (invt(l)% pa% shortname), &
                     ' : unknown range : ', invt(l)% ti% range
              end select
            endif
          end do
        endif
      end do diff
      if (table(i)% found) cycle
      !-------------------
      ! take sum of fields
      !-------------------
      j = 1
      last_time(j) = table(i)% ver_time
      ix(j)        = 0
      do
        ix(j) = ix(j) + 1
        if (ix(j) > n) then
          j = j - 1
          if (j < 1) exit
          cycle
        endif

        if ((invt(ix(j))% pa% iname     == table(i)% name  .or.   &
             invt(ix(j))% pa% shortname == table(i)% shortnames(1)) .and. &
             invt(ix(j))% ti% ver_time  == last_time(j)             .and. &
             invt(ix(j))% ti% range     == WMO5_RANGE                     ) then

          next_time = last_time(j) - invt(ix(j))% ti% rng_time
          if (next_time == last_time(j)) cycle
          if (next_time == table(i)% rng_start) then
            table(i)% found = .true.
            allocate (table(i)% ptr (lb(1):ub(1), lb(2):ub(2), lb(4):ub(4)))
            table(i)% ptr = ptr (:,:,ix(1),:)
            do l = 2, j
              select case (table(i)% erange)
              case (R_MIN)
                table(i)% ptr = min (table(i)% ptr, ptr (:,:,ix(l),:))
              case (R_MAX)
                table(i)% ptr = max (table(i)% ptr, ptr (:,:,ix(l),:))
              case default
                call finish('grib_trange','invalid erange')
              end select
            end do
          else if  (next_time > table(i)% rng_start) then
            if (j < size(ix)) then
              j = j + 1
              last_time(j) = next_time
              ix(j)        = 0
              cycle
            else
              cycle
            endif
          else
            cycle
          endif
        endif
      end do

     end do

     !----------------------------
     ! insert in atmospheric state
     !----------------------------
     ntr = size (state(1,1)% m(1)% i% tr)
     deallocate (ptr)
     allocate (ptr (lb(1):ub(1), lb(2):ub(2), 1:ntr, lb(4):ub(4)))
     do l = 1, size (state,2)
      !-------------------------------
      ! print out info on observations
      !-------------------------------
      if (verb .and. dace% lpio .and. k==1) then
        write(6,*)
        write(6,*) '  Observations with time range indicator at ',&
                      cyyyymmddhh(state(k,l)% time),' :',           &
                      count (table(:)% ver_time == state(k,l)% time)
        write(6,*)
        write(6,*) '  no.  varno   range    shortnames        resolved'
        do j = 1, size (table)
          i = table(j)% i      ! sort for nice printout
          if (table(i)% ver_time /= state   (k,l)% time) cycle
          write (6,'(i6,2x,a8,i5,2x,2x,2a12,l2,1x,a)') j,&
                 name_value(varno, table(i)% varno),     &
                                   table(i)% rng_h,      &
                                   table(i)% shortnames, &
                                   table(i)% found,      &
            trim(text_value(varno, table(i)% varno))
        end do
        write(6,*)
      endif
      !-------
      ! insert
      !-------
      do j = 1, size (obs_grib)
        itr = 0
        if (obs_grib(j)% name == '') cycle
        do i = 1, size (table)
          if ( .not.                         table(i)% found         ) cycle
          if ((obs_grib(j)% name          /= table(i)% name) .and.   &
              (obs_grib(j)% shortnames(1) /= table(i)% shortnames(1))) cycle
          if ( state (k,l)% time          /= table(i)% ver_time      ) cycle
          itr = itr + 1
          if (itr > ntr) call finish('grib_trange','itr > ntr')
          ptr (:,:,itr,:) = table(i)% ptr
          rng     (itr)   = table(i)% rng_h * 60
        end do
        if (itr > 0) then
          call allocate (state(k,l), obs_grib(j)% name, trange=rng(1:itr), idx=idx)
#ifdef _CRAYFTN
!DIR$ IVDEP
#endif
          state(k,l)% m(idx)% ptr = ptr (:,:,1:itr,:)
        endif
      end do
     end do
     deallocate (ptr)

     !--------
     ! cleanup
     !--------
     if (dace% lpio .and. k==1) then
       write(6,*) ' cleanup'
       write(6,*)
     endif
     do i = 1, size (table)
       if (table(i)% found) deallocate (table(i)% ptr)
     end do

    end do ! ensemble member loop over k

99  continue
    !--------------
    ! final cleanup
    !--------------
    do k = 1, nens
       if (associated (files(k)% invt)) deallocate (files(k)% invt)
    end do
    deallocate (files)

  contains

    subroutine same_order (inv, ref, i)
      !---------------------------------------------
      ! Rearrange inventory inv to same order as ref
      !---------------------------------------------
      type(t_inventory), intent(inout) :: inv(:) ! inventory to sort
      type(t_inventory), intent(in)    :: ref(:) ! reference inventory
      integer,           intent(out)   :: i      ! Unmatched entry in reference
      !----------------
      ! Local variables
      !----------------
      integer           :: j, k, n         ! Indices
      type(t_inventory) :: new(size (inv))

      i = 0
      n = size (inv)
      if (n /= size (ref)) return

l1:   do j = 1, n
         do k = 1, n
            if ( inv(k)% pa == ref(j)% pa .and. &
                 inv(k)% ti == ref(j)% ti .and. &
                 inv(k)% lv == ref(j)% lv       ) then
               new(j) = inv(k)
               cycle l1
            end if
         end do
         !if(dace% lpio) print*,"### same_order: unmatched position", j
         i = j       ! Return unmatched position in reference
         return
      end do l1
      inv(:) = new(:)

    end subroutine same_order

  end subroutine grib_trange_ens

!===============================================================================

  subroutine apply_trange (spot, obs, y,  vmax_10m, tmin_2m, tmax_2m,   &
                                          tot_prec, aswdir_s, aswdifd_s )
  !-----------------------------------------------------------------
  ! apply observation operator for variables valid for a time range.
  ! derived for COSMO verification with MEC
  ! to be used with ICON as well
  !-----------------------------------------------------------------
  type(t_spot)      ,intent(in)    :: spot         ! report header meta data
  type(t_obs_block) ,intent(in)    :: obs          ! observation data
  real(wp)          ,intent(inout) :: y        (:) ! model equivalent
  real(wp)          ,intent(in)    :: vmax_10m (:) ! wind gust
  real(wp)          ,intent(in)    :: tmin_2m  (:) ! minimum temperature (K)
  real(wp)          ,intent(in)    :: tmax_2m  (:) ! maximum Temperature (K)
  real(wp)          ,intent(in)    :: tot_prec (:) !
  real(wp)          ,intent(in)    :: aswdir_s (:) !
  real(wp)          ,intent(in)    :: aswdifd_s(:) !

    integer  :: i, j ! observation index
    integer  :: k    ! time range index
    integer  :: im   ! time range (minutes)

    !-----------------------------------------------------------------------------
    !   height correction for 2m temperature
    !  vip (3)  =  t_2m  -  zlapse *(hsurf - zstalt)
    ! (alternative lapse rates : constant : -.0065 (e.g for height diff > 500m, or
    !  (Damrath): noon: -.0110 ; midnight : -.0045  --> with temporal interpolat.)
    !      ( for  Tmax: -.0130 ; for Tmin : -.0020 )
    !-----------------------------------------------------------------------------
    real(wp) ,parameter :: zlapse      = -.0065_wp
    real(wp) ,parameter :: zlapse_tmin = -.0020_wp
    real(wp) ,parameter :: zlapse_tmax = -.0130_wp
    real(wp)            :: dz
    dz = 0._wp
    if (spot% z > -900._wp .and. spot% gp_bg /= rvind) &
     dz = spot% z - (spot% gp_bg / gacc)

    do i = 1, spot% o% n
      j = spot% o% i + i
      if (obs% o% body(j)% lev_typ == VN_TRTR) then
        im = nint (obs% o% olev(j) * 60._wp)
        do k = 1, nt
          if (tr(k) == im) then
            select case (obs% o% varno (j))
            case (VN_TMAX)
              if (tmax_2m  (k) >= 0._wp) y (i) = tmax_2m  (k) + dz * zlapse_tmax
            case (VN_TMIN)
              if (tmin_2m  (k) >= 0._wp) y (i) = tmin_2m  (k) + dz * zlapse_tmin
            case (VN_GUST)
              if (vmax_10m (k) >= 0._wp) y (i) = vmax_10m (k)
            case (VN_RR)
              if (tot_prec (k) >= 0._wp) y (i) = tot_prec (k)
            case (VN_RAD_DF)
              if (aswdifd_s(k) >= 0._wp) y (i) = im * 60_wp *    aswdifd_s(k)
            case (VN_RAD_GL)                   ! im * 60 : average->accum
              if (aswdir_s (k) >= 0._wp .and.                                &
                  aswdifd_s(k) >= 0._wp) y (i) =                             &
                                    im * 60_wp * (aswdir_s (k) + aswdifd_s(k))
            end select
            exit
          endif
        end do
      endif
    end do
  end subroutine apply_trange
!==============================================================================
!--------------------------------------------------------
! define MPI allgatherv routine for derived type t_trange
!--------------------------------------------------------
#define VECTOR
#define DERIVED type(t_trange)
#define p_allgather_DERIVED allgatherv_tr
#undef MPI_TYPE
#include "p_allgather.incf"
!==============================================================================
end module mo_obs_trange
