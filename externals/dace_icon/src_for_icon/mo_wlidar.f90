!
!+ ADM-Aeolus observation operator (wind lidar)
!
MODULE mo_wlidar
!
! Description:
!   ADM-Aeolus observation operator (wind lidar)
!
! Current Code Owner: DWD, Alexander Cress
!    phone: +49 69 8062 2716
!    fax:   +49 69 8062 3721
!    email: alexander.cress@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Andreas Rhodin  DWD  2018  original code
!==============================================================================

  !=============
  ! modules used
  !=============
  !-------------------------
  ! general purpose routines
  !-------------------------
  use mo_exception,  only: finish           ! abort routine
  use mo_kind,       only: wp, sp, i8       ! kind parameters
  use mo_mpi_dace,   only: dace             ! MPI group info
  use mo_time,       only: cyyyymmddhhmmss  ! convert time to string
  use mo_physics,    only: d2r,            &! pi /180
                           gacc             ! gravity acceleration
  !----------------------------
  ! acess observation data type
  !----------------------------
  use mo_t_datum,    only: t_datum          ! report body derived type
  use mo_t_obs,      only: t_obs,          &! observation derived type
                           t_spot,         &! observation report derived type
                           t_head,         &! observation header derived type
                           source,         &! list of report source files
                           set_xuv,        &! set unit vectors, zenith angle
                           new_spot,       &! reserve memory
                           new_obs,        &! reserve memory
                           new_int,        &! reserve memory for interp.space
                           set_int_insitu, &! set interpolation space
                           set_vqc_insitu, &! subroutine to set VQC bounds
                           shrink_report,  &! remove passive observations
                           ITY_ICOL,       &! horizontally interpolated input
                           CHR_LIN,        &! H is linear
                           TSK_INIT,       &! task: initialisation
                           TSK_READ,       &!       read observations
                           TSK_SETUP_COLS, &!       specify required columns
                           TSK_SET_CHR,    &!       set characteristics
                           TSK_SETUP_FUL0, &!       setup size of
                           TSK_SETUP_FULL, &!       description of PSAS-space
                           TSK_K,          &!       set up linear operator
                           TSK_Y,          &!       nonlinear operator
                           TSK_R,          &!       observational errors
                           TSK_SHRINK,     &!       release unused observations
                           z2p_hlos         !       derive pressure from height ?
  use mo_wlidar_bc,  only: read_wlidar_bc_nml
  use mo_obs_set,    only: t_obs_block      ! observation derived type
  use mo_obs_tables, only: rept_use,       &! report type usage table
                           decr_rpt_use,   &! degrade status of report
                           check_report_1, &! basic checks on reports
                           check_report_0   ! basic checks on reports
  use mo_fdbk_tables,only: OT_WLIDAR,      &! wind lidar obstape code
                           VN_HLOS,        &! HorizontalLineOfSight wind code
                           VN_P,           &! pressure (level) code
                           VN_HEIGHT        ! height   (level) code
  use mo_t_use,      only: t_use,          &! data type to hold state
                           use_0,          &! default values of type use
                           STAT_ACTIVE,    &! flag for active    observation
                           STAT_DISMISS,   &! flag for dismissed observation
                           CHK_NONE,       &! flag for no check
                           CHK_INSDAT,     &! flag for insufficient data
                           CHK_DOMAIN       ! flag for out of domain
  use mo_obs_err,    only: obs_err          ! get observation error from table
  !----------------------------
  ! acess atmospheric data type
  !----------------------------
  use mo_atm_state,  only: t_atm            ! atm. state data type
  use mo_t_col,      only: t_cols,         &! model columns data type
                           COL_UV,         &! specification of fields
                           COL_GEO          !
  !------------------------------------
  ! Interface to interpolation routines
  !------------------------------------
  use mo_grid_intpol,only: idx_init         ! get grid-indices & interp.coeff.
  !-----------------------------------
  ! acess vector and matrix data types
  !-----------------------------------
  use mo_dec_matrix, only: t_vector_segm    ! vector segment data type
  !--------------------------------------
  ! obs data header info read from NetCDF
  !--------------------------------------
  use mo_head_netcdf,only: ncid,          &! NetCDF file Id
!                          imissing,      &! NetCDF _FillValue for integer
                           rmissing,      &! NetCDF _FillValue for reals
                           s1cent,        &! data centre
                           stime,         &! header observation time (section1)
                           db_time,       &! data bank time
                           s1cents,       &! data sub centre
!                          mlah,          &! latitude
!                          mloh,          &! longitude
                           istidn,        &! here: satellite identifier
                           obs_time,      &! body observation time
                           get_int,       &! read integer from NetCDF file
                           get_real        ! read real    from NetCDF file
  use mo_usstd,     only:  p_h_usstd       ! pressure from geopotential height

  implicit none

!------------------------------------------------------------------------------
  !================
  ! public entities
  !================
  private
  public :: process_wlidar     ! general purpose WLIDAR processing routine
  public :: read_wlidar_netcdf ! read WLIDAR observation from netCDF file

!------------------------------------------------------------------------------
  !=================
  ! module variables
  !=================

!==============================================================================
contains
!==============================================================================

  subroutine process_wlidar (task, spot, obs, atm, cols, xi, y, Jo, Jo_atm, &
                            state)
  integer            ,intent(in)             :: task    ! what to do
  type(t_spot)       ,intent(inout),optional :: spot    ! SPOT observations
  type(t_obs_block)  ,intent(inout),optional :: obs     ! observation data type
  type(t_atm)        ,intent(in)             :: atm     ! atmospheric state
  type(t_cols)       ,intent(in)   ,optional :: cols    ! model columns
  type(t_vector_segm),intent(in)   ,optional :: xi      ! interpolated values
  type(t_vector_segm),intent(inout),optional :: y       ! observed quantity
  real(wp)           ,intent(inout),optional :: Jo      ! obs. cost funct. Jo
  type(t_atm)        ,intent(inout),optional :: Jo_atm  ! gradient:d Jo/d atm
  integer            ,intent(in)   ,optional :: state   ! status flag

    !================
    ! local variables
    !================
    integer  :: tsk     ! task (local copy)
    integer  :: i, n    ! observation index and len
    integer  :: ii      ! input index
    integer  :: j, k, l ! indices
    integer(i8) :: icol    ! columns to retrieve
    integer  :: ierr    ! error return parameter
    logical  :: change  ! report has changed ?
!   integer       ,pointer :: ty (:)        ! unpacked levels
!   real(sp)      ,pointer :: acc (:)       ! unpacked accuracy
    real(wp) :: H(2)    ! temporary Jacoby matrix
    real(wp) :: a       ! azimuth
    real(wp) :: s, c    ! sin, cos of azimuth
    real(wp) :: u, v    ! wind components
    real(wp) :: w       ! model equivalent to observation
    real(wp) :: z       ! height(gpm) for p.-interp.
    real(wp) :: wzp     ! weight for p.-interpolation
    real(wp) :: err     ! Observation error estimate
!   real(sp) :: plev    ! pressure from height
!   Test auf wp gesetzt
    real(wp) :: plev    ! pressure from height

    !======================
    ! executable statements
    !======================
    tsk = task
    if(tsk==0) return

    !==================================================
    ! skip tasks not required for this observation type
    !==================================================
    tsk = task
    tsk = iand (task, not (&
          TSK_READ))        ! Input is read by separate routine
    if (tsk == 0) return

    !===============
    ! initialisation
    !===============
    if (iand (TSK_INIT,tsk) /= 0) then
      tsk = tsk - TSK_INIT
      call read_wlidar_bc_nml
      if (tsk == 0) return
    endif

    !=============================================
    ! TSK_SET_CHR: set observation characteristics
    !=============================================
    if (iand (TSK_SET_CHR,tsk) /= 0) then
      spot% int_type  = ITY_ICOL   ! expect horizontally interpolated input
      spot% cost      = 1._wp      ! operator is cheap
      spot% nr        = spot% o% n ! diagonal R
      spot% char = CHR_LIN         ! this is a nonlinear operator
      tsk = tsk - TSK_SET_CHR
      if (tsk == 0) return
    endif

    !=========================================
    ! TSK_SETUP_COLS: specify required columns
    !=========================================
    if (iand (TSK_SETUP_COLS,tsk) /= 0) then
      i = spot% o%i  +1
      n = spot% o%n+i-1
      icol = COL_UV
      if (obs% o% body (spot%o%i+1)% plev < 0._wp) &
        icol = icol + COL_GEO
      call idx_init (      &
            spot% col% c,  &! <-  column descriptor
            spot% col% h,  &!  -> interpolation coefficients
            obs% o% mc,    &! <-> model column descriptors
            icol,          &! <-  fields required
            0,             &! <-  tracers required
            atm% grid,     &! <-  model grid
            spot% i_time,  &! <-  time slot
            spot% w_time   )! <-  time interpolation weight
      if (spot% col% h% imc(1,1) == 0) call decr_rpt_use (spot, CHK_DOMAIN)
      tsk = tsk - TSK_SETUP_COLS
      if (tsk == 0) return
    endif

    !===========================
    ! tsk == TSK_SETUP_FUL0:
    ! set up interpolation space
    ! part 1: reserve memory
    !===========================
    if (iand (TSK_SETUP_FUL0,tsk) /= 0) then
      call new_int (obs% o, spot, spot% o% n * 2)
      tsk = tsk - TSK_SETUP_FUL0
      if (tsk == 0) return
    endif

    !==============================================================
    ! tsk == TSK_SETUP_FULL:
    ! part 2: interpolate pressure from height using the bg-profile
    ! setup description of interpolation-space
    !==============================================================
    if (iand (TSK_SETUP_FULL,tsk) /= 0) then
      if (dace% pe == obs% o% pe .and. &
          z2p_hlos == 2                ) then
        j = sum (maxloc (spot% col% h% w))
        k = spot% col% h% imc (j,1)
!NEC$ nomove
        do i = spot%o%i+1, spot%o%i+spot%o%n
          !-----------------------------------------------------------
          ! interpolate pressure if height is the independent quantity
          !-----------------------------------------------------------
          if (obs% o% body (i)% plev < 0._wp) then
            n = size(cols% col(k)% geo)
            plev    = -1._wp
            z       = obs% o% olev (i) * gacc
            if (z <= cols% col(k)% geo(n)) then
              plev = exp (cols% col(k)% p(n))
            elseif (z >= cols% col(k)% geo(1)) then
              plev = exp (cols% col(k)% p(1))
            else
              do l = 2, size(cols% col(k)% geo)
                if (z > cols% col(k)% geo(l)) then
                  wzp     = (z                      - cols% col(k)% geo(l)) / &
                            (cols% col(k)% geo(l-1) - cols% col(k)% geo(l))
                  plev = exp (      wzp  * cols% col(k)% p(l-1) + &
                             (1._wp-wzp) * cols% col(k)% p(l  )   )
                  exit
                endif
              end do
            endif
!           print*,'Interpolate pressure: ', obs% o% olev (i), z, plev
            obs% o% body (i)% plev    = plev
            spot% ps                  = plev
          endif
        end do
      endif
      !---------------------------
      ! set up interpolation space
      !---------------------------
      call set_int_insitu (spot, obs% o)

      tsk = tsk - TSK_SETUP_FULL
      if (tsk == 0) return
    endif

    !===========================================
    ! tsk == TSK_SETUP_FULL:
    ! setup description of PSAS-space
    ! observed values were set up while reading.
    !===========================================
    if (iand (TSK_SETUP_FULL,tsk) /= 0) then
      tsk = tsk - TSK_SETUP_FULL
      if (tsk == 0) return
    endif

    !=========================================
    ! run operator or set up H (Jacoby-matrix)
    !=========================================
    if (iand (TSK_K+TSK_Y,tsk) /= 0) then
      if (spot% pe_eval == dace% pe) then
        !-------------------------------
        ! get indices, check consistence
        !-------------------------------
        i = spot%o% i + 1
        n = spot%o% n
        if (spot%o% n /= 1) call finish ('process_wlidar','o%n/=1')
        if (spot%i% n /= 2) call finish ('process_wlidar','i%n/=2')
        if (obs% o% varno(i) /= VN_HLOS) call finish ('process_wlidar','HLOS')
        ii = spot%i% i + 1
        !-------------------------------
        ! evaluate operator and Jacobian
        !-------------------------------
        a  = obs% o% body(i)% obs_par(1)
        c  = cos (d2r * a)
        s  = sin (d2r * a)
        u  = xi% x(ii  )
        v  = xi% x(ii+1)
        w  = - u * s - v * c
        H  = [-s, -c]
        !--------------------------------
        ! store results, forward operator
        !--------------------------------
        if (iand (TSK_Y,tsk) /= 0) then
          if (present(y)) y%x (i) = w
        endif
        !------------------------
        ! store results, Jacobian
        !------------------------
        if (iand (TSK_K,tsk) /= 0) then
          if (present(y)) y%x (i) = w
          k = obs% H% ia (ii)
          do j=1,spot% i% n                  ! columns
            obs% H% ia (spot% i% i +j) = k   ! column index
            l = (j - 1) * n + 1
            obs% H% packed (k) = H(l)        ! coefficient
            obs% H% ja     (k) = i           ! row index
            k = k + 1
          end do
          obs% H% ia (ii + spot% i% n) = k   ! column index
          obs% yi% x (i)               = w   ! result
        endif
      endif
      if (iand (TSK_Y,tsk) /= 0) tsk = tsk - TSK_Y
      if (iand (TSK_K,tsk) /= 0) tsk = tsk - TSK_K
      if (tsk == 0) return
    endif

    !===============================================
    ! tsk == TSK_R
    ! set up R (observation error covariance matrix)
    ! specify observation errors via namelist
    !===============================================
     if (iand (TSK_R,tsk) /= 0) then
      if (dace% pe == obs% o% pe) then
       n = spot% o% n
      !================================
      ! setup HLOS observational errors
      !================================
       k  = obs% R% ia (spot% o% i+1)
       ii = spot% o% i+1
      obs% R% ia (ii) = k
      select case (obs% o% varno (spot%o% i+1))
      case (VN_HLOS)
        obs% R% packed(k) = obs_err (OT_WLIDAR,                   &
                                     spot%hd% codetype,           &
                                     VN_HLOS,                     &
                                     dble(obs%o% body (spot%o% i+1)% plev),   &
                                     obs%o% body (spot%o% i+1)% o,&
                              ierr = ierr                         ) **2
        if (ierr /= 0) obs% R% packed(k) = 1._wp
      case default
        call finish('process_wlidar (TSK_R)','invalid observation type')
      end select
      obs% R% ja (k) = ii
      k = k + 1
      obs% R% ia (spot% o% i + n + 1) = k
      !=================
      ! setup vqc bounds
      !=================
      call set_vqc_insitu (spot, obs% o)
     endif
     tsk = tsk - TSK_R
     if (tsk == 0) return
    endif

    !==================================================================
    ! tsk == TSK_SHRINK:
    ! release unused observations in the report
    !==================================================================
    if (iand (TSK_SHRINK,tsk) /= 0) then
      call shrink_report (spot, obs%o, state, change)
      if (change) call set_int_insitu (spot, obs% o)
      tsk = tsk - TSK_SHRINK
      if (tsk == 0) return
    endif

    !===========
    ! left tasks
    !===========
    if (tsk /= 0) then
      if (dace% lpio) write (6,*) 'process_wlidar: unknown task',tsk
      call finish ('process_wlidar','unknown task')
    endif

  end subroutine process_wlidar

!==============================================================================

  subroutine read_wlidar_netcdf (ifile, i_source, obs, head, lkeep, nkeep)
  !================================================================
  ! Read HLOS WLIDAR observations from netCDF (converted BUFR data)
  !=================================================================
  integer       ,intent(in)           :: ifile    ! # of netCDF file read
  integer       ,intent(inout)        :: i_source ! # of records in source-file
  type (t_obs)  ,intent(inout)        :: obs      ! observations data type
  type (t_head) ,intent(in)           :: head     ! header data already encoded
  logical       ,intent(out)          :: lkeep    ! accept observation ?
  integer       ,intent(out)          :: nkeep    ! number of accepted obsvs.

    !================
    ! local variables
    !================
    type (t_use)          :: use        ! status variable
    type (t_head)         :: hd         ! report header
    type (t_spot)         :: spt0, spti ! report meta data
    type (t_spot) ,save   :: empty      !
    integer               :: len_report ! number of reports from file
    integer               :: n          ! number of reports to process
    integer               :: is, isum   ! report index variable
    integer               :: entry      ! total report count
    logical               :: lk         ! flag to keep this report
    real(sp) ,allocatable :: nhlsw (:)  ! Horizontal Line Of Sight wind
    real(sp) ,allocatable :: nhoil (:)  ! Horizontal observation integration length
    real(sp) ,allocatable :: neshl (:)  ! error estimate of HLOS
    real(sp) ,allocatable :: mppp  (:)  ! pressure
    real(sp) ,allocatable :: mda1  (:)  ! bearing or azimuth
    real(sp) ,allocatable :: mde1  (:)  ! elevation
    real(sp) ,allocatable :: mhhir1(:)  ! height at centre-of-gravitiy
    real(sp) ,allocatable :: mlah1 (:)  ! latitude at centre-of-gravity
    real(sp) ,allocatable :: mloh1 (:)  ! longitude at centre-of-gravity
    integer  ,allocatable :: mcofl (:)  ! confidence flag
    integer  ,allocatable :: mlclt (:)  ! Lidar L2b classification type
    integer  ,allocatable :: mrech (:)  ! satellite reciever

    lkeep = .false.
    nkeep = 0
    isum  = 0

    !------------------------
    ! get dimension of fields
    !------------------------
    if (.not. allocated (istidn)) call finish ('read_wlidar_netcdf','istidn')
    len_report = size (istidn)
    n          = min (len_report, rept_use(OT_WLIDAR)% max_proc)

!   write(6,*) ' BIN im Einlesen der Wind Lidar Daten '

    !----------------
    ! allocate arrays
    !----------------
    allocate (mhhir1 (len_report))
    allocate (mlah1  (len_report))
    allocate (mloh1  (len_report))
    allocate (nhlsw  (len_report))
    allocate (nhoil  (len_report))
    allocate (neshl  (len_report))
    allocate (mppp   (len_report))
    allocate (mda1   (len_report))
    allocate (mde1   (len_report))
    allocate (mcofl  (len_report))
    allocate (mlclt  (len_report))
    allocate (mrech  (len_report))

    !--------------------------------
    ! read variables from NetCDF file
    !--------------------------------
    call get_real (mhhir1, 'MHHIR1',-999._sp)
    call get_real (mlah1,  'MLAH1', -999._sp)
    call get_real (mloh1,  'MLOH1', -999._sp)
    call get_real (nhlsw,  'NHLSW', -999._sp)
    call get_real (nhoil,  'NHOIL', -999._sp)
    call get_real (neshl,  'NESHL', -999._sp)
    call get_real (mppp ,  'MPPP' , -999._sp)
    call get_real (mda1 ,  'MDA1' , -999._sp)
    call get_real (mde1 ,  'MDE1' , -999._sp)
    call get_int  (mrech,  'MRECH', -1      )
    call get_int  (mcofl,  'MCOFL', -1      )
    call get_int  (mlclt,  'MLCLT', -1      )

    !-------------------------------
    ! preset total number of reports
    !-------------------------------
    entry = sum (source(1:ifile-1)% entries)

    !------------------
    ! loop over reports
    !------------------
    do is = 1, n

      entry   = entry    + 1

      !-----------------------
      ! initialize usage flags
      !-----------------------
      use = use_0

      !--------------------
      ! define head section
      !--------------------
      hd            = head
      hd% time      = stime   (is)
      hd% db_time   = db_time (is)
      hd% source    = ifile
      hd% record    = is
      hd% id        = entry
      hd% center    = s1cent  (is)
      hd% subcenter = s1cents (is)

!     if (is < 100) then
      if (.FALSE. ) then
        write (6,'()')
        write (6,'( 8(a16, i8,/),   &
                  & 2(a16, a ,/),   &
                  & 6(a16, i8,/) )' )                      &
          'pe='         ,dace% pe,                         &
          'head is='    ,is,                               &
          'obstype='    , hd% obstype  ,                   &
          'dbkz='       , hd% dbkz     ,                   &
          'modtype='    , hd% modtype  ,                   &
          'buf_type='   , hd% buf_type ,                   &
          'buf_subtype=', hd% buf_subtype,                 &
          'codetype='   , hd% codetype ,                   &
          'time='       , cyyyymmddhhmmss (hd% time)   ,   &
          'db_time='    , cyyyymmddhhmmss (hd% db_time),   &
          'dbk='        , hd% idbk     ,                   &
          'source='     , hd% source   ,                   &
          'record='     , hd% record   ,                   &
          'id='         , hd% id       ,                   &
          'center='     , hd% center   ,                   &
          'subcenter='  , hd% subcenter
      endif

      !--------------------------------------------
      ! perform simple generic check on report type
      !--------------------------------------------
      call check_report_0 (use, hd, 1)
      if (use% state <= STAT_DISMISS) cycle

      !------------------
      ! create new report
      !------------------
      spt0             = empty
      spt0% use        = use
      spt0% hd         = hd
      spti             = spt0

      !--------------------------
      ! check for sufficient data
      !--------------------------
      if (nhlsw (is) == -999._sp .or. &
          neshl (is) == -999._sp .or. &
          mda1  (is) == -999._sp .or. &
          mde1  (is) == -999._sp .or. &
          mhhir1(is) == -999._sp .or. &
          mrech (is) == -1       .or. &
          mcofl (is) == -1       .or. &
          mlah1 (is) == rmissing .or. &
          mloh1 (is) == rmissing      ) then
        call decr_rpt_use (spti, CHK_INSDAT, comment='read_wlidar_netcdf')
        cycle
      endif

      if (mppp  (is) == -999._sp .and. z2p_hlos <= 0) then
        call decr_rpt_use (spti, CHK_INSDAT, comment='pressure not reported in BUFR')
        cycle
      endif

      if (mcofl (is) == 1 ) then
         call decr_rpt_use (spti, CHK_INSDAT, comment='Confidence flag is set to invalid')
        cycle
      endif

      if (mhhir1(is) < 0 ) then
         call decr_rpt_use (spti, CHK_INSDAT, comment='Height below zero')
         cycle
      endif
!
! The HLOS geometric heights were systematically 250 m too low due to a known error in the
! LOS pointing knowledge determined during the CP (star-tracker calibration issue) and
! therefore all heihgts are increased 250 m. Valid till February 2019.

!      mhhir1(is) = mhhir1(is) + 250.0_sp

      if (mrech(is) == 1 .AND. mlclt(is) == 1) then
         call decr_rpt_use (spti, CHK_INSDAT, comment='Rayleigh cloudy winds')
         cycle
      endif

      if (mrech(is) == 0 .AND. mlclt(is) == 0) then
         call decr_rpt_use (spti, CHK_INSDAT, comment='Mie clear winds')
         cycle
      endif

      if (mhhir1(is) <  1000.0_sp) then
         call decr_rpt_use (spti, CHK_INSDAT, comment='HLOS winds unter 1000 m')
         cycle
      endif

      if (mrech(is) == 1 .AND. mlclt(is) == 0 .AND. nhoil(is) < 60._sp ) then
         call decr_rpt_use (spti, CHK_INSDAT, comment='Rayleigh clear hlos winds with NHOIL less than 60')
         cycle
      endif

      if (mrech(is) == 1 .AND. mlclt(is) == 0 .AND. neshl(is) > 6._sp ) then
         call decr_rpt_use (spti, CHK_INSDAT, comment='Rayleigh clear hlos error estimate threshold')
         cycle
      endif

      if (mrech(is) == 0 .AND. mlclt(is) == 1 .AND. nhoil(is) < 5._sp ) then
         call decr_rpt_use (spti, CHK_INSDAT, comment='Mie cloudy hlos winds with NHOIL less than 5')
         cycle
      endif

      if (mrech(is) == 0 .AND. mlclt(is) == 1 .AND. neshl(is) > 4._sp ) then
         call decr_rpt_use (spti, CHK_INSDAT, comment='Mie cloudy hlos error estimate threshold')
         cycle
      endif
      if (mrech(is) == 0 .AND. mhhir1(is) > 15000._sp ) then
         call decr_rpt_use (spti, CHK_INSDAT, comment='Mie cloudy hlos winds above 15000 m')
         cycle
      endif

!     PRINT*,'Reviever choice: ', mrech(is), mlclt(is)
      isum = isum + 1

      !------------------------------
      ! process report header entries
      !------------------------------
      spti% col% c% dlat = mlah1    (is)
      spti% col% c% dlon = mloh1    (is)
      spti% actual_time  = obs_time (is)
      spti% ident        = istidn   (is)
      spti% stret        = mrech    (is)
      spti% col% nlev    = 1
      call set_xuv (spti)
      !------------------------
      ! set center / processing
      !------------------------
      spti% statid = 'ADM'

!     if (isum < 100) then
      if (.FALSE. ) then
        write (6,'()')
        write (6,'(   a20, i6  ,a, /, &
                  &   a20, i6  ,   /, &
                  & 2(a20,f8.3 ,   /),&
                  &   a20, a   ,   / ,&
                  &   a20, a   ,   / ,&
                  & 6(a20, i5      /))' )                    &
              'pe=',dace% pe,'  spti ',                      &
              'spti% corme        = ', spti% corme        ,  &
              'spti% col% c% dlat = ', spti% col% c% dlat ,  &
              'spti% col% c% dlon = ', spti% col% c% dlon ,  &
              'spti% actual_time  = ', cyyyymmddhhmmss (spti% actual_time),  &
              'spti% statid       = ', spti% statid       ,  &
              'spti% ident        = ', spti% ident        ,  &
              'spti% stret        = ', spti% stret        ,  &
              'spti% col% nlev    = ', spti% col% nlev
      endif

      !------------------------------------------------------------------
      ! derive pressure level from height
      ! z2p_hlos = -1 ! for testing: use pressure as nominal coordinate
      ! z2p_hlos =  0 ! do not derive pressure level, take from BUFR
      ! z2p_hlos =  1 ! derive pressure level from US standard atmosphere
      ! z2p_hlos =  2 ! as 1) later re-calculated from background profile
      !------------------------------------------------------------------
      select case (z2p_hlos)
      case (1)
        mppp(is)  = p_h_usstd (real(mhhir1(is),wp))
      case (2)
        mppp(is)  = -999._wp  ! postpone until bg profile is present
      end select

      !----------------
      ! standard checks
      !----------------
      call check_report_1 (spti)
      lk = spti% use% state > STAT_DISMISS
      if (lk) then
        call check_store_wlidar (spti, obs, lk,                           &
             nhlsw(is), neshl(is), mppp(is), mhhir1(is), mda1(is), mde1(is), mcofl(is))
        if (lk) nkeep = nkeep + 1
      endif
    end do

    write(6,*) 'read_wlidar_netcdf: keep', nkeep, 'observations out of', n

  end subroutine read_wlidar_netcdf

!==============================================================================

  subroutine check_store_wlidar (spot, obs, lkeep,                            &
                                 nhlsw, neshl, mppp, mhhir1, mda1, mde1, mcofl)
  type(t_spot)    ,intent(inout)  :: spot    ! report derived type
  type(t_obs)     ,intent(inout)  :: obs     ! observations derived type
  logical         ,intent(out)    :: lkeep   ! keep or reject observation
  real(sp)        ,intent(in)     :: nhlsw   ! Horizontal Line Of Sight wind
  real(sp)        ,intent(in)     :: neshl   ! error estimate of HLOS
  real(sp)        ,intent(in)     :: mppp    ! pressure
  real(sp)        ,intent(in)     :: mhhir1  ! height
  real(sp)        ,intent(in)     :: mda1    ! bearing or azimuth
  real(sp)        ,intent(in)     :: mde1    ! elevation
  integer         ,intent(in)     :: mcofl   ! confidence flag

    type(t_spot),pointer :: spt          ! temporary
    type(t_datum)        :: bod          ! body derived type
    integer              :: id           ! observation id
    integer              :: i1, in       ! index range
    integer              :: no           ! number of observations in the report

    lkeep = .true.

!   !----------------------
!   ! evaluate quality flag
!   !----------------------
!   call decr_rpt_use (spot, CHK_NOTUSED, STAT_PASSIVE, &
!                                         comment='....')
!   !------------------------------
!   ! exit if no valid data present
!   !------------------------------
!   if (.........) then
!     lkeep = .false.
!     call decr_rpt_use (spot, CHK_INSDAT, comment='....')
!     return
!   endif
!   !--------------------------
!   ! report selection (filter)
!   !--------------------------
!  if (...) call decr_rpt_use (spot, CHK_BLACKLIST)

    if (spot% use% state > STAT_DISMISS) then
      !--------------------------------------------------
      ! new report header entry in DACE  observation list
      !--------------------------------------------------
      call new_spot (obs,1, set_id=.true.)
      spt => obs% spot (obs% n_spot)
      id      = spt% id
      spt     = spot
      spt% id = id
      no      = 1
      !------------------
      ! fill in body info
      !------------------
      bod % use % state = STAT_ACTIVE
      bod % use % check = CHK_NONE
      bod % mn          = 'hlos'
      call new_obs (obs, no, spot=spt)
      i1 = spt% o% i+1
      in = spt% o% i + spt% o% n
      obs % varno (i1 : in)             = VN_HLOS   ! HLOS code
      obs %  olev (i1 : in)             = mhhir1    ! height
      obs %  body (i1 : in)             = bod
      obs %  body (i1 : in)% o          = nhlsw     ! HLOS wind
      obs %  body (i1 : in)% plev       = mppp      ! pressure
      obs %  body (i1 : in)% eo         = 1._sp     ! observation error
      obs %  body (i1 : in)% ac         = neshl     ! error estimate
      obs %  body (i1 : in)% obs_par(1) = mda1      ! azimuth
      obs %  body (i1 : in)% obs_par(2) = mde1      ! elevation
!     obs %  body (i1 : in)% lev_sig    =           ! instrument id ?
      obs %  body (i1 : in)% lev_typ    = VN_HEIGHT ! height code
!     obs %  body (i1 : in)% pcc        =           ! quality measure
      if (z2p_hlos == -1) then
        obs %  olev (i1 : in)           = mppp      ! pressure
        obs %  body (i1 : in)% lev_typ  = VN_P      ! pressure code
      endif
    else
      lkeep = .false.
    endif

  end subroutine check_store_wlidar

!==============================================================================
end module mo_wlidar
