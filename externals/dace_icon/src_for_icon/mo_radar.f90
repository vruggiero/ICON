!
!+ Radar observation handling in DACE
!
MODULE mo_radar
!
! Description:
!   Radar observation handling in LETKF.
!   The link to the radar operator (subroutine process_radar) is
!   located in module 'mo_radar_obs'.
!
! Current Code Owner: DWD, Elisabeth Bauernschubert
!    phone: +49 69 8062 2715
!    fax:   +49 69 8062 3721
!    email: elisabeth.bauernschubert@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_19        2012-04-16 Andreas Rhodin
!  module for radar observation handling
! V1_20        2012-06-18 Andreas Rhodin
!  cleanup: remove unused variables
! V1_22        2013-02-13 Andreas Rhodin
!  dealiasing of radar radial winds; write observations back to feedback files
! V1_26        2013/06/27 Andreas Rhodin
!  LETKF RADAR operator: correctly check consistency of fof input files
! V1_28        2014/02/26 Andreas Rhodin
!  comments changed
! V1_29        2014/04/02 Andreas Rhodin
!  option to constrain radar reflectivities to minimum value
!  use VN_RREFL instead of VN_REFL for LETKF radar reflectivity observation
! V1_31        2014-08-21 Hendrik Reich
!  bug fix in constrain_refl
! V1_47        2016-06-06 Andreas Rhodin
!  move routines split_radar, join_radar to new module mo_split_obs, rename
! V1_51        2017-02-24 Andreas Rhodin
!  subroutine thin_radar: be more verbose, remove invalid radial winds
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Andreas Rhodin            DWD  2012  original source
! Elisabeth Bauernschubert  DWD  2018
!==============================================================================

  !=============
  ! modules used
  !=============
  !-------------------------
  ! general purpose routines
  !-------------------------
  use mo_exception,   only: finish,      &! abort routine
                            message       ! write warning
  use mo_kind,        only: wp, sp        ! working precision kind parameter
  use mo_mpi_dace,    only: dace,        &! MPI group info
                            p_bcast       ! broadcast routine
  use mo_namelist,    only: position_nml,&! position namelist
                            nnml,        &! namelist Fortran unit number
                            POSITIONED    ! ok    code from position_nml
! use mo_time,        only: chhmm         ! derive string from time
  use mo_physics,     only: d2r,         &! conversion: degree  -> radians
                            r2d,         &! conversion: degree <-  radians
                            rearth,      &! earth radius
                            pi            ! 3.1415....
  use mo_t_netcdf,    only: stanc,       &! NetCDF start  parameter
                            counc,       &! NetCDF count  parameter
                            strnc,       &! NetCDF stride parameter
                            get_var       ! read variable
  use mo_dec_matrix,  only: t_vector,    &! vector derived type
                            t_bvector,   &! boolean vector derived type
                            t_vector_segm ! vector segment data type
  !-----------------------------
  ! access atmospheric data type
  !-----------------------------
  use mo_atm_state,   only: t_atm         ! atm. state data type
  !-----------------------------
  ! access observation data type
  !-----------------------------
  use mo_obs_set,     only: t_obs_block   ! obs data type
  use mo_t_obs,       only: t_obs,       &! observation derived type
                            t_spot,      &! report derived type
                            read_cdfin,  &! flag to read COSMO observations
                            TSK_INIT,    &! FLAGS: initialize module
                            TSK_READ,    &!  read observations
                            new_par,     &! reserve memory for parameters
                            set_xuv       ! set components of t_spot
  use mo_t_datum,     only: t_datum,     &! observation body derived type
                            rvind         ! invalid value
  use mo_t_use,       only: CHK_NONE,    &! no check applied
!                           CHK_THIN,    &! thinning check value
                            CHK_INSDAT,  &! insufficient/inconsistent data
                            CHK_NOTUSED, &! observation not used flag value
                            CHK_OPERATOR,&! Operator not applicable
                            CHK_FG,      &! obs - fg error check
                            STAT_FORGET, &!
                            STAT_DISMISS,&! status    flag value
                            STAT_OBS_ONLY,&!
                            STAT_PAS_REJ ,&! status   flag value
                            STAT_REJECTED,&! status   flag value
                            STAT_ACTIVE, &! status    flag value
                            stat_mnem,   &! get mnemonic of state
                            decr_use      ! decrease the state of an observation
  use mo_obs_tables,  only: rept_use,    &! report type usage table
                            decr_rpt_use  ! change use-flags of report
  use mo_t_col,       only: t_cols        ! model columns data type
  use mo_fdbk_tables, only: init_fdbk_tables,&! initialise the tables
                            OT_RADAR,    &! radar observation type flag value
                            VN_RADVEL,   &! radial velocity        flag value
                            VN_RREFL,    &! radar reflectivity     flag value
                            VN_REFL       !       reflectivity     flag value
  !---------------
  ! radar operator
  !---------------
  use mo_radar_obs,   only: init_radar_obs, &! initialise mo_cosmo
                            read_radar_obs, &! read COSMO data
                            radar2dace       ! fill radar data into t_obs
  implicit none

!------------------------------------------------------------------------------
  !================
  ! public entities
  !================
  private
  public :: process_radar    ! specific RADAR operator processing routine
  public :: t_radar          ! radar obstype specific table entry
  public :: load, store      ! load, store t_radar in t_obs
  public :: thin_radar       ! thinning of radar observations
  public :: read_nml_radar   ! read namelist /RADAR_OBS/
  public :: read_fdbk_radar  ! read radar specific entries from feedback file
  public :: radar_tables     ! prepare radartables for writing to feedback file
  public :: shrink_par_radar ! adapt table 'par' to reduced set of observations
  public :: alias_radar      ! de-alias radial wind observations
  public :: check_radar      ! check model equivalents of radial wind obs.
  public :: constrain_refl   ! constrain reflectivity (obs,fg)

  public :: luse_eo_factor   ! use eo from input feedback file as factor
                             ! and multiply it with eo given by dace namelist

!------------------------------------------------------------------------------
  !=================
  ! module variables
  !=================

  !---------------------
  ! Namelist /RADAR_OBS/
  !---------------------
  integer  :: use_refl   = STAT_ACTIVE  ! radar reflectivity usage flag
  integer  :: use_radvel = STAT_ACTIVE  ! radial velocity    usage flag
  integer  :: iprintout  =   0          ! steering of printout
  logical  :: luse_eo_factor = .false.  ! use eo from input feedback file as factor and
                                        ! multiply it with eo given by dace namelist
  logical  :: dealias_fg = .true.       ! dealias radial wind (by first guess)
  real(wp) :: chk_alias  = 2._wp        ! check dealiasing (compare to spread)
  real(wp) :: ofg_alias  = 0.7_wp       ! check dealiasing (compare to o-fg)
  real(wp) :: min_refl   = -999._wp     ! constrain reflectivity (obs,fg)
  !------------------------------------------------
  ! +++ CURRENTLY not used: min/max_range/xdist +++
  !------------------------------------------------
! real(wp) :: min_range  =   0._wp     ! minimum range                  (km)
! real(wp) :: max_range  = 999._wp     ! maximum range                  (km)
! real(wp) :: min_hdist  =   0._wp     ! minimum hor.  distance between rays
! real(wp) :: min_vdist  =   0._wp     ! minimum vert. distance between rays
! real(wp) :: min_rdist  =   0._wp     ! minimum distance within a ray  (km)

  namelist /RADAR_OBS/ use_refl, use_radvel, iprintout,      &
                       luse_eo_factor, dealias_fg, chk_alias, &
                       ofg_alias, min_refl
!                      min_range, max_range, min_hdist, min_vdist, min_rdist
  !------------------------------------------------------
  ! radar specific observation table entry (for each ray)
  ! to be stored in obs% par
  !------------------------------------------------------
  type t_radar
    real(sp) :: radar_azimuth      ! azimuth of ray         (degree)
    real(sp) :: radar_elevation    ! elevation of ray       (degree)
    integer  :: radar_nrange       ! number of bins per ray
    real(sp) :: radar_range_start  ! distance of first bin       (m)
    real(sp) :: radar_drange       ! distance between bins       (m)
    integer  :: radar_nbody        ! number of corresponding body entries
  end type t_radar

  integer ,save :: ray_int_size  = 0 ! length of data type t_radar
  integer ,save :: ray_byte_size = 0 ! length of data type t_radar

  !-----------
  ! interfaces
  !-----------
  interface load
    module procedure load_radar
  end interface load

  interface store
    module procedure store_radar
  end interface store


!==============================================================================
contains
!==============================================================================
  subroutine process_radar (task, spot, obs, atm, cols, xi, y, Jo, Jo_atm, &
                            state)
  integer            ,intent(in)             :: task    ! what to do
  type(t_spot)       ,intent(inout),optional :: spot    ! SPOT observations
  type(t_obs_block)  ,intent(inout),optional :: obs     ! observation data type
  type(t_atm)        ,intent(in)             :: atm     ! atmospheric state
  type(t_cols)       ,intent(in)   ,optional :: cols    ! model columns
  type(t_vector_segm),intent(in)   ,optional :: xi      ! interpolated values
  type(t_vector_segm),intent(inout),optional :: y       ! model equivalents
  real(wp)           ,intent(inout),optional :: Jo      ! obs. cost funct. Jo
  type(t_atm)        ,intent(inout),optional :: Jo_atm  ! gradient:d Jo/d atm
  integer            ,intent(in)   ,optional :: state   ! status flag
  !-----------------------------------------------------------------------
  ! This subroutine is called from various points in the assimilation code
  ! in order to perform specific tasks for the RADAR operator.
  !------------------------------------------------------------------------

    integer                   :: tsk          ! task, local copy

    !==============================
    ! observation non_specific part (called once)
    !==============================
    if (rept_use(OT_RADAR)% use (CHK_NONE) <= STAT_FORGET) return
    if (rept_use(OT_RADAR)% init           == 0          ) return
    tsk = task

    !------------------------------------------
    ! tsk == TSK_INIT:
    ! module initialisation, namelist read etc.
    !------------------------------------------
    if (iand (TSK_INIT,tsk) /= 0) then
      if (rept_use(OT_RADAR)% init >= 3) call init_radar_obs (atm% grid)
      tsk=tsk-TSK_INIT
    endif
    if (tsk==0) return

    !-----------------------------
    ! tsk == TSK_READ:
    ! read RADAR observation files
    !-----------------------------
    if (iand (TSK_READ,tsk) /= 0) then
      if (read_cdfin) then
        !--------------
        ! read the file
        !--------------
        call init_fdbk_tables
        call read_radar_obs ()
        call radar2dace ()
      endif
      tsk=tsk-TSK_READ
    endif
    if (tsk==0) return

!     !----------------------
!     ! tsk == TSK_Y:
!     ! run operator EMVORADO
!     ! +++ work in progress +++
!     !----------------------
!     if (iand (TSK_Y,tsk) /= 0) then
!       call compute_radar_mod
!       tsk=tsk-TSK_Y
!     endif
!     if (tsk==0) return

!     !-----------
!     ! left tasks
!     !-----------
!     if (tsk /= 0) then
!       if (dace% lpio) write (6,*) 'process_radar: unknown task',tsk
!       call finish ('process_radar','unknown task')
!     endif

  end subroutine process_radar
!==============================================================================
  subroutine set_size
  !-----------------------------------------------
  ! store sizes of derived data type T_RADAR
  ! (in termes of size of component OBS% PAR) into
  ! private module variables RAY_INT_SIZE
  !-----------------------------------------------
    type (t_radar) :: ray
    type (t_obs)   :: obs
    if (ray_int_size == 0) then
      ray_int_size  = size (transfer (ray, obs% par))
      ray_byte_size = size (transfer (ray, (/' '/) ))
    endif
  end subroutine set_size
!------------------------------------------------------------------------------
  subroutine load_radar (obs, spot, rays)
  type (t_obs)   ,intent(in)  :: obs      ! data of all observations
  type (t_spot)  ,intent(in)  :: spot     ! meta data of this observation
  type (t_radar) ,pointer     :: rays (:) ! radar meta data
  !------------------------------------------------------------------
  ! Load the data from component PAR of OBS from position provided by
  ! SPOT. Store into RAYS. allocate RAYS with size required.
  !------------------------------------------------------------------
    integer ,pointer :: par (:)
    integer          :: m,n

    if (ray_int_size == 0) call set_size
    n = spot% s% n
    m = n * ray_int_size
    allocate (rays (n))
    par  => obs % par (spot% p% i+1 : spot% p% i + spot% p% n)
    rays =  transfer (par, rays)

  end subroutine load_radar
!------------------------------------------------------------------------------
  subroutine store_radar (obs, spot, rays)
  type (t_obs)   ,intent(inout) :: obs      ! data of all observations
  type (t_spot)  ,intent(inout) :: spot     ! meta data of this observation
  type (t_radar) ,intent(in)    :: rays (:) ! occultation data
  !-----------------------------------------------------------------------
  ! Store the data from variables RAYS in the component PAR of
  ! OBS at position provided by SPOT. Allocate memory for PAR if required.
  !-----------------------------------------------------------------------
    integer ,pointer :: par (:)
    integer          :: n, m

    if (sum(rays%radar_nbody) /= spot%o%n)           &
      call finish ('store_radar','inconsistent table')

    if (ray_int_size == 0) call set_size
    n = spot% s% n
    m = n * ray_int_size
    if (n /= size(rays)) call finish('store_radar','n /= size(rays)')
    if (spot% p% i < 0)  call new_par (obs, m, spot=spot)
    if (m < spot% p% n)  spot% p% n = m
    if (m > spot% p% n)  call finish('store_radar','m > spot% p% n')
    par => obs % par (spot% p% i+1 : spot% p% i + spot% p% n)
    par = transfer (rays, par)

  end subroutine store_radar
!==============================================================================
  subroutine thin_radar (obs)
  !-------------------------------
  ! thinning of radar observations
  !-------------------------------
  type (t_obs) ,intent(inout) :: obs ! observation container

    integer :: is         ! report      index
    integer :: i          ! observation index
    integer :: i0, i1, in ! observation (range) index bounds
    integer :: n          ! number of observations in record
    integer :: nw_invalid ! number of invalid wind observations
    integer :: nw_notused ! number of wind observations not used
    integer :: nr_notused ! number of reflectivities not used
    integer :: nu_notused ! number of unknown observations not used

!!$    integer :: i_str      ! stride to use vor range thinning
!!$    integer :: i_min      ! minimum range index
!!$    integer :: i_max      ! maximum range index
!!$    integer :: i          ! range index
!!$    integer :: j          ! azimut index
!!$    integer :: k          ! elevation index
!!$    integer :: jk         ! combined elevation*azimut index
!!$    integer :: ir (nra)   ! actual ranges in report
!!$    logical :: la (nra)   ! active ranges in report
!!$    integer :: si         ! spec_index(first obs.) - 1
!!$
!!$    i_str = nint (min_rdist * 1000._wp / dra)
!!$    i_min = nint (min_range * 1000._wp / dra)
!!$    i_max = nint (max_range * 1000._wp / dra)

!!$    jk = 0
    !------------------------------------------------------
    ! loop over reports, currently one azimut and elevation
    !------------------------------------------------------
    nw_invalid = 0
    nw_notused = 0
    nr_notused = 0
    nu_notused = 0
    do is = 1, obs% n_spot
      if (obs% spot(is)% hd% obstype /= OT_RADAR) cycle
      !--------------------------------------------------------------
      ! set index bounds for obs. in report, currently complete range
      !--------------------------------------------------------------
      n  = obs% spot(is)% o%n
      i0 = obs% spot(is)% o%i
      i1 = obs% spot(is)% o%i + 1
      in = obs% spot(is)% o%i + n
      !--------------------------------------------------------------
      ! check for observed quantity, currently all the same in report
      !--------------------------------------------------------------
      do i = i1, in
        select case (obs% varno (i1))
        case (VN_RADVEL)
          if (abs(obs% body(i)% o) >= 999._sp) then
            obs% body(i)% o = rvind
            if (obs% body(i)% use% state > STAT_DISMISS) nw_invalid = nw_invalid + 1
            call decr_use (obs% body(i)% use, STAT_DISMISS, CHK_INSDAT)
          endif
          if (obs% body(i)% use% state > use_radvel) then
            if (obs% body(i)% use% state >= STAT_ACTIVE) nw_notused = nw_notused + 1
            call decr_use (obs% body(i)% use, use_radvel, CHK_NOTUSED)
          endif
        case (VN_RREFL, VN_REFL)
          obs% varno (i) = VN_RREFL                  !+++ patch for old data +++!
          if (obs% body(i)% use% state > use_refl) then
            if (obs% body(i)% use% state >= STAT_ACTIVE) nr_notused = nr_notused + 1
            call decr_use (obs% body(i)% use, use_refl, CHK_NOTUSED)
          endif
        case default
          if (obs% body(i)% use% state > STAT_DISMISS) then
            nu_notused = nu_notused + 1
            call decr_use (obs% body(i)% use, STAT_DISMISS, CHK_NOTUSED)
          endif
        end select
!!$   if (obs% spot(is)% use% state <= STAT_DISMISS) cycle
!!$      !----------------------------------------------------------------------
!!$      ! check for valid spec_idex: spec_index = i + (j-1)*nra + (k-1)*nra*naz
!!$      ! may be incorrect due to short integer overflow in NetCDF
!!$      !----------------------------------------------------------------------
!!$      si = obs% body(i1)% spec_index-1
!!$      if (si >= 0) then
!!$        ir(1) = mod(si , nra) + 1
!!$        jk    =     si / nra
!!$      else
!!$        ir(1) = 1
!!$        jk    = jk + 1
!!$      endif
!!$      j = mod(jk , naz) + 1
!!$      k =     jk / naz  + 1
!!$      do i = 2, n
!!$        if (obs% body(i0+i)% spec_index > 0) then
!!$          ir(i) = mod((obs% body(i0+i)% spec_index-1) , nra) + 1
!!$        else
!!$          ir(i) = ir(i-1) + 1
!!$        endif
!!$      end do
!!$      !--------------------------------------------------
!!$      ! temporary storage of indices  in header and  body
!!$      !--------------------------------------------------
!!$      obs% body(i1:in)% lev_sig = ir(1:n)  ! range index
!!$      obs% spot(is)   % stzen   = ele (k)
!!$      obs% spot(is)   % phase   = nint (100 * (az_start + (j-1)*daz)) !   azimut
!!$      !---------
!!$      ! thinning
!!$      !---------
!!$      la (1:n) = .true.
!!$      do i = 1, n
!!$        if (ir(i) < i_min)            la(i) = .false.
!!$        if (ir(i) > i_max)            la(i) = .false.
!!$        if (0 < i_str) then
!!$          if (mod (ir(i),i_str) /= 0) la(i) = .false.
!!$        endif
!!$      end do
!!$      !---------
!!$      ! printout
!!$      !---------
!!$      if (iprintout > 0) then
!!$        where (.not.la (1:n)) ir(1:n) = - ir(1:n)
!!$        write(6,'(a,a,a,3i8,2i4,a,200i5)')                                 &
!!$          'RADAR ', obs% spot(is)% statid, chhmm(obs% spot(is)% hd% time), &
!!$          is, jk+1, si, j, k, ' --', ir(1:n)
!!$      endif
!!$      !---------------------
!!$      ! dismiss observations
!!$      !---------------------
!!$      if (any (la(1:n))) then
!!$        do i = 1, n
!!$          if (.not. la (i)) &
!!$            call decr_use (obs% body(i0+i)% use, STAT_DISMISS, CHK_THIN)
!!$        end do
!!$      else
!!$        call decr_rpt_use (obs% spot(is), CHK_THIN, STAT_DISMISS)
!!$      endif
      end do
    end do

    if (dace% lpio) then
      if (nu_notused > 0 .or. &
          nr_notused > 0 .or. &
          nw_notused > 0 .or. &
          nw_invalid > 0      ) then
        write(6,'(a)') repeat('-',79)
        write(6,'()')
        write(6,'(a)') '  thinning of radar observations'
        write(6,'()')
!!$     write(6,'(a,i6,4x,a)') 'i_str    =',i_str,'stride used for thinning'
!!$     write(6,'(a,i6,4x,a)') 'i_min    =',i_min,'minimum range index'
!!$     write(6,'(a,i6,4x,a)') 'i_max    =',i_max,'maximum range index'
!!$     write(6,'()')
        write(6,'(a,i10,1x,a)') '    unknown observations       :',nu_notused,stat_mnem(STAT_DISMISS)
        write(6,'(a,i10,1x,a)') '    reflectivities    not used :',nr_notused,stat_mnem(use_refl)
        write(6,'(a,i10,1x,a)') '    wind observations not used :',nw_notused,stat_mnem(use_radvel)
        write(6,'(a,i10,1x,a)') '    wind observations invalid  :',nw_invalid,stat_mnem(STAT_DISMISS)
        write(6,'()')
      endif
    endif
  end subroutine thin_radar
!==============================================================================
  subroutine shrink_par_radar (obs, mask)
  !-------------------------------------------------
  ! adapt table 'par' to reduced set of observations
  !-------------------------------------------------
  type (t_obs) ,intent(inout) :: obs      ! observation container
  logical      ,intent(in)    :: mask (:) ! indicates observations to keep

    target                 :: mask
    integer                :: is       ! report index
    logical       ,pointer :: m (:)    ! pointer to spot subset of mask
    type(t_spot)  ,pointer :: s        ! pointer to report header
    type(t_radar) ,pointer :: r (:)    ! pointer to table 'par'
    integer                :: i        ! ray index
    integer                :: i0, in   ! observation index bounds
    integer                :: n, l     ! observation index bounds

    do is = 1, obs% n_spot
      if (obs% spot(is)% hd% obstype /= OT_RADAR) cycle
      s => obs% spot(is)
      m => mask (s%o%i+1:s%o%i+s%o%n)
      n =  count (m)
      !----------
      ! no change
      !----------
      if (n == s%o%n) cycle
      !----------------
      ! no body entries
      !----------------
      if (n == 0    ) then
        s%s%n =  0
        s%p%n =  0
        s%p%i = -1
      else
        !-----------
        ! read table
        !-----------
        call load  (obs, s, r)
        !--------------
        ! count entries
        !--------------
        in = 0
        do i = 1, s%s%n
          i0 = in
          in = i0+r(i)% radar_nbody
          r(i)% radar_nbody = count (m(i0+1:in))
        end do
        !------------
        ! write table
        !------------
        l = s%o%n
        s%o%n = n               ! preliminary entry for consistency test
        call store (obs, s, r)
        s%o%n = l
        deallocate (r)
      endif
    end do

  end subroutine shrink_par_radar
!==============================================================================
  subroutine read_nml_radar
  !--------------------------
  ! read namelist /RADAR_OBS/
  !--------------------------
    integer :: ierr
    !-------------
    ! set defaults
    !-------------
!   min_range  =   0._wp      ! minimum range                       (km)
!   max_range  = 999._wp      ! maximum range                       (km)
!   min_hdist  =   0._wp      ! minimum horizontal distance between rays
!   min_vdist  =   0._wp      ! minimum vertical   distance between rays
!   min_rdist  =   0._wp      ! minimum distance within a ray       (km)
    use_refl   = STAT_ACTIVE  ! radar reflectivity usage flag
    use_radvel = STAT_ACTIVE  ! radial velocity    usage flag
    iprintout  =   0          ! steering of printout
    luse_eo_factor = .false.  ! eo from fdbk as relative factor
    dealias_fg = .true.       ! dealias radial wind (by first guess)
    chk_alias  = 2._wp        ! check dealiasing (compare to spread)
    ofg_alias  = 0.7_wp       ! check dealiasing (compare to o-fg)
    min_refl   = -999._wp     ! constrain reflectivity (obs,fg)
    !--------------
    ! read namelist
    !--------------
    if (dace% lpio) then
      call position_nml ('RADAR_OBS', status=ierr)
      select case (ierr)
      case (POSITIONED)
#if defined(__ibm__)
        read (nnml ,nml=RADAR_OBS, iostat=ierr)
        if(ierr/=0) call finish('read_nml_radar','ERROR in namelist /RADAR_OBS/')
#else
        read (nnml ,nml=RADAR_OBS)
#endif
      end select
    endif
    !------------------------
    ! broadcast to other PE's
    !------------------------
!   call p_bcast (min_range  ,dace% pio)
!   call p_bcast (max_range  ,dace% pio)
!   call p_bcast (min_hdist  ,dace% pio)
!   call p_bcast (min_vdist  ,dace% pio)
!   call p_bcast (min_rdist  ,dace% pio)
    call p_bcast (use_refl   ,dace% pio)
    call p_bcast (use_radvel ,dace% pio)
    call p_bcast (iprintout  ,dace% pio)
    call p_bcast (luse_eo_factor ,dace% pio)
    call p_bcast (dealias_fg ,dace% pio)
    call p_bcast (chk_alias  ,dace% pio)
    call p_bcast (ofg_alias  ,dace% pio)
    call p_bcast (min_refl   ,dace% pio)
    !---------------------------
    ! print namelist /RADAR_OBS/
    !---------------------------
    if (dace% lpio) then
      write(6,'(a)') repeat('-',79)
      write(6,'()')
      write(6,'(a)') '  namelist /RADAR_OBS/'
      write(6,'()')
!     write(6,'(a,f10.3,a)') 'min_range =',min_range ,'minimum range      (km)'
!     write(6,'(a,f10.3,a)') 'max_range =',max_range ,'maximum range      (km)'
!     write(6,'(a,f10.3,a)') 'min_hdist =',min_hdist ,'minimum hor. distance'
!     write(6,'(a,f10.3,a)') 'min_vdist =',min_vdist ,'minimum ver. distance'
!     write(6,'(a,f10.3,a)') 'min_rdist =',min_rdist ,'minimum distance in ray'
      write(6,'(a,i6,4x,a)') 'use_refl  =',use_refl  ,'reflectivity usage flag'
      write(6,'(a,i6,4x,a)') 'use_radvel=',use_radvel,'radial velocity usage'
      write(6,'(a,i6,4x,a)') 'iprintout =',iprintout ,'printout flag'
      write(6,'(a,l1,5x,a)') 'luse_eo_factor=',luse_eo_factor,'eo from fdbk as relative factor'
      write(6,'(a,l1,9x,a)') 'dealias_fg=',dealias_fg,'dealias wind using fg'
      write(6,'(a,f10.3,a)') 'chk_alias =',chk_alias, 'check dealiasing (spread)'
      write(6,'(a,f10.3,a)') 'ofg_alias =',ofg_alias, 'check dealiasing (spread)'
      write(6,'(a,f10.3,a)') 'min_refl  =',min_refl,  'constrain reflectivity'
      write(6,'()')
    endif
  end subroutine read_nml_radar
!==============================================================================
  subroutine constrain_refl (o, H_det, H_x, det_run)
  type (t_obs)    ,intent(inout) :: o    (:)  ! observation
  type (t_vector) ,intent(inout) :: H_det     ! first guess
  type (t_vector) ,intent(inout) :: H_x  (:)  ! ensemble fg
  logical         ,intent(in)    :: det_run   ! deterministic run present
  !-----------------------------------------------------
  ! constrain reflectivity (observation and first guess)
  ! to minimum value (min_refl)
  !-----------------------------------------------------

    integer  :: ib             ! box index
    integer  :: is             ! report index
    integer  :: i              ! observation index
    integer  :: k              ! ensemble index
    integer  :: i1, in         ! observation index bounds

    !-----------------------------
    ! loop over radar reflectivity
    !-----------------------------
    if (min_refl <= -999._wp ) return
    do ib = 1, size(o)
      if (o(ib)% pe /= dace% pe) cycle
      do is = 1, o(ib)% n_spot
        if (o(ib)% spot(is)% hd% obstype /= OT_RADAR) cycle
        i1 = o(ib)% spot(is)% o%i + 1
        in = o(ib)% spot(is)% o%i + o(ib)% spot(is)% o%n
        do i = i1, in
          if (o(ib)% varno(i) /= VN_RREFL) cycle
          !-----------------------
          ! constrain reflectivity
          !-----------------------
          o(ib)% body(i)% o     = max (o(ib)% body(i)% o,   real (min_refl,sp))
          do k = 1, size (H_x)
            H_x(k)% s(ib)% x(i) = max (H_x(k)% s(ib)% x(i), min_refl)
          end do
          if (det_run) then
            H_det% s(ib)% x(i)  = max (H_det% s(ib)% x(i),  min_refl)
          endif
        end do
      end do
    end do

  end subroutine constrain_refl
!==============================================================================
  subroutine check_radar (o, H_det, H_x, det_run)
    type (t_obs)    ,intent(inout) :: o    (:)  ! observation
    type (t_vector) ,intent(in)    :: H_det     ! first guess
    type (t_vector) ,intent(in)    :: H_x  (:)  ! ensemble fg
    logical         ,intent(in)    :: det_run   ! deterministic run present

    integer  :: ib             ! box index
    integer  :: is             ! report index
    integer  :: i              ! observation index
    integer  :: k              ! ensemble index
    integer  :: i1, in         ! observation index bounds
    integer  :: counter        ! invalid model equivalent spot counter

    counter = 0
    do ib=1, size(o)
      do is = 1, o(ib)% n_spot
        if (o(ib)% spot(is)% hd% obstype /= OT_RADAR) cycle
        i1 = o(ib)% spot(is)% o%i + 1
        in = o(ib)% spot(is)% o%i + o(ib)% spot(is)% o%n
        do i = i1, in
          ! Check all radar obs. (radial wind and reflectivity)
!         if (o(ib)% varno(i) /= VN_RADVEL) cycle
          ! Set observation to STAT_OBS_ONLY if one of the model equivalents is invalid:
          do k = 1, size (H_x)
            if (abs(H_x(k)% s(ib)% x(i))   > 999._wp .and. &
                o(ib)% body(i)% use% state > STAT_OBS_ONLY) then
              counter = counter + 1
              call decr_use (o(ib)% body(i)% use, STAT_OBS_ONLY, CHK_OPERATOR)
              exit
            end if
          end do
        end do
      end do
    end do

    if (counter > 0 ) then
      write (*, '(a,i7,a)') 'check_radar(): ', counter, &
           ' model equivalents with fill_values found where the corresponding&
           & obs were flagged as ACTIVE or PASSIVE!'
    end if

  end subroutine check_radar
!==============================================================================
  subroutine alias_radar (o, H_det, H_x, det_run)
  type (t_obs)    ,intent(inout) :: o    (:)  ! observation
  type (t_vector) ,intent(in)    :: H_det     ! first guess
  type (t_vector) ,intent(in)    :: H_x  (:)  ! ensemble fg
  logical         ,intent(in)    :: det_run   ! deterministic run present
  !---------------------------------------------------------------
  ! de-alias the radial wind observation
  ! assumes: model fg is de-aliased, observation is not de-aliased
  !---------------------------------------------------------------

    integer  :: ib             ! box index
    integer  :: is             ! report index
    integer  :: i              ! observation index
    integer  :: k              ! ensemble index
    integer  :: i1, in         ! observation index bounds
    integer  :: ie, id         ! nyquist frequency index
    real(wp) :: min_o,  max_o  ! min/max observed value
    real(wp) :: min_fg, max_fg ! min/max of first guess ensemble
    real(wp) :: obs            ! observed value
    real(wp) :: mean_fg        ! mean of first guess ensemble
    real(wp) :: vnyquist       ! nyquist wind speed
    real(wp) :: offset         ! correction term
    real(wp) :: cobs           ! corrected observation

    !-----------------------------
    ! loop over radar radial winds
    !-----------------------------
    if (.not. dealias_fg) return
    do ib = 1, size(o)
      if (o(ib)% pe /= dace% pe) cycle
      do is = 1, o(ib)% n_spot
        if (o(ib)% spot(is)% hd% obstype /= OT_RADAR) cycle
        vnyquist = o(ib)% spot(is)% params(1)
        i1 = o(ib)% spot(is)% o%i + 1
        in = o(ib)% spot(is)% o%i + o(ib)% spot(is)% o%n
        min_o =  huge (min_o)
        max_o = -huge (min_o)
        do i = i1, in
          if (o(ib)% varno(i) /= VN_RADVEL) cycle
          select case (o(ib)% body(i)% use% state)
          case (:STAT_OBS_ONLY, STAT_PAS_REJ, STAT_REJECTED)
             cycle                ! do not de_alias
          end select
          !-----------------------------------
          ! check for valid Nyquist wind speed
          !-----------------------------------
          if (vnyquist <= 0._wp) then
             call decr_use (o(ib)% body(i)% use, STAT_REJECTED, CHK_FG)
             cycle                ! bad, do not de_alias
          end if
          !------------------------------------
          ! check if the observation is aliased
          !------------------------------------
          obs = o(ib)% body(i)% o
          min_o = min (min_o, obs)
          max_o = max (max_o, obs)
!         if (abs (obs) > vnyquist) then
!           call finish ('alias_radar','abs(obs) > vnyquist')
!         endif
          !----------------------------------
          ! de-aliasing for the ensemble mean
          !----------------------------------
          min_fg  = H_x(1)% s(ib)% x(i)
          max_fg  = min_fg
          mean_fg = min_fg
          do k = 2, size (H_x)
            min_fg  = min (min_fg, H_x(k)% s(ib)% x(i))
            max_fg  = max (max_fg, H_x(k)% s(ib)% x(i))
            mean_fg = mean_fg +    H_x(k)% s(ib)% x(i)
          end do
          mean_fg = mean_fg / size (H_x)
          ie      = nint ((mean_fg - obs) / (2._wp * vnyquist))
          !--------------------------------------------------------------
          ! dismiss report if spread is too large
          ! or offset for ensemble and deterministic run are inconsistent
          !--------------------------------------------------------------
          if (max_fg - min_fg > chk_alias * vnyquist)                &
            call decr_use (o(ib)% body(i)% use, STAT_REJECTED, CHK_FG)
          if (det_run) then
            id = nint( (H_det% s(ib)% x(i) - obs) / (2._wp*vnyquist))
            if (ie /= id)                                              &
              call decr_use (o(ib)% body(i)% use, STAT_REJECTED, CHK_FG)
          endif
          !------------------
          ! actually de-alias
          !------------------
          offset  = vnyquist * 2 * ie
          cobs    = obs + offset
          o(ib)% body(i)% o = cobs
          if (abs(mean_fg - cobs) > vnyquist)                &
            call finish ('alias_radar','abs(o-fg) > vnyquist')
          !----------------------------------------------------------------
          ! dismiss report if obs - fg is too large in relation to vnyquist
          !----------------------------------------------------------------
          if (abs (mean_fg - cobs) > ofg_alias * vnyquist)           &
            call decr_use (o(ib)% body(i)% use, STAT_REJECTED, CHK_FG)
          if (det_run) then
            if (abs (H_det% s(ib)% x(i) - cobs) > ofg_alias * vnyquist)&
              call decr_use (o(ib)% body(i)% use, STAT_REJECTED, CHK_FG)
          endif
        end do
        !-------------------------------------------------
        ! issue WARNING if the observation was not aliased
        !-------------------------------------------------
        if (min_o < - vnyquist .or. &
            max_o >   vnyquist      ) then
          write(6,*) 'WARNING alias_radar: abs(obs) > vnyquist: min,max,nyq=',&
            o(ib)% spot(is)% statid, min_o, max_o, vnyquist
        endif
      end do
    end do
  end subroutine alias_radar
!==============================================================================
  subroutine read_fdbk_radar (o, n_radar)
  type (t_obs) ,intent(inout) :: o          ! observation data type variable
  integer      ,intent(in)    :: n_radar    ! total number of table entries
  !-----------------------------------------------
  ! read radar specific entries from feedback file
  !-----------------------------------------------

    type (t_radar) ,allocatable ,target :: rt (:) ! RADAR specific table
    type (t_radar) ,pointer             :: r  (:) ! pointer to table
    type (t_spot)  ,pointer             :: s      ! pointer to report
    integer                             :: i      ! report index
!   integer                             :: s1, s2 ! index bounds
    logical                             :: ok     ! true for consistent data
    allocate (rt (n_radar))

    !----------------------------
    ! read radar specific entries
    !----------------------------
    stanc = 1
    counc = n_radar
    strnc = 1
    call get_var (rt% radar_azimuth     ,'radar_azimuth'    )
    call get_var (rt% radar_elevation   ,'radar_elevation'  )
    call get_var (rt% radar_nrange      ,'radar_nrange'     )
    call get_var (rt% radar_range_start ,'radar_range_start')
    call get_var (rt% radar_drange      ,'radar_drange'     )
    call get_var (rt% radar_nbody       ,'radar_nbody'      )

    !------------------------------
    ! store in component t_obs% par
    !------------------------------
    do i = 1, o% n_spot
      if  (o% spot(i)% hd% obstype /= OT_RADAR) cycle
      s => o% spot(i)
      r => rt (s%s%i : s%s%i + s%s%n - 1)
      call store_radar (o, s, r)
      call set_ray_coordinates (s, o%body, o%olev, r, ok)
      if (.not.ok) call decr_rpt_use (s, CHK_INSDAT, STAT_DISMISS)
      !--------------------------------------------------------
      ! misuse station zenith angle to store elevation
      ! (currently the same for all observations in the report)
      !--------------------------------------------------------
      if (all (r% radar_elevation == r(1)% radar_elevation)) then
        s% stzen = r(1)% radar_elevation
      endif
    end do

    deallocate (rt)

  end subroutine read_fdbk_radar
!------------------------------------------------------------------------------
  subroutine radar_tables  (radar, spt, obs, ix, ib, n, mask)
  type (t_radar)  ,pointer    :: radar (:) ! table to be allocated and filled
  type (t_spot)   ,intent(in) :: spt   (:) ! observation headers
  type (t_obs)    ,intent(in) :: obs   (:) ! observation 'boxes'
  integer         ,intent(in) :: ix    (:) ! offsets for individual reports
  integer         ,intent(in) :: ib    (:) ! observation box index
  integer         ,intent(in) :: n         ! total # of radar table entries
  type (t_bvector),intent(in) :: mask      ! mask for valid observations
  !---------------------------------------------------------------
  ! prepare the radar specific tables for writing to feedback file
  !---------------------------------------------------------------

    integer                 :: i      ! report index
    integer                 :: k      ! ray index
    integer                 :: j      ! observation index
    integer                 :: jk     ! number of observations per ray
    type (t_radar) ,pointer :: r (:)  ! radar operator table for 1 report
    logical        ,pointer :: m (:)  ! mask

    !---------------
    ! allocate table
    !---------------
    allocate (radar (n))
    if (n == 0) return

    !-----------------------------
    ! decode radar specific tables
    !-----------------------------
    do i = 1, size (spt)
      if (spt(i)% hd% obstype /= OT_RADAR) cycle
      if (spt(i)% s% n        == 0       ) cycle
      call load (obs(ib(i)), spt(i), r)
      if (sum (r(:)% radar_nbody) /= spt(i)%o%n)     &
      call finish('radar_tables','inconsistent table')
      !-----------------------------------------------------
      ! make radar_nbody consistent with active observations
      !-----------------------------------------------------
      m  => mask% s(ib(i))% x (spt(i)%o%i+1:spt(i)%o%i+spt(i)%o%n)
      j = 0
      do k = 1, size(r)
        jk = r(k)% radar_nbody
        r(k)% radar_nbody = count (m (j+1:j+jk))
        j = j + jk
      end do
      !-----------------------
      ! combine to final table
      !-----------------------
      radar (ix(i)+1 : ix(i)+spt(i)% s% n) = r
      deallocate (r)
    end do

  end subroutine radar_tables
!==============================================================================
  subroutine set_ray_coordinates (s, b, olev, rs, ok)
  type (t_spot)  ,intent(inout) :: s       ! report variable
  type (t_datum) ,intent(inout) :: b   (:) ! observation body variable
  real(wp)       ,intent(inout) :: olev(:) ! level of observation
  type (t_radar) ,intent(in)    :: rs  (:) ! radar rays variable
  logical        ,intent(out)   :: ok      ! true for consistent data

    !----------------
    ! local variables
    !----------------
    target                 :: rs
    integer                :: ir        ! ray index
!   integer                :: nr        ! observations in ray
    integer                :: io,i0, in ! observation index bounds
    type(t_radar) ,pointer :: r         ! pointer to table
    !--------------------------------
    ! parameters to rad2geo_const_vec
    !--------------------------------
    real(wp)              :: latsta     ! geographical latitude of radar station
    real(wp)              :: lonsta     ! geographical longitude of radar station
    real(wp)              :: altsta     ! altitude of radar station asl (m)
    integer               :: nra        ! number of ranges
    integer               :: naz = 1    ! number of azimuths
    integer               :: nel = 1    ! number of elevations
    real(wp) ,allocatable :: ra (:)     ! array of ranges
    real(wp)              :: az (1)     ! array of azimuths
    real(wp)              :: el (1)     ! array of elevations
    !----------------------------------
    ! parameters from rad2geo_const_vec
    !----------------------------------
    real(wp) ,allocatable :: lat (:)    ! geographical latitude
    real(wp) ,allocatable :: lon (:)    ! geographical longitude
    real(wp) ,allocatable :: h   (:)    ! height above mean sea level

    nra = maxval (rs% radar_nbody)
    allocate (ra (nra))
    allocate (lat(nra))
    allocate (lon(nra))
    allocate (h  (nra))

    !------------------
    ! consistency check
    !------------------
    if (dace% lpio .and. iprintout > 0) then
      write(6,*) '--------------------------------------------------------------------------------------'
      write(6,*) 'i, azimuth, elevation, nrange, range_start, drange, nbody, 1st, lst, nobs - spec_index'
    endif
    io = s% o% i
    i0 = 0
    ok = .true.
    do ir = 1, s% s% n
      r => rs (ir)
      in = i0 + r% radar_nbody
      !------------------
      ! consistency check
      !------------------
      if (                  r% radar_nbody <  0                              ) ok = .false.
      if (any(b(io+i0+2:io+in)% spec_index <= b(io+i0+1:io+in-1)% spec_index)) ok = .false.
      if (any(b(io+i0+1:io+in)% spec_index <  1                             )) ok = .false.
      if (any(b(io+i0+1:io+in)% spec_index >  r% radar_nrange               )) ok = .false.
      !---------
      ! printout
      !---------
      if (dace% lpio .and. iprintout > 0) then
        write(6,*) ir, r% radar_azimuth, r% radar_elevation, r% radar_nrange, &
                       r% radar_range_start, r% radar_drange, r% radar_nbody, &
                       i0+1,in,s% o% n,' - ',b(io+i0+1:io+in)% spec_index
      endif
      !------------------------
      ! restore ray coordinates
      !------------------------
      if (ok) then
        latsta    = s%col%c% dlat
        lonsta    = s%col%c% dlon
        altsta    = s      % z
        nra       = r% radar_nbody
        ra (:nra) = r% radar_range_start + r% radar_drange * (b(io+i0+1:io+in)% spec_index -1)
        az        = r% radar_azimuth
        el        = r% radar_elevation
        call rad2geo_const_vec (latsta,lonsta,altsta,naz,nra,nel,ra,az,el,lat,lon,h)
        b(io+i0+1:io+in)% lat        = lat(:nra)
        b(io+i0+1:io+in)% lon        = lon(:nra)
!        b(io+i0+1:io+in)% l2c        = ra (:nra)
        b(io+i0+1:io+in)% obs_par(1) = r% radar_azimuth
        b(io+i0+1:io+in)% obs_par(2) = r% radar_elevation
        !---------
        ! printout
        !---------
        if (dace% lpio .and. iprintout > 0) then
          write(6,*) ir, r% radar_azimuth, r% radar_elevation, r% radar_nrange, &
                         r% radar_range_start, r% radar_drange, r% radar_nbody, &
                         i0+1,in,s% o% n,' - ',real(olev(io+i0+1:io+in),sp)
          write(6,*) ir, r% radar_azimuth, r% radar_elevation, r% radar_nrange, &
                     r% radar_range_start, r% radar_drange, r% radar_nbody,     &
                     i0+1,in,s% o% n,' - ',real(h(:nra),sp)
        endif
      endif
      i0 = in
    end do
    if (in /= s% o% n) ok = .false.

    if (.not. ok) call message ('set_ray_coordinates','inconsistent rays: '//s%statid)

  end subroutine set_ray_coordinates
!------------------------------------------------------------------------------
  !==============================================================================
  !+ Module procedure from src_radar for the coordinate transformation of radar
  !  radar coordinates to geographical coordinates for constant beam propagation
  !------------------------------------------------------------------------------

  SUBROUTINE rad2geo_const_vec (latsta,lonsta,altsta,naz,nra,nel,ra,az,el,lat,lon,h)

    !------------------------------------------------------------------------------
    !
    ! Description: Coordinate transformation of radar coordinates
    !              (range, azimuth and elevation) to
    !              geographical coordinates (lat, lon and height asl)
    !
    ! Method:      Assumption: 4/3 earth model
    !              Formulas with reference to  Appendices B and D of Dissertation
    !              of Uli Blahak
    !
    !------------------------------------------------------------------------------
    !
    ! Subroutine / Function arguments
    ! Scalar arguments with intent(in):

    !------------------------------------------------------------------------------
    ! Parameter list:
    REAL (KIND=wp), INTENT (IN)      ::        &
         latsta,  & ! geographical latitude of radar station
         lonsta,  & ! geographical longitude of radar station
         altsta     ! altitude of radar station asl (m)

    INTEGER, INTENT (IN)::        &
         nra,naz,nel! number of ranges, azimuths and elevations (array dimensions)

    REAL (KIND=wp), INTENT (IN)      ::        &
         ra(nra), & ! array of ranges
         az(naz), & ! array of azimuths
         el(nel)    ! array of elevations

    REAL (KIND=wp), INTENT (OUT)     ::        &
         lat(naz*nra*nel),& ! geographical latitude
         lon(naz*nra*nel),& ! geographical longitude
         h(naz*nra*nel)     ! height above mean sea level
    !
    ! Local scalars:
!   CHARACTER (LEN=32) yzroutine

    INTEGER              :: ira,iaz,iel,irp
    REAL    (KIND=wp)    :: h1,s1,d1,d2,d3,re


    !- End of header
    !==============================================================================

    !------------------------------------------------------------------------------
    !- Begin SUBROUTINE rad2geo_const_vec
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    ! Section 1: Initializations
    !------------------------------------------------------------------------------

!   yzroutine(:) = ' '
!   yzroutine = 'rad2geo_const_vec'

    re = 4.0_wp/3.0_wp*(rearth + altsta)

    !if (my_cart_id == 0) then
    !   write(my_cart_id+700,'(2(A7),3(A10,"(",A3,")"),3A12') 'latsta','lonsta','az','m','ra','n','el','o','lon','lat','alt'
    !end if

    DO ira = 1,nra
      DO iel = 1,nel
        DO iaz = 1,naz
          ! ZY >> use new index irp to replace m,n,o indices
          ! Determine the index of the radar point, based on az,ra,el
          irp = iaz + (ira - 1)*naz  + (iel - 1)*naz*nra
          ! ZY <<

          ! ZY >> Vectorized Version
          ! calculate height of radar beam above surface (formula B.8)
          h1 = SQRT(re**2 + ra(ira)**2 + 2*re*ra(ira)*SIN(el(iel)*d2r)) - re
          ! calculate length of circle at height of radar station (formula B.9)
          s1 = re*ASIN(ra(ira)*COS(el(iel)*d2r)/(re+h1))
          ! calculate dummy variable
          d1 = s1/(rearth+altsta)

          ! calculate height above mean sea level (only dependent on range and elevation!)
          ! (formula B.14)

          h(irp) =  altsta + h1

          ! compute geographical latitude (dependent on range, azimuth and elevation)
          ! (formula D.4)
          d2 = SIN(latsta*d2r)*COS(d1) + COS(latsta*d2r)*SIN(d1)*COS(az(iaz)*d2r)
          ! restrict argument of asin to [-1,1]
          d2 = MIN(MAX(d2,-1.0_wp),1.0_wp)
          lat(irp) = ASIN(d2)

          ! compute geographical longitude (dependent on range, azimuth and elevation)
          ! (formula D.5)
          d3 =(COS(d1)-SIN(latsta*d2r)*SIN(lat(irp)))/(COS(latsta*d2r)*COS(lat(irp)))
          ! restrict argument of acos to [-1,1]
          d3 = MIN(MAX(d3,-1.0_wp),1.0_wp)
          lon(irp) = lonsta*d2r + SIGN(1.0_wp,(pi - az(iaz)*d2r)) * ACOS(d3)
          ! ZY <<

          !if (my_cart_id == 0) then
          !   write(my_cart_id+700,'(2(F7.2),3(F10.2,"(",I3,")"),3F12.5)') latsta, lonsta, az(iaz),iaz, &
          !         ra(ira),ira,el(iel),iel,lon(iaz,ira,iel)*r2d, lat(iaz,ira,iel)*r2d,h(ira,iel)
          !end if

        END DO
      END DO
    END DO

    lat = lat * r2d
    lon = lon * r2d

  END SUBROUTINE rad2geo_const_vec
!==============================================================================
end module mo_radar
