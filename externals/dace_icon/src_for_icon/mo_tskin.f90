!
!+ provides skin temperature retrieval options
!
MODULE mo_tskin
!
! Description:
!   skin temperature: dynamic retrieval
!   email:  mahdiyeh.mousavi@dwd.de
!
!   Current Code Owner: DWD, Mahdiyeh Mousavi
!
! History:
! Version      Date       Name
! ------------ ---------- ----
!
! Code Description:
! Language: Fortran 2003
! Software Standards:
!
!==============================================================================
  !=============
  ! modules used
  !=============
  !-------------------------
  ! general purpose routines
  !-------------------------
  use mo_kind,              only: wp                        ! working precision kind parameter

  use mo_exception,         only: finish                    ! abort in case of error

  use mo_mpi_dace,          only: dace,                    &! MPI group info
                                  p_sum,                   &!
                                  p_max

  use mo_namelist,          only: position_nml,            &! position namelist
                                  nnml,                    &! namelist Fortran unit number
                                  POSITIONED                ! ok    code from position_nml
  use mo_instrid,           only: ir_instr

  !-----------------------
  ! access to observations
  !-----------------------
  use mo_rad,               only: m_chan,                  &! max number of channels
                                  rad_set,                 &! radiance meta data
                                  n_set,                   &!
                                  t_rad_set,               &! rad_set type
                                  n_styp,                  &!
                                  t_tskin_opt,             &!
                                  assign_tskin_opt,        &!
                                  instr_type,              &!
                                  ITYP_IR,                 &!
                                  construct, destruct,     &!
                                  OPTV_CLD_FRC

  use mo_t_obs,             only: t_obs,                   &! observation derived type
                                  t_spot,                  &! spot derived type
                                  ldeb,usd,dpref            ! debug selected spot
  use mo_fdbk_tables,       only: OT_RAD,                  &! Radiances report type ID
                                  TF_TSKIN,                &
                                  TF_TSKIN_DYNRET,         &
                                  TF_TSKIN_FAILED

  !-------------------------
  ! rttov related modules
  !-------------------------
  use mo_t_tovs,            only: t_tovs,                  &! observation operator specific type
                                  t_tovs_instr,            &! information on instruments in t_tovs
                                  get_tovs_rs,             &! get rad_set and instrument info for t_tovs
                                  store,                   &
                                  TTOVS_FLAG, TTOVS_CI,    &
                                  TTOVS_BASE

  use mo_rttov,             only: t_rttov_prof,            &! arguments to rttov
                                  call_rttov                ! call to rttovs interface subroutines

  use mo_rtifc,             only: rts_name,                &!
                                  rts_land,rts_sea,rts_ice  !

  implicit none

  !================
  ! public entities
  !================
  private

  public :: set_tskin
  public :: read_nml_tskin
  public :: print_tskin_stat
  public :: TSKIN_CALC


  !=================
  ! Parameters
  !=================
  integer, parameter :: m_tskin            = 50      ! Max number of tskin to be set via namelist
  integer, parameter :: MODE_TSKIN_DEFAULT = 0       ! Code for default tskin
  integer, parameter :: MODE_TSKIN_DYNRET  = 1       ! Code for dynamical retrieval of tskin

  integer, parameter :: TSKIN_SUCCESS = 2**TF_TSKIN(0)
  integer, parameter :: TSKIN_CALC    = TSKIN_SUCCESS + &
                                        2**TF_TSKIN(1)


  !=================
  ! Module variables
  !=================
  type(t_tskin_opt),  target :: tskin_set(m_tskin)
  integer                    :: n_nml           = 0
  integer                    :: id              = 0
  logical                    :: mask_pr(m_chan) = .false. ! used in set_tskin
  logical                    :: mask_dr(m_chan) = .false. ! used in set_tskin

  ! Statistics
  type t_tskin_stat
    integer                  :: id       =  -1     ! same as t_tskin_opt ID
    integer                  :: satid    =  -1
    integer                  :: instr    =  -1
    integer                  :: styp     =  -1
    character(len=120)       :: descr    =  ''
    integer                  :: n        =   0
    integer                  :: n_f      =   0     ! failed
    integer                  :: n_mix    =   0     ! mixed surfaces
    integer                  :: n_mix_f  =   0     ! mixed surfaces failed
    real(wp)                 :: tskin    =   0._wp
    real(wp)                 :: dt_mn    =   0._wp
    real(wp)                 :: dt_s2    =   0._wp
  end type t_tskin_stat
  type(t_tskin_stat), pointer :: stat(:)  => NULL()
  integer                     :: nstat    =  0
  logical                     :: l_pr     =  .true.
  
contains

  !----------------------------------------------------------------------------
  ! Reads /TSKIN/ namelist into tskin_set
  !----------------------------------------------------------------------------
  subroutine read_nml_tskin_sub
    character(len=*),  parameter :: proc = 'read_nml_tskin_sub'
    type(t_tskin_opt), pointer   :: tso   => NULL()
    integer                      :: i_read
    logical                      :: first
    integer                      :: ierr_pos
    integer                      :: ityp, k, n_cdyn
    character(len=300)           :: msg
    character(len=80)            :: hint
    integer                      :: sfchl(4)
    integer                      :: i, ns
    ! Namelist content
    integer                      :: inst
    integer                      :: prio
    integer                      :: styp(n_styp)
    integer                      :: mode
    integer                      :: cdyn (m_chan)
    character(len=120)           :: descr


    namelist /TSKIN/ inst, prio, styp, mode, cdyn

    write(*,*) 'READ TSKIN namelists:'
    first  = .true.
    i_read = 0
    n_nml  = 1
    loop_tskin: do
      call position_nml ('TSKIN' ,lrewind=first ,status=ierr_pos)
      if (ierr_pos == POSITIONED) then
        i_read = i_read + 1
        call init
        write(hint, '("TSKIN namelist number ",I2,1x)') i_read
        read (nnml ,nml=TSKIN)
        call process_nml
      else
        exit loop_tskin
      end if
      first = .false.
    end do loop_tskin

    n_nml = n_nml - 1

  contains

    subroutine init
      inst         = -1
      prio         =  0
      styp(:)      = -1
      mode         = -1
      cdyn(:)      = -1
      descr        = ''
    end subroutine init


    subroutine process_nml

      if (inst < 0) call finish(proc, 'Missing instrument '//trim(hint))

      if (descr == '') then
        ityp = 0
        if (ir_instr(inst)) ityp = ibset(ityp, ITYP_IR)
        descr = trim(hint)//' '//c_tskin(mode, ityp, 0)
      end if
      hint = trim(hint)//' "'//trim(descr)//'"'
      n_cdyn = count(cdyn > 0)

      ! Check if channels chosen for retrival are valid
      if (n_cdyn > 0) then
        if (.not. all(cdyn(1:n_cdyn) > 0)) call finish(proc, trim(hint)//': invalid cdyn')
        ! Check if availability of cdyn channels agrees with mode
        if (mode == MODE_TSKIN_DEFAULT) call finish(proc,  trim(hint)//': mode not compatible with cdyn')
      end if

      ! Check surface type validity and assign proper stypes
      if (any(styp(:) < -1 .or. styp(:) > n_styp)) then
        write(msg,*) trim(hint)//': invalid styp=',styp(:)
        call finish('read_nml_tskin',trim(msg))
      end if
      ns = count(styp >= 0)
      if (ns == 0) then
        styp = (/ 0, 1, 2/)
        ns   = 3
      end if
      styp(1:ns) = pack(styp, mask=(styp >= 0))
      ! Check mode
      select case (mode)
      case (MODE_TSKIN_DYNRET)
        ! Check if instrument is supported
        if (.not. any(inst == (/21,16,44,56/)) ) then
          call finish('read_nml_tskin','Tskin retrieval not supported for this instrument')
        end if
      case default
        write(msg, *) 'Default tskin.'
      end select

      ! Fill corresponding tskin_set element
      if (n_nml > size(tskin_set)) call  finish(proc, 'Too many TSKIN namelists. Recompile with &
           &increased m_tskin')
      tso => tskin_set(n_nml)
      tso%nml          = n_nml
      tso%inst         = inst
      tso%prio         = prio
      tso%ns           = ns
      tso%styp         = styp
      tso%descr        = trim(descr)
      tso%mode         = mode
      tso%n_cdyn       = n_cdyn

      if (mode == MODE_TSKIN_DYNRET) then
        allocate(tso%cdyn(n_cdyn))
        do k = 1, n_cdyn
          tso%cdyn(k) = cdyn(k)
        end do
      end if

      n_nml = n_nml + 1
      ! printout
      write(*,*)
      write(*,'(3x,a)') trim(hint)
      write(*,'(5x,a,1X,a)')                     ' descr        = ', trim(tso%descr)
      write(*,'(5x,a,1X,i2)')                    ' inst         = ', tso%inst
      write(*,'(5x,a,1X,i2)')                    ' prio         = ', tso%prio
      write(*,'(5x,a,10(1X,i2))')                ' styp         = ', tso%styp(1:tso%ns)
      write(*,'(5x,a,1X,i2)')                    ' mode         = ', tso%mode
      if (tso%mode == MODE_TSKIN_DYNRET .and. tso%n_cdyn > 0 ) then
        write(*,'(5x,a,1X,i2)')                  ' ncdyn        = ', tso%n_cdyn
        write(*,'(5x,a,*(8(1x,i5,:),/,20x))')    ' cdyn         = ', tso%cdyn(1:tso%n_cdyn)
      end if

    end subroutine process_nml

  end subroutine read_nml_tskin_sub


  subroutine read_nml_tskin
    character(len=*),  parameter   :: proc = 'read_nml_tskin'
    type(t_rad_set),   pointer     :: rs
    type(t_tskin_opt), pointer     :: tso => null()
    type(t_tskin_opt)              :: aux
    integer,           allocatable :: instr_nml(:), n_prio_instr(:), n_set_instr(:), i0(:)
    integer                        :: i, ii, j, k, instr, iprio, is, ityp
    integer                        :: n0, n1
    integer                        :: n_tskin_opt
    integer                        :: n_instr_nml
    integer                        :: n_prio
    integer                        :: styp
    logical                        :: l_swap
    ! printout
    integer                        :: n_instr_print
    integer                        :: instr_print(100)

    ! read /TSKIN/ namelists
    call read_nml_tskin_sub

    !-----------------
    ! Prepare tskin_set
    !-----------------
    ! Sort for inst and prio
    !================ when there is only one namelist id does not enter this part =====
    do i = n_nml-1,1,-1
      do k = 1,i
        l_swap = .false.
        if (tskin_set(k)%inst >  tskin_set(k+1)%inst) then
          l_swap = .true.
        elseif (tskin_set(k)%inst == tskin_set(k+1)%inst) then
          if (tskin_set(k)%prio > tskin_set(k+1)%prio) then
            l_swap = .true.
          end if
        end if
        if (l_swap) then
          aux            = tskin_set(k)
          tskin_set(k)   = tskin_set(k+1)
          tskin_set(k+1) = aux
        end if
      end do
    end do

    ! count number of different sets and priorities per instrument
    allocate(instr_nml(n_nml), n_prio_instr(n_nml), n_set_instr(n_nml), i0(n_nml))
    instr_nml(:)    = -1
    n_prio_instr(:) = -1
    n_set_instr(:)  = -1
    i0(:)           = -1
    n_instr_nml = 0
    do i = 1, n_nml
      ii = -1
      if (i > 1) then
        if (tskin_set(i)%inst == tskin_set(i-1)%inst) ii = n_instr_nml
      end if
      if (ii < 0) then
        n_instr_nml = n_instr_nml + 1
        ii          = n_instr_nml
        instr_nml   (ii) = tskin_set(i)%inst
        i0          (ii) = i
        n_set_instr (ii) = 0
        n_prio_instr(ii) = 0
      end if
      n_set_instr(ii)  = n_set_instr(ii) + 1

      if (all(tskin_set(i)%prio /= tskin_set(i0(ii):i-1)%prio)) then
        n_prio_instr(ii) = n_prio_instr(ii) + 1
      end if
    end do


    !----------------------------------
    ! cycle over satellites/instruments
    !----------------------------------
    n_instr_print = 0
    id = 0
    do i = 1, size(rad_set)
      rs => rad_set(i)
      ! only
      if (any(rs%grid == tskin_set(1:n_nml)%inst) ) then
        ! Derive dimensions, allocate
        n_prio = 0
        n_tskin_opt = 0
        do instr = 1, rs%n_instr
          ityp = instr_type(rs, instr)
          ! Add default
          if (btest(ityp, ITYP_IR)) n_tskin_opt = n_tskin_opt + 1
          call find_instr(ii)
          if (ii > 0) then
            n_tskin_opt = n_tskin_opt + n_set_instr(ii)
            n_prio = max(n_prio, n_prio_instr(ii)+ 1)
          end if
        end do
        n_prio = max(n_prio, 1)
        allocate(rs%tskin_opt(n_tskin_opt), rs%tskin_ind(n_prio, 0:n_styp-1, rs%n_chan))
        rs%tskin_ind = -1


        n_tskin_opt = 0
        do instr = 1, rs%n_instr
          if (any(rs%instr(instr) == tskin_set(:)%inst) ) then
            n0 = rs%o_ch_i(instr) + 1
            n1 = rs%o_ch_i(instr) + rs%n_ch_i(instr)
            ityp = instr_type(rs, instr)
            ! ! Add defaults
            if (btest(ityp, ITYP_IR)) then
              n_tskin_opt = n_tskin_opt + 1
              tso => rs%tskin_opt(n_tskin_opt)
              call construct(tso)
              tso%inst  = rs%instr(instr)
              tso%mode  = MODE_TSKIN_DEFAULT
              tso%ns    = 3 ! number of surface types
              tso%styp  = (/ 0, 1, 2/)
              tso%descr = 'Default method: '//c_tskin(tso%mode, ityp, 0)
              tso%id    = id ; id = id + 1
              rs%tskin_ind(:,:,n0:n1) = n_tskin_opt
            end if
            ! Add options from namelist
            call find_instr(ii)
            if (ii > 0) then
              do j = i0(ii), i0(ii) + n_set_instr(ii) - 1
                n_tskin_opt = n_tskin_opt + 1
                tso => rs%tskin_opt(n_tskin_opt)
                call assign_tskin_opt(tso, tskin_set(j))
                tso%id = id ; id = id + 1
                iprio = tskin_set(j)%prio + 1

                do is = 1, tso%ns
                  styp  = tso%styp(is)
                  ! all channels
                  rs%tskin_ind(iprio,styp,n0:n1) = n_tskin_opt
                end do ! do is=1 on surface type
              end do ! do j on priority
            end if ! end if ii > 0
            ! Set tskin_ind to invalid for non-required entries
            do iprio = n_prio, 2, -1
              where(rs%tskin_ind(iprio,:,:) == rs%tskin_ind(iprio-1,:,:)) rs%tskin_ind(iprio,:,:) = -1
            end do

            ! Print resulting tskin retrieval rules for instrument
            if (.not.any(instr_print(1:n_instr_print) == rs% instr(instr))) then
              write(*,*)
              write(*,'(1x,"Skin temperature retrieval method for instrument ",I3.3)') rs% instr(instr)
              call print_instr_opt(rs,instr)
              n_instr_print = n_instr_print + 1
              instr_print(n_instr_print) = rs% instr(instr)
            end if

          end if ! end if rs%instr(instr)
        end do ! end do instr

        rs%n_tskin_opt = n_tskin_opt

      end if ! if rs%grid
    end do ! end do i rad_set

    call destruct(tskin_set)

  contains

    subroutine find_instr(ii)
      integer, intent(out) :: ii
      integer :: j
      ii = -1
      do j = 1, n_instr_nml
        if (rs%instr(instr) == instr_nml(j)) then
          ii = j
          exit
        end if
      end do
    end subroutine find_instr

  end subroutine read_nml_tskin


  subroutine print_instr_opt(rs,instr)
    type(t_rad_set),   intent(in)  :: rs
    integer,           intent(in)  :: instr
    type(t_tskin_opt), pointer     :: tso => null()
    integer,           parameter   :: ls = 8
    integer,           parameter   :: ds = ls + 1
    character(len=400)             :: msg
    character(len=ls)              :: str = ''
    integer                        :: cdyn, n_prio, ityp
    integer                        :: n0, n1
    integer                        :: j,k,ic,its,iprio

    n_prio = size(rs%tskin_ind,1)
    msg = 'chan   |'
    do j = 1, n_styp
      msg = trim(msg)//rts_name(j-1) ; k = ls+j*n_prio*ds ; msg(k:k) = '|'
    end do

    write(*,'(3x,A)') trim(msg)
    if (n_prio > 1) then
      msg = '       |'
      do j = 1, n_styp
        msg = trim(msg)//repeat(' ',ls)//'|'
        do k = 2, n_prio
          write(msg(len_trim(msg)+1:),'("fallbck",I1)') k-1
          msg = trim(msg)//repeat(' ',8-ls)//'|'
        end do
      end do
      write(*,'(3x,A)') trim(msg)
    end if

    n0 = rs%o_ch_i(instr) + 1
    n1 = rs%o_ch_i(instr) + rs%n_ch_i(instr)

    do ic = n0, n1
    ityp = instr_type(rs,instr,ic)
      msg = ''
      cdyn = 0
      do j = 0, n_styp-1
        do iprio = 1, n_prio
          its = rs%tskin_ind(iprio, j, ic)
          if (its > 0) then
            tso => rs%tskin_opt(its)
            if  (tso%mode == MODE_TSKIN_DYNRET ) then
              if (any(rs%chan(ic) == tso%cdyn )) cdyn = rs%chan(ic)
            end if
            str = c_tskin(tso%mode,ityp,cdyn) !rs%chan(ic)
          else
            str = repeat('-',ls)
          end if
          msg = trim(msg)//str//'|'
        end do
      end do
      write(*,'(3x,I7.7,"|", A)') rs%chan(ic), trim(msg)
    end do


  end subroutine print_instr_opt


  function c_tskin(mode, ityp, cdyn) result(c)
  character(len=8)      :: c
    integer, intent(in) :: mode
    integer, intent(in) :: ityp
    integer, intent(in) :: cdyn

    if (btest(ityp,ITYP_IR)) then
      select case(mode)
      case(MODE_TSKIN_DEFAULT)
        c = 'MODELTS'
      case(MODE_TSKIN_DYNRET)
        c = 'DR'
        if (cdyn > 0) write(c(3:8),'(A)') "--USED"
      case default
        c = '??????'
      end select
    end if
  end function c_tskin

  !----------------------------------------------------------------------------
  ! computes tskin and sets it in ttovs and var, optionally stores it in obs%par
  !----------------------------------------------------------------------------
  subroutine set_tskin(spot, ttovs, rtp, obs, stat)
    type (t_spot),       intent(inout)           :: spot    ! spot
    type (t_tovs),       intent(inout), target   :: ttovs   ! t_tovs variable
    type (t_rttov_prof), intent(inout)           :: rtp     ! variables passed to rttov
    type (t_obs),        intent(inout)           :: obs     ! observations, needed for GRODY
    integer,             intent(out)             :: stat    ! error status

    character(len=*), parameter   :: proc = 'set_tskin'
    character(len=300)            :: msg
    type(t_rad_set),  pointer     :: rs
    type(t_tskin_opt), pointer    :: tso
    type(t_tovs_instr)            :: ti
    integer                       :: i_tskin(m_chan)
    integer                       :: i, j, i1, in, n_chan, iat, ityp, styp, ic
    integer                       :: chan, ii,  mod
    integer                       :: i_prio, m_prio
    integer                       :: ns, is
    logical                       :: msk_dyn(1:n_styp)
    integer,  pointer             :: fl(:)
    integer                       :: i_rs
    logical                       :: ld
    logical                       :: l_rawbt = .true.

    stat = 0
    i_tskin(:) = -1
    call init_stat

    ttovs%ts  = spot%ts_bg

    if (iand(ttovs%init, TTOVS_CI) == 0) call finish(proc, 'ttovs%ci not available')
    call get_tovs_rs(ttovs, rs=rs, ti=ti, i_rs=i_rs)

    if (.not. associated(rs%tskin_opt)) then
      call store(obs, spot, ttovs, tovs_io=TTOVS_BASE)
      return
    end if
    ! TODO: does not work for GEO:
    ! if (.not.btest(rs%gopts%opt_vars, OPTV_CLD_FRC)) call finish(proc,&
    !      'cloud fraction not available.')

    fl => ttovs%flag

    ld = ldeb(spot)

    if (ld) write(usd,*) dpref,proc,' satid/instr',spot%hd%satid, spot%hd%grid_id
    if (ld) write(usd,*) dpref,proc,' cld_cov',ttovs%cloud_imag,rs%max_tskin_cld_cov
    !if cloud cover is above threshold or not available, do not retrieve ts
    if ( ttovs%cloud_imag > rs%max_tskin_cld_cov .or. ttovs%cloud_imag < 0.) return

    n_chan = ttovs%nchan
    ns     = ttovs%ns
    m_prio = size(rs%tskin_ind,1)

    ! surface type check
    do i = 1, ns
      styp = ttovs%rt_stype(i)
      if (styp < lbound(rs%tskin_ind,2) .or. styp > ubound(rs%tskin_ind,2)) then
        write(0,*) proc//': sat/grid=',rs%satid,rs%grid
        write(0,*) proc//': styp=',styp
        write(0,*) proc//': bounds of tskin_ind=',lbound(rs%tskin_ind,2),ubound(rs%tskin_ind,2)
        call finish(proc, 'invalid surface type')
      end if
    end do

    ! ld = .true. ; write(dpref,*) spot%hd%id
    if (ld) write(usd,*) dpref,proc,' start',ns,ttovs%rt_stype(1:ns),ttovs%w_styp(1:ns)
    if (ld) write(usd,*) dpref,proc,' loop_len',ti%n_instr,ns,m_prio
    ! !-----------------------------------
    ! ! loop over instruments, select IR
    ! !-----------------------------------
    instr_loop: do i = 1, ti% n_instr                 ! do i (loop over instruments)
      ii = ti%ii(i)
      i1 = ti% o_ch_i(i) + 1
      in = ti% o_ch_i(i) + ti% n_ch_i(i)
      mask_pr(  1:n_chan) = .false.
      mask_dr(  1:n_chan) = .false.
      if (ld) write(usd,*) dpref,proc,' instr_loop',i,ii,i1,in

      styp_loop: do is = 1, ns
        if (ld) write(usd,*) dpref,proc,' styp_loop',is
        styp = ttovs%rt_stype(is)

        if (ld) write(usd,*) dpref,proc,' styp',styp,ttovs%w_styp(is),is,&
        ' stlsf', spot%stlsf, 'e_bg_ts', ttovs%e_bg_ts

        prio_loop: do i_prio = 1, m_prio

          i_tskin(i1:in) = rs%tskin_ind(i_prio, styp, ttovs%ci(i1:in))
          if (any(i_tskin(i1:in) > ubound(rs%tskin_opt,1))) call finish(proc,'invalid i_tskin')

          if (ld) then
            do j = i1,in
              if (i_tskin(j) > 0) then
                write(usd,*) dpref,proc,' i_tskin',j,i_tskin(j),fl(j), &
                     iand(fl(j), TSKIN_SUCCESS),rs%tskin_opt(i_tskin(j))%mode
              else
                write(usd,*) dpref,proc,' i_tskin',j,i_tskin(j),fl(j), &
                     iand(fl(j), TSKIN_SUCCESS)
              end if
            end do
          end if
          ! when 'TSKIN_SUCCESS' is not set in fl and i_tskin is larger
          ! than zero, then mask_pr is set to true
          mask_pr(i1:in) = (iand(fl(i1:in), TSKIN_SUCCESS) == 0) .and. &
                           (i_tskin(i1:in) > 0)
          if (ld) write(usd,*) dpref,proc,' mask_pr',count(mask_pr(i1:in))
          if (count(mask_pr(i1:in)) <= 0) then
            exit prio_loop
          end if

          mask_dr(i1:in) = .false.

          do j = i1, in

            if (.not.mask_pr(j)) cycle ! mask_pr==F-> cycle, mask_pr==T->continue

            tso => rs%tskin_opt(i_tskin(j))
            mod = tso%mode

            if (ld) write(usd,*) dpref,proc,' chan',i,i_prio,j,mod

            if (mod == MODE_TSKIN_DYNRET) then
              ic = -2
              if (tso%n_cdyn > 0) then
                if ( any(tso%cdyn(:) == rs%chan(ttovs%ci(j)))) then
                  ic = 0
                else
                  ic = -1
                end if
              else
                call finish(proc, 'cdyn not provided')
              end if
            end if
            select case(mod)
            case(MODE_TSKIN_DYNRET)
              if (ic <= -2) then
                write(0,*) 'spot/instr/chan',spot%hd%id,ii,rs%instr(ii),j,rs%chan(ttovs%ci(j))
                if (associated(tso%cdyn)) write(0,*) 'tso%decr ',trim(tso%descr),' tso%cdyn ', tso%cdyn
                call finish(proc, 'Did not find cdyn channel in rs%chan')
              elseif(ic == -1) then
                mask_dr(j) = .false.
              else
                mask_dr(j) = .true.
              end if

            end select
          end do

          ! Retrieve skin temperature
          if (ld) write(usd,*) dpref,proc,' mask_dr',count(mask_dr(i1:in))

          if (count(mask_dr(i1:in)) > 0) then
            ! Comment by RF: please note, that in the current implementation trace gases are not taken
            !                into account in call_rttov - although the "main" RTTOV calls in
            !                process_tovs_mult@mo_tovs are using actual trace gas fields. This might
            !                degrade the quality of the Tskin retrieval.
            !                So we should use a channel without O3 dependency like e.g. IASI 1027.
            if (l_rawbt) then
              call call_rttov('t', ttovs, rtp, stat, msk=mask_dr(1:n_chan), styp=styp, ldebug=ld,&
                  obs_bt=real(obs% body(spot%o%i+1 : spot%o%i+ttovs%nchan)%o ,kind=wp), spt_hd_id=spot%hd%id)
            else
              call call_rttov('t', ttovs, rtp, stat, msk=mask_dr(1:n_chan), styp=styp, ldebug=ld,&
                  obs_bt=real(obs% body(spot%o%i+1 : spot%o%i+ttovs%nchan)%o - &
                  obs% body(spot%o%i+1 : spot%o%i+ttovs%nchan)%bc ,kind=wp), spt_hd_id=spot%hd%id)
            end if ! end if l_rawbt

            ! Nothing to distribute, only Set TF_TSKIN_* flag
            ! If retrieval is carried out successfully or not
            ! set TSKIN_* flag accordingly
            if (rtp%drtsfl) then
              chan_loop: do j = i1, in                       ! do j (loop over channels)
                if (.not.mask_pr(j)) cycle
                tso => rs%tskin_opt(i_tskin(j))
                mod = tso%mode
                if (ld) write(usd,*) dpref,proc,' tskin ',spot%hd%id,j,rtp%ssv(1,:)

                if (mod == MODE_TSKIN_DYNRET) then
                  fl(j) = ibset(fl(j), TF_TSKIN_DYNRET)
                  fl(j) = ibclr(fl(j), TF_TSKIN_FAILED)
                  ! save tskin - spot is not useful here since it changes in psas interpolation
                  ttovs% ts = rtp% ssv(1,1)
                  ! to be used for setting e_bg_ts_land
                  ttovs% l_ts_dr = .true.
                end if
                if (ld) write(usd,*) dpref,proc,' tskin fl dyret ',spot%hd%id,j,fl(j)
              end do chan_loop
            else
              chan_loop2: do j = i1, in                       ! do j (loop over channels)
                if (.not.mask_pr(j)) cycle
                tso => rs%tskin_opt(i_tskin(j))
                mod = tso%mode
                if (ld) write(usd,*) dpref,proc,' tskin ',spot%hd%id,j,rtp%ssv(1,:)

                if (mod == MODE_TSKIN_DYNRET) then
                  fl(j) = ibclr(fl(j), TF_TSKIN_DYNRET)
                  fl(j) = ibset(fl(j), TF_TSKIN_FAILED)
                  ! save tskin - spot is not useful here since it changes in psas interpolation
                  !ttovs% ts = rtp% ssv(1,1)
                end if
                if (ld) write(usd,*) dpref,proc,' tskin fl fail ',spot%hd%id,j,fl(j)
              end do chan_loop2
            end if ! if derive tskin flag is true

            call add_stat(tso%id, styp, .not.rtp%drtsfl, rtp%ssv(1,1), spot%ts_bg, mix=(ns>1))

          end if ! end if mask_dr
          if (stat > 0) then
            call finish(proc,'call_rttov("t",..) failed')
          end if
        end do prio_loop
        if (ld) write(usd,*) dpref,proc,' styp_loop done',is,ns

      end do styp_loop

    end do instr_loop

    call store(obs, spot, ttovs, tovs_io=TTOVS_FLAG+TTOVS_BASE)

  end subroutine set_tskin

  !--------------------
  ! Statistics routines
  !--------------------

  subroutine init_stat
    type(t_rad_set),    pointer :: rs
    type(t_tskin_opt),  pointer :: tso
    type(t_tskin_stat), pointer :: st
    integer :: iset, iopt, is, nch, i, j, ioff, n

    if (associated(stat)) return

    nstat = 0
    do iset = 1, n_set
      rs => rad_set(iset)
      do iopt = 1, rs%n_tskin_opt
        tso => rs%tskin_opt(iopt)
        nstat = nstat + tso%ns
      end do
    end do
    allocate(stat(nstat))

    n = 0
    do iset = 1, n_set
      rs => rad_set(iset)
      do iopt = 1, rs%n_tskin_opt
        tso => rs%tskin_opt(iopt)
        do is = 1, tso%ns
          n = n + 1
          st => stat(n)
          st%satid = rs%satid
          st%id    = tso%id
          st%instr = tso%inst
          st%styp  = tso%styp(is)
          st%descr = trim(tso%descr)
        end do
      end do
    end do

  end subroutine init_stat


  subroutine add_stat(id, styp, fail, ts_ret, ts_bg, mix)
    integer,  intent(in) :: id
    integer,  intent(in) :: styp
    logical,  intent(in) :: fail
    real(wp), intent(in) :: ts_ret
    real(wp), intent(in) :: ts_bg
    logical,  intent(in) :: mix

    type(t_tskin_stat), pointer :: st
    integer,            save    :: ii = 1
    integer                     :: i
    logical                     :: lf

    if (stat(ii)%id /= id .or. stat(ii)%styp /= styp) then
      lf = .false.
      do i = ii+1,nstat
        if (stat(i)%id == id .and. stat(i)%styp == styp) then
          ii = i
          lf = .true.
        end if
      end do
      if (.not.lf) then
        do i = 1,nstat
          if (stat(i)%id == id .and. stat(i)%styp == styp) then
            ii = i
            lf = .true.
          end if
        end do
      end if
      if (.not.lf) then
        write(0,*) 'id,styp',id,styp
        call finish('add_stat', 'Did not find tskin_stat entry')
      end if
    end if
    st => stat(ii)
    st%n = st%n + 1
    if (fail) then
      st%n_f = st%n_f + 1
    else
      st%tskin = st%tskin + ts_ret
      st%dt_mn = st%dt_mn + (ts_ret - ts_bg)
      st%dt_s2 = st%dt_s2 + (ts_ret - ts_bg)**2
    end if
    if (mix) then
      st%n_mix = st%n_mix + 1
      if (fail) st%n_mix_f = st%n_mix_f + 1
    end if

  end subroutine add_stat

  
  subroutine print_tskin_stat
    type(t_tskin_stat), pointer :: st
    character(len=300)          :: msg
    real(wp)                    :: tskin, dt_mn, dt_s2, dt_sd, fr_mix, fr_fail
    integer                     :: n_omit
    integer                     :: i, j, ic, nch
    integer                     :: mstat

    if (.not.l_pr) return ! Not really necessary, but avoids waiting for MPI below
    mstat = p_max(nstat)
    if (mstat <= 0) return

    if (nstat <= 0) call init_stat
    ! summarize over all PEs
    do i = 1, nstat
      st => stat(i)
      st%n       = p_sum(st%n      )
      st%n_f     = p_sum(st%n_f    )
      st%n_mix   = p_sum(st%n_mix  )
      st%n_mix_f = p_sum(st%n_mix_f)
      st%tskin   = p_sum(st%tskin  )
      st%dt_mn   = p_sum(st%dt_mn  )
      st%dt_s2   = p_sum(st%dt_s2  )
    end do

    if (dace%lpio) then
      write(*,*)
      write(*,*) 'tskin retrieval statistics:'
      write(*,*) '(*_m: #mixed surfaces, *_f: #method failed)'
      do i = 1, nstat
        st => stat(i)
        if (st%n > 0) then
          write(*,*)
          write(*,'(3x,A," satid=",I3," instr=",I3," styp=",I1)') &
               trim(st%descr), st%satid, st%instr, st%styp
          write(*,'(5x,1x,"       n",1x,"     n_f",1x,"   n/n_f",1x,"mn(ts)",1x,"mn(dts)",1x,"sd(dts)",1x,&
               &"     n_m","    n_m_f",1x,"n_m/n_m_f",1x,"n_m_f/n_f")')
          n_omit   = 0
          tskin = -1._wp
          if (st%n > st%n_f) then
            tskin = st%tskin / (st%n-st%n_f)
            dt_mn = st%dt_mn / (st%n-st%n_f)
            dt_s2 = st%dt_s2 / (st%n-st%n_f)
            dt_sd = dt_s2 - dt_mn**2
            if (dt_sd >= 0._wp) then
              dt_sd = sqrt(dt_sd)
            else
              dt_sd = 0._wp
            end if
          end if
          fr_mix= 0._wp
          if (st%n_mix > 0) fr_mix=(100._wp*st%n_mix_f)/st%n_mix
          fr_fail= 0._wp
          if (st%n_f > 0) fr_fail=(100._wp*st%n_mix_f)/st%n_f
          write(msg,'(2(1x,I8),1x,F7.3,"%",1x,F6.2,2(1x,F7.3),2(1x,I8),2(2x,F7.3,"%"))') &
               st%n,st%n_f,(100._wp*st%n_f)/real(st%n,wp),tskin,dt_mn,dt_sd,&
               st%n_mix,st%n_mix_f,fr_mix,fr_fail
          write(*,'(5x,A)') trim(msg)
        end if
      end do
    end if
    
    deallocate(stat)
    nstat = 0
    l_pr  = .false.

  end subroutine print_tskin_stat


END MODULE mo_tskin
