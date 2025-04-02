!
!+ MW emissivity models
!
MODULE mo_emis
!
! Description:
!   MW emissivity models: sea_model, atlases, dynamic retrieval
!
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
  use mo_kind,              only: wp, sp                    ! working precision kind parameter

  use mo_exception,         only: finish                    ! abort in case of error

  use mo_mpi_dace,          only: dace,                    &! MPI group info
                                  p_sum,                   &!
                                  p_max

  use mo_namelist,          only: position_nml,            &! position namelist
                                  nnml,                    &! namelist Fortran unit number
                                  POSITIONED                ! ok    code from position_nml
  use mo_instrid,           only: mw_instr, ir_instr, vis_instr

  use mo_p_output,          only: add_line

  !-----------------------
  ! access to observations
  !-----------------------
  use mo_rad,               only: m_chan,                  &! max number of channels
                                  rad_set,                 &! radiance meta data
                                  n_set,                   &! number of valid rad_sets
                                  set_indx,                &! get sat/instr/chan set index
                                  print_rad_set,           &! print rad_set for debugging
                                  t_rad_set,               &! rad_set type
                                  warning,                 &!
                                  n_styp,                  &!
                                  t_emis_opt,              &!
                                  assign_emis_opt,         &!
                                  instr_type,              &!
                                  ITYP_MW,ITYP_IR,ITYP_VIS,&
                                  construct, destruct

  use mo_t_obs,             only: t_obs,                   &! observation derived type
                                  t_spot,                  &! spot derived type
                                  ldeb,usd,dpref            ! debug selected spot
  use mo_fdbk_tables,       only: OT_RAD,                  &! Radiances report type ID
                                  TF_EMIS,                 &
                                  TF_EMIS_SEA_MOD,         &
                                  TF_EMIS_DYNRET,          &
                                  TF_EMIS_GRODY,           &
                                  TF_EMIS_TLSM,            &
                                  TF_EMIS_CNRM,            &
                                  TF_EMIS_CAMEL07,         &
                                  TF_EMIS_CAMELCL,         &
                                  TF_EMIS_UWIR,            &
                                  TF_EMIS_FAILED,          &
                                  TF_REFL_BRDF

  !-------------------------
  ! rttov related modules
  !-------------------------
  use mo_t_tovs,            only: t_tovs,                  &! observation operator specific type
                                  t_tovs_instr,            &! information on instruments in t_tovs
                                  get_tovs_rs,             &! get rad_set and instrument info for t_tovs
                                  store,                   &
                                  TTOVS_FLAG, TTOVS_CI,    &
                                  TTOVS_SPEC

  use mo_rttov,             only: t_rttov_prof,            &! arguments to rttov
                                  call_rttov                ! call to rttovs interface subroutines

  use mo_t_use,             only: STAT_REJECTED

  !-------------------------
  ! cloud detection modules
  !-------------------------
  use mo_cloud_indices,     only: mw_emiss,                & ! Subroutine that determines Grody emissivity
                                  sfchl_amsua,             & ! surface channels AMSU-A
                                  sfchl_atms                 ! surface channels ATMS

  use mo_rtifc,             only: rts_name,                &!
                                  rts_land,rts_sea,rts_ice  !

  implicit none

  !================
  ! public entities
  !================
  private

  public :: set_emis
  public :: read_nml_emis
  public :: update_emis_opt
  public :: print_emis_stat
  public :: dynret_w_pref
  public :: dynret_avg
  public :: m_emis
  public :: MODE_ATLAS
  public :: ATLS_TLSM
  public :: ATLS_CNRM
  public :: ATLS_UWIR
  public :: ATLS_CML
  public :: ATLS_CML_CL
  public :: ATLS_BRDF
  public :: EMIS_CALC
  public :: EMIS_SUCCESS


  !=================
  ! Parameters
  !=================
  integer, parameter :: m_emis           = 50      ! Max number of mw_emi_set to be set via namelist

  ! Codes for the different emissivity models
  integer, parameter :: MODE_SEA_MODEL   = 0       ! Code for sea model (e.g. FASTEM)
  integer, parameter :: MODE_DYNRET      = 1       ! Code for dynamic retrieval mode
  integer, parameter :: MODE_GRODY       = 2       ! Code for Grody mode
  integer, parameter :: MODE_ATLAS       = 3       ! Code for atlas mode
  integer, parameter :: ATLS_TLSM        = 1       ! Code for TELSEM atlas (MW)
  integer, parameter :: ATLS_CNRM        = 2       ! Code for CNRM atlas (MW)
  integer, parameter :: ATLS_UWIR        = 1       ! Code for UWIR atlas (MW)
  integer, parameter :: ATLS_CML         = 2       ! Code for CAMEL2007 atlas (IR)
  integer, parameter :: ATLS_CML_CL      = 3       ! Code for CAMEL clim. atlas (IR)
  integer, parameter :: ATLS_BRDF        = 4       ! Code for BRDF atlas (VIS)

  integer, parameter :: EMIS_SUCCESS = 2**TF_EMIS(0) + &
                                       2**TF_EMIS(1) + &
                                       2**TF_EMIS(2) + &
                                       2**TF_EMIS(3) + &
                                       2**TF_EMIS(4) + &
                                       2**TF_EMIS(5) + &
                                       2**TF_EMIS(6)
  integer, parameter :: EMIS_CALC    = EMIS_SUCCESS + &
                                       2**TF_EMIS(7)


  !=================
  ! Module variables
  !=================
  integer                    :: n_nml = 0
  integer, save              :: id    = 0
  type(t_emis_opt),  target  :: emis_set(m_emis)
  logical                    :: mask_pr(    m_chan) ! used in set_emis
  logical                    :: mask_dr(    m_chan) ! used in set_emis
  logical                    :: mask_at(1:3,m_chan) ! used in set_emis
  logical                    :: mask_br(    m_chan) ! used in set_emis
  logical                    :: mask_ss(    m_chan)

  ! Statistics
  type t_emis_stat
    integer                  :: id       =  -1     ! same as t_emis_opt ID
    integer                  :: satid    =  -1
    integer                  :: instr    =  -1
    integer                  :: prio     =  -1
    integer                  :: styp     =  -1
    character(len=120)       :: descr    =  ''
    integer                  :: n_chan   =  0
    integer,         pointer :: chan (:) => NULL()
    integer,         pointer :: n    (:) => NULL()
    integer,         pointer :: nfail(:) => NULL()
    integer,         pointer :: n_mix(:) => NULL()
    integer,         pointer :: n_mix_fail(:) => NULL()
    real(wp),        pointer :: emiss(:) => NULL()
  end type t_emis_stat
  type(t_emis_stat), pointer :: stat(:)  => NULL()
  integer                    :: nstat    =  0
  type t_styp_stat
    integer                  :: ist      = 0  ! Bits of surface
    integer                  :: ns       = 0  ! number od surface types
    integer                  :: n        = 0  ! number of spots with this styp combination
    real(wp)                 :: w(0:n_styp-1)= 0._wp
  end type t_styp_stat
  type(t_styp_stat), pointer :: stst(:,:)=> null()
  integer                    :: nstst    = 0
  logical                    :: l_pr     =  .true.

  ! namelist settings
  real(wp)                   :: dynret_w_pref = 0.1_wp
  integer                    :: dynret_avg    = 0

contains

  !----------------------------------------------------------------------------
  ! Reads /EMIS/ namelist into emis_set
  !----------------------------------------------------------------------------
  subroutine read_nml_emis_sub
    character(len=17), parameter :: proc = 'read_nml_emis_sub'
    type(t_emis_opt),  pointer   :: eo   => NULL()
    integer                      :: i_read
    logical                      :: first
    integer                      :: ierr_pos
    integer                      :: ityp, n_chan, k
    character(len=300)           :: msg
    character(len=80)            :: hint
    integer                      :: sfchl(4)
    integer                      :: i, ns
    ! Namelist content
    integer                      :: inst
    integer                      :: prio
    integer                      :: styp(n_styp)
    integer                      :: mode
    integer                      :: chan (m_chan)
    integer                      :: atls (m_chan)
    integer                      :: atlas
    integer                      :: cdyn (m_chan)
    logical                      :: angcorr
    real(wp)                     :: max_dst
    character(len=120)           :: descr

    ! The MW_EMIS namelist are obsolete. Since they were operational for some time,
    ! they are kept for backwards compatibility

    namelist /MW_EMIS/ inst, prio, styp, mode, chan, cdyn, atls
    namelist /EMIS/    inst, prio, styp, mode, chan, cdyn, atlas, angcorr, max_dst, descr

    write(*,*) 'READ EMIS namelists:'
    first  = .true.
    i_read = 0
    n_nml  = 1
    loop_emis: do
      call position_nml ('EMIS' ,lrewind=first ,status=ierr_pos)
      if (ierr_pos == POSITIONED) then
        i_read = i_read + 1
        call init
        write(hint, '("EMIS namelist number ",I2,1x)') i_read
        read (nnml ,nml=EMIS)
        call process_nml
      else
        exit loop_emis
      end if
      first = .false.
    end do loop_emis

    first  = .true.
    i_read = 0
    loop_mw_emis: do
      call position_nml ('MW_EMIS' ,lrewind=first ,status=ierr_pos)
      if (ierr_pos == POSITIONED) then
        i_read = i_read + 1
        call init
        write(hint, '("MW_EMIS namelist number ",I2,1x)') i_read
        read (nnml ,nml=MW_EMIS)
        call process_nml
      else
        exit loop_mw_emis
      end if
      first = .false.
    end do loop_mw_emis

    n_nml = n_nml - 1

  contains

    subroutine init
      inst         = -1
      prio         =  0
      styp(:)      = -1
      mode         =  0
      chan(:)      =  0
      atls(:)      = -1
      atlas        = -1
      cdyn(:)      = -1
      angcorr      = .true.
      max_dst      = 0._wp
      descr        = ''
    end subroutine init


    subroutine process_nml

      if (inst < 0) call finish(proc, 'Missing instrument '//trim(hint))
      if (atlas < 0) atlas = atls(1)
      i = maxval(atls)
      if (any(atls > 0 .and. atls /= i)) call finish(proc, 'Multiple atlases in EMIS not supported')
      if (descr == '') then
        ityp = 0
        if (mw_instr(inst)) ityp = ibset(ityp, ITYP_MW)
        if (ir_instr(inst)) ityp = ibset(ityp, ITYP_IR)
        descr = trim(hint)//' '//c_emis(mode, ityp, atlas, 0)
      end if
      hint = trim(hint)//' "'//trim(descr)//'"'
      n_chan = count(chan > 0)

      ! Consistency checks
      if (n_chan > 0) then
        if (.not. all(chan(1:n_chan) > 0)) call finish(proc, trim(hint)//': invalid chan')
      end if
      if (any(styp(:) < -1 .or. styp(:) > n_styp)) then
        write(msg,*) trim(hint)//': invalid styp=',styp(:)
        call finish('read_nml_emis',trim(msg))
      end if
      ns = count(styp >= 0)
      if (ns == 0) then
        styp = (/ 0, 1, 2/)
        ns   = 3
      end if
      styp(1:ns) = pack(styp, mask=(styp >= 0))
      select case (mode)
      ! case (MODE_SEA_MODEL)
      !   cycle loop_nml
      case (MODE_ATLAS)
        if (atlas < 1 .or. atlas > 3) then
          call finish(proc, trim(hint)//' has mode 1 but no valid atlas identifier')
        else if (atlas == 2 .and. mw_instr(inst) .and. &
             .not. any(inst == (/ 3, 4, 10, 15, 19 /)) ) then
          call finish(proc, trim(hint)//' uses CNRM atlas for an instrument different &
               &than AMSU-A/-B, ATMS, MHS or SSMI/S')
        end if
      case (MODE_GRODY)
        if ( inst /= 3 .and. inst /= 19 ) call finish(proc, trim(hint)//&
             ' uses Grody method for an instrument different than AMSU-A or ATMS')
        ! Check window channels (needed) present
        if ( inst == 3  ) sfchl = sfchl_amsua
        if ( inst == 19 ) sfchl = sfchl_atms
        do k = 1, size(sfchl)
          if (.not.any(chan(1:n_chan) == sfchl(k) )) call finish(proc, trim(hint)//&
               ' uses Grody method but surface channel is not available.')
        end do
      case (MODE_DYNRET)
        ! Check all channels in cdyn are in chan
        do k = 1, n_chan
          if (.not.any( chan(1:n_chan) == cdyn(k)) .and. cdyn(k) /= 0) call finish(proc, &
               trim(hint)//' uses dynamical retrievals with a missmatch between chan(:) and cdyn(:)')
        end do
      case default
        call finish('read_nml_emis', trim(hint)//': unknown mode')
      end select

      ! Fill corresponding emis_set element
      if (n_nml > size(emis_set)) call  finish(proc, 'Too many EMIS namelists. Recompile with &
           &increased m_emis')
      eo => emis_set(n_nml)
      eo%nml          = n_nml
      eo%inst         = inst
      eo%prio         = prio
      eo%ns           = ns
      eo%styp         = styp
      eo%descr        = trim(descr)
      eo%mode         = mode
      eo%atls         = atlas
      eo%n_chan       = n_chan
      eo%angcorr      = angcorr
      eo%max_dst      = max_dst
      allocate(eo%chan (n_chan))
      if (mode == MODE_DYNRET) allocate(eo%cdyn(n_chan))
      do k= 1, n_chan
        eo%chan(k) = chan(k)
        if (mode == MODE_DYNRET) eo%cdyn(k) = cdyn(k)
      end do
      n_nml = n_nml + 1
      ! printout
      write(*,*)
      write(*,'(3x,a)') trim(hint)
      write(*,'(5x,a,1X,a)')                     ' descr        = ', trim(eo%descr)
      write(*,'(5x,a,1X,i2)')                    ' inst         = ', eo%inst
      write(*,'(5x,a,1X,i2)')                    ' prio         = ', eo%prio
      write(*,'(5x,a,10(1X,i2))')                ' styp         = ', eo%styp(1:eo%ns)
      write(*,'(5x,a,1X,i2)')                    ' mode         = ', eo%mode
      write(*,'(5x,a,1X,i2)')                    ' atlas        = ', eo%atls
      if (eo%n_chan > 0) then
        write(*,'(5x,a,*(8(1x,i5,:),/,20x))')    ' chan         = ', eo%chan(1:eo%n_chan)
        if (eo%mode == MODE_DYNRET) then
          write(*,'(5x,a,*(8(1x,i5,:),/,20x))')  ' cdyn         = ', eo%cdyn(1:eo%n_chan)
        end if
      end if
      write(*,'(5x,a,l1)')                       ' angcorr      = ', eo%angcorr
      write(*,'(5x,a,f8.3," km")')               ' max_dst      = ', eo%max_dst
    end subroutine process_nml

  end subroutine read_nml_emis_sub


  subroutine read_nml_emis
    character(len=13), parameter   :: proc = 'read_nml_emis'
    type(t_rad_set),   pointer     :: rs
    type(t_emis_opt),  pointer     :: eo => null()
    type(t_emis_opt)               :: aux
    integer,           allocatable :: instr_nml(:), n_prio_instr(:), n_set_instr(:), i_vis(:), i0(:)
    integer                        :: i, ii, j, k, instr, iprio, is, ityp
    integer                        :: n0, n1
    integer                        :: n_emis_opt
    integer                        :: n_instr_nml
    integer                        :: n_prio, n_vis
    integer                        :: styp
    logical                        :: l_swap
    ! printout
    integer                        :: n_instr_print
    integer                        :: instr_print(100)

    ! read /EMIS/ namelists
    call read_nml_emis_sub

    !-----------------
    ! Prepare emis_set
    !-----------------
    ! Sort for inst and prio
    do i = n_nml-1,1,-1
      do k = 1,i
        l_swap = .false.
        if (emis_set(k)%inst >  emis_set(k+1)%inst) then
          l_swap = .true.
        elseif (emis_set(k)%inst == emis_set(k+1)%inst) then
          if (emis_set(k)%prio > emis_set(k+1)%prio) then
            l_swap = .true.
          end if
        end if
        if (l_swap) then
          aux           = emis_set(k)
          emis_set(k)   = emis_set(k+1)
          emis_set(k+1) = aux
        end if
      end do
    end do

    ! count number of different sets and priorities per instrument
    allocate(instr_nml(n_nml), n_prio_instr(n_nml), n_set_instr(n_nml), i0(n_nml))
    n_instr_nml = 0
    do i = 1, n_nml
      ii = -1
      if (i > 1) then
        if (emis_set(i)%inst == emis_set(i-1)%inst) ii = n_instr_nml
      end if
      if (ii < 0) then
        n_instr_nml = n_instr_nml + 1
        ii          = n_instr_nml
        instr_nml   (ii) = emis_set(i)%inst
        i0          (ii) = i
        n_set_instr (ii) = 0
        n_prio_instr(ii) = 0
      end if
      n_set_instr(ii)  = n_set_instr(ii) + 1
      if (all(emis_set(i)%prio /= emis_set(i0(ii):i-1)%prio)) n_prio_instr(ii) = n_prio_instr(ii) + 1
    end do

    !----------------------------------
    ! cycle over satellites/instruments
    !----------------------------------
    n_instr_print = 0
    id = 0
    do i = 1, size(rad_set)
      rs => rad_set(i)
      ! Derive dimensions, allocate
      n_prio = 0
      n_emis_opt = 0
      do instr = 1, rs%n_instr
        ityp = instr_type(rs, instr)
        ! Add default
        if (btest(ityp, ITYP_MW) .or. btest(ityp, ITYP_IR)) n_emis_opt = n_emis_opt + 1
        if (btest(ityp, ITYP_VIS)                         ) n_emis_opt = n_emis_opt + 1
        !  (for ITYP_VIS we reserve space for the BRDF atlas)
        call find_instr(ii)
        if (ii > 0) then
          n_emis_opt = n_emis_opt + n_set_instr(ii)
          n_prio = max(n_prio, n_prio_instr(ii)+1)
        end if
      end do
      n_prio = max(n_prio, 1)
      allocate(rs%emis_opt(n_emis_opt), rs%emis_ind(n_prio, 0:n_styp-1, rs%n_chan))

      n_emis_opt = 0
      do instr = 1, rs%n_instr
        n0 = rs%o_ch_i(instr) + 1
        n1 = rs%o_ch_i(instr) + rs%n_ch_i(instr)
        ityp = instr_type(rs, instr)
        ! Add defaults
        if (btest(ityp, ITYP_MW) .or. btest(ityp, ITYP_IR)) then
          n_emis_opt = n_emis_opt + 1
          eo => rs%emis_opt(n_emis_opt)
          call construct(eo)
          eo%inst  = rs%instr(instr)
          eo%mode  = MODE_SEA_MODEL
          eo%ns    = 3
          eo%styp  = (/ 0, 1, 2/)
          eo%descr = 'Default method: '//c_emis(eo%mode, ityp, 0, 0)
          eo%id    = id ; id = id + 1
          rs%emis_ind(:,:,n0:n1) = n_emis_opt
        end if
        if (btest(ityp, ITYP_VIS) .and. associated(rs%iflag)) then
          ! default for vis chans: BRDF atlas
          ! This would work only if iflag is available, i.e. after reading the data.
          ! Since we have to call this routine before reading the data, iflag is evaluated
          ! in update_emis_opt
          allocate(i_vis(rs%n_ch_i(instr)))
          where (rs%iflag(n0:n1) == 1)
            i_vis = (/ (j, j=n0,n1) /)
          elsewhere
            i_vis = 0
          end where
          n_vis = count(i_vis > 0)
          if (n_vis > 0) then
            n_emis_opt = n_emis_opt + 1
            eo => rs%emis_opt(n_emis_opt)
            call construct(eo)
            eo%inst  = rs%instr(instr)
            eo%mode  = MODE_ATLAS
            eo%atls  = ATLS_BRDF
            eo%ns    = 3
            eo%styp  = (/ 0, 1, 2/)
            eo%descr = 'Default method: BRDF atlas'
            eo%id    = id ; id = id + 1
            i_vis(1:n_vis) = pack(i_vis, mask=(i_vis > 0))
            rs%emis_ind(:,:,i_vis(1:n_vis)) = n_emis_opt
          end if
          deallocate(i_vis)
        end if
        ! Add options from namelist
        call find_instr(ii)
        if (ii > 0) then
          do j = i0(ii), i0(ii) + n_set_instr(ii) - 1
            n_emis_opt = n_emis_opt + 1
            eo => rs%emis_opt(n_emis_opt)
            call assign_emis_opt(eo, emis_set(j))
            eo%id = id ; id = id + 1
            iprio = emis_set(j)%prio + 1
            do is = 1, eo%ns
              styp  = eo%styp(is)
              if (eo%n_chan > 0) then
                ! selected channels
                do k = 1, eo%n_chan
                  where(rs%chan(n0:n1) == eo%chan(k))
                    rs%emis_ind(iprio,styp,n0:n1) = n_emis_opt
                  end where
                end do
              else
                ! all channels
                rs%emis_ind(iprio,styp,n0:n1) = n_emis_opt
              end if
            end do
          end do
        end if
        ! Set emis_ind to invalid for non-required entries
        do iprio = n_prio, 2, -1
          where(rs%emis_ind(iprio,:,:) == rs%emis_ind(iprio-1,:,:)) rs%emis_ind(iprio,:,:) = -1
        end do

        ! Print resulting emissivity retrieval rules for instrument
        if (.not.any(instr_print(1:n_instr_print) == rs% instr(instr))) then
          write(*,*)
          write(*,'(1x,"Emissivity calculation methods for instrument ",I3.3)') rs% instr(instr)
          call print_instr_opt(rs,instr)
          n_instr_print = n_instr_print + 1
          instr_print(n_instr_print) = rs% instr(instr)
        end if
      end do
      rs%n_emis_opt = n_emis_opt
    end do

    call destruct(emis_set)

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

  end subroutine read_nml_emis


  subroutine update_emis_opt
    ! Purpose: we call read_nml_emis before the data is read, i.e. in read_nml_emis it is
    ! not known, which channel is a VIS channel. To determine the VIS channels on the basis
    ! of the channel numbers only will fail, because of different channel numbering conventions
    ! in different datasets. So, unfortunately, we have to set the defaults for VIS channels here.
    ! It's really not elegant to have an extra routine for this purpose!
    ! Unlike read_nml_emis, this is to be called on ALL PEs
    character(len=13), parameter   :: proc = 'read_nml_emis'
    type(t_rad_set),   pointer     :: rs
    type(t_emis_opt),  pointer     :: eo => null()
    integer                        :: i, j, k, instr, iopt, ii, ityp
    integer                        :: nprio, n0, n1, n_vis
    integer,           allocatable :: i_vis(:)
    integer                        :: n_instr_print
    integer                        :: instr_print(100)

    n_instr_print = 0
    do i = 1, n_set
      rs => rad_set(i)
      do instr = 1, rs%n_instr
        n0 = rs%o_ch_i(instr) + 1
        n1 = rs%o_ch_i(instr) + rs%n_ch_i(instr)
        if (associated(rs%iflag)) then
          ityp = instr_type(rs, instr)
          if (btest(ityp, ITYP_VIS)) then
            ! default for vis chans: BRDF atlas
            allocate(i_vis(rs%n_ch_i(instr)))
            where (rs%iflag(n0:n1) == 1)
              i_vis = (/ (j, j=n0,n1) /)
            elsewhere
              i_vis = 0
            end where
            n_vis = count(i_vis > 0)
            if (n_vis > 0) then
              iopt = -1
              do j = 1, rs%n_emis_opt
                eo => rs%emis_opt(j)
                if ( eo%inst  == rs%instr(instr) .and. &
                     eo%mode  == MODE_ATLAS      .and. &
                     eo%atls  == ATLS_BRDF ) then
                  iopt = j
                  exit
                end if
              end do
              if (iopt <= 0) then
                if (rs%n_emis_opt >= size(rs%emis_opt)) &
                     call finish(proc, 'No space in emis_opt for BRDF default')
                rs%n_emis_opt = rs%n_emis_opt + 1
                iopt = rs%n_emis_opt
                eo => rs%emis_opt(iopt)
                call construct(eo)
                eo%inst  = rs%instr(instr)
                eo%mode  = MODE_ATLAS
                eo%atls  = ATLS_BRDF
                eo%ns    = 3
                eo%styp  = (/ 0, 1, 2/)
                eo%descr = 'Default method: BRDF atlas'
                eo%id    = id ; id = id + 1
              end if
              i_vis(1:n_vis) = pack(i_vis, mask=(i_vis > 0))
              nprio = size(rs%emis_ind,1)
              do j = 1, n_vis
                do k = 0, n_styp-1
                  ii = rs%emis_ind(1,k,i_vis(j))
                  if (rs%emis_opt(ii)%nml < 0) then ! default option, not modified by nml
                    rs%emis_ind(2:nprio,k,i_vis(j)) = rs%emis_ind(1:nprio-1,k,i_vis(j))
                    rs%emis_ind(1,k,i_vis(j))= iopt
                  end if
                end do
              end do
              ! Print resulting emissivity retrieval rules for instrument
              if (dace%lpio.and..not.any(instr_print(1:n_instr_print) == rs% instr(instr))) then
                write(*,*)
                write(*,'(1x,"Updated emissivity calculation methods for instrument ",I3.3)') rs% instr(instr)
                call print_instr_opt(rs,instr)
                n_instr_print = n_instr_print + 1
                instr_print(n_instr_print) = rs% instr(instr)
              end if
            end if
            deallocate(i_vis)
          end if
        end if
      end do
    end do
  end subroutine update_emis_opt


  subroutine print_instr_opt(rs,instr)
    type(t_rad_set),   intent(in)  :: rs
    integer,           intent(in)  :: instr
    type(t_emis_opt),  pointer     :: eo => null()
    integer,           parameter   :: ls = 6
    integer,           parameter   :: ds = ls + 1
    character(len=300)             :: msg, msg_prev
    character(len=ls)              :: str = ''
    integer                        :: cdyn, n_omit, n_prio, ityp
    integer                        :: n0, n1
    integer                        :: j,k,ic,ie,iprio

    n_prio = size(rs%emis_ind,1)

    msg = 'chan |'
    do j = 1, n_styp
      msg = trim(msg)//rts_name(j-1) ; k = ls+j*n_prio*ds ; msg(k:k) = '|'
    end do
      ! msg = trim(msg)//'land'   ; k = ls+1*n_prio*ds ; msg(k:k) = '|'
      ! msg = trim(msg)//'sea'    ; k = ls+2*n_prio*ds ; msg(k:k) = '|'
      ! msg = trim(msg)//'seaice' ; k = ls+3*n_prio*ds ; msg(k:k) = '|'
    write(*,'(3x,A)') trim(msg)
    if (n_prio > 1) then
      msg = '     |'
      do j = 1, n_styp
        msg = trim(msg)//repeat(' ',ls)//'|'
        do k = 2, n_prio
          write(msg(len_trim(msg)+1:),'("fallb",I1)') k-1
          msg = trim(msg)//repeat(' ',6-ls)//'|'
        end do
      end do
      write(*,'(3x,A)') trim(msg)
    end if
    msg_prev = ''
    n_omit = 0
    n0 = rs%o_ch_i(instr) + 1
    n1 = rs%o_ch_i(instr) + rs%n_ch_i(instr)
    do ic = n0, n1
      ityp = instr_type(rs,instr,ic)
      msg = ''
      do j = 0, n_styp-1
        do iprio = 1, n_prio
          ie = rs%emis_ind(iprio, j, ic)
          if (ie > 0) then
            eo => rs%emis_opt(ie)
            if (eo%mode == MODE_ATLAS .or. eo%mode == MODE_DYNRET) then
              do k = 1, eo%n_chan
                if (eo%chan(k) == rs%chan(ic)) then
                  if (associated(eo%cdyn)) cdyn = eo%cdyn(k)
                  exit
                end if
              end do
            end if
            str = c_emis(eo%mode,ityp,eo%atls,cdyn)
          else
            str = repeat('-',ls)
          end if
          msg = trim(msg)//str//'|'
        end do
      end do
      if (msg /= msg_prev .or. ic == n1) then
        write(*,'(3x,I5.5,"|",A)') rs%chan(ic),trim(msg)
        n_omit = 0
      else
        n_omit = n_omit + 1
        if (n_omit == 1) write(*,'(3x," ... |",A)') trim(msg)
      end if
      msg_prev = trim(msg)
    end do
  end subroutine print_instr_opt


  function c_emis(mode, ityp, atls, cdyn) result(c)
  character(len=6)      :: c
    integer, intent(in) :: mode
    integer, intent(in) :: ityp
    integer, intent(in) :: atls
    integer, intent(in) :: cdyn

    if (btest(ityp,ITYP_MW)) then
      select case(mode)
      case(MODE_SEA_MODEL)
        c = 'MWSEA'
      case(MODE_DYNRET)
        write(c,'("DYNR",I2.2)') cdyn
      case(MODE_GRODY )
        c = 'GRODY '
      case(MODE_ATLAS )
        select case(atls)
        case(ATLS_TLSM)
          c = 'TELSEM'
        case(ATLS_CNRM)
          c = 'CNRM  '
        case default
          c = '??????'
        end select
      case default
        c = '??????'
      end select
    elseif (btest(ityp,ITYP_IR)) then
      select case(mode)
      case(MODE_SEA_MODEL)
        c = 'IRSEA '
      case(MODE_DYNRET)
        c = 'DYNRET'
        if (cdyn > 0) write(c(5:6),'(I2.2)') cdyn
      case(MODE_ATLAS )
        select case(atls)
        case(ATLS_CML_CL)
          c = 'CMLCL '
        case(ATLS_CML)
          c = 'CML07 '
        case(ATLS_UWIR)
          c = 'UWIR  '
        case default
          c = '??????'
        end select
      case default
        c = '??????'
      end select
    elseif (btest(ityp,ITYP_VIS)) then
      select case(mode)
      case(MODE_ATLAS )
        select case(atls)
        case(ATLS_BRDF)
          c = 'BRDF  '
        case default
          c = '??????'
        end select
      case default
        c = '??????'
      end select
    end if
  end function c_emis


  !----------------------------------------------------------------------------
  ! computes emissivity and sets it in ttovs and var, optionally stores it in obs%par
  !----------------------------------------------------------------------------
  subroutine set_emis(spot, ttovs, rtp, obs, stat)
    type (t_spot),       intent(inout)           :: spot    ! spot
    type (t_tovs),       intent(inout), target   :: ttovs   ! t_tovs variable
    type (t_rttov_prof), intent(inout)           :: rtp     ! variables passed to rttov
    type (t_obs),        intent(inout)           :: obs     ! observations, needed for GRODY
    integer,             intent(out)             :: stat    ! error status

    character(len=8), parameter   :: proc = 'set_emis'
    character(len=300)            :: msg
    type(t_rad_set),  pointer     :: rs
    type(t_emis_opt), pointer     :: eo
    type(t_tovs_instr)            :: ti
    logical                       :: l_grody
    real(wp),         allocatable :: e_grody(:)
    integer                       :: i_emis(m_chan)
    integer                       :: n_failed, aux, i, j, i1, in, k, n_chan, iat, ityp, styp, ic
    integer                       :: chan, ii,  mod
    integer                       :: i_prio, m_prio
    integer                       :: tfe
    logical                       :: lfail
    ! specularity calc.
    real(wp)                      :: fr,snf
    integer                       :: nspec
    ! surface type averaging
    real(wp), parameter           :: eps_wdyn = epsilon(1._wp)
    integer                       :: ns, is
    logical                       :: l_avg
    logical                       :: msk_dyn(1:n_styp)
    real(wp)                      :: w_dynret, ws, w, fac_dyn
    real(wp)                      :: max_dst(3)
    integer,  pointer             :: fl(:)
    integer,  allocatable         :: fl_  (:,:)
    real(wp), allocatable         :: emis_(:,:)
    ! Surface type statistics
    integer                       :: i_rs
    logical                       :: ld

    stat = 0

    call init_stat

    if (iand(ttovs%init, TTOVS_CI) == 0) call finish(proc, 'ttovs%ci not available')
    call get_tovs_rs(ttovs, rs=rs, ti=ti, i_rs=i_rs)

    call add_stst(i_rs,ttovs)

    n_chan = ttovs%nchan
    ns     = ttovs%ns
    m_prio = size(rs%emis_ind,1)

    do i = 1, ns
      styp = ttovs%rt_stype(i)
      if (styp < lbound(rs%emis_ind,2) .or. styp > ubound(rs%emis_ind,2)) then
        write(0,*) proc//': sat/grid=',rs%satid,rs%grid
        write(0,*) proc//': styp=',styp
        write(0,*) proc//': bounds of emis_ind=',lbound(rs%emis_ind,2),ubound(rs%emis_ind,2)
        call finish(proc, 'invalid surface type')
      end if
    end do

    l_avg = (ns > 1)

    if (l_avg) then
      allocate(fl(n_chan))
      allocate(fl_(ns,n_chan),emis_(ns,n_chan))
    else
      fl => ttovs%flag
    end if

    nspec = 0

    ld = ldeb(spot)
    if (ld) write(usd,*) dpref,proc,' start',ns,ttovs%rt_stype(1:ns),ttovs%w_styp(1:ns)
    if (ld) write(usd,*) dpref,proc,' loop_len',ti%n_instr,ns,m_prio,l_avg
    !-----------------------------------
    ! loop over instruments, select MW
    !-----------------------------------
    instr_loop: do i = 1, ti% n_instr                 ! do i (loop over instruments)
      ii = ti%ii(i)
      i1 = ti% o_ch_i(i) + 1
      in = ti% o_ch_i(i) + ti% n_ch_i(i)
      mask_pr(  1:n_chan) = .false.
      mask_dr(  1:n_chan) = .false.
      mask_at(:,1:n_chan) = .false.
      mask_br(  1:n_chan) = .false.
      mask_ss(  1:n_chan) = .false.

      if (l_avg) then
        emis_(:,:) = 0._wp
        fl_  (:,:) = 0
      end if

      ! Specularity calculation
      if (rs% iopts(ii)% do_lambertian) then
        nspec = nspec + 1
        if (nspec > ttovs%nspec) call finish(proc, 'ttovs%nspec too small')
        ttovs% spec(nspec) = 0._wp
        snf = min(max(rtp%snf,0._wp),1._wp)
        do is = 1, ttovs%ns
          styp = ttovs%rt_stype(is)
          select case(styp)
          case(rts_land,rts_ice)
            fr = ttovs%w_styp(is) * (1._wp - snf)
            ttovs% spec(nspec) = ttovs% spec(nspec) + fr * rs%iopts(ii)%specularity(styp)
            if (snf > 0._wp) then
              fr = ttovs%w_styp(is) * snf
              ttovs% spec(nspec) = ttovs% spec(nspec) + fr * rs%iopts(ii)%specularity_snow
            end if
          case(rts_sea)
            ttovs% spec(nspec) = ttovs% spec(nspec) + ttovs%w_styp(is) * rs%iopts(ii)%specularity(styp)
          end select
        end do
        if (ld) write(usd,*) dpref,proc,' spec',ttovs% spec(nspec)
      end if


      styp_loop: do is = 1, ns
        styp = ttovs%rt_stype(is)
        if (ld) write(usd,*) dpref,proc,' styp',styp,ttovs%w_styp(is),is

        if (l_avg) then
          fl(i1:in) = 0
        end if

        n_failed = 0

        prio_loop: do i_prio = 1, m_prio

          if (ld) write(usd,*) dpref,proc,' prio:',i_prio

          i_emis(i1:in) = rs%emis_ind(i_prio, styp, ttovs%ci(i1:in))


          if (ld) then
            do j = i1,in
              if (i_emis(j) > 0) then
                write(usd,*) dpref,proc,' i_emis',j,i_emis(j),fl(j), &
                     iand(fl(j), EMIS_SUCCESS),rs%emis_opt(i_emis(j))%mode
              else
                write(usd,*) dpref,proc,' i_emis',j,i_emis(j),fl(j), &
                     iand(fl(j), EMIS_SUCCESS)
              end if
            end do
          end if

          mask_pr(i1:in) = (iand(fl(i1:in), EMIS_SUCCESS) == 0) .and. &
                           (i_emis(i1:in) > 0)

          if (ld) write(usd,*) dpref,proc,' mask_pr',count(mask_pr(i1:in))
          if (count(mask_pr(i1:in)) <= 0) exit prio_loop

          mask_dr(  i1:in) = .false.
          mask_at(:,i1:in) = .false.
          mask_br(  i1:in) = .false.
          mask_ss(  i1:in) = .false.
          l_grody          = .false.
          max_dst          = 0._wp
          do j = i1, in
            if (.not.mask_pr(j)) cycle
            rtp% emis(1,j) = -1._wp
            ityp = instr_type(rs,ii,chan=j)
            eo => rs%emis_opt(i_emis(j))
            mod = eo%mode
            if (mod == MODE_DYNRET .or. mod == MODE_ATLAS) then
              ic = -1
              if (eo%n_chan > 0) then
                do k = 1, eo%n_chan
                  if (eo%chan(k) == rs%chan(ttovs%ci(j))) then
                    ic = k
                    exit
                  end if
                end do
              end if
            end if
            if (ld) write(usd,*) dpref,proc,' chan',i,i_prio,j,mod,eo%n_chan,ic,ityp
            select case(mod)
            case(MODE_DYNRET)
              if (ic <= 0) then
                write(0,*) 'spot/instr/chan',spot%hd%id,ii,rs%instr(ii),j,rs%chan(ttovs%ci(j))
                if (associated(eo%chan)) write(0,*) 'eo',eo%n_chan,eo%chan
                if (associated(eo%chan)) write(0,*) 'eo%decr',trim(eo%descr)
                call finish(proc, 'Did not find channel in eo%chan')
              end if
              mask_dr(j) = (eo%cdyn(ic) == 0)
            case(MODE_ATLAS)
              if (btest(ityp, ITYP_VIS)) then
                mask_br(j) = (eo%atls==ATLS_BRDF)
              else
                mask_at(eo%atls,j) = (eo%n_chan == 0 .or. ic > 0)
                max_dst(eo%atls) = eo%max_dst
              end if
            case(MODE_GRODY)
              l_grody = .true.
            case(MODE_SEA_MODEL)
              if (l_avg) mask_ss(j) = .true.
            end select
          end do

          ! Retrieve emissivities
          if (ld) write(usd,*) dpref,proc,' mask_dr',count(mask_dr(i1:in))
          if (count(mask_dr(i1:in)) > 0) &
               call call_rttov('e', ttovs, rtp, stat, msk=mask_dr(1:n_chan), styp=styp, ldebug=ld,&
                               obs_bt=real(obs% body(spot%o%i+1 : spot%o%i+ttovs%nchan)%o,kind=wp) )
          do iat = 1, 3
            if (ld) write(usd,*) dpref,proc,' mask_at',iat,count(mask_at(iat,i1:in))
            if (count(mask_at(iat,i1:in)) > 0) &
                 call call_rttov('a', ttovs, rtp, stat, msk=mask_at(iat,1:n_chan), atlas=iat, &
                                 styp=styp, max_dst=max_dst(iat), ldebug=ld)
          end do
          if (ld) write(usd,*) dpref,proc,' mask_br',count(mask_br(i1:in))
          if (count(mask_br(i1:in)) > 0) &
               call call_rttov('b', ttovs, rtp, stat, msk=mask_br(1:n_chan), atlas=iat, &
                               styp=styp, ldebug=ld)
          if (ld) write(usd,*) dpref,proc,' mask_ss',count(mask_ss(i1:in))
          if (count(mask_ss(i1:in)) > 0) then
            call call_rttov('s', ttovs, rtp, stat, msk=mask_ss(1:n_chan), styp=styp, &
                 ldebug=ld, ierr_no_stop=(/8/)) ! 8: ERR_NO_RTTOV_LIB (happens for version before 13.2)
            if (stat == 8) call call_rttov('f', ttovs, rtp, stat, msk=mask_ss(1:n_chan), &
                                           styp=styp, ldebug=ld)
          end if
          if (l_grody) then
            allocate(e_grody(n_chan))
            call emis_grody(rs% instr(ii), rs% chan(ttovs%ci(i1:in)), &
                 real(obs% body(spot%o%i+1 : spot%o%i+ttovs%nchan)%o, kind=wp), spot=spot, emi=e_grody, stat=stat)
            if (stat > 0 ) then
              call finish(proc,'emis_grody failed')
            else if (stat < 0) then
              write(msg, '(A,I2,A,F6.2,A,F7.2,A,I2,A)') 'Grody method can not be used for instrument ', &
                   rs%instr(ii), ' , lat ', spot%col%c%dlat, " and lon ", spot%col%c%dlon, &
                   " because channel ", -stat, "is missing."
              call warning(msg)
              stat = 0
            end if
          end if

          ! Distribute results to the right channels according to namelist
          chan_loop: do j = i1, in                       ! do j (loop over channels)
            if (.not.mask_pr(j)) cycle
            chan = rs%chan(ttovs%ci(j))
            eo => rs%emis_opt(i_emis(j))
            mod = eo%mode
            ityp = instr_type(rs, i, chan=j)

            select case (mod)
            case (MODE_SEA_MODEL)
              if (.not.l_avg) rtp% emis(1,j) = -1._wp
              tfe = TF_EMIS_SEA_MOD
            case (MODE_DYNRET)
              if (mask_dr(j)) then
                ! nothing to do
              else
                aux = -1
                do k = 1, eo%n_chan
                  if (eo%chan(k) == chan) then
                    aux = eo%cdyn(k)
                    exit
                  end if
                end do
                if (aux <= 0) then
                  write(0,*) 'spot',spot%hd%id
                  write(0,*) 'rs',rs%satid,rs%grid,rs%instr(ii),chan
                  write(0,*) 'eo%chan',eo%chan(1:eo%n_chan)
                  call finish(proc,'failed to find cdyn entry')
                end if
                do k=i1,in
                  if (rs% chan(ttovs%ci(k)) == aux) then
                    rtp% emis(1,j) = rtp% emis(1,k)
                    exit
                  end if
                end do
              end if
              tfe = TF_EMIS_DYNRET
            case (MODE_GRODY)
              rtp%emis(1,j) = e_grody(j-i1+1)
              tfe = TF_EMIS_GRODY
            case (MODE_ATLAS)
              ! Nothing to distribute, the selection of right channels happens in call_rttov
              ! (x% emis(1,:) = unpack(emissiv(:,1), msk, x% emis(1,:)), see call_rttov)
              if (btest(ityp,ITYP_MW)) then
                if (mask_at(ATLS_TLSM  ,j)) tfe = TF_EMIS_TLSM
                if (mask_at(ATLS_CNRM  ,j)) tfe = TF_EMIS_CNRM
              end if
              if (btest(ityp,ITYP_IR)) then
                if (mask_at(ATLS_UWIR  ,j)) tfe = TF_EMIS_UWIR
                if (mask_at(ATLS_CML   ,j)) tfe = TF_EMIS_CAMEL07
                if (mask_at(ATLS_CML_CL,j)) tfe = TF_EMIS_CAMELCL
              end if
              if (btest(ityp,ITYP_VIS)) then
                if (mask_br(j)) tfe = TF_REFL_BRDF
              end if
            end select

            if (ld) write(usd,*) dpref,proc,' emis',j,chan,rtp%emis(1,j)

            ! NOTE on TF_EMIS_FAILED:
            ! If TF_EMIS_FAILED is set in the final tovs%flag then the observation will
            ! be rejected in mo_fg_checks.f90. We would like to reject observations if
            ! "non-sea-emissivity-methods" fail and we have to use the sea emissivity
            ! model (MODE_SEA_MODEL) as a fallback over non-sea terrain (land/ice).
            ! Therefore, the TF_EMIS_FAILED bit (coming from previous methods) is kept in
            ! the case mod==MODE_SEA_MODEL below in contrast to the other models.
            lfail = .false.
            if (mod == MODE_SEA_MODEL) then
              fl(j) = ibset(fl(j), tfe)
              if (l_avg) then
                if (rtp% emis(1,j) < 0._wp .or. rtp% emis(1,j) > 1._wp) then
                  ! write(0,*) 'bad emis from sea_model',spot%hd%id,j,rtp%emis(1,j), '--> restrict to range 0..1'
                  rtp% emis(1,j) = max(min(rtp% emis(1,j), 1._wp), 0._wp)
                  fl(j) = ibset(fl(j), TF_EMIS_FAILED)
                  lfail = .true.
                  n_failed = n_failed + 1
                end if
              end if
            elseif (rtp% emis(1,j) > 0._wp .and. rtp% emis(1,j) < 1._wp) then
              fl(j) = ibset(fl(j), tfe)
              fl(j) = ibclr(fl(j), TF_EMIS_FAILED) ! Unset TF_EMIS_FAILED from previous method
            else
              fl(j) = ibset(fl(j), TF_EMIS_FAILED)
              lfail = .true.
              n_failed = n_failed + 1
            end if
            if (ld) write(usd,*) dpref,proc,' emis tf',j,chan,tfe,fl(j),lfail

            call add_stat(eo%id, styp, chan, lfail, rtp%emis(1,j), mix=l_avg)

          end do chan_loop

          if (n_failed == 0) EXIT prio_loop

        end do prio_loop

        if (l_avg) then
          emis_(is,i1:in) = rtp%emis(1,i1:in)
          fl_  (is,i1:in) = fl(i1:in)
          if (ld) then
            do j = i1,in
              write(usd,*) dpref,proc,' styp_emis',j,rtp%emis(1,j),ttovs%w_styp(is),fl(j)
            end do
          end if
        end if
        if (any(iand(fl(i1:in),EMIS_CALC) == 0)) then
          write(0,*) proc//' failed for spot=',spot%hd%id,'instr=',i,'styp=',styp,is
          do j = i1,in
            if (iand(ttovs%flag(j),EMIS_CALC) == 0) write(0,*) proc//' failed for instr',i,&
                 'chan',j,rs%chan(ttovs%ci(j))
          end do
          call finish(proc,'Failed to set emissivity')
        end if

      end do styp_loop

      if (l_avg) then
        do j = i1,in
          if (ld) write(usd,*) dpref,proc,' avg',rs%instr(ti%ii(i)),rs%chan(ttovs%ci(j)),&
               ns,ttovs%rt_stype(1:ns),ttovs%w_styp(1:ns)
          msk_dyn(1:ns) = (fl_(1:ns,j)==2**TF_EMIS_DYNRET)
          w_dynret = sum(ttovs%w_styp(1:ns), mask=msk_dyn(1:ns))
          if (ld) write(usd,*) dpref,proc,' avg dyn',msk_dyn(1:ns),w_dynret
          if (w_dynret > dynret_w_pref) then
            ! If there is a significant part of "dynamical retrieval surface(s)", we use ONLY the
            ! dynamical retrieval, since averaging dynamical retrieval results with other
            ! emissivities (atlas, FASTEM, ...) does not make physical sense.
            ! However, nevertheless it might be useful to mix dynamically retrieved emissivity
            ! (representing the whole FOV) with the value of a physical emissivity model
            ! (representing e.g. a sea surface part of the FOV). The sea surface value might be
            ! interpreted as an a priori first guess and the dynamically retrieved value as an
            ! observation. The larger the sea fraction the more we should trust in the first guess.
            ! (David Duncan, personal communication ITSC-24). This idea might be formuated as a
            ! simple weighted mean (as in the "else" code chunk). This is achieved by setting
            ! dynret_w_pref >= 1.
            do is = 1, ns
              if (msk_dyn(is)) then
                rtp%emis(1,j) = emis_(is,j)
                ttovs%flag(j) = ior(ttovs%flag(j), fl_(is,j))
                if (ld) write(usd,*) dpref,proc,' avg dyn res',is,rtp%emis(1,j),fl_(is,j)
                exit
              end if
            end do
          else
            rtp%emis(1,j) = 0._wp
            ws            = 0._wp
            do is = 1, ns
              if (dynret_avg > 0) then
                if (w_dynret > eps_wdyn .and. dynret_w_pref > eps_wdyn) then
                  ! with the aid of fac_dyn we modify the weights of the different surfaces
                  ! such that DYNRET-surfaces have a weight of 1. for w_dynret == dynret_w_pref.
                  fac_dyn = 1._wp/dynret_w_pref
                  if (.not.msk_dyn(is)) then
                    fac_dyn = (1._wp - fac_dyn * w_dynret) / (1._wp - w_dynret)
                  end if
                else
                  fac_dyn = 1._wp
                end if
              else
                if (msk_dyn(is)) then
                  fac_dyn = 0._wp
                else
                  fac_dyn = 1._wp
                end if
              end if
              w = fac_dyn * ttovs%w_styp(is)
              if (w > 0._wp) then
                rtp%emis(1,j) = rtp%emis(1,j) + w * emis_(is,j)
                ws            = ws            + w
                ttovs%flag(j) = ior(ttovs%flag(j), fl_(is,j))
              end if
              if (ld) write(usd,*) dpref,proc,' avg sum',ns,ttovs%rt_stype(is),emis_(is,j),&
                   fl_(is,j),ttovs%w_styp(is),w,ws,fac_dyn
            end do
            if (ws > 0._wp) then
              rtp%emis(1,j) = rtp%emis(1,j)/ws
            else
              call finish(proc,'invalid emissivity averaging weights')
            end if
            if (ld) write(usd,*) dpref,proc,' avg res',rtp%emis(1,j),ws
          end if
          if (ld) write(usd,*) dpref,proc,' final_emis',j,rtp%emis(1,j)
        end do
      end if
      if (any(rtp%emis(1,i1:in) > 1._wp)) then
        write(0,*) 'spot',spot%hd%id,ns,ttovs%rt_stype(1:ns),ttovs%w_styp(1:ns)
        do j = i1,in
          if (rtp%emis(1,j) >= 1._wp) then
            if (obs% body(spot%o%i+j)%use%state <= STAT_REJECTED ) then
              write(0,'("bad_emis=",G13.6," --> set to -1!  instr=",I3," chan=",I5," lat=",F7.2," lon=",F7.2)') &
                   rtp%emis(1,j),rs%instr(ii),rs%chan(ttovs%ci(j)),spot%col%c%dlat,spot%col%c%dlon
              rtp%emis(1,j) = -1._wp
            else
              call finish(proc,'bad emissivity (>1.)')
            end if
          end if
        end do
      end if

    end do instr_loop

    if (ld) then
      do j = 1, n_chan
        write(usd,*) dpref,proc,' final',j,rtp%emis(1,j),ttovs%flag(j),&
             (/ (btest(ttovs%flag(j),TF_EMIS(i)), i=lbound(TF_EMIS,1),ubound(TF_EMIS,1)) /)
      end do
    end if

    call store(obs, spot, ttovs, tovs_io=TTOVS_FLAG+TTOVS_SPEC)

    if (l_avg) deallocate(fl)

  end subroutine set_emis
!------------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  ! Computes Grody emissivity
  !----------------------------------------------------------------------------
  subroutine emis_grody(ins, cha, obt, spot, emi, stat)
    integer     ,           intent(in)   :: ins       !instrument rttov id
    integer     ,           intent(in)   :: cha(:)    ! channel numbers
    real(wp)    ,           intent(in)   :: obt(:)    ! observed brightness temperatures
    type(t_spot),           intent(in)   :: spot      ! spot info
    real(wp)    ,           intent(out)  :: emi(:)    ! emissivity
    integer     , optional, intent(out)  :: stat      ! error status

    character(len=100)     :: msg
    real(wp)               :: e_wndw(4)
    integer                :: id_sfc
    integer                :: ierr
    integer                :: k

    !---------------------------------------------------------------
    ! Check that it is for AMSU-A or ATMS instruments
    !---------------------------------------------------------------
    if ( ins /= 3 .and. ins /= 19 ) then
      write(msg, '(A,I2)') 'Grody MW emissivity can not be used for instrument', ins
      call warning(msg)
      return
    end if

    !---------------------------------------------------------------
    ! Check for mismatched array sizes
    !---------------------------------------------------------------
    if ( size(cha) /= size(emi) .or. size(cha) /= size(obt)) then
      call finish('emis_grody','Missmatch in input arrays')
    end if

    !---------------------------------------------------------------
    ! Fetch window emissivities
    !---------------------------------------------------------------
    call mw_emiss(ins, cha, obt, spot% stzen, spot% sl_bg, spot% ts_bg, ierr, id_sfc=id_sfc, e_ret=e_wndw)
    if (present(stat)) stat = ierr

    !---------------------------------------------------------------
    ! Set the emissivities for all channels
    !---------------------------------------------------------------
    do k = 1, size(cha)
      select case (cha(k))
      case(1:3)
        emi(k) = e_wndw(k)
      case(4:14)
        emi(k) = e_wndw(3)
      case(15)
        if (ins == 19) then   !atms
          emi(k) = e_wndw(3)
        else                  !amsua
          emi(k) = e_wndw(4)
        end if
      case(16)
        emi(k) = e_wndw(4)
      case default ! humidity sounding channels
        emi(k) = 0._wp
      end select
    end do

  end subroutine emis_grody


  elemental subroutine destruct_emis_stat(st)
    type(t_emis_stat), intent(inout) :: st
    if (associated(st%chan )) deallocate(st%chan )
    if (associated(st%n    )) deallocate(st%n    )
    if (associated(st%nfail)) deallocate(st%nfail)
    if (associated(st%n_mix)) deallocate(st%n_mix)
    if (associated(st%n_mix_fail)) deallocate(st%n_mix_fail)
    if (associated(st%emiss)) deallocate(st%emiss)
  end subroutine destruct_emis_stat


  subroutine init_stat
    type(t_rad_set),   pointer :: rs
    type(t_emis_opt),  pointer :: eo
    type(t_emis_stat), pointer :: st
    integer :: iset, iopt, is, nch, i, j, ioff, n

    if (associated(stat)) return

    nstat = 0
    do iset = 1, n_set
      rs => rad_set(iset)
      do iopt = 1, rs%n_emis_opt
        eo => rs%emis_opt(iopt)
        nstat = nstat + eo%ns
      end do
    end do
    allocate(stat(nstat))

    n = 0
    do iset = 1, n_set
      rs => rad_set(iset)
      do iopt = 1, rs%n_emis_opt
        eo => rs%emis_opt(iopt)
        do is = 1, eo%ns
          n = n + 1
          st => stat(n)
          st%satid = rs%satid
          st%id    = eo%id
          st%instr = eo%inst
          st%prio  = eo%prio
          st%styp  = eo%styp(is)
          st%descr  = trim(eo%descr)
          if (eo%n_chan > 0) then
            nch = eo%n_chan
          else
            ioff = -1
            nch  =  0
            do i = 1, rs%n_instr
              if (rs%instr(i) == eo%inst) then
                if (ioff < 0) ioff = rs%o_ch_i(i)
                nch = nch + rs%n_ch_i(i)
              endif
            end do
          end if
          st%n_chan = nch
          allocate(st%chan(nch), st%n(nch), st%nfail(nch), st%n_mix(nch), st%n_mix_fail(nch), st%emiss(nch))
          st%n = 0 ; st%nfail = 0 ; st%n_mix = 0 ; st%n_mix_fail = 0 ; st%emiss = 0._wp
          if (eo%n_chan > 0) then
            st%chan(1:nch) = eo%chan(1:nch)
          else
            st%chan(1:nch) = rs%chan(ioff+1:ioff+nch)
          end if
        end do
      end do
    end do

    nstst = 0
    do i = 0,n_styp-1
      nstst = nstst + 2**i
    end do
    allocate(stst(n_set,nstst))
    do i = 1, nstst
      stst(:,i)%ist = i
      n = 0
      do j = 0,n_styp-1
        if (btest(i,j)) n = n + 1
      end do
      stst(:,i)%ns  = n
    end do

  end subroutine init_stat


  subroutine add_stat(id, styp, chan, fail, emis, mix)
    integer,  intent(in) :: id
    integer,  intent(in) :: styp
    integer,  intent(in) :: chan
    logical,  intent(in) :: fail
    real(wp), intent(in) :: emis
    logical,  intent(in) :: mix

    type(t_emis_stat), pointer :: st
    integer,           save    :: ii = 1
    integer                    :: i, ic
    logical                    :: lf

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
        call finish('add_stat', 'Did not find emis_stat entry')
      end if
    end if
    st => stat(ii)
    lf = .false.
    do ic = 1, st%n_chan
      if (st%chan(ic) == chan) then ; lf = .true. ; exit ; endif
    end do
    if (.not.lf) then
      write(0,*) 'id,styp,chan',id,styp,chan
      call finish('add_stat', 'Did not find channel')
    end if
    st%n(ic) = st%n(ic) + 1
    if (fail) then
      st%nfail(ic) = st%nfail(ic) + 1
    else
      st%emiss(ic) = st%emiss(ic) + emis
    end if
    if (mix) then
      st%n_mix(ic) = st%n_mix(ic) + 1
      if (fail) st%n_mix_fail(ic) = st%n_mix_fail(ic) + 1
    end if

  end subroutine add_stat

  subroutine add_stst(i_rs,ttovs)
    integer,      intent(in) :: i_rs
    type(t_tovs), intent(in) :: ttovs

    type(t_styp_stat), pointer :: s
    integer                    :: is,i,ist,styp

    if (i_rs <= 0 .or. i_rs > n_set) call finish('add_stst','invalid i_rs')

    ist = 0
    do is = 1, ttovs%ns
      ist = ibset(ist, ttovs%rt_stype(is))
    end do
    s => null()
    do i = 1,nstst
      if (stst(i_rs,i)%ist == ist) then
        s => stst(i_rs,i)
        exit
      end if
    end do
    if (.not.associated(s)) call finish('add_stst','failed to determine entry in stst')

    s%n = s%n + 1
    do is = 1, ttovs%ns
      styp = ttovs%rt_stype(is)
      s%w(styp) = s%w(styp) + ttovs%w_styp(is)
    end do

  end subroutine add_stst


  subroutine print_emis_stat
    type(t_emis_stat), pointer :: st
    character(len=300)         :: msg, msg_prev
    real(wp)                   :: emiss, fr_mix, fr_fail
    integer                    :: n_omit
    integer                    :: i, j, ic, nch
    integer                    :: mstat
    ! surface type statistics
    type(t_styp_stat), pointer :: sts
    type(t_styp_stat)          :: stst_
    integer                    :: i_rs
    logical                    :: l_swap, lstyp(n_styp)

    if (.not.l_pr) return ! Not really necessary, but avoids waiting for MPI below
    mstat = p_max(nstat)
    if (mstat <= 0) return

    if (nstat <= 0) call init_stat
    ! summarize over all PEs
    do i = 1, nstat
      st => stat(i)
      st%n    (1:st%n_chan) = p_sum(st%n    (1:st%n_chan))
      st%nfail(1:st%n_chan) = p_sum(st%nfail(1:st%n_chan))
      st%n_mix(1:st%n_chan) = p_sum(st%n_mix(1:st%n_chan))
      st%n_mix_fail(1:st%n_chan) = p_sum(st%n_mix_fail(1:st%n_chan))
      st%emiss(1:st%n_chan) = p_sum(st%emiss(1:st%n_chan))
    end do
    nstst = p_max(nstst)
    do i_rs = 1, n_set
      do i = 1, nstst
        sts => stst(i_rs,i)
        sts%n    = p_sum(sts%n)
        sts%w(:) = p_sum(sts%w(:))
      end do
    end do

    if (dace%lpio) then
      write(*,*)
      write(*,*) 'surface type statistics (TOVS):'
      do i_rs = 1, n_set
        ! Sort statistics
        do i = nstst-1,1,-1
          do j = 1, i
            l_swap = .false.
            if (stst(i_rs,j)%ns > stst(i_rs,j+1)%ns) then
              l_swap = .true.
            elseif (stst(i_rs,j)%ns == stst(i_rs,j+1)%ns) then
              l_swap = stst(i_rs,j)%ist > stst(i_rs,j+1)%ist
            end if
            if (l_swap) then
              stst_          = stst(i_rs,j  )
              stst(i_rs,j  ) = stst(i_rs,j+1)
              stst(i_rs,j+1) = stst_
            end if
          end do
        end do
        ! print statistics
        write(*,*)
        write(*,'(3x,"satid=",I3," grid=",I3)') rad_set(i_rs)%satid,rad_set(i_rs)%grid
        do i = 1, nstst
          sts => stst(i_rs,i)
          if (sts%n > 0) then
            msg = ''
            do j = 0, n_styp-1
              if (btest(sts%ist,j)) then
                write(msg(j*7+1:),'(A6)') rts_name(j)
              else
                write(msg(j*7+1:),'("------")')
              end if
            end do
            write(*,'(5x,A,2x,I8,3(1x,F8.6))') trim(msg),sts%n,sts%w(:)/sts%n
          end if
        end do
      end do


      write(*,*)
      write(*,*) 'emissivity statistics:'
      write(*,*) '(*_m: #mixed surfaces, *_f: #method failed)'
      do i = 1, nstat
        st => stat(i)
        nch = st%n_chan
        if (sum(st%n(1:nch)) > 0) then
          write(*,*)
          write(*,'(3x,A," satid=",I3," instr=",I3," styp=",I1)') &
               trim(st%descr), st%satid, st%instr, st%styp
          write(*,'(5x," chan",1x,"       n",1x,"     n_f",1x,"   n/n_f",1x," emiss",1x,&
               &"     n_m","    n_m_f",1x,"n_m/n_m_f",1x,"n_m_f/n_f")')
          n_omit   = 0
          msg_prev = ''
          do ic = 1, nch
            if (st%n(ic) > 0) then
              emiss = -1._wp
              if (st%n(ic) > st%nfail(ic)) emiss=st%emiss(ic)/(st%n(ic) - st%nfail(ic))
              fr_mix= 0._wp
              if (st%n_mix(ic) > 0) fr_mix=(100._wp*st%n_mix_fail(ic))/st%n_mix(ic)
              fr_fail= 0._wp
              if (st%nfail(ic) > 0) fr_fail=(100._wp*st%n_mix_fail(ic))/st%nfail(ic)
              write(msg,'(2(1x,I8),1x,F7.3,"%",1x,F6.3,2(1x,I8),2(2x,F7.3,"%"))') &
                   st%n(ic),st%nfail(ic),(100._wp*st%nfail(ic))/real(st%n(ic),wp),emiss,&
                   st%n_mix(ic),st%n_mix_fail(ic),fr_mix,fr_fail
              if (msg /= msg_prev .or. ic == nch) then
                write(*,'(5x,I5,A)') st%chan(ic),trim(msg)
                n_omit = 0
              else
                n_omit = n_omit + 1
                if (n_omit == 1) write(*,'(5x,"  ...",A)') trim(msg)
              end if
              msg_prev = trim(msg)
            end if
          end do
        end if
      end do
    end if

    call destruct_emis_stat(stat)
    deallocate(stat)
    nstat = 0
    l_pr  = .false.

  end subroutine print_emis_stat



END MODULE mo_emis
