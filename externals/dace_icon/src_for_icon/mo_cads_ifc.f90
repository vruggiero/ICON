!
!+ Interface to the CADS (CLoud and Aerosol Detection Software) by NWPSAF.
!
MODULE mo_cads_ifc
!
! Description:
! Interface to the CADS (CLoud and Aerosol Detection Software) by NWPSAF.
! This module allows to call different version of CADS.
! CADS contains the well-known McNally-Watts cloud detection scheme.
!
! Current Maintainer: DWD, Robin Faulwetter
!    phone: +49 69 8062 2746
!    fax:   +49 69 8062 3721
!    email: robin.faulwetter@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_1         2020/12/15 Robin Faulwetter
!  intial version
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Robin Faulwetter  DWD
! Olaf Stiller      DWD
!=======================================================================

  !=============
  ! Modules used
  !=============

  use mo_kind,              only: wp, sp
  use mo_exception,         only: finish
  use mo_mpi_dace,          only: dace,                            &
                                  p_bcast
  use mo_namelist,          only: position_nml,                    &! position namelist
                                  nnml,                            &! namelist Fortran unit
                                  POSITIONED                        ! ok code from position_nml
  use mo_fortran_units,     only: get_unit_number,                 &! reserve a unit number
                                  return_unit_number                ! release the un
  use mo_rad,               only: t_rad_set,                       &
                                  rad_set,                         &
                                  n_set,                           &
                                  lev2p
  use mo_t_obs,             only: usd
  use mo_t_tovs,            only: t_tovs,                          &
                                  t_tovs_instr,                    &
                                  get_tovs_rs,                     &
                                  get_im_ch_ind,                   &
                                  IMCL_FRAC, IMCL_MEAN, IMCL_STDV, &
                                  IMCH_MEAN_O, IMCH_MEAN_B,        &
                                  IMCH_STDV, IMCH_EMIS, TTOVS_IM
  use mo_cloud_ir,          only: cloud_detect_setup_hss,          &!
                                  cloud_detect_hss,                &!
                                  l_debug_2 => ldebug_ass
  use mo_cads31,            only: cads_setup_cloud_31,             &!
                                  cads_detect_cloud_31,            &!
                                  cads_setup_aerosol_31,           &!
                                  cads_detect_aerosol_31,          &!
                                  cads_setup_land_sensitivity_31,  &!
                                  cads_detect_land_sensitivity_31, &!
                                  cads_setup_trace_gas_31,         &!
                                  cads_detect_trace_gas_31,        &!
                                  p_bcast_cloud_31,                &!
                                  p_bcast_aerosol_31,              &!
                                  p_bcast_land_sens_31,            &!
                                  p_bcast_trace_gas_31,            &!
                                  set_cloud_detect_chans_31,       &!
                                  check_sensor_id_31,              &!
!                                  get_cloud_opt_31,                &!
                                  l_debug_31 => l_debug,           &!
                                  dwd_sort_31 => dwd_sort,         &!
                                  l_rank_p_31 => l_rank_p,         &!
                                  t_aer_out_31 => t_aer_out,       &!
                                  t_cld_out_31 => t_cld_out,       &!
                                  usd_31 => usd

  implicit none
  private

  public :: cads_init_ifc
  public :: cads_ifc
  public :: toplev
  public :: bottomlev
  public :: MNW_CLOUD
  public :: MNW_AEROSOL
  public :: MNW_LAND
  public :: MNW_TRACEGAS
  public :: MNW_SURF
  public :: CLOUDY
  public :: CLOUDFREE
  public :: l_aerosol
  public :: l_land_sens
  public :: l_trace_gas
  public :: use_im_frac
  public :: IMFRC_FAR
  public :: IMFRC_CADS

  character(len=20) :: version   = '2_0'
  integer           :: version_major = -1
  integer           :: version_minor = -1
  ! Remarks on different versions:
  ! 1.) In particular with l2c_type=2 (default), but also with l2c_type=3 we have
  !     different results because of the different sorting routines used in V2 and V3
  !     (V2: sortrx(DWD), V3: CADS_Detect_Cloud_Heapsort). The sorting routines behave
  !     differently for channels that have exactly the same height. This is usual
  !     for l2c_type=2, but it might also happen for l2c_type=3 for the window
  !     channels. In V3 the DWD sorting routine might be used by setting dwd_sort=T.
  ! 2.) The cross-band-flagging uses different heights:
  !     V2  : flagging below i__firstCloudyChannel
  !     V3.1: flagging below i__LastClearChannel
  ! 3.) The default for n__bandtouse (band to use in crossband-flagging) differs:
  !     V2  : n__bandtouse = 1 2 1 4 5 , i.e. no crossband-flagging für bands 2,4,5
  !     V3.1: n__bandtouse = 1 1 1 1 1 , i.e. crossband-flagging for all bands with band 1
  ! 4.) For flagging those channels, that were not supplied to CADS, the following
  !     routines require the cloud (top) level for each band. In version 2 this is
  !     calculated after the CADS-call by using the flags given by CADS. For each
  !     band the level of the HIGHEST CLOUDY channel IN THE BAND is used as the
  !     cloud level. This is very suboptimal, since this might be very different
  !     (and much lower) compared to the cloud level in CADS, which is the level of
  !     the LOWEST CLOUDFREE channel in the BAND USED FOR CROSS-BAND-FLAGGING.
  !     In V3.1 the cloud level given by CADS is used.
  !     For backwards-compatibility reasons the problematic behaviour is kept for V2.

  logical           :: init       = .false.

  integer, parameter :: mx_sensor = 200
  integer, parameter :: n_sensor   = 9
  integer, parameter :: sensors(n_sensor) = (/11,16,27,28,57,59,94,97,98/)

  ! Bits to control use of cloud fraction from imager:
  integer, parameter :: IMFRC_FAR  = 0
  integer, parameter :: IMFRC_CADS = 1

  logical, save      :: l_first = .true.

  ! Namelist variables
  real(kind=wp)     :: toplev                        = 0._wp
  real(kind=wp)     :: bottomlev                     = 1.e20_wp
  logical           :: use_all_chans     (mx_sensor) = .true.
  logical           :: l_aerosol         (mx_sensor) = .false.
  logical           :: l_land_sens       (mx_sensor) = .false.
  logical           :: l_trace_gas       (mx_sensor) = .false.
  logical           :: use_bcor_cloud    (mx_sensor) = .true.
  logical           :: use_bcor_aerosol  (mx_sensor) = .true.
  logical           :: use_bcor_trace_gas(mx_sensor) = .true.
  logical           :: dwd_sort                      = .false. ! Use dwd sort routine
  integer           :: use_im_frac       (mx_sensor) = 0
  integer           :: surf_check        (mx_sensor) = 0
  integer           :: si_surf_check     (mx_sensor) = 0
  real(wp)          :: si_thresh         (mx_sensor) = huge(0._wp)
  real(wp)          :: imager_cld_thresh (mx_sensor) = 5._wp   ! {%}
  real(wp)          :: imager_psurf_frac (mx_sensor) = 0.1_wp
  integer           :: l2c_aer                       = 0
  logical           :: write_nhgt                    = .false.
  logical           :: aer_out                       = .false.
  logical           :: cld_out                       = .false.


  NAMELIST /CADS/ version, toplev, bottomlev, use_all_chans,    &
       l_aerosol, l_land_sens, l_trace_gas, use_bcor_cloud,     &
       use_bcor_aerosol, use_bcor_trace_gas, dwd_sort,          &
       use_im_frac, imager_cld_thresh, imager_psurf_frac,       &
       l2c_aer, write_nhgt, aer_out, cld_out,                   &
       surf_check, si_surf_check, si_thresh


  ! Flags set in this routine
  ! Bits
  integer, parameter :: MNW_CLOUD    = 0
  integer, parameter :: MNW_AEROSOL  = 1
  integer, parameter :: MNW_LAND     = 2
  integer, parameter :: MNW_TRACEGAS = 3
  integer, parameter :: MNW_SURF     = 4
  ! flag values
  integer, parameter :: CLOUDFREE = 0
  integer, parameter :: CLOUDY    = 2 ** MNW_CLOUD
  ! surf check bits
  integer, parameter :: SC_LAND   = 0
  integer, parameter :: SC_SEA    = 1
  integer, parameter :: SC_ICE    = 2

contains

  subroutine cads_init_ifc
    character(len=300) :: msg = ''
    integer :: i, j, stat
    integer :: l

    ! required for "use_all_chans" option
    integer :: isens
    integer, parameter :: mx_bands = 8
    integer, parameter :: mx_chans = 8461
    integer            :: nbands   = 0
    integer            :: nchans(mx_bands) = 0
    integer            :: chans(mx_chans, mx_bands) = 0
    integer            :: istart, iend, ii, iset
    type(t_rad_set), pointer :: rs => null()

    if (init) call finish('cads_init_ifc','double call to cads_init_ifc')

    ! Read CADS namelist
    if (dace%lpio) then
       call position_nml('CADS',lrewind=.true., status=stat)
       if (stat == POSITIONED) read(nnml, nml=CADS)
       write(*,*)
       write(*,*) 'CADS namelist'
       write(*,'(3x,"version (namelist)  = ",A)') trim(version)
    end if
    call p_bcast(version,           dace%pio)
    call p_bcast(toplev,            dace%pio)
    call p_bcast(bottomlev,         dace%pio)
    call p_bcast(use_all_chans,     dace%pio)
    call p_bcast(l_aerosol,         dace%pio)
    call p_bcast(l_land_sens,       dace%pio)
    call p_bcast(l_trace_gas,       dace%pio)
    call p_bcast(use_bcor_cloud,    dace%pio)
    call p_bcast(use_bcor_aerosol,  dace%pio)
    call p_bcast(use_bcor_trace_gas,dace%pio)
    call p_bcast(dwd_sort,          dace%pio)
    call p_bcast(surf_check,        dace%pio)
    call p_bcast(si_surf_check,     dace%pio)
    call p_bcast(si_thresh,         dace%pio)
    call p_bcast(use_im_frac,       dace%pio)
    call p_bcast(imager_cld_thresh, dace%pio)
    call p_bcast(imager_psurf_frac, dace%pio)
    call p_bcast(l2c_aer,           dace%pio)
    call p_bcast(write_nhgt,        dace%pio)
    call p_bcast(aer_out,           dace%pio)
    call p_bcast(cld_out,           dace%pio)

    ! Determine major/minor version from version
    l = len_trim(version)
    do i = 1, l
      if (verify(version(i:i), '0123456789') /= 0) exit
    end do
    i = min(i-1,l)
    read(version(1:i),*,iostat=stat) version_major
    if (stat /= 0) call finish('cads_init_ifc', &
         'Failed to get major version from version="'//trim(version)//'"')
    i = i+2
    if (i <= l) then
       read(version(i:),*,iostat=stat) version_minor
       if (stat /= 0) call finish('cads_init_ifc', &
            'Failed to get minor version from version="'//trim(version)//'"')
    else
       version_minor = 1
    end if
    write(msg,*) version_major
    version = adjustl(msg)
    write(msg,*) version_minor
    version = trim(version)//'_'//adjustl(msg)
    if (dace%lpio) then
      ! Print info from CADS namelist
       write(*,'(3x,"version (processed)  = ",A)') trim(version)
       write(*,'(3x,"toplev               = ",F8.3," [hPa]")') toplev
       write(*,'(3x,"bottomlev            = ",F8.3," [hPa]")') bottomlev
       write(*,'(3x,"dwd_sort             = ",L1)') dwd_sort
       write(*,'(3x,"l2c_aer              = ",I1)') l2c_aer
       write(*,'(3x,"write_nhgt           = ",L1)') write_nhgt
       write(*,'(3x,"aer_out              = ",L1)') aer_out
       write(*,'(3x,"cld_out              = ",L1)') cld_out
       do ii = 1, n_sensor
         i = sensors(ii)
         write(*,'(3x,"instrument ",I3,":")') i
         write(*,'(5x,"use_all_chans      = ",L1)')   use_all_chans     (i)
         write(*,'(5x,"l_aerosol          = ",L1)')   l_aerosol         (i)
         write(*,'(5x,"l_land_sens        = ",L1)')   l_land_sens       (i)
         write(*,'(5x,"l_trace_gas        = ",L1)')   l_trace_gas       (i)
         write(*,'(5x,"use_bcor_cloud     = ",L1)')   use_bcor_cloud    (i)
         write(*,'(5x,"use_bcor_aerosol   = ",L1)')   use_bcor_aerosol  (i)
         write(*,'(5x,"use_bcor_trace_gas = ",L1)')   use_bcor_trace_gas(i)
         write(*,'(5x,"use_bcor_trace_gas = ",L1)')   use_bcor_trace_gas(i)
         write(*,'(5x,"surf_check         = ",I1)')   surf_check        (i)
         write(*,'(5x,"si_surf_check      = ",I1)')   si_surf_check     (i)
         write(*,'(5x,"si_thresh          = ",F6.3)') si_thresh         (i)
         write(*,'(5x,"use_im_frac        = ",I1)')   use_im_frac       (i)
         write(*,'(5x,"imager_cld_thresh  = ",F6.2)') imager_cld_thresh (i)
         write(*,'(5x,"imager_psurf_frac  = ",F6.2)') imager_psurf_frac (i)
       end do
    end if


    select case(version_major)
    case(2)
       ! CADS version 1 modified by DWD
      if (any(l_land_sens)) call finish('cads_init_ifc','l_land_sens not implemented for &
           &CADS version < 3')
      if (any(l_trace_gas)) call finish('cads_init_ifc','l_trace_gas not implemented for &
           &CADS version < 3')
      call cloud_detect_setup_hss(stat, l_aerosol=l_aerosol)
      if (stat /= 0) call finish('cads_init_ifc','cloud_detect_setup_hss failed.')

    case(3)
      select case(version_minor)
      case(1)
        if (dace%lpio) then
          call cads_setup_cloud_31
          if (stat /= 0) call finish('cads_init_ifc','cads_setup_cloud (3.1) failed.')
          if (any(l_aerosol)) then
            call cads_setup_aerosol_31
            if (stat /= 0) call finish('cads_init_ifc','cads_setup_aerosol (3.1) failed.')
          end if
          if (any(l_land_sens)) then
            call cads_setup_land_sensitivity_31
            if (stat /= 0) call finish('cads_init_ifc','cads_setup_land_sensitivity (3.1) failed.')
          end if
          if (any(l_trace_gas)) then
            call cads_setup_trace_gas_31
            if (stat /= 0) call finish('cads_init_ifc','cads_setup_trace_gas (3.1) failed.')
          end if
        end if
        call p_bcast_cloud_31(dace%pio)
        if (any(l_aerosol  )) call p_bcast_aerosol_31  (dace%pio)
        if (any(l_land_sens)) call p_bcast_land_sens_31(dace%pio)
        if (any(l_trace_gas)) call p_bcast_trace_gas_31(dace%pio)
        dwd_sort_31 = dwd_sort
        l_rank_p_31 = (l2c_aer >= 3)
        usd_31      = usd

      case default
        write(msg,*) version_minor
        msg = 'CADS minor version '//adjustl(trim(msg))//' not implemented'
        call finish('cads_init_ifc',trim(msg))
      end select
    case default
      write(msg,*) version_major
      msg = 'CADS major version '//adjustl(trim(msg))//' not implemented'
      call finish('cads_init_ifc',trim(msg))
    end select

    if (version_major /= 2) then
      ! Modify channels to be used
      do i = 1, size(sensors)
        isens = sensors(i)
        if (use_all_chans(isens)) then
          ! Get Band information from TOVS_OBS_CHAN namelists (stored in rad_set)
          do iset = 1, n_set
            rs => rad_set(iset)
            if (rs% id < 0) cycle
            if (any(rs% instr(1:rs%n_instr) == isens)) then
              do ii = 1, rs%n_instr
                if (rs% instr(ii) == isens) exit
              end do
              istart = rs% o_ch_i(ii) + 1
              iend   = rs% o_ch_i(ii) + rs% n_ch_i(ii)
              nbands = maxval(rs%band(istart:iend))
              do j = 1, nbands
                nchans(j) = count(rs%band(istart:iend) == j)
                chans(1:nchans(j), j) = pack(rs% chan(istart:iend), &
                     rs%band(istart:iend) == j)
              end do
              select case(version)
              case('3_1')
                call set_cloud_detect_chans_31(isens, nbands, nchans, chans)
              case default
                call finish('cads_init_ifc','set_cloud_detect_chans not implemented for given CADS version.')
              end select
              exit
            end if
          end do ! rad_set
        end if ! use_all chans
      end do ! sensors
    end if

    init = .true.

  end subroutine cads_init_ifc


  subroutine cads_ifc(rttov_id, nchans, chans, band, b, o_bc, o, chan2lev, plev,  &
       ttovs, instr, l_im_frac, psurf, flag, cloud_lev, lat, lon, l2c_max, si_max,&
       status, lprint, aer_type, im_flag)
    integer,       intent(in)            :: rttov_id
    integer,       intent(in)            :: nchans
    integer,       intent(in)            :: chans(:)
    integer,       intent(in)            :: band(:)
    real(kind=wp), intent(in)            :: b(:)
    real(kind=wp), intent(in),  target   :: o_bc(:)
    real(kind=wp), intent(in),  target   :: o(:)
    real(kind=wp), intent(in),  target   :: chan2lev(:)
    real(kind=wp), intent(in)            :: plev(:)
    type(t_tovs),  intent(inout)         :: ttovs
    integer,       intent(in)            :: instr
    logical,       intent(in)            :: l_im_frac
    real(kind=wp), intent(in)            :: psurf
    integer,       intent(out)           :: flag(:)
    real(kind=wp), intent(out)           :: cloud_lev(:)
    real(kind=wp), intent(in)            :: lat ! only for debugging (aer_out)
    real(kind=wp), intent(in)            :: lon ! only for debugging (aer_out)
    real(kind=wp), intent(out)           :: l2c_max
    real(kind=wp), intent(out)           :: si_max
    integer,       intent(out)           :: status
    logical,       intent(in),  optional :: lprint
    integer,       intent(out), optional :: aer_type
    integer,       intent(out), optional :: im_flag

    character(len=8),    parameter :: proc = 'cads_ifc'
    real(kind=wp),       parameter :: cl_huge = huge(cl_huge)

    type(t_tovs_instr)             :: ti
    type(t_rad_set),     pointer   :: rs => null()
    character(len=300)             :: msg = ''
    real(kind=wp),       pointer   :: o_p  (:) => null()
    real(kind=wp),       pointer   :: l2c_p(:) => null()
    real(kind=wp)                  :: p, cloud_lev_min
    integer                        :: nch, nlev
    integer                        :: ind(size(chans))
    integer                        :: flag_aux(size(flag))
    integer                        :: aerosol_type
    integer                        :: i_toplev, i_bottomlev, i, ic
    logical                        :: lpr
    real(kind=wp)                  :: land_frac
    ! Writing of cads_nhgt.agr (write_nhgt)
    integer,             parameter :: np = 200
    real(wp)                       :: p_min, p_max, lp, d_lp, nhgt, pd(1)
    integer                        :: iunit
    ! Writing of debug output
    type(t_aer_out_31)             :: ae_out
    type(t_cld_out_31)             :: cl_out
    character(len=30)              :: fname
    character(len=6)               :: position
    ! surf_check
    logical                        :: l_sc_l2c, l_sc_si
    logical,           allocatable :: mask(:)
    ! imager aided cloud detection
    integer                        :: im_flag_
    logical                        :: l_im
    integer                        :: i_imch_v, ii, ich, imcl
    integer                        :: n_im_ch,  n_im_cl, nc
    integer                        :: im_chan   (ttovs%n_im_ch)
    real(kind=wp)                  :: im_cl_frac(ttovs%n_im_cl)
    real(kind=wp)                  :: im_cl_mean(ttovs%n_im_ch, ttovs%n_im_cl)
    real(kind=wp)                  :: im_cl_stdv(ttovs%n_im_ch, ttovs%n_im_cl)
    real(kind=wp)                  :: im_o_stdv (ttovs%n_im_ch)
    real(kind=wp)                  :: im_o_mean (ttovs%n_im_ch)
    real(kind=wp)                  :: im_fg     (ttovs%n_im_ch)
    real(kind=wp)                  :: im_emis   (ttovs%n_im_ch)
    real(kind=wp)                  :: sfrac

    status = 0

    nlev = size(plev)
    if (size(chans   ) /= nchans) call finish(proc, 'invalid size of array "chans"')
    if (size(b       ) /= nchans) call finish(proc, 'invalid size of array "b"')
    if (size(o_bc    ) /= nchans) call finish(proc, 'invalid size of array "o_bc"')
    if (size(o       ) /= nchans) call finish(proc, 'invalid size of array "o"')
    if (size(chan2lev) /= nchans) call finish(proc, 'invalid size of array "chan2lev"')
    if (size(flag    ) /= nchans) call finish(proc, 'invalid size of array "flag"')

    if (present(lprint)) then
      lpr = lprint
    else
      lpr = .false.
    end if
    ! lpr = .true.
    ! usd = 0

    if (l_first) then
      position = 'rewind'
    else
      position = 'append'
    end if

    if (.not.init) call finish(proc, 'mo_cads_ifc is not initialized. Call cads_init_ifc first.')

    i_topLev    = minloc(abs(plev - toplev   *100._wp),dim=1) ! RTTOV level that is closest to toplev
    i_bottomLev = minloc(abs(plev - bottomlev*100._wp),dim=1) ! RTTOV level that is closest to bottomlev

    land_frac = real(ttovs%land_frac, kind=wp)

    cloud_lev(:) = cl_huge

    l2c_max = huge(l2c_max)
    si_max  = huge(si_max )
    ind = (/ (i, i=1, nchans) /)
    nch = nchans
    l_sc_l2c = (ttovs%rt_stype(1) == 0 .and. btest(surf_check   (rttov_id), SC_LAND)) .or. &
               (ttovs%rt_stype(1) == 1 .and. btest(surf_check   (rttov_id), SC_SEA )) .or. &
               (ttovs%rt_stype(1) == 2 .and. btest(surf_check   (rttov_id), SC_ICE ))
    l_sc_si  = (ttovs%rt_stype(1) == 0 .and. btest(si_surf_check(rttov_id), SC_LAND)) .or. &
               (ttovs%rt_stype(1) == 1 .and. btest(si_surf_check(rttov_id), SC_SEA )) .or. &
               (ttovs%rt_stype(1) == 2 .and. btest(si_surf_check(rttov_id), SC_ICE ))
    rs => null()
    if (l_sc_l2c .or. l_sc_si) then
      allocate(mask(nchans))
      mask = .true.
      if (l_sc_l2c) then
        do i = 2, nlev
          if (lpr) write(usd,*) i,plev(i),psurf
          if (plev(i) >= psurf) then
            l2c_max = i-1 + (psurf - plev(i-1))/(plev(i) - plev(i-1))
            if (lpr) write(usd,*) 'l2c_max',l2c_max
            exit
          end if
        end do
        l2c_max = min(l2c_max, real(nlev-1,wp))  ! nlev-1 is crucial for calc. on model levels, where
                                                 ! l2c is always above surface.
        where (chan2lev(:) >= l2c_max) mask(:) = .false.
      end if
      if (l_sc_si) then
        if (.not.associated(ttovs%sinfl)) call finish(proc, 'ttovs%surf_infl not available.&
             & si_surf_check requires surf_infl_mode to be set.')
        call get_tovs_rs(ttovs, rs=rs, ti=ti)
        do i = ti%o_ch_i(instr)+1, ti%o_ch_i(instr)+ti%n_ch_i(instr)
          if (ttovs%sinfl(i) > si_thresh(rttov_id)) then
            where(chans(:) == rs%chan(ttovs%ci(i))) mask = .false.
          end if
        end do
        si_max = si_thresh(rttov_id)
      end if

      nch = count(mask)
      if (lpr) write(usd,*) 'cads_use',nch,nchans,l2c_max,psurf,ttovs%rt_stype,lat,lon,nlev,plev(1),plev(nlev)
      ind(1:nch) = pack(ind, mask=mask)
    end if

    n_im_ch =  ttovs%n_im_ch
    n_im_cl =  ttovs%n_im_cl
    l_im = (iand(ttovs%init,TTOVS_IM) == TTOVS_IM) .and. (min(n_im_cl,n_im_ch) > 0)
    if (l_im) then
      if (lpr) write(usd,*) 'imager n_im*:',n_im_ch,n_im_cl
      im_chan(:) = 0
      nc         = 0
      do i = 1, ttovs%n_im_cl_v
        imcl = mod(ttovs%im_cl_v(i), 100)
        ich  = ttovs%im_cl_v(i) / 100
        if (ich > 0) then
#if defined(__GFORTRAN__) && (__GNUC__ <= 8)    /* FINDLOC is Fortran 2008.  */
          do ii = 1, nc
             if (im_chan(ii) == ich) exit
          end do
          if (ii > nc) ii = 0
#else
          ii   = findloc(im_chan(1:nc), ich, dim=1)
#endif
          if (ii <= 0) then
            nc = nc + 1
            if (nc > size(im_chan)) call finish(proc, 'inconsistent im_chan and ttovs%im_cl_v.')
            im_chan(nc) = ich
            ii  = nc
          end if
          if (imcl == IMCL_MEAN) then
            im_cl_mean(ii, 1:n_im_cl) = ttovs%im_cl(i,1:n_im_cl)
          elseif (imcl == IMCL_STDV) then
            im_cl_stdv(ii, 1:n_im_cl) = ttovs%im_cl(i,1:n_im_cl)
          end if
        else if (imcl == IMCL_FRAC) then
          im_cl_frac(    1:n_im_cl) = ttovs%im_cl(i,1:n_im_cl)
        end if
      end do
      call get_im_ch_ind(ttovs, IMCH_MEAN_B, i_imch_v)
      im_fg(1:n_im_ch) = ttovs%im_ch(i_imch_v,1:n_im_ch)
      call get_im_ch_ind(ttovs, IMCH_EMIS, i_imch_v)
      im_emis(1:n_im_ch) = ttovs%im_ch(i_imch_v,1:n_im_ch)
      sfrac = sum(im_cl_frac(1:n_im_cl))
      if (lpr) write(usd,*) 'imager cl frac ',im_cl_frac(1:n_im_cl),'sum=',sfrac
      ! Calculate Overall stdv (and mean)
!NEC$ nomove
      do i = 1, n_im_ch
        if (lpr) then
          write(usd,*) 'imager cl mean',i,im_chan(i),im_cl_mean(i,1:n_im_cl)
          write(usd,*) 'imager cl stdv',i,im_chan(i),im_cl_stdv(i,1:n_im_cl)
        end if
        if (sfrac > 1.E-10) then
          im_o_mean(i) = dot_product(im_cl_frac(1:n_im_cl),im_cl_mean(i,1:n_im_cl)) / sfrac
          im_o_stdv(i) = dot_product(im_cl_frac(1:n_im_cl),(im_cl_mean(i,1:n_im_cl) + im_cl_stdv(i,1:n_im_cl))**2) + &
                         dot_product(im_cl_frac(1:n_im_cl),(im_cl_mean(i,1:n_im_cl) - im_cl_stdv(i,1:n_im_cl))**2)
          im_o_stdv(i) = im_o_stdv(i) / (2.*sfrac)
          im_o_stdv(i) = im_o_stdv(i) - im_o_mean(i)**2
          if (im_o_stdv(i) > 1.E-20_wp) then
            im_o_stdv(i) = sqrt(im_o_stdv(i))
          else
            im_o_stdv(i) = 0._wp
          end if
        else
          im_o_mean(i) = 0._wp
          im_o_stdv(i) = 0._wp
        end if
        if (lpr) then
          write(usd,*) 'imager obs mean/stdv',i,im_chan(i),im_o_mean(i),im_o_stdv(i)
          write(usd,*) 'imager fg           ',i,im_chan(i),im_fg(i)
        end if
      end do
      call get_im_ch_ind(ttovs, IMCH_MEAN_O, i_imch_v)
      ttovs%im_ch(i_imch_v,1:n_im_ch) = im_o_mean(1:n_im_ch)
      call get_im_ch_ind(ttovs, IMCH_STDV  , i_imch_v)
      ttovs%im_ch(i_imch_v,1:n_im_ch) = im_o_stdv(1:n_im_ch)
    end if

    im_flag_ = 0
    if (btest(use_im_frac(rttov_id), IMFRC_CADS)) then
      if (l_im_frac) then
        if (ttovs%cloud_imag > imager_cld_thresh(rttov_id)) im_flag_ = ibset(im_flag_, 3) ! bits 0..2 by CADS*_Detect_Cloud_Imager
        if (lpr) write(usd,*) 'imager frac flag ',ttovs%cloud_imag,imager_cld_thresh(rttov_id),im_flag_
      else
        call finish(proc,'imager cloud fraction required by CADS but not available.')
      end if
    end if

    select case(version_major)
    case(2)
      l_debug_2 = lpr ! Does not work properly
      call cloud_detect_hss(rttov_id,             &
                            nch,                  &
                            chans(ind(1:nch)),    &
                            b(ind(1:nch)),        &
                            o_bc(ind(1:nch)),     &
                            chan2lev(ind(1:nch)), &
                            flag(1:nch),          &
                            i_toplev,             &
                            i_bottomlev,          &
                            status)
      flag(ind(1:nch)) = flag(1:nch)
      if (allocated(mask)) where(.not.mask(1:nchans)) flag = ibset(flag, MNW_SURF)
      do i = 1, nchans
        if (.not.btest(flag(i), MNW_CLOUD)) cycle
        if (band(i) <= 0 .or. band(i) > ubound(cloud_lev,1)) then
          write(0,*) i,'chan',chans(i),'band',band(i),'ubound(cloud_lev)',ubound(cloud_lev,1)
          call finish(proc,'invalid band or size of cloud_lev')
        end if
        cloud_lev(band(i)) = min(cloud_lev(band(i)), chan2lev(i))
      end do
    case(3)
      select case(version_minor)
      case(1)
        flag = 0

        if (.not.check_sensor_id_31(rttov_id, typ='cloud')) &
             call finish('cads_ifc '//trim(version), 'invalid sensor for cloud detection')
        if (use_bcor_cloud(rttov_id)) then
          o_p => o_bc
        else
          o_p => o
        end if
        l_debug_31 = lpr
        cl_out%l_out = cld_out
        call cads_detect_cloud_31(rttov_id,             &
                                  nch,                  &
                                  chans(ind(1:nch)),    &
                                  i_toplev,             &
                                  i_bottomlev,          &
                                  n_im_ch,              &
                                  im_chan,              &
                                  n_im_cl,              &
                                  flag_aux(1:nch),      &
                                  o_p(ind(1:nch)),      &
                                  b(ind(1:nch)),        &
                                  chan2lev(ind(1:nch)), &
                                  im_cl_frac,           &
                                  im_cl_mean,           &
                                  im_o_stdv,            &
                                  im_fg,                &
                                  Cloud_lev,            &
                                  K__Imager_Flag=im_flag_,&
                                  cld_out=cl_out)
        where(flag_aux(1:nch) /= 0) flag(ind(1:nch)) = ibset(flag(ind(1:nch)), MNW_CLOUD)
        if (allocated(mask)) where(.not.mask(1:nchans)) flag = ibset(flag, MNW_SURF)
        if (cld_out) then
          iunit = get_unit_number()
          write(fname, '("cld_out_",I4.4,".dat")') dace%pe
          open(iunit, file=fname, position=position)
          write(iunit,'(I3,2(1x,F9.4),1x,4(1x,I4))') &
               rttov_id, lat, lon, cl_out%i__chan_low, cl_out%i__chan_high, cl_out%i__scenario_index, cl_out%i__start_channel
          close(iunit)
          if (l_im) then
            write(fname, '("cld_out_imag_",I4.4,".dat")') dace%pe
            open(iunit, file=fname, position=position)
            write(iunit,'(I3,2(1x,F9.4),1x,100(1x,E13.6))') &
                 rttov_id, lat, lon, im_o_mean(1:n_im_ch), im_o_mean(1:n_im_ch)-im_fg(1:n_im_ch), &
                 im_emis(1:n_im_ch),cl_out%Z__Wsqdev, &
                 im_o_stdv(1:n_im_ch), cl_out%Z__intercluster_max, cl_out%Z__intercluster_sqdev, sfrac
            close(iunit)
          end if
          call return_unit_number(iunit)
        end if

        if (l_aerosol(rttov_id)) then
          if (.not.check_sensor_id_31(rttov_id, typ='aerosol')) &
               call finish('cads_ifc '//trim(version), 'invalid sensor for aerosol detection')
          if (use_bcor_aerosol(rttov_id)) then
            o_p => o_bc
          else
            o_p => o
          end if
          select case(l2c_aer)
          case(0)
            l2c_p => chan2lev
          case(1:)
            allocate(l2c_p(nchans))
            l2c_p = chan2lev
            call lev2p(plev, l2c_p)
            if (l2c_aer == 1) l2c_p = log(l2c_p)
          end select
          ae_out%l_out = aer_out
          call cads_detect_aerosol_31(rttov_id,     &
                                      nchans,       &
                                      chans,        &
                                      aerosol_type, &
                                      flag_aux,     &
                                      land_frac,    & ! TODO: land fraction
                                      o_p,          &
                                      l2c_p,        &
                                      aer_out=ae_out)
          where(flag_aux /= 0) flag = ibset(flag, MNW_AEROSOL)
          if (present(aer_type)) aer_type = aerosol_type
          if (l2c_aer > 0) deallocate(l2c_p)
          if (aer_out) then
            write(fname, '("aer_out_",I4.4,".dat")') dace%pe
            iunit = get_unit_number()
            open(iunit, file=fname, position=position)
            write(iunit,'(I3,2(1x,F9.4),1x,I1,1x,F6.3,3(1x,E13.6))') &
                 rttov_id, lat, lon, ae_out%aer_type, ae_out%hgt_thresh, &
                 ae_out%aod_other, ae_out%aod_dust, ae_out%aod_ash
            close(iunit)
            call return_unit_number(iunit)
          end if
        end if

        if (l_land_sens(rttov_id)) then
          if (.not.check_sensor_id_31(rttov_id, typ='land')) &
               call finish('cads_ifc '//trim(version), 'invalid sensor for land detection')
          call cads_detect_land_sensitivity_31(rttov_id,  &
                                               nchans,    &
                                               land_frac, & ! TODO: land fraction
                                               chan2lev,  &
                                               flag_aux)
          where(flag_aux /= 0) flag = ibset(flag, MNW_LAND)
        end if

        if (l_trace_gas(rttov_id)) then
          if (.not.check_sensor_id_31(rttov_id, typ='cloud')) &
               call finish('cads_ifc '//trim(version), 'invalid sensor for trace gas detection')
          if (use_bcor_trace_gas(rttov_id)) then
            o_p => o_bc
          else
            o_p => o
          end if
          call cads_detect_trace_gas_31(rttov_id, &
                                        nchans,   &
                                        chans,    &
                                        o_p,      &
                                        b,        &
                                        flag_aux)
          where(flag_aux /= 0) flag = ibset(flag, MNW_TRACEGAS)
        end if
      case default
        call finish(proc,'Version '//trim(version)//' not implemented')
      end select
    end select

    ! make sure, that all bands have a value for cloud_lev
    cloud_lev_min = minval(cloud_lev(:))
    where(cloud_lev(:) == cl_huge) cloud_lev(:) = cloud_lev_min

    ! Imager aided "false alarm" detection
    if (btest(use_im_frac(rttov_id),IMFRC_FAR)) then
      if (l_im_frac) then
        if (ttovs%soza <= 80._sp) then ! Use imager only during day
          if (ttovs%cloud_imag >= 0._sp .and. ttovs%cloud_imag < imager_cld_thresh(rttov_id)) then
            ! 1. Determine the bottom (imager_psurf_frac)% of the surface pressure
            p = (1._wp - imager_psurf_frac(rttov_id)) * psurf
            ! 2. Identify at which or between which two RTTOV levels this lies.
            !    Choose the upper RTTOV level (higher in the atmosphere) as a threshold
            do ic = ubound(plev,1), lbound(plev,1), -1
              if (plev(ic) < p) exit
            end do
            ! 3. Declare channels with a cloud_lev below this (lower in the atm.)
            !    as cloudfree
            do i = 1, nchans
              if (cloud_lev(band(i)) >= ic .and. cloud_lev(band(i)) /= cl_huge) then
                flag(i) = ibclr(flag(i), MNW_CLOUD)
              end if
            end do
            cloud_lev = huge(cloud_lev) ! indicate "no clouds" to calling routine
          end if
        end if
      else
        call finish(proc,'imager cloud fraction required by CADS (FAR) but not available.')
      end if
    end if

    if (present(im_flag)) im_flag = im_flag_

    if (lpr) then
!    if (.true.) then
      write(usd,*)
      write(usd,*) 'cads input/output'
      write(usd,*) 'instrument:',rttov_id
      write(usd,*) 'nchans    :',nchans
      write(usd,*) 'cloud_lev :',cloud_lev
      write(usd,*) 'lat/lon   :',lat,lon
      write(usd,*) 'psurf     :',psurf
      if (present(im_flag)) write(usd,*) 'im_flag   :',im_flag
      allocate(l2c_p(nchans))
      l2c_p = chan2lev
      call lev2p(plev, l2c_p)
      do i = 1, nchans
        write(msg,'(I5,1x,I5,3x,F6.2,2x,F10.2,2x,3(1x,F7.3),1x,I9,4(1x,L1))') i, chans(i),chan2lev(i), &
             l2c_p(i),b(i),o_bc(i),o(i), flag(i),btest(flag(i),MNW_CLOUD),btest(flag(i),MNW_AEROSOL),&
             btest(flag(i),MNW_LAND),btest(flag(i),MNW_TRACEGAS)
        if (allocated(mask)) write(msg(80:),'(L1)') mask(i)
        write(usd,*) trim(msg)
      end do
      do i = 1, nlev
        write(usd,*) 'lev/p',i,plev(i)
      end do
      deallocate(l2c_p)
    end if

    if (write_nhgt .and. lpr) then
      pd(1)  = minval(chan2lev)
      call lev2p(plev, pd) ; p_min = pd(1)
      pd(1) = maxval(chan2lev)
      call lev2p(plev, pd) ; p_max = pd(1)
      p_max = min(p_max, 101325._wp)
      d_lp = (log(p_max) - log(p_min)) / (np - 1._wp)
      !print*,'nhgt',p_min,p_max,log(p_min),log(p_max)
      iunit = get_unit_number()
      open(iunit, file='cads_nhgt.agr')
      write(iunit,'("@target G0.S0")')
      write(iunit,'("@type xy")')
      write(iunit,'("@s0 legend  ""ECMWF 137 model levels""")')
      do i = 1, np
        lp = log(p_min) + (i-1) * d_lp
        p = exp(lp)
        !print*,'nhgt',i,lp,log(p_min),log(p_max),p
        call nhgt2p(nhgt, p, p_min, p_max, lback=.true.)
        write(iunit,'(f11.4,1x,f7.5)') p*0.01_wp, nhgt
      end do
      write(iunit,'("&")')
      write(iunit,'("@target G0.S1")')
      write(iunit,'("@type xy")')
      write(iunit,'("@s1 legend  ""DWD ",I3," levels""")') size(plev)
      lp = log(p_min)
      do i = 1, np
        p = exp(lp)
        call nhgt2p(nhgt, p, p_min, p_max, pl=plev, lback=.true.)
        write(iunit,'(f11.4,1x,f7.5)') p*0.01_wp, nhgt
        lp = lp + d_lp
      end do
      close(iunit)
      call return_unit_number(iunit)
    end if

    l_first = .false.

  end subroutine cads_ifc

  include "nhgt2p.incf"


end MODULE mo_cads_ifc
