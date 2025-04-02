!
!+ Interface to the RTTOV fast radiative transfer model package
!
MODULE mo_rttov
!
! Description:
!   - Definition of the datatype 't_rad' to hold observations of
!     top of the atmosphere radiances.
!   - Definition of datatype 't_rttov_p' and 't_rttov_s' for atmospheric
!     profiles and single level data passed to rttov.
!   - Routines to read radiance observations from a NetCDF file.
!   - Interface to the RTTOV fast radiative transfer model package.
!
! Current Maintainer: DWD, Robin Faulwetter
!    phone: +49 69 8062 2746
!    fax:   +49 69 8062 3721
!    email: robin.faulwetter@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_1         2008/11/05 Andreas Rhodin
!  First operational 3D-Var release
! V1_4         2009/03/26 Andreas Rhodin
!  Cope with fillvalues in 1dvar feedback file
! V1_7         2009/08/24 Andreas Rhodin
!  Prepare for NOAA-19
! V1_8         2009/12/09 Detlef Pingel
!  Interface to RTTOV9
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Detlef Pingel
!  changes for HIRS, IASI, SATPP
! V1_19        2012-04-16 Harald Anlauf
!  bugfixes: RTTOV9 adjoint code (conversion ppmv ->  kg/kg)
!            for IASI + minimal channel selection
!            for optional arguments
!            fix disambiguation bug for rttov_init_min=.true.
!  read coefficients on 1PE
!  remove support for RTTOV-6
!  changes for RTTOV10                            (A.Messer)
! V1_22        2013-02-13 Robin Faulwetter
!  Implementation of vectorized K-mode.
!  loading of new style SATPP files. (Superobing).
!  preparations for RTTOV-10 51 Layer coefficients.
!  Provide satellite/sun azimut angle to rttov for FASTEM.
! V1_23        2013-03-26 Robin Faulwetter
!  Implemented processing of CrIS data
! V1_26        2013/06/27 Robin Faulwetter
!  Introduced a check on the influence of the surface onto radiances.
!  Introduced USE_MWSURF bit. Corrected the usage of other USE_* bits.
! V1_27        2013-11-08 Andreas Rhodin
!  fix bug: x% tb(i1:i2) = emissiv(i1:i2,1)
! V1_28        2014/02/26 Robin Faulwetter
!  Enable compilation without RTTOV7
! V1_29        2014/04/02 Andreas Rhodin
!  introduce sink variable snow-fraction for IR emissivities
! V1_31        2014-08-21 Harald Anlauf
!  enable 53/54 levels with RTTOV10
!  preparations for emissivity adjoint calculations
! V1_35        2014-11-07 Robin Faulwetter
!  Enabled processing of chinese satellite radiances
! V1_47        2016-06-06 Robin Faulwetter
!  Many improvements for radiances.
! V1_48        2016-10-06 Robin Faulwetter
!  Implemented RTTOV12
! V1_51        2017-02-24 Harald Anlauf
!  Cleanup of (deprecated) RTTOV7 stuff
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Andreas Rhodin    DWD  2001  original code
!                        2003  adapted to new file format / 3D-Var
! Detlef Pingel     DWD  2004  change to RTTOV version 7
!                        2004  write artificial data
! Christina Koepken DWD  2005  corrections
! Christina Koepken DWD  2007  introduce NOAA18,METOP
! Christina Koepken DWD  2008  introduce AMSU-B,MHS
!==============================================================================
#if defined(__ICON__) && !defined(__USE_RTTOV)
#undef _RTTOV_VERSION
#endif
  !=============
  ! modules used
  !=============
  !------------------------
  ! general purpose modules
  !------------------------
  use mo_kind,          only: wp,                &! working precision kind
                              sp                  ! single  precision kind
  use mo_exception,     only: finish              ! abort on error condition
  use mo_mpi_dace,      only: dace,              &! MPI group info
                              self                ! MPI group info (COMM_SELF)
  use mo_dace_string,   only: char2               ! transform integer to string
  use mo_physics,       only: RDRD,              &! gas const. (dry air)/gas const. (water)
                              ppmv_moist_q,      &! ppmv_moist from spec.hum.
                              ppmv_dry_q,        &! ppmv_moist from spec.hum.
                              dppmv_moist_dq,    &! d(ppmv_moist)/dq
                              dppmv_dry_dq,      &! d(ppmv_dry)/dq
                              q_ppmv_dry,        &! spec.hum. from ppmv_dry
                              q_ppmv_dry_trg,    &! tracegas spec.conc. from ppmv_dry
                              ppmv_moist_trg      ! tracegas ppmv_moist from spec.conc.
#if (_RTTOV_VERSION >= 12)
  use rttov_const,      only: gas_unit_specconc,   &!
                              gas_unit_ppmvdry,    &!
                              gas_unit_ppmv
#endif

#if (_RTTOV_VERSION > 0)
  use rttov_const,      only: qmin_rttov => qmin,  &! minimum allowed water vapour
                              qmax_rttov => qmax,  &! maximum allowed water vapour
                              tmin_rttov => tmin,  &! minimum allowed temperature
                              tmax_rttov => tmax    ! maximum allowed temperature

#endif
  use mo_rtifc,         only: rtifc_check_config,  &! Check RTTOV version, level number, set nlevs_top
                              rtifc_set_opts,      &! set RTTOV options
                              rtifc_get_opts,      &! get RTTOV options
                              rtifc_init,          &! initializes the rttov modules
                              rtifc_coef_prop,     &! get RTTOV coefficient properties
                              rtifc_fill_input,    &! fills variable part for call of rttov
                              rtifc_direct,        &! renews instr. dependent part of rttov_direct
                              rtifc_k,             &! renews instr. dependent part of rttov_k
                              rtifc_errmsg,        &! rttov interface error messages
                              rtifc_version,       &! rttov library version string
                              rtifc_vers,          &! rtifc branch (RTTOV version, it is compiled for)
                              read1pe,             &! Read/distrib. coeffs.
                              nlevs_top,           &! Additional rttov levels
                              NO_ERROR, WARN_RTTOV,&! exit codes
                              qmin_ifc,            &! minimum humidity
                              qmax_ifc              ! maximum humidity

#if (_RTTOV_VERSION >= 12)
  use mo_rtifc,         only: default_gas_units,   &!
                              rtifc_emis_atlas,    &! mw emissivity from atlas
                              rtifc_emis_retrieve, &! mw emissivity dynamic retrieval
                              rtifc_emis_sea,      &! mw emissivity sea model
                              rtifc_tskin_retrieve,&
                              rtifc_brdf_atlas
#endif
  use mo_rad,           only: t_rad_set,           &!
                              destruct,            &!
                              m_chan,              &!
                              rad_set,             &!
                              n_set,               &!
                              lev2chan
  use mo_satid,         only: satid_bufr2rttov
  use mo_t_tovs,        only: t_tovs,              &!
                              t_tovs_instr,        &!
                              get_tovs_rs
  use mo_t_obs,         only: usd,dpref

  implicit none
!------------------------------------------------------------------------------
  !================
  ! public entities
  !================
  private
  !---------------------
  ! rttov version to use
  !---------------------
  public :: version
  !---------------------------
  ! RTTOV interface data types
  !---------------------------
  public :: t_rttov_prof    ! data type holding RTTOV input profile information
  public :: construct       ! construct t_rttov_var
  public :: destruct        ! destruct t_rttov_var, t_rad, t_rad_set
  !-------------------------
  ! RTTOV interface routines
  !-------------------------
  public :: set_rttov_vers  ! set RTTOV version to be used by this module
  public :: clean_mo_rttov  ! deallocate module variables
  public :: call_rttvi      ! initialize rttov package
  public :: call_rttov      ! call rttov routine
  public :: rttov_bounds    ! apply rttov limits to variables
  public :: q2rt_humi       ! conversion from spec.humi. to RTTOV input humidity unit
  public :: rt_trg          ! conversion of ozone units for RTTOV
  public :: get_preslev     ! get pressure levels
  !-------------------------
  ! Options
  !-------------------------
  public :: rt_humi         ! Humidity input switch
  public :: rt_levs_exact   ! Use exact RTTOV levels (from RTTOV directly)
  public :: rt_hard_limits  ! check on RTTOV hard limits
  !------------------------------------------------------------------
  ! rttov hard coded dimensions (no constants, set by set_rttov_vers)
  !------------------------------------------------------------------
  public :: jplev           ! no. of pressure levels
  public :: jpnav           ! no. of profile variables
  public :: jpnsav          ! no. of surface air variables
  public :: jpnssv          ! no. of skin variables
  public :: jpnsv           ! no. of solar variables
  public :: jpncv           ! no. of cloud variables
  public :: jpnsat          ! max. no. of ?
  !----------------------------------------------------------------
  ! hard coded dimensions for single level and dummy sink variables
  !----------------------------------------------------------------
  public :: nsv              ! # of single level input variables
  public :: nbs              ! # of bound constrained sink variables
  !-----------------------------
  ! min, max, reference profiles
  !-----------------------------
  public :: qmin, qmax      ! hard limits for Q
  public :: tmin, tmax      ! hard limits for T
  public :: preslev         ! pressure levels for rttov
  public :: preshPa         ! pressure levels in hPa
  public :: lnp             ! ln (preslev)
  !-------------------------
  ! values returned by RTTVI
  !-------------------------
  public :: rttvi_called
!------------------------------------------------------------------------------
  !================================
  ! Module variables and data types
  !================================

  integer, parameter :: jpnsat_def = 50
  integer, parameter :: jpnav_def  = 2
  integer, parameter :: jpnsav_def = 5
  integer, parameter :: jpnssv_def = 1
  integer, parameter :: jpncv_def  = 2
  ! hard coded dimensions for single level and dummy sink variables
  integer ,parameter :: nsv = 5 ! # of single level input variables
  integer ,parameter :: nbs = 3 ! # of bound constrained sink variables

  ! RTTOV
  integer, save      :: rt_humi        = 0 ! Switch for humidity input to RTTOV
                                           ! 0: humidity in ppmv over moist air with wrong conversion spec.humi.->ppmv
                                           ! 1: humidity in ppmv over moist air with fixed conversion spec.humi.->ppmv
                                           ! 2: humisity in ppmv over dry air (for use with new rttov11 coefficient files)
  logical, save      :: rt_levs_exact   = .true.
  integer, save      :: rt_hard_limits  = 2
  real(wp),save      :: eps_bound       = 0.0001_wp   ! modifier for RTTOV hard limits
  logical ,save      :: rttvi_called    = .false. ! rttov routine initialized?

  !--------------------------------------------
  ! Quantities set by subroutine set_rttov_vers
  !--------------------------------------------
  integer :: version = -999  ! RTTOV version to use
  integer :: jplev   = -999  ! no. of pressure levels
  integer :: jpnav   = -999  ! no. of profile variables
  integer :: jpnsav  = -999  ! no. of surface air variables
  integer :: jpnssv  = -999  ! no. of skin variables
  integer :: jpnsv   = -999  ! no. of solar variables
  integer :: jpncv   = -999  ! no. of cloud variables
  integer :: jpnsat  = -999  ! max. no. of ?
  !----------------------------------------------------------------------------


  !=====================================================================
  ! data type holding variable input/output to/from rttov/rttovad/rttovk
  !=====================================================================
  type t_rttov_prof
    real(wp) ,pointer :: av (:,:,:)=>NULL()  ! (jplev,jpnav ) atmospheric profile vars.
    real(wp) ,pointer :: sav  (:,:)=>NULL()  ! (jpnsav)       surface air variables
    real(wp) ,pointer :: ssv  (:,:)=>NULL()  ! (jpnssv)       surface skin variables
    real(wp) ,pointer :: cv   (:,:)=>NULL()  ! (jpncv )       cloud variables
    real(wp) ,pointer :: sv   (:,:)=>NULL()  ! (jpnsv )       solar variables
    real(wp) ,pointer :: emis (:,:)=>NULL()  ! (npc, nchan )  surface emissivities
    real(wp)          :: hsurf     = 0._wp   !                surface height
    real(wp)          :: snf       = 0._wp   !                snow fraction (for RTOVP.nc)
!    integer           :: stype     = 0       !                surface type
    integer           :: i_o3      = -1      ! index of ozone in av array
    integer           :: i_co2     = -1      ! index of co2 in av array
    logical           :: drtsfl    = .false. ! derived skin temperature flag
  end type t_rttov_prof

!------------------------------------------------------------------------------
  integer              :: nProf            ! number of profiles in one RTTOV call
  integer              :: err_rttov               ! RTTOV-interface error code
  !-------------------------------
  ! parameters returned from RTTVI
  !-------------------------------
  real(wp) ,allocatable :: preslev (:) !(jplev) pressure levels for rttov (Pa)
  real(wp) ,allocatable :: preshPa (:) !(jplev)    pressure levels in hPa
  real(wp) ,allocatable, target :: lnp (:) !(jplev) ln (preslev)
  real(wp)              :: tmin        ! min temp
  real(wp)              :: tmax        ! max temp
  real(wp)              :: qmin        ! min q
  real(wp)              :: qmax        ! max q
  !---------------
  ! NetCDF indices
  !---------------
  integer            :: status      ! NetCDF return variable from functions

!==============================================================================
  !-----------
  ! Interfaces
  !-----------

  interface construct
    module procedure construct_rttov_prof
  end interface construct

  interface destruct
    module procedure destruct_rttov_prof
  end interface destruct

!==============================================================================
contains
!==============================================================================
  subroutine set_rttov_vers (vers,levels)
  integer, intent(in) :: vers
  integer, intent(in) :: levels

  !--------------------------------
  ! set RTTOV version number
  ! set version dependent constants
  ! allocate some arrays
  !--------------------------------
    if (version > 0) then
      if (version /= vers) call finish ('set_rttov_vers',&
           'attempt to reset RTTOV version!')
      return
    end if

    !-------------------------
    ! set RTTOV version number
    !-------------------------
    version       = vers
    jplev         = levels

    select case (version)
    !-----------------------------------
    ! set constants from RTTOV version 7
    !-----------------------------------
    case (7)
      call finish ('set_rttov_vers', 'RTTOV7 not implemented anymore.')
    case (9)
      call finish ('set_rttov_vers', 'RTTOV9 not implemented anymore.')
    case (10,12,13)
      jpnav   = jpnav_def    ! no. of profile variables
      jpnsav  = jpnsav_def   ! no. of surface air variables
      jpnssv  = jpnssv_def   ! no. of skin variables
      jpncv   = jpncv_def    ! no. of cloud variables
      jpnsat  = jpnsat_def   ! max. no. of sensors to be used
      jpnsv   = 2            ! no of solar variables - lat/lon
    case default
      call finish ('set_rttov_vers','RTTOV version not supported: '//char2(version))
    end select

    call rtifc_check_config(version, levels, status)
    if (status /= NO_ERROR) then
       call finish('set_rttov_vers', &
                   trim(rtifc_errmsg(status))//' '//char2(version)//&
                   ' - binary was compiled for version '//&
                   char2(rtifc_vers))
    end if

    !--------------------------------
    ! allocate some arrays
    !--------------------------------
    allocate (preslev (jplev))       ! pressure levels for rttov (Pa)
    allocate (preshPa (jplev))       ! pressure levels in hPa
    allocate (lnp     (jplev))       ! ln (preslev)

  end subroutine set_rttov_vers
!------------------------------------------------------------------------------
  subroutine clean_mo_rttov
  !------------------
  ! deallocate arrays
  !------------------
    if (version <=0) return
    deallocate (preslev   ) ! 43 pressure levels for rttov (Pa)
    deallocate (preshPa   ) !    pressure levels in hPa
    deallocate (lnp       ) ! ln (preslev)

  end subroutine clean_mo_rttov
!------------------------------------------------------------------------------
  subroutine construct_rttov_prof (var, nch, jakobi, npc, n_av, nlev)
  type (t_rttov_prof),intent(out) :: var    ! rttov variables data type
  integer            ,intent(in)  :: nch    ! number of channels
  integer            ,intent(in)  :: jakobi ! 0:forward,1:jacobian,2:kch=2
  integer,  optional ,intent(in)  :: npc    ! # of PC (for jakobi==2)
  integer,  optional ,intent(in)  :: n_av   ! # of profile variables
  integer,  optional ,intent(in)  :: nlev   ! # of profile variables

  integer :: nav, nlev_
  !---------------------------------------------------------------------------
  ! allocate components of derived type t_rttov_prof :
  ! jakobi = 0: for forward model
  !          1: for Jakobi matrix
  !          2:     Jakobi matrix for 3DVAR to RTTOV input variable conversion
  !---------------------------------------------------------------------------
    if (present(n_av)) then
      nav = n_av
    else
      nav = jpnav
    end if
    if (present(nlev)) then
      nlev_ = nlev
    else
      nlev_ = jplev
    end if

    select case (jakobi)

    case (0)                                   ! forward calculations:
      allocate (var% av  (nlev_ ,nav    ,  1)) ! atmospheric profile vars.
      allocate (var% sav        (jpnsav ,  1)) ! surface air variables
      allocate (var% ssv        (jpnssv ,  1)) ! surface skin variables
      allocate (var% sv         (jpnsv  ,  1)) ! solar variables
      allocate (var% cv         (jpncv  ,  1)) ! cloud variables
      allocate (var% emis       (    1  ,nch)) ! surface emissivities

    case (1)                                   ! Jakobi matrix: d tb / d ...
      allocate (var% av  (nlev_ ,nav    ,nch)) ! atmospheric profile vars.
      allocate (var% sav        (jpnsav ,nch)) ! surface air variables
      allocate (var% ssv        (jpnssv ,nch)) ! surface skin variables
      allocate (var% sv         (jpnsv  ,nch)) ! solar variables
      allocate (var% cv         (jpncv  ,nch)) ! cloud variables
      allocate (var% emis       (    1  ,nch)) ! surface emissivities

    case (2)                                   ! Jakobi: d 3dvar / d rttov
      if(.not.present(npc)) call finish('construct_rttov_prof','npc')
      allocate (var% av  (nlev_ ,nav    ,  2)) ! atmospheric profile vars.
      allocate (var% sav        (jpnsav ,  2)) ! surface air variables
      allocate (var% ssv        (jpnssv ,  1)) ! surface skin variables
      allocate (var% sv         (jpnsv  ,  1)) ! solar variables
      allocate (var% cv         (jpncv  ,  1)) ! cloud variables
      allocate (var% emis       (  npc  ,nch)) ! surface emissivities

    case default
      call finish('construct_rttov_prof','jakobi not in valid range (0:2)')
    end select

  end subroutine construct_rttov_prof
!------------------------------------------------------------------------------
  subroutine destruct_rttov_prof (var)
  type (t_rttov_prof) ,intent(inout) :: var
    if (.not. associated (var% av)) return
    deallocate (var% av  )
    deallocate (var% sav )
    deallocate (var% ssv )
    deallocate (var% sv  )
    deallocate (var% cv  )
    deallocate (var% emis)
  end subroutine destruct_rttov_prof
! !------------------------------------------------------------------------------
  subroutine call_rttvi (data_path)
    character(len=*)  ,intent(in) :: data_path  ! Path to coefficient files
    !--------------------------------------------
    ! call the rttov initialization routine RTTVI
    !--------------------------------------------
    integer              :: n_instr                 !
    integer              :: rt_instr(3,jpnsat)      ! (1,1:nInstr): rttov platform   id's;
                                                    ! (2,1:nInstr): rttov satellite  id's;
                                                    ! (3,1:nInstr): rttov instrument id's;
                                                    ! nInstr: summarized number of all instruments
                                                    ! on all satellites which should be processed;
    integer              :: numchans   (jpnsat)     ! number of chans per instrument
    integer              :: channels(m_chan,jpnsat) ! list of channels to extract (channels,1:instr)
                                                    ! for each instrument
    integer              :: chan_idx(m_chan,jpnsat) ! Channel index list
    integer              :: iopts(jpnsat)
    integer              :: mx_chans

    integer, parameter      :: invalid = -huge(0)
    integer                 :: i,j,k            ! indices
    integer                 :: ioff,nch
    integer                 :: nl
    integer                 :: mpi_comm_type    ! MPI communicator for rtifc_init
    integer                 :: n_proc           ! # procs. in group
    integer                 :: my_proc_id       ! rank of current processor
    integer                 :: io_proc_id       ! rank of I/O processor
    type(t_rad_set),pointer :: s                ! dataset description
    integer                 :: instr, iset      ! loop indices
    integer                 :: platf, sat       ! auxiliary variables
    character(len=5)        :: cchan (jpnsat)   ! channel for printout
    logical                 :: l_new

    !--------------------------------------------------
    ! currently rtifc_init is called with fixed setup:
    ! NOAA-15, NOAA-16, NOAA-18, METOP-2, AQUA, NOAA-19
    !--------------------------------------------------

    ! TODO: simplify: maybe, combine set_rttov_vers and call_rttvi
    if (version < 0) call finish('call_rttvi', 'set_rttov_vers not called so far.')

#if (_RTTOV_VERSION > 0)

    ! initialize parameter arrays
    rt_instr = 0
    numchans = 0
    channels = invalid

    i = 0
    set_loop: do iset = 1, n_set
      s => rad_set(iset)
      call satid_bufr2rttov(s%satid, sat, platf)
      if (sat < 0 .or. platf < 0) call finish('call_rttvi', 'unknown satellite')
      s% rttov_satid = sat
      s% platform    = platf

      instr_loop: do instr = 1, s% n_instr
        ioff = s%o_ch_i(instr)
        nch  = s%n_ch_i(instr)
        l_new = .true.
        if (s%rttov_indx(instr) > 0 .and. associated(s%chidx)) then
          l_new = any(s%chidx(ioff+1:ioff+nch) <= 0)
        end if
        if (l_new) then
          i = i + 1
          if (i > jpnsat) call finish('call_rttvi', 'jpnsat too small')
          if (s%rttov_indx(instr) <= 0) call finish('call_rttvi','invalid rttov options &
               &index, maybe rtifc_set_opts(new=T) was not called.')
          iopts     (       i) = s%rttov_indx(instr)
          rt_instr  (1    , i) = platf
          rt_instr  (2    , i) = sat
          rt_instr  (3    , i) = s% instr(instr)
          numchans  (       i) = nch
          channels  (1:nch, i) = s%chan(ioff+1:ioff+nch)
        end if
      end do instr_loop
      if (s% im_rttov_indx > 0) then
        i = i+ 1
        if (i > jpnsat) call finish('call_rttvi', 'jpnsat too small')
        nch                  = s% n_im_chan
        iopts     (       i) = s% im_rttov_indx
        rt_instr  (1    , i) = platf
        rt_instr  (2    , i) = sat
        rt_instr  (3    , i) = s% im_rttov_id
        numchans  (       i) = nch
        channels  (1:nch, i) = s% im_chan(1:s%n_im_chan)
      end if
    end do set_loop
    n_instr = i
    mx_chans = maxval (numchans   (1:n_instr))
    if (n_instr <= 0) RETURN

    !-----------
    ! call RTTVI
    !-----------
    !------------
    ! print input
    !------------
    if (dace% lpio) then
       write(6,'(a)')repeat('-',79)
       write(6,*)
       write(6,'(a,I2,a)'  ) 'Input to '//rtifc_version()//':'
       write(6,*)
       write(6,'(a,(100i5))') 'n_instr    :',n_instr
       write(6,'(a,(100i5))') 'platform   :',rt_instr(1,1:n_instr)
       write(6,'(a,(100i5))') 'satellite  :',rt_instr(2,1:n_instr)
       write(6,'(a,(100i5))') 'instrument :',rt_instr(3,1:n_instr)
       write(6,'(a,(100i5))') 'numchans   :',numchans(  1:n_instr)
       write(6,*)
       do j = 1, mx_chans
         do i = 1, n_instr
           if (channels(j,i) /= invalid) then
             write(cchan(i),'(i5)') channels (j,i)
           else
             cchan(i) = ''
           end if
         end do
         write(6,'(a,(100a5))') 'channel    :', cchan(1:n_instr)
       end do
       write(6,*)
    endif

    mpi_comm_type = self% comm  ! MPI communicator for rtifc_init
    n_proc        = 1           ! # processors in comm. group
    my_proc_id    = 0           ! rank of current processor
    io_proc_id    = 0           ! rank of processor performing the I/O

    !-------------------------------------
    ! Read coeffs on I/O pe and distribute
    !-------------------------------------
    if (read1pe) then
      mpi_comm_type = dace% comm
      n_proc        = dace% npe
      my_proc_id    = dace% pe
      io_proc_id    = dace% pio
    end if

    ! call rtifc_init:
    call rtifc_init(                         &
         rt_instr     (1:3       ,1:n_instr),& ! <- info on processed instruments
         channels     (1:mx_chans,1:n_instr),& ! <- info on processed channels/instrument
         numchans     (           1:n_instr),& ! <- number of channels/instr. defined in channels
         chan_idx     (1:mx_chans,1:n_instr),&
         iopts        (           1:n_instr),& ! <- RTTOV options
         my_proc_id,                         & ! <- rank of current processor
         n_proc,                             & ! <- number of procs. in comm. group
         io_proc_id,                         & ! <- rank of processor performing the I/O
         mpi_comm_type,                      & ! <- MPI communicator
         status        = err_rttov,          &
         path_coefs    = data_path           ) ! <- path to coefficient files
    if (err_rttov /= 0) then
       if (dace% lpio) print*,trim (rtifc_errmsg (err_rttov, .true.))
       call finish('call_rttvi (rtifc_init)', trim (rtifc_errmsg (err_rttov, .true.)))
    endif

    ! Store channel indices s%chidx
    ! and get number of rttov levels
    i = 0
    do iset = 1, n_set
      s => rad_set(iset)
      if (.not.associated(s%chidx)) then
        allocate(s%chidx(s%n_chan))
        s%chidx = -1
      end if
      do instr = 1, s%n_instr
        ioff = s%o_ch_i(instr)
        nch  = s%n_ch_i(instr)
        i = i+1
        s%chidx(ioff+1:ioff+nch) = chan_idx(1:nch, i)
        call rtifc_coef_prop(s%rttov_indx(instr), nlevs=nl)
        s%iopts(instr)%rt_nlevs = nl
        if (nl > jplev+1 .or. nl < jplev) &
             call rtifc_set_opts(s%rttov_indx(instr), addinterp=.true.)
      end do
      if (s% im_rttov_indx > 0) then
        i = i+ 1
        s% im_chan_indx(1:s%n_im_chan) = chan_idx(1:s%n_im_chan, i)
      end if
    end do

    ! Get rttov pressure levels
    preslev = -1._wp
    if (rt_levs_exact) then
      do k = 1, n_instr
        call rtifc_coef_prop(iopts(k), nlevs=nl)
        if (jplev <= nl .and. jplev >= nl-1) then
          call rtifc_coef_prop(iopts(k), preslev=preslev)
          exit
        end if
      end do
    end if
    if (any(preslev < 0._wp)) then
      call get_preslev(jplev, preslev)
    end if

    preshPa = preslev
    preslev = preslev * 100._wp  ! hPa -> Pa
    lnp     = log (preslev)

    !---------------------------
    ! set dummy qmin,qmax values
    !---------------------------
    if (rt_hard_limits <= 1) then
      ! Old, wrong version: use (wrong) q limits in wrong unit [ppmv_dry] from RTTOV as bounds for qv [kg/kg]
      ! Kept only for backwards compatibilty
      qmin    = 1.e-10_wp
      qmax    = 1._wp
      qmin_ifc = qmin_rttov * 2
    else
      qmin    = q_ppmv_dry(qmin_rttov) * (1._wp + eps_bound)
      qmax    = q_ppmv_dry(qmax_rttov) * (1._wp - eps_bound)
      qmin_ifc = q2rt_humi(qmin)
      qmax_ifc = q2rt_humi(qmax)
    end if
    tmin = tmin_rttov
    tmax = tmax_rttov
    if (dace% lpio .and. rt_hard_limits > 0) then
       write(6,*)
       write(6,'(A20,4(2x,A13))') 'RTTOV hard limits:','tmin','tmax','qmin','qmax'
       write(6,'(20x,4(2x,E13.6))') tmin,tmax,qmin,qmax
     end if

#endif

    rttvi_called = .true.

  end subroutine call_rttvi


  subroutine get_preslev(nlev, p)
    integer,  intent(in) :: nlev
    real(wp), intent(out) :: p(:)

    if (size(p) < nlev) call finish('get_preslev', 'array p too small')

    select case (nlev)
    case (43) ! RTTOV9 levels and RTTOV10 levels (in RTTOV-9 compat mode,
              ! using RTTOV9 coefficients)
      p =                                                                 &
       (/           0.10,   0.29,   0.69,   1.42,  2.611,  4.407,   6.95, &
           10.37,  14.81,  20.40,  27.26,  35.51,  45.29,  56.73,  69.97, &
           85.18, 102.05, 122.04, 143.84, 167.95, 194.36, 222.94, 253.71, &
          286.60, 321.50, 358.28, 396.81, 436.95, 478.54, 521.46, 565.54, &
          610.60, 656.43, 702.73, 749.12, 795.09, 839.95, 882.80, 922.46, &
          957.44, 985.88,1005.43,1013.25  /)
    case (44) ! RTTOV10 levels (using RTTOV9 coefficients)
      p =                                                                 &
       (/  0.005,   0.10,   0.29,   0.69,   1.42,  2.611,  4.407,   6.95, &
           10.37,  14.81,  20.40,  27.26,  35.51,  45.29,  56.73,  69.97, &
           85.18, 102.05, 122.04, 143.84, 167.95, 194.36, 222.94, 253.71, &
          286.60, 321.50, 358.28, 396.81, 436.95, 478.54, 521.46, 565.54, &
          610.60, 656.43, 702.73, 749.12, 795.09, 839.95, 882.80, 922.46, &
          957.44, 985.88,1005.43,1013.25  /)
    case (50) ! RTTOV-10 levels (in RTTOV-9 compat mode)
      p =                                                                  &
       (/             0.01,    0.10,    0.20,     0.50,    0.80,    1.20,  &
             1.60,    2.20,    2.70,    3.50,     4.20,    5.00,    6.95,  &
            10.37,   14.81,   20.40,   27.26,    35.51,   45.29,   56.73,  &
            69.97,   85.18,  102.05,  122.04,   143.84,  167.95,  194.36,  &
           222.94,  253.71,  286.60,  321.50,   358.38,  396.81,  436.95,  &
           478.54,  521.46,  565.54,  610.60,   656.43,  702.73,  749.12,  &
           795.09,  839.95,  882.80,  922.46,   957.44,  985.88, 1005.43,  &
          1025.00, 1050.00 /)
    case (51) ! RTTOV-10 levels
      p =                                                                  &
       (/    0.005,   0.01,    0.10,    0.20,     0.50,    0.80,    1.20,  &
             1.60,    2.20,    2.70,    3.50,     4.20,    5.00,    6.95,  &
            10.37,   14.81,   20.40,   27.26,    35.51,   45.29,   56.73,  &
            69.97,   85.18,  102.05,  122.04,   143.84,  167.95,  194.36,  &
           222.94,  253.71,  286.60,  321.50,   358.38,  396.81,  436.95,  &
           478.54,  521.46,  565.54,  610.60,   656.43,  702.73,  749.12,  &
           795.09,  839.95,  882.80,  922.46,   957.44,  985.88, 1005.43,  &
          1025.00, 1050.00 /)
    case (53) ! RTTOV-10 levels (in RTTOV-9 compat mode)
      p =                                                                  &
       (/             0.0131,    0.0304,    0.0644,    0.1263,    0.2324,  &
           0.4052,    0.6749,    1.0801,    1.6691,    2.5011,    3.6462,  &
           5.1864,    7.2150,    9.8368,   13.1672,   17.3308,   22.4601,  &
          28.6937,   36.1735,   45.0430,   55.4433,   67.5109,   81.3744,  &
          97.1505,  114.9420,  134.8320,  156.8850,  181.1390,  207.6090,  &
         236.2780,  267.1010,  300.0000,  334.8650,  371.5530,  409.8890,  &
         449.6680,  490.6520,  532.5770,  575.1540,  618.0710,  660.9960,  &
         703.5860,  745.4840,  786.3280,  825.7550,  863.4050,  898.9280,  &
         931.9850,  962.2590,  989.4510, 1013.2900, 1033.5400, 1050.0000  /)
    case (54) ! RTTOV-10 levels
      p =                                                                  &
       (/  0.0050,    0.0131,    0.0304,    0.0644,    0.1263,    0.2324,  & ! 1- 6
           0.4052,    0.6749,    1.0801,    1.6691,    2.5011,    3.6462,  & ! 7-12
           5.1864,    7.2150,    9.8368,   13.1672,   17.3308,   22.4601,  & !13-18
          28.6937,   36.1735,   45.0430,   55.4433,   67.5109,   81.3744,  & !19-24
          97.1505,  114.9420,  134.8320,  156.8850,  181.1390,  207.6090,  & !25-30
         236.2780,  267.1010,  300.0000,  334.8650,  371.5530,  409.8890,  & !31-36
         449.6680,  490.6520,  532.5770,  575.1540,  618.0710,  660.9960,  & !37-42
         703.5860,  745.4840,  786.3280,  825.7550,  863.4050,  898.9280,  & !43-48
         931.9850,  962.2590,  989.4510, 1013.2900, 1033.5400, 1050.0000  /) !48-54
    case (101) ! RTTOV-10 levels
      p =                                                                  &
       (/  0.0050,    0.0161,    0.0384,    0.0769,    0.1370,    0.2244,    0.3454,    0.5064, &
           0.7140,    0.9753,    1.2972,    1.6872,    2.1526,    2.7009,    3.3398,    4.0077, &
           4.9204,    5.8776,    6.9567,    8.1655,    9.5119,   11.0038,   12.6492,   14.4559, &
          16.4318,   18.5847,   20.9224,   23.4526,   26.1829,   29.1210,   32.2744,   35.6505, &
          39.2566,   43.1001,   47.1882,   51.5278,   56.1260,   60.9895,   66.1253,   71.5398, &
          77.2396,   83.2310,   89.5204,   96.1138,  103.0170,  110.2370,  117.7780,  125.6460, &
         133.8460,  142.3850,  151.2660,  160.4960,  170.0780,  180.0180,  190.3200,  200.9890, &
         212.0280,  223.4420,  235.2340,  247.4080,  259.9690,  272.9190,  286.2620,  300.0000, &
         314.1370,  328.6750,  343.6180,  358.9660,  374.7240,  390.8930,  407.4740,  424.4700, &
         441.8820,  459.7120,  477.9610,  496.6300,  515.7200,  535.2320,  555.1670,  575.5250, &
         596.3060,  617.5110,  639.1400,  661.1920,  683.6670,  706.5650,  729.8860,  753.6280, &
         777.7900,  802.3710,  827.3710,  852.7880,  878.6200,  904.8660,  931.5240,  958.5910, &
         986.0670, 1013.9500, 1042.2300, 1070.9200, 1100.0000   /)


    case default
      call finish('call_rttvi', 'failed to determine pressure levels')
    end select
  end subroutine get_preslev


  ! TODO: remove this routine. This is partly implemented in mo_rtifc.
  subroutine rttov_bounds (x)
    type(t_rttov_prof) ,intent(inout)  :: x                      ! vars. to check
    !-------------------------------------------------------------------
    ! Check if the RTTOV input variables (x) are within specified bounds.
    !-------------------------------------------------------------------
    integer :: k
    if (rt_hard_limits <= 0) return
    if (x% sav(1,1) < tmin) then
      x% sav(1,1) = tmin
    endif
    if (x% sav(1,1) > tmax) then
      x% sav(1,1) = tmax
    endif
    if (x% sav(2,1) < qmin) then
      x% sav(2,1) = qmin
    endif
    if (x% sav(2,1) > qmax) then
      x% sav(2,1) = qmax
    endif

    do k=1,size(x% av,1)
      if (x% av (k,1,1) < tmin) then
        x% av (k,1,1) = tmin
      endif
      if (x% av (k,1,1) > tmax) then
        x%  av (k,1,1) = tmax
      endif
      if (x% av (k,2,1) < qmin) then
        x%  av (k,2,1) = qmin
      endif
      if (x% av (k,2,1) > qmax) then
        x%  av (k,2,1) = qmax
      endif
    end do
  end subroutine rttov_bounds
!==============================================================================
  elemental function q2rt_humi(q) result(rth)
    real(wp)             :: rth
    real(wp), intent(in) :: q

#if (_RTTOV_VERSION >= 12)
    select case(default_gas_units)
    case(gas_unit_specconc)
      ! Nothing to do
      rth = q
    case(gas_unit_ppmv)
      ! ppmv over moist air
      rth = ppmv_moist_q(q)
    case(gas_unit_ppmvdry)
      ! ppmv over dry air (required for RTTOV10 with RTTOV11 coefficient files)
      rth = ppmv_dry_q(q)
    end select
#else
    if (rt_humi == 1) then
      ! ppmv over moist air
      rth = ppmv_moist_q(q)
    elseif (rt_humi == 2) then
      ! ppmv over dry air (required for RTTOV10 with RTTOV11 coefficient files)
      rth = ppmv_dry_q(q)
    else
      ! ppmv over moist air with wrong conversion (which was routinely used for a long time)
      rth = 1.d6/RDRD * q
    end if
#endif
  end function q2rt_humi
!==============================================================================
  elemental function rt_trg(trg,q,gas_id) result(trg_rt)
    real(wp)             :: trg_rt
    real(wp), intent(in) :: trg
    real(wp), intent(in) :: q
    integer,  intent(in) :: gas_id ! gas_id_* from rttov_const

#if (_RTTOV_VERSION >= 12)
    ! select case(default_gas_units)
    ! case(gas_unit_specconc)
    !   ! Nothing to do
    !   trg_rt = trg
    ! case(gas_unit_ppmv)
    !   ! ppmv over moist air
    !   trg_rt = ppmv_moist_trg(trg,q,gas_id)
    ! case(gas_unit_ppmvdry)
    !   ! ppmv over dry air
    !   trg_rt = ppmv_dry_trg(trg,q,gas_id)
    ! end select

    ! Input is in ppmv_dry
    select case(default_gas_units)
    case(gas_unit_specconc)
      ! Nothing to do
      trg_rt = q_ppmv_dry_trg(trg,q,gas_id)
    case(gas_unit_ppmv)
      ! ppmv over moist air
      trg_rt = q_ppmv_dry_trg(trg,q,gas_id)    ! ppmv_dry -> specconc
      trg_rt = ppmv_moist_trg(trg_rt,q,gas_id) ! specconc -> ppmv_moist
    case(gas_unit_ppmvdry)
      ! ppmv over dry air
      trg_rt = trg
    end select
#else
    trg_rt = -HUGE(0._wp)        ! undefined
#endif

  end function rt_trg

!==============================================================================
  subroutine call_rttov (ckey, tovs, x, ifail, msk, tb, x_a, ldebug, obs_bt, atlas, styp, &
                         max_dst, ierr_no_stop, spt_hd_id)
    character,         intent(in)                      :: ckey     ! f=forward a=adjoint k=jacobi
    type(t_tovs),      intent(in)                      :: tovs      ! constant input
    type(t_rttov_prof),intent(inout)                   :: x          ! variables to/from rttov
    integer,           intent(inout)                   :: ifail
    logical,           intent(in),    optional         :: msk(:)
    real(wp),          intent(out),   optional         :: tb(:)
    type(t_rttov_prof),intent(inout), optional, target :: x_a        ! Jacobi matrices
    logical,           intent(in),    optional         :: ldebug
    real(wp),          intent(in),    optional         :: obs_bt(:)
    integer,           intent(in),    optional         :: atlas
    integer,           intent(in),    optional         :: styp
    real(wp),          intent(in),    optional         :: max_dst
    integer,           intent(in),    optional         :: ierr_no_stop(:)
    integer,           intent(in),    optional         :: spt_hd_id

    ! meta information
    type(t_rad_set),     pointer     :: rs => null()
    type(t_tovs_instr)               :: ti

    ! Indices
    integer                          :: iopt                ! options index
    integer                          :: ksen                ! sensor    index
    integer                          :: chidx(m_chan)       ! Channel indices
    integer                          :: k,nc,i1,i2          ! indices
    real(wp),            allocatable :: q2m_rttv    (:)     ! tmp. 2m humidity, unit rttv
    integer ,            allocatable :: lProfs      (    :) ! rttov profile index
    real(wp),            allocatable :: press       (:,  :) ! used pressure grid (profile)
    real(wp),            allocatable :: humi_rttv   (:,:  ) ! tmp. humidity profile, unit ppmv
    real(wp),            allocatable :: emissiv_k   (  :,:) ! k-matrix of surface emissivities
    real(wp),            allocatable :: emissiv     (  :,:) ! surface emissivities
    real(wp),            allocatable :: dr_tskin    (    :) ! surface derived skin temperature
    integer ,            allocatable :: emissiv_flag(  :,:) ! surface emissivities flags #LB
    real(wp),            allocatable :: temp_k      (:,:,:) ! gradient matrix temperature profiles
    real(wp),            allocatable :: humi_k      (:,:,:) ! gradient matrix humidity profiles
    real(wp),            allocatable :: t2m_k       (  :,:) ! gradient of 2m temperature
    real(wp),            allocatable :: psurf_k     (  :,:) ! gradient of surface pressure
    real(wp),            allocatable :: q2m_k       (  :,:) ! gradient of 2m humidity
    real(wp),            allocatable :: u10m_k      (  :,:) ! gradient of 10m wv
    real(wp),            allocatable :: v10m_k      (  :,:) ! gradient of 10m wv
    real(wp),            allocatable :: sTemp_k     (  :,:) ! gradient of surface skin temperature
    real(wp),            allocatable :: ctp_k       (  :,:) ! gradient resp. cloud top press.
    real(wp),            allocatable :: cfraction_k (  :,:) ! gradient resp. cloud fraction
    real(wp),            allocatable :: T_b         (  :,:) ! brightness temperature
    real(wp),            allocatable :: radTotal    (  :,:) ! total radiance
    real(wp),            allocatable :: radOverc    (:,:,:) ! overcast radiance
    real(sp),            allocatable :: levChans    (  :,:) ! assignment of channels to levels
    real(wp),            allocatable :: spec        (  :,:)
    integer,             allocatable :: channum     (:)
    integer                          :: ipr
    integer                          :: styp_
    logical                          :: ldeb,lkeep_ozone_data,lkeep_co2_data
    logical                          :: lstop

    ifail = 0

#if (_RTTOV_VERSION <= 0)
    return
#endif

    if (any(ckey == (/'e','a','b','s','t'/))) then
      if (.not.present(msk)) &
           call finish('call_rttov','mask missing for ckey='//ckey)
    elseif (ckey == 'k') then
      if (present(msk)) &
           call finish('call_rttov','mask with wrong ckey')
    end if
    if (ckey == 'a' .and. .not.present(atlas)) &
           call finish('call_rttov','ckey=a without atlas ID.')



    if (present(ldebug)) then
      ldeb = ldebug
    else
      ldeb = .false.
    end if
    !---------------------------------------------------
    ! 1. Setup: check arguments, array allocations
    !---------------------------------------------------
    if (all(ckey /= (/'e', 'a', 'b', 's', 't'/))) call rttov_bounds (x)

    call get_tovs_rs(tovs, rs=rs, ti=ti)

    !-----------------------------------------------
    ! 2. Call RTTOV
    !-----------------------------------------------
    ! 2.1 General choice: calculate TBs
    !-----------------------------------------------
    !----------------------------------------------------------------------
    ! 2.2.1 Loop over all possible sensors
    !       Setup of channels, setup of current satid, sensor, surface params
    !----------------------------------------------------------------------

     ! so far, pass only 1 profile to RTTOV
     nProf = 1

     ! allocate some temporary arrays for rtifc_fill_input:
     allocate(press    (tovs%nlev    ,nProf))
     allocate(humi_rttv(tovs%nlev    ,nProf ))
     allocate(q2m_rttv (              nProf ))

     ! fill temporary arrays for rtifc_fill_input:
     if (tovs%i_p > 0) then
       do k = 1, nProf
         press(:,k) = tovs%av(:,tovs%i_p) * 0.01_wp ! fill pressure profiles for all nProf profiles
       enddo
     else
       do k = 1, nProf
         press(:,k) = preshPa(:) ! fill pressure profiles for all nProf profiles
       enddo
     end if

     ! convert units of humidity
     humi_rttv = q2rt_humi(x%av (:,2,:))
     q2m_rttv  = q2rt_humi(x%sav(  2,:))

     if (present(styp)) then
       styp_ = styp
     else
       styp_ = tovs%rt_stype(1)
     end if

     ! fill the rttov profile structure array
     ! (for processing more than one profile at once)
     call rtifc_fill_input (        &
          status     = err_rttov,   & ! -->  exit status
          iopts      = rs%rttov_indx(ti%ii(1:ti%n_instr)), &!
          press      = press(:,:),  & ! <--  used pressure grid (profile) [hPa]
          temp       = x%av(:,1,:), & ! <--  temperature profiles         [K]
          humi       = humi_rttv,   & ! <--  water vapour profile      [ppmv]
          t2m        = x%sav(1,:),  & ! <--  2m temperature               [K]
          q2m        = q2m_rttv,    & ! <--  2m humidity               [ppmv]
          psurf      = x%sav(3,:),  & ! <--  surface pressure           [hPa]
          hsurf      = (/x%hsurf/), & ! <--  surface height (elevation)  [km]
          u10m       = x%sav(4,:),  & ! <--  10m U wind component       [m/s]
          v10m       = x%sav(5,:),  & ! <--  10m V wind component       [m/s]
          stemp      = x%ssv(1,:),  & ! <--  surface skin temperature     [K]
          stype      = (/styp_/),   & ! <--  surface type: 0=land, 1=sea, 2=sea-ice
          lat        = x%sv(1,:),   & ! <--  latitude of ground point   [deg]
          lon        = x%sv(2,:),   & ! <--  longitude of ground point   [deg]
          sat_zen    = (/real(tovs%saza,wp)/),   & ! <--  satellite zenith angle     [deg]
          ctp        = x%cv (1,:),  & ! <--  cloud top pressure         [hPa]
          cfraction  = x%cv (2,:),  & ! <--  cloud fraction            [0..1]
          sat_azi    = (/real(tovs%boa,wp)/),    & ! <--  satellite azimuth angle     [deg]
          sun_zen    = (/real(tovs%soza,wp)/),   & ! <--  sun zenith angle           [deg]
          sun_azi    = (/real(tovs%soa,wp)/),    & ! <--  sun azimuth angle           [deg]
          snw_frc    = (/x%snf/) )

         lstop = (err_rttov /= NO_ERROR .and. err_rttov /= WARN_RTTOV)
         if (lstop .and. present(ierr_no_stop)) lstop = all(err_rttov /= ierr_no_stop)
         if (lstop) call finish('call_rttov (rtifc_fill_input)',   &
                      trim (rtifc_errmsg (err_rttov, .true.)))

     do ksen = 1, ti%n_instr

       i1     = ti% o_ch_i(ksen) + 1
       i2     = ti% o_ch_i(ksen) + ti%n_ch_i(ksen)

       ! number of channels and their position i1,i2
       ! in arrays for current sensor:
       if (present(msk)) then
         nc = count(msk(i1:i2))
       else
         nc = ti% n_ch_i(ksen)
       end if

       if (nc == 0) cycle

       iopt = rs%rttov_indx(ti%ii(ksen))

       allocate(emissiv(nc, nProf)) ! calculated emissivities
       if (present(msk)) then
         chidx  (1:nc)   = pack(rs%chidx(tovs%ci(i1:i2)), mask=msk(i1:i2))
         emissiv(1:nc,1) = pack(x%emis(1,i1:i2), mask=msk(i1:i2))
       else
         chidx(1:nc)     = rs%chidx(tovs%ci(i1:i2))
         emissiv(1:nc,1) = x%emis(1,i1:i2)
       end if
       if (rs%iopts(ti%ii(ksen))%do_lambertian) then
         allocate(spec(nc, nProf))
         k     = count(rs%iopts(1:ti%ii(ksen))%do_lambertian)
         spec  = min(max(real(tovs%spec(k),wp),0._wp),1._wp)
       else
         allocate(spec(0,0))
       end if

       allocate(lProfs (nc       )) ! index list of profiles
       lProfs=1

       ! Temporarily disable O3, CO2, since we do not have the input fields here
       ! (It does not affect the emissivity calcs. done with this module)
       call rtifc_get_opts(iopt=iopt, ozone_data=lkeep_ozone_data, co2_data=lkeep_co2_data)
       call rtifc_set_opts(iopt=iopt, ozone_data=.false.,          co2_data=.false.       )

      !--------------
      ! 2.2.3 Do RT
      !--------------
      select case (ckey)
      case('e','a','b','s','t')
         !--------------------------------------------------
         ! MW emissivity (dynamical retrievals and atlases)
         !-------------------------------------------------
#if (_RTTOV_VERSION >= 12)
         select case(ckey)
         case('e')
            call rtifc_emis_retrieve(iopt,lProfs, chidx(1:nc), pack(obs_bt(:), mask=msk), &
                                     spec, emissiv(:,1), err_rttov, pe=dace%pe, ldeb=ldeb)
        case('t')
            allocate(channum(nc),dr_tskin(nprof))
            channum(1:nc)   = pack(rs%chan(tovs%ci(i1:i2)), mask=msk(i1:i2))
            !print *, 'call_rttov channum ', channum(1:nc)
            call rtifc_tskin_retrieve(iopt,lProfs, channum, chidx(1:nc), pack(obs_bt(:), mask=msk), &
                                    spec, emissiv(1:nc,1), dr_tskin(:), err_rttov, x%drtsfl,        &
                                    pe=dace%pe, ldeb=ldeb, spt_hd_id=spt_hd_id)
            if (ldeb) write(usd,*) 'debug_spot call_rttov tskin_mdl ', x%ssv(1,:), ' tskin_dr ' , dr_tskin(:)
            if (all(dr_tskin(:) > 0)) x%ssv(1,:) = dr_tskin(:)
            deallocate(channum,dr_tskin)
         case('a')
           call  rtifc_emis_atlas(iopt, atlas, lProfs, chidx(1:nc), emissiv(:,1), err_rttov,&
                                  ldebug=ldebug, max_dst=max_dst)
         case('s')
           call  rtifc_emis_sea(iopt, lProfs, chidx(1:nc), emissiv(:,1), err_rttov,&
                                  ldeb=ldebug)
         case('b')
           allocate(emissiv_flag(nc, nProf)) ! emissivity flags
           call  rtifc_brdf_atlas(iopt, lProfs, chidx(1:nc), emissiv(:,1), err_rttov, emissiv_flag(:,1))
           if (ldeb) write(usd,*) 'debug_spot call_rttov emis_flag',emissiv_flag(:,1)
           where(emissiv_flag(:,1) >= 3 .and. emissiv_flag(:,1) <= 6) emissiv(:,1) = -1._wp
           deallocate(emissiv_flag)
         end select
#else
         call finish('call_rttov', 'RTTOV version does not provide emissivity/brdf atlases')
#endif
         if (ldeb) write(usd,*) 'debug_spot call_rttov emis',err_rttov,emissiv(:,1)

         x% emis(1,:) = unpack(emissiv(:,1), msk, x% emis(1,:))

         lstop = (err_rttov /= NO_ERROR .and. err_rttov /= WARN_RTTOV)
         if (lstop .and. present(ierr_no_stop)) lstop = all(err_rttov /= ierr_no_stop)
         if (lstop) call finish('call_rttov', trim (rtifc_errmsg (err_rttov)))

      case ('f')
      !--------------
      ! forward model
      !--------------
         ipr=0 !; if (dace% pe==3) ipr=1
         allocate(T_b    (nc, nProf)) ! calculated brightness temperatures
         ! call RTTOV direct:
         call rtifc_direct (                  &
              iopt       = iopt,              & ! <--  rttov instrument index
              lProfs     = lProfs,            & ! <--  list of profile indices
              chans      = chidx ( 1:nc),     & ! <--  list of channel indices
              emissiv    = emissiv( 1:nc,:),  & ! <--> emissivities -
              T_b        = T_b,               & !  --> calculated brightness
              status     = err_rttov,         & !  --> exit status
              iprint     = (/ipr/))

         ! store value of temporary variables:
         if ( present(tb) )   tb      (i1:i2)  = T_b     (1:nc,1)
         if (present(msk)) then
           x% emis(1,:) = unpack(emissiv(:,1), msk, x% emis(1,:))
         else
           x% emis (1,i1:i2)  = emissiv (1:nc,1)
         end if
         if (ldeb) write(usd,*) dpref,'call_rttov emis_out',x%emis(1,:)

         ! initialize variables/gradients unchanged by RTTOV
         ! x%av  ( :,3:,      1) = 0._wp

         lstop = (err_rttov /= NO_ERROR .and. err_rttov /= WARN_RTTOV)
         if (lstop .and. present(ierr_no_stop)) lstop = all(err_rttov /= ierr_no_stop)
         if (lstop) call finish('call_rttov(rtifc_direct)',trim(rtifc_errmsg(err_rttov)))

      case ('k')
        if (rs%gopts%lev_mode > 0) call finish('call_rttov', 'lev_mode > 0 currently not supported')

       ! allocate temporary variables for rtifc_k:
       allocate(emissiv_k(              i1:i2, nProf)) ! k-matrix of surface emissivities
       allocate(temp_k   (size(preshPa),i1:i2, nProf)) ! gradient matrix t
       allocate(humi_k   (size(preshPa),i1:i2, nProf)) ! gradient matrix wv
       allocate(t2m_k    (              i1:i2, nProf)) ! gradient of 2m t
       allocate(q2m_k    (              i1:i2, nProf)) ! gradient of 2m wv
       allocate(sTemp_k  (              i1:i2, nProf)) ! gradient of surf. skin t
       allocate(ctp_k    (              i1:i2, nProf)) ! gradient resp. cloud top pressure
       allocate(cfraction_k (           i1:i2, nProf)) ! gradient resp. cloud fraction
       allocate(psurf_k  (              i1:i2, nProf)) ! gradient of surf. pressure
       allocate(u10m_k   (              i1:i2, nProf)) ! gradient of 10m wv
       allocate(v10m_k   (              i1:i2, nProf)) ! gradient of 10m wv
       allocate(radTotal (              i1:i2, nProf)) ! total radiance    (-> cld. detect.)
       allocate(levChans (              i1:i2, nProf)) ! vert. localization of channels
       allocate(radOverc (jplev+nlevs_top-1,i1:i2, nProf))       ! overcast radiance (-> cld. detect.)
       allocate(T_b    (nc, nProf)) ! calculated brightness temperatures

       ! initialize inout variable emissiv_k
       emissiv_k = 0._wp

       ipr=0 !; if (dace% pe==0) ipr=1
       ! call RTTOV K:

       call rtifc_k (                         &
            iopt       = iopt,                & ! <--  rttov instrument index
            lProfs     = lProfs,              & ! <--  list of profile indices
            chans      = chidx     ( 1:nc),   & ! <--  list of channel indices
            emissiv    = emissiv   (i1:i2,:), & ! <--> emissivities -
            emissiv_k  = emissiv_k (i1:i2,:), & !  --> k-matrix of surface emissivities
            temp_k     = temp_k,              & !  --> gradient matrix t
            humi_k     = humi_k,              & !  --> gradient matrix wv
            ctp_k      = ctp_k,               & !  --> gradient resp. cloud top pressure
            cfraction_k= cfraction_k,         & !  --> gradient resp. cloud fraction
            t2m_k      = t2m_k,               & !  --> gradient of 2m t
            q2m_k      = q2m_k,               & !  --> gradient of 2m wv
            sTemp_k    = sTemp_k,             & !  --> gradient of surf. skin t
            T_b        = T_b,                 & !  --> calculated brightness temperatures
            status     = err_rttov,           & !  --> exit status
            psurf_k    = psurf_k,             & !  --> gradient of surf. pressure
            u10m_k     = u10m_k,              & !  --> gradient of 10m wv
            v10m_k     = v10m_k,              & !  --> gradient of 10m wv
            radTotal   = radTotal,            & !  --> total radiance (for cld. detect.)
            radOvercast= radOverc,            & !  --> overcast radiance (for cld. detect.)
            iprint     = (/ipr/))


#if (_RTTOV_VERSION >= 12)
       select case(default_gas_units)
       case(gas_unit_specconc)
         ! Nothing to do
       case(gas_unit_ppmv)
         ! ppmv over moist air
         humi_k(:,:,1) = humi_k(:,:,1) * dppmv_moist_dq(x%av (:,2,:))
         q2m_k   (:,1) = q2m_k   (:,1) * dppmv_moist_dq(x%sav(  2,:))
       case(gas_unit_ppmvdry)
         ! ppmv over dry air
         humi_k(:,:,1) = humi_k(:,:,1) * dppmv_dry_dq(x%av (:,2,:))
         q2m_k   (:,1) = q2m_k   (:,1) * dppmv_dry_dq(x%sav(  2,:))
       end select
#else
       if (rt_humi == 1) then
         ! ppmv over moist air
         humi_k(:,:,1) = humi_k(:,:,1) * dppmv_moist_dq(x%av (:,2,:))
         q2m_k   (:,1) = q2m_k   (:,1) * dppmv_moist_dq(x%sav(  2,:))
       elseif (rt_humi == 2) then
         ! ppmv over dry air (required for RTTOV10 with RTTOV11 coefficient files)
         humi_k(:,:,1) = humi_k(:,:,1) * dppmv_dry_dq(x%av (:,2,:))
         q2m_k   (:,1) = q2m_k   (:,1) * dppmv_dry_dq(x%sav(  2,:))
       else
         ! ppmv over moist air with wrong conversion (which was routinely used for a long time)
         humi_k(:,:,1) = humi_k(:,:,1) * 1.d6/RDRD
         q2m_k   (:,1) = q2m_k   (:,1) * 1.d6/RDRD
       end if
#endif


       ! store values of temporary varibles:
       x_a%emis(  1,i1:i2)  = emissiv_k(  :,1)
       x_a%av  (:,1,i1:i2)  = temp_k   (:,:,1)
       x_a%av  (:,2,i1:i2)  = humi_k   (:,:,1)
       x_a%sav (  1,i1:i2)  = t2m_k    (  :,1)
       x_a%sav (  2,i1:i2)  = q2m_k    (  :,1)
       x_a%sav (  3,i1:i2)  = psurf_k  (  :,1)
       x_a%sav (  4,i1:i2)  = u10m_k   (  :,1)
       x_a%sav (  5,i1:i2)  = v10m_k   (  :,1)
       x_a%ssv (  1,i1:i2)  = sTemp_k  (  :,1)
       x_a%cv  (  1,i1:i2)  = ctp_k    (  :,1)
       x_a%cv  (  2,i1:i2)  = cfraction_k(:,1)

       if ( present(tb) ) tb(i1:i2) = T_b(:,1)
       x%emis (  1,i1:i2)  = emissiv  (:,1)

       ! initialize variables/gradients unchanged by RTTOV
       ! x_a%ssv(   2:,i1:i2  ) = 0._wp

!        ! do level-to-channel assignment:
       ! TODO
       ! Currently, just one instrument is supprted. In a next step individual options
       ! can be provided for each instrument in each SatPP file and a loop over all
       ! instruments is required:
       ! rs%iopts(1)%l2c_abs_lim = ...iopts(instument)...
       if (rs%iopts(ti%ii(ksen))%l2c_type > 0) then
         if (rs%iopts(1)%l2c_use_rad) then
           call lev2chan(radTotal, radOverc,spread((/.true./),1,nc), levChans, &
                         rs%iopts(ti%ii(ksen))%l2c_type,                       &
                         rs%iopts(1)%l2c_rel_lim, rs%iopts(1)%l2c_abs_lim)
         else
           call lev2chan(T_b, radOverc,spread((/.true./),1,nc), levChans,      &
                         rs%iopts(ti%ii(ksen))%l2c_type,                       &
                         rs%iopts(1)%l2c_rel_lim, rs%iopts(1)%l2c_abs_lim)
         end if
!            if (jplev /= size(radOverc,1)) then
!                ! RTTOV10 is not running in rttov9 compat mode
!                ! We have to do something on this, however, for a first shot
!                ! just shift levChans up by one (radOverc has dimension
!                ! layers = levels - 1)
!                levChans(:,:) = levChans(:,:) + 1._sp
!            endif
         tovs%l2c(i1:i2) = real(levChans(i1:i2,1), kind=sp)

         if (ldeb) then
           do k = i1,i2
             write(usd,*) 'debug_spot store l2c:',k,rs%chan(tovs%ci(k)),tovs%l2c(k)
           end do
         end if
       endif

       if (ldeb) then
         do k = i1,i2
           write(usd,*) 'debug_spot rttov result',k,rs%chan(tovs%ci(k)),T_b(k,1)
         end do
       end if

       ! deallocate temporary variables for rtifc_k:
       deallocate (radTotal, radOverc, levChans)
       deallocate (emissiv_k, temp_k, humi_k, t2m_k, q2m_k, sTemp_k, &
                   psurf_k, u10m_k, v10m_k, ctp_k, cfraction_k       )

       lstop = (err_rttov /= NO_ERROR .and. err_rttov /= WARN_RTTOV)
       if (lstop .and. present(ierr_no_stop)) lstop = all(err_rttov /= ierr_no_stop)
       if (lstop) call finish('call_rttov (rtifc_k)',trim(rtifc_errmsg(err_rttov)))

      case default
        call finish('call_rttov:','invalid ckey='//ckey)
      end select

      call rtifc_set_opts(iopt=iopt, ozone_data=lkeep_ozone_data, co2_data=lkeep_co2_data)

      ! deallocate variables for
      ! rtifc_direct and rtifc_k:
      deallocate(lProfs, emissiv)

    enddo

    ! deallocate auxiliary quantities
    ! for rtifc_fill_input:
    deallocate(press,humi_rttv,q2m_rttv)

    !-------------------
    ! 2.3 copy output arrays
    !-------------------
    if (err_rttov /= NO_ERROR .and. err_rttov /= WARN_RTTOV) then
      ifail = err_rttov
    else
      ifail = 0
    end if

  end subroutine call_rttov
!==============================================================================
end module mo_rttov
