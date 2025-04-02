!
!+ Apply cloud detection for radiance data
!
MODULE mo_cloud
!
! Description:
!   Apply cloud detection for radiance data
!   3dvar top level cloud detection module
!
! Current Maintainer: DWD, Robin Faulwetter
!    phone: +49 69 8062 2746
!    fax:   +49 69 8062 3721
!    email: robin.faulwetter@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_13        2011/11/01 Detlef Pingel
!  3dvar interface to module mo_clouddetection
! V1_19        2012-04-16 Andreas Rhodin
!  cloud_detect: cleanup and bug fixes for IASI (band_size decr_use)
!  revised IASI diagnostics in feedback file
!  use flags flg_cld/flg_cld_.. for steering of AMSU-A cloud flags
! V1_20        2012-06-18 Andreas Rhodin
!  cleanup
! V1_21        2012-09-04 Andreas Messer
!  add support for ATMS
! V1_22        2013-02-13 Robin Faulwetter
!  implementation of vectorized K-mode; Add support for ATMS
! V1_23        2013-03-26 Robin Faulwetter
!  Implemented processing of CrIS data.
!  Improve bias correction modules: t_decay = 0. : infinite memory
!                                   t_decay < 0. : static bias correction.
! V1_26        2013/06/27 Robin Faulwetter
!  Introduced a check on the influence of the surface onto radiances.
!  Introduced USE_MWSURF bit.
!  Corrected the usage of other USE_* bits.
! V1_27        2013-11-08 Robin Faulwetter
!  Rewrite of radiance flagging (for assimilation of radiances over land)
! V1_28        2014/02/26 Robin Faulwetter
!  Bugfixes in MW cloud check; obs - fg based cloud checks for AMSUA and ATMS
! V1_35        2014-11-07 Robin Faulwetter
!  Enabled processing of chinese satellite radiances.
!  estimate cloud top height dummy sink variable first guess from MNW (E.Lange)
!  Fix for missing initialization (H.Anlauf)
! V1_43        2015-08-19 Robin Faulwetter
!  features for assimilation of high peaking radiances over land/clouds
! V1_47        2016-06-06 Robin Faulwetter
!  Many improvements for radiances
! V1_48        2016-10-06 Robin Faulwetter
!  reorganize interface to subroutine amsub_cloud_check
! V1_49        2016-10-25 Robin Faulwetter
!  Option to use mapped amsub 89ghz channel instead of amsua 89ghz channel
!  for amsua cloud detection
! V1_50        2017-01-09 Robin Faulwetter
!  Restructured/unified cloud detection for radiances;
!  Added features to rttov12; Improvements in rttov handling.
! V1_51        2017-02-24 Harald Anlauf
!  Fix non-conforming format
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Andreas Rhodin   DWD  2010  original source
! Detlef Pingel    DWD  2010
! Robin Faulwetter DWD
!=====================================================================

  !-------------
  ! Modules used
  !-------------
  use mo_kind,          only: wp                   ! working precision kind
  use mo_exception,     only: finish               ! abort routine
  use mo_mpi_dace,      only: dace,               &! MPI group info
                              p_bcast,            &! broadcast routine
                              p_sum
  use mo_dace_string,   only: intstr2array,       &
                              tolower
  use mo_namelist,      only: position_nml,       &! position namelist
                              nnml,               &! namelist Fortran unit
                              POSITIONED           ! ok code from position_nml
  use mo_obs_set,       only: t_obs_set            ! obs. data derived type
  use mo_t_obs,         only: t_spot,             &! report/fov derived type
                              t_ilev,             &! channel/instrument pair
                              operator(==),       &! comparison of t_ilev
                              debug_spot,         &! debug selected spot
                              spt_hd_debug, ldeb, &
                              usd, dpref
  use mo_fdbk_tables,   only: OT_RAD,             &! Radiances report type ID
                              SUR_MWSURF,         &! surftype: bad surf by MW sounder
                              SUR_LAND,           &! surftype: land
                              SUR_SEA,            &! surftype: sea
                              TF_SURF_RETR,       &! tovs_flag: wrong surface type retrieved
                              TF_SURF_INFL,       &! tovs_flag: surface influence detected
                              TF_CLOUD,           &! tovs_flag: cloud
                              TF_AEROSOL,         &! tovs_flag: aerosol
                              TF_DESERT_DUST,     &! tovs_flag: desert dust
                              TF_VOLCANIC_ASH,    &! tovs_flag: volcanic ash
                              TF_TRACE_GAS         ! tovs_flag: trace gas
  use mo_dec_matrix,    only: t_vector             ! vector data type
  use mo_t_use,         only: CHK_CLOUD,          &! reject due to cloud flag
                              CHK_SURF,           &! reject due to surface type
                              CHK_INSDAT,         &! reject due to insufficient data
                              STAT_PAS_REJ,       &! report is monitored + rejected
                              STAT_REJECTED
  use mo_tovs,          only: decr_tovs_use,      &
                              clchk_amsub,        &
                              param_file_amsub,   &
                              amsub_ec_bnd,       &
                              amsub_delta_ch20_ch18,&
                              p_bcast
  use mo_t_tovs,        only: t_tovs,             &! radiance data derived type
                              t_tovs_instr,       &! info on instruments in t_tovs
                              load,               &! load  t_tovs
                              store,              &! store t_tovs
                              TTOVS_BASE,         &! flag for loading/storing t_tovs
                              TTOVS_IM,           &! flag for loading/storing t_tovs
                              TTOVS_IMCH,         &! flag for loading/storing t_tovs
                              TTOVS_CLDLEV,       &! flag for loading/storing t_tovs
                              destruct             ! deallocate t_tovs components
  use mo_rad,           only: USE_CLOUDDET,       &
                              USE_CLEAR,          &
                              USE_MWSURF,         &
                              USE_MWICE,          &
                              USE_QCNOTUSE,       &! do not use, e.g. for cloud det.
                              USE_LAND,           &
                              t_rad_set,          &
                              rad_set,            &
                              n_set,              &
                              back_nml,           &
                              construct,          &
                              destruct,           &
                              assignment(=),      &
                              m_instr,            &
                              m_bd,               &
                              m_oe,               &
                              lev2p,              &
                              amsua_chan_id,      &
                              amsub_chan_id,      &
                              OPTV_CLD_FRC         ! imager cloud fraction availability
  use mo_rttov,         only: preslev
  use mo_range_fparse,  only: t_range_fparse,     &
                              construct,          &
                              destruct,           &
                              p_bcast,            &
                              init
  use mo_tovs_prof,     only: c_var
  use mo_cloud_ir,      only: hirs_cloud_check, &
                              read_nml_hirs
  use mo_cloud_params,  only: param_file_scatt,   &
                              param_file_si_amsua,&
                              param_amsub,        &
                              read_thresh_params, &
                              read_param_amsub,   &
                              get_threshold,      &
                              get_threshold_buehler,&
                              inv_thr,            &
                              p_bcast,            &
                              dlat_thresh_def,    &
                              dsz_thresh_def,     &
                              buehler_version
  use mo_cloud_indices, only: mw_emiss,           &
                              si_amsua,           &!
                              lwp_KaSiRu94,       &!
                              LWP_QinZou2016,     &!
                              TPW_QinZou2016,     &!
                              L_index_QinZou2016, &
                              inv_ind
  use mo_cads_ifc,      only: toplev,             &
                              bottomlev,          &
                              cads_init_ifc,      &
                              cads_ifc,           &
                              CLOUDFREE,          &
                              CLOUDY,             &
                              MNW_CLOUD,          &
                              MNW_AEROSOL,        &
                              MNW_LAND,           &
                              MNW_TRACEGAS,       &
                              MNW_SURF,           &
                              l_aerosol,          &
                              l_land_sens,        &
                              l_trace_gas,        &
                              use_im_frac,        &
                              IMFRC_CADS

  implicit none

  !----------------
  ! Public entities
  !----------------
  private
  public :: cloud_detect            ! apply cloud detection routines
  public :: cloud_detect_init       ! setup cloud detection
  public :: calc_tovs_cloud         ! Calc. cloud for use in RTTOV
  public :: PH_ABC                  ! cloud detection after bias correction
  public :: PH_BFG                  ! cloud detection before first guess checks
  public :: PH_ATH                  ! cloud detection after thinning


  ! Information on phase of processing
  integer, parameter :: n_phases = 4
  integer, parameter :: PH_CC    = 0 ! cloud detection in calc_tovs_cloud
  integer, parameter :: PH_ABC   = 1 ! cloud detection after bias correction
  integer, parameter :: PH_BFG   = 2
  integer, parameter :: PH_ATH   = 3
  integer, parameter :: phases(0:n_phases-1) = (/PH_CC, PH_ABC, PH_BFG, PH_ATH /)

  !---------------------
  ! check parameters
  !---------------------
  type t_chk_par
    character(len=3000)        :: thresh_str = ''
    type(t_range_fparse)       :: threshold
    integer                    :: ch(3)      = -1
    integer                    :: instr(3)   = -1
    logical                    :: use_rawbt  = .TRUE.
  end type t_chk_par
  type(t_range_fparse),  save  :: empty_trf

  ! Information on a cloud check
  integer,           parameter :: CLD_CHK_SEA   = 0
  integer,           parameter :: CLD_CHK_LAND  = 1
  integer,           parameter :: CLD_CHK_SEAICE= 2
  integer,           parameter :: notappl_prio  = 3
  integer,           parameter :: notappl_instr = 4
  character(len=20), parameter :: c_notappl(CLD_CHK_SEA:notappl_instr) = (/ &
       'sea surface         ', 'land surface        ', &
       'no sea ice surface  ', 'priority            ', &
       'instr. not avail.   ' /)

  integer,           parameter  :: len_fnct      = 80
  integer,           parameter  :: nsat          = 20
  integer,           parameter  :: ninstr        = 3
  integer,           parameter  :: ngrid         = 5

  type t_cld_check
    integer                     :: instr(ninstr) = -1    ! instruments
    integer                     :: satids(nsat)  = -1    ! restriction to satids
    integer                     :: grids(ngrid)  = -1    ! restriction to grids
    integer                     :: phase         = -99   ! processing phase
    integer                     :: priority      = -99   ! priority
    character(len=80)           :: description   = ''    ! description
    character(len=80)           :: chain         = ''    ! identifies a "chain" of checks
    character(len=len_fnct)     :: fnct          = ''    ! function name
    integer                     :: flags         = ibset(0,CLD_CHK_LAND) + ibset(0,CLD_CHK_SEA)
                                                         ! flags that control the test
    character(len=1000)         :: flagchan      = 'instr'
    integer                     :: tf_cloud_ind(ninstr) = -1
    ! parameters for special functions
    type(t_chk_par)             :: par                   ! parameters for cloud check
  end type t_cld_check
  type(t_cld_check),  save      :: empty_chk

  ! Function parameters
  character(len=len_fnct), parameter :: FNCT_EMISS          = 'emiss_cloud_check'
  character(len=len_fnct), parameter :: FNCT_MNW            = 'mnw'
  character(len=len_fnct), parameter :: FNCT_ALWAYS         = 'always_cloudy'
  character(len=len_fnct), parameter :: FNCT_HIRS           = 'hirs_cloud_check'
  character(len=len_fnct), parameter :: FNCT_POL_DIFF       = 'polarization_diff'
  character(len=len_fnct), parameter :: FNCT_LWP            = 'lwp_approx'
  character(len=len_fnct), parameter :: FNCT_SCATIND        = 'scattering_index'
  character(len=len_fnct), parameter :: FNCT_FGDEP          = 'fgdep'
  character(len=len_fnct), parameter :: FNCT_SI_AMSUA       = 'si_amsua'
  character(len=len_fnct), parameter :: FNCT_AMSUB_LWP      = 'lwp_qinzou2016'
  character(len=len_fnct), parameter :: FNCT_AMSUB_TPW      = 'tpw_qinzou2016'
  character(len=len_fnct), parameter :: FNCT_AMSUB_LINDEX   = 'l_index_qinzou2016'
  character(len=len_fnct), parameter :: FNCT_AMSUB_LWP_S    = 'lwp_qinzou2016_s'
  character(len=len_fnct), parameter :: FNCT_AMSUB_TPW_S    = 'tpw_qinzou2016_s'
  character(len=len_fnct), parameter :: FNCT_AMSUB_BUEHLER  = 'amsub_buehler'
  character(len=len_fnct), parameter :: FNCT_AMSUB_3_THRESH = 'amsub_3_thresh'
  character(len=len_fnct), parameter :: FNCT_AMSUB_EC1      = 'amsub_ec1'
  character(len=len_fnct), parameter :: FNCT_AMSUB_5_3      = 'amsub_5_3'
  character(len=len_fnct), parameter :: FNCT_AMSUB_EC2      = 'amsub_ec2'

  ! Error handling
  integer, parameter :: ERR_INV_SURF       =  1           ! config file for amsub/mhs c.d. not read
  integer, parameter :: ERR_INV_THR        =  2           ! config file for amsub/mhs c.d. not read
  integer, parameter :: ERR_ALL_MISSING    =  3           ! All channels missing
  integer, parameter :: ERR_PHYS           =  4           ! config file for amsub/mhs c.d. not rea
  integer, parameter :: ERR_LAT            =  5           ! Invalid latitude
  integer, parameter :: ERR_CADS           =  6           ! cads_ifc failed
  integer, parameter :: ERR_MW_EMISS       =  7           ! cads_ifc failed
  integer, parameter :: ERR_CF_AMSUB       =  8           ! config file for amsub/mhs c.d. not read
  integer, parameter :: ERR_NO_BUEHLER_PAR =  9           ! did not find parameters for mod. Buehler
  integer, parameter :: ERR_NO_HIRS_PARAM  = 10           ! HIRS cloud detection parameters missing

  type t_err
    integer                     :: n     = 0
    integer                     :: code  = 0
    integer                     :: ichk  = -1
    character(len=80)           :: descr = ''
  end type t_err

  integer,          parameter   :: im_fl_nb = 4  ! 0..3, 0..2 by CADS*_Detect_Cloud_Imager, 3 by dace (see IMFRC_CADS bit)
  integer,          parameter   :: im_fl_u  = 2**im_fl_nb-1
  character(len=7), parameter   :: im_fl_name(im_fl_nb) = (/'fg     ', 'intercl', 'homog  ', 'cldfrac'/)

  type t_cld_tmp
    ! type to store temporary stuff on cloud detection
    integer                     :: n_fovs            =  0         ! #FOVs
    integer,        allocatable :: chk       (:)                  ! #cloud checks on grid    (dim:n_chk)
    logical,        allocatable :: chk_done  (:)                  ! channel mask for flagging
    integer                     :: n_chk                          ! #checks on grid
    integer,        allocatable :: n_data    (:)                  ! #Tests                   (dim:n_chk)
    integer,        allocatable :: n_notappl (:,:)                ! #Not applied tests       (dim:n_chk,nflags)
    integer,        allocatable :: n_applied (:)                  ! #Tests applied           (dim:n_chk)
    integer,        allocatable :: n_failed  (:)                  ! #Tests that failed       (dim:n_chk)
    integer,        allocatable :: n_cld_chan(:,:)                ! #cloudy data             (dim:n_chan,n_chk)
    integer,        allocatable :: n_srf_chan(:,:)                ! #data with "bad" surface (dim:n_chan,n_chk)
    integer,        allocatable :: n_aer_chan(:,:)                ! #data with aerosol       (dim:n_chan,n_chk)
    integer,        allocatable :: n_trg_chan(:,:)                ! #data with aerosol       (dim:n_chan,n_chk)
    integer,        allocatable :: n_use_mhs1(:)                  ! #data MHS1 instead of AMSUA15       (n_chk)
    type(t_ilev),   allocatable :: chi       (:)                  ! channel identification
    type(t_err),    pointer     :: err       (:)     => null()    ! error message info
    integer                     :: nerr              =  0         ! number of errors that occurred
    integer                     :: n_im_flag(0:im_fl_u) = 0       ! number of FOVs, colocated imager flag
    logical,        allocatable :: flagmask  (:,:)                ! channel mask for flagging
  end type t_cld_tmp


  integer,           parameter      :: n_surf    =  12            ! surface types from AMSU-A cloud check

  type(t_cld_check), save, pointer  :: check(:) => null()         ! global info in all cloud checks
  type(t_cld_check), save, pointer  :: chk      => null()         ! auxiliary pointer to cloud check info
  integer,           save           :: n_chk    = 0               ! number of active cloud checks

  !-------------------
  ! Namelist variables
  !-------------------
  integer,           save           :: verbose                    ! verbosity level (0: silent, >0: talkative)
  logical,           save           :: do_clddetect_irs           ! do IRS    cloud detection in 3dvar
  logical,           save           :: do_clddetect_amsua         ! do AMSU-A cloud detection in 3dvar
  logical,           save           :: do_clddetect_iasi          ! do IASI   cloud detection in 3dvar
  logical,           save           :: do_clddetect_hirs          ! do HIRS   cloud detection in 3dvar
  logical,           save           :: do_clddetect_atms          ! do ATMS   cloud detection in 3dvar
  logical,           save           :: do_clddetect_atms_amsub    ! do AMSU-B like cloud detection for ATMS humidity chan.
  logical,           save           :: do_clddetect_cris          ! do CrIS   cloud detection in 3dvar
  logical,           save           :: do_clddetect_mwts          ! do MWTS   cloud detection in 3dvar
  logical,           save           :: do_clddetect_gmi           ! do GMI    cloud detection in 3dvar
  logical,           save           :: do_clddetect_ssmis         ! do SSMI/S cloud detection in 3dvar
  logical,           save           :: do_clddetect_saphir        ! do SAPHIR cloud detection in 3dvar
  logical,           save           :: do_clddetect_amsr2         ! do AMSR2  cloud detection in 3dvar
  logical,           save           :: do_clddetect_mwhs2         ! do MWHS2  cloud detection in 3dvar
  logical,           save           :: do_clddetect_mhs           ! do MHS    cloud detection in 3dvar
  logical,           save           :: do_clddetect_mwts3         ! do MWTS-3 cloud detection in 3dvar
  logical,           save, target   :: use_amsua_surftype(n_surf) ! use surface types 1-12 from amsua cloud check
  logical,           save, target   :: use_atms_surftype (n_surf) ! use surface types 1-12 from atms  cloud check
  logical,           save, target   :: use_mwts3_surftype(n_surf) ! use surface types 1-12 from mwts3 cloud check
  logical,           save, target   :: use_amsua_ice     (n_surf) ! surface types considered as ice for MWICE flag
  logical,           save, target   :: use_atms_ice      (n_surf) ! surface types considered as ice for MWICE flag
  logical,           save, target   :: use_mwts3_ice     (n_surf) ! surface types considered as ice for MWICE flag
!  real(kind=wp),     save           :: topLev      =  0._wp       ! cloud search top level (hPa)
!  real(kind=wp),     save           :: bottomLev   =  1.E20_wp    ! cloud search bottom level (hPa)

  namelist /TOVS_CLOUD/ verbose, do_clddetect_amsua, do_clddetect_atms,        &
       do_clddetect_atms_amsub, do_clddetect_iasi, do_clddetect_hirs,          &
       do_clddetect_cris, do_clddetect_mwts, do_clddetect_gmi,                 &
       do_clddetect_ssmis, do_clddetect_saphir, do_clddetect_amsr2,            &
       do_clddetect_mwhs2, do_clddetect_mhs, do_clddetect_mwts3,               &
       toplev, bottomlev, use_amsua_surftype, use_atms_surftype,               &
       use_mwts3_surftype,use_amsua_ice,use_atms_ice, use_mwts3_ice,           &
       clchk_amsub, param_file_amsub, amsub_ec_bnd, amsub_delta_ch20_ch18,     &
       dlat_thresh_def, dsz_thresh_def, param_file_scatt, param_file_si_amsua, &
       do_clddetect_irs, buehler_version

  interface p_bcast
    module procedure bcast_cld_check
  end interface

  interface construct
    module procedure construct_chk_par
    module procedure construct_cld_check
  end interface construct


contains

  !-----------------------------------------------------------------
  ! Wrapper for cloud detection routines
  !-----------------------------------------------------------------
  subroutine cloud_detect(obs, y, phase, mask_rs)
    type (t_obs_set) ,intent(inout)        :: obs   ! observation
    type (t_vector)  ,intent(in)           :: y     ! background
    integer          ,intent(in)           :: phase ! phase of processing
    logical          ,intent(in), optional :: mask_rs(:) !

    character(len=*), parameter         :: proc = 'cloud_detect'
    integer                             :: ib                   ! observation 'box'      index
    integer                             :: is                   ! observation report/fov index
    type(t_spot),    pointer            :: s                    ! pointer to report/fov data
    integer                             :: ierr                 ! Error flag
    integer                             :: ierr_max             ! Maximum error within a priority level
    type(t_cld_check),pointer           :: chk2 => null()       ! auxiliary pointer to cloud check
    integer                             :: ichk                 ! check loop index
    integer                             :: nchk                 ! total number of checks
    integer                             :: o_ch,n_ch            ! channel indices
    integer                             :: nc                   ! number of useable channels (mask==.true.)
    integer                             :: i,j,k,ii             ! indices
    integer                             :: i1, in               ! indices
    integer                             :: iset                 ! rad_set stuff
    integer                             :: tovs_flag
    logical,         target             :: l_check_amsua_surf   ! Check estimated surface type for AMSUA
    logical,         target             :: l_check_atms_surf    ! Check estimated surface type for AMSUA
    logical,         target             :: l_check_mwts3_surf   ! Check estimated surface type for MWTS3
    logical                             :: l_apply              ! apply test
    logical,         allocatable        :: mask (:)
    logical,         allocatable        :: flagmask(:)          ! mask for flagging channels
    logical                             :: clear
    character(50)                       :: decr_hint  = ''
    ! t_tovs and corresponding rad_set
    type(t_tovs),    target             :: tt                   ! radiance data type
    type(t_tovs_instr), target          :: ti                   ! info on instruments in t_tovs
    type(t_rad_set), pointer            :: rs => null()
    integer,         pointer            :: n_instr => null()
    integer,         pointer            :: n_chan  => null()
    integer,         pointer            :: o_ch_i(:) => null()
    integer,         pointer            :: n_ch_i(:) => null()
    integer,         pointer            :: n_lev
    integer                             :: instr(m_instr)
    integer,         allocatable        :: chan(:)
    integer,         allocatable        :: flag(:)
    integer,         allocatable        :: band(:)
    logical                             :: ld

    ! statistics of results
    type(t_cld_tmp), allocatable,target :: cld_tmp(:)
    type(t_cld_tmp), pointer            :: ctmp => null()

    if (.not.any(phases(:) == phase)) then
      call finish(proc,'called with invalid phase')
    end if

    if (n_set <= 0) return

    l_check_amsua_surf = .not.all(use_amsua_surftype(1:n_surf)) .or. &
                         .not.all(use_amsua_ice     (1:n_surf))
    l_check_atms_surf  = .not.all(use_atms_surftype (1:n_surf)) .or. &
                         .not.all(use_atms_ice      (1:n_surf))
    l_check_mwts3_surf = .not.all(use_mwts3_surftype(1:n_surf)) .or. &
                         .not.all(use_mwts3_ice     (1:n_surf))

    call construct_cld_tmp
    nchk = sum(cld_tmp(1:n_set)%n_chk)

    if (nchk > 0) then

      !---------------------------------
      ! apply cloud detection algorithms
      !---------------------------------
      do ib = 1,size(obs% o)                  ! loop over observation 'boxes'
        if (obs% o(ib)%pe /= dace%pe) cycle
        do is=1,obs% o(ib)% n_spot            ! loop over reports/fovs
          s => obs% o(ib)%spot(is)
          if (s% hd% obstype /= OT_RAD) cycle ! only treat radiances

          ld = ldeb(s)
          ! load t_tovs, and derive some stuff from corresponding rad_set
          ! TODO: think about partial read
          call load(obs% o(ib),s, tt, i_rs=iset, rs=rs, ti=ti)

          if (present(mask_rs)) then
            if (.not.mask_rs(iset)) cycle
          end if

          n_instr => ti%n_instr
          n_chan  => tt%nchan
          n_lev   => tt%nlev
          o_ch_i  => ti%o_ch_i
          n_ch_i  => ti%n_ch_i
          instr(1:n_instr) = rs%instr(ti%ii(1:n_instr))
          allocate(chan(n_chan), flag(n_chan), band(n_chan))
          chan = rs%chan(tt%ci(1:n_chan))
          flag = rs%flag(tt%ci(1:n_chan))
          band = rs%band(tt%ci(1:n_chan))

          if (iset <= 0 .or. iset > n_set) call finish(proc,'failed to find rad_set entry')
          ctmp => cld_tmp(iset)

          ctmp%n_fovs = ctmp%n_fovs + 1
          check_loop_1: do ichk = 1, ctmp%n_chk
            if (ctmp%chk(ichk) > 0 .and. ctmp%chk(ichk) <= n_chk) then
              chk => check(ctmp%chk(ichk))
            else
              call finish(proc,'invalid ctmp%chk entry (1)')
            end if
            do i = 1, n_instr
              if (any(chk%instr(:) == instr(i))) cycle check_loop_1
            end do
            ctmp%n_notappl(ichk,notappl_instr) = ctmp%n_notappl(ichk,notappl_instr) + 1
          end do check_loop_1

          ctmp%chk_done(:) = .false.

          do i = 1, n_instr
            !-------------------------------
            ! get offset and channel numbers
            !-------------------------------
            o_ch = o_ch_i(i)
            n_ch = n_ch_i(i)
            if (o_ch < 0 .or. o_ch + n_ch > s% o% n) &
                 call finish(proc,'channel offset mismatch')
            i1 = s%o%i+o_ch+1
            in = s%o%i+o_ch+n_ch

            if (allocated(mask)) deallocate(mask)
            allocate(mask(n_ch))
            mask = .not.btest(flag(o_ch+1:o_ch+n_ch),USE_QCNOTUSE)
            nc = count(mask)

            ierr_max = 0

            check_loop: do ichk = 1, ctmp%n_chk
              chk => check(ctmp%chk(ichk))

              if (ctmp%chk_done(ichk)               ) CYCLE check_loop
              if (.not.any(chk%instr(:) == instr(i))) CYCLE check_loop

              !!!! TODO: same surface test! -> not s%stlsf and tt%rt_stype.
              !!!! If possible, all should be tt%rt_stype!!!!
              l_apply = .true.
              if (.not.btest(chk%flags,CLD_CHK_LAND)) then ! do not apply above land
                if (btest(s%stlsf, SUR_LAND)) then
                  ctmp%n_notappl(ichk,CLD_CHK_LAND) = ctmp%n_notappl(ichk,CLD_CHK_LAND) + 1
                  l_apply = .false.
                end if
              end if
              if (.not.btest(chk%flags,CLD_CHK_SEA)) then  ! do not apply above sea
                if (btest(s%stlsf, SUR_SEA) .and. .not.btest(s%stlsf, SUR_LAND)) then
                  ctmp%n_notappl(ichk,CLD_CHK_SEA) = ctmp%n_notappl(ichk,CLD_CHK_SEA) + 1
                  l_apply = .false.
                end if
              end if
              if (btest(chk%flags,CLD_CHK_SEAICE)) then     ! apply above sea ice
                if (tt%rt_stype(1) == 2) then
                   l_apply = .true.
                else
                   ctmp%n_notappl(ichk,CLD_CHK_SEAICE) = ctmp%n_notappl(ichk,CLD_CHK_SEAICE) + 1
                end if
              end if

              if (l_apply) then
                ctmp%n_data(ichk) = ctmp%n_data(ichk) + 1
                ctmp%chk_done(ichk) = .true.

                ! generate flagmask for individual spot on the basis of the general flagmask
                allocate(flagmask(n_chan))
                k = 1
                ii = 1
                do j = 1, n_chan
                  do while (j > o_ch_i(ii) + n_ch_i(ii))
                    ii = ii + 1
                  end do
                  do while (ctmp%chi(k)%instr /= instr(ii) .or. ctmp%chi(k)%value /= chan(j))
                    k = k + 1
                    if (k > ubound(ctmp%chi,1)) &
                         call finish('cloud_check','Failed to set up flagmask')
                  end do
                  flagmask(j) = ctmp%flagmask(k,ichk)
                end do

                tovs_flag = -1
                do j = 1, ninstr
                  if (chk%instr(j) == instr(i)) then
                    tovs_flag = chk%tf_cloud_ind(j)
                  end if
                end do
                if (tovs_flag < 1 .or. tovs_flag > size(TF_CLOUD)) then
                  call finish(proc,'invalid chk%tf_cloud_ind')
                end if
                tovs_flag = TF_CLOUD(tovs_flag)

                if (ld) write(usd,*) dpref,trim(chk%fnct),tovs_flag
                ! Apply check
                select case(chk%fnct)
                case(FNCT_EMISS)
                  call emiss(ierr)
                  call store(obs% o(ib),s, tt, tovs_io=TTOVS_BASE)
                case(FNCT_AMSUB_BUEHLER, FNCT_AMSUB_3_THRESH, FNCT_AMSUB_EC1,   FNCT_AMSUB_5_3,&
                     FNCT_AMSUB_EC2,     FNCT_AMSUB_LWP,      FNCT_AMSUB_LWP_S, FNCT_AMSUB_TPW,&
                     FNCT_AMSUB_TPW_S,   FNCT_AMSUB_LINDEX)
                  call amsub(ierr)
                case(FNCT_ALWAYS)
                  call always_cloudy(ierr)
                case(FNCT_POL_DIFF)
                  call polarization_diff(ierr)
                case(FNCT_LWP)
                  call lwp_approx(ierr)
                case(FNCT_SCATIND)
                  call scatindex(ierr)
                case(FNCT_SI_AMSUA)
                  call si_amsua_wrap(ierr)
                case(FNCT_FGDEP)
                  call fgdep(ierr)
                case(FNCT_MNW)
                  call mnw_wrap(ierr)
                case(FNCT_HIRS)
                  call hirs_wrap(ierr)
                case default
                  call finish(proc,'function "'//trim(chk%fnct)//'" is not implemented.')
                end select

                deallocate(flagmask)
              elseif (.not.btest(s%stlsf, SUR_LAND).and.&
                      .not.btest(s%stlsf, SUR_SEA)) then
                ! Use fallback for invalid surface flags
                ierr = ERR_INV_SURF
                call add_err(ierr, ichk, descr='invalid surface type')
              else
                ierr = 0
              end if

              ierr_max = max(ierr_max, abs(ierr))

              call debug_spot(sp=s, obs=obs%o(ib), flags=flag, hint=proc//'/'//trim(chk%fnct))

              ! Check for lower priority checks
              lower_prio: do j = ichk+1, ctmp%n_chk
                chk2 => check(ctmp%chk(j))
                if (.not.any(chk2%instr(:) == instr(i))) EXIT  lower_prio
                if (ctmp%chk_done(j)                   ) CYCLE lower_prio
                if (chk2%chain /= chk%chain            ) CYCLE lower_prio
                if (chk2%priority == chk%priority) then
                  ! same priority level
                  EXIT lower_prio
                else
                  ! next priority level
                  if (ierr_max == 0) then
                    ! Deactivate following checks with lower priority
                    ctmp%n_notappl(j,notappl_prio) = ctmp%n_notappl(j,notappl_prio) + 1
                    ctmp%chk_done (j) = .true.
                  else
                    ! Next priority level
                    ierr_max = 0
                    EXIT lower_prio
                  end if
                end if
              end do lower_prio

            end do check_loop

          end do
          call destruct (tt)
          deallocate(chan, flag, band)
        end do
      end do
    end if

    call print_stat

    call destruct_cld_tmp

  contains

    ! Microwave emissivity check
    ! TODO: SSMIS
    subroutine emiss(ierr)
      integer, intent(out)  :: ierr

      real(wp)              :: t_b_bc(4) ! obs bright.temp., bias-corr.
      real(wp)              :: tskin     ! surf. skin temp.
      real(wp)              :: landfrac  ! land fraction
      integer               :: surftype  ! estimated surface type
      integer               :: chans(4)
      integer               :: j, k
      integer               :: i_body(4)   ! indices of AMSUA 1-3,15
      integer               :: i_amsua
      logical               :: mask_(n_chan), l_amsub
      logical,  pointer     :: l_surf       => null()
      logical,  pointer     :: use_surf (:) => null()
      logical,  pointer     :: use_ice(:) => null()

      ierr = 0
      if (instr(i) /= 3 .and. instr(i) /= 19 .and. instr(i) /= 132) &
           call finish(proc//'/emiss','invalid instrument')

      i_body(:) = -1
      chans (:) = -1
      do j=1,n_instr
         do k = o_ch_i(j)+1,o_ch_i(j)+n_ch_i(j)
            if (.not.btest(flag(k),USE_QCNOTUSE)) then
               i_amsua = amsua_chan_id(instr(j), chan(k))
               if (i_amsua >= 1 .and. i_amsua <=3) then
                  if (i_body(i_amsua) < 0 .or. instr(j) == instr(i)) then
                    i_body(i_amsua) = s%o%i + k
                    chans (i_amsua) = chan(k)
                  end if
               elseif (i_amsua == 15) then
                  l_amsub = any(instr(j) == (/4,15,10,34,71,73/))
                  if (i_body(4) < 0 .or. instr(j) == instr(i) .or. &
                      (rs%gopts%use_amsub_89ghz .and. l_amsub)) then
                     i_body(4) = s%o%i + k
                     chans (4) = chan(k)
                     if (l_amsub) then
                       ctmp%n_use_mhs1(ichk) = ctmp%n_use_mhs1(ichk) + 1
                       select case(instr(i))
                       case(3)   ; chans(4) = 15
                       case(19)  ; chans(4) = 16
                       case(132) ; chans(4) = 18
                       end select
                     end if
                  end if
               end if
            end if
         end do
      end do

      do j = 1, 4
        if (i_body(j) <= 0) then
          select case(j)
          case(1:3)
            ierr = -j
          case(4)
            ierr = -15
          end select
        end if
      end do

      if (ierr == 0) then
        t_b_bc(:) = obs%o(ib)%body (i_body(:))%o

        if (ld) write(usd,*) dpref,trim(FNCT_EMISS),'t_b_bc',t_b_bc(:),chans(:)

        if (tt%ts > 0) then
          tskin = tt%ts
        else
          tskin = s% ts_bg
        end if
        ! TODO: use t_tovs%rt_stype
        landfrac = s% sl_bg

        call mw_emiss(instr(i), &
                      chans,    &
                      t_b_bc,   &
                      s%stzen,  &
                      landfrac, &
                      tskin,    &
                      ierr,     &
                      surftype, &
                      ld=ld)
        clear = all(surftype /= (/10,11/) )   ! cloud/rain over water
        if (ierr > 0) ierr = ERR_MW_EMISS

      else
        clear = .false.
      end if

      if (ierr == 0) then
        ctmp%n_applied(ichk) = ctmp%n_applied(ichk) + 1
      else
        call add_err(ierr,ichk, n=j)
        return
      end if

      if (.not. clear) then
        call decr_tovs_use(s, obs% o(ib), CHK_CLOUD, tovs=tt, bit=USE_CLEAR,&
             chan_mask=flagmask(:), hint_debug='cld_detect '//trim(FNCT_EMISS),&
             tovs_flag=tovs_flag)
        call add_stat('cld', 1, n_chan, bit=USE_CLEAR, chan_mask=flagmask(:))
      end if

      !------------------------------
      ! check for surface type to use
      !------------------------------
      if (instr(i) == 3) then
        l_surf   => l_check_amsua_surf
        use_surf => use_amsua_surftype
        use_ice  => use_amsua_ice
      elseif (instr(i) == 19) then
        l_surf   => l_check_atms_surf
        use_surf => use_atms_surftype
        use_ice  => use_atms_ice
      elseif (instr(i) == 132) then
        l_surf   => l_check_mwts3_surf
        use_surf => use_mwts3_surftype
        use_ice  => use_mwts3_ice
      end if

      if (l_surf) then
        if ((surftype > 0).and.(surftype <= n_surf)) then
          if (.not.use_surf(surftype)) then
            mask_(:) = .true.
            ! Do not flag obs. that are very likely over ice if USE_MWICE is set
            if (use_ice(surftype) .and. (s%fi_bg > 0._wp .or. s%ts_bg <= 272._wp)) then
              where(btest(rs%flag(tt%ci(1:n_chan)), USE_MWICE)) &
                   mask_(1:n_chan) = .false.
            end if
            if (ld) write(usd,*) dpref,trim(FNCT_EMISS),' mw_ice',&
                 surftype,use_ice(surftype),s%fi_bg,s%ts_bg,mask_(1:n_chan)
            call decr_tovs_use(s, obs% o(ib), CHK_SURF, tovs=tt, tovs_flag=TF_SURF_RETR, &
                 &not_bit=USE_MWSURF, chan_mask=mask_,                                   &
                 &hint_debug='cld_detect '//trim(FNCT_EMISS)//' surf')
            s% stlsf = ibset(s% stlsf, SUR_MWSURF)
            call add_stat('srf',1,n_chan,not_bit=USE_MWSURF)
          endif
          tt%mw_stype = surftype
        endif
      endif

    end subroutine emiss

    ! AMSU-B cloud checks
    subroutine amsub(ierr)
      integer, intent(out)    :: ierr

      real(wp)              :: o    (5)  ! observed (bias-corrected) BT
      real(wp)              :: fc   (5)  ! First Guess BT
      real(wp)              :: lat       ! latitude
      real(wp)              :: thr, x
      real(wp)              :: hl_max, bt_h_min
      integer               :: j, k, i_amsub
      integer               :: ibody(5)
      logical               :: l_sea, l_high, l_low_high

      ierr    = 0

      ibody(:) = -1
      do j=1,n_instr
         do k = o_ch_i(j)+1,o_ch_i(j)+n_ch_i(j)
            if (.not.btest(flag(k),USE_QCNOTUSE)) then
               i_amsub = amsub_chan_id(instr(j), chan(k))
               if (i_amsub > 0) then
                  if (ibody(i_amsub) < 0 .or. instr(j) == instr(i)) then
                     ibody(i_amsub) = s%o%i + k
                  end if
               end if
            end if
         end do
      end do

      if (.not. chk%par%use_rawbt) then
        where(ibody(:) > 0) o (:) = obs%o(ib)%body(ibody(:))%o
      else
        where(ibody(:) > 0) o (:) = obs%o(ib)%body(ibody(:))%o - obs%o(ib)%body(ibody(:))%bc
      end if
      where(ibody(:) > 0)   fc(:) = y%s(ib)%x(ibody(:))

      clear = .false.

#define CHKCHAN(ch) if (ibody(ch) <= 0) ierr = -ch

      select case(chk%fnct)
      case(FNCT_AMSUB_LWP,FNCT_AMSUB_LWP_S,FNCT_AMSUB_TPW,FNCT_AMSUB_TPW_S,FNCT_AMSUB_LINDEX)
        ! see:
        ! Qin Zhengkun and Zou Xiaolei, 2016: Development and initial assessment of a new land index
        ! for microwave humidity sounder cloud detection. J. Meteor. Res., 30(1), 012-037,
        ! doi: 10.1007/s13351-016-5076-4
        if (chk%par%threshold%func == '') &
             call finish(proc//'/amsub(QinZou)','missing parameters')
        thr = get_threshold(chk%par%threshold, s)
        if (thr /= inv_thr) then
          if (chk%fnct == FNCT_AMSUB_LWP_S .or. chk%fnct == FNCT_AMSUB_TPW_S) then
            l_sea = .true.
          else
            l_sea = btest(s%stlsf, SUR_SEA).and..not.btest(s%stlsf, SUR_LAND)
          end if
          ! TODO: obsolete? get_tovs_var might be used for index calc.
          select case(chk%fnct)
          case (FNCT_AMSUB_LWP, FNCT_AMSUB_LWP_S)
            CHKCHAN(1)
            CHKCHAN(2)
            if (ierr == 0) x = LWP_QinZou2016(o(1:2), fc(1:2), l_sea)
          case (FNCT_AMSUB_TPW, FNCT_AMSUB_TPW_S)
            CHKCHAN(1)
            CHKCHAN(2)
            if (ierr == 0) x = TPW_QinZou2016(o(1:2), fc(1:2), l_sea)
          case (FNCT_AMSUB_LINDEX)
            do j = 1, 5
              CHKCHAN(j)
            end do
            if (ierr == 0) x = L_index_QinZou2016(o(1:5))
          end select
          if (x == inv_ind) ierr = ERR_PHYS
          if (ierr == 0) then
            clear = (x < thr)
          else
            clear = .false.
          end if
        else
          ierr = ERR_INV_THR
        end if
      case (FNCT_AMSUB_BUEHLER, FNCT_AMSUB_3_THRESH, FNCT_AMSUB_5_3)
        l_high     = any(chk%fnct == (/FNCT_AMSUB_BUEHLER,FNCT_AMSUB_3_THRESH/))
        l_low_high = any(chk%fnct == (/FNCT_AMSUB_BUEHLER,FNCT_AMSUB_5_3     /))
        lat = s%col%c%dlat
        if (lat < -90._wp .or. lat > 90._wp) ierr = ERR_LAT
        if (l_low_high) then
          CHKCHAN(3)
          CHKCHAN(5)
        elseif(l_high) then
          CHKCHAN(3)
        end if
        if (ierr == 0) then
          call get_threshold_buehler(s%stzen, lat, o(3), rs%satid, instr(i), bt_h_min, hl_max, ierr)
          if (ld) write(usd,*) dpref,chk%fnct,'get_threshold_buehler', bt_h_min, hl_max,ierr
          if ( ierr == 0 .and. bt_h_min /= inv_thr .and. hl_max /= inv_thr .and. &
               (buehler_version == 0 .or. bt_h_min > 0._wp)) then
            clear = .true.
            if (l_high     .and. o(3)      <  bt_h_min) clear = .false.
            if (l_low_high .and. o(5)-o(3) <= hl_max  ) clear = .false.
          elseif (ierr == 2) then
            ierr = ERR_CF_AMSUB
          elseif (ierr > 0) then
            ierr = ERR_NO_BUEHLER_PAR
          else
            ierr = ERR_INV_THR
          end if
          if (ld) write(usd,*) dpref,chk%fnct,clear,ierr
        end if
      case (FNCT_AMSUB_EC1)
        CHKCHAN(2)
        if (ierr == 0) then
          clear = o(2) - fc(2) > param_amsub% ec_bnd(1) .and.  &
                  o(2) - fc(2) < param_amsub% ec_bnd(2)
        end if
      case (FNCT_AMSUB_EC2)
        CHKCHAN(1)
        if (ierr == 0) then
          clear = o(1) - fc(1) > param_amsub% ec_bnd(1) .and.  &
                  o(1) - fc(1) < param_amsub% ec_bnd(2)
        end if
      case default
        call finish('amsub@'//proc, 'function "'//trim(chk%fnct)//'" not implemented.')
      end select

      if (ierr == 0) then
        ctmp%n_applied(ichk) = ctmp%n_applied(ichk) + 1
        if (.not.clear) then
          call decr_tovs_use(s, obs%o(ib), CHK_CLOUD, tovs=tt, bit=USE_CLEAR, &
               chan_mask=flagmask(:), hint_debug='cld_check '//trim(chk%fnct),&
               tovs_flag=tovs_flag)
          call add_stat('cld', 1, n_chan, bit=USE_CLEAR, chan_mask=flagmask(:))
        end if
      else
        call add_err(ierr,ichk)
      end if

    end subroutine amsub


    subroutine hirs_wrap(ierr)
      integer, intent(out)    :: ierr

      integer  ,allocatable :: chans (:)
      real(wp) ,allocatable :: o     (:)  ! observed (bias-corrected) BT
      real(wp) ,allocatable :: fc    (:)  ! First Guess BT
      real(wp)              :: ts
      integer  ,allocatable :: cloudy(:)
      integer  ,allocatable :: state (:)

      allocate(chans(nc), o(nc), fc(nc), cloudy(nc), state(nc))
      chans(:) = pack(chan(o_ch+1:o_ch+n_ch),          mask=mask)
      o    (:) = pack(obs%o(ib)%body(i1:in)%o,         mask=mask)
      fc   (:) = pack(  y%s(ib)%x(i1:in),              mask=mask)
      state(:) = pack(obs%o(ib)%body(i1:in)%use%state, mask=mask)

      ts = s%ts_bg  ! This differs with the old version, where run_rttov is used
                    ! var_f% ssv  (1,1) = spot% ts_bg &
                    !                     + sigma_var_tskin * xi% x(spot%i%i+npv+1) ! tsurf
      call hirs_cloud_check(rs%satid, chans, o, fc, ts, state, cloudy, ierr, debug=ld)

      if (ierr == 0) then
        ctmp%n_applied(ichk) = ctmp%n_applied(ichk) + 1
      else
        if (ierr > 0) ierr = ERR_NO_HIRS_PARAM
        call add_err(ierr,ichk)
      end if

      ! Do the flagging always, also if test failed, i.e. no fallback cloud detection method for HIRS
      call decr_tovs_use(s, obs%o(ib), CHK_CLOUD, tovs=tt, bit=USE_CLEAR, instr=instr(i), &
           chan_mask=unpack(cloudy==1, mask, .true.), hint_debug='hirs_cld_check',&
           tovs_flag=tovs_flag)
      call add_stat('cld',o_ch+1,o_ch+n_ch,bit=USE_CLEAR, chan_mask=unpack(cloudy==1, mask, .true.))
      if (any(cloudy == 3)) then
        call decr_tovs_use(s, obs%o(ib), CHK_INSDAT, state=STAT_REJECTED, tovs=tt, bit=USE_CLEAR, instr=instr(i), &
             chan_mask=unpack(cloudy==3 .and. state > STAT_PAS_REJ .and. state/= STAT_REJECTED, mask, .true.),      &
             hint_debug='hirs_cld_check', tovs_flag=tovs_flag)
      end if

    end subroutine hirs_wrap


    ! McNally-Watts cloud detection
    subroutine mnw_wrap(ierr)
      integer, intent(out) :: ierr
      !-----------------------------------------------------------------
      ! IR cloud detection routine to be called with fg and obs,
      ! applies McNally-Watts cloud detection algorithm to observations,
      ! in addition application of the IASI level 2 cloud flags
      ! documentation of the method:
      ! http://research.metoffice.gov.uk/research/interproj/nwpsaf/IR_aerosol_cloud_detect/user_manual/index.html
      !-----------------------------------------------------------------
      integer                :: chans            (n_ch) ! all channels
      integer                :: mnw_flags        (n_ch) ! flags from MNW cld. detection
      real(wp)               :: level_chans      (n_ch) ! channels used for cld. detect.
      real(wp), allocatable  :: cloud_lev        (:)    ! cloud top level for each band
      real(wp)               :: t_b_fg           (n_ch) ! fg  bright.temp.
      real(wp)               :: t_b_bc           (n_ch) ! obs bright.temp., bias-corr.
      real(wp)               :: t_b_obs          (n_ch) ! obs bright.temp.
      real(wp)               :: p                (n_lev)
      real(wp)               :: l2c_max                 ! maximum allowed l2c (from cads_ifc surf_check)
      real(wp)               :: si_max                  ! maximum allowed sinfl
      real(wp)               :: ct(1)                   ! cloud top
      logical                :: fmask            (n_ch)
      character(len=100)     :: msg
      character(len=80)      :: logname
      character(len=8)       :: instr_name
      integer                :: iband
      integer                :: n_ch_used            ! no.of channels for cld.detect.
      integer                :: ioff
      integer                :: aer_type
      integer                :: tf
      integer                :: k
      integer                :: im_flag
      integer                :: m_band
      integer                :: tt_io
      logical                :: l_aer, l_im_cl, l_im_fr

      ierr = 100

      write(instr_name, '("instr",I3.2)') instr(i)
      ! sounder channels used for cloud detection:
      where (band(o_ch+1:o_ch+n_ch)/=0                   .and. &
             btest(flag(o_ch+1:o_ch+n_ch), USE_CLOUDDET) .and. &
             mask(1:n_ch))
        fmask(:) = .true.
      elsewhere
        fmask(:) = .false.
      end where
      n_ch_used = count(fmask)

      if (dace% lpio) then
        if (n_ch_used <= 0) write(6,*) '*** WARNING: no channels &
             &used in '//trim(instr_name)//' clouddetection.'
      end if

      chans     = chan(o_ch+1:o_ch+n_ch)
      t_b_fg (:) = y% s(ib)% x    (s%o%i+o_ch+1:s%o%i+o_ch+n_ch)
      t_b_bc (:) = obs%o(ib)%body (s%o%i+o_ch+1:s%o%i+o_ch+n_ch)%o
      t_b_obs(:) = t_b_bc(:) - obs%o(ib)%body (s%o%i+o_ch+1:s%o%i+o_ch+n_ch)%bc
      if (tt%i_p > 0) then
        p(1:n_lev) = tt%av(1:n_lev,tt%i_p)
      else
        p(1:n_lev) = preslev(1:n_lev)
      end if

      ! Calculate offset in tt%l2c
      ioff = sum(n_ch_i(1:i-1), mask=rs% iopts(1:i-1)% l2c_type > 0)

      if_used_chans: if (n_ch_used > 0) then
        !--------------------------------------------
        ! first guess and bias-corrected observations
        !--------------------------------------------
        level_chans(1:n_ch) = real(tt%l2c(ioff+1:ioff+n_ch), kind=wp)
        if (any(level_chans(1:n_ch) <= 0._wp)) then
          do k = 1, n_ch
            write(0,*) k,level_chans(k)
          end do
          write(msg,'("levels <= 0. in mnw cloud detection for spot ",I10)') s%hd%id
          call finish("mnw_wrap", trim(msg))
        end if
        m_band = maxval(band(o_ch+1:o_ch+n_ch))
        allocate(cloud_lev(m_band))
        ! Check for imager use
        l_im_fr = btest(rs%gopts%opt_vars, OPTV_CLD_FRC) .and. btest(use_im_frac(instr(i)),IMFRC_CADS)
        l_im_cl = (iand(tt%init,TTOVS_IM) == TTOVS_IM .and. min(tt%n_im_cl,tt%n_im_ch) > 0)
        !------------------------------------------------------
        ! call subroutine to perform a cloud detection
        ! for the hyperspectral sounders (AIRS, IASI, AND CRIS)
        !------------------------------------------------------
        call cads_ifc(instr(i),                                 &
                      n_ch_used,                                &
                      pack(chans,                  mask=fmask), &
                      pack(band(o_ch+1:o_ch+n_ch), mask=fmask), &
                      pack(t_b_fg (:),             mask=fmask), &
                      pack(t_b_bc (:),             mask=fmask), &
                      pack(t_b_obs(:),             mask=fmask), &
                      pack(level_chans(:),         mask=fmask), &
                      p(1:n_lev),                               &
                      tt,                                       &
                      i,                                        &
                      btest(rs%gopts%opt_vars, OPTV_CLD_FRC),   &
                      s%ps_bg,                                  &
                      mnw_flags(1:n_ch_used),                   &
                      cloud_lev,                                &
                      s%col%c%dlat,                             &
                      s%col%c%dlon,                             &
                      l2c_max,                                  &
                      si_max,                                   &
                      ierr,                                     &
                      aer_type=aer_type,                        &
                      im_flag=im_flag,                          &
                      lprint=any(spt_hd_debug==s%hd%id)         )

        mnw_flags = unpack(mnw_flags, fmask, (/(CLOUDFREE,k=1,n_ch)/))

        if (ierr == 0) then
          ctmp%n_applied(ichk) = ctmp%n_applied(ichk) + 1
          if (ld) write(usd,*) dpref,'imager flag',im_flag,(im_flag >= 8)
          tt_io = 0
          if (l_im_cl) then
            tt%cld_flg = ior(tt%cld_flg, im_flag)
            tt_io = TTOVS_BASE+TTOVS_IMCH
          elseif (l_im_fr) then
            tt%cld_flg = ior(tt%cld_flg, im_flag)
            tt_io = TTOVS_BASE
          end if
          if (rs%iopts(i)%l_cldlev) then
            tt%cldlev(1:m_band)  = cloud_lev(1:m_band) ! typecast wp->sp
            tt_io = tt_io + TTOVS_CLDLEV
            if (ld) write(usd,*) dpref,'cloud_lev',m_band,tt%cldlev(1:m_band)
          end if
          if (ld) write(usd,*) dpref,'store_tovs',tt_io
          if (tt_io /= 0) then
            call store(obs% o(ib),s, tt, tovs_io=tt_io)
          end if
        else
          ierr = ERR_CADS
        end if

      else
        mnw_flags(:) = CLOUDY
        im_flag      = 0
        ierr         = ERR_ALL_MISSING
      end if if_used_chans

      if (ierr /= 0) then
        call add_err(ierr,ichk)
        return
      end if

      !--------------------------------------------------
      ! apply cloud flags,
      ! print resulting cloud-free channels (diagnostics)
      !--------------------------------------------------
      if (verbose > 0 .and. n_ch_used > 0) then
        write(logname,'(A,I3.3)') trim(instr_name)//'_cloud_detection',dace% pe
        open(10,file=logname,POSITION='APPEND')
        write(10,*) trim(instr_name)//' cloud detection',dace% pe,s%hd%id
        write(10,*) 'spot',s%hd%id,s%id,s%col%c%dlat,s%col%c%dlon
        write(10,*) 'number of channels', n_ch_used
        write(10,*) 'cloud free channels:'
        write(10,*) '      ch/ind    ch/no    fg      obs_bc    fg-obs_bc     ch/height'
        write(10,*) '-------------------------------------------------------------------'
        do k=1, n_ch
          if (band(o_ch+k)    == 0 .or. &
              mnw_flags(k) /= CLOUDFREE) cycle
          write(10,'(i10,i11,4f10.3)') &
               k,chans(k),t_b_fg(k),t_b_bc(k),t_b_fg(k)-t_b_bc(k),tt%l2c(ioff+k)
        enddo
      endif

      ! apply McNally-Watts cloud flag to individual channels:
      if (verbose > 0 ) then
        write(10,*) 'cloudy channels:'
        write(10,*) '      ch/ind    ch/no    fg      obs_bc    fg-obs_bc     ch/height'
        write(10,*) '-------------------------------------------------------------------'
        if (n_ch_used <= 0) write(*,*) 'all channels'
      endif

      ! Flag channels that were not used in the cloud detection
      l_aer = any(btest(mnw_flags(1:n_ch),MNW_AEROSOL))
      if (verbose > 0 ) then
         write(10,*) 'aerosol',l_aer
         write(10,*) 'cloudlev',cloud_lev
         write(10,*) 'channels not used in cloud detection:'
         write(10,*) '      ch/ind    ch/no    fg      obs_bc    fg-obs_bc     ch/height      cloudy'
         write(10,*) '-------------------------------------------------------------------------------'
      endif
      do k=1, n_ch
         iband = band(o_ch+k)
         if (btest(flag(o_ch+k), USE_CLOUDDET)) cycle
         mnw_flags(k) = CLOUDY
         if (iband > 0) then
            if (tt%l2c(ioff+k) < cloud_lev(iband)) mnw_flags(k) = CLOUDFREE
         end if
         if (l_aer) mnw_flags(k) = ibset(mnw_flags(k), MNW_AEROSOL)
         if (l2c_max < huge(l2c_max)) then
           if (tt%l2c(ioff+k) >= l2c_max) mnw_flags(k) = ibset(mnw_flags(k), MNW_SURF)
         end if
         if (si_max < huge(si_max)) then
           if (tt%sinfl(ioff+k) >= si_max) mnw_flags(k) = ibset(mnw_flags(k), MNW_SURF)
         end if

         if (verbose > 0 .and. n_ch_used > 0) write(10,'(i10,i11,4f10.3,2i11)') &
              k,chans(k),t_b_fg(k),t_b_bc(k), t_b_fg(k)-t_b_bc(k),tt%l2c(ioff+k),&
              mnw_flags(k)
      end do

      if (rs%iopts(ti%ii(i))%cloud_mode == -1 .and. any(btest(mnw_flags(1:n_ch),MNW_CLOUD))) then
        ct(1) = minval(cloud_lev)
        call lev2p(p, ct(1:1))
        tt%cloud_top = ct(1)
        if (ld) write(usd,*) dpref,trim(FNCT_MNW),' cloud_top',ct(1),minval(cloud_lev)
        call store(obs% o(ib),s, tt, tovs_io=TTOVS_BASE)
      end if

      ! Flag because of cloud
      fmask = btest(mnw_flags(1:n_ch),MNW_CLOUD).and.flagmask(o_ch+1:o_ch+n_ch)
      call decr_tovs_use(s, obs% o(ib), CHK_CLOUD, tovs=tt, bit=USE_CLEAR, instr=instr(i), &
           chan_mask=fmask, hint_debug='cld_detect '//trim(FNCT_MNW), tovs_flag=tovs_flag)
      call add_stat('cld',o_ch+1,o_ch+n_ch,bit=USE_CLEAR, chan_mask=fmask)
      if (l_aerosol(instr(i))) then
        ! Flag because of aerosols
        tf = ibset(0, TF_AEROSOL)
        select case(aer_type)
        case(1)
          tf = ibset(tf,TF_DESERT_DUST)
        case(2)
          tf = ibset(tf,TF_VOLCANIC_ASH)
        end select
        fmask = btest(mnw_flags(1:n_ch),MNW_AEROSOL).and.flagmask(o_ch+1:o_ch+n_ch)
        call decr_tovs_use(s, obs% o(ib), CHK_CLOUD, tovs=tt, bit=USE_CLEAR, instr=instr(i), &
             chan_mask=fmask, hint_debug='cld_detect '//trim(FNCT_MNW)//' aerosol', &
             tovs_flag_val=tf)
        call add_stat('aer',o_ch+1,o_ch+n_ch,bit=USE_CLEAR, chan_mask=fmask)
      end if
      if (l_land_sens(instr(i))) then
        ! Flag because of land sensitivity
        fmask = btest(mnw_flags(1:n_ch),MNW_LAND).and.flagmask(o_ch+1:o_ch+n_ch)
        call decr_tovs_use(s, obs% o(ib), CHK_SURF, tovs=tt, not_bit=USE_LAND, instr=instr(i), &
             chan_mask=fmask, hint_debug='cld_detect '//trim(FNCT_MNW)//' land_sens', &
             tovs_flag=TF_SURF_INFL)
        call add_stat('srf',o_ch+1,o_ch+n_ch,not_bit=USE_LAND, chan_mask=fmask)
      end if
      if (l_trace_gas(instr(i))) then
        ! Flag because of aerosols
        fmask = btest(mnw_flags(1:n_ch),MNW_TRACEGAS).and.flagmask(o_ch+1:o_ch+n_ch)
        call decr_tovs_use(s, obs% o(ib), CHK_CLOUD, tovs=tt, bit=USE_CLEAR, instr=instr(i), &
             chan_mask=fmask, hint_debug='cld_detect '//trim(FNCT_MNW)//' trace_gas', &
             tovs_flag=TF_TRACE_GAS)
        call add_stat('trg',o_ch+1,o_ch+n_ch,bit=USE_CLEAR, chan_mask=fmask)
      end if
      if (l2c_max < huge(l2c_max) .or. si_max < huge(si_max)) then
        ! Flag because of land sensitivity
        fmask = btest(mnw_flags(1:n_ch),MNW_SURF).and.flagmask(o_ch+1:o_ch+n_ch)
        call decr_tovs_use(s, obs% o(ib), CHK_SURF, tovs=tt, instr=instr(i), &
             chan_mask=fmask, hint_debug='cld_detect '//trim(FNCT_MNW)//' surf_check', &
             tovs_flag=TF_CLOUD(size(TF_CLOUD)))
        call add_stat('srf',o_ch+1,o_ch+n_ch, chan_mask=fmask)
      end if

      if (l_im_cl .or. l_im_fr) then
        if (im_flag >= 0 .and. im_flag <= im_fl_u) ctmp%n_im_flag(im_flag) = ctmp%n_im_flag(im_flag) + 1
      end if

      ! Keep information on cloudiness of channels, that shall be used also in the cloudy case.
      ! These channels are not flagged as cloudy by the above command. This information
      ! might be evaluated by "calc_tovs_quality".
      if (rs%gopts%quality_mode == 6) then
        if (any(.not.btest(flag(o_ch+1:o_ch+n_ch), USE_CLEAR) .and. &
             mnw_flags(1:n_ch) /= CLOUDFREE)) then
          s%pcc = 0
        else
          s%pcc = 100
        end if
        if (verbose > 0) write(10,*) 'pcc',s%pcc
      end if

      if (verbose > 0 .and. n_ch_used > 0) close(10)

    end subroutine mnw_wrap

    subroutine get_chan(ichan, found, o, fc, other_instrs, use_rawbt, ibody)
      integer,  intent(in)            :: ichan
      logical,  intent(out)           :: found
      real(wp), intent(out), optional :: o
      real(wp), intent(out), optional :: fc
      integer,  intent(in),  optional :: other_instrs(:)
      logical,  intent(in),  optional :: use_rawbt
      integer,  intent(out), optional :: ibody

      logical              :: l_rawbt = .false.
      integer              :: k, j
      integer              :: i1_other, n_instrs_other, n_ch_other, o_ch_other
!     logical, allocatable :: mask_other(:)

      if (present(other_instrs)) then
         n_instrs_other = count(other_instrs >= 0)
      else
         n_instrs_other = 0
      end if
      if (present(use_rawbt)) then
         l_rawbt = use_rawbt
      else
         l_rawbt = .false.
      end if

      found = .false.
      if (present(ibody)) ibody = -1
      if (n_instrs_other == 0) then
        do k = 1,n_ch
          if( mask(k) ) then
            if(chan(o_ch+k) == ichan) then
              if (present(o)) then
                if (l_rawbt) then
                  o = obs%o(ib)%body(i1+k-1)%o - obs%o(ib)%body(i1+k-1)%bc
                else
                  o = obs%o(ib)%body(i1+k-1)%o
                end if
              end if
              if (present(fc)) fc = y%s(ib)%x(i1+k-1)
              if (present(ibody)) ibody = i1+k-1
              found = .true.
              exit
            end if
          end if
        end do
      else
        do j=1,n_instr
          if( j /= i .and. any(instr(j) == other_instrs(:)) ) then
            o_ch_other = o_ch_i(j)
            n_ch_other = n_ch_i(j)
            i1_other   = s%o%i + o_ch_other + 1
            do k = 1,n_ch_other
              if(.not.btest(flag(o_ch_other+k),USE_QCNOTUSE)) then
                if(chan(o_ch_other+k) == ichan) then
                  if (present(o)) then
                    if (l_rawbt) then
                      o = obs%o(ib)%body(i1_other+k-1)%o - obs%o(ib)%body(i1_other+k-1)%bc
                    else
                      o = obs%o(ib)%body(i1_other+k-1)%o
                    end if
                  end if
                  if (present(fc)) fc = y%s(ib)%x(i1_other+k-1)
                  if (present(ibody)) ibody = i1_other+k-1
                  found = .true.
                  exit
                end if
              end if
            end do
          end if
        end do
      end if

    end subroutine get_chan


    ! Polarization difference check
    subroutine polarization_diff(ierr)
      integer, intent(out)  :: ierr

      real(wp)                        :: o  (2)      ! observed (bias-corrected) BT ( 37 GHz V, 37 GHz H )
      real(wp)                        :: fc (2)      ! first guess BT               ( 37 GHz V, 37 GHz H )
      real(wp)                        :: quot, thr
      integer                         :: j
      logical                         :: found(2)
      real(wp), parameter             :: small = 0.01_wp    ! below obs.err./fg.err.

      ! Initialization
      if ((chk%par%ch(1) <= 0) .or. (chk%par%ch(2) <= 0) .or. (chk%par%threshold%func == '')) &
           call finish(proc//'/polarization_diff','missing parameters')
      clear     = .true.
      ierr      = 0
      decr_hint = 'cld_check '//trim(chk%fnct)

      do j = 1, 2
        call get_chan(chk%par%ch(j), found(j), o=o(j), fc=fc(j), &
                      other_instrs=(/chk%par%instr(j)/), use_rawbt=chk%par%use_rawbt)
        if (.not.found(j) .and. chk%par%instr(j) < 0) &
             call get_chan(chk%par%ch(j), found(j), o=o(j), fc=fc(j), &
                           other_instrs=(/instr(i)/), use_rawbt=chk%par%use_rawbt)
      end do

      ! Check for cloudy spot
      if(all(found)) then
        clear = .false.
        if ( abs(fc(1) - fc(2)) > small ) then
          quot = ( o(1) - o(2) ) / ( fc(1) - fc(2) )
          thr = get_threshold(chk%par%threshold, s)
          if (thr /= inv_thr) then
            clear = (quot > thr)
            ctmp%n_applied(ichk) = ctmp%n_applied(ichk) + 1
            ! chan_mask: Kanaele die nicht geflaggt werden sollen auf FALSE setzen
          else
            ierr = ERR_INV_THR
          end if
        else
          ierr = ERR_PHYS
        end if
        if( .not.clear ) then
          call decr_tovs_use(s, obs%o(ib), CHK_CLOUD, tovs=tt, bit=USE_CLEAR, &
               chan_mask=flagmask(:), hint_debug=decr_hint, tovs_flag=tovs_flag)
          call add_stat('cld', 1, n_chan, bit=USE_CLEAR, chan_mask=flagmask(:))
        end if
        if (ld) write(usd,*) dpref,trim(chk%fnct),' o=',o,' fc=',fc,clear,ierr
      else
        do j = 1, 2
          if (.not.found(j)) ierr = -chk%par%ch(j)
        end do
      end if
      if (ierr /= 0) call add_err(ierr,ichk)
      !
    end subroutine polarization_diff


    ! Liquid water path approximation check
    subroutine lwp_approx(ierr)
      integer, intent(out)  :: ierr

      real(wp)       :: o(2)  ! observed (bias-corrected) BT ( 22 GHz V, 37 GHz V )
      real(wp)       :: lwp, thr
      integer        :: j
      logical        :: found(2)

      ! Initialization
      if ((chk%par%ch(1) <= 0) .or. (chk%par%ch(2) <= 0) .or. (chk%par%threshold%func == '')) &
           call finish(proc//'/lwp_diff','missing parameters')
      clear     = .true.
      ierr      = 0
      decr_hint = 'cld_check '//trim(chk%fnct)

      do j = 1, 2
        call get_chan(chk%par%ch(j), found(j), o=o(j), other_instrs=(/chk%par%instr(j)/), &
                      use_rawbt=chk%par%use_rawbt)
        if (.not.found(j) .and. chk%par%instr(j) < 0) &
             call get_chan(chk%par%ch(j), found(j), o=o(j), other_instrs=(/instr(i)/), &
                           use_rawbt=chk%par%use_rawbt)
      end do

      ! Check for cloudy spot
      if(all(found)) then
        lwp = lwp_KaSiRu94(o)
        thr = get_threshold(chk%par%threshold, s)
        if (thr /= inv_thr) then
          clear = (lwp < thr)
          ctmp%n_applied(ichk) = ctmp%n_applied(ichk) + 1
          if( .not.clear ) then
            call decr_tovs_use(s, obs%o(ib), CHK_CLOUD, tovs=tt, bit=USE_CLEAR, &
                 chan_mask=flagmask(:), hint_debug=decr_hint, tovs_flag=tovs_flag)
            call add_stat('cld',1, n_chan, bit=USE_CLEAR, chan_mask=flagmask(:))
          end if
        else
          clear = .false.
          ierr  = ERR_INV_THR
        end if
      else
        do j = 1, 2
          if (.not.found(j)) ierr = -chk%par%ch(j)
        end do
      end if
      if (ierr /= 0) call add_err(ierr,ichk)
      !
    end subroutine lwp_approx


    ! Scattering index cloud check
    subroutine scatindex(ierr)
      integer, intent(out)  :: ierr

      real(wp) :: o(2)       ! observed BT ( 150 GHz, 90 GHz )
      real(wp) :: thr
      integer  :: j
      logical  :: found(2)

      ! Initialization
      if ((chk%par%ch(1) <= 0) .or. (chk%par%ch(2) <= 0)) &
           call finish(proc//'/scatindex','missing parameters')

      clear     = .true.
      ierr      = 0
      decr_hint = 'cld_check '//trim(chk%fnct)
      do j = 1, 2
        call get_chan(chk%par%ch(j), found(j), o=o(j), other_instrs=(/chk%par%instr(j)/), &
             use_rawbt=chk%par%use_rawbt)
        ! if an instrument has different resolutions, for example for horizontal and vertical polarized
        ! channels, we are 'mapping' one resolution onto the other resolution of the same instrument
        ! so that we have two instruments in the satpp-file with the same instr-id
        ! then the following if-loop is necessary (for example for amsr2)
        if (.not.found(j) .and. chk%par%instr(j) < 0) &
             call get_chan(chk%par%ch(j), found(j), o=o(j), other_instrs=(/instr(i)/), &
                           use_rawbt=chk%par%use_rawbt)
        if (.not.found(j) .and. instr(i) == 3 .and. chk%par%ch(j) == 15 .and. rs%gopts%use_amsub_89ghz) then
           call get_chan(1, found(j), o=o(j), use_rawbt=chk%par%use_rawbt, other_instrs=(/4,15/))
           if (found(j)) ctmp%n_use_mhs1(ichk) = ctmp%n_use_mhs1(ichk) + 1
        end if
      end do

      ! Check for cloudy spot
      if (all(found)) then
        if (chk%par%threshold%func == '') then
          thr = get_threshold('sc', instr(i), rs%satid, s)
        else
          thr = get_threshold(chk%par%threshold, s)
        end if
        if (thr /= inv_thr) then
          clear = (o(2) - o(1) < thr)
          ctmp%n_applied(ichk) = ctmp%n_applied(ichk) + 1
          if( .not.clear ) then
            call decr_tovs_use(s, obs%o(ib), CHK_CLOUD, tovs=tt, bit=USE_CLEAR, &
                 chan_mask=flagmask(:), hint_debug=decr_hint, tovs_flag=tovs_flag)
            call add_stat('cld',1, n_chan, bit=USE_CLEAR, chan_mask=flagmask(:))
          end if
        else
          ierr = ERR_INV_THR
        end if
      else
        do j = 1, 2
          if (.not.found(j)) ierr = -chk%par%ch(j)
        end do
      end if
      if (ierr /= 0) call add_err(ierr,ichk)
      !
    end subroutine scatindex


    subroutine si_amsua_wrap(ierr)
      integer, intent(out)  :: ierr

      real(wp) :: o(3)       ! observed BT ( 150 GHz, 90 GHz )
      real(wp) :: thr
      integer  :: j
      logical  :: found(3)

      ! Initialization
      if ((chk%par%ch(1) <= 0) .or. (chk%par%ch(2) <= 0) .or. (chk%par%ch(3) <= 0)) &
           call finish(proc//'/si_amsua_wrap','missing parameters')

      found(1:3) = .false.
      clear      = .true.
      ierr       = 0
      decr_hint  = 'cld_check '//trim(chk%fnct)

      do j = 1, 3
        call get_chan(chk%par%ch(j), found(j), o=o(j), other_instrs=(/chk%par%instr(j)/), &
                      use_rawbt=chk%par%use_rawbt)
        if (.not.found(j) .and. chk%par%instr(j) < 0) &
             call get_chan(chk%par%ch(j), found(j), o=o(j), other_instrs=(/instr(i)/), &
                           use_rawbt=chk%par%use_rawbt)
        if (.not.found(j) .and. instr(i) == 3 .and. chk%par%ch(j) == 15 .and. rs%gopts%use_amsub_89ghz) then
          call get_chan(1, found(j), o=o(j), use_rawbt=chk%par%use_rawbt, other_instrs=(/4,15/))
          if (found(j)) ctmp%n_use_mhs1(ichk) = ctmp%n_use_mhs1(ichk) + 1
        end if
      end do

      ! Check for cloudy spot
      if (all(found(1:3))) then
        if (chk%par%threshold%func == '') then
          thr = get_threshold('sia', instr(i), rs%satid, s)
        else
          thr = get_threshold(chk%par%threshold, s)
        end if
        if (thr /= inv_thr) then
          clear = (si_amsua(o) < thr)
          ctmp%n_applied(ichk) = ctmp%n_applied(ichk) + 1
          if( .not.clear ) then
            call decr_tovs_use(s, obs%o(ib), CHK_CLOUD, tovs=tt, bit=USE_CLEAR, &
                 chan_mask=flagmask(:), hint_debug=decr_hint)
            call add_stat('cld',1, n_chan, bit=USE_CLEAR, chan_mask=flagmask(:))
          end if
        else
          ierr = ERR_INV_THR
        end if
      else
        do j = 1, 3
          if (.not.found(j)) ierr = -chk%par%ch(j)
        end do
      end if
      if (ierr /= 0) call add_err(ierr,ichk)
      !
    end subroutine si_amsua_wrap


    subroutine fgdep(ierr)
      integer, intent(out)  :: ierr

      real(wp) :: o  (1)      ! observed BT
      real(wp) :: fc (1)      ! first guess BT
      real(wp) :: thr
      logical  :: found1

      ! Initialization
      if ((chk%par%ch(1) <= 0) .or. (chk%par%threshold%func == '')) &
           call finish(proc//'/fgdep','missing parameters')
      found1    = .false.
      clear     = .true.
      ierr      = 0
      decr_hint = 'cld_check '//trim(chk%fnct)

      ! Find data
      call get_chan(chk%par%ch(1), found1, o=o(1), fc=fc(1), other_instrs=(/chk%par%instr(1)/), &
                    use_rawbt=chk%par%use_rawbt)
      if (.not.found1 .and. chk%par%instr(1) < 0) &
           call get_chan(chk%par%ch(1), found1, o=o(1), fc=fc(1), other_instrs=(/instr(i)/), &
                         use_rawbt=chk%par%use_rawbt)

      ! Check for cloudy spot
      if( found1 ) then
        ! call fgdep_cld_check(o, fc, chk%par%threshold, s, clear, ierr)
        thr = get_threshold(chk%par%threshold, s)
        if (thr /= inv_thr) then
          clear =  (abs(o(1) - fc(1)) < thr)
          ctmp%n_applied(ichk) = ctmp%n_applied(ichk) + 1
          if( .not.clear ) then
            call decr_tovs_use(s, obs%o(ib), CHK_CLOUD, tovs=tt, bit=USE_CLEAR, &
                 chan_mask=flagmask(:), hint_debug=decr_hint, tovs_flag=tovs_flag)
            call add_stat('cld', 1, n_chan, bit=USE_CLEAR, chan_mask=flagmask(:))
          end if
        else
          ierr = ERR_INV_THR
        end if
      else
        if (.not.found1) ierr = -chk%par%ch(1)
      end if
      if (ierr /= 0) call add_err(ierr,ichk)
      !
    end subroutine fgdep


    ! Flag always as cloudy (useful as fallback)
    subroutine always_cloudy(ierr)
      integer, intent(out) :: ierr
      ierr = 0
      ctmp%n_applied(ichk) = ctmp%n_applied(ichk) + 1
      call decr_tovs_use(s, obs% o(ib), CHK_CLOUD, tovs=tt, bit=USE_CLEAR, &
                         chan_mask=flagmask(:), hint_debug='cld_detect '//trim(FNCT_ALWAYS),&
                         tovs_flag=tovs_flag)
      call add_stat('cld', 1, n_chan, bit=USE_CLEAR, chan_mask=flagmask(:))
    end subroutine always_cloudy


    subroutine construct_cld_tmp
      character(len=300)   :: fstr = ''
      integer, allocatable :: ichan(:)
      integer              :: iset, i, j
      integer              :: nc, nchk, i_s, i_e, instr, n, stat, n_loop,m, n_chain
      logical              :: use_chk(n_chk), l_swap
      integer              :: chk_instr(n_chk), chk_aux(n_chk)

      allocate(cld_tmp(n_set))
      do iset = 1, n_set
        rs  => rad_set(iset)
        ctmp => cld_tmp (iset)
        nchk = 0
        use_chk(:) = .false.
        do i = 1, rs%n_instr
          do j = 1, n_chk
            if ( any(check(j)%instr == rs%instr(i)) .and. &
                 check(j)%phase == phase) then
              if (any(check(j)%satids >= 0)) then
                if (.not.any(check(j)%satids == rs%satid)) cycle
              end if
              if (any(check(j)%grids >= 0)) then
                if (.not.any(check(j)%grids == rs%grid)) cycle
              end if
              use_chk(j)   = .true.
              chk_instr(j) = rs%instr(i)
            end if
          end do
        end do
        nchk = count(use_chk(:))
        allocate(ctmp%chk(nchk))
        ctmp%chk (1:nchk) = pack( (/(j,j=1,n_chk)/), mask=use_chk(1:n_chk))
        chk_instr(1:nchk) = pack(chk_instr(1:n_chk), mask=use_chk(1:n_chk))

        ! Sort for increasing instruments (for beauty) and decreasing priority (important)
        do i = nchk-1,1,-1
          do j = 1, i
            if (chk_instr(j+1) == chk_instr(j)) then
              l_swap = check(ctmp%chk(j+1))%priority > check(ctmp%chk(j))%priority
            else
              l_swap = chk_instr(j+1) < chk_instr(j)
            end if
            if (l_swap) then
              m = chk_instr(j+1) ; chk_instr(j+1) = chk_instr(j) ; chk_instr(j) = m
              m = ctmp%chk (j+1) ; ctmp%chk (j+1) = ctmp%chk (j) ; ctmp%chk (j) = m
            end if
          end do
        end do
        ! Make sure, that checks with the same "chain" are following each other
        i = 1
        do while (i <= nchk)
          chk => check(ctmp%chk(i))
          if (chk%chain /= '') then
            n_chain = count(check(ctmp%chk(1:nchk))%chain == chk%chain)
            if (n_chain > 1) then
              chk_aux(i+1:nchk) = ctmp%chk(i+1:nchk)
              ctmp%chk(i+1:i+n_chain-1) = pack(chk_aux(i+1:nchk), &
                   mask=check(chk_aux(i+1:nchk))%chain == chk%chain)
              ctmp%chk(i+n_chain:nchk)  = pack(chk_aux(i+1:nchk), &
                   mask=check(chk_aux(i+1:nchk))%chain /= chk%chain)
              i = i + n_chain - 1
            end if
          end if
          i = i + 1
        end do

        ctmp%n_chk             = nchk
        ctmp%n_fovs            = 0
        ctmp%nerr              = 0
        ctmp%n_im_flag         = 0
        nc = rs%n_chan
        allocate(ctmp%chk_done   (   nchk                          ))
        allocate(ctmp%n_data     (   nchk                          )) ; ctmp%n_data     = 0
        allocate(ctmp%n_notappl  (   nchk,CLD_CHK_SEA:notappl_instr)) ; ctmp%n_notappl  = 0
        allocate(ctmp%n_applied  (   nchk                          )) ; ctmp%n_applied  = 0
        allocate(ctmp%n_failed   (   nchk                          )) ; ctmp%n_failed   = 0
        allocate(ctmp%n_cld_chan (nc,nchk                          )) ; ctmp%n_cld_chan = 0
        allocate(ctmp%n_srf_chan (nc,nchk                          )) ; ctmp%n_srf_chan = 0
        allocate(ctmp%n_aer_chan (nc,nchk                          )) ; ctmp%n_aer_chan = 0
        allocate(ctmp%n_trg_chan (nc,nchk                          )) ; ctmp%n_trg_chan = 0
        allocate(ctmp%n_use_mhs1 (   nchk                          )) ; ctmp%n_use_mhs1 = 0
        allocate(ctmp%flagmask   (nc,nchk                          )) ; ctmp%flagmask   = .false.
        allocate(ctmp%chi        (nc                               ))
        do i = 1, rs%n_instr
          ctmp%chi(rs%o_ch_i(i)+1:rs%o_ch_i(i)+rs%n_ch_i(i))%instr = rs%instr(i)
        end do
        ctmp%chi(1:nc)%value = rs%chan(1:nc)
        ! Convert t_cld_check%flagchan into t_cld_tmp%flagmask
        allocate(ichan(nc))
        do j = 1, nchk
          chk => check(ctmp%chk(j))
          chk%flagchan = trim(adjustl(tolower(chk%flagchan)))
          if (chk%flagchan == 'all') then
            ctmp%flagmask(:,j) = .true.
          elseif (chk%flagchan == 'instr') then
            do i = 1,nc
              if (any(ctmp%chi(i)%instr == chk%instr(:))) then
                ctmp%flagmask(i,j) = .true.
              else
                ctmp%flagmask(i,j) = .false.
              end if
            end do
          else
            i_s = 1
            i_e = 0
            n_loop = 0
            do while (i_s <= len_trim(chk%flagchan) .and. n_loop <=100)
              n_loop = n_loop + 1
              i = index(chk%flagchan(i_s:), 'instr=')
              if (i <= 0) then
                instr = -1
                i_e = len_trim(chk%flagchan)
              else
                i_s = i_s + 6
                i = scan(chk%flagchan(i_s:), '0123456789')
                i_s = i_s + i - 1
                read(chk%flagchan(i_s:),*,iostat=stat) instr
                i = scan(chk%flagchan(i_s:), ' ,;')
                if (i == 0) then
                  i_s = len_trim(chk%flagchan) + 1
                else
                  i_s = i_s + i
                end if
                i_e = index(chk%flagchan(i_s:), 'instr=')
                if (i_e == 0) then
                  i_e = len_trim(chk%flagchan)
                else
                  i_e = i_s + i_e - 2
                end if
              end if
              if (instr > 0 .and. all(ctmp%chi(:)%instr /= instr)) then
                ! Avoid wrong size of ichan (in intstr2array) for multiple-instrument rule in flagchan
                i_s = i_e + 1
                cycle
              end if
              if (i_e >= i_s) then
                if (adjustl(chk%flagchan(i_s:i_e)) == 'all') then
                  ichan = ctmp%chi(:)%value
                  n = nc
                else
                  call intstr2array(chk%flagchan(i_s:i_e),' ,;',n=n,iarr=ichan, status=stat, fstr=fstr)
                  if (stat == 2) then
                    write(0,'("WARNING construct_cld_tmp/intstr2array (satid ",I3.3,&
                         &" grid=",I3.3,"): ",A)') rs%satid, rs%grid, trim(fstr)
                  elseif (stat /= 0) then
                    call finish('construct_cld_tmp/intstr2array',trim(fstr))
                  end if
                end if
              elseif (instr >= 0) then
                ichan = ctmp%chi(:)%value
                n = nc
              end if
              do i = 1, nc
                if (instr >= 0) then
                  if (ctmp%chi(i)%instr /= instr) cycle
                end if
                if (any(ctmp%chi(i)%value == ichan(1:n))) ctmp%flagmask(i,j) = .true.
              end do
              i_s = i_e + 1
            end do
          end if
        end do
        deallocate(ichan)
      end do
    end subroutine construct_cld_tmp

    subroutine destruct_cld_tmp
      integer :: iset
      do iset = 1, size(cld_tmp)
        ctmp => cld_tmp(iset)
        if (allocated (ctmp%chk       )) deallocate(ctmp%chk       )
        if (allocated (ctmp%chk_done  )) deallocate(ctmp%chk_done  )
        if (allocated (ctmp%n_data    )) deallocate(ctmp%n_data    )
        if (allocated (ctmp%n_notappl )) deallocate(ctmp%n_notappl )
        if (allocated (ctmp%n_applied )) deallocate(ctmp%n_applied )
        if (allocated (ctmp%n_failed  )) deallocate(ctmp%n_failed  )
        if (allocated (ctmp%n_cld_chan)) deallocate(ctmp%n_cld_chan)
        if (allocated (ctmp%n_srf_chan)) deallocate(ctmp%n_srf_chan)
        if (allocated (ctmp%n_aer_chan)) deallocate(ctmp%n_aer_chan)
        if (allocated (ctmp%n_trg_chan)) deallocate(ctmp%n_trg_chan)
        if (allocated (ctmp%n_use_mhs1)) deallocate(ctmp%n_use_mhs1)
        if (allocated (ctmp%chi       )) deallocate(ctmp%chi       )
        if (associated(ctmp%err       )) deallocate(ctmp%err       )
        if (allocated (ctmp%flagmask  )) deallocate(ctmp%flagmask  )
      end do
    end subroutine destruct_cld_tmp

    subroutine add_err(code,ichk,descr,n)
      integer,          intent(in)            :: code
      integer,          intent(in)            :: ichk
      character(len=*), intent(in),  optional :: descr
      integer,          intent(out), optional :: n

      type(t_err), pointer :: earr(:) => null()
      integer              :: ind, i
      character(len=80)    :: descr_loc = ''

      ctmp%n_failed(ichk) = ctmp%n_failed(ichk) + 1

      ind = -1
      do i = 1, ctmp%nerr
        if ( ctmp%err(i)%code  == code  .and. &
             ctmp%err(i)%ichk  == ichk ) then
          ind = i
          exit
        end if
      end do
      if (ind <=0) then
        allocate(earr(ctmp%nerr+1))
        if  (ctmp%nerr > 0) earr(1:ctmp%nerr) = ctmp%err(1:ctmp%nerr)
        if (associated(ctmp%err)) deallocate(ctmp%err)
        ctmp%err => earr
        ctmp%nerr = ctmp%nerr + 1
        ind = ctmp%nerr
        descr_loc = ''
        if (present(descr))  descr_loc = descr(1:min(len_trim(descr),len(descr_loc)))
        if (descr_loc == '') call cloud_detect_errmsg(code, descr_loc)
        ctmp%err(ind) = t_err(0, code, ichk, descr_loc)
      end if
      ctmp%err(ind)%n = ctmp%err(ind)%n + 1
      if (present(n)) n = ctmp%err(ind)%n

    end subroutine add_err

    subroutine add_stat(stat,is,ie,bit,not_bit,chan_mask)
      character(len=*), intent(in)           :: stat
      integer,          intent(in)           :: is  ! start channel
      integer,          intent(in)           :: ie  ! end channel
      integer,          intent(in), optional :: bit
      integer,          intent(in), optional :: not_bit
      logical,          intent(in), optional :: chan_mask(:)

      integer :: j, k, ii

      ii = 1
      do j = is, ie
        do while (j > o_ch_i(ii) + n_ch_i(ii) .and. ii < n_instr)
          ii = ii + 1
        end do
        if (present(bit)) then
          if (.not.btest(flag(j), bit)) cycle
        end if
        if (present(not_bit)) then
          if (btest(flag(j), not_bit)) cycle
        end if
        if (present(chan_mask)) then
          if (.not.chan_mask(j-is+1)) cycle
        end if
        do k = j, size(ctmp%chi)
          if (t_ilev(chan(j),instr(ii)) == ctmp%chi(k)) then
            select case(stat)
            case('cld')
              ctmp%n_cld_chan(k,ichk) = ctmp%n_cld_chan(k,ichk) + 1
            case('srf')
              ctmp%n_srf_chan(k,ichk) = ctmp%n_srf_chan(k,ichk) + 1
            case('aer')
              ctmp%n_aer_chan(k,ichk) = ctmp%n_aer_chan(k,ichk) + 1
            case('trg')
              ctmp%n_trg_chan(k,ichk) = ctmp%n_trg_chan(k,ichk) + 1
            case default
              call finish(proc//'/add_stat','invalid stat='//trim(stat))
            end select
            exit
          end if
        end do
      end do
    end subroutine add_stat

    subroutine print_stat
      type(t_err), pointer :: err(:) => null()
      type(t_err), pointer :: err_snd(:) => null()
      type(t_err), target  :: err_dum(0)
      character(len=300)   :: msg, msg_
      integer              :: iset, ichk, i, j, n, nch
      integer              :: nerr, ierr
      logical              :: l_srf, l_aer, l_trg
      if (dace% lpio) then
        write(6,*)
        write(6,'(1x,A,I1," (",A,")")') 'Cloud detection statistics phase ',phase,get_c_phase(phase)
        write(6,*) '------------------------------------------------'
      end if
      do iset = 1, size(cld_tmp)
        ctmp => cld_tmp(iset)
        rs  => rad_set(iset)
        ctmp%n_fovs            = p_sum(ctmp%n_fovs    )
        ctmp%n_data            = p_sum(ctmp%n_data    )
        ctmp%n_notappl         = p_sum(ctmp%n_notappl )
        ctmp%n_applied         = p_sum(ctmp%n_applied )
        ctmp%n_failed          = p_sum(ctmp%n_failed  )
        ctmp%n_cld_chan        = p_sum(ctmp%n_cld_chan)
        ctmp%n_srf_chan        = p_sum(ctmp%n_srf_chan)
        ctmp%n_aer_chan        = p_sum(ctmp%n_aer_chan)
        ctmp%n_trg_chan        = p_sum(ctmp%n_trg_chan)
        ctmp%n_use_mhs1        = p_sum(ctmp%n_use_mhs1)
        ctmp%n_im_flag         = p_sum(ctmp%n_im_flag )
        ! collect error messages on p_io
        nerr = p_sum(ctmp%nerr)
        if (nerr > 0) then
          if (.not.dace% lpio) nerr = 0
          allocate(err(nerr))
          if (ctmp%nerr > 0) then
            err_snd => ctmp%err(1:ctmp%nerr)
          else
            err_snd => err_dum
          end if
          call p_gather_err(err_snd, err(1:nerr), dace% pio)
          if (dace% lpio) then
            do i = 1, nerr
              do j = i+1,nerr
                if (err(i)%code==err(j)%code .and. err(i)%ichk==err(j)%ichk) then
                  err(i)%n = err(i)%n + err(j)%n
                  err(j)%code = 0
                end if
              end do
            end do
            i = count(err(1:nerr)%code /= 0)
            err(1:i) = pack(err(1:nerr), mask=err(1:nerr)%code/=0)
            if (associated(ctmp%err)) deallocate(ctmp%err)
            ctmp%err  => err
            ctmp%nerr =  i
            err => null()
          else
            deallocate(err)
          end if
        end if

        if (dace% lpio) then
          if (any(ctmp%n_data > 0)) then
            write(6,'(3x,"satellite ",I3.3,"  grid ",I3.2)') rs%satid, rs%grid
            write(6,'(5x,"#FOVs",T45,": ",I9)') ctmp%n_fovs
            if (ctmp%n_fovs <= 0) cycle
            ! Checks
            do ichk = 1, ctmp%n_chk
              if (ctmp%chk(ichk) <= 0 .or. ctmp%chk(ichk) > n_chk) CYCLE
              chk => check(ctmp%chk(ichk))
              write(6,'(5x,I1,". check: ",A,"  (instr=",I3.3,", priority=",I2,")")') &
                   ichk, trim(chk%description), chk%instr(1), chk%priority  !TODO: all instruments
              do i = notappl_instr, CLD_CHK_SEA, -1
                if (ctmp%n_notappl(ichk,i) > 0) &
                     write(6,'(7x,"#Not applied (",A,")",T45,":",2x,I9,2x,F7.3,"%")') trim(c_notappl(i)), &
                     ctmp%n_notappl(ichk,i),(100.*ctmp%n_notappl(ichk,i))/ctmp%n_fovs
              end do
              write(6,'(7x,"#Tests  (1)",T45,":",2x,I9,2x,F7.3,"%")') ctmp%n_data(ichk),(100.*ctmp%n_data(ichk))/ctmp%n_fovs
              ! Statistics on test
              if (ctmp%n_data(ichk) > 0) then
                write(6,'(7x,"#Successful tests  (2)",T45,":",2x,I9,2x,F7.3,"%")') &
                     ctmp%n_applied(ichk),(100.*ctmp%n_applied(ichk))/ctmp%n_data(ichk)
                write(6,'(7x,"#Failed tests",T45,":",2x,I9,2x,F7.3,"%")') &
                     ctmp%n_failed(ichk),(100.*ctmp%n_failed(ichk))/ctmp%n_data(ichk)
                do ierr = 1, ctmp%nerr
                  if (ctmp%err(ierr)%ichk == ichk) then
                    write(6,'(7x,"#Failed tests with errcode=",I4.3,T45,":",2x,I9,2x,F7.3,"%  (",A,")")') &
                         ctmp%err(ierr)%code,ctmp%err(ierr)%n,(100.*ctmp%err(ierr)%n)/ctmp%n_data(ichk),&
                         trim(ctmp%err(ierr)%descr)
                  end if
                end do
                if (any(chk%instr == 3) .and. ctmp%n_use_mhs1(ichk) > 0) then
                  write(6,'(7x,"#use mhs1 instead of amsua15",T45,":",2x,I9,2x,F7.3,"%")') &
                       ctmp%n_use_mhs1(ichk),(100.*ctmp%n_use_mhs1(ichk))/ctmp%n_data(ichk)
                end if
                if (chk%fnct == FNCT_MNW .and. any(ctmp%n_im_flag > 0)) then
                  n = sum(ctmp%n_im_flag(1:))
                  write(6,'(7x,"#coloc. imag. flag > 0",T45,":",2x,I9,2x,F7.3,"%")') &
                       n,(100.*n)/ctmp%n_data(ichk)
                  do i = 0, im_fl_nb-1
                    write(6,'(7x,"#coloc. imag. bit ",I1," (",I1,") meaning",T45,":",2x,A)') &
                         i,2**i,im_fl_name(i+1)
                  end do
                  write(6,'(27x,"bits 0 1 2 3")')
                  do i = 0, im_fl_u
                    write(6,'(7x,"#coloc. imag. flag=",I2,4x,A,T45,":",2x,I9,2x,F7.3,"%")') &
                         i,im_fl_str(i),ctmp%n_im_flag(i),(100.*ctmp%n_im_flag(i))/ctmp%n_data(ichk)
                  end do
                  do i = 0, im_fl_nb-1
                    n = 0
                    do j = 0, im_fl_u
                      if (btest(j,i)) n = n + ctmp%n_im_flag(j)
                    end do
                    write(6,'(7x,"#coloc. imag. bit ",I1," (",I1,") ",A," total",T45,":",2x,I9,2x,F7.3,"%")') &
                         i,2**i,im_fl_name(i+1),n,(100.*n)/ctmp%n_data(ichk)
                  end do
                end if
                ! Flags
                if (any(ctmp%n_cld_chan(:,ichk) > 0) .or. any(ctmp%n_srf_chan(:,ichk) > 0)) then
                  ! header line for channel statistics
                  if (ctmp%n_applied(ichk) > 0) then
                    write(msg,'(7x,"          instr   chan",T45,"       cloud     %FOVs %tests(2)")')
                    n = ctmp%n_applied(ichk)
                  else
                    write(msg,'(7x,"          instr   chan",T45,"       cloud     %FOVs %tests(1)")')
                    n = ctmp%n_data(ichk)
                  end if
                  nch = size(ctmp%chi)
                  l_srf = any(ctmp%n_srf_chan(1:nch,ichk) > 0)
                  if (l_srf) msg = trim(msg)//'    surface     %FOVs %tests(1)'
                  l_aer = any(ctmp%n_aer_chan(1:nch,ichk) > 0)
                  if (l_aer) msg = trim(msg)//'    aerosol     %FOVs %tests(1)'
                  l_trg = any(ctmp%n_trg_chan(1:nch,ichk) > 0)
                  if (l_trg) msg = trim(msg)//'   tracegas     %FOVs %tests(1)'
                  write(6,*) trim(msg)
                  ! Statistics for each channel
                  do i = 1, nch
                    if ( ctmp%n_cld_chan(i,ichk) > 0 .or. &
                         ctmp%n_srf_chan(i,ichk) > 0 .or. &
                         ctmp%n_aer_chan(i,ichk) > 0.and.l_aer .or. &
                         ctmp%n_trg_chan(i,ichk) > 0.and.l_trg) then
                      write(msg,'(7x,"#flagged    ",I3.3,2x,I5.5,T45,":",(2x,I9,2(2x,F7.3,"%")))') &
                           ctmp%chi(i)%instr, ctmp%chi(i)%value,                                   &
                           ctmp%n_cld_chan(i,ichk),(100.*ctmp%n_cld_chan(i,ichk))/ctmp%n_fovs,     &
                                                   (100.*ctmp%n_cld_chan(i,ichk))/n
                      if (l_srf) then
                        write(msg_,'(2x,I9,2(2x,F7.3,"%"))') &
                             ctmp%n_srf_chan(i,ichk),(100.*ctmp%n_srf_chan(i,ichk))/ctmp%n_fovs, &
                                                     (100.*ctmp%n_srf_chan(i,ichk))/n
                        msg = trim(msg)//trim(msg_)
                      end if
                      if (l_aer) then
                        write(msg_,'(2x,I9,2(2x,F7.3,"%"))') &
                             ctmp%n_aer_chan(i,ichk),(100.*ctmp%n_aer_chan(i,ichk))/ctmp%n_fovs, &
                                                     (100.*ctmp%n_aer_chan(i,ichk))/n
                        msg = trim(msg)//trim(msg_)
                      end if
                      if (l_trg) then
                        write(msg_,'(2x,I9,2(2x,F7.3,"%"))') &
                             ctmp%n_trg_chan(i,ichk),(100.*ctmp%n_trg_chan(i,ichk))/ctmp%n_fovs, &
                                                     (100.*ctmp%n_trg_chan(i,ichk))/n
                        msg = trim(msg)//trim(msg_)
                      end if
                      write(6,*) trim(msg)
                    end if
                  end do
                end if
              end if
            end do
            write(6,*) ''
          end if
        end if
      end do
    end subroutine print_stat

    function im_fl_str(im_fl) result(s)
      character(len=2*im_fl_nb-1) :: s
      integer, intent(in) :: im_fl
      integer :: i, is
      s = repeat(' ',2*im_fl_nb-1)
      do i = 0, im_fl_nb-1
        is = (i+1)*2-1
        if (btest(im_fl,i)) write(s(is:is),'(I1)') i
      end do
    end function im_fl_str

  end subroutine cloud_detect


  subroutine cloud_detect_init
  !-----------------------------------------------
  ! read namelists /TOVS_CLOUD/ and /TOVS_CLOUD_CHECK/
  !-----------------------------------------------
    integer                 :: ierr
    integer                 :: i, j
    type(t_cld_check)       :: cdum
    logical                 :: l_swap

    ! Reading of TOVS_CLOUD_CHECK namelists
    integer                 :: inml, ios
    character(len=80)       :: msg = ''
    type(t_cld_check)       :: cld_chk_read
    character(len=len_fnct) :: amsub_fnct

    ! TOVS_CLOUD_CHECK namelist variables
    type t_fg_chk  ! deprecated: just for backwards compatibility with old namelists
      real(wp)              :: threshold = -huge(0._wp)
      integer               :: chan      = -1
    end type t_fg_chk
    integer                 :: instr(ninstr) = -1
    integer                 :: satids(nsat)  = -1
    integer                 :: grids(ngrid)  = -1
    integer                 :: phase         = -99
    character(len=10)       :: c_phase       = ''
    integer                 :: priority      = -99
    character(len=80)       :: description   = ''
    character(len=80)       :: chain         = ''
    character(len=len_fnct) :: fnct          = ''
    integer                 :: flags         = ibset(0,CLD_CHK_LAND) + ibset(0,CLD_CHK_SEA)
    character(len=1000)     :: flagchan      = ''
    integer                 :: tf_cloud_ind(ninstr)= -1
    ! parameters for special functions
    type(t_chk_par)         :: par_
    type(t_fg_chk)          :: fg  ! deprecated: just for backwards compatibility with old namelists
                                   ! the content is stored into t_chk_par
    ! Attach tf_cloud_ind to each test
    logical                 :: lused(size(TF_CLOUD))
    integer                 :: instrs(100)
    integer                 :: n_instr, ins, ii, ibit


    !---------------------------
    ! 1. read namelist TOVS_CLOUD
    !---------------------------
    ! 1.a set defaults
    !------------------
    verbose                  = 0        ! verbosity level
    do_clddetect_irs         = .false.  ! do IRS cloud detection in 3dvar
    do_clddetect_amsua       = .true.   ! do AMSU-A cloud detection in 3dvar
    do_clddetect_atms        = .true.   ! do ATMS cloud detection in 3dvar
    do_clddetect_atms_amsub  = .false.  ! do AMSU-B-like cloud detection in 3dvar for ATMS humidity channels
    do_clddetect_iasi        = .true.   ! do IASI   cloud detection in 3dvar
    do_clddetect_hirs        = .true.   ! do HIRS   cloud detection in 3dvar
    do_clddetect_cris        = .true.   ! do CrIS   cloud detection in 3dvar
    do_clddetect_mwts        = .false.  ! do MWTS   cloud detection in 3dvar
    do_clddetect_gmi         = .false.  ! do GMI    cloud detection in 3dvar
    do_clddetect_ssmis       = .false.  ! do SSMI/S cloud detection in 3dvar
    do_clddetect_amsr2       = .false.  ! do AMSR2  cloud detection in 3dvar
    do_clddetect_saphir      = .false.  ! do SAPHIR cloud detection in 3dvar
    do_clddetect_mwhs2       = .false.  ! do MWHS2  cloud detection in 3dvar
    do_clddetect_mhs         = .false.  ! do MHS    cloud detection in 3dvar
    do_clddetect_mwts3       = .false.  ! do MWTS-3 cloud detection in 3dvar
    use_amsua_surftype       = .true.   ! use surface types 1-12 from amsua cloud check
    use_atms_surftype        = .true.   ! use surface types 1-12 from amsua cloud check
    use_mwts3_surftype       = .true.   ! use surface types 1-12 from mwts3 cloud check
    use_amsua_ice            = .false.  ! use surface types 1-12 from amsua cloud check
    use_atms_ice             = .false.  ! use surface types 1-12 from amsua cloud check
    use_mwts3_ice            = .false.  ! use surface types 1-12 from mwts3 cloud check
    toplev                   = 0._wp    ! cloud search top level
    bottomlev                = 1.E20_wp ! cloud search bottom level
    param_file_scatt         = ''
    param_file_si_amsua      = ''
    ! clchk_amsub, param_fil_amsub, amsub_ec_bnd, amsub_delta_ch20_ch18 are not initialized,
    ! because they are still part of the /TOVS_OBS/ namelist for compatibility reasons

    !------------------
    ! 1.b read namelist
    !------------------
    if (dace% lpio) then
      call position_nml ('TOVS_CLOUD', status=ierr)
      select case (ierr)
      case (POSITIONED)
#if defined(__ibm__)
        read (nnml ,nml=TOVS_CLOUD ,iostat=ios)
        if (ios/=0) call finish ('cloud_detect_init',             &
                                  'ERROR in namelist /TOVS_CLOUD/')
#else
        read (nnml ,nml=TOVS_CLOUD)
#endif
      end select
    endif

    !--------------------------
    ! 1.c broadcast namelist
    !--------------------------
    call p_bcast (ierr                    ,dace% pio)
    call p_bcast (verbose                 ,dace% pio)
    call p_bcast (do_clddetect_irs        ,dace% pio)
    call p_bcast (do_clddetect_amsua      ,dace% pio)
    call p_bcast (do_clddetect_atms       ,dace% pio)
    call p_bcast (do_clddetect_atms_amsub ,dace% pio)
    call p_bcast (do_clddetect_iasi       ,dace% pio)
    call p_bcast (do_clddetect_hirs       ,dace% pio)
    call p_bcast (do_clddetect_cris       ,dace% pio)
    call p_bcast (do_clddetect_mwts       ,dace% pio)
    call p_bcast (do_clddetect_gmi        ,dace% pio)
    call p_bcast (do_clddetect_ssmis      ,dace% pio)
    call p_bcast (do_clddetect_saphir     ,dace% pio)
    call p_bcast (do_clddetect_mwhs2      ,dace% pio)
    call p_bcast (do_clddetect_mhs        ,dace% pio)
    call p_bcast (do_clddetect_amsr2      ,dace% pio)
    call p_bcast (do_clddetect_mwts3      ,dace% pio)
    call p_bcast (use_amsua_surftype      ,dace% pio)
    call p_bcast (use_atms_surftype       ,dace% pio)
    call p_bcast (use_mwts3_surftype      ,dace% pio)
    call p_bcast (use_amsua_ice           ,dace% pio)
    call p_bcast (use_atms_ice            ,dace% pio)
    call p_bcast (use_mwts3_ice           ,dace% pio)
    call p_bcast (topLev                  ,dace% pio)
    call p_bcast (bottomLev               ,dace% pio)
    call p_bcast (clchk_amsub             ,dace% pio)
    call p_bcast (param_file_amsub        ,dace% pio)
    call p_bcast (amsub_ec_bnd            ,dace% pio)
    call p_bcast (amsub_delta_ch20_ch18   ,dace% pio)
    call p_bcast (dlat_thresh_def         ,dace% pio)
    call p_bcast (dsz_thresh_def          ,dace% pio)
    call p_bcast (param_file_scatt        ,dace% pio)
    call p_bcast (param_file_si_amsua     ,dace% pio)
    call p_bcast (buehler_version         ,dace% pio)
    !------------------------------------
    ! printout of some namelist variables
    !------------------------------------
    if (dace% lpio) then
      write(6,'(a)') repeat('-',79)
      write(6,'( )')
      if (ierr == POSITIONED) then
        write(6,'(a)') '  namelist /TOVS_CLOUD/ read'
      else
        write(6,'(a)') '  namelist /TOVS_CLOUD/ not present, using defaults'
      end if
      write(6,'( )')
      write(6,'( )')
      write(6,'(a, i6)')         '   verbose                  =',verbose
      write(6,'(a,5x,l1)')       '   do_clddetect_irs         =',do_clddetect_irs
      write(6,'(a,5x,l1)')       '   do_clddetect_amsua       =',do_clddetect_amsua
      write(6,'(a,5x,l1)')       '   do_clddetect_atms        =',do_clddetect_atms
      write(6,'(a,5x,l1)')       '   do_clddetect_atms_amsub  =',do_clddetect_atms_amsub
      write(6,'(a,5x,l1)')       '   do_clddetect_iasi        =',do_clddetect_iasi
      write(6,'(a,5x,l1)')       '   do_clddetect_cris        =',do_clddetect_cris
      write(6,'(a,5x,l1)')       '   do_clddetect_hirs        =',do_clddetect_hirs
      write(6,'(a,5x,l1)')       '   do_clddetect_mwts        =',do_clddetect_mwts
      write(6,'(a,5x,l1)')       '   do_clddetect_gmi         =',do_clddetect_gmi
      write(6,'(a,5x,l1)')       '   do_clddetect_ssmis       =',do_clddetect_ssmis
      write(6,'(a,5x,l1)')       '   do_clddetect_saphir      =',do_clddetect_saphir
      write(6,'(a,5x,l1)')       '   do_clddetect_mwhs2       =',do_clddetect_mwhs2
      write(6,'(a,5x,l1)')       '   do_clddetect_mhs         =',do_clddetect_mhs
      write(6,'(a,5x,l1)')       '   do_clddetect_amsr2       =',do_clddetect_amsr2
      write(6,'(a,5x,l1)')       '   do_clddetect_mwts3       =',do_clddetect_mwts3
      write(6,'(a)')             '   ! MW Emissivity cloud check'
      write(6,'(a,5x,20(l1,1x))')'   use_amsua_surftype       =',use_amsua_surftype(1:n_surf)
      write(6,'(a,5x,20(l1,1x))')'   use_atms_surftype        =',use_atms_surftype(1:n_surf)
      write(6,'(a,5x,20(l1,1x))')'   use_mwts3_surftype       =',use_mwts3_surftype(1:n_surf)
      write(6,'(a,5x,20(l1,1x))')'   use_amsua_ice            =',use_amsua_ice(1:n_surf)
      write(6,'(a,5x,20(l1,1x))')'   use_atms_ice             =',use_atms_ice(1:n_surf)
      write(6,'(a,5x,20(l1,1x))')'   use_mwts3_ice            =',use_mwts3_ice(1:n_surf)
      write(6,'(a)')             '   ! McNally-Watts cloud check'
      write(6,'(a, E13.6)')      '   topLev                   =',topLev
      write(6,'(a, E13.6)')      '   bottomLev                =',bottomLev
      write(6,'(a)')             '   ! AMSUB-like cloud checks'
      write(6,'(a, i6)')         '   clchk_amsub              =',clchk_amsub
      write(6,'(a, a)')          '   param_file_amsub         =',trim(param_file_amsub)
      write(6,'(a, 2(1x,F7.3))') '   amsub_ec_bnd             =',amsub_ec_bnd
      write(6,'(a,  (1x,F7.3))') '   amsub_delta_ch20_ch18    =',amsub_delta_ch20_ch18
      write(6,'(a,  (1x,F7.4))') '   dlat_thresh_def          =',dlat_thresh_def
      write(6,'(a,  (1x,F7.4))') '   dsz_thresh_def           =',dsz_thresh_def
      write(6,'(a, i6)')         '   buehler_version          =',buehler_version
      write(6,'(a)')             '   ! scattering index test'
      write(6,'(a, a)')          '   param_file_scatt         =',trim(param_file_scatt)
      write(6,'(a, a)')          '   param_file_si_amsua      =',trim(param_file_si_amsua)
      write(6,'( )')
    endif

    !-------------------------------
    ! 2. Set up default cloud checks
    !-------------------------------
    n_chk = 0
    if (do_clddetect_irs    ) n_chk = n_chk + 2
    if (do_clddetect_amsua  ) n_chk = n_chk + 2
    if (do_clddetect_atms   ) n_chk = n_chk + 2
    if (do_clddetect_iasi   ) n_chk = n_chk + 2
    if (do_clddetect_cris   ) n_chk = n_chk + 4
    if (do_clddetect_mwts   ) n_chk = n_chk + 2
    if (do_clddetect_gmi    ) n_chk = n_chk + 4
    if (do_clddetect_ssmis  ) n_chk = n_chk + 4
    if (do_clddetect_saphir ) n_chk = n_chk + 2
    if (do_clddetect_mwhs2  ) n_chk = n_chk + 3
    if (do_clddetect_mhs    ) n_chk = n_chk + 9
    if (do_clddetect_amsr2  ) n_chk = n_chk + 2
    if (abs(clchk_amsub) > 0) n_chk = n_chk + 2
    if (do_clddetect_hirs   ) n_chk = n_chk + 1
    if (do_clddetect_mwts3  ) n_chk = n_chk + 2

    allocate(check(n_chk))

    ! Phase "after bias correction" defaults
    n_chk = 0
    if (do_clddetect_irs) then
      call add_check(instr=(/57/),                 phase=PH_ABC,   priority= 0,             &
                     description='McNally-Watts IRS cloud check',  fnct=trim(FNCT_MNW))
      call add_check(instr=(/57/),                 phase=PH_ABC,   priority=-1,             &
                     description='Always cloudy fallback',         fnct=trim(FNCT_ALWAYS))
    end if
    if (do_clddetect_amsua) then
      call add_check(instr=(/ 3/),                 phase=PH_ABC,   priority= 0,             &
                     description='AMSU-A emissivity cloud check' , fnct=trim(FNCT_EMISS),   &
                     flags=ibset(0,CLD_CHK_SEA))
      call add_check(instr=(/ 3/),                 phase=PH_ABC,   priority=-1,             &
                     description='Always cloudy fallback',         fnct=trim(FNCT_ALWAYS))
    end if
    if (do_clddetect_atms) then
      if (do_clddetect_atms_amsub) then
        flagchan = '1-16'
      else
        flagchan = empty_chk%flagchan
      end if
      call add_check(instr=(/19/),                 phase=PH_ABC,   priority= 0,             &
                     description='ATMS emissivity cloud check',    fnct=trim(FNCT_EMISS),   &
                     flags=ibset(0,CLD_CHK_SEA),                   flagchan=trim(flagchan))
      call add_check(instr=(/19/),                 phase=PH_ABC,   priority=-1,             &
                     description='Always cloudy fallback',         fnct=trim(FNCT_ALWAYS),  &
                     flagchan=trim(flagchan))
    end if
    if (do_clddetect_iasi) then
      call add_check(instr=(/16/),                 phase=PH_ABC,   priority= 0,             &
                     description='McNally-Watts IASI cloud check', fnct=trim(FNCT_MNW))
      call add_check(instr=(/16/),                 phase=PH_ABC,   priority=-1,             &
                     description='Always cloudy fallback',         fnct=trim(FNCT_ALWAYS))
    end if
    if (do_clddetect_cris) then
      call add_check(instr=(/27/),                 phase=PH_ABC,   priority= 0,             &
                     description='McNally-Watts CrIS cloud check', fnct=trim(FNCT_MNW))
      call add_check(instr=(/27/),                 phase=PH_ABC,   priority=-1,             &
                     description='Always cloudy fallback',         fnct=trim(FNCT_ALWAYS))
      call add_check(instr=(/28/),                 phase=PH_ABC,   priority= 0,             &
                     description='McNally-Watts CrIS-FSR cloud check', fnct=trim(FNCT_MNW))
      call add_check(instr=(/28/),                 phase=PH_ABC,   priority=-1,             &
                     description='Always cloudy fallback',         fnct=trim(FNCT_ALWAYS))
    end if
    if (do_clddetect_mwts) then
      call add_check(instr=(/72/),                 phase=PH_ABC,   priority= 0,             &
                     description='MWTS first guess cloud check',   fnct=trim(FNCT_FGDEP))
      call add_check(instr=(/72/),                 phase=PH_ABC,   priority=-1,             &
                     description='Always cloudy fallback',         fnct=trim(FNCT_ALWAYS))
    end if
    if (do_clddetect_gmi) then
      call add_check(instr=(/71/),                 phase=PH_ABC,   priority= 0,             &
                     description='GMI polarization difference',    fnct=trim(FNCT_POL_DIFF),&
                     par=t_chk_par('0.95',empty_trf,(/6,7,-1/),(/-1,-1,-1/),.FALSE.))
      call add_check(instr=(/71/),                 phase=PH_ABC,   priority= 0,             &
                     description='GMI scattering index',           fnct=trim(FNCT_SCATIND), &
                     par=t_chk_par('-17.0',empty_trf,(/11,9,-1/),(/-1,-1,-1/),.TRUE.),      &
                     flagchan='8-13')
      call add_check(instr=(/71/),                 phase=PH_ABC,   priority= 0,             &
                     description='GMI fg dep. (166 GHz)',          fnct=trim(FNCT_FGDEP),   &
                     par=t_chk_par('2.5',empty_trf,(/11,-1,-1/),(/-1,-1,-1/),.FALSE.),      &
                     flagchan='8-13')
      call add_check(instr=(/71/),                 phase=PH_ABC,   priority=-1,             &
                     description='Always cloudy fallback',         fnct=trim(FNCT_ALWAYS))
    end if
    if (do_clddetect_ssmis) then
      call add_check(instr=(/10/),                 phase=PH_ABC,   priority= 0,             &
                     description='SSMI/S polarization difference', fnct=trim(FNCT_POL_DIFF),&
                     par=t_chk_par('0.9',empty_trf,(/16,15,-1/),(/-1,-1,-1/),.TRUE.),       &
                     flagchan='12-18')
      call add_check(instr=(/10/),                 phase=PH_ABC,   priority= 0,             &
                     description='SSMI/S scattering index',        fnct=trim(FNCT_SCATIND), &
                     par=t_chk_par('-6.0',empty_trf,(/8,17,-1/),(/-1,-1,-1/),.FALSE.),      &
                     flagchan='8-11')
      call add_check(instr=(/10/),                 phase=PH_ABC,   priority= 0,             &
                     description='SSMI/S fg dep. (150 GHz)',       fnct=trim(FNCT_FGDEP),   &
                     par=t_chk_par('4.0',empty_trf,(/8,-1,-1/),(/-1,-1,-1/),.FALSE.),       &
                     flagchan='8-11')
      call add_check(instr=(/10/),                 phase=PH_ABC,   priority=-1,             &
                     description='Always cloudy fallback',         fnct=trim(FNCT_ALWAYS))
    end if
    if (do_clddetect_saphir) then
      call add_check(instr=(/34/),                 phase=PH_ABC,   priority= 0,             &
                     description='SAPHIR fg dep. (183 pm 11 GHz)', fnct=trim(FNCT_FGDEP),   &
                     par=t_chk_par('0.5',empty_trf,(/6,-1,-1/),(/-1,-1,-1/),.FALSE.))
      call add_check(instr=(/34/),                 phase=PH_ABC,   priority=-1,             &
                     description='Always cloudy fallback',         fnct=trim(FNCT_ALWAYS))
    end if
    if (do_clddetect_mwhs2) then
      call add_check(instr=(/73/),                 phase=PH_ABC,   priority=0,              &
                     description='MWHS2 scattering index',         fnct=trim(FNCT_SCATIND), &
                     par=t_chk_par('-20.0',empty_trf,(/10,1,-1/),(/-1,-1,-1/),.TRUE.))
      call add_check(instr=(/73/),                 phase=PH_ABC,   priority=0,              &
                     description='MWHS2 fg dep. (150 GHz)',        fnct=trim(FNCT_FGDEP),   &
                     par=t_chk_par('2.0',empty_trf,(/10,-1,-1/),(/-1,-1,-1/),.TRUE.))
      call add_check(instr=(/73/),                 phase=PH_ABC,   priority=-1,             &
                     description='Always cloudy fallback',         fnct=trim(FNCT_ALWAYS))
    end if
    if (do_clddetect_mhs) then
       ! METOP-B
      call add_check(instr=(/15/), satids=(/3/),   phase=PH_ABC,   priority=0,              &
                     description='MHS scattering index',           fnct=trim(FNCT_SCATIND), &
                     par=t_chk_par('-20.0',empty_trf,(/2,1,-1/),(/-1,-1,-1/),.TRUE.))
      call add_check(instr=(/15/), satids=(/3/),   phase=PH_ABC,   priority=0,              &
                     description='MHS fg dep. (157 GHz)',          fnct=trim(FNCT_FGDEP),   &
                     par=t_chk_par('1.0',empty_trf,(/2,-1,-1/),(/-1,-1,-1/),.FALSE.))
       ! METOP-A
      call add_check(instr=(/15/), satids=(/4/),   phase=PH_ABC,   priority=0,              &
                     description='MHS scattering index',           fnct=trim(FNCT_SCATIND), &
                     par=t_chk_par('-2.0',empty_trf,(/2,1,-1/),(/-1,-1,-1/),.TRUE.))
      call add_check(instr=(/15/), satids=(/4/),   phase=PH_ABC,   priority=0,              &
                     description='MHS fg dep. (157 GHz)',          fnct=trim(FNCT_FGDEP),   &
                     par=t_chk_par('3.5',empty_trf,(/2,-1,-1/),(/-1,-1,-1/),.FALSE.))
       ! NOAA-18
      call add_check(instr=(/15/), satids=(/209/), phase=PH_ABC,   priority=0,              &
                     description='MHS scattering index',           fnct=trim(FNCT_SCATIND), &
                     par=t_chk_par('-20.0',empty_trf,(/2,1,-1/),(/-1,-1,-1/),.TRUE.))
      call add_check(instr=(/15/), satids=(/209/), phase=PH_ABC,   priority=0,              &
                     description='MHS fg dep. (157 GHz)',          fnct=trim(FNCT_FGDEP),   &
                     par=t_chk_par('4.0',empty_trf,(/2,-1,-1/),(/-1,-1,-1/),.FALSE.))
       ! NOAA-19
      call add_check(instr=(/15/), satids=(/223/), phase=PH_ABC,   priority=0,              &
                     description='MHS scattering index',           fnct=trim(FNCT_SCATIND), &
                     par=t_chk_par('-20.0',empty_trf,(/2,1,-1/),(/-1,-1,-1/),.TRUE.))
      call add_check(instr=(/15/), satids=(/223/), phase=PH_ABC,   priority=0,              &
                     description='MHS fg dep. (157 GHz)',          fnct=trim(FNCT_FGDEP),   &
                     par=t_chk_par('4.0',empty_trf,(/2,-1,-1/),(/-1,-1,-1/),.FALSE.))
      ! all satellites
      call add_check(instr=(/15/),                 phase=PH_ABC,   priority=-1,             &
                     description='Always cloudy fallback',         fnct=trim(FNCT_ALWAYS))
    end if
    if (do_clddetect_amsr2) then
      call add_check(instr=(/63/),                 phase=PH_ABC,   priority= 0,             &
                     description='AMSR2 polarization difference',  fnct=trim(FNCT_POL_DIFF),&
                     par=t_chk_par('0.95',empty_trf,(/11,12,-1/),(/-1,-1,-1/),.TRUE.))
      call add_check(instr=(/63/),                 phase=PH_ABC,   priority=-1,             &
                     description='Always cloudy fallback',         fnct=trim(FNCT_ALWAYS))
    end if
    if (do_clddetect_mwts3) then
      call add_check(instr=(/132/),                phase=PH_ABC,   priority= 0,             &
                     description='MWTS3 emissivity cloud check' , fnct=trim(FNCT_EMISS),    &
                     flags=ibset(0,CLD_CHK_SEA))
      call add_check(instr=(/132/),                phase=PH_ABC,   priority=-1,             &
                     description='Always cloudy fallback',         fnct=trim(FNCT_ALWAYS))
    end if

    ! Phase "before first guess checks" defaults
    select case(clchk_amsub)
    case(-1)
      description = 'AMSU-B always cloudy'
      amsub_fnct = FNCT_ALWAYS
    case(0)
      ! No cloud check for amsub
      amsub_fnct = ''
    case(1)
      description = 'AMSU-B Buehler(2007) cloud check'
      amsub_fnct = FNCT_AMSUB_BUEHLER
    case(2)
      description = 'AMSU-B ECMWF type I check'
      amsub_fnct = FNCT_AMSUB_EC1
    case(3)
      description = 'AMSU-B chan5-chan3 check'
      amsub_fnct = FNCT_AMSUB_5_3
    case(4)
      description = 'AMSU-B ECMWF type II check'
      amsub_fnct = FNCT_AMSUB_EC2
    case default
      call finish('cloud_detect_init','invalid value for clchk_amsub')
    end select
    if (amsub_fnct /= '') then
      i = 2
      instr(1:i) = (/ 4, 15/)
      if (do_clddetect_atms_amsub) then
        i = i + 1
        instr(i) = 19
      end if
      call add_check(instr=instr(1:i), phase=PH_BFG, priority= 0,          &
                     description=trim(description), fnct=trim(amsub_fnct), &
                     par=t_chk_par('',empty_trf,(/-1,-1,-1/),(/-1,-1,-1/),.FALSE.),     &
                     flagchan='instr=4 1-5 instr=15 1-5 instr=19 17-22')
      call add_check(instr=instr(1:i), phase=PH_BFG, priority= -1, &
                     description='Always cloudy fallback', fnct=trim(FNCT_ALWAYS), &
                     flagchan='instr=4 1-5 instr=15 1-5 instr=19 17-22')
    end if

    ! Phase "after thinning" defaults
    if (do_clddetect_hirs) then !.and. biascor_mode /= BC_NOBC for full backwards compatibility
      call add_check(instr=(/0/), phase=PH_ATH, priority= 0, description='HIRS cloud check', &
                     fnct=trim(FNCT_HIRS))
      ! The "Always cloudy" fallback is not called here, since the hirs cloud check has
      ! it's own channel-dependent fallback
    end if

    !-----------------------------------
    ! 3. Read TOVS_CLOUD_CHECK namelists
    !-----------------------------------
    inml  = 1
    do
      if (dace% lpio) then
        call position_nml ('TOVS_CLOUD_CHECK' ,lrewind=(inml==1) ,status=ierr)
        select case (ierr)
        case (POSITIONED)
          instr        = empty_chk%instr
          satids       = empty_chk%satids
          grids        = empty_chk%grids
          phase        = empty_chk%phase
          c_phase      = ''
          priority     = empty_chk%priority
          description  = empty_chk%description
          chain        = empty_chk%chain
          fnct         = empty_chk%fnct
          flags        = empty_chk%flags
          flagchan     = empty_chk%flagchan
          tf_cloud_ind = empty_chk%tf_cloud_ind
          fg%threshold = -huge(0._wp)
          fg%chan      = -1
          call construct(par_)
#if 1  /* was: defined (__GFORTRAN__) || defined (_CRAYFTN) */
          call read_nml_1
          if (ios /= 0) THEN
            call back_nml(nnml, ios)
            if (ios == 0) call read_nml_2
          end if
#else  /* ??? */
          call read_nml_2
          if (ios /= 0) THEN
            call back_nml(nnml, ios)
            if (ios == 0) call read_nml_1
          end if
#endif
          if (ios /= 0) then
            write(msg, '("Failed to read TOVS_CLOUD_CHECK namelist number ",I2)') inml
            call finish('cloud_detect_init',trim(msg))
          end if

          if (any(tf_cloud_ind > size(TF_CLOUD))) then
            write(msg, '("Failed to read TOVS_CLOUD_CHECK namelist number ",I2,&
                 &" tf_cloud_ind too large: ",10(1x,I2))') inml,tf_cloud_ind,size(TF_CLOUD)
            call finish('cloud_detect_init',trim(msg))
          end if

          ! deprecated: just for backwards compatibility with old fg_cld_check
          if (tolower(fnct) == 'fg_cld_check') then
            fnct = FNCT_FGDEP
            par_%use_rawbt = .false.
          end if

          if (phase <= 0 .and. c_phase /= '') then
            ! determine phase from c_phase
            do i = 1, n_phases
              if (trim(get_c_phase(i)) == trim(c_phase)) phase = i
            end do
          end if
          ! remove double entries from "instr"
          j = 1
          do i = 2, size(instr)
            if (.not.any(instr(1:i-1) == instr(i))) then
              j = j + 1
              if (i /= j) instr(j) = instr(i)
            else
              instr(i) = -1
            end if
          end do

          cld_chk_read = t_cld_check(instr       , &
                                     satids      , &
                                     grids       , &
                                     phase       , &
                                     priority    , &
                                     description , &
                                     chain       , &
                                     fnct        , &
                                     flags       , &
                                     flagchan    , &
                                     tf_cloud_ind, &
                                     par_        )

          call destruct(par_%threshold)
       end select
      end if
      call p_bcast (ierr, dace% pio)
      if (ierr /= POSITIONED) exit

      call p_bcast(cld_chk_read, dace% pio)
      call add_check(new_chk=cld_chk_read)

      inml = inml + 1
    end do

    ! Sort according to phase/instrument/priority
    do i = n_chk - 1, 1, -1
      do j = 1, i
        if (check(j)%phase == check(j+1)%phase) then
          if (check(j)%priority == check(j+1)%priority) then
            l_swap = check(j)%instr(1) > check(j+1)%instr(1)
          else
            l_swap = check(j)%priority < check(j+1)%priority
          end if
        else
          l_swap = check(j)%phase > check(j+1)%phase
        end if
        if (l_swap) then
          cdum       = check(j)
          check(j)   = check(j+1)
          check(j+1) = cdum
        end if
      end do
    end do

    ! Attach a tf_cloud_ind bit to each test
    !   Determine instruments
    n_instr = 0
    do i = 1, n_chk
      do j = 1, ninstr
        ins = check(i)%instr(j)
        if (ins >= 0) then
          if (.not.any(instrs(1:n_instr) == ins)) then
            n_instr = n_instr + 1
            instrs(n_instr) = ins
          end if
        end if
      end do
    end do
    do i = 1, n_instr
      ins = instrs(i)
      lused(:) = .false.
      ! Determine checks with a tf_cloud_ind given by the namelist
      do j = 1, n_chk
        do ii = 1, ninstr
          if (check(j)%instr(ii) == ins) then
            if (check(j)%tf_cloud_ind(ii) > 0) lused(check(j)%tf_cloud_ind(ii)) = .true.
          end if
        end do
      end do
      ibit = 1
      do j = 1, n_chk
        do ii = 1, ninstr
          if (check(j)%instr(ii) == ins) then
            if (check(j)%tf_cloud_ind(ii) <= 0) then
              do while (lused(ibit))
                ibit = ibit + 1
                if (ibit > size(TF_CLOUD)) then
                  lused(:) = .false.
                  ibit = 1
                end if
              end do
              check(j)%tf_cloud_ind(ii) = ibit
              lused(ibit) = .true.
            end if
          end if
        end do
      end do
    end do

    ! print checks
    if (dace% lpio) then
      do i = 1, n_chk
        call print_chk(check(i))
      end do
    end if

    !---------------------------------------
    ! 4. Initialize cloud detection routines
    !---------------------------------------
    if (n_chk > 0) then
      ! McNally-Watts cloud detection
      if (any(check(1:n_chk)%fnct == FNCT_MNW)) call cads_init_ifc

      ! AMSU-B cloud check parameters
      call read_param_amsub (param_file_amsub)
      param_amsub% check           = amsub_fnct
      param_amsub% ec_bnd          = amsub_ec_bnd
      param_amsub% delta_ch20_ch18 = amsub_delta_ch20_ch18
      call p_bcast(param_amsub ,dace% pio)

      ! Scattering test thresholds
      call read_thresh_params(ios)
      if (ios/=0) call finish ('cloud_detect_init','cannot read "'//trim(param_file_scatt)//&
           &'" or "'//trim(param_file_si_amsua)//'"')

      if (do_clddetect_hirs) call read_nml_hirs

    end if

  contains

    subroutine read_nml_1
      type t_par1
        real(wp) :: threshold  = -huge(0._wp)
        integer  :: ch1        = -1
        integer  :: ch2        = -1
        integer  :: ch(3)      = -1
        integer  :: instr(3)   = -1
        logical  :: use_rawbt  = .TRUE.
      end type t_par1
      type(t_par1) :: par
      type(t_par1) :: empty_par
      namelist /TOVS_CLOUD_CHECK/ instr, satids, grids, phase, c_phase, priority, description, &
           chain, fnct, flags, flagchan, tf_cloud_ind, fg, par

      par = empty_par
      read (nnml ,nml=TOVS_CLOUD_CHECK, iostat=ios)
      if (ios /= 0) RETURN

      call construct(par_)
      par_%ch(1)            = par%ch1
      par_%ch(2)            = par%ch2
      where(par%ch   (1:3) /= -1) par_%ch   (1:3) = par%ch   (1:3)
      where(par%instr(1:3) /= -1) par_%instr(1:3) = par%instr(1:3)
      par_%use_rawbt        = par%use_rawbt
      if (par%threshold /= -huge(0._wp)) then
        write(par_%thresh_str,'(e20.13)') par%threshold
      else if (fg%threshold /= -huge(0._wp)) then
        ! Just for backwards compatibilty
        write(par_%thresh_str,'(e20.13)') fg%threshold
        par_%ch(1)              = fg%chan
      end if
    end subroutine read_nml_1

    subroutine read_nml_2
      type t_par2
        character(len=3000) :: threshold  = ''
        integer             :: ch1        = -1
        integer             :: ch2        = -1
        integer             :: ch(3)      = -1
        integer             :: instr(3)   = -1
        logical             :: use_rawbt  = .TRUE.
      end type t_par2
      type(t_par2) :: par
      type(t_par2) :: empty_par
      namelist /TOVS_CLOUD_CHECK/ instr, satids, grids, phase, c_phase, priority, description, &
           chain, fnct, flags, flagchan, tf_cloud_ind, par

      par = empty_par
      read (nnml ,nml=TOVS_CLOUD_CHECK, iostat=ios)
      if (ios /= 0) RETURN

      call construct(par_)
      par_%thresh_str  = trim(par%threshold)
      par_%ch(1)       = par%ch1
      par_%ch(2)       = par%ch2
      where(par%ch   (1:3) /= -1) par_%ch   (1:3) = par%ch   (1:3)
      where(par%instr(1:3) /= -1) par_%instr(1:3) = par%instr(1:3)
      par_%use_rawbt = par%use_rawbt
    end subroutine read_nml_2

    subroutine add_check(instr, satids, grids, phase, priority, description, chain, &
                         fnct, flags, flagchan, tf_cloud_ind, par, new_chk)
      integer,           intent(in), optional :: instr(:)
      integer,           intent(in), optional :: satids(:)
      integer,           intent(in), optional :: grids(:)
      integer,           intent(in), optional :: phase
      integer,           intent(in), optional :: priority
      character(len=*),  intent(in), optional :: description
      character(len=*),  intent(in), optional :: chain
      character(len=*),  intent(in), optional :: fnct
      integer,           intent(in), optional :: flags
      character(len=*),  intent(in), optional :: flagchan
      integer,           intent(in), optional :: tf_cloud_ind(ninstr)
      type(t_chk_par),   intent(in), optional :: par
      type(t_cld_check), intent(in), optional :: new_chk

      type(t_cld_check), pointer  :: check_(:)
      type(t_cld_check)           :: chk_loc
      integer                     :: i, j, n_chk_old, ni, ni_, ns, ns_, ng, ng_
      integer                     :: stat

      if ( present(instr).and.present(phase).and.present(priority).and.&
           present(description).and.present(fnct)) then
        ni = min(count(instr >= 0), size(chk_loc%instr))
        chk_loc%instr(1:ni) = instr(1:ni)
        chk_loc%phase       = phase
        chk_loc%priority    = priority
        chk_loc%description = description
        chk_loc%fnct        = fnct
        call construct(chk_loc%par)
      elseif (present(new_chk)) then
        chk_loc             = new_chk
      else
        call finish('cloud_detect_init','invalid call to add_check')
      end if
      ni = count(chk_loc%instr >= 0)
      if (present(par         )) chk_loc%par          = par
      if (present(flags       )) chk_loc%flags        = flags
      if (present(flagchan    )) chk_loc%flagchan     = flagchan
      if (present(chain       )) chk_loc%chain        = chain
      if (present(tf_cloud_ind)) chk_loc%tf_cloud_ind = tf_cloud_ind
      if (present(satids)) then
         ns = size(satids)
         chk_loc%satids(1:ns) = satids
      end if
      ns = count(chk_loc%satids >= 0)
      if (present(grids)) then
         ng = size(grids)
         chk_loc%grids(1:ng) = grids
      end if
      ng = count(chk_loc%grids >= 0)
      chk => null()

      chk_loop: do i = 1, n_chk
        if (trim(check(i)%description) /= trim(chk_loc%description)) cycle
        ni_ = count(check(i)%instr >= 0)
        if (ni /= ni_) cycle
        do j = 1, ni
          if (.not.any(check(i)%instr(1:ni_) == chk_loc%instr(j))) cycle chk_loop
        end do
        ns_ = count(check(i)%satids >= 0)
        if (ns /= ns_) cycle
        do j = 1, ns
          if (.not.any(check(i)%satids(1:ns_) == chk_loc%satids(j))) cycle chk_loop
        end do
        ng_ = count(check(i)%grids >= 0)
        if (ng /= ng_) cycle
        do j = 1, ng
          if (.not.any(check(i)%grids(1:ng_) == chk_loc%grids(j))) cycle chk_loop
        end do
        chk => check(i)
        exit
      end do chk_loop

      if (.not.associated(chk)) then
        ! Add new check
        n_chk = n_chk + 1
        if (associated(check)) then
          n_chk_old = size(check)
        else
          n_chk_old = 0
        end if
        if (n_chk > n_chk_old) then
          allocate(check_(n_chk))
          if (n_chk_old > 0) then
            ! We do a simple assignment(=) here, although check%threshold is a derived
            ! type that contains pointers. But that is okay here,
            check_(1:n_chk_old) = check(1:n_chk_old)
            deallocate(check)
          end if
          check => check_
        end if
        chk => check(n_chk)
        call construct(chk)
        if (chk_loc%fnct == '') call finish('cloud_detect_init/add_check',&
             'missing function name (namelist entry "fnct") for "'//trim(chk_loc%description)//'"')
        if (.not.any(phases(:) == chk_loc%phase) .and. dace% lpio) &
             write(6,*) '*** WARNING: added cloud check, which is not active (unknown phase)'
      end if

      chk%instr(1:ni) = chk_loc%instr(1:ni)
      if (any(chk_loc%satids       /= empty_chk%satids)       ) chk%satids         = chk_loc%satids
      if (any(chk_loc%grids        /= empty_chk%grids)        ) chk%grids          = chk_loc%grids
      if (chk_loc%phase            /= empty_chk%phase         ) chk%phase          = chk_loc%phase
      if (chk_loc%priority         /= empty_chk%priority      ) chk%priority       = chk_loc%priority
      if (chk_loc%description      /= empty_chk%description   ) chk%description    = chk_loc%description
      if (chk_loc%chain            /= empty_chk%chain         ) chk%chain          = chk_loc%chain
      if (chk_loc%fnct             /= empty_chk%fnct          ) chk%fnct           = chk_loc%fnct
      if (chk_loc%flags            /= empty_chk%flags         ) chk%flags          = chk_loc%flags
      if (chk_loc%flagchan         /= empty_chk%flagchan      ) chk%flagchan       = chk_loc%flagchan
      if (any(chk_loc%tf_cloud_ind /= empty_chk%tf_cloud_ind) ) chk%tf_cloud_ind   = chk_loc%tf_cloud_ind
      if (chk_loc%par%thresh_str   /= empty_chk%par%thresh_str) then
        chk%par%thresh_str = chk_loc%par%thresh_str
        chk%par%threshold  = chk_loc%par%threshold
        call init(chk%par%thresh_str, c_var, chk%par%threshold, &
             status=stat, used_vnames=.true.)
        if (stat /= 0) call finish('cloud_detect_init','failed to interpret threshold "'//&
             trim(cld_chk_read%par%thresh_str)//'" for check "'//trim(chk%description)//'"')
      end if
      where (chk_loc%par%ch   (1:3) /= empty_chk%par%ch   (1:3)) chk%par%ch   (1:3) = chk_loc%par%ch   (1:3)
      where (chk_loc%par%instr(1:3) /= empty_chk%par%instr(1:3)) chk%par%instr(1:3) = chk_loc%par%instr(1:3)
      if (chk_loc%par%use_rawbt  .neqv. empty_chk%par%use_rawbt ) &
           chk%par%use_rawbt  = chk_loc%par%use_rawbt

    end subroutine add_check

  end subroutine cloud_detect_init
  !=========================================================================
  subroutine print_chk(chk)
    type(t_cld_check), intent(in) :: chk
    integer                       :: i
    write(6,*)
    write(6,*) 'TOVS_CLOUD_CHECK'
    write(6,'(3x,"instr               : ",20(1x,I3.3))')    pack(chk%instr,  mask=chk%instr >=0)
    if (any(chk%satids >=0)) &
         write(6,'(3x,"satids              : ",20(1x,I3.3))')    pack(chk%satids, mask=chk%satids>=0)
    if (any(chk%grids >=0)) &
         write(6,'(3x,"grids               : ",20(1x,I3.3))')    pack(chk%grids,  mask=chk%grids>=0)
    write(6,'(3x,"phase               : ",I3,"  (",A,")")') chk%phase, trim(get_c_phase(chk%phase))
    write(6,'(3x,"priority            : ",I3  )')           chk%priority
    write(6,'(3x,"description         : ",A   )')           trim(chk%description)
    if (chk%chain /= '') write(6,'(3x,"chain               : ",A   )')           trim(chk%chain)
    write(6,'(3x,"fnct                : ",A   )')           trim(chk%fnct)
    write(6,'(3x,"flags               : ",I8  )')           chk%flags
    write(6,'(3x,"flagchan            : ",A   )')           trim(chk%flagchan)
    write(6,'(3x,"tf_cloud_ind        : ",10(I2,1x))')      pack(chk%tf_cloud_ind, mask=chk%instr >=0)
    write(6,'(3x," -> bit in tovs_flag: ",10(I2,1x))')      TF_CLOUD(pack(chk%tf_cloud_ind, mask=chk%instr >=0))
    if (chk%fnct == FNCT_POL_DIFF .or. chk%fnct == FNCT_LWP   .or. &
        chk%fnct == FNCT_SCATIND  .or. chk%fnct == FNCT_FGDEP .or. &
        chk%fnct == FNCT_SI_AMSUA) then
      if (chk%fnct == FNCT_SCATIND .and. chk%par%thresh_str == '') then
        write(6,'(3x,"par%threshold       = ",A)')           ''''' (taken from parameter file)'
      else
        write(6,'(3x,"par%threshold       = ",A)')           trim(chk%par%thresh_str)
      end if
      write(6,'(3x,"par%ch(1)           = ",I5)')          chk%par%ch(1)
      write(6,'(3x,"par%ch(2)           = ",I5)')          chk%par%ch(2)
      if (chk%fnct == FNCT_SI_AMSUA) &
           write(6,'(3x,"par%ch(3)           = ",I5)')     chk%par%ch(3)
      do i = 1, size(chk%par%instr)
        if (chk%par%instr(i) >= 0) &
             write(6,'(3x,"par%instr(",I1,")        = ",I5)') i,chk%par%instr(i)
      end do
      write(6,'(3x,"par%use_rawbt       = ",L5)')          chk%par%use_rawbt
    elseif(chk%fnct == FNCT_AMSUB_LWP.or. chk%fnct == FNCT_AMSUB_TPW) then
      write(6,'(3x,"par%threshold       = ",A)')           trim(chk%par%thresh_str)
    elseif(chk%fnct == FNCT_AMSUB_BUEHLER) then
      write(6,'(3x,"par%use_rawbt       = ",L5)')          chk%par%use_rawbt
    end if

  end subroutine print_chk
  !=========================================================================

  character(len=10) function get_c_phase(phase)
    integer, intent(in) :: phase

    select case(phase)
    case(PH_CC)
      get_c_phase = 'cloud_calc'
    case(PH_ABC)
      get_c_phase = 'after_bc'
    case(PH_BFG)
      get_c_phase = 'before_fg'
    case(PH_ATH)
      get_c_phase = 'after_thin'
    case default
      get_c_phase = 'NOT ACTIVE'
    end select

  end function get_c_phase



  elemental subroutine construct_chk_par(x)
    type(t_chk_par), intent(out) :: x
    call construct(x%threshold)
  end subroutine construct_chk_par

  elemental subroutine construct_cld_check(x)
    type(t_cld_check), intent(out) :: x
    call construct(x%par)
  end subroutine construct_cld_check

  subroutine bcast_cld_check(c, source, comm)
    type (t_cld_check)  ,intent(inout) :: c
    integer             ,intent(in)    :: source
    integer   ,optional ,intent(in)    :: comm

    integer :: lcom
    integer :: count = 0

    lcom = dace% comm ;if(present(comm)) lcom = comm

    if (count == 0) count = size(transfer(c,(/' '/)))

    if (dace% pe /= source) call construct (c)
    call p_bcast_derivedtype (c, count, source, lcom)
    call p_bcast(c% par% threshold, source, lcom)

  end subroutine bcast_cld_check


  subroutine calc_tovs_cloud(obs, atm, xi, y)
    use mo_atm_state,     only: t_atm
    use mo_dec_matrix,    only: t_vector
    use mo_radbias_3dv,   only: radbias_bfg
    use mo_t_obs,         only: TSK_Y,po_context,POC_PRE
    use mo_tovs,          only: process_tovs_mult

    type(t_obs_set), intent(inout), target   :: obs    ! observations
    type(t_atm),     intent(in),    optional :: atm    ! atmospheric state
    type(t_vector),  intent(in),    optional :: xi     ! background (interp.space)
    type(t_vector),  intent(inout), optional :: y      ! result     (observ.space)

    character(len=m_oe)      :: str
    type(t_rad_set), pointer :: rs => null()
    integer                  :: i, j
    integer                  :: instr, ni
    logical                  :: mask_rs(n_set)
    logical                  :: mask_instr(n_set, m_instr)
    logical                  :: mask_chk(n_chk)

    mask_rs    = .false.
    mask_instr = .false.
    do i = 1, n_set
      rs => rad_set(i)
      ni = rs%n_instr
      !mask_rs(i) = any(rs%iopts(:)%cloud_mode < 0)
      mask_instr(i,1:ni) = rs%iopts(1:ni)%l_cldlev
      mask_rs   (i     ) = any(mask_instr(i,1:ni))
    end do
    if (.not.any(mask_rs)) return

    if (dace%lpio) print*,'Calculate TOVS clouds'

    po_context = POC_PRE
    call process_tovs_mult(TSK_Y, obs, atm, xi, y, mask_rs=mask_rs)

    call radbias_bfg(obs, y, xi, mask_rs=mask_rs)

    ! Check whether cloud checks are activated for this purpose. If not, then
    ! force all cloud checks for the selected instruments (with mask_rs=T) to
    ! be used here (set phase=PH_CC)
    do i = 1, n_set
      if (.not.mask_rs(i)) cycle
      rs => rad_set(i)
      do instr = 1, rs%n_instr
        if (.not.mask_instr(i,instr)) cycle
        mask_chk = .false.
        do j = 1, n_chk
          chk => check(j)
          if ( any(check(j)%instr == rs%instr(instr))) then
            if (any(check(j)%satids >= 0)) then
              if (.not.any(check(j)%satids == rs%satid)) cycle
            end if
            if (any(check(j)%grids >= 0)) then
              if (.not.any(check(j)%grids == rs%grid)) cycle
            end if
            mask_chk(j) = .true.
          end if
        end do
        if (.not.any(check(1:n_chk)%phase == PH_CC .and. mask_chk(1:n_chk))) then
          if (dace%lpio) print*,'Force cloud check phase '//get_c_phase(PH_CC)//':'
          do j = 1, n_chk
            if (mask_chk(j) .and. .not.check(j)%phase == PH_CC) then
              check(j)%phase = PH_CC
              if (dace%lpio) write(*,'(3x,I2,1x,A)') j,trim(check(j)%description)
            end if
          end do
        end if
      end do
    end do

    call cloud_detect(obs, y, PH_CC, mask_rs=mask_rs)

  end subroutine calc_tovs_cloud



  subroutine cloud_detect_errmsg(ierr, errmsg)
    integer,          intent(in)  :: ierr
    character(len=*), intent(out) :: errmsg

    character(len=300) :: msg = ''
    integer            :: l

    select case(ierr)
    case(:-1)
      write(msg, '("Channel ",I2," missing")') -ierr
    case(ERR_INV_SURF     )
      msg = 'Invalid surface type'
    case(ERR_INV_THR      )
      msg = 'Invalid threshold, get_threshold failed'
    case(ERR_ALL_MISSING  )
      msg = 'All channels missing'
    case(ERR_PHYS         )
      msg = 'physical constraint(s) not fulfilled'
    case(ERR_LAT          )
      msg = 'Invalid latitude'
    case(ERR_CADS         )
      msg = 'CADS failed (call cads_ifc)'
    case(ERR_MW_EMISS     )
      msg = 'mw_emiss failed'
    case(ERR_CF_AMSUB     )
      msg = 'Config file for amsub/mhs Buehler cloud detection not read (param_file_amsub)'
    case(ERR_NO_BUEHLER_PAR)
      msg = 'Did not find parameters for modified Buehler cloud detection.'
    case(ERR_NO_HIRS_PARAM)
      msg = 'Missing HIRS cloud detection parameters (&HIRS_CLOUD_CHECK namelist)'
    case default
      msg = 'Unknown error'
    end select

    l = min(len(errmsg), len_trim(msg))
    errmsg = msg(1:l)

  end subroutine cloud_detect_errmsg

!==============================================================================
#define DERIVED type(t_err)
#undef  MPI_TYPE
#define p_gather_DERIVED p_gather_err
#include "p_gather_derived.incf"
#undef  DERIVED
!==============================================================================


end module mo_cloud
