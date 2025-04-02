!
!+ Routines and types for radiance data processing
!
module mo_rad
!
! Description:
!   Interface from 3D-Var to radiance data processing
!   to be exchanged with COSMO
!
! Current Code Owner: DWD, Robin Faulwetter
!    phone: +49 69 8062 2746
!    fax:   +49 69 8062 3721
!    email:  robin.faulwetter@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_22        2013-02-13 Robin Faulwetter
!  gather routines and types for radiance data processing
! V1_23        2013-03-26 Robin Faulwetter
!  Implemented processing of CrIS data
! V1_26        2013/06/27 Robin Faulwetter
!  Generalize the 'chan'-entry in TOVS_OBS_CHAN namelists:
!    Now, a string instead of an integer might be used as flag-value.
!    The string might be: a logical array (e.g. 'T F T ...'),
!                      or a sequence of bit numbers (e.g. '0 5 10 ...'),
!                      or a sequence of bit names (e.g. 'sea clear',
!                                                       'sea,land,passive').
!    Bit ranges like '0-10', 'sea-passive' can be used for the flag value.
!  Introduced a check on the influence of the surface onto radiances.
!  Introduced USE_MWSURF bit. Corrected the usage of other USE_* bits.
! V1_27        2013-11-08 Robin Faulwetter
!  Rewrite of radiance flagging, Implemented FOV-dependent obserrors
! V1_28        2014/02/26 Andreas Rhodin
!  changes for 1dvar-mode, IR emissivity, cloud-top-height/fraction
! V1_29        2014/04/02 Andreas Rhodin
!  reduce default memory requirement for rttov_mult_prof on scalar machines
! V1_31        2014-08-21 Robin Faulwetter
!  Unify mo_rad with COSMO. Improve mo_rttov_ifc. New write_rttov_prof routine
! V1_35        2014-11-07 Robin Faulwetter
!  adaptions for COSMO, CrIS, chinese satellites
! V1_42        2015-06-08 Robin Faulwetter
!  Unify modules with COSMO; pass op_na flag from RTTOV to 3dvar
! V1_43        2015-08-19 Robin Faulwetter
!  features required for high peaking radiances over land/clouds
! V1_45        2015-12-15 Harald Anlauf
!  write_rtovp: replace bubblesort by quicksort
! V1_47        2016-06-06 Robin Faulwetter
!  Many improvements for radiances.
!  specify NF90_NETCDF4 when creating monRTOVP.nc file, option to write B_ii
! V1_48        2016-10-06 Robin Faulwetter
!  Implemented RTTOV12
! V1_49        2016-10-25 Robin Faulwetter
!  Option to use mapped amsub 89ghz channel instead of amsua 89ghz channel
!  for amsua cloud detection
! V1_50        2017-01-09 Robin Faulwetter
!  Added features to rttov12
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
!==============================================================================

!-------------
! Modules used
!-------------

  use mo_kind,         only: wp, sp, dp, i8     ! precision kind parameters

  use mo_satid,        only: satname,                         &
                             sat_name2id => satid,            &
                             satid_bufr2rttov

  use mo_mpi_dace,     only: p_bcast,                         &
                             dace

  use mo_exception,    only: finish

  use mo_dace_string,  only: uppercase => toupper,            &
                             tolower,                         &
                             split2

  use mo_instrid,      only: instr_rttov,                     &
                             mw_instr,                        &
                             ir_instr,                        &
                             vis_instr
#ifdef HRIT
  use hrit,            only: hrit_open, hrit_get_data,        &
                             hrit_close, hrit_silent,         &
                             hrit_talk, hrit_debug
#endif
  use iso_fortran_env, only: stderr => error_unit,            &
                             stdout => output_unit
  use netcdf,          only: nf90_strerror,                   &
                             nf90_open,                       &
                             nf90_create,                     &
                             nf90_close,                      &
                             nf90_def_dim,                    &
                             nf90_def_var,                    &
                             nf90_put_var,                    &
                             nf90_put_att,                    &
                             nf90_enddef,                     &
                             nf90_inq_dimid,                  &
                             nf90_inquire_dimension,          &
                             nf90_inq_varid,                  &
                             nf90_inquire_variable,           &
                             nf90_inq_attname,                &
                             nf90_get_att,                    &
                             nf90_get_var,                    &
!                            NF90_64BIT_OFFSET,               &
                             NF90_NETCDF4,                    &
                             NF_CLOBBER => nf90_clobber,      &
                             NF_NOERR   => nf90_noerr,        &
                             NF_NOWRITE => nf90_nowrite,      &
                             NF_INT     => nf90_INT,          &
                             NF_CHAR    => nf90_CHAR,         &
                             NF_FLOAT   => nf90_FLOAT,        &
                             NF_DOUBLE  => nf90_DOUBLE,       &
                             NF_FILL_FLOAT => nf90_fill_float,&
                             NF_FILL_INT => nf90_fill_int,    &
                             nf90_fill_float,                 &
                             nf90_fill_int,                   &
                             NF_GLOBAL  => nf90_global
  use mo_instrid,      only: hss_instr

  use mo_algorithms,   only: index              ! Derive permutation for sorting

  use mo_range_fparse, only: t_range_fparse,                  &
                             p_bcast,                         &
                             destruct,                        &
                             assignment(=)

  implicit none


  !----------------
  ! Public entities
  !----------------
  private
  !-------------------------
  ! derived type definitions
  !-------------------------
  public :: t_radv            ! type to hold radiance data
  public :: t_rad_set         ! type to hold radiance meta data
  public :: t_rad_gopt        ! type to hold common options for all instruments
  public :: t_rad_iopt        ! type to hold instrument specific options
  public :: t_emis_opt        ! type to hold emissivity retrieval options
  public :: t_tskin_opt       ! type to hold tskin retrieval options
  !----------
  ! constants
  !----------
  public :: m_instr           ! max. number of instruments in set
  public :: m_chan            ! max. number of channels in namelist
  public :: m_rad_set         ! max. number of radiance sat. datasets
  public :: m_bd              ! max. number of bands
  public :: m_oe              ! max. obserr string length
  public :: USE_SEA           ! usage flag  0    (1) Sea
  public :: USE_SEAICE        !             1    (2) SeaIce
  public :: USE_LAND          !             2    (4) Land
  public :: USE_HIGHLAND      !             3    (8) Highland
  public :: USE_MISMATCH      !             4   (16) Surface Mismatch
  public :: USE_CLEAR         !             5   (32) Clear
  public :: USE_QCNOTUSE      !             6   (64) Do not use (e.g. for QC ...)
  public :: USE_CLOUD2        !             7  (128) CloudFlag2
  public :: USE_CLOUD3        !             8  (256) CloudFlag3
  public :: USE_PASS_EKF      !             9  (512) Passive in EnKF, not in psas
  public :: USE_PASSIVE       !            10 (1024) Passive in psas and EnKF
  public :: USE_CLOUDDET      !            11 (2048) use channel for IR cloud detection
  public :: USE_BCOR          !            12 (4096) force to keep the channel, since it is required in biascor
  public :: USE_MWSURF        !            13
  public :: USE_MWICE         !            19
  public :: USE_SURFINFL      !            14
  public :: USE_MINTSURF      !            15
  public :: USE_L2C_CLD       !            16
  public :: USE_L2C_SURF      !            17
  public :: USE_BLK_PP        !            18
  public :: WRM_OPEN          ! parameter for write_rtovp
  public :: WRM_WRITE         ! parameter for write_rtovp
  public :: WRM_CLOSE         ! parameter for write_rtovp
  public :: f_ins             ! temporary. 'level'=chan+f_ins*instrument
  public :: rinvalid          ! value for real invalid value
  public :: default_flags     ! default flags (used if NO tovs_obs_chan nml is available)
  public :: n_styp
  public :: ITYP_MW
  public :: ITYP_IR
  public :: ITYP_VIS
  public :: OPTV_L2C
  public :: OPTV_NWC_FLG
  public :: OPTV_ORB_PH
  public :: OPTV_INS_TMP
  public :: OPTV_CLD_FRC
  public :: OPTV_CLD_FLG

  !----------
  ! variables
  !----------
  public :: n_set             ! number of radiance sat. datasets used
  public :: rad_set           ! radiances dataset descriptions
  ! public :: ltest             ! auxiliary debug/test switch
  !-----------------------
  ! routines and operators
  !-----------------------
  public :: construct         ! set components of type t_rad_set
  public :: destruct          ! deallocate t_radv/t_rad_set content
  public :: assignment (=)    ! t_rad_set = t_rad_set
  public :: chan_indx         ! get channel index in sat/instr/chan set
  public :: set_indx          ! get sat/instr/chan set index
  public :: instr_chan        ! get instrument from channel index
  public :: read_satpp_feedbk ! read observation data files from SAT_PP
  public :: print_rad_set     ! print content of t_rad_set
  public :: print_radv        ! print content of t_radv
  public :: link_rad          ! associate t_radv structure with t_rad_set
  public :: reduced_rad_set   !
  public :: read_tovs_obs_chan_nml ! read namelist /TOVS_OBS_CHAN/
  public :: sat_sun_azimuth_angle  ! derive angle between satellite and sun
  public :: get_satangles     ! calculate satellite zenith and azimuth angle
  public :: thin_superob_radv ! superobbing routine
  public :: avg_angle         ! Averaging of angles
  public :: warning
  public :: back_nml           ! go to start of last namelist
  public :: sat_poly2fparse    ! Conversion of old obsolete "sat_poly" string to new fparser input
  public :: p_bcast
  public :: lev2chan
  public :: lev2p
  public :: instr_type
  public :: assign_emis_opt
  public :: assign_tskin_opt
  public :: amsua_chan_id
  public :: amsub_chan_id
  public :: irxsr_chan_id
  public :: seviri_chan_id
  !--------------------------------------------------
  ! parameters may be modified in namelist /TOVS_OBS/
  !--------------------------------------------------
  public :: n_max_calc_k_def  ! maximum number of simultaneous forward  calculations
  public :: n_max_prof_k_def  ! maximum number of simultaneous K-matrix calculations
  public :: n_max_calc_y_def  ! maximum number of profiles in simultaneous forward  calculations
  public :: n_max_prof_y_def  ! maximum number of profiles in simultaneous K-matrix calculations
  public :: quality_mode_def  !
  public :: l2c_type_def      !
  public :: l2c_rel_lim_def   !
  public :: l2c_abs_lim_def   !
  public :: l2c_use_rad_def   !
  public :: l2c_max_def       !
  public :: l2c_max_trop_def  !
  public :: l2c_max_midlat_def!
  public :: l2c_max_polar_def !
  public :: surf_infl_mode_def!
  public :: max_surf_infl_def !
  public :: d_stemp_def       !
  public :: d_emiss_def       !
  public :: l_max_surf_hgt_def!
  public :: pp_flags_instr_def!
  public :: pp_flags_def      !
  public :: use_amsub_89ghz_def!
  public :: usd


#ifdef HRIT
  public :: read_hrit
#endif
  public :: write_rtovp
  public :: t_nc_attr
  public :: MP_NETCDF
  public :: MP_PSAS_ANA
  public :: MP_H
  public :: MP_B
  public :: MP_NOSORT
  public :: err_msg

  !-----------------
  ! Other parameters
  !-----------------
  integer , parameter :: f_ins    = 100    ! 'level'=chan+f_ins*instrument
  real(wp), parameter :: rinvalid = -9999._wp
  integer,  parameter :: m_instr  =   4             ! max. number of instruments in set
  integer,  parameter :: m_chan   = 600             ! max. number of channels in namelist
  integer,  parameter :: m_rad_set= 60              ! maximum number of radiance datasets
  integer,  parameter :: m_oe     = 300             ! maximum length of obserr string
  integer,  parameter :: m_bd     =  10             ! max. number of bands.
                                                    ! (not relevant for hyperspec. sounders,
                                                    !  but for band-"misuse" in surfinfl checks)
  integer,  parameter :: m_im_ch  =   2             ! Maximum number of imager channels
  integer,  parameter :: n_styp    = 3              ! Number of different surface types to consider
                                                    ! seperately for the mw emissivity model

  character(len=20) :: name_use_bit(0:31)
  logical           :: obsolete_bit(0:31)
#define DEF_BIT(bitname,bit,name,lobs) integer,parameter::bitname=bit;data name_use_bit(bit)/name/;data obsolete_bit(bit)/lobs/

  DEF_BIT(USE_SEA     ,  0, 'SEA'     , .false.) !    (1)    Sea
  DEF_BIT(USE_SEAICE  ,  1, 'SEAICE'  , .false.) !    (2)    SeaIce
  DEF_BIT(USE_LAND    ,  2, 'LAND'    , .false.) !    (4)    Land
  DEF_BIT(USE_HIGHLAND,  3, 'HIGHLAND', .false.) !    (8)    Highland
  DEF_BIT(USE_MISMATCH,  4, 'MISMATCH', .false.) !   (16)    Surface Mismatch
  DEF_BIT(USE_CLEAR   ,  5, 'CLEAR'   , .false.) !   (32)    Clear
  DEF_BIT(USE_QCNOTUSE,  6, 'QCNOTUSE', .false.) !   (64)    Do not use in QC, useful for passive chans that are used in QC
  DEF_BIT(USE_CLOUD2  ,  7, 'CLOUD2'  , .true. ) !  (128)    CloudFlag2
  DEF_BIT(USE_CLOUD3  ,  8, 'CLOUD3'  , .true. ) !  (256)    CloudFlag3
  DEF_BIT(USE_PASS_EKF,  9, 'PASS_EKF', .false.) !  (512)    Passive in EnKF, not in psas
  DEF_BIT(USE_PASSIVE , 10, 'PASSIVE' , .false.) ! (1024)    Passive in psas and EnKF
  DEF_BIT(USE_CLOUDDET, 11, 'CLOUDDET', .false.) ! (2048)    use channel for IR cloud detection
  DEF_BIT(USE_BCOR    , 12, 'BCOR'    , .false.) ! (4096)    force to keep the channel, since it is required in biascor
  DEF_BIT(USE_MWSURF  , 13, 'MWSURF'  , .false.) ! (8192)    use also if surftype derived from mw is not selected (use_a*_surftype)
  DEF_BIT(USE_SURFINFL, 14, 'SURFINFL', .false.) !(16384)    use if the influence of the surface is lower than a threshold
  DEF_BIT(USE_MINTSURF, 15, 'MINTSURF', .false.) !(32768)    use also if tsurf is smaller than min_tsurf
  DEF_BIT(USE_L2C_CLD , 16, 'L2C_CLD' , .false.) !(65536)    do not flag for clouds if l2c is above threshold
  DEF_BIT(USE_L2C_SURF, 17, 'L2C_SURF', .false.) !(131072)   do not flag for surface disturbances if l2c is above threshold
  DEF_BIT(USE_BLK_PP  , 18, 'BLK_PP'  , .false.) !(262144)   do not flag for surface disturbances if l2c is above threshold
  DEF_BIT(USE_MWICE   , 19, 'MWICE'   , .false.) !           use also if surftype derived from mw is ice (use_a*_ice+fr_ice,ts_bg)

  data name_use_bit(20:31) / '','','','','','','','','','','','' /

  ! default flags, that are used (by setup_rad_set), if NO tovs_obs_chan namelist is
  ! available (e.g. mec).
  integer, parameter :: default_flags = 2**USE_SEA      + 2**USE_SEAICE + 2**USE_LAND + &
                                        2**USE_HIGHLAND + 2**USE_CLEAR


  !--------------
  ! surface flags
  !--------------

  real(kind=dp), parameter :: pi                 = 3.1415926535897931_dp
  real(kind=dp), parameter :: deg2rad            = pi/180._dp
  real(kind=dp), parameter :: rad2deg            = 180._dp/pi

  ! Codes for different instrument/channel types
  integer, parameter :: ITYP_MW          = 1
  integer, parameter :: ITYP_IR          = 2
  integer, parameter :: ITYP_VIS         = 3

  ! Bits for optional variables (set in t_rad_gopt%opt_vars)
  ! Please, make sure that this is consistent with c_optv in mo_fdbk_tables.f90
  integer, parameter :: OPTV_L2C         = 0
  integer, parameter :: OPTV_NWC_FLG     = 1
  integer, parameter :: OPTV_ORB_PH      = 2
  integer, parameter :: OPTV_INS_TMP     = 3
  integer, parameter :: OPTV_CLD_FRC     = 4
  integer, parameter :: OPTV_CLD_FLG     = 5

  !-----------------------------
  ! t_rad_gopt - general options
  !-----------------------------
  type t_rad_gopt
    ! --- RTTOV call ---
    integer           :: n_max_calc_y        =  -1          ! maximum number of simultaneous forward  calculations
    integer           :: n_max_calc_k        =  -1          ! maximum number of simultaneous K-matrix calculations
    integer           :: n_max_prof_y        =  -1          ! maximum number of profiles in simultaneous forward  calc.
    integer           :: n_max_prof_k        =  -1          ! maximum number of profiles in simultaneous K-matrix calc.
                                                            ! NOTE: n_max_prof affects amount of memory required by 3dvar
                                                            !       n_max_calc affects amount of memory required by RTTOV
    ! --- Thinning ---
    integer           :: quality_mode        =  -1          ! Option for calculation of artificial quality
                                                            ! (useful for thinning)
    ! --- PP_FLAGS from sat_pp ---
    integer           :: pp_flags_instr      = 0            ! instrument index for reading PP_FLAGS
    integer           :: pp_flags            = 0            ! PP_FLAGS to be used for blacklisting
    ! --- cloud check for AMSUA
    logical           :: use_amsub_89ghz     = .false.      ! use 89GHz AMSU-B channel instead of AMSU-A
                                                            ! in AMSUA-A cloud check
    ! --- surface influence test ---
    logical           :: l_max_surf_hgt      = .false.      ! Read MAX_SURFACE_HEIGHT instead of SURFACE_HEIGHT
    ! --- superobbing
    integer           :: thin_superob_box_lines   = 0
    integer           :: thin_superob_box_fovs    = 0
    real(wp)          :: max_size_sobox      = 50._wp       ! maximum distance between pixels in a superob box [km] (security check)
    integer           :: max_timediff_sobox  = 15           ! maximum time diff. between pixels in a superob box [min] (security check)
    integer           :: thin_superob_mode   = 1            ! 1 thin, 2 superob
    ! --- provisional until all-sky solution with channel groups existst
    integer           :: lev_mode            = 0
    ! --- optional variables (read and write)
    integer           :: opt_vars            = 0
    ! --- use colocated imager info (relevant for hyperspectral cloud det.)
    integer           :: use_imager          = 0
    integer           :: im_atlas_id         = -1
    integer           :: im_tskin_opt(3)     = (/-1,-1,0/)  ! (instr,chan,styp-bits) for Tskin retrieval
    real(wp)          :: im_tskin_cfrac      = 0._wp
  end type t_rad_gopt


  !-----------------------------------------
  ! t_rad_iopt - instrument specific options
  !-----------------------------------------
  type t_rad_iopt
    ! --- level2channel assignment and tests ---
    integer           :: l2c_type            =  -1          ! Option to calculate assigned level
    real(wp)          :: l2c_rel_lim         =  huge(1._wp) ! relative limit used in level to channel assignment
    real(wp)          :: l2c_abs_lim         =  huge(1._wp) ! absolute limit used in level to channel assignment
    logical           :: l2c_use_rad         =  .true.      ! Use radiances or BTs in level to channel assignment
    real(wp)          :: l2c_max             = -1._wp       ! maximum allowed level in l2c check (default)
    real(wp)          :: l2c_max_trop        = -1._wp       ! maximum allowed level in l2c check tropics
    real(wp)          :: l2c_max_midlat      = -1._wp       ! maximum allowed level in l2c check midlatitudes
    real(wp)          :: l2c_max_polar       = -1._wp       ! maximum allowed level in l2c check polar regions
    ! --- surface influence test ---
    integer           :: surf_infl_mode(m_bd)= 0            ! surface influence test mode
    real(wp)          :: max_surf_infl(m_bd) = huge(1._wp)  ! Maximum influence by surface
    real(wp)          :: d_stemp(m_bd)       = 0._wp        ! delta(t_skin)
    real(wp)          :: d_emiss(m_bd)       = 0._wp        ! delta(emissivity)

    integer           :: cloud_mode          = 0            ! provisional until all-sky solution with channel groups exists
    logical           :: do_lambertian       = .false.
    real(wp)          :: specularity(0:n_styp-1) = 1._wp
    real(wp)          :: specularity_snow    = 1._wp
    integer           :: rad_out             = 0            ! see OUT_* flags in mo_rtifc
    integer           :: use_o3              = 0
    integer           :: use_co2             = 0
    integer           :: rt_nlevs            = 0            ! number of RTTOV levels
    integer           :: rt_mw_emis_mod      = -1
    ! --- (pre-)calculate and store cldlev (with McNally-Watts)
    ! depends on cloud_mode and obserr-model
    logical           :: l_cldlev            = .false.
  end type t_rad_iopt

  type t_emis_opt
    integer           :: id                  = -1
    integer           :: nml                 = -1
    integer           :: inst                = -1           ! RTTOV instrument ID
    integer           :: prio                =  0           ! Priority (0...), higher number are used as fallback
    integer           :: styp(n_styp)        = -1           ! surface type (0=land,1=sea,2=sea-ice,-1=all)
    integer           :: ns                  =  0
    character(len=120):: descr               =  ''
    integer           :: mode                =  -1
    integer           :: atls                =  -1
    integer           :: n_chan              =  0
    integer, pointer  :: chan(:)             => NULL()
    integer, pointer  :: cdyn(:)             => NULL()
    logical           :: angcorr             = .true.
    real(wp)          :: max_dst             = 0._wp
  end type t_emis_opt

  type t_bc_aux
    integer           :: n_tr                = 0             ! number of transmission factors
    real(wp), pointer :: p_tr(:,:)           => null()       ! pressure values for layers
    integer,  pointer :: ib_ch(:)            => null()       ! bitmask for each channel
    integer,  pointer :: type(:)             => null()       ! type of predictor
  end type t_bc_aux

  type t_tskin_opt
    integer           :: id                  = -1
    integer           :: nml                 = -1
    integer           :: inst                = -1           ! RTTOV instrument ID
    integer           :: prio                =  0           ! Priority (0...), higher number are used as fallback
    integer           :: styp(n_styp)        = -1           ! surface type (0=land,1=sea,2=sea-ice,-1=all)
    integer           :: ns                  =  0
    character(len=120):: descr               =  ''
    integer           :: mode                =  -1
    integer           :: n_cdyn              =  0
    integer, pointer  :: cdyn(:)             => NULL()
  end type t_tskin_opt

  !-----------------------------------
  ! t_rad_gopt and t_rad_iopt defaults
  !-----------------------------------
  ! default memory requirements for rttov_mult_prof=T
#ifdef __NEC__
 ! on vector machine with lots of memory
 integer              :: n_max_calc_k_def   =   30000 ! maximum number of simultaneous forward  calculations
 integer              :: n_max_prof_k_def   =   10000 ! maximum number of simultaneous K-matrix calculations
 integer              :: n_max_calc_y_def   = 4000000 ! maximum number of profiles in simultaneous forward  calculations
 integer              :: n_max_prof_y_def   =  300000 ! maximum number of profiles in simultaneous K-matrix calculations
                                                      ! NOTE: n_max_prof affects amount of memory for 3dvar-specific arrays
                                                      !       n_max_calc affects amount of memory for RTTOV
#else
  ! no vector machine
  integer              :: n_max_calc_k_def   =  1000   ! maximum number of simultaneous forward  calculations
  integer              :: n_max_prof_k_def   =  3000   ! maximum number of simultaneous K-matrix calculations
  integer              :: n_max_calc_y_def   = 40000   ! maximum number of profiles in simultaneous forward  calculations
  integer              :: n_max_prof_y_def   =  3000   ! maximum number of profiles in simultaneous K-matrix calculations
#endif
  ! --- Thinning ---
  integer              :: quality_mode_def   = 0       ! Option for calculation of artificial quality (useful for thinning)
  ! --- level2channel assignment and tests ---
  integer              :: l2c_type_def       = 2       ! Option to calculate assigned level
  real(wp)             :: l2c_rel_lim_def    = 0.01_wp ! relative limit used in level to channel assignment
  real(wp)             :: l2c_abs_lim_def    = -1._wp  ! absolute limit used in level to channel assignment
  logical              :: l2c_use_rad_def    = .true.  ! Use radiances or BTs in level to channel assignment
  real(wp)             :: l2c_max_def        = -1._wp  ! maximum allowed level in l2c check (default)
  real(wp)             :: l2c_max_trop_def   = -1._wp  ! maximum allowed level in l2c check
  real(wp)             :: l2c_max_midlat_def = -1._wp  ! maximum allowed level in l2c check
  real(wp)             :: l2c_max_polar_def  = -1._wp  ! maximum allowed level in l2c check
  ! --- surface influence test ---
  integer              :: surf_infl_mode_def = 0       ! surface influence test mode
  real(wp)             :: max_surf_infl_def  = 0.1_wp  ! Maximum influence by surface
  real(wp)             :: d_stemp_def        = 3._wp   ! delta(t_skin)
  real(wp)             :: d_emiss_def        = 0.1_wp  ! delta(emissivity)
  logical              :: l_max_surf_hgt_def = .false. ! Read MAX_SURFACE_HEIGHT instead of SURFACE_HEIGHT
  ! --- PP_FLAGS from sat_pp ---
  integer              :: pp_flags_instr_def  = 0      ! instrument index for reading PP_FLAGS
  integer              :: pp_flags_def        = 0      ! PP_FLAGS to be used for blacklisting
  ! --- cloud check for AMSUA
  logical              :: use_amsub_89ghz_def = .false.! Use radiances or BTs in level to channel assignment

  !-----------------------------------------------------
  ! flag monitor_prof, mon_ana_prof
  !-----------------------------------------------------
  integer, parameter :: MP_NETCDF   =  2 ! write NetCDF file
  integer, parameter :: MP_PSAS_ANA =  4 ! write info on PSAS analysis
  integer, parameter :: MP_H        =  8 ! write H matrices
  integer, parameter :: MP_B        = 16 ! write B in interpolation space
  integer, parameter :: MP_NOSORT   = 32 ! skip ordering in *RTOVP.nc
  !-----------------------------------------------------
  ! Flags to be set in "mode" in write_rtovp
  !-----------------------------------------------------
  integer, parameter :: WRM_OPEN      = 1
  integer, parameter :: WRM_WRITE     = 2
  integer, parameter :: WRM_CLOSE     = 4

  !-----------------------------------------------------------------------
  ! specify the set of instruments on a satellite which have their fov
  ! at the same location and are handled jointly by bias correction
  ! and observation operators.
  !-------------------------------------------------------------------
  type t_rad_set
    integer           :: id                  = -1      ! ID to identify the set fastly
    character(len=4)  :: source              = ''      ! source of t_rad_set, e.g. 'nml'
    integer           :: satid               = -1      ! Satellite ID
    integer           :: rttov_satid         =  0      ! RTTOV platform ID
    integer           :: platform            =  0      ! RTTOV platform ID
    integer           :: n_chan              =  0      ! total number of channels
    integer           :: n_instr             =  0      ! number of instruments in set
    integer           :: n_sens              =  0      ! number of sensors in set
    integer           :: n_fov               =  0      ! number of Fields Of View
    integer           :: n_lev               =  0      ! Number of PSAS levels (interp. space)
                                                       ! (currently only required in write_rtovp, where
                                                       !  n_lev should be known before t_radv%n_lev
                                                       !  is available)
!   integer           :: n_pc_emiss          =  0      ! Number of emissivity PCs
    integer           :: grid                =  0      ! RTTOV-ID of target instrument
    integer           :: grid_wmo            =  0      ! WMO ID of target instrument
    integer           :: flag_instr          = -1      ! sensor selection based on instrument
    integer           :: select_sens         = -1      ! sensor selection
    type(t_rad_gopt)  :: gopts                         ! general options
    type(t_rad_iopt)  :: iopts     (m_instr)           ! instrument specific options
    integer           :: instr     (m_instr) = -1      ! instrument(s) (rttov ID)
    integer           :: instr_wmo (m_instr) = -1      ! instrument(s) (WMO ID)
    integer           :: n_ch_i    (m_instr) =  0      ! number of channels per instr.
    integer           :: o_ch_i    (m_instr) =  0      ! offset of channels for instr.
    integer           :: rttov_indx(m_instr) = -1      ! index for instrument in RTTOV setup
    integer  ,pointer :: sensor_instr(:)     => NULL() ! mapping of sensor index to instrument (rttov_id)
    integer  ,pointer :: chan (:)            => NULL() ! list of channels
    integer  ,pointer :: chidx(:)            => NULL() ! channel indices for RTTOV
    integer  ,pointer :: iflag(:)            => NULL() ! internal flags
    integer  ,pointer :: band (:)            => NULL() ! number of band for each channel
    integer  ,pointer :: flag (:)            => NULL() ! flag from TOVS_OBS_CHAN namelist (see USE_* above)
    !----------------------------------------
    ! IR principal component emissivity model
    !----------------------------------------
    real(wp) ,pointer :: wavenum   (:)       => NULL() ! wave number
    real(wp) ,pointer :: emis_land (:)       => NULL() ! mean emissivity over land
    real(wp) ,pointer :: emis_snow (:)       => NULL() ! mean emissivity over snow
    real(wp) ,pointer :: emis_sice (:)       => NULL() ! mean emissivity over sea ice
    real(wp) ,pointer :: emis_pc (:,:)       => NULL() ! emissivity principal component vectors
    real(wp) ,pointer :: w1        (:)       => NULL() ! spectral interpolation coefficient
    integer  ,pointer :: i1        (:)       => NULL() ! spectral interpolation index
    !-------------------
    ! Emissivity options
    !-------------------
    type(t_emis_opt),     pointer :: emis_opt(:) => NULL() ! Emissivity options
    integer                       :: n_emis_opt          =  0      ! number of emissivity options
    integer,              pointer :: emis_ind(:,:,:)     => NULL() ! indices in emis_opt for each channel, surftyp
                                                       ! and fallback level
    !-------------------
    ! Tskin options
    !-------------------
    type(t_tskin_opt),    pointer :: tskin_opt(:)         => NULL() ! tskin options
    integer                       :: n_tskin_opt          =  0      ! number of tskin options
    integer,              pointer :: tskin_ind(:,:,:)     => NULL() ! indices in tskin_opt for each channel,
                                                                    ! surftyp and fallback level
    real(sp)                      :: max_tskin_cld_cov    = 0.0_sp
    character(len=10)             :: chan_cld_cov         = ''      ! 
    !---------------------------
    ! observation-error handling
    !---------------------------
    real(sp),             pointer :: var  (:)            => NULL() ! obserror variance
    character(len=m_oe),  pointer :: oe_str(:)           => NULL()
    type(t_range_fparse), pointer :: oe_trf(:)           => NULL()
    !-------------------------------
    ! Biascorrection auxiliary stuff
    !-------------------------------
    type(t_bc_aux)                :: bc
    !-------------------------------------
    ! description of colocated imager info
    !-------------------------------------
    integer                       :: n_im_chan       =  0
    integer                       :: n_im_cluster    =  0
    integer                       :: im_rttov_id     = -1
    integer                       :: im_rttov_indx   = -1
    integer                       :: im_chan(m_im_ch)=  0          ! imager channel number
    integer                       :: im_chan_indx(m_im_ch) = 0     ! imager channel RTTOV index
  end type t_rad_set


  !--------------------------------------------------------
  ! Type to keep radiance data (suitable for vectorization)
  !--------------------------------------------------------
  type t_radv
     !---------------
     ! identification
     !---------------
     character(len=128) :: filename        = ''     ! 1DVAR feedback file name
     integer            :: file_id         = -1     ! file id
     integer            :: model_date      = 0      ! Model date  YYYYMMDDHH
     integer            :: n_rec           = 0      ! Number of records
     integer            :: n_lev           = 0      ! Number of levels
     integer            :: n_pc_emiss      = 0      ! Number of emissivity PCs
     type (t_rad_set)   :: i                        ! set of instruments & channels used
     character(len=3)   :: p_unit          = 'hpa'  ! unit of pressure variables
     !---------------------------------------------------------
     ! position in 3dvar derived type, input file or model grid
     !---------------------------------------------------------
     integer   ,pointer :: pe          (:) =>NULL() ! Origin PE
     integer   ,pointer :: ind         (:) =>NULL() ! Index on origin PE
     integer   ,pointer :: i_box       (:) =>NULL() ! box index (3dvar) or zonal gridpoint index (COSMO)
     integer   ,pointer :: j_box       (:) =>NULL() ! box index (3dvar) or meridional gridpoint index (COSMO)
     integer   ,pointer :: i_reprt     (:) =>NULL() ! report index
     integer   ,pointer :: i_body    (:,:) =>NULL() ! body   index
     integer   ,pointer :: obsnum      (:) =>NULL() ! obs. number from input file
     integer            :: ideb(5)         = -1     ! spots with debug output
     !--------------------------------------
     ! description of array spec (see below)
     !--------------------------------------
     integer            :: n_trg           =  0     ! number of tracegases
     integer            :: i_o3            = -1     ! index of ozone profile
     integer            :: i_co2           = -1     ! index of co2 profile
     ! ... further tracegases
     !---------
     ! location
     !---------
     integer   ,pointer :: date     (:)    =>NULL() ! Date of measurement YYYYMMDD
     integer   ,pointer :: time     (:)    =>NULL() ! Time of measurement hhmmss
     integer   ,pointer :: date_d   (:)    =>NULL() ! Decoding date of measurement YYYYMMDD
     integer   ,pointer :: time_d   (:)    =>NULL() ! Decoding time of measurement hhmmss
     integer   ,pointer :: ntstep   (:)    =>NULL() ! model time step
     real(wp)  ,pointer :: dlat     (:)    =>NULL() ! latitude
     real(wp)  ,pointer :: dlon     (:)    =>NULL() ! longitude
     integer   ,pointer :: fov      (:)    =>NULL() ! Measurement field of view
     integer   ,pointer :: scanl    (:)    =>NULL() ! Scanline number
     real(wp)  ,pointer :: stzen    (:)    =>NULL() ! Satellite zenith angle
     real(wp)  ,pointer :: stazi    (:)    =>NULL() ! Satellite azimuth angle
     real(wp)  ,pointer :: sunzen   (:)    =>NULL() ! Sun zenith angle
     real(wp)  ,pointer :: sunazi   (:)    =>NULL() ! Sun azimuth angle
     real(wp)  ,pointer :: orb_phase(:)    =>NULL() ! orbit phase
     real(wp)  ,pointer :: instr_temp(:,:) =>NULL() ! instrument temp.
     real(wp)  ,pointer :: landfr   (:,:)  =>NULL() ! Land fraction
     integer   ,pointer :: stype    (:,:)  =>NULL() ! Surface type
     real(wp)  ,pointer :: shgt     (:,:)  =>NULL() ! Topography
     real(wp)  ,pointer :: specularity(:,:)=>NULL() ! specularity
     !-----------
     ! origin
     !-----------
     integer   ,pointer :: center     (:)  =>NULL() ! Originating center
     integer   ,pointer :: subcenter  (:)  =>NULL() ! Originating subcenter
     !-----------
     ! predictors
     !-----------
     real(wp)  ,pointer :: pred    (:,:,:) =>NULL() ! predictors
     real(wp)  ,pointer :: tr      (:,:,:) =>NULL()
     real(wp)  ,pointer :: plev      (:,:) =>NULL()
     !----------------
     ! quality control
     !----------------
     integer   ,pointer :: mdlsfc      (:) =>NULL() ! model land/ice/snow flags
     integer   ,pointer :: op_na       (:) =>NULL() ! operator not applicable
     logical   ,pointer :: not_rej   (:,:) =>NULL() ! 3DVAR valid data
     integer   ,pointer :: cloudy    (:,:) =>NULL() ! 1:cloudy 2:free 3:invalid
     integer   ,pointer :: state     (:,:) =>NULL() ! 3DVAR status
     integer   ,pointer :: flags     (:,:) =>NULL() ! 3DVAR flags
     logical   ,pointer :: valid     (:,:) =>NULL() ! whether channels are useful
     real(wp)  ,pointer :: sinfl     (:,:) =>NULL() ! surface influence on observation
     real(wp)  ,pointer :: plevel    (:,:) =>NULL() ! assigned height level !RF:NEW
     integer   ,pointer :: pp_flags  (:,:) =>NULL() ! sat_pp flags
     integer   ,pointer :: nwc_flg     (:) =>NULL() ! nwc flag
     integer   ,pointer :: r_state     (:) =>NULL() ! nwc flag
     !---------------------------------------------------------------
     ! quantities in observation space (brightness temperatures)
     ! ECMWF convention for sign of bcor_ (differs from 3dvar):
     !   bt_bcor = bt_obs - bcor_
     !---------------------------------------------------------------
     real(wp)  ,pointer :: bt_fg     (:,:) =>NULL() ! First guess
     real(wp)  ,pointer :: bt_obs    (:,:) =>NULL() ! Measured bright.temp.
     real(wp)  ,pointer :: bt_bcor   (:,:) =>NULL() ! Bias corrected bright.temp.
     real(wp)  ,pointer :: bcor_     (:,:) =>NULL() ! Bias correction, ECMWF sign
     real(wp)  ,pointer :: bt_fg_cs  (:,:) =>NULL() ! First guess clear sky bright.temp.
     real(wp)  ,pointer :: rad_fg    (:,:) =>NULL() ! First guess all sky radiance
     real(wp)  ,pointer :: rad_fg_cs (:,:) =>NULL() ! First guess clear sky radiance
     real(wp)  ,pointer :: emiss     (:,:) =>NULL() ! emissivity
     real(wp)  ,pointer :: emis_pc   (:,:) =>NULL() ! emissivity
     real(wp)  ,pointer :: e_fg      (:,:) =>NULL() ! background error
     real(wp)  ,pointer :: e_obs     (:,:) =>NULL() ! observation error
     !------------------------------------------------------------------
     ! Additional imager information, colocated with original instrument
     !------------------------------------------------------------------
     real(wp)  ,pointer :: im_cl_mean(:,:) =>NULL() ! (imager) cluster mean
     real(wp)  ,pointer :: im_cl_stdv(:,:) =>NULL() ! (imager) cluster stddev
     integer   ,pointer :: im_cl_frac(:,:) =>NULL() ! (imager) cluster fraction
     integer   ,pointer :: im_chan   (:)   =>NULL() ! (imager) channel number
     !---------
     ! profiles
     !---------
     real(wp)  ,pointer :: p         (:,:) =>NULL() ! pressure
     real(wp)  ,pointer :: hh        (:,:) =>NULL() ! height of layers
     real(wp)  ,pointer :: h_fg      (:,:) =>NULL() ! First guess height
     real(wp)  ,pointer :: t_fg      (:,:) =>NULL() ! First guess temperature
     real(wp)  ,pointer :: q_fg      (:,:) =>NULL() ! First guess spec.humidity
     real(wp)  ,pointer :: t2m         (:) =>NULL() ! First guess 2m temp.
     real(wp)  ,pointer :: q2m         (:) =>NULL() ! First guess 2m humidity
     real(wp)  ,pointer :: ps_fg       (:) =>NULL() ! First guess surface press.
     real(wp)  ,pointer :: gp_fg       (:) =>NULL() ! First guess surface geopotential
     real(wp)  ,pointer :: ts_fg       (:) =>NULL() ! First guess surface temp.
     real(wp)  ,pointer :: tsm_fg      (:) =>NULL() ! Model first guess surface temp.
     real(wp)  ,pointer :: u10_fg      (:) =>NULL() ! First guess 10 m zonal wind.
     real(wp)  ,pointer :: v10_fg      (:) =>NULL() ! First guess 10 m meridional wind.
     real(wp)  ,pointer :: v10_abs_fg  (:) =>NULL() ! First guess 10 m abs. wind speed.
     real(wp)  ,pointer :: snw_frc     (:) =>NULL() ! First guess cloud fraction
     real(wp)  ,pointer :: cld_top     (:) =>NULL() ! First guess cloud top height
     real(wp)  ,pointer :: cfraction   (:) =>NULL() ! First guess simple cloud fraction
     real(wp)  ,pointer :: clwde     (:,:) =>NULL() ! First guess water eff. diameter profile
     real(wp)  ,pointer :: icede     (:,:) =>NULL() ! First guess ice eff. diameter profile
     real(wp)  ,pointer :: cfrac     (:,:) =>NULL() ! First guess cloud fraction profile
     real(wp)  ,pointer :: cld_fg  (:,:,:) =>NULL() ! First guess cloud water content
     real(wp)  ,pointer :: cld_frc   (:,:) =>NULL() ! Cloud fraction from 1D-Var
     real(wp)  ,pointer :: clw       (:,:) =>NULL() ! First guess cloud liquid water
     real(wp)  ,pointer :: trg     (:,:,:) =>NULL() ! First guess of tracegases
     !----------------
     ! sens. functions
     !----------------
     real(wp)  ,pointer :: H_t     (:,:,:) =>NULL() ! H temperature (chan,lev,prof)
     real(wp)  ,pointer :: H_q     (:,:,:) =>NULL() ! H humidity    (chan,lev,prof)
     real(wp)  ,pointer :: H_ts      (:,:) =>NULL() ! H temperature (chan,prof)
     real(wp)  ,pointer :: H_ps      (:,:) =>NULL() ! H humidity    (chan,prof)
     real(wp)  ,pointer :: H_em_pc (:,:,:) =>NULL() ! H emissivity  (chan,pc,prof)
     real(wp)  ,pointer :: K_t     (:,:,:) =>NULL() ! K temperature (chan,lev,prof)
     real(wp)  ,pointer :: K_q     (:,:,:) =>NULL() ! K humidity    (chan,lev,prof)
     real(wp)  ,pointer :: K_ts      (:,:) =>NULL() ! K temperature (chan,prof)
     real(wp)  ,pointer :: K_ps      (:,:) =>NULL() ! K humidity    (chan,prof)
     !-----------------
     ! background error
     !-----------------
     real(wp)  ,pointer :: t_eb      (:,:) =>NULL() ! Bg. error temperature
     real(wp)  ,pointer :: q_eb      (:,:) =>NULL() ! Bg. error spec.humidity
     real(wp)  ,pointer :: ps_eb       (:) =>NULL() ! Bg. error surface press.
     real(wp)  ,pointer :: ts_eb       (:) =>NULL() ! Bg. error surface temp.
     real(wp)  ,pointer :: cld_top_eb  (:) =>NULL() ! Bg. error cloud top height
     real(wp)  ,pointer :: cld_frc_eb(:,:) =>NULL() ! Bg. error cloud fraction
     real(wp)  ,pointer :: snw_frc_eb  (:) =>NULL() ! Bg. error snow fraction
     real(wp)  ,pointer :: emis_pc_eb(:,:) =>NULL() ! Bg. error emissivity PC
     !-----------------------------
     ! background error covariances
     !-----------------------------
      real(wp)  ,pointer :: B_tt   (:,:,:) =>NULL() ! B temp.-temp. (lev,lev,prof)
      real(wp)  ,pointer :: B_qq   (:,:,:) =>NULL() ! B humi.-humi. (lev,lev,prof)
      real(wp)  ,pointer :: B_tq   (:,:,:) =>NULL() ! B temp.-humi. (lev,lev,prof)
     !-----------------------------
     ! error covariances in obs space
     !-----------------------------
      real(wp)  ,pointer :: R      (:,:,:) =>NULL() ! R     (chan,chan,prof)
      real(wp)  ,pointer :: HBHR   (:,:,:) =>NULL() ! HBH+R (chan,chan,prof)
  end type t_radv

  !----------------------------
  ! attributes for *RTOVP files
  !----------------------------
  type t_nc_attr
    character(len=80) :: name = ''
    character(len=80) :: c    = ''
    integer           :: i    = NF_FILL_INT
  end type t_nc_attr

  !--------------------------------------------------
  ! Mapping between channels of different instruments
  !--------------------------------------------------
  integer, save, target :: amsua_amsua_chans(15) = (/ 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15/)
  integer, save, target :: mhs_amsua_chans  (15) = (/-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 1/)
  integer, save, target :: atms_amsua_chans (15) = (/ 1, 2, 3, 5, 6, 7, 8, 9,10,11,12,13,14,15,16/)
  integer, save, target :: mwts3_amsua_chans(15) = (/ 1, 2, 3, 5, 7, 9,10,11,12,13,14,15,16,17,-1/)
  integer, save, target :: mwhs2_amsua_chans(15) = (/-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 1/)
  integer, save, target :: amsub_amsub_chans (5) = (/ 1, 2, 3, 4, 5/)
  integer, save, target :: atms_amsub_chans  (5) = (/16,17,22,20,18/)
  integer, save, target :: saphir_amsub_chans(5) = (/-1,-1, 2, 3, 5/)
  integer, save, target :: mwhs2_amsub_chans (5) = (/ 1,10,11,13,15/)
  integer, save, target :: ssmis_amsub_chans (5) = (/17, 8,11,10, 9/)
  integer, save, target :: gmi_amsub_chans   (5) = (/ 8,10,-1,12,13/)
  integer, save, target :: seviri_irxsr_chans(7) = (/ 1, 2, 3, 4, 5, 7, 8/) ! 3.9,6.2 ,7.3,8.7,9.7,12.0,13.4
  integer, save, target :: mviri_irxsr_chans (7) = (/-1, 1,-1,-1,-1,-1,-1/) !
  integer, save, target :: abi_irxsr_chans   (7) = (/ 1, 2, 4, 5, 6, 9,10/) ! 3.9,6.15,7.3,8.4,9.6,12.3,13.3
  integer, save, target :: goesim_irxsr_chans(7) = (/ 1, 2,-1,-1,-1,-1, 4/) ! 3.9,6.15,7.3,8.4,9.6,12.3,13.3
  integer, save, target :: fci_irxsr_chans   (7) = (/ 1, 2, 3, 4, 5, 7, 8/)
  integer, save, target :: sev_seviri_chans (11) = (/ 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11/)
  integer, save, target :: fci_seviri_chans (11) = (/ 3, 4, 7, 9,10,11,12,13,14,15,16/)
  
  !=== define module variables:
  type (t_rad_set), save, target :: rad_set(m_rad_set)     ! Radiances dataset descriptions
  integer,          save         :: n_set        = 0       ! Number of useful radiance sat. datasets
  character(len=3000)            :: msg          = ''
  character(len=300), save       :: err_msg      = ''
  integer                        :: usd          = 6
  ! logical                        :: ltest        = .false. ! auxiliary switch for new features/debugging

  !-----------
  ! interfaces
  !-----------
  interface construct
    module procedure construct_rad_set
    module procedure construct_rad_set_1dvar
    module procedure construct_radv
    module procedure construct_rad_gopts
    module procedure construct_rad_iopts
    module procedure construct_emis_opt
    module procedure construct_tskin_opt
  end interface construct

  interface destruct
    module procedure destruct_rad_set
    module procedure destruct_radv
    module procedure destruct_emis_opt
    module procedure destruct_tskin_opt
 end interface

  interface assignment (=)
    module procedure assign_rad_set
    module procedure assign_emis_opt
    module procedure assign_tskin_opt
  end interface assignment (=)

  interface p_bcast
    module procedure bcast_rad_set
    module procedure bcast_emis_opt
    module procedure bcast_tskin_opt
  end interface

#if (defined (__GFORTRAN__) && (__GNUC__ >= 10)) || defined (NAGFOR)
  interface
     subroutine p_bcast_derivedtype (buffer, count, source, comm)
       type(*) ,intent(inout)     :: buffer       ! variable to bcast
       integer ,intent(in)        :: count        ! len(byte) of variable
       integer ,intent(in)        :: source       ! source processor index
       integer ,intent(in)        :: comm         ! communicator
     end subroutine p_bcast_derivedtype
  end interface
#endif

contains
!------------------------------------------------------------------------------

  subroutine read_satpp_feedbk(file, r, valid, r_inv, i_inv, lprint, lread, istart, iend, &
                               status, pe)
  !----------------------------------------
  ! read observation data files from SAT_PP
  !----------------------------------------
  character(len=*), intent(in)            :: file    ! File to read
  type(t_radv),     intent(out), target   :: r       ! t_radv to be filled
  logical,          intent(out)           :: valid   ! Whether the file is valid
  real(wp),         intent(in),  optional :: r_inv   ! float/double fillvalue
  integer,          intent(in),  optional :: i_inv   ! integer fillvalue
  logical,          intent(in),  optional :: lprint  ! print or not
  logical,          intent(in),  optional :: lread   ! read the file or infer the number of obs
  integer,          intent(in) , optional :: istart  ! first record to be read
  integer,          intent(in) , optional :: iend    ! last record to be read
  integer,          intent(out), optional :: status  ! exit status
  integer,          intent(in) , optional :: pe

    type t_nc_dim
      character(len=80) :: name
      integer           :: ncid = -1
      integer           :: len  = 0
    end type t_nc_dim

    type(t_nc_dim),   target      :: d_meas  = t_nc_dim('nMeas'             , -1, 0)
    type(t_nc_dim),   target      :: d_instr = t_nc_dim('nInstr'            , -1, 0)
    type(t_nc_dim),   target      :: d_sens  = t_nc_dim('nSens'             , -1, 0)
    type(t_nc_dim),   target      :: d_chan  = t_nc_dim('nChannel'          , -1, 0)
    type(t_nc_dim),   target      :: d_chcl  = t_nc_dim('nImagerChanCluster', -1, 0)
    type(t_nc_dim),   target      :: d_cl    = t_nc_dim('nImagerCluster'    , -1, 0)

    character(len=3), parameter   :: cind = '  '
    character(len=80)             :: hgt_vname = ''
    type(t_rad_set),  pointer     :: s => null()
    logical                       :: l_print
    logical                       :: l_read
    integer,          allocatable :: chan_blkl(:,:)
    integer                       :: instr_sc (m_instr)
    integer                       :: instr_ec (m_instr)
    integer,          pointer     :: n_instr => null()
    integer,          pointer     :: n_chan  => null()
    integer,          pointer     :: n_meas  => null()
    integer                       :: nf_stat, stat, n_cc
    integer                       :: ncid
    integer                       :: i_e
    integer                       :: i_s
    integer                       :: i, j, instr
    integer                       :: satpp_version
    integer                       :: n_dummy
    integer,          allocatable :: start_sensor(:), end_sensor(:)
    integer,          allocatable :: mapping(:)
    real(wp)                      :: rfill
    integer                       :: ifill
    integer                       :: pe_loc

    integer,          allocatable :: idummy(:,:)

    pe_loc = -1
    if (present(pe)) pe_loc=pe
    err_msg = ' file "'//trim(file)//'"'

    if (present(status)) status = 1
    valid = .false.

    if (present(lprint)) then ; l_print=lprint ; else ; l_print=.true. ; endif
    if (present(lread )) then ; l_read =lread  ; else ; l_read =.true. ; endif

    s       => r% i
    n_meas  => r%n_rec
    n_instr => s%n_instr
    n_chan  => s%n_chan


    if (l_print) then
      print*, repeat ('-',79)
      if (l_read) then
        print*,cind//'read sat_pp file '//trim(file)
      else
        print*,cind//'get number of records in sat_pp file '//trim(file)
      end if
    end if

    ! Open the file
    nf_stat = nf90_open(file, NF_NOWRITE, ncid)
    if (nf_stat /= NF_NOERR) then
      call nf_err('nf90_open')
      err_msg = 'Failed to open'//trim(err_msg)
      return
    endif

    ! get the version of the file
    call detect_satpp_version(satpp_version,stat)
    if (stat /= 0) then
       err_msg = 'Failed to detect version of'//trim(err_msg)
       return
    end if

#define CHECK_ERR(text) if(nf_stat/=nf_noerr)then;err_msg='infer dim '//trim(text)//' in'//trim(err_msg);return;endif
    ! Get dimension lengths
    call get_dim(d_meas)
    CHECK_ERR("d_meas")
    call get_dim(d_instr)
    CHECK_ERR("d_instr")
    call get_dim(d_Sens)
    CHECK_ERR("d_sens")
    call get_dim(d_chan)
    CHECK_ERR("d_chan")
#undef CHECK_ERR

    n_meas   = d_meas% len
    n_instr  = d_instr%len
    n_chan   = d_chan% len
    s%n_sens = d_sens% len

!   write (*,*) cind,'sat_pp version ', satpp_version, ' detected'

    select case (satpp_version)
      case (0)
        n_dummy  = s%n_instr
      case (1)
        n_dummy  = s%n_sens
      case default
        call Error('Unknown satpp_version','INPUT')
        err_msg = 'Unknown satpp_version of'//trim(err_msg)
        return
    end select

    ! Check whether file is empty
    if ((n_meas <= 0) .or. (n_instr <= 0) .or. (n_chan <= 0)) then
      call Error('Empty input file '//trim(file), 'INPUT')
      call nf_close
      if (nf_stat /= NF_NOERR) then
         err_msg = 'Failed to close'//trim(err_msg)
         return
      end if
      if (present(status)) status = 0
      err_msg = 'Zero dimenasion length in'//trim(err_msg)
      return
    end if

    ! Return if reading is not requested
    if (.not.l_read) then
      call nf_close
      if (nf_stat /= NF_NOERR) then
         err_msg = 'Failed to close'//trim(err_msg)
         return
      end if
      valid = .true.
      if (present(status)) status = 0
      return
    end if

    ! Check whether there are too many instruments in the input file
    if (n_instr > m_instr) then
      call Error('Too many instruments in '//trim(file), 'INPUT')
      err_msg = 'Too many instruments in'//trim(err_msg)
      call nf_close
    end if

    ! Define start and end record
    if (present(istart)) then
      i_s = max(1, istart)
    else
      i_s = 1
    end if
    if (present(iend)) then
      i_e = min(iend, n_meas)
    else
      i_e= n_meas
    end if
    if (i_e < i_s) then
      call Error('End record is below start record.', 'INTERNAL')
      err_msg = 'End record is below start record:'//trim(err_msg)
      return
    end if
    n_meas = i_e - i_s + 1

    r% filename = trim(file)

#define CHECK_ERR(text) if(stat/=0)then;err_msg='Failed to read '//trim(text)//' in'//trim(err_msg);return;endif
    ! Set up the description of the dataset (the t_rad_set structure)
    call get_var('SATID'      ,.true., stat, idata0=s% satid)
    CHECK_ERR('SATID')
    call get_var('GRID'       ,.true., stat, idata0=s% grid)
    CHECK_ERR('GRID')
    call get_var('RTTOV_ID'   ,.true., stat, idata1=s% instr(1:n_instr))
    CHECK_ERR('RTTOV_ID')
    call get_var('INSTRUMENT' ,.true., stat, idata1=s% instr_wmo(1:n_instr))
    CHECK_ERR('INSTRUMENT')
    ! Temporary workaround for bug in sat_pp instrument-ID determination
    do i = 1, n_instr
      instr = instr_rttov(s% instr(i), s% satid)
      if (instr > -1) s% instr_wmo(i) = instr
    end do
    ! GRID in the input file is the number of the target instrument. Here it should
    ! be the RTTOV ID:
    if ((s%grid > 0) .and. (s%grid <= n_instr)) then
      s%grid_wmo = s%instr_wmo(s%grid)
      s%grid     = s%instr(s%grid)
    else
      call Error('Invalid grid in file '//trim(file), 'INPUT')
      call nf_close
      err_msg = 'Invalid grid in'//trim(err_msg)
      return
    end if

    call get_var('INSTR_SC'   ,.true., stat, idata1=instr_sc(1:n_instr))
    CHECK_ERR('INSTR_SC')
    call get_var('INSTR_EC'   ,.true., stat, idata1=instr_ec(1:n_instr))
    CHECK_ERR('INSTR_EC')

    allocate(s% chan  (1:n_chan))
    call get_var('CHANNEL'    ,.true., stat, idata1=s% chan(1:n_chan))
    if (stat /= 0) then
       deallocate(s% chan)
       CHECK_ERR('CHANNEL')
    endif

    ! nonobligatory data
    allocate(s% iflag(1:n_chan))
    call get_var('VIS'    ,.false., stat, idata1=s% iflag(1:n_chan))
    IF (stat /= 0) then
      deallocate(s% iflag)
    end IF


    s% o_ch_i(1:n_instr) = instr_sc(1:n_instr) - 1
    s% n_ch_i(1:n_instr) = instr_ec(1:n_instr) - instr_sc(1:n_instr) + 1

    ! Read the data
    ! Macros
#define CHECK_ERR_EXIT(var,text) if(stat/=0)then;deallocate(var);err_msg='Var '//text//' in'//trim(err_msg);return;endif
#undef CHECK_ERR
#define CHECK_ERR(var)      if (stat /= 0) deallocate(var)

    i = set_indx(rad_set, satid=s%satid, grid=s%grid)
    if (i < 0) then
       err_msg = 'No tovs_obs_chan namelist corresponding to'//trim(err_msg)
       if (pe_loc==1) then
          print*,'ERROR in reading'
          print*,'err_msg: '//trim(err_msg)
       end if
       return
    end if

    s%gopts = rad_set(i)%gopts
    s%iopts = rad_set(i)%iopts

    s%n_im_cluster = 0
    s%n_im_chan    = 0
    n_cc           = 0
    if (s%gopts%use_imager > 0) then
      call get_dim(d_chcl, l_opt=.true.)
      call get_dim(d_cl, l_opt=.true.)
      s%n_im_cluster = d_cl  %len
      n_cc           = d_chcl%len
      if (s%n_im_cluster > 0) s%n_im_chan  =  n_cc / s%n_im_cluster
      if (s%n_im_chan > m_im_ch) then
        err_msg = 'too many imager channels in'//trim(err_msg)//'. Recompile with larger m_im_ch.'
        return
      end if
    end if

    ! obligatory data
    allocate(r% date  (1:n_meas))
    call get_var('DATE'            ,.true. , stat, idata1=r% date  (1:n_meas))
    CHECK_ERR_EXIT(r%date,'DATE')

    allocate(r% time  (1:n_meas))
    call get_var('TIME'            ,.true. , stat, idata1=r% time  (1:n_meas))
    CHECK_ERR_EXIT(r%time,'TIME')

    allocate(r% dlat  (1:n_meas))
    call get_var('LAT'             ,.true. , stat, rdata1=r% dlat  (1:n_meas))
    CHECK_ERR_EXIT(r%dlat,'LAT')

    allocate(r% dlon  (1:n_meas))
    call get_var('LON'             ,.true. , stat, rdata1=r% dlon  (1:n_meas))
    CHECK_ERR_EXIT(r%dlon,'LON')

    allocate(r% fov   (1:n_meas))
    call get_var('FOV'             ,.true. , stat, idata1=r% fov   (1:n_meas))
    CHECK_ERR_EXIT(r%fov,'FOV')
    WHERE(r% fov(1:n_meas) /= ifill) r% fov(1:n_meas) = r% fov(1:n_meas) - 1
                                            ! sat_pp convention: start with 1
                                            ! convention here  : start with 0

    allocate(r% scanl (1:n_meas))
    call get_var('SCANLINE'        ,.false., stat, idata1=r% scanl (1:n_meas))
    CHECK_ERR_EXIT(r%scanl,'SCANLINE')

    allocate(r% stzen (1:n_meas))
    call get_var('SAT_ZENITH'      ,.true. , stat, rdata1=r% stzen (1:n_meas))
    CHECK_ERR_EXIT(r%stzen,'SAT_ZENITH')

    allocate(r% stazi (1:n_meas))
    call get_var('SAT_AZIMUTH'     ,.true. , stat, rdata1=r% stazi (1:n_meas))
    CHECK_ERR_EXIT(r%stazi,'SAT_AZIMUTH')

    allocate(r% landfr(1:s%n_sens, 1:n_meas))
    call get_var('LAND_FRACTION'   ,.true. , stat, rdata2=r% landfr(1:n_dummy,1:n_meas))
    CHECK_ERR_EXIT(r%landfr,'LAND_FRACTION')

    allocate(r% stype(1:s%n_sens, 1:n_meas))
    call get_var('SURFACE_TYPE'    ,.true. , stat, idata2=r% stype (1:n_dummy,1:n_meas))
    CHECK_ERR_EXIT(r%stype,'SURFACE_TYPE')

    hgt_vname = 'SURFACE_HEIGHT'
    IF (s%gopts%l_max_surf_hgt) hgt_vname = 'MAX_'//trim(hgt_vname)
    allocate(r% shgt(1:s%n_sens, 1:n_meas))
    call get_var(hgt_vname         ,.true. , stat, rdata2=r% shgt  (1:n_dummy,1:n_meas))
    CHECK_ERR_EXIT(r%shgt,'*SURFACE_HEIGHT')

    allocate(r% bt_obs(1:n_chan,  1:n_meas))
    call get_var('BT_OBS'          ,.true. , stat, rdata2=r% bt_obs(1:n_chan,1:n_meas))
    CHECK_ERR_EXIT(r%bt_obs,'BT_OBS')

    allocate(chan_blkl(1:n_chan,  1:n_meas))
    call get_var('CHAN_BLACKLIST'  ,.true. , stat, idata2=chan_blkl(1:n_chan,1:n_meas))
    CHECK_ERR_EXIT(chan_blkl,'CHAN_BLACKLIST')
    allocate(r% valid(1:n_chan,  1:n_meas))
    where(chan_blkl(:,:) == 1)
      r% valid(:,:) = .false.
    elsewhere
      r% valid(:,:) = .true.
    end where

    allocate(start_sensor(1:n_instr))
    call get_var('START_SENSOR', .true. , stat, idata1=start_sensor(:))
    CHECK_ERR_EXIT(start_sensor,'START_SENSOR')

    allocate(end_sensor(1:n_instr))
    call get_var('END_SENSOR', .true. , stat, idata1=end_sensor(:))
    CHECK_ERR_EXIT(end_sensor,'END_SENSOR')

    !-------------------
    ! nonobligatory data
    !-------------------
    allocate(r% sunzen(1:n_meas))
    call get_var('SOLAR_ZENITH'    ,.false., stat, rdata1=r% sunzen(1:n_meas))
    CHECK_ERR(r%sunzen)

    allocate(r% sunazi(1:n_meas))
    call get_var('SOLAR_AZIMUTH'   ,.false., stat, rdata1=r% sunazi(1:n_meas))
    if (stat /= 0) then ! For some reason this field was renamed in newer files
      call get_var('MDAS' ,.false., stat, rdata1=r% sunazi(1:n_meas))
    endif
    CHECK_ERR(r%sunazi)

    allocate(r% orb_phase(1:n_meas))
    call get_var('ORBIT_PHASE'    ,.false., stat, rdata1=r% orb_phase(1:n_meas))
    CHECK_ERR(r%orb_phase)

    allocate(r% instr_temp(1:s%n_sens,1:n_meas))
    call get_var('INSTR_TEMP'    ,.false., stat, rdata2=r% instr_temp(1:n_dummy, 1:n_meas))
    CHECK_ERR(r%instr_temp)

    allocate(r% obsnum(1:n_meas))
    call get_var('OBSNUM'          ,.false., stat, idata1=r% obsnum(1:n_meas))
    CHECK_ERR(r%obsnum)

    select case (satpp_version)
      case (0)
        allocate(r% center(1:n_meas))
        call get_var('CENTRE_1B'       ,.false., stat, idata1=r% center(1:n_meas))
        CHECK_ERR(r%center)

        allocate(r% subcenter(1:n_meas))
        call get_var('SUBCENTRE_1B'    ,.false., stat, idata1=r% subcenter(1:n_meas))
        CHECK_ERR(r%subcenter)
      case (1)
        allocate(r% center(1:n_meas))
        call get_var('CENTRE_1B'       ,.false., stat, idata0=r% center(1))
        CHECK_ERR(r%center)

        allocate(r% subcenter(1:n_meas))
        call get_var('SUBCENTRE_1B'    ,.false., stat, idata0=r% subcenter(1))
        CHECK_ERR(r%subcenter)
      case default
        call Error('Unknown satpp_version','INPUT')
        return
    end select

    allocate(r% date_d(1:n_meas))
    call get_var('SECTION2_DECODING_DATE' ,.false., stat, idata1=r% date_d(1:n_meas))
    CHECK_ERR(r%date_d)

    allocate(r% time_d(1:n_meas))
    call get_var('SECTION2_DECODING_TIME' ,.false., stat, idata1=r% time_d(1:n_meas))
    CHECK_ERR(r%time_d)

    allocate(r% pp_flags(1:n_instr,1:n_meas))
    call get_var('PP_FLAGS' ,.false., stat, idata2=r% pp_flags(1:n_instr,1:n_meas))
    CHECK_ERR(r% pp_flags)

    allocate(idummy(1:n_instr,1:n_meas),r% cld_frc(1:n_instr,1:n_meas))
    call get_var('CLOUD_COVER' ,.false., stat, idata2=idummy)
    if (stat == 0) then
      r% cld_frc = idummy
    else
      deallocate(idummy, r%cld_frc)
      allocate(idummy(1:n_chan,1:n_meas),r% cld_frc(1:n_chan,1:n_meas))
      call get_var('CLOUD_COVER' ,.false., stat, idata2=idummy)
      if (stat == 0) then
        r% cld_frc = idummy
      else
        deallocate(r%cld_frc)
      end if
    end if
    deallocate(idummy)

    !---------------
    ! Optional stuff
    !---------------
    allocate(r% nwc_flg(1:n_meas))
    call get_var('NWC_FLAG' ,.false., stat, idata1=r% nwc_flg(1:n_meas))
    CHECK_ERR(r%nwc_flg)

    if (n_cc > 0 .and. s%n_im_cluster > 0) then
      call get_var('IMAGER_RTTOV_ID',.false., stat, idata0=s% im_rttov_id)
#define CHECK_ERR_EXIT(var,text) if(stat/=0)then;deallocate(var);err_msg='Read '//text//' in'//trim(err_msg);s%n_im_chan=0;return;endif
      if (stat == 0) then
        allocate(r% im_chan(1:n_cc))
        call get_var('IMAGER_CHAN'   ,.true., stat, idata1=r% im_chan(:))
        CHECK_ERR_EXIT(r% im_chan,'IMAGER_CHAN')
        if (s%im_rttov_id == 5 .and. maxval(r%im_chan(1:n_cc)) > 3) then
          !!! WARNING: bad hack of channel numbers/indices
          r%im_chan(1:n_cc) = r%im_chan(1:n_cc) - 3
          if (any(r%im_chan(1:n_cc) <= 0)) then
            err_msg = 'incompatible avhrr  channel numbers'
            return
          end if
        end if
        i = 0
        do j = 1, n_cc
          if (.not.any(s%im_chan(1:i) == r%im_chan(j))) then
            i = i + 1
            if (i > s%n_im_chan) then
              err_msg = 'imager channel numbers inconsistent (1)'
              return
            end if
            s%im_chan(i) = r%im_chan(j)
          end if
        end do
        if (i < s%n_im_chan) then
          err_msg = 'imager channel numbers inconsistent (2)'
          return
        end if
        allocate(r% im_cl_mean(1:n_cc,1:n_meas))
        call get_var('IMAGER_CL_MEAN',.true., stat, rdata2=r% im_cl_mean(:,:))
        CHECK_ERR_EXIT(r% im_cl_mean,'IMAGER_CL_MEAN')
        allocate(r% im_cl_stdv(1:n_cc,1:n_meas))
        call get_var('IMAGER_CL_STDV',.true., stat, rdata2=r% im_cl_stdv(:,:))
        CHECK_ERR_EXIT(r% im_cl_stdv,'IMAGER_CL_STDV')
        allocate(r% im_cl_frac(1:s%n_im_cluster,1:n_meas))
        call get_var('IMAGER_CL_FRACTION',.true., stat, idata2=r% im_cl_frac(:,:))
        CHECK_ERR_EXIT(r% im_cl_frac,'IMAGER_CL_FRAC')
      end if
    end if

    ! Close the file
    call nf_close
    if (nf_stat /= NF_NOERR) then
       err_msg = 'Failed to close'//trim(err_msg)
       return
    end if

    ! Prepare instrument mapping for sensors
    allocate(mapping(s%n_sens))
    do i = 1, n_instr
      mapping(start_sensor(i):end_sensor(i))        = i
    end do

    allocate(s%sensor_instr(s%n_sens))
    s%sensor_instr(:) = s%instr(mapping(:))

    ! Adapt to different sat_pp versions
    select case (satpp_version)
      case (0)
        ! some fields in old satpp are stored in dimension
        ! n_instr
        if (n_instr < s%n_sens) then
          r%landfr(:,:) = r%landfr(mapping(:),:)
          r%stype(:,:)  = r%stype(mapping(:),:)
          r%shgt(:,:)   = r%shgt(mapping(:),:)
        endif
      case (1)
        ! center and subcenter are no longer an array
        if (associated(r%center))    r%center(2:n_meas)    = r%center(1)
        if (associated(r%subcenter)) r%subcenter(2:n_meas) = r%subcenter(1)
      case default
        call Error('Unknown satpp_version','INPUT')
        err_msg = 'Unknown satpp_version in'//trim(err_msg)
        return
    end select

    if (present(status)) status = 0
    valid = .true.

  contains

    subroutine get_var(name, required, stat, idata0, idata1, idata2, rdata1, rdata2)
      character(len=*), intent(in)            :: name
      logical,          intent(in)            :: required
      integer,          intent(out)           :: stat
      integer,          intent(out), optional :: idata0
      integer,          intent(out), optional :: idata1(:)
      integer,          intent(out), optional :: idata2(:,:)
      real(wp),         intent(out), optional :: rdata1(:)
      real(wp),         intent(out), optional :: rdata2(:,:)

      character(len=180)      :: attname
      type(t_nc_dim), pointer :: d => null()
      integer                 :: shape_arg(2)
      integer                 :: dimids(2)
      integer                 :: start(2)
      integer                 :: icount(2)
      integer                 :: n_dims, n_dims_nc, n_atts
      integer                 :: type, type_arg
      integer                 :: ncid_var
      integer                 :: i

      stat = 1

      ! Get info about the arguments
      if     (present(idata0)) then
        n_dims = 0
        type_arg = NF_INT
      elseif (present(idata1)) then
        n_dims = 1
        shape_arg(1:1) = shape(idata1)
        type_arg = NF_INT
      elseif (present(idata2)) then
        n_dims = 2
        shape_arg(1:2) = shape(idata2)
        type_arg = NF_INT
      elseif (present(rdata1)) then
        n_dims = 1
        shape_arg(1:1) = shape(rdata1)
        type_arg = NF_FLOAT
      elseif (present(rdata2)) then
        n_dims = 2
        shape_arg(1:2) = shape(rdata2)
        type_arg = NF_FLOAT
      else
        stat = -99
        call Error('invalid arguments to get_var', 'INTERNAL')
        return
      end if

      ! inquire the variable
      nf_stat = nf90_inq_varid(ncid, name, ncid_var)
      if (nf_stat /= NF_NOERR) then
        if (required) call nf_err('nf90_inq_varid', 'Failed to get NetCDF-ID of variable '//trim(name))
        return
      end if
      nf_stat = nf90_inquire_variable(ncid, ncid_var, &
           &xtype=type, ndims=n_dims_nc, dimids=dimids, natts=n_atts)
      if (nf_stat /= NF_NOERR) then
        if (required) call nf_err('nf90_inquire_varible', 'Failed to inquire variable '//trim(name))
        return
      end if

      ! Check whether arguments are in agreement with the variable
      if (n_dims /= n_dims_nc) then
        if (required) call error('Number of dimensions of variables '//trim(name)//&
             &' is not as expected.', 'INPUT')
        return
      end if
      if ( ((type_arg == NF_INT)  .and.(type /= NF_INT)) .or. &
           ((type_arg == NF_FLOAT).and.(type /= NF_FLOAT).and.(type /= NF_DOUBLE)) ) then
        if (required) call error('Type of variable '//trim(name)//' is not as expected', 'INPUT')
        return
      end if

      ! Get Fillvalue
      ifill = NF_FILL_INT
      rfill = NF_FILL_FLOAT
      do i = 1, n_atts
        nf_stat = nf90_inq_attname(ncid, ncid_var, i, attname)
        if (nf_stat /= NF_NOERR) then
          if (required) call nf_err('nf90_inq_attname', &
               &'Failed to inquire attributes of variable '//trim(name))
          return
        end if
        if (attname == '_FillValue') then
          select case (type)
          case (NF_INT)
            nf_stat = nf90_get_att(ncid, ncid_var, attname, ifill)
          case (NF_FLOAT, NF_DOUBLE)
            nf_stat = nf90_get_att(ncid, ncid_var, attname, rfill)
          end select
          if (nf_stat /= NF_NOERR) then
            if (required) call nf_err('nf90_get_att', &
                 &'Failed to get attribute "_FillValue" of variable '//trim(name))
            return
          end if
          exit
        end if
      end do

      ! check whether dimensions are in agreement and set up start and count arrays
      do i = 1, n_dims
        d => find_dim(dimids(i))
        if (.not.associated(d)) then
          if (required) call Error('Variable '//trim(name)//' has invalid dimensions.', 'INPUT')
          return
        end if
        if (d% ncid == d_meas% ncid) then
          start(i)  = i_s
          icount(i) = i_e - i_s + 1
        else
          start(i)  = 1
          icount(i) = d%len
        end if
        if (icount(i) > shape_arg(i)) then
          if (required) call Error('Dimension lengths of variable '//trim(name)//&
               &' not as expected.', 'INPUT')
          return
        end if
      end do

      ! Read the variable
      if     (present(idata0)) then
        nf_stat = nf90_get_var(ncid, ncid_var, idata0)
        if (present(i_inv) .and. (idata0==ifill)) idata0 = i_inv
      elseif (present(idata1)) then
        nf_stat = nf90_get_var(ncid, ncid_var, idata1(1:icount(1)), &
             &start=start(1:1), count=icount(1:1))
        if (present(i_inv)) &
             where(idata1(1:icount(1)) == ifill) idata1(1:icount(1)) = i_inv
      elseif (present(idata2)) then
        nf_stat = nf90_get_var(ncid, ncid_var, idata2(1:icount(1),1:icount(2)), &
             &start=start(1:2), count=icount(1:2))
        if (present(i_inv)) &
             where(idata2(1:icount(1),1:icount(2)) == ifill) &
             idata2(1:icount(1),1:icount(2)) = i_inv
      elseif (present(rdata1)) then
        nf_stat = nf90_get_var(ncid, ncid_var, rdata1(1:icount(1)), &
             &start=start(1:1), count=icount(1:1))
         if (present(r_inv)) &
             where(rdata1(1:icount(1)) == rfill) rdata1(1:icount(1)) = r_inv
      elseif (present(rdata2)) then
        nf_stat = nf90_get_var(ncid, ncid_var, rdata2(1:icount(1),1:icount(2)), &
             &start=start(1:2), count=icount(1:2))
        if (present(r_inv)) &
             where(rdata2(1:icount(1),1:icount(2)) == rfill) &
             rdata2(1:icount(1),1:icount(2)) = r_inv
      end if
      if (nf_stat /= 0) then
        if (required) call nf_err('nf90_get_var', 'Failed to read variable '//trim(name))
        return
      end if

      stat = 0

    end subroutine get_var


    ! The following function shall determine the 'version' of the sat_pp
    ! File. At the time superobbing was invented, the dimension of several
    ! fields within the sat_pp file was changed

    subroutine detect_satpp_version(version,stat)
      integer, intent(out) :: stat
      integer, intent(out) :: version

      integer                 :: ncid_var,ncid_dim
      integer                 :: dimids(2)

      stat    = 1
      version = 0

      ! inquire the variable
      nf_stat = nf90_inq_varid(ncid, 'SENSOR', ncid_var)
      if (nf_stat /= NF_NOERR) then
        call nf_err('nf90_inq_varid', 'Failed to get NetCDF-ID of variable '//trim('SENSOR'))
        return
      end if

      nf_stat = nf90_inquire_variable(ncid, ncid_var, dimids=dimids)
      if (nf_stat /= NF_NOERR) then
        call nf_err('nf90_inquire_varible', 'Failed to inquire variable '//trim('SENSOR'))
        return
      end if

      nf_stat = nf90_inq_dimid(ncid, 'nSens', ncid_dim)
      if (nf_stat /= NF_NOERR) then
        call nf_err('nf90_inq_imid', 'Failed to inquire dimension id '//trim('nSens'))
        return
      end if

      if (ncid_dim == dimids(1)) then
        version = 1
      endif

      stat = 0
    end subroutine detect_satpp_version

    function find_dim(id)
      integer,        intent(in) :: id
      type(t_nc_dim), pointer    :: find_dim

      if     (id == d_meas%ncid) then
        find_dim => d_meas
      elseif (id == d_instr%ncid) then
        find_dim => d_instr
      elseif (id == d_sens%ncid) then
        find_dim => d_sens
      elseif (id == d_chan%ncid) then
        find_dim => d_chan
      elseif (id == d_chan%ncid) then
        find_dim => d_chan
      elseif (id == d_chcl%ncid) then
        find_dim => d_chcl
      elseif (id == d_cl%ncid) then
        find_dim => d_cl
      else
        find_dim => null()
      end if

    end function find_dim

    subroutine nf_close
      nf_stat = nf90_close(ncid)
      if (nf_stat /= NF_NOERR) call nf_err('nf90_close')
    end subroutine nf_close

    subroutine get_dim(dim, l_opt)
      type(t_nc_dim), intent(inout)          :: dim
      logical,        intent(in),   optional :: l_opt
      logical :: lopt

      if (present(l_opt)) then
        lopt = l_opt
      else
        lopt = .false.
      end if
      nf_stat = nf90_inq_dimid(ncid, dim%name, dim%ncid)
      if (nf_stat /= NF_NOERR) then
        if (.not.lopt) then
          call nf_err('nf90_inq_dim')
          call warning('Dimension '//trim(dim%name)//' not found in input file')
        else
          nf_stat = 0
        end if
        return
      endif
      nf_stat = nf90_inquire_dimension(ncid, dim%ncid, len=dim%len)
      if (nf_stat /= NF_NOERR) then
        call nf_err('nf90_inquire_dimension')
        call warning('Failed to get length of dimension '//trim(dim%name))
        return
      endif
      if (l_print) print '(3x,A,T40,A,I7)',cind//'dimension '//trim(dim%name),' : ',dim%len

    end subroutine get_dim

    subroutine nf_err(proc, msg)
      character(len=*), intent(in)           :: proc
      character(len=*), intent(in), optional :: msg

      call nf_error(nf_stat, proc, msg)
      if (present(status)) then
        status = nf_stat
      else
        call finish(trim(proc)//'@read_satpp_feedbk', trim(msg))
      end if

    end subroutine nf_err

    subroutine error(msg, type)
      character(len=*), intent(in) :: msg
      character(len=*), intent(in) :: type

      call warning('*** '//trim(type)//' ERROR : '//trim(msg), nowarn=.true.)
      if (present(status)) then
        status = nf_stat
      else
        call finish('read_satpp_feedbk', trim(type)//': '//trim(msg))
      end if

    end subroutine error

  end subroutine read_satpp_feedbk
!------------------------------------------------------------------------------

  function chan_indx (instr, ichan, set, ii)
  !--------------------------------------------------------
  ! figure out suitable channel index in sat/instr/chan set
  ! return -1 in case of no match
  !--------------------------------------------------------
  integer         ,intent(in) :: instr      ! instrument number
  integer         ,intent(in) :: ichan      ! channel number
  type(t_rad_set) ,intent(in) :: set        ! set of instruments/channels
  integer,optional,intent(out):: ii         ! instrument index
  integer                     :: chan_indx  ! channel index in set

    integer :: i  ! instrument index
    integer :: j  ! channel index loop variable
    integer :: ic ! channel number (instrument number stripped)
    !---------------------------------------------------
    ! strip instrument number as currently used in 3dvar
    !---------------------------------------------------
    select case (instr)
    case (3,4,15,19)
      ic = mod (ichan, f_ins)
    case default
      ic = ichan
    end select
    !-------------------
    ! find channel index
    !-------------------
    chan_indx = - 1
    if (present(ii)) ii = -1
    if (ic <= 0) return
    do i = 1, set% n_instr
      if (instr == set% instr (i)) then
        if (present(ii)) ii = i
        do j = set% o_ch_i(i)+1, set% o_ch_i(i) + set% n_ch_i(i)
          if (set% chan(j) == ichan) then
            chan_indx = j
            return
          endif
        enddo
      endif
    end do
  end function chan_indx
!------------------------------------------------------------------------------

  function set_indx (set, satid, grid)
  !---------------------------------------------
  ! figure out suitable sat/instr/chan set index
  ! return -1 in case of no match
  !----------------------------------------------
  integer                               :: set_indx  ! array element index in set(:)
  type(t_rad_set) ,intent(in)           :: set(:)    ! set of instruments/channels
  integer         ,intent(in), optional :: satid     ! satellite id
  integer         ,intent(in), optional :: grid      ! grid ID

    integer :: i
    set_indx = -1
    if (.not.present(satid) .and. .not.present(grid)) return
    do i = 1, size (set)
      if (present(satid)) then
        if (satid /= set(i)% satid) cycle
      end if
      if (present(grid)) then
        if (grid /= set(i)% grid) cycle
      end if
      set_indx = i
      return
    end do
  end function set_indx
!------------------------------------------------------------------------------

  function instr_chan (ichan, set) result (instr)
  integer         ,intent(in) :: ichan  ! channel index
  type(t_rad_set) ,intent(in) :: set    ! set of instruments
  integer                     :: instr  ! rttov instrument id
  !----------------------------------
  ! get instrument from channel index
  !----------------------------------
    integer :: i
    instr = - 1
    do i = 1, set% n_instr
      if (ichan <= set% o_ch_i(i) + set% n_ch_i(i)) then
        instr = set% instr(i)
        return
      endif
    end do
  end function instr_chan
!------------------------------------------------------------------------------

  elemental subroutine assign_rad_set (y, x)
    !===============================================
    ! assignment: t_rad_set = t_rad_set
    ! purpose: correctly allocate pointer components
    !===============================================
    type(t_rad_set) ,intent(inout) :: y
    type(t_rad_set) ,intent(in)    :: x

    integer :: i

    y% id            = x% id
    y% source        = x% source
    y% satid         = x% satid
    y% rttov_satid   = x% rttov_satid
    y% platform      = x% platform
    y% n_chan        = x% n_chan
    y% n_instr       = x% n_instr
    y% n_sens        = x% n_sens
    y% n_fov         = x% n_fov
    y% n_lev         = x% n_lev
    y% grid          = x% grid
    y% grid_wmo      = x% grid_wmo
    y% flag_instr    = x% flag_instr
    y% select_sens   = x% select_sens
!    y% n_pc_emiss    = x% n_pc_emiss
    y% instr         = x% instr
    y% instr_wmo     = x% instr_wmo
    y% n_ch_i        = x% n_ch_i
    y% o_ch_i        = x% o_ch_i
    y% rttov_indx    = x% rttov_indx
    y% gopts         = x% gopts
    y% iopts         = x% iopts
    if (associated(y% chan )) deallocate(y% chan )
    if (associated(x% chan )) then
      allocate(y% chan (size(x% chan )))
      y% chan        = x% chan
    end if
    if (associated(y% chidx)) deallocate(y% chidx)
    if (associated(x% chidx)) then
      allocate(y% chidx(size(x% chidx)))
      y% chidx       = x% chidx
    end if
    if (associated(y% iflag)) deallocate(y% iflag)
    if (associated(x% iflag)) then
      allocate(y% iflag(size(x% iflag)))
      y% iflag       = x% iflag
    end if
    if (associated(y% band )) deallocate(y% band )
    if (associated(x% band )) then
      allocate(y% band (size(x% band )))
      y% band        = x% band
    end if
    if (associated(y% flag )) deallocate(y% flag )
    if (associated(x% flag )) then
      allocate(y% flag (size(x% flag )))
      y% flag        = x% flag
    end if
    if (associated(y% var  )) deallocate(y% var  )
    if (associated(x% var  )) then
      allocate(y% var  (size(x% var  )))
      y% var         = x% var
    end if
    if (associated(y%sensor_instr)) deallocate(y%sensor_instr)
    if (associated(x%sensor_instr)) then
      allocate(y% sensor_instr  (size(x% sensor_instr)))
      y% sensor_instr(:) = x% sensor_instr(:)
    endif
    if (associated(y% oe_str)) deallocate(y% oe_str)
    if (associated(x% oe_str)) then
      allocate(y% oe_str(size(x% oe_str)))
      do i = 1, size(x% oe_str)
        y% oe_str(i) = x% oe_str(i)
      end do
    end if
    if (associated(y% oe_trf)) deallocate(y% oe_trf)
    if (associated(x% oe_trf)) then
      allocate(y% oe_trf(size(x% oe_trf)))
      do i = 1, size(x% oe_trf)
        y% oe_trf(i) = x% oe_trf(i)
      end do
    end if
    if (associated(y% wavenum)) deallocate(y% wavenum)
    if (associated(x% wavenum)) then
      allocate(y% wavenum(size(x% wavenum)))
      y% wavenum     = x% wavenum
    end if
    if (associated(y% emis_land)) deallocate(y% emis_land)
    if (associated(x% emis_land)) then
      allocate(y% emis_land(size(x% emis_land)))
      y% emis_land     = x% emis_land
    end if
    if (associated(y% emis_snow)) deallocate(y% emis_snow)
    if (associated(x% emis_snow)) then
      allocate(y% emis_snow(size(x% emis_snow)))
      y% emis_snow     = x% emis_snow
    end if
    if (associated(y% emis_sice)) deallocate(y% emis_sice)
    if (associated(x% emis_sice)) then
      allocate(y% emis_sice(size(x% emis_sice)))
      y% emis_sice     = x% emis_sice
    end if
    if (associated(y% emis_pc)) deallocate(y% emis_pc)
    if (associated(x% emis_pc)) then
      allocate(y% emis_pc(size(x% emis_pc,1),size(x% emis_pc,2)))
      y% emis_pc     = x% emis_pc
    end if
    if (associated(y% i1)) deallocate(y% i1)
    if (associated(x% i1)) then
      allocate(y% i1(size(x% i1)))
      y% i1     = x% i1
    end if
    if (associated(y% w1)) deallocate(y% w1)
    if (associated(x% w1)) then
      allocate(y% w1(size(x% w1)))
      y% w1     = x% w1
    end if
    if (associated(y% emis_ind)) deallocate(y% emis_ind)
    if (associated(x% emis_ind)) then
      allocate(y% emis_ind(size(x% emis_ind,1),0:size(x% emis_ind,2)-1,size(x% emis_ind,3)))
      y% emis_ind     = x% emis_ind
    end if
    y% n_emis_opt     = x% n_emis_opt
    if (associated(y% emis_opt)) deallocate(y% emis_opt)
    if (associated(x% emis_opt)) then
      allocate(y% emis_opt(size(x% emis_opt,1)))
      y% emis_opt     = x% emis_opt
    end if
    if (associated(y% bc% p_tr)) deallocate(y% bc% p_tr)
    if (associated(x% bc% p_tr)) then
      allocate(y% bc% p_tr(size(x% bc% p_tr,1),size(x% bc% p_tr,2)))
      y% bc% p_tr     = x% bc% p_tr
    end if
    if (associated(y% bc% ib_ch)) deallocate(y% bc% ib_ch)
    if (associated(x% bc% ib_ch)) then
      allocate(y% bc% ib_ch(size(x% bc% ib_ch,1)))
      y% bc% ib_ch     = x% bc% ib_ch
    end if
    if (associated(y% bc% type)) deallocate(y% bc% type)
    if (associated(x% bc% type)) then
      allocate(y% bc% type(size(x% bc% type,1)))
      y% bc% type     = x% bc% type
    end if
    if (associated(y% tskin_ind)) deallocate(y% tskin_ind)
    if (associated(x% tskin_ind)) then
      allocate(y% tskin_ind(size(x% tskin_ind,1),0:size(x% tskin_ind,2)-1,size(x% tskin_ind,3)))
      y% tskin_ind     = x% tskin_ind
    end if
    y% n_tskin_opt     = x% n_tskin_opt
    y% max_tskin_cld_cov = x% max_tskin_cld_cov
    y% chan_cld_cov      = x% chan_cld_cov
    if (associated(y% tskin_opt)) deallocate(y% tskin_opt)
    if (associated(x% tskin_opt)) then
      allocate(y% tskin_opt(size(x% tskin_opt,1)))
      y% tskin_opt     = x% tskin_opt
    end if
    y% n_im_chan     = x% n_im_chan
    y% n_im_cluster  = x% n_im_cluster
    y% im_rttov_id   = x% im_rttov_id
    y% im_chan       = x% im_chan
    y% im_chan_indx  = x% im_chan_indx

  end subroutine assign_rad_set


  elemental subroutine destruct_emis_opt(eo)
    type (t_emis_opt) ,intent(inout) :: eo
    if (associated(eo%chan)) deallocate(eo%chan)
    if (associated(eo%cdyn)) deallocate(eo%cdyn)
    call construct_emis_opt(eo)
  end subroutine destruct_emis_opt

    elemental subroutine destruct_tskin_opt(tso)
    type (t_tskin_opt) ,intent(inout) :: tso
    if (associated(tso%cdyn)) deallocate(tso%cdyn)
    call construct_tskin_opt(tso)
  end subroutine destruct_tskin_opt

  elemental subroutine construct_emis_opt(eo)
    type (t_emis_opt) ,intent(out) :: eo
  end subroutine construct_emis_opt

  elemental subroutine construct_tskin_opt(tso)
    type (t_tskin_opt) ,intent(out) :: tso
  end subroutine construct_tskin_opt

  elemental subroutine assign_emis_opt(y, x)
    type(t_emis_opt) ,intent(inout) :: y
    type(t_emis_opt) ,intent(in)    :: x

    call destruct_emis_opt(y)
    y%nml     = x%nml
    y%id      = x%id
    y%inst    = x%inst
    y%prio    = x%prio
    y%styp    = x%styp
    y%ns      = x%ns
    y%descr   = x%descr
    y%mode    = x%mode
    y%atls    = x%atls
    y%n_chan  = x%n_chan
    y%angcorr = x%angcorr
    y%max_dst = x%max_dst
    if (x%n_chan > 0) then
      allocate(y%chan(y%n_chan))
      y%chan(1:y%n_chan) = x%chan(1:x%n_chan)
      if (associated(x%cdyn)) then
        allocate(y%cdyn(y%n_chan))
        y%cdyn(1:y%n_chan) = x%cdyn(1:x%n_chan)
      end if
    end if

  end subroutine assign_emis_opt

  elemental subroutine assign_tskin_opt(y, x)
    type(t_tskin_opt) ,intent(inout) :: y
    type(t_tskin_opt) ,intent(in)    :: x

    call destruct_tskin_opt(y)
    y%nml     = x%nml
    y%id      = x%id
    y%inst    = x%inst
    y%prio    = x%prio
    y%styp    = x%styp
    y%ns      = x%ns
    y%descr   = x%descr
    y%mode    = x%mode
    y%n_cdyn  = x%n_cdyn
    if (x%n_cdyn > 0) then
      allocate(y%cdyn(y%n_cdyn))
      y%cdyn(1:y%n_cdyn) = x%cdyn(1:x%n_cdyn)
    end if

  end subroutine assign_tskin_opt

!------------------------------------------------------------------------------
  function reduced_rad_set (set, mask, oerr_par) result (r)
  type(t_rad_set)                       :: r
  type(t_rad_set), intent(in)           :: set
  logical,         intent(in)           :: mask(:)
  logical,         intent(in), optional :: oerr_par

    integer                     :: i, j, n_count
    logical                     :: l_oerr_par

    if (present(oerr_par)) then
      l_oerr_par = oerr_par
    else
      l_oerr_par = .true.
    end if

    if (size(mask) < set% n_chan) return

    r% id            = set% id
    r% source        = set% source
    r% satid         = set% satid
    r% rttov_satid   = set% rttov_satid
    r% platform      = set% platform
    r% grid          = set% grid
    r% grid_wmo      = set% grid_wmo
    r% flag_instr    = set% flag_instr
    r% select_sens   = set% select_sens
    r% n_fov         = set% n_fov
!    r% n_pc_emiss    = set% n_pc_emiss
    r% rttov_indx    = set% rttov_indx
    r% gopts         = set% gopts
    r% iopts         = set% iopts
    r% n_chan        = count(mask(1:set% n_chan))
    r% n_instr       = 0
    r% n_sens        = 0
    r% max_tskin_cld_cov  = set% max_tskin_cld_cov
    r% chan_cld_cov       = set% chan_cld_cov
    if (r% n_chan <= 0) return
    do i = 1, set% n_instr
      n_count = count(mask(set%o_ch_i(i)+1:set%o_ch_i(i)+set%n_ch_i(i)))
      if (n_count > 0) then
        r% n_instr = r% n_instr + 1
        r% instr         (r% n_instr) = set% instr     (i)
        r% instr_wmo     (r% n_instr) = set% instr_wmo (i)
        r% rttov_indx    (r% n_instr) = set% rttov_indx(i)
        r% n_ch_i        (r% n_instr) = n_count
        if (r% n_instr > 1) then
          r% o_ch_i      (r%n_instr)  = r% o_ch_i(r%n_instr-1) + r% n_ch_i(r%n_instr-1)
        else
          r% o_ch_i      (r%n_instr)  = 0
        end if
        r% iopts(r%n_instr)% l2c_type = set% iopts(i)% l2c_type
      end if
    end do
    if (associated(set% chan)) then
      allocate(r% chan (r%n_chan))
      r% chan  = pack(set% chan (1:set%n_chan), mask(1:set%n_chan))
    end if
    if (associated(set% chidx)) then
      allocate(r% chidx(r%n_chan))
      r% chidx = pack(set% chidx(1:set%n_chan), mask(1:set%n_chan))
    end if
    if (associated(set% iflag)) then
      allocate(r% iflag(r%n_chan))
      r% iflag = pack(set% iflag(1:set%n_chan), mask(1:set%n_chan))
    end if
    if (associated(set% band )) then
      allocate(r% band (r%n_chan))
      r% band  = pack(set% band (1:set%n_chan), mask(1:set%n_chan))
    end if
    if (associated(set% flag )) then
      allocate(r% flag (r%n_chan))
      r% flag  = pack(set% flag (1:set%n_chan), mask(1:set%n_chan))
    end if
    if (associated(set% var  )) then
      allocate(r% var  (r%n_chan))
      r% var   = pack(set% var  (1:set%n_chan), mask(1:set%n_chan))
    end if
    if (l_oerr_par) then
      if (associated(set% oe_str)) then
        allocate(r% oe_str(r%n_chan))
        r% oe_str = pack(set% oe_str(1:set%n_chan), mask(1:set%n_chan))
      end if
    end if
    if (associated(set% sensor_instr  )) then
      do i = 1, r%n_instr
        r%n_sens =r%n_sens + count(set%sensor_instr(:) == r%instr(i))
      end do
      allocate(r%sensor_instr(r%n_sens))
      j = 0
      do i = 1,set%n_sens
        if (any(set%sensor_instr(i) == r%instr(1:r%n_instr))) then
          j = j + 1
          r%sensor_instr(j) = set%sensor_instr(i)
        endif
      end do
    endif

  end function reduced_rad_set
!------------------------------------------------------------------------------

  subroutine construct_rad_set (set, satid, n_fov, grid, instr_1, chan_1, &
                                                         instr_2, chan_2, &
                                                         instr_3, chan_3  )
  !=======================================================
  ! preset components of type t_rad_set
  ! by passing list of instruments and associated channels
  !=======================================================
  type (t_rad_set)  ,intent(out) :: set        ! derived type variable to init.
  integer           ,intent(in)  :: satid     ! satellite id
  integer           ,intent(in)  :: n_fov      ! number of Fields Of View
  integer           ,intent(in)  :: grid       ! interpolation grid instrument
  integer           ,intent(in)  :: instr_1    ! 1st instrument used
  integer           ,intent(in)  :: chan_1 (:) ! channels used by 1st instr.
  integer ,optional ,intent(in)  :: instr_2    ! 2nd instrument used
  integer ,optional ,intent(in)  :: chan_2 (:) ! channels used by 2nd instr.
  integer ,optional ,intent(in)  :: instr_3    ! 3rd instrument used
  integer ,optional ,intent(in)  :: chan_3 (:) ! channels used by 3rd instr.

    set% id           = -1
    set% satid        = satid
    set% n_fov        = n_fov
    set% n_instr      = 1
    set% instr (1)    = instr_1
    set% grid         = grid
    set% o_ch_i(1)    = set% n_chan
    set% n_ch_i(1)    = size (chan_1)
    set% n_chan       = set% n_chan + set% n_ch_i(1)
!    set% n_pc_emiss   = -1
    call construct(set%gopts)
    call construct(set%iopts)
    if (present (instr_2)) then
      set% n_instr   = 2
      set% instr (2) = instr_2
      set% o_ch_i(2) = set% n_chan
      set% n_ch_i(2) = size (chan_2)
      set% n_chan    = set% n_chan + set% n_ch_i(2)
    endif
    if (present (instr_3)) then
      set% n_instr   = 3
      set% instr (3) = instr_3
      set% o_ch_i(3) = set% n_chan
      set% n_ch_i(3) = size (chan_3)
      set% n_chan    = set% n_chan + set% n_ch_i(3)
    endif
    allocate (set% chan (set% n_chan))
    set% chan (set% o_ch_i(1)+1:                     &
               set% o_ch_i(1)+set% n_ch_i(1)) = chan_1
    if (present (chan_2))                           &
      set% chan (set% o_ch_i(2)+1:                     &
                 set% o_ch_i(2)+set% n_ch_i(2)) = chan_2
    if (present (chan_3))                           &
      set% chan (set% o_ch_i(3)+1:                     &
                 set% o_ch_i(3)+set% n_ch_i(3)) = chan_3

  end subroutine construct_rad_set
!------------------------------------------------------------------------------
  subroutine construct_rad_gopts(x)
    type(t_rad_gopt), intent(out) :: x
  end subroutine construct_rad_gopts

  subroutine construct_rad_iopts(x)
    type(t_rad_iopt), intent(out) :: x(:)
  end subroutine construct_rad_iopts
!------------------------------------------------------------------------------
  subroutine construct_rad_set_1dvar (set, satid)
  !====================================
  ! preset components of type t_rad_set
  ! (setup as used in 1dvar)
  !====================================
  type (t_rad_set) ,intent(out) :: set
  integer          ,intent(in)  :: satid

    integer :: i

  integer, parameter  ::  n_fov_noaa  = 56  ! for NOAA satellites
  integer, parameter  ::  n_fov_aqua  = 30  ! for Aqua satellite

    set% satid = satid
    select case (satid)
    case (784)
      call construct (set, satid, n_fov_aqua, &! Aqua
                        3,  3, (/(i,i=1,15)/))  ! AMSUA
    case (209,223,004)
      call construct (set, satid, n_fov_noaa, &! NOAA 18,19, Metop
                        0,  0, (/(i,i=1,20)/), &! HIRS
                            3, (/(i,i=1,15)/), &! AMSUA
                           15, (/(i,i=1, 5)/))  ! MHS
    case (206,207,208)
      call construct (set, satid, n_fov_noaa, &! NOAA 15,16,17
                        0,  0, (/(i,i=1,20)/), &! HIRS
                            3, (/(i,i=1,15)/), &! AMSUA
                            4, (/(i,i=1, 5)/))  ! AMSUB
    end select

  end subroutine construct_rad_set_1dvar
!------------------------------------------------------------------------------
  elemental subroutine destruct_rad_set (rad_set)
  type (t_rad_set) ,intent(inout) :: rad_set
    if (associated (rad_set% chan        )) deallocate (rad_set% chan        )
    if (associated (rad_set% chidx       )) deallocate (rad_set% chidx       )
    if (associated (rad_set% iflag       )) deallocate (rad_set% iflag       )
    if (associated (rad_set% band        )) deallocate (rad_set% band        )
    if (associated (rad_set% flag        )) deallocate (rad_set% flag        )
    if (associated (rad_set% var         )) deallocate (rad_set% var         )
    if (associated (rad_set% sensor_instr)) deallocate (rad_set% sensor_instr)
    if (associated (rad_set% wavenum     )) deallocate (rad_set% wavenum     )
    if (associated (rad_set% emis_land   )) deallocate (rad_set% emis_land   )
    if (associated (rad_set% emis_snow   )) deallocate (rad_set% emis_snow   )
    if (associated (rad_set% emis_sice   )) deallocate (rad_set% emis_sice   )
    if (associated (rad_set% emis_pc     )) deallocate (rad_set% emis_pc     )
    if (associated (rad_set% i1          )) deallocate (rad_set% i1          )
    if (associated (rad_set% w1          )) deallocate (rad_set% w1          )
    if (associated (rad_set% emis_ind    )) deallocate (rad_set% emis_ind    )
    if (associated (rad_set% emis_opt    )) deallocate (rad_set% emis_opt    )
    if (associated (rad_set% tskin_ind   )) deallocate (rad_set% tskin_ind   )
    if (associated (rad_set% tskin_opt   )) deallocate (rad_set% tskin_opt   )
    if (associated (rad_set% oe_str      )) deallocate (rad_set% oe_str      )
    if (associated (rad_set% oe_trf      )) then
      call destruct(rad_set% oe_trf(:))
      deallocate(rad_set% oe_trf)
    end if
    if (associated (rad_set% bc% p_tr    )) deallocate (rad_set% bc% p_tr    )
    if (associated (rad_set% bc% ib_ch   )) deallocate (rad_set% bc% ib_ch   )
    if (associated (rad_set% bc% type    )) deallocate (rad_set% bc% type    )
    call construct_rs          (rad_set)

  contains

    elemental subroutine construct_rs (rs)
      type(t_rad_set), intent(out) :: rs
    end subroutine construct_rs

  end subroutine destruct_rad_set

!------------------------------------------------------------------------------
  !> \todo this is quite specific for COSMO. Implement a more general code.
  subroutine construct_radv(radv , nobs, stat, tmpl)
    !------------------------------------------------------------------------------
    !
    ! Description:
    !   This subroutine allocates the memory for rad derived type
    !------------------------------------------------------------------------------

    ! Parameters
    type(t_radv), intent(out)          :: radv  ! t_radv to be constructed
    integer,      intent (in)          :: nobs  ! number of observables
    integer,      intent(out)          :: stat  ! error status variable
    type(t_radv), intent(in), optional :: tmpl  ! template
    ! Local variables
    integer                   :: n_instr, n_chan, n_lev, n_sens, n_cld, n_pc, n_trg, n_im

    !-------------------------------------------------------------------------------
    !- Section 1: Allocate Memory
    !-------------------------------------------------------------------------------
    allocate(radv%ntstep(nobs), &
             radv%i_box(nobs),  &
             radv%j_box(nobs),  &
             radv%date(nobs),   &
             radv%time(nobs),   &
             radv%dlat(nobs),   &
             radv%dlon(nobs),   &
             radv%fov (nobs),   &
             radv%stzen (nobs), &
             radv%stazi (nobs), &
             stat = stat)

    if (present(tmpl) .and. stat == 0) then
      n_instr = tmpl%i%n_instr
      n_chan  = tmpl%i%n_chan
      n_sens  = tmpl%i%n_sens
      n_lev   = tmpl%n_lev
      n_pc    = tmpl%n_pc_emiss
      n_trg   = tmpl%n_trg

      allocate(radv%landfr(n_sens,nobs), &
               radv%stype(n_sens,nobs),  &
               radv%shgt(n_sens,nobs),   &
               radv%valid(n_chan,nobs),  &
               radv%bt_obs(n_chan,nobs), &
               radv%sinfl(n_chan,nobs),  &
               stat = stat)

      if (associated(tmpl%i_reprt)    .and. stat == 0) allocate(radv%i_reprt(nobs),         stat = stat)
      if (associated(tmpl%obsnum)     .and. stat == 0) allocate(radv%obsnum(nobs),          stat = stat)
      if (associated(tmpl%date_d)     .and. stat == 0) allocate(radv%date_d(nobs),          stat = stat)
      if (associated(tmpl%time_d)     .and. stat == 0) allocate(radv%time_d(nobs),          stat = stat)
      if (associated(tmpl%scanl)      .and. stat == 0) allocate(radv%scanl(nobs),           stat = stat)
      if (associated(tmpl%sunzen)     .and. stat == 0) allocate(radv%sunzen(nobs),          stat = stat)
      if (associated(tmpl%sunazi)     .and. stat == 0) allocate(radv%sunazi(nobs),          stat = stat)
      if (associated(tmpl%orb_phase)  .and. stat == 0) allocate(radv%orb_phase(nobs),       stat = stat)
      if (associated(tmpl%instr_temp) .and. stat == 0) allocate(radv%instr_temp(n_sens, nobs),stat = stat)
      if (associated(tmpl%center)     .and. stat == 0) allocate(radv%center(nobs),          stat = stat)
      if (associated(tmpl%subcenter)  .and. stat == 0) allocate(radv%subcenter(nobs),       stat = stat)
      if (associated(tmpl%mdlsfc)     .and. stat == 0) allocate(radv%mdlsfc(nobs),          stat = stat)
      if (associated(tmpl%op_na )     .and. stat == 0) allocate(radv%op_na(nobs),           stat = stat)
      if (associated(tmpl%nwc_flg)    .and. stat == 0) allocate(radv%nwc_flg(nobs),         stat = stat)
      if (associated(tmpl%r_state)    .and. stat == 0) allocate(radv%r_state(nobs),         stat = stat)
      if (associated(tmpl%specularity).and. stat == 0) allocate(radv%specularity(n_chan,nobs), stat = stat)
      if (associated(tmpl%state)      .and. stat == 0) allocate(radv%state(n_chan,nobs),    stat = stat)
      if (associated(tmpl%flags)      .and. stat == 0) allocate(radv%flags(n_chan,nobs),    stat = stat)
      if (associated(tmpl%bt_fg)      .and. stat == 0) allocate(radv%bt_fg(n_chan,nobs),    stat = stat)
      if (associated(tmpl%bt_bcor)    .and. stat == 0) allocate(radv%bt_bcor(n_chan,nobs),  stat = stat)
      if (associated(tmpl%bcor_)      .and. stat == 0) allocate(radv%bcor_(n_chan,nobs),    stat = stat)
      if (associated(tmpl%bt_fg_cs)   .and. stat == 0) allocate(radv%bt_fg_cs(n_chan,nobs), stat = stat)
      if (associated(tmpl%rad_fg)     .and. stat == 0) allocate(radv%rad_fg(n_chan,nobs),   stat = stat)
      if (associated(tmpl%rad_fg_cs)  .and. stat == 0) allocate(radv%rad_fg_cs(n_chan,nobs),stat = stat)
      if (associated(tmpl%p)          .and. stat == 0) allocate(radv%p(n_lev, nobs),        stat = stat)
      if (associated(tmpl%hh)         .and. stat == 0) allocate(radv%hh(n_lev, nobs),       stat = stat)
      if (associated(tmpl%h_fg)       .and. stat == 0) allocate(radv%h_fg(n_instr,nobs),    stat = stat)
      if (associated(tmpl%t_fg)       .and. stat == 0) allocate(radv%t_fg(n_lev,nobs),      stat = stat)
      if (associated(tmpl%q_fg)       .and. stat == 0) allocate(radv%q_fg(n_lev,nobs),      stat = stat)
      if (associated(tmpl%t2m)        .and. stat == 0) allocate(radv%t2m(nobs),             stat = stat)
      if (associated(tmpl%q2m)        .and. stat == 0) allocate(radv%q2m(nobs),             stat = stat)
      if (associated(tmpl%ps_fg)      .and. stat == 0) allocate(radv%ps_fg(nobs),           stat = stat)
      if (associated(tmpl%gp_fg)      .and. stat == 0) allocate(radv%gp_fg(nobs),           stat = stat)
      if (associated(tmpl%ts_fg)      .and. stat == 0) allocate(radv%ts_fg(nobs),           stat = stat)
      if (associated(tmpl%tsm_fg)     .and. stat == 0) allocate(radv%tsm_fg(nobs),          stat = stat)
      if (associated(tmpl%u10_fg)     .and. stat == 0) allocate(radv%u10_fg(nobs),          stat = stat)
      if (associated(tmpl%v10_fg)     .and. stat == 0) allocate(radv%v10_fg(nobs),          stat = stat)
      if (associated(tmpl%v10_abs_fg) .and. stat == 0) allocate(radv%v10_abs_fg(nobs),      stat = stat)
      if (associated(tmpl%snw_frc)    .and. stat == 0) allocate(radv%snw_frc(nobs),         stat = stat)
      if (associated(tmpl%cld_top)    .and. stat == 0) allocate(radv%cld_top(nobs),         stat = stat)
      if (associated(tmpl%cfraction)  .and. stat == 0) allocate(radv%cfraction(nobs),       stat = stat)
      if (associated(tmpl%clwde)      .and. stat == 0) allocate(radv%clwde(n_lev-1,nobs),   stat = stat)
      if (associated(tmpl%icede)      .and. stat == 0) allocate(radv%icede(n_lev-1,nobs),   stat = stat)
      if (associated(tmpl%cfrac)      .and. stat == 0) allocate(radv%cfrac(n_lev-1,nobs),   stat = stat)
      if (associated(tmpl%cld_fg)     .and. stat == 0) then
        n_cld = size(tmpl%cld_fg, 1)
        allocate(radv%cld_fg(n_cld,n_lev-1,nobs), stat = stat)
      end if
      if (associated(tmpl%cld_frc)   .and. stat == 0) allocate(radv%cld_frc(n_lev-1,nobs),  stat = stat)
      if (associated(tmpl%clw    )   .and. stat == 0) allocate(radv%clw    (n_lev-1,nobs),  stat = stat)
      if (associated(tmpl%trg)       .and. stat == 0) allocate(radv%trg      (n_trg , n_lev, nobs), stat = stat)
      if (associated(tmpl%H_t)       .and. stat == 0) allocate(radv%H_t      (n_chan, n_lev, nobs), stat = stat)
      if (associated(tmpl%H_q)       .and. stat == 0) allocate(radv%H_q      (n_chan, n_lev, nobs), stat = stat)
      if (associated(tmpl%H_ts)      .and. stat == 0) allocate(radv%H_ts     (n_chan,        nobs), stat = stat)
      if (associated(tmpl%H_ps)      .and. stat == 0) allocate(radv%H_ps     (n_chan,        nobs), stat = stat)
      if (associated(tmpl%H_em_pc)   .and. stat == 0) allocate(radv%H_em_pc  (n_chan, n_pc,  nobs), stat = stat)
      if (associated(tmpl%K_t)       .and. stat == 0) allocate(radv%K_t      (n_chan, n_lev, nobs), stat = stat)
      if (associated(tmpl%K_q)       .and. stat == 0) allocate(radv%K_q      (n_chan, n_lev, nobs), stat = stat)
      if (associated(tmpl%K_ts)      .and. stat == 0) allocate(radv%K_ts     (n_chan,        nobs), stat = stat)
      if (associated(tmpl%K_ps)      .and. stat == 0) allocate(radv%K_ps     (n_chan,        nobs), stat = stat)
      if (associated(tmpl%B_tt)      .and. stat == 0) allocate(radv%B_tt     (n_lev,  n_lev, nobs), stat = stat)
      if (associated(tmpl%B_qq)      .and. stat == 0) allocate(radv%B_qq     (n_lev,  n_lev, nobs), stat = stat)
      if (associated(tmpl%B_tq)      .and. stat == 0) allocate(radv%B_tq     (n_lev,  n_lev, nobs), stat = stat)
      if (associated(tmpl%R)         .and. stat == 0) allocate(radv%R        (n_chan,n_chan, nobs), stat = stat)
      if (associated(tmpl%HBHR)      .and. stat == 0) allocate(radv%HBHR     (n_chan,n_chan, nobs), stat = stat)
      if (associated(tmpl%t_eb)      .and. stat == 0) allocate(radv%t_eb     (n_lev,         nobs), stat = stat)
      if (associated(tmpl%q_eb)      .and. stat == 0) allocate(radv%q_eb     (n_lev,         nobs), stat = stat)
      if (associated(tmpl%ps_eb)     .and. stat == 0) allocate(radv%ps_eb    (               nobs), stat = stat)
      if (associated(tmpl%ts_eb)     .and. stat == 0) allocate(radv%ts_eb    (               nobs), stat = stat)
      if (associated(tmpl%cld_top_eb).and. stat == 0) allocate(radv%cld_top_eb(              nobs), stat = stat)
      if (associated(tmpl%cld_frc_eb).and. stat == 0) allocate(radv%cld_frc_eb(n_lev,        nobs), stat = stat)
      if (associated(tmpl%snw_frc_eb).and. stat == 0) allocate(radv%snw_frc_eb(              nobs), stat = stat)
      if (associated(tmpl%im_cl_mean) .and. stat == 0) then
        n_im = size(tmpl%im_cl_mean,1)
        allocate(radv%im_cl_mean(n_im, nobs), stat=stat)
      end if
      if (associated(tmpl%im_cl_stdv) .and. stat == 0) then
        n_im = size(tmpl%im_cl_stdv,1)
        allocate(radv%im_cl_stdv(n_im, nobs), stat=stat)
      end if
      if (associated(tmpl%im_cl_frac) .and. stat == 0) then
        n_im = size(tmpl%im_cl_frac,1)
        allocate(radv%im_cl_frac(n_im, nobs), stat=stat)
      end if
      if (associated(tmpl%im_chan) .and. stat == 0) then
        n_im = size(tmpl%im_chan,1)
        allocate(radv%im_chan(n_im), stat=stat)
      end if

    endif
    !-------------------------------------------------------------------------------
    !- Section 2: Assign values
    !-------------------------------------------------------------------------------
    if (stat == 0) then
      radv%n_rec      = nobs

      if (present(tmpl)) then
        radv%filename   = tmpl%filename
        radv%n_lev      = tmpl%n_lev
        radv%model_date = tmpl%model_date
        radv%i          = tmpl%i
      endif
    endif
    !------------------------------------------------------------------------------
    !- end of the subroutine
    !------------------------------------------------------------------------------
  end subroutine construct_radv

!------------------------------------------------------------------------------

  elemental subroutine destruct_radv (rad)
  !===============================
  ! free components of type t_radv
  !===============================
  type (t_radv) ,intent(inout) :: rad

    rad% filename   = ''
    rad% model_date = 0
    rad% n_rec      = 0
    rad% n_lev      = 0

    call destruct_rad_set (rad% i)

    !---------
    ! position
    !---------
    if (associated (rad% pe        )) deallocate (rad% pe        )
    if (associated (rad% ind       )) deallocate (rad% ind       )
    if (associated (rad% i_box     )) deallocate (rad% i_box     )
    if (associated (rad% j_box     )) deallocate (rad% j_box     )
    if (associated (rad% i_reprt   )) deallocate (rad% i_reprt   )
    if (associated (rad% i_body    )) deallocate (rad% i_body    )
    if (associated (rad% obsnum    )) deallocate (rad% obsnum    )
    !---------
    ! location
    !---------
    if (associated (rad% date      )) deallocate (rad% date      )
    if (associated (rad% time      )) deallocate (rad% time      )
    if (associated (rad% date_d    )) deallocate (rad% date_d    )
    if (associated (rad% time_d    )) deallocate (rad% time_d    )
    if (associated (rad% ntstep    )) deallocate (rad% ntstep    )
    if (associated (rad% dlon      )) deallocate (rad% dlon      )
    if (associated (rad% dlat      )) deallocate (rad% dlat      )
    if (associated (rad% fov       )) deallocate (rad% fov       )
    if (associated (rad% scanl     )) deallocate (rad% scanl     )
    if (associated (rad% stzen     )) deallocate (rad% stzen     )
    if (associated (rad% stazi     )) deallocate (rad% stazi     )
    if (associated (rad% sunzen    )) deallocate (rad% sunzen    )
    if (associated (rad% sunazi    )) deallocate (rad% sunazi    )
    if (associated (rad% orb_phase )) deallocate (rad% orb_phase )
    if (associated (rad% instr_temp)) deallocate (rad% instr_temp)
    if (associated (rad% landfr    )) deallocate (rad% landfr    )
    if (associated (rad% stype     )) deallocate (rad% stype     )
    if (associated (rad% shgt      )) deallocate (rad% shgt      )
    if (associated (rad% specularity)) deallocate (rad% specularity)
    !-------
    ! origin
    !-------
    if (associated (rad% center    )) deallocate (rad% center    )
    if (associated (rad% subcenter )) deallocate (rad% subcenter )
    !-----------
    ! predictors
    !-----------
    if (associated (rad% pred      )) deallocate (rad% pred      )
    if (associated (rad% tr        )) deallocate (rad% tr        )
    if (associated (rad% plev      )) deallocate (rad% plev      )
    !----------------
    ! quality control
    !----------------
    if (associated (rad% mdlsfc    )) deallocate (rad% mdlsfc    )
    if (associated (rad% op_na     )) deallocate (rad% op_na     )
    if (associated (rad% not_rej   )) deallocate (rad% not_rej   )
    if (associated (rad% cloudy    )) deallocate (rad% cloudy    )
    if (associated (rad% state     )) deallocate (rad% state     )
    if (associated (rad% flags     )) deallocate (rad% flags     )
    if (associated (rad% valid     )) deallocate (rad% valid     )
    if (associated (rad% sinfl     )) deallocate (rad% sinfl     )
    if (associated (rad% plevel    )) deallocate (rad% plevel    )
    if (associated (rad% pp_flags  )) deallocate (rad% pp_flags )
    if (associated (rad% nwc_flg   )) deallocate (rad% nwc_flg  )
    if (associated (rad% r_state   )) deallocate (rad% r_state  )
     !--------------------------------
    ! quantities in observation space
    !--------------------------------
    if (associated (rad% bt_fg     )) deallocate (rad% bt_fg     )
    if (associated (rad% bt_obs    )) deallocate (rad% bt_obs    )
    if (associated (rad% bt_bcor   )) deallocate (rad% bt_bcor   )
    if (associated (rad% bt_fg_cs  )) deallocate (rad% bt_fg_cs  )
    if (associated (rad% rad_fg    )) deallocate (rad% rad_fg    )
    if (associated (rad% rad_fg_cs )) deallocate (rad% rad_fg_cs )
    if (associated (rad% bcor_     )) deallocate (rad% bcor_     )
    if (associated (rad% emiss     )) deallocate (rad% emiss     )
    if (associated (rad% emis_pc   )) deallocate (rad% emis_pc   )
    if (associated (rad% e_fg      )) deallocate (rad% e_fg      )
    if (associated (rad% e_obs     )) deallocate (rad% e_obs     )
    !----------------------
    ! colocated imager info
    !----------------------
    if (associated (rad% im_cl_mean)) deallocate (rad% im_cl_mean)
    if (associated (rad% im_cl_stdv)) deallocate (rad% im_cl_stdv)
    if (associated (rad% im_cl_frac)) deallocate (rad% im_cl_frac)
    if (associated (rad% im_chan   )) deallocate (rad% im_chan   )
    !---------
    ! profiles
    !---------
    if (associated (rad% p         )) deallocate (rad% p         )
    if (associated (rad% hh        )) deallocate (rad% hh        )
    if (associated (rad% h_fg      )) deallocate (rad% h_fg      )
    if (associated (rad% t_fg      )) deallocate (rad% t_fg      )
    if (associated (rad% q_fg      )) deallocate (rad% q_fg      )
    if (associated (rad% trg       )) deallocate (rad% trg       )
    !------------------
    ! single level data
    !------------------
    if (associated (rad% t2m       )) deallocate (rad% t2m       )
    if (associated (rad% q2m       )) deallocate (rad% q2m       )
    if (associated (rad% ps_fg     )) deallocate (rad% ps_fg     )
    if (associated (rad% gp_fg     )) deallocate (rad% gp_fg     )
    if (associated (rad% ts_fg     )) deallocate (rad% ts_fg     )
    if (associated (rad% tsm_fg    )) deallocate (rad% tsm_fg    )
    if (associated (rad% u10_fg    )) deallocate (rad% u10_fg    )
    if (associated (rad% v10_fg    )) deallocate (rad% v10_fg    )
    if (associated (rad% v10_abs_fg)) deallocate (rad% v10_abs_fg)
    if (associated (rad% snw_frc   )) deallocate (rad% snw_frc   )
    if (associated (rad% cld_top   )) deallocate (rad% cld_top   )
    if (associated (rad% cld_fg    )) deallocate (rad% cld_fg    )
    if (associated (rad% cld_frc   )) deallocate (rad% cld_frc   )
    if (associated (rad% H_t       )) deallocate (rad% H_t       )
    if (associated (rad% H_q       )) deallocate (rad% H_q       )
    if (associated (rad% H_ts      )) deallocate (rad% H_ts      )
    if (associated (rad% H_ps      )) deallocate (rad% H_ps      )
    if (associated (rad% H_em_pc   )) deallocate (rad% H_em_pc   )
    if (associated (rad% K_t       )) deallocate (rad% K_t       )
    if (associated (rad% K_q       )) deallocate (rad% K_q       )
    if (associated (rad% K_ts      )) deallocate (rad% K_ts      )
    if (associated (rad% K_ps      )) deallocate (rad% K_ps      )
    if (associated (rad% B_tt      )) deallocate (rad% B_tt      )
    if (associated (rad% B_qq      )) deallocate (rad% B_qq      )
    if (associated (rad% B_tq      )) deallocate (rad% B_tq      )
    if (associated (rad% R         )) deallocate (rad% R         )
    if (associated (rad% HBHR      )) deallocate (rad% HBHR      )
    if (associated (rad% t_eb      )) deallocate (rad% t_eb      )
    if (associated (rad% q_eb      )) deallocate (rad% q_eb      )
    if (associated (rad% ts_eb     )) deallocate (rad% ts_eb     )
    if (associated (rad% ps_eb     )) deallocate (rad% ps_eb     )
    if (associated (rad% cld_top_eb)) deallocate (rad% cld_top_eb)
    if (associated (rad% cfraction )) deallocate (rad% cfraction )
    if (associated (rad% clwde     )) deallocate (rad% clwde     )
    if (associated (rad% icede     )) deallocate (rad% icede     )
    if (associated (rad% cfrac     )) deallocate (rad% cfrac     )
    if (associated (rad% clw       )) deallocate (rad% clw       )
    if (associated (rad% cld_frc_eb)) deallocate (rad% cld_frc_eb)
    if (associated (rad% snw_frc_eb)) deallocate (rad% snw_frc_eb)
    if (associated (rad% emis_pc_eb)) deallocate (rad% emis_pc_eb)

  end subroutine destruct_radv
!------------------------------------------------------------------------------

  subroutine print_rad_set(set, hint, header, indent, unit)
  !----------------------------
  ! printout of t_rad_set entry
  !----------------------------
  type(t_rad_set),  intent(in)           :: set
  integer,          intent(in), optional :: hint   ! integer to print in head line
  character(len=*), intent(in), optional :: header ! string to print in head line
  integer,          intent(in), optional :: indent ! number of indent blanks
  integer,          intent(in), optional :: unit   ! output unit (default: stdout)

    character(len=1000) :: str
    integer             :: i, j, l, ind, instr, iunit, mx_len, mx_bit, n_si

    if (present(unit)) then
      iunit = unit
    else
      iunit = stdout
    end if

    if (present(indent)) then
      ind = indent
    else
      ind = 0
    end if
    if (present(header)) then
      write(str,*) repeat(' ',ind)//trim(header)
      if (present(hint)) write(str(len_trim(str)+3:),*) hint
    elseif (present(hint)) then
      write(str,*) repeat(' ',ind)//'t_rad_set ',hint
    else
      write(str,*) repeat(' ',ind)//'t_rad_set'
    end if
    write(iunit,*) trim(str)
    write(iunit,*) repeat(' ',ind+2)//'id             =', set%id
    write(iunit,*) repeat(' ',ind+2)//'source         = ', trim(set%source)
    write(iunit,*) repeat(' ',ind+2)//'satid          =', set%satid
    write(iunit,*) repeat(' ',ind+2)//'rttov_satid    =', set%rttov_satid
    write(iunit,*) repeat(' ',ind+2)//'platform       =', set%platform
    write(iunit,*) repeat(' ',ind+2)//'grid           =', set%grid
    write(iunit,*) repeat(' ',ind+2)//'grid_wmo       =', set%grid_wmo
    write(iunit,*) repeat(' ',ind+2)//'flag_instr     =', set%flag_instr
    write(iunit,*) repeat(' ',ind+2)//'select_sens    =', set%select_sens
    write(iunit,*) repeat(' ',ind+2)//'n_instr        =', set%n_instr
    write(iunit,*) repeat(' ',ind+2)//'n_chan         =', set%n_chan
    write(iunit,*) repeat(' ',ind+2)//'n_fov          =', set%n_fov
!    write(iunit,*) repeat(' ',ind+2)//'n_pc_emiss     =', set%n_pc_emiss
    write(iunit,*) repeat(' ',ind+2)//'instr          =', set%instr     (1:set%n_instr)
    write(iunit,*) repeat(' ',ind+2)//'instr_wmo      =', set%instr_wmo (1:set%n_instr)
    write(iunit,*) repeat(' ',ind+2)//'n_ch_i         =', set%n_ch_i    (1:set%n_instr)
    write(iunit,*) repeat(' ',ind+2)//'o_ch_i         =', set%o_ch_i    (1:set%n_instr)
    write(iunit,*) repeat(' ',ind+2)//'rttov_indx     =', set%rttov_indx(1:set%n_instr)
    ! colocated imager
    write(iunit,*) repeat(' ',ind+2)//'im_rttov_id    =', set%im_rttov_id
    write(iunit,*) repeat(' ',ind+2)//'n_im_chan      =', set%n_im_chan
    write(iunit,*) repeat(' ',ind+2)//'n_im_cluster   =', set%n_im_cluster
    write(iunit,*) repeat(' ',ind+2)//'im_chan        =', set%im_chan(1:set%n_im_chan)
    write(iunit,*) repeat(' ',ind+2)//'im_chan_indx   =', set%im_chan_indx(1:set%n_im_chan)
    ! general options
    write(iunit,*) repeat(' ',ind+2)//'n_max_calc_k   =', set%gopts%n_max_calc_k
    write(iunit,*) repeat(' ',ind+2)//'n_max_calc_y   =', set%gopts%n_max_calc_y
    write(iunit,*) repeat(' ',ind+2)//'n_max_prof_k   =', set%gopts%n_max_prof_k
    write(iunit,*) repeat(' ',ind+2)//'n_max_prof_y   =', set%gopts%n_max_prof_y
    write(iunit,*) repeat(' ',ind+2)//'quality_mode   =', set%gopts%quality_mode
    write(iunit,*) repeat(' ',ind+2)//'pp_flags_instr =', set%gopts%pp_flags_instr
    write(iunit,*) repeat(' ',ind+2)//'pp_flags       =', set%gopts%pp_flags
    write(iunit,*) repeat(' ',ind+2)//'use_amsub_89ghz=', set%gopts%use_amsub_89ghz
    write(iunit,*) repeat(' ',ind+2)//'l_max_surf_hgt =', set%gopts%l_max_surf_hgt
    write(iunit,*) repeat(' ',ind+2)//'thin_superob_box_lines=', set%gopts%thin_superob_box_lines
    write(iunit,*) repeat(' ',ind+2)//'thin_superob_box_fovs =', set%gopts%thin_superob_box_fovs
    write(iunit,*) repeat(' ',ind+2)//'max_size_sobox        =', set%gopts%max_size_sobox
    write(iunit,*) repeat(' ',ind+2)//'max_timediff_sobox    =', set%gopts%max_timediff_sobox
    write(iunit,*) repeat(' ',ind+2)//'thin_superob_mode     =', set%gopts%thin_superob_mode
    write(iunit,*) repeat(' ',ind+2)//'lev_mode       =', set%gopts%lev_mode
    write(iunit,*) repeat(' ',ind+2)//'opt_vars       =', set%gopts%opt_vars
    write(iunit,*) repeat(' ',ind+2)//'use_imager     =', set%gopts%use_imager
    write(iunit,*) repeat(' ',ind+2)//'im_atlas_id    =', set%gopts%im_atlas_id
    write(iunit,*) repeat(' ',ind+2)//'im_tskin_opt   =', set%gopts%im_tskin_opt
    write(iunit,*) repeat(' ',ind+2)//'im_tskin_cfrac =', set%gopts%im_tskin_cfrac
    ! Emissivity
    write(iunit,*) repeat(' ',ind+2)//'n_emis_opts    =', set%n_emis_opt
    do i = 1, set%n_emis_opt
      write(iunit,*) repeat(' ',ind+2)//'emis_opts      =', trim(set%emis_opt(i)%descr),set%emis_opt(i)%inst
    end do
    ! Tskin
    write(iunit,*) repeat(' ',ind+2)//'max_tskin_cld_cov =', set%max_tskin_cld_cov
    write(iunit,*) repeat(' ',ind+2)//'chan_cld_cov      =', set%chan_cld_cov
    write(iunit,*) repeat(' ',ind+2)//'n_tskin_opts      =', set%n_tskin_opt
    do i = 1, set%n_tskin_opt
      write(iunit,*) repeat(' ',ind+2)//'tskin_opts      =', trim(set%tskin_opt(i)%descr),set%tskin_opt(i)%inst
    end do
    ! instrument specific options, arrays
    write(iunit,*) repeat(' ',ind+2)//'l2c_type       =', set%iopts%l2c_type
    write(iunit,*) repeat(' ',ind+2)//'l2c_rel_lim    =', set%iopts%l2c_rel_lim
    write(iunit,*) repeat(' ',ind+2)//'l2c_abs_lim    =', set%iopts%l2c_abs_lim
    write(iunit,*) repeat(' ',ind+2)//'l2c_use_rad    =', set%iopts%l2c_use_rad
    write(iunit,*) repeat(' ',ind+2)//'l2c_max        =', set%iopts%l2c_max
    write(iunit,*) repeat(' ',ind+2)//'l2c_max_trop   =', set%iopts%l2c_max_trop
    write(iunit,*) repeat(' ',ind+2)//'l2c_max_midlat =', set%iopts%l2c_max_midlat
    write(iunit,*) repeat(' ',ind+2)//'l2c_max_polar  =', set%iopts%l2c_max_polar
    do i = 1, set%n_instr
      write(iunit,*) repeat(' ',ind+2)//'instrument ',i,':'
      n_si = 0
      do j = 1, m_bd
        if (set%iopts(i)%surf_infl_mode(j) > 0) then
          n_si = j
        else
          exit
        end if
      end do
      if (associated(set%band)) then
        n_si = min(n_si, maxval(set%band(1:set%n_chan)))
        n_si = max(n_si, 1)
      end if
      write(iunit,'(a,*(5(1x,i9,:),/,21x))') repeat(' ',ind+5)//'surf_infl_mode =', set%iopts(i)%surf_infl_mode(1:n_si)
      write(iunit,'(a,*(5(1x,g9.2,:),/,21x))') repeat(' ',ind+5)//'max_surf_infl  =', set%iopts(i)%max_surf_infl(1:n_si)
      write(iunit,'(a,*(5(1x,g9.2,:),/,21x))') repeat(' ',ind+5)//'d_stemp        =', set%iopts(i)%d_stemp(1:n_si)
      write(iunit,'(a,*(5(1x,g9.2,:),/,21x))') repeat(' ',ind+5)//'d_emiss        =', set%iopts(i)%d_emiss(1:n_si)
      write(iunit,*) repeat(' ',ind+2)//'do_lambertian  =', set%iopts(i)%do_lambertian
      write(iunit,'(1x,a,*(5(1x,g9.2,:),/,21x))') repeat(' ',ind+2)//'specularity    =', set%iopts(i)%specularity(0:n_styp-1)
      write(iunit,*) repeat(' ',ind+2)//'specularity_snow=', set%iopts(i)%specularity_snow
    end do
    write(iunit,*) repeat(' ',ind+2)//'cloud_mode     =', set%iopts%cloud_mode
    write(iunit,*) repeat(' ',ind+2)//'rad_out        =', set%iopts%rad_out
    write(iunit,*) repeat(' ',ind+2)//'use_o3         =', set%iopts%use_o3
    write(iunit,*) repeat(' ',ind+2)//'use_co2        =', set%iopts%use_co2
    write(iunit,*) repeat(' ',ind+2)//'rt_mw_emis_mod =', set%iopts%rt_mw_emis_mod
    if (associated(set%chan)) then
      do mx_bit = 31, 0, -1
        if (name_use_bit(mx_bit) /= '') exit
      end do
      mx_len = 0
      do i = 0, mx_bit
        mx_len = max(len_trim(name_use_bit(i)), mx_len)
      end do

      do i = mx_len, 1, -1
        str = '|'//repeat(' |',mx_bit+1)
        do j = 0, mx_bit
          l = len_trim(name_use_bit(j))
          if (l >= i) str((j+1)*2:(j+1)*2) = name_use_bit(j)(l-i+1:l-i+1)
        end do
        if (i > 1) then
          write(iunit,'(1x,A)') repeat(' ',ind+2+34)//trim(str)
        else
          write(iunit,'(1x,A,2x,A6,2x,A6,2x,A6,2x,A6,2x,A,2x,A)') &
               &repeat(' ',ind+2),' index', ' instr', '  chan', '  band', &
               &trim(str), '  var(obserr)'
        end if
      end do
      msg = repeat(' ',34)//'|'//repeat(' |',mx_bit+1)
      write(iunit,*) repeat(' ',ind+2)//trim(msg)
      instr = 1
      do i = 1, size(set%chan)
        if (i > set%o_ch_i(instr)+set%n_ch_i(instr) .and. instr<set%n_instr) instr=instr+1
        msg = ''
        write(msg( 5: 8), '(I4)') i
        write(msg(11:16), '(I6)') set%instr(instr)
        write(msg(19:24), '(I6)') set%chan(i)
        if (associated(set%band )) write(msg(27:32), '(I6)'   ) set%band (i)
!        if (associated(set%var  )) write(msg(35:47), '(F13.7)') set%var  (i)
        if (associated(set%flag )) then
          msg = msg(1:34)//'|'//repeat(' |',mx_bit+1)
          do j = 0, mx_bit
            if (btest(set%flag(i), j)) msg(34+(j+1)*2:34+(j+1)*2) = 'y'
          end do
          do j = 0, mx_bit
            if (obsolete_bit(j)) msg(34+(j+1)*2:34+(j+1)*2) = '-'
          end do
        end if
        l = len_trim(msg)
        if (associated(set% var)) write(msg(l+3:),'(F13.7)') set%var(i)
        if (associated(set% oe_str)) write(msg(l+3:),'(A)') trim(set% oe_str(i))
        write(iunit,*) repeat(' ',ind+2)//trim(msg)
      end do
    end if

  end subroutine print_rad_set
!------------------------------------------------------------------------------
  subroutine print_radv(r, hint, header, indent, unit)
    type(t_radv),     intent(in)           :: r
    integer,          intent(in), optional :: hint
    character(len=*), intent(in), optional :: header
    integer,          intent(in), optional :: indent
    integer,          intent(in), optional :: unit

!   character(len=60) :: msg
    integer :: ind, iunit

    if (present(unit)) then
      iunit = unit
    else
      iunit = stdout
    end if

    if (present(indent)) then
       ind = indent
    else
       ind = 0
    end if
    if (present(header)) then
      write(iunit,*) repeat(' ',ind)//trim(header)
    elseif (present(hint)) then
      write(iunit,*) repeat(' ',ind)//'t_radv ',hint
    else
      write(iunit,*) repeat(' ',ind)//'t_radv'
    end if
    call print_rad_set(r%i, hint, header, indent+2, unit)
    write(iunit,*) repeat(' ',ind+2)//'filename     =', r%filename
    write(iunit,*) repeat(' ',ind+2)//'file_id      =', r%file_id
    write(iunit,*) repeat(' ',ind+2)//'model_date   =', r%model_date
    write(iunit,*) repeat(' ',ind+2)//'n_rec        =', r%n_rec
    write(iunit,*) repeat(' ',ind+2)//'n_lev        =', r%n_lev
    write(iunit,*) repeat(' ',ind+2)//'p_unit       =', r%p_unit

#define PR_ARR(var) if (associated(r%var)) then; call print_var('var',shape(r%var)); else; call print_var('var'); endif
    PR_ARR(pe)
    PR_ARR(j_box)
    PR_ARR(i_reprt)
    PR_ARR(i_body)
    PR_ARR(obsnum)
    PR_ARR(date)
    PR_ARR(time)
    PR_ARR(date_d)
    PR_ARR(time_d)
    PR_ARR(ntstep)
    PR_ARR(dlat)
    PR_ARR(dlon)
    PR_ARR(fov)
    PR_ARR(scanl)
    PR_ARR(stzen)
    PR_ARR(stazi)
    PR_ARR(sunzen)
    PR_ARR(sunazi)
    PR_ARR(orb_phase)
    PR_ARR(instr_temp)
    PR_ARR(landfr)
    PR_ARR(stype)
    PR_ARR(shgt)
    PR_ARR(center)
    PR_ARR(subcenter)
    PR_ARR(pred)
    PR_ARR(tr)
    PR_ARR(mdlsfc)
    PR_ARR(op_na)
    PR_ARR(not_rej)
    PR_ARR(cloudy)
    PR_ARR(state)
    PR_ARR(flags)
    PR_ARR(valid)
    PR_ARR(sinfl)
    PR_ARR(plevel)
    PR_ARR(pp_flags)
    PR_ARR(bt_fg)
    PR_ARR(bt_obs)
    PR_ARR(bt_bcor)
    PR_ARR(bt_fg_cs)
    PR_ARR(rad_fg)
    PR_ARR(rad_fg_cs)
    PR_ARR(bcor_)
    PR_ARR(emiss)
    PR_ARR(emis_pc)
    PR_ARR(e_fg)
    PR_ARR(e_obs)
    PR_ARR(p)
    PR_ARR(hh)
    PR_ARR(h_fg)
    PR_ARR(t_fg)
    PR_ARR(q_fg)
    PR_ARR(t2m)
    PR_ARR(q2m)
    PR_ARR(ps_fg)
    PR_ARR(gp_fg)
    PR_ARR(ts_fg)
    PR_ARR(tsm_fg)
    PR_ARR(u10_fg)
    PR_ARR(v10_fg)
    PR_ARR(v10_abs_fg)
    PR_ARR(snw_frc)
    PR_ARR(cld_top)
    PR_ARR(cfraction)
    PR_ARR(clwde)
    PR_ARR(icede)
    PR_ARR(cfrac)
    PR_ARR(cld_fg)
    PR_ARR(cld_frc)
    PR_ARR(clw)
    PR_ARR(trg)
    PR_ARR(H_t)
    PR_ARR(H_q)
    PR_ARR(H_ts)
    PR_ARR(H_ps)
    PR_ARR(H_em_pc)
    PR_ARR(K_t)
    PR_ARR(K_q)
    PR_ARR(K_ts)
    PR_ARR(K_ps)
    PR_ARR(B_tt)
    PR_ARR(B_qq)
    PR_ARR(B_tq)
    PR_ARR(R  )
    PR_ARR(HBHR)
    PR_ARR(t_eb)
    PR_ARR(q_eb)
    PR_ARR(ps_eb)
    PR_ARR(ts_eb)
    PR_ARR(cld_top_eb)
    PR_ARR(cld_frc_eb)
    PR_ARR(snw_frc_eb)
    PR_ARR(emis_pc_eb)
    PR_ARR(im_cl_mean)
    PR_ARR(im_cl_stdv)
    PR_ARR(im_cl_frac)
    PR_ARR(im_chan)

  contains

    subroutine print_var(name, sh)
      character(len=*), intent(in) :: name
      integer,          intent(in), optional :: sh(:)
      character(len=120) :: line

      if (present(sh)) then
        write(line,'(A,T14,": ",10(1x,I9))') trim(name), sh
      else
        write(line,'(A,T14,": ",A)') trim(name), 'not allocated'
      end if
      write(iunit,'(1x,A,A)') repeat(' ',ind+2), trim(line)

    end subroutine print_var

  end subroutine print_radv
!------------------------------------------------------------------------------

  !> Read a TOVS_OBS_CHAN namelist
  subroutine read_tovs_obs_chan_nml(nnml, s, status, i_read)
    integer,         intent(in)            :: nnml   !< unit for reading the namelist
                                                     !! it is assumed that the connected
                                                     !! file is positioned at a namelist
    type(t_rad_set), intent(out)           :: s      !< The set that is to be filled
    integer,         intent(in),  optional :: i_read !< Number of the namelist (for
                                                     !! error output.
    integer,         intent(out), optional :: status !< Exit status
                                                     !! ==0 :  okay
                                                     !! > 0 : non-critical error
                                                     !! < 0 : critical error

    character(len=40), parameter :: proc = 'read_tovs_obs_chan_nml'
    character(len=80)    :: hint
    integer              :: instr_list(m_instr)
    integer              :: i, j, k, n, stat
    integer              :: n_chan
    integer              :: n_instr
    integer              :: ioff
    integer              :: itype
    logical              :: mask(m_chan)
    logical              :: l_stat
    logical              :: l_new

    !----------------------------
    ! Namelist TOVS_OBS_CHAN
    !----------------------------
    character(len=8)     :: c_satellite          ! satellite name
    integer              :: satid                ! WMO satellite ID
    integer              :: grid                 ! RTTOV ID of target grid
    integer              :: flag_instr           ! sensor selection
    integer              :: select_sens          ! sensor selection
    real(sp)             :: max_tskin_cld_cov
    character(len=10)    :: chan_cld_cov
    ! --- general options which will be stored in type t_rad_gopt ---
    integer              :: n_max_calc_y         ! maximum number of simultaneous forward  calculations
    integer              :: n_max_calc_k         ! maximum number of simultaneous K-matrix calculations
    integer              :: n_max_prof_y         ! maximum number of profiles in simultaneous forward  calculations
    integer              :: n_max_prof_k         ! maximum number of profiles in simultaneous K-matrix calculations
                                                 ! NOTE: n_max_prof affects the amount of memory required by 3dvar-specific array
                                                 !       n_max_calc affects the amount of memory required by RTTOV
    integer              :: quality_mode         ! Option for calculation of artificial quality (useful for thinning)
    integer              :: pp_flags_instr       ! instrument index for reading PP_FLAGS
    integer              :: pp_flags             ! PP_FLAGS to be used for blacklisting
    logical              :: use_amsub_89ghz      ! use AMSU-B 89GHz channel instead of AMSU-A 89GHz in AMSUA-A cloud check
    logical              :: l_max_surf_hgt       ! Read MAX_SURFACE_HEIGHT instead of SURFACE_HEIGHT
    integer              :: thin_superob_box_lines ! Superobbing box extent in line direction
    integer              :: thin_superob_box_fovs  ! Superobbing box extent in fov direction
    integer              :: thin_superob_mode    ! 1 for thinning, 2 for superob
    real(wp)             :: max_size_sobox       ! maximal  side length of superob-box in km
    integer              :: max_timediff_sobox   ! maximal observation time difference of pixels within a superob-box
    integer              :: lev_mode             ! provisional

    ! --- instrument specific options which will be stored in array of type t_rad_iopt ---
    integer              :: l2c_type       (m_instr) ! Option for level to channel assignment
    real(wp)             :: l2c_rel_lim    (m_instr) ! relative limit used in level to channel assignment
    real(wp)             :: l2c_abs_lim    (m_instr) ! relative limit used in level to channel assignment
    logical              :: l2c_use_rad    (m_instr) ! use radiances or BTs in level to channel assignment
    real(wp)             :: l2c_max        (m_instr) ! maximum allowed level in l2c check (default)
    real(wp)             :: l2c_max_trop   (m_instr) ! maximum allowed level in l2c check tropics
    real(wp)             :: l2c_max_midlat (m_instr) ! maximum allowed level in l2c check midlatitudes
    real(wp)             :: l2c_max_polar  (m_instr) ! maximum allowed level in l2c check polar regions
    integer              :: surf_infl_mode (m_instr,m_bd) ! surface influence test mode
    real(wp)             :: max_surf_infl  (m_instr,m_bd) ! Maximum influence by surface
    real(wp)             :: d_stemp        (m_instr,m_bd) ! delta(t_skin)
    real(wp)             :: d_emiss        (m_instr,m_bd) ! delta(emissivity)
    integer              :: cloud_mode     (m_instr) ! provisional
    integer              :: rad_out        (m_instr) ! provisional
    logical              :: do_lambertian  (m_instr) ! provisional
    real(wp)             :: specularity    (0:n_styp-1)  ! provisional
    real(wp)             :: specularity_snow             ! provisional
    integer              :: use_o3         (m_instr) ! use o3 in RTTOV
    integer              :: use_co2        (m_instr) ! use co2 in RTTOV
    integer              :: rt_mw_emis_mod (m_instr) ! emissivity model version
    ! --- channel specific options ---
    integer              :: ichan     (m_chan)   ! channel numbers
    integer              :: instr     (m_chan)   ! RTTOV IDs of instruments
    integer              :: flag      (m_chan)   ! channel usage flags
    integer              :: band      (m_chan)   ! channel band
    character(len=m_oe)  :: oe_str    (m_chan)   ! obserror parameterization
    ! --- colocated imager options
    integer              :: use_imager           ! use colocated imager data
    integer              :: im_atlas_id          ! emissivity atlas
    integer              :: im_tskin_opt(3)      ! (instr,chan,stype_bits) for tskin retrieval
    real(wp)             :: im_tskin_cfrac       ! cloud cover threshold for imager tskin retr.

    l_stat = present(status)
    if (l_stat) status = 0

    if (present(i_read)) then
      write(hint, '(" number ",I2)') i_read
    else
      hint = ''
    end if

    call init
    ichan(:)        =    0  ! channel numbers
    instr(:)        =   -1  ! instrument default
    flag(:)         =   -1  ! channel usage flag
    band(:)         =   -1  ! band default
    oe_str(:)       =   ''  ! obserror
!    call construct(oe_par(:))

#define ERR(STAT) if (l_stat) then; status=STAT; return; else; call finish(trim(proc),trim(msg)); endif

    !----------------------------------------------------------------
    ! Try 3 different variants of namelist input one after the other.
    ! Further tryals if error is idicated by iostat return variable.
    ! Some errors are not properly detected by some compilers,
    ! thus the sequence of tests depends on the compiler.
    !
    ! Variants:
    ! read_nml_type_1: chan = real real real        real        real
    ! read_nml_type_2: chan = real real char(len=*) real        real
    ! read_nml_type_3: chan = real real char(len=*) char(len=*) real
    ! hybrid type_2/3: chan = real real char(len=*) 'real'      real
    !----------------------------------------------------------------
#ifdef __SX__
    call read_nml_type_3
    if (stat/=0) then
      call back_nml(nnml, stat)
      if (stat==0) then
        call read_nml_type_2
        if (stat/=0) then
          call back_nml(nnml, stat)
          if (stat==0) call read_nml_type_1
#elif  defined (__GFORTRAN__) || defined (_CRAYFTN)
! CRAY compiler:
! 1.) A hybrid type_2/3 namelist is read in as type_2 but the fourth entry 'real' is
!     interpreted wrongly (always equal to 0.601E-153).
!     Therefore read_nml_type_3 is called before read_nml_type_2
! 2.) The __SX__ version above works correctly as well
!
! GNU compiler:
! 1.) A type_1 namelist is accepted by read_nml_type_3. This should not happen, since the
!     numbers read in for entry 3 are misinterpreted later.
!     Therefore read_nml_type_1 should be called before read_nml_type_3
! 2.) A hybrid type_2/3 namelist is read correctly as type_3.
! 3.) The #else version below works as well.
    call read_nml_type_1
    if (stat/=0) then
      call back_nml(nnml, stat)
      if (stat==0) then
        call read_nml_type_3
        if (stat/=0) then
          call back_nml(nnml, stat)
          if (stat==0) call read_nml_type_2
#else
! INTEL compiler:
! 1.) A type_1 namelist is accepted by read_nml_type_3. This should not happen, since the
!     numbers read in for entry 3 are misinterpreted later.
!     Therefore read_nml_type_1 should be called before read_nml_type_3
! 2.) Namelist type_2 is accepted by read_nml_type_3 and interpreted correctly.
! 3.) The GNU/CRAY version above works as well.
    call read_nml_type_1
    if (stat/=0) then
      call back_nml(nnml, stat)
      if (stat==0) then
        call read_nml_type_2
        if (stat/=0) then
          call back_nml(nnml, stat)
          if (stat==0) call read_nml_type_3
#endif
          if (stat/=0) then
            write(msg, *) 'Failed to read TOVS_OBS_CHAN namelist'//trim(hint),&
                 ' stat=', stat
            call warning(msg)
            call destruct(s)
            ERR(1)
          else
            print*,'namelist type: ',itype
          end if
        end if
      end if
    end if

    if (satid <= 0) satid = sat_name2id(c_satellite)
    if (satid <= 0) then
      msg = 'Unknown satellite "'//trim(c_satellite)//'" in &
           &TOVS_OBS_CHAN namelist'//trim(hint)
      call warning(trim(msg))
      call destruct(s)
      ERR(2)
    endif

    !------------------------------------------
    ! default grid: HIRS(0), AMSU-A(3) for AQUA
    !------------------------------------------
    if (grid < 0) then
      if (satid == 784) then
        grid = 3
      else
        grid = 0
      end if
    endif

    ! Count valid channels
    n_chan = count(ichan(:) > 0)
    if (n_chan == 0) then
      msg = 'No valid channels in TOVS_OBS_CHAN namelist'//trim(hint)
      call warning(trim(msg))
      call destruct(s)
      ERR(3)
    end if

    ! Count instruments
    n_instr = 0
    do i = 1, size(instr)
      if (ichan(i) > 0) then
        l_new = (n_instr == 0)
        if (.not.l_new) then
          l_new = .not.any(instr_list(1:n_instr) == instr(i))
        endif
        if (l_new) then
          if (n_instr < m_instr) then
            n_instr = n_instr + 1
            instr_list(n_instr) = instr(i)
          else
            msg = 'Too many instruments in TOVS_OBS_CHAN namelist'//trim(hint)
            call warning(trim(msg))
            call destruct(s)
            ERR(4)
          end if
        end if
      end if
    end do
    ! Check validity of flag_instr
    if (flag_instr > 0) then
      k = -1
      do j = 1, n_instr
        if (instr_list(j) == flag_instr) then
          k = j
          exit
        end if
      end do
      if (k < 0) then
        msg = 'Invalid flag_instr in TOVS_OBS_CHAN namelist'//trim(hint)
        call warning(trim(msg))
        call destruct(s)
        ERR(5)
      end if
    elseif (flag_instr < -3) then
      msg = 'Invalid flag_instr in TOVS_OBS_CHAN namelist'//trim(hint)
      call warning(trim(msg))
      call destruct(s)
      ERR(6)
    end if
    ! Check validity of max_tskin_cld_cov
    if (max_tskin_cld_cov < 0 .and. max_tskin_cld_cov > 100) then
      msg = 'Invalid max_tskin_cld_cov in TOVS_OBS_CHAN namelist, set to 25'//trim(hint)
      call warning(trim(msg))
      max_tskin_cld_cov = 0.0_sp
    end if
    ! Check validity of grid
    if (grid == -1) grid = 0 ! default HIRS
    k = -1
    do j = 1, n_instr
      if (instr_list(j) == grid) then
        k = j
        exit
      end if
    end do
    if (k < 0) then
      msg = 'Invalid grid in TOVS_OBS_CHAN namelist'//trim(hint)
      call warning(trim(msg))
      call destruct(s)
      ERR(7)
    end if
    ! prepare l2c check
    do j=1, m_instr
       if (l2c_max_trop(j) < 0._wp) then
          if (l2c_max(j) > 0._wp) then
             l2c_max_trop(j) = l2c_max(j)
          else
             l2c_max_trop(j) = l2c_max_trop_def
          end if
       end if
       if (l2c_max_midlat(j) < 0._wp) then
          if (l2c_max(j) > 0._wp) then
             l2c_max_midlat(j) = l2c_max(j)
          else
             l2c_max_midlat(j) = l2c_max_midlat_def
          end if
       end if
       if (l2c_max_polar(j) < 0._wp) then
          if (l2c_max(j) > 0._wp) then
             l2c_max_polar(j) = l2c_max(j)
          else
             l2c_max_polar(j) = l2c_max_polar_def
          end if
       end if
    end do

    s% satid                     = satid
    s% grid                      = grid
    s% flag_instr                = flag_instr
    s% select_sens               = select_sens
    s% n_instr                   = n_instr
    s% n_chan                    = n_chan
    s% instr                     = instr_list
    s% max_tskin_cld_cov         = max_tskin_cld_cov
    s% chan_cld_cov              = chan_cld_cov
    ! general options
    s% gopts% n_max_calc_y       = n_max_calc_y
    s% gopts% n_max_calc_k       = n_max_calc_k
    s% gopts% n_max_prof_y       = n_max_prof_y
    s% gopts% n_max_prof_k       = n_max_prof_k
    s% gopts% quality_mode       = quality_mode
    s% gopts% pp_flags_instr     = pp_flags_instr
    s% gopts% pp_flags           = pp_flags
    s% gopts% use_amsub_89ghz    = use_amsub_89ghz
    s% gopts% l_max_surf_hgt     = l_max_surf_hgt
    ! superobbing options
    s% gopts% thin_superob_box_lines  = thin_superob_box_lines
    s% gopts% thin_superob_box_fovs   = thin_superob_box_fovs
    s% gopts% max_size_sobox     = max_size_sobox  ! maximal size in km of pixels in a superob box (security check)
    s% gopts% max_timediff_sobox = max_timediff_sobox  ! in min
    s% gopts% lev_mode           = lev_mode
    s% gopts% thin_superob_mode  = thin_superob_mode
    s% gopts% use_imager         = use_imager
    s% gopts% im_atlas_id        = im_atlas_id
    s% gopts% im_tskin_opt       = im_tskin_opt
    s% gopts% im_tskin_cfrac     = im_tskin_cfrac

    ! instrument specific options, arrays
    s% iopts% l2c_rel_lim    = l2c_rel_lim
    s% iopts% l2c_abs_lim    = l2c_abs_lim
    s% iopts% l2c_use_rad    = l2c_use_rad
    s% iopts% l2c_max        = l2c_max
    s% iopts% l2c_max_trop   = l2c_max_trop
    s% iopts% l2c_max_midlat = l2c_max_midlat
    s% iopts% l2c_max_polar  = l2c_max_polar
    s% iopts% cloud_mode     = cloud_mode
    s% iopts% rad_out        = rad_out
    s% iopts% do_lambertian  = do_lambertian
    s% iopts% use_o3         = use_o3
    s% iopts% use_co2        = use_co2
    s% iopts% rt_mw_emis_mod = rt_mw_emis_mod

    call satid_bufr2rttov(s% satid, s% rttov_satid, platform=s% platform)
    s% instr_wmo = instr_rttov (s% instr, s% satid)

    allocate(s%chan  (1:n_chan))
    allocate(s%band  (1:n_chan))
    allocate(s%flag  (1:n_chan))
    allocate(s%oe_str(1:n_chan))
!    allocate(s%oe_par(1:n_chan))

    loop_instr: do i = 1, n_instr

      s% iopts(i)% surf_infl_mode(:) = surf_infl_mode(i,:)
      s% iopts(i)% max_surf_infl (:) = max_surf_infl (i,:)
      s% iopts(i)% d_stemp       (:) = d_stemp       (i,:)
      s% iopts(i)% d_emiss       (:) = d_emiss       (i,:)
      s% iopts(i)% specularity   (:) = specularity   (  :)
      s% iopts(i)% specularity_snow  = specularity_snow

      ! Determine namelist entries that belong to this instrument
      mask(:) = (ichan(:) > 0) .and. (instr(:) == s%instr(i))

      ! Check for valid obserror specification
!      if (any(mask(:) .and. oe_par(:)%n(1) < 0)) then
      if (any(mask(:) .and. oe_str(:) == '')) then
        msg = 'There are channels with invalid obserror in namelist &
             &TOVS_OBS_CHAN namelist'//trim(hint)
        call warning(trim(msg))
        call destruct(s)
        ERR(8)
      end if

      ! Check for valid band
      if (hss_instr(s%instr(i))) then
        if (any(mask(:) .and. band(:) <= 0)) then
          call warning('There are channels with invalid band in namelist &
               &TOVS_OBS_CHAN namelist number'//trim(hint))
          call destruct(s)
          ERR(9)
        end if
        if (l2c_type(i) == -999) then
          ! Set default for l2c_type
          s% iopts(i)% l2c_type = l2c_type_def
        else
          s% iopts(i)% l2c_type = 0
        end if
      else
        s% iopts(i)% l2c_type = 0
      end if
      if (s%iopts(i)% l2c_type > 0) s%gopts%opt_vars = ibset(s%gopts%opt_vars, OPTV_L2C)

      if (s%iopts(i)%cloud_mode >= 2) s%gopts%lev_mode = 1

      n = count(mask(:))
      if (i > 1) then
        ioff = s%o_ch_i(i-1) + s%n_ch_i(i-1)
      else
        ioff = 0
      end if
      s% n_ch_i(i) = n
      s% o_ch_i(i) = ioff

      s% chan  (ioff+1:ioff+n) = pack(ichan  (:), mask(:))
      s% band  (ioff+1:ioff+n) = pack(band   (:), mask(:))
      s% flag  (ioff+1:ioff+n) = pack(flag   (:), mask(:))
      s% oe_str(ioff+1:ioff+n) = pack(oe_str (:), mask(:))
!      s% oe_par(ioff+1:ioff+n) = pack(oe_par (:), mask(:))

    end do loop_instr

    s% n_sens = s%n_instr ! will be overwritten by link_rad in var3d

    do i = 1, s%n_chan
      s%oe_str(i) = trim(sat_poly2fparse(s%oe_str(i)))
    end do

    s% source = 'nml'

#undef ERR

  contains

    subroutine read_nml_type_1
    !-------------------------------------------------------------
    ! variant 1 of namelist input: chan = real real real real real
    !-------------------------------------------------------------
      real(kind=sp) :: chan(5, m_chan)  ! instr, chan, flag, var, band
      integer       :: i

      namelist /TOVS_OBS_CHAN/ c_satellite,        &! satellite name
                               satid,              &! WMO satellite ID
                               grid,               &! ID (RTTOV) of target grid
                               flag_instr,         &! ID (RTTOV) of instr. to get flags from
                               select_sens,        &! sensor to get flags from
                               n_max_calc_y,       &! maximum number of simultaneous forward  calculations
                               n_max_calc_k,       &! maximum number of simultaneous K-matrix calculations
                               n_max_prof_y,       &! maximum number of profiles in simultaneous forward  calculations
                               n_max_prof_k,       &! maximum number of profiles in simultaneous K-matrix calculations
                                                    ! NOTE: n_max_prof affects the amount of memory required by 3dvar-specific arrays
                                                    !       n_max_calc affects the amount of memory required by RTTOV
                               quality_mode,       &! Option for calculation of artificial quality (useful for thinning)
                               l2c_type,           &! Option for level to channel assignment
                               l2c_rel_lim,        &! relative limit used in level to channel assignment
                               l2c_abs_lim,        &! absolute limit used in level to channel assignment
                               l2c_use_rad,        &! use radiances or BTs in level to channel assignment
                               l2c_max,            &! maximum allowed level in l2c check (default)
                               l2c_max_trop,       &! maximum allowed level in l2c check tropics
                               l2c_max_midlat,     &! maximum allowed level in l2c check midlatitudes
                               l2c_max_polar,      &! maximum allowed level in l2c check polar regions
                               surf_infl_mode,     &! surface influence test mode
                               max_surf_infl,      &! Maximum influence by surface
                               d_stemp,            &! delta(t_skin)
                               d_emiss,            &! delta(emissivity)
                               l_max_surf_hgt,     &! Read MAX_SURFACE_HEIGHT instead of SURFACE_HEIGHT
                               pp_flags_instr,     &!
                               pp_flags,           &!
                               use_amsub_89ghz,    &!
                               cloud_mode,         &!
                               rad_out,            &!
                               do_lambertian,      &!
                               specularity,        &!
                               specularity_snow,   &!
                               lev_mode,           &!
                               use_o3,             &!
                               use_co2,            &!
                               use_imager,         &!
                               im_atlas_id,        &!
                               im_tskin_opt,       &!
                               im_tskin_cfrac,     &!
                               rt_mw_emis_mod,     &!
                               thin_superob_box_lines,  &!
                               thin_superob_box_fovs,   &!
                               thin_superob_mode,  &!
                               max_size_sobox,     &! maximum dist. between pixels in superob. box
                               max_timediff_sobox, &! maximum allowed time difference of superobbed pixels
                               chan,               &! instr, chan, flag, var, band
                               max_tskin_cld_cov,  &! allowed cloud fraction in tskin retrieval
                               chan_cld_cov         ! select channel for multi-channel cloud_cover input

      call init
      chan(1,:)       = -1._sp
      chan(2,:)       = -1._sp
      chan(3,:)       = -1._sp
      chan(4,:)       = -99._sp
      chan(5,:)       = -1._sp

      read (nnml ,nml=TOVS_OBS_CHAN, iostat=stat)
      if (stat /= 0) return

      instr(:) = nint(chan(1,:))
      ichan(:) = nint(chan(2,:))
      flag (:) = nint(chan(3,:))
      band (:) = nint(chan(5,:))
      do i = 1, size(chan,2)
        write(oe_str(i),'(e14.7)') chan(4,i)
      end do

      itype = 1

    end subroutine read_nml_type_1


    subroutine read_nml_type_2
    !--------------------------------------------------------------------
    ! variant 2 of namelist input: chan = real real char(len=*) real real
    !--------------------------------------------------------------------

      type t_chan
         real(wp)           :: instr       = -1._wp
         real(wp)           :: ichan       = -1._wp
         character(len=300) :: flag        = ''
         real(wp)           :: var         = -99._wp
         real(wp)           :: band        = -1._wp
      end type t_chan

      type(t_chan)      :: chan(m_chan)
      integer           :: i

      namelist /TOVS_OBS_CHAN/ c_satellite,        &! satellite name
                               satid,              &! WMO satellite ID
                               grid,               &! ID (RTTOV) of target grid
                               flag_instr,         &! ID (RTTOV) of instr. to get flags from
                               select_sens,        &! sensor to get flags from
                               n_max_calc_y,       &! maximum number of simultaneous forward  calculations
                               n_max_calc_k,       &! maximum number of simultaneous K-matrix calculations
                               n_max_prof_y,       &! maximum number of profiles in simultaneous forward  calculations
                               n_max_prof_k,       &! maximum number of profiles in simultaneous K-matrix calculations
                                                    ! NOTE: n_max_prof affects the amount of memory required by 3dvar-specific arrays
                                                    !       n_max_calc affects the amount of memory required by RTTOV
                               quality_mode,       &! Option for calculation of artificial quality (useful for thinning)
                               l2c_type,           &! Option for level to channel assignment
                               l2c_rel_lim,        &! relative limit used in level to channel assignment
                               l2c_abs_lim,        &! absolute limit used in level to channel assignment
                               l2c_use_rad,        &! use radiances or BTs in level to channel assignment
                               l2c_max,            &! maximum allowed level in l2c check (default)
                               l2c_max_trop,       &! maximum allowed level in l2c check tropics
                               l2c_max_midlat,     &! maximum allowed level in l2c check midlatitudes
                               l2c_max_polar,      &! maximum allowed level in l2c check polar regions
                               surf_infl_mode,     &! surface influence test mode
                               max_surf_infl,      &! Maximum influence by surface
                               d_stemp,            &! delta(t_skin)
                               d_emiss,            &! delta(emissivity)
                               l_max_surf_hgt,     &! Read MAX_SURFACE_HEIGHT instead of SURFACE_HEIGHT
                               pp_flags_instr,     &!
                               pp_flags,           &!
                               use_amsub_89ghz,    &!
                               cloud_mode,         &!
                               rad_out,            &!
                               do_lambertian,      &!
                               specularity,        &!
                               specularity_snow,   &!
                               lev_mode,           &!
                               use_o3,             &!
                               use_co2,            &!
                               use_imager,         &!
                               im_atlas_id,        &!
                               im_tskin_opt,       &!
                               im_tskin_cfrac,     &!
                               rt_mw_emis_mod,     &!
                               thin_superob_box_lines,  &!
                               thin_superob_box_fovs,   &!
                               thin_superob_mode,  &
                               max_size_sobox,     &! maximum side length of superobbing box in km
                               max_timediff_sobox, &! maximum allowed time difference of superobbed pixels
                               chan,               &! new convention
                               max_tskin_cld_cov,  &! allowed cloud fraction in tskin retrieval
                               chan_cld_cov         ! select channel for multi-channel cloud_cover input

      call init
      chan(:)         = t_chan(-1._wp,-1._wp,'',-99._wp,-1._wp)

      read (nnml ,nml=TOVS_OBS_CHAN, iostat=stat)
      if (stat /= 0) return

      chan_loop: do i = 1, size(chan)
         if (chan(i)% ichan == 0) cycle

         instr(i) = chan(i)% instr
         ichan(i) = chan(i)% ichan
         flag (i) = get_flags(chan(i)% flag)
         band (i) = chan(i)% band
         write(oe_str(i),'(e14.7)') chan(i)% var

      end do chan_loop

      itype = 2

    end subroutine read_nml_type_2

    subroutine read_nml_type_3
    !---------------------------------------------------------------------------
    ! variant 3 of namelist input: chan = real real char(len=*) char(len=*) real
    !---------------------------------------------------------------------------

      type t_chan
         real(wp)            :: instr       = 0._wp
         real(wp)            :: ichan       = 0._wp
         character(len=300)  :: flag        = ''
         character(len=m_oe) :: obserr_par  = '-99.'
         real(wp)            :: band        = -1._wp
      end type t_chan

      type(t_chan)      :: chan(m_chan)
      integer           :: i

      namelist /TOVS_OBS_CHAN/ c_satellite,        &! satellite name
                               satid,              &! WMO satellite ID
                               grid,               &! ID (RTTOV) of target grid
                               flag_instr,         &! ID (RTTOV) of instr. to get flags from
                               select_sens,        &! sensor to get flags from
                               n_max_calc_y,       &! maximum number of simultaneous forward  calculations
                               n_max_calc_k,       &! maximum number of simultaneous K-matrix calculations
                               n_max_prof_y,       &! maximum number of profiles in simultaneous forward  calculations
                               n_max_prof_k,       &! maximum number of profiles in simultaneous K-matrix calculations
                                                    ! NOTE: n_max_prof affects the amount of memory required by 3dvar-specific arrays
                                                    !       n_max_calc affects the amount of memory required by RTTOV
                               quality_mode,       &! Option for calculation of artificial quality (useful for thinning)
                               l2c_type,           &! Option for level to channel assignment
                               l2c_rel_lim,        &! relative limit used in level to channel assignment
                               l2c_abs_lim,        &! absolute limit used in level to channel assignment
                               l2c_use_rad,        &! use radiances or BTs in level to channel assignment
                               l2c_max,            &! maximum allowed level in l2c check (default)
                               l2c_max_trop,       &! maximum allowed level in l2c check tropics
                               l2c_max_midlat,     &! maximum allowed level in l2c check midlatitudes
                               l2c_max_polar,      &! maximum allowed level in l2c check polar regions
                               surf_infl_mode,     &! surface influence test mode
                               max_surf_infl,      &! Maximum influence by surface
                               d_stemp,            &! delta(t_skin)
                               d_emiss,            &! delta(emissivity)
                               l_max_surf_hgt,     &! Read MAX_SURFACE_HEIGHT instead of SURFACE_HEIGHT
                               pp_flags_instr,     &!
                               pp_flags,           &!
                               use_amsub_89ghz,    &!
                               cloud_mode,         &!
                               rad_out,            &!
                               do_lambertian,      &!
                               specularity,        &!
                               specularity_snow,   &!
                               lev_mode,           &!
                               use_o3,             &!
                               use_co2,            &!
                               use_imager,         &!
                               im_atlas_id,        &!
                               im_tskin_opt,       &!
                               im_tskin_cfrac,     &!
                               rt_mw_emis_mod,     &!
                               thin_superob_box_lines,  &!
                               thin_superob_box_fovs,   &!
                               thin_superob_mode,  &!
                               max_size_sobox,     &! maximum side length of superobbing box in km
                               max_timediff_sobox, &! maximum allowed time difference of superobbed pixels
                               chan,               &! new convention
                               max_tskin_cld_cov,  &! allowed cloud fraction in tskin retrieval
                               chan_cld_cov         ! select channel for multi-channel cloud_cover input

      call init
      chan(:)         = t_chan(-1._wp,-1._wp,'','',-1._wp)

      read (nnml ,nml=TOVS_OBS_CHAN, iostat=stat)
      if (stat /= 0) return

      chan_loop: do i = 1, size(chan)
         if (chan(i)% ichan == 0) cycle

         instr(i)  = chan(i)% instr
         ichan(i)  = chan(i)% ichan
         flag (i)  = get_flags(chan(i)% flag)
         oe_str(i) = trim(chan(i)% obserr_par)
         band (i)  = chan(i)% band

      end do chan_loop

      itype = 3

    end subroutine read_nml_type_3

    integer function get_flags(str)
      character(len=*), intent(in) :: str

      integer           :: j, k, istat, n_str, m, bit_start, bit_end
      character(len=20) :: substr(32)
      logical           :: l_flag
      get_flags = 0

      call split2(str, sep=' ,;', array=substr, n=n_str)
      if (n_str >= 1) then
        read(substr(1),*,iostat=istat) l_flag
        if (istat == 0) then
          ! Interpret as a logical array
          do j = 1, n_str
            read(substr(j),*,iostat=istat) l_flag
            if (istat==0 .and. l_flag) get_flags = ibset(get_flags, j-1)
          end do
        else
          ! Interpret as a sequence of bit identifiers
          do j = 1, n_str
            substr(j) = uppercase(substr(j))
            m = scan(substr(j), '-')
            if (m > 0) then
              if (m > 1) then
                bit_start = max(get_bit_num(substr(j)(1:m-1)), 0)
              else
                bit_start = 0
              end if
              if (m < len_trim(substr(j))) then
                bit_end = max(get_bit_num(substr(j)(m+1:)), 0)
              else
                bit_end = 31
              end if
              do k = bit_start, bit_end
                get_flags = ibset(get_flags, k)
              end do
            else
              k = get_bit_num(substr(j))
              if (k >= 0) then
                if (k > 31 .and. n_str == 1) then
                  get_flags = k
                elseif (k <= 31) then
                  get_flags = ibset(get_flags, k)
                else
                  write(0,*) 'k=',k
                  call error('invalid flag (1): "'//trim(substr(j))//'"')
                end if
              elseif (substr(j) == 'L2C') then
                  get_flags = ibset(get_flags, USE_L2C_CLD)
                  get_flags = ibset(get_flags, USE_L2C_SURF)
              else
                write(0,*) 'k=',k
                call error('invalid flag (2): "'//trim(substr(j))//'"')
              end if
            end if
          end do
        end if
      end if
    end function get_flags

    subroutine init
      flag_instr         = -1    ! where to get quality flags from
                                 ! (for mapped, i.e. multi-instrument FOVs only):
                                 ! >   =0: get flag from instrument with this RTTOV ID
                                 !  -1: get flag from grid instrument
                                 !  -2: worst flag of all requested instruments
                                 !      (watch out for meaning of 'worst'   !)
                                 !  -3: as for -1, but flag as 'mismatch'
                                 !      if inconsistent for diff. requested instruments
      select_sens        = -1
      c_satellite        = ''    ! satellite name
      satid              = -1    ! WMO satellite ID
      grid               = -1    ! ID of target grid
      n_max_calc_y       = n_max_calc_y_def
      n_max_calc_k       = n_max_calc_k_def
      n_max_prof_y       = n_max_prof_y_def
      n_max_prof_k       = n_max_prof_k_def
      quality_mode       = quality_mode_def
      l2c_type           = -999
      l2c_rel_lim        = l2c_rel_lim_def
      l2c_abs_lim        = l2c_abs_lim_def
      l2c_use_rad        = l2c_use_rad_def
      l2c_max            = -1._wp
      l2c_max_trop       = -1._wp
      l2c_max_midlat     = -1._wp
      l2c_max_polar      = -1._wp
      surf_infl_mode     = surf_infl_mode_def
      max_surf_infl      = max_surf_infl_def
      d_stemp            = d_stemp_def
      d_emiss            = d_emiss_def
      l_max_surf_hgt     = l_max_surf_hgt_def
      pp_flags_instr     = pp_flags_instr_def
      pp_flags           = pp_flags_def
      use_amsub_89ghz    = use_amsub_89ghz_def
      cloud_mode         = 0
      rad_out            = 0
      do_lambertian      = .false.
      specularity        = 1._wp
      specularity_snow   = 1._wp
      lev_mode           = 0
      use_o3             = 0
      use_co2            = 0
      use_imager         = 0
      im_atlas_id        = -1
      im_tskin_opt       = (/-1,-1,0/)
      im_tskin_cfrac     = 0._wp
      rt_mw_emis_mod     = -1
      thin_superob_box_lines  = 0
      thin_superob_box_fovs   = 0
      thin_superob_mode  = 1
      max_size_sobox     = 50._wp
      max_timediff_sobox = 15
      max_tskin_cld_cov  = 0.0_sp ! allowed cloud fraction in tskin retrieval
      chan_cld_cov       = ''
    end subroutine init

!------------------------------------------------------------------------------

    integer function get_bit_num(str)
      character(len=*), intent(in) :: str
      integer :: k, stat
      get_bit_num = -1
      read(str,*,iostat=stat) k
      if (stat == 0) then
         if (k >= 0 .and. k <=31) get_bit_num = k
      else
         do k = 0, 31
            if (trim(str) == trim(name_use_bit(k))) then
               get_bit_num = k
               return
            end if
         end do
      end if
    end function get_bit_num

    subroutine error(msg)
      character(len=*), intent(in) :: msg
      call warning('*** read_tovs_obs_chan_nml : '//trim(msg), nowarn=.true.)
      if (present(status)) then
        status = 13
      else
        call finish('read_tovs_obs_chan_nml', trim(msg))
      end if
    end subroutine error

  end subroutine read_tovs_obs_chan_nml


  function sat_poly2fparse(str) result(s)
    ! Convert string from old obsolete "sat_poly" format to new fparser input.
    ! This is necessary for backwards compatibilty of TOVS_OBS_CHAN namelists
    character(len=500) :: s
    character(len=*), intent(in) :: str

    character(len=15), parameter :: numbers='0123456789eE-+.'
    character(len=30)            :: pstr(100)
    character(len=11)            :: pow
    integer                      :: i, n, stat
    real(wp)                     :: x

    s = trim(str)
    if (verify(s, numbers//' ') > 0 .or. len_trim(s) == 0) return ! Contains other characters than numbers and blanks

    call split2(s, array=pstr, n=n, status=stat)
    if (stat /= 0) return

    do i = 1, n
      read(pstr(i),*,iostat=stat) x
      if (stat /= 0) then
        ! Not a real number, is not a "sat_poly"-string, most likely a string that should
        ! be passed to mo_fparser or mo_range_fparser
        s = trim(str)
        return
      end if
      select case(i)
      case(1)
        s = trim(pstr(i))
      case(2:)
        if (verify(pstr(i)(1:1), '+-') > 0) s = trim(s)//'+'
        write(pow,'("*sat_zen^",I2)') i-1
        s = trim(s)//trim(pstr(i))//trim(pow)
      end select
    end do

  end function sat_poly2fparse
!------------------------------------------------------------------------------
  !> Go back to the last line beginning with '&'
  subroutine back_nml(iunit, stat)
    integer, intent(in)  :: iunit
    integer, intent(out) :: stat

    character(len=120) :: line
    integer            :: j

    backspace (iunit)
    do j = 1, 10000
       read(iunit,*,iostat=stat) line
       if (stat /= 0) return
       line=adjustl(line)
       backspace (iunit)
       if (line(1:1) == '&') then
          stat = 0
          return
       end if
       backspace (iunit)
    end do
    stat = -1

  end subroutine back_nml

!------------------------------------------------------------------------------

  !> Associate a t_radv structure (e.g. read from a sat.data file) with one of the
  !! t_rad_set structures in rad_set (e.g. read from a namelist).
  !! \todo Add instrument to t_radv, if it is required for the bias correction.
  subroutine link_rad(r, status, i_set)
    type(t_radv), intent(inout), target   :: r
    integer,      intent(out),   optional :: status
    integer,      intent(out),   optional :: i_set

    character(len=8), parameter :: proc = 'link_rad'
    character(len=180)       :: descr = ''
    character(len=180)       :: msg   = ''
    type(t_rad_set), pointer :: s  => NULL()
    type(t_rad_set), pointer :: rs => NULL()
    integer                  :: i, j, k, iset, is, ie, ic
    integer                  :: ioff_r, ioff_s, stat
    integer,         pointer :: n_r, n_s
    logical                  :: l_stat

    l_stat = present(status)
    if (l_stat) status = 0

#define ERR(STAT) call warning(msg); if(l_stat)then; status=STAT; return; else; call finish(proc,trim(msg)); endif

    s => r% i
    write(descr, '(A,I3.3,A,I2.2,A)') '(file "'//trim(r% filename)//&
         &'", satid ',s%satid,', grid ',s%grid,')'

    ! Find the correct rad_set entry for the dataset
    iset = set_indx(rad_set(:), satid=s% satid, grid=s% grid)
    if (present(i_set)) i_set = iset
    if (iset <= 0) then
      msg = 'failed to link radiance data '//trim(descr)//&
           &' to a "rad_set" entry.'
      call warning(trim(msg))
      ERR(1)
    end if

    rs => rad_set(iset)

    ! Compare the instruments
    ! Is there a rad_set entry for each instrument in r?
    do i = 1, s% n_instr
      if (.not.any(rs% instr(1:rs%n_instr) == s%instr(i))) then
        write(msg, '(1x,A,I2.2,A)', iostat=stat) 'instrument ',s%instr(i), &
             &trim(descr)//' is not contained in the corresponding rad_set entry.'
        ERR(2)
      end if
    end do
    ! Is each instrument in rad_set contained in the dataset r?
    i = 1
    do while (i <= rs% n_instr)
      if (.not.any(s% instr(1:rs%n_instr) == rs%instr(i))) then
        write(msg, '(1x,A,I2.2,A)', iostat=stat) 'instrument ',rs%instr(i), &
             &' (in rad_set) is not contained in the corresponding dataset '//trim(descr)
        call warning(msg)
        ioff_r = rs% o_ch_i(i)
        n_r    => rs% n_ch_i(i)
        if (any(btest(rs% flag(ioff_r+1:ioff_r+n_r), USE_BCOR))) then
          write(msg, '(1x,A,I2.2,A)', iostat=stat) 'instrument ',rs%instr(i), &
               &' is required for the bias-correction. Since it is not contained in &
               &the corresponding dataset '//trim(descr)//' the dataset is skipped.'
          ERR(3)
        end if
        ! Delete instrument from the rad_set entry
        call del_instr_rad_set(rs, i) !rs% n_instr is decreased by 1 here
        if (rs% n_instr <= 0) then
          msg = 'no instrument contained in rad_set.'
          ERR(4)
        end if
      else
        i = i + 1
      end if
    end do

    if (rs%n_instr < s%n_instr) then
       ! e.g. GMI: here we should process the instrument as two instruments with the same ID
       ! -> split s (from TOVS_OBS_CHAN namelist) into two instruments (with the same ID)
       is = 1
       ie = is
       do while (ie < s%n_instr)
          do ie = is, s%n_instr
             if (ie >= s%n_instr) EXIT
             if (s%instr(ie+1) /= s%instr(is)) EXIT
          end do
          if (ie > is) then
             ! Found double instrument
             do j = 1, rs%n_instr
                if (rs%instr(j) == s%instr(is)) EXIT
             end do
             j = min(j, rs%n_instr)
             if (rs%instr(j) /= s%instr(is)) then
                write(msg,'("instrument ",I3.3," is not contained in corresponding TOVS_OBS_CHAN namelist ",A)') &
                     s%instr(is),trim(descr)
                ERR(10)
             end if
             ic = s%chan(s%o_ch_i(is) + s%n_ch_i(is))
             do k = 1, rs%n_ch_i(j)
                if (rs%chan(rs%o_ch_i(j) + k) == ic) EXIT
             end do
             k = min(k,rs%n_ch_i(j))
             if (rs%chan(rs%o_ch_i(j) + k) /= ic) then
                write(msg,'("channel ",I4.4," of instrument ",I3.3," is not contained in &
                     &corresponding TOVS_OBS_CHAN namelist ",A)') s%instr(is), ic, trim(descr)
                ERR(11)
             end if
             rs%instr (j+1:rs%n_instr+1) = rs%instr (j  :rs%n_instr)
             rs%o_ch_i(j+2:rs%n_instr+1) = rs%o_ch_i(j+1:rs%n_instr)
             rs%n_ch_i(j+2:rs%n_instr+1) = rs%n_ch_i(j+1:rs%n_instr)
             rs%n_ch_i(j+1) = rs%n_ch_i(j) - k
             rs%n_ch_i(j  ) = k
             rs%o_ch_i(j+1) = rs%o_ch_i(j) + k
             rs%n_instr = rs%n_instr + 1
          end if
          is = ie + 1
       end do
    end if

    if (s%n_instr /= rs%n_instr) then
      msg = 'number of instruments in rad_set and dataset '//trim(descr)//' is not equal.'
      ERR(5)
    end if
    if (any(s% instr(1:s%n_instr) /= rs% instr(1:rs%n_instr))) then
      msg = 'the instruments in rad_set and dataset '//trim(descr)//' are not equal.'
      ERR(6)
    end if

    ! Compare the channels for each instrument. Reorder or delete them if necessary.
    ! For each channel in r there must be an entry in rad_set (because the obs error and
    ! the use-flag is required).
    instr_loop: do i = 1, rs% n_instr
      ioff_r =  rs% o_ch_i(i)
      ioff_s =  s % o_ch_i(i)
      n_r    => rs% n_ch_i(i)
      n_s    => s % n_ch_i(i)
      ! Check whether there is an entry in rad_set for each channel in the dataset
      do j = ioff_s + 1, ioff_s + n_s
        if (.not. any(rs%chan(ioff_r+1:ioff_r+n_r) == s% chan(j))) then
          write(msg, '(1x,A,I4.4,A,I2.2)', iostat=stat) 'channel ', s% chan(j), &
               &', instrument ',s%instr(i)
          msg = trim(msg)//' '//trim(descr)//' is not contained in the &
               &corresponding rad_set entry.'
          ERR(7)
        end if
      end do
      ! Check whether channels in rad_set can be deleted.
      j = ioff_r + 1
      do while (j <= ioff_r + n_r)
        if (.not. any(s%chan(ioff_s+1:ioff_s+n_s) == rs% chan(j))) then
          call del_chan_rad_set(rs, j, i)
        else
          j = j + 1
        end if
      end do
      ! Now the numbers of channels should be equal in rad_set and r.
      if (n_s /= n_r) then
        write(msg, '(1x,A,I2.2)', iostat=stat) 'numbers of channels in &
             &rad_set and dataset '//trim(descr)//' do not agree for instrument ',&
             &s%instr(i)
        ERR(8)
      end if
      ! Reorder channels in rad_set according to the ordering in r.
      do j = 1, n_s
        if (rs% chan(ioff_r+j) /= s% chan(ioff_s+j)) then
          do k = 1, n_r
            if (rs% chan(ioff_r+k) == s% chan(ioff_s+j)) exit
          end do
          call swap_chan(rs, ioff_r+k, ioff_r+j)
        end if
      end do
      ! Check whether a complete channel is to be ignored
      j = 1
      do while (j <= n_r)
        if (rs% flag(ioff_r+j) == 0) then
          call del_chan_rad_set(rs, ioff_r+j, i)
          call del_chan_radv   (r,  ioff_s+j, i)
        else
          j = j + 1
        end if
      end do
      ! evaluste options that might require cldlev (pre-)calculation
      call eval_cldlev_opt(rs, i)
    end do instr_loop
    rs% n_chan = sum(rs% n_ch_i(1:rs% n_instr))
    if (any(rs% chan(1:rs% n_chan) /= s% chan(1:rs% n_chan))) then
      msg = 'some channels in rad_set and dataset '//trim(descr)//' do not agree.'
      ERR(9)
    end if

    ! Indicate that the dataset is to be used furthermore
    s% id  = iset
    rs% id = s% id
    rs% source = s% source

    ! Copy info
    rs% instr_wmo(:) = s% instr_wmo(:)
    rs% grid_wmo     = s% grid_wmo
    rs% n_sens       = s% n_sens
    allocate(rs%sensor_instr(rs%n_sens))
    if (associated(s%sensor_instr)) rs%sensor_instr(:) = s%sensor_instr(:)
    if (associated(s%iflag)) then
      allocate(rs%iflag(rs%n_chan))
      rs% iflag(:) = s% iflag(:)
    end if
    if (any(rs%iopts(1:rs%n_instr)% l2c_type > 0)) &
                                   rs%gopts%opt_vars = ibset(rs%gopts%opt_vars, OPTV_L2C    )
    if (associated(r% orb_phase )) rs%gopts%opt_vars = ibset(rs%gopts%opt_vars, OPTV_ORB_PH )
    if (associated(r% instr_temp)) rs%gopts%opt_vars = ibset(rs%gopts%opt_vars, OPTV_INS_TMP)
    if (associated(r% nwc_flg   )) rs%gopts%opt_vars = ibset(rs%gopts%opt_vars, OPTV_NWC_FLG)
    if (associated(r% cld_frc   )) rs%gopts%opt_vars = ibset(rs%gopts%opt_vars, OPTV_CLD_FRC)
    rs% n_im_chan    = s% n_im_chan
    rs% n_im_cluster = s% n_im_cluster
    rs% im_rttov_id  = s% im_rttov_id
    rs% im_rttov_indx= s% im_rttov_indx
    rs% im_chan      = s% im_chan
    rs% im_chan_indx = s% im_chan_indx
    if (rs%n_im_chan > 0 .and. rs%n_im_cluster > 0 .and. rs%im_rttov_id >= 0) &
         rs%gopts%opt_vars = ibset(rs%gopts%opt_vars, OPTV_CLD_FLG)

    s = rs

#undef ERR

  end subroutine link_rad


  subroutine eval_cldlev_opt(rs,instr)
    use mo_fparse, only : replace_shortcuts, lt
    type(t_rad_set), intent(inout), target :: rs
    integer,         intent(in)            :: instr
    type(t_rad_iopt),pointer :: iopt
    character(len=lt)        :: str
    integer                  :: j
    iopt => rs%iopts(instr)
    iopt%l_cldlev = (iopt%cloud_mode < 0)
    if (.not.iopt%l_cldlev) then
      do j = rs%o_ch_i(instr)+1,rs%o_ch_i(instr)+rs%n_ch_i(instr)
        str = trim(tolower(rs%oe_str(j)))
        call replace_shortcuts(str)
        if (check_word(trim(str),'cldlev') .or. check_word(trim(str),'p_cldlev')) then
          iopt%l_cldlev = .true.
          exit
        end if
      end do
    end if
  end subroutine eval_cldlev_opt


  function check_word(str, ss) result(lf)
    logical :: lf
    character(len=*), intent(in) :: str
    character(len=*), intent(in) :: ss
    integer :: l, ls, k
    lf = .false.
    l  = len_trim(str)
    ls = len_trim(ss)
    k = index(str, ss)
    if (k > 0) then
      lf = .true.
      if (k > 1) then
        if (index('abcdefghijklmnopqrstuvwxyz_0123456789',str(k-1:k-1)) > 0) lf = .false.
      end if
      if (k+ls <= l) then
        if (index('abcdefghijklmnopqrstuvwxyz_0123456789',str(k+ls:k+ls)) > 0) lf = .false.
      end if
    end if
  end function check_word
!------------------------------------------------------------------------------

  subroutine swap_chan(s, i, j)
  type(t_rad_set), intent(inout) :: s
  integer,         intent(in)    :: i, j
    integer                 :: i_aux
    real(sp)                :: s_aux
    real(wp)                :: r_aux
    real(wp), allocatable   :: r_pc_aux(:)
    integer,  allocatable   :: i_mw_aux(:,:)
    character(len=m_oe)     :: oe_aux
    type(t_range_fparse)    :: oe_trf
    if (associated(s%chan )) then
      i_aux            = s% chan(i)
      s% chan(i)       = s% chan(j)
      s% chan(j)       = i_aux
    end if
    if (associated(s%chidx)) then
      i_aux            = s% chidx(i)
      s% chidx(i)      = s% chidx(j)
      s% chidx(j)      = i_aux
    end if
     if (associated(s%iflag)) then
      i_aux            = s% iflag(i)
      s% iflag(i)      = s% iflag(j)
      s% iflag(j)      = i_aux
    end if
    if (associated(s%band )) then
      i_aux            = s% band (i)
      s% band (i)      = s% band (j)
      s% band (j)      = i_aux
    end if
    if (associated(s%flag )) then
      i_aux            = s% flag (i)
      s% flag (i)      = s% flag (j)
      s% flag (j)      = i_aux
    end if
    if (associated(s%wavenum )) then
      r_aux            = s% wavenum (i)
      s% wavenum (i)   = s% wavenum (j)
      s% wavenum (j)   = r_aux
    end if
    if (associated(s%emis_land )) then
      r_aux            = s% emis_land (i)
      s% emis_land (i) = s% emis_land (j)
      s% emis_land (j) = r_aux
    end if
    if (associated(s%emis_snow )) then
      r_aux            = s% emis_snow (i)
      s% emis_snow (i) = s% emis_snow (j)
      s% emis_snow (j) = r_aux
    end if
    if (associated(s%emis_sice )) then
      r_aux            = s% emis_sice (i)
      s% emis_sice (i) = s% emis_sice (j)
      s% emis_sice (j) = r_aux
    end if
    if (associated(s%emis_pc   )) then
      allocate(r_pc_aux(size(s% emis_pc, 1)))
      r_pc_aux   (:)   = s% emis_pc   (:,i)
      s% emis_pc (:,i) = s% emis_pc   (:,j)
      s% emis_pc (:,j) = r_pc_aux     (:)
      deallocate(r_pc_aux)
    end if
    if (associated(s%i1 )) then
      i_aux     = s% i1 (i)
      s% i1 (i) = s% i1 (j)
      s% i1 (j) = i_aux
    end if
    if (associated(s%w1 )) then
      r_aux     = s% w1 (i)
      s% w1 (i) = s% w1 (j)
      s% w1 (j) = r_aux
    end if
    if (associated(s%emis_ind )) then
      allocate(i_mw_aux(size(s% emis_ind,1), 0:size(s% emis_ind,2)-1))
      i_mw_aux    (:,:)  = s% emis_ind (:,:,i)
      s% emis_ind(:,:,i) = s% emis_ind (:,:,j)
      s% emis_ind(:,:,j) = i_mw_aux     (:,:)
      deallocate(i_mw_aux)
    end if
    if (associated(s%tskin_ind )) then
      allocate(i_mw_aux(size(s% tskin_ind,1), 0:size(s% tskin_ind,2)-1))
      i_mw_aux    (:,:)  = s% tskin_ind (:,:,i)
      s% tskin_ind(:,:,i) = s% tskin_ind (:,:,j)
      s% tskin_ind(:,:,j) = i_mw_aux     (:,:)
      deallocate(i_mw_aux)
    end if
    if (associated(s%var  )) then
      s_aux            = s% var  (i)
      s% var  (i)      = s% var  (j)
      s% var  (j)      = s_aux
    end if
    if (associated(s%oe_str)) then
      oe_aux           = s% oe_str(i)
      s% oe_str(i)     = s% oe_str(j)
      s% oe_str(j)     = oe_aux
    end if
    if (associated(s%oe_trf)) then
      oe_trf           = s% oe_trf(i)
      s% oe_trf(i)     = s% oe_trf(j)
      s% oe_trf(j)     = oe_trf
    end if
    if (associated(s%bc% ib_ch )) then
      i_aux            = s% bc% ib_ch (i)
      s% bc% ib_ch (i) = s% bc% ib_ch (j)
      s% bc% ib_ch (j) = i_aux
    end if

  end subroutine swap_chan
!------------------------------------------------------------------------------

  !> Delete an instrument from t_rad_set
  subroutine del_instr_rad_set(s, instr)
    type(t_rad_set), intent(inout) :: s
    integer,         intent(in)    :: instr

    integer,             pointer :: idum  (:)     => NULL()
    integer,             pointer :: idum3 (:,:,:) => NULL()
    real(sp),            pointer :: sdum  (:)     => NULL()
    character(len=m_oe), pointer :: odum  (:)     => NULL()
    type(t_range_fparse),pointer :: trfdum(:)     => NULL()
    real(wp),            pointer :: wdum  (:)     => NULL()
    real(wp),            pointer :: wdum2 (:,:)   => NULL()
    integer                      :: ioff, n, ni, nch

    ni   = s%n_instr
    ioff = s% o_ch_i(instr)
    n    = s% n_ch_i(instr)
    nch  = s% n_chan - n

    if (ni > instr) then
      s% iopts     (instr:ni-1) = s% iopts     (instr+1:ni)
      s% instr     (instr:ni-1) = s% instr     (instr+1:ni)
      s% instr_wmo (instr:ni-1) = s% instr_wmo (instr+1:ni)
      s% n_ch_i    (instr:ni-1) = s% n_ch_i    (instr+1:ni)
      s% o_ch_i    (instr:ni-1) = s% o_ch_i    (instr+1:ni) - n
      s% rttov_indx(instr:ni-1) = s% rttov_indx(instr+1:ni)
    end if
    s%n_instr = s%n_instr - 1

    if ((s% n_instr <= 0) .or. (nch <= 0)) then
      if (associated(s% chan     )) deallocate(s% chan     )
      if (associated(s% chidx    )) deallocate(s% chidx    )
      if (associated(s% iflag    )) deallocate(s% iflag    )
      if (associated(s% band     )) deallocate(s% band     )
      if (associated(s% flag     )) deallocate(s% flag     )
      if (associated(s% wavenum  )) deallocate(s% wavenum  )
      if (associated(s% emis_land)) deallocate(s% emis_land)
      if (associated(s% emis_snow)) deallocate(s% emis_snow)
      if (associated(s% emis_sice)) deallocate(s% emis_sice)
      if (associated(s% emis_pc  )) deallocate(s% emis_pc  )
      if (associated(s% i1       )) deallocate(s% i1       )
      if (associated(s% w1       )) deallocate(s% w1       )
      if (associated(s% emis_ind )) deallocate(s% emis_ind )
      if (associated(s% tskin_ind)) deallocate(s% tskin_ind)
      if (associated(s% var      )) deallocate(s% var      )
      if (associated(s% oe_str   )) deallocate(s% oe_str   )
      if (associated(s% bc% p_tr )) deallocate(s% bc% p_tr )
      if (associated(s% bc% ib_ch)) deallocate(s% bc% ib_ch)
      if (associated(s% bc% type )) deallocate(s% bc% type )
!      call destruct(s% oe_par)
      return
    end if

    if (associated(s% chan )) then
      allocate(idum(nch))
      if (ioff > 0)             idum(1:ioff)     = s% chan(1:ioff)
      if (ioff + n < s% n_chan) idum(ioff+1:nch) = s% chan(ioff+n+1:s% n_chan)
      deallocate(s% chan)
      s% chan => idum
    end if
    if (associated(s% chidx )) then
      allocate(idum(nch))
      if (ioff > 0)             idum(1:ioff)     = s% chidx(1:ioff)
      if (ioff + n < s% n_chan) idum(ioff+1:nch) = s% chidx(ioff+n+1:s% n_chan)
      deallocate(s% chidx)
      s% chidx => idum
    end if
    if (associated(s% iflag )) then
      allocate(idum(nch))
      if (ioff > 0)             idum(1:ioff)     = s% iflag(1:ioff)
      if (ioff + n < s% n_chan) idum(ioff+1:nch) = s% iflag(ioff+n+1:s% n_chan)
      deallocate(s% iflag)
      s% iflag => idum
    end if
    if (associated(s% band )) then
      allocate(idum(nch))
      if (ioff > 0)             idum(1:ioff)     = s% band(1:ioff)
      if (ioff + n < s% n_chan) idum(ioff+1:nch) = s% band(ioff+n+1:s% n_chan)
      deallocate(s% band)
      s% band => idum
    end if
    if (associated(s% flag )) then
      allocate(idum(nch))
      if (ioff > 0)             idum(1:ioff)     = s% flag(1:ioff)
      if (ioff + n < s% n_chan) idum(ioff+1:nch) = s% flag(ioff+n+1:s% n_chan)
      deallocate(s% flag)
      s% flag => idum
    end if
    if (associated(s% wavenum  )) then
      allocate(wdum(nch))
      if (ioff > 0)             wdum(1:ioff)     = s% wavenum (1:ioff)
      if (ioff + n < s% n_chan) wdum(ioff+1:nch) = s% wavenum (ioff+n+1:s% n_chan)
      deallocate(s% wavenum )
      s% wavenum  => wdum
    end if
    if (associated(s% emis_land)) then
      allocate(wdum(nch))
      if (ioff > 0)             wdum(1:ioff)     = s% emis_land(1:ioff)
      if (ioff + n < s% n_chan) wdum(ioff+1:nch) = s% emis_land(ioff+n+1:s% n_chan)
      deallocate(s% emis_land )
      s% emis_land=> wdum
    end if
    if (associated(s% emis_snow)) then
      allocate(wdum(nch))
      if (ioff > 0)             wdum(1:ioff)     = s% emis_snow(1:ioff)
      if (ioff + n < s% n_chan) wdum(ioff+1:nch) = s% emis_snow(ioff+n+1:s% n_chan)
      deallocate(s% emis_snow )
      s% emis_snow=> wdum
    end if
    if (associated(s% emis_sice)) then
      allocate(wdum(nch))
      if (ioff > 0)             wdum(1:ioff)     = s% emis_sice(1:ioff)
      if (ioff + n < s% n_chan) wdum(ioff+1:nch) = s% emis_sice(ioff+n+1:s% n_chan)
      deallocate(s% emis_sice )
      s% emis_sice=> wdum
    end if
    if (associated(s% emis_pc  )) then
      allocate(wdum2(size(s%emis_pc,1),nch))
      if (ioff > 0)             wdum2(:,1:ioff)     = s% emis_pc  (:,1:ioff)
      if (ioff + n < s% n_chan) wdum2(:,ioff+1:nch) = s% emis_pc  (:,ioff+n+1:s% n_chan)
      deallocate(s% emis_pc )
      s% emis_pc  => wdum2
    end if
    if (associated(s% i1)) then
      allocate(idum(nch))
      if (ioff > 0)             idum(1:ioff)     = s% i1(1:ioff)
      if (ioff + n < s% n_chan) idum(ioff+1:nch) = s% i1(ioff+n+1:s% n_chan)
      deallocate(s% i1 )
      s% i1=> idum
    end if
    if (associated(s% w1)) then
      allocate(wdum(nch))
      if (ioff > 0)             wdum(1:ioff)     = s% w1(1:ioff)
      if (ioff + n < s% n_chan) wdum(ioff+1:nch) = s% w1(ioff+n+1:s% n_chan)
      deallocate(s% w1 )
      s% w1=> wdum
    end if
    if (associated(s% emis_ind)) then
      allocate(idum3(size(s%emis_ind,1),0:size(s%emis_ind,2)-1,nch))
      if (ioff > 0)             idum3(:,:,1:ioff)     = s% emis_ind  (:,:,1:ioff)
      if (ioff + n < s% n_chan) idum3(:,:,ioff+1:nch) = s% emis_ind  (:,:,ioff+n+1:s% n_chan)
      deallocate(s% emis_ind )
      s% emis_ind => idum3
    end if
    if (associated(s% tskin_ind)) then
      allocate(idum3(size(s%tskin_ind,1),0:size(s%tskin_ind,2)-1,nch))
      if (ioff > 0)             idum3(:,:,1:ioff)     = s% tskin_ind  (:,:,1:ioff)
      if (ioff + n < s% n_chan) idum3(:,:,ioff+1:nch) = s% tskin_ind  (:,:,ioff+n+1:s% n_chan)
      deallocate(s% tskin_ind )
      s% tskin_ind => idum3
    end if
    if (associated(s% var  )) then
      allocate(sdum(nch))
      if (ioff > 0)             sdum(1:ioff)     = s% var (1:ioff)
      if (ioff + n < s% n_chan) sdum(ioff+1:nch) = s% var (ioff+n+1:s% n_chan)
      deallocate(s% var )
      s% var  => sdum
    end if
    if (associated(s% oe_str)) then
      allocate(odum(nch))
      if (ioff > 0)             odum(1:ioff)     = s% oe_str (1:ioff)
      if (ioff + n < s% n_chan) odum(ioff+1:nch) = s% oe_str (ioff+n+1:s% n_chan)
      deallocate(s% oe_str )
      s% oe_str  => odum
    end if
    if (associated(s% oe_trf)) then
      allocate(trfdum(nch))
      if (ioff > 0)             trfdum(1:ioff)     = s% oe_trf (1:ioff)
      if (ioff + n < s% n_chan) trfdum(ioff+1:nch) = s% oe_trf (ioff+n+1:s% n_chan)
      call destruct(s% oe_trf )
      s% oe_trf  => trfdum
    end if
    if (associated(s% bc%ib_ch )) then
      allocate(idum(nch))
      if (ioff > 0)             idum(1:ioff)     = s% bc%ib_ch(1:ioff)
      if (ioff + n < s% n_chan) idum(ioff+1:nch) = s% bc%ib_ch(ioff+n+1:s% n_chan)
      deallocate(s% bc%ib_ch)
      s% bc%ib_ch => idum
    end if

    s% n_chan = nch

  end subroutine del_instr_rad_set
!------------------------------------------------------------------------------

  !> Delete a channel from t_rad_set
  subroutine del_chan_rad_set(s, ic, instr)
    type(t_rad_set), intent(inout) :: s
    integer,         intent(in)    :: ic      !< channel to be deleted.
    integer,         intent(in)    :: instr   !< instrument index

    integer,             pointer :: idum     (:) => NULL()
    integer,             pointer :: idum3(:,:,:) => NULL()
    real(sp),            pointer :: sdum     (:) => NULL()
    character(len=m_oe), pointer :: odum     (:) => NULL()
    type(t_range_fparse),pointer :: trfdum   (:) => NULL()
    real(wp),            pointer :: wdum     (:) => NULL()
    real(wp),            pointer :: wdum2  (:,:) => NULL()
    integer                      :: ni, nch

    ni  = s%n_instr
    nch = s% n_chan - 1
    s% n_ch_i(instr) = s% n_ch_i(instr) - 1
    if (ni > instr) s% o_ch_i(instr+1:ni) = s% o_ch_i(instr+1:ni) - 1

    if (associated(s% chan )) then
      allocate(idum(nch))
      if (ic > 1)         idum(1:ic-1) = s% chan(1:ic-1)
      if (ic < s% n_chan) idum(ic:nch) = s% chan(ic+1:s% n_chan)
      deallocate(s% chan)
      s% chan => idum
    end if
    if (associated(s% chidx)) then
      allocate(idum(nch))
      if (ic > 1)         idum(1:ic-1) = s% chidx(1:ic-1)
      if (ic < s% n_chan) idum(ic:nch) = s% chidx(ic+1:s% n_chan)
      deallocate(s% chidx)
      s% chidx=> idum
    end if
    if (associated(s% iflag)) then
      allocate(idum(nch))
      if (ic > 1)         idum(1:ic-1) = s% iflag(1:ic-1)
      if (ic < s% n_chan) idum(ic:nch) = s% iflag(ic+1:s% n_chan)
      deallocate(s% iflag)
      s% iflag=> idum
    end if
    if (associated(s% flag )) then
      allocate(idum(nch))
      if (ic > 1)         idum(1:ic-1) = s% flag(1:ic-1)
      if (ic < s% n_chan) idum(ic:nch) = s% flag(ic+1:s% n_chan)
      deallocate(s% flag)
      s% flag => idum
    end if
    if (associated(s% band )) then
      allocate(idum(nch))
      if (ic > 1)         idum(1:ic-1) = s% band(1:ic-1)
      if (ic < s% n_chan) idum(ic:nch) = s% band(ic+1:s% n_chan)
      deallocate(s% band)
      s% band => idum
    end if
    if (associated(s% emis_land)) then
      allocate(wdum(nch))
      if (ic > 1)         wdum(1:ic-1) = s% emis_land(1:ic-1)
      if (ic < s% n_chan) wdum(ic:nch) = s% emis_land(ic+1:s% n_chan)
      deallocate(s% emis_land)
      s% emis_land => wdum
    end if
    if (associated(s% emis_snow)) then
      allocate(wdum(nch))
      if (ic > 1)         wdum(1:ic-1) = s% emis_snow(1:ic-1)
      if (ic < s% n_chan) wdum(ic:nch) = s% emis_snow(ic+1:s% n_chan)
      deallocate(s% emis_snow)
      s% emis_snow => wdum
    end if
    if (associated(s% emis_sice)) then
      allocate(wdum(nch))
      if (ic > 1)         wdum(1:ic-1) = s% emis_sice(1:ic-1)
      if (ic < s% n_chan) wdum(ic:nch) = s% emis_sice(ic+1:s% n_chan)
      deallocate(s% emis_sice)
      s% emis_sice => wdum
    end if
    if (associated(s% emis_pc  )) then
      allocate(wdum2(size(s%emis_pc,1),nch))
      if (ic > 1)         wdum2(:,1:ic-1) = s% emis_pc  (:,1:ic-1)
      if (ic < s% n_chan) wdum2(:,ic:nch) = s% emis_pc  (:,ic+1:s% n_chan)
      deallocate(s% emis_pc)
      s% emis_pc   => wdum2
    end if
    if (associated(s% i1)) then
      allocate(idum(nch))
      if (ic > 1)         idum(1:ic-1) = s% i1(1:ic-1)
      if (ic < s% n_chan) idum(ic:nch) = s% i1(ic+1:s% n_chan)
      deallocate(s% i1)
      s% i1 => idum
    end if
    if (associated(s% w1)) then
      allocate(wdum(nch))
      if (ic > 1)         wdum(1:ic-1) = s% w1(1:ic-1)
      if (ic < s% n_chan) wdum(ic:nch) = s% w1(ic+1:s% n_chan)
      deallocate(s% w1)
      s% w1 => wdum
    end if
    if (associated(s% emis_ind)) then
      allocate(idum3(size(s%emis_ind,1),0:size(s%emis_ind,2)-1,nch))
      if (ic > 1)         idum3(:,:,1:ic-1) = s% emis_ind(:,:,1:ic-1)
      if (ic < s% n_chan) idum3(:,:,ic:nch) = s% emis_ind(:,:,ic+1:s% n_chan)
      deallocate(s% emis_ind)
      s% emis_ind => idum3
    end if
    if (associated(s% tskin_ind)) then
      allocate(idum3(size(s%tskin_ind,1),0:size(s%tskin_ind,2)-1,nch))
      if (ic > 1)         idum3(:,:,1:ic-1) = s% tskin_ind(:,:,1:ic-1)
      if (ic < s% n_chan) idum3(:,:,ic:nch) = s% tskin_ind(:,:,ic+1:s% n_chan)
      deallocate(s% tskin_ind)
      s% tskin_ind => idum3
    end if
    if (associated(s% var  )) then
      allocate(sdum(nch))
      if (ic > 1)         sdum(1:ic-1) = s% var (1:ic-1)
      if (ic < s% n_chan) sdum(ic:nch) = s% var (ic+1:s% n_chan)
      deallocate(s% var )
      s% var  => sdum
    end if
    if (associated(s% oe_str  )) then
      allocate(odum(nch))
      if (ic > 1)         odum(1:ic-1) = s% oe_str (1:ic-1)
      if (ic < s% n_chan) odum(ic:nch) = s% oe_str (ic+1:s% n_chan)
      deallocate(s% oe_str )
      s% oe_str  => odum
    end if
    if (associated(s% oe_trf  )) then
      allocate(trfdum(nch))
      if (ic > 1)         trfdum(1:ic-1) = s% oe_trf (1:ic-1)
      if (ic < s% n_chan) trfdum(ic:nch) = s% oe_trf (ic+1:s% n_chan)
      deallocate(s% oe_trf )
      s% oe_trf  => trfdum
    end if
    if (associated(s% bc%ib_ch)) then
      allocate(idum(nch))
      if (ic > 1)         idum(1:ic-1) = s% bc%ib_ch(1:ic-1)
      if (ic < s% n_chan) idum(ic:nch) = s% bc%ib_ch(ic+1:s% n_chan)
      deallocate(s% bc%ib_ch)
      s% bc%ib_ch=> idum
    end if

    s% n_chan = nch

  end subroutine del_chan_rad_set
!------------------------------------------------------------------------------

  subroutine del_chan_radv(r, ichan, instr)
    type(t_radv), intent(inout) :: r
    integer,      intent(in)    :: ichan
    integer,      intent(in)    :: instr

    integer :: nch

    call del_chan_rad_set(r%i, ichan, instr)
    nch = r%i%n_chan
    if (ichan <= nch) then
      if (associated(r%not_rej  )) r%not_rej  (ichan:nch,:)   = r%not_rej  (ichan+1:nch+1,:)
      if (associated(r%cloudy   )) r%cloudy   (ichan:nch,:)   = r%cloudy   (ichan+1:nch+1,:)
      if (associated(r%state    )) r%state    (ichan:nch,:)   = r%state    (ichan+1:nch+1,:)
      if (associated(r%flags    )) r%flags    (ichan:nch,:)   = r%flags    (ichan+1:nch+1,:)
      if (associated(r%valid    )) r%valid    (ichan:nch,:)   = r%valid    (ichan+1:nch+1,:)
      if (associated(r%sinfl    )) r%sinfl    (ichan:nch,:)   = r%sinfl    (ichan+1:nch+1,:)
      if (associated(r%plevel   )) r%plevel   (ichan:nch,:)   = r%plevel   (ichan+1:nch+1,:)
      if (associated(r%bt_fg    )) r%bt_fg    (ichan:nch,:)   = r%bt_fg    (ichan+1:nch+1,:)
      if (associated(r%bt_obs   )) r%bt_obs   (ichan:nch,:)   = r%bt_obs   (ichan+1:nch+1,:)
      if (associated(r%bt_bcor  )) r%bt_bcor  (ichan:nch,:)   = r%bt_bcor  (ichan+1:nch+1,:)
      if (associated(r%bt_fg_cs )) r%bt_fg_cs (ichan:nch,:)   = r%bt_fg_cs (ichan+1:nch+1,:)
      if (associated(r%rad_fg   )) r%rad_fg   (ichan:nch,:)   = r%rad_fg   (ichan+1:nch+1,:)
      if (associated(r%rad_fg_cs)) r%rad_fg_cs(ichan:nch,:)   = r%rad_fg_cs(ichan+1:nch+1,:)
      if (associated(r%bcor_    )) r%bcor_    (ichan:nch,:)   = r%bcor_    (ichan+1:nch+1,:)
      if (associated(r%emiss    )) r%emiss    (ichan:nch,:)   = r%emiss    (ichan+1:nch+1,:)
      if (associated(r%H_t      )) r%H_t      (ichan:nch,:,:) = r%H_t      (ichan+1:nch+1,:,:)
      if (associated(r%H_t      )) r%H_q      (ichan:nch,:,:) = r%H_q      (ichan+1:nch+1,:,:)
      if (associated(r%H_ts     )) r%H_ts     (ichan:nch,:)   = r%H_ts     (ichan+1:nch+1,:)
      if (associated(r%H_ps     )) r%H_ps     (ichan:nch,:)   = r%H_ps     (ichan+1:nch+1,:)
      if (associated(r%H_em_pc  )) r%H_em_pc  (ichan:nch,:,:) = r%H_em_pc  (ichan+1:nch+1,:,:)
      if (associated(r%K_t      )) r%K_t      (ichan:nch,:,:) = r%K_t      (ichan+1:nch+1,:,:)
      if (associated(r%K_t      )) r%K_q      (ichan:nch,:,:) = r%K_q      (ichan+1:nch+1,:,:)
      if (associated(r%K_ts     )) r%K_ts     (ichan:nch,:)   = r%K_ts     (ichan+1:nch+1,:)
      if (associated(r%K_ps     )) r%K_ps     (ichan:nch,:)   = r%K_ps     (ichan+1:nch+1,:)
    end if

  end subroutine del_chan_radv
!------------------------------------------------------------------------------
  !> Computation of sat_zenith and sat_azimuth angles. Compare GetSatAngles.c in
  !! NWCSAF software
  elemental subroutine get_satangles(lat, lon, sat_lon, sat_hgt, sat_zenith, sat_azimuth, equator_radius, n_pol_radius)
    real(kind=wp),      intent(in)           :: lat            ! latitude  of observation
    real(kind=wp),      intent(in)           :: lon            ! longitude of observation
    real(kind=wp),      intent(in)           :: sat_lon        ! longitude of satellite
    real(kind=wp),      intent(in)           :: sat_hgt        ! height of satellite
    real(kind=wp),      intent(out)          :: sat_zenith
    real(kind=wp),      intent(out)          :: sat_azimuth
    real(kind=wp),      intent(in), optional :: equator_radius ! height of satellite
    real(kind=wp),      intent(in), optional :: n_pol_radius   ! height of satellite

    real(kind=dp), parameter :: equator_radius_def = 6378.169_wp
    real(kind=dp), parameter :: n_pol_radius_def   = 6356.5838_wp
    real(kind=dp), parameter :: s_pol_radius_def   = 6356.5838_wp
    real(kind=wp), parameter :: rfill              = NF90_FILL_FLOAT

    real(kind=dp)            :: sat(3)                 !< cartesian coordinates of satellite
    real(kind=dp)            :: cosf, sinf, cosl, sinl !< sin/cos of lat/lon
    real(kind=dp)            :: r_e                    !< earth radius at (lat,lon)
    real(kind=dp)            :: earth_normal(3)        !< normalized cart. coord. of (lat,lon)
    real(kind=dp)            :: earth(3)               !< cart. coord. of (lat,lon)
    real(kind=dp)            :: earth_sat(3)           !< vector earth - satellite
    real(kind=dp)            :: mod
    real(kind=dp)            :: cos_sat
    real(kind=dp)            :: cosll, sinll, cosa, sina
    real(kind=dp)            :: c, dazi
    real(kind=dp)            :: equ_rad, npo_rad

    if (lon == sat_lon .and. lat == 0._wp) then
       sat_zenith  = 0._wp
       sat_azimuth = 0._wp
       return
    else
       sat_zenith  = rfill
       sat_azimuth = rfill
    end if

    if (lat == rfill .or. lon == rfill) return

    if (present(equator_radius)) then
      equ_rad = equator_radius
    else
      equ_rad = equator_radius_def
    end if
    if (present(n_pol_radius)) then
      npo_rad = n_pol_radius
    else
      npo_rad = n_pol_radius_def
    end if

    sat(1) = (sat_hgt + equ_rad) * dcos(sat_lon * deg2rad)
    sat(2) = (sat_hgt + equ_rad) * dsin(sat_lon * deg2rad)
    sat(3) = 0._dp

    cosf = dcos(lat * deg2rad)
    sinf = dsin(lat * deg2rad)
    cosl = dcos(lon * deg2rad)
    sinl = dsin(lon * deg2rad)

    r_e = npo_rad / dsqrt(1._dp - ((equ_rad**2 - npo_rad**2) / &
         (equ_rad**2) * cosf**2))

    earth_normal(1) = cosf * cosl
    earth_normal(2) = cosf * sinl
    earth_normal(3) = sinf

    earth = r_e * earth_normal

    earth_sat = sat - earth

    mod = dsqrt(earth_sat(1)**2 + earth_sat(2)**2 + earth_sat(3)**2)
    if (mod == 0._dp) return

    earth_sat = earth_sat / mod

    cos_sat = earth_normal(1) * earth_sat(1) + &
              earth_normal(2) * earth_sat(2) + &
              earth_normal(3) * earth_sat(3)

    sat_zenith = real(dacos(cos_sat) * rad2deg, kind=wp)

    ! azimuth angle
    cosll = dcos((lon - sat_lon) * deg2rad)
    sinll = dsin((lon - sat_lon) * deg2rad)
    cosa  = cosf * cosll
    sina  = dsin(dacos(cosa))

    if (sina == 0._dp) return

    c = datan2(sinll/sina, -sinf*cosll/sina)
    if (c < 0._dp) c = c + 2._dp * pi
    dazi = 2._dp * pi - c
    if (dazi > 2._dp * pi) dazi = dazi - 2._dp * pi

    sat_azimuth = real(dazi * rad2deg, kind=wp)

  end subroutine get_satangles
!------------------------------------------------------------------------------

  real(wp) function sat_sun_azimuth_angle(stzen, stazi, sunzen, sunazi)
  !------------------------------------------------------
  ! derive angle between satellite zenith/azimuth and sun
  !------------------------------------------------------
  real(wp), intent(in) :: stzen    ! satellite zenith angle
  real(wp), intent(in) :: stazi    ! satellite azimuth angle
  real(wp), intent(in) :: sunzen   ! sun zenith angle
  real(wp), intent(in) :: sunazi   ! sun azimuth angle

    real(wp) :: pi
    real(wp) :: qpi180
    real(wp) :: q180pi
    real(wp) :: x1, x2, y1, y2, z1, z2

    if (  (stzen  /= NF_FILL_FLOAT).and.&
         &(sunzen /= NF_FILL_FLOAT).and.&
         &(stazi  /= NF_FILL_FLOAT).and.&
         &(sunazi /= NF_FILL_FLOAT)) then
      pi     = acos(-1._wp)
      qpi180 = pi/180._wp
      q180pi = 180._wp/pi
      z1     = cos (stzen  * qpi180)
      z2     = cos (sunzen * qpi180)
      y1     = cos (stazi  * qpi180)
      y2     = cos (sunazi * qpi180)
      x1     = sqrt (1._wp - y1**2)
      x2     = sqrt (1._wp - y2**2)
      sat_sun_azimuth_angle = q180pi * &
           acos (z1*z2+sqrt((1._wp-z1**2)*(1._wp-z2**2))*(x1*x2+y1*y2))
    else
      sat_sun_azimuth_angle = NF_FILL_FLOAT
    endif

  end function sat_sun_azimuth_angle
!------------------------------------------------------------------------------
  subroutine warning(text, nowarn)
    character(len=*), intent(in)           :: text
    logical,          intent(in), optional :: nowarn

#ifdef __SX__
    integer :: mxlen = 133
#else
    integer :: mxlen = 9999
#endif
    character(len=11), parameter :: pref = '*** WARNING'
    integer                      :: l, lp, is, ie, stat
    logical                      :: warn

    warn = .true.
    IF (present(nowarn)) warn = .not.nowarn

    l = len_trim(text)
    IF (warn) THEN
      lp = len_trim(pref) + 1
    ELSE
      lp = 0
    END IF
    is = 1
    DO
      ie = min(l, is-1+mxlen-lp)
      IF (ie < is) EXIT
      IF (warn) THEN
        write(0,'(A)',iostat=stat) trim(pref)//' '//text(is:ie)
      ELSE
        write(0,'(A)',iostat=stat) text(is:ie)
      END IF
      is = ie + 1
      IF (is > l) EXIT
    END DO

  end subroutine warning
!------------------------------------------------------------------------------
#ifdef HRIT
  subroutine read_hrit(fnames, r, lprint, lread, istart, iend, qual_mask,       &
       &qual_geo_mask, qual_rad_mask, minlat, maxlat, minlon, maxlon, time_mode,&
       &zenmax, thin_col, thin_line, status, unit, pe)
    character(len=*), intent(in)            :: fnames(:)
    type(t_radv),     intent(out)           :: r
    logical,          intent(in),  optional :: lprint
    logical,          intent(in),  optional :: lread
    integer,          intent(in),  optional :: istart
    integer,          intent(in),  optional :: iend
    integer,          intent(in),  optional :: qual_mask
    integer,          intent(in),  optional :: qual_geo_mask
    integer,          intent(in),  optional :: qual_rad_mask
    real,             intent(in),  optional :: minlat
    real,             intent(in),  optional :: maxlat
    real,             intent(in),  optional :: minlon
    real,             intent(in),  optional :: maxlon
    real,             intent(in),  optional :: zenmax
    integer,          intent(in),  optional :: thin_col
    integer,          intent(in),  optional :: thin_line
    character(len=*), intent(in),  optional :: time_mode
    integer,          intent(out), optional :: status
    integer,          intent(in),  optional :: unit ! stdout unit
    integer,          intent(in),  optional :: pe

    character(len=9), parameter :: proc = 'read_hrit'
    character(len=180)          :: msg  = ''
    integer, allocatable :: satid(:)
    integer, allocatable :: qual(:)
    integer, allocatable :: qual_rad(:)
    integer, allocatable :: qual_geo(:)
    integer              :: stat
    integer              :: verbosity, nfovs, nchan, nsat, i
    integer              :: is, ie
    integer              :: ibit
    logical              :: lpr, lrd, lst

    if (present(lprint)) then
       lpr = lprint
    else
       lpr = .true.
    end if
    if (lpr) then
       verbosity = hrit_talk
    else
       verbosity = hrit_silent
    end if
    lst = present(status)

#define ERR(STAT) if (lst) then; call warning(trim(msg)); status=STAT; return; else; call finish(proc,trim(msg)); endif

    call hrit_open(fnames, stat, nfovs=nfovs, nchan=nchan, nsat=nsat,           &
         verbosity=verbosity, minlat=minlat, maxlat=maxlat, minlon=minlon,      &
         maxlon=maxlon, zenmax=zenmax, thin_line=thin_line, thin_col=thin_col,  &
         print_unit=unit, lskip=.true., pe=pe)
    if (stat /= 0 .and. stat /= 25) then
      write(msg,*) 'hrit_open failed, stat=',stat
      ERR(stat)
    end if
    if (nsat > 1) then
      msg = 'not more than one satellite supported.'
      ERR(11)
    end if

    r%n_rec        = nfovs
    r%i%n_fov      = nfovs
    r%i%n_chan     = nchan
    r%i%n_instr    = 1
    r%i%n_sens     = 1
!    r%i%n_pc_emiss = -1

    if (present(lread)) then
       lrd = lread
    else
       lrd = .true.
    end if
    if (lrd) then
       is = 1
       if (present(istart)) is  = max(istart,1)
       ie = nfovs
       if (present(iend)) then
          if (iend > 0) then
             ie = min(iend,ie)
          end if
       end if
       nfovs  = ie - is + 1

       r%filename       = 'HRIT'
       r%n_rec          = nfovs

       if (nfovs <= 0) return

       allocate(satid(nfovs), r%time(nfovs), r%date(nfovs), r%scanl(nfovs),     &
            r%fov(nfovs), qual(nfovs), qual_rad(nfovs), qual_geo(nfovs),        &
            r%dlat(nfovs), r%dlon(nfovs), r%center(nfovs),                      &
            r%stazi(nfovs), r%stzen(nfovs), r%i%chan(nchan), r%obsnum(nfovs),   &
            r%valid(nchan,nfovs), r%bt_obs(nchan,nfovs), stat=stat)
       if (stat /= 0) then
         msg = 'allocation failed.'
         ERR(12)
       end if
       r%valid(:,:) = .true.

       call hrit_get_data(stat, satid=satid, time=r%time, date=r%date,          &
            line=r%scanl, column=r%fov, quality=qual, quality_rad=qual_rad,     &
            quality_geo=qual_geo, lat=r%dlat, lon=r%dlon, bt=r%bt_obs,          &
            sat_zenith=r%stzen, sat_azimuth=r%stazi, channel=r%i%chan,          &
            istart=is, iend=ie, time_mode=time_mode, pe=pe)
       if (stat /= 0) then
         msg = 'hrit_get_data failed.'
         ERR(13)
       end if

       r%i%satid        = satid(1)
       r%i%grid         =  21 ! SEVIRI
       r%i%grid_wmo     = 207 ! SEVIRI
       r%i%instr(1)     = r%i%grid
       r%i%instr_wmo(1) = r%i%grid_wmo
       r%i%n_ch_i(1)    = nchan
       r%i%o_ch_i(1)    = 0

       r%obsnum(1:nfovs) = (/ (i, i=istart, iend) /)
       r%center(1:nfovs) = 254

       if (present(qual_mask)) then
          do ibit = 0, 4
             if (.not. btest(qual_mask, ibit)) then
                forall(i=1:nfovs, qual(i) == ibit)
                   r%valid(i,:) = .false.
                end forall
             end if
          end do
       end if
       if (present(qual_geo_mask)) then
          do ibit = 0, 4
             if (.not. btest(qual_geo_mask, ibit)) then
                forall(i=1:nfovs, qual_geo(i) == ibit)
                   r%valid(i,:) = .false.
                end forall
             end if
          end do
       end if
       if (present(qual_rad_mask)) then
          do ibit = 0, 4
             if (.not. btest(qual_rad_mask, ibit)) then
                forall(i=1:nfovs, qual_rad(i) == ibit)
                   r%valid(i,:) = .false.
                end forall
             end if
          end do
       end if

       call hrit_close(stat, pe=pe)
       if (stat /= 0) then
         msg = 'hrit_close failed.'
         ERR(14)
       end if

    end if

#undef(ERR)

  end subroutine read_hrit

#endif

!------------------------------------------------------------------------------


  !> TODO: write spec array
  subroutine write_rtovp(r_arr, m_prof, mode, status, name, n_prof, n_par,     &
        n_lev, n_chan, n_pc, offset, par, cloud_mask, date, pe, lsort, lH, lK, &
        lB, lHBH, attrs, keep_sorting)
    type(t_radv),     intent(in),  target   :: r_arr(:)   ! data stucture
    integer,          intent(in)            :: m_prof ! Mode: see MP_* parameters
    integer,          intent(in),  optional :: mode   ! Mode: see M_* parameters
    integer,          intent(out), optional :: status ! exit status
    character(len=*), intent(in),  optional :: name   ! Name of file
    integer,          intent(in),  optional :: n_prof ! Total number of profiles in final file
    integer,          intent(in),  optional :: n_par  ! Total number of params in final file
    integer,          intent(in),  optional :: n_lev  ! Total number of levels in final file
    integer,          intent(in),  optional :: n_chan ! Total number of channels in final file
    integer,          intent(in),  optional :: n_pc   ! Total number of emissivity PCs in final file
    integer,          intent(in),  optional :: offset ! Profile offset
    character(len=*), intent(in),  optional :: par    ! Parameter-variable value
    integer,          intent(in),  optional :: cloud_mask(:) ! Mask for cloud type
    character(len=*), intent(in),  optional :: date
    integer,          intent(in),  optional :: pe
    logical,          intent(in),  optional :: lsort
    logical,          intent(in),  optional :: lH
    logical,          intent(in),  optional :: lK
    logical,          intent(in),  optional :: lB
    logical,          intent(in),  optional :: lHBH
    type(t_nc_attr),  intent(in),  optional :: attrs(:)
    target attrs
    logical,          intent(in),  optional :: keep_sorting

    integer, parameter :: mx_par = 8

    type t_nc_var
      integer           :: id
      integer           :: type
      character(len=20) :: name
      character(len=10) :: unit
      character(len=50) :: long_name
      logical           :: lpar
    end type t_nc_var

    type t_nc_info
      character(len=180)      :: fname     = ''
      integer                 :: set_id    = -1
      integer                 :: id        = -1
      integer                 :: d_prof    = -1
      integer                 :: d_lev     = -1
      integer                 :: d_len8    = -1
      integer                 :: d_len16   = -1
      integer                 :: d_par     = -1
      integer                 :: d_chan    = -1
      integer                 :: d_cld     = -1
      integer                 :: d_pc      = -1
      type(t_nc_var)          :: v_par     = t_nc_var(-1, NF_CHAR,  'parameter', '',       'parameter (fg,ana,..)',    .true. )
      type(t_nc_var)          :: v_id      = t_nc_var(-1, NF_INT,   'id',        '',       'observation id',           .false.)
      type(t_nc_var)          :: v_lat     = t_nc_var(-1, NF_FLOAT, 'lat',       'degree', 'latitude of observation',  .false.)
      type(t_nc_var)          :: v_lon     = t_nc_var(-1, NF_FLOAT, 'lon',       'degree', 'longitude of observation', .false.)
      type(t_nc_var)          :: v_time    = t_nc_var(-1, NF_INT,   'time',      'min',    'observation minus reference time',&
                                                                                                                       .false.)
      type(t_nc_var)          :: v_statid  = t_nc_var(-1, NF_CHAR,  'statid',    '',       'station id in character form',    &
                                                                                                                       .false.)
      type(t_nc_var)          :: v_level   = t_nc_var(-1, NF_FLOAT, 'level',     'Pa',     'RTTOV pressure levels',    .false.)
      type(t_nc_var)          :: v_hh      = t_nc_var(-1, NF_FLOAT, 'h',         'm',      'COSMO layer boundaries',   .false.)
      type(t_nc_var)          :: v_ts      = t_nc_var(-1, NF_FLOAT, 'ts',        'K',      'surface temperature',      .true. )
      type(t_nc_var)          :: v_ps      = t_nc_var(-1, NF_FLOAT, 'ps',        'Pa',     'surface pressure',         .true. )
      type(t_nc_var)          :: v_v10m    = t_nc_var(-1, NF_FLOAT, 'v10m',      'm/s',    '10 m wind speed',          .true. )
      type(t_nc_var)          :: v_cldtop  = t_nc_var(-1, NF_FLOAT, 'cldtop',    'Pa',     'cloud top pressure',       .true. )
      type(t_nc_var)          :: v_cldfrac = t_nc_var(-1, NF_FLOAT, 'cldfrac',   '',       'cloud fraction profile',   .false. )
      type(t_nc_var)          :: v_cld     = t_nc_var(-1, NF_FLOAT, 'cld',       'g/m**3', 'cloud water content',      .true. )
      type(t_nc_var)          :: v_o3      = t_nc_var(-1, NF_FLOAT, 'o3',        'kg/kg',  'O3 conc. (dry air)',       .false.)
      type(t_nc_var)          :: v_co2     = t_nc_var(-1, NF_FLOAT, 'co2',       'kg/kg',  'CO2 conc. (dry air)',      .false.)
      type(t_nc_var)          :: v_snwfrac = t_nc_var(-1, NF_FLOAT, 'snowfrac',  '',       'snow fraction',            .true. )
      type(t_nc_var)          :: v_t       = t_nc_var(-1, NF_FLOAT, 't',         'K',      'temperature',              .true. )
      type(t_nc_var)          :: v_q       = t_nc_var(-1, NF_FLOAT, 'q',         '',       'specific humidity',        .true. )
      type(t_nc_var)          :: v_emiss   = t_nc_var(-1, NF_FLOAT, 'emissivity','',       'emissivity',               .true. )
      type(t_nc_var)          :: v_emis_pc = t_nc_var(-1, NF_FLOAT, 'emis_pc',   '',       'emissivity principal component &
                                                                                                       &coefficients', .true. )
      type(t_nc_var)          :: v_chan    = t_nc_var(-1, NF_INT,   'chan',      '',       'channel number',           .false.)
      type(t_nc_var)          :: v_instr   = t_nc_var(-1, NF_INT,   'instr',     '',       'WMO instrument number',    .false.)
      type(t_nc_var)          :: v_H_t     = t_nc_var(-1, NF_FLOAT, 'H_t',       'K',      'H temperature',            .false.)
      type(t_nc_var)          :: v_H_ts    = t_nc_var(-1, NF_FLOAT, 'H_ts',      'K',      'H surface temperature',    .false.)
      type(t_nc_var)          :: v_H_q     = t_nc_var(-1, NF_FLOAT, 'H_q',       '',       'H specific humidity',      .false.)
      type(t_nc_var)          :: v_H_ps    = t_nc_var(-1, NF_FLOAT, 'H_ps',      'Pa',     'H surface pressure',       .false.)
      type(t_nc_var)          :: v_H_em_pc = t_nc_var(-1, NF_FLOAT, 'H_emis_pc', '',       'H emissivity',             .false.)
      type(t_nc_var)          :: v_K_t     = t_nc_var(-1, NF_FLOAT, 'K_t',       'K',      'K temperature',            .false.)
      type(t_nc_var)          :: v_K_ts    = t_nc_var(-1, NF_FLOAT, 'K_ts',      'K',      'K surface temperature',    .false.)
      type(t_nc_var)          :: v_K_q     = t_nc_var(-1, NF_FLOAT, 'K_q',       '',       'K specific humidity',      .false.)
      type(t_nc_var)          :: v_K_ps    = t_nc_var(-1, NF_FLOAT, 'K_ps',      'Pa',     'K surface pressure',       .false.)
      type(t_nc_var)          :: v_B_tt    = t_nc_var(-1, NF_FLOAT, 'B_tt',      'K**2',   'B temperature',            .false.)
      type(t_nc_var)          :: v_B_qq    = t_nc_var(-1, NF_FLOAT, 'B_qq',      '',       'B humidity',               .false.)
      type(t_nc_var)          :: v_B_tq    = t_nc_var(-1, NF_FLOAT, 'B_tq',      'K',      'B temperature-humidity',   .false.)
      type(t_nc_var)          :: v_R       = t_nc_var(-1, NF_FLOAT, 'R  ',       '',       'obs eror covariances  ',   .false.)
      type(t_nc_var)          :: v_HBHR    = t_nc_var(-1, NF_FLOAT, 'HBHR',      '',       'Bg+obs-error-covs obs-space',.false.)
      type(t_nc_var)          :: v_tsm     = t_nc_var(-1, NF_FLOAT, 'tsm',       'K',      'model surface temperature',.false. )
      ! --- params ---
      integer                 :: npar           =  0
      character(len=16)       :: params(mx_par) =  ''
      ! --- other stuff ---
      integer,        pointer :: cld_mask(:)    => null()
      logical                 :: lsorted        =  .false.
      integer,        pointer :: ipos(:)        => null()
    end type t_nc_info

    ! Variables fort sorted output
    logical, pointer :: lsorted
    type t_sort
      integer :: iprof  ! index in loop over sets and profiles
      integer :: reprt  ! Report ID
      integer :: ibox   ! box
    end type t_sort
    type(t_sort), allocatable :: s(:)
    integer,    pointer       :: ipos(:) => null()

    type t_i1
      integer,  pointer     :: d(:)  => null()
      logical               :: alloc =  .false.
      integer               :: fill  =  NF_FILL_INT
    end type t_i1
    type(t_i1)              :: i1(size(r_arr))
    integer,    allocatable :: i1dat(:)
    type t_i2
      integer,  pointer     :: d(:,:) => null()
      logical               :: alloc  =  .false.
      integer               :: fill   =  NF_FILL_INT
    end type t_i2
    type(t_i2)              :: i2(size(r_arr))
    integer,    allocatable :: i2dat(:,:)
    type t_r1
      real(wp), pointer     :: d(:)  => null()
      logical               :: alloc =  .false.
      real(wp)              :: fill  =  NF_FILL_FLOAT
    end type t_r1
    type(t_r1)              :: r1(size(r_arr))
    real(wp),   allocatable :: r1dat(:)
    type t_r2
      real(wp), pointer     :: d(:,:) => null()
      logical               :: alloc  =  .false.
      real(wp)              :: fill   =  NF_FILL_FLOAT
    end type t_r2
    type(t_r2)              :: r2(size(r_arr))
    real(wp),   allocatable :: r2dat(:,:)
    type t_r3
      real(wp), pointer     :: d(:,:,:) => null()
      logical               :: alloc    =  .false.
      real(wp)              :: fill     =  NF_FILL_FLOAT
    end type t_r3
    type(t_r3)              :: r3(size(r_arr))
    real(wp),   allocatable :: r3dat(:,:,:)

    type(t_nc_info), pointer, save :: nc_arr(:) => null()
    type(t_nc_info), pointer, save :: nc_tmp(:) => null()
    type(t_nc_info), pointer, save :: nc => null()
    integer,                  save :: n_nc = 0

    type(t_radv),      pointer     :: r   => null()
    character(len=80), parameter   :: proc='write_rtovp'
    type(t_nc_attr),   pointer     :: a
    character(len=180)             :: fname
    character(len=180)             :: msg
    character(len=16)              :: param = ''
    integer                        :: np(size(r_arr))
    integer                        :: idum(m_chan)
    integer                        :: stat
    integer                        :: m, i, j, istart, ipar, iset
    integer                        :: nset
    integer                        :: nprof, nlev, nchan, npar, ncld, nch, npc
    integer                        :: lev_mode
    integer                        :: pe_loc
    logical                        :: l_open, l_write, l_close
    logical                        :: l_cld, l_H, l_K, l_B, l_HBH, l_pc, l_keep_sorting
    logical                        :: l_o3, l_co2

    status = 0

    if (iand(m_prof,MP_NETCDF)==0) return

    if (present(pe)) then
      pe_loc = pe
    else
      pe_loc = -1
    end if

    if (present(mode)) then
      m = mode
    else
      m = WRM_OPEN + WRM_WRITE + WRM_CLOSE
    end if
    l_open  = iand(m, WRM_OPEN)  /= 0
    l_write = iand(m, WRM_WRITE) /= 0
    l_close = iand(m, WRM_CLOSE) /= 0

    if (present(lH)) then ; l_H=lH ; else ; l_H=.false. ; endif
    if (present(lK)) then ; l_K=lK ; else ; l_K=.false. ; endif
    if (present(lB)) then ; l_B=lB ; else ; l_B=.false. ; endif
    if (present(lHBH)) then ; l_HBH=lHBH ; else ; l_HBH=.false. ; endif
    if (present(keep_sorting)) then ; l_keep_sorting=keep_sorting ; else ; l_keep_sorting=.false. ; endif

    if (present(name)) then
      fname = trim(name)
    else
      fname = 'monRTOVP.nc'
    end if

    nset = size(r_arr)

    nc => null()
    do i = 1, n_nc
!      if (nc_arr(i)%fname == fname .or. nc_arr(i)%set_id == r_arr(:)%i%id) then
      if (nc_arr(i)%fname == fname) then
        nc => nc_arr(i)
        exit
      end if
    end do
    if (.not.associated(nc)) then
      allocate(nc_tmp(n_nc+1))
      if (associated(nc_arr)) then
        nc_tmp(1:n_nc) = nc_arr(1:n_nc)
        deallocate(nc_arr)
      end if
      nc_arr => nc_tmp
      n_nc = n_nc + 1
      nc => nc_arr(n_nc)
      nc%fname = fname
      if (nset > 0) nc%set_id = r_arr(1)%i%id
    end if

    np = r_arr(:)%n_rec
    if (present(n_prof)) then
      nprof = n_prof
    else
      nprof = sum(np(:))
    end if
!    if (nprof <= 0) return
    if (present(n_lev)) then
      nlev = n_lev
    else
      nlev  = maxval(r_arr(:)%i%n_lev)
    end if
    if (present(n_chan)) then
      nchan = n_chan
    else
      nchan = maxval(r_arr(:)%i%n_chan)
    end if
    if (present(n_pc)) then
      npc = n_pc
    else
      npc = maxval(r_arr(:)%n_pc_emiss)
    end if
    l_pc = (npc > 0)

    l_cld = present(cloud_mask)
    iset = -1
    do i = 1, nset
      if (associated(r_arr(i)%cld_fg)) then
        iset = i
        l_cld = .true.
        exit
      end if
    end do
    if (l_cld) then
      if (.not.associated(nc%cld_mask)) then
        if (present(cloud_mask)) then
          ncld = size(cloud_mask)
          allocate(nc%cld_mask(ncld))
          nc%cld_mask(:) = cloud_mask(:)
        else
          ncld = size(r_arr(iset)%cld_fg,1)
          allocate(nc%cld_mask(ncld))
          nc%cld_mask(:) = (/ (i, i=1, ncld)/)
        end if
      else
        ncld = size(nc%cld_mask)
      end if
    end if

    l_o3  = .false.
    l_co2 = .false.
    do i = 1, nset
      r => r_arr(i)
      l_o3  = l_o3  .or. any(r%i%iopts(1:r%i%n_instr)%use_o3  > 0)
      l_co2 = l_co2 .or. any(r%i%iopts(1:r%i%n_instr)%use_co2 > 0)
    end do

    lev_mode = r_arr(1)%i%gopts%lev_mode

    lsorted => nc%lsorted
    ipos    => nc%ipos
    if (present(lsort)) lsorted=lsort

#undef  NF_ERR
#define NF_ERR(nproc, text) if (stat/=0) then ; call nf_err(nproc,text) ; return ; endif

    if (l_open) then
      if (present(n_par)) then
        npar = n_par
      else
        npar = 1
      end if
      stat = nf90_create(fname,                       &
                         ior(NF_CLOBBER,              &
                             NF90_NETCDF4     ), nc%id)
!                            NF90_64BIT_OFFSET), nc%id)
      NF_ERR('nf90_create',trim(fname))
      if (present(date)) then
        stat = nf90_put_att(nc%id, NF_GLOBAL, 'analysis_time', trim(date))
        NF_ERR('nf90_put_att',date)
      end if
      ! Dimensions
      stat = nf90_def_dim(nc%id, 'profile'  ,nprof ,nc%d_prof )
      NF_ERR('nf90_def_dim','profile')
      stat = nf90_def_dim(nc%id, 'level'    ,nlev  ,nc%d_lev  )
      NF_ERR('nf90_def_dim','level')
      stat = nf90_def_dim(nc%id, 'len8'     ,8     ,nc%d_len8 )
      NF_ERR('nf90_def_dim','len8')
      stat = nf90_def_dim(nc%id, 'len16'    ,16    ,nc%d_len16)
      NF_ERR('nf90_def_dim','len16')
      stat = nf90_def_dim(nc%id, 'parameter',npar  ,nc%d_par  )
      NF_ERR('nf90_def_dim','parameter')
      stat = nf90_def_dim(nc%id, 'channel'  ,nchan ,nc%d_chan )
      NF_ERR('nf90_def_dim','channel')
!!$      if (l_cld) then
!!$        stat = nf90_def_dim(nc%id, 'cloud_type', ncld,nc%d_cld )
!!$        NF_ERR('nf90_def_dim','cloud_type')
!!$      end if
      if (l_pc) then
        stat = nf90_def_dim(nc%id, 'pc_coeff', npc,nc%d_pc)
        NF_ERR('nf90_def_dim','pc_coeff')
      end if
      ! Variables
#define DEF_VAR(var, dims) call def_var(var, dims) ; if (stat /= 0) then ; return ; endif
      DEF_VAR(nc%v_par,    (/nc%d_len16, nc%d_par/))
      DEF_VAR(nc%v_id,     (/nc%d_prof/))
      DEF_VAR(nc%v_lat,    (/nc%d_prof/))
      DEF_VAR(nc%v_lon,    (/nc%d_prof/))
      DEF_VAR(nc%v_time,   (/nc%d_prof/))
      DEF_VAR(nc%v_statid, (/nc%d_len8, nc%d_prof/))
      if (lev_mode <= 0) then
        DEF_VAR(nc%v_level,  (/nc%d_lev/))
      else
        DEF_VAR(nc%v_level,  (/nc%d_lev, nc%d_prof/))
      end if
      DEF_VAR(nc%v_ts,     (/nc%d_par, nc%d_prof/))
      DEF_VAR(nc%v_tsm,    (/nc%d_prof/))
      DEF_VAR(nc%v_ps,     (/nc%d_par, nc%d_prof/))
      DEF_VAR(nc%v_v10m,   (/nc%d_par, nc%d_prof/))
      DEF_VAR(nc%v_cldtop, (/nc%d_par, nc%d_prof/))
 !     DEF_VAR(nc%v_cfract, (/nc%d_par, nc%d_prof/))
      if (l_cld) then
     !   DEF_VAR(nc%v_cld_ind, (/nc%d_cld/))
     !   DEF_VAR(nc%v_cld,     (/nc%d_cld, nc%d_lev, nc%d_par, nc%d_prof/))
     !   DEF_VAR(nc%v_cldfrac, (/          nc%d_lev, nc%d_par, nc%d_prof/))
      else
     !   DEF_VAR(nc%v_cldfrac, (/nc%d_par, nc%d_prof/))
      end if
      if (l_o3 ) then
        DEF_VAR(nc%v_o3,   (/nc%d_lev, nc%d_prof/))
      end if
      if (l_co2) then
        DEF_VAR(nc%v_co2,  (/nc%d_lev, nc%d_prof/))
      end if
      DEF_VAR(nc%v_snwfrac,(/nc%d_par, nc%d_prof/))
      DEF_VAR(nc%v_t,      (/nc%d_lev, nc%d_par, nc%d_prof/))
      DEF_VAR(nc%v_q,      (/nc%d_lev, nc%d_par, nc%d_prof/))
      if (l_pc) then
        DEF_VAR(nc%v_emiss,   (/nc%d_chan,nc%d_par, nc%d_prof/))
        DEF_VAR(nc%v_emis_pc, (/nc%d_pc  ,nc%d_par, nc%d_prof/))
      end if
      if (l_H .or. l_K .or. l_pc .or. l_HBH) then
        DEF_VAR(nc%v_chan,    (/nc%d_chan,           nc%d_prof/))
        DEF_VAR(nc%v_instr,   (/nc%d_chan,           nc%d_prof/))
      end if
      if (l_H) then
        DEF_VAR(nc%v_H_ts,    (/nc%d_chan,           nc%d_prof/))
        DEF_VAR(nc%v_H_ps,    (/nc%d_chan,           nc%d_prof/))
        if (l_pc) then
          DEF_VAR(nc%v_H_em_pc,(/nc%d_chan, nc%d_pc, nc%d_prof/))
        endif
      end if
      if (l_K) then
        DEF_VAR(nc%v_K_ts,    (/nc%d_chan,           nc%d_prof/))
        DEF_VAR(nc%v_K_ps,    (/nc%d_chan,           nc%d_prof/))
      end if
      if (l_H) then
        DEF_VAR(nc%v_H_t,     (/nc%d_chan, nc%d_lev, nc%d_prof/))
        DEF_VAR(nc%v_H_q,     (/nc%d_chan, nc%d_lev, nc%d_prof/))
      end if
      if (l_K) then
        DEF_VAR(nc%v_K_t,     (/nc%d_chan, nc%d_lev, nc%d_prof/))
        DEF_VAR(nc%v_K_q,     (/nc%d_chan, nc%d_lev, nc%d_prof/))
      end if
      if (l_B) then
        DEF_VAR(nc%v_B_tt,    (/nc%d_lev,  nc%d_lev, nc%d_prof/))
        DEF_VAR(nc%v_B_qq,    (/nc%d_lev,  nc%d_lev, nc%d_prof/))
        DEF_VAR(nc%v_B_tq,    (/nc%d_lev,  nc%d_lev, nc%d_prof/))
      end if
      if (l_HBH) then
        DEF_VAR(nc%v_R  ,     (/nc%d_chan, nc%d_chan,nc%d_prof/))
        DEF_VAR(nc%v_HBHR,    (/nc%d_chan, nc%d_chan,nc%d_prof/))
      end if

      if (present(attrs)) then
        do i = 1, size(attrs)
          a => attrs(i)
          if (a%name /= '') then
            if (a%c /= '') then
              stat = nf90_put_att(nc%id, NF_GLOBAL, trim(a%name), trim(a%c))
              NF_ERR('nf90_put_att','')
            elseif (a%i /= NF_FILL_INT) then
              stat = nf90_put_att(nc%id, NF_GLOBAL, trim(a%name), a%i)
              NF_ERR('nf90_put_att','')
            end if
          end if
        end do
      end if

      stat = nf90_enddef(nc%id)
      NF_ERR('nf90_enddef','')

      nc%npar = 0

      if (associated(ipos) .and. .not.l_keep_sorting) deallocate(ipos)
    end if

    if (l_write .and. nprof > 0) then
      if (nc%id < 0) then
        msg = 'output file was not opened before writing.'
        call warning(trim(msg))
        if (present(status)) then
          status = stat
          return
        else
          call finish(trim(proc), trim(msg))
        end if
      end if

      if (present(offset)) then
        istart = offset + 1
      else
        istart = 1
      end if
      if (present(par)) then
        param = par(1:min(16,len_trim(par)))
      else
        param = 'fg'
      end if
      ipar = -1
      do i = 1, nc%npar
        if (param == nc%params(i)) then
          ipar = i
          exit
        end if
      end do
      if (ipar < 0) then
        if (nc%npar >= mx_par) then
          msg = 'too many parameters in write_rtovp'
          call warning(trim(msg))
          if (present(status)) then
            status = stat
            return
          else
            call finish(trim(proc), trim(msg))
          end if
        end if
        nc%npar = nc%npar + 1
        ipar = nc%npar
        nc%params(ipar) = param
        stat = nf90_put_var(nc%id, nc%v_par%id,  param, start=(/1,ipar/), count=(/len(param),1/))
        NF_ERR('nf90_put_var','parameter')
      end if

      ! Prepare sorted output
      if (lsorted .and. .not.associated(ipos)) then
        do iset = 1, nset
          if (np(iset) > 0 .and. .not.associated(r_arr(iset)%i_reprt)) lsorted = .false.
        end do
        if (.not.lsorted) then
          write(stderr, '(A)') '*** WARNING: sorted output into '//&
               &trim(nc%fname)//' is not possible because of missing metadata.'
        else
          allocate(s(nprof), nc%ipos(nprof))
          ipos => nc%ipos
          i = 0
          do iset = 1, nset
            if (np(iset) > 0) then
              s(i+1:i+np(iset))% iprof = (/ (j, j=i+1,i+np(iset)) /)
              s(i+1:i+np(iset))% reprt = r_arr(iset)%i_reprt(1:np(iset))
              s(i+1:i+np(iset))% ibox  = r_arr(iset)%i_box(1:np(iset))
              i = i + np(iset)
            end if
          end do
          ! Sort according to report index
          call sort_reports (s(1:nprof))
          ! Fill ipos array
          ipos(s(1:nprof)%iprof) = (/ (i, i=1,nprof) /)
          deallocate(s)
        end if
      end if

      ! --- 0d data
      if (ipar==1 .and. l_cld) then
!        stat = nf90_put_var(nc%id, nc%v_cld_ind%id, nc%cld_mask(:))
!        NF_ERR('nf90_put_var','cld_ind')
      end if

      ! --- 1d integer data
      allocate(i1dat(nprof))
#define WR_VAR(NAM,NCV) do i=1,nset; if (associated(r_arr(i)%NAM)) i1(i)%d=>r_arr(i)%NAM; enddo; call write_var_i1(nc%NCV)
      if (ipar == 1) then
        WR_VAR(obsnum,v_id)
        WR_VAR(time,v_time)
      end if
      deallocate(i1dat)

      ! --- 1d real data
      allocate(r1dat(nprof))
#undef WR_VAR
#define WR_VAR(NAM,NCV) do i=1,nset; if (associated(r_arr(i)%NAM)) r1(i)%d=>r_arr(i)%NAM; enddo; call write_var_r1(nc%NCV)
      if (ipar == 1) then
        WR_VAR(dlat,v_lat)
        WR_VAR(dlon,v_lon)
        WR_VAR(tsm_fg, v_tsm)
      end if
      WR_VAR(ts_fg,v_ts)
      WR_VAR(ps_fg,v_ps)
      do i=1,nset
        r => r_arr(i)
        if (associated(r%v10_abs_fg)) then
          r1(i)%d => r%v10_abs_fg
        else if (associated(r%v10_fg) .and. associated(r%u10_fg)) then
          allocate(r1(i)%d(np(i)))
          r1(i)%alloc = .true.
          r1(i)%d(1:np(i)) = ws(r%u10_fg(1:np(i)),r%v10_fg(1:np(i)))
        end if
      enddo
      call write_var_r1(nc%v_v10m)
      WR_VAR(cld_top,v_cldtop)
      if (.not.l_cld) then
        do i=1,nset
          if (associated(r_arr(i)%cld_frc)) then
            r1(i)%d => r_arr(i)%cld_frc(1,:)
          end if
        end do
!        call write_var_r1(nc%v_cldfrac)
      end if
      WR_VAR(snw_frc,v_snwfrac)
      deallocate(r1dat)

      ! --- 1d char data
      if (ipar == 1) call write_statid

      ! --- 2d real data
      allocate(r2dat(nlev,nprof))
#undef WR_VAR
#define WR_VAR(NAM,NCV) do i=1,nset; if (associated(r_arr(i)%NAM)) r2(i)%d=>r_arr(i)%NAM; enddo; call write_var_r2(nc%NCV)
      WR_VAR(t_fg,v_t)
      WR_VAR(q_fg,v_q)
      if (ipar == 1) then
        if (lev_mode <= 0) then
          do i = 1, nset
            if (associated(r_arr(i)%p)) then
              stat = nf90_put_var(nc%id, nc%v_level%id, r_arr(i)%p(:,1))
              NF_ERR('nf90_put_var','statid')
              exit
            end if
          end do
        else
          WR_VAR(p,v_level)
        end if
      end if
#define WR_TRG(NAM,ISP,NCV) do i=1,nset; if (r_arr(i)%ISP > 0) r2(i)%d=>r_arr(i)%trg(r_arr(i)%ISP,:,:); enddo; call write_var_r2(nc%NCV)
      if (l_o3 .and. trim(param) == 'fg') then
        WR_TRG(o3,i_o3,v_o3)
      end if
      if (l_co2 .and. trim(param) == 'fg') then
        WR_TRG(co2,i_co2,v_co2)
      end if
      deallocate(r2dat)
      if (l_pc) then
        allocate(r2dat(npc,nprof))
        WR_VAR(emis_pc,v_emis_pc)
        deallocate(r2dat)
      end if
      if (l_cld) then
        allocate(r2dat(nlev-1,nprof))
   !     WR_VAR(cld_frc,v_cldfrac)
        deallocate(r2dat)
      end if

      ! --- 3d real data
#undef WR_VAR
#define WR_VAR(NAM,NCV) do i=1,nset; if (associated(r_arr(i)%NAM)) r3(i)%d=>r_arr(i)%NAM; enddo; call write_var_r3_cld(nc%NCV)
      if (l_cld) then
        allocate(r3dat(ncld,nlev-1,nprof))
   !     WR_VAR(cld_fg,v_cld)
        deallocate(r3dat)
      end if
      if (l_H .or. l_K .or. l_pc .or. l_HBH) then
        ! 2dim int
        allocate(i2dat(nchan,nprof))
        do i=1,nset
          nch = r_arr(i)% i% n_chan
          allocate(i2(i)%d(nch,np(i)))
          i2(i)%alloc = .true.
          i2(i)%d(1:nch,1:np(i)) = spread(r_arr(i)% i% chan, 2, np(i))
        end do
        call write_var_i2_chan(nc%v_chan)
        do i=1,nset
          nch = r_arr(i)%i%n_chan
          allocate(i2(i)%d(nch,np(i)))
          do j=1,r_arr(i)%i%n_instr
            idum(r_arr(i)%i%o_ch_i(j)+1:r_arr(i)%i%o_ch_i(j)+r_arr(i)%i%n_ch_i(j)) = &
                 r_arr(i)%i%instr_wmo(j)
          end do
          i2(i)%alloc = .true.
          i2(i)%d(1:nch,1:np(i)) = spread(idum(1:nch), 2, np(i))
        end do
        call write_var_i2_chan(nc%v_instr)
        deallocate(i2dat)
      end if

      if (l_H .or. l_K .or. l_B .or. l_HBH) then
        ! 3dim real
#undef WR_VAR
#define WR_VAR(NAM,NCV) do i=1,nset; if (associated(r_arr(i)%NAM)) r3(i)%d=>r_arr(i)%NAM; enddo; call write_var_r3_chan(nc%NCV)
        allocate(r3dat(nchan,nlev,nprof))
        if (l_K) then
          WR_VAR(K_t,v_K_t)
          WR_VAR(K_q,v_K_q)
        end if
        if (l_H) then
          WR_VAR(H_t,v_H_t)
          WR_VAR(H_q,v_H_q)
          if (l_pc) then
            deallocate(r3dat)
            allocate(r3dat(nchan,npc,nprof))
            WR_VAR(H_em_pc,v_H_em_pc)
            call write_var_r3_chan(nc%v_H_em_pc)
          endif
        end if
        deallocate(r3dat)
        if (l_B) then
#undef WR_VAR
#define WR_VAR(NAM,NCV) do i=1,nset; if (associated(r_arr(i)%NAM)) r3(i)%d=>r_arr(i)%NAM; enddo; call write_var_r3(nc%NCV)
          allocate(r3dat(nlev,nlev,nprof))
          WR_VAR(B_tt,v_B_tt)
          WR_VAR(B_qq,v_B_qq)
          WR_VAR(B_tq,v_B_tq)
          deallocate(r3dat)
        end if
        if (l_HBH) then
#undef WR_VAR
#define WR_VAR(NAM,NCV) do i=1,nset; if (associated(r_arr(i)%NAM)) r3(i)%d=>r_arr(i)%NAM; enddo; call write_var_r3_chan2(nc%NCV)
          allocate(r3dat(nchan,nchan,nprof))
          WR_VAR(HBHR,v_HBHR)
          WR_VAR(R,v_R)
          deallocate(r3dat)
        end if
      end if
      if (l_H .or. l_K .or. l_pc) then
        ! 2dim real
#undef WR_VAR
#define WR_VAR(NAM,NCV) do i=1,nset; if (associated(r_arr(i)%NAM)) r2(i)%d=>r_arr(i)%NAM; enddo; call write_var_r2_chan(nc%NCV)
        allocate(r2dat(nchan,nprof))
        if (l_H) then
          WR_VAR(H_ts,v_H_ts)
          WR_VAR(H_ps,v_H_ps)
        end if
        if (l_K) then
          WR_VAR(K_ts,v_K_ts)
          WR_VAR(K_ps,v_K_ps)
        end if
        if (l_pc) then
          WR_VAR(emiss,v_emiss)
        end if
        deallocate(r2dat)
#undef WR_VAR
      end if

    end if

    if (.not.l_keep_sorting .and. associated(ipos)) deallocate(ipos)

    if (l_close) then
      stat = nf90_close(nc%id)
      NF_ERR('nf90_close',trim(fname))
      nc%set_id = -999
      if (all(nc_arr(1:n_nc)%set_id == -999)) then
        n_nc = 0
        do i = 1, n_nc
          if (associated(nc_arr(i)%cld_mask)) deallocate(nc_arr(i)%cld_mask)
        end do
        deallocate(nc_arr)
      end if
      nc => null()
    end if

  contains

#define TYPE i1
#define write_var_TYPE write_var_i1
#define TYPEdat i1dat
#include "write_rtovp_var.incf"
#undef  TYPE
#undef  write_var_TYPE
#undef  TYPEdat
    !--------------------------------------------------------------------------
#define TYPE r1
#define write_var_TYPE write_var_r1
#define TYPEdat r1dat
#include "write_rtovp_var.incf"
#undef  TYPE
#undef  write_var_TYPE
#undef  TYPEdat
    !--------------------------------------------------------------------------
#define TYPE r2
#define write_var_TYPE write_var_r2
#define TYPEdat r2dat
#define DIMS2
#include "write_rtovp_var.incf"
#undef  DIMS2
#undef  TYPE
#undef  write_var_TYPE
#undef  TYPEdat
    !--------------------------------------------------------------------------
#define TYPE r3
#define write_var_TYPE write_var_r3
#define TYPEdat r3dat
#define DIMS3
#include "write_rtovp_var.incf"
#undef  DIMS3
#undef  TYPE
#undef  write_var_TYPE
#undef  TYPEdat
    !--------------------------------------------------------------------------
#define TYPE r3
#define write_var_TYPE write_var_r3_cld
#define TYPEdat r3dat
#define DIMS3
#define DAT_MASK nc%cld_mask(:),:,
#include "write_rtovp_var.incf"
#undef  DIMS3
#undef  TYPE
#undef  write_var_TYPE
#undef  TYPEdat
    !--------------------------------------------------------------------------
#define TYPE i2
#define write_var_TYPE_chan write_var_i2_chan
#define TYPEdat i2dat
#define DIMCH
#include "write_rtovp_var_chan2.incf"
#undef  DIMCH
#undef  TYPE
#undef  write_var_TYPE_chan
#undef  TYPEdat
    !--------------------------------------------------------------------------
#define TYPE r2
#define write_var_TYPE_chan write_var_r2_chan
#define TYPEdat r2dat
#include "write_rtovp_var_chan.incf"
#undef  TYPE
#undef  write_var_TYPE_chan
#undef  TYPEdat
    !--------------------------------------------------------------------------
#define TYPE r3
#define write_var_TYPE_chan write_var_r3_chan
#define TYPEdat r3dat
#define DIMS3
#include "write_rtovp_var_chan.incf"
#undef  DIMS3
#undef  TYPE
#undef  write_var_TYPE_chan
#undef  TYPEdat
    !--------------------------------------------------------------------------
#define TYPE r3
#define write_var_TYPE_chan write_var_r3_chan2
#define TYPEdat r3dat
#define DIMS3
#define DIMCH
#include "write_rtovp_var_chan2.incf"
#undef  DIMCH
#undef  DIMS3
#undef  TYPE
#undef  write_var_TYPE_chan
#undef  TYPEdat

    subroutine write_statid
      integer :: i, n, iset
      character(len=8) :: c1dat(nprof)
      i = 0
      do iset = 1, nset
        n = r_arr(iset)% n_rec
        if (lsorted) then
          c1dat(ipos(i+1:i+n)) = satname(r_arr(iset)%i%satid)
        else
          c1dat(i+1:i+n) = satname(r_arr(iset)%i%satid)
        end if
        i = i + n
      end do
      stat = nf90_put_var(nc%id, nc%v_statid%id, c1dat(1:nprof), start=(/1, istart/), count=(/8, nprof/))
      NF_ERR('nf90_put_var','statid')
    end subroutine write_statid


    subroutine def_var(v, dims)
      type(t_nc_var), intent(inout) :: v
      integer,        intent(in)    :: dims(:)
#undef NF_ERR
#define NF_ERR(nproc, text) if (stat/=0) then ; call nf_err(nproc, text) ; return ; endif
      stat = nf90_def_var(nc%id, trim(v%name) , v%type, dims, v%id)
      NF_ERR('nf90_def_var',trim(v%name))
      if (v%unit /= '') then
        stat = nf90_put_att(nc%id, v%id, 'units', v%unit)
        NF_ERR('nf90_put_att',trim(v%name)//': units')
      end if
      if (v%long_name /= '') then
        stat = nf90_put_att(nc%id, v%id, 'long_name', v%long_name)
        NF_ERR('nf90_put_att',trim(v%name)//': long_name')
      end if
      if (v%type == NF_FLOAT) then
        stat = nf90_put_att(nc%id, v%id, '_FillValue', NF_FILL_FLOAT)
        NF_ERR('nf90_put_att',trim(v%name)//': _FillValue')
      elseif (v%type == NF_INT) then
        stat = nf90_put_att(nc%id, v%id, '_FillValue', NF_FILL_INT)
        NF_ERR('nf90_put_att',trim(v%name)//': _FillValue')
      end if
    end subroutine def_var

    subroutine nf_err(nproc, msg)
      character(len=*), intent(in)           :: nproc
      character(len=*), intent(in), optional :: msg

      if (stat /= 0) then
        call nf_error(stat, trim(proc)//', '//trim(nproc), msg)
        if (present(status)) then
          status = stat
        else
          call finish(trim(proc), trim(nproc))
        end if
      end if

    end subroutine nf_err

    elemental real(wp) function ws(u,v)
      real(wp), intent(in) :: u
      real(wp), intent(in) :: v
      ws = sqrt(u**2 + v**2)
    end function ws

    subroutine sort_reports (s)
      type(t_sort), intent(inout) :: s(:)
      !------------------------------
      ! Use efficient sorting routine
      !------------------------------
      integer     :: idx(size (s))      ! Permutation
      integer(i8) :: key(size (s))      ! Sort key

      !-----------------------------------
      ! Sort according to ibox, then reprt
      !-----------------------------------
      key  = s(:)% reprt + s(:)% ibox * 2_i8 ** 32
      idx  = index (key)
      s(:) = s(idx(:))

    end subroutine sort_reports

  end subroutine write_rtovp


  subroutine lev2chan(tot, overc, valid, l2c, l2c_type, rel_lim, abs_lim, ipr)
    real(wp), intent(in)            :: tot(:,:)      ! (chan, prof)
    real(wp), intent(in)            :: overc(:,:,:)  ! (lev, chan, prof)
    logical,  intent(in)            :: valid(:,:)    ! mask
    real(sp), intent(out)           :: l2c(:,:)      ! result: assigned level
    integer,  intent(in)            :: l2c_type
    real(wp), intent(in),  optional :: rel_lim
    real(wp), intent(in),  optional :: abs_lim
    integer,  intent(in),  optional :: ipr(:)

    real(wp),parameter :: cloudlimit = 0.01_wp

    real(wp) :: tmp(size(tot,1),size(tot,2)), tmp_aux (size(tot,1),size(tot,2))
    real(wp) :: r_lim, a_lim, fac
    integer  :: n_chan, n_prof
    integer  :: n_lev
    integer  :: i,k,l
    integer  :: ideb, idc
    !integer :: ii(size(tot,1)), nc_, ii_, ip, j

    r_lim = cloudlimit
    if (present(rel_lim)) then
      if (rel_lim > 0._wp) r_lim = rel_lim
    end if
    a_lim = -1._wp
    if (present(abs_lim)) then
      if (abs_lim > 0._wp) then
        a_lim = abs_lim
        r_lim = -1._wp
      end if
    end if

    n_chan = size(tot,1)
    n_prof = size(tot,2)
    n_lev  = size(overc,1)


    if (present(ipr)) then
      do i = 1, size(ipr)
        ideb = ipr(i)
        if (ideb <= 0 .or. ideb > n_prof) cycle
        do idc = 1, n_chan
          write(usd,*) 'debug_spot lev2chan tot',ideb,idc,tot(idc,ideb),r_lim*tot(idc,ideb)
          do k = 1, n_lev
            write(usd,*) 'debug_spot lev2chan overc',ideb,idc,k,overc(k,idc,ideb),tot(idc,ideb)
          end do
        end do
      end do
    end if

    select case(l2c_type)
    case(1)
      ! level to channel assignment as proposed by the manual of the
      ! cloud detection software:
      l2c(:,:) = -1._sp
      if (r_lim > 0._wp) then
        do l = 1, n_lev
          do i = 1, n_prof
            do k = 1, n_chan
              if (valid(k, i).and.(l2c(k,i) < 0._sp)) then
                if (abs( (tot(k,i) - overc(l, k, i)) / tot(k,i)) < r_lim) &
                     l2c(k,i) = real(l, kind=sp)
              end if
            end do
          end do
        end do
      else !absolute limit
        do l = 1, n_lev
          do i = 1, n_prof
            do k = 1, n_chan
              if (valid(k, i).and.(l2c(k,i) < 0._sp)) then
                if (abs(tot(k,i) - overc(l, k, i)) < a_lim) &
                     l2c(k,i) = real(l, kind=sp)
              end if
            end do
          end do
        end do
      end if
      where (valid(:,:) .and. l2c(:,:) < 0._sp) &
        l2c(:,:) = real(n_lev+1, sp)                ! Surface level

    case(2)

      ! level to channel assignment as proposed by the manual of the
      ! cloud detection software but with the modification that the level
      ! is not assigned if the scaled difference is once more larger than
      ! r_lim below.
      ! This is needed since due to the fact that on the top of the model
      ! region an ecmwf short term forecast is put the scaled difference
      ! is a little bit strange in the transition region:

      l2c(:,:) = real(n_lev+1, sp)
      if (r_lim > 0._wp) then
        where(valid(:,:))
          tmp_aux(:,:) = abs(r_lim * tot(:,:))
        elsewhere
          tmp_aux(:,:) = 0._wp
          tmp    (:,:) = 0._wp
        end where

        do l = n_lev, 1, -1
          where(valid(:,:)) tmp(:,:) = abs(tot(:,:) - overc(l,:,:))
          where(valid(:,:) .and. l2c(:,:)==real(l+1, sp) .and. tmp(:,:) < tmp_aux(:,:))
            l2c(:,:) = real(l, sp)
          end where
        end do
      else !absolute limit
        where(.not.valid(:,:))
          tmp    (:,:) = 0._wp
        end where

        do l = n_lev, 1, -1
          where(valid(:,:)) tmp(:,:) = abs(tot(:,:) - overc(l,:,:))
          where(valid(:,:) .and. l2c(:,:)==real(l+1, sp) .and. tmp(:,:) < a_lim)
            l2c(:,:) = real(l, sp)
          end where
        end do
      end if

      ! More comprehensible but slower code:
      ! l2c(:,:) = -1._sp
      ! do l = 1, n_lev
      !   do i = 1, n_prof
      !     do k = 1, n_chan
      !       if ((l2c(k,i) < 0._sp) .and.                                    &
      !           all(abs((tot(k,i) - overc(l,k,i)) / tot(k,i)) < r_lim))&
      !         l2c(k,i) = real(l, sp)
      !     end do
      !   end do
      ! end do
      ! where (valid(:,:) .and. l2c(:,:) < 0._sp) &
      !   l2c(:,:) = real(n_lev+1, sp)               ! Surface level

    case(3)
      ! As case 2, but the result is a real number, i.e. the associated level
      ! is somewhere between two RTTOV levels.

      if (r_lim > 0._wp) then
        where (valid(:,:) .and. tot(:,:) > 0._wp)
          tmp(:,:) = abs((tot(:,:) - overc(n_lev,:,:)) / tot(:,:))
        elsewhere
          l2c(:,:) = real(n_lev+1, sp)
          tmp(:,:) = 0._wp
        end where

        fac = maxval(tmp(:,:), mask=valid(:,:))
        if (fac > 0._wp) then
          fac = 1._wp/fac
        else
          fac = 1._wp
        end if

        where (valid(:,:))
          where (tmp(:,:) >= r_lim)
            ! Channel influenced by cloud in lowest level -> below surface level
            ! We add tmp-r_lim in order to have a (hopefully) unique order of the
            ! surface channels
            l2c(:,:) = real(n_lev, sp) + tmp(:,:) * fac
          elsewhere
            l2c(:,:) = 1._sp
            tmp_aux (:,:) = tmp(:,:)
          end where
        end where

        do l = n_lev-1, 1, -1
          where (valid(:,:) .and. tot(:,:) > 0._wp)
            where (l2c(:,:) == 1._sp)
              tmp(:,:) = abs((tot(:,:) - overc(l,:,:)) / tot(:,:))
              where (tmp(:,:) >= r_lim)
                l2c(:,:) = real(( (tmp(:,:) - r_lim       ) * (l+1) +     &
                                  (r_lim    - tmp_aux(:,:)) * l       ) / &
                                (tmp(:,:) - tmp_aux (:,:)), kind=sp)
              end where
              tmp_aux (:,:) = tmp(:,:)
            end where
          end where
        end do

      else !absolute limit
        where (valid(:,:))
          tmp(:,:) = abs(tot(:,:) - overc(n_lev,:,:))
        elsewhere
          l2c(:,:) = real(n_lev+1, sp)
          tmp(:,:) = 0._wp
        end where

        where(valid(:,:))
          where (tmp(:,:) >= a_lim)
            ! Channel influenced by cloud in lowest level -> surface channel
            l2c(:,:) = real(n_lev, sp) + (tmp(:,:) - a_lim) /tot(:,:)
          elsewhere
            l2c(:,:) = 1._sp
            tmp_aux (:,:) = tmp(:,:)
          end where
        end where

        do l = n_lev-1, 1, -1
          where (valid(:,:))
            where (l2c(:,:) == 1._sp)
              tmp(:,:) = abs(tot(:,:) - overc(l,:,:))
              where (tmp(:,:) >= a_lim)
                l2c(:,:) = real(( (tmp(:,:) - a_lim       ) * (l+1) +     &
                                  (a_lim    - tmp_aux(:,:)) * l       ) / &
                                (tmp(:,:) - tmp_aux (:,:)), kind=sp)
              end where
              tmp_aux (:,:) = tmp(:,:)
            end where
          end where
        end do
      end if

      ! ! Test ambiguity in result
      ! do ip = 1, n_prof
      !   ii = (/ (i, i=1,n_chan) /)
      !   nc_ = count(valid(:,ip))
      !   ii(1:nc_) = pack(ii, mask=valid(:,ip))
      !   do i = nc_-1,1,-1
      !     do j = 1, i
      !       if (l2c(ii(j),ip) > l2c(ii(j+1),ip)) then
      !         ii_ = ii(j)
      !         ii(j) = ii(j+1)
      !         ii(j+1) = ii_
      !       end if
      !     end do
      !   end do
      !   do i = 1, nc_-1
      !     if (l2c(ii(i),ip) == l2c(ii(i+1),ip)) then
      !       write(0,*) 'Equal values in l2c',ip,i,ii(i),l2c(ii(i),ip),i+1,ii(i+1),l2c(ii(i+1),ip)
      !       write(0,*) 'bt',tot(ii(i),ip),overc(n_lev,ii(i),ip),tot(ii(i+1),ip),overc(n_lev,ii(i+1),ip)
      !     end if
      !   end do
      ! end do


    case default
      l2c(:,:) = -1._sp
    end select

    if (present(ipr)) then
      do i = 1, size(ipr)
        ideb = ipr(i)
        if (ideb <= 0 .or. ideb > n_prof) cycle
        do idc = 1, n_chan
          write(usd,*) 'debug_spot lev2chan result',idc,l2c(idc,ideb)
        end do
      end do
    end if

  end subroutine lev2chan


  subroutine lev2p(p, l2p)
    real(wp), intent(in)    :: p(:)
    real(wp), intent(inout) :: l2p(:)

    real(wp) :: frac
    integer  :: nl
    integer  :: l0
    integer  :: i

    nl = size(p)
    do i = 1,size(l2p)
      if (l2p(i) < 1) then
        l2p(i) = log(p(1)) + (l2p(i) - 1) * (log(p(2)) - log(p(1)))
      elseif (l2p(i) >= nl) then
        l2p(i) = min(l2p(i),3._wp*nl) ! Avoid extreme values
        l2p(i) = log(p(nl)) + (l2p(i) - nl) * (log(p(nl)) - log(p(nl-1)))
      else
        l0 = int(l2p(i))
        frac = l2p(i) - l0
        l2p(i) = frac * log(p(l0+1)) + (1._wp -frac) * log(p(l0))
      end if
      l2p(i) = exp(l2p(i))
    end do

  end subroutine lev2p


  !> \todo Does not work in COSMO
  subroutine nf_error(stat, proc, msg)
    integer,          intent(in)           :: stat
    character(len=*), intent(in)           :: proc
    character(len=*), intent(in), optional :: msg

    write(stderr, '(A)') '*** NetCDF error (routine '//trim(proc)//') : '//&
         &nf90_strerror(stat)
    if (present(msg)) write(stderr, '(A)') '*** I/O error '//trim(msg)

  end subroutine nf_error


  function instr_type(rs, instr, chan, ci) result(typ)
    integer                               :: typ
    type(t_rad_set), intent(in)           :: rs
    integer,         intent(in)           :: instr
    integer,         intent(in), optional :: chan
    integer,         intent(in), optional :: ci(:)
    integer :: n0, n1
    typ = 0
    if (mw_instr(rs%instr(instr))) typ = ibset(typ, ITYP_MW)
    if (ir_instr(rs%instr(instr))) typ = ibset(typ, ITYP_IR)
    if (associated(rs%iflag)) then
      if (present(chan)) then
        if (rs%iflag(chan) == 1) typ = ibset(0, ITYP_VIS)
      elseif (present(ci)) then
        if (all(rs%iflag(ci) == 1)) then
          typ = ibset(0, ITYP_VIS)
        elseif (any(rs%iflag(ci) == 1)) then
          typ = ibset(typ, ITYP_VIS)
        end if
      else
        n0 = rs%o_ch_i(instr) + 1
        n1 = rs%o_ch_i(instr) + rs%n_ch_i(instr)
        if (all(rs%iflag(n0:n1) == 1)) then
          typ = ibset(0, ITYP_VIS)
        elseif (any(rs%iflag(n0:n1) == 1)) then
          typ = ibset(typ, ITYP_VIS)
        end if
      end if
    elseif (.not.present(chan).and.vis_instr(rs%instr(instr))) then
      typ = ibset(typ, ITYP_VIS)
    end if
    if (typ == 0) then
      write(0,*) 'instr',instr,rs%instr(instr),rs%satid
      call finish('instr_type', 'failed to determine instrument type')
    end if
  end function instr_type

  !> Merge observations in t_radv (see also explanation of superob_tovs@mo_tovs.f90)
  subroutine thin_superob_radv(r, ind_, id, ibase, spot_mask, status)
    type(t_radv), intent(inout)           :: r
    integer,      intent(in)              :: ind_(:)  ! indices of obs to be merged
    integer,      intent(in)              :: id(:)    ! obs-IDs of obs to be merged
    integer,      intent(in)              :: ibase(:) ! "ibase" array from superob_tovs
                                                      ! see explanantion there
    logical,      intent(out)             :: spot_mask(:) ! Mask of "base" obs, i.e. obs
                                                          ! to be kept later
    integer,      intent(out)  , optional :: status

    integer, allocatable :: ib(:)     ! "base" FOVs of the so-boxes
    integer, allocatable :: nobs(:)   ! number of FOVs per so-box
    integer, allocatable :: ii(:)     ! so-box number that corresponds to entries in ibase/id/ind_
    integer, allocatable :: ind(:,:)  ! for all obs in all so-boxes the indices in t_radv
    integer, allocatable :: iso(:)    ! indices of obs to be merged
    logical, allocatable :: mask(:)   ! auxiliary array for merging

    integer :: i,j,k,n, m
    integer :: nbox, ibox, mx_nobs
    integer :: nv
    logical :: ldeb

    status = 0

    n = size(ibase)

    ! 1. Reorganize ibase/id arrays into more convenient arrays
    ! Set up number of boxes ("nbox"), box size ("nobs"), ibase value of each box ("ib")
    ! and box for each obs ("ii")
    allocate(ib(n), nobs(n), ii(n))
    nbox = 0
    do i = 1, n
      ibox = 0
      do j = 1, nbox
        if (ib(j) == abs(ibase(i))) then
          ibox = j
          exit
        end if
      end do
      if (ibox <= 0) then
        nbox = nbox + 1
        nobs(nbox) = 0
        ib(nbox) = abs(ibase(i))
        ibox = nbox
      end if
      nobs(ibox) = nobs(ibox) + 1
      ii(i) = ibox
    end do
    mx_nobs = maxval(nobs(1:nbox))
    ! Set up the 2d index array
    allocate(ind(nbox,mx_nobs))
    nobs = 0
    ib = 0
    do i = 1, n
      ldeb = any(r%ideb(:) == ibase(i))
      ibox = ii(i)
      nobs(ibox) = nobs(ibox) + 1
      ! Get index of obs in t_radv
      if (ind_(i) > 0) then
        ind(ibox, nobs(ibox)) = ind_(i)
        ! Safeguard: check obsnum
        if (r% obsnum(ind_(i)) /=  id(i)) then
          write(0,*) 'inconsistent indices in thin_superob_radv: i=',i,' ibox=',ibox,&
               ' nobs=',nobs(ibox),' ind=',ind_(i),' id=',id(i),' obsnum=',r% obsnum(ind_(i))
          status = 1
          return
        end if
      else
        ! Obs is not originating from this PE and was sent from another PE
        ! So we have to search for matching obsnum
        k = 0
        do j = 1, r%n_rec
          if (r%obsnum(j) == id(i)) then
            k = j
            exit
          end if
        end do
        if (k <= 0) then
          write(0,*) 'Did not find obsnum in thin_superob_radv: i=',i,' ibox=',ibox,&
               ' nobs=',nobs(ibox),' id=',id(i)
          status = 2
          return
        end if
        ind(ibox, nobs(ibox)) = k
      end if
      ! Store index of base FOV of the box in ib
      if (ibase(i) == ind_(i)) ib(ibox) = ibase(i)
    end do

    ! 2. Superobbing, i.e. combine quantities of FOVs to be superobbed
    allocate(iso(mx_nobs), mask(mx_nobs))
    spot_mask(:) = .false.
    do ibox = 1, nbox
      if (ib(ibox) <= 0) then
        write(0,*) 'no base FOV found for box',ibox,nobs(ibox),ind(ibox,1:nobs(ibox))
        status = 3
        return
      end if

      ! TODO: safeguard check: distances between FOVs

      i = ib(ibox)

      ldeb = any(r%ideb(:) == i)

      m = nobs(ibox)
      iso(1:m) = ind(ibox,1:m)
      spot_mask(i) = .true.
      !write(0,*) 'box',ibox,i,m,iso(1:m),r%n_rec
      if (ldeb) write(usd,*) 'debug_spot superob',i,m,iso(1:m)
      if (ldeb) write(usd,*) 'debug_spot superob ids',i,m,r%obsnum(iso(1:m))
      if (r%i% gopts% thin_superob_mode == 1 ) then
        ! Just thinning, nothing to do here
      else if (r%i% gopts% thin_superob_mode == 2 ) then
        if (associated(r%dlat    )) call angle_1     (r%dlat    )
        if (associated(r%dlon    )) call angle_1     (r%dlon    )

        if (associated(r%time    )) call i_avg_1     (r%time    )
        if (associated(r%nwc_flg )) call ior_1       (r%nwc_flg )
        if (associated(r%scanl   )) call i_avg_1     (r%scanl   ) ! TODO: do we want to do this?
        if (associated(r%fov     )) call i_avg_1     (r%fov     ) ! TODO: do we want to do this?
        if (associated(r%mdlsfc  )) call ior_1       (r%mdlsfc  )
        if (associated(r%stzen   )) call r_absavg_1  (r%stzen   )
        if (associated(r%stazi   )) call angle_1     (r%stazi   )
        if (associated(r%sunzen  )) call r_absavg_1  (r%sunzen  )
        if (associated(r%sunazi  )) call angle_1     (r%sunazi  )
        if (associated(r%orb_phase))call angle_1     (r%orb_phase)
        if (associated(r%instr_temp))call r_avg_2    (r%instr_temp)

        if (associated(r%bt_obs  )) call r_avg_2     (r%bt_obs  )
        if (associated(r%bt_bcor )) call r_avg_2     (r%bt_bcor )
        if (associated(r%bcor_   )) call r_avg_2     (r%bcor_   )
        if (associated(r%bt_fg   )) call r_avg_2     (r%bt_fg   )
        if (associated(r%bt_fg_cs)) call r_avg_2     (r%bt_fg_cs)
        if (associated(r%state   )) call state_2     (r%state   )

        if (associated(r%flags   )) call ior_2       (r%flags   )
        if (associated(r%emiss   )) call r_avg_2     (r%emiss   )
        if (associated(r%shgt    )) call r_avg_2     (r%shgt    )
        if (associated(r%stype   )) call stype_2     (r%stype   )
        if (associated(r%r_state )) call state_1     (r%r_state )
      else

          write(0,*) 'Namelist variable thin_superob_mode has invalid value: ', r%i%gopts%thin_superob_mode
          status = 4
          return

      end if

      ! Please, feel free to add further variables ...

    end do

  contains

    subroutine r_absavg_1(x)
      real(kind=wp), intent(inout) :: x(:)
      nv = count(x(iso(1:m)) /= rinvalid)
      if (nv > 0) then
        x(i) = sum(abs(x(iso(1:m))), mask=(x(iso(1:m))/=rinvalid)) / (1.*nv)
      else
        x(i) = rinvalid
      end if
    end subroutine r_absavg_1

    subroutine r_avg_1(x)
      real(kind=wp), intent(inout) :: x(:)
      nv = count(x(iso(1:m)) /= rinvalid)
      if (nv > 0) then
        x(i) = sum(x(iso(1:m)), mask=(x(iso(1:m))/=rinvalid)) / (1.*nv)
      else
        x(i) = rinvalid
      end if
    end subroutine r_avg_1

    subroutine r_avg_2(x)
      real(kind=wp), intent(inout) :: x(:,:)
      integer :: k
      do k = lbound(x,1),ubound(x,1)
        call r_avg_1(x(k,:))
      end do
    end subroutine r_avg_2


    subroutine i_avg_1(x)
      integer, intent(inout) :: x(:)
      nv = m
      if (nv > 0) then
        x(i) = int(sum(x(iso(1:m)), mask=(x(iso(1:m))/=rinvalid)) / (1.*nv))
      else
        x(i) = 0
      end if
    end subroutine i_avg_1

    subroutine angle_1(x)
      real(kind=wp), intent(inout) :: x(:)

      x(i) = avg_angle(x(iso(1:m)))

    end subroutine angle_1


    subroutine state_1(x)
      integer, intent(inout) :: x(:)
      ! Since STAT_PASSIVE and STAT_REJECTED art badly chosen, we have to modify
      ! the state variable
      integer, parameter :: st_map(0:15) = (/0,1,2,3,4,5,6,8,7,9,10,11,12,13,14,-1/)
      integer :: k, st, st_save

      st = 14 !ACCEPTED
      nv = 0
      st_save=x(i)
      do k = 1,m
        if (x(iso(k)) >= 0 .and. x(iso(k)) <= 15) then
          if (st_map(x(iso(k))) < st_map(st)) st = x(iso(k))
          nv = nv + 1
        end if
      end do
      if (nv > 0) then
        x(i) = st
      else
        x(i) = 4 ! STAT_DISMISS
      end if
      ! If a channel was not available in some spots, it might happen, that status=STAT_DEFAULT
      ! Since this state is not accepted later, we set a minimum state here (STAT_DISMISS)
      x(i) = max(x(i),4)
      !Workaround to avoid that single channel gets STAT_DISMISS in spot which is
      !not STAT_DISMISS
      if (x(i) == 4 .and. st_save > 4) then
        if (st_save == 6 .or. st_save == 7) then
          x(i) = 6
        else
          x(i) = 8
        end if
      end if

    end subroutine state_1


    subroutine state_2(x)
      integer, intent(inout) :: x(:,:)
      integer :: k
      do k = lbound(x,1),ubound(x,1)
        call state_1(x(k,:))
      end do
    end subroutine state_2


    subroutine ior_1(x)
      integer, intent(inout) :: x(:)
      integer :: k
      x(i) = x(iso(1))
      do k = 2,m
        x(i) = ior(x(i), x(iso(k)))
      end do
    end subroutine ior_1

    subroutine ior_2(x)
      integer, intent(inout) :: x(:,:)
      integer :: k
      do k = lbound(x,1),ubound(x,1)
        call ior_1(x(k,:))
      end do
    end subroutine ior_2

    subroutine stype_2(x)
      integer, intent(inout) :: x(:,:)
      integer :: k
      do k = 1, size(x,1)
        if (all(x(k,iso(1:m)) == 0)) then
          x(k,i) = 0
        elseif (all(x(k,iso(1:m)) == 2)) then
          x(k,i) = 2
        else
          x(k,i) = 1
        end if
      end do
    end subroutine stype_2

  end subroutine thin_superob_radv


  function avg_angle(x) result(a)
    real(kind=wp)             :: a
    real(kind=wp), intent(in) :: x(:)
    real(kind=wp) :: vx, vy
    integer :: i, nv
    ! Average 2d vector
    nv = 0
    vx = 0._wp
    vy = 0._wp
    do i = 1, size(x)
      if (x(i) /= rinvalid) then
        vx = vx + dcos(x(i) * deg2rad)
        vy = vy + dsin(x(i) * deg2rad)
        nv = nv + 1
      end if
    end do
    if (nv > 0) then
      vx = vx / (1._wp*nv)
      vy = vy / (1._wp*nv)
      a = datan2(vy,vx) * rad2deg
    else
      a = rinvalid
    end if
  end function avg_angle


  function amsua_chan_id(instr, ch) result(id)
    integer             :: id
    integer, intent(in) :: instr
    integer, intent(in) :: ch

    integer, pointer :: map(:)
    integer          :: i

    id = -1
    map => null()
    select case(instr)
    case(3)
      map => amsua_amsua_chans
    case(4,15)
      map => mhs_amsua_chans
    case(19)
      map => atms_amsua_chans
    case(73)
      map => mwhs2_amsua_chans
    case(132)
      map => mwts3_amsua_chans
    case default
      return
    end select
    do i = 1, size(map)
      if (map(i) == ch) then
        id = i
        return
      end if
    end do
  end function amsua_chan_id

  function amsub_chan_id(instr, ch) result(id)
    integer             :: id
    integer, intent(in) :: instr
    integer, intent(in) :: ch

    integer, pointer :: map(:)
    integer          :: i

    id = -1
    map => null()
    select case (instr)
    case(4,15)
       map => amsub_amsub_chans
    case(10)
       map => ssmis_amsub_chans
    case(19)
       map => atms_amsub_chans
    case(34)
       map => saphir_amsub_chans
    case(71)
       map => gmi_amsub_chans
    case(73)
       map => mwhs2_amsub_chans
    case default
       return
    end select

    do i = 1, size(map)
      if (map(i) == ch) then
        id = i
        return
      end if
    end do
  end function amsub_chan_id

  function irxsr_chan_id(instr, ch) result(id)
    integer             :: id
    integer, intent(in) :: instr
    integer, intent(in) :: ch

    integer, pointer :: map(:)
    integer          :: i

    id = -1
    map => null()
    select case (instr)
    case(20)
       map => mviri_irxsr_chans
    case(21)
       map => seviri_irxsr_chans
    case(22)
       map => goesim_irxsr_chans
    case(44,56)
       map => abi_irxsr_chans
    case(61)
       map => fci_irxsr_chans
    case default
       return
    end select

    do i = 1, size(map)
      if (map(i) == ch) then
        id = i
        return
      end if
    end do
  end function irxsr_chan_id

  function seviri_chan_id(instr, ch) result(id)
    integer             :: id
    integer, intent(in) :: instr
    integer, intent(in) :: ch

    integer, pointer :: map(:)
    integer          :: i

    id = -1
    map => null()
    select case (instr)
    ! case(20)
    !    map => mviri_seviri_chans
    case(21)
       map => sev_seviri_chans
    ! case(22)
    !    map => goesim_seviri_chans
    ! case(44,56)
    !    map => abi_seviri_chans
    case(61)
       map => fci_seviri_chans
    case default
       return
    end select

    do i = 1, size(map)
      if (map(i) == ch) then
        id = i
        return
      end if
    end do
  end function seviri_chan_id


  subroutine bcast_rad_set (s, source, comm)
    type (t_rad_set)      ,intent(inout) :: s
    integer               ,intent(in)    :: source
    integer   ,optional   ,intent(in)    :: comm

    integer :: lcom
    integer :: count = 0
    integer :: i
    integer :: n_pc, n_styp, n_prio, n_emis_opt, n_pr, n_tskin_opt
    integer :: m_styp, m_prio

    lcom = dace% comm ; if (present(comm)) lcom = comm

    if (count == 0) count = size(transfer(s,(/' '/)))

    if (dace% pe /= source) call destruct (s)
    call p_bcast_derivedtype (s, count, source, lcom)

    if (associated(s% emis_pc)) then
      if (dace% pe == source) n_pc = size(s% emis_pc,1)
      call p_bcast(n_pc, source, lcom)
    end if
    if (associated(s% emis_ind)) then
      if (dace% pe == source) then
        n_prio = size(s% emis_ind,1)
        n_styp = size(s% emis_ind,2)
      end if
      call p_bcast(n_prio, source, lcom)
      call p_bcast(n_styp, source, lcom)
    end if
    if (associated(s% tskin_ind)) then
      if (dace% pe == source) then
        n_prio = size(s% tskin_ind,1)
        n_styp = size(s% tskin_ind,2)
      end if
      call p_bcast(n_prio, source, lcom)
      call p_bcast(n_styp, source, lcom)
    end if
    if (associated(s% emis_opt)) then
      if (dace% pe == source) n_emis_opt = size(s% emis_opt)
      call p_bcast(n_emis_opt, source, lcom)
      ! It is crucial to allocate emis_opt with the given size, and not with
      ! s%n_emis_opt, since the size might be intenionally larger than s%n_emis_opt in order
      ! to add further entries later in update_emis_opt
    end if
    if (associated(s% tskin_opt)) then
      if (dace% pe == source) n_tskin_opt = size(s% tskin_opt)
      call p_bcast(n_tskin_opt, source, lcom)
      ! It is crucial to allocate tskin_opt with the given size, and not with
      ! s%n_tskin_opt, since the size might be intenionally larger than s%n_tskin_opt in order
      ! to add further entries later in update_tskin_opt
    end if

    if (dace% pe /= source) then
      if (associated(s% sensor_instr)) allocate(s% sensor_instr(s%n_sens      ))
      if (associated(s% chan        )) allocate(s% chan        (s%n_chan      ))
      if (associated(s% chidx       )) allocate(s% chidx       (s%n_chan      ))
      if (associated(s% iflag       )) allocate(s% iflag       (s%n_chan      ))
      if (associated(s% band        )) allocate(s% band        (s%n_chan      ))
      if (associated(s% flag        )) allocate(s% flag        (s%n_chan      ))
      if (associated(s% wavenum     )) allocate(s% wavenum     (s%n_chan      ))
      if (associated(s% emis_land   )) allocate(s% emis_land   (s%n_chan      ))
      if (associated(s% emis_snow   )) allocate(s% emis_snow   (s%n_chan      ))
      if (associated(s% emis_sice   )) allocate(s% emis_sice   (s%n_chan      ))
      if (associated(s% emis_pc     )) allocate(s% emis_pc     (n_pc, s%n_chan))
      if (associated(s% emis_opt    )) allocate(s% emis_opt    (n_emis_opt    ))
      if (associated(s% tskin_opt   )) allocate(s% tskin_opt   (n_tskin_opt    ))
      if (associated(s% emis_ind    )) allocate(s% emis_ind    (n_prio,0:n_styp-1, s%n_chan))
      if (associated(s% tskin_ind   )) allocate(s% tskin_ind   (n_prio,0:n_styp-1, s%n_chan))
      if (associated(s% var         )) allocate(s% var         (s%n_chan      ))
      if (associated(s% oe_str      )) allocate(s% oe_str      (s%n_chan      ))
      if (associated(s% oe_trf      )) allocate(s% oe_trf      (s%n_chan      ))
      if (associated(s% bc% p_tr    )) allocate(s% bc% p_tr    (s%bc%n_tr, 2  ))
      if (associated(s% bc% ib_ch   )) allocate(s% bc% ib_ch   (s%n_chan      ))
      if (associated(s% bc% type    )) allocate(s% bc% type    (s%n_chan      ))
    end if

    if (associated (s% sensor_instr)) call p_bcast (s% sensor_instr(1:s%n_sens),  source, lcom)
    if (associated (s% chan        )) call p_bcast (s% chan        (1:s%n_chan),  source, lcom)
    if (associated (s% chidx       )) call p_bcast (s% chidx       (1:s%n_chan),  source, lcom)
    if (associated (s% iflag       )) call p_bcast (s% iflag       (1:s%n_chan),  source, lcom)
    if (associated (s% band        )) call p_bcast (s% band        (1:s%n_chan),  source, lcom)
    if (associated (s% flag        )) call p_bcast (s% flag        (1:s%n_chan),  source, lcom)
    if (associated (s% wavenum     )) call p_bcast (s% wavenum     (1:s%n_chan),  source, lcom)
    if (associated (s% emis_land   )) call p_bcast (s% emis_land   (1:s%n_chan),  source, lcom)
    if (associated (s% emis_snow   )) call p_bcast (s% emis_snow   (1:s%n_chan),  source, lcom)
    if (associated (s% emis_sice   )) call p_bcast (s% emis_sice   (1:s%n_chan),  source, lcom)
    if (associated (s% emis_pc     )) call p_bcast (s% emis_pc     (1:n_pc,1:s%n_chan),  source, lcom)
    if (associated (s% emis_opt    )) then
      do i = 1, s%n_emis_opt
        call p_bcast (s% emis_opt(i),source, lcom)
      end do
    end if
    if (associated (s% tskin_opt    )) then
      do i = 1, s%n_tskin_opt
        call p_bcast (s% tskin_opt(i),source, lcom)
      end do
    end if
    if (associated (s% emis_ind    )) call p_bcast (s% emis_ind    (1:n_prio,0:n_styp-1,1:s%n_chan),  source, lcom)
    if (associated (s% tskin_ind   )) call p_bcast (s% tskin_ind    (1:n_prio,0:n_styp-1,1:s%n_chan),  source, lcom)
    if (associated (s% var         )) call p_bcast (s% var         (1:s%n_chan),  source, lcom)
!   if (associated (s% oe_str      )) call p_bcast (s% oe_str      (1:s%n_chan),  source, lcom)
    if (associated (s% oe_str)) then
      do i = 1, s%n_chan
        call p_bcast (s% oe_str(i),  source, lcom)
      end do
    end if
    if (associated (s% bc% ib_ch   )) call p_bcast (s% bc% ib_ch   (1:s%n_chan),    source, lcom)
    if (associated (s% bc% p_tr    )) call p_bcast (s% bc% p_tr    (1:s%bc%n_tr,:), source, lcom)
    if (associated (s% bc% type    )) call p_bcast (s% bc% type    (1:s%bc%n_tr  ), source, lcom)

  end subroutine bcast_rad_set


  subroutine bcast_emis_opt (eo, source, comm)
    type (t_emis_opt)     ,intent(inout) :: eo
    integer               ,intent(in)    :: source
    integer   ,optional   ,intent(in)    :: comm

    integer :: lcom
    integer :: count = 0

    lcom = dace% comm ; if (present(comm)) lcom = comm

    if (count == 0) count = size(transfer(eo,(/' '/)))

    if (dace% pe /= source) call destruct (eo)
    call p_bcast_derivedtype (eo, count, source, lcom)

    if (dace% pe /= source) then
      if (associated(eo%chan)) allocate(eo%chan(eo%n_chan))
      if (associated(eo%cdyn)) allocate(eo%cdyn(eo%n_chan))
    end if

    if (associated(eo%chan)) call p_bcast(eo%chan(1:eo%n_chan), source, lcom)
    if (associated(eo%cdyn)) call p_bcast(eo%cdyn(1:eo%n_chan), source, lcom)

  end subroutine bcast_emis_opt


  subroutine bcast_tskin_opt (tso, source, comm)
    type (t_tskin_opt)    ,intent(inout) :: tso
    integer               ,intent(in)    :: source
    integer   ,optional   ,intent(in)    :: comm

    integer :: lcom
    integer :: count = 0

    lcom = dace% comm ; if (present(comm)) lcom = comm

    if (count == 0) count = size(transfer(tso,(/' '/)))

    if (dace% pe /= source) call destruct (tso)
    call p_bcast_derivedtype (tso, count, source, lcom)

    if (dace% pe /= source) then
      if (associated(tso%cdyn)) allocate(tso%cdyn(tso%n_cdyn))
    end if

    if (associated(tso%cdyn)) call p_bcast(tso%cdyn(1:tso%n_cdyn), source, lcom)

  end subroutine bcast_tskin_opt



end module mo_rad
