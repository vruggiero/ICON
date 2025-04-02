!
!+ Definition of generic observation data type T_OBS and level operations
!
MODULE mo_t_obs
!
! Description:
!   This module MO_T_OBS defines the generic observation data type T_OBS
!   and low level operations on this data type. This module will be used
!   by the specific observation modules MO_SYNOP, MO_TEMP,.. . Finally
!   these modules are used by MO_OBS which provides the high level
!   interface to the observation data type and procedures.
!
! Current Maintainer: DWD, Harald Anlauf
!    phone: +49 69 8062 4941
!    fax:   +49 69 8062 3721
!    email: harald.anlauf@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_1         2008/11/05 Andreas Rhodin
!  First operational 3D-Var release
! V1_4         2009/03/26 Andreas Rhodin
!  changes for verification mode
! V1_5         2009/05/25 Andreas Rhodin
!  subroutine release_mem_1: new optional parameter "sort"
! V1_8         2009/12/09 Andreas Rhodin
!  new components: t_spot% sozen             ( solar zenith angle )
!                          mon_file, mon_rec ( position in monitoring file )
!                          prc1d             ( 1dvar processing flag )
! V1_9         2010/04/20 Andreas Rhodin
!  shrink_report: check for components being associated
!  change spot% phase from i1 to i2 to be able to hold GPSRO PCD flags
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Detlef Pingel
!  changes for IASI, SATPP
! V1_19        2012-04-16 Andreas Rhodin
!  option to specify separate input directory for observations
!  new subroutines: set_lluv, shrink_reports, join_obs, empty_source
! V1_20        2012-06-18 Andreas Rhodin
!  allow up to 100 observation input files
! V1_22        2013-02-13 Andreas Rhodin
!  changes for STD, RADAR and COSMO observation operators
! V1_23        2013-03-26 Andreas Rhodin
!  added RADAR operator
! V1_26        2013/06/27 Andreas Rhodin
!  option to fix misbehaviour of check_no_obs
!  (dont set passive observations to passive_rejected)
! V1_27        2013-11-08 Andreas Rhodin
!  implementation of Variational Bias Correction (VarBC)
!  set flag fix_no_obs to true by default (changes passive observation status)
! V1_28        2014/02/26 Andreas Rhodin
!  preparations for bound constrained sink variables
! V1_29        2014/04/02 Andreas Rhodin
!  transfer first guess from monitoring to analysis pass
! V1_31        2014-08-21 Andreas Rhodin
!  consider snow fraction (derived from snow height) in IR emissivity model.
!  modified rttov specific vertical interpolation (nwv_rad=3 in /observations/)
!  change defaults in namelist /observations/ : read_bufr = F, read_NetCDF = F
! V1_35        2014-11-07 Andreas Rhodin
!  revise setting of mdlsfc (model surface flag) for sea-ice and snow
!  new optional parameters to subroutine shrink_report (re-ordering of levels)
! V1_37        2014-12-23 Harald Anlauf
!  Change ijdp, index_x from i2 to int
! V1_42        2015-06-08 Andreas Rhodin
!  implement temporal interpolation for COSMO MEC
! V1_43        2015-08-19 Andreas Rhodin
!  changes for MEC and MISR aircraft observations
! V1_44        2015-09-30 Harald Anlauf
!  preparations for Jason-2
!  move type t_bufr_inv to mo_t_obs to simplify dependencies
! V1_46        2016-02-05 Harald Anlauf
!  Fix DRIBU station id; extend SYNOP/DRIBU code for monitoring of T2m
! V1_47        2016-06-06 Andreas Rhodin
!  add obstype as selection criterium in namelist /rules/
! V1_48        2016-10-06 Robin Faulwetter
!  Implemented level-based thinning
!  handling of passive variables (for MEC)
! V1_49        2016-10-25 Harald Anlauf
!  t_spot: generalize bc_airep -> bc_index
! V1_50        2017-01-09 Andreas Rhodin
!  option for higher order interpolation in apply_L_m; 10 character 'statid'
! V1_51        2017-02-24 Andreas Rhodin
!  change meaning of sink variable mode (1/2)
!  changes to pass observations from COSMO fof-files through the Var-scheme
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Andreas Rhodin MPIfM/DWD  2000-2008
! Gerhard Paul   DWD        2008       modifications for NetCDF interface
! Harald Anlauf  DWD        2008       modifications for SX8
!==============================================================================
#include "tr15581.incf"
  !=============
  ! modules used
  !=============
  use mo_kind,       only: wp,             &! working precision kind parameter
                           sp,             &! single  precision kind parameter
                           i1, i2, i8       ! integer precision kind parameter
  use mo_exception,  only: finish,         &! exit on error
                           message          ! warning message
  use mo_time,       only: t_time,         &! time/date information data type
                           solar            ! calculate solar zenith angle
  use mo_mpi_dace,   only: dace,           &! MPI communication info
                           p_bcast,        &! generic MPI_BCAST  routine
                           p_allgather      ! generic MPI_ALLGATHER routine
  use mo_namelist,   only: position_nml,   &! routine to position nml group
                           nnml,           &! namelist fortran unit number
                           POSITIONED       ! position_nml: OK    return flag
  use mo_dace_string,only: char3            ! conversion: int -> char(len=3)
  use mo_vqc,        only: svqc,           &! default   a priory probability
                           vqc_form         !
  use mo_obs_rules,  only: get_rule,       &! routine to get a rule
                           t_set,          &! result data type
                           iud,            &! undefined integer value
                           rud,            &! undefined real value
                           t_ilev,         &! instrument "level", e.g. channel
                           empty_ilev,     &! empty t_ilev
                           operator(==)     ! t_ilev comparison
  use mo_t_use,      only: t_use            ! status flag data type
  use mo_t_datum,    only: t_datum,        &! observed datum data type
                           rvind,          &! invalid value
                           qbit_conv        ! qbit convention in BUFR/NetCDF
  use mo_sink,       only: t_sink,         &! sink variable meta data
                           set_v_sink       ! set sink variable meta data
  use mo_wigos,      only: t_wsi,          &! WIGOS station identifier data type
                           read_wigos_nml   ! Read WIGOS namelist
  use mo_fdbk_tables,only: init_fdbk_tables,    &! initialise the table (fill in table entries)
                           VN_U,VN_V,VN_W, VN_Z,&! variable numbers
                           VN_T, VN_RH, VN_Q,   &
                           VN_VT, VN_RAWBT,     &
                           VN_BENDANG, VN_REFR, &
                           VN_FF, VN_T2M, VN_PS,&
                           VN_RH2M, VN_HEIGHT,  &
                           VN_U10M, VN_V10M,    &
                           VN_NH, VN_N_L, VN_N, &
                           VN_P, VN_TD2M, VN_DD,&
                           VN_N_M, VN_N_H,      &
                           VN_GUST, VN_RR,      &
                           VN_TMIN, VN_TMAX,    &
                           VN_PTEND,VN_VV,      &
                           VN_RAD_GL,VN_RAD_DF, &
                           VN_RAD_LW, VN_SDEPTH,&
                           VN_TSEA, VN_HLOS,    &
                           VN_VGUST, VN_TURB,   &
                           VN_CEIL, VN_WW,      &
                           VN_GCLG, VN_ICLG,    &
                           VN_ZPD,  VN_MIXR,    &
                           MS_LAND, MS_SEA,     &! model surface characteristics
                           MS_ICE, MS_NO_ICE,   &
                           MS_SNOW, MS_NO_SNOW, &
                           n_ot,                &! number of observation types
              obstype_t => obstype,             &! observation type table
                           OT_RAD,              &! RAD observation type
                           varno                 ! 'varno' table
  use mo_t_table    ,only: name_value       ! find name of table entry
  use mo_run_params ,only: obsinput,       &! default input-data path
                           path_file        ! concatenate path/filename
  use mo_algorithms ,only: index            ! sort
  use mo_instrid,    only: rttov_instr      ! RTTOV number from instrument id
  use mo_dace_string,only: split            ! put words of string into array
#if (_RTTOV_VERSION >= 12)
  use rttov_const,  only : interp_rochon
#endif


  implicit none
  !================
  ! public entities
  !================
  private
  !----------------------
  ! data type definitions
  !----------------------
  public :: t_obs        ! generic observation             data type
  public :: t_spot       ! component of t_obs
  public :: t_datum      ! component of t_obs
  public :: t_mcols      ! component of t_obs
  public :: t_mcol       ! component of t_mcols
  public :: t_head       ! component of t_spot
  public :: t_index      ! component of t_spot
  public :: t_box        ! component of t_spot
  public :: t_lev        ! component of t_spot
  public :: t_vic        ! vertical   interpolation coefficient data type
  public :: t_hic1       ! horizontal interpolation coefficient data type (1d)
  public :: t_hic        ! horizontal interpolation coefficient data type (2d)
  public :: t_coord      ! coordinates
  public :: t_icol       ! component of t_spot
  public :: t_imcol      ! component of t_spot
  public :: t_obsq       ! observed quantity table entry
  public :: add_source   ! add source-file to list
  public :: empty_source ! clear source-file list
  public :: source       ! list of source-files
  public :: t_source     ! derived type of 'source'
  public :: n_source     ! number of source files
  public :: m_source     ! max number of source files
  public :: n_filetype   ! number of files of specific type
  public :: m_filetype   ! number of different file types
  public :: p_bcast      ! broadcast t_source
  public :: t_flgbits    ! satpp bit field type
  public :: mwv          ! max.no.coefficients for vert.interpol.
  public :: t_bufr_inv   ! observation inventory data type
  public :: bufr_inv     ! observation inventory data
  public :: bufr_inv0    ! empty observation inventory template
  public :: empty_spot   ! empty spot
  public :: npmax        ! size of some arrays in t_hic

  !-----------
  ! operations
  !-----------
  public :: construct      ! constructor module procedures
  public :: destruct       ! destructor  module procedures
  public :: release_mem    ! release unused memory in t_obs
  public :: join_obs       ! join observations from multiple boxes
  public :: gather_obs     ! gather observations from processors
  public :: scatter_obs    ! scatter observations to  processors
  public :: reorder_obs    ! order observations according to PEs
  public :: filter_obs     ! filter observations for valid PEs
  public :: mean_x_obs     ! set mean coordinates of observ. in box
  public :: new_spot       ! reserve memory for new spot
  public :: unique_spot    ! ensure unique obs id in parallel environment
  public :: unitvector     ! return unit vector on the sphere
  public :: set_xuv        ! set components of t_icol data type
  public :: set_lluv       ! set lat,lon,du,dv from unit vector
  public :: new_par        ! reserve memory for new parameters
  public :: new_obs        ! reserve memory for new observations
  public :: new_int        ! reserve memory for physical space
  public :: new_sink       ! reserve memory for sink variable
  public :: col2lev        ! set component 'levs' of derived type t_obs
  public :: forget_intp    ! remove any information on the interpolation space
  public :: set_vqc_insitu ! set bounds for variational quality control routine
  public :: set_int_insitu ! set interpolation space
  public :: shrink_report  ! remove passive observations in report
  public :: shrink_reports ! remove passive observations in reports
  public :: set_inbx       ! set component inbx in t_spot
  public :: mdlsfc_frac    ! model surface flag from land/ice/snow fraction
  public :: frac_mdlsfc    ! land/ice/snow fraction from model surface flag
  public :: set_sink       ! define background values for sink variables
  public :: stat_hash      ! hash function for 'statid'

  !----------
  ! constants
  !----------
  public :: invalid                              ! invalid observation value
  public :: SYNOP, TEMP, TOVS, GPSRO, AMV,      &! observation operator ident.
            AIREP, SATEM, GPSGB, COSMO, RADAR,  &!
            SOIL, WLIDAR, SCATT, MWR, OP_X       !
  public :: OBS_TV, OBS_RH, OBS_U, OBS_V, OBS_H,&! observation type identifier
            OBS_HS, OBS_BT, OBS_BDA, OBS_Q,     &!
            OBS_CHI, OBS_DUM, OBS_RFR, OBS_T,   &!
            OBS_CTR, OBS_FF, OBS_X, OBS_HLOS,   &!
            OBS_ZPD, OBS_MIXR                    !
  public :: INT_TV, INT_RH, INT_U, INT_V, INT_H,&! interpolation type identifier
            INT_PSI, INT_CHI                     !
  public :: TSK_INIT, TSK_READ                   ! task identifier
  public :: TSK_SET_CHR, TSK_SETUP_COLS, TSK_SETUP_FULL, TSK_SETUP_FUL0
  public :: TSK_R, TSK_Y, TSK_YH, TSK_K, TSK_XK, TSK_H, TSK_SETUP_OP
  public :: TSK_SHRINK
  public :: BUF_P, BUF_GP, BUF_T, BUF_TD, BUF_FF, BUF_DD ! BUFR codes
  public :: ITY_ICOL, ITY_ICOLS, ITY_MCOLS
  public :: CHR_ID,CHR_LIN,CHR_NONL,CHR_INV,CHR_EXP,CHR_REQ,CHR_EVAL,CHR_COM
  public :: FT_BUFR, FT_NETCDF, FT_ODB, FT_SATPP, FT_ROPIC,             &
            FT_UNKNOWN, FT_MISSING, FT_FEEDBACK, FT_SPECIFIC, FT_EMPTY, &
            FT_RAW_FEEDBACK, FT_CDFIN
  public :: mstep     ! number of predefined feedback file types (cof, mon, ekf, ver)
  public :: pref_step ! prefix for feedback files
  public :: comm_step ! comment for feedback files
  !--------------------------------------------------------
  ! byte sizes of observation data type and its components,
  ! routine to set byte sizes
  !--------------------------------------------------------
  public :: obs_bytes, spot_bytes, body_bytes, int_bytes, wp_bytes
  public :: set_byte_size
  !--------------------------------------------------------
  ! provide names (character strings) for integer constants
  !--------------------------------------------------------
  public :: rct_name ! name of report (code) type  ('TEMP','SYNOP',..)
  public :: rct_key  ! key  of report (code) type  ('TEMP','SYNOP',..)
  public :: oq_name  ! name of observed quantity ('t','rh',..)
  public :: oq_key   ! key  of observed quantity ('t','rh',..)
  public :: oq_varno ! varno of observed quantity (written to feedback file)
  public :: varno_oq ! varno to observed quantity (written to feedback file)
  public :: ft_name  ! name of file type ('BUFR', 'NETCDF' ..)
  public :: tsk_name ! name of task
  public :: obsq     ! table of observed quantities (OBS_H,_TV,_U,..)
  !------------------------
  ! namelist /OBSERVATIONS/
  !------------------------
  public :: read_bufr    ! read BUFR files
  public :: read_NetCDF  ! read NetCDF files
  public :: read_cdfin   ! read COSMO observation files
  public :: dace_op      ! use dace obs. operators (MEC)
  public :: n_dace_op    ! length of non-empty entries in dace_op
  public :: obs_local    ! keep observations on gridpoint PE (for COSMO MEC)
  public :: obs_path     ! path to read BUFR files
  public :: bufr_files   ! BUFR   files
  public :: obs_files    ! NetCDF files
  public :: fdb_files    ! Feedback (input) files
  public :: read_obs_nml ! subroutine to read namelist OBSERVATIONS
  public :: bufr_verb    ! verbose flag for BUFR decoding routines
  public :: bufr_pause   ! wait (for carriage return) after each record read
  public :: netcdf_verb  ! verbosity level of NetCDF decoding routines
  public :: derive_dbkz  ! derive suitable DBKZ if not present
  public :: fix_no_obs   ! set.false. to revert no_obs-check for compatibility
  public :: int_vh       ! interpolation: 1. vertical --> 2. horizontal
  public :: int_nn       ! horizontal interpolation: nearest neighbour only
  public :: int_rad_hum  ! humidity variable for vertical interp. for radiances
  public :: int_rad_hum_pblend ! for int_rad_hum==2: lower/upper height for blending
                               ! RH-interpolated and Q-interpolated profile
  public :: vint_lin_t   ! linear vertical interpolation for temperature
  public :: vint_lin_z   ! linear vertical interpolation for geopotential
  public :: vint_lin_uv  ! linear vertical interpolation for wind
  public :: vint_lin_tov ! linear vertical interpolation for RTTOV operator
  public :: pcc_conv     ! per cent confidence usage: 0=none 1=pegasus 2=sky
  public :: corme_conv   ! correction message handling
  public :: nwv_rad      ! radiances vert.intpol.flag -1:as usual, 3:Rochon
  public :: vint_rttov   ! if nwv_rad==4 (RTTOV interp.) this is the interpolation mode
  public :: vint_rttov_b ! interpolation mode for rttov-interp. of B
  public :: ndv_rad      ! vertical differentiation flag
  public :: z2p_amv      ! AMV: derive p. from height 0:no; 1:US std atm; 2:bg
  public :: z2p_hlos     ! WLIDAR: derive p. from height 0:no; 1:US std atm; 2:bg
  public :: ptop_lapse   ! pressure bounds for derivation ..
  public :: pbot_lapse   ! of lapse rate (for SYNOP t2m extrapolation)
  public :: monitor_ff   ! monitor wind speed (in addition to wind components)
  public :: monitor_dd   ! monitor wind direction (in addition to u, v)
  public :: fdbk_split   ! split feedback files
  !--------------------------------------
  ! Instrument "level" information
  !--------------------------------------
  public :: t_ilev
  public :: empty_ilev
  public :: operator(==)
  public :: match_ilev
  !----------------------------
  ! Context of process_obs call
  !----------------------------
  public :: po_context   ! context of call to process_obs, e.g. fg, ana, outer loop number
  public :: po_ilns      ! context of call to process_obs: call number in linesearch
  public :: po_lmon      ! context of call to process_obs: monitoring call?
  public :: POC_PRE      ! process_obs context: pre-calculation (before first guess)
  public :: POC_FG       ! process_obs context: first guess
  public :: POC_ANA      ! process_obs context: analysis
  public :: POC_TOP      ! process_obs context: test obs. operator
  public :: POC_NN       ! process_obs context: not known
  !-----------------------
  ! Debug selected spot(s)
  !-----------------------
  public :: spt_hd_debug
  public :: spt_debug
  public :: usd
  public :: debug_spot
  public :: ldeb
  public :: dpref
  public :: ldeb_all
  public :: ldeb_spot
!------------------------------------------------------------------------------
  !==========
  ! constants
  !==========
! real(wp) ,parameter :: invalid = -999._wp
  real(wp) ,parameter :: invalid = rvind
  integer  ,parameter :: npmax   = 16 ! max. number of neighbour gridpoints
  integer  ,parameter :: mstep   =  4 ! # of feedback file types (cof, mon, ekf, ver)
  character(len=3)  ,parameter :: pref_step (mstep) = &
       (/'mon',                    'cof',                    'ekf',                    'ver'                    /)
  character(len=22) ,parameter :: comm_step (mstep) = &
       (/'monitoring            ', 'analysis              ', 'Ensemble Kalman Filter', 'verification          ' /)
  !-------------------------------------------------------
  ! byte sizes of observation data type and its components
  !-------------------------------------------------------
  integer :: obs_bytes  = 0
  integer :: spot_bytes = 0
  integer :: body_bytes = 0
  integer :: int_bytes  = 0
  integer :: wp_bytes   = 0
  !------------------------------------------------------------------------------
  ! Module Types
  ! These codes are used to identify the modules (program code) which is used to
  ! handle specific observation types. Currently certain observation types may be
  ! handled by different alternative modules dependent on the application
  ! (regional data assimilation, COSMO / ICON-LAM or global data assimilation,
  ! ICON). The choice can be made at runtime. These codes are only used intenally
  ! in the data assimilation code.
  !------------------------------------------------------------------------------
  integer ,parameter :: SYNOP  =     1 ! SYNOP, BUOY, ...
  integer ,parameter :: TEMP   =     2 ! TEMP, PILOT, Wind-profiler, ..
  integer ,parameter :: TOVS   =     4 ! Radiances from Sattelites
  integer ,parameter :: GPSRO  =     8 ! GNSS radio occultations
  integer ,parameter :: AMV    =    16 ! Satellite winds
  integer ,parameter :: PSTEMP =    32 ! pseudo-temp (diagnostics only)
  integer ,parameter :: AIREP  =    64 ! Aircraft observations
  integer ,parameter :: SATEM  =   128 ! Satellite derived profiles
  integer ,parameter :: GPSGB  =   256 ! Ground based GPS observations
  integer ,parameter :: COSMO  =   512 ! COSMO conventional data
  integer ,parameter :: RADAR  =  1024 ! Volume radar operator
  integer ,parameter :: SOIL   =  2048 ! ASCAT soil moisture
  integer ,parameter :: WLIDAR =  4096 ! Wind Lidar
  integer ,parameter :: SCATT  =  8192 ! Scatterometer observations
  integer ,parameter :: MWR    = 16384 ! Microwave radiometer (ground-based)
  integer ,parameter :: OP_X   = 32768 ! next (unused) module type
  !--------------------
  ! observed quantities
  !--------------------
  integer ,parameter :: OBS_H   =     1 ! geopotential
  integer ,parameter :: OBS_T   =     2 ! temperature
  integer ,parameter :: OBS_RH  =     4 ! relative humidity
  integer ,parameter :: OBS_U   =     8 ! wind component u
  integer ,parameter :: OBS_V   =    16 ! wind component v
  integer ,parameter :: OBS_HS  =    32 ! surface geopotential
  integer ,parameter :: OBS_BT  =    64 ! brightness temperature
  integer ,parameter :: OBS_BDA =   128 ! bending angle
  integer ,parameter :: OBS_Q   =   256 ! specific humidity
  integer ,parameter :: OBS_DUM =   512 ! dummy (sink) variable
  integer ,parameter :: OBS_CHI =  1024 ! velocity potential
  integer ,parameter :: OBS_TV  =  2048 ! virtual temperature
  integer ,parameter :: OBS_RFR =  4096 ! refractivity
  integer ,parameter :: OBS_CTR =  8192 ! additional control variable
  integer ,parameter :: OBS_FF  = 16384 ! wind speed
  integer ,parameter :: OBS_DRH = 32768 ! RH dummy (currently not used)
  integer ,parameter :: OBS_HLOS= 65536 ! hor. line-of-sight wind
  integer ,parameter :: OBS_ZPD =131072 ! slant/zenith path delay
  integer ,parameter :: OBS_MIXR=262144 ! mixing ratio
  integer ,parameter :: OBS_X   =524288 ! next (unused) value
  !------------------------
  ! interpolated quantities
  !------------------------
  integer ,parameter :: INT_H   =    1 ! geopotential
  integer ,parameter :: INT_RH  =    2 ! relative humidity
  integer ,parameter :: INT_CHI =    3 ! velocity potential
  integer ,parameter :: INT_TV  =    4 ! virtual temperature
  integer ,parameter :: INT_PSI =    5 ! stream   function
  integer ,parameter :: INT_U   =    6 ! wind component u
  integer ,parameter :: INT_V   =    7 ! wind component v
  integer ,parameter :: INT_X   =    8 ! next (unused) observation type

  type t_obsq
    character(len=8)  :: name ! mnemonic
    integer           :: key  ! key
    character(len=32) :: desc ! description
  end type t_obsq

  !----------------------
  ! observation data type
  !----------------------
  type (t_obsq) ,parameter :: obsq(18) =                        &
    (/ t_obsq ('z'     ,OBS_H   ,'geopotential'),               &
       t_obsq ('t'     ,OBS_T   ,'temperature'),                &
       t_obsq ('rh'    ,OBS_RH  ,'relative humidity'),          &
       t_obsq ('u'     ,OBS_U   ,'wind u component'),           &
       t_obsq ('v'     ,OBS_V   ,'wind v component'),           &
       t_obsq ('zs'    ,OBS_HS  ,'geopotential at the surface'),&
       t_obsq ('rawbt' ,OBS_BT  ,'brightness temperature'),     &
       t_obsq ('bangl' ,OBS_BDA ,'bending angle'),              &
       t_obsq ('q'     ,OBS_Q   ,'specific humidity'),          &
       t_obsq ('sink'  ,OBS_DUM ,'dummy sink variable'),        &
       t_obsq ('chi'   ,OBS_CHI ,'velocity potential'),         &
       t_obsq ('vt'    ,OBS_TV  ,'virtual temperature'),        &
       t_obsq ('rfr'   ,OBS_RFR ,'refractivity'),               &
       t_obsq ('cntrl' ,OBS_CTR ,'additional control variable'),&
       t_obsq ('ff'    ,OBS_FF  ,'wind speed')                 ,&
       t_obsq ('hlos'  ,OBS_HLOS,'hor. line-of-sight wind')    ,&
       t_obsq ('zpd'   ,OBS_ZPD ,'zenith path delay')          ,&
       t_obsq ('mixr'  ,OBS_MIXR,'mixing ratio')               /)
  !----------------------------------------------------------------------
  ! 'interpolation type' of the observation operator:
  ! ITY_ICOL:  The observation operator depends on exactly one horizontal
  !            location.
  ! ITY_ICOLS: The observation operator depends on several horizontal
  !            locations.  (Not yet implemented)
  ! ITY_MCOLS: The observation operator depends on several model columns.
  !----------------------------------------------------------------------
  integer ,parameter :: ITY_ICOL  = 1 ! interpolation at one column
  integer ,parameter :: ITY_ICOLS = 2 ! interpolation at multiple columns
  integer ,parameter :: ITY_MCOLS = 3 ! dependence on multiple model columns
  !-----------
  ! BUFR codes
  !-----------
  integer ,parameter :: BUF_P      =  7004 ! Pressure              code
  integer ,parameter :: BUF_GP     = 10003 ! Geopotential          code
  integer ,parameter :: BUF_T      = 12001 ! Temperature           code
  integer ,parameter :: BUF_TD     = 12003 ! Dew Point Temperature code
  integer ,parameter :: BUF_FF     = 11002 ! Wind speed            code
  integer ,parameter :: BUF_DD     = 11001 ! Wind direction        code
  !--------------------------------------------------------------
  ! flags for tasks to be performed by 'setup_obs', 'process_obs'
  !--------------------------------------------------------------
  integer ,parameter :: TSK_INIT       =     1 ! initialize modules
  integer ,parameter :: TSK_READ       =     2 ! read observations
  integer ,parameter :: TSK_SET_CHR    =     4 ! set observ. characteristics
  integer ,parameter :: TSK_SHRINK     =     8 ! release unused obs. in report
  integer ,parameter :: TSK_SETUP_COLS =    16 ! setup columns
  integer ,parameter :: TSK_SETUP_FUL0 =    32 ! setup interpolation space
  integer ,parameter :: TSK_SETUP_FULL =    64 ! setup description of PSAS-space
  integer ,parameter :: TSK_R          =   128 ! setup observational error
  integer ,parameter :: TSK_Y          =   256 ! run forward operator
  integer ,parameter :: TSK_YH         =   512 ! run linear or forward operator
  integer ,parameter :: TSK_H          =  1024 ! run tangent linear operator
  integer ,parameter :: TSK_K          =  2048 ! evaluate linear operator
  integer ,parameter :: TSK_XK         =  4096 ! exchange linear operator
  integer ,parameter :: TSK_SETUP_OP   =  8192 ! set up obs. operator
! integer ,parameter :: TSK_FEEDBACK   =  8192 ! write feedback
!                                   2147483647 ! huge(1)
  !-----------------------------------------
  ! characteristics of observation operators
  !-----------------------------------------
  integer ,parameter :: CHR_ID         =   1 ! identity operation
  integer ,parameter :: CHR_LIN        =   2 ! linear
  integer ,parameter :: CHR_NONL       =   4 ! nonlinear
  integer ,parameter :: CHR_INV        =   8 ! invertible
  integer ,parameter :: CHR_EXP        =  16 ! expensive
  integer ,parameter :: CHR_REQ        =  32 ! required on this PE
  integer ,parameter :: CHR_EVAL       =  64 ! evaluated on this PE
  integer ,parameter :: CHR_COM        = 128 ! R, H are communicated over PE's

  !-----------------
  ! input file types
  !-----------------
  integer ,parameter :: FT_BUFR         =  1 ! BUFR   file
  integer ,parameter :: FT_NETCDF       =  2 ! NETCDF file (converted from BUFR)
  integer ,parameter :: FT_ODB          =  3 ! read from ODB (not supported)
  integer ,parameter :: FT_SATPP        =  4 ! file from SATPP preprocessing
  integer ,parameter :: FT_ROPIC        =  5 ! ROPIC (radio occultation) file
  integer ,parameter :: FT_FEEDBACK     =  6 ! NetCDF feedback file (fof etc.)
  integer ,parameter :: FT_SPECIFIC     =  7 ! observation type specific
  integer ,parameter :: FT_UNKNOWN      =  8 ! unknown file format
  integer ,parameter :: FT_MISSING      =  9 ! missing file
  integer ,parameter :: FT_EMPTY        = 10 ! empty file
  integer ,parameter :: FT_RAW_FEEDBACK = 11 ! NetCDF raw feedback file (rof etc.)
  integer ,parameter :: FT_CDFIN        = 12 ! NetCDF file (for use with COSMO operators)
  integer ,parameter :: m_filetype      = 12 ! number of different file types

  !----------------------------
  ! context of process_obs call
  !----------------------------
  integer ,parameter :: POC_FG      =  0 ! run obs.operators on first guess
  integer ,parameter :: POC_ANA     = -1 ! run obs.operators on analysis
  integer ,parameter :: POC_TOP     = -2 ! test obs.operator
  integer ,parameter :: POC_PRE     = -3 ! run obs.operators preliminarily on first guess
  integer ,parameter :: POC_NN      =-99 ! not known
                                         ! values > 0 are the outer loop number

!==============================================================================
  !=================
  ! type definitions
  !=================

  !-----------------------
  ! index to another array
  !-----------------------
  type t_index
#ifdef USE_SEQUENCE
    SEQUENCE
#endif
    integer        :: i       = -1      ! index
    integer        :: n       =  0      ! size
  end type t_index

  !------------------------------------
  ! vertical interpolation coefficients
  !------------------------------------
  integer ,parameter :: mwv = 9 ! number of points for vert. interpolation
  type t_vic
#ifdef USE_SEQUENCE
    SEQUENCE
#endif
    real(wp)      :: wh (mwv) ! weights for vertical interpol. (height)
    real(wp)      :: wt (mwv) ! weights for vertical interpol. (temperature)
    type(t_index) :: g        ! index   for vertical interpol. grid
    real(wp)      :: ezn      ! d e_h / d logp  / e_t
    real(wp)      :: exzn     ! e_h             / e_t
    type(t_index) :: i        ! index to 'interpolation space'
  end type t_vic

  !--------------------------------------
  ! horizontal interpolation coefficients
  !--------------------------------------
  integer ,parameter :: mwh = 2         ! no.points for horizontal interpol.
  type t_hic1
#ifdef USE_SEQUENCE
    SEQUENCE
#endif
    real(wp)         :: w (mwh) = 0._wp ! weights for horizontal interpol.
    integer          :: ix      = 0     ! index   for horizontal interpol.
  end type t_hic1

  !----------------------------------------------------------------------------
  !------------------------------
  ! generic observation data type
  !------------------------------
  !
  ! The derived type 't_obs' holds all observation specific
  ! information. In general an array with elements of type 't_obs' is
  ! used, with most of the components of each array element being
  ! allocated only on one processor.
  !
  ! The component 'spot' holds the header data for the 'n_spot' reports
  ! stored in an array element.
  !
  ! The total number of observations stored is 'n_obs'. Arrays of this
  ! size are used for the observed values ('obs') the type of the observed
  ! quantity ('t_obs') as well as forecast ('fc'),analysis ('ana') and
  ! quality control flags ('qcf'). Arrays may be allocated as required.
  !
  ! The number of model parameters interpolated to the locations of
  ! observations (or the locations required by the observation operators)
  ! is 'n_int'. The type of the parameters is given by 't_int', their
  ! height by 'lev'.
  !
  ! 'n_par' is the size of the array 't_par', which is used to store
  ! parameters specific to individual report types. The content of this
  ! array is interpreted by the respective observation operators.
  !

  !------------------------------------------------------------
  ! model column descriptor
  !
  ! describes a model column required to evaluate observations.
  !------------------------------------------------------------
  type t_mcol
#ifdef USE_SEQUENCE
    SEQUENCE
#endif
    integer     :: ijdtp(5) ! model grid indices and processor: i1,i2,id,t,pe
    integer     :: icol     !
    integer(i8) :: iatm     ! ored bit field for atmospheric parameter ids
    integer     :: itrac    ! ored bit field for tracers
  end type t_mcol
  !--------------------------------
  ! set of model column descriptors
  !--------------------------------
  type t_mcols
#ifdef USE_SEQUENCE
    SEQUENCE
#endif
    integer      ,pointer  :: idx (:,:,:,:)=>NULL() ! flag field for usage
    type(t_mcol) ,pointer  :: c   (:)      =>NULL() ! columns used(i1,i2,id)
    integer                :: n                     ! number of columns used
    integer                :: pe                    ! observation PE
  end type t_mcols
  !
  ! generic observation data type
  !
  type t_obs
#ifdef USE_SEQUENCE
    SEQUENCE
#endif
    !---------------------------------------------------------------------
    ! characteristics of this (type t_obs) array element (observation box)
    !---------------------------------------------------------------------
    integer                :: ibox      =  0      ! box index
    integer                :: pe        = -1      ! preconditioning PE
    integer                :: n_spot    =  0      ! number of observations
    integer                :: n_spt     =  0      ! number of spots
    real(wp)               :: bcc   (3) =  0._wp  ! box center coordinates
    logical                :: bcast     = .true.  ! broadcast results
    !--------------------------------------
    ! characteristics of each report or FoV
    !--------------------------------------
    type (t_spot) _POINTER :: spot  (:) => NULL() ! header information
    integer                :: n_par     =  0      ! number of parameters
    integer       _POINTER :: par   (:) => NULL() ! parameter types
    !-------------------------------------------
    ! characteristics of each single observation
    !-------------------------------------------
    integer                :: n_obs     =  0      ! size of observational space
    type(t_datum) _POINTER :: body  (:) => NULL() ! datum meta information
    integer       _POINTER :: varno (:) => NULL() ! observation type
    real(wp)      _POINTER :: bger  (:) => NULL() ! background error
    real(wp)      _POINTER :: s_vqc (:) => NULL() ! VQC stdev. for baddata
    integer       _POINTER :: f_vqc (:) => NULL() ! VQC formulation
    real(wp)      _POINTER :: olev  (:) => NULL() ! levels in observation space

    integer                :: n_lev     =  0      ! number of levels(int.space)
    type(t_lev)   _POINTER :: levs(:)   => NULL() ! level information
    !-----------------------------------------------
    ! characteristics of the interpolated quantities
    !-----------------------------------------------
    integer                :: n_int     =  0      ! size of interpolated space
    integer       _POINTER :: t_int (:) => NULL() ! types in interpolated space
    real(wp)      _POINTER :: lev   (:) => NULL() ! levels
    real(wp)      _POINTER :: bgeri (:) => NULL() ! background error
    real(wp)      _POINTER :: bgi   (:) => NULL() ! background
    !---------------------------
    ! (nonlinear) sink variables
    !---------------------------
    integer                :: n_sink    =  0      ! number of sink variables
    type(t_sink)  _POINTER :: sink  (:) => NULL() ! sink variables
    !----------------------------------
    ! VBC (Variational bias correction)
    !----------------------------------
    integer                :: n_bcp     =  0      ! size of 'bcpred'
    real(wp)      _POINTER :: bcpred(:) => NULL() ! predictor values
    !-----------------------------------------
    ! model columns required for interpolation
    !-----------------------------------------
    type(t_mcols)          :: mc                  ! neighbour model columns
    integer                :: n_time    = 1       ! number of time slices
  end type t_obs
  !
  ! The header (meta) information applicable to each report type is stored
  ! in the derived type 't_spot'.
  !
  type t_box
#ifdef USE_SEQUENCE
    SEQUENCE
#endif
    integer(i2)    :: pe      = -1      ! processor index
    integer(i2)    :: box     = -1      ! box index
    integer        :: spot    = -1      ! spot index
  end type t_box
  !
  ! coordinates
  !
  type t_coord
#ifdef USE_SEQUENCE
    SEQUENCE
#endif
    real(wp) :: dlon        = invalid ! longitude [degree east]
    real(wp) :: dlat        = invalid ! latitude  [degree north]
    real(wp) :: x  (3)      = invalid ! unit vector on the sphere
    real(wp) :: du (3)      = invalid ! unit vector to the east
    real(wp) :: dv (3)      = invalid ! unit vector to the north
  end type t_coord
  !
  ! horizontal interpolation coefficients
  !
  type t_hic
#ifdef USE_SEQUENCE
    SEQUENCE
#endif
    integer       :: imc (npmax, 2) = 0       ! indices to model columns
    real(wp)      :: w   (npmax)    = 0._wp   ! weights for model columns
    real(wp)      :: w1  (npmax)    = 0._wp   ! gradient ..
    real(wp)      :: w2  (npmax)    = 0._wp   ! .. or sin/cosine
    type(t_index) :: l                        ! index to levels
    integer       :: ijdp(4)        = -1      ! nearest gridpoint: i,j,k,pe
    integer       :: iw12           = 0       ! w1,w2 usage: 0 = none
  end type t_hic                              !              1 = gradient (1/m)
  !                                           ! not implem.: 2 = sin/cosine
  ! interpolation point (column) descriptor
  !
  type t_icol
#ifdef USE_SEQUENCE
    SEQUENCE
#endif
    integer       :: nlev        = 0       ! number of levels
    type(t_coord) :: c                     ! coordinates
    ! integer       :: np          = 0       ! number of parameters
    ! integer(i8)   :: ipar        = 0_i8    ! parameter bit flags
    integer       :: nt          = 0       ! number of tracers
    integer       :: itrac       = 0       ! tracer bit flag
    type(t_hic)   :: h                     ! horizontal interpolation coefs.
  end type t_icol
  !
  ! satpp bit field type
  !
  type t_flgbits
    integer          :: ntest   = 0        ! number of tests/bits
    integer, pointer :: test (:)=>NULL()   ! array of test IDs
    integer, pointer :: btfld(:)=>NULL()   ! bit field (instruments)
    integer, pointer :: instr(:)=>NULL()   ! instrument IDs (RTTOV)
  end type t_flgbits
  !
  ! information available without decoding the body of a (BUFR) report
  !
  type t_head
#ifdef USE_SEQUENCE
    SEQUENCE
#endif
    integer      :: obstype     =  0      ! observation type as used in OI/CMA
    integer      :: codetype    =  0      ! code type as used in CMA
    integer      :: buf_type    = -1      ! BUFR message type
    integer      :: buf_subtype = -1      ! BUFR message subtype
    integer      :: modtype     =  0      ! observation (code) type
    integer      :: dbkz        = -1      ! data bank 'Kennzahl'
    integer      :: idbk        =  1      ! index to dbkz entry
    type(t_time) :: time                  ! time
    type(t_time) :: db_time               ! data bank time
    integer      :: source      =  0      ! source-file of Report
    integer      :: record      =  0      ! record in source-file
    integer      :: subset      =  0      ! subset in BUFR-file
    integer      :: id          =  0      ! observation ID from source/subset
    !
    ! position in monitoring file
    !
    integer      :: mon_file    = -1      ! monitoring file number
    integer      :: mon_rec     = -1      ! monitoring file record
    !
    ! new for BUFR2NetCDF / new feedback file format
    !
    integer(i2)  :: center      = -1      ! originating center
    integer(i2)  :: subcenter   = -1      !         sub-center
    integer(i2)  :: satid       = -1      ! satellite id
    integer(i2)  :: grid_id     =  0      ! RTTOV-ID of target grid
  end type t_head

  type t_spot
#ifdef USE_SEQUENCE
    SEQUENCE
#endif
    !
    ! identification
    !
    type(t_head)     :: hd
    integer          :: id         = 0       ! unique identification (internal)
    integer          :: is         = 0       ! spot index
    character(len=10):: statid     = ''      ! station name (ships call sign)
    integer(i8)      :: stat_hash  = 0       ! "hash value" of station name
    type(t_wsi)      :: wsi                  ! WIGOS station identifier
    integer          :: ident      = 0       ! numeric station id
    integer          :: sttyp      = -1      ! station instrument, sat. sensor, radiosonde typ
    integer          :: stret      = -1      ! station retrieval type; rs solar and IR radiation correction
    integer          :: stlsf      =  0      ! station land/sea flag, AMV:land/sea qualifier
    integer          :: soiltype   = -1      ! soil type of soil retrieval observations
    integer          :: stclf      = -1      ! station cloud flag
    integer          :: meas_type  = -1      ! type of measuring equipment used; (0 02 003 WMO Table B)
    integer          :: tracking   = -1      ! Tracking technique used;          (0 02 014 WMO Table B)
                                             ! GPSRO: GNSS satellite classif.    (0 02 020 WMO Table B)
    real(wp)         :: stzen      = invalid ! station/satellite zenith  angle
    real(wp)         :: stazi      = invalid ! station/satellite azimuth angle
    real(wp)         :: sozen      = invalid ! solar zenith angle
    real(wp)         :: soazi      = invalid ! solar azimuth angle
    real(wp)         :: z          = -999._wp! station height        [m]
    real(wp)         :: ps         = invalid ! pressure at station height [Pa]
!   real(wp)         :: zs         = invalid ! surface height        [m] (10001)
!   real(wp)         :: zg         = invalid ! geopotential height   [gpm?]
    integer(i2)      :: phase      = -1      ! aircraft phase, field of view, GPSRO PCD flags
    integer(i1)      :: pcc        = -1      ! per cent confidence [0..100]
    integer(i1)      :: ncot       = -1      ! nominal confidence threshold [0..100]
    integer          :: sender_id  = -1      ! Transmitter ID      (0 01 150 WMO Table B)
    type(t_time)     :: actual_time          ! observation time
    integer          :: spec_rfl   =  0      ! observation type specific report flags
    integer          :: center_id  = -1      ! processing center ID, e.g. GPSGB
    !
    ! status flags
    !
    type(t_use)     :: use                  ! 'check' and 'status' flags
    character(len=8):: comment    = ''      ! comment to 'check' and 'status'
    integer         :: corme      = 0       ! corrected message
    !
    ! technical parameters
    !
    integer         :: int_type   = 0       ! interpolation type
    integer         :: char       = 0       ! characteristics
    integer         :: pe_eval    = -1      ! PE the operator is evaluated on
    integer         :: n_dest     = 0       ! number of PEs H and R are send to
    integer,pointer :: pe_dest(:) =>NULL()  ! PEs H and R are send to
    type (t_box)    :: srbx                 ! send/receive box actually used
    type (t_box)    :: inbx                 ! input             in this pe/box
    type (t_box)    :: fgbx                 ! first guess check in this pe/box
    type (t_box)    :: pre1                 ! preconditioning   in this pe/box
    real(wp)        :: cost       = invalid ! cost to evaluate operator
    !
    ! background model values at location of the observation
    !
    real(wp)        :: ps_bg      = invalid ! background surface  pressure
    real(wp)        :: gp_bg      = invalid ! background surface  geopotential
    real(wp)        :: sl_bg      = invalid ! background sea-land mask [0=sea]
    real(wp)        :: fi_bg      = invalid ! background sea sea ice fraction
    real(wp)        :: tl_bg      = invalid ! background temp. at lowest level
    real(wp)        :: ts_bg      = invalid ! background temp. at surface
    real(wp)        :: pz_bg      = invalid ! background factor ps(hPa)/z(gpm)
    real(wp)        :: hs_bg      = invalid ! background snow height (m)
    real(wp)        :: dtdz_bg    = invalid ! background dt/dz (850..700h Pa)
    real(wp)        :: u10bg      = invalid ! background 10m wind
    real(wp)        :: v10bg      = invalid ! background 10m wind
!   real(wp)        :: t2mbg      = invalid ! background  2m temperature
    real(wp)        :: z0_bg      = invalid ! background roughless length (m)
    real(wp)        :: ssd_bg     = invalid ! background SSO standard dev.(m)
    real(wp)        :: tw_bg      = invalid ! background sea/water temperature
    integer         :: mdlsfc     = 0       ! model surface flag
    !
    ! horizontal location
    !
    integer         :: n_spt      = 1       ! number of spots
    type(t_icol)    :: col                  ! column descriptor
    type(t_hic1)    :: hic                  ! horizontal interp.coef (bg error)
    integer         :: i_time     = 1       ! time slice
    real(wp)        :: w_time     = 0._wp   ! time slice fraction
    real(wp)        :: emod       = 0._wp   ! bg error modification factor
    !
    ! model columns for operator type ITY_MCOLS
    !
    integer               :: mke      = 0      ! number of model levels
    type(t_imcol),pointer :: imcol(:) =>NULL() ! model column indices
    !
    ! References to other components of t_obs
    !
    integer         :: nr        = 0     ! size of R matrix block
    type (t_index)  :: o                 ! observation space info
    type (t_index)  :: i                 ! interpolated space info
    type (t_index)  :: l                 ! interpolated space, sorted by level
    type (t_index)  :: p                 ! parameter info
    type (t_index)  :: s                 ! observation type specific table
    type (t_index)  :: d                 ! dummy sink variable
    !----------------------------------------------------------
    ! VBC specific entries
    ! bcp  pointer to predictors   (in t_obs, fof specific)
    ! ibcb pointer to coefficients (in t_vqc, channel specific)
    !----------------------------------------------------------
    type (t_index)  :: bcp               ! bias correction predictors
    integer         :: ibcb              ! bias correction coefficient box
    integer         :: ibcfov            ! bias correction FOV index
    integer         :: bc_index          ! conv.bias correction entry index
    !
    ! Observation-specific forward model definition
    !
    integer         :: fwmodel   = 0     ! Forward model id
    real(sp)        :: params(2) = 0._sp ! Model-specific parameter set
    !
    ! COMET specific component
    !
!c  real(wp),pointer:: p_wf(:) =>NULL()  ! Jakobian pressure levels
  end type t_spot
  !
  ! The components 'o,i,c,p,l' of type 't_index' provide the size
  ! (component 'n') and position (offset 'i') of the information related
  ! to this report stored in the arrays 't_obs,obs,fc,ana,qcf' ('o'),
  ! 't_int' ('i'), 't_ctr' ('c') and 'par' ('p').
  !
  type t_lev
#ifdef USE_SEQUENCE
    SEQUENCE
#endif
    integer       :: its  ! variable types (iored) at this level
    real(wp)      :: z    ! level (ln p)
    type(t_vic)   :: vi   ! vertical interpolation coefficients
    type(t_index) :: i    ! index to interpolation space
    type(t_index) :: o    ! index to observation   space
    integer       :: is   ! spot index
  end type t_lev
  !
  ! temporary
  !
  type t_imcol
#ifdef USE_SEQUENCE
    SEQUENCE
#endif
    integer       :: pe    = -1    !
    integer       :: imc(2)=  0    ! index to model column
    integer(i8)   :: iatm  =  0_i8 ! ored bit field: atmospheric parameter ids
    integer       :: natm  =  0    ! number of atmospheric parameters
    integer       :: nsl   =  0    ! number of single level parameters
    integer       :: itrac =  0    ! ored bit field for tracers
    type(t_coord) :: c             ! coordinates: lat,lon, x,du,dv
    type(t_hic1)  :: hi            ! horizontal interp.coef (bg error)
    real(wp)      :: emod  = 0._wp ! bg error modification factor
  end type t_imcol

  !----------------------------------------------------------------------------
  !----------------------------
  ! list of Report source files
  !----------------------------
  type t_source
    character(len=128) :: file      = ''    ! filename
    character(len=128) :: path      = ''    ! path
    integer            :: obstype   = 0     ! use specific input-routine;else 0
    integer            :: codetype  = -1    ! CMA codetype
    integer            :: filetype  = 0     ! 1:BUFR, 2:NetCDF, 3:ODB, ...
    logical            :: used      =.false.! used on this PE
    logical            :: complete  =.false.! read complete file at once
    integer            :: pe        = -1    ! processor (for complete==T)
    integer            :: reports   = 0     ! number of reports
    integer            :: subsets   = 0     ! number of subsets
    integer            :: entries   = 0     ! total number of entries
    integer            :: nobs      = 0     ! total number of observations
    integer            :: cost      = 0     ! cost per entry (read & decode)
    integer            :: format    = 0     ! file format:
                                            ! (0:old, 2:new)
!   character(10)      :: model        = '' ! model creating feedback file
    real(sp)           :: resolution(2)= 0. ! grid resolution in feedback file
    integer            :: domain(2)    = -1 ! hor.domain size parameters nx,ny
  end type t_source
  integer ,parameter  :: m_source                = 150 ! max number allowed
  type(t_source),save :: source (m_source)
  type(t_source),save :: source_0
  integer             :: n_source                = 0 ! actual number
  integer             :: n_filetype (m_filetype) = 0 ! numver of files of type
  !----------------------------------------------------------------------------

  !-----------------------------------------------------------
  ! observation input file inventory data type definition
  ! holds information on distribution of observations over PEs
  ! seperate list items are hold for each observation type
  !-----------------------------------------------------------
  type t_bufr_inv
    !----------------------------------------------
    ! distribution of observations over input files
    !----------------------------------------------
    logical :: file    (m_source) = .false. ! file n holds obstype
    integer :: subseto (m_source) = 0       ! offset for sub-sets in file
    !-----------------------------------------------
    ! total number of BUFR records and subsets
    ! a subset corresponds to one observation report
    !-----------------------------------------------
    integer :: nrec               = 0       ! number of records
    integer :: nsubset            = 0       ! number of sub-sets
    !--------------------------------------------
    ! subsequent entries are different on each PE
    !--------------------------------------------
    integer :: subsets            = 0       ! no of sub-sets to skip on this PE
    integer :: subsetl            = 0       ! last sub-set to read on this PE
  end type t_bufr_inv
  !-----------------------------------------
  ! observation input file inventory
  ! seperate entry for each observation type
  !-----------------------------------------
  type(t_bufr_inv) ,save :: bufr_inv0       ! empty bufr_inv
  type(t_bufr_inv) ,save :: bufr_inv (n_ot) ! observation file inventory
  target                 :: bufr_inv
  !----------------------------------------------------------------------------
  !===============
  ! more constants
  !===============
  type (t_obs)   ,save :: empty_obs
  type (t_spot)  ,save :: empty_spot
  type (t_datum) ,save :: empty_body
  !----------------------------------------------------------------------------
  !===========
  ! interfaces
  !===========
!------------------------------------------------------------------------------
  interface construct
    module procedure construct_ob
    module procedure construct_obs
  end interface construct
!------------------------------------------------------------------------------
  interface destruct
    module procedure destruct_ob
    module procedure destruct_obs
    module procedure destruct_spots
    module procedure destruct_mc
    module procedure destruct_mcs
  end interface destruct
!------------------------------------------------------------------------------
  interface set_xuv
    module procedure set_xuv_icol
    module procedure set_xuv_imcol
    module procedure set_xuv_coord
    module procedure set_xuv_spot
  end interface
!------------------------------------------------------------------------------
  interface release_mem
    module procedure release_mem_1
    module procedure release_mem_0
  end interface release_mem
!------------------------------------------------------------------------------
  interface p_bcast
    module procedure p_bcast_source
    module procedure p_bcast_bufr_inv
    module procedure p_bcast_bufr_inv_1d
  end interface p_bcast
!------------------------------------------------------------------------------
  interface debug_spot
    module procedure debug_spot_sp
    module procedure debug_spot_obs
  end interface debug_spot
!------------------------------------------------------------------------------
  !========================
  ! public module variables
  !========================
  integer  ,save      :: po_context = -99
  integer  ,save      :: po_ilns    = -99
  logical  ,save      :: po_lmon    = .false.

!------------------------------------------------------------------------------
  !=========================
  ! private module variables
  !=========================
  real(wp) ,parameter :: pi      = 3.141592653589793238_wp

  integer  ,save      :: spot_id =      0
  integer  ,parameter :: n_spot  =   1000
  integer  ,parameter :: n_obs   = 100000
  integer  ,parameter :: n_par   =    500
  integer  ,parameter :: n_int   =    500
  integer  ,parameter :: n_sink  =    100
  !=======================
  ! Debug selected spot(s)
  !=======================
  integer, parameter :: nsd             =  5
  integer            :: spt_debug(nsd)  = -1   ! debug_spot is active for spot% id == spt_debug
  integer            :: spt_hd_debug(nsd)= -1  ! debug_spot is active for spot% hd %id == spt_hd_debug
  integer            :: usd             =  6   ! unit for spot_debug output
  character(len=24)  :: dpref           = ''   ! prefix for current spot for debug output
  integer            :: i_shd           = -1
  integer            :: i_sd            = -1
  logical            :: ldeb_all        = .false.
  logical            :: ldeb_spot       = .false.
!------------------------------------------------------------------------------
  !=========
  ! namelist
  !=========
  logical           :: read_bufr     = .false. ! read BUFR   files
  logical           :: read_NetCDF   = .false. ! read NetCDF files
  logical           :: read_cdfin    = .false. ! read COSMO observation files
  character(len=512):: dace_obs_op   = ''      ! use dace obs. operators in MEC (string)
  character(len=18) :: dace_op(n_ot) = ''      ! use dace obs. operators (array)
  integer           :: n_dace_op     = -1      ! length of non-empty entries in dace_op
  logical           :: obs_local     = .false. ! keep obs on gridpoint PE (MEC)
  character(len=128):: obs_path      = ''      ! path to read observation files
  character(len=64) :: bufr_files(63)= ''      ! names of BUFR files
  character(len=64) :: obs_files(120)= ''      ! observation file names (NetCDF)
  character(len=64) :: fdb_files (63)= ''      ! feedback (input) files
  integer           :: bufr_verb     = 0       ! verbose flag for BUFR decoding
  logical           :: bufr_pause    = .false.
  integer           :: netcdf_verb   = 0       ! verbosity level of NetCDF decoding
  logical           :: derive_dbkz   = .false. ! derive DBKZ if not present
  logical           :: fix_no_obs    = .true.  ! set F to revert no_obs-check
  logical           :: int_vh        = .true.  ! interpolation: 1.vert.,2.hor.
  logical           :: int_nn        = .false. ! hor. interp.: nearest neighb.
  integer           :: int_rad_hum   = 0       !  humidity variable for vertical interp. for radiances
                                               ! 0: RH, 1:QV
  real(kind=wp)     :: int_rad_hum_pblend(2)   ! for int_rad_hum==2: lower/upper height for blending
                                               ! RH-interpolated and Q-interpolated profile
                                               ! p>int_rad_hum_pblend(1) : only RH-interpolated profile
                                               ! int_rad_hum_pblend(1)>p>int_rad_hum_pblend(2) : blending
                                               ! p<int_rad_hum_pblend(2) : only Q-interpolated profile
  logical           :: vint_lin_t    = .false. ! linear vert.intp. for temp.
  logical           :: vint_lin_z    = .false. ! linear vert.intp. for geop.
  logical           :: vint_lin_uv   = .false. ! linear vert.intp. for wind
  logical           :: vint_lin_tov  = .false. ! linear vert.intp. for RTTOV
  logical           :: monitor_ff    = .false. ! monitor wind speed
  logical           :: monitor_dd    = .false. ! monitor wind direction
  integer           :: nwv_rad       = -1      ! radiance vertical interp. flag
  integer           :: vint_rttov    = -1      ! RTTOV interpolation mode
  integer           :: vint_rttov_b  = -1      ! RTTOV interpolation mode
  integer           :: ndv_rad       =  3      ! vertical differentiation flag
  integer           :: z2p_amv       =  0      ! AMV: derive pressure from height
  integer           :: z2p_hlos      =  0      ! AMV: derive pressure from height
  real(wp)          :: ptop_lapse    =  700._wp! pressure bounds for derivation
  real(wp)          :: pbot_lapse    =  850._wp! of lapse rate (t2m extrapolat.)
  integer           :: pcc_conv      =  1      ! per cent confidence usage flag
  integer           :: fdbk_split(mstep,n_ot)=0! split feedback files
                                               !   1: separate file for each satellite
                                               !   2: separate file for each isntrument grid

  integer           :: corme_conv    = -1      ! correction message handling:
                                               ! -1 NetCDF=off BUFR=off
                                               !  0 NetCDF=on  BUFR=off
                                               !  1 NetCDF=on  BUFR=on
                                               !  2 handle corme, ABCD separately
  namelist /OBSERVATIONS/ obs_path,    &
                          bufr_files,  &
                          obs_files,   &
                          fdb_files,   &
                          bufr_verb,   &
                          bufr_pause,  &
                          netcdf_verb, &
                          derive_dbkz, &
                          fix_no_obs,  &
                          int_vh,      &
                          int_nn,      &
                          int_rad_hum, &
                          int_rad_hum_pblend, &
                          vint_lin_t,  &
                          vint_lin_z,  &
                          vint_lin_uv, &
                          vint_lin_tov,&
                          monitor_ff,  &
                          monitor_dd,  &
                          nwv_rad,     &
                          vint_rttov,  &
                          vint_rttov_b,&
                          ndv_rad,     &
                          z2p_amv,     &
                          read_bufr,   &
                          read_NetCDF, &
                          read_cdfin,  &
                          dace_obs_op, &
                          obs_local,   &
                          qbit_conv,   &
                          pcc_conv,    &
                          corme_conv,  &
                          ptop_lapse,  &
                          pbot_lapse,  &
                          fdbk_split,  &
                          z2p_hlos,    &
                          spt_hd_debug,&
                          spt_debug,   &
                          usd
!==============================================================================
contains
!==============================================================================
  function ft_name (ft)
  integer ,intent(in) :: ft
  character (len=8)   :: ft_name
    select case (ft)
    case(FT_BUFR)
      ft_name = 'BUFR'
    case(FT_NETCDF)
      ft_name = 'NETCDF'
    case(FT_CDFIN)
      ft_name = 'CDFIN'
    case(FT_ODB)
      ft_name = 'ODB'
    case(FT_SATPP)
      ft_name = 'SATPP'
    case(FT_ROPIC)
      ft_name = 'ROPIC'
    case(FT_FEEDBACK)
      ft_name = 'FEEDBACK'
    case(FT_RAW_FEEDBACK)
      ft_name = 'RAW_FEEDBACK'
    case(FT_UNKNOWN)
      ft_name = 'UNKNOWN'
    case(FT_SPECIFIC)
      ft_name = 'SPECIFIC'
    case(FT_MISSING)
      ft_name = 'MISSING'
    case(FT_EMPTY)
      ft_name = 'EMPTY'
    case default
      call finish ('ft_name','invalid file type')
    end select
  end function ft_name
!------------------------------------------------------------------------------
  function rct_name (ot)
  integer ,intent(in) :: ot
  character (len=6)   :: rct_name
    select case (ot)
    case (SYNOP)
      rct_name = 'SYNOP'
    case (TEMP)
      rct_name = 'TEMP'
    case (TOVS)
      rct_name = 'TOVS'
    case (GPSRO)
      rct_name = 'GPSRO'
    case (AMV)
      rct_name = 'AMV'
    case (PSTEMP)
      rct_name = 'PSTEMP'
    case (AIREP)
      rct_name = 'AIREP'
    case (SATEM)
      rct_name = 'SATEM'
    case (GPSGB)
      rct_name = 'GPSGB'
    case (COSMO)
      rct_name = 'COSMO'
    case (RADAR)
      rct_name = 'RADAR'
    case (SOIL)
      rct_name = 'SOIL'
    case (WLIDAR)
      rct_name = 'WLIDAR'
    case (SCATT)
      rct_name = 'SCATT'
    case default
      write(0,*) 'rct_name: invalid observation type ',ot
      call finish ('rct_name','invalid observation type')
    end select
  end function rct_name
!------------------------------------------------------------------------------
  function rct_key (rct_name)
  integer                        :: rct_key
  character (len=*) ,intent(in)  :: rct_name
    select case (rct_name)
    case ('SYNOP')
      rct_key = SYNOP
    case ('TEMP')
      rct_key = TEMP
    case ('TOVS')
      rct_key = TOVS
    case ('GPSRO')
      rct_key = GPSRO
    case ('AMV')
      rct_key = AMV
    case ('PSTEMP')
      rct_key = PSTEMP
    case ('AIREP')
      rct_key = AIREP
    case ('SATEM')
      rct_key = SATEM
    case ('GPSGB')
      rct_key = GPSGB
    case ('COSMO')
      rct_key = COSMO
    case ('RADAR')
      rct_key = RADAR
    case ('SOIL')
      rct_key = SOIL
    case ('WLIDAR')
      rct_key = WLIDAR
    case ('SCATT')
      rct_key = SCATT
    case default
      call finish ('rct_key','invalid observation type: '//rct_name)
    end select
  end function rct_key
!------------------------------------------------------------------------------
  function oq_name (oq)
  integer ,intent(in) :: oq
  character (len=3)   :: oq_name
    select case (oq)
    case (OBS_H)     ! geopotential
      oq_name = 'h'
    case (OBS_TV)    ! virtual temperature
      oq_name = 'tv'
    case (OBS_RH)    ! relative humidity
      oq_name = 'rh'
    case (OBS_U)     ! wind component u
      oq_name = 'u'
    case (OBS_V)     ! wind component v
      oq_name = 'v'
    case (OBS_HS)    ! surface geopotential
      oq_name = 'hs'
    case (OBS_BT)    ! brightness temperature
      oq_name = 'bt'
    case (OBS_BDA)   ! bending angle
      oq_name = 'bda'
    case (OBS_Q)     ! specific humidity
      oq_name = 'q'
    case (OBS_DUM)   ! dummy sink variable
      oq_name = 'du'
    case (OBS_CTR)   ! additional control variable
      oq_name = 'ctr'
    case (OBS_CHI)   ! velocity potential
      oq_name = 'chi'
   case (OBS_RFR)    ! refractivity
      oq_name = 'rfr'
    case (OBS_T)     ! temperature
      oq_name = 't'
    case (OBS_FF)    ! wind speed
      oq_name = 'ff'
    case (OBS_HLOS)  ! hor. line-of-sight wind (projected component)
      oq_name = 'hls'
    case (OBS_ZPD)   ! slant/zenith path delay
      oq_name = 'zpd'
    case (OBS_MIXR)   ! slant/zenith path delay
      oq_name = 'mixr'
    case (OBS_DRH)   ! relative humidity, dummy (sink) variable
      oq_name = 'drh'
    case default
      write(0,*)   'oq_name:  invalid observed quantity:',oq
      call finish ('oq_name','invalid observed quantity')
    end select
  end function oq_name
!------------------------------------------------------------------------------
  elemental function oq_varno (oq)
  integer ,intent(in) :: oq
  integer             :: oq_varno
    select case (oq)
    case (OBS_H)     ! geopotential
      oq_varno = VN_Z
    case (OBS_TV)    ! virtual temperature
      oq_varno = VN_VT
    case (OBS_RH)    ! relative humidity
      oq_varno = VN_RH
    case (OBS_U)     ! wind component u
      oq_varno = VN_U
    case (OBS_V)     ! wind component v
      oq_varno = VN_V
    case (OBS_HS)    ! surface geopotential
      oq_varno = -1
    case (OBS_BT)    ! brightness temperature
      oq_varno = VN_RAWBT
    case (OBS_BDA)   ! bending angle
      oq_varno = VN_BENDANG
    case (OBS_Q)     ! specific humidity
      oq_varno = VN_Q
    case (OBS_DUM,OBS_DRH)   ! dummy sink variable
      oq_varno = -1
    case (OBS_CTR)   ! additional control variable
      oq_varno = -1
    case (OBS_CHI)   ! velocity potential
      oq_varno = -1
    case (OBS_T)     ! temperature
      oq_varno = VN_T
    case (OBS_RFR)   ! refractivity
      oq_varno = VN_REFR
    case (OBS_FF)    ! wind speed
      oq_varno = VN_FF
    case (OBS_HLOS)  ! hor. line-of-sight wind (projected component)
      oq_varno = VN_HLOS
    case (OBS_ZPD)   ! slant/zenith path delay
      oq_varno = VN_ZPD
    case (OBS_MIXR)  ! mixing ratio
      oq_varno = VN_MIXR
    case default
      oq_varno = -1
    end select
  end function oq_varno
!------------------------------------------------------------------------------
  elemental function varno_oq (varno)
  integer ,intent(in) :: varno
  integer             :: varno_oq
    select case (varno)
    case (VN_Z)         ! geopotential
      varno_oq = OBS_H
    case (VN_VT)        ! virtual temperature
      varno_oq = OBS_TV
    case (VN_RH)        ! relative humidity
      varno_oq = OBS_RH
    case (VN_U)         ! wind component u
      varno_oq = OBS_U
    case (VN_V)         ! wind component v
      varno_oq = OBS_V
    case (VN_RAWBT)     ! brightness temperature
      varno_oq = OBS_BT
    case (VN_BENDANG)   ! bending angle
      varno_oq = OBS_BDA
    case (VN_Q)         ! specific humidity
      varno_oq = OBS_Q
    case (VN_T)         ! temperature
      varno_oq = OBS_T
    case (VN_T2M)       ! 2m temperature
      varno_oq = OBS_T
    case (VN_REFR)      ! refractivity
      varno_oq = OBS_RFR
    case (VN_FF)        ! wind speed
      varno_oq = OBS_FF
    case (VN_HLOS)      ! hor. line-of-sight wind (projected component)
      varno_oq = OBS_HLOS
    case (VN_ZPD)       ! slant/zenith path delay
      varno_oq = OBS_ZPD
    case (VN_MIXR)      ! mixing ratio
      varno_oq = OBS_MIXR
    case default
      varno_oq = -1
    end select
  end function varno_oq
!------------------------------------------------------------------------------
  function oq_key (oq_name)
  integer                        :: oq_key
  character (len=3) ,intent(in)  :: oq_name
    select case (oq_name)
    case ('h')       ! geopotential
      oq_key = OBS_H
    case ('tv')      ! virtual temperature
      oq_key = OBS_TV
    case ('rh')      ! relative humidity
      oq_key = OBS_RH
    case ('u')       ! wind component u
      oq_key = OBS_U
    case ('v')       ! wind component v
      oq_key = OBS_V
    case ('hs')      ! surface geopotential
      oq_key = OBS_HS
    case ('bt')      ! brightness temperature
      oq_key = OBS_BT
    case ('bda')     ! bending angle
      oq_key = OBS_BDA
    case ('q')       ! specific humidity
      oq_key = OBS_Q
    case ('du')      ! dummy sink variable
      oq_key = OBS_DUM
    case ('ctr')     ! additional control variable
      oq_key = OBS_CTR
    case ('chi')     ! velocity potential
      oq_key = OBS_CHI
    case ('rfr')     ! refractivity
      oq_key = OBS_RFR
    case ('t')       ! temperature
      oq_key = OBS_T
    case ('ff')      ! wind speed
      oq_key = OBS_FF
    case ('hls')     ! hor. line-of-sight wind (projected component)
      oq_key = OBS_HLOS
    case ('zpd')     ! zenith path delay
      oq_key = OBS_ZPD
    case ('mixr')    ! mixing ratio
      oq_key = OBS_MIXR
    case default
      call finish ('oq_key','invalid observed quantity: '//oq_name)
    end select
  end function oq_key
!------------------------------------------------------------------------------
  function tsk_name (task)
  integer ,intent(in) :: task
  character (len=16)  :: tsk_name
    select case (task)
    case (TSK_INIT)
      tsk_name = 'TSK_INIT'
    case (TSK_READ)
      tsk_name = 'TSK_READ'
    case (TSK_SET_CHR)
      tsk_name = 'TSK_SET_CHR'
    case (TSK_SHRINK)
      tsk_name = 'TSK_SHRINK'
    case (TSK_SETUP_COLS)
      tsk_name = 'TSK_SETUP_COLS'
    case (TSK_SETUP_FUL0)
      tsk_name = 'TSK_SETUP_FUL0'
    case (TSK_SETUP_FULL)
      tsk_name = 'TSK_SETUP_FULL'
    case (TSK_R)
      tsk_name = 'TSK_R'
    case (TSK_Y)
      tsk_name = 'TSK_Y'
    case (TSK_YH)
      tsk_name = 'TSK_YH'
    case (TSK_H)
      tsk_name = 'TSK_H'
    case (TSK_K)
      tsk_name = 'TSK_K'
    case (TSK_XK)
      tsk_name = 'TSK_XK'
    case default
      write(tsk_name,'(a4,i8.8,4x)') 'TSK_',task
    end select
  end function tsk_name
!==============================================================================
  pure function unitvector (dlat, dlon) result (x)
    real(wp), intent(in) :: dlat, dlon  ! Latitude, longitude [degree]
    real(wp)             :: x(3)
    !---------------------------------
    ! return unit vector on the sphere
    !---------------------------------
    real(wp) ,parameter :: pdr = pi / 180._wp
    real(wp)            :: sinlat, coslat, sinlon, coslon
    sinlat = sin (pdr * dlat)
    coslat = cos (pdr * dlat)
    sinlon = sin (pdr * dlon)
    coslon = cos (pdr * dlon)
    x(1)   = coslon * coslat ! x
    x(2)   = sinlon * coslat ! y
    x(3)   =          sinlat ! z
  end function unitvector
!------------------------------------------------------------------------------
  elemental subroutine set_xuv_coord (col)
  type (t_coord) ,intent(inout) :: col
  !----------------------------------
  ! set up unit vectors on the sphere
  !----------------------------------
    real(wp) ,parameter :: pdr = pi / 180._wp
    real(wp)            :: sinlat, coslat, sinlon, coslon
    sinlat     =  sin (pdr * col% dlat)
    coslat     =  cos (pdr * col% dlat)
    sinlon     =  sin (pdr * col% dlon)
    coslon     =  cos (pdr * col% dlon)
    col% x (1) =  coslon * coslat ! x
    col% x (2) =  sinlon * coslat ! y
    col% x (3) =           sinlat ! z
    col% du(1) = -sinlon          ! u - vector
    col% du(2) =  coslon
    col% du(3) =  0._wp
    col% dv(1) = -coslon * sinlat ! v - vector
    col% dv(2) = -sinlon * sinlat
    col% dv(3) =           coslat
  end subroutine set_xuv_coord
!------------------------------------------------------------------------------
  elemental subroutine set_lluv (col)
  type (t_coord) ,intent(inout) :: col
  !---------------------------------
  ! set up lat lon and du, dv from x
  !---------------------------------
    real(wp) ,parameter :: rdp = 180._wp / pi
    real(wp)            :: sinlat, coslat, sinlon, coslon
!   sinlat      =  sin (pdr * col% dlat)
!   coslat      =  cos (pdr * col% dlat)
!   sinlon      =  sin (pdr * col% dlon)
!   coslon      =  cos (pdr * col% dlon)
!   col% x (1) =  coslon * coslat ! x
!   col% x (2) =  sinlon * coslat ! y
!   col% x (3) =           sinlat ! z
    sinlat =       col% x (3)
    coslat = sqrt (col% x (1) ** 2 + col% x (2) ** 2)
    coslon =       col% x (1) / coslat
    sinlon =       col% x (2) / coslat
    col% dlat  = rdp * atan2 (sinlat, coslat)
    col% dlon  = rdp * atan2 (sinlon, coslon)
    col% du(1) = -sinlon          ! u - vector
    col% du(2) =  coslon
    col% du(3) =  0._wp
    col% dv(1) = -coslon * sinlat ! v - vector
    col% dv(2) = -sinlon * sinlat
    col% dv(3) =           coslat
  end subroutine set_lluv
!------------------------------------------------------------------------------
  elemental subroutine set_xuv_icol (col)
  type (t_icol) ,intent(inout) :: col
    call set_xuv_coord (col% c)
  end subroutine set_xuv_icol
!------------------------------------------------------------------------------
  elemental subroutine set_xuv_imcol (col)
  type (t_imcol) ,intent(inout) :: col
    call set_xuv_coord (col% c)
  end subroutine set_xuv_imcol
!------------------------------------------------------------------------------
  elemental subroutine set_xuv_spot (s)
  type (t_spot) ,intent(inout) :: s
    !-------------------------------
    ! set unit vectors on the sphere
    !-------------------------------
    real(wp) :: zenith,  azimuth
    call set_xuv_coord (s% col% c)
    !---------------------------------------
    ! set solar zenith angle if not yet done
    !---------------------------------------
    if (s% sozen == invalid .or. s% soazi == invalid) then
      call solar (s% actual_time, s% col% c% dlat, s% col% c% dlon, &
                  azimuth, zenith                                   )
      if (s% sozen == invalid) s% sozen = zenith
      if (s% soazi == invalid) s% soazi = azimuth
    endif
  end subroutine set_xuv_spot
!==============================================================================
  elemental function mdlsfc_frac (spt) result (mdlsfc)
  !------------------------------------------------------------
  ! derive model surface flag from 3dvar land/ice/snow fraction
  !------------------------------------------------------------
  type (t_spot) ,intent(in) :: spt    ! 3dvar report meta data
  integer                   :: mdlsfc ! model surface flag

    mdlsfc                         = 0
    if   (spt% sl_bg  < 1._wp) then
                               mdlsfc = ibset (mdlsfc, MS_SEA    )
      if (spt% fi_bg  > 0._wp) mdlsfc = ibset (mdlsfc, MS_ICE    )
      if (spt% fi_bg  < 1._wp) mdlsfc = ibset (mdlsfc, MS_NO_ICE )
    endif
    if   (spt% sl_bg  > 0._wp) then
                               mdlsfc = ibset (mdlsfc, MS_LAND   )
!+++  fractional snow cover not considered so far +++
      if (spt% hs_bg  > 0._wp) mdlsfc = ibset (mdlsfc, MS_SNOW   )
      if (spt% hs_bg == 0._wp) mdlsfc = ibset (mdlsfc, MS_NO_SNOW)
    endif

  end function mdlsfc_frac
!------------------------------------------------------------------------------
  elemental subroutine frac_mdlsfc (spt, mdlsfc)
  !------------------------------------------------------------
  ! derive 3dvar land/ice/snow fraction from model surface flag
  !------------------------------------------------------------
  type (t_spot) ,intent(inout) :: spt    ! 3dvar report meta data
  integer       ,intent(in)    :: mdlsfc ! model surface flag

    spt%                                  sl_bg =              0.5_wp
    if   (btest (mdlsfc, MS_LAND))   spt% sl_bg = spt% sl_bg + 0.5_wp
    if   (btest (mdlsfc, MS_SEA ))   spt% sl_bg = spt% sl_bg - 0.5_wp

    spt%                                  fi_bg = 0.0_wp
    if   (btest (mdlsfc, MS_ICE )) then
                                     spt% fi_bg = 1.0_wp
      if (btest (mdlsfc, MS_NO_ICE)) spt% fi_bg = 0.5_wp
    endif

!   snow not considered so far !!

  end subroutine frac_mdlsfc
!==============================================================================
  subroutine col2lev (obs)
  type (t_obs) ,intent (inout) :: obs
  !=========================================
  ! compress column observation meta data to
  !          level  observation meta data
  !=========================================
    integer               :: i          ! index: spots
    integer               :: k          ! index: observations
    integer               :: l          ! index: levels
    real(wp)              :: z          ! vertical coordinate
    integer               :: n_lev      ! number of levels in obs data type
    integer               :: n_lev_spot ! number of levels per spot
    type(t_spot) ,pointer :: s
    !-----------------------------------
    ! first loop: count number of levels
    !-----------------------------------
    n_lev = 0
    do i = 1, obs% n_spot
      s => obs% spot(i)
      n_lev_spot = 0
      z = -huge(z)
      s% l% i = n_lev
      do k = s% i% i + 1, s% i% i + s% i% n
        if (z /= obs% lev(k)) then
          n_lev      = n_lev      + 1
          n_lev_spot = n_lev_spot + 1
          z = obs% lev(k)
        endif
      end do
      s% l% n = n_lev_spot
    end do
    !------------------------------------------
    ! allocate level component with proper size
    !------------------------------------------
    if (associated (obs% levs)) deallocate (obs% levs)
    allocate (obs% levs (n_lev))
    obs% n_lev = n_lev
    !----------------------------------------------
    ! second loop: derive and store level meta data
    !----------------------------------------------
    l = 0
    do i = 1, obs% n_spot
      s => obs% spot(i)
      z = -huge(z)
      do k = s% i% i + 1, s% i% i + s% i% n
        if (z /= obs% lev(k)) then
          l = l + 1
          z = obs% lev(k)
          obs% levs(l)% its  = 0   ! observation identifiers, 'iored'
          obs% levs(l)% z    = z   ! transformed coordinate (ln(p))
          obs% levs(l)% i% i = k-1 ! index to 'int'  component
          obs% levs(l)% i% n = 1   ! no.observations per level
          obs% levs(l)% is   = i   ! spot index
        else
          obs% levs(l)% i% n = obs% levs(l)% i% n + 1
        endif
        obs% levs(l)% its = ior(obs% levs(l)% its, obs% t_int (k))
      end do
    end do
  end subroutine col2lev
!==============================================================================
  subroutine set_byte_size
    if (obs_bytes == 0) then
      obs_bytes  = size (transfer (empty_obs  ,(/' '/)))
      spot_bytes = size (transfer (empty_spot ,(/' '/)))
      body_bytes = size (transfer (empty_body ,(/' '/)))
      int_bytes  = size (transfer (1          ,(/' '/)))
      wp_bytes   = size (transfer (1._wp      ,(/' '/)))

!     if (dace% lpio) then
!        write(0,*)
!        write(0,*) "set_byte_size: obs_bytes  =", obs_bytes
!        write(0,*) "set_byte_size: spot_bytes =", spot_bytes
!        write(0,*) "set_byte_size: body_bytes =", body_bytes
!        write(0,*) "set_byte_size: int_bytes  =", int_bytes
!        write(0,*) "set_byte_size: wp_bytes   =", wp_bytes
!        write(0,*)
!        call flush (0)
!     end if
    endif
  end subroutine set_byte_size
!------------------------------------------------------------------------------
  subroutine mean_x_obs (obs)
  type(t_obs) ,intent(inout) :: obs (:)
  !--------------------------------------------------
  ! determine mean coordinates of observations in box
  !--------------------------------------------------
    target                :: obs
    type (t_obs) ,pointer :: o
    integer               :: ib, is
    real(wp)              :: xlen

    do ib = 1, size (obs)
      o => obs(ib)
      o% bcc = 0._wp
      do is = 1, o% n_spot
        o% bcc = o% bcc + o% spot(is)% o% n * o% spot(is)% col% c% x
      end do
      xlen = sqrt (sum (o% bcc**2))
      if (xlen > 0._wp) o% bcc = o% bcc / xlen
    end do
  end subroutine mean_x_obs
!------------------------------------------------------------------------------
  subroutine reorder_obs (src, dst, ibx)
  type(t_obs) ,intent(inout)        :: src      ! source
  type(t_obs) ,intent(inout)        :: dst      ! destination
  type(t_box) ,intent(in) ,optional :: ibx (:)  ! destination boxes and PEs
  !----------------------------------------------------
  ! sort observations according to their destination PE
  !----------------------------------------------------

    type(t_obs) :: d(1)
    integer     :: ib, ns

    if (present(ibx)) src% spot% srbx = ibx
    d(1) = dst
    do ib=1, dace% npe
      ns = src% n_spot
      call scatter_obs (src,d, mask=src%spot(:ns)%srbx% pe==ib-1, allpe=.true.)
    end do
    dst = d(1)
  end subroutine reorder_obs
!------------------------------------------------------------------------------
  subroutine filter_obs (src, dst, ibx)
  type(t_obs) ,intent(inout)        :: src      ! source
  type(t_obs) ,intent(inout)        :: dst      ! destination
  type(t_box) ,intent(in) ,optional :: ibx (:)  ! destination boxes and PEs
  !--------------------------------------------------------------------
  ! filter observations according to valid destination PE indices (>=0)
  !--------------------------------------------------------------------

    type(t_obs) :: d(1)

    if (present(ibx)) src% spot% srbx = ibx
    d(1) = dst
    call scatter_obs (src, d, mask=src%spot%srbx% pe>=0, allpe=.true.)
    dst = d(1)
  end subroutine filter_obs
!------------------------------------------------------------------------------
  subroutine scatter_obs (src, dst, ibx, mask, allpe)
  type(t_obs) ,intent(inout)        :: src
  type(t_obs) ,intent(inout)        :: dst  (:)
  type(t_box) ,intent(in) ,optional :: ibx  (:)
  logical     ,intent(in) ,optional :: mask (:)
  logical     ,intent(in) ,optional :: allpe

    integer              :: i, j, l, n, it, if, si
    logical              :: msk (src% n_spot)
    logical              :: thispe
    type(t_obs) ,pointer :: d
    type(t_obs)          :: old
    integer              :: nd , no
    target               :: dst
    integer              :: n_new

    if (present(ibx)) src% spot% srbx = ibx
    thispe = .true. ;if(present(allpe)) thispe = .not.allpe
    n = src% n_spot
    !----------------------------------------
    ! loop over elements of destination array
    !----------------------------------------
    dst% bcast  = .not. thispe
    dst% n_time = src% n_time
    do l = 1, size (dst)

      if (present(mask)) then
        msk = mask
      else
        msk = (src% spot(1:n)% srbx% box == l)
      end if

      old   =  dst(l)
      d     => dst(l)
      n_new =  count(msk)
      d% ibox   = l
      d% n_spot = d% n_spot + n_new
      d% n_obs  = d% n_obs  + sum  (src% spot(1:n)% o% n   ,msk)
      d% n_int  = d% n_int  + sum  (src% spot(1:n)% i% n   ,msk)
!!!   d% n_ctr  = d% n_ctr  + sum  (src% spot(1:n)% c% n   ,msk)
      d% n_par  = d% n_par  + sum  (src% spot(1:n)% p% n   ,msk)
      d% n_sink = d% n_sink + sum  (src% spot(1:n)% d% n   ,msk)
      d% n_bcp  = d% n_bcp  + sum  (src% spot(1:n)% bcp% n ,msk)

      !----------------------------------
      ! skip remaining part if allpe == F
      !----------------------------------
      if (thispe) then
        msk = (msk .and. src% spot(1:n)% srbx% pe == dace% pe)
        if (count(msk) /= n_new) then
          if (count(msk)==0) cycle
          call finish ('scatter_obs','count(msk) /= n_new')
        end if
      end if
      !--------------------------------------------------------------
      ! (re)allocate components with correct size (observation space)
      !--------------------------------------------------------------
      nd = d%   n_obs
      no = old% n_obs
      if (d% n_spot /= old% n_spot .or. .not. associated (old% spot)) then
        allocate (d% spot (d% n_spot))
        if (associated(old% spot)) then
          d%spot(:old% n_spot) = old% spot(:old% n_spot)
          deallocate (old% spot)
        endif
      else
        d% spot => old% spot
      end if
      if (associated (src% body )) then
        if (nd /= no) then
          allocate (d% body  (nd))
          if (associated (old%body)) then
!$omp parallel workshare
            d% body(:no) = old% body(:no)
!$omp end parallel workshare
            deallocate (old% body)
          endif
        else
          d% body => old% body
        end if
      endif
      if (associated (src% varno)) then
        if (nd /= no) then
          allocate (d% varno  (nd))
          if (associated (old%varno)) then
            d% varno(:no) = old% varno(:no)
            deallocate (old% varno)
          endif
        else
          d% varno => old% varno
        end if
      endif
      if (associated (src% olev )) then
        if (nd /= no) then
          allocate (d% olev  (nd))
          if (associated (old%olev)) then
            d% olev(:no) = old% olev(:no)
            deallocate (old% olev)
          endif
        else
          d% olev => old% olev
        end if
      endif
      if (associated (src% bger  )) then
        if (nd /= no) then
          allocate (d% bger  (nd))
          if (associated (old%bger)) then
            d% bger(:no) = old% bger(:no)
            deallocate (old% bger)
          endif
        else
          d% bger => old% bger
        end if
      endif

      !--------------------------------------------
      ! copy component elements (observation space)
      !--------------------------------------------
      it = old% n_obs
      j  = old% n_spot
      do i=1, n
        if(.not.msk(i)) cycle
        j  = j + 1
        si = src% spot(i)% o% n
        if = src% spot(i)% o% i
        d% spot(j)       = src% spot(i)
        d% spot(j)% o% i = it
        d% spot(j)% o% n = si
        d% spot(j)% pe_dest => NULL()
        d% spot(j)% imcol   => NULL()
        if(associated(src% body))d% body(it+1:it+si)=src% body(if+1:if+si)
        if(associated(src%varno))d%varno(it+1:it+si)=src%varno(if+1:if+si)
        if(associated(src% olev))d% olev(it+1:it+si)=src% olev(if+1:if+si)
        if(associated(src% bger))d% bger(it+1:it+si)=src% bger(if+1:if+si)
        it = it + si
      end do

      !-----------------------------
      ! same for interpolation space
      !-----------------------------
      nd = d%   n_int
      no = old% n_int
      if (associated (src% t_int)) then
        if (nd /= no) then
          allocate (d% t_int (nd))
          allocate (d% lev   (nd))
          if (associated (old%t_int)) then
            d% t_int(:no) = old% t_int(:no)
            deallocate (old% t_int)
          endif
          if (associated (old%lev)) then
            d% lev(:no) = old% lev(:no)
            deallocate (old% lev)
          endif
        else
          d% t_int => old% t_int
          d% lev   => old% lev
        end if
      endif

      if (associated (src% bgeri  )) then
        if (nd /= no) then
          allocate (d% bgeri  (nd))
          if (associated (old%bgeri)) then
            d% bgeri(:no) = old% bgeri(:no)
            deallocate (old% bgeri)
          endif
        else
          d% bgeri => old% bgeri
        end if
      endif

      if (associated (src% bgi  )) then
        if (nd /= no) then
          allocate (d% bgi  (nd))
          if (associated (old%bgi)) then
            d% bgi(:no) = old% bgi(:no)
            deallocate (old% bgi)
          endif
        else
          d% bgi => old% bgi
        end if
      endif

      it = old% n_int
      j  = old% n_spot
      do i=1, n
        if(.not.msk(i)) cycle
        j = j + 1
        si = src% spot(i)% i% n
        if = src% spot(i)% i% i
        d% spot(j)% i% i = it
        d% spot(j)% i% n = si
        if(associated(src% t_int)) then
          d% t_int(it+1:it+si)= src% t_int (if+1:if+si)
          d% lev  (it+1:it+si)= src% lev   (if+1:if+si)
        endif
        if(associated(src% bgeri))d% bgeri(it+1:it+si)= src% bgeri (if+1:if+si)
        if(associated(src% bgi  ))d% bgi  (it+1:it+si)= src% bgi   (if+1:if+si)
        it = it + si
      end do

      !------------------------------------
      ! same for bias correction predictors
      !------------------------------------
      if (associated (src% bcpred)) then
        nd = d%   n_bcp
        no = old% n_bcp
        if (nd /= no) then
          allocate (d% bcpred (nd))
          if (associated (old%bcpred)) then
            d% bcpred(:no) = old% bcpred(:no)
            deallocate (old% bcpred)
          endif
        else
          d% bcpred => old% bcpred
        end if
        it = old% n_bcp
        j  = old% n_spot
        do i=1, n
          if(.not.msk(i)) cycle
          j = j + 1
          si = src% spot(i)% bcp% n
          if = src% spot(i)% bcp% i
          d% spot(j)% bcp% i = it
          d% spot(j)% bcp% n = si
          if (si > 0) &
          d% bcpred(it+1:it+si) = src% bcpred(if+1:if+si)
          it = it + si
        end do
      endif

      !----------------------------------
      ! same for operator specific params
      !----------------------------------
      if (associated (src% par  )) then
        nd = d%   n_par
        no = old% n_par
        if (nd /= no) then
          allocate (d% par (nd))
          if (associated (old%par)) then
!!$omp parallel workshare
            d% par(:no) = old% par(:no)
!!$omp end parallel workshare
            deallocate (old% par)
          endif
        else
          d% par => old% par
        end if
        it = old% n_par
        j  = old% n_spot
        do i=1, n
          if(.not.msk(i)) cycle
          j = j + 1
          si = src% spot(i)% p% n
          if = src% spot(i)% p% i
          d% spot(j)% p% i = it
          d% spot(j)% p% n = si
          if (si > 0) &! <-- for Intel compiler runtime bounds/pointer checks
          d% par   (it+1:it+si) = src% par   (if+1:if+si)
          it = it + si
        end do
      endif

      !------------------------
      ! same for sink variables
      !------------------------
      if (associated (src% sink  )) then
        nd = d%   n_sink
        no = old% n_sink
        if (nd /= no) then
          allocate (d% sink (nd))
          if (associated (old%sink)) then
            d% sink(:no) = old% sink(:no)
            deallocate (old% sink)
          endif
        else
          d% sink => old% sink
        end if
        it = old% n_sink
        j  = old% n_spot
        do i=1, n
          if(.not.msk(i)) cycle
          j = j + 1
          si = src% spot(i)% d% n
          if = src% spot(i)% d% i
          d% spot(j)% d% i = it
          d% spot(j)% d% n = si
          if (si > 0) &! silence valgrind
          d% sink   (it+1:it+si) = src% sink   (if+1:if+si)
          it = it + si
        end do
      endif

    end do

  end subroutine scatter_obs
!------------------------------------------------------------------------------
  subroutine gather_obs (src, dst, isrc)
  type(t_obs) ,intent(in)    :: src(:)
  type(t_obs) ,intent(inout) :: dst
  integer     ,intent(in)    :: isrc(:)

    integer :: i, l, n, it, if, si
    logical :: mask(dst% n_spot)

    type (t_index) :: io, ii, ip

    do l = 1, size (src)
      mask = isrc(1:dst% n_spot) == l

      n = dst% n_spot

      if (.not.associated(dst% spot)) allocate (dst% spot (dst% n_spot))
      if = 1
      do i=1,n
        if(.not.mask(i)) cycle
        io = dst% spot(i)% o
        ii = dst% spot(i)% i
        ip = dst% spot(i)% p
        dst% spot(i) = src(l)% spot(if)
        dst% spot(i)% o = io
        dst% spot(i)% i = ii
        dst% spot(i)% p = ip
        if = it + 1
      end do

      if (associated (src(l)% body ) .and..not. associated (dst% body ))      &
        allocate (dst% body  (dst% n_obs))
      if (associated (src(l)% varno) .and..not. associated (dst% varno))      &
        allocate (dst% varno (dst% n_obs))
      if (associated (src(l)% olev ) .and..not. associated (dst% olev))       &
        allocate (dst% olev  (dst% n_obs))
      if (associated (src(l)% bger ) .and..not. associated (dst% bger))       &
        allocate (dst% bger   (dst% n_obs))
      do i=1, dst% n_spot
        if(.not.mask(i)) cycle
        si = src(l)% spot(i)% o% n
        if = src(l)% spot(i)% o% i
        it = dst  % spot(i)% o% i
        dst% spot(i)% o% i = it
        if (associated(src(l)%body )) &
          dst%body (it+1:it+si) = src(l)%body (if+1:if+si)
        if (associated(src(l)%varno)) &
          dst%varno(it+1:it+si) = src(l)%varno(if+1:if+si)
        if (associated(src(l)%  olev)) &
          dst% olev(it+1:it+si) = src(l)% olev(if+1:if+si)
        if (associated(src(l)%  bger)) &
          dst% bger(it+1:it+si) = src(l)% bger(if+1:if+si)
      end do

      if (associated (src(l)% bgeri) .and..not. associated (dst% bgeri)) &
        allocate (dst% bgeri  (dst% n_int))
      if (associated (src(l)% bgi  ) .and..not. associated (dst% bgi  )) &
        allocate (dst% bgi    (dst% n_int))
      if (associated (src(l)% t_int)) then
        do i=1, dst% n_spot
          if(.not.mask(i)) cycle
          si = src(l)% spot(i)% i% n
          if = src(l)% spot(i)% i% i
          it = dst  % spot(i)% i% i
          dst% spot(i)% i% i = it
          if (associated (src(l)% t_int)) &
            dst% t_int (it+1:it+si) = src(l)% t_int (if+1:if+si)
          if (associated (src(l)% lev)) &
            dst% lev   (it+1:it+si) = src(l)% lev   (if+1:if+si)
          if (associated (src(l)% bgeri)) &
            dst% bgeri (it+1:it+si) = src(l)% bgeri (if+1:if+si)
          if (associated (src(l)% bgi  )) &
            dst% bgi   (it+1:it+si) = src(l)% bgi   (if+1:if+si)
        end do
      endif

!!!      if (associated (src(l)% t_ctr)) then
!!!        do i=1, dst% n_spot
!!!          if(.not.mask(i)) cycle
!!!          si = src(l)% spot(i)% c% n
!!!          if = src(l)% spot(i)% c% i
!!!          it = dst  % spot(i)% c% i
!!!          dst% spot(i)% c% i = it
!!!          dst% t_ctr (it+1:it+si) = src(l)% t_ctr (if+1:if+si)
!!!        end do
!!!      endif

      if (associated (src(l)% par  )) then
        do i=1, dst% n_spot
          if(.not.mask(i)) cycle
          si = src(l)% spot(i)% p% n
          if = src(l)% spot(i)% p% i
          it = dst  % spot(i)% p% i
          dst% spot(i)% p% i = it
          dst% par   (it+1:it+si) = src(l)% par   (if+1:if+si)
        end do
      endif
    end do
  end subroutine gather_obs
!------------------------------------------------------------------------------
  subroutine join_obs (src, dst)
  !----------------------------------------------------------------
  ! join observations from multiple boxes (elements of type(t_obs))
  ! into one box (element)
  !----------------------------------------------------------------
  type(t_obs) ,intent(in)  :: src(:)
  type(t_obs) ,intent(out) :: dst

    integer :: ib, is, ip, io, il, ii, n

    !-----------------------
    ! return if src is empty
    !-----------------------
    if (size(src) == 0) return

    !---------------------
    ! copy some parameters
    !---------------------
    dst% pe = src(1)% pe

    !-----------------------------------
    ! set array sizes in destination box
    !-----------------------------------
    dst% n_spot = sum (src% n_spot)
    dst% n_spt  = sum (src% n_spt)
    dst% n_par  = sum (src% n_par)
    dst% n_obs  = sum (src% n_obs)
    dst% n_lev  = sum (src% n_lev)
    dst% n_int  = sum (src% n_int)

    !---------
    ! allocate
    !---------
    if (associated (src(1)% spot )) allocate (dst% spot  (dst% n_spot))
    if (dst% n_par /= 0           ) allocate (dst% par   (dst% n_par ))
    if (associated (src(1)% body )) allocate (dst% body  (dst% n_obs ))
    if (associated (src(1)% varno)) allocate (dst% varno (dst% n_obs ))
    if (associated (src(1)% bger )) allocate (dst% bger  (dst% n_obs ))
    if (associated (src(1)% s_vqc)) allocate (dst% s_vqc (dst% n_obs ))
    if (associated (src(1)% f_vqc)) allocate (dst% f_vqc (dst% n_obs ))
    if (associated (src(1)% olev )) allocate (dst% olev  (dst% n_obs ))
    if (associated (src(1)% levs )) allocate (dst% levs  (dst% n_lev ))
    if (associated (src(1)% t_int)) allocate (dst% t_int (dst% n_int ))
    if (associated (src(1)% lev  )) allocate (dst% lev   (dst% n_int ))
    if (associated (src(1)% bgeri)) allocate (dst% bgeri (dst% n_int ))
    if (associated (src(1)% bgi  )) allocate (dst% bgi   (dst% n_int ))

    !-----
    ! copy
    !-----
    is = 0
    ip = 0
    io = 0
    il = 0
    ii = 0
    do ib = 1, size (src)
      n  = src(ib)% n_spot

      if (associated (dst% spot )) dst% spot  (is+1:is+n) = src(ib)% spot (1:n)
      is = is + n
      n  = src(ib)% n_par

      if (n > 0)                   dst% par   (ip+1:ip+n) = src(ib)% par  (1:n)
      ip = ip + n
      n  = src(ib)% n_obs

      if (n > 0) then ! <-- for Intel compiler runtime bounds/pointer checks
      if (associated (dst% body )) dst% body  (io+1:io+n) = src(ib)% body (1:n)
      if (associated (dst% varno)) dst% varno (io+1:io+n) = src(ib)% varno(1:n)
      if (associated (dst% bger )) dst% bger  (io+1:io+n) = src(ib)% bger (1:n)
      if (associated (dst% s_vqc)) dst% s_vqc (io+1:io+n) = src(ib)% s_vqc(1:n)
      if (associated (dst% f_vqc)) dst% f_vqc (io+1:io+n) = src(ib)% f_vqc(1:n)
      if (associated (dst% olev )) dst% olev  (io+1:io+n) = src(ib)% olev (1:n)
      end if
      io = io + n
      n  = src(ib)% n_lev

      if (associated (dst% levs )) dst% levs  (il+1:il+n) = src(ib)% levs (1:n)
      il = il + n
      n  = src(ib)% n_int

      if (n > 0) then ! <-- for Intel compiler runtime bounds/pointer checks
      if (associated (dst% t_int)) dst% t_int (ii+1:ii+n) = src(ib)% t_int(1:n)
      if (associated (dst% lev  )) dst% lev   (ii+1:ii+n) = src(ib)% lev  (1:n)
      if (associated (dst% bgeri)) dst% bgeri (ii+1:ii+n) = src(ib)% bgeri(1:n)
      if (associated (dst% bgi  )) dst% bgi   (ii+1:ii+n) = src(ib)% bgi  (1:n)
      end if
      ii = ii + n
    end do

    !-----------------
    ! set index arrays
    !-----------------
    ip = 0
    io = 0
    il = 0
    ii = 0
    do is = 1, dst% n_spot
      dst% spot(is)% o% i = io; io = io + dst% spot(is)% o% n
      dst% spot(is)% p% i = ip; ip = ip + dst% spot(is)% p% n
      dst% spot(is)% l% i = il; il = il + dst% spot(is)% l% n
      dst% spot(is)% i% i = ii; ii = ii + dst% spot(is)% i% n
    end do

  end subroutine join_obs
!==============================================================================
  subroutine destruct_spots (spots)
  type (t_spot) ,intent(inout) :: spots(:)
    integer :: i
    do i=1,size(spots)
      if (associated(spots(i)% imcol))   deallocate (spots(i)% imcol)
      if (associated(spots(i)% pe_dest)) deallocate (spots(i)% pe_dest)
    end do
  end subroutine destruct_spots
!------------------------------------------------------------------------------
  subroutine destruct_obs (obs)
  type (t_obs) ,intent(inout) :: obs (:)
    integer :: i
    do i=1,size(obs)
      call destruct (obs(i))
    end do
  end subroutine destruct_obs
!------------------------------------------------------------------------------
  subroutine destruct_ob (obs)
  type (t_obs) ,intent(inout) :: obs
    type (t_obs) :: empty
    if (associated (obs% spot )) then
      call destruct (obs% spot)
      deallocate    (obs% spot)
    endif
    call destruct (obs% mc)
    if (associated (obs% body  )) deallocate (obs% body  )
    if (associated (obs% varno )) deallocate (obs% varno )
    if (associated (obs% bger  )) deallocate (obs% bger  )
    if (associated (obs% s_vqc )) deallocate (obs% s_vqc )
    if (associated (obs% f_vqc )) deallocate (obs% f_vqc )
    if (associated (obs% olev  )) deallocate (obs% olev  )
    if (associated (obs% t_int )) deallocate (obs% t_int )
    if (associated (obs% bgeri )) deallocate (obs% bgeri )
    if (associated (obs% bgi   )) deallocate (obs% bgi   )
    if (associated (obs% levs  )) deallocate (obs% levs  )
    if (associated (obs% lev   )) deallocate (obs% lev   )
!!! if (associated (obs% t_ctr )) deallocate (obs% t_ctr )
    if (associated (obs% par   )) deallocate (obs% par   )
    if (associated (obs% bcpred)) deallocate (obs% bcpred)
    if (associated (obs% sink  )) deallocate (obs% sink  )
    obs = empty
  end subroutine destruct_ob
!------------------------------------------------------------------------------
  subroutine destruct_mcs (mc)
  type (t_mcols) ,intent(inout) :: mc (:)
    integer :: i
    do i=1,size(mc)
      call destruct (mc(i))
    end do
  end subroutine destruct_mcs
!------------------------------------------------------------------------------
  subroutine destruct_mc  (mc)
  type (t_mcols) ,intent(inout) :: mc
    if (associated (mc% c  )) deallocate (mc% c  )
    if (associated (mc% idx)) deallocate (mc% idx)
  end subroutine destruct_mc
!==============================================================================
  subroutine construct_ob (obs)
  type (t_obs) ,intent(out) :: obs
    allocate(obs% spot(0))
  end subroutine construct_ob
!------------------------------------------------------------------------------
  subroutine construct_obs (obs)
  type (t_obs) ,intent(out) :: obs (:)
    integer :: i
    do i=1,size(obs)
      allocate(obs(i)% spot(0))
    end do
  end subroutine construct_obs
!==============================================================================
  subroutine unique_spot (obs)
  type (t_obs) ,intent(inout) :: obs
  !-------------------------------------------------
  ! Set unique spot_ids in the parallel environment.
  ! To be called after the last call to 'new_spot'.
  !-------------------------------------------------
    integer :: n_spots (0:dace% npe)
    integer :: pe

    if (spot_id<0) call finish('uniqe_spot','called twice')
    call p_allgather (spot_id, n_spots(1:))
    n_spots(0) = 0
    do pe = 1, dace% npe
      n_spots(pe) = n_spots(pe) + n_spots(pe-1)
    end do
    obs% spot(1:obs% n_spot)% id = obs% spot(1:obs% n_spot)% id + n_spots(dace% pe)
!   spot_id = -1  ! set to -1 to check for multiple calls
    spot_id =  0
  end subroutine unique_spot
!------------------------------------------------------------------------------
  subroutine new_spot (obs, n, set_id)
  type (t_obs) ,intent(inout) :: obs    ! container for observations
  integer      ,intent(in)    :: n      ! number of reports to add
  logical      ,intent(in)    :: set_id ! flag to set unique report id
  !---------------------------------------------------------------
  ! reserve memory for a new report
  ! meta data only, subsequently call new_spot for the data itself
  !---------------------------------------------------------------
    type (t_spot) ,pointer :: tmp (:)
    integer                :: i
    if (.not. associated(obs% spot)) allocate (obs% spot (n_spot))
    if (obs% n_spot + n > size(obs% spot)) then
      tmp => obs% spot
      if (n==1) then
        allocate (obs% spot ((obs% n_spot + n) * 2 ))
      else
        allocate (obs% spot ((obs% n_spot + n)     ))
      endif
      obs% spot (:obs% n_spot) = tmp (:obs% n_spot)
      deallocate (tmp)
    endif
    if (set_id) then
      if (spot_id >= 0) then
        do i=1,n
          spot_id = spot_id + 1
          obs% spot(obs% n_spot+i)% id = spot_id
        end do
      endif
    endif
    obs% n_spot = obs% n_spot + n
  end subroutine new_spot
!------------------------------------------------------------------------------
  subroutine release_mem_1 (obs, sort)
  type (t_obs) ,intent(inout)           :: obs (:) ! observation data type
  logical      ,intent(in)    ,optional :: sort    ! flag to reorder observs.
    integer :: i
    do i=1,size(obs)
      call release_mem (obs(i), sort=sort)
    end do
  end subroutine release_mem_1
!------------------------------------------------------------------------------
  subroutine release_mem_0 (obs, keep, sort)
  type (t_obs) ,intent(inout)           :: obs     ! observation data type
  logical      ,intent(in)    ,optional :: keep(:) ! mask for reports to keep
  logical      ,intent(in)    ,optional :: sort    ! flag to reorder observs.
  !-----------------------------------------
  ! release unused memory, update obs% types
  !-----------------------------------------
!   integer                     :: i
    type (t_index) ,allocatable :: ix  (:)  ! index array
    integer        ,allocatable :: key (:)  ! sort index
    integer        ,allocatable :: ipm (:)  ! permutation array
    logical                     :: changed  ! .true. if order is changed
    logical                     :: so

    if (obs% n_spot == 0) then

      if(associated(obs% body )) deallocate (obs% body )
      if(associated(obs% varno)) deallocate (obs% varno)
      if(associated(obs% olev )) deallocate (obs% olev )
      if(associated(obs% bger )) deallocate (obs% bger )
      if(associated(obs% bgeri)) deallocate (obs% bgeri)
      if(associated(obs% bgi  )) deallocate (obs% bgi  )
      if(associated(obs% t_int)) deallocate (obs% t_int)
      if(associated(obs% lev  )) deallocate (obs% lev  )
      if(associated(obs% par  )) deallocate (obs% par  )
      if(associated(obs% sink )) deallocate (obs% sink )
      if(associated(obs% spot )) deallocate (obs% spot )
      allocate     (obs% spot(0))

    else

      if (.not.associated(obs% spot)) allocate (obs% spot(0))
      call release_spot (obs% spot, obs% n_spot, keep)

      !----------------------------------
      ! sort observations according to id
      !----------------------------------
      so = .true.; if(present(sort)) so = sort
      if (so) then
        allocate (key (obs% n_spot))
        allocate (ipm (obs% n_spot))
        key = obs% spot% hd% id        ! set sort key
        ipm = index (key)              ! index array
        obs% spot = obs% spot (ipm)    ! reorder observations within box
        deallocate (key)
        deallocate (ipm)
      endif

      allocate (ix  (obs% n_spot))
      ix = obs% spot% o
      call reorder (obs% spot% o, ix, changed, obs% n_obs)
      call release_body (obs%  body, obs% n_obs, obs% spot% o, ix, changed)
      call release_int  (obs% varno, obs% n_obs, obs% spot% o, ix, changed)
      call release_real (obs%  olev, obs% n_obs, obs% spot% o, ix, changed)
      call release_real (obs%  bger, obs% n_obs, obs% spot% o, ix, changed)
      if ( size(obs% bger) == 0 ) deallocate (obs% bger)

      ix = obs% spot% i
      call reorder (obs% spot% i, ix, changed, obs% n_int)
      call release_int  (obs% t_int, obs% n_int, obs% spot% i, ix, changed)
      call release_real (obs%   lev, obs% n_int, obs% spot% i, ix, changed)
      call release_real (obs% bgeri, obs% n_int, obs% spot% i, ix, changed)
      call release_real (obs% bgi  , obs% n_int, obs% spot% i, ix, changed)
      if ( size(obs% bgeri) == 0 ) deallocate (obs% bgeri)
      if ( size(obs% bgi  ) == 0 ) deallocate (obs% bgi  )

      ix = obs% spot% p
      call reorder (obs% spot% p, ix, changed, obs% n_par)
      call release_int  (obs%   par, obs% n_par, obs% spot% p, ix, changed)

      ix = obs% spot% d
      call reorder (obs% spot% d, ix, changed, obs% n_sink)
      call release_sink (obs%  sink, obs% n_sink, obs% spot% d, ix, changed)

      deallocate (ix)

    endif

  contains
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    subroutine release_int (a, n, ixn, ixo, c)
    integer        ,pointer    :: a (:)   ! integer array to rearrange
    integer        ,intent(in) :: n       ! final size of i
    type (t_index) ,intent(in) :: ixn (:) ! new indices
    type (t_index) ,intent(in) :: ixo (:) ! old indices
    logical        ,intent(in) :: c       ! changed
      integer ,pointer :: t (:)
      integer          :: j, m
      m = size (ixn)
      if (associated (a)) then
        if (size(a) > n .or. c) then
          allocate (t (n))
          if (c) then
            do j = 1, m
              t (ixn(j)% i+1 : ixn(j)% i+ixn(j)%n) = &
              a (ixo(j)% i+1 : ixo(j)% i+ixo(j)%n)
            end do
          else
            t = a (:n)
          endif
          deallocate (a)
          a => t
        endif
      endif
    end subroutine release_int
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    subroutine release_real (a, n, ixn, ixo, c)
    real(wp)       ,pointer    :: a (:)   ! real array to rearrange
    integer        ,intent(in) :: n       ! final size of i
    type (t_index) ,intent(in) :: ixn (:) ! new indices
    type (t_index) ,intent(in) :: ixo (:) ! old indices
    logical        ,intent(in) :: c       ! changed
      real(wp) ,pointer :: t (:)
      integer           :: j, m
      m = size (ixn)
      if (associated (a)) then
        if (size(a) > n .or. c) then
          allocate (t (n))
          if (c) then
            do j = 1, m
              t (ixn(j)% i+1 : ixn(j)% i+ixn(j)%n) = &
              a (ixo(j)% i+1 : ixo(j)% i+ixo(j)%n)
            end do
          else
            t = a (:n)
          endif
          deallocate (a)
          a => t
        endif
      else
        allocate (a(0))
      endif
    end subroutine release_real
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    subroutine release_body (a, n, ixn, ixo, c)
    type (t_datum) ,pointer    :: a (:)   ! derived type array to rearrange
    integer        ,intent(in) :: n       ! final size of i
    type (t_index) ,intent(in) :: ixn (:) ! new indices
    type (t_index) ,intent(in) :: ixo (:) ! old indices
    logical        ,intent(in) :: c       ! changed
      type (t_datum) ,pointer :: t (:)
      integer                 :: j, m
      m = size (ixn)
      if (associated (a)) then
        if (size(a) > n .or. c) then
          allocate (t (n))
          if (c) then
            do j = 1, m
              t (ixn(j)% i+1 : ixn(j)% i+ixn(j)%n) = &
              a (ixo(j)% i+1 : ixo(j)% i+ixo(j)%n)
            end do
          else
            t = a (:n)
          endif
          deallocate (a)
          a => t
        endif
      else
        allocate (a(0))
      endif
    end subroutine release_body
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    subroutine release_sink (a, n, ixn, ixo, c)
    type (t_sink)  ,pointer    :: a (:)   ! derived type array to rearrange
    integer        ,intent(in) :: n       ! final size of i
    type (t_index) ,intent(in) :: ixn (:) ! new indices
    type (t_index) ,intent(in) :: ixo (:) ! old indices
    logical        ,intent(in) :: c       ! changed
      type (t_sink) ,pointer :: t (:)
      integer                :: j, m
      m = size (ixn)
      if (associated (a)) then
        if (size(a) > n .or. c) then
          allocate (t (n))
          if (c) then
            do j = 1, m
              t (ixn(j)% i+1 : ixn(j)% i+ixn(j)%n) = &
              a (ixo(j)% i+1 : ixo(j)% i+ixo(j)%n)
            end do
          else
            t = a (:n)
          endif
          deallocate (a)
          a => t
        endif
      else
        allocate (a(0))
      endif
    end subroutine release_sink
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    subroutine release_spot (spot, n_spot, keep)
    type (t_spot) ,pointer              :: spot (:) ! report meta data
    integer       ,intent(inout)        :: n_spot   ! no of valid reports
    logical       ,intent(in) ,optional :: keep (:) ! mask for reports to keep

      integer                :: i, j, n_old
      type (t_spot) ,pointer :: old (:)

      n_old = n_spot
      if (present(keep)) n_spot = count (keep(1:n_old))
      if (associated (spot)) then
        if (size(spot) > n_spot) then
          old => spot
          allocate (spot (n_spot))
          if (present (keep)) then
            j = 0
            do i = 1, n_old
              if (.not. keep(i)) cycle
              j = j + 1
              spot(j) = old(i)
            end do
          else
            spot = old (:n_spot)
          endif
          deallocate (old)
        endif
      endif
    end subroutine release_spot
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    subroutine reorder (new, old, changed, n)
    type (t_index) ,intent(out) :: new (:)
    type (t_index) ,intent(in)  :: old (:)
    logical        ,intent(out) :: changed
    integer        ,intent(out) :: n
      integer :: i,j,m
      changed = .false.
      m       = size (old)
      i       = 0
      do j = 1, m
        if (old(j)% i < 0) then
          new(j) = old(j)
          cycle
        endif
        changed = changed .or. old(j)% i /= i
        new(j)% i =         i
        new(j)% n = old(j)% n
        i         = old(j)% n + i
      end do
      n = i
    end subroutine reorder
  end subroutine release_mem_0
!------------------------------------------------------------------------------
  subroutine new_par (obs, n, spot)
  type (t_obs)  ,intent(inout) :: obs   ! observation data type
  integer       ,intent(in)    :: n     ! words required
  type (t_spot) ,intent(inout) :: spot  ! specific observation info
  !-------------------------------------
  ! reserve memory for 'obs% par'-buffer
  !-------------------------------------
    integer ,pointer   :: tmp (:)
    integer            :: nr
    integer            :: new_n
    !----------------------
    ! set indices in 'spot'
    !----------------------
    if (spot% p% i < 0) then
      spot% p% i = obs% n_par
    else
      if (spot% p% i+spot% p% n > obs% n_par) &
        call finish('new_par','inconsistent indices')
      if (spot% p% n < n) spot% p% i = obs% n_par
    endif
    spot% p% n = n
    nr = spot% p% i + spot% p% n - obs% n_par
    !---------------------------------------------------
    ! at first time allocate with standard buffer length
    !---------------------------------------------------
    if (.not. associated(obs% par)) allocate (obs% par (n_par))
    !----------------------------
    ! increase buffer if required
    !----------------------------
    if (nr <= 0) return
    if (obs% n_par + nr > size(obs% par)) then
      new_n = nint((obs% n_par + nr) * 1.5)
      tmp => obs% par
      allocate (obs% par (new_n))
      obs% par (:obs% n_par) = tmp(:obs% n_par)
      deallocate (tmp)
    endif
    !------------
    ! set counter
    !------------
    obs% n_par  = obs% n_par + nr
  end subroutine new_par
!------------------------------------------------------------------------------
  subroutine new_int (obs, spot, n)
  type (t_obs)      ,intent(inout) :: obs   ! observation data type
  type (t_spot)     ,intent(inout) :: spot  ! specific observation info
  integer           ,intent(in)    :: n     ! words required
  !-------------------------------------
  ! reserve memory for 'obs% phs'-buffer
  !-------------------------------------
    real(wp)     ,pointer :: tmp  (:) ! temporary for reallocation
    integer      ,pointer :: itmp (:) ! temporary for reallocation
!   type(t_sink) ,pointer :: snk  (:) ! temporary for reallocation
    integer               :: nr       ! number of missing entries
    integer               :: new_n    ! new number of entries
    !----------------------
    ! set indices in 'spot'
    !----------------------
    if (spot% i% i < 0) then
      spot% i% i = obs% n_int
    else
      if (spot% i% i+spot% i% n > obs% n_int) &
        call finish('new_int','inconsistent indices')
      if (spot% i% n < n) spot% i% i = obs% n_int
    endif
    spot% i% n = n
    nr = spot% i% i + spot% i% n - obs% n_int
    !---------------------------------------------------
    ! at first time allocate with standard buffer length
    !---------------------------------------------------
    if (.not. associated(obs% t_int)) then
      allocate (obs% t_int (n_int))
      obs% t_int = 0
      allocate (obs% lev   (n_int))
      obs% lev = - huge(1._wp)
    endif
    !----------------------------
    ! increase buffer if required
    !----------------------------
    if (nr <= 0) return
    if (obs% n_int + nr > size(obs% t_int)) then
      new_n = nint((obs% n_int + nr) * 1.5)
      itmp => obs% t_int
      tmp  => obs%   lev
      allocate (obs% t_int (new_n))
      allocate (obs%   lev (new_n))
      obs% t_int (:obs% n_int)    = itmp(:obs% n_int)
      obs%   lev (:obs% n_int)    =  tmp(:obs% n_int)
      obs% t_int ( obs% n_int+1:) = 0
      obs%   lev ( obs% n_int+1:) = 0._wp
      deallocate (itmp)
      deallocate ( tmp)
    endif
    !------------
    ! set counter
    !------------
    obs% n_int  = obs% n_int + nr
  end subroutine new_int
!------------------------------------------------------------------------------
  subroutine new_sink (obs, spot, n)
  type (t_obs)      ,intent(inout) :: obs   ! observation data type
  type (t_spot)     ,intent(inout) :: spot  ! specific observation info
  integer           ,intent(in)    :: n     ! words required
  !-------------------------------------------
  ! reserve memory for sink variable meta data
  !-------------------------------------------
    type(t_sink) ,pointer :: snk  (:) ! temporary for reallocation
    integer               :: nr       ! number of missing entries
    integer               :: new_n    ! new number of entries
    !----------------------
    ! set indices in 'spot'
    !----------------------
    if (spot% d% i < 0) then
      spot% d% i = obs% n_sink
    else
      if (spot% d% i+spot% d% n > obs% n_sink) &
        call finish('new_int(sink)','inconsistent indices')
      if (spot% d% n < n) spot% d% i = obs% n_sink
    endif
    spot% d% n = n
    nr = spot% d% i + spot% d% n - obs% n_sink
    !---------------------------------------------------
    ! at first time allocate with standard buffer length
    !---------------------------------------------------
    if (.not. associated(obs% sink)) then
      allocate (obs% sink (n_sink))
    endif
    !----------------------------
    ! increase buffer if required
    !----------------------------
    if (nr <= 0) return
    if (obs% n_sink + nr > size(obs% sink)) then
      new_n = nint((obs% n_sink + nr) * 1.5)
      snk  => obs% sink
      allocate (obs% sink (new_n))
      obs% sink  (:obs% n_sink)   = snk(:obs% n_sink)
      deallocate (snk)
    endif
    !------------
    ! set counter
    !------------
    obs% n_sink  = obs% n_sink + nr
  end subroutine new_sink
!------------------------------------------------------------------------------
  subroutine new_obs (obs, n, spot)
  type(t_obs)  ,intent(inout)           :: obs  ! observations derived type
  integer      ,intent(in)              :: n    ! no. observations to allocate
  type(t_spot) ,intent(inout) ,optional :: spot ! report derived type
  !---------------------------------------------------------------------
  ! Reserve memory for a new report.
  ! 'n' is the amount of memory required for the report.
  ! The pointer index 'spot% o' and the counter 'obs% n_obs' is updated.
  ! The following components are reallocated if memory is short:
  !   obs% obs
  !   obs% t_obs
  !   obs% body
  !   obs% olev
  ! Actually 1.5 times the total amount of memory required is allocated
  ! so that some spare memory is available for future requests.
  ! If 'spot' is missing exactly obs% n_obs + n words are allocated
  ! without updating 'obs% n_obs' for future requests.
  !---------------------------------------------------------------------
    real(wp)      ,pointer ::   tmp (:)
    integer       ,pointer :: t_tmp (:)
    type(t_datum) ,pointer :: b_tmp (:)
    integer                :: nr        ! no. additional observations required
    integer                :: new_n

    !----------------------
    ! set indices in 'spot'
    !----------------------
    if (present (spot)) then
      if (spot% o% i < 0) then
        spot% o% i = obs% n_obs
      else
        if (spot% o% i+spot% o% n > obs% n_obs) &
          call finish('new_obs','inconsistent indices')
        if (spot% o% n < n) spot% o% i = obs% n_obs
      endif
      spot% o% n = n
      nr    = spot% o% i + spot% o% n - obs% n_obs
      new_n = nint((obs% n_obs + nr) * 1.5)
    else
      nr    = n
      new_n = obs% n_obs + nr
    endif
    !---------------------------------------------------
    ! at first time allocate with standard buffer length
    !---------------------------------------------------
    if (.not. associated(obs%  body)) allocate (obs%  body (n_obs))
    if (.not. associated(obs% varno)) allocate (obs% varno (n_obs))
    if (.not. associated(obs%  olev)) then
      allocate (obs%  olev (n_obs))
      obs% olev = 0._wp
    endif
    !----------------------------
    ! increase buffer if required
    !----------------------------
    if (nr <= 0) return
    if (obs% n_obs + nr > size(obs% varno)) then
      tmp => obs% olev
      allocate (obs% olev (new_n))
      obs% olev (1:obs% n_obs ) = tmp(1:obs% n_obs)
      obs% olev (  obs% n_obs+1:) = 0._wp
      deallocate (tmp)
      b_tmp => obs% body
      allocate (obs% body  (new_n))
      obs% body  (1:obs% n_obs) = b_tmp(1:obs% n_obs)
      deallocate (b_tmp)
      t_tmp => obs% varno
      allocate (obs% varno (new_n))
      obs% varno (1:obs% n_obs) = t_tmp(1:obs% n_obs)
      deallocate (t_tmp)
    endif
    !------------
    ! set counter
    !------------
    if (present (spot)) obs% n_obs  = obs% n_obs + nr
  end subroutine new_obs
!------------------------------------------------------------------------------
  subroutine read_obs_nml
  !-----------------------------
  ! read namelist /OBSERVATIONS/
  !-----------------------------
    integer :: ierr, i, j, ftype
    logical :: exists
    logical :: obs_nml_read = .false.   ! namelist not yet read
    character(len=255) :: filename

    !------------------------------
    ! don't read the namelist twice
    !------------------------------
    if (obs_nml_read) then
      if(dace% lpio) then
        write(6,'()')
        write(6,'(a)') repeat ('-',79)
        write(6,'(a)') ' namelist /OBSERVATIONS/ skipped (already read)'
        write(6,'()')
      endif
      return
    endif
    obs_nml_read = .true.
    !------------------
    ! initialise tables
    !------------------
    call init_fdbk_tables
    !-------------
    ! set defaults
    !-------------
    obs_path    = obsinput ! path to read BUFR files
    bufr_files  = ''       ! BUFR   file names
    obs_files   = ''       ! NetCDF file names
    fdb_files   = ''       ! feedback file names
    bufr_verb   = 0        ! verbose flag for buffer decoding routines
    bufr_pause  = .false.  ! wait (for carriage return) after each record
    netcdf_verb = 0        ! verbosity level of NetCDF decoding routines
    derive_dbkz = .false.  ! derive suitable DBKZ if not present
    fix_no_obs  = .true.   ! set to .false. to revert no_obs-check
    int_vh      = .true.   ! interpolation: 1.vert.,2.hor.
    int_nn      = .false.  ! horizontal interpolation: nearest neighbour
    int_rad_hum = 0        ! humidity variable for vertical interp. for radiances
    int_rad_hum_pblend = (/300._wp,100._wp/) ! Blending levels for humidity interpolation for vertical interp. of radiances
    vint_lin_t  = .false.  ! linear vertical interpolation for temperature
    vint_lin_z  = .false.  ! linear vertical interpolation for geopotential
    vint_lin_uv = .false.  ! linear vertical interpolation for wind
    vint_lin_tov= .false.  ! linear vertical interpolation for RTTOV
#if (_RTTOV_VERSION >= 12)
    vint_rttov  = interp_rochon
    vint_rttov_b= interp_rochon
#endif
    monitor_ff  = .false.  ! monitor wind speed     (in addition to u, v)
    monitor_dd  = .false.  ! monitor wind direction (in addition to u, v)
    read_bufr   = .false.  ! read observations from BUFR   files
    read_NetCDF = .false.  ! read observations from NetCDF files
    read_cdfin  = .false.  ! read observations with COSMO routines
    dace_obs_op = ''       ! use dace obs. operators (MEC)
    obs_local   = .false.  ! keep observations on gridpoint PE (for COSMO MEC)
    fdbk_split(:,:) = 0    ! split feedback files
    !--------------
    ! read namelist
    !--------------
    if (dace% lpio) then
      call position_nml ('OBSERVATIONS',status=ierr)
      select case (ierr)
      case (POSITIONED)
#if defined(__ibm__)
        read (nnml ,nml=OBSERVATIONS, iostat=ierr)
        if (ierr/=0) call finish ('nml_run_flags',                  &
                                  'ERROR in namelist /OBSERVATIONS/')
#else
        read (nnml ,nml=OBSERVATIONS)
#endif
      end select
    endif

    !-----------------------------
    ! derive array from dace_obs_op
    !-----------------------------
    if (dace% lpio) then
      call split(dace_op,dace_obs_op,n_dace_op)
    end if

    !----------
    ! broadcast
    !----------
    call p_bcast (obs_path    ,dace% pio)
    call p_bcast (bufr_files  ,dace% pio)
    call p_bcast (obs_files   ,dace% pio)
    call p_bcast (fdb_files   ,dace% pio)
    call p_bcast (bufr_verb   ,dace% pio)
    call p_bcast (bufr_pause  ,dace% pio)
    call p_bcast (netcdf_verb ,dace% pio)
    call p_bcast (derive_dbkz ,dace% pio)
    call p_bcast (fix_no_obs  ,dace% pio)
    call p_bcast (int_vh      ,dace% pio)
    call p_bcast (int_nn      ,dace% pio)
    call p_bcast (int_rad_hum ,dace% pio)
    call p_bcast (int_rad_hum_pblend ,dace% pio)
    call p_bcast (vint_lin_t  ,dace% pio)
    call p_bcast (vint_lin_z  ,dace% pio)
    call p_bcast (vint_lin_uv ,dace% pio)
    call p_bcast (vint_lin_tov,dace% pio)
    call p_bcast (monitor_ff  ,dace% pio)
    call p_bcast (monitor_dd  ,dace% pio)
    call p_bcast (nwv_rad     ,dace% pio)
    call p_bcast (ndv_rad     ,dace% pio)
    call p_bcast (vint_rttov  ,dace% pio)
    call p_bcast (vint_rttov_b,dace% pio)
    call p_bcast (z2p_amv     ,dace% pio)
    call p_bcast (read_bufr   ,dace% pio)
    call p_bcast (read_NetCDF ,dace% pio)
    call p_bcast (read_cdfin  ,dace% pio)
    call p_bcast (dace_op     ,dace% pio)
    call p_bcast (n_dace_op   ,dace% pio)
    call p_bcast (obs_local   ,dace% pio)
    call p_bcast (qbit_conv   ,dace% pio)
    call p_bcast (pcc_conv    ,dace% pio)
    call p_bcast (corme_conv  ,dace% pio)
    call p_bcast (ptop_lapse  ,dace% pio)
    call p_bcast (pbot_lapse  ,dace% pio)
    call p_bcast (fdbk_split  ,dace% pio)
    call p_bcast (z2p_hlos    ,dace% pio)
    call p_bcast (spt_hd_debug,dace% pio)
    call p_bcast (spt_debug   ,dace% pio)
    call p_bcast (usd         ,dace% pio)
    do i=1,size(obs_files)
      if (obs_files(i) /=' ') then
        if (dace% lpio) then
           !-------------------------------
           ! Check existence of NetCDF file
           !-------------------------------
           filename = path_file (obs_path, obs_files(i))
           inquire (file=trim (filename), exist=exists)
           if (.not. exists) then
              call message ("read_obs_nml", &
                            "WARNING: file not found: " // trim (filename))
           end if
        end if
        call p_bcast (exists, dace% pio)
        ftype = FT_MISSING; if (exists) ftype = FT_NETCDF
        call add_source (obs_path, obs_files(i), ftype, complete=.true.)
      endif
    end do

    do i=1,size(fdb_files)
      if (fdb_files(i) /=' ') then
        if (dace% lpio) then
           !---------------------------------
           ! Check existence of feedback file
           !---------------------------------
           filename = path_file (obs_path, fdb_files(i))
           inquire (file=trim (filename), exist=exists)
           if (.not. exists) then
              call message ("read_obs_nml", &
                            "WARNING: file not found: " // trim (filename))
           end if
        end if
        call p_bcast (exists, dace% pio)
        ftype = FT_MISSING; if (exists) ftype = FT_FEEDBACK
        call add_source (obs_path, fdb_files(i), ftype, complete=.true.)
      endif
    end do

    do i=1,size(bufr_files)
      if (bufr_files(i) /=' ') then
        if (dace% lpio) then
           !-----------------------------
           ! Check existence of BUFR file
           !-----------------------------
           filename = path_file (obs_path, bufr_files(i))
           inquire (file=trim (filename), exist=exists)
           if (.not. exists) then
              call message ("read_obs_nml", &
                            "WARNING: file not found: " // trim (filename))
           end if
        end if
        call p_bcast (exists, dace% pio)
        ftype = FT_MISSING; if (exists) ftype = FT_BUFR
        call add_source (obs_path, bufr_files(i), ftype)
      endif
    end do
    if (all (bufr_files == '')) read_bufr   = .false.
    if (all (obs_files  == '')) read_NetCDF = .false.

    ldeb_spot = any(spt_hd_debug > 0 .or. spt_debug > 0)
    !---------
    ! printout
    !---------
    if (dace% lpio) then
      write(6,'()')
      write(6,'(a)') repeat ('-',79)
      write(6,'(a)') ' namelist /OBSERVATIONS/ read'
      write(6,'()')
      write(6,'(a,a)')    '  obs_path        = ', trim(obs_path)
!  character(len=64) :: bufr_files(63)= ''      ! names of BUFR files
!  character(len=64) :: obs_files (63)= ''      ! observation file names (NetCDF)
!  character(len=64) :: fdb_files (63)= ''      ! feedback (input) files
!     write(6,'(a,l1)')   '  obs_nml_read    = ', obs_nml_read
      write(6,'(a,l1)')   '  read_bufr       = ', read_bufr
      write(6,'(a,l1)')   '  read_NetCDF     = ', read_NetCDF
      write(6,'(a,l1)')   '  read_cdfin      = ', read_cdfin
      write(6,'(a,a)')    '  dace_obs_op     = ', trim(dace_obs_op)
      write(6,'(a,l1)')   '  obs_local       = ', obs_local
      write(6,'(a,i1)')   '  bufr_verb       = ', bufr_verb
      write(6,'(a,l1)')   '  bufr_pause      = ', bufr_pause
      write(6,'(a,i1)')   '  netcdf_verb     = ', netcdf_verb
      write(6,'(a,l1)')   '  derive_dbkz     = ', derive_dbkz
      write(6,'(a,l1)')   '  fix_no_obs      = ', fix_no_obs
      write(6,'(a,l1)')   '  int_vh          = ', int_vh
      write(6,'(a,l1)')   '  int_nn          = ', int_nn
      write(6,'(a,i0)')   '  int_rad_hum     = ', int_rad_hum
      write(6,'(a,2(1x,F5.0))')   '  int_rad_hum_pblend= ', int_rad_hum_pblend
      write(6,'(a,l1)')   '  vint_lin_t      = ', vint_lin_t
      write(6,'(a,l1)')   '  vint_lin_z      = ', vint_lin_z
      write(6,'(a,l1)')   '  vint_lin_uv     = ', vint_lin_uv
      write(6,'(a,l1)')   '  vint_lin_tov    = ', vint_lin_tov
      write(6,'(a,l1)')   '  monitor_ff      = ', monitor_ff
      write(6,'(a,l1)')   '  monitor_dd      = ', monitor_dd
      write(6,'(a,i0)')   '  nwv_rad         = ', nwv_rad
      write(6,'(a,i0)')   '  ndv_rad         = ', ndv_rad
      write(6,'(a,i0)')   '  vint_rttov      = ', vint_rttov
      write(6,'(a,i0)')   '  vint_rttov_b    = ', vint_rttov_b
      write(6,'(a,i0)')   '  z2p_amv         = ', z2p_amv
      write(6,'(a,i0)')   '  z2p_hlos        = ', z2p_hlos
      write(6,'(a,f4.0)') '  ptop_lapse      = ', ptop_lapse
      write(6,'(a,f4.0)') '  pbot_lapse      = ', pbot_lapse
      write(6,'(a,i0)')   '  qbit_conv       = ', qbit_conv
      write(6,'(a,i0)')   '  pcc_conv        = ', pcc_conv
      write(6,'(a,i0)')   '  corme_conv      = ', corme_conv
      write(6,'(A,20(1x,A6))')   '                   ',(trim(name_value(obstype_t,i)),i=1,obstype_t%n)
      do j = 1, mstep
        write(6,'(A,20(1x,I6))') '  fdbk_split('//trim(pref_step(j))//') =',(fdbk_split(j,i),i=1,obstype_t%n)
      end do
    endif

    call read_wigos_nml ()
  end subroutine read_obs_nml
!==============================================================================
  subroutine set_vqc_insitu (spot, obs)
  type (t_spot) ,intent(in)    :: spot ! observation meta information
  type (t_obs)  ,intent(inout) :: obs  ! observation data type
  !-----------------------------------------------------------------------
  ! Set the variational quality control bound (multiple of the
  ! observational error) for in situ observations.  For a given
  ! observation (TEMP,SYNOP) described by 'spot' the corresponding
  ! elements of the component field 'obs% s_vqc' are set according to the
  ! specifications in namelist /RULES/.
  !----------------------------------------------------------------------
    type (t_set) :: t
    type (t_set) :: gp
    type (t_set) :: uv
    type (t_set) :: q
    type (t_set) :: p
    real(sp)     :: sgm_vq
    integer      :: frm_vq      ! var. quality control formulation
    integer      :: i, l
    !------------------------------------------------------------------
    ! make sure the array to hold the a priory probability is allocated
    !------------------------------------------------------------------
    if (.not. associated (obs% s_vqc)) then
      allocate (obs% s_vqc (obs% n_obs))
      obs% s_vqc = svqc
    endif
    if (.not. associated (obs% f_vqc)) then
      allocate (obs% f_vqc (obs% n_obs))
      obs% f_vqc = vqc_form
    endif
    !---------------------------------------
    ! get the settings from namelist /RULES/
    !---------------------------------------
    call get_rule (type     = spot% hd% modtype,  &! <- module type
                   obstype  = spot% hd% obstype,  &! <- observation type
                   codetype = spot% hd% codetype, &! <- code type
                   bf_type  = iud,                &! <- no BUFR     type
                   bf_subt  = iud,                &! <- no BUFR  subtype
                   db_kz    = iud,                &! <- no Datenbankkennzahl
                   stname   = '',                 &! <- no Station name
                   lat      = spot% col% c% dlat, &! <- latitude
                   lon      = spot% col% c% dlon, &! <- longitude
                   t        = t,                  &! -> info on temperature
                   gp       = gp,                 &! -> info on geopotential
                   uv       = uv,                 &! -> info on wind
                   q        = q,                  &! -> info on humidity
                   p        = p)                   ! -> info on pressure
    !---------------------------------------------
    ! set the appropriate elements of 'obs% s_vqc'
    !---------------------------------------------
    do i=1,spot% o% n
      l = spot% o% i + i
      sgm_vq = rud
      frm_vq = iud
      select case (obs% varno(l))
      case (VN_Z, VN_HEIGHT)              ! geopotential, height
        sgm_vq = gp% sgm_vq
        frm_vq = gp% frm_vq
      case (VN_PS, VN_P)                  ! surface pressure
        sgm_vq = p% sgm_vq
        frm_vq = p% frm_vq
      case (VN_VT)                        ! (virtual) temperature
        sgm_vq = t% sgm_vq
        frm_vq = t% frm_vq
      case (VN_RH, VN_RH2M)               ! relative humidity
        sgm_vq = q% sgm_vq
        frm_vq = q% frm_vq
      case (VN_U, VN_V, VN_U10M, VN_V10M) ! wind component
        sgm_vq = uv% sgm_vq
        frm_vq = uv% frm_vq
      case (VN_FF, VN_HLOS)               ! wind speed, line of sight wind
        sgm_vq = uv% sgm_vq
        frm_vq = uv% frm_vq
      case (VN_T, VN_T2M)                 ! temperature
        sgm_vq = t% sgm_vq
        frm_vq = t% frm_vq
      case (VN_Q)                         ! specific humidity
        sgm_vq = q% sgm_vq
        frm_vq = q% frm_vq
      case (VN_N, VN_NH, VN_N_L, VN_N_M, VN_N_H, VN_CEIL, VN_W, &
            VN_TD2M, VN_GUST, VN_RR, VN_TMIN, VN_TMAX, VN_PTEND,&
            VN_DD, VN_VV, VN_RAD_GL, VN_RAD_DF, VN_RAD_LW,      &
            VN_SDEPTH, VN_TSEA, VN_WW, VN_GCLG, VN_ICLG         )
        !-------------------------
        ! passive, not assimilated
        !-------------------------
        sgm_vq = 100._sp
        frm_vq = 0
      case default
        call finish('set_vqc_bound','unknown observation: '//char3(obs% varno(l)))
      end select
      if (sgm_vq == rud) then
        obs% s_vqc (l) = svqc
      else
        obs% s_vqc (l) = sgm_vq
      endif
      if (frm_vq == iud) then
        obs% f_vqc (l) = vqc_form
      else
        obs% f_vqc (l) = frm_vq
      endif
    end do

  end subroutine set_vqc_insitu
!==============================================================================
  subroutine set_int_insitu (spt, obs)
  type (t_spot) ,intent(inout) :: spt
  type (t_obs)  ,intent(inout) :: obs
  !----------------------------------------------
  ! set interpolation space for in-situ-variables
  ! (TEMP, PILOT, AIREP,..)
  !----------------------------------------------
    integer  :: oi, ii, i, ni
    real(wp) :: olev  ! level of observation
    real(wp) :: plev  ! pressure level
    logical  :: tq, uv
    if (spt% o% i <  0) return
    if (spt% o% n == 0) then
      spt% i% i = -1
      spt% i% n =  0
    else
      oi = spt% o% i
      !--------------------------
      ! count params in int-space
      !--------------------------
      ni   = 0
      olev = -huge(olev)
      do i = 1, spt% o% n
        if (olev /= obs% olev (oi+i)) then
          olev = obs% olev (oi+i)
          tq   = .false.
          uv   = .false.
        endif
        select case (obs% varno(oi+i))
        case (VN_RH, VN_VT, VN_Q, VN_T)
          if (tq) cycle
          tq = .true.
          ni = ni + 2
        case (VN_U, VN_V, VN_FF, VN_DD, VN_HLOS)
          if (uv) cycle
          uv = .true.
          ni = ni + 2
        case (VN_Z)
          ni = ni + 1
        case default
        end select
      end do
      !------------------------
      ! set interpolation space
      !------------------------
      call new_int (obs, spt, ni)
      ii   = spt% i% i
      olev = -huge(olev)
      ni  = 0
      do i = 1, spt% o% n
        if (olev /= obs% olev (oi+i)) then
          olev = obs% olev (oi+i)
          tq   = .false.
          uv   = .false.
          select case (obs% body(oi+i)% lev_typ)
          case (VN_P)
            plev = obs% olev (oi+i)
          case default
            plev = obs% body(oi+i)% plev
          end select
        endif
        select case (obs% varno(oi+i))
        case (VN_RH, VN_VT, VN_Q, VN_T)
          !--------------------------------
          ! for T, Q we need both tv and rh
          !--------------------------------
          if (plev <= 0.) then
            write(0,*)'obstype, i varno lev_typ olev plev =',             &
            spt% hd% obstype, i, obs% varno(oi+i),                        &
            obs% body(oi+i)% lev_typ,obs% olev (oi+i),obs% body(oi+i)% plev
            call finish ('set_int_insitu','plev <= 0')
          endif
          if (tq) cycle
          tq = .true.
          obs% t_int(ii+ni+1) = OBS_TV
          obs% t_int(ii+ni+2) = OBS_RH
          obs% lev  (ii+ni+1) = log (plev)
          obs% lev  (ii+ni+2) = obs% lev  (ii+ni+1)
          ni = ni + 2
        case (VN_U, VN_V, VN_FF, VN_DD, VN_HLOS)
          !---------------------
          ! we need both u and v
          !---------------------
          if (plev <= 0.) then
            write(0,*)'Wind obstype, i varno lev_typ olev plev =',             &
            spt% hd% obstype, i, obs% varno(oi+i),                        &
            obs% body(oi+i)% lev_typ,obs% olev (oi+i),obs% body(oi+i)% plev
            call finish ('set_int_insitu','plev <= 0')
          endif
          if (uv) cycle
          uv = .true.
          obs% t_int(ii+ni+1) = OBS_U
          obs% t_int(ii+ni+2) = OBS_V
          obs% lev  (ii+ni+1) = log (plev)
          obs% lev  (ii+ni+2) = obs% lev  (ii+ni+1)
          ni = ni + 2
        case (VN_Z)
          !--------------------------------------
          ! for Z just interpolate the components
          !--------------------------------------
          if (plev <= 0.) then
            write(0,*)'obstype, i varno lev_typ olev plev =',             &
            spt% hd% obstype, i, obs% varno(oi+i),                        &
            obs% body(oi+i)% lev_typ,obs% olev (oi+i),obs% body(oi+i)% plev
            call finish ('set_int_insitu','plev <= 0')
          endif
          ni = ni + 1
          obs% t_int(ii+ni) = varno_oq (obs% varno(oi+i))
          obs% lev  (ii+ni) = log (plev)
        case (VN_U10M,VN_V10M,VN_T2M,VN_RH2M,VN_NH,VN_N_L,VN_N_M,VN_N_H, &
              VN_HEIGHT,VN_W,VN_PS,VN_VGUST,VN_TURB)
          !---------------
          !not implemented
          !---------------
        case default
          call finish ('set_int_insitu',                               &
                       'varno = '//name_value (varno, obs% varno(oi+i)))
        end select
      end do
    endif
  end subroutine set_int_insitu
!==============================================================================
  subroutine shrink_reports (obs, state)
  type(t_obs) ,intent(inout) :: obs     ! observation data
  integer     ,intent(in)    :: state   ! status to keep
    integer :: is
    logical :: change
    do is = 1, obs% n_spot
      call shrink_report (obs% spot(is) , obs, state, change)
    end do
  end subroutine shrink_reports
!------------------------------------------------------------------------------
  subroutine shrink_report (spt, obs, state, change, mask, lindx, nl)
  !---------------------------------------------------------------
  ! keep only observations with status >= state
  ! only if change==.true. the optional return arguments are valid
  !                        and pointers will be allocated
  !---------------------------------------------------------------
  type(t_spot) ,intent(inout)          :: spt      ! report data
  type(t_obs)  ,intent(inout)          :: obs      ! observation data
  integer      ,intent(in)             :: state    ! status to keep
  logical      ,intent(out)            :: change   ! report was changed
  logical      ,pointer     ,optional  :: mask (:) ! true for observations kept
  integer      ,pointer     ,optional  :: lindx(:) ! index for level selection
  integer                   ,optional  :: nl       ! number of levels kept

    integer                    :: i, j, k, l, m ! indices
    type(t_datum) ,pointer     :: bd (:)     ! pointer to report body data
    integer                    :: cnt        ! counter
    real(wp)                   :: llev, mlev
    integer       ,allocatable :: li (:)

    change = .false.
    if (spt% use% state <  state) return
    if (spt% o%   n     <= 0            ) return
    bd => obs% body(spt% o% i+1: spt% o% i+spt% o% n)
    cnt = count (bd% use% state >= state)
    if (cnt /= spt% o% n) then
      change = .true.
      if (present(mask)) then
        allocate (mask (spt% o% n))
        mask = .false.
      endif
      if (present(lindx)) then
        allocate (li (spt% col% nlev))
      endif
      llev = huge(1._wp)
      mlev = huge(1._wp)
      l    = 0
      m    = 0
      j    = spt% o% i
      do i = 1, spt% o% n
        k = i + spt% o% i
        if (obs% olev(k) /= llev) then
          l    = l + 1
          llev = obs% olev(k)
        endif
        if (bd(i)% use% state >= state) then
          if (obs% olev(k) /= mlev) then
            m     = m + 1
            mlev  = obs% olev(k)
            if (present (lindx)) li(m) = l
          endif
          if (present(mask)) mask(i) = .true.
          j = j + 1
                                      obs% body  (j) = obs% body  (k)
                                      obs% varno (j) = obs% varno (k)
!                                     obs% qcf   (j) = obs% qcf   (k)
          if(associated(obs% bger ))  obs% bger  (j) = obs% bger  (k)
          if(associated(obs% s_vqc))  obs% s_vqc (j) = obs% s_vqc (k)
          if(associated(obs% f_vqc))  obs% f_vqc (j) = obs% f_vqc (k)
                                      obs% olev  (j) = obs% olev  (k)
        end if
      end do
      spt% o% n = cnt
      if (present (nl)   ) nl = m
      if (present (lindx)) then
        if (l /= spt% col% nlev) call finish ('shrink_report','l /= spt% col% nlev')
        allocate (lindx (m))
        lindx = li (1:m)
      endif
    end if
  end subroutine shrink_report
!==============================================================================
  subroutine set_inbx (obs)
  type(t_obs)  ,intent(inout)          :: obs(:)   ! observation data
  !-----------------------------
  ! set component inbx in t_spot
  !-----------------------------
    integer :: ib, is
    do ib = 1, size (obs)
      do is = 1, obs(ib)% n_spot
        obs(ib)% spot(is)% inbx% pe   = dace% pe
        obs(ib)% spot(is)% inbx% box  = ib
        obs(ib)% spot(is)% inbx% spot = is
      end do
    end do
  end subroutine set_inbx
!==============================================================================
  subroutine forget_intp (obs)
  type (t_obs) ,intent(inout) :: obs (:)
  !------------------------------------------------------------------------
  ! remove any information on the interpolation space in derived type t_obs
  ! to be called before t_obs is associated with a new background model
  !------------------------------------------------------------------------
    integer        :: ib, is
    type (t_index) :: x ! empty

    do ib = 1, size(obs)

      obs(ib)% n_lev   = 0
      obs(ib)% n_int   = 0
      obs(ib)% n_spt   = 0
      obs(ib)% mc% n   = 0
      if (associated (obs(ib)% levs   )) deallocate (obs(ib)% levs   )
      if (associated (obs(ib)% t_int  )) deallocate (obs(ib)% t_int  )
      if (associated (obs(ib)% lev    )) deallocate (obs(ib)% lev    )
      if (associated (obs(ib)% bgeri  )) deallocate (obs(ib)% bgeri  )
      if (associated (obs(ib)% bgi    )) deallocate (obs(ib)% bgi    )

      if (associated (obs(ib)% mc% idx)) deallocate (obs(ib)% mc% idx)

      if (associated (obs(ib)% mc% c  )) deallocate (obs(ib)% mc% c  )
      do is = 1, obs(ib)% n_spot

        obs(ib)% spot(is)% n_spt = 1
        obs(ib)% spot(is)% mke   = 0
        obs(ib)% spot(is)% i     = x
        obs(ib)% spot(is)% l     = x

        if (associated (obs(ib)% spot(is)% imcol)) &
            deallocate (obs(ib)% spot(is)% imcol)

      end do
    end do

  end subroutine forget_intp
!==============================================================================
  subroutine set_sink (spot, obs, isink, bg, err,      &
                       iobs, min, max, lbnd, ubnd, mode)
  type(t_spot) ,intent(in)           :: spot  ! report / FOV
  type(t_obs)  ,intent(inout)        :: obs   ! observation meta data
  integer      ,intent(in)           :: isink ! sink variable index
  real(wp)     ,intent(in)           :: bg    ! background value to set
  real(wp)     ,intent(in)           :: err   ! background error to set
  integer      ,intent(in) ,optional :: iobs  ! reference to observ. in report
  real(wp)     ,intent(in) ,optional :: min   ! lower bound
  real(wp)     ,intent(in) ,optional :: max   ! upper bound
  real(wp)     ,intent(in) ,optional :: lbnd  ! nonlinear range
  real(wp)     ,intent(in) ,optional :: ubnd  ! nonlinear range
  integer      ,intent(in) ,optional :: mode  ! 0:no, 1:quadratic, 2:log trnsf.
  !--------------------------------------------
  ! define background values for sink variables
  !--------------------------------------------
    integer  :: i ! index
    integer  :: m ! mode, set from 1 to 2 if bg == min or bg == max
    real(wp) :: x ! background, constrained to min/max value

    if (isink > spot% d% n) call finish('set_sink','isink > d% n')
    !-------------------
    ! check bounds, mode
    !-------------------
    i = spot% d% i + isink
    m = obs% sink(i)% mode; if (present(mode)) m = mode
    x = bg
    if (present(min)) then
      if (x <= min) then
        x = min
        if (m == 1) m = 2
      endif
    endif
    if (present(max)) then
      if (x >= max) then
        x = max
        if (m == 1) m = 2
      endif
    endif
!   if (err == 0._wp) m = 0
    obs%                     sink(i)% bg   = x
    obs%                     sink(i)% err  = err
    if (present (iobs)) obs% sink(i)% iobs = iobs
    if (present (min )) obs% sink(i)% min  = min
    if (present (max )) obs% sink(i)% max  = max
    if (present (lbnd)) obs% sink(i)% lbnd = lbnd
    if (present (ubnd)) obs% sink(i)% ubnd = ubnd
    if (present (mode)) obs% sink(i)% mode = m
    call set_v_sink (obs% sink(i))
  end subroutine set_sink
!==============================================================================
  subroutine add_source (path, file, filetype, obstype, codetype, &
                         complete, entries, nobs, file_format     )
  !------------------------
  ! Add source file to list
  !------------------------
  character(len=*)  ,intent(in) :: path       ! file path name
  character(len=*)  ,intent(in) :: file       ! file name
  integer           ,intent(in) :: filetype   ! 1:BUFR, 2:NetCDF, ...
  integer ,optional ,intent(in) :: obstype    ! use obst.specific input
  integer ,optional ,intent(in) :: codetype   ! assume CMA codetype
  logical ,optional ,intent(in) :: complete   ! read complete file on pe
  integer ,optional ,intent(in) :: entries    ! number of entries in file
  integer ,optional ,intent(in) :: nobs       ! number of observations
  integer ,optional ,intent(in) :: file_format! format (old/intermediate/new)
    integer :: n
    n_source              = n_source              + 1
    n_filetype (filetype) = n_filetype (filetype) + 1
    n = n_source
    if (n > m_source) call finish ('add_source','n_source > m_source')
    source(n)           = source_0
    source(n)% path     = path
    source(n)% file     = file
    source(n)% filetype = filetype
    if (present (obstype))      source(n)% obstype  = obstype
    if (present (codetype))     source(n)% codetype = codetype
    if (present (complete))     source(n)% complete = complete
    if (present (entries))      source(n)% entries  = entries
    if (present (nobs))         source(n)% nobs     = nobs
    if (present (file_format))  source(n)% format   = file_format
    call p_bcast        (n_source,   dace% pio)
    call p_bcast        (n_filetype, dace% pio)
    call p_bcast_source (source,     dace% pio)
  end subroutine add_source
!------------------------------------------------------------------------------
  subroutine empty_source ()
  !-----------------------
  ! Clear source file list
  !-----------------------
    n_source   = 0
    n_filetype = 0
  end subroutine empty_source
!==============================================================================
  subroutine match_ilev(tl, n_lev, instr, levs, si, oi, ind, iatms)
    type(t_ilev)  ,intent(in)                    :: tl(:)
    integer       ,intent(out)                   :: n_lev
    integer       ,intent(in)  ,optional ,target :: instr
    integer       ,intent(in)  ,optional ,target :: levs(:)
    type(t_spot)  ,intent(in)  ,optional ,target :: si
    type(t_obs)   ,intent(in)  ,optional         :: oi
    integer       ,intent(out) ,optional         :: ind(:)
    integer       ,intent(out) ,optional         :: iatms ! 1 for ATMS temp chans, 2 for ATMS hum. chans, 0 otherwise

    integer ,pointer    :: ii    => null()
    integer ,pointer    :: l(:) => null()
    integer             :: i, j, m, is, ie, nl, ni, mx_lev, mn_lev, id

    if (present(instr)) then
      ii => instr
      id =  rttov_instr(ii) ! Assume OT_RAD
    elseif (present(si)) then
      ii => si% sttyp
      id =  ii
      if (si% hd% obstype == OT_RAD) id = rttov_instr(id, int(si% hd% satid))
    else
      call finish('match_ilev','Either instr or si parameter required.')
    end if

    if (present(levs)) then
      l => levs
      is =  1
      m = size(l)
    elseif (present(si) .and. present(oi)) then
      if (si% hd% obstype == OT_RAD) then
        ! Consider only instrument ii (usually the grid instrument)
        is = 0
        ie = -1
        do i = si%o%i+1, si%o%i + si%o%n
          if (ii == oi%body(i)%lev_sig) then
            if (is == 0) is = i
            ie = i
          end if
        end do
      else
        is = si%o%i + 1
        ie = si%o%i + si%o%n
      end if
      if (ie >= is) then
        m = ie - is + 1
        allocate(l(m))
        l = int(oi% olev(is:ie))
      else
        m = 0
        l => null()
      end if
    else
      call finish('match_ilev','Either levs or si and oi parameters required.')
    end if

    nl    = size(tl)
    n_lev = m
    if (present(iatms)) then
      mx_lev = -1
      mn_lev = huge(1)
    end if
    if (present(ind)) then
      ni = size(ind)
      if (ni < n_lev) call error_loc('1')
      ind(1:m) = (/ (j, j=1,m) /)
    end if

    if (any(tl(:)%instr == id)) then
      n_lev = 0
      do i = 1, nl
        if (tl(i)%instr == id) then
          do j = 1, m
            if (l(j) == tl(i)%value) then
              n_lev = n_lev + 1
              if (present(ind)) then
                if (ni < n_lev) call error_loc('2')
                ind(n_lev) = j
              end if
              if (present(iatms)) then
                mn_lev = min(mn_lev, l(j))
                mx_lev = max(mx_lev, l(j))
              end if
              exit
            end if
          end do
        end if
      end do
    elseif (any(tl(:)%instr < 0)) then
      n_lev = 0
      do i = 1, nl
        do j = 1, m
          if (l(j) == tl(i)%value) then
            n_lev = n_lev + 1
            if (present(ind)) then
              if (ni < n_lev) call error_loc('3')
              ind(n_lev) = j
            end if
            if (present(iatms)) then
              mn_lev = min(mn_lev, l(j))
              mx_lev = max(mx_lev, l(j))
            end if
            exit
          end if
        end do
      end do
    elseif (present(iatms)) then
      if (associated(l)) then
        mn_lev = minval(l(:))
        mx_lev = maxval(l(:))
      else
        mn_lev = 0
        mx_lev = 0
      end if
    end if

    if (present(iatms)) then
      iatms = 0
      if (id == 19) then
        if (mx_lev <= 16) then
          iatms = 1
        elseif (mn_lev >= 15) then
          iatms = 2
        end if
      end if
    end if

    if (present(ind)) ind(1:n_lev) = ind(1:n_lev) + is - 1
    if (.not.present(levs) .and. present(si) .and. present(oi)) then
      if (associated(l)) deallocate(l)
    end if

  contains

    subroutine error_loc(hint)
      character, intent(in) :: hint
      character(len=1000)   :: msg
      msg = ''
      write(msg, '(" id=",I3," ni=",I9," n_lev=",I9," present args:",6(1x,L1))') &
           id,ni,n_lev, present(instr), present(levs), present(si), present(oi), &
           present(ind), present(iatms)
      msg = 'index array too small ('//trim(hint)//'):'//trim(msg)
      call finish('match_ilev',trim(msg))
    end subroutine error_loc

  end subroutine match_ilev
!------------------------------------------------------------------------------
  elemental function stat_hash (statid)
  !---------------------------------------------------
  ! hash function for character(len=10) -> integer(i8)
  ! assumes printable uppercase ASCII characters
  !---------------------------------------------------
  character(len=10) ,intent(in) :: statid
  integer (i8)                  :: stat_hash

    integer(i1) :: istat1 (10)
    integer(i8) :: t
    integer     :: i

    istat1 = transfer (statid, istat1)
    istat1 = ibclr(istat1, 7)

    stat_hash = 0_i8
    do i = 1, size(istat1)
      select case (istat1(i)/32)
      case (0)
        istat1(i) = 0_i1                        ! non printable characters
      case (1,2)
        istat1(i) = istat1(i) - 32_i1           ! uppercase characters, numbers
      case (3)
        istat1(i) = istat1(i) - 64_i1           ! lowercase -> uppercase
      end select
      t = istat1(i)
      stat_hash = stat_hash + ishft(t,(i-1)*6)  ! 6 bits/character left
    end do

  end function stat_hash
!==============================================================================

  function ldeb(spot)
    logical :: ldeb
    type(t_spot), intent(in) :: spot
    integer                  :: i
    ldeb = .false.
    if (.not.ldeb_spot) return
    do i = 1,nsd
      if (spot%hd%id == spt_hd_debug(i)) then
        if (i /= i_shd) then ! update dpref
          write(dpref,'("debug_spot ",I12)') spot%hd%id
          i_shd = i
        end if
        ldeb = .true.
        exit
      end if
      if (spot%id == spt_debug(i)) then
        if (i /= i_sd) then ! update dpref
          write(dpref,'("debug_spot ",I12)') spot%id
          i_sd = i
        end if
        ldeb = .true.
        exit
      end if
      if (ldeb_all) then
        write(dpref,'("debug_spot ",I12)') spot%hd%id
        ldeb = .true.
        exit
      end if
    end do
  end function ldeb


  subroutine debug_spot_obs(obs, flags, hint, msg)
    type(t_obs),      intent(in)           :: obs
    integer,          intent(in), optional :: flags(:)
    character(len=*), intent(in), optional :: hint
    character(len=*), intent(in), optional :: msg
    integer :: is
    if (.not.ldeb_spot) return
    do is = 1, obs%n_spot
      call debug_spot(obs% spot(is), obs=obs, flags=flags, hint=hint, msg=msg)
    end do
  end subroutine debug_spot_obs


  subroutine debug_spot_sp(sp, obs, flags, hint, msg)
    type(t_spot),     intent(in)                   :: sp
    type(t_obs),      intent(in), optional, target :: obs
    integer,          intent(in), optional         :: flags(:)
    character(len=*), intent(in), optional         :: hint
    character(len=*), intent(in), optional         :: msg
    type(t_obs), pointer :: obs_
    character(len=300)   :: prefix  =  ''
    integer              :: i

    if (.not.ldeb_spot) return
    if (.not.ldeb(sp)) return
    obs_ => null()
    if (present(obs)) obs_ => obs

    if (present(hint)) then
      prefix = dpref//trim(hint)
    else
      prefix = dpref
    end if
    write(usd,'(1x,A,1x,I3,2(1x,F8.3),1x,I3,1x,I3)') trim(prefix),&
         sp%hd%source, sp%col%c%dlat, sp%col%c%dlon, sp%pcc,sp%sttyp
    write(usd,*) trim(prefix)//' report use',sp%use%state,&
         sp%use%check,sp%use%flags
    if (present(msg)) then
      write(usd,*) trim(prefix)//' '//trim(msg)
    else
      if (associated(obs_)) then
        if (associated(obs_%body)) then
          if (present(flags)) then
            do i = 1, sp%o%n
              write(usd,'(1x,A,1x,I5,1x,I9,2(1x,I2),1x,I9,2(1x,I4))') &
                   trim(prefix),i,flags(i),                           &
                   obs_%body(sp%o%i+i)%use%state,                     &
                   obs_%body(sp%o%i+i)%use%check,                     &
                   obs_%body(sp%o%i+i)%use%flags,                     &
                   obs_%body(sp%o%i+i)%lev_sig,                       &
                   nint(obs_%olev(sp%o%i+i))
            end do
          else
            do i = 1, sp%o%n
              write(usd,'(1x,A,1x,I5,10x,2(1x,I2),1x,I9,2(1x,I4))') &
                   trim(prefix),i,                                  &
                   obs_%body(sp%o%i+i)%use%state,                   &
                   obs_%body(sp%o%i+i)%use%check,                   &
                   obs_%body(sp%o%i+i)%use%flags,                   &
                   obs_%body(sp%o%i+i)%lev_sig,                     &
                   nint(obs_%olev(sp%o%i+i))
            end do
          end if
        end if
      end if
    end if
  end subroutine debug_spot_sp

!------------------------------------------------------------------------------
! Broadcast source file list
!------------------------------------------------------------------------------
#define VECTOR
#define DERIVED type(t_source),dimension(:)
#define p_bcast_DERIVED p_bcast_source
#undef  MPI_TYPE
#include "p_bcast.incf"
!==============================================================================
#undef VECTOR
#undef DERIVED
#undef p_bcast_DERIVED
#define DERIVED type(t_bufr_inv)
#define p_bcast_DERIVED p_bcast_bufr_inv
#undef MPI_TYPE
#include "p_bcast.incf"
!------------------------------------------------------------------------------
#define VECTOR
#undef  DERIVED
#define DERIVED type(t_bufr_inv),DIMENSION(:)
#undef  p_bcast_DERIVED
#define p_bcast_DERIVED p_bcast_bufr_inv_1d
#undef MPI_TYPE
#include "p_bcast.incf"
!==============================================================================
end module mo_t_obs
