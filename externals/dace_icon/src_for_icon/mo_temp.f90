!
!+ Profile (TEMP, PILOT) observation operator specific routines
!
MODULE mo_temp
!
! Description:
!   Profile (TEMP, PILOT) observation operator specific routines.
!
! Current Maintainer: DWD, Harald Anlauf, Alexander Cress
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
!  Changes for verification mode
! V1_5         2009/05/25 Andreas Rhodin
!  modified TEMP/PILOT level selection (for monitoring)
! V1_6         2009/06/10 Andreas Rhodin
!  set u,v quality flags consistently with ff,dd
! V1_7         2009/08/24 Andreas Rhodin
!  fix humidity significant levels
! V1_8         2009/12/09 Andreas Rhodin
!  Changes for GME30L60 (model top at 5hPa)
!  Changes for COSMO    (pass pressure coordinates to interpolation routine)
! V1_9         2010/04/20 Andreas Rhodin
!  TSK_SHRINK in subroutines process: pass parameter 'state' to shrink_report
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_11        2010/06/16 Harald Anlauf
!  read_temp_netcdf: reduce verbosity for netcdf_verb==0
! V1_13        2011/11/01 Harald Anlauf
!  revised merge; optimisations for NEC SX-9
! V1_15        2011/12/06 Andreas Rhodin
!  replace unused status flag value STAT_USED by STAT_OBS_ONLY
! V1_18        2012/01/06 Harald Anlauf
!  Revised level selection and correction message processing (bugfix)
! V1_19        2012-04-16 Andreas Rhodin
!  subroutine correct_report: declare dummy variable 'btemp' as a pointer
! V1_20        2012-06-18 Harald Anlauf
!  properly initialize len_level(0),len_report(0)
! V1_22        2013-02-13 Harald Anlauf
!  bugfixes for merge of twin levels
!  changes for COSMO multilevel data and wind profilers
! V1_26        2013/06/27 Andreas Rhodin
!  detailed diagnostic in case of 'trh_tvgh' failure
! V1_27        2013-11-08 Harald Anlauf
!  Dismiss TEMPs with lat=0,lon=0
! V1_28        2014/02/26 Alexander Cress
!  split drifting temps
! V1_31        2014-08-21 Andreas Rhodin
!  changes for ECMWF TEMP (BUFR2)NetCDF input
!  revise level significance handling (for new BUFR templates)
! V1_35        2014-11-07 Alexander Cress
!  changes for TEMP BUFR reports
! V1_40        2015-02-27 Harald Anlauf
!  Prepare for Australian wind profiler (A.Cress)
! V1_42        2015-06-08 Andreas Rhodin
!  implement temporal interpolation for COSMO MEC
! V1_44        2015-09-30 Harald Anlauf
!  merge_report: handle cases where RS type is coded only in part B
!                increase tolerance to delta(lon,lat)=0.1 deg
!  split_report: bug fix (A.Rhodin)
! V1_46        2016-02-05 Andreas Rhodin
!  base decisions on new flag 'vct', not 'ivctype'
! V1_47        2016-06-06 Harald Anlauf
!  Recognize Australian wind profilers 95xxx
! V1_48        2016-10-06 Andreas Rhodin
!  passive monitoring of FF,DD; ignore NH,N_L,RH2M,T2M,U10M,V10M,HEIGHT
! V1_50        2017-01-09 Harald Anlauf
!  high-resolution BUFR TEMPs; additional VAD wind profilers (A.Cress)
! V1_51        2017-02-24 Harald Anlauf
!  check_store_temp: check for bogus/multiple surface levels
!  BUFR TEMPs: handle missing drift information at top of ascent
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Andreas Rhodin  MPIfM  2000  original code
! Andreas Rhodin  DWD    2001  changes for PSAS
! Oliver Schmid   DWD    2005  new obs data type
! Gerhard Paul    DWD    2008  input in NetCDF format
! Harald Anlauf   DWD    2008  fixes and optimizations for NEC SX8
! Alexander Cress DWD    2011  extensions for wind profilers
!==============================================================================

!-------------------------------------------------
! uncomment to exit in case of unknown BUFR codes:
!
!#define CHECKCODES
!-------------------------------------------------

  !=============
  ! modules used
  !=============
  !-------------------------
  ! general purpose routines
  !-------------------------
  use mo_exception,   only: finish,      &! abort routine
                            message       ! issue a warning
  use mo_kind,        only: wp,          &! working precision
                            sp,          &! single  precision
                            i2,          &! 2-byte integer
                            i8            ! 8-byte integer
  use mo_physics,     only: d2r,         &! conversion factor degree -> radians
                            gacc,        &! gravity acceleration
                            t0c,         &! 273.15 K
                            fd_uv         ! calculate ff,dd from u,v
  use mo_cntrlvar,    only: trh_tvgh      ! generalized humidity transform
  use mo_mpi_dace,    only: dace,        &! MPI group info
                            p_bcast       ! broadcast routine
  use mo_namelist,    only: position_nml,&! position namelist
                            nnml,        &! namelist Fortran unit number
                            POSITIONED    ! ok    code from position_nml
  use mo_fdbk_tables, only: LS_SURFACE,     &! level significance id
                            LS_STANDARD,    &!
                            LS_TROPO,       &!
                            LS_MAX,LS_SIGN, &!
                            LS_SUPEROBS,    &! superobservation (layer average)
                            OT_TEMP,        &! observation type id
                            OT_PILOT,       &!
                            OC_SCADA,       &! codetype: SCADA report
                            varno            ! 'varno' table
  use mo_t_table,     only: name_value       ! find name of table entry
  !-----------------------------
  ! access atmospheric data type
  !-----------------------------
  use mo_atm_grid,  only: t_grid,         &! type definition for grid info
                          VCT_P_HYB        ! GME vertical coordinate
  use mo_atm_state, only: t_atm            ! atm. state data type
  use mo_time,      only: init_time,      &! initialise time data type
                          operator(<),    &! compare times
                          operator(>=),   &! compare times
                          operator(==),   &! compare times
                          operator(/=),   &! compare times
                          operator(-),    &! subtract times
                          cyyyymmdd,      &! time -> string conversion routine
                          chhmmss,        &! time -> string conversion routine
                          minutes,        &! time differnce -> real
                          print,          &! generic print routine
                          cyyyymmddhhmm,  &! derive string from time
                          cyyyymmddhhmmss  ! derive string from time
  use mo_physics,   only: es_t,           &! sat. vap. pressure <- temperature
                          esw_t,          &! sat. vap. pressure over water
                          esw_t_hardy,    &!   " (ITS90 formulation, Hardy)
                          rearth           ! earth radius
!                         rh_q,           &! relative humidity  <- specific h.
!                         rh_q_adj         ! adjoint: rel.hum.  <- specific h.
  use mo_usstd,     only: p_h_usstd        ! derive p from h(gpm) US std.atm.
  !-----------------------------
  ! access observation data type
  !-----------------------------
  use mo_t_datum,   only: t_datum,        &! data type for one observed datum
                          rvind,          &! missing value indicator (real)
                          inv_datum,      &! constant for invalid datum
                          set_datum,      &! set t_datum% o  (observed value)
!                         set_qbits,      &! set t_datum% qc (quality bits)
                          cmp_datum,      &! compare datum
                          merge_datum,    &! merge datum
                          print,          &! generic print routine
                          QC_OK,          &!  quality flag: OK
                          QC_MISS,        &!                missing
                          QC_NOUSE,       &!                do not use
                          QC_INCON,       &!                inconsistent
                          QC_CLIM,        &!                climatological range
                          SRC_DER, SRC_DOK ! derived flag (3D-Var)
  use mo_t_use,     only: decr_use,       &! decrease the state of a datum
                          t_use,          &! status variable data type
                          use_0,          &! default values of type use
                          chk,            &! list of mnemonics for CHK_..
                          STAT_PASSIVE,   &!
                          STAT_OBS_ONLY,  &!
                          STAT_FORGET,    &!
                          STAT_DISMISS,   &!
                          CHK_INSDAT,     &!
                          CHK_NOIMPL,     &!
                          CHK_NOTUSED,    &!
                          CHK_CORR,       &!
                          CHK_CORRERR,    &!
                          CHK_MERGE,      &!
                          CHK_REDUNDANT,  &!
                          CHK_DBLERR,     &!
                          CHK_NONE,       &!
                          CHK_THIN,       &!
                          CHK_HEIGHT,     &!
                          CHK_DOMAIN       !
  use mo_t_obs,     only: t_obs,          &!
                          t_spot,         &!
                          t_head,         &! observation data type
                          derive_dbkz,    &! derive DBKZ if not present
                          new_spot,       &! reserve memory
                          new_obs,        &! reserve memory
                          new_par,        &! reserve memory
                          set_int_insitu, &! set interpolation space
                          set_xuv,        &! set unit vectors, zenith angle
                          unitvector,     &! get unit vector on the sphere
                          set_vqc_insitu, &! subroutine to set VQC bounds
                          shrink_report,  &! release unused obs. in report
                          source,         &! list   of Report source files
                          corme_conv,     &! corme handling (TEMP/BUFR)
                          monitor_ff,     &! flag to monitor wind speed
                          monitor_dd,     &! flag to monitor wind direction
                          TSK_INIT,       &! FLAGS: initialize module
                          COSMO,          &! module flag value
                          TEMP,           &! module flag value
                          TSK_READ,       &!  read observations
                          TSK_SET_CHR,    &!  set observation characteristics
                          TSK_SETUP_COLS, &!  determine model columns required
                          TSK_SETUP_FULL, &!  setup description of PSAS-space
                          TSK_SETUP_FUL0, &!  setup interpolation space
                          TSK_SHRINK,     &!  release unused obs. in report
                          TSK_K,          &!          evaluate linear operator
                          TSK_Y,          &!          evaluate nonlinear oper.
                          TSK_R,          &!   setup observational errors
                          CHR_ID,         &! H is the identity operator
                          CHR_NONL,       &! H is nonlinear
                          CHR_INV,        &! H is invertible
                          CHR_EXP,        &! H is 'expensive'
                          ITY_ICOL,       &! interpolation type: column
                          netcdf_verb      ! verbosity of NetCDF decoding
  use mo_wigos,     only: operator (==)    ! Same WIGOS station ID?
  use mo_fdbk_tables,only:VN_U, VN_V,     &! wind component code
                          VN_FF, VN_DD,   &! wind direction,speed code
                          VN_W,           &!
                          VN_U10M,VN_V10M,&!
                          VN_RH, VN_RH2M, &! rh             code
                          VN_T,  VN_T2M,  &! T              code
                          VN_Z,           &! geopotential   code
                          VN_HEIGHT,      &!
                          VN_NH,          &! cloud cover
                          VN_P,           &!         pressure
                          VN_PS,          &! surface pressure
                          VN_N_L, VN_N_M   ! cloud cover
  use mo_wmo_tables,only: WMO0_ECMWF       ! generating center
  use mo_obs_tables,only: decr_rpt_use,   &!
                          rept_stat,      &! report statistics table
                          idb_dbk,        &! index in table rept_stat
!                         rept_char,      &! observation type characteristics
                          check_report_0, &! init. flags, standard checks
                          check_report_1, &!
                          rept_use,       &!
                          obstyp,         &!
                          update_statistics!
  use mo_t_col,     only: t_cols,         &! model columns data type
                          set_lev_insitu, &! observation pres.levels from cols
                          COL_TV,         &! column flag: virt.temp.
                          COL_RH,         &!              rel. hum .
                          COL_UV,         &!              hor.wind comp.
                          COL_GEOH,       &!              geop.h. half levs.
                          COL_GEO,        &!              geop.h. full levs.
                          COL_PH           !              press.  half levs.
  use mo_bufr_dwd,  only: t_bufr,         &! BUFR record data type
                          inv_bufr,       &! indicator for invalid value
                          ty_loop_cnt,    &! indicator for loop counter
                          bufr_get_character,  &! decode BUFR character entry
#ifdef  CHECKCODES
                          bufr_get_entry_texts,&!
                          bufr_get_entry_units,&!
#endif
                          bufr_print_sections, &!
                          bufr_print_subset     !
  use mo_obs_set,   only: t_obs_block      ! obs data type
  use mo_obstypes,  only: t_obsid,        &! observation id table entry
                          obstype_dbkz,   &! derive obsids from dbkz
                          obstype_bufr     ! derive obsids from bufr type
  !--------------------------------------
  ! obs data header info read from NetCDF
  !--------------------------------------
  use mo_head_netcdf,only:ncid,           &! NetCDF file id
                          dimids_max,     &! max number of NetCDF dimension ids
                          imissing,       &! NetCDF _FillValue for integer
                          rmissing,       &! NetCDF _FillValue for reals
                          s2ikz,          &! DWD-internal classifier
                          s1cat,          &! data category
                          s1catls,        &! local data sub category
                          s1cent,         &! data centre
                          stime,          &! header observation time (section1)
                          db_time,        &! data bank time
                          s1cents,        &! data sub centre
                          istidn,         &! WMO numeric station number combined
                          s1updat,        &! update sequence no.
                          mlah,           &! latitude
                          mloh,           &! longitude
                          obs_time,       &! body observation time
                          ystidn,         &! any type of station identifier as variable
                          lwsi,           &! WIGOS station identifier valid
                          wsi,            &! WIGOS station id stored representation
                          wsihash,        &! Station id hash
                          get_real         ! read real variable from NetCDF file
  !---------------------
  ! netCDF f90 interface
  !---------------------
  use netcdf,       only: nf90_Inquire,          &!
                          nf90_Inquire_Dimension,&!
                          nf90_Inquire_Variable, &!
                          nf90_inq_varid,        &!
                          nf90_get_var,          &!
                          nf90_get_att,          &!
                          NF90_MAX_NAME,         &!
                          NF90_NOERR              !
  !------------------------
  ! access matrix data type
  !------------------------
  use mo_dec_matrix,only: t_vector_segm    ! vector segment data type

  !----------------------------------
  ! set observation error covariances
  !----------------------------------
  use mo_obs_err,   only: temp_obs_err
  !------------------------------------
  ! Interface to interpolation routines
  !------------------------------------
  use mo_grid_intpol,only: idx_init, &
                           Grid_Indices
  implicit none
!------------------------------------------------------------------------------
  !================
  ! public entities
  !================
  private

  !----------------------------------------
  ! general operations on TEMP observations
  !----------------------------------------
  public :: process_temp     ! calculate costfunction, derivatives, ...
  public ::   print_temp     ! TEMP print routine
  public :: check_dbl_pilot  ! dismiss PILOTs collocated with TEMPs
  public :: merge_report     ! merge parts A,B,C,D
  public :: split_report     ! split report if it drifts out of grid-cell
  public :: select_levels    ! select levels to use
  public :: reduce_profiles  ! reduce profile information (average/thin)
! public :: check_temp_cons  ! check TEMP report internal consistency
  public :: read_fdbk_temp   ! fix level count after fof-input.
  !---------------------------------------------------------------
  ! specific data type and routines passed to BUFR reading routine
  !---------------------------------------------------------------
  public ::           t_temp ! data type holding levels
  public :: check_store_temp ! check and store data into generic data type
  public :: read_temp_bufr   ! read TEMP observation from BUFR record
  public :: read_temp_netcdf ! read TEMP observation from netCDF file
  public :: read_temp_nml    ! read namelist related to TEMP observation
  !---------------------
  ! (namelist) variables
  !---------------------
  public :: fgchk_inv_t      ! stability-dep. factor for fg-check limit for T
  public :: fgchk_inv_q      ! stability-dep. factor for fg-check limit for RH
!------------------------------------------------------------------------------
  !================================
  ! Module variables and data types
  !================================

  !-----------------------------------------------------------------
  ! temporary missing value for level_sig in feedback file (set msb)
  !-----------------------------------------------------------------
  integer,     parameter :: LS_MISSING  = digits (inv_datum % lev_sig)
  integer,     parameter :: LS_OTHER    = LS_MISSING - 1
  integer(i2), parameter :: LEV_MISSING = ibset (0_i2, LS_MISSING)

  !-------------------------
  ! private module variables
  !-------------------------
  integer ,save      :: temp_int_size =    0 ! integer size of type t_temp
  integer ,parameter :: inv_secs      = -999 ! invalid time offset [s]

  !-------------------------------
  ! special data base 'Kennzahlen'
  !-------------------------------
  integer :: kz_pilot (12) = (/508,509,510,511,512,513,514,515,548,10553,10600,584/)
  integer :: kz_temp   (9) = (/520,521,522,523,10520,10521,10526,10527,10574/)
  integer :: kz_tship  (9) = (/776,777,778,779,10776,10777,10782,10783,10785/)
  integer :: kz_pship  (8) = (/764,765,766,767,768,769,770,771/)
  integer :: kz_tmobil(11) = (/516,517,518,519,524,525,526,527,10516,10517,10570/)
  integer :: kz_tdrop  (5) = (/780,781,782,783,10780/)
  !-------------------------
  ! Corresponding codetypes:
  !-------------------------
  ! PILOT     : 32
  ! PILOT SHIP: 33
  ! WIND Prof.: 34, 132, ...
  ! TEMP      : 35,  109 (BUFR), 231 (BUFR, descent)
  ! TEMP SHIP : 36,  111 (BUFR)
  ! TEMP MOBIL: 37
  ! TEMP DROP : 135, 230 (BUFR)
  !------------------
  ! special codetypes
  !------------------
  integer :: ct_wprof (6) = (/132,133,134,135,136,137/)

  !-----------------------------------------------------------------------------
  ! Vertical Sounding Significance 0 08 001 (7 bits)           (BUFR Bit number)
  !-----------------------------------------------------------------------------
  integer, parameter :: lev_other = 128 ! Anything else       (3dvar convention)
  integer, parameter :: lev_surf  =  64 ! Surface                            (1)
  integer, parameter :: lev_std   =  32 ! Standard level                     (2)
  integer, parameter :: lev_tropo =  16 ! Tropopause level                   (3)
  integer, parameter :: lev_max_v =   8 ! Maximum wind level                 (4)
  integer, parameter :: lev_sig_t =   4 ! Significant level, temperature/rh  (5)
  integer, parameter :: lev_sig_v =   2 ! Significant level, wind            (6)
  integer, parameter :: lev_miss  =   1 ! Missing value                  (All 7)
  integer, parameter :: lev_all   = 255 ! All
  !-----------------------------------------
  ! fortran bit count for significant levels
  !-----------------------------------------
  integer, parameter :: lbit_other  = 7
  integer, parameter :: lbit_surf   = 6
  integer, parameter :: lbit_std    = 5
  integer, parameter :: lbit_tropo  = 4
  integer, parameter :: lbit_max_v  = 3
  integer, parameter :: lbit_sig_t  = 2
  integer, parameter :: lbit_sig_v  = 1
  integer, parameter :: lbit_miss   = 0

  !--------------------------------------------------------------------------------
  ! EXTENDED Vertical Sounding Significance 0 08 042 (18 bits)    (BUFR Bit number)
  !--------------------------------------------------------------------------------
  integer, parameter :: lev_surf_18  = 2**17   ! Surface                        (1)
  integer, parameter :: lev_std_18   = 2**16   ! Standard level                 (2)
  integer, parameter :: lev_tropo_18 = 2**15   ! Tropopause level               (3)
  integer, parameter :: lev_max_v_18 = 2**14   ! Maximum wind level             (4)
  integer, parameter :: lev_sig_t_18 = 2**13   ! Significant level, temperature (5)
  integer, parameter :: lev_sig_h_18 = 2**12   ! Significant level, humidity    (6)
  integer, parameter :: lev_sig_v_18 = 2**11   ! Significant level, wind        (7)
  integer, parameter :: lev_regn_18  = 2**3    ! Determ. by regional decision  (15)
  integer, parameter :: lev_miss_18  = 1       ! Missing value             (All 18)
  integer, parameter :: lev_all_18   = 2**18-1 ! All


  !------------------
  ! Namelist TEMP_OBS
  !------------------
  real(wp)  :: top_p       =   0._wp ! top level for radiosonde profiles  (hPa)
  real(wp)  :: top_q       = 300._wp ! top level for humidity observations(hPa)
  integer   :: top_q_mode  = 0       ! humidity obs. handling above tropopause
                                     ! = 0 : use    TD/RH up to top_q (always)
                                     ! = 1 : ignore TD/RH      above tropopause
                                     ! > 1 : ignore TD/RH at & above tropopause
  integer   :: satvp_form  = 1       ! Formula for sat. vapor pres. over water
                                     ! 1: Magnus-Tetens, 2: Hardy
  logical   :: prt_data    = .false.
  real(wp)  :: vmainpl(24) = -1._wp  ! main pressure level values [hPa]
  real(wp)  :: wmainpl(30) = -1._wp  ! pressure level values for profilers
  real(wp)  :: dz_coloc    = -1._wp  ! vert. separation for collocations (Pa)
  real(wp)  :: dzv_coloc   = -1._wp  ! vert. separation for zero wind
  real(wp)  :: dx_coloc    =  1._wp  ! hori. separation for collocations (km)
  logical   :: rtl_coloc   = .true.  ! remove temp levels in collocations
  integer   :: active_t    = lev_std ! + lev_sig_t + lev_tropo ! active  levels
  integer   :: active_h    = lev_surf
  integer   :: active_q    = lev_std ! + lev_sig_t
  integer   :: active_v    = lev_std ! + lev_sig_v + lev_max_v
  integer   :: active_v_wp = lev_std ! for wind profilers
  integer   :: passive_t   = lev_all                           ! passive levels
  integer   :: passive_h   = lev_all
  integer   :: passive_q   = lev_all
  integer   :: passive_v   = lev_all
  logical   :: split       = .false. ! split report if drifting out of gridcell
  logical   :: std_t0miss  = .true.  ! use standard levels with zero/miss. time
  logical   :: rm_t0c_ff0  = .true.  ! dismiss levels with apparent miss.values
  integer   :: corme_mode  = 0       ! Mode for correction message handling
  integer   :: check_cons  = 0       ! Check consistency of report
  integer   :: thin_bufr(2)= 0       ! thinning interval for BUFR TEMPs (s)
  logical   :: layer_avg   = .false. ! apply layer averaging (superobbing)?
  integer   :: supob_thr(2) = 300    ! levels in report required for superobbing
                                     ! (ascending, descending)
  integer   :: prefer       = 1      ! level use preference: 1=std, 2=layer-avg.
  real(wp)  :: avl_levs(15) = 0._wp  ! averaging layers: level definition list
  real(wp)  :: avl_incr(15) = 0._wp  ! averaging layers: level increments
  real(wp)  :: dp_asurf     = 5._wp  ! dist. of lowest layer above surface [hPa]
  real(wp)  :: fgchk_inv_t  = 0._wp  ! stability/inversion-dep. enhancement of
                                     !   fg-check limit for T (recommended: .8)
  real(wp)  :: fgchk_inv_q  = 0._wp  ! stability/inversion-dep. enhancement of
                                     !   fg-check limit for RH (recommended: .6)
                                     !   (== 0. means: no enhancement,
                                     !    valid range: 0. <= fgchk_inv_x <= 1.)
  namelist  /TEMP_OBS/ active_t, active_q, active_v, active_h,     &
                       passive_t, passive_q, passive_v, passive_h, &
                       vmainpl, wmainpl, top_p, top_q, prt_data,   &
                       dz_coloc, dzv_coloc, dx_coloc, rtl_coloc,   &
                       split, corme_mode, check_cons, thin_bufr,   &
                       layer_avg, supob_thr, prefer, avl_levs,     &
                       avl_incr, dp_asurf, active_v_wp,            &
                       fgchk_inv_t, fgchk_inv_q, top_q_mode,       &
                       std_t0miss, rm_t0c_ff0, satvp_form
  target :: vmainpl, wmainpl

  !-------------------------------------------------------------------
  ! almost private temp datatype (only passed to BUFR reading routine)
  ! holds one level of TEMP observations
  !-------------------------------------------------------------------
  type t_temp
    type(t_datum) :: p               ! pressure
    type(t_datum) :: t               ! temperature
    type(t_datum) :: td              ! dew point temperature
    type(t_datum) :: rh              ! relative humidity (derived)
    type(t_datum) :: ff              ! wind direction
    type(t_datum) :: dd              ! wind speed
    type(t_datum) :: uu              ! wind speed
    type(t_datum) :: vv              ! wind speed
    type(t_datum) :: gp              ! geopotential
    real          :: dlat = -999.    ! latitude  (high accuracy) [deg]
    real          :: dlon = -999.    ! longitude (high accuracy) [deg]
    integer       :: secs = 0        ! time since launch (seconds)
    integer       :: vss  = lev_miss ! vertical sounding significance
    real          :: height = -999.  ! profiler, no pressure level reported [m]
  end type t_temp

  !-------------------------------------------------------------------
  ! general logical
  !-------------------------------------------------------------------
  logical       :: lpr_temp = .false. ! print temp/pilot reports from netcdf
!==============================================================================
contains
!==============================================================================
  subroutine process_temp (task, spot, obs, atm, cols, xi, y, Jo, Jo_atm, &
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
    integer                :: tsk         ! task (local copy)
    !-----------------------------------
    ! arguments to interpolation routine
    !-----------------------------------
    integer                :: i,j,k,l,n,ii,io,iu,ju
    type(t_temp)  ,pointer :: o  (:)      ! observations
    type(t_datum) ,pointer :: ob (:)      ! observation body
    real(wp)      ,pointer :: lv (:)      ! levels (ln p)
    integer       ,pointer :: ty (:)      ! observation types
    real(wp)  ,allocatable :: e  (:)      ! observational errors
    logical                :: change
    real(wp)  ,allocatable :: Hnew(:)     ! Jacobi matrix
    integer                :: jtv, jrh    ! indices
    real(wp)  ,allocatable :: yn  (:)

    real(wp)               :: p           ! pressure
    real(wp)               :: rh          ! relative humidity
    real(wp)               :: t           ! temperature
    real(wp)               :: gh          ! generalized humidity
    real(wp)               :: tv          ! virtual temperature
    real(wp)               :: dt_tv, dt_gh, drh_tv, drh_gh ! gradients
    real(wp)               :: ff, dd      ! wind speed, direction
    real(wp)               :: ff_x, ff_y, dd_x, dd_y ! gradients
    integer                :: ifail       ! error code
    integer(i8)            :: col_geop    ! p-lev for geop.height
    integer                :: nlev
    integer       ,pointer :: lindx (:)
    real(wp)               :: ol          ! temporary: obs.level

    !======================
    ! executable statements
    !======================
    tsk = task
    if(tsk==0) return
    !----------------------------------------------------
    ! In dependence of the value of flag TASK
    ! proceed as follows:
    !
    !  TSK_INIT       ! initialize modules
    !  TSK_READ       ! read observations
    !  TSK_SET_CHR    ! set observation characteristics
    !  TSK_SETUP_COLS ! setup columns
    !  TSK_SETUP_FUL0 ! setup interpolation space
    !  TSK_SETUP_FULL ! setup description of PSAS-space
    !  TSK_R          ! setup observational error
    !  TSK_Y          ! run forward operator
    !  TSK_H          ! run tangent linear operator
    !  TSK_K          ! evaluate linear operator
    !----------------------------------------------------
    !==================================================
    ! skip tasks not required for this observation type
    !==================================================
    tsk = iand (task, not (&
          TSK_READ))        ! BUFR file is read in module mo_obs
    if (tsk == 0) return
    !===============
    ! initialisation
    !===============
    if (iand (TSK_INIT,tsk) /= 0) then
      call read_temp_nml
      tsk = tsk - TSK_INIT
      if (tsk == 0) return
    endif
    !=============================================
    ! TSK_SET_CHR: set observation characteristics
    !=============================================
    if (iand (TSK_SET_CHR,tsk) /= 0) then
      tsk = tsk - TSK_SET_CHR
      spot% int_type  = ITY_ICOL
      spot% cost      = spot% col% nlev
      spot% nr        = spot% o%n
!     if(count(obs%o% varno (spot%o%i+1:spot%o%i+spot%o%n) == VN_Z) > 1) &
!       spot% nr      = spot% o%n * spot% o%n    ! geop. correlations only
      spot% char      = CHR_ID
      if (any (obs%o% varno (spot%o%i+1:spot%o%i+spot%o%n) == VN_T   &
          .or. obs%o% varno (spot%o%i+1:spot%o%i+spot%o%n) == VN_RH))&
        spot% char = CHR_NONL + CHR_INV + CHR_EXP
      if (any (obs%o% varno (spot%o%i+1:spot%o%i+spot%o%n) == VN_FF  &
          .or. obs%o% varno (spot%o%i+1:spot%o%i+spot%o%n) == VN_DD))&
        spot% char = CHR_NONL + CHR_EXP
      if (tsk == 0) return
    endif
    !==================================================================
    ! tsk == TSK_SHRINK:
    ! release unused observations in the report
    !==================================================================
    if (iand (TSK_SHRINK,tsk) /= 0) then
      call shrink_report (spot, obs%o, state, change, lindx=lindx, nl=nlev)
      if (change) then
        call set_int_insitu (spot, obs% o)
        if (spot% p% n > 0) then
          call load_temp  (obs% o, spot, o)
          spot% col% nlev = nlev
          call store_temp (obs% o, spot, o(lindx))
          deallocate (o)
        else
          spot% col% nlev = nlev
        endif
        deallocate (lindx)
      endif
      tsk = tsk - TSK_SHRINK
      if (tsk == 0) return
    endif
    !========================
    ! determine model columns
    !========================
    if (iand (TSK_SETUP_COLS,tsk) /= 0) then

      select case (atm% grid% vct)
      case (VCT_P_HYB)
        col_geop = COL_GEOH+COL_PH ! GME/HRM: get geop.height from half levels
      case default
        col_geop = COL_GEO         ! COSMO/ICON:              from full levels
      end select

      call idx_init (                  &
            spot% col% c,              &! <-  column descriptor
            spot% col% h,              &!  -> interpolation coefs.
            obs% o% mc,                &! <-> model column descriptors
            COL_TV + COL_RH + COL_UV + &! <-  fields required
            col_geop,                  &!
            0,                         &! <-  tracers required
            atm% grid,                 &! <-  model grid
            spot% i_time,              &! <-  time slot
            spot% w_time               )! <-  time interpolation weight

      if (spot% col% h% imc(1,1) == 0) call decr_rpt_use (spot, CHK_DOMAIN)

      tsk = tsk - TSK_SETUP_COLS
      if (tsk == 0) return
    endif
    !===============
    ! PSAS-space
    !===============
    if (iand (TSK_SETUP_FUL0,tsk) /= 0) then
      !-----------------------------------------------------------
      ! interpolate pressure if height is the independent quantity
      !-----------------------------------------------------------
      if (dace% pe == obs% o% pe .and. spot% hd% obstype == OT_PILOT) then
         call set_lev_insitu (spot, obs% o, cols)
      end if
      !---------------------------
      ! set up interpolation space
      !---------------------------
      call set_int_insitu (spot, obs% o)
      tsk = tsk - TSK_SETUP_FUL0
      if (tsk == 0) return
    endif
    !===========================================
    ! setup description of PSAS-space
    ! observed values were set up while reading.
    !===========================================
    if (iand (TSK_SETUP_FULL,tsk) /= 0) then
      tsk = tsk - TSK_SETUP_FULL
      if (tsk == 0) goto 888
    endif

    !===============================================
    ! tsk == TSK_R
    ! set up R (observation error covariance matrix)
    !===============================================
    if (iand (TSK_R,tsk) /= 0) then
      n = spot% o% n
      i = spot% o% i
      ob => obs% o% body  (i + 1 : i + n)
      ty => obs% o% varno (i + 1 : i + n)
      lv => obs% o% olev  (i + 1 : i + n)
      if (obs% o% pe == dace% pe) then
        !==========================================
        ! setup TEMP observational error covariance
        !==========================================
        allocate (e (n))
        if (spot% p% n > 0) then
          call load_temp (obs% o, spot, o)
          call temp_obs_err (obs,spot,e,ty,lv,ob%o, o%p%o, o%t%o, lcov=.TRUE.)
          deallocate (o)
        else
          e = ob% eo
        endif

        k = obs% R% ia (i + 1)
        do j=1,n
          obs% R% ia (i + j) = k
          obs% R% packed (k) = e(j) ** 2
          obs% R% ja (k) = i + j
          k = k + 1
        end do
        obs% R% ia (i + n + 1) = k

        deallocate (e)
        !=================
        ! setup vqc bounds
        !=================
        call set_vqc_insitu (spot, obs% o)
      endif
      tsk = tsk - TSK_R
      if (tsk == 0) goto 888
    endif

    !=================================================
    ! evaluate observation operator, set up H operator
    !=================================================
    if (iand (TSK_K+TSK_Y,tsk) /= 0) then
      if (spot% pe_eval == dace% pe) then
        n = spot% o% n
        if (iand (TSK_K,tsk) /= 0) then
          allocate (Hnew(spot% o% n*spot% i% n))
          allocate (yn  (spot% o% n))
          Hnew = huge(0._wp)
        endif
        j  = 1                  ! Index of output variable(s)
        ii = spot%i% i + 1      ! Index of input  variable(s)
        ol = -huge(ol)
        do i = 1, spot%o% n
          io = spot%o% i + i
          if (ol /= obs% o% olev(io)) then
            ol = obs% o% olev(io)
            p  = obs% o% body(io)% plev
            tv = -1._wp
            ff = -1._wp
            iu = -1             ! Index of first wind input  component (u)
            ju = -1             ! Index of first wind output component
          endif
          select case (obs% o% varno(io))
          case (VN_Z)
            if (present(y)) y%x (io) = xi%x (ii)
            if (iand (TSK_K,tsk) /= 0) then
              yn(i) = xi%x (ii)
              l     = (j - 1) * n + i
              Hnew (l) = 1._wp
            endif
            ii = ii + 1         ! Update input  index
            j  = j  + 1         ! Update output index

          case (VN_U, VN_V)
            if (iu == -1) then
              iu = ii           ! Remember index of input  u component
              ju = j            ! Remember index of output variable u
              ii = ii + 2       ! u,v (input)  always needed  in pairs
              j  = j  + 2       ! U,V (output) always treated in pairs
            endif
            if (obs% o% varno(io) == VN_U) then
              if (present(y)) y%x (io) = xi%x (iu)
              if (iand (TSK_K,tsk) /= 0) then
                yn(i) = xi%x (iu)
                Hnew ((ju - 1) * n + i) = 1._wp
              endif
            endif
            if (obs% o% varno(io) == VN_V) then
              if (present(y)) y%x (io) = xi%x (iu+1)
              if (iand (TSK_K,tsk) /= 0) then
                yn(i) = xi%x (iu+1)
                Hnew ((ju    ) * n + i) = 1._wp
              endif
            endif

          case (VN_FF, VN_DD)
            if (iu == -1) then
              iu = ii           ! Remember index of input  u component
              ju = j            ! Remember index of output variable
              ii = ii + 2       ! u,v (input) always needed in pairs
              j  = j  + 2       ! FF and/or DD adjoint calc. treated at once
            endif
            if (ff < 0._wp) then
              call fd_uv (ff, dd, xi% x(iu), xi% x(iu+1), &
                          ff_x, ff_y, dd_x, dd_y          )
            endif
            if (obs% o% varno(io) == VN_FF) then
              if (present(y)) y%x (io) = ff
              if (iand (TSK_K,tsk) /= 0) then
                yn(i) = ff
                Hnew ((ju - 1) * n + i) = ff_x
                Hnew ((ju    ) * n + i) = ff_y
              endif
            endif
            if (obs% o% varno(io) == VN_DD) then
              if (present(y)) y%x (io) = dd
              if (iand (TSK_K,tsk) /= 0) then
                yn(i) = dd
                Hnew ((ju - 1) * n + i) = dd_x
                Hnew ((ju    ) * n + i) = dd_y
              endif
            endif

          case (VN_T, VN_RH)
            if (tv < 0._wp) then
              jtv = j
              jrh = j+1
              tv = xi%x (ii)
              gh = xi%x (ii+1)
              call trh_tvgh (t,                            &!  -> t
                             rh,                           &!  -> rh
                             tv,                           &! <-  tv
                             gh,                           &! <-  gh
                             p,                            &! <-  p
                             ifail,                        &! -> error code
                             dt_tv, dt_gh, drh_tv, drh_gh)  !  -> gradient
              if (ifail < 0) then
                write(6,*)'***************************************************'
                write(6,*)'process_temp (TSK_K): trh_tvgh failed:',ifail
                write(6,*)'station : ', trim (spot% statid)
                write(6,*)'location: lat,lon,z =', &
                     [ real :: spot%col%c% dlat, spot%col%c% dlon, spot% z ]
                write(6,*)'i,j,ii  =',i,j,ii
                write(6,*)'tv,gh,p =',tv,gh,p
                write(6,*)'t,rh    =',t,rh
                write(6,*)'ni,no,  =',spot%i% n, spot%o% n
                write(6,*)'xi      =',xi%x (spot%i%i+1:spot%i%i+spot%i%n)
                write(6,*)'t_int   =',obs%o%t_int(spot%i%i+1:spot%i%i+spot%i%n)
                write(6,*)'varno   =',obs%o%varno(spot%o%i+1:spot%o%i+spot%o%n)
                write(6,*)'olev    =',obs%o%olev (spot%o%i+1:spot%o%i+spot%o%n)
                write(6,*)'***************************************************'
                call finish ('process_temp (TSK_K)','trh_tvgh failed')
              endif
              ii = ii + 2       ! tv,gh (input) always needed in pairs
              j  = j  + 2       ! T/RH (output) adjoint calc. treated at once
            endif
            if (obs% o% varno(io) == VN_T) then
              if (present(y)) y%x (io) = t
              if (iand (TSK_K,tsk) /= 0) then
                yn(i) = t
                Hnew ((jtv - 1) * n + i) = dt_tv
                Hnew ((jrh - 1) * n + i) = dt_gh
              endif
            else if (obs% o% varno(io) == VN_RH) then
              if (present(y)) y%x (io) = rh
              if (iand (TSK_K,tsk) /= 0) then
                yn(i) = rh
                Hnew ((jtv - 1) * n + i) = drh_tv
                Hnew ((jrh - 1) * n + i) = drh_gh
              endif
            endif
          case default
            if (present(y))            y%x (io) = rvind
            if (iand (TSK_K,tsk) /= 0) yn  (i)  = rvind
!           write (0,*) 'process_temp, TSK_K: invalid obstype',&
!                        obs% o% varno(io)
!           call finish('process_temp, TSK_K','invalid observation type')
          end select
        end do
        if (iand (TSK_K,tsk) /= 0) then
          k = obs% H% ia (spot% i% i + 1)
          do j=1,spot% i% n                             ! columns
            obs% H% ia (spot% i% i +j) = k              ! column index
            do i=1,spot% o% n                           ! rows
              l = (j - 1) * n + i
              if (Hnew(l) /= huge(0._wp)) then
                obs% H% packed (k) = Hnew(l)            ! coefficient
                obs% H% ja (k) = spot% o% i + i         ! row index
                k = k + 1
              endif
            end do
          end do
          obs% H% ia (spot% i% i + spot% i% n + 1) = k  ! column index
!         obs% xi% x (spot% i% i+1:spot% i% i+spot% i% n) = &
!              xi% x (spot% i% i+1:spot% i% i+spot% i% n)
          obs% yi% x (spot% o% i+1:spot% o% i+spot% o% n) = yn (1:spot%o% n)
          deallocate (Hnew)
          deallocate (yn)
        endif
      endif
      if (iand (TSK_K,tsk) /= 0) tsk = tsk - TSK_K
      if (iand (TSK_Y,tsk) /= 0) tsk = tsk - TSK_Y
      if (tsk == 0) goto 888
    endif
    !=============
    ! unknown task
    !=============
    if (tsk /= 0) then
      write (0,*)  'process_temp:  unknown task',tsk
      call finish ('process_temp','unknown task')
    endif
    !========================
    ! release local variables
    !========================
888 continue
  end subroutine process_temp
!------------------------------------------------------------------------------
  subroutine print_temp (spot ,obs)
  type(t_spot) ,intent(in) :: spot
  type(t_obs)  ,intent(in) :: obs
    character(len=3)       :: c
!   type (t_temp) ,pointer :: o(:)
!   integer                :: k
!++++++++++++++++++++++
! Currently not working
!++++++++++++++++++++++
    c = '____'
    print *
    print *,'type (t_spot)'
    print *,c,'id          :',spot% id
    print *,c,'type        : TEMP'
    print *,c,'ident       :',spot% ident
    print *,c,'header_time :',cyyyymmdd(spot% hd% time), &
                              chhmmss  (spot% hd% time)
    print *,c,'actual_time :',cyyyymmdd(spot% actual_time), &
                              chhmmss  (spot% actual_time)
!   print *,c,'decode_time :',cyyyymmdd(spot% decode_time), &
!                             chhmmss  (spot% decode_time)
    print *,c,'col% dlon   :',spot% col% c% dlon
    print *,c,'col% dlat   :',spot% col% c% dlat
    print *,c,'z           :',spot% z
    print *,c,'col% x      :',spot% col% c% x
    print *,c,'cost        :',spot% cost
    print *,c,'col% nlev   :',spot% col% nlev
    print *,c,'o           :',spot% o
    print *,c,'i           :',spot% i
    print *,c,'p           :',spot% p
    if (spot% o% n > 0) then
      write(6,'(a)') &
      '        p         gp        t         td        ff        dd'
!     if (associated (obs% obs)) then
!       call load_temp (obs ,spot ,o)
!!!        call unpack_temp (spot ,obs, obs% obs, o)
!       c = 'obs_'
!       do k=1,spot% col% nlev
!         write(6,'(a,f10.0,f10.0,f10.2,f10.2,f10.2,f10.2)') &
!                   c, o(k)%p,o(k)%gp, o(k)%t,o(k)%td,o(k)%ff,o(k)%dd
!       end do
!       deallocate (o)
!     endif
!     if (associated (obs% ana)) then
!!!        call unpack_spot (spot ,obs, obs% ana, o)
!       c = 'ana_'
!       do k=1,spot% col% nlev
!         write(6,'(a,f10.0,f10.0,f10.2,f10.2,f10.2,f10.2)') &
!                   c, o(k)%p,o(k)%gp, o(k)%t,o(k)%td,o(k)%ff,o(k)%dd
!       end do
!       deallocate (o)
!     endif
    endif
  end subroutine print_temp
!==============================================================================
  subroutine read_temp_bufr (bufr, spt, obs, lkeep, cc)

  type (t_bufr) ,intent(inout)        :: bufr  ! BUFR record to decode
  type (t_spot) ,intent(inout)        :: spt   ! meta information to set
  type (t_obs)  ,intent(inout)        :: obs   ! observations data type to set
  logical       ,intent(out)          :: lkeep ! accept observation ?
  integer       ,intent(in) ,optional :: cc    ! part of year ccyy

    integer                :: nlev     ! number of levels
    type (t_temp) ,pointer :: btemp(:) !
    integer                :: i
    logical                :: unused   ! true for unknown BUFR entries
    integer                :: novalid  ! number of valid observations in BUFR
    !-----------------------------------------
    ! quantities derived from the BUFR message
    !-----------------------------------------
    integer         :: ival  ! value decoded from BUFR (integer)
    real(sp)        :: rval  ! value decoded from BUFR (real)
    character(len=8):: ymnem ! value decoded from BUFR (string)
    integer         :: itype ! type of BUFR message
    integer         :: yyyy,mo,dd,hh,mi ! actual time read

    !----------------
    ! index variables
    !----------------
    integer       :: is               ! sub-set    index
    integer       :: ie               ! entry in sub-set
    integer       :: id               ! descriptor index

    yyyy  = 0      ! date of observation
    mo    = 0      ! "
    dd    = 0      ! "
    hh    = 0      ! "
    mi    = 0      ! "

    !-------------------------------------------
    ! estimate number of levels, allocate memory
    !-------------------------------------------
    lkeep = .false.
    if (bufr% sec3% num_subsets /= 1) then
      call decr_rpt_use (spt, CHK_NOIMPL, &
                         comment='read_temp_bufr: cannot handle subsets')
      return
    endif
    is = 1
    nlev = 0
    do ie=1,bufr% nbufdat(is)
      id = bufr% idescidx(is,ie)
      itype = bufr% itype(id)
      if (itype == ty_loop_cnt) nlev = nlev + 1
    end do
    if(nlev < 1) then
      call decr_rpt_use (spt, CHK_INSDAT, comment='read_temp_bufr: no. levels < 1')
      return
    endif

    !------------------------------------
    ! correction message handling for SKY
    !------------------------------------
    if (corme_conv >= 1) spt% corme = bufr% sec1% update
    !---------------------
    ! loop over data, copy
    !---------------------
    is = 1
    allocate (btemp (nlev))
    call construct_temp (btemp)
    nlev   = 0
    lkeep  = .true.
    unused = .false.
    do ie=1,bufr% nbufdat(is)
      !---------------------------------------
      ! decode real/character datum, mnemonics
      !---------------------------------------
      ival  = bufr% ibufdat (is,ie)
      id    = bufr% idescidx(is,ie)
      itype = bufr% itype(id)
      ymnem = bufr% ymnem(id)
      rval  = rvind
      if (bufr% is_char(id) == 0) then
        IF (ival /= inv_bufr) &
          rval = ival * bufr% scale(id)
      endif

      if (itype == ty_loop_cnt) then
        nlev = nlev + 1
        spt% col% nlev = spt% col% nlev + 1
      endif

      if(ival /= inv_bufr) then
        select case (ymnem)
        !----------------------
        ! station values follow
        !----------------------
        case ('YDDDD') ! SHIP OR MOBILE LAND STATION IDENTIFIER
          call bufr_get_character (bufr, ival, spt% statid)
        case ('YXXNN') ! aircraft flight number (for drop sondes)
          if (spt% statid=='') call bufr_get_character (bufr,ival,spt% statid)
!       case ('YCCCC') ! ICAO LOCATION INDICATOR
!         call bufr_get_character (bufr, ival, spt% statid)
!       case ('MABNN') ! BUOY/PLATFORM IDENTIFIER
!         write(spt%  statid,'("BP_",i5.5)') ival
        case ('MII')   ! WMO   block number
          spt% ident  = 1000      * int(rval)
        case ('NIII')  ! WMO station number
          spt% ident  = spt% ident + int(rval)
        case ('MOGC')  ! generating centre
          if(rval/=rvind) spt% hd% center = int(rval)
        case ('MGAPP') ! generating application
        case ('NRARA') ! radiosonde type
          if(rval/=rvind) spt% sttyp = int(rval)
        case ('MSR')   ! radiosonde computational method
          if(rval/=rvind) spt% stret = int(rval)
        case ('MJJJ')  ! year
          yyyy = ival
          if(yyyy<100 .and. present(cc)) yyyy = yyyy + cc * 100
        case ('MMM')   ! month
          mo   = ival
        case ('MYY')   ! day
          dd   = ival
        case ('MGG')   ! hour
          hh   = ival
        case ('NGG')   ! minute
          mi   = ival
        case ('MLAH')  ! latitude  (high accuracy)                    [deg]
          spt% col% c% dlat = rval
        case ('MLOH')  ! longitude (high accuracy)                    [deg]
          spt% col% c% dlon = rval
        case ('MLALA')  ! latitude  (coarse accuracy)                 [deg]
          if (spt% col% c% dlat == rvind) spt% col% c% dlat = rval
        case ('MLOLO')  ! longitude (coarse accuracy)                 [deg]
          if (spt% col% c% dlon == rvind) spt% col% c% dlon = rval
        case ('MHP')   ! height of station                            [m]
          spt%   z     = rval
        case ('MCORME')! cor message mark                      [CODE_TABLE]
          if (corme_conv >= 1) spt% corme = ival
        !----------------------------------
        ! values for pressure levels follow
        !----------------------------------
                       !-------------------------
        case ('MPN')   ! pressure (vert.location)             [Pa]
                       !-------------------------
          call set_datum (btemp(nlev)% p ,rval ,spt% corme)
                       !-----------------------------
        case ('MHN')   ! geopotential (vert.location)         [m^2/s^2]
                       !-----------------------------
          call set_datum (btemp(nlev)% gp ,rval ,spt% corme)
          btemp(nlev)% gp% o = btemp(nlev)% gp% o / gacc     ! -> gpm
                       !-------------
        case ('NHNHN') ! geopotential                         [m^2/s^2]
                       !-------------
          call set_datum (btemp(nlev)% gp ,rval ,spt% corme)
          btemp(nlev)% gp% o = btemp(nlev)% gp% o / gacc     ! -> gpm
                       !---------------------
        case ('MTN')   ! dry bulb temperature                 [K]
                       !---------------------
          call set_datum (btemp(nlev)% t ,rval ,spt% corme)
                       !----------------------
        case ('MTDN')  ! dew point temperature                [K]
                       !----------------------
          call set_datum (btemp(nlev)% td ,rval ,spt% corme)
                       !---------------
        case ('NDNDN') ! WIND DIRECTION                       [degree]
                       !---------------
          call set_datum (btemp(nlev)% dd ,rval ,spt% corme)
                       !-----------
        case ('NFNFN') ! WIND SPEED                           [m/s]
                       !-----------
          call set_datum (btemp(nlev)% ff ,rval ,spt% corme)
                       !-------------------------------
        case ('MVTSO') ! VERTICAL SOUNDING SIGNIFICANCE
                       !-------------------------------
          btemp(nlev)% vss = ival
          call ffvss_from_vss (btemp(nlev))  ! to feedback file convention
        !--------------
        ! BUFR4 entries
        !--------------
                       !----------------------------------------
        case ('MEVSS') ! EXTENDED VERTICAL SOUNDING SIGNIFICANCE  [flag_table]
                       !----------------------------------------
          if (iand (ival, lev_miss_18 ) /= 0) ival = lev_miss_18
          btemp(nlev)% vss = ival
          call ffvss_from_evss (btemp(nlev))     ! to feedback file convention
          btemp(nlev)% vss = vss_from_evss(ival) ! translate to old convention
                                                 ! .. for further processing
                       !--------------------
        case ('NHHH')  ! GEOPOTENTIAL HEIGHT                 [gpm]
                       !--------------------
          call set_datum (btemp(nlev)% gp ,rval ,spt% corme)
                       !--------------------
        case ('NHHHN') ! GEOPOTENTIAL HEIGHT                 [gpm]
                       !--------------------
          call set_datum (btemp(nlev)% gp ,rval ,spt% corme)
                       !---------------------------------
        case ('MTDBT') ! TEMPERATURE/DRY BULB TEMPERATURE    [K]
                       !---------------------------------
          call set_datum (btemp(nlev)% t ,rval ,spt% corme)
                       !----------------------
        case ('MTDNH') ! DEW-POINT TEMPERATURE               [K]
                       !----------------------
          call set_datum (btemp(nlev)% td ,rval ,spt% corme)
                       !------------------------------------
        case ('NSASA') ! TRACKING TECHNIQUE,STATUS OF SYSTEM [CODE_TABLE]
                       !------------------------------------
          spt% tracking = ival
                       !---------------------------------
        case ('NA4')   ! TYPE OF MEASURING EQUIPMENT USED    [CODE_TABLE]
                       !---------------------------------
          spt% meas_type = ival
                         !---------------------------
        case ('Loop001') ! Start of Loop (wind shear)
                         !---------------------------
          exit
#ifdef  CHECKCODES
        !--------------------------------
        ! ignore for unused codes (PILOT)
        !--------------------------------
        case ('MOBITQ')  ! OVERALL QUALITY BITS                    [CODE_TABLE]
        case ('MADDF')   ! ASSOCIATED FIELD SIGNIFICANCE           [CODE_TABLE]
        case ('MA1')     ! WMO REGION NUMBER                       [CODE_TABLE]
        case ('MHHA')    ! HEIGHT OF LAND SURFACE                  [M]
        case ('MADDF0')  ! ASSOCIATED FIELD SIGNIFICANCE           [CODE_TABLE]
        case ('Loop000') ! Start of Loop
        case ('Lcnt000') ! Loop Counter
        case ('MDREP')   ! DELAYED DESCRIPTOR REPLICATION FACTOR   [NUMERIC]
        case ('NVBVB')   ! ABSOLUT WIND SHEAR IN 1 KM LAYER BELOW  [M/S]
        case ('NVAVA')   ! ABSOLUT WIND SHEAR IN 1 KM LAYER ABOVE  [M/S]
        case ('NFNFNQ')  ! Q-BITS FOR FOLLOWING VALUE              [CODE TABLE]
        case ('NDNDNQ')  ! Q-BITS FOR FOLLOWING VALUE              [CODE TABLE]
        !-------------------------------
        ! ignore for unused codes (TEMP)
        !-------------------------------
        case ('MGGACT')  ! ACTUAL HOUR OF OBSERVATION              [HOUR]
        case ('NGGACT')  ! ACTUAL MINUTE OF OBSERVATION            [MINUTE]
        case ('NSR')     ! SOLAR AND INFRARED RADIATION CORRECTION [CODE_TABLE]
        case ('MTNQ')    ! Q-BITS FOR FOLLOWING VALUE              [CODE TABLE]
        case ('MTDNQ')   ! Q-BITS FOR FOLLOWING VALUE              [CODE TABLE]
        case ('MHNQ')    ! Q-BITS FOR FOLLOWING VALUE              [CODE TABLE]
        !--------------------------------------
        ! ignore for unused codes (PILOT BUFR4)
        !--------------------------------------
        case ('MTISI')   ! TIME SIGNIFICANCE                       [CODE_TABLE]
        case ('MSEC')    ! SECOND                                  [SECOND]
        case ('MHOSNN')  ! HEIGHT OF STATION GROUND ABOVE MEAN SEA [M]
        case ('MHOBNN')  ! HEIGHT OF BAROMETER ABOVE MEAN SEA LEVEL[M]
        case ('MEDRE')   ! EXTENDED DELAYED DESCRIPT.REPLIC.FACTOR [NUMERIC]
!       case ('Loop001') ! START OF LOOP
        case ('Lcnt001') ! Loop Counter                            []
        case ('MEVSS0')  ! EXTENDED VERTICAL SOUNDING SIGNIFICANCE [FLAG_TABLE]
        case ('Loop002') ! Start of Loop                           []
        case ('MDREP0')  ! DELAYED DESCRIPTOR REPLICATION FACTOR   [NUMERIC]
        case ('NHHH0')   ! GEOPOTENTIAL HEIGHT                     [GPM]
        case ('Lcnt002') ! Loop Counter                            []
        case ('YSUPL')   ! CHARACTERS                              [CCITT IA5]
        !-------------------------------------
        ! ignore for unused codes (TEMP BUFR4)
        !-------------------------------------
        case ('MVTSU') ! VERTICAL SIGNIFICANCE (SURFACE OBSERV.)   [CODE_TABLE]
        case ('MNH')   ! CLOUD AMOUNT                              [CODE_TABLE]
        case ('NH')    ! HEIGHT OF BASE OF CLOUD                   [M]
        case ('MCC')   ! CLOUD TYPE                                [CODE_TABLE]
        case ('MCC0')  ! CLOUD TYPE                                [CODE_TABLE]
        case ('MCC1')  ! CLOUD TYPE                                [CODE_TABLE]
        case ('MPN0')  ! PRESSURE (VERT.LOCATION)                  [PA]
        case ('MTN00') ! SEA/WATER TEMPERATURE                     [K]
        !--------------------
        ! print unknown codes
        !--------------------
        case default
          call bufr_get_entry_texts (bufr)
          call bufr_get_entry_units (bufr)
          write(0,*) bufr% itype(id),bufr% ifxy(id)
          write(0,'(a,a,a)') ymnem, bufr% ytext(id),      &
                          '['//trim(bufr% yunit(id))//']'
          unused = .true.
#endif
        end select
      endif
    end do
    !--------------------------------------
    ! in case of unknown codes:
    ! print suspicious BUFR record and exit
    !--------------------------------------
    if (unused) then
      call bufr_print_sections (bufr)
      call bufr_print_subset   (bufr, is)
      call finish ('read_temp_bufr','code(s) not implemented')
    endif

    !-----------------------------------------------
    ! if station name is not set, use station number
    !-----------------------------------------------
    if (spt% statid =='') write(spt% statid,'(i5.5)') spt% ident
    where (iand (btemp% vss, lev_miss) /=0 ) btemp% vss = lev_miss

    !----------------
    ! standard checks
    !----------------
    lkeep = .true.
    call init_time (spt% actual_time, yyyy, mo, dd, hh, mi)
    call check_report_1 (spt)
    if (spt% use% state <= STAT_DISMISS) lkeep = .false.

    if (lkeep) then

      !--------------------------------------------
      ! if not present, derive pressure from height
      ! (US Standard Atmosphere)
      !--------------------------------------------
      if (rept_use(spt%hd% obstype)% deriv_p /= 0) then
        do i=1,size(btemp)
          if (btemp(i)% p % qc /= QC_OK .and. btemp(i)% gp% qc == QC_OK .and. &
              btemp(i)% gp% o  >  0.                                  ) then
             btemp(i)% p % o   = p_h_usstd (real(btemp(i)% gp% o,wp))
             btemp(i)% p % qc  = QC_OK
             btemp(i)% p % src = SRC_DER
             btemp(i)% gp% qc  = QC_NOUSE
          endif
        end do
      endif

      !-------------------------
      ! check for double entries
      !-------------------------
      if (corme_conv >= 2) then
        call correct_report (obs, spt, btemp, lkeep, obs% n_spot)
      else
        call correct_merge  (obs, spt, btemp, lkeep, obs% n_spot)
      endif

    endif

    if (lkeep) then
      call check_store_temp (btemp, spt, obs, novalid)
      lkeep = novalid > 0
    endif

   if (prt_data) then
      print *
      print *, spt% statid
      print *, spt% col% c% dlat, spt% col% c% dlon
      print *, spt% hd% buf_type, spt% hd% buf_subtype
      call print (spt% hd% time)
      call print (spt%     actual_time)
      call print (spt% hd% db_time)
      print *, spt% hd% dbkz
      call print_temp_data (btemp)
   endif

    deallocate (btemp)

  end subroutine read_temp_bufr
!==============================================================================
  subroutine reduce_profiles (obs)
    !----------------------------------------------------
    ! Reduce profile information:
    ! - derive layer averages from full TEMP profiles
    !   or interpolated significant levels (TODO)
    ! - thin excessive data from high-resolution profiles
    !----------------------------------------------------
    type(t_obs) ,intent(inout) :: obs   ! observations data type

    integer               :: is       ! report index
    integer               :: i0,il,in ! level  index bounds
    integer               :: j, k     ! indices
    integer               :: kmax, n  ! level count estimates
    real(wp), allocatable :: paux(:)  ! auxiliary grid for superobbing
    real(wp), allocatable :: ptmp(:)  ! temporary

    if (.not. (layer_avg .or. thin_bufr(1) > 0)) return

    if (layer_avg) then
       !----------------------------------------
       ! prepare for layer averaging/superobbing
       !----------------------------------------
       n = count (avl_levs > 0._wp)
       if (any (avl_levs(1:n-1) <= avl_levs(2:n)) .or.       &
           any (avl_levs(1:n)   <= 0._wp        ) .or.       &
           any (avl_incr(1:n-1) <= 0._wp        ) .or. n <= 1) then
          call finish ("reduce_profiles","bad avl_levs or avl_incr")
       end if
       kmax = nint (sum ((avl_levs(1:n-1)-avl_levs(2:n))/avl_incr(1:n-1))) + n
       allocate (ptmp(kmax))
       ptmp    = 0._wp
       ptmp(1) = avl_levs(1)                      ! first avg.layer bound
       j       = 1
       do k = 2, kmax
          if (avl_incr(j) <= 0._wp) then          ! no increments left
             kmax = k
             exit
          end if
          ptmp(k) = ptmp(k-1) - avl_incr(j)       ! next layer bound
          if (ptmp(k) <= 0._wp) then              ! p>0 for ln(p) to be real
             ptmp(k) = 0._wp
             kmax    = k - 1
             exit
          end if
          if (ptmp(k) <= avl_levs(j+1)) j = j + 1 ! take next increment
       end do
       allocate (paux(kmax))
       paux(1:kmax) = ptmp(1:kmax)
       deallocate (ptmp)
    end if

    !------------------
    ! loop over reports
    !------------------
    do is = 1, obs% n_spot
      if (obs% spot(is)% hd% modtype == COSMO       ) cycle
      if (obs% spot(is)% use% state  <= STAT_DISMISS) cycle
      select case (obs% spot(is)% hd% obstype)
      case (OT_TEMP, OT_PILOT)
        i0 = obs% spot(is)% o% i
        in = i0 + obs% spot(is)% o% n
        il = i0 + 1
        if (in == i0)  cycle
        !----------------------------------------
        ! superobbing of high-resolution profiles
        !----------------------------------------
        if (layer_avg) then
           !-------------------------------------
           ! currently implemented only for TEMPs
           !-------------------------------------
           select case (obs% spot(is)% hd% obstype)
           case (OT_TEMP)
              call average_layers (obs, obs% spot(is), paux, is)
           end select
        end if
        !-------------------------------------------------------
        ! thinning of high-resolution data with time information
        !-------------------------------------------------------
        if (thin_bufr(1) > 0) then
           call thin_profile (obs, obs% spot(is), is)
        end if
      end select
    end do
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  contains
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    subroutine average_layers (obs, spot, pref, is)
      type(t_obs)  ,intent(inout) :: obs      ! data of all observations
      type(t_spot) ,intent(inout) :: spot     ! meta data of this observation
      real(wp)     ,intent(in)    :: pref(:)  ! interpolation reference levels
      integer      ,intent(in)    :: is       ! report index
      !----------------
      ! Local variables
      !----------------
      integer     :: i                   ! level index
      integer     :: i0,il,in            ! level index bounds
      integer     :: j, jj,k             ! indices
      integer     :: k0, ktop            ! level indices
      integer     :: n                   ! No. interpolation levels
      integer     :: nlev                ! No. observed levels
      integer     :: n_new               ! No. combined levels
      real(wp)    :: d                   ! pressure difference
      real(wp)    :: ptop                ! Top of profile
      real(wp)    :: psur, pbot          ! Surface pressure, level above surf.
      real(wp)    :: p1, p2              ! Temporaries
      real(wp)    :: ptol                ! Tolerance for difference
      real(wp)    :: w1, w2              ! Interpolation weights
      real(wp)    :: w, sumw             ! Integration weights
      real(wp)    :: lat1, lat2          ! Latitude  range
      real(wp)    :: lon1, lon2          ! Longitude range
      real(wp), allocatable :: pobs(:)   ! Observed pressure levels
      real(wp), allocatable :: pint(:)   ! pressure grid for interpolation
      real(wp), allocatable :: lnpint(:) ! Vertical  coord. for interpolation
      real(wp), allocatable :: lnpobs(:) ! Reference coord. for interpolation
      type(t_temp) ,pointer :: temp(:)   ! TEMP data
      type(t_temp) ,pointer :: tint(:)   ! Interpolated values
      type(t_temp) ,pointer :: t_av(:)   ! Layer averages
      type(t_temp) ,pointer :: tnew(:)   ! Resulting TEMP profile
      type(t_temp)          :: ttmp      ! Temporary
      real(wp), parameter   :: tol1 = 0.5_wp    ! Adjustment top/bottom level
      real(wp), parameter   :: tol2 = 0.1_wp    ! Adjustment other levels

      if (.not. layer_avg)  return
      nlev = spot%col% nlev
      select case (spot% hd% codetype)
      case default
         if (nlev < supob_thr(1)) return
      case (230)  ! TEMP DROP BUFR: deferred until below
      end select
      !-----------------------------------------
      ! Check drift for excessive/missing values
      !-----------------------------------------
      i0   = spot% o% i
      in   = i0 + spot% o% n
      il   = i0 + 1
      lat1 = minval (obs% body(il:in)% lat)
      lat2 = maxval (obs% body(il:in)% lat)
      lon1 = minval (obs% body(il:in)% lon)
      lon2 = maxval (obs% body(il:in)% lon)
      if (abs (lat1-lat2) >  10._wp .or. &
          abs (lon1-lon2) > 400._wp .or. &
                     lat2 == rvind  .or. &
                     lon2 == rvind       ) return
      !--------------------------------------
      ! Load TEMP specific observational data
      !--------------------------------------
      call load_temp (obs, spot, temp)
      allocate (pobs(nlev))
      pobs(:) = temp(:)% p% o
      !--------------------------------------------------------
      ! TEMP DROP BUFR: use threshold depending on top pressure
      !--------------------------------------------------------
      if (spot% hd% codetype == 230) then
         if (nlev < supob_thr(2) * (1._wp - minval (pobs) / 101325._wp)) then
            deallocate (pobs, temp)
            return
         end if
      end if
      !-------------------------
      ! Search for surface level
      !-------------------------
      psur = 999999._wp
      do i = il, in
         if (btest (obs% body(i)% lev_sig, LS_SURFACE)) then
            psur = obs% olev(i)
            exit
         end if
      end do
      !-----------------------------------------------
      ! Search for lowest reported level above surface
      ! which has not just extrapolated geopotential.
      ! Account for minimal distance to surface
      !-----------------------------------------------
      do n = 1, nlev
         pbot = pobs(n)
         if (temp(n)% t% qc == QC_OK .or. temp(n)% ff% qc == QC_OK) exit
      end do
      if (psur < 999999._wp) then
         do i = il, in
            if (.not. btest (obs% body(i)% lev_sig, LS_SURFACE) &
                .and.        obs% olev(i) < psur                ) then
               pbot = obs% olev(i)
               exit
            end if
         end do
      end if
      pbot = min (pbot, psur - dp_asurf)
      !--------------------------------------------
      ! Determine top level for layer calculations:
      ! must provide at least valid T or FF,
      ! geopotential may have been extrapolated!
      !--------------------------------------------
      ptop = 999999._wp
      do n = nlev, 1, -1
         if (temp(n)% t% qc == QC_OK .or. temp(n)% ff% qc == QC_OK) then
            ptop = pobs(n)
            exit
         end if
      end do
      if (ptop == 999999._wp) then
         deallocate (temp)
         return
      end if
      !-------------------------------------------------------------------
      ! Determine number of reference levels to use for layer calculations
      !-------------------------------------------------------------------
      ktop = count (pref >= ptop)
      if (ktop < size (pref)) then
         ptol = (pref(ktop) - pref(ktop+1)) * tol1
         if (abs (ptop - pref(ktop+1)) < ptol) ktop = ktop + 1
      end if
      k0 = minloc (abs (pref - pbot), dim=1)    ! Lowest p-ref level
      n  = ktop - k0 + 1
      if (prt_data) then
         print *, "average_layer: ", trim (spot% statid)
         print *, "pbot,ptop,ok,n(ref)", real (pbot), real (ptop), pbot>ptop, n
      end if
      !---------------------------------------------
      ! We need at least n=2 levels (for n-1 layers)
      !---------------------------------------------
      if (n <= 1) then
         deallocate (temp)
         return
      end if
      allocate (pint(n))
      pint(1:n) = pref(k0:k0+n-1)
      !----------------------------------------------------------
      ! Allow adjustment of interpolation levels within tolerance
      ! in order to avoid excessive interpolation
      !----------------------------------------------------------
      if (     pint(1) > pbot        .or.           &
          abs (pint(1) - pbot) < tol1 * avl_incr(1) ) then
!print *, "pint(1) -> pbot", pint(1), pbot
         pint(1) = pbot
      end if
      if (     pint(n) < ptop                      .or.     &
          abs (pint(n) - ptop) < tol1 * (pint(n-1)-pint(n)) ) then
!print *, "pint(n) -> ptop", pint(n), ptop
         pint(n) = ptop
      end if
      do k = 1, n
         if (k < n) then
            ptol = (pint(k) - pint(k+1)) * tol2
         else
            ptol = (pint(k-1) - pint(k)) * tol2
         end if
         i = minloc (abs (pint(k) - pobs(:)), dim=1)
         d =              pint(k) - pobs(i)
         if (d /= 0._wp .and. abs (d) < ptol .and. (k > 1 .or. d > 0._wp)) then
            i = minloc (abs (pint(k) - pobs(:)), dim=1) - lbound (pobs,1) + 1
!print *, "pint(k) -> pobs(i)", pint(k), pobs(i)
            pint(k) = pobs(i)
         end if
      end do
      !----------------------------
      ! Derive interpolation levels
      !----------------------------
      allocate (lnpobs(nlev))
      allocate (lnpint(n))
      allocate (tint  (n))
      lnpobs(:) = log (pobs(:))
      lnpint(:) = log (pint(:))
      do k = 1, n
         i = minloc (abs (pint(k) - pobs(:)), dim=1)
         d =              pint(k) - pobs(i)
         !---------------------------------------
         ! Copy data, metadata from nearest level
         !---------------------------------------
         tint(k) = temp(i)
         if (abs (d) > 0.01_wp) then
!print *, "### k,pint,d=", k, pint(k), d
            if (d > 0._wp) then
               j = i-1                  ! Interpolate from below
            else
               j = i+1                  ! Interpolate from above
            end if
            w1 = (lnpint(k) - lnpobs(j)) / (lnpobs(i) - lnpobs(j))
            w2 = 1._wp - w1
            ! Pressure
            if (temp(i)% p% qc == QC_OK .and. &
                temp(j)% p% qc == QC_OK       ) then
               tint(k)% p% o   = exp (w1 * lnpobs(i) + w2 * lnpobs(j))
               tint(k)% p% src = SRC_DER
            else
               tint(k)% p% qc  = QC_NOUSE
               write(0,*) "statid,i,j,qc(i),qc(j)=", spot% statid, &
                    i, j, temp(i)% p% qc,temp(j)% p% qc
               call finish ("average_layers", "Should not happen: p%qc /= 0")
            end if
            ! Temperature
            if (temp(i)% t% qc == QC_OK .and. &
                temp(j)% t% qc == QC_OK       ) then
               tint(k)% t% o   = w1 * temp(i)% t% o + w2 * temp(j)% t% o
               tint(k)% t% src = SRC_DER
            else
               tint(k)% t% qc  = QC_NOUSE
            end if
            ! Rel. humidity
            if (temp(i)% rh% qc == QC_OK .and. &
                temp(j)% rh% qc == QC_OK       ) then
               tint(k)% rh% o   = w1 * temp(i)% rh% o + w2 * temp(j)% rh% o
               tint(k)% rh% src = SRC_DER
            else
               tint(k)% rh% qc  = QC_NOUSE
            end if
            ! Wind
            if (temp(i)% uu% qc == QC_OK .and. &
                temp(j)% uu% qc == QC_OK       ) then
               tint(k)% uu% o   = w1 * temp(i)% uu% o + w2 * temp(j)% uu% o
               tint(k)% vv% o   = w1 * temp(i)% vv% o + w2 * temp(j)% vv% o
               tint(k)% uu% src = SRC_DER
               tint(k)% vv% src = SRC_DER
            else
               tint(k)% uu% qc  = QC_NOUSE
               tint(k)% vv% qc  = QC_NOUSE
            end if
            ! Geopotential
            if (temp(i)% gp% qc == QC_OK .and. &
                temp(j)% gp% qc == QC_OK       ) then
               tint(k)% gp% o   = w1 * temp(i)% gp% o + w2 * temp(j)% gp% o
               tint(k)% gp% src = SRC_DER
            else
               tint(k)% gp% qc  = QC_NOUSE
            end if
            ! Coordinates and time offset (approximated)
            tint(k)% dlat = real (w1 * temp(i)% dlat + w2 * temp(j)% dlat)
            tint(k)% dlon = real (w1 * temp(i)% dlon + w2 * temp(j)% dlon)
            tint(k)% secs = nint (w1 * temp(i)% secs + w2 * temp(j)% secs)
!print*,"++",pint(k),i,j, temp(i)% p% o, temp(j)% p% o
         end if
      end do
      !-------------------------------------------------------------
      ! Vertical coordinate for averaging: prefer height over log(p)
      !-------------------------------------------------------------
!     i0 = minloc (abs (temp(:)% p% o - psur), dim=1)
!     if (all (temp(i0:)% gp% qc == QC_OK) .and. &
!         all (tint(  :)% gp% qc == QC_OK)       ) then
!        lnpobs(:) = temp(:)% gp% o
!        lnpint(:) = tint(:)% gp% o
!     end if
      !----------------------------
      ! Derive (n-1) layer averages
      !----------------------------
      allocate (t_av(n))
      t_av(n)% p% o = 0._sp     ! Auxiliary layer above top
      do k = 1, n-1
         p1 = pint(k)
         p2 = pint(k+1)
         i  = maxloc (pobs, mask=(pobs <= p1), dim=1)
         j  = minloc (pobs, mask=(pobs >= p2), dim=1)
         call init (ttmp, sumw)
         if (i <= j) then
!print*, "::", k,i,j,p1,p2
            w1    = abs (lnpint(k) - lnpobs(i))
            if (w1 > 0._wp) then
               call add (ttmp, tint(k), w1, sumw)
               call add (ttmp, temp(i), w1, sumw)
            end if
            do jj = i, j-1
               w  = abs (lnpobs(jj) - lnpobs(jj+1))
!print*, ":: w=", w
               call add (ttmp, temp(jj),   w, sumw)
               call add (ttmp, temp(jj+1), w, sumw)
            end do
            w2   = abs (lnpobs(j) - lnpint(k+1))
            if (w2 > 0._wp) then
               call add (ttmp, temp(j),   w2, sumw)
               call add (ttmp, tint(k+1), w2, sumw)
            end if
            ttmp% secs = (temp(i)% secs + temp(j)% secs) / 2
         else
            w = abs (lnpint(k) - lnpint(k+1))
!print*, "::intonly:", k,i,j,p1,p2,w
            call add (ttmp, tint(k),   w, sumw)
            call add (ttmp, tint(k+1), w, sumw)
            ttmp% secs = (tint(k)% secs + tint(k+1)% secs) / 2
         end if
         if (sumw > 0._wp) then
            ! Normalize
            call normalize (ttmp, sumw)
            call set_ff_dd (ttmp)
            t_av(k) = ttmp
!print *, "av:", ttmp%p%o, ttmp%t%o, ttmp%rh%o*100, ttmp%ff%o, ttmp%dd%o, ttmp%gp%o
         end if
      end do
      deallocate (tint)
      deallocate (pint, lnpobs, lnpint)
      !-------------------------------------------
      ! Merge observed profile with layer averages
      !-------------------------------------------
      allocate (tnew(nlev+n-1))
      i = 1
      j = 1
      do k = 1, nlev + n-1
         if (i <= nlev .and. (j >= n .or. temp(i)% p% o >= t_av(j)% p% o)) then
            tnew(k) = temp(i)
            i = i + 1           !min (i + 1, nlev)
         else if (j < n) then
            tnew(k) = t_av(j)
            j = j + 1           !min (j + 1, n-1)
         else
            write(0,*) k,i,j,temp(i)% p% o, t_av(j)% p% o, temp(i)% p% o == t_av(j)% p% o
            call finish ("average_layers", "cannot happen!")
         end if
      end do
      deallocate (temp, t_av)
      call check_store_temp (tnew, spot, obs, n_new, repl=is)
      deallocate (tnew)
    end subroutine average_layers
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    subroutine init (t, sumw)
      type(t_temp), intent(out) :: t
      real(wp),     intent(out) :: sumw

      t% p % o   = 0._sp
      t% t % o   = 0._sp
      t% rh% o   = 0._sp
      t% uu% o   = 0._sp
      t% vv% o   = 0._sp
      t% gp% o   = 0._sp
      t% p % qc  = QC_OK
      t% t % qc  = QC_OK
      t% rh% qc  = QC_OK
      t% uu% qc  = QC_OK
      t% vv% qc  = QC_OK
      t% gp% qc  = QC_OK
      t% p % src = SRC_DER
      t% t % src = SRC_DER
      t% rh% src = SRC_DER
      t% uu% src = SRC_DER
      t% vv% src = SRC_DER
      t% gp% src = SRC_DER
      t% p % lev_sig = ibset (0_i2, LS_SUPEROBS)
      t% t % lev_sig = ibset (0_i2, LS_SUPEROBS)
      t% rh% lev_sig = ibset (0_i2, LS_SUPEROBS)
      t% uu% lev_sig = ibset (0_i2, LS_SUPEROBS)
      t% vv% lev_sig = ibset (0_i2, LS_SUPEROBS)
      t% gp% lev_sig = ibset (0_i2, LS_SUPEROBS)
      t% dlat    = 0._sp
      t% dlon    = 0._sp
      t% secs    = 0
      sumw       = 0._wp
    end subroutine init
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    subroutine add (t1, t2, w, sumw)
      type(t_temp), intent(inout) :: t1
      type(t_temp), intent(in)    :: t2
      real(wp),     intent(in)    :: w
      real(wp),     intent(inout) :: sumw

      sumw        = sumw + w
      t1% p % qc  = max (t1% p % qc, t2% p % qc)
      t1% t % qc  = max (t1% t % qc, t2% t % qc)
      t1% rh% qc  = max (t1% rh% qc, t2% rh% qc)
      t1% uu% qc  = max (t1% uu% qc, t2% uu% qc)
      t1% vv% qc  = max (t1% vv% qc, t2% vv% qc)
      t1% gp% qc  = max (t1% gp% qc, t2% gp% qc)
      if (t1% p % qc == 0) t1% p % o = t1% p % o + w * log(real(t2% p % o,kind=wp))
      if (t1% t % qc == 0) t1% t % o = t1% t % o + w *      t2% t % o
      if (t1% rh% qc == 0) t1% rh% o = t1% rh% o + w *      t2% rh% o
      if (t1% uu% qc == 0) t1% uu% o = t1% uu% o + w *      t2% uu% o
      if (t1% vv% qc == 0) t1% vv% o = t1% vv% o + w *      t2% vv% o
      if (t1% gp% qc == 0) t1% gp% o = t1% gp% o + w *      t2% gp% o
      t1% dlat                       = t1% dlat  + w *      t2% dlat
      t1% dlon                       = t1% dlon  + w *      t2% dlon
    end subroutine add
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    subroutine normalize (t, sumw)
      type(t_temp), intent(inout) :: t
      real(wp),     intent(in)    :: sumw

      real(wp) :: nrm
      nrm = 1 / sumw
      if (t% p % qc == 0) t% p % o = exp (t% p % o * nrm)
      if (t% t % qc == 0) t% t % o =      t% t % o * nrm
      if (t% rh% qc == 0) t% rh% o =      t% rh% o * nrm
      if (t% uu% qc == 0) t% uu% o =      t% uu% o * nrm
      if (t% vv% qc == 0) t% vv% o =      t% vv% o * nrm
      if (t% gp% qc == 0) t% gp% o =      t% gp% o * nrm
      t% dlat                      =      t% dlat  * nrm
      t% dlon                      =      t% dlon  * nrm
    end subroutine normalize
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    subroutine set_ff_dd (t)
      type(t_temp), intent(inout) :: t
      real(wp) :: ff, dd
      if (t% uu% qc == 0 .and. t% vv% qc == 0) then
         call fd_uv (ff, dd, real (t% uu% o,wp), real (t% vv% o, wp))
         t% ff% o       = ff
         t% dd% o       = dd
         t% ff% qc      = QC_OK
         t% dd% qc      = QC_OK
         t% ff% src     = SRC_DER
         t% dd% src     = SRC_DER
         t% ff% lev_sig = t% uu% lev_sig
         t% dd% lev_sig = t% uu% lev_sig
      end if
    end subroutine set_ff_dd
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    subroutine thin_profile (obs, spot, is)
      type(t_obs)  ,intent(inout) :: obs      ! data of all observations
      type(t_spot) ,intent(inout) :: spot     ! meta data of this observation
      integer      ,intent(in)    :: is       ! report index

      integer               :: i0,il,in ! level index bounds
      integer               :: i,im,k   ! level indices
      integer               :: t,t0,t1  ! time "levels"
      integer               :: nlev     ! old no. levels
      integer               :: nnew     ! new no. levels
      integer               :: nval     ! no. valid observations
      integer               :: ti       ! thinning interval index
      logical               :: ldrop    ! dropsonde?
      logical,  allocatable :: keep(:)  ! keep level?
      type(t_temp) ,pointer :: temp(:)  ! TEMP specific data type
      type(t_temp) ,pointer :: tnew(:)  ! TEMP specific data type
      !------------------------------------------------------
      ! Mask for all "classical" level significances to keep,
      ! including superobservations
      !------------------------------------------------------
      integer(i2),parameter :: mask_sig = ibset (0_i2, LS_SURFACE)  &!
                                        + ibset (0_i2, LS_STANDARD) &!
                                        + ibset (0_i2, LS_SIGN)     &!
                                        + ibset (0_i2, LS_TROPO)    &!
                                        + ibset (0_i2, LS_MAX)      &!
                                        + ibset (0_i2, LS_SUPEROBS)  !

      nlev = spot% col% nlev
      if (nlev <= 1) return

      i0 = spot% o% i
      in = i0 + spot% o% n
      il = i0 + 1
      if (all (obs% body(il:in)% secs == 0)) return

      call load_temp (obs, spot, temp)
      allocate (keep(nlev))
      keep(:) = .false.
      !--------------------------------------------------
      ! Dropsonde detection, taking into account that the
      ! lowest and highest reported levels may have been
      ! extrapolated, and having no valid time
      !--------------------------------------------------
      do i = nlev, 1, -1
         t1 = temp(i)% secs
         if (t1 /= 0) exit
      end do
      do i = 1, nlev
         t0 = temp(i)% secs
         if (t0 /= 0) exit
      end do
      i0    = i
      ldrop = t0 > t1
      if (ldrop) temp(:)% secs = - temp(:)% secs    ! Reverse times
      if (ldrop .and. temp(1)% secs == 0) then
         i = i0                                     ! Skip extrapolated levels
      else
         i = 1                                      ! Thin from beginning
      end if
      t0 = temp(i)% secs
      ti = 1
      do
         !------------------------------
         ! Switch thinning above 100 hPa
         !------------------------------
         if (temp(i)% p% o <= 10000._sp) ti = 2
         !------------------
         ! Next time "level"
         !------------------
         t1   = t0 + thin_bufr(ti)
         !---------------------------------------------
         ! Search for next level with time offset >= t1
         !---------------------------------------------
         im   = nlev
         do k = i+1, nlev
            if (temp(k)% secs >= t1) then
               im = k-1
               exit
            end if
         end do
         !------------------------------------
         ! Apply thinning, keep all marked obs
         !------------------------------------
         t    = temp(i)% secs
         do k = i, im
            keep(k) = temp(k)% secs                   <= t    .or. &
                      temp(k)% secs                   >= t1   .or. &
                iand (temp(k)% p % lev_sig, mask_sig) /= 0_i2 .or. &
                iand (temp(k)% t % lev_sig, mask_sig) /= 0_i2 .or. &
                iand (temp(k)% rh% lev_sig, mask_sig) /= 0_i2 .or. &
                iand (temp(k)% uu% lev_sig, mask_sig) /= 0_i2 .or. &
                iand (temp(k)% gp% lev_sig, mask_sig) /= 0_i2
         end do
         i = im + 1
         if (i > nlev) exit
         t0 = t1
      end do
      if (ldrop) temp(:)% secs = - temp(:)% secs        ! Restore times

      nnew = count (keep)
      if (nnew > 0) then
         !----------------------
         ! Store selected levels
         !----------------------
         allocate (tnew(nnew))
         tnew(1:nnew) = pack (temp, mask=keep)
         call check_store_temp (tnew, spot, obs, nval, repl=is)
         deallocate (tnew)
      else
         call decr_rpt_use (spot, CHK_THIN)
      end if
      deallocate (temp)

    end subroutine thin_profile
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine reduce_profiles
!==============================================================================
  subroutine select_levels (obs)
  !--------------------------------------------------
  ! select the levels to be assimilated
  ! based on level significance and namelist settings
  ! disabled for COSMO multilevel data
  !--------------------------------------------------
  type (t_obs) ,intent(inout) :: obs   ! observations data type

    integer     :: is       ! report index
    integer     :: i        ! level  index
    integer     :: i0,il,in ! level  index bounds
    real(wp)    :: olev     ! pressure level of observation
    integer(i2) :: lev_sig  ! level significance to test
    integer(i2) :: act_t, act_q, act_v, act_h, act_v_wp, act_v_
    real(wp) ,pointer :: vmpl (:)
    logical     :: lavg     ! prefer layer averages
    logical     :: junior   ! lower-ranking ("junior") level?

    !----------------------------------------------------
    ! derive Feedback-File vertical sounding significance
    !   from               vertical sounding significance
    !----------------------------------------------------
    call active_levels (act_t    ,active_t    ,layer_avg)
    call active_levels (act_q    ,active_q    ,layer_avg)
    call active_levels (act_v_   ,active_v    ,layer_avg)
    call active_levels (act_v_wp ,active_v_wp ,layer_avg)
    !-------------------------------------------------
    ! Beware to not use layer average for geopotential
    !-------------------------------------------------
    call active_levels (act_h, active_h)

    !------------------
    ! loop over reports
    !------------------
    do is = 1, obs% n_spot
      if (obs% spot(is)% hd% modtype == COSMO       ) cycle
      if (obs% spot(is)% use% state  <= STAT_DISMISS) cycle
      select case (obs% spot(is)% hd% obstype)
      case (OT_TEMP, OT_PILOT)
        i0 = obs% spot(is)% o% i
        in = i0 + obs% spot(is)% o% n
        il = i0 + 1
        if (in == i0)  cycle
        !----------------------------------------------------------
        ! list of main pressure levels in addition to those in BUFR
        !----------------------------------------------------------
        vmpl  => vmainpl
        act_v =  act_v_
        !----------------------------------
        ! different list for wind profilers
        !----------------------------------
        if (    obs% spot(is)% hd% obstype  == OT_PILOT .and. &
            any(obs% spot(is)% hd% codetype == ct_wprof)    ) then
          vmpl  => wmainpl
          act_v =  act_v_wp
        endif
        call mark_p_levels ()
        call mark_z_levels ()
        !---------------------
        ! select active levels
        !---------------------
        if (layer_avg .and. prefer == 2) then
           lavg = any (btest (obs% body(il:in)% lev_sig, LS_SUPEROBS))
        else
           lavg = .false.
        end if
        olev = 999999._wp
        do i = obs% spot(is)% o% i + obs% spot(is)% o% n, &
               obs% spot(is)% o% i + 1, -1
          !--------------------------------------------------------
          ! reset other level significance values for surface level
          ! do not use levels below surface (except geopotential)
          !--------------------------------------------------------
          if (btest (obs% body(i)% lev_sig, LS_SURFACE)) then
            olev =   obs% olev(i)
            obs% body(i)% lev_sig = ibset (0_i2, LS_SURFACE)
          endif
          if (obs% olev(i) > olev) then
            if (btest (obs% body(i)% lev_sig, LS_STANDARD) &
                .and.  obs% varno(i) == VN_Z               ) then
              call decr_use (obs% body(i)% use, STAT_PASSIVE, check=CHK_HEIGHT)
            else
              call decr_use (obs% body(i)% use,               check=CHK_HEIGHT)
            end if
          end if
          !----------------------------------------------
          ! select levels as requested by active_.. flags
          !----------------------------------------------
          lev_sig = obs% body(i)% lev_sig
          select case (lev_sig)
          case (:-1)
            lev_sig = LEV_MISSING
          case (0)
            lev_sig = ibset (0_i2, LS_OTHER)
          end select
          junior = lavg .neqv. btest (lev_sig, LS_SUPEROBS)

          select case (obs% varno(i))
          case (VN_T)
            if (iand (lev_sig, act_t) == 0 .or. junior)               &
              call decr_use (obs% body(i)% use, STAT_PASSIVE, CHK_THIN)
          case (VN_Z)
            if (iand (lev_sig, act_h) == 0)                           &
              call decr_use (obs% body(i)% use, STAT_PASSIVE, CHK_THIN)
          case (VN_RH)
            if (iand (lev_sig, act_q) == 0 .or. junior)               &
              call decr_use (obs% body(i)% use, STAT_PASSIVE, CHK_THIN)
          case (VN_U, VN_V )
            if (iand (lev_sig, act_v) == 0 .or. junior)               &
              call decr_use (obs% body(i)% use, STAT_PASSIVE, CHK_THIN)
          case (VN_FF,VN_DD)
            select case (obs% spot(is)% hd% codetype)
            case default    ! redundant with U,V
              call decr_use (obs% body(i)% use, STAT_PASSIVE, CHK_NOTUSED)
            case (OC_SCADA) ! SCADA (FF only, accept)
!             call decr_use (obs% body(i)% use, STAT_PASSIVE, CHK_NOTUSED)
            end select
          case (VN_NH,VN_N_L,VN_N_M,VN_RH2M,VN_T2M,VN_U10M,VN_V10M,VN_HEIGHT,VN_W)
            !-----------------------------------------------
            ! these quantities are currently not assimilated
            !-----------------------------------------------
            call decr_use (obs% body(i)% use, STAT_OBS_ONLY, CHK_NOIMPL)
          case (VN_PS)
            call decr_use (obs% body(i)% use, STAT_DISMISS,  CHK_NOIMPL)
          case default
            call finish ('select_levels',                                     &
                         'invalid variable '//name_value(varno, obs% varno(i)))
          end select
        end do
      end select

      !---------------------------------
      ! thinning of main pressure levels
      !---------------------------------
      !-------------------------------
      ! thinning of significant levels
      !-------------------------------
      !----------------------
      ! add artificial levels
      !----------------------
    end do
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  contains
    subroutine mark_p_levels ()
      !--------------------------------------------
      ! Pressure level selection for TEMP and PILOT
      !--------------------------------------------
      integer  :: i, im, ii(1)
      real(wp) :: olev     ! pressure level of observation
      real(wp) :: dlev     ! difference to mean pressure level
      real(wp) :: d        ! min. distance to main plev for 1st and last level

      if (any (obs% body(il:in)% lev_typ == VN_HEIGHT)) return
      !------------------------------------------------
      ! mark main pressure levels according to the list
      !------------------------------------------------
      do i=1,size(vmpl)
        if (vmpl(i)<=0._wp) cycle
        ii = minloc (abs(vmpl(i)-obs% olev(il:in)))
        im = i0 + ii(1)
        !-------------------------------------------------------
        ! do not mark surface levels
        ! nor first level if distance to main plev > 5%
        !-------------------------------------------------------
        if (btest (obs% body(im)% lev_sig, LS_SURFACE))      cycle
        olev = obs% olev(im)
        dlev = olev - vmpl(i)
        d    = 0.05_wp * olev
        if (olev == obs% olev(il) .and. dlev < -d)           cycle
        !---------------------------------------------------------------
        ! Is this level really the one closest to the current main plev?
        !---------------------------------------------------------------
        if (abs (dlev) > 0._wp) then
           if (i /= minloc (abs (vmpl(:) - olev), dim=1))    cycle
        end if
        !----------------------------------------------------------
        ! set level-significance for all observations at this level
        !----------------------------------------------------------
        do im = i0 + ii(1), in
          if (olev /= obs% olev(im)) exit
          obs% body(im)% lev_sig = ibclr (ibset (obs% body(im)% lev_sig,&
                                                 LS_STANDARD          ),&
                                                 LS_MISSING             )
        end do
      end do
    end subroutine mark_p_levels

    subroutine mark_z_levels ()
      integer :: i
      if (obs% spot(is)% hd% obstype     /= OT_PILOT)   return
      if (any (obs% body(il:in)% lev_typ /= VN_HEIGHT)) return
      !------------------------------------------------
      ! TODO: level selection for profilers.
      ! Set all levels to "standard" for now (testing).
      !------------------------------------------------
      do i = il, in
         obs% body(i)% lev_sig = ibclr (ibset (obs% body(i)% lev_sig,&
                                               LS_STANDARD         ),&
                                               LS_MISSING            )
      end do
    end subroutine mark_z_levels

    subroutine active_levels (active, levels, superobs)
    integer(i2)       ,intent(out) :: active
    integer           ,intent(in)  :: levels
    logical ,optional ,intent(in)  :: superobs
    !----------------------------------------------------
    ! derive Feedback-File vertical sounding significance
    !   from               vertical sounding significance
    !----------------------------------------------------
      active = 0_i2
      if (iand (levels, lev_other) /= 0) active = ibset (active, LS_OTHER)
      if (iand (levels, lev_surf ) /= 0) active = ibset (active, LS_SURFACE)
      if (iand (levels, lev_std  ) /= 0) active = ibset (active, LS_STANDARD)
      if (iand (levels, lev_tropo) /= 0) active = ibset (active, LS_TROPO)
      if (iand (levels, lev_max_v) /= 0) active = ibset (active, LS_MAX)
      if (iand (levels, lev_sig_t) /= 0) active = ibset (active, LS_SIGN)
      if (iand (levels, lev_sig_v) /= 0) active = ibset (active, LS_SIGN)
      if (iand (levels, lev_miss ) /= 0) active = ibset (active, LS_MISSING)
      if (present (superobs)) then
         if (superobs)                   active = ibset (active, LS_SUPEROBS)
      end if
    end subroutine active_levels
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine select_levels
!==============================================================================
  subroutine check_dbl_pilot (obs)
  type (t_obs)  ,intent(inout)        :: obs   ! observations data type
  !-------------------------------------
  ! check for collocated TEMP and PILOTs
  !-------------------------------------
    integer               :: i, j             !       spot index
    integer               :: it               ! TEMP  spot index
    integer               :: ip               ! PILOT spot index
    integer               :: itx(obs% n_spot) ! index array for TEMP  indices
    integer               :: ipx(obs% n_spot) ! index array for PILOT indices
    integer               :: nt               ! number of TEMP  reports
    integer               :: np               ! number of PILOT reports
    integer               :: i1t, int         ! TEMP level index bounds
    integer               :: i1p, inp         ! TEMP level index bounds
    integer               :: nwt, nwp         ! number of wind obsv.in TEMP/P.
    integer               :: nwr              ! number of levels removed
    integer               :: lt, lp, l        ! TEMP, PILOT level index
    real(wp)              :: dzmin, dz        ! pressure   difference (Pa)
    real(wp)              :: dx               ! horizontal difference (km)
    real(wp)              :: dv               ! wind vector difference
    type(t_spot) ,pointer :: sptt, sptp       ! pointer to meta data
    character(len=24)     :: action
    logical               :: rmlev
    integer               :: nffdd            ! # of monitoring vars (ff/dd)
    !--------------------------------------------
    ! preselect PILOTs (no wind profilers), TEMPs
    !--------------------------------------------
    np = 0
    nt = 0
    do ip = 1, obs% n_spot
      it = ip
      select case (obs% spot(ip)% hd% obstype)
      case (OT_PILOT)
        if (any(obs% spot(ip)% hd% codetype == ct_wprof)) cycle
        np = np + 1
        ipx(np) = ip
      case (OT_TEMP)
        nt = nt + 1
        itx(nt) = it
      end select
    end do
    if (np == 0) return
    if (nt == 0) return

    nffdd = count ([monitor_ff, monitor_dd])
    !----------------
    ! loop over TEMPs
    !----------------
    do j = 1, nt
      it = itx(j)
      if (obs% spot(it)% use% state   <= STAT_DISMISS) cycle
      sptt => obs% spot(it)
      i1t = sptt% o% i + 1
      int = sptt% o% i + sptt% o% n
      nwt = count (obs% varno(i1t:int)             == VN_U .and.  &
                   obs% body (i1t:int)% use% state >  STAT_DISMISS )
      if (nwt == 0) cycle
      !-----------------
      ! loop over PILOTs
      !-----------------
      do i = 1, np
        ip = ipx(i)
        !-------------------------------------
        ! check for collocated TEMP and PILOTs
        !-------------------------------------
        if (obs% spot(ip)% use% state   <= STAT_DISMISS)              cycle
        dx = (rearth/1000._wp) * &
             sqrt(sum( (sptt% col%c% x - obs% spot(ip)% col%c% x)**2))
        if (dx > dx_coloc)                                            cycle
        sptp => obs% spot(ip)
        i1p = sptp% o% i + 1
        inp = sptp% o% i + sptp% o% n
        nwp = count (obs% varno(i1p:inp)             == VN_U .and.  &
                     obs% body (i1p:inp)% use% state >  STAT_DISMISS )
        !-------------------------------------------------
        ! check for collocated observations level by level
        !-------------------------------------------------
        nwr = 0
        if (dz_coloc >= 0._wp) then
          if (rtl_coloc) then
            action =' remove levels in TEMP'
          else
            action =' remove levels in PILOT'
          endif
          !-----------------
          ! PILOT level loop
          !-----------------
          do lp = i1p, inp
            if (obs% varno(lp)             /= VN_U)         cycle
            if (obs% body (lp)% use% state <= STAT_DISMISS) cycle
            !----------------
            ! TEMP level loop
            !----------------
            dzmin = huge(dzmin)
            lt = 0
            do l = i1t, int
              if (obs% varno(l) /= VN_U) cycle
              dz = abs (obs% olev(lp) - obs% olev (l))
              if (dz < dzmin) then
                lt    = l
                dzmin = dz
                dv    = sqrt (  (obs% body(lp  )%o - obs% body(lt  )%o)**2 &
                              + (obs% body(lp+1)%o - obs% body(lt+1)%o)**2 )
              endif
            end do
            rmlev = ( dzmin <=  dz_coloc .or.              &
                     (dzmin <= dzv_coloc .and. dv == 0._wp))

            if (dzmin /= huge(dzmin)) then
              write (6,'(a,1x,2(a,2f9.3," - "),f8.3,2l2,2i5,3f12.3)') '#COLOC',&
              sptt% statid, sptt% col%c%dlat, sptt% col%c%dlon,                &
              sptp% statid, sptp% col%c%dlat, sptp% col%c%dlon,                &
              dzmin, rmlev, obs% body(lt)% use% state > STAT_DISMISS,          &
              lp-i1p,lt-i1t,obs% olev(lp), obs% olev (lt), dv
            endif

            if (rmlev) then
              !------------------------------------------------
              ! remove collocated observations from TEMP report
              !------------------------------------------------
              if (obs% body(lt  )% use% state > STAT_DISMISS) then
                if (rtl_coloc) then
                  call decr_use (obs% body(lt  )% use, STAT_DISMISS, CHK_REDUNDANT)
                  call decr_use (obs% body(lt+1)% use, STAT_DISMISS, CHK_REDUNDANT)
                  if (nffdd >= 1 .and. lt+2 <= int)                               &
                    call decr_use(obs%body(lt+2)% use, STAT_DISMISS, CHK_REDUNDANT)
                  if (nffdd == 2 .and. lt+3 <= int)                               &
                    call decr_use(obs%body(lt+3)% use, STAT_DISMISS, CHK_REDUNDANT)
                else
                  call decr_use (obs% body(lp  )% use, STAT_DISMISS, CHK_REDUNDANT)
                  call decr_use (obs% body(lp+1)% use, STAT_DISMISS, CHK_REDUNDANT)
                  if (nffdd >= 1 .and. lp+2 <= inp)                               &
                    call decr_use(obs%body(lp+2)% use, STAT_DISMISS, CHK_REDUNDANT)
                  if (nffdd == 2 .and. lp+3 <= inp)                               &
                    call decr_use(obs%body(lp+3)% use, STAT_DISMISS, CHK_REDUNDANT)
                endif
                nwr = nwr + 1
              endif
            endif

          end do
        else
          call decr_rpt_use (sptp ,CHK_REDUNDANT,comment='collocated TEMP/PILOT')
          action =' remove PILOT report'
        endif
        !---------
        ! printout
        !---------
        if (nwr == nwp .and. .not. rtl_coloc) then
          call decr_rpt_use (sptp ,CHK_REDUNDANT,comment='collocated TEMP/PILOT')
          action =' remove PILOT report'
        endif
        write (6,'(a,2i7,1x,2(a,2f9.3," - "),f8.3,3i4,a)') '#COLOC',&
          it,ip,                                                    &
          sptt% statid, sptt% col%c%dlat, sptt% col%c%dlon,         &
          sptp% statid, sptp% col%c%dlat, sptp% col%c%dlon,         &
          dx, nwt, nwp, nwr, action
      end do
    end do
  end subroutine check_dbl_pilot
!==============================================================================
  !----------------------------------------------------------------------
  ! Routine passed to BUFR reading routine
  ! insert a new TEMP or PILOT report into the observation data structure
  !----------------------------------------------------------------------
  subroutine check_store_temp (temp, spot, obs, novalid, repl, pass)
  type(t_temp) ,pointer              :: temp(:) ! Temp specific data type
  type(t_spot) ,intent(inout)        :: spot    ! report meta data
  type(t_obs)  ,intent(inout)        :: obs     ! observation data structure
  integer      ,intent(out)          :: novalid ! number of valid observations
  integer      ,intent(in) ,optional :: repl    ! observation to replace
  integer      ,intent(in) ,optional :: pass    ! checking pass

    type(t_temp) ,pointer :: t(:)
    type(t_spot) ,pointer :: s        ! new report meta data
    type(t_spot)          :: s2       ! backup copy of empty new report
    integer               :: typ     (7*size(temp))
    real(sp)              :: lev     (7*size(temp))
    real(wp)              :: olev    (7*size(temp))
    real(sp)              :: lev_lat (7*size(temp))
    real(sp)              :: lev_lon (7*size(temp))
    integer               :: lev_tim (7*size(temp))
    type(t_datum)         :: bod     (7*size(temp)) ! body entry
    logical               :: mask    (  size(temp))
    logical               :: luse    (  size(temp)) ! use this level?
    integer               :: k, nk, np, k0, i
    integer               :: no       ! number of observations in report
    integer               :: nlev     ! number of levels       in report
    real(sp)              :: p        ! pressure of current level
    real(sp)              :: ps       ! surface pressure
    real(sp)              :: ptp      ! tropopause level pressure in report
    character(len=32)     :: com      ! comment line
    logical               :: lkeep    ! keep report
    integer               :: vss,vss0 ! vertical sounding significance (masked)
    logical               :: hres     ! high-resolution BUFR TEMP?
    integer               :: minsec   ! minimum time offset
    integer               :: maxsec   ! maximum time offset
    real(sp)              :: zav      ! average geopotential of set of levels
    real(sp)              :: gp       ! temporary geopotential
    real(sp)              :: plev     ! pressure level (if p is used)
    real(sp)              :: rlev     ! reported level (p or z)
    integer               :: k1       ! temporary index
    logical               :: ldown    ! downsonde?
    integer               :: t0,t1    ! time "levels"
    integer               :: lpass    ! local copy of pass
    integer               :: levtype  ! report level type
    real(wp)              :: ff, dd   ! temporaries: wind speed, direction
    real(sp),parameter    :: t0c_sp   = real (t0c,sp) ! 273.15K
    real(wp),parameter    :: cut      = cos (5 * d2r) ! cos(5 deg)
    real(wp)              :: ref_uvec(3)
    real(wp)              :: tmp_uvec(3)
    !-----------------------------------------
    ! Mask for "classical" level significances
    !-----------------------------------------
    integer, parameter    :: vss_mask = lev_surf  + lev_std   + lev_tropo &
                                      + lev_sig_t + lev_max_v + lev_sig_v
    !--------------------------------------------------------------------
    ! Level significances conditionally accepted for zero or missing time
    !--------------------------------------------------------------------
    integer, parameter    :: vss_t0miss = lev_std + lev_tropo

    com  = "check_store_temp: "
    nlev = size (temp)
    !---------------------------
    ! determine surface pressure
    !---------------------------
    ps = huge(ps)
    k0 = 0
    do k = 1, nlev
      if (      temp(k)% p% qc          == QC_OK .and. &
          iand (temp(k)% vss, lev_surf) /= 0           ) then
        ps = temp(k)% p% o
        k0 = k
        exit
      endif
    end do
    !-----------------------------------------
    ! Check for bogus/multiple surface levels,
    ! as seen on a Korean radiosonde.
    ! (Should we dismiss the entire ascent?)
    !-----------------------------------------
    if (k0 > 0) then
       do k = k0+1, size (temp)
          if (iand (temp(k)% vss, lev_surf) /= 0) then
             temp(k)% p% qc =  QC_NOUSE
!print *, "### remove bogus surface level: ", spot% statid, ":", temp(k)% p% o
          endif
       end do
    end if

    !-------------------------------------------------
    ! Dismiss upper-air levels with T=TD=0 C, FF=0 m/s
    !-------------------------------------------------
    if (rm_t0c_ff0) then
       do k = 1, nlev
          if (iand (temp(k)% vss, lev_surf) == 0        .and. &
               abs (temp(k)% t % o - t0c_sp) < 1.e-4_sp .and. &
               abs (temp(k)% td% o - t0c_sp) < 1.e-4_sp .and. &
                    temp(k)% ff% o          == 0._sp          ) then
             temp(k)% p% qc = QC_NOUSE
!print *, "### TEMP with bad missing values: ", spot% statid, ":", temp(k)% p% o
          end if
       end do
    end if

    !----------------------------------------------------------
    ! Treat as high-resolution TEMP for duplicate-level checks?
    ! Apply stricter high-resolution checks only for pass <= 1.
    !----------------------------------------------------------
    minsec = minval (temp(:)% secs)
    maxsec = maxval (temp(:)% secs)
    hres   = maxsec > 0
    lpass  = 1; if (present (pass)) lpass = pass
    if (lpass > 1) hres = .false.
    !--------------------------------------------------
    ! Downsonde detection, taking into account that the
    ! lowest and highest reported levels may have been
    ! extrapolated, and having no valid time
    !--------------------------------------------------
    ldown  = .false.
    if (hres) then
      do k = 1, nlev
         if (temp(k)% p% qc /= QC_OK) cycle
         t0 = temp(k)% secs
         if (t0 /= 0) exit
      end do
      do k = nlev, 1, -1
         if (temp(k)% p% qc /= QC_OK) cycle
         t1 = temp(k)% secs
         if (t1 /= 0) exit
      end do
      ldown = t0 > t1
    end if

    !-----------------------------------------------------------------------
    ! Consider "surface levels" of downsondes as generally suspicious.
    ! They are extrapolated (c.f. regulation: B/C26.6.1) or sometimes bogus.
    !-----------------------------------------------------------------------
    if (ldown .and. lpass == 1 .and. ps < huge(ps)) then
       do k = 1, nlev
          if (      temp(k)% p% qc          == QC_OK .and. &
              iand (temp(k)% vss, lev_surf) /= 0           ) then
             temp(k)% p% qc = QC_NOUSE
!print *, "### bogus downsonde surface level: ", spot% statid, ":", temp(k)%p% o
          end if
       end do
       ps = huge(ps)
    end if

    !------------------------------------
    ! determine tropopause level pressure
    ! (only for first checking pass)
    !------------------------------------
    ptp = -1._sp
    if (top_q_mode > 0 .and. lpass <= 1) then
      do k = 1, nlev
        if (      temp(k)% p% qc           == QC_OK .and. &
            iand (temp(k)% vss, lev_tropo) /= 0           ) then
          ptp = temp(k)% p% o
          exit
        endif
      end do
      if (ptp > 0._sp) then
        ! top_q_mode = 1: use humidity at tropopause and below
        !            > 1: use humidity only below tropopause
        if (top_q_mode > 1) ptp = nearest (ptp, +1._sp)
      end if
    end if

    !---------------------------------------------------------------
    ! Determine level type used for profiler from first valid level,
    ! but require pressure level for TEMP and classic PILOT reports
    !---------------------------------------------------------------
    levtype = VN_P

    if (     spot% hd% obstype  == OT_PILOT .and. &
        any (spot% hd% codetype == ct_wprof)      ) then
      lvtyp: block
        !----------------------------------------------
        ! First check for any level with valid pressure
        !----------------------------------------------
        do k = 1, nlev
           if (temp(k)% p%  qc == QC_OK) then
              levtype = temp(k)% p%  lev_typ
              if (levtype == -1) levtype = VN_P
              exit lvtyp
           end if
        end do
        !-------------------------------------------
        ! Then check for any level with valid height
        !-------------------------------------------
        if (any (temp(1:nlev)% height > -999._sp)) then
           levtype = VN_HEIGHT
           exit lvtyp
        end if
        !-----------------------------------------------------------
        ! Finally check for any level with valid geopotential height
        !-----------------------------------------------------------
        do k = 1, nlev
           if (temp(k)% gp% qc == QC_OK) then
              levtype = temp(k)% gp% lev_typ
              if (levtype == -1) levtype = VN_P
              exit lvtyp
           end if
        end do
      end block lvtyp

      if (.not. any (levtype == [VN_P,VN_HEIGHT])) then
         write(0,*) "check_store_temp: statid,codetype,levtype= ", &
              spot% statid, spot% hd% codetype, levtype
         levtype = -1
      end if
    end if

    !--------------------------------------------------
    ! loop over levels
    ! check for invalid levels first, merge twin levels
    !--------------------------------------------------
    p          = huge(p)
    k0         = -1
    t0         = 0
    if (ldown) then
       t0      = maxsec
    else
       t0      = minsec
    end if
    do k = 1 ,size(temp)
      if (temp(k)% p% qc /= QC_OK) cycle
      !---------------------------------------
      ! only accept decreasing pressure values
      !             p < ps
      !---------------------------------------
      if (lpr_temp .or. prt_data) then
         if (temp(k)% p% o > ps) then
           print '(a,i3,a,i3,a,2f12.3)' ,'pe=',dace% pe,' check_store_temp  k=',k,&
                 '  p > psurf:', temp(k)% p% o, ps
         endif
      endif
      if (temp(k)% p% o <=      0._sp) temp(k)% p% qc =  QC_INCON ! check for ..
      if (temp(k)% p% o >  120000._sp) temp(k)% p% qc =  QC_INCON ! .. invalid values
      if (temp(k)% p% o >  p         ) temp(k)% p% qc =  QC_INCON ! monotonous?
      if (temp(k)% p% o >  ps        ) then
         !-----------------------------------------------------
         ! c.f. WMO Regulations for TEMP reporting, B/C25.8.2.2
         !-----------------------------------------------------
         if (temp(k)% t % qc == QC_OK) temp(k)% t % qc = QC_NOUSE ! or QC_MISS
         if (temp(k)% td% qc == QC_OK) temp(k)% td% qc = QC_NOUSE
         if (temp(k)% ff% qc == QC_OK) temp(k)% ff% qc = QC_NOUSE
         if (temp(k)% dd% qc == QC_OK) temp(k)% dd% qc = QC_NOUSE
         !---------------------------------------------------
         ! remove non-standard extrapolations of geopotential
         !---------------------------------------------------
         if (hres .and. iand (temp(k)% vss, lev_std) == 0)      &
                                       temp(k)% p% qc =  QC_INCON
      end if
      if (temp(k)% p% o < top_p .and.  temp(k)% p% qc == QC_OK) &
                                       temp(k)% p% qc =  QC_NOUSE
      !----------------------------------------------------------------
      ! time offset must be monotonous, except for extrapolated levels,
      ! (c.f. B/C25.8.2.4), and excepting artificial levels.
      !----------------------------------------------------------------
      if (hres .and. .not. btest (temp(k)% p% lev_sig, LS_SUPEROBS)) then
         if (ldown) then
            if (temp(k)% secs > t0)    temp(k)% p% qc =  QC_INCON
         else
            if (      temp(k)% secs             == 0 .and.           &
                iand (temp(k)% vss, vss_t0miss) /= 0 .and. std_t0miss) then
               ! use standard levels with zero or missing time
            else
               if (   temp(k)% secs < t0          .and.            &
                  (k < nlev .or. iand (temp(k)% vss, lev_std +     &
                                                     lev_other) == 0)) then
                                       temp(k)% p% qc =  QC_INCON
               end if
            end if
         end if
      end if
      !------------------------------------------------------------------------
      ! High-resolution downsondes: skip levels without valid drift information
      !------------------------------------------------------------------------
      if (hres .and. ldown) then
         if (temp(k)% dlat == rvind .or. &
             temp(k)% dlon == rvind      ) temp(k)% p% qc = QC_NOUSE
      end if
      !-------------------------------------------------------------
      ! High-resolution profiles: skip levels near effective 0E/0N:
      ! bad drift data, e.g. Modem M20, TEMP 04339 on 2023-04-25 23Z
      !-------------------------------------------------------------
      if (hres) then
         if (abs (temp(k)% dlat) < 1.e-3_sp .and. &
             abs (temp(k)% dlon) < 1.e-3_sp       ) temp(k)% p% qc = QC_NOUSE
      end if
      !------------------------------------------------------------
      ! High-resolution profiles: discard levels with very unlikely
      ! or obviously wrong drift data
      !------------------------------------------------------------
      if (hres .and. temp(k)% p% qc == QC_OK) then
         if (       temp(k)% dlat                       /= rvind  .and. &
              (abs (temp(k)% dlat)                      >  90._wp .or.  &
               abs (temp(k)% dlat - spot% col% c% dlat) >   5._wp      )) then
            temp(k)% p% qc = QC_NOUSE
         end if
         if (       temp(k)% dlon                       /= rvind  .and. &
               abs (temp(k)% dlon - spot% col% c% dlon) >  10._wp       ) then
            !------------------------------------------------------------------
            ! Check great-circle distance on unit-sphere for near-polar ascents
            !------------------------------------------------------------------
            ref_uvec = unitvector (dlat=spot% col% c% dlat, &
                                   dlon=spot% col% c% dlon  )
            tmp_uvec = unitvector (dlat=real (temp(k)% dlat,wp), &
                                   dlon=real (temp(k)% dlon,wp)  )
            if (sum (tmp_uvec * ref_uvec) < cut) temp(k)% p% qc = QC_NOUSE
         end if
      end if
      !------------------
      ! merge twin levels
      !------------------
      if (temp(k)% p% o == p    .and.  temp(k)% p% qc == QC_OK) then
        !-----------------------------------------------------------------
        ! For high-resolution radiosondes we apply the following policy
        ! to multiple input levels with same the coded pressure value:
        ! - level significance 0 is assumed as raw data, values /= 0
        !   indicate derived (e.g. standard or significant) levels.
        ! - same time stamps indicate redundant data.
        ! - we want to keep the data with highest vertical significance.
        !   (careful to ignore the "lev_other" bit in comparisons!)
        ! - for different time stamps we choose a "central" representative
        !   level below to prevent selection bias.  The decision will be
        !   based on geopotential height, which has higher effective
        !   resolution in the BUFR as pressure where this matters.
        !-----------------------------------------------------------------
        if (hres) then
          vss0 = iand (temp(k0)% vss, vss_mask)
          vss  = iand (temp(k) % vss, vss_mask)
          if (vss0 /= 0 .or. vss /= 0 .or. temp(k)% secs == temp(k0)% secs) then
             if (vss0 >= vss) then
                temp(k )% p% qc = QC_INCON   ! Mark for later removal
             else
                temp(k0)% p% qc = QC_INCON   ! Mark for later removal
             end if
          else
             ! vss=0 for both, different times: defer duplicate resolution
          end if
        else
          call merge_level (temp(k), temp(k0))
          temp(k0)% p% qc = QC_INCON      ! Mark merged part for later removal
        end if

!!! Consistency check of level significances:
!print*,"after merge_level:p,vss,sig(p,t,td,rh,ff,v):",p,int(temp(k)%vss,i2),temp(k)%p%lev_sig,&
!temp(k)%t%lev_sig,temp(k)%td%lev_sig,temp(k)%rh%lev_sig,temp(k)%ff%lev_sig,temp(k)%uu%lev_sig
      endif

      if(temp(k)% p% qc == QC_OK) then
        p  = temp(k)% p% o
        k0 = k
        if (.not. btest (temp(k)% p% lev_sig, LS_SUPEROBS) &
            .and.        temp(k)% secs /= 0                ) t0 = temp(k)% secs
      endif
    end do
    !-------------------------------------------------------------
    ! Resolve duplicate pressure levels for high-resolution sondes
    !-------------------------------------------------------------
    if (hres) then
       do k = size (temp), 2, -1
          if (temp(k)% p% qc == QC_OK) then
             p   = temp(k)% p% o
             vss = iand (temp(k)% vss, vss_mask)
             if (vss /= 0) then
                !--------------------------------------------------
                ! Discard all other levels with same coded pressure
                !--------------------------------------------------
                do k0 = k-1, 1, -1
                   if (temp(k0)% p% o  /= p)     exit
                   if (temp(k0)% p% qc == QC_OK) temp(k0)% p% qc = QC_INCON
                end do
             else
                !-------------------------------------------------
                ! Search for lowest level with same coded pressure
                !-------------------------------------------------
                k1 = k
                do k0 = k-1, 1, -1
                   if (temp(k0)% p% o  /= p)     exit
                   k1 = k0
                end do
                if (k1 < k) then
                   mask(k1:k) = temp(k1:k)% p% qc == QC_OK
                   zav = (maxval (temp(k1:k)% gp% o, mask=mask(k1:k)) +        &
                          minval (temp(k1:k)% gp% o, mask=mask(k1:k)) ) * 0.5_sp
                   k0  =  minloc (abs (temp(k1:k)% gp% o - zav),               &
                                  mask=mask(k1:k), dim=1) + k1-1
                   if (k0 >= k1 .and. k0 <= k) then
                      temp(k1:k0-1)% p% qc = QC_INCON
                      temp(k0+1:k) % p% qc = QC_INCON
                   end if
                end if
             end if
          end if
       end do
    end if
    !------------------------------------------
    ! Check plausibility of geopotential height
    ! for profilers it must be ascending
    !------------------------------------------
    do k = 1, size (temp)
      if (temp(k)% gp% qc == QC_OK .and. (temp(k)% gp% o < -500._sp .or. &
                                          temp(k)% gp% o > 6.0e4_sp     )) then
          temp(k)% gp% qc =  QC_CLIM  ! Geopot. height out of expected range
      end if
    end do
    if (spot% hd% obstype == OT_PILOT .and. levtype == VN_HEIGHT) then
       gp = - huge (gp)
       do k = 1, size (temp)
          if (temp(k)% gp% qc /= QC_OK) cycle
          if (temp(k)% gp% o < gp) then
             temp(k)% gp% qc = QC_INCON         ! not monotonous
          else
             gp = temp(k)% gp% o
          end if
       end do
       gp = - huge (gp)
       do k = 1, size (temp)
          if (temp(k)% height == -999._sp) cycle
          if (temp(k)% height < gp) then
             temp(k)% gp% qc = QC_INCON         ! not monotonous
             temp(k)% height = -999._sp
          else
             gp = temp(k)% height
          end if
       end do
    end if
    !-------------------------------------------
    ! count number of degrees of freedom / level
    !-------------------------------------------
    luse = .false.
    no = 0
    k0 = 0
    do k = 1 ,size(temp)
      nk = 0
      select case (levtype)
      case (VN_P, -1)
         if (temp(k)% p% qc == QC_OK) then
            plev = temp(k)% p% o
            rlev = plev
            levtype = VN_P
            luse(k) = .true.
         end if
      case (VN_HEIGHT)
         if (temp(k)% height >  -999.) then
            plev = inv_datum% plev  ! -1.
            rlev = temp(k)% height  ! [gpm]
            luse(k) = .true.
         end if
!     case default
!        write(0,*) "??? statid,qc: ", spot% statid, temp(k)% p% qc, &
!             temp(k)% gp% qc, temp(k)% height, levtype
      end select
      if (luse(k)) then
        !------------------------------------------------------------------
        ! Check drift information, remember last valid horizontal position.
        ! Do not use wind data where drift information was missing.
        !------------------------------------------------------------------
        if (temp(k)% dlat /= rvind .and. temp(k)% dlon /= rvind) then
           k0 = k
        else
           if (      temp(k)% secs             == 0 .and.           &
               iand (temp(k)% vss, vss_t0miss) /= 0 .and. std_t0miss) then
              !--------------------------------------------------------------
              ! Exception: allow standard levels for badly formatted reports.
              ! We only trust the interpolation, but not the location,
              ! unless it is an extrapolated level at the top.
              !--------------------------------------------------------------
              if (k0 > 0 .and. k == size (temp)) then
                 temp(k)% dlat   = temp(k0)% dlat
                 temp(k)% dlon   = temp(k0)% dlon
              end if
           else if (k0 > 0) then
              temp(k)% dlat   = temp(k0)% dlat
              temp(k)% dlon   = temp(k0)% dlon
              temp(k)% ff% qc = QC_NOUSE
              temp(k)% dd% qc = QC_NOUSE
           end if
        end if
        !--------------------
        ! geopotential height
        !--------------------
        if (temp(k)% gp% qc == QC_OK .and. &
            spot% hd% obstype /= OT_PILOT ) then
          if (iand(passive_h, temp(k)% vss) /=0 ) then
            no         = no + 1
            nk         = nk + 1
            typ    (no)= VN_Z
            lev    (no)= plev
            olev   (no)= rlev
            lev_lat(no)= temp(k)% dlat
            lev_lon(no)= temp(k)% dlon
            lev_tim(no)= temp(k)% secs
            bod    (no)= temp(k)% gp
          endif
        endif
        !------------
        ! temperature
        !------------
        if (temp(k)% t% qc == QC_OK .and. (temp(k)% t% o < 100._sp .or. &
                                           temp(k)% t% o > 400._sp     )) then
            temp(k)% t% qc =  QC_CLIM   ! Temperature out of expected range
        end if
        if (temp(k)% t% qc == QC_OK) then
          if (iand(passive_t, temp(k)% vss) /=0 ) then
            no         = no + 1
            nk         = nk + 1
            typ    (no)= VN_T
            lev    (no)= plev
            olev   (no)= rlev
            lev_lat(no)= temp(k)% dlat
            lev_lon(no)= temp(k)% dlon
            lev_tim(no)= temp(k)% secs
            bod    (no)= temp(k)% t
          endif
        endif
        !---------
        ! humidity
        !---------
        if (temp(k)% td% qc == QC_OK .and. (temp(k)% td% o < 100._sp .or. &
                                            temp(k)% td% o > 400._sp     )) then
            temp(k)% td% qc =  QC_CLIM      ! Dew-point temperature out of range
        end if
        if (temp(k)% td% qc == QC_OK .and.  &
            temp(k)% t%  qc == QC_OK .and.  &
            temp(k)% td% o  >  temp(k)% t% o) then
            temp(k)% td% qc =  QC_INCON     ! Invalid dew-point temperature
        end if
        if (temp(k)% p% o < ptp) then       ! Humidity data above tropopause
           if (temp(k)% td% qc == QC_OK) temp(k)% td% qc = QC_NOUSE
           if (temp(k)% rh% qc == QC_OK) temp(k)% rh% qc = QC_NOUSE
        end if
        !-------------------------------------------------------------
        ! Derive from dewpoint temperature if rel. humidity is missing
        !-------------------------------------------------------------
        if (temp(k)% rh% qc == QC_MISS .and.                        &
            temp(k)% td% qc == QC_OK   .and. temp(k)% t% qc == QC_OK) then
          if (iand(temp(k)% vss, passive_q) /= 0        .and. &
                   temp(k)% p% o            >= 100._wp * top_q) then
            temp(k)% rh% src = SRC_DOK
            temp(k)% rh% qc  = temp(k)% td% qc
            select case (satvp_form)
            case (1) ! Magnus-Tetens
              temp(k)% rh% o = esw_t (real(temp(k)% td% o, wp)) &
                             / esw_t (real(temp(k)% t%  o, wp))
            case (2) ! Hardy: used by Vaisala and other manufacturers
              temp(k)% rh% o = esw_t_hardy (real(temp(k)% td% o, wp)) &
                             / esw_t_hardy (real(temp(k)% t%  o, wp))
            end select
            !-----------------------------------------------------
            ! Set significance of derived quantities from observed
            !-----------------------------------------------------
            temp(k)% rh% lev_sig = temp(k)% td% lev_sig
          endif
        end if
        if (temp(k)% rh% qc == QC_OK .and. temp(k)% rh% o <  0._sp) then
            temp(k)% rh% qc =  QC_MISS
        end if
        if (temp(k)% rh% qc == QC_OK .and.    &
            temp(k)% p % o  >= 100._wp * top_q) then
          no         = no + 1
          nk         = nk + 1
          typ    (no)= VN_RH
          lev    (no)= plev
          olev   (no)= rlev
          lev_lat(no)= temp(k)% dlat
          lev_lon(no)= temp(k)% dlon
          lev_tim(no)= temp(k)% secs
          bod    (no)= temp(k)% rh
        endif
        !-----
        ! wind
        !-----
        if (temp(k)% ff% qc == QC_OK .and. temp(k)% dd% qc == QC_OK) then
          temp(k)% uu% src = SRC_DOK
          temp(k)% uu% qc  = QC_OK
          temp(k)% uu% o   = temp(k)% ff% o * (-sin (d2r * temp(k)% dd% o))
          temp(k)% vv% src = SRC_DOK
          temp(k)% vv% qc  = QC_OK
          temp(k)% vv% o   = temp(k)% ff% o * (-cos (d2r * temp(k)% dd% o))
          !-----------------------------------------------------
          ! Set significance of derived quantities from observed
          !-----------------------------------------------------
          temp(k)% uu% lev_sig = temp(k)% ff% lev_sig
          temp(k)% vv% lev_sig = temp(k)% ff% lev_sig
        else if (spot% hd% codetype == OC_SCADA) then
          !-------------------------------------------
          ! Wind speed-only measurements from "towers"
          !-------------------------------------------
        else if (temp(k)% uu% qc == QC_OK .and. temp(k)% vv% qc == QC_OK) then
          !------------------------------------------------------
          ! Set wind speed, direction from wind vector components
          ! set significance of derived quantities from observed
          !------------------------------------------------------
          call fd_uv (ff, dd, real(temp(k)% uu% o,wp), real(temp(k)% vv% o,wp))
          temp(k)% ff% o       = ff
          temp(k)% dd% o       = dd
          temp(k)% ff% src     = SRC_DOK
          temp(k)% dd% src     = SRC_DOK
          temp(k)% ff% qc      = QC_OK
          temp(k)% dd% qc      = QC_OK
          temp(k)% ff% lev_sig = temp(k)% uu% lev_sig
          temp(k)% dd% lev_sig = temp(k)% uu% lev_sig
        else
          !----------------------------------------------
          ! set u,v quality flags consistently with ff,dd
          !----------------------------------------------
          if(temp(k)% ff% qc == QC_OK .or.     &
             temp(k)% ff% qc == temp(k)% dd% qc) then
             temp(k)% uu% qc =  temp(k)% dd% qc
             temp(k)% vv% qc =  temp(k)% dd% qc
          else if (temp(k)% dd% qc == QC_OK) then
             temp(k)% uu% qc =  temp(k)% ff% qc
             temp(k)% vv% qc =  temp(k)% ff% qc
          else
             temp(k)% uu% qc =  QC_INCON
             temp(k)% vv% qc =  QC_INCON
          endif
        endif
        if (iand(temp(k)% vss, passive_v) /= 0) then
          if (temp(k)% uu% qc == QC_OK) then
            no         = no + 1
            nk         = nk + 1
            typ    (no)= VN_U
            lev    (no)= plev
            olev   (no)= rlev
            lev_lat(no)= temp(k)% dlat
            lev_lon(no)= temp(k)% dlon
            lev_tim(no)= temp(k)% secs
            bod    (no)= temp(k)% uu
          endif
          if (temp(k)% vv% qc == QC_OK) then
            no         = no + 1
            nk         = nk + 1
            typ    (no)= VN_V
            lev    (no)= plev
            olev   (no)= rlev
            lev_lat(no)= temp(k)% dlat
            lev_lon(no)= temp(k)% dlon
            lev_tim(no)= temp(k)% secs
            bod    (no)= temp(k)% vv
          endif
          if (temp(k)% ff% qc == QC_OK) then
            if ((temp(k)% dd% qc    == QC_OK .and. monitor_ff) .or. &
                 spot% hd% codetype == OC_SCADA                     ) then
              no         = no + 1
              nk         = nk + 1
              typ    (no)= VN_FF
              lev    (no)= plev
              olev   (no)= rlev
              lev_lat(no)= temp(k)% dlat
              lev_lon(no)= temp(k)% dlon
              lev_tim(no)= temp(k)% secs
              bod    (no)= temp(k)% ff
            endif
          endif
          if (temp(k)% dd% qc == QC_OK) then
            if (temp(k)% ff% qc == QC_OK .and. monitor_dd) then
              no         = no + 1
              nk         = nk + 1
              typ    (no)= VN_DD
              lev    (no)= plev
              olev   (no)= rlev
              lev_lat(no)= temp(k)% dlat
              lev_lon(no)= temp(k)% dlon
              lev_tim(no)= temp(k)% secs
              bod    (no)= temp(k)% dd
            endif
          endif
        endif
      endif
      !----------------------------------------------------------
      ! skip level if there is no valid observation at this level
      !----------------------------------------------------------
      if (nk == 0) luse(k) = .false.
    end do

    !------------------------------
    ! shrink buffer to minimum size
    !------------------------------
    np = count (luse)
    if (np/=size(temp)) then
      t => temp
      allocate (temp(np))
      temp = pack (t, luse)
      deallocate (t)
    endif
    !---------------------------
    ! check for valid station Id
    !---------------------------
    if (spot% statid == '00000' .or. spot% statid == '') then
       no  = 0
       com = "check_store_temp: station id"
    end if
    !-----------------------------------------
    ! Drop bogus reports at exactly (0°N,0°E).
    !-----------------------------------------
    if (spot% col% c% dlat == 0._wp .and. &
        spot% col% c% dlon == 0._wp) then
       no  = 0
       com = "check_store_temp: location"
    end if
    !--------------------
    ! store data into OBS
    !--------------------
    novalid = no
    lkeep = no > 0
    if (lkeep) then
      if (present (repl)) then
        s => obs% spot (repl)
      else
        call new_spot (obs, 1, set_id=.true.)
        s => obs% spot (obs% n_spot)
      endif
      s2           = s
      s            = spot
      s% id        = s2% id
      s% o         = s2% o
      s% p         = s2% p
      s% col% nlev = np
      if (temp_int_size == 0) call set_size
      call new_obs (obs, no, s)
      obs% varno (s%o%i+1:s%o%i+s%o%n)          = typ    (1:no)
      obs% body  (s%o%i+1:s%o%i+s%o%n)          = bod    (1:no)
      obs% body  (s%o%i+1:s%o%i+s%o%n)% plev    = lev    (1:no)
      obs% olev  (s%o%i+1:s%o%i+s%o%n)          = olev   (1:no)
      obs% body  (s%o%i+1:s%o%i+s%o%n)% lev_typ = levtype
      obs% body  (s%o%i+1:s%o%i+s%o%n)% lat     = lev_lat(1:no)
      obs% body  (s%o%i+1:s%o%i+s%o%n)% lon     = lev_lon(1:no)
      obs% body  (s%o%i+1:s%o%i+s%o%n)% secs    = lev_tim(1:no)
      call store_temp (obs, s, temp)
      call set_xuv (s)
!     if (s% hd% codetype == OC_SCADA) then
!        ! SCADA (FF only)
!        call decr_rpt_use (s, CHK_NOTUSED, STAT_PASSIVE)
!     end if
      !--------------------------------------
      ! certain quantities are always passive
      !--------------------------------------
      do i = 1, s% o% n
        select case (obs% varno(s% o% i+i))
        case (VN_FF, VN_DD)
          if (s% hd% codetype /= OC_SCADA) then
            ! When u,v are available, set ff,dd to passive
            call decr_use (obs% body(s% o% i+i)% use, STAT_PASSIVE, CHK_NOTUSED)
          end if
        end select
      end do
    else
      call decr_rpt_use (spot ,CHK_INSDAT, STAT_DISMISS, comment=trim (com))
    endif
   contains
!------------------------------------------------------------------------------
    subroutine print_line
      integer :: ip, it, itd, iv, igp
      ip  = count (temp% p  %qc == QC_OK)
      it  = count (temp% t  %qc == QC_OK)
      itd = count (temp% td %qc == QC_OK .and. temp% t  %qc == QC_OK)
      iv  = count (temp% ff %qc == QC_OK .and. temp% dd %qc == QC_OK)
      igp = count (temp% gp %qc == QC_OK)
      write(6,                                                         &
        '(a,":pe,stid,lon,lat,kz,n:p,t,td,v,gp=",i4,i7,2f7.2,1x,6i3)') &
        obstyp(spot% hd% obstype)% name,                               &
        dace% pe, spot% ident, spot% col% c% dlon, spot% col% c% dlat, &
        spot% hd% dbkz, ip, it, itd, iv, igp
    end subroutine print_line
  end subroutine check_store_temp
!==============================================================================
  !---------------------------
  ! Private auxiliary routines
  !---------------------------
!------------------------------------------------------------------------------
  subroutine set_size
  !-------------------------------------------------------
  ! store sizes of derived data types T_TEMP
  ! (in termes of size of component OBS% PAR) into
  ! private module variables TEMP_INT_SIZE
  !-------------------------------------------------------
    type (t_temp)  :: temp
    type (t_obs)   :: obs
    if (temp_int_size == 0) &
      temp_int_size = size (transfer (temp, obs% par))
  end subroutine set_size
!------------------------------------------------------------------------------
  subroutine load_temp (obs, spot, temp)
  type (t_obs)  ,intent(in)  :: obs      ! data of all observations
  type (t_spot) ,intent(in)  :: spot     ! meta data of this observation
  type (t_temp) ,pointer     :: temp (:) ! TEMP level information
  !------------------------------------------------------------------------
  ! Load the data from components PAR, OBS of OBS from position provided by
  ! SPOT. Store into TEMP. allocate TEMP with size required.
  !------------------------------------------------------------------------
    allocate (temp (spot% col% nlev))
    if (temp_int_size == 0) call set_size
      if (spot% p% n <= 0) call finish('load_temp','spot% p% n <= 0')
      temp = transfer (obs% par (spot% p% i+1 : spot% p% i + spot% p% n), temp)
  end subroutine load_temp
!------------------------------------------------------------------------------
  subroutine store_temp  (obs, spot, temp)
  type (t_obs)  ,intent(inout) :: obs      ! data of all observations
  type (t_spot) ,intent(inout) :: spot     ! meta data of this observation
  type (t_temp) ,intent(in)    :: temp (:) ! TEMP level information
  !-----------------------------------------------------------------------
  ! Store the data from variable TEMP in the component PAR of
  ! OBS at position provided by SPOT. Allocate memory for PAR if required.
  !-----------------------------------------------------------------------
    integer :: par(1)  !+++ work around bug in NEC sxf90/2.0rev360-ia64
    if (temp_int_size == 0) call set_size
    call new_par (obs, spot% col% nlev * temp_int_size, spot=spot)
    obs % par (spot% p% i+1 : spot% p% i + spot% p% n) = &
      transfer(temp, par)
  end subroutine store_temp
!------------------------------------------------------------------------------
  subroutine read_fdbk_temp (o)
  !-----------------------------------------------------------------
  ! fix level count after observation was read from feedback file.
  ! - Occasionally fof-files contain surface level data
  !   at the same level as regular data.
  ! - Occasionally fof-files contain standard pressure level data
  !   at the same level as regular data.
  ! This data stems from incorrect BUFR reports and is removed here.
  !-----------------------------------------------------------------
  type (t_obs) ,intent(inout) :: o          ! observation data type variable

    integer ,parameter :: mm = 10 ! max # of standard pressure levels
    integer            :: i       ! report index
    integer            :: j, k    ! observation body indices
    integer            :: l       ! level count without separating surface data
    integer            :: ls      ! level count with    separating surface data
    integer            :: lm      ! level count with    main pressure level
    real(wp)           :: pstd(mm)! double standard pressure levels
    real(wp)           :: olev    ! current level
    integer            :: lsigs   !
    integer            :: lsigm   !

    !-----------------------
    ! loop over TEMP reports
    !-----------------------
    do i = 1, o% n_spot
      if (o% spot(i)% hd% obstype /= OT_TEMP) cycle
      !-------------------------------------
      ! loop over observations in the report
      !-------------------------------------
      l    = 0             ! level count
      ls   = 0             ! dispensable surface  level count
      lm   = 0             ! dispensable standard level count
      pstd = -huge(olev)   ! list of dispensable standard levels
      olev = -huge(olev)   ! current level
      do j = 1, o% spot(i)% o% n
        k = j + o% spot(i)% o% i
        !-----------------------
        ! detect new level value
        !-----------------------
        if (o% olev(k) /= olev) then
          l     = l + 1
          lsigs = 0      ! no double surface layer
          lsigm = 0      ! no double standard layer
          olev = o% olev(k)
        endif
        !---------------------------------
        ! detect for double surface layers
        !---------------------------------
        if (btest (o% body(k)% lev_sig, LS_SURFACE)) then
          select case (lsigs)
          case (0)
            lsigs = 1        ! this level is surface layer
          case (2)
            lsigs = 3        ! this level is surface + other layer
            ls    = ls + 1
          end select
        else
          select case (lsigs)
          case (0)
            lsigs = 2        ! this level is other layer
          case (1)
            lsigs = 3        ! this level is surface + other layer
            ls    = ls + 1
          end select
        endif
        !------------------------------------
        ! test for double main pressure level
        !------------------------------------
        if (btest (o% body(k)% lev_sig, LS_STANDARD)) then
          select case (lsigm)
          case (0)
            lsigm = 1       ! this level is standard layer
          case (2)
            lsigm = 3       ! this level is standard + other layer
            lm    = lm + 1
            if (lm>mm) call finish('read_fdbk_temp','lm>mm')
            pstd(lm) = olev
          end select
        else
          select case (lsigm)
          case (0)
            lsigm = 2       ! this level is other layer
          case (1)
            lsigm = 3       ! this level is standard + other layer
            lm    = lm + 1
            if (lm>mm) call finish('read_fdbk_temp','lm>mm')
            pstd(lm) = olev
          end select
        endif
      end do
      !-----------------------------------------------------------------
      ! do nothing if level count from fof corresponds with that in DACE
      !-----------------------------------------------------------------
      if (l == o% spot(i)% col% nlev) cycle
      !------------------------------------------------------------
      ! remove surface observations and adjust level count from fof
      !------------------------------------------------------------
      if (ls > 1) call finish('read_fdbk_temp','more than 1 surface layer')
      if (ls == 1) then
        do j = 1, o% spot(i)% o% n
          k = j + o% spot(i)% o% i
          if (btest (o% body(k)% lev_sig, LS_SURFACE)) then
            if (o% body(k)% use% state > STAT_DISMISS) &
                o% body(k)% use% state = STAT_DISMISS
          endif
        end do
        o% spot(i)% col% nlev = o% spot(i)% col% nlev - ls
      endif
      !-----------------------------------------------------------------
      ! do nothing if level count from fof corresponds with that in DACE
      !-----------------------------------------------------------------
      if (l == o% spot(i)% col% nlev) cycle
      !-----------------------------------------------------------------
      ! remove non-standard observations and adjust level count from fof
      !-----------------------------------------------------------------
      if (lm > 0) then
        do j = 1, o% spot(i)% o% n
          k = j + o% spot(i)% o% i
          if (any (o% olev(k) == pstd (:))) then
            if (.not. btest (o% body(k)% lev_sig, LS_STANDARD)) then
              if (o% body(k)% use% state > STAT_DISMISS) &
                  o% body(k)% use% state = STAT_DISMISS
            endif
          endif
        end do
        o% spot(i)% col% nlev = o% spot(i)% col% nlev - lm
      endif
      !-----------------------------------------------------------------
      ! do nothing if level count from fof corresponds with that in DACE
      !-----------------------------------------------------------------
      if (l == o% spot(i)% col% nlev) cycle
      !--------------------------------------------------------
      ! level count still inconsistent, remove the whole report
      !--------------------------------------------------------
      write(6,*)'read_fdbk_temp cannot fix: statid,l,ls,lm,nlev: ', &
                 o% spot(i)% statid,l,ls,o% spot(i)% col% nlev
      write(0,*)'read_fdbk_temp cannot fix: statid,l,ls,lm,nlev: ', &
                 o% spot(i)% statid,l,ls,o% spot(i)% col% nlev
!     call finish ('read_fdbk_temp','cannot fix: l,m,nlev')
      if (o% spot(i)% use% state > STAT_DISMISS) &
          o% spot(i)% use% state = STAT_DISMISS
      where (o% body(o% spot(i)% o% i + 1:&
                     o% spot(i)% o% i +   &
                     o% spot(i)% o% n     )% use% state > STAT_DISMISS) &
             o% body(o% spot(i)% o% i + 1:&
                     o% spot(i)% o% i +   &
                     o% spot(i)% o% n     )% use% state = STAT_DISMISS
      o% spot(i)% col% nlev = l
    end do

  end subroutine read_fdbk_temp
!------------------------------------------------------------------------------
  pure subroutine construct_temp (temp)
  type (t_temp) ,intent(out) :: temp (:)
!#if defined(__SX__)
!   ! Default initialisation does not work with sxf90 rev.360 (SX-6)
!   temp = t_temp (inv_datum,inv_datum,inv_datum,inv_datum,inv_datum,&
!                  inv_datum,inv_datum,inv_datum,inv_datum,          &
!                  -999.,-999.,lev_miss)
!#endif
    temp(:)% p  %mn = 'p'
    temp(:)% t  %mn = 't'
    temp(:)% td %mn = 'td'
    temp(:)% rh %mn = 'rh'
    temp(:)% ff %mn = 'ff'
    temp(:)% dd %mn = 'dd'
    temp(:)% uu %mn = 'u'
    temp(:)% vv %mn = 'v'
    temp(:)% gp %mn = 'gp'
  end subroutine construct_temp
!------------------------------------------------------------------------------
  subroutine merge_temp (dst, src)
  type (t_temp) ,pointer    :: dst (:)
  type (t_temp) ,intent(in) :: src (:)
  !-----------------------------------
  ! merge TEMP reports 'src' and 'dst'
  !-----------------------------------
    integer :: i1                     ! level index for     'dst'
    integer :: i2                     ! level index for     'src'
    integer :: n1                     ! number of levels in 'dst'
    integer :: n2                     ! number of levels in 'src'
    integer :: n                      ! number of levels in merged report
    type (t_temp) ,pointer :: new (:) ! temporary for merged report
    n1 = size (dst)
    n2 = size (src)
    allocate (new (n1+n2))
    i1 = 1
    i2 = 1
    n  = 0
    !----------------------------------------------
    ! loop over sorted (descending) pressure levels
    ! find corresponding levels in 'src' and 'dst'
    !----------------------------------------------
    do
      !-----------------------------------
      ! exit if both reports are processed
      !-----------------------------------
      if (i1 > n1 .and. i2 > n2) exit
      !--------------------
      ! skip invalid levels
      !--------------------
      if (i1 <= n1) then
        if (dst(i1)% p% qc /= QC_OK ) then
          i1 = i1 + 1
          cycle
        endif
      endif
      if (i2 <= n2) then
        if (src(i2)% p% qc /= QC_OK) then
          i2 = i2 + 1
          cycle
        endif
      endif
      !-----------------------------------
      ! process new level in merged report
      !-----------------------------------
      n = n + 1
      if (i2 > n2) then
        !------------------------------------------
        ! no more entries in 'src', take fron 'dst'
        !------------------------------------------
        new (n) = dst(i1)
        i1 = i1 + 1
        cycle
      endif
      if (i1 > n1) then
        !------------------------------------------
        ! no more entries in 'dst', take fron 'src'
        !------------------------------------------
        new (n) = src(i2)
        i2 = i2 + 1
        cycle
      endif
      if (dst(i1)% p% o > src(i2)% p% o) then
        !---------------------------------
        ! next descending level from 'dst'
        !---------------------------------
        new (n) = dst(i1)
        i1 = i1 + 1
        cycle
      endif
      if (src(i2)% p% o > dst(i1)% p% o) then
        !---------------------------------
        ! next descending level from 'src'
        !---------------------------------
        new (n) = src(i2)
        i2 = i2 + 1
        cycle
      endif
      !--------------------------
      ! 2 matching levels, merge!
      !--------------------------
      if (dst(i1)% p% o /= src(i2)% p% o) then
        call finish ('merge_temp','pressure level mismatch')
      endif
      new (n) = dst(i1)
      call merge_level (new(n), src(i2))
      i1 = i1 + 1
      i2 = i2 + 1
    end do
    deallocate (dst)
    allocate (dst (n))
    dst = new (:n)
    deallocate (new)
  end subroutine merge_temp
!------------------------------------------------------------------------------
  subroutine merge_level (dst, src)
  type (t_temp) ,intent(inout) :: dst
  type (t_temp) ,intent(in)    :: src

!logical::diff,m12,m21
!diff=.false.
!call cmp_level(dst, src, diff,m12,m21)
!if(diff) print *,'dddddddddddddddddddddddddddddd'
!if(diff) call print_temp_data((/dst/))
!if(diff) print *,'ssssssssssssssssssssssssssssss'
!if(diff) call print_temp_data((/src/))
    call merge_datum (dst% p   ,src% p   )
    call merge_datum (dst% t   ,src% t   )
    call merge_datum (dst% td  ,src% td  )
    call merge_datum (dst% rh  ,src% rh  ) ! (derived from t,td)
    call merge_datum (dst% ff  ,src% ff  )
    call merge_datum (dst% dd  ,src% dd  )
!call print (src% uu, 'src%uu (in)')
!call print (dst% uu, 'dst%uu (in)')
    call merge_datum (dst% uu  ,src% uu  ) ! (derived from ff,dd)
!call print (dst% uu, 'dst%uu (out)')
    call merge_datum (dst% vv  ,src% vv  )
    call merge_datum (dst% gp  ,src% gp  )
    dst% vss =   ior (dst% vss ,src% vss )
!if(diff) print *,'dddddddddddddddddddddddddddddd'
!if(diff) call print_temp_data((/dst/))
!if(diff) print *,'ffffffffffffffffffffffffffffff'
!!if(diff) read(5,*)
  end subroutine merge_level
!------------------------------------------------------------------------------
  subroutine cmp_temp (t1, t2, diff, m12, m21, z2)
  !------------------
  ! compare two TEMPs
  !------------------
  type (t_temp) ,intent(in)  :: t1(:), t2(:)
  logical       ,intent(out) :: diff  ! TEMPs differ in some observed value
  logical       ,intent(out) :: m12   ! TEMP 1 extends TEMP 2
  logical       ,intent(out) :: m21   ! TEMP 2 extends TEMP 1
  logical       ,intent(out) :: z2    ! no information in TEMP 2
    integer :: i1, i2, n1, n2

    diff = .false.
    m12  = .false.
    m21  = .false.
    z2   = .true.
    n1 = size (t1)
    n2 = size (t2)
    i1 = 1
    i2 = 1
    do
      !-----------------------
      ! exit if no more levels
      !-----------------------
      if (i1 > n1 .and. i2 > n2) exit
      !---------------------------------------
      ! skip invalid pressure levels in TEMP 1
      !---------------------------------------
      if (i1 <= n1) then
       if (t1(i1)% p% qc /= QC_OK) then
        i1 = i1 + 1
        cycle
       endif
      endif
      !---------------------------------------
      ! skip invalid pressure levels in TEMP 2
      !---------------------------------------
      if (i2 <= n2) then
       if (t2(i2)% p% qc /= QC_OK) then
        i2 = i2 + 1
        cycle
       endif
      endif
      !---------------------------
      ! only TEMP 1 or TEMP 2 left
      !---------------------------
      if (i1 > n1) then
        m21 = .true.
        i2 = i2 + 1
        cycle
      endif
      if (i2 > n2) then
        m12 = .true.
        i1 = i1 + 1
        cycle
      endif
      if (t1(i1)% p% o > t2(i2)% p% o) then
        m12 = .true.
        i1  = i1 + 1
        cycle
      endif
      if (t2(i2)% p% o > t1(i1)% p% o) then
        m21 = .true.
        z2  = .false.
        i2  = i2 + 1
        cycle
      endif
      if (t1(i1)% p% o /= t2(i2)% p% o) then
        call finish ('cmp_temp','pressure level mismatch')
      endif
      call cmp_level (t1(i1), t2(i2), diff, m12, m21)
      z2  = .false.
      i1 = i1 + 1
      i2 = i2 + 1
    end do
  end subroutine cmp_temp
!------------------------------------------------------------------------------
  subroutine cmp_level (l1, l2, diff, m12, m21)
  type (t_temp) ,intent(in)    :: l1, l2
  logical       ,intent(inout) :: diff
  logical       ,intent(inout) :: m12, m21
    call cmp_datum (l1% p    ,l2% p    ,diff ,m12 ,m21)
    call cmp_datum (l1% t    ,l2% t    ,diff ,m12 ,m21)
    call cmp_datum (l1% td   ,l2% td   ,diff ,m12 ,m21)
    call cmp_datum (l1% rh   ,l2% rh   ,diff ,m12 ,m21) ! (derived from t,td)
    call cmp_datum (l1% ff   ,l2% ff   ,diff ,m12 ,m21)
    call cmp_datum (l1% dd   ,l2% dd   ,diff ,m12 ,m21)
    call cmp_datum (l1% uu   ,l2% uu   ,diff ,m12 ,m21)
    call cmp_datum (l1% vv   ,l2% vv   ,diff ,m12 ,m21)
    call cmp_datum (l1% gp   ,l2% gp   ,diff ,m12 ,m21)
  end subroutine cmp_level
!------------------------------------------------------------------------------
  subroutine print_temp_datum (t)
  type (t_temp) ,intent(in) :: t
    call print (t% p  ,'pressure')
    call print (t% t  ,'temperature')
    call print (t% td ,'dew point temperature')
    call print (t% rh ,'relative humidity')
    call print (t% ff ,'wind speed')
    call print (t% dd ,'wind direction')
    call print (t% uu ,'wind component')
    call print (t% vv ,'wind component')
    call print (t% gp ,'geopotential')
    write(6,'(5x,a,f15.5,9x,a)') 'hgt',t% height,'height    [m]'
    write(6,'(5x,a,f15.5,9x,a)') 'lat',t% dlat,'latitude  [deg]'
    write(6,'(5x,a,f15.5,9x,a)') 'lon',t% dlon,'longitude [deg]'
    write(6,'(5x,a,9x,i6,9x,a)') 'sec',t% secs,'time since launch [sec]'
    write(6,'(5x,a,9x,i6,9x,a)') 'vss',t% vss ,'vertical sounding significance'
  end subroutine print_temp_datum
!------------------------------------------------------------------------------
  subroutine print_temp_data (t)
  type (t_temp) ,intent(in) :: t (:)
    integer :: i
    write(6,*)
    do i=1,size (t)
      call print_temp_datum (t(i))
    end do
  end subroutine print_temp_data
!------------------------------------------------------------------------------
  subroutine read_temp_nml
  !-------------------------
  ! read namelist /TEMP_OBS/
  !-------------------------
    integer :: ierr
    logical :: temp_obs_read = .false.
    !------------------------------
    ! don't read the namelist twice
    !------------------------------
    if (temp_obs_read) then
      if(dace% lpio) then
        write(6,'(a,/)') repeat('-',79)
        write(6,'(a,/)') '  Namelist /TEMP_OBS/ already read'
      endif
      return
    end if
    temp_obs_read = .true.
    !-------------
    ! set defaults
    !-------------
    top_p      =   0._wp   ! top level for radiosonde profiles   (hPa)
    top_q      = 300._wp   ! top level for humidity observations (hPa)
    top_q_mode = 0         ! humidity obs. handling above tropopause
    satvp_form = 1         ! Formula for sat. vapor pres. over water
    active_t   = lev_std   ! + lev_sig_t + lev_tropo
    active_q   = lev_std   ! + lev_sig_t
    active_v   = lev_std   ! + lev_sig_v + lev_max_v
    active_v_wp= lev_std   ! for wind profilers
    active_h   = lev_surf
    passive_t  = lev_all
    passive_q  = lev_all
    passive_v  = lev_all
    passive_h  = lev_all
    vmainpl    = -1._wp
    wmainpl    = -1._wp
    prt_data   = .false.
!   prt_data   = .true.
    dz_coloc   = -1._wp    ! vert. separation for collocations (Pa)
   dzv_coloc   = -1._wp    ! vert. separation for zero wind
    dx_coloc   =  1._wp    ! hori. separation for collocations (km)
    rtl_coloc  = .true.    ! remove temp levels in collocations
    split      = .false.   ! split report if drifting out of grid-cell
    std_t0miss = .true.    ! use standard levels with zero/miss. time
    rm_t0c_ff0 = .true.    ! dismiss levels with apparent missing values
    corme_mode = 0         ! Mode for correction message handling
    check_cons = 0         ! Check consistency of report
    thin_bufr  = 0         ! thinning interval for BUFR TEMPs (s)
    layer_avg  = .false.   ! apply layer averaging (superobbing)?
    supob_thr  = [300,300] ! levels in report required for superobbing
                           ! (ascending, descending)
    prefer     = 1         ! level use preference: 1=std, 2=layer-avg.
    avl_levs   = 0._wp     ! averaging layers: level definition list
    avl_incr   = 0._wp     ! averaging layers: level increments
    dp_asurf   = 5._wp     ! distance of lowest layer above surface [hPa]
    fgchk_inv_t= 0._wp     ! stability/inversion-dep. enhancement of
                           !   fg-check limit for T (recommended: .8)
    fgchk_inv_q= 0._wp     ! stability/inversion-dep. enhancement of
                           !   fg-check limit for RH (recommended: .6)
                           !   (== 0. means: no enhancement,
                           !    valid range: 0. <= fgchk_inv_x <= 1.)
    !---------------------------------
    ! read namelist, consistency check
    !---------------------------------
    if (dace% lpio) then
      call position_nml ('TEMP_OBS', status=ierr)
      select case (ierr)
      case (POSITIONED)
#if defined(__ibm__)
        read (nnml ,nml=TEMP_OBS, iostat=ierr)
        if(ierr/=0) call finish('read_temp_nml','ERROR in namelist /TEMP_OBS/')
#else
        read (nnml ,nml=TEMP_OBS)
#endif
        if (vmainpl(1) < 0._wp) then
          vmainpl (1:18) = (/1000.,925.,850.,700.,500.,400., &
                              300.,250.,200.,150.,100., 70., &
                               50., 30., 20., 10.,  7.,  5.  /)
        endif
        vmainpl = vmainpl * 100._wp  ! hPa -> Pa
        if (wmainpl(1) < 0._wp) then
          wmainpl (1:21) = (/1000.,950.,925.,900.,850.,800.,750., &
                              700.,650.,600.,550.,500.,450.,400., &
                              350.,300.,250.,200.,150.,100., 50.  /)
        endif
        wmainpl = wmainpl * 100._wp  ! hPa -> Pa
        if (all (avl_levs == 0._wp .and. avl_incr == 0._wp)) then
           avl_levs(1:7) = [ 1050, 100, 70, 20, 10, 7, 5 ]
           avl_incr(1:7) = [   25,  15, 10,  5,  3, 2, 0 ]
        end if
        avl_levs = avl_levs * 100._wp  ! hPa -> Pa
        avl_incr = avl_incr * 100._wp  ! hPa -> Pa
        dp_asurf = dp_asurf * 100._wp  ! hPa -> Pa
        top_p    = top_p    * 100._wp  ! hPa -> Pa
        fgchk_inv_t = min( max( fgchk_inv_t, 0._wp ), 1._wp )
        fgchk_inv_q = min( max( fgchk_inv_q, 0._wp ), 1._wp )
      case default
        write(6,'(/,a)') '  namelist /TEMP_OBS/ not present, using defaults'
      end select
    endif
    !------------------------
    ! broadcast to other PE's
    !------------------------
    call p_bcast (top_p,       dace% pio)
    call p_bcast (top_q,       dace% pio)
    call p_bcast (top_q_mode,  dace% pio)
    call p_bcast (active_t,    dace% pio)
    call p_bcast (active_q,    dace% pio)
    call p_bcast (active_v,    dace% pio)
    call p_bcast (active_v_wp, dace% pio)
    call p_bcast (active_h,    dace% pio)
    call p_bcast (passive_t,   dace% pio)
    call p_bcast (passive_q,   dace% pio)
    call p_bcast (passive_v,   dace% pio)
    call p_bcast (passive_h,   dace% pio)
    call p_bcast (vmainpl,     dace% pio)
    call p_bcast (wmainpl,     dace% pio)
    call p_bcast (dz_coloc,    dace% pio)
    call p_bcast (dzv_coloc,   dace% pio)
    call p_bcast (dx_coloc,    dace% pio)
    call p_bcast (rtl_coloc,   dace% pio)
    call p_bcast (split,       dace% pio)
    call p_bcast (std_t0miss,  dace% pio)
    call p_bcast (rm_t0c_ff0,  dace% pio)
    call p_bcast (prt_data,    dace% pio)
    call p_bcast (corme_mode,  dace% pio)
    call p_bcast (check_cons,  dace% pio)
    call p_bcast (thin_bufr,   dace% pio)
    call p_bcast (layer_avg,   dace% pio)
    call p_bcast (supob_thr,   dace% pio)
    call p_bcast (prefer,      dace% pio)
    call p_bcast (avl_levs,    dace% pio)
    call p_bcast (avl_incr,    dace% pio)
    call p_bcast (dp_asurf,    dace% pio)
    call p_bcast (fgchk_inv_t, dace% pio)
    call p_bcast (fgchk_inv_q, dace% pio)
    call p_bcast (satvp_form,  dace% pio)
    !---------------------------------------
    ! consistently set level selection flags
    !---------------------------------------
    passive_t = ior (passive_t, active_t)
    passive_h = ior (passive_h, active_h)
    passive_q = ior (passive_q, active_q)
    passive_v = ior (passive_v, active_v)
    !--------------
    ! sanity checks
    !--------------
    thin_bufr = max (thin_bufr, 0)
    if (thin_bufr(2) == 0) thin_bufr(2) = thin_bufr(1)
    if (prefer <= 0 .or. prefer > 2)                     &
         call finish ('read_temp_nml','prefer not 1 or 2')
    if (satvp_form < 1 .or. satvp_form > 2)                  &
         call finish ('read_temp_nml','satvp_form not 1 or 2')
    !--------------------------
    ! print namelist /TEMP_OBS/
    !--------------------------
    if (dace% lpio) then
      write(6,'(a)') repeat('-',79)
      write(6,'()')
      write(6,'(a)') '  namelist /TEMP_OBS/'
      write(6,'()')
      write(6,'(a,100f8.0)') 'vmainpl    =',pack(vmainpl,vmainpl>=0._wp)/100.
      write(6,'(a,100f8.0)') 'wmainpl    =',pack(wmainpl,wmainpl>=0._wp)/100.
      write(6,'(a,2i5,a)')   'thin_bufr  =',thin_bufr,'  thinning interval for BUFR TEMPs (s)'
      write(6,'(a,f10.2,a)') 'top_p      =',top_p/100,'  top level for radiosonde profiles   (hPa)'
      write(6,'(a,f10.2,a)') 'top_q      =',top_q,    '  top level for humidity observations (hPa)'
      write(6,'(a,i8,4x,a)') 'top_q_mode =',top_q_mode, 'humidity obs. handling above tropopause'
      write(6,'(a,i8,4x,a)') 'satvp_form =',satvp_form, 'saturation vapor pressure formula'
      write(6,'(a,f10.2,a)') 'dz_coloc   =',dz_coloc, '  vert. separation for collocations (Pa)'
      write(6,'(a,f10.2,a)') 'dzv_coloc  =',dz_coloc, '  vert. separation for zero wind'
      write(6,'(a,f10.2,a)') 'dx_coloc   =',dx_coloc, '  hori. separation for collocations (km)'
      write(6,'(a,l8,4x,a)') 'rtl_coloc  =',rtl_coloc,  'remove TEMP levels in collocations'
      write(6,'(a,l8,4x,a)') 'split      =',split,      'split TEMPS drifting out of grid-cell'
      write(6,'(a,l8,4x,a)') 'std_t0miss =',std_t0miss, 'use std. levels with zero/miss. time'
      write(6,'(a,l8,4x,a)') 'rm_t0c_ff0 =',rm_t0c_ff0, 'dismiss levels with apparent miss. values'
      write(6,'(a,i8,4x,a)') 'corme_mode =',corme_mode, 'mode for correction message handling'
      write(6,'(a,i8,4x,a)') 'check_cons =',check_cons, 'internal consistency check for TEMPs'
      write(6,'(a,i8,4x,a)') 'active_t   =',active_t,   'active  temperature  level selection'
      write(6,'(a,i8,4x,a)') 'active_h   =',active_h,   'active  height       level selection'
      write(6,'(a,i8,4x,a)') 'active_q   =',active_q,   'active  humidity     level selection'
      write(6,'(a,i8,4x,a)') 'active_v   =',active_v,   'active  wind         level selection'
      write(6,'(a,i8,4x,a)') 'active_v_wp=',active_v_wp,'active  windprofiler level selection'
      write(6,'(a,i8,4x,a)') 'passive_t  =',passive_t,  'passive temperature  level selection'
      write(6,'(a,i8,4x,a)') 'passive_h  =',passive_h,  'passive height       level selection'
      write(6,'(a,i8,4x,a)') 'passive_q  =',passive_q,  'passive humidity     level selection'
      write(6,'(a,i8,4x,a)') 'passive_v  =',passive_v,  'passive wind         level selection'
      write(6,'(a,l8,4x,a)') 'layer_avg  =',layer_avg,  'apply layer averaging / superobbing'
      write(6,'(a,2i5,2x,a)')'supob_thr  =',supob_thr,  'no. levels required for superobbing'
      write(6,'(a,i8,4x,a)') 'prefer     =',prefer,     'level preference (1=std, 2=layer-avg.)'
      write(6,'(a,99f8.1)')  'avl_levs   = ',pack (avl_levs, avl_levs>0._wp)/100._wp
      write(6,'(a,99f8.1)')  'avl_incr   = ',pack (avl_incr, avl_incr>0._wp)/100._wp
      write(6,'(a,f9.2,a)')  'dp_asurf   = ',dp_asurf/100._wp,                                &
                                                    '  dist. lowest layer above surface (hPa)'
      write(6,'(a,f10.2,a)') 'fgchk_inv_t=',fgchk_inv_t,                                      &
                                            '  stability-dep. factor for fg-check limit for T'
      write(6,'(a,f10.2,a)') 'fgchk_inv_q=',fgchk_inv_q,                                      &
                                            '  stability-dep. factor for fg-check limit for RH'
      write(6,'()')
    endif
  end subroutine read_temp_nml
!==============================================================================
  subroutine read_temp_netcdf (ifile, i_source, obs, lkeep, nkeep)
  integer       ,intent(in)           :: ifile     ! Number of netCDF file read
  integer       ,intent(inout)        :: i_source  ! number of records in source-file
  type (t_obs)  ,intent(inout)        :: obs       ! observations data type to set
  logical       ,intent(out)          :: lkeep     ! accept observation ?
  integer       ,intent(out)          :: nkeep     ! number of accepted obsvs.

  !============================================================================
  ! Read TEMP/PILOT/WINDPROFILER observations from netCDF (converted BUFR data)
  !============================================================================
  !-----------------------------------------
  ! quantities derived from the BUFR message
  !-----------------------------------------
  integer ,allocatable      :: nrara   (:)    ! RADIOSONDE TYPE/SYSTEM
  integer ,allocatable      :: nsr     (:)    ! SOLAR AND INFRARED RADIATION CORRECTION
  integer ,allocatable      :: nsasa   (:)    ! TRACKING TECHNIQUE,STATUS OF SYSTEM
  integer ,allocatable      :: na4     (:)    ! TYPE OF MEASURING EQUIPMENT USED
  integer ,allocatable      :: mtisi   (:)    ! TIME SIGNIFICANCE
  real    ,allocatable      :: mhosnn  (:)    ! HEIGHT OF STATION GROUND ABOVE MEAN SEA
  real    ,allocatable      :: mhobnn  (:)    ! HEIGHT OF BAROMETER ABOVE MEAN SEA LEVEL
  integer ,allocatable      :: mh      (:)    ! HEIGHT [ m ]
  integer ,allocatable      :: mhp     (:)    ! HEIGHT OF STATION [ m ]
  integer ,allocatable      :: mseqm   (:)    ! STATION ELEVATION QUALITY MARK
  integer ,allocatable      :: mvtsu   (:)    ! VERTICAL SIGNIFICANCE (SURFACE OBSERV.)
  integer ,allocatable      :: mnh     (:)    ! CLOUD AMOUNT
  real    ,allocatable      :: nh      (:)    ! HEIGHT OF BASE OF CLOUD
  integer ,allocatable      :: mcc     (:)    ! CODE_TABLE FOR CLOUD TYPE
  integer ,allocatable      :: mcc0    (:)    ! CODE_TABLE FOR CLOUD TYPE
  integer ,allocatable      :: mcc1    (:)    ! CODE_TABLE FOR CLOUD TYPE
  integer ,allocatable      :: mvtsu0  (:)    ! VERTICAL SIGNIFICANCE
  real    ,allocatable      :: mtn00   (:)    ! SEA/WATER TEMPERATURE
! Extensions for windprofiler
!
! integer ,allocatable      :: mtyoan (:)     ! TYPE OF ANTENNA
! real    ,allocatable      :: m3dbbw (:)     ! 3 DB BEAMWIDTH
! real    ,allocatable      :: mmfreq (:)     ! MEAN FREQUENCY
! integer ,allocatable      :: mrglen (:)     ! RANGE-GATE LENGTH
! integer ,allocatable      :: mmsest (:)     ! MEAN SPEED ESTIMATION
! integer ,allocatable      :: mwce   (:)     ! WIND COMPUTATION ENHANCEMENT
! integer ,allocatable      :: mwps   (:)     ! NOAA WIND PROFILER, SUBMODE INFORMATION
! integer ,allocatable      :: msetp  (:)     ! TIME PERIOD OR DISPLACEMENT
! loop 000
  integer ,allocatable      :: medre   (:)    ! EXTENDED DELAYED DESCRIPT.REPLIC.FACTOR
  integer ,allocatable      :: mdrep   (:)    ! DELAYED DESCRIPT.REPLIC.FACTOR
  integer ,allocatable      :: nltpd   (:,:)  ! LONG TIME PERIOD OR DISPLACEMENT
  integer ,allocatable      :: mevss   (:,:)  ! EXTENDED VERTICAL SOUNDING SIGNIFICANCE
  integer ,allocatable      :: mvss    (:,:)  !          VERTICAL SOUNDING SIGNIFICANCE
  real    ,allocatable      :: mpn     (:,:)  ! PRESSURE (VERT.LOCATION)
  real    ,allocatable      :: nhhhn   (:,:)  ! GEOPOTENTIAL HEIGHT [ GPM ]; PILOT NHHH
  real    ,allocatable      :: mladh   (:,:)  ! LATITUDE DISPLACEMENT  (HIGH ACCURACY)
  real    ,allocatable      :: mlodh   (:,:)  ! LONGITUDE DISPLACEMENT (HIGH ACCURACY)
  real    ,allocatable      :: mtdbt   (:,:)  ! TEMPERATURE/DRY BULB TEMPERATURE
  real    ,allocatable      :: mtdnh   (:,:)  ! DEW-POINT TEMPERATURE
  real    ,allocatable      :: ndndn   (:,:)  ! WIND DIRECTION
  real    ,allocatable      :: nfnfn   (:,:)  ! WIND SPEED [ M/S ]
  real    ,allocatable      :: mu      (:,:)  ! u-component [ M/S ]
  real    ,allocatable      :: mv      (:,:)  ! v-component [ M/S ]
! integer ,allocatable      :: mquinz  (:,:)  ! WINDPROFILER QUALITY INFORMATION
! integer ,allocatable      :: mhg     (:,:)  ! WINDPROFILER HEIGHT
! integer ,allocatable      :: maddf   (:,:)  ! ASSOCIATED FIELD SIGNIFICANCE
! real    ,allocatable      :: mwmps   (:,:)  ! W-COMPONENT IN M/S
  real    ,allocatable      :: nsinor  (:,:)  ! SIGNAL TO NOISE RATIO
! integer ,allocatable      :: nwpm    (:,:)  ! WIND PROFILER, MODE INFORMATION
! integer ,allocatable      :: nwpq    (:,:)  ! NOAA WIND PROFILER, QUAL. CONTR. RESULTS
  real    ,allocatable      :: nstdff  (:,:)  ! STANDARD DEVIATION WIND SPEED
! real    ,allocatable      :: nstdfv  (:,:)  ! STANDARD DEVIATION VERTICAL WIND SPEED

! integer ,allocatable      :: mdrep   (:)    ! DELAYED DESCRIPTOR REPLICATION FACTOR
! loop 001
! integer ,allocatable      :: nltpd0  (:,:)  ! LONG TIME PERIOD OR DISPLACEMENT
! integer ,allocatable      :: mevss0  (:,:)  ! EXTENDED VERTICAL SOUNDING SIGNIFICANCE
! real    ,allocatable      :: mpn0    (:,:)  ! PRESSURE (VERT.LOCATION)
! real    ,allocatable      :: nhhh0   (:,:)  ! GEOPOTENTIAL HEIGHT [ GPM ]             ; PILOT
! real    ,allocatable      :: mladh0  (:,:)  ! LATITUDE DISPLACEMENT (HIGH ACCURACY)
! real    ,allocatable      :: mlodh0  (:,:)  ! LONGITUDE DISPLACEMENT (HIGH ACCURACY)
! real    ,allocatable      :: nvbvb   (:,:)  ! ABSOLUT WIND SHEAR IN 1 KM LAYER BELOW
! real    ,allocatable      :: nvava   (:,:)  ! ABSOLUT WIND SHEAR IN 1 KM LAYER ABOVE

  !---------------------------------------
  ! NetCDF variable IDs for body   section
  !---------------------------------------
  ! variable ID's in NetCDF file for
  !    expansion  in NetCDF file BUFR- data section4

  integer              :: varid_NRARA    ! RADIOSONDE TYPE/SYSTEM
  integer              :: varid_NSR      ! SOLAR AND INFRARED RADIATION CORRECTION
  integer              :: varid_NSASA    ! TRACKING TECHNIQUE,STATUS OF SYSTEM
  integer              :: varid_NA4      ! TYPE OF MEASURING EQUIPMENT USED
  integer              :: varid_MTISI    ! TIME SIGNIFICANCE
  integer              :: varid_MHOSNN   ! HEIGHT OF STATION GROUND ABOVE MEAN SEA
  integer              :: varid_MHOBNN   ! HEIGHT OF BAROMETER ABOVE MEAN SEA LEVEL
  integer              :: varid_mh       ! HEIGHT [ m ]
  integer              :: varid_MHP      ! HEIGHT OF STATION [ m ]
  integer              :: varid_MSEQM    ! STATION ELEVATION QUALITY MARK
  integer              :: varid_MVTSU    ! VERTICAL SIGNIFICANCE (SURFACE OBSERV.)
  integer              :: varid_MNH      ! CLOUD AMOUNT
  integer              :: varid_NH       ! HEIGHT OF BASE OF CLOUD
  integer              :: varid_MCC      ! CODE_TABLE FOR CLOUD TYPE
  integer              :: varid_MCC0     ! CODE_TABLE FOR CLOUD TYPE
  integer              :: varid_MCC1     ! CODE_TABLE FOR CLOUD TYPE
  integer              :: varid_MVTSU0   ! VERTICAL SIGNIFICANCE (SURFACE OBSERV.)
  integer              :: varid_MTN00    ! SEA/WATER TEMPERATURE
! Extensions for windprofiler
!
! integer              :: varid_MTYOAN   ! TYPE OF ANTENNA
! integer              :: varid_M3DBBW   ! 3 DB BEAMWIDTH
! integer              :: varid_MMFREQ   ! MEAN FREQUENCY
! integer              :: varid_MRGLEN   ! RANGE-GATE LENGTH
! integer              :: varid_MMSEST   ! MEAN SPEED ESTIMATION
! integer              :: varid_MWCE     ! WIND COMPUTATION ENHANCEMENT
! integer              :: varid_MWPS     ! NOAA WIND PROFILER, SUBMODE INFORMATION
! integer              :: varid_MSETP    ! TIME PERIOD OR DISPLACEMENT
! loop 000
  integer              :: varid_MEDRE    ! EXTENDED DELAYED DESCRIPT.REPLIC.FACTOR
  integer              :: varid_MDREP    ! DELAYED DESCRIPT.REPLIC.FACTOR
  integer              :: varid_NLTPD    ! LONG TIME PERIOD OR DISPLACEMENT
  integer              :: varid_MEVSS    ! EXTENDED VERTICAL SOUNDING SIGNIFICANCE
  integer              :: varid_MPN      ! PRESSURE (VERT.LOCATION) [Pa]
  integer              :: varid_NHHHN    ! GEOPOTENTIAL HEIGHT [ GPM ]; PILOT NHHH
  integer              :: varid_MLADH    ! LATITUDE DISPLACEMENT (HIGH ACCURACY)
  integer              :: varid_MLODH    ! LONGITUDE DISPLACEMENT(HIGH  ACCURACY)
  integer              :: varid_MTDBT    ! TEMPERATURE/DRY BULB TEMPERATURE
  integer              :: varid_MTDNH    ! DEW-POINT TEMPERATURE
  integer              :: varid_NDNDN    ! WIND DIRECTION
  integer              :: varid_NFNFN    ! WIND SPEED [ M/S ]
  integer              :: varid_MU       ! u-component [ M/S ]
  integer              :: varid_MV       ! v-component [ M/S ]
! integer              :: varid_MQINZ    ! WINDPROFILER QUALITY INFORMATION
  integer              :: varid_MHG      ! WINDPROFILER HEIGHT
  integer              :: varid_MHOSEN   ! SCADA Wind power turbine HEIGHT
! integer              :: varid_MADDF    ! ASSOCIATED FIELD SIGNIFICANCE
! integer              :: varid_MWMPS    ! W-COMPONENT IN M/S
  integer              :: varid_NSINOR   ! SIGNAL TO NOISE RATIO
! integer              :: varid_NWPM     ! WIND PROFILER, MODE INFORMATION
! integer              :: varid_NWPQ     ! NOAA WIND PROFILER, QUAL. CONTR. RESULTS
  integer              :: varid_NSTDFF   ! STANDARD DEVIATION WIND SPEED
! integer              :: varid_NSTDFV   ! STANDARD DEVIATION VERTICAL WIND SPEED


! integer              :: varid_MDREP    ! DELAYED DESCRIPTOR REPLICATION FACTOR
! loop 001
! integer              :: varid_NLTPD0   ! LONG TIME PERIOD OR DISPLACEMENT
! integer              :: varid_MEVSS0   ! EXTENDED VERTICAL SOUNDING SIGNIFICANCE
! integer              :: varid_MPN0     ! PRESSURE (VERT.LOCATION)
! integer              :: varid_NHHH0    ! GEOPOTENTIAL HEIGHT [ GPM ]             ; PILOT
! integer              :: varid_MLADH0   ! LATITUDE DISPLACEMENT (HIGH ACCURACY)
! integer              :: varid_MLODH0   ! LONGITUDE DISPLACEMENT (HIGH ACCURACY)
! integer              :: varid_NVBVB    ! ABSOLUT WIND SHEAR IN 1 KM LAYER BELOW
! integer              :: varid_NVAVA    ! ABSOLUT WIND SHEAR IN 1 KM LAYER ABOVE

  !================
  ! local variables
  !================
  integer                :: nlev     ! number of levels
  type (t_temp) ,pointer :: btemp(:) !
  integer                :: novalid  ! number of observations actually stored

  type (t_spot)        :: spt0, spti ! observation meta data
  type (t_spot) ,save  :: empty      !
  type (t_use)         :: use        ! status variable


  type (t_head)        :: head         ! meta information

  integer              :: bufr_type    ! BUFR message type    read
  integer              :: bufr_subtype ! BUFR message subtype read
  integer              :: centre       ! generating centre

  integer              :: obstype      ! observation type
  integer              :: report_subt  ! observation subtype (Datenbankkennz.)
  integer              :: report_subti ! observation subtype index
  type (t_obsid)       :: obsid        ! observation id table entry

  ! variable for     NetCDF file concerning unlimited dimension, Maximum number of dimensions,
  !                                         number of attributes
  integer              :: status         ! NetCDF status variable
! integer              :: ncdims         ! NetCDF number of dimensions defined in this NetCDF file
! character (LEN=40)   :: yname          ! NetCDF dimension name
  integer              :: ncvars         ! NetCDF number of variables  defined in this NetCDF file
! integer              :: ncatts         ! NetCDF number of attributes defined in this NetCDF file
! integer              :: unlim_dimid    ! NetCDF id for of unlimited dimension defined in this NetCDF file
  integer              :: numDims        ! NetCDF number of dimensions for individual variable in NetCDF variable
  integer              :: numAtts        ! NetCDF number of attributes for individual variable in NetCDF variable
  integer              :: dimid_level    ! NetCDF dimension id for  Levels
  integer              :: dimid_report   !                     for  Reports

  integer              :: len_level      ! number of vertical levels in NetCDF-File 1.type of field
  integer              :: len_report     ! number of reports in NetCDF file         1.type of field
  integer              :: len_level0     ! number of vertical levels in NetCDF-File 2.type of field
  integer              :: len_report0    ! number of reports in NetCDF file         2.type of field
  integer              :: j              ! loop index
  integer              :: ilev           ! loop index vertical levels
  integer              :: nc1            ! first  dimension for netcdf 1.type of field getvar in start / count
  integer              :: nc2            ! second dimension for netcdf 1.type of field getvar in start / count
! integer              :: nc3            ! first  dimension for netcdf 2.type of field getvar in start / count
! integer              :: nc4            ! second dimension for netcdf 2.type of field getvar in start / count
  integer              :: entry1,entry   ! position in source file (subset)
  logical              :: lk             ! lkeep, local copy
  character(len=5)     :: dim1_length    ! variable length descriptor
  logical              :: l_newbufr      ! "new" BUFR format?
  integer              :: dimids (dimids_max)

  integer ,allocatable :: ifield  (:,:)  !
  real    ,allocatable :: rfield  (:,:)  !
! integer ,allocatable :: ifield0 (:,:)  !
! real    ,allocatable :: rfield0 (:,:)  !
  integer ,allocatable :: ifield1 (:)    !
  real    ,allocatable :: rfield1 (:)    !

  !----------------
  ! index variables
  !----------------
! integer       :: nreport          ! number of observations (default)

  integer       :: is               ! sub-set    index
! integer       :: ie               ! entry in sub-set
! integer       :: id               ! descriptor index
  integer       :: i                ! spot index
! integer       :: j1               ! spot index

  logical       :: lpr_extd         ! extended  printing of temps
  integer       :: npr_extd         ! number of extended  printing of temps

! logical for meteorological variables(over all reports)
  logical         :: l_mpn   , l_height, l_wind  , l_mtdbt, l_mtdnh,   &
                     l_mladh , l_mlodh , l_mevss , l_nltpd,            &
                     l_ndndn , l_nfnfn , l_nhhhn , l_nhhh, l_mhg, l_nstdff, &
                     l_mhosen, l_mu    , l_mv    , l_uv
! logical for technological variables(over all reports)
  logical         :: l_nrara , l_nsr , l_nsasa, l_na4   , l_mtisi, l_mhosnn,&
                     l_mhobnn, l_mh  , l_mseqm, l_mvtsu , l_mnh  , l_nh    ,&
                     l_mcc   , l_mcc0, l_mcc1 , l_mvtsu0, l_mtn00, l_medre ,&
                     l_mhp   , l_mdrep, l_nsinor

  character(NF90_MAX_NAME)   :: yname_v  ! NetCDF variable name

  !----------------------------------------------
  ! Misc. local variables
  ! used by contained subroutine 'vector_indices'
  !----------------------------------------------
  integer              :: n_idx
  integer, allocatable :: idx(:)         ! Auxiliary arrays, used in
  logical, allocatable :: mask(:)        ! checks for double entries
!------------------------------------------------------------------------------
  lpr_temp = .false.; if (netcdf_verb > 1) lpr_temp = .true.
  lpr_extd = .true.
  npr_extd =   2
! if( lpr_extd) npr_extd  =   10
  if( lpr_extd) npr_extd  =  300

  !------------------------------
  ! get default number of reports
  !------------------------------
  nkeep       = 0
  !---------------------------
  ! get variables ids and name
  !---------------------------
  status = nf90_Inquire (ncid, nVariables=ncvars)
  if ( lpr_temp .and. lpr_extd ) then
    do j = 1 , ncvars
    status = nf90_Inquire_Variable(ncid, j,name=yname_v )
    write (6,'(a,i3,a,a)') &
         'nf90_Inquire_Variable(',j,')  Variable name(o): ',trim(yname_v)
    end do
  endif
  !------------------------
  ! get dimension of fields
  ! from wind direction or
  ! from u-component
  !------------------------
  nc1 = 0
  nc2 = 0
  len_level   = 0
  len_report  = 0
  dim1_length = ""

  status = nf90_inq_varid (ncid, 'NDNDN', varid_NDNDN)
  if (status /= NF90_NOERR) then
    status = nf90_inq_varid (ncid, 'MU', varid_NDNDN)
  end if
  if (status == NF90_NOERR) then
    status = nf90_Inquire_Variable (ncid, varid_NDNDN,                         &
                                    ndims=numDims, dimids=dimids, natts=numAtts)
    status = nf90_get_att (ncid, varid_NDNDN, 'dim1_length', dim1_length)
  else
    !------------------------------------------------
    ! Fallback: just get number of reports from MEDRE
    !------------------------------------------------
    status = nf90_inq_varid (ncid, 'MEDRE', varid_MEDRE)
    if (status == NF90_NOERR) dim1_length = 'MEDRE'
    status = nf90_Inquire_Variable  (ncid, varid_MEDRE,          &
                                     ndims=numDims, dimids=dimids)
    status = nf90_Inquire_Dimension (ncid, dimids(1), len=len_report)
    nc2 = len_report
  end if
  if (status /= NF90_NOERR) then
    call finish('read_temp_netcdf','NDNDN/MU/MEDRE all not present')
  end if

  if(numDims >= 2) then
    dimid_level  = dimids(1)
    dimid_report = dimids(2)
    status = nf90_Inquire_Dimension(ncid, dimid_level,  len=len_level)
    status = nf90_Inquire_Dimension(ncid, dimid_report, len=len_report)
    nc1       = len_level
    nc2       = len_report
  endif

! status = nf90_inq_varid (ncid, 'NVBVB' ,  varid_NVBVB)
! status = nf90_Inquire_Variable(ncid, varid_NVBVB, ndims=numDims, dimids=dimids, natts=numAtts)

  len_level0  = 0
  len_report0 = 0

  if(numDims >= 2) then
    dimid_level  = dimids(1)
    dimid_report = dimids(2)
    status = nf90_Inquire_Dimension(ncid, dimid_level,  len=len_level0)
    status = nf90_Inquire_Dimension(ncid, dimid_report, len=len_report0)
  endif

  if ( lpr_temp  ) then
    write (6,'()')
    write (6,'(a,i3,a)' )   'pe=',dace% pe, ' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write (6,'(a,i3,a)' )   'pe=',dace% pe, ' BEGINN  MO_TEMP.F90     !!!!!!!!!!!!!!!!!!!!!!!'
    write (6,'(a,i3,a)' )   'pe=',dace% pe, ' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write (6,'(a,i4,/,6x,a,i4,/,6x,a,i4,/,6x,a,10i4)')   'pe=',dace% pe,            &
                            ' varid_NDNDN number of dimensions: ',numDims,          &
                            ' varid_NDNDN number of attributes: ',numAtts,          &
                            ' varid_NDNDN ids    of dimensions: ',dimids(1:numDims)
    write (6,'(a,i3,a,i8,/ , &
             & 1x  ,a,i8,/ , &
             & 1x  ,a,i10,a,i10) ') &
                        'pe=',dace% pe,  ' Reports in BUFR File:',len_report ,      &
                                         ' Levels  in BUFR File:',len_level  ,      &
                                         ' nc1=',nc1,' nc2=',nc2
  endif
  !----------------------------------------
  ! define number of reports in source-file
  !----------------------------------------
  i_source = len_report

  allocate (nrara   (len_report))
  allocate (nsr     (len_report))
  allocate (nsasa   (len_report))
  allocate (na4     (len_report))
  allocate (mtisi   (len_report))
  allocate (mhosnn  (len_report))
  allocate (mhobnn  (len_report))
  allocate (mh      (len_report))
  allocate (mhp     (len_report))
  allocate (mseqm   (len_report))
  allocate (mvtsu   (len_report))
  allocate (mnh     (len_report))
  allocate (nh      (len_report))
  allocate (mcc     (len_report))
  allocate (mcc0    (len_report))
  allocate (mcc1    (len_report))
  allocate (mvtsu0  (len_report))
  allocate (mtn00   (len_report))
! allocate (mtyoan  (len_report))
! allocate (m3dbbw  (len_report))
! allocate (mmfreq  (len_report))
! allocate (mrglen  (len_report))
! allocate (mmsest  (len_report))
! allocate (mwce    (len_report))
! allocate (mwps    (len_report))
! allocate (msetp   (len_report))
! loop 000
  allocate (medre   (len_report))
  allocate (mdrep   (len_report))

  allocate (nltpd   (len_level,len_report))
  allocate (mevss   (len_level,len_report))
  allocate (mvss    (len_level,len_report))
  allocate (mpn     (len_level,len_report))
  allocate (nhhhn   (len_level,len_report))
  allocate (mladh   (len_level,len_report))
  allocate (mlodh   (len_level,len_report))
  allocate (mtdbt   (len_level,len_report))
  allocate (mtdnh   (len_level,len_report))
  allocate (ndndn   (len_level,len_report))
  allocate (nfnfn   (len_level,len_report))
  allocate (mu      (len_level,len_report))
  allocate (mv      (len_level,len_report))
! allocate (mquinz  (len_level,len_report))
! allocate (mhg     (len_level,len_report))
! allocate (maddf   (len_level,len_report))
! allocate (mwmps   (len_level,len_report))
  allocate (nsinor  (len_level,len_report))
! allocate (nwpm    (len_level,len_report))
! allocate (nwpq    (len_level,len_report))
  allocate (nstdff  (len_level,len_report))
! allocate (nstdfv  (len_level,len_report))

! loop 001
! allocate (mdrep   (len_report0))

! allocate (nltpd0  (len_level0,len_report0))
! allocate (mevss0  (len_level0,len_report0))
! allocate (mpn0    (len_level0,len_report0))
! allocate (nhhh0   (len_level0,len_report0))
! allocate (mladh0  (len_level0,len_report0))
! allocate (mlodh0  (len_level0,len_report0))
! allocate (nvbvb   (len_level0,len_report0))
! allocate (nvava   (len_level0,len_report0))


  allocate (ifield  (len_level ,len_report))
  allocate (rfield  (len_level ,len_report))
! allocate (ifield0 (len_level0,len_report0))
! allocate (rfield0 (len_level0,len_report0))
  allocate (ifield1 (           len_report ))
  allocate (rfield1 (           len_report ))

  !-------------------------------
  ! get technological information
  !-------------------------------

  !-----------------------
  ! RADIOSONDE TYPE/SYSTEM
  !-----------------------
  nrara = -1
  status = nf90_inq_varid (ncid, 'NRARA'  ,  varid_NRARA )
  l_nrara = .FALSE.
  if (status == nf90_noerr) then
    l_nrara = .TRUE.
    status = nf90_get_var (ncid, varid_NRARA , ifield1)

!   check for missing values
    if ( all( ifield1 == imissing ))  then
       l_nrara = .FALSE.
    else
       nrara = ifield1
       where ( ifield1 == imissing ) nrara = -1
    endif
  endif

  !-----------------------------------------
  ! SOLAR AND INFRARED RADIATION CORRECTION
  ! TRACKING TECHNIQUE,STATUS OF SYSTEM
  !-----------------------------------------
  nsr   = -1
  status = nf90_inq_varid (ncid, 'NSR'  ,  varid_NSR )
  l_nsr = .FALSE.
  if (status == nf90_noerr) then
    l_nsr = .TRUE.
    status = nf90_get_var (ncid, varid_NSR , ifield1)

!   check for missing values
    if ( all( ifield1 == imissing ))  then
       l_nsr = .FALSE.
    else
       nsr = ifield1
       where ( ifield1 == imissing ) nsr = -1
    endif
  endif

  nsasa   = -1
  status = nf90_inq_varid (ncid, 'NSASA'  ,  varid_NSASA )
  l_nsasa = .FALSE.
  if (status == nf90_noerr) then
    l_nsasa = .TRUE.
    status = nf90_get_var (ncid, varid_NSASA , ifield1)

!   check for missing values
    if ( all( ifield1 == imissing ))  then
       l_nsasa = .FALSE.
    else
       nsasa = ifield1
       where ( ifield1 == imissing ) nsasa = -1
    endif
  endif
  !---------------------------------
  ! TYPE OF MEASURING EQUIPMENT USED
  ! TIME SIGNIFICANCE
  !---------------------------------
  na4   = -1
  status = nf90_inq_varid (ncid, 'NA4'  ,  varid_NA4 )
  l_na4 = .FALSE.
  if (status == nf90_noerr) then
    l_na4 = .TRUE.
    status = nf90_get_var (ncid, varid_NA4 , ifield1)

!   check for missing values
    if ( all( ifield1 == imissing ))  then
       l_na4 = .FALSE.
    else
       na4 = ifield1
       where ( ifield1 == imissing ) na4 = -1
    endif
  endif

  mtisi   = -1
  status = nf90_inq_varid (ncid, 'MTISI'  ,  varid_MTISI )
  l_mtisi = .FALSE.
  if (status == nf90_noerr) then
    l_mtisi = .TRUE.
    status = nf90_get_var (ncid, varid_MTISI , ifield1)

!   check for missing values
    if ( all( ifield1 == imissing ))  then
       l_mtisi = .FALSE.
    else
       mtisi = ifield1
       where ( ifield1 == imissing ) mtisi = -1
    endif
  endif

  !-----------------------------------------
  ! HEIGHT OF STATION GROUND ABOVE MEAN SEA
  ! HEIGHT OF BAROMETER ABOVE MEAN SEA LEVEL
  !-----------------------------------------
  mhosnn   = -999.
  status = nf90_inq_varid (ncid, 'MHOSNN'  ,  varid_MHOSNN )
  l_mhosnn = .FALSE.
  if (status == nf90_noerr) then
    l_mhosnn = .TRUE.
    status = nf90_get_var (ncid, varid_MHOSNN , rfield1)

!   check for missing values
    if ( all( rfield1 == rmissing ))  then
       l_mhosnn = .FALSE.
    else
       mhosnn = rfield1
       where ( rfield1 == rmissing ) mhosnn = -999.
    endif
  endif

  mhobnn   = -999.
  status = nf90_inq_varid (ncid, 'MHOBNN'  ,  varid_MHOBNN )
  l_mhobnn = .FALSE.
  if (status == nf90_noerr) then
    l_mhobnn = .TRUE.
    status = nf90_get_var (ncid, varid_MHOBNN , rfield1)

!   check for missing values
    if ( all( rfield1 == rmissing ))  then
       l_mhobnn = .FALSE.
    else
       mhobnn = rfield1
       where ( rfield1 == rmissing ) mhobnn = -999.
    endif
  endif

  !--------------------------------
  ! HEIGHT [ m ]
  ! STATION ELEVATION QUALITY MARK
  !--------------------------------

  if (s2ikz(1) /= 10553 .and. s2ikz(1) /= 10600) then
    mh   = -999
    status = nf90_inq_varid (ncid, 'MH'  ,  varid_MH )
    l_mh = .FALSE.
    if (status == nf90_noerr) then
      l_mh = .TRUE.
      status = nf90_get_var (ncid, varid_MH , ifield1)

!   check for missing values
      if ( all( ifield1 == imissing ))  then
        l_mh = .FALSE.
      else
        mh = ifield1
        where ( ifield1 == imissing ) mh = -999
      endif
    endif
  endif

  mhp  = -999
  status = nf90_inq_varid (ncid, 'MHP'  ,  varid_MHP )
  l_mhp = .FALSE.
  if (status == nf90_noerr) then
    l_mhp = .TRUE.
    status = nf90_get_var (ncid, varid_MHP , ifield1)

!   check for missing values
    if ( all( ifield1 == imissing ))  then
      l_mhp = .FALSE.
    else
      mhp = ifield1
      where ( ifield1 == imissing ) mhp = -999
    endif
  endif

  mseqm   = -1
  status = nf90_inq_varid (ncid, 'MSEQM'  ,  varid_MSEQM )
  l_mseqm = .FALSE.
  if (status == nf90_noerr) then
    l_mseqm = .TRUE.
    status = nf90_get_var (ncid, varid_MSEQM , ifield1)

!   check for missing values
    if ( all( ifield1 == imissing ))  then
       l_mseqm = .FALSE.
    else
       mseqm = ifield1
       where ( ifield1 == imissing ) mseqm = -1
    endif
  endif
  !----------------------------------------
  ! VERTICAL SIGNIFICANCE (SURFACE OBSERV.)
  ! CLOUD AMOUNT
  !----------------------------------------

  mvtsu   = -1
  status = nf90_inq_varid (ncid, 'MVTSU'  ,  varid_MVTSU )
  l_mvtsu = .FALSE.
  if (status == nf90_noerr) then
    l_mvtsu = .TRUE.
    status = nf90_get_var (ncid, varid_MVTSU , ifield1)

!   check for missing values
    if ( all( ifield1 == imissing ))  then
       l_mvtsu = .FALSE.
    else
       mvtsu = ifield1
       where ( ifield1 == imissing ) mvtsu = -1
    endif
  endif

  mnh   = -1
  status = nf90_inq_varid (ncid, 'MNH'  ,  varid_MNH )
  l_mnh = .FALSE.
  if (status == nf90_noerr) then
    l_mnh = .TRUE.
    status = nf90_get_var (ncid, varid_MNH , ifield1)

!   check for missing values
    if ( all( ifield1 == imissing ))  then
       l_mnh = .FALSE.
    else
       mnh = ifield1
       where ( ifield1 == imissing ) mnh = -1
    endif
  endif

  !--------------------------
  ! HEIGHT OF BASE OF CLOUD
  ! CODE_TABLE FOR CLOUD TYPE
  !--------------------------

  nh   = -1
  status = nf90_inq_varid (ncid, 'NH'  ,  varid_NH )
  l_nh = .FALSE.
  if (status == nf90_noerr) then
    l_nh = .TRUE.
    status = nf90_get_var (ncid, varid_NH , rfield1)

!   check for missing values
    if ( all( rfield1 == rmissing ))  then
       l_nh = .FALSE.
    else
       nh = rfield1
       where ( rfield1 == rmissing ) nh = -1.
    endif
  endif

  mcc   = -1
  status = nf90_inq_varid (ncid, 'MCC'  ,  varid_MCC )
  l_mcc = .FALSE.
  if (status == nf90_noerr) then
    l_mcc = .TRUE.
    status = nf90_get_var (ncid, varid_MCC , ifield1)

!   check for missing values
    if ( all( ifield1 == imissing ))  then
       l_mcc = .FALSE.
    else
       mcc = ifield1
       where ( ifield1 == imissing ) mcc = -1
    endif
  endif

  mcc0   = -1
  status = nf90_inq_varid (ncid, 'MCC0'  ,  varid_MCC0 )
  l_mcc0 = .FALSE.
  if (status == nf90_noerr) then
    l_mcc0 = .TRUE.
    status = nf90_get_var (ncid, varid_MCC0 , ifield1)

!   check for missing values
    if ( all( ifield1 == imissing ))  then
       l_mcc0 = .FALSE.
    else
       mcc0 = ifield1
       where ( ifield1 == imissing ) mcc0 = -1
    endif
  endif

  mcc1   = -1
  status = nf90_inq_varid (ncid, 'MCC1'  ,  varid_MCC1 )
  l_mcc1 = .FALSE.
  if (status == nf90_noerr) then
    l_mcc1 = .TRUE.
    status = nf90_get_var (ncid, varid_MCC1 , ifield1)

!   check for missing values
    if ( all( ifield1 == imissing ))  then
       l_mcc1 = .FALSE.
    else
       mcc1 = ifield1
       where ( ifield1 == imissing ) mcc1 = -1
    endif
  endif
  !----------------------
  ! VERTICAL SIGNIFICANCE
  ! SEA/WATER TEMPERATURE
  !----------------------
  mvtsu0   = -1
  status = nf90_inq_varid (ncid, 'MVTSU0'  ,  varid_MVTSU0 )
  l_mvtsu0 = .FALSE.
  if (status == nf90_noerr) then
    l_mvtsu0 = .TRUE.
    status = nf90_get_var (ncid, varid_MVTSU0 , ifield1)

!   check for missing values
    if ( all( ifield1 == imissing ))  then
       l_mvtsu0 = .FALSE.
    else
       mvtsu0 = ifield1
       where ( ifield1 == imissing ) mvtsu0 = -1
    endif
  endif

  mtn00   = -1
  status = nf90_inq_varid (ncid, 'MTN00'  ,  varid_MTN00 )
  l_mtn00 = .FALSE.
  if (status == nf90_noerr) then
    l_mtn00 = .TRUE.
    status = nf90_get_var (ncid, varid_MTN00 , rfield1)

!   check for missing values
    if ( all( rfield1 == rmissing ))  then
       l_mtn00 = .FALSE.
    else
       mtn00 = rfield1
       where ( rfield1 == rmissing ) mtn00 = -1.
    endif
  endif


! loop 000
  !----------------------------------------
  ! EXTENDED DELAYED DESCRIPT.REPLIC.FACTOR
  !----------------------------------------
  medre   = -1
  status = nf90_inq_varid (ncid, 'MEDRE'  ,  varid_MEDRE )
  l_medre = .FALSE.
  if (status == nf90_noerr) then
    l_medre = .TRUE.
    status = nf90_get_var (ncid, varid_MEDRE , ifield1)

!   check for missing values
    if ( all( ifield1 == imissing ))  then
       l_medre = .FALSE.
    else
       medre = ifield1
       where ( ifield1 == imissing ) medre = -1
    endif
  endif

  !----------------------------------------
  ! DELAYED DESCRIPT.REPLIC.FACTOR
  !----------------------------------------
  mdrep   = -1
  status = nf90_inq_varid (ncid, 'MDREP'  ,  varid_MDREP )
  l_mdrep = .FALSE.
  if (status == nf90_noerr) then
    l_mdrep = .TRUE.
    status = nf90_get_var (ncid, varid_MDREP , ifield1)

!   check for missing values
    if ( all( ifield1 == imissing ))  then
       l_mdrep = .FALSE.
    else
       mdrep = ifield1
       where ( ifield1 == imissing ) mdrep = -1
    endif
  endif

  !----------------------------------------
  ! LONG TIME PERIOD OR DISPLACEMENT
  ! EXTENDED VERTICAL SOUNDING SIGNIFICANCE
  !----------------------------------------
  nltpd   = inv_secs
  status  = nf90_inq_varid (ncid, 'NLTPD'  ,  varid_NLTPD )
  l_nltpd = .FALSE.
  if (status == nf90_noerr) then
    l_nltpd = .TRUE.
    status = nf90_get_var (ncid, varid_NLTPD , ifield)

!   check for missing values
    if ( all( ifield == imissing ))  then
       l_nltpd = .FALSE.
    else
       nltpd = ifield
       where ( ifield == imissing ) nltpd = inv_secs
    endif
  endif

! initialise as in structure t_temp
! mevss   = -1
  mevss   = lev_miss_18
  mvss    = lev_miss

  status = nf90_inq_varid (ncid, 'MEVSS'  ,  varid_MEVSS )
  l_mevss = .FALSE.
  if (status == nf90_noerr) then
    l_mevss = .TRUE.
    status = nf90_get_var (ncid, varid_MEVSS , ifield)

!   check for missing values
    if ( all( ifield == imissing ))  then
       l_mevss = .FALSE.
    else
       mevss = ifield
       where ( ifield == imissing ) mevss = lev_miss_18
    endif

!   check for all bits set (undefined)
    where (iand (mevss , lev_miss_18) /=0 ) mevss = lev_miss_18

    !--------------------------------------------
    ! transform to vertical sounding significance
    !--------------------------------------------
    mvss = vss_from_evss (mevss)

  endif

  !------------------------------
  ! PRESSURE (VERT.LOCATION) [Pa]
  !------------------------------
  mpn    = rvind
  l_mpn  = .FALSE.
  status = nf90_inq_varid (ncid, 'MPN', varid_MPN)
  if (status == nf90_noerr) then
     call get_real (mpn, 'MPN', rvind)
     l_mpn = .not. all (mpn == rvind)
  endif
  !---------------------------------------------------------
  ! GEOPOTENTIAL HEIGHT [ GPM ]; PILOT NHHH; WINDPROFILER MH
  !---------------------------------------------------------
  nhhhn    = rvind
  l_nhhhn  = .FALSE.
  l_nhhh   = .FALSE.
  l_height = .FALSE.
  l_mhg    = .FALSE.
  l_mhosen = .FALSE.
  status   = nf90_inq_varid (ncid, 'NHHHN', varid_NHHHN)
  if (status == nf90_noerr) then
     call get_real (nhhhn, 'NHHHN', rvind)
     l_nhhhn = .not. all (nhhhn == rvind)
  else
     status = nf90_inq_varid (ncid, 'NHHH', varid_NHHHN)
     if (status == nf90_noerr) then
       l_nhhh   = .TRUE.
       status = nf90_get_var (ncid, varid_NHHHN   , ifield )
   !   check for missing values
       if ( all( ifield == imissing ))  then
          l_nhhh   = .FALSE.
       else
          nhhhn = ifield
          where ( ifield == imissing ) nhhhn = rvind
      endif
    endif
  endif

  if (any (s2ikz(1) == [ 10553, 10600, 548 ])) then
    status = nf90_inq_varid (ncid, 'MH'  ,  varid_MHG  )
    if (status == nf90_noerr) then
      l_mhg = .TRUE.
      status = nf90_get_var (ncid, varid_MHG , ifield)

!     check for missing values
      if ( all( ifield == rmissing ))  then
        l_mhg = .FALSE.
      else
        nhhhn = ifield
        where ( ifield == imissing ) nhhhn = rvind
      endif
    endif
  endif

  if (s2ikz(1) == 584 ) then
    status = nf90_inq_varid (ncid, 'MHOSEN'  ,  varid_MHOSEN  )
    if (status == nf90_noerr) then
      l_mhosen = .TRUE.
      status = nf90_get_var (ncid, varid_MHOSEN , rfield)

!     check for missing values
      if ( all( rfield == rmissing ))  then
        l_mhosen = .FALSE.
      else
        nhhhn = rfield
        where ( ifield == imissing ) nhhhn = rvind
      endif
    endif
  endif

  if ( .not. (l_nhhhn .or. l_nhhh .or. l_mhg .or. l_mhosen)) then
    !--------------------------------------------
    ! geopotential height for ECMWF BUFR: 'NHNHN'
    !--------------------------------------------
    nhhhn    = rvind
    status = nf90_inq_varid (ncid, 'NHNHN'  ,  varid_NHHHN )
    l_nhhhn  = .FALSE.
    if (status == nf90_noerr) then
      l_nhhhn = .TRUE.
      status = nf90_get_var (ncid, varid_NHHHN , rfield)
      !-------------------------
      ! check for missing values
      !-------------------------
      if ( all( rfield == rmissing ))  then
         l_nhhhn = .FALSE.
      else
         where ( rfield == rmissing )
           nhhhn = rvind
         elsewhere
           nhhhn = rfield / gacc     ! -> gpm
         endwhere
      endif
    endif
  endif

  if ( l_nhhhn .or. l_nhhh .or. l_mhg .or. l_mhosen) l_height = .TRUE.
  !---------------------------------------
  ! LATITUDE DISPLACEMENT  (HIGH ACCURACY)
  ! LONGITUDE DISPLACEMENT (HIGH ACCURACY)
  !---------------------------------------
  mladh   = rvind
  status = nf90_inq_varid (ncid, 'MLADH'  ,  varid_MLADH )
  l_mladh = .FALSE.
  if (status == nf90_noerr) then
    l_mladh = .TRUE.
    status = nf90_get_var (ncid, varid_MLADH , rfield)

!   check for missing values
    if ( all( rfield == rmissing ))  then
       l_mladh = .FALSE.
    else
       mladh = rfield
       where ( rfield == rmissing ) mladh = rvind
    endif
  endif

  mlodh   = rvind
  status = nf90_inq_varid (ncid, 'MLODH'  ,  varid_MLODH )
  l_mlodh = .FALSE.
  if (status == nf90_noerr) then
    l_mlodh = .TRUE.
    status = nf90_get_var (ncid, varid_MLODH , rfield)

!   check for missing values
    if ( all( rfield == rmissing ))  then
       l_mlodh = .FALSE.
    else
       mlodh = rfield
       where ( rfield == rmissing ) mlodh = rvind
    endif
  endif

  !---------------------------------
  ! TEMPERATURE/DRY BULB TEMPERATURE
  !---------------------------------
  mtdbt   = rvind
  l_mtdbt = .FALSE.
  status = nf90_inq_varid (ncid, 'MTDBT'  ,  varid_MTDBT )
  if (status /= nf90_noerr) then
    !-----------------------
    ! try MTN for ecmwf BUFR
    !-----------------------
    status = nf90_inq_varid (ncid, 'MTN'  ,  varid_MTDBT )
  endif
  if (status == nf90_noerr) then
    l_mtdbt = .TRUE.
    status = nf90_get_var (ncid, varid_MTDBT , rfield)
    !-------------------------
    ! check for missing values
    !-------------------------
    if ( all( rfield == rmissing ))  then
       l_mtdbt = .FALSE.
    else
       mtdbt = rfield
       where ( rfield == rmissing ) mtdbt = rvind
    endif
  endif

  !----------------------
  ! DEW-POINT TEMPERATURE
  !----------------------
  mtdnh   = rvind
  l_mtdnh = .FALSE.
  status = nf90_inq_varid (ncid, 'MTDNH'  ,  varid_MTDNH )
  if (status /= nf90_noerr) then
    !------------------------
    ! try MTDN for ecmwf BUFR
    !------------------------
    status = nf90_inq_varid (ncid, 'MTDN'  ,  varid_MTDNH )
  endif
  if (status == nf90_noerr) then
    l_mtdnh = .TRUE.
    status = nf90_get_var (ncid, varid_MTDNH , rfield)
    !-------------------------
    ! check for missing values
    !-------------------------
    if ( all( rfield == rmissing ))  then
       l_mtdnh = .FALSE.
    else
       mtdnh = rfield
       where ( rfield == rmissing ) mtdnh = rvind
    endif
  endif

  !---------------------------------------
  ! WIND DIRECTION
  ! WIND SPEED [ M/S ]
  !---------------------------------------
  ndndn   = rvind
  status = nf90_inq_varid (ncid, 'NDNDN'  ,  varid_NDNDN )
  l_ndndn = .FALSE.
  if (status == nf90_noerr) then
    l_ndndn = .TRUE.
    status = nf90_get_var (ncid, varid_NDNDN , ifield)

!   check for missing values
    if ( all( ifield == imissing ))  then
       l_ndndn = .FALSE.
    else
       ndndn = ifield
       where ( ifield == imissing ) ndndn = rvind
    endif
  endif

  nfnfn   = rvind
  status = nf90_inq_varid (ncid, 'NFNFN'  ,  varid_NFNFN )
  l_nfnfn = .FALSE.
  if (status == nf90_noerr) then
    l_nfnfn = .TRUE.
    status = nf90_get_var (ncid, varid_NFNFN , rfield)

!   check for missing values
    if ( all( rfield == rmissing ))  then
       l_nfnfn = .FALSE.
    else
       nfnfn = rfield
       where ( rfield == rmissing ) nfnfn = rvind
    endif
  endif

  l_wind  = ( l_ndndn  .and.  l_nfnfn )

  !----------------------------------------
  ! horizontal wind components (e.g. lidar)
  !----------------------------------------
  l_mu = .false.
  l_mv = .false.
  mu   = rvind
  mv   = rvind
  status = nf90_inq_varid (ncid, 'MU', varid_MU)
  if (status == NF90_NOERR) then
    status = nf90_get_var (ncid, varid_MU, rfield)
    l_mu = (status == NF90_NOERR)

!   check for missing values
    if (all (rfield == rmissing)) then
       l_mu = .false.
    else
       mu = rfield
       where (rfield == rmissing) mu = rvind
    endif
  endif
  status = nf90_inq_varid (ncid, 'MV', varid_MV)
  if (status == NF90_NOERR) then
    status = nf90_get_var (ncid, varid_MV, rfield)
    l_mv = (status == NF90_NOERR)

!   check for missing values
    if (all (rfield == rmissing)) then
       l_mv = .false.
    else
       mv = rfield
       where (rfield == rmissing) mv = rvind
    endif
  endif
  l_uv = (l_mu .and. l_mv)

  !----------------------------------------------
  ! Spezielle Felder fuer Windprofiler Monitoring
  ! 1. Standard Deviation WIND SPEED [ M/S ]
  !----------------------------------------------

  nstdff  = rvind
  status = nf90_inq_varid (ncid, 'NSTDFF' ,  varid_NSTDFF )
  l_nstdff = .FALSE.
  if (status == nf90_noerr) then
    l_nstdff = .TRUE.
    status = nf90_get_var (ncid, varid_NSTDFF , rfield)

!   check for missing values
    if ( all( rfield == rmissing ))  then
       l_nstdff = .FALSE.
    else
       nstdff = rfield
       where ( rfield == rmissing ) nstdff = rvind
    endif
  endif

  !--------------------------------
  ! 2. Signal to Noise Ratio [ DB ]
  !--------------------------------
  nsinor   = rvind
  status = nf90_inq_varid (ncid, 'NSINOR'  ,  varid_NSINOR )
  l_nsinor = .FALSE.
  if (status == nf90_noerr) then
    l_nsinor = .TRUE.
    status = nf90_get_var (ncid, varid_NSINOR , ifield)

!   check for missing values
    if ( all( ifield == imissing ))  then
       l_nsinor = .FALSE.
    else
       nsinor = ifield
       where ( ifield == imissing ) nsinor = rvind
    endif
  endif

  !----------------------------------
  ! list of defined variables in file
  !----------------------------------
  if (netcdf_verb > 0) then
! meteorological part
    write (6,'(a,i3,1x,a,a,/,2x,12l8)') 'pe=',dace% pe ,                     &
     'l_mpn l_height l_wind   l_uv  l_mtdbt l_mtdnh l_mevss l_nltpd l_mladh',&
     'l_mlodh l_nstdff l_nsinor',                                            &
      l_mpn,l_height,l_wind,  l_uv ,l_mtdbt,l_mtdnh,l_mevss,l_nltpd,l_mladh, &
      l_mlodh,l_nstdff,l_nsinor
! technological / cloud part
    write (6,'(a,i3,1x,a,a,/,2x,17l8)') 'pe=',dace% pe ,                      &
     'l_nrara l_nsr  l_nsasa  l_na4  l_mtisi l_mhosnn l_mhobnn l_mh  ',       &
     'l_mseqm l_mvtsu l_mnh   l_nh   l_mcc   l_mcc0  l_mcc1 l_mvtsu0 l_mtn00',&
      l_nrara,l_nsr, l_nsasa, l_na4, l_mtisi,l_mhosnn,l_mhobnn,l_mh,          &
      l_mseqm,l_mvtsu,l_mnh,  l_nh,  l_mcc,  l_mcc0, l_mcc1,l_mvtsu0,l_mtn00
  end if

  !-------------------------------
  ! preset total number of reports
  !-------------------------------
  entry   = sum (source(1:ifile-1)% entries)

  !------------------
  ! loop over reports
  !------------------
! original: is number of subsets
! netcdf  : is number of report
  do is = 1, len_report

  entry1  = entry    + 1
  entry   = entry    + 1

  !---------------
  ! initialize use
  !---------------
  use = use_0

  !--------------------
  ! define head section
  !--------------------
  !=======================================
  ! derive observation type specifications
  !=======================================
  report_subt  = s2ikz(is)
  bufr_type    = s1cat(is)
  bufr_subtype = s1catls(is)
  centre       = s1cent(is)

  if (report_subt >= 0) then
    !---------------------------------
    ! derive information from DWD dbkz
    !---------------------------------
    obsid      = obstype_dbkz (report_subt)
  else
    !---------------------------------------------------------------
    ! or from bufr_type, bufr_subtype specified by generating center
    !---------------------------------------------------------------
    obsid      = obstype_bufr (bufr_type, bufr_subtype, centre)
    !------------------------
    ! optionally set DWD dbkz
    !------------------------
    if (derive_dbkz) report_subt   = obsid% dbkz
  endif
  !-------------------------------------------------------
  ! set CMA obstype, BUFR type, BUFR subtype, if not given
  !-------------------------------------------------------
  obstype                          = obsid% obstype
  if (bufr_type   <0) bufr_type    = obsid% bufrtype
  if (bufr_subtype<1 .and. obsid% centre == WMO0_ECMWF) &
                      bufr_subtype = obsid% subtype
  if (obstype < 0) cycle
  report_subti  = idb_dbk (report_subt, obstype)

  head% obstype     = obstype
  head% dbkz        = report_subt
! head% modtype     = rept_char(obstype)% mod
  head% modtype     = TEMP
  head% buf_type    = bufr_type
  head% buf_subtype = bufr_subtype
  head% codetype    = obsid% codetype
  head% time        = stime(is)
  head% db_time     = db_time(is)
  head% idbk        = report_subti
  head% source      = ifile
  head% record      = is
  head% id          = entry1
  head% center      = s1cent(is)
  head% subcenter   = s1cents(is)

  !------------------------------------------------------
  ! set different code types for different wind profilers
  !------------------------------------------------------
  if (head% dbkz == 10553 .and. head% codetype == 132) then
    select case (istidn(is)/1000)
    case (01:17)
      head% codetype = 132 ! European  wind profiler
    case (45:47)
      head% codetype = 134 ! Japanese  wind profiler
    case (48)
      head% codetype = 132 ! Singapore wind profiler
    case (61)
      head% codetype = 132 ! French wind profiler in La Reunion/Madagaskar
    case (70:74)
      head% codetype = 136 !  wind/profiler/rass report (USA)
    case (81)
      head% codetype = 132 ! French wind profiler in French Guiana
    case (89)
      head% codetype = 132 ! Australian wind profiler (Davis Station/Antarctic)
    case (91:95)
      head% codetype = 132 !  wind/profiler/rass report (Australia)
    case default
                     ! 133 ! European sodar/rass report
                     ! 137 ! radar VAD wind profile report
      write (0,*) &
        'read_temp_netcdf: skipping unknown wind profiler codetype=',&
        head% codetype, 'dbkz=', head% dbkz, 'statid =',istidn(is)
      cycle
    end select
  endif

  !------------------------------------------------------
  ! set different code types for different vad wind profiles
  !------------------------------------------------------
  if (head% dbkz == 10600 .and. head% codetype == 137) then
    select case (istidn(is)/1000)
    case (01:14)
      head% codetype = 137 !  vad wind profile report (Europe)
    case (60:62)
      head% codetype = 137 !  vad wind profile report (Northern Africa)
    case (70:74)
      head% codetype = 137 !  vad wind profile report (Canada)
    case (91:95)
      head% codetype = 137 !  vad wind profile report (Australia)
    case default
                     ! 133 ! European sodar/rass report
      write (0,*) &
        'read_temp_netcdf: skipping unknown wind profiler codetype=',&
        head% codetype, 'dbkz=', head% dbkz, 'statid =',istidn(is)
      cycle
    end select
  endif


  if ( lpr_temp .and. is < npr_extd ) then
    write (6,'()')
    write (6,'( 8(a16, i6  ,/),   &
              & 2(a16, a   ,/),   &
              & 6(a16, i6  ,/) )' )                      &
      'pe='         ,dace% pe,                           &
      'head is='    ,is,                                 &
      'obstype='    , head% obstype  ,                   &
      'dbkz='       , head% dbkz     ,                   &
      'modtype='    , head% modtype  ,                   &
      'buf_type='   , head% buf_type ,                   &
      'buf_subtype=', head% buf_subtype,                 &
      'codetype='   , head% codetype   ,                 &
      'time='       , cyyyymmddhhmmss (head% time)   ,   &
      'db_time='    , cyyyymmddhhmmss (head% db_time),   &
      'dbk='        , head% idbk     ,                   &
      'source='     , head% source   ,                   &
      'record='     , head% record   ,                   &
      'id='         , head% id       ,                   &
      'center='     , head% center    ,                  &
      'subcenter='  , head% subcenter
  endif

  !------------------------------------------------
  ! exclude windprofiler without measurement type
  !------------------------------------------------
  if (na4(is) == -1 .and. (head% dbkz == 10553 .or. head% dbkz == 10600)) cycle

  !--------------------------------------------
  ! perform simple generic check on report type
  !--------------------------------------------
  call check_report_0 (use, head, 1)
  if (use% state <= STAT_DISMISS) cycle

  !------------------
  ! create new report
  !------------------
  spt0             = empty
  spt0% use        = use
  spt0% hd         = head
  spti             = spt0
! spti% hd% subset = is

  !-----------------------
  ! check number of levels
  !-----------------------
  nlev = 0
  select case (dim1_length)
  case ('MEDRE')
    nlev = medre(is)
  case ('MDREP')
    nlev = mdrep(is)
  case default
    call finish ('read_temp_netcdf','invalid dim1_length: '//dim1_length)
  end select

  if(nlev < 1) then
    call decr_rpt_use (spti, CHK_INSDAT, comment='read_temp_netcdf: no. levels < 1')
    if ( lpr_temp ) then
       print * , 'pe=',dace% pe,' read_temp_netcdf  after decr_rpt_use CHK_INSDAT NO levels is=',is
    endif
    cycle
  endif

  !----------------------------------------
  ! allocate memory
  ! construct 'empty' Report data structure
  !----------------------------------------
  allocate (btemp (nlev))
  call construct_temp (btemp)

  !------------------------------
  ! process the following entries
  !------------------------------
  spti% corme        = max ( s1updat(is), 0)
  spti% col% c% dlat = mlah(is)
  spti% col% c% dlon = mloh(is)
  spti% z            = mhobnn(is)
  spti% actual_time  = obs_time(is)
  spti% statid       = ystidn(is)
  spti% ident        = istidn(is)
  spti% col% nlev    = nlev
  spti% sttyp        = nrara(is)                ! radiosonde type/system
  spti% stret        = nsr  (is)                ! solar and infrared radiation correction
  spti% meas_type    = na4  (is)                ! type of measuring equipment used
  spti% tracking     = nsasa(is)                ! tracking technique,status of system
  !--------------------------------------------------------------------
  ! Map station id to integer "hash value" for vectorizable comparisons
  ! Be careful about reports with invalid/non-standard WIGOS station ID
  ! that may report redundantly using traditional station identifiers.
  !--------------------------------------------------------------------
  if (lwsi(is)) then
     spti% wsi       = wsi    (is)
  end if
  if (lwsi(is) .and. spti% statid(1:1) == "_") then
     spti% stat_hash = wsihash(is)
  else
     spti% stat_hash = transfer (spti% statid, spti% stat_hash)
  end if

  if ( lpr_temp  .and. is < npr_extd ) then
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
          'spti% col% nlev    = ', spti% col% nlev    ,  &
          'spti% sttyp        = ', spti% sttyp        ,  &
          'spti% stret        = ', spti% stret        ,  &
          'spti% meas_type    = ', spti% meas_type    ,  &
          'spti% tracking     = ', spti% tracking
  endif

! default: -999._wp
! height of barometer above mean sea level
  if ( mhobnn(is) /= -999.) then
    spti%  z       = mhobnn(is)
  else if ( mhosnn(is) /= -999. ) then
! alternative: height of station ground above mean sea  [m]
    spti%  z          = mhosnn(is)
! alternative: height of station [m]
  else if ( mhp(is) /= -999 ) then
    spti%  z          = mhp(is)
  endif

! mtisi    ! time significance
! mh       ! height [ m ]
! mseqm    ! station elevation quality mark(mobil)
! mvtsu    ! vertical significance (surface observ.)
! mnh      ! cloud amount
! nh       ! height of base of cloud
! mcc      ! code_table for cloud type
! mcc0     ! code_table for cloud type
! mcc1     ! code_table for cloud type
! mvtsu0   ! vertical significance (surface observ.)
! mtn00    ! sea/water temperature

  !------------------
  ! loop over levels
  !------------------
! do ilev = 1,medre(is)
  do ilev = 1,spti% col% nlev
    ! long time period or displacement   [s]
    btemp(ilev)% secs = merge (nltpd(ilev,is), 0, nltpd(ilev,is) /= inv_secs)

    !  pressure (vert.location)          [Pa]
    call set_datum (btemp(ilev)% p  ,mpn(ilev,is)   ,spti% corme)

    !  nhhhn geopotential height [ gpm ]; pilot nhhh
    call set_datum (btemp(ilev)% gp ,nhhhn(ilev,is) ,spti% corme)

    ! temperature/dry bulb temperature   [K]
    call set_datum (btemp(ilev)% t  ,mtdbt(ilev,is) ,spti% corme)

    ! dew-point temperature              [K]
    call set_datum (btemp(ilev)% td ,mtdnh(ilev,is) ,spti% corme)

    ! wind direction                     [degree]
    call set_datum (btemp(ilev)% dd ,ndndn(ilev,is) ,spti% corme)

    ! wind speed                         [m/s]
    call set_datum (btemp(ilev)% ff ,nfnfn(ilev,is) ,spti% corme)

    ! horizontal wind components         [m/s]
    if (l_uv) then
      call set_datum (btemp(ilev)% uu ,mu(ilev,is)  ,spti% corme)
      call set_datum (btemp(ilev)% vv ,mv(ilev,is)  ,spti% corme)
    end if

    ! extended vertical sounding significance
    btemp(ilev)% vss = mevss (ilev,is)

    !--------------------------------------------------------------
    !   Windprofiler Quality Information written into feedback file
    !   1. Signal to Noise [DB]
    !--------------------------------------------------------------

    if (nsinor(ilev,is) /= rvind ) then
       btemp(ilev)%uu%pcc = int (nsinor(ilev,is))
       btemp(ilev)%vv%pcc = int (nsinor(ilev,is))
    endif
    !------------------------------
    !   2. Standard Deviation [m/s]
    !------------------------------
    if (nstdff(ilev,is) /= rvind ) then
       btemp(ilev)%uu%ac = nstdff(ilev,is)
       btemp(ilev)%vv%ac = nstdff(ilev,is)
    endif

    !--------------------------------------
    ! translate to feedback-file convention
    !--------------------------------------
    call ffvss_from_evss (btemp(ilev))
    !---------------------------------------------------------
    ! translate from extended vertical sounding significance
    !             to          vertical sounding significance
    ! needed for subsequent subroutines i. e. check_store_temp
    !----------------------------------------------------------
    btemp(ilev)% vss = mvss (ilev,is)

    if (l_mladh) then
       ! latitude displacement (high accuracy)
       btemp(ilev)% dlat = mlah(is) + mladh(ilev,is)
       ! longitude displacement (high accuracy)
       btemp(ilev)% dlon = mloh(is) + mlodh(ilev,is)
    else
       ! latitude  (high accuracy)
       btemp(ilev)% dlat = mlah(is)
       ! longitude (high accuracy)
       btemp(ilev)% dlon = mloh(is)
    end if

    ! loop 001 : not stored yet
    !  mdrep    ! delayed descriptor replication factor
    !  nltpd0   ! long time period or displacement
    !  mevss0   ! extended vertical sounding significance
    !  mpn0     ! pressure (vert.location)
    !  nhhh0    ! geopotential height [ gpm ]             ; pilot
    !  mladh0   ! latitude displacement (high accuracy)
    !  mlodh0   ! longitude displacement (high accuracy)
    !  nvbvb    ! absolut wind shear in 1 km layer below
    !  nvava    ! absolut wind shear in 1 km layer above

    !----------------------------------------------------------------------------
    ! WMO table 0 08 001 (7 bits)
    ! Vertical Sounding Significance                         (BUFR Bit number)
    !----------------------------------------------------------------------------
    !                     lev_surf  =  64 ! Surface                           (1)
    !                     lev_std   =  32 ! Standard level                    (2)
    !                     lev_tropo =  16 ! Tropopause level                  (3)
    !                     lev_max_v =   8 ! Maximum wind level                (4)
    !                     lev_sig_t =   4 ! Significant level, temperature/rh (5)
    !                     lev_sig_v =   2 ! Significant level, wind           (6)
    !                     lev_miss  =   1 ! Missing value                 (All 7)
    !----------------------------------------------------------------------------


    !--------------------------------------------------------------------------
    ! WMO table 0 08 042 (18 bits)
    ! Extended vertical sounding significance
    !
    ! bit no       meaning                 (bit beginning on the left hand side
    !   1  2**17   surface                 (  1   2  3  4 ..             18
    !   2     16   Standard level          ( 18  17 16 15 ..              1
    !   3     15   Tropopause level
    !   4     14   Maximum wind level
    !   5     13   Significant temperature level
    !   6     12   Significant humidity level
    !   7     11   Significant wind level
    !   8     10   Beginning of missing temperature data
    !   9      9   End of missing temperature data
    !  10      8   Beginning of missing humidity data
    !  11      7   End of missing humidity data
    !  12      6   Beginning of missing wind data
    !  13      5   End of missing wind data
    !  14      4   Top of wind sounding
    !  15      3   Level determined by regional decision
    !  16      2   Reserved
    !  17      1   Pressure level originally indicated by height as the
    !              vertical coordinate
    !  18 all  0   Missing value
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    ! vertical sounding significance for FEEDBackfile
    !
    ! bit no       meaning                (bit beginning on the RIGHT hand side
    !                                      Fortran convention
    !   0          Surface                (  8   7   6   5  4   3   2   1   0
    !   1          Standard level
    !   2          Tropopause level
    !   3          Maximum wind level
    !   4          Significant temperature  or wind level
    !   6          Missing significance (in BUFR)
    !
    ! significance is defined per level and variable
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    ! WMO table 0 02 003
    ! Type of measuring equipment used  (na4)
    !
    ! Code figure  meaning
    !    0         Pressure Instrument associated with wind measuring equipment
    !    1         Optical theodolite
    !    2         Radio theodolite
    !    3         Radar
    !    4         VLF Omega
    !    5         Loran C
    !    6         Wind profiler
    !    7         Satellite navigation
    !    8         Radio-Acoustic Sounding System (RASS)
    !    9         Sodar
    ! 10-13        Reserved
    !   14         Pressure instrument associated with wind measuring equipment
    !              but pressure element failed during ascent
    !   15         Missing value
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    ! WMO COMMON CODE TABLE C-7: Tracking technique/status of system used
    !  Code Table 3872 - sasa for alphanumeric code (nsasa)
    !  Code Table 0 02 014 in BUFR
    !
    ! Code figure  meaning
    !    0         No windfinding
    !    1         Automatic with auxiliary optical direction finding
    !    2         Automatic with auxiliary radio direction finding
    !    3         Automatic with auxiliary ranging
    !    4         Not used
    !    5         Automatic with multiple VLF-Omega signals
    !    6         Automatic cross chain Loran-C
    !    7         Automatic with auxiliary wind profiler
    !    8         Automatic satellite navigation
    ! 9-18         Reserved
    !   ...
    !  127         Missing value
    !--------------------------------------------------------------------------

    !------------------------
    ! end of loop over levels
    !------------------------
  end do

  ! !-----------------------------------------------
  ! ! if station name is not set, use station number
  ! !-----------------------------------------------
  ! if (spti% statid =='') write(spti% statid,'(i5.5)') spt% ident

  where (iand (btemp% vss, lev_miss) /=0 ) btemp% vss = lev_miss

  !----------------
  ! standard checks
  !----------------
  call check_report_1 (spti)
  lk = spti% use% state > STAT_DISMISS

  if (lk) then
    !--------------------------------------------
    ! if not present, derive pressure from height
    ! (US Standard Atmosphere)
    !--------------------------------------------
    if (rept_use(spti%hd% obstype)% deriv_p /= 0) then
      do i=1,size(btemp)
        if (btemp(i)% p % qc /= QC_OK .and. btemp(i)% gp% qc == QC_OK .and. &
            btemp(i)% gp% o  >  0.                                  ) then
            btemp(i)% p % o   = p_h_usstd (real(btemp(i)% gp% o,wp))
            btemp(i)% p % qc  = QC_OK
            btemp(i)% p % src = SRC_DER
            btemp(i)% gp% qc  = QC_NOUSE
        endif
      end do
    endif

    !----------------------------------------------------------------
    ! For profilers not reporting pressure but height, use the latter
    ! if derivation of pressure from height is not enabled.
    !----------------------------------------------------------------
    if (     spti%hd% obstype            == OT_PILOT  .and. &
        any (spti%hd% codetype           == ct_wprof) .and. &
             rept_use(OT_PILOT)% deriv_p == 0               ) then
      do i = 1, size (btemp)
        if (btemp(i)% gp% qc == QC_OK     .and. &
            btemp(i)% p % qc /= QC_OK     .and. &
            btemp(i)% gp% o  >   -999._sp .and. &
            btemp(i)% gp% o  <  99999._sp       ) then
           btemp(i)% height = btemp(i)% gp% o
        end if
      end do
    end if

    !++++++++++++++++++++++++++++++++++++
    ! inhibit correction message handling
    !++++++++++++++++++++++++++++++++++++
    if (corme_conv == -1) spti% corme = 0

    !-------------------------
    ! check for double entries
    !-------------------------

    !------------------------------------------------------
    ! exclude "new" BUFR data from standard duplicate check
    !------------------------------------------------------
    select case (spti% hd% dbkz)
    case (10553,10600,10520,10521,10526,10527,10776,10777,10780,10782,10783, &
          10574,10785,10516,10517,10570,584)
       l_newbufr = .true.
    case default
       l_newbufr = .false.
    end select

    if (.not. l_newbufr) then

#if defined (__SX__)
      !------------------------------------------------------
      ! Vectorized preselection of observations to compare to
      !------------------------------------------------------
      call vector_indices     ! sets n_idx, 'idx'
      !------------------------------------------
      ! Now we walk only through the list of hits
      !------------------------------------------
      if (corme_conv >= 2) then
        call correct_report (obs, spti, btemp, lk, n_idx, idx=idx)
      else
        call correct_merge  (obs, spti, btemp, lk, n_idx, idx=idx)
      endif
#else
      if (corme_conv >= 2) then
        call correct_report (obs, spti, btemp, lk, obs% n_spot)
      else
        call correct_merge  (obs, spti, btemp, lk, obs% n_spot)
      endif
#endif
    endif
  endif

  if (lk) then
!   print *
!   print *, "Before check_store_temp:"
!   call print_temp_data (btemp)
!   print *, "------------------------"
    call check_store_temp (btemp, spti, obs, novalid)
    nkeep = nkeep + 1
  endif

  if (prt_data) then
    print *, 'pe=',dace% pe,' mo_temp: station-id lat lon typ subtyp'
    print '(1x,a,2f12.5,2i8)', spti% statid,     &
         spti% col% c% dlat, spti% col% c% dlon, &
         spti% hd% buf_type, spti% hd% buf_subtype
    print *, 'pe=',dace% pe,' mo_temp: station-time actual time  db-time'
    call print (spti% hd% time)
    call print (spti%     actual_time)
    call print (spti% hd% db_time)
    print *, 'pe=',dace% pe,' mo_temp: dbkz =',spti% hd% dbkz, ' nkeep =', nkeep
    call print_temp_data (btemp)
  endif

  deallocate (btemp)
  !-------------------------
  ! end of loop over reports
  !-------------------------
  end do
! if any temp/pilot stored
  lkeep = nkeep > 0

  !-----------
  ! deallocate
  !-----------

   if (allocated  (nrara   )) deallocate (nrara   )
   if (allocated  (nsr     )) deallocate (nsr     )
   if (allocated  (nsasa   )) deallocate (nsasa   )
   if (allocated  (na4     )) deallocate (na4     )
   if (allocated  (mtisi   )) deallocate (mtisi   )
   if (allocated  (mhosnn  )) deallocate (mhosnn  )
   if (allocated  (mhobnn  )) deallocate (mhobnn  )
   if (allocated  (mh      )) deallocate (mh      )
   if (allocated  (mseqm   )) deallocate (mseqm   )
   if (allocated  (mvtsu   )) deallocate (mvtsu   )
   if (allocated  (mnh     )) deallocate (mnh     )
   if (allocated  (nh      )) deallocate (nh      )
   if (allocated  (mcc     )) deallocate (mcc     )
   if (allocated  (mcc0    )) deallocate (mcc0    )
   if (allocated  (mcc1    )) deallocate (mcc1    )
   if (allocated  (mvtsu0  )) deallocate (mvtsu0  )
   if (allocated  (mtn00   )) deallocate (mtn00   )

! loop 000
   if (allocated  (medre   )) deallocate (medre   )

   if (allocated  (nltpd   )) deallocate (nltpd   )
   if (allocated  (mevss   )) deallocate (mevss   )
   if (allocated  (mvss    )) deallocate (mvss    )
   if (allocated  (mpn     )) deallocate (mpn     )
   if (allocated  (nhhhn   )) deallocate (nhhhn   )
   if (allocated  (mladh   )) deallocate (mladh   )
   if (allocated  (mlodh   )) deallocate (mlodh   )
   if (allocated  (mtdbt   )) deallocate (mtdbt   )
   if (allocated  (mtdnh   )) deallocate (mtdnh   )
   if (allocated  (ndndn   )) deallocate (ndndn   )
   if (allocated  (nfnfn   )) deallocate (nfnfn   )
   if (allocated  (mu      )) deallocate (mu      )
   if (allocated  (mv      )) deallocate (mv      )

! loop 001
!  if (allocated  (mdrep   )) deallocate (mdrep   )

!  if (allocated  (nltpd0  )) deallocate (nltpd0  )
!  if (allocated  (mevss0  )) deallocate (mevss0  )
!  if (allocated  (mpn0    )) deallocate (mpn0    )
!  if (allocated  (nhhh0   )) deallocate (nhhh0   )
!  if (allocated  (mladh0  )) deallocate (mladh0  )
!  if (allocated  (mlodh0  )) deallocate (mlodh0  )
!  if (allocated  (nvbvb   )) deallocate (nvbvb   )
!  if (allocated  (nvava   )) deallocate (nvava   )


  if (allocated  (ifield  )) deallocate  (ifield  )
  if (allocated  (rfield  )) deallocate  (rfield  )
! if (allocated  (ifield0 )) deallocate  (ifield0 )
! if (allocated  (rfield0 )) deallocate  (rfield0 )
  if (allocated  (ifield1 )) deallocate  (ifield1 )
  if (allocated  (rfield1 )) deallocate  (rfield1 )

  contains
!------------------------------------------------------------------------------
    subroutine vector_indices
    !----------------------------------------------------------
    ! Calculate index array
    ! for vectorized preselection of observations to compare to
    !----------------------------------------------------------
      integer :: i, k, n
      !--------------------------------------------------------------------
      ! Map station id to integer "hash value" for vectorizable comparisons
      !--------------------------------------------------------------------
!     spti% stat_hash = transfer (spti% statid, spti% stat_hash)
      !------------------------------------------------------
      ! Vectorized preselection of observations to compare to
      !------------------------------------------------------
      n = obs% n_spot
      k = 0; if (allocated (mask)) k = size (mask)
      if (n > k) then
        if (allocated (mask)) deallocate (mask)
        allocate (mask(max (n, 2*k)))
      end if
      n_idx = 0
      do i = 1, n
        mask(i) = obs% spot(i)% stat_hash   == spti% stat_hash   .and. &
                  obs% spot(i)% hd% obstype == spti% hd% obstype
        if (mask(i)) n_idx = n_idx + 1
      end do
      if (n_idx > 0) then
        k = 0; if (allocated (idx)) k = size (idx)
        if (n_idx > k) then
          if (allocated (idx)) deallocate (idx)
          allocate (idx(n_idx))
        end if
        idx(1:n_idx) = pack ( (/ (k, k=1,n) /), mask(1:n))
      end if
    end subroutine vector_indices
!------------------------------------------------------------------------------
  end subroutine read_temp_netcdf
!==============================================================================
  subroutine correct_report (obs, spti, btemp, lkeep, n, idx)
  type (t_obs)      ,intent(inout) :: obs      ! observation container
  type (t_spot)     ,intent(inout) :: spti     ! new observation to check
  type (t_temp)     ,pointer       :: btemp(:) ! TEMP/PILOT specific info
  logical           ,intent(inout) :: lkeep    ! accept observation ?
  integer           ,intent(in)    :: n        ! number of previous obsv.
  integer ,optional ,intent(in)    :: idx(:)   ! preselected observations
  !---------------------------------------------------------------
  ! Check for double entries.
  ! The current report ('spti', 'btemp') is checked versus the
  ! previous observations (contained in 'obs') for double entries.
  ! for better performance on vector processors relevant entries
  ! may be preselected by the index array 'idx'.
  ! This routine only handles correction messages and double reports for
  ! the same dbkz. Merges of parts A,B,C,D are handled in subroutine
  ! merge_report.
  !---------------------------------------------------------------

    integer               :: i, k            ! indices
    character(len=32)     :: com             ! comment line
    logical               :: diff,m12,m21,z2 ! flags from cmp_temp
    type(t_temp) ,pointer :: t1 (:)          ! TEMP/PILOT specific info
    integer               :: action          ! what to do
    integer               :: chk_            ! check (reason for rejection)
    integer               :: nonew           ! number of observations stored
    integer               :: noold           ! number of observations stored
    logical               :: same_time       ! same reference or launch time?

    !------------------------------------
    ! loop over previously stored reports
    !------------------------------------
    do k = 1, n
      i = k; if(present (idx)) i = idx(k)
      !----------------------------------------------------
      ! check for corresponding station name, subtype, time
      !----------------------------------------------------
      if   (obs% spot(i)%        statid  == spti%        statid  .and. &
            obs% spot(i)%    hd% obstype == spti%    hd% obstype .and. &
            obs% spot(i)%    hd% dbkz    == spti%    hd% dbkz    .and. &
!           obs% spot(i)%    hd% time    == spti%    hd% time    .and. &
            obs% spot(i)% col%c% dlon    == spti% col%c% dlon    .and. &
            obs% spot(i)% col%c% dlat    == spti% col%c% dlat          ) then

       same_time = obs% spot(i)% hd% time    == spti% hd% time    ! Reference
       if (same_time .and. corme_mode >= 3)                     & ! +
       same_time = obs% spot(i)% actual_time == spti% actual_time ! Launch time

       if (same_time) then
        call load_temp (obs, obs% spot(i), t1)
        call cmp_temp (t1, btemp, diff, m12, m21, z2)
!write(*,*) "after cmp_temp: diff, m12, m21, z2 =", diff, m12, m21, z2
        !-----------------------------------------------------------------
        ! decide what to do, set flag 'action' :
        ! 1=keep old report, 2=keep new, 3=keep both, 4=merge,-=suspicious
        !-----------------------------------------------------------------
        action = 0
        chk_   = CHK_NONE
        write (com,'(a,2i8)') 'correct_report:', obs% spot(i)% hd% id, &
                                                 spti%         hd% id
        !---------------------------
        ! Handle correction messages
        !---------------------------
        if (spti% corme > 0 .or. obs% spot(i)% corme > 0) then
         select case (corme_mode)
         case (0)
          !------------
          ! Old version
          !------------
          if     (spti% corme - obs% spot(i)% corme == 1) then ! 2 corrects 1
            action = 4                                         ! merge
            chk_   = CHK_CORR
          elseif (obs% spot(i)% corme - spti% corme == 1) then ! 1 corrects 2
            action = 4                                         ! merge
            chk_   = CHK_CORR
          else
            action = -5                                        ! dismiss both
            chk_   = CHK_CORRERR
          endif
         case default                                          ! corme_mode /= 0
          !-------------------------------------------------------------
          ! Allow for gaps in sequence numbers; merge only if
          ! "correction message" is more than 20% shorter but consistent
          !-------------------------------------------------------------
          if     (spti% corme > obs% spot(i)% corme)     then  ! 2 corrects 1
            if (diff .or. (size(btemp) >= size(t1)*4/5)) then
              action = 2                                       ! replace (1 by 2)
              chk_   = CHK_CORR
            else
              action = 4                                       ! merge
              chk_   = CHK_CORR
            end if
          elseif (spti% corme < obs% spot(i)% corme)     then  ! 1 corrects 2
            if (diff .or. (size(t1) >= size(btemp)*4/5)) then
              action = 1                                       ! dismiss new
              chk_   = CHK_CORR
            else
              action = 4                                       ! merge
              chk_   = CHK_CORR
            end if
          else
            action = -5                                        ! dismiss both
            chk_   = CHK_CORRERR
          endif
         end select
        endif

        if (action == 0 .and. .not. diff) then                 ! no substantial difference
          if (m21 .and. .not. m12) then
            action = -2                                        ! keep new
            chk_   = CHK_REDUNDANT
          else if (m12 .and. .not. m21)   then
            action = -1                                        ! keep old
            chk_   = CHK_REDUNDANT
          else if (.not.m12 .and. .not. m21) then              ! no difference
            action = 1                                         ! keep old
            chk_   = CHK_REDUNDANT
          else
            action = -4                                        ! merge
            chk_   = CHK_MERGE
!           action = -5                                        ! dismiss both
!           chk_   = CHK_DBLERR
!           action = -1                                        ! keep old
!           chk_   = CHK_DBLERR
          endif
        endif
        if (action == 0) then                                  ! substantial difference
          action = -5                                          ! dismiss both
          chk_   = CHK_DBLERR
        endif

        if (obs% spot(i)% actual_time /= spti% actual_time) &  ! different actual time
          action = - abs (action)

        !---------
        ! printout
        !---------
        write(6,&
'(a,i4,1x,a,i5," - ",i3,i8,i3,i8," - ",2(a,1x,a,l2," - "),a,1x,a," - ",4l2,2i3,l2,i3,1x,a)')&
          '#CORME',spti% hd% obstype,spti% statid, spti% hd% dbkz,                      &
          spti% hd% source, spti% hd% record,                                           &
          obs% spot(i)% hd% source, obs% spot(i)% hd% record,                           &
          chhmmss(spti% actual_time),chhmmss(obs% spot(i)% actual_time),                &
          spti% actual_time==obs% spot(i)% actual_time,                                 &
          chhmmss(spti% hd%    time),chhmmss(obs% spot(i)% hd%    time),                &
          spti% hd%    time==obs% spot(i)% hd%    time,                                 &
          chhmmss(spti% hd% db_time),chhmmss(obs% spot(i)% hd% db_time),                &
          diff, m12, m21, z2,obs% spot(i)% corme, spti% corme, spti% corme > 0,         &
          action, chk(chk_)% mnem
        nonew = -999
        noold = obs% spot(i)%o%n

        !-----------------------------------------------------------------
        ! actually correct reports according to flag 'action'
        ! 1=keep old report, 2=keep new, 3=keep both, 4=merge,-=suspicious
        !-----------------------------------------------------------------
        select case (abs(action))
        case (1)
          !---------------------------------
          ! keep old report: dismiss new one
          !---------------------------------
          call decr_rpt_use (spti         ,chk_ ,comment=com)
        case (2)
          !-------------------------------------------------
          ! keep new report: dismiss old, replace by new one
          !-------------------------------------------------
          call decr_rpt_use (obs% spot(i) ,chk_ ,comment=com)
          call check_store_temp (btemp, spti, obs, nonew, repl = i)
        case (3)
          !---------------------------------------------
          ! keep both reports: currently not implemented
          !---------------------------------------------
          call finish('correct_report','dont keep both')
        case (4)
          !----------------------------------------------
          ! merge reports: dismiss old, replace by merged
          !----------------------------------------------
          call merge_temp (t1, btemp)
          spti% corme = max (spti% corme, obs% spot(i)% corme)
          call decr_rpt_use (obs% spot(i) ,chk_ ,comment=com)
          call check_store_temp (t1, spti, obs, nonew, repl = i)
          if (prt_data) then
             print '(a,a,a)', "### Corrected TEMP ", spti% statid, ":"
             call print_temp_data (t1)
          end if
        case (5)
          !---------------------
          ! dismiss both reports
          !---------------------
          call decr_rpt_use (spti         ,chk_ ,comment=com)
          call decr_rpt_use (obs% spot(i) ,chk_ ,comment=com)
          select case (spti% hd% dbkz)
          !----------------------------------------
          ! print message for DOUBLE TEMP SHIP BUFR
          !----------------------------------------
          case (10520, 10521, 10526, 10527, 10776, 10777, 10782, 10783, &
                10574, 10785, 10516, 10517, 10570)
          write (6,*)'### DOUBLE TEMP BUFR rejected: ',                   &
            spti% statid, spti% hd% dbkz, cyyyymmddhhmm  (spti% hd% time),&
            ' ',  cyyyymmddhhmmss(spti% actual_time), size(t1), size(btemp)
          end select
        case default
          call finish ('correct_report','invalid action')
        end select

        if (nonew /= -999) write(6,'(a,i4,1x,a,i5," - ",i3,i8,a,2i6)') &
          '#CORME',spti% hd% obstype,spti% statid,spti% hd% dbkz,      &
          spti% hd% source, spti% hd% record,                          &
          ' old,new: ',noold, nonew
        lkeep = .false.
        deallocate (t1)
        exit

       end if ! same_time
      endif
    end do
  end subroutine correct_report
!------------------------------------------------------------------------------
  subroutine correct_merge  (obs, spti, btemp, lkeep, n, idx)
  type (t_obs)      ,intent(inout) :: obs      ! observation container
  type (t_spot)     ,intent(inout) :: spti     ! new observation to check
  type (t_temp)     ,intent(in)    :: btemp(:) ! TEMP/PILOT specific info
  logical           ,intent(inout) :: lkeep    ! accept observation ?
  integer           ,intent(in)    :: n        ! number of previous obsv.
  integer ,optional ,intent(in)    :: idx(:)   ! preselected observations
  !------------------------------------------------------------------
  ! Check for double entries.
  ! Old Version: Correction messages and sections A,B,C,D are handled
  ! in the same loop.
  ! The current report ('spti', 'btemp') is checked versus the
  ! previous observations (contained in 'obs') for double entries.
  ! for petter performance on vector processors relevant entries
  ! may be preselected by the index array 'idx'.
  ! New version: cf. subroutines correct_report, merge_report
  !------------------------------------------------------------------

    integer               :: i, k            ! indices
    character(len=32)     :: com             ! comment line
    logical               :: diff,m12,m21,z2 ! flags from cmp_temp
    type(t_temp) ,pointer :: t1 (:)          ! TEMP/PILOT specific info
    integer               :: novalid         ! number of valid observations

    !------------------------------------
    ! loop over previously stored reports
    !------------------------------------
    do k = 1, n
      i = k; if(present (idx)) i = idx(k)
      !----------------------------------------------------
      ! check for corresponding station name, subtype, time
      !----------------------------------------------------
      if   (obs% spot(i)%     statid  == spti% statid      .and. &
            obs% spot(i)% hd% obstype == spti% hd% obstype .and. &
      ((any(obs% spot(i)% hd% dbkz    == kz_pilot)         .and. &
                 any(spti% hd% dbkz   == kz_pilot))        .or.  &
       (any(obs% spot(i)% hd% dbkz    == kz_temp )         .and. &
                 any(spti% hd% dbkz   == kz_temp ))        .or.  &
       (any(obs% spot(i)% hd% dbkz    == kz_tship)         .and. &
                 any(spti% hd% dbkz   == kz_tship))        .or.  &
       (any(obs% spot(i)% hd% dbkz    == kz_tmobil)        .and. &
                 any(spti% hd% dbkz   == kz_tmobil))       .or.  &
       (any(obs% spot(i)% hd% dbkz    == kz_tdrop)         .and. &
                 any(spti% hd% dbkz   == kz_tdrop))        .or.  &
       (any(obs% spot(i)% hd% dbkz    == kz_pship)         .and. &
                 any(spti% hd% dbkz   == kz_pship)))       .and. &
            obs% spot(i)% actual_time == spti% actual_time) then

        call load_temp (obs, obs% spot(i), t1)
        call cmp_temp (t1, btemp, diff, m12, m21, z2)
        !----------------
        ! correct reports
        !----------------
        if (spti% corme > 0) then
          write (com,'("corme,dbkz:",i2,i5,",new:",i2,i5)')&
            obs% spot(i)% corme, obs% spot(i)% hd% dbkz,   &
            spti% corme, spti% hd% dbkz
          call decr_rpt_use (obs% spot(i), CHK_CORR, comment=com)
          if (diff .or. m21) THEN
            call merge_temp (t1, btemp)
            call check_store_temp (t1, spti, obs, novalid, repl = i)
          endif
        else
          !--------------
          ! merge reports
          !--------------
          if (diff .or. m21) then
          !----------------------------------------------
          ! the new report extends the old one or differs
          !----------------------------------------------
            call merge_temp (t1, btemp)
            write (com,'("corme,dbkz:",i2,i5,",new:",i2,i5)')&
              obs% spot(i)% corme, obs% spot(i)% hd% dbkz,   &
              spti% corme, spti% hd% dbkz
            call decr_rpt_use (obs% spot(i), CHK_MERGE, comment=com)
            call check_store_temp (t1, spti, obs, novalid, repl = i)
          else
            !-------------------------------------------
            ! the new report does not extend the old one
            !-------------------------------------------
            write (com,'("corme,dbkz:",i2,i5,",new:",i2,i5)')&
              obs% spot(i)% corme, obs% spot(i)% hd% dbkz,   &
              spti% corme, spti% hd% dbkz
            if (z2) then
              call decr_rpt_use (spti, CHK_INSDAT, comment=com)
            else
              call decr_rpt_use (spti, CHK_REDUNDANT, comment=com)
            endif
          endif
        endif
        lkeep = .false.
        deallocate (t1)
        exit
      end if
    end do
  end subroutine correct_merge
!------------------------------------------------------------------------------
  subroutine merge_report  (obs)
  type (t_obs)      ,intent(inout) :: obs      ! observation container
  !-------------------------------------------------------------------
  ! Check for double entries.
  ! New Version: handling of sections A,B,C,D.
  ! Corresponding TEMP or PILOT reports are merged.  Correction messages
  ! and double reports should be handled by subroutine correct_report.
  !-------------------------------------------------------------------
    integer                :: i1, k1, j1        ! indices
    integer                :: i2, k2            ! indices
    type (t_spot) ,pointer :: spt2, spt1        ! observation to check
    type (t_spot)          :: spt               ! observation to check
    type (t_temp) ,pointer :: t1 (:)            ! TEMP/PILOT specific info
    type (t_temp) ,pointer :: t2 (:)            ! TEMP/PILOT specific info
    integer                :: nonew             ! number of observation stored
    character(len=40)      :: com               ! comment line
    logical                :: changed
    logical                :: first
    logical                :: mask(obs% n_spot) ! Mask of eligible reports
    integer                :: nrep, n_idx       ! No. of reports to process
    integer   ,allocatable :: idx  (:)          ! Index pointer to reports
    integer   ,allocatable :: idx1 (:)          ! Index pointer to reports
    integer   ,allocatable :: idx2 (:)          ! Index pointer to reports
    integer   ,allocatable :: prio (:)          ! priority
    logical                :: same_time         ! same reference or launch time?
    logical                :: same_id           ! same station id?
    logical                :: ldown             ! downsondes?
    logical                :: lmerge            ! merge reports
    logical   ,allocatable :: down (:)          ! downsonde specific checks?
    real(wp)               :: dt                ! launch time difference

    if (corme_conv < 2) return

    !----------------------------------------
    ! Set up index vector to eligible reports
    !   TEMP or (PILOT but not wind profiler)
    !----------------------------------------
    do k2 = 1, obs% n_spot
       mask(k2) =      obs% spot(k2)% hd%  obstype == OT_TEMP  .or.  &
                  (    obs% spot(k2)% hd%  obstype == OT_PILOT .and. &
                   all(obs% spot(k2)% hd% codetype /= ct_wprof)      )
    end do
    nrep = count (mask)
    allocate (idx(nrep))
    idx(:nrep) = pack ((/ (k2, k2=1,obs% n_spot) /), mask)
    allocate (idx1 (nrep))
    allocate (idx2 (nrep))
    allocate (prio (nrep))
    allocate (down (nrep))
    down = .false.
    !----------------------------------------
    ! Hash station names for fast comparisons
    !----------------------------------------
    do k2 = 1, nrep
       i2 = idx(k2)
       spt2 => obs% spot(i2)
       down(k2) = any (spt2% hd% dbkz == kz_tdrop)
       if (down(k2)) then
!!!       if (.not. spt2% wsi% valid) & !!! TODO: WSI for dropsondes?
          spt2% stat_hash = transfer ("DROP    "  , spt2% stat_hash)
!      else
!         spt2% stat_hash = transfer (spt2% statid, spt2% stat_hash)
       end if
       select case (obs% spot(i2)% hd% dbkz)
       case (1:999)               ! old format A/B/C/D to merge
         prio (k2) = 1
       case (10520, 10776, 10780) ! new TEMP/SHIP/DROP BUFR format, thinned
         prio (k2) = 4
       case (10521, 10777)        ! new TEMP/SHIP BUFR format (up to 100 hPa)
         prio (k2) = 2
       case (10526, 10782, 10574, 10785, 10516, 10570) ! TEMP BUFR, full res.
         prio (k2) = 5
       case (10527, 10783, 10517) ! new TEMP BUFR format (up to 100 hPa)
         prio (k2) = 3
       case default               ! unknown
         prio (k2) = 0
       end select
    end do
    !------------------
    ! loop over reports
    !------------------
    do k2 = 1, nrep
      i2   = idx(k2)
      spt2 => obs% spot(i2)
      if (spt2% use% state <= STAT_DISMISS) cycle
      if (spt2% p% n       <= 0           ) cycle
      call load_temp (obs, spt2, t2)
      changed = .false.
      !--------------------------------------------
      ! Setup index vector for detailed comparisons
      !--------------------------------------------
      n_idx = 0
      do k1 = k2+1, nrep
         i1 = idx(k1)
         if (obs%spot(i1)%      stat_hash == spt2% stat_hash .and. &
             obs%spot(i1)%    hd% obstype == spt2% hd% obstype     ) then
            n_idx = n_idx + 1
            idx1(n_idx) = idx(k1)
            idx2(n_idx) =     k1
         end if
      end do
      !-------------------------------
      ! loop over reports stored later
      !-------------------------------
      first = .true.
      do k1 = 1,n_idx
        i1 = idx1(k1)
        j1 = idx2(k1)
        spt1 => obs% spot(i1)
        !-------------------------
        ! ignore dismissed reports
        !-------------------------
        if (spt1% use% state <= STAT_DISMISS) cycle
        !---------------------------------------------------------
        ! check for corresponding station name, subtype, time.
        ! coordinates may differ due to different accuracy stored,
        ! time may differ for BUFR format compared with other
        ! (dbkz>=10000 provides seconds, not yet used here.)
        ! Dropsonde IDs may have aircraft registration or "DROP".
        !---------------------------------------------------------
        ldown   = down(j1) .and. down(k2)
        same_id = ldown .or. (spt1% statid == spt2% statid)
        !----------------------------------------------------------------------
        ! "Same station" check for stations with valid, hashed WIGOS station ID
        !----------------------------------------------------------------------
        if (same_id .and. spt1% wsi% valid .and. spt1% statid(1:1) == "_") then
           same_id = spt1% wsi == spt2% wsi
        end if

        if    (same_id                                              .and. &
               spt1%    hd% obstype == spt2%    hd% obstype         .and. &
          abs (spt1% col%c% dlon    -  spt2% col%c% dlon) < 0.10_wp .and. &! tship
          abs (spt1% col%c% dlat    -  spt2% col%c% dlat) < 0.10_wp .and. &! tship
        ((any (spt1%    hd% dbkz    == kz_pilot)                    .and. &
          any (spt2%    hd% dbkz    == kz_pilot))                   .or.  &
         (any (spt1%    hd% dbkz    == kz_temp )                    .and. &
          any (spt2%    hd% dbkz    == kz_temp ))                   .or.  &
         (any (spt1%    hd% dbkz    == kz_tship)                    .and. &
          any (spt2%    hd% dbkz    == kz_tship))                   .or.  &
         (any (spt1%    hd% dbkz    == kz_tmobil)                   .and. &
          any (spt2%    hd% dbkz    == kz_tmobil))                  .or.  &
         (any (spt1%    hd% dbkz    == kz_tdrop)                    .and. &
          any (spt2%    hd% dbkz    == kz_tdrop))                   .or.  &
         (any (spt1%    hd% dbkz    == kz_pship)                    .and. &
          any (spt2%    hd% dbkz    == kz_pship)))                        ) then

         same_time = spt1% hd% time    == spt2% hd% time          ! Reference
         if (ldown) then
            dt        = abs (minutes (spt1% actual_time - spt2% actual_time))
            same_time = dt < 2._wp
            lmerge    = same_time .or. prio(j1)      /=  prio(k2)
         else
            if (same_time .and. corme_mode >= 4)             &    ! +
            same_time = spt1% actual_time == spt2% actual_time    ! Launch time
            lmerge    = same_time .or. prio(j1) > 1 .or. prio(k2) > 1
         end if

         if (lmerge) then

          !--------------------------------------------------
          ! Fix issues with RS type coded only in TEMP part B
          ! fix similar issues for measurement equipment type
          ! for solar and infrared radiation correction, and
          ! for tracking technique.
          !--------------------------------------------------
          if   (spt1% hd% obstype == OT_TEMP) then
            if (spt1% sttyp * spt2% sttyp < 0) then
             !print*,"statid,spt1,spt2:",spt1%statid,spt1%hd%dbkz,spt1%sttyp,spt2%hd%dbkz,spt2%sttyp
              if (spt1% sttyp < 0) then
                 spt1% sttyp = spt2% sttyp
              else
                 spt2% sttyp = spt1% sttyp
              end if
            end if
            if (spt1% meas_type * spt2% meas_type < 0) then
              if (spt1% meas_type < 0) then
                 spt1% meas_type = spt2% meas_type
              else
                 spt2% meas_type = spt1% meas_type
              end if
            end if
            if (spt1% stret * spt2% stret < 0) then
              if (spt1% stret < 0) then
                 spt1% stret = spt2% stret
              else
                 spt2% stret = spt1% stret
              end if
            end if
            if (spt1% tracking * spt2% tracking < 0) then
              if (spt1% tracking < 0) then
                 spt1% tracking = spt2% tracking
              else
                 spt2% tracking = spt1% tracking
              end if
            end if
          end if
          !----------------------------------
          ! handling for new TEMP BUFR format
          !----------------------------------
          if (prio(j1) == 1 .and. prio(k2) == 1) then
            ! merge parts A/B/C/D
          else if (prio(j1) == 0 .or. prio(k2) == 0) then
            write (0,*) 'priorities:',prio(j1),prio(k2),      &
                        ', dbkz:',spt1% hd% dbkz,spt2% hd% dbkz
            call finish ('merge_report','prio(j1) == 0 .or. prio(k2) == 0')
          else if (prio(j1) > prio(k2)) then

            write (com,'("corme,prio,dbkz:",2i2,i6," - ",2i2,i6)') &
              spt2% corme, prio(k2), spt2% hd% dbkz,               &
              spt1% corme, prio(j1), spt1% hd% dbkz

            call decr_rpt_use (spt2, CHK_REDUNDANT, comment=com)    ! keep last
            exit
          else if (prio(j1) < prio(k2)) then

            write (com,'("corme,prio,dbkz:",2i2,i6," - ",2i2,i6)') &
              spt1% corme, prio(j1), spt1% hd% dbkz,               &
              spt2% corme, prio(k2), spt2% hd% dbkz

            call decr_rpt_use (spt1, CHK_REDUNDANT, comment=com)    ! keep first
            cycle
          else if (corme_mode <= 1) then

            write (com,'("corme,prio,dbkz:",2i2,i6," - ",2i2,i6)') &
              spt2% corme, prio(k2), spt2% hd% dbkz,               &
              spt1% corme, prio(j1), spt1% hd% dbkz

            call decr_rpt_use (spt2, CHK_REDUNDANT, comment=com)    ! keep last
            exit

          else
            ! goto merge
          endif

          !---------
          ! printout
          !---------

          if (first) write(6,'(a)') '#MERGE'
          first = .false.
          write(6,&
'(a,i4,1x,a,a," - ",i3,i8,i3,i8," - ",2(a,1x,a,l2," - "),a,1x,a," - ",2i6," - ",2i6)')&
            '#MERGE',spt2% hd% obstype,spt2% statid, spt1% statid,          &
            spt2% hd% source, spt2% hd% record,                             &
            spt1% hd% source, spt1% hd% record,                             &
            chhmmss(spt2% actual_time),chhmmss(spt1% actual_time),          &
            spt2% actual_time==spt1% actual_time,                           &
            chhmmss(spt2% hd%    time),chhmmss(spt1% hd%    time),          &
            spt2% hd%    time==spt1% hd%    time,                           &
            chhmmss(spt2% hd% db_time),chhmmss(spt1% hd% db_time),          &
            spt2% hd% dbkz, spt1% hd% dbkz, spt2% o% n, spt1% o% n

          !--------------
          ! merge reports
          !--------------
          call load_temp (obs, spt1, t1)
          call merge_temp (t2, t1)
          changed = .true.

          write (com,'("corme,prio,dbkz:",2i2,i6," - ",2i2,i6)') &
            spt1% corme, prio(j1), spt1% hd% dbkz,               &
            spt2% corme, prio(k2), spt2% hd% dbkz

          if      (spt2% hd% dbkz <= spt1% hd% dbkz) then
            call decr_rpt_use (spt1, CHK_MERGE, comment=com)
          else
            call decr_rpt_use (spt2, CHK_MERGE, comment=com)
            spt2 => spt1
            i2   =    i1
          endif
          deallocate (t1)

         end if ! same_time .or. prio()...
        end if
      end do
      if (changed) then
        spt = spt2
        call check_store_temp (t2, spt, obs, nonew, repl = i2)
        write(6,'(a,i4,1x,a,a,2i6)') &
          '#MERGE',spt2% hd% obstype,spt2% statid,'  new: ',nonew, spt2% o% n
        if (prt_data) then
           print '(a,a,a)', "### Merged TEMP ", spt2% statid, ":"
           call print_temp_data (t2)
        end if
      endif
      deallocate (t2)
    end do

  end subroutine merge_report

!------------------------------------------------------------------------------
  subroutine split_report (obs, grid)
  type (t_obs)  ,intent(inout) :: obs      ! observation container
  type (t_grid) ,intent(in)    :: grid     ! model grid (for grid spacing)
  !----------------------------------------------------------
  ! splits TEMP-reports to account for drift out of grid cell
  !----------------------------------------------------------

    integer     ,parameter    :: npmax = 4    ! max.# of neighbour gridpoints
    type(t_spot),pointer      :: spt2         ! reference to 'old' report
    type(t_spot)              :: s            ! local copy of spt2
    type(t_spot)              :: sptnew       ! new report to store
    type(t_temp),pointer      :: t1 (:)       ! new TEMP specific obs. data
    type(t_temp),pointer      :: t2 (:)       ! old TEMP specific obs. data
    type(t_datum),allocatable :: bdy(:)       ! report body
    integer                   :: i            ! observation index within report
    integer                   :: i1           ! first index within report
    integer                   :: in           ! last  index within report
    integer                   :: is           ! report index within box
    integer                   :: ns           ! new report index within box
    integer                   :: idx1(npmax,4)! neighbour grid point indices
    integer                   :: idx2(npmax,4)! neighbour grid point indices
    real(wp)                  :: w(npmax)     ! interpol. weights (not used)
    real(wp)                  :: lat1, lon1   ! coordinates of current report
    real(wp)                  :: lat2, lon2   ! coordinates of current obs.
    integer                   :: no1          ! observation indices (in t_temp)
    integer                   :: no2          ! observation indices (in t_temp)
    integer                   :: nonew        ! number of obs. stored
    integer                   :: noold        ! number of obs. before storage
    integer                   :: np           ! number of neighbours
    integer                   :: np1          ! number of neighbours
    logical                   :: first        ! flag for first report
    integer                   :: rtype        ! report type
    integer                   :: ctype        ! code type
    integer                   :: idbk         ! data base id slot
    integer                   :: status       ! report status

    !--------------------
    ! check namelist flag
    !--------------------
    if (.not. split) return

    !--------------------------------
    ! loop over TEMP or PILOT reports
    !--------------------------------
    do is = 1, obs% n_spot
      rtype  = obs% spot(is)% hd % obstype
      ctype  = obs% spot(is)% hd % codetype
      idbk   = obs% spot(is)% hd % idbk
      status = obs% spot(is)% use% state
      if (.not. (rtype == OT_TEMP  .or.  &
                (rtype == OT_PILOT .and. &
             all(ctype /= ct_wprof) ) )  ) cycle
      if (status <= STAT_DISMISS)          cycle
      spt2  => obs% spot(is)
      s     =  spt2
      first = .true.

      !------------------------------------
      ! skip report if displacement is zero
      !------------------------------------
      i1    = s% o% i + 1
      in    = s% o% i + s% o% n
      lat1  = obs% body (i1)% lat
      lon1  = obs% body (i1)% lon
      lat2  = obs% body (in)% lat
      lon2  = obs% body (in)% lon
      if (lat1 == lat2 .and. lon1 == lon2) cycle

      if (prt_data) print *, "split_report: ", trim (spt2% statid)
      !-------------------------------------------
      ! skip report if any displacement is invalid
      !-------------------------------------------
      if (any (obs% body (i1:in)% lat == rvind)) cycle
      if (any (obs% body (i1:in)% lon == rvind)) cycle

      !-------------------------------------------
      ! find neighbour grid-points for first entry
      ! split only if lowest level inside domain
      !-------------------------------------------
      call Grid_Indices &
           (lon1,       & ! <-- geodetic longitude
            lat1,       & ! <-- geodetic latitude
            grid,       & ! <-- grid data type
            idx1,       & ! --> Grid point indices [Point, index]
            w,          & ! --> Weight
            np1)          ! --> number of points returned

      if (np1 == 0) then
         if (prt_data) print *, "split_report: ", trim (spt2% statid), &
              " no splitting since first point outside domain"
         cycle
      end if

      !---------------------------------------------
      ! find neighbour grid-points for last entry.
      ! if level in same horizontal grid-cell: cycle
      ! if outside domain, check for usable segments
      !---------------------------------------------
      call Grid_Indices &
          (lon2,        & ! <-- geodetic longitude
           lat2,        & ! <-- geodetic latitude
           grid,        & ! <-- grid data type
           idx2,        & ! --> Grid point indices [Point, index]
           w,           & ! --> Weight
           np           ) ! --> number of points returned

      if (np /= 0 .and. all(idx1(1:np,:) == idx2(1:np,:))) then
         if (prt_data) print *, "split_report: ", trim (spt2% statid), &
              " no splitting since start and endpoint in same grid-cell"
         cycle
      end if

      !---------------------------------------------------
      ! if TEMP specific parameters are not present: cycle
      !---------------------------------------------------
      if (spt2% p% n <= 0) then
        call message ('split_report',                                         &
                      'TEMP specific parameters missing: '//trim(spt2% statid))
        cycle
      endif

      !----------------------
      ! keep old body entries
      !----------------------
      allocate  (bdy  (i1:in))
      bdy = obs% body (i1:in)

      !--------------------------------------
      ! load TEMP specific observational data
      !--------------------------------------
      call load_temp (obs, spt2, t2)
      no2   = 1
      no1   = 1

      if (prt_data) then
         print *, "split_report: consider splitting report ", &
              trim (spt2% statid), " with", size (t2), "levels"
      end if
      !-----------------------------
      ! loop over subsequent entries
      !-----------------------------
      do i = i1+1, in+1
        if (i <= in) then
          !------------------
          ! same level: cycle
          !------------------
          if (obs% olev(i-1) == obs% olev(i)) cycle
          !-----------------------------------
          ! new level: check for new grid-cell
          !-----------------------------------
          !---------------------------------------------
          ! find neighbour grid-points for current level
          !---------------------------------------------
          lat2 = obs% body (i)% lat
          lon2 = obs% body (i)% lon
          call Grid_Indices &
              (lon2,        & ! <-- geodetic longitude
               lat2,        & ! <-- geodetic latitude
               grid,        & ! <-- grid data type
               idx2,        & ! --> Grid point indices [Point, index]
               w,           & ! --> Weight
               np           ) ! --> number of points returned

          !------------------------------------------
          ! level in same horizontal grid-cell: cycle
          !------------------------------------------
          if (np > 0 .and. np1 > 0 .and. all(idx1(1:np,:) == idx2(1:np,:))) then
            no2 = no2 + 1
            cycle
          endif
          !-------------------------------------------------------------
          ! consecutive levels outside domain: gather in single segment;
          ! domain boundary crossed if np==0 && np1>0 || np1==0 && np>0
          !-------------------------------------------------------------
          if (np == 0 .and. np1 == 0) then
            no2 = no2 + 1
            cycle
          end if
        endif

        if (prt_data) print *, "split_report: splitting levels", no1,"...",no2
        !-----------------------------------------------
        ! level in new report (or end of report) : split
        !-----------------------------------------------
        t1 => t2(no1:no2)
        sptnew = s
        sptnew% col% c% dlon = obs% body (i1)% lon
        sptnew% col% c% dlat = obs% body (i1)% lat

        ! safeguard against displacement data out of domain
        if (np1 == 0) call decr_rpt_use (sptnew, CHK_DOMAIN)

        if (first) then
          call check_store_temp (t1, sptnew, obs, nonew, repl=is, pass=2)
          ns = is
        else
          call check_store_temp (t1, sptnew, obs, nonew         , pass=2)
          ns = obs% n_spot
        endif

        !------------------------------
        ! check report size consistency
        !------------------------------
        noold = i - i1
        if (nonew /= noold) then
          write(*,*) "split_report: statid, nonew /= noold : ", &
               s% statid, nonew, noold
          call print_temp_data (t1)
          call finish('split_report','nonew /= noold, statid= ' // s% statid)
        endif

        !----------------------
        ! keep old status flags
        !----------------------
        obs% body (obs% spot(ns)% o% i + 1 :                  &
                   obs% spot(ns)% o% i + obs% spot(ns)% o% n) &
          =  bdy  (i1:i1-1+nonew)

        !-------------------------
        ! update report statistics
        !-------------------------
        if (.not.first) then
          call update_statistics (rept_stat (0,   rtype), &
                                  rept_stat (idbk,rtype), &
                                  status,             1   )
          call update_statistics (rept_stat (0,   rtype), &
                                  rept_stat (idbk,rtype), &
                                  STAT_FORGET,       -1   )
        endif

        !----------------
        ! update counters
        !----------------
        no2   = no2 + 1
        no1   = no2
        i1    = i
        idx1  = idx2
        np1   = np      ! remember whether inside/outside domain
        first = .false.
      end do

      !------------------------------
      ! deallocate TEMP specific data
      !------------------------------
      deallocate (t2)
      deallocate(bdy)
    end do

  end subroutine split_report
!------------------------------------------------------------------------------
  elemental function vss_from_evss (ievss) result (ivss)
  integer ,intent(in) :: ievss
  integer             :: ivss
  !-----------------------------------------------
  ! derive          vertical sounding significance
  !   from Extended vertical sounding significance
  !-----------------------------------------------
    ivss   = 0
    if   (iand (ievss, lev_surf_18 ) /= 0)  ivss = ibset (ivss , lbit_surf )
    if   (iand (ievss, lev_std_18  ) /= 0)  ivss = ibset (ivss , lbit_std  )
    if   (iand (ievss, lev_tropo_18) /= 0)  ivss = ibset (ivss , lbit_tropo)
    if   (iand (ievss, lev_max_v_18) /= 0)  ivss = ibset (ivss , lbit_max_v)
    if   (iand (ievss, lev_sig_t_18) /= 0)  ivss = ibset (ivss , lbit_sig_t)
    if   (iand (ievss, lev_sig_h_18) /= 0)  ivss = ibset (ivss , lbit_sig_t)
    if   (iand (ievss, lev_sig_v_18) /= 0)  ivss = ibset (ivss , lbit_sig_v)
    if   (iand (ievss, lev_regn_18 ) /= 0)  ivss = ibset (ivss , lbit_other)
    if   (iand (ievss, lev_miss_18 ) /= 0)  ivss = ibset (ivss , lbit_miss )
    if   (      ivss                 == 0)  ivss = lev_other

  end function vss_from_evss
!------------------------------------------------------------------------------
  pure subroutine ffvss_from_evss (temp)
  type (t_temp) ,intent(inout) :: temp
  !----------------------------------------------------
  ! derive Feedback-File vertical sounding significance
  !   from Extended      vertical sounding significance
  !----------------------------------------------------

    !---------------
    ! missing values
    !---------------
    if (iand (temp% vss, lev_miss_18) /= 0) then
      temp% p % lev_sig = LEV_MISSING
      temp% t % lev_sig = LEV_MISSING
      temp% td% lev_sig = LEV_MISSING
      temp% ff% lev_sig = LEV_MISSING
      goto 99
    endif

    !-------------
    ! other values
    !-------------
    temp% p  % lev_sig = 0_i2
    temp% t  % lev_sig = 0_i2
    temp% td % lev_sig = 0_i2
    temp% ff % lev_sig = 0_i2

    if (iand (temp% vss, lev_surf_18 ) /= 0) then
              temp% p % lev_sig = ibset (temp% p % lev_sig, LS_SURFACE)
              temp% t % lev_sig = ibset (temp% t % lev_sig, LS_SURFACE)
              temp% td% lev_sig = ibset (temp% td% lev_sig, LS_SURFACE)
              temp% ff% lev_sig = ibset (temp% ff% lev_sig, LS_SURFACE)
    endif
    if (iand (temp% vss, lev_std_18  ) /= 0) then
              temp% p % lev_sig = ibset (temp% p % lev_sig, LS_STANDARD)
              temp% t % lev_sig = ibset (temp% t % lev_sig, LS_STANDARD)
              temp% td% lev_sig = ibset (temp% td% lev_sig, LS_STANDARD)
              temp% ff% lev_sig = ibset (temp% ff% lev_sig, LS_STANDARD)
    endif
    if (iand (temp% vss, lev_tropo_18) /= 0) then
              temp% p % lev_sig = ibset (temp% p % lev_sig, LS_TROPO)
              temp% t % lev_sig = ibset (temp% t % lev_sig, LS_TROPO)
              temp% td% lev_sig = ibset (temp% td% lev_sig, LS_TROPO)
              temp% ff% lev_sig = ibset (temp% ff% lev_sig, LS_TROPO)
    endif
    if (iand (temp% vss, lev_max_v_18) /= 0) then
              temp% ff% lev_sig = ibset (temp% ff% lev_sig, LS_MAX)
    endif
    if (iand (temp% vss, lev_sig_t_18) /= 0) then
              temp% t % lev_sig = ibset (temp% t % lev_sig, LS_SIGN)
    endif
    if (iand (temp% vss, lev_sig_h_18) /= 0) then
              temp% td% lev_sig = ibset (temp% td% lev_sig, LS_SIGN)
    endif
    if (iand (temp% vss, lev_sig_v_18) /= 0) then
              temp% ff% lev_sig = ibset (temp% ff% lev_sig, LS_SIGN)
    endif

    !-----------------------------
    ! flags for derived quantities
    !-----------------------------
99  temp% rh % lev_sig= temp% td% lev_sig
    temp% dd % lev_sig= temp% ff% lev_sig
    temp% uu % lev_sig= temp% ff% lev_sig
    temp% vv % lev_sig= temp% ff% lev_sig
    temp% gp % lev_sig= temp% p % lev_sig
  end subroutine ffvss_from_evss
!------------------------------------------------------------------------------
  pure subroutine ffvss_from_vss (temp)
  type (t_temp) ,intent(inout) :: temp
  !----------------------------------------------------
  ! derive Feedback-File vertical sounding significance
  !   from               vertical sounding significance
  !----------------------------------------------------
    temp% p  % lev_sig = 0_i2
    temp% t  % lev_sig = 0_i2
    temp% td % lev_sig = 0_i2
    temp% ff % lev_sig = 0_i2
    temp% rh % lev_sig = 0_i2
    temp% dd % lev_sig = 0_i2
    temp% uu % lev_sig = 0_i2
    temp% vv % lev_sig = 0_i2
    temp% gp % lev_sig = 0_i2
    if (iand (temp% vss, lev_miss ) /= 0) then
      temp% vss = lev_miss
      return
    endif
    if (iand (temp% vss, lev_surf ) /= 0) then
      temp% p % lev_sig = ibset (temp% p % lev_sig, LS_SURFACE)
      temp% t % lev_sig = ibset (temp% t % lev_sig, LS_SURFACE)
      temp% td% lev_sig = ibset (temp% td% lev_sig, LS_SURFACE)
      temp% ff% lev_sig = ibset (temp% ff% lev_sig, LS_SURFACE)
    endif
    if (iand (temp% vss, lev_std  ) /= 0) then
      temp% p % lev_sig = ibset (temp% p % lev_sig, LS_STANDARD)
      temp% t % lev_sig = ibset (temp% t % lev_sig, LS_STANDARD)
      temp% td% lev_sig = ibset (temp% td% lev_sig, LS_STANDARD)
      temp% ff% lev_sig = ibset (temp% ff% lev_sig, LS_STANDARD)
    endif
    if (iand (temp% vss, lev_tropo) /= 0) then
      temp% p % lev_sig = ibset (temp% p % lev_sig, LS_TROPO)
      temp% t % lev_sig = ibset (temp% t % lev_sig, LS_TROPO)
      temp% td% lev_sig = ibset (temp% td% lev_sig, LS_TROPO)
      temp% ff% lev_sig = ibset (temp% ff% lev_sig, LS_TROPO)
    endif
    if (iand (temp% vss, lev_max_v) /= 0) then
      temp% ff% lev_sig = ibset (temp% ff% lev_sig, LS_MAX)
    endif
    if (iand (temp% vss, lev_sig_t) /= 0) then
      temp% t % lev_sig = ibset (temp% t % lev_sig, LS_SIGN)
      temp% td% lev_sig = ibset (temp% td% lev_sig, LS_SIGN)
    endif
    if (iand (temp% vss, lev_sig_v) /= 0) then
      temp% ff% lev_sig = ibset (temp% ff% lev_sig, LS_SIGN)
    endif
    temp% rh % lev_sig= temp% td% lev_sig
    temp% dd % lev_sig= temp% ff% lev_sig
    temp% uu % lev_sig= temp% ff% lev_sig
    temp% vv % lev_sig= temp% ff% lev_sig
    temp% gp % lev_sig= temp% p % lev_sig
  end subroutine ffvss_from_vss
!------------------------------------------------------------------------------
! subroutine check_temp_cons (obs)
!   type(t_obs), intent(inout) :: obs         ! observation container
    !--------------------------------------------------
    ! Check consistency of observational data in report
    !--------------------------------------------------
! end subroutine check_temp_cons
!==============================================================================
end module mo_temp
