!
!+ Routines specific to the aircraft observation operator
!
MODULE mo_airep
!
! Description:
!   Routines specific to the aircraft observation operator
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
!  changes for verification
! V1_5         2009/05/25 Andreas Rhodin
!  Changes for reading observations from feedback-file
! V1_7         2009/08/24 Harald Anlauf
!  Remove unused variables
! V1_8         2009/12/09 Harald Anlauf
!  namelist airep_obs: add chk_phase, chk_rollangle
!                      for checks of report quality
!  check_store_airep : set reports on latitude 0 to PASSIVE
!                      drop bogus reports at (0N,0E)
! V1_9         2010/04/20 Andreas Rhodin
!  TSK_SHRINK in subroutines process: pass parameter 'state' to shrink_report
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_11        2010/06/16 Harald Anlauf
!  read_airep_netcdf: reduce verbosity of output for netcdf_verb==0
! V1_13        2011/11/01 Andreas Rhodin
!  adaptions for bias correction
! V1_15        2011/12/06 Andreas Rhodin
!  option to specify codetypes for observation error tables
! V1_20        2012-06-18 Harald Anlauf
!  properly initialize len_level(0),len_report(0)
! V1_22        2013-02-13 Alexander Cress
!  re-derive phase (ascend/descend/level) for bias correction
! V1_26        2013/06/27 Andreas Rhodin
!  fix fault caused by zero Jakobians drh/dgh (for version=2 in /cntrvar/)
! V1_29        2014/04/02 Harald Anlauf
!  Handle MODE-S data, part 2: properly deal with humidity undefined data
! V1_41        2015-03-05 Andreas Rhodin
!  fix for zero pressure level reported by aircrafts
! V1_42        2015-06-08 Andreas Rhodin
!  implement temporal interpolation for COSMO MEC
! V1_47        2016-06-06 Andreas Rhodin
!  add obstype as selection criterium in namelist /rules/
! V1_48        2016-10-06 Andreas Rhodin
!  implement FF,DD passive monitoring for aircraft observations
! V1_49        2016-10-25 Harald Anlauf
!  fix uninitialized nfnfn bug
! V1_51        2017-02-24 Andreas Rhodin
!  TSK_R: check for invalid observations, multilevel reports
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Andreas Rhodin  DWD  2004-2008  original source
! Oliver Schmid   DWD  2005/06    new obs data type
! Harald Anlauf   DWD  2008       optimizations for SX8
! Gerhard Paul    DWD  2008       NetCDF input
!-----------------------------------------------------------

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
  use mo_exception, only: finish           ! abort routine
  use mo_kind,      only: wp, sp, i1, i8   ! kind parameters
  use mo_mpi_dace,  only: dace,           &! MPI group info
                          p_bcast          ! generic MPI bcast routine
  use mo_obs_set,   only: t_obs_block      ! observation data type
  use mo_usstd,     only: p_h_usstd        ! p from h(gpm) US std.atm.
  use mo_namelist,  only: position_nml,   &! position namelist
                          nnml,           &! namelist Fortran unit number
                          POSITIONED       ! ok    code from position_nml
  use mo_run_params,only: flag_biasc_airep ! aircraft bias correction
  !-----------------------------
  ! access atmospheric data type
  !-----------------------------
  use mo_atm_state, only: t_atm            ! atm. state data type
  use mo_t_col,     only: t_cols,         &! model columns data type
                          COL_UV,         &! specification of fields in column
                          COL_TV, COL_RH   !
  use mo_physics,   only: d2r,            &! pi/180.
                          r2d,            &! 180./pi
                          rhw_m,          &! rel.hum.(water) from mixing ratio
                          rhw_m_hardy,    &!  " (ITS90 formulation, Hardy)
                          fd_uv            ! calculate ff,dd from u,v
  use mo_cntrlvar,  only: trh_tvgh         ! generalized humidity transform
  !-----------------------------
  ! access observation data type
  !-----------------------------
  use mo_t_obs,     only: t_obs,          &!
                          t_spot,         &!
                          t_head,         &! observation data type
                          derive_dbkz,    &! derive DBKZ if not present
                          new_spot,       &! reserve memory
                          new_obs,        &! reserve memory
                          set_xuv,        &! set unit vectors, solar zenith
                          invalid,        &! invalid observation value
                          set_vqc_insitu, &! subroutine to set VQC bounds
                          set_int_insitu, &! set interpolation space
                          shrink_report,  &! remove passive observations
                          source,         &! list   of Report source files
                          monitor_ff,     &! flag to monitor wind speed
                          monitor_dd,     &! flag to monitor wind direction
                          AIREP,          &! flag to process by this module
                          TSK_INIT,       &!  initialisation
                          TSK_READ,       &!  read observations
                          TSK_SET_CHR,    &!  set observation characteristics
                          TSK_SETUP_COLS, &!  determine model columns required
                          TSK_SETUP_FUL0, &!  setup interpolation space
                          TSK_SETUP_FULL, &!  setup description of PSAS-space
                          TSK_SHRINK,     &!  release unused obs. in report
                          TSK_K,          &!  set up linear operator
                          TSK_Y,          &!  evaluate nonlinear operator
                          TSK_R,          &!   setup observational errors
                          tsk_name,       &! derive mnemonic for TSK_.. flags
                          CHR_ID,         &! H is the identity operator
                          CHR_NONL,       &! H is nonlinear
                          CHR_INV,        &! H is invertible
                          CHR_EXP,        &! H is 'expensive'
                          ITY_ICOL,       &! interpolation type: column
                          netcdf_verb      ! verbosity of NetCDF decoding
  use mo_fdbk_tables,only:VN_U,           &!          wind component code
                          VN_V,           &!                         code
                          VN_FF, VN_DD,   &! wind direction, speed   code
                          VN_RH,          &!          rh             code
                          VN_T,           &!          T              code
                          VN_VGUST,       &!          VGUST          code
                          VN_TURB,        &!          TURB           code
                          VN_FLEV,        &! nominal flight level code
                          VN_P,           &! pressure (level) code
                          OT_AIREP         ! observation type id
  use mo_wmo_tables,only: WMO0_ECMWF       ! generating center
  use mo_t_use,     only: decr_use,       &! decrease the state of a datum
                          t_use,          &! status variable data type
                          use_0,          &! default values of type use
                          STAT_DISMISS,   &!
                          STAT_PASSIVE,   &!
                          CHK_INSDAT,     &!
                          CHK_CORR,       &!
                          CHK_DOMAIN,     &!
                          CHK_CORRERR,    &!
                          CHK_QI,         &! quality index
                          CHK_NOTUSED      !
  use mo_obs_tables,only: check_report_0, &! init. flags, standard checks
                          check_report_1, &! standard checks
                          rept_use,       &! report type usage table
                          decr_rpt_use,   &! change use-flags of report
                          idb_dbk          ! index in table rept_stat
!                         rept_char        ! observation type characteristics
  use mo_obs_err,   only: obs_err          ! get observation error from table
  use mo_time,      only: init_time,      &! initialise time data type
                          cyyyymmddhhmm,  &! derive string from time
                          cyyyymmddhhmmss,&! derive string from time
                          operator (==)    ! compare times
  use mo_t_datum,   only: t_datum,        &! data type for one observed datum
                          rvind,          &! missing value indicator (real)
!#if defined(__SX__)
!                         inv_datum,      &! invalid datum
!#endif
                          SRC_DER,        &! derived quantity flag value
                          QC_OK,          &! QC OK      flag value
                          QC_NOUSE,       &! DO NOT USE flag value
                          QC_MISS,        &! missing or invalid value
                          QC_CLIM,        &! out of climatological range
                          set_datum,      &! set t_datum% o  (observed value)
                          set_qbits        ! set t_datum% qc (quality bits)
  use mo_bufr_dwd,  only: t_bufr,         &! BUFR record data type
                          bufr_get_character,   &!
#ifdef CHECKCODES
                          bufr_get_entry_texts, &!
                          bufr_get_entry_units, &!
#endif
                          bufr_print_sections,  &!
                          bufr_print_subset,    &!
                          inv_bufr         ! indicator for invalid value
  use mo_obs_rules, only: get_rule,       &! routine to get a rule
                          iud,            &! undefined integer value
                          rud,            &! undefined real    value
                          t_set            ! result data type

  use mo_obstypes,  only: t_obsid,        &! observation id table entry
                          obstype_dbkz,   &! derive obsids from dbkz
                          obstype_bufr     ! derive obsids from bufr type
  use mo_conv_bc   ,only: biascor_mode,   &! mode used for bias correction
                          t_decay,        &! biasc. accum. decay time (days)
                          force_level,    &! force phase to level above
                          fr_land_bc,     &! land fraction req. for bc stats
                          n_required,     &! number of entries required for bc
                          bc_fallback,    &! action if biasc-file not present
                          scan_tracks,    &! scan tracks, set phase
                          time_sep,       &! time separation (s)
                          rate_asc,       &! ascend  rate required
                          rate_desc,      &! descend rate required
                          rate_lev,       &! level   rate required
                          BC_NOBC, BC_FG   ! values for 'biascor_mode'
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
                          s1updat,        &! update sequence no.
                          mlah,           &! latitude
                          mloh,           &! longitude
                          obs_time,       &! body observation time
                          ystidn           ! any type of station identifier as variable
  !---------------------
  ! netCDF f90 interface
  !---------------------
  use netcdf,       only: nf90_Inquire_Dimension,&!
                          nf90_Inquire_Variable, &!
                          nf90_inq_varid,        &!
                          nf90_get_var,          &!
                          NF90_FLOAT,            &!
                          NF90_INT,              &!
                          NF90_NOERR              !
  !------------------------------------
  ! Interface to interpolation routines
  !------------------------------------
  use mo_grid_intpol,only: idx_init
  !------------------------
  ! access matrix data type
  !------------------------
  use mo_dec_matrix,only: t_vector_segm    ! vector segment data type
  implicit none
!------------------------------------------------------------------------------
  !================
  ! public entities
  !================
  private
  public :: process_airep     ! general purpose AIREP processing routine
  public :: t_atm             ! argument data type to process_airep
  public :: t_cols            ! argument data type to process_airep
  public :: read_airep_bufr   ! read AIREP observation from BUFR record
  public :: read_airep_netcdf ! read AIREP observation from netCDF file
  public :: check_store_airep ! store accepted data in observation data type
  public :: t_airep           ! type for temporary storage of AIREP observats.
  public :: read_airep_nml    ! read namelist /AIREP_OBS/
!------------------------------------------------------------------------------
  !======================
  ! data type definitions
  !======================
  !---------------------------------------------------
  ! type t_airep: temporarily store AIREP observations
  !---------------------------------------------------
  type t_airep
    type (t_datum) :: total ! quality control flag
    type (t_datum) :: p     ! pressure level             [Pa]
    type (t_datum) :: z     ! altitude                   [m]
    type (t_datum) :: t     ! temperature                [K]
    type (t_datum) :: ff    ! wind speed                 [m/s]
    type (t_datum) :: dd    ! wind direction             [degree]
    type (t_datum) :: rh    ! relative humidity          [%]
    type (t_datum) :: mr    ! mixing ratio               [kg/kg]
  end type t_airep

  !-------------------
  ! Namelist AIREP_OBS
  !-------------------
  logical  :: black_t_nop   = .false. ! blacklist t-obs. if pressure is missing
  logical  :: use_regnum    = .true.  ! use registration number, not flight no.
  logical  :: use_tailnum   = .true.  ! use aircraft tail number, not flight no.
  logical  :: require_q     = .false. ! require q for sufficient data
  logical  :: reject_q_0    = .true.  ! reject zero value humidity observations
  logical  :: bug_multilev  = .false. ! emulate bug reading multilevel reports?
  integer  :: chk_phase     = 0       ! check for flight phase
  integer  :: chk_rollangle = 0       ! check for roll-angle quality
  integer  :: flev_mode     = 0       ! flight level reporting: 1=p,2=flev
  integer  :: satvp_form    = 1       ! Formula for sat. vapor pres. over water
                                      ! 1: Magnus-Tetens, 2: Hardy
  real(wp) :: top_q         = 250._wp ! top level for humidity observations(hPa)
  namelist  /AIREP_OBS/ black_t_nop, use_regnum,use_tailnum,&
                        require_q, reject_q_0, bug_multilev,&
                        chk_phase, chk_rollangle, flev_mode,&
                        biascor_mode, t_decay, force_level, &! for bias corr.
                        fr_land_bc, n_required, bc_fallback,&!
                        scan_tracks, time_sep,              &!
                        rate_asc, rate_desc, rate_lev,      &!
                        satvp_form, top_q                    !

contains
!==============================================================================
  subroutine read_airep_nml
  !--------------------------
  ! read namelist /AIREP_OBS/
  !--------------------------
    integer       :: ierr
    logical ,save :: first = .true.
    !------------------------
    ! read namelist only once
    !------------------------
    if (.not. first) return
    first = .false.
    !-------------
    ! set defaults
    !-------------
    black_t_nop   = .false. ! blacklist t-obs. if pressure is missing
    use_regnum    = .true.  ! use registration number, not flight no.
    use_tailnum   = .true.  ! use aircraft tail number, not flight no.
    require_q     = .false. ! require q for sufficient data
    reject_q_0    = .true.  ! reject zero value humidity observations
    bug_multilev  = .false. ! emulate bug reading multilevel reports?
    chk_phase     = 0       ! check for flight phase
    chk_rollangle = 0       ! check for roll-angle quality
    flev_mode     = 0       ! flight level reporting preference
    satvp_form    = 1       ! Formula for sat. vapor pres. over water
    top_q         = 250._wp ! top level for humidity observations (hPa)
    !--------------------------------------------
    ! set defaults depending on 'ga3_biasc_airep'
    !--------------------------------------------
    select case (flag_biasc_airep)
    case (-1)
      biascor_mode = BC_NOBC
      bc_fallback  = .false.
    case ( 0)
      biascor_mode = - BC_FG
      bc_fallback  = .true.
    case ( 1)
      biascor_mode = - BC_FG
      bc_fallback  = .false.
    end select
    !--------------
    ! read namelist
    !--------------
    if (dace% lpio) then
      call position_nml ('AIREP_OBS', status=ierr)
      select case (ierr)
      case (POSITIONED)
#if defined(__ibm__)
        read (nnml ,nml=AIREP_OBS, iostat=ierr)
        if (ierr/=0) call finish ('read_airep_nml',              &
                                  'ERROR in namelist /AIREP_OBS/')
#else
        read (nnml ,nml=AIREP_OBS)
#endif
      end select
      !-----------------------------------------------------
      ! adjust 'biascor_mode' depending on 'ga3_biasc_airep'
      !-----------------------------------------------------
      if (biascor_mode < 0) then
        select case (flag_biasc_airep)
        case (-1)
          biascor_mode = BC_NOBC
        case (0:1)
          biascor_mode = - biascor_mode
        end select
      endif

    endif
    if (force_level < 1000._wp) force_level = force_level * 100._wp ! hPa -> Pa
    if (top_q       < 1000._wp) top_q       = top_q       * 100._wp ! hPa -> Pa
    call p_bcast (black_t_nop,   dace% pio)
    call p_bcast (use_regnum,    dace% pio)
    call p_bcast (use_tailnum,   dace% pio)
    call p_bcast (require_q,     dace% pio)
    call p_bcast (reject_q_0,    dace% pio)
    call p_bcast (bug_multilev,  dace% pio)
    call p_bcast (chk_phase,     dace% pio)
    call p_bcast (chk_rollangle, dace% pio)
    call p_bcast (flev_mode,     dace% pio)
    call p_bcast (satvp_form,    dace% pio)
    call p_bcast (biascor_mode,  dace% pio)
    call p_bcast (t_decay,       dace% pio)
    call p_bcast (force_level,   dace% pio)
    call p_bcast (fr_land_bc,    dace% pio)
    call p_bcast (n_required,    dace% pio)
    call p_bcast (bc_fallback,   dace% pio)
    call p_bcast (scan_tracks,   dace% pio)
    call p_bcast (time_sep,      dace% pio)
    call p_bcast (rate_asc,      dace% pio)
    call p_bcast (rate_desc,     dace% pio)
    call p_bcast (rate_lev,      dace% pio)
    call p_bcast (top_q,         dace% pio)
    !---------------------------
    ! Print namelist /AIREP_OBS/
    !---------------------------
    if (dace% lpio) then
      write(6,'(a)') repeat('-',79)
      write(6,'()')
      write(6,'(a)')      ' Namelist /AIREP_OBS/:'
      write(6,'()')
      write(6,'(a,l6)'  ) ' black_t_nop       = ', black_t_nop
      write(6,'(a,l6)'  ) ' use_regnum        = ', use_regnum
      write(6,'(a,l6)'  ) ' use_tailnum       = ', use_tailnum
      write(6,'(a,l6)'  ) ' require_q         = ', require_q
      write(6,'(a,l6)'  ) ' reject_q_0        = ', reject_q_0
      write(6,'(a,l6)'  ) ' bug_multilev      = ', bug_multilev
      write(6,'(a,i6)'  ) ' chk_phase         = ', chk_phase
      write(6,'(a,i6)'  ) ' chk_rollangle     = ', chk_rollangle
      write(6,'(a,i6)'  ) ' flev_mode         = ', flev_mode
      write(6,'(a,i6)'  ) ' satvp_form        = ', satvp_form
      write(6,'(a,i6)'  ) ' biascor_mode      = ', biascor_mode
      write(6,'(a,l6)'  ) ' bc_fallback       = ', bc_fallback
      write(6,'(a,l6)'  ) ' scan_tracks       = ', scan_tracks
      write(6,'(a,i6)'  ) ' time_sep          = ', time_sep
      write(6,'(a,f6.0)') ' rate_asc          = ', rate_asc
      write(6,'(a,f6.0)') ' rate_desc         = ', rate_desc
      write(6,'(a,f6.0)') ' rate_lev          = ', rate_lev
      write(6,'(a,f6.0)') ' t_decay           = ', t_decay
      write(6,'(a,f6.0)') ' force_level (hPa) = ', force_level / 100._wp
      write(6,'(a,f6.2)') ' fr_land_bc        = ', fr_land_bc
      write(6,'(a,i6)'  ) ' n_required        = ', n_required
      write(6,'(a,f6.0)') ' top_q       (hPa) = ', top_q       / 100._wp
      write(6,'()')
    endif
    if (flev_mode < 0 .or. flev_mode > 2) &
         call finish ("read_airep_nml","invalid flev_mode")
    if (satvp_form < 1 .or. satvp_form > 2)                   &
         call finish ('read_airep_nml','satvp_form not 1 or 2')
  end subroutine read_airep_nml
!==============================================================================
  subroutine process_airep (task, spot, obs, atm, cols, xi, y, Jo, Jo_atm,&
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
    integer                :: i,j,k,n,l
    integer                :: tsk        ! task (local copy)
    real(sp)      ,pointer :: lv (:)     ! unpacked levels
    real(sp)      ,pointer :: o  (:)     ! observed value
    real(sp)               :: lv0        ! last level processed
    integer       ,pointer :: ty (:)     ! unpacked types
    real(wp)               :: p          ! pressure
    real(wp)               :: rh         ! relative humidity
    real(wp)               :: t          ! temperature
    real(wp)               :: gh         ! generalized humidity
    real(wp)               :: tv         ! virtual temperature
    real(wp)               :: dt_tv, dt_gh, drh_tv, drh_gh ! gradients
    real(wp)               :: uu, vv     ! wind components
    real(wp)               :: ff, dd     ! wind speed, direction
    real(wp)               :: ff_x, ff_y, dd_x, dd_y ! gradients
    integer                :: ifail      ! error code
    integer(i8)            :: colx       ! column fields required
    integer                :: ii, io     ! indices
    real(wp)  ,allocatable :: Hnew(:)    !
    integer                :: jtv, jrh   ! indices
    integer                :: ju, iu     ! indices
    real(wp)               :: yn(8)
    integer                :: i0,iii     !+++ index variables for NEC SX
    logical                :: change     ! argument from shrink_report
    real(wp)               :: e_ff, e_dd ! wind observation error
    real(wp)               :: e_o        ! observation error
    !======================
    ! executable statements
    !======================
    tsk = task
    if(tsk==0) return
    !----------------------------------------------------
    ! In dependence of the value of flag TASK
    ! proceed as follows:
    !  TSK_INIT       =     1 ! initialize modules
    !  TSK_READ       =     2 ! read observations
    !  TSK_SET_CHR    =     4 ! set observation characteristics
    !  TSK_SHRINK     =     8 ! release unused obs. in report
    !  TSK_SETUP_COLS =    16 ! setup columns
    !  TSK_SETUP_FUL0 =    32 ! setup interpolation space
    !  TSK_SETUP_FULL =    64 ! setup description of PSAS-space
    !  TSK_R          =   128 ! setup observational error
    !  TSK_Y          =   256 ! run forward operator
    !  TSK_H          =  1024 ! run tangent linear operator
    !  TSK_K          =  2048 ! evaluate linear operator
    !==================================================
    ! skip tasks not required for this observation type
    !==================================================
    tsk = iand (task, not (&
          TSK_READ         &! BUFR file is read in module mo_obs
    ))
    if (tsk == 0) return
    !===============
    ! initialisation
    !===============
    if (iand (TSK_INIT,tsk) /= 0) then
      call read_airep_nml
      tsk = tsk - TSK_INIT
      if (tsk == 0) return
    endif

    !=============================================
    ! TSK_SET_CHR: set observation characteristics
    !=============================================
    if (iand (TSK_SET_CHR,tsk) /= 0) then
      spot% int_type  = ITY_ICOL
      spot% cost      = 1._wp
      spot% nr        = spot% o%n
      spot% char      = CHR_ID
      if (any (obs%o% varno (spot%o%i+1:spot%o%i+spot%o%n) == VN_T   &
          .or. obs%o% varno (spot%o%i+1:spot%o%i+spot%o%n) == VN_RH))&
        spot% char    = CHR_NONL + CHR_INV + CHR_EXP
      if (any (obs%o% varno (spot%o%i+1:spot%o%i+spot%o%n) == VN_FF  &
          .or. obs%o% varno (spot%o%i+1:spot%o%i+spot%o%n) == VN_DD))&
        spot% char = CHR_NONL + CHR_EXP
      tsk = tsk - TSK_SET_CHR
      if (tsk == 0) return
    endif

    !===========
    ! PSAS-space
    !===========
    if (iand (TSK_SETUP_FUL0,tsk) /= 0) then
      call set_int_insitu (spot, obs% o)
      tsk = tsk - TSK_SETUP_FUL0
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

    !==================================================================
    ! tsk == TSK_SETUP_COLS:
    ! determine the model column indices required by AIREP observations
    ! specify the required fields (wind components)
    !==================================================================
    if (iand (TSK_SETUP_COLS,tsk) /= 0) then

      colx = 0_i8
      i = spot% o%i  +1
      n = spot% o%n+i-1
      if (any(obs% o% varno(i:n) == VN_U))      colx =        COL_UV
      if (any(obs% o% varno(i:n) == VN_T)  .or. &
          any(obs% o% varno(i:n) == VN_RH))     colx = colx + COL_TV + COL_RH
      call idx_init (     &
            spot% col% c, &! <-  column descriptor
            spot% col% h, &!  -> interpolation coefficients
            obs% o% mc,   &! <-> model column descriptors
            colx,         &! <-  fields required
            0,            &! <-  tracers required
            atm% grid,    &! <-  model grid
            spot% i_time, &! <-  time slot
            spot% w_time  )! <-  time interpolation weight

      if (spot% col% h% imc(1,1) == 0) call decr_rpt_use (spot, CHK_DOMAIN)

      tsk = tsk - TSK_SETUP_COLS
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

    !===============================================
    ! tsk == TSK_R
    ! set up R (observation error covariance matrix)
    !===============================================
    if (iand (TSK_R,tsk) /= 0) then
     if (dace% pe == obs% o% pe) then
      n = spot% o% n
      ty => obs% o% varno (spot% o% i + 1 : spot% o% i + n)
!     lv => obs% o% olev  (spot% o% i + 1 : spot% o% i + n)
      lv => obs% o% body  (spot% o% i + 1 : spot% o% i + n)% plev
      o  => obs% o% body  (spot% o% i + 1 : spot% o% i + n)% o
      !=================================
      ! setup AIREP observational errors
      !=================================
      k    = obs% R% ia (spot% o% i+1)
      i0   = spot% o% i  !+++ intermediate index variables for NEC SX
      lv0  = -huge(lv0)
      do i=1,n
        iii=i0+i     !+++ intermediate index variables for NEC SX
        obs% R% ia (iii) = k
        select case (ty(i))
        case (VN_U, VN_V, VN_FF, VN_DD)
          if (lv(i) /= lv0) then
            e_ff = obs_err (OT_AIREP, spot%hd% codetype, &
                            VN_U,  real (lv(i),wp), 0._sp)
            ff   = 0._wp
            lv0  = lv(i)
          endif
          if (ty(i) == VN_U) then
            if (obs% o% body(iii  )% o /= rvind .and. &
                obs% o% body(iii+1)% o /= rvind       ) then
              ff = sqrt (obs% o% body(iii  )% o **2 &
                 +       obs% o% body(iii+1)% o **2 )
            endif
          endif
          if (ty(i) == VN_FF) ff =       obs% o% body(iii  )% o
          if (ty(i) == VN_DD) then
            if (ff > 0._wp) then
              e_dd = r2d * e_ff / ff
            else
              e_dd = 360._wp
            endif
            e_dd = min (e_dd, 360._wp)
            obs% o% body(iii)% eo = e_dd
            obs% R% packed(k)     = e_dd ** 2
          else
            obs% o% body(iii)% eo = e_ff
            obs% R% packed(k)     = e_ff ** 2
          endif
        case (VN_T)
          e_o = obs_err (OT_AIREP, spot%hd% codetype, VN_T,  real (lv(i), wp), o(i))
          obs% o% body(iii)% eo = e_o
          obs% R% packed(k)     = e_o **2
        case (VN_RH)
          e_o = obs_err (OT_AIREP, spot%hd% codetype, VN_RH, real (lv(i), wp), o(i))
          obs% o% body(iii)% eo = e_o
          obs% R% packed(k)     = e_o **2
        case default
          call finish('process_airep (TSK_R)','invalid observation type')
        end select
        obs% R% ja (k) = iii
         k = k + 1
      end do
      obs% R% ia (spot% o% i + n + 1) = k
      !=================
      ! setup vqc bounds
      !=================
      call set_vqc_insitu (spot, obs% o)
     endif
     tsk = tsk - TSK_R
     if (tsk == 0) return
    endif

    !=========================
    ! set up H (Jakoby-matrix)
    !=========================
    if (iand (TSK_K,tsk) /= 0) then
      if (spot% pe_eval == dace% pe) then
        n = spot% o% n
        allocate (Hnew(spot% o% n * spot% i% n))
        Hnew = huge(0._wp)
        j  = 1
        ii = spot%i% i + 1
        tv = -1._wp
        ff = -1._wp
        do i = 1, spot%o% n
          l=(j - 1) * n + i
          io = spot%o% i + i
          select case (obs% o% varno(io))
          case (VN_U)
            ju    = j
            iu    = ii
            uu    = xi%x (iu)
            yn(i) = uu
            if (present(y)) y%x (io) = uu
            Hnew (l) = 1._wp
            j  = j  + 2
            ii = ii + 2
          case (VN_V)
            vv    = xi%x (iu+1)
            yn(i) = vv
            if (present(y)) y%x (io) = vv
            l = (ju) * n + i
            Hnew (l) = 1._wp
          case (VN_FF)
            if (ff < 0._wp) call fd_uv (ff, dd, uu, vv, ff_x, ff_y, dd_x, dd_y)
            yn(i) = ff
            if (present(y)) y%x (io) = ff
            Hnew ((ju - 1) * n + i) = ff_x
            Hnew ((ju    ) * n + i) = ff_y
          case (VN_DD)
            if (ff < 0._wp) call fd_uv (ff, dd, uu, vv, ff_x, ff_y, dd_x, dd_y)
            yn(i) = dd
            if (present(y)) y%x (io) = dd
            Hnew ((ju - 1) * n + i) = dd_x
            Hnew ((ju    ) * n + i) = dd_y
          case (VN_T, VN_RH)
            if (tv < 0._wp) then
              jtv = j
              jrh = j+1
              tv = xi%x (ii)
              gh = xi%x (ii+1)
              p  = exp (obs% o% lev (ii))
              call trh_tvgh (t,                            &!  -> t
                             rh,                           &!  -> rh
                             tv,                           &! <-  tv
                             gh,                           &! <-  gh
                             p,                            &! <-  p
                             ifail,                        &! -> error code
                             dt_tv, dt_gh, drh_tv, drh_gh)  !  -> gradient
              if (ifail < 0) &
                call finish ('process_airep (TSK_K)','trh_tvgh failed')
              ii = ii + 2
              j  = j  + 2
            endif
            if (obs% o% varno(io) == VN_T) then
              yn(i) = t
              if (present(y)) y%x (io) = t
              Hnew ((jtv - 1) * n + i) = dt_tv
              Hnew ((jrh - 1) * n + i) = dt_gh
            else if (obs% o% varno(io) == VN_RH) then
              yn(i) = rh
              if (present(y)) y%x (io) = rh
              Hnew ((jtv - 1) * n + i) = drh_tv
              Hnew ((jrh - 1) * n + i) = drh_gh
            endif
          case default
            write (0,*) 'process_airep, TSK_K: invalid obstype',&
                         obs% o% varno(io)
            call finish('process_airep, TSK_K','invalid observation type')
          end select
        end do
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
        deallocate (Hnew)
!       obs% xi% x (spot% i% i+1:spot% i% i+spot% i% n) = &
!            xi% x (spot% i% i+1:spot% i% i+spot% i% n)
        obs% yi% x (spot% o% i+1:spot% o% i+spot% o% n) = yn (1:spot%o% n)
      endif
      tsk = tsk - TSK_K
      if (tsk == 0) return
    endif

    !=====================================
    ! TSK_Y: evaluate observation operator
    !=====================================
    if (iand (TSK_Y,tsk) /= 0) then
      if (spot% pe_eval == dace% pe) then
        if (spot%o% n > size (yn))                                      &
             call finish ("process_airep","TSK_Y: spot%o% n > size (yn)")
        ii = spot%i% i + 1
        tv = -1._wp
        ff = -1._wp
        do i = 1, spot%o% n
          io = spot%o% i + i
          select case (obs% o% varno(io))
          case (VN_U)
            iu    = ii
            uu    = xi%x (iu)
            yn(i) = uu
            if (present(y)) y%x (io) = uu
            ii = ii + 2
          case (VN_V)
            vv    = xi%x (iu+1)
            yn(i) = vv
            if (present(y)) y%x (io) = vv
          case (VN_FF)
            if (ff < 0._wp) call fd_uv (ff, dd, uu, vv)
            yn(i) = ff
            if (present(y)) y%x (io) = ff
          case (VN_DD)
            if (ff < 0._wp) call fd_uv (ff, dd, uu, vv)
            yn(i) = dd
            if (present(y)) y%x (io) = dd
          case (VN_T, VN_RH)
            if (tv < 0._wp) then
              tv = xi%x (ii)
              gh = xi%x (ii+1)
              p  = exp (obs% o% lev (ii))
              call trh_tvgh (t,                            &!  -> t
                             rh,                           &!  -> rh
                             tv,                           &! <-  tv
                             gh,                           &! <-  gh
                             p,                            &! <-  p
                             ifail                       )  ! -> error code
              if (ifail < 0) &
                call finish ('process_airep (TSK_Y)','trh_tvgh failed')
              ii = ii + 2
            endif
            if (obs% o% varno(io) == VN_T) then
              yn(i) = t
              if (present(y)) y%x (io) = t
            else if (obs% o% varno(io) == VN_RH) then
              yn(i) = rh
              if (present(y)) y%x (io) = rh
            endif
          case (VN_VGUST,VN_TURB)
            ! Currently not implemented
            yn(i) = rvind
            if (present(y)) y%x (io) = rvind
          case default
            write (0,*) 'process_airep, TSK_Y: invalid obstype',&
                         obs% o% varno(io)
            call finish('process_airep, TSK_Y','invalid observation type')
          end select
        end do
      endif
      tsk = tsk - TSK_Y
      if (tsk == 0) return
    endif

    !===============================
    ! abort if invalid task is given
    !===============================
    call finish('process_airep','unknown task: '//tsk_name(tsk))

  end subroutine process_airep
!==============================================================================
  subroutine read_airep_bufr (bufr, spt, obs, lkeep, nkeep, cc)
  type (t_bufr) ,intent(inOUT)        :: bufr  ! BUFR record to decode
  type (t_spot) ,intent(inout)        :: spt   ! meta information to set
  type (t_obs)  ,intent(inout)        :: obs   ! observations data type to set
  logical       ,intent(out)          :: lkeep ! accept observation ?
  integer       ,intent(out)          :: nkeep ! number of accepted obsvs.
  integer       ,intent(in) ,optional :: cc    ! part of year ccyy
  !=======================================
  ! Read AIREP observations from BUFR data
  !=======================================
    !-----------------------------------------
    ! quantities derived from the BUFR message
    !-----------------------------------------
    integer         :: ival  ! value decoded from BUFR (integer)
    real(sp)        :: rval  ! value decoded from BUFR (real)
    character(len=8):: ymnem ! value decoded from BUFR (string)
!   integer         :: itype ! type of BUFR message
    integer         :: yyyy,mo,dd,hh,mi,ss ! actual time read

    !----------------
    ! index variables
    !----------------
    integer         :: is         ! sub-set    index
    integer         :: ie         ! entry in sub-set
    integer         :: id         ! descriptor index
    integer         :: i          ! spot index

    type (t_airep)  :: s, s0      ! AIREP observation read
    type (t_spot)   :: spt0, spti ! observation meta data
    logical         :: lbadqi     ! bad quality
    logical         :: lrepl      ! report replaced by correction message
    integer         :: unused     ! no. mnemonics not yet handled

    unused          = 0
    lkeep           = .false.
    nkeep           = 0
    spt0            = spt
    spt0% hd% id    = spt% hd% id - 1
    spt0% phase     = 7           ! default: missing value (phase of flight)
    s0% ff% lev_sig = 3           ! default: missing value (roll angle quality)
    !---------------------
    ! loop over data, copy
    !---------------------
    do is = 1, bufr% sec3% num_subsets
      spti             = spt0
      spti% hd% id     = spt0% hd% id + is
      spti% hd% subset = is
      ss               = 0
      s                = s0
      !----------------------------------------------------------
      ! construct 'empty' Report data structure
      ! create additional 'reports' for subsequent subset members
      !----------------------------------------------------------
      call construct_airep (s)
      if (is > 1) call check_report_0 (spti% use, spti% hd, 1)
      readdata: do ie=1,bufr% nbufdat(is)
        !---------------------------------------
        ! decode real/character datum, mnemonics
        !---------------------------------------
        ival  = bufr% ibufdat (is,ie)
        id    = bufr% idescidx(is,ie)
        if (bufr% is_char(id) == 0 .and. ival == inv_bufr) cycle
!       itype = bufr% itype(id)
        ymnem = bufr% ymnem(id)
        rval  = rvind
        if (bufr% is_char(id) == 0) then
          IF (ival /= inv_bufr) &
            rval = ival * bufr% scale(id)
        else
!         call bufr_get_character (bufr, ival, text)
        endif

        select case (ymnem)
        !------------------------------
        ! process the following entries
        !------------------------------
        case ('MCORME')
          spti% corme = ival
        case ('MLAH')   ! latitude (high accuracy)
          spti% col% c% dlat = rval
        case ('MLOH')   ! longitude (high accuracy)
          spti% col% c% dlon = rval
        case ('MLALA')  ! latitude (coarse accuracy)
          if (spti% col% c% dlat == rvind) spti% col% c% dlat = rval
        case ('MLOLO')  ! longitude (coarse accuracy)
          if (spti% col% c% dlon == rvind) spti% col% c% dlon = rval
        case ('MJJJ')   ! year
          yyyy = ival
          if(yyyy<100 .and. present(cc)) yyyy = yyyy + cc * 100
        case ('MMM')    ! month
          mo   = ival
        case ('MYY')    ! day
          dd   = ival
!       case ('MYYQ')   ! Q-BITS FOR FOLLOWING VALUE

        case ('MGG')    ! hour
          hh   = ival
!       case ('MGGQ')   ! Q-BITS FOR FOLLOWING VALUE

        case ('NGG')    ! minute
          mi   = ival
        case ('MSEC')   ! second
          ss   = ival
        case ('YXXNN')  ! aircraft flight number
          if(spti% statid=='') call bufr_get_character (bufr,ival,spti% statid)
        case ('YAIRN')  ! AIRCRAFT REGISTRATION NUMBER
          if (use_regnum) call bufr_get_character (bufr,ival,spti% statid)
        case ('YATNO')  ! AIRCRAFT TAIL  NUMBER
          if (use_tailnum) call bufr_get_character (bufr,ival,spti% statid)
        case ('MPN')    ! pressure (vert.location)
          call set_datum (s% p  ,rval ,spti% corme)
        case ('NDNDN')  ! wind direction
          call set_datum (s% dd ,rval ,spti% corme)
        case ('NDNDNQ')  ! Q-BITS FOR FOLLOWING VALUE
          call set_qbits (s% dd ,ival)
        case ('MIAA')   ! INDICATED AIRCRAFT ALTITUDE        [m]
          call set_datum (s% z  ,rval ,spti% corme)
        case ('MHHH')   ! HEIGHT OR ALTITUDE (VERT.LOCATION) [m]
          call set_datum (s% z  ,rval ,spti% corme)
        case ('NFNFN')  ! wind speed
          call set_datum (s% ff ,rval ,spti% corme)
        case ('NFNFNQ')  ! Q-BITS FOR FOLLOWING VALUE (wind speed)
          call set_qbits (s% ff ,ival)
        case ('MTN')    ! TEMPERATURE/DRY BULB TEMPERATURE
          call set_datum (s% t  ,rval ,spti% corme)
        case ('MUUU')   ! RELATIVE HUMIDITY                  [%]
          if (rval /= rvind) &
          call set_datum (s% rh ,rval / 100._sp ,spti% corme)
          if (reject_q_0 .and. s% rh% o == 0.) s% rh% qc = QC_NOUSE
        case ('MUUUO')  ! RELATIVE HUMIDITY                  [%]
          if (rval /= rvind) &
          call set_datum (s% rh ,rval / 100._sp ,spti% corme)
          if (reject_q_0 .and. s% rh% o == 0.) s% rh% qc = QC_NOUSE
        case ('MMIXR')  ! MIXING RATIO
          call set_datum (s% mr ,rval ,spti% corme)
          if (reject_q_0 .and. s% mr% o == 0.) s% mr% qc = QC_NOUSE
        case ('MMRQ')    ! MIXING RATIO QUALITY
          select case (ival)
          case (0,1,2)
          case default
            s% rh% qc = QC_NOUSE
            s% mr% qc = QC_NOUSE
          end select
        case ('MPHAI')  ! PHASE OF AIRCRAFT FLIGHT
          spti% phase = ival
        case ('MQARA')  ! AIRCRAFT ROLL ANGLE QUALITY
          s% ff% lev_sig = ival
        !-----------------------------------------------------------
        ! Meaning of some BUFR code table entries:
        !
        ! MPHAI (8004) : PHASE OF AIRCRAFT FLIGHT
        ! Code figure
        !   0 -   1   Reserved
        !      2      Unsteady (UNS)
        !      3      Level flight, routine observation (LVR)
        !      4      Level flight, highest wind encountered (LVW)
        !      5      Ascending (ASC)
        !      6      Descending (DES)
        !      7      Missing value
        !
        ! MQARA (2064) : AIRCRAFT ROLL ANGLE QUALITY
        ! Code figure
        !      0      Good
        !      1      Bad
        !      2      Reserved
        !      3      Missing value
        ! Note:       Bad is currently defined as a roll angle > 5
        !                 degrees from vertical
        !-----------------------------------------------------------
        !------------------
        ! SKY BUFR4 entries
        !------------------
        case ('NFLEV')  ! flight level                             [M]
          call set_datum (s% z  ,rval ,spti% corme)
        case ('MTDBT')  ! temperature/dry bulb temperature         [K]
          call set_datum (s% t  ,rval ,spti% corme)
        case ('NDEPF')  ! detailed phase of flight                 [CODE_TABLE]
          if ( spti% phase == 7 ) then
!           if mphai and ndepf are defined; mphai is sufficient
            spti% phase = ival
            !-------------------------------------
            !   transfer from table  8009 to 8004
            !-------------------------------------
            if( ival == 15                  ) spti% phase = 7
            if( ival == 14  .or. ival == 12  .or.      &
                ival == 10  .or. ival ==  8 ) spti% phase = 2
            if( ival == 13  .or. ival == 11 ) spti% phase = 6
            if( ival ==  9  .or. ival ==  7 ) spti% phase = 5
            if( ival ==  1  .or. ival ==  0 ) spti% phase = 2
          endif

#ifndef CHECKCODES
        case default
          !-----------------------------
          ! ignore the following entries
          !-----------------------------
#else
        !---------------------------------------------
        ! ignore the following entries in AMDAR, AIREP
        !---------------------------------------------
        case ('MOBITQ')  ! OVERALL QUALITY BITS
        case ('MADDF')   ! ASSOCIATED FIELD SIGNIFICANCE
        case ('YCCCC')   ! ICAO LOCATION INDICATOR
        case ('YDDDD')   ! SHIP OR MOBILE LAND STATION IDENTIFIER
        case ('MADDF0')  ! ASSOCIATED FIELD SIGNIFICANCE
        case ('Loop000') ! Start of Loop - 101000
        case ('MDREP')   ! DELAYED DESCRIPTOR REPLICATION FACTOR
        case ('Lcnt000') ! Loop Counter
        case ('YSUPL')   ! 008 CHARACTERS
        !------------------------------------------------------
        ! in addition ignore the following entries in ACARS_USA
        !------------------------------------------------------
        case ('NIX')     ! TYPE OF STATION
        case ('NIW')     ! TYPE OF INSTRUMENT.FOR WIND MEASUREMENT
        case ('MPOTO')   ! PRECISION OF TEMPERATURE OBSERVATION
        case ('NADRS')   ! TYPE OF AIRCRAFT DATA RELAY SYSTEM
        case ('NOSLL')   ! ORIGINAL SPECIF. OF LATITUDE/LONGITUDE
        case ('YAGRS')   ! ACARS GROUND RECEIVING STATION
        case ('NGGTI')   ! TIME INCREMENT
        case ('NGGTM')   ! DURAT.OF TIME RELAT.TO FOLLOWING VALUE
        case ('011235')  ! Unknown Descriptor
        case ('MAIV')    ! ACARS INTERPOLATED VALUES
        !-----------------------------------------------------
        ! in addition ignore the following entries in ACARS_EU
        !-----------------------------------------------------
        case ('NS1')     ! AIRCRAFT NAVIGATIONAL SYSTEM
        !-----------------------------------------------------
        ! in addition ignore the following entries in ACARS_LH
        !-----------------------------------------------------
        case ('MSREP')   ! SHORT DELAYED DESCRIPTOR REPLIC. FACTOR
        case ('MSREP0':'MSREP4')
        case ('Loop001':'Loop005')  ! Start of Loop
        case ('Lcnt001':'Lcnt005')  ! Loop Counter
        !----------------------------------------------
        ! since 2007062318 ignore the following entries
        !----------------------------------------------
        !--------------
        ! BUFR4 entries
        !--------------
        case ('MDREP0') ! DELAYED DESCRIPTOR REPLICATION FACTOR    [NUMERIC]
        case ('MB')     ! DEGREE OF TURBULENCE                     [CODE_TABLE]
        case ('NMDEWX') ! MAXIMUM DERIVED EQUIVALENT VERTICAL GUST [M/S]
        !------------------------
        ! check for unknown codes
        !------------------------
        case default
          call bufr_get_entry_texts (bufr)
          call bufr_get_entry_units (bufr)
          write(0,'(a,a,a)') ymnem, bufr% ytext(id),      &
                        '['//trim(bufr% yunit(id))//']'
          unused = unused + 1
          if (unused > 50) exit readdata
#endif
        end select
      end do readdata
      !--------------------------------------
      ! in case of unknown codes:
      ! print suspicious BUFR record and exit
      !--------------------------------------
      if (unused /=0) then
        call bufr_print_sections (bufr)
        call bufr_print_subset   (bufr, is)
        call finish ('read_airep_bufr','code(s) not implemented')
      endif

      !------------------------------------------
      ! Check flight phase, quality of roll-angle
      !------------------------------------------
      lbadqi = .false.
      select case (chk_phase)
      case (0)
         ! Accept any
      case (1)
         ! Reject 'unsteady'
         if (spti% phase    == 2) lbadqi = .true.
      end select
      select case (chk_rollangle)
      case (0)
         ! Accept any
      case (1)
         ! Reject 'bad'
         if (s% ff% lev_sig == 1) lbadqi = .true.
      case (2)
         ! Require 'good'
         if (s% ff% lev_sig /= 0) lbadqi = .true.
      end select
      if (lbadqi) call decr_rpt_use (spti, CHK_QI)
      !----------------
      ! standard checks
      !----------------
      lkeep    = .true.
      spti% ps = s% p% o
      call init_time (spti% actual_time, yyyy, mo, dd, hh, mi, ss)
      call check_report_1 (spti)
      if (spti% use% state <= STAT_DISMISS) lkeep        = .false.
      if (spti% statid     == ' '         ) spti% statid = 'AIREP'

      if (lkeep) then
        !-------------------------------------------------------
        ! check for double entries in case of correction message
        !-------------------------------------------------------
       if (spti% corme/=0) then
         write(6,*) '---------------------------------------------------------'
         write(6,*) 'AIREP correction message:'
         write(6,*) '  statid,corme,file,rec,dbkz,spot,time,lat,lon,p=',   &
           spti% statid, spti% corme, spti% hd% source, spti% hd% record,  &
           spti% hd% dbkz,obs% n_spot+1,cyyyymmddhhmm(spti%  actual_time), &
           real(spti% col% c% dlon), real(spti% col% c% dlat), real(s% p% o)
         lrepl = .false.
         write(6,*) '  replace messages:'
         do i=obs% n_spot,1,-1
           if (      obs% spot(i)% statid      == spti% statid      &
               .and. obs% spot(i)% actual_time == spti% actual_time &
               .and. obs% spot(i)% hd% dbkz    == spti% hd% dbkz    &
               .and. obs% spot(i)% corme       == spti% corme - 1   &
                                                                    ) then
             lrepl = .true.
             call decr_rpt_use (obs% spot(i), CHK_CORR)
             write(6,*) '  statid,corme,file,rec,dbkz,spot,time,lat,lon,p=',&
               obs% spot(i)% statid, obs% spot(i)% corme,                   &
               obs% spot(i)% hd% source, obs% spot(i)% hd% record,          &
               obs% spot(i)% hd% dbkz,obs% n_spot+1,                        &
               cyyyymmddhhmm(obs% spot(i)%  actual_time),                   &
               real(obs% spot(i)% col% c% dlon),                            &
               real(obs% spot(i)% col% c% dlat),                            &
               real(obs%  olev (obs% spot(i)%o%i+1))
           endif
         end do
         if (.not.lrepl) then
           lkeep = .false.
           call decr_rpt_use (spti ,CHK_CORRERR ,STAT_DISMISS)
           write(6,*) '  dismiss this and all messages:'
           do i=obs% n_spot,1,-1
             if (      obs% spot(i)% statid      == spti% statid      &
                 .and. obs% spot(i)% hd% dbkz    == spti% hd% dbkz    &
                                                                    ) then
               call decr_rpt_use (obs% spot(i) ,CHK_CORRERR ,STAT_DISMISS)
               write(6,*) '  statid,corme,file,rec,dbkz,spot,time,lat,lon,p=',&
                 obs% spot(i)% statid, obs% spot(i)% corme,                   &
                 obs% spot(i)% hd% source, obs% spot(i)% hd% record,          &
                 obs% spot(i)% hd% dbkz,obs% n_spot+1,                        &
                 cyyyymmddhhmm(obs% spot(i)%  actual_time),                   &
                 real(obs% spot(i)% col% c% dlon),                            &
                 real(obs% spot(i)% col% c% dlat),                            &
                 real(obs%  olev (obs% spot(i)%o%i+1))
             endif
           end do
         endif
       endif
       !--------------------------------------------------------
       ! insert in list if most important information is present
       !--------------------------------------------------------
       if (lkeep) call check_store_airep (s, spti, obs, lkeep)
       if (lkeep) then
         nkeep = nkeep + 1
       endif
      endif
    end do
    lkeep = nkeep > 0
  contains
  !----------------------------------------------------------------------------
    subroutine andbits (target, source, mask, shift)
    integer(i1) ,intent(inout) :: target
    integer     ,intent(in)    :: source
    integer     ,intent(in)    :: mask
    integer     ,intent(in)    :: shift
!#if defined(__SX__)
!     integer(i1) :: k              ! Temporary for workaround of bug on SX-6
!     k = int (source+mask, i1)
!     target = iand (target, ishftc (k,shift))
!#else
      target = iand (target, ishftc (int(source+mask,i1),shift))
!#endif
    end subroutine andbits
  end subroutine read_airep_bufr
!==============================================================================
  subroutine check_store_airep (airep, spot, obs, lkeep, repl)
  type(t_airep),intent(inout)        :: airep ! AIREP level information
  type(t_spot) ,intent(inout)        :: spot  ! meta data of this observation
  type(t_obs)  ,intent(inout)        :: obs   ! data of all observations
  logical      ,intent(out)          :: lkeep ! flag: keep or reject
  integer      ,intent(in) ,optional :: repl  ! observation to replace
  !==================================================
  ! check data for consistency
  ! store accepted data in observation data structure
  !==================================================

    type (t_spot) ,pointer :: s
    integer                :: no, io
    integer                :: id
    type (t_set)           :: use_t, use_q, use_uv
    logical                :: use_flev
    integer                :: stat_chkqi           ! status after QI check

    lkeep = .false.

    !---------------------------
    ! get usage flags for t,uv,q
    !---------------------------
    call get_rule (type     = spot% hd% modtype,  &! <- module      type
                   obstype  = spot% hd% obstype,  &! <- observation type
                   codetype = spot% hd% codetype, &! <- code type
                   bf_type  = iud,                &! <- no BUFR     type
                   bf_subt  = iud,                &! <- no BUFR  subtype
                   db_kz    = iud,                &! <- no Datenbankkennzahl
                   stname   = '',                 &! <- no Station Name
                   lat      = spot% col% c% dlat, &! <- latitude
                   lon      = spot% col% c% dlon, &! <- longitude
                   t        = use_t,              &! ->
                   uv       = use_uv,             &! ->
                   q        = use_q)               ! ->
    no = 0

    !-------------------------------------------
    ! derive pressure from height if not present
    !-------------------------------------------
    if (airep% p% qc /= QC_OK .and. airep% z% qc == QC_OK .and. &
        rept_use(OT_AIREP)% deriv_p /= 0                  .and. &
        airep% z% o  >  0.                                      ) then
      airep% p% o   = p_h_usstd (real(airep% z% o,wp))
      airep% p% qc  = QC_OK
      airep% p% src = SRC_DER
    endif

    !---------------------------------
    ! check for invalid pressure value
    !---------------------------------
    if (airep% p% o <= 0._sp) airep%  p% qc = QC_MISS

    !--------------------------------------
    ! Sanity check for reported temperature
    !--------------------------------------
    if (airep% t%  qc == QC_OK .and. (airep% t%  o < 100._sp .or. &
                                      airep% t%  o > 400._sp     )) then
        airep% t%  qc =  QC_CLIM  ! Temperature out of expected range
    end if
    if (airep% rh% qc == QC_OK .and.  airep% rh% o < 0._sp) then
        airep% rh% qc =  QC_MISS  ! Unphysical rel.humidity
    end if

    !---------------------------------
    ! check for presence of valid data
    !---------------------------------
    if   (spot% col% c% dlon /= invalid      .and. &
          spot% col% c% dlat /= invalid      .and. &
          airep%  p% qc      == 0)            then
      if (use_uv   % use     >  STAT_DISMISS .and. &
          airep% ff% qc      == 0 .and. airep% dd% qc == 0) then
        no = no + 2
        if (monitor_ff) no = no + 1
        if (monitor_dd) no = no + 1
      endif
      if (use_t    % use     >  STAT_DISMISS .and. &
          airep%  t% qc      == 0)                          no = no + 1
      if (use_q    % use     >  STAT_DISMISS .and. &
          airep% rh% qc      == 0            .and. &
          airep%  p% o       >= top_q      ) then
        no = no + 1
      else if (use_q% use    >  STAT_DISMISS .and. &
          airep% mr% qc      == 0            .and. &
          airep% mr% o       >= 0._sp        .and. &
          airep%  t% qc      == 0            .and. &
          airep%  p% o       >= top_q      ) then
        no = no + 1
        select case (satvp_form)
        case (1) ! Magnus-Tetens
          airep% rh% o = rhw_m (real(airep% mr% o,wp), &
                                real(airep% t % o,wp), &
                                real(airep% p % o,wp))
        case (2) ! Hardy
          airep% rh% o = rhw_m_hardy (real(airep% mr% o,wp), &
                                      real(airep% t % o,wp), &
                                      real(airep% p % o,wp))
        end select
        airep% rh% qc  = QC_OK
        airep% rh% src = SRC_DER
      endif
    endif
    !--------------------------------------------
    ! Drop bogus reports at exactly (0°N,0°E),
    ! set potentially bad ones on 0°N to passive.
    ! (Mostly ACARS data from US aircraft.)
    !--------------------------------------------
    if   (spot% col% c% dlat == 0._wp) then
      if (spot% col% c% dlon == 0._wp) then
        call decr_rpt_use (spot, CHK_NOTUSED, STAT_DISMISS)
      else
        call decr_rpt_use (spot, CHK_NOTUSED, STAT_PASSIVE)
      end if
    end if

    !--------------------------------------------------------
    ! check for presence of humidity observation if requested
    !--------------------------------------------------------
    if (require_q .and. airep% rh% qc /= QC_OK) no = 0

    !-------------------------------------
    ! Ensure consistent roll angle quality
    !-------------------------------------
    airep% dd% lev_sig = airep% ff% lev_sig
    airep% t % lev_sig = airep% ff% lev_sig
    airep% rh% lev_sig = airep% ff% lev_sig

    !--------------------------------------
    ! store data into observation data type
    !--------------------------------------
    if (no > 0) then
      if (present (repl)) then
        s => obs% spot(repl)
      else
        call new_spot (obs, 1, set_id=.true.)
        s => obs% spot(obs% n_spot)
      end if
      id    = s% id
      s     = spot
      s% id = id
!     s% int_type  = ITY_ICOL
      s% col% nlev = 1
!     s% cost      = 1._wp
!     s% char      = CHR_ID
      s% ps        = airep% p% o
      call new_obs (obs, no, s)
!     s% nr        = no
      call set_xuv (s)
      select case (flev_mode)
      case default
         use_flev = .false.                     ! never report flight level
      case (1)
         use_flev = airep% p% src == SRC_DER    ! prefer origin of p reported
      case (2)
         use_flev = airep% z% qc  == QC_OK      ! prefer flight level
      end select
      stat_chkqi  = rept_use(OT_AIREP)% use(CHK_QI)   ! status after QI check
      io = 0
      !-----
      ! wind
      !-----
      if (use_uv   % use >  STAT_DISMISS .and. &
          airep% ff% qc == 0 .and. airep% dd% qc == 0) then
        obs% varno (s%o%i+1)       = VN_U
        obs%  body (s%o%i+1)       = airep% ff
        obs%  body (s%o%i+1) %o    = airep% ff% o * (-sin (d2r * airep% dd% o))
        obs%  body (s%o%i+1) %plev = airep% p% o
        obs% varno (s%o%i+2)       = VN_V
        obs%  body (s%o%i+2)       = airep% ff
        obs%  body (s%o%i+2) %o    = airep% ff% o * (-cos (d2r * airep% dd% o))
        obs%  body (s%o%i+2) %plev = airep% p% o
        if (use_flev) then
         obs% body (s%o%i+1) %lev_typ = VN_FLEV
         obs% body (s%o%i+2) %lev_typ = VN_FLEV
         obs% olev (s%o%i+1)          = airep% z% o
         obs% olev (s%o%i+2)          = airep% z% o
        else
         obs% body (s%o%i+1) %lev_typ = VN_P
         obs% body (s%o%i+2) %lev_typ = VN_P
         obs% olev (s%o%i+1)          = airep% p% o
         obs% olev (s%o%i+2)          = airep% p% o
        end if
        if (chk_rollangle == 3 .and. airep% ff% lev_sig == 1) then
           call decr_use (obs% body(s%o%i+1)% use, stat_chkqi, check=CHK_QI)
           call decr_use (obs% body(s%o%i+2)% use, stat_chkqi, check=CHK_QI)
        end if
        io = io + 2
        if (monitor_ff) then
          obs% varno (s%o%i+io+1)          = VN_FF
          obs%  body (s%o%i+io+1)          = airep% ff
          obs%  body (s%o%i+io+1) %plev    = airep% p% o
          if (use_flev) then
           obs% body (s%o%i+io+1) %lev_typ = VN_FLEV
           obs% olev (s%o%i+io+1)          = airep% z% o
          else
           obs% body (s%o%i+io+1) %lev_typ = VN_P
           obs% olev (s%o%i+io+1)          = airep% p% o
          end if
          call decr_use (obs% body (s%o%i+io+1)% use, STAT_PASSIVE, check=CHK_NOTUSED)
          if (chk_rollangle == 3 .and. airep% ff% lev_sig == 1) then
             call decr_use (obs% body(s%o%i+io+1)% use, stat_chkqi, check=CHK_QI)
          end if
          io = io + 1
        endif
        if (monitor_dd) then
          obs% varno (s%o%i+io+1)          = VN_DD
          obs%  body (s%o%i+io+1)          = airep% dd
          obs%  body (s%o%i+io+1) %plev    = airep% p% o
          if (use_flev) then
           obs% body (s%o%i+io+1) %lev_typ = VN_FLEV
           obs% olev (s%o%i+io+1)          = airep% z% o
          else
           obs% body (s%o%i+io+1) %lev_typ = VN_P
           obs% olev (s%o%i+io+1)          = airep% p% o
          end if
          call decr_use (obs% body (s%o%i+io+1)% use, STAT_PASSIVE, check=CHK_NOTUSED)
          if (chk_rollangle == 3 .and. airep% ff% lev_sig == 1) then
             call decr_use (obs% body(s%o%i+io+1)% use, stat_chkqi, check=CHK_QI)
          end if
          io = io + 1
        endif
      endif
      !------------
      ! temperature
      !------------
      if (use_t   % use >  STAT_DISMISS .and. &
          airep% t% qc == 0) then
        obs% varno (s%o%i+io+1)          = VN_T
        obs%  body (s%o%i+io+1)          = airep% t
        obs%  body (s%o%i+io+1) %plev    = airep% p% o
        if (use_flev) then
         obs% body (s%o%i+io+1) %lev_typ = VN_FLEV
         obs% olev (s%o%i+io+1)          = airep% z% o
        else
         obs% body (s%o%i+io+1) %lev_typ = VN_P
         obs% olev (s%o%i+io+1)          = airep% p% o
        end if
!       s% char = CHR_NONL + CHR_INV + CHR_EXP
!!      s% char = CHR_LIN  + CHR_INV + CHR_EXP
        if (black_t_nop .and. airep% p% src == SRC_DER) &
          call decr_use (obs% body (s%o%i+io+1)% use, check = CHK_INSDAT)
        io = io + 1
      endif
      !---------
      ! humidity
      !---------
      if (use_q    % use >  STAT_DISMISS .and. &
          airep% rh% qc  == QC_OK        .and. &
          airep%  p% o   >= top_q              ) then
        obs% varno (s%o%i+io+1)          = VN_RH
        obs%  body (s%o%i+io+1)          = airep% rh
        obs%  body (s%o%i+io+1) %plev    = airep% p% o
        if (use_flev) then
         obs% body (s%o%i+io+1) %lev_typ = VN_FLEV
         obs% olev (s%o%i+io+1)          = airep% z% o
        else
         obs% body (s%o%i+io+1) %lev_typ = VN_P
         obs% olev (s%o%i+io+1)          = airep% p% o
        end if
!       s% char = CHR_NONL + CHR_INV + CHR_EXP
!!      s% char = CHR_LIN  + CHR_INV + CHR_EXP
        io = io + 1
      endif
!     call set_int_insitu (s, obs)
      lkeep = .true.
    else
      call decr_rpt_use (spot ,CHK_INSDAT, STAT_DISMISS)
    endif
  end subroutine check_store_airep
!==============================================================================
  pure subroutine construct_airep (s)
  type (t_airep) ,intent(out) :: s
!#if defined(__SX__)
!   ! Default initialisation does not work with sxf90 rev.360 (SX-6)
!   s = t_airep (inv_datum,inv_datum,inv_datum,inv_datum,&
!                inv_datum,inv_datum,inv_datum,inv_datum)
!#endif
    s% total %mn =''
    s% p     %mn ='p'
    s% z     %mn ='z'
    s% t     %mn ='t'
    s% ff    %mn ='ff'
    s% dd    %mn ='dd'
    s% rh    %mn ='rh'
    s% mr    %mn ='mr'
  end subroutine construct_airep
!==============================================================================
  subroutine read_airep_netcdf (ifile, i_source, obs, rec1, recl, lkeep, nkeep)
  integer       ,intent(in)           :: ifile     ! Number of netCDF file read
  integer       ,intent(inout)        :: i_source  ! number of records in source-file
  type (t_obs)  ,intent(inout)        :: obs       ! observations data type to set
  integer       ,intent(in)           :: rec1      ! first record to read
  integer       ,intent(in)           :: recl      ! last  record to read
  logical       ,intent(out)          :: lkeep     ! accept observation ?
  integer       ,intent(out)          :: nkeep     ! number of accepted obsvs.

  !==========================================================
  ! Read AIREP observations from netCDF (converted BUFR data)
  !==========================================================
  !-----------------------------------------
  ! quantities derived from the BUFR message
  !-----------------------------------------
  real    ,allocatable      :: mpn     (:,:)  !
  real    ,allocatable      :: nflev   (:,:)  ! integer in BUFR
  real    ,allocatable      :: ndndn   (:,:)  ! integer in BUFR
  real    ,allocatable      :: nfnfn   (:,:)  !
  real    ,allocatable      :: mtdbt   (:,:)  !
  real    ,allocatable      :: mmixr   (:,:)  !
  real    ,allocatable      :: muuu    (:,:)  !
  real    ,allocatable      :: muuuo   (:,:)  !
  real    ,allocatable      :: mpoto   (:  )  !
  real    ,allocatable      :: mlah0   (:,:)  !
  real    ,allocatable      :: mloh0   (:,:)  !
  integer ,allocatable      :: ndepf   (:  )  !
  integer ,allocatable      :: mphai   (:  )  !
  integer ,allocatable      :: mqara   (:,:)  !
  integer ,allocatable      :: mmrq    (:  )  !
  integer ,allocatable      :: muuuoq  (:  )  !
  integer ,allocatable      :: mdrep   (:  )  !

  !---------------------------------------
  ! NetCDF variable IDs for body   section
  !---------------------------------------
  ! variable ID's in NetCDF file for
  !    expansion  in NetCDF file BUFR- data section4

  integer              :: varid_MPN      ! PRESSURE (VERT.LOCATION)          [Pa]
  integer              :: varid_NFLEV    ! FLIGHT LEVEL(i); MIAA(i); MHHH(r) [m]
                                         ! INDICATED AIRCRAFT ALTITUDE       [m]
                                         ! HEIGHT OR ALTITUDE (VERT.LOCATION)[m]
  integer              :: varid_NDNDN    ! WIND DIRECTION                    [degree]
  integer              :: varid_NFNFN    ! WIND SPEED                        [m/s] (real)
  integer              :: varid_MTDBT    ! TEMPERATURE/DRY BULB TEMPERATURE  [K]   (real)
! integer              :: varid_MTN      ! TEMPERATURE/DRY BULB TEMPERATURE  [K]   (real) (old formats)
  integer              :: varid_MMIXR    ! MIXING RATIO                      [KG/KG] (real)
  integer              :: varid_MUUU     ! RELATIVE HUMIDITY                 [%]
  integer              :: varid_MUUUO    ! RELATIVE HUMIDITY TAMDAR          [%]
  integer              :: varid_MPOTO    ! PRECISION OF TEMPERATURE OBSERVATION [K]
  integer              :: varid_MLAH0    ! latitude(level)
  integer              :: varid_MLOH0    ! longitude(level)
  integer              :: varid_NDEPF    ! DETAILED PHASE OF FLIGHT
  integer              :: varid_MPHAI    ! PHASE OF AIRCRAFT FLIGHT
  integer              :: varid_MQARA    ! AIRCRAFT ROLL ANGLE QUALITY
  integer              :: varid_MMRQ     ! MIXING RATIO QUALITY
  integer              :: varid_MUUUOQ   ! RELATIVE HUMIDITY QUALITY TAMDAR
  integer              :: varid_MDREP    ! DELAYED DESCRIPTOR REPLICATION FACTOR [NUMERIC]

  !================
  ! local variables
  !================

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
! integer              :: ncvars         ! NetCDF number of variables  defined in this NetCDF file
! integer              :: ncatts         ! NetCDF number of attributes defined in this NetCDF file
! integer              :: unlim_dimid    ! NetCDF id for of unlimited dimension defined in this NetCDF file
  integer              :: numDims        ! NetCDF number of dimensions for individual variable in NetCDF variable
  integer              :: numAtts        ! NetCDF number of attributes for individual variable in NetCDF variable
  integer              :: dimid_level    ! NetCDF dimension id for  Levels
! integer              :: dimid_report   !                     for  Reports

  integer              :: len_level      ! number of vertical levels in NetCDF-File
  integer              :: len_report     ! number of reports in NetCDF file
! integer              :: j              ! loop index
  integer              :: ilev           ! loop index vertical levels
  integer              :: nc1            ! first  dimension for netcdf getvar in start / count
  integer              :: nc2            ! second dimension for netcdf getvar in start / count
  integer              :: start1(1)      ! start indices netcdf getvar (1-d)
  integer              :: start2(2)      ! start indices netcdf getvar (2-d)
  integer              :: count2(2)      ! counts for netcdf getvar (2-d)
  integer              :: entry1,entry   ! position in source file (subset)
  integer              :: xtype          ! NetCDF variable type

  integer ,allocatable :: ifield  (:,:)  !
  real    ,allocatable :: rfield  (:,:)  !

  integer ,allocatable :: ifield1 (:  )  !
  real    ,allocatable :: rfield1 (:  )  !

  character(len=14)    :: cdate          ! Character form of date

  integer              :: dimids (dimids_max)

  !----------------
  ! index variables
  !----------------
  integer         :: nreport        ! number of observations (default)

  integer         :: is             ! sub-set    index
! integer         :: ie             ! entry in sub-set
! integer         :: id             ! descriptor index
! integer         :: i              ! spot index

  type (t_airep)      :: s          ! AIREP observation read
  type (t_spot)       :: spt0, spti ! observation meta data
  type (t_spot) ,save :: empty      !
  type (t_use)        :: use        ! status variable

  logical             :: lbadqi     ! bad quality
! logical             :: lrepl      ! report replaced by correction message
  logical             :: lpr_airep  ! airep reports from netcdf are printed
  logical             :: lpr_extd   ! extended  printing of aireps
  integer             :: npr_extd   ! number of extended printing of aireps
! integer             :: unused     ! no. mnemonics not yet handled

! logical for meteorological variables(over all reports)
  logical         :: l_press, l_height, l_height_r, l_wind, l_mtdbt, l_mtn, &
                     l_muuu,  l_mmixr, l_mpoto , l_mlah0, l_mloh0,          &
                     l_ndepf, l_mphai, l_mqara, l_mmrq, l_muuuo, l_muuuoq
! logical for technological variables(over all reports)
  logical         :: l_mdrep
!------------------------------------------------------------------------------
  lpr_airep = .false.; if (netcdf_verb > 1) lpr_airep = .true.
  lpr_extd  = .true.
  npr_extd  =   2
  if( lpr_extd) npr_extd  = 100
  !------------------------------
  ! get default number of reports
  !------------------------------
  nreport     = recl - rec1 + 1

! unused      = 0
  lkeep       = .false.
  nkeep       = 0
  !------------------------
  ! get dimension of fields
  !------------------------
  status = nf90_inq_varid (ncid, 'NDNDN' ,  varid_NDNDN)
  status = nf90_Inquire_Variable(ncid, varid_NDNDN, ndims=numDims, dimids=dimids, natts=numAtts)

  nc1 = 0
  nc2 = 0
  len_level   = 0
  len_report  = nreport

  if(numDims == 1) then
!   dimid_report = dimids(1)
!   status = nf90_Inquire_Dimension(ncid, dimid_report, len=len_report)
    len_level = 1
    nc1       = len_report
    nc2       = 1
    start1    = [rec1]
    start2    = [rec1,1]
    count2    = [nc1,nc2]
  else if(numDims >= 2) then
    dimid_level  = dimids(1)
!   dimid_report = dimids(2)
    status = nf90_Inquire_Dimension(ncid, dimid_level,  len=len_level)
!   status = nf90_Inquire_Dimension(ncid, dimid_report, len=len_report)
    nc1       = len_level
    nc2       = len_report
    start1    =   [rec1]
    start2    = [1,rec1]
    count2    = [nc1,nc2]

!   1.case: level > 1  rep  = 1
!   2.case: level = 1  rep  > 1
!   3.case: level > 1  rep  > 1
!   if ( len_level > 1  .and. len_report == 1 ) then
!     nc1       = len_level
!     nc2       = len_report
!   else if ( len_level == 1  .and. len_report >  1 ) then
!     nc1       = len_level
!     nc2       = len_report
!   else if ( len_level >  1  .and. len_report >  1 ) then
!     nc1       = len_level
!     nc2       = len_report
!   endif

  endif

  if ( lpr_airep ) then
    write (6,'()')
    write (6,'(a,i3,a)' )   'pe=',dace% pe, ' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write (6,'(a,i3,a)' )   'pe=',dace% pe, ' BEGINN  MO_AIREP.F90    !!!!!!!!!!!!!!!!!!!!!!!'
    write (6,'(a,i3,a)' )   'pe=',dace% pe, ' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write (6,'(a,i4,/,6x,a,i4,/,6x,a,i4,/,6x,a,10i4)')   'pe=',dace% pe,            &
                            ' varid_NDNDN number of dimensions: ',numDims,          &
                            ' varid_NDNDN number of attributes: ',numAtts,          &
                            ' varid_NDNDN ids    of dimensions: ',dimids(1:numDims)
    write (6,'(a,i3,a,i8,/ , &
             & 1x  ,a,i8,/ , &
             & 1x  ,a,i10,a,i10,1x,a,i10,a,i10) ')                              &
                            'pe=',dace% pe,' Reports in BUFR File:',len_report ,&
                                           ' Levels  in BUFR File:',len_level  ,&
                                           ' nc1=',nc1,' nc2=',nc2             ,&
                                           ' rec1=',rec1,' recl=',recl
  endif
  !----------------------------------------
  ! define number of reports in source-file
  !----------------------------------------
  i_source = len_report

  allocate (mpn     (len_level,len_report))
  allocate (nflev   (len_level,len_report))
  allocate (ndndn   (len_level,len_report))
  allocate (nfnfn   (len_level,len_report))
  allocate (mtdbt   (len_level,len_report))
  allocate (mmixr   (len_level,len_report))
  allocate (muuu    (len_level,len_report))
  allocate (muuuo   (len_level,len_report))
  allocate (muuuoq  (          len_report))

  allocate (ndepf   (          len_report))
  allocate (mphai   (          len_report))
  allocate (mqara   (len_level,len_report))
  allocate (mmrq    (          len_report))
  allocate (mpoto   (          len_report))
  allocate (mdrep   (          len_report))
  allocate (mlah0   (len_level,len_report))
  allocate (mloh0   (len_level,len_report))

  allocate (ifield  (len_level,len_report))
  allocate (rfield  (len_level,len_report))
  allocate (ifield1 (          len_report))
  allocate (rfield1 (          len_report))

  !-------------------------------
  ! get meteorological information
  !-------------------------------

  !---------
  ! pressure
  !---------
  mpn = rvind
  l_press = .FALSE.
  status = nf90_inq_varid (ncid, 'MPN'  ,  varid_mpn )
  if (status == nf90_noerr) then
    l_press = .TRUE.

    status = nf90_get_var (ncid, varid_mpn, rfield(:,:), start=start2, count=count2)

!   check for missing values
    if ( all( rfield == rmissing ))  then
       l_press = .FALSE.
    else
       mpn = rfield
       where ( rfield == rmissing ) mpn = rvind
    endif
  endif

  !-----------
  ! height [m]
  !-----------
  nflev = rvind
  l_height   = .FALSE.
  l_height_r = .FALSE.
  status = nf90_inq_varid (ncid, 'NFLEV',  varid_nflev)
  if (status == nf90_noerr) then
    l_height = .TRUE.
  else
    status = nf90_inq_varid (ncid, 'MIAA' ,  varid_nflev)
    if (status == nf90_noerr) then
      l_height = .TRUE.
    else
      status = nf90_inq_varid (ncid, 'MHHH' ,  varid_nflev)
      if (status == nf90_noerr) then
        l_height_r = .TRUE.
      endif
    endif
  endif
  if ( l_height ) then

! ---------------------
! standard
! ---------------------
! also for 2- dim
!!if(numDims == 1) then
    status = nf90_get_var (ncid, varid_nflev, ifield(:,:), start=start2, count=count2)
!!endif
    if ( all( ifield == imissing ))  then
       l_height = .FALSE.
    else
       where ( ifield == imissing ) ifield = iud
       rfield = ifield
       where ( rfield == rud ) rfield = rvind
       nflev = rfield
    endif
  else if ( l_height_r ) then
    status = nf90_get_var (ncid, varid_nflev, rfield(:,:), start=start2, count=count2)
    if ( all( rfield == rmissing ))  then
       l_height = .FALSE.
    else
       where ( rfield == rmissing ) rfield = rvind
       nflev = rfield
    endif
  endif

  !---------
  ! wind
  !---------
  ndndn = rvind
  nfnfn = rvind
  l_wind = .FALSE.
  status = nf90_inq_varid (ncid, 'NDNDN',  varid_ndndn)
  if (status == nf90_noerr) then
    l_wind   = .TRUE.
    status = nf90_get_var (ncid, varid_ndndn, ifield(:,:), start=start2, count=count2)
    if ( all( ifield == imissing ))  then
       l_wind = .FALSE.
    else
       where ( ifield == imissing ) ifield = iud
       rfield = ifield
       where ( rfield == rud ) rfield = rvind
       ndndn = rfield
    endif
    if ( l_wind ) then
      status = nf90_inq_varid (ncid, 'NFNFN',  varid_nfnfn)
      if (status == nf90_noerr) then
      status = nf90_get_var (ncid, varid_nfnfn, rfield(:,:), start=start2, count=count2)
        if ( all( rfield == rmissing ))  then
           l_wind = .FALSE.
        else
           where ( rfield == rmissing ) rfield = rvind
           nfnfn = rfield
        endif
      endif
    endif
  endif

  !------------
  ! temperature
  !------------
  mtdbt = rvind
  l_mtdbt = .FALSE.
  l_mtn   = .FALSE.

  status = nf90_inq_varid (ncid, 'MTDBT'  ,  varid_mtdbt )
  if (status == nf90_noerr) then
    l_mtdbt = .TRUE.
    status = nf90_get_var (ncid, varid_mtdbt, rfield(:,:), start=start2, count=count2)
!   check for missing values
    if ( all( rfield == rmissing ))  then
       l_mtdbt = .FALSE.
    else
       mtdbt = rfield
       where ( rfield == rmissing ) mtdbt = rvind
    endif

  else
     status = nf90_inq_varid (ncid, 'MTN'    ,  varid_mtdbt   )
     if (status == nf90_noerr) then
       l_mtn   = .TRUE.
       status = nf90_get_var (ncid, varid_mtdbt, rfield(:,:), start=start2, count=count2)
   !   check for missing values
       if ( all( rfield == rmissing ))  then
          l_mtn   = .FALSE.
       else
          mtdbt = rfield
          where ( rfield == rmissing ) mtdbt = rvind
      endif
    endif

  endif

  !-------------
  ! rel.humidity
  !-------------
  muuu = rvind
  l_muuu = .FALSE.
  status = nf90_inq_varid (ncid, 'MUUU'  ,  varid_muuu )
  if (status == nf90_noerr) then
    l_muuu = .TRUE.
    status = nf90_Inquire_Variable (ncid, varid_muuu, xtype=xtype)
    if (status /= nf90_noerr) &
         call finish('read_airep_netcdf',&
                     'Inquire_Variable(MUUU) failed for '//trim (source(ifile)% file))
    !---------------------------------------
    ! type-specific check for missing values
    !---------------------------------------
    select case (xtype)
    case (NF90_INT)
       if (netcdf_verb > 0)                                       &
       write(6,'(a,i3,3a)') 'pe=',dace% pe,' read_airep_netcdf: ',&
            trim (source(ifile)% file)," : MUUU is integer"
       status = nf90_get_var (ncid, varid_muuu, ifield(:,:), start=start2, count=count2)
       if ( all ( ifield == imissing )) then
          l_muuu = .FALSE.
       else
          where ( ifield == imissing ) ifield = iud
          rfield = ifield
          where ( rfield == rud      ) rfield = rvind
          muuu = rfield
       endif
    case (NF90_FLOAT)
       if (netcdf_verb > 0)                                       &
       write(6,'(a,i3,3a)') 'pe=',dace% pe,' read_airep_netcdf: ',&
            trim (source(ifile)% file)," : MUUU is float"
       status = nf90_get_var (ncid, varid_muuu, rfield(:,:), start=start2, count=count2)
       if ( all ( rfield == rmissing )) then
          l_muuu = .FALSE.
       else
          where ( rfield == rmissing ) rfield = rvind
          muuu = rfield
       endif
    case default
       write(6,*) 'pe=',dace% pe,"read_airep_netcdf: xtype =", xtype
       call finish('read_airep_netcdf','invalid type for MUUU in '//trim (source(ifile)% file))
    end select
  endif

  !---------------------
  ! rel. humidity TAMDAR
  !---------------------
  muuuo = rvind
  l_muuuo = .FALSE.
  status = nf90_inq_varid (ncid, 'MUUUO'  ,  varid_muuuo )
  if (status == nf90_noerr) then
    l_muuuo = .TRUE.
    status = nf90_Inquire_Variable (ncid, varid_muuuo, xtype=xtype)
    if (status /= nf90_noerr) &
         call finish('read_airep_netcdf',&
                     'Inquire_Variable(MUUUO) failed for '//trim (source(ifile)% file))
    !---------------------------------------
    ! type-specific check for missing values
    !---------------------------------------
    select case (xtype)
    case (NF90_INT)
       if (netcdf_verb > 0)                                       &
       write(6,'(a,i3,3a)') 'pe=',dace% pe,' read_airep_netcdf: ',&
            trim (source(ifile)% file)," : MUUUO is integer"
       status = nf90_get_var (ncid, varid_muuuo, ifield(:,:), start=start2, count=count2)
       if ( all ( ifield == imissing )) then
          l_muuuo = .FALSE.
       else
          where ( ifield == imissing ) ifield = iud
          rfield = ifield
          where ( rfield == rud      ) rfield = rvind
          muuuo = rfield
       endif
    case (NF90_FLOAT)
       if (netcdf_verb > 0)                                       &
       write(6,'(a,i3,3a)') 'pe=',dace% pe,' read_airep_netcdf: ',&
            trim (source(ifile)% file)," : MUUUO is float"
       status = nf90_get_var (ncid, varid_muuuo, rfield(:,:), start=start2, count=count2)
       if ( all ( rfield == rmissing )) then
          l_muuuo = .FALSE.
       else
          where ( rfield == rmissing ) rfield = rvind
          muuuo = rfield
       endif
    case default
       write(6,*) 'pe=',dace% pe,"read_airep_netcdf: xtype =", xtype
       call finish('read_airep_netcdf','invalid type for MUUUO in '//trim (source(ifile)% file))
    end select
  endif

  !-------------
  ! mix.ratio
  !-------------
  mmixr = rvind
  l_mmixr = .FALSE.
  status = nf90_inq_varid (ncid, 'MMIXR'  ,  varid_mmixr )
  if (status == nf90_noerr) then
    l_mmixr = .TRUE.
    status = nf90_get_var (ncid, varid_mmixr  , rfield(:,:), start=start2, count=count2)
!   check for missing values
    if ( all( rfield == rmissing ))  then
       l_mmixr = .FALSE.
    else
       mmixr = rfield
       where ( rfield == rmissing ) mmixr = rvind
    endif
  endif

  !-------------------------------------
  ! precision of temperature observation
  !-------------------------------------
  mpoto   = rvind
  l_mpoto = .FALSE.
  status = nf90_inq_varid (ncid, 'MPOTO'  ,  varid_mpoto )
  if (status == nf90_noerr) then
    l_mpoto = .TRUE.
    status = nf90_get_var (ncid, varid_mpoto, rfield1(:), start=start1)
!   check for missing values
    if ( all( rfield1 == rmissing ))  then
       l_mpoto = .FALSE.
    else
       mpoto = rfield1
       where ( rfield1 == rmissing ) mpoto = rvind
    endif
  endif

  !-------------------------------------
  ! level dependent hor.coordinates
  !-------------------------------------
  l_mdrep = .false.
  if (len_level > 1) then
    status = nf90_inq_varid (ncid, 'MDREP', varid_mdrep)
    if (status == nf90_noerr) then
       l_mdrep = .true.
       status = nf90_get_var (ncid, varid_mdrep, ifield1(:), start=start1)
!      check for missing values
       if (all (ifield1 == imissing)) then
          l_mdrep = .false.
       else
          mdrep = ifield1
          where (ifield1 == imissing) mdrep = len_level     ! or 0 or 1?
          mdrep = min (mdrep, len_level)
       endif
    end if
  end if
  if (.not. l_mdrep) mdrep = len_level

  mlah0 = rvind
  l_mlah0 = .FALSE.
  l_mloh0 = .FALSE.

  if (l_mdrep) then
    status = nf90_inq_varid (ncid, 'MLAH0', varid_mlah0)
    if (status == nf90_noerr) then
      l_mlah0 = .TRUE.
      status = nf90_get_var (ncid, varid_mlah0, rfield(:,:), start=start2, count=count2)
!     check for missing values
      if ( all( rfield == rmissing ))  then
         l_mlah0 = .FALSE.
      else
         mlah0 = rfield
         where ( rfield == rmissing ) mlah0 = rvind
      endif
    endif
  endif
  if (l_mlah0) then
    status = nf90_inq_varid (ncid, 'MLOH0', varid_mloh0)
    if (status == nf90_noerr) then
      l_mloh0 = .TRUE.
      status = nf90_get_var (ncid, varid_mloh0, rfield(:,:), start=start2, count=count2)
!     check for missing values
      if ( all( rfield == rmissing ))  then
         l_mloh0 = .FALSE.
      else
         mloh0 = rfield
         where ( rfield == rmissing ) mloh0 = rvind
      endif
    else
      l_mlah0 = .FALSE.
      l_mloh0 = .FALSE.
    endif
  endif
  ! Pretend mlah0/mloh0 were missing
  if (bug_multilev) l_mlah0 = .FALSE.

  !-------------------------------------
  ! phase of aircraft flight
  ! aircraft roll angle quality
  ! detailed phase of aircraft flight
  !-------------------------------------
! missing in table
  ndepf = 15
  mphai = 7
  mqara = 3
  l_ndepf = .FALSE.
  l_mphai = .FALSE.

  status = nf90_inq_varid (ncid, 'NDEPF'  ,  varid_ndepf )
  if (status == nf90_noerr) then
    l_ndepf = .TRUE.
    status = nf90_get_var (ncid, varid_ndepf, ifield1(:), start=start1)
!   check for missing values
    if ( all( ifield1 == imissing ))  then
       l_ndepf = .FALSE.
    else
       ndepf = ifield1
!      where ( ifield1 == imissing ) ndepf = -1
       where ( ifield1 == imissing ) ndepf = 15
    endif

!-------------------------------------
!   transfer from table  8009 to 8004
!-------------------------------------
    mphai = ndepf
    where( ndepf == 15                   ) mphai = 7
    where( ndepf == 14  .or. ndepf == 12  .or.      &
           ndepf == 10  .or. ndepf ==  8 ) mphai = 2
    where( ndepf == 13  .or. ndepf == 11 ) mphai = 6
    where( ndepf ==  9  .or. ndepf ==  7 ) mphai = 5
    where( ndepf ==  1  .or. ndepf ==  0 ) mphai = 2
  else
    status = nf90_inq_varid (ncid, 'MPHAI'  ,  varid_mphai )
    if (status == nf90_noerr) then
       l_mphai = .TRUE.
       status = nf90_get_var (ncid, varid_mphai, ifield1(:), start=start1)
!      check for missing values
       if ( all( ifield1 == imissing ))  then
          l_mphai = .FALSE.
       else
          mphai = ifield1
!         where ( ifield1 == imissing ) mphai = -1
          where ( ifield1 == imissing ) mphai =  7
       endif
    endif
  endif

  status = nf90_inq_varid (ncid, 'MQARA'  ,  varid_mqara )
  if (status == nf90_noerr) then
    l_mqara = .TRUE.
    status = nf90_get_var (ncid, varid_mqara, ifield(:,:), start=start2, count=count2)
!   check for missing values
    if ( all( ifield == imissing ))  then
       l_mqara = .FALSE.
    else
       mqara = ifield
!      where ( ifield == imissing ) mqara = -1
       where ( ifield == imissing ) mqara =  3
    endif
  endif

  !-------------------------------------
  ! quality of humidity
  !-------------------------------------
  l_mmrq = .FALSE.
  mmrq   = -1
  status = nf90_inq_varid (ncid, 'MMRQ'  ,  varid_mmrq )
  if (status == nf90_noerr) then
    l_mmrq  = .TRUE.
    status = nf90_get_var (ncid, varid_MMRQ, ifield1(:), start=start1)
!   check for missing values
    if ( all( ifield1 == imissing ))  then
       l_mmrq = .FALSE.
    else
       mmrq = ifield1
       where ( ifield1 == imissing ) mmrq = -1
    endif
  endif

  !-------------------------------------
  ! quality of humidity TAMDAR
  !-------------------------------------
  l_muuuoq = .FALSE.
  muuuoq   = -1
  status = nf90_inq_varid (ncid, 'MUUUOQ'  ,  varid_muuuoq )
  if (status == nf90_noerr) then
    l_muuuoq  = .TRUE.
    status = nf90_get_var (ncid, varid_MUUUOQ, ifield1(:), start=start1)
!   check for missing values
    if ( all( ifield1 == imissing ))  then
       l_muuuoq = .FALSE.
    else
       muuuoq = ifield1
       where ( ifield1 == imissing ) muuuoq = -1
    endif
  endif

  !----------------------------------
  ! list of defined variables in file
  !----------------------------------
  if (netcdf_verb > 0) then
    write (6,'(a,i3,1x,a,a,/,5x,17l8)') 'pe=',dace% pe ,                 &
     'l_press l_height l_height_r l_wind l_mtdbt l_mtn  l_muuu l_muuuo ',&
     'l_mmixr',                                                          &
      l_press,l_height,l_height_r,l_wind,l_mtdbt,l_mtn, l_muuu,l_muuuo,  &
      l_mmixr
    write (6,'(a,i3,1x,a,a,/,5x,17l8)') 'pe=',dace% pe ,                 &
     'l_mpoto  l_mdrep  l_mlah0 l_mloh0 l_ndepf l_mphai l_mqara l_mmrq ',&
     'l_muuuoq',                                                         &
      l_mpoto, l_mdrep, l_mlah0,l_mloh0,l_ndepf,l_mphai,l_mqara,l_mmrq,  &
      l_muuuoq
  end if

  !-------------------------------
  ! preset total number of reports
  !-------------------------------
  entry   = sum (source(1:ifile-1)% entries) + rec1 - 1

  !------------------
  ! loop over reports
  !------------------
! original: is number of subsets
! netcdf  : is number of report
  do is = 1, len_report

    entry1  = entry    + 1
!   entry   = entry    + nsubset
    entry   = entry    + 1

   !------------------
   ! loop over levels
   !------------------
!  do ilev = 1,len_level
   do ilev = 1,mdrep(is)
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
! test
    if (lpr_airep) then
       write (6,'( 8(a, i6 ))')                                    &
            'pe=',dace% pe,' ifile ',ifile,' is=',is,' ilev ',ilev,&
            ' report_subt ',report_subt,' bufr_type ',bufr_type,   &
            ' bufr_subtype ',bufr_subtype, ' centre ',centre
    end if
! test

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
!   head% modtype     = rept_char(obstype)% mod
    head% modtype     = AIREP
    head% buf_type    = bufr_type
    head% buf_subtype = bufr_subtype
    head% codetype    = obsid% codetype
    head% time        = stime(is)
    head% db_time     = db_time(is)
    head% idbk        = report_subti
    head% source      = ifile
    head% record      = is + rec1 - 1
    head% id          = entry1
    head% center      = s1cent(is)
    head% subcenter   = s1cents(is)

    if ( lpr_airep .and. (is + ilev) < npr_extd ) then
      write (6,'()')
      write (6,'( 9(a16, i6  ,/),   &
                & 2(a16, a   ,/),   &
                & 6(a16, i6  ,/) )' )                       &
         'pe='         ,dace% pe,                           &
         'head is='    ,is,                                 &
         'ilev='       ,ilev,                               &
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

    !--------------------------------------------
    ! perform simple generic check on report type
    !--------------------------------------------
    call check_report_0 (use, head, 1)
    if (use% state <= STAT_DISMISS) then
      lkeep = .false.
! test
      if (lpr_airep) then
        print *,'pe=',dace% pe,' after check_report_0 STAT_DISMISS lkeep=',lkeep
        print *,'pe=',dace% pe,' after check_report_0 STAT_DISMISS  use =',use
        print *,'pe=',dace% pe,' after check_report_0 head% obstype     =',head% obstype
      end if
! test

      cycle
    endif

    !------------------
    ! create new report
    !------------------
    spt0        = empty
    spt0% use   = use
    spt0% hd    = head
    spt0% phase = 31

    spti = spt0

    ! Treat level index like a subset, even if it technically isn't
    if (len_level > 1) spti% hd% subset = ilev
    !----------------------------------------------------------
    ! construct 'empty' Report data structure
    ! create additional 'reports' for subsequent subset members
    !----------------------------------------------------------
    call construct_airep (s)

    !------------------------------
    ! process the following entries
    !------------------------------
    spti% corme        = max ( s1updat(is), 0)
    if ( len_level == 1 ) then
      spti% col% c% dlat =   mlah(is)
      spti% col% c% dlon =   mloh(is)
    else if ( len_level >  1 .and. .not.l_mlah0) then
      spti% col% c% dlat =   mlah(is)
      spti% col% c% dlon =   mloh(is)
    else if ( len_level >  1 .and.  l_mlah0) then
      spti% col% c% dlat =   mlah0(ilev,is)
      spti% col% c% dlon =   mloh0(ilev,is)
    endif
    spti% actual_time  =     obs_time(is)
    spti% statid       =     ystidn(is)

    if ( lpr_airep .and. (is + ilev) < npr_extd ) then
      cdate = cyyyymmddhhmmss (spti% actual_time)
      write (6,'()')
      write (6,'(   a21, i6  ,a, /, &
                &   a21, i6  ,   /, &
                & 2(a21,f8.3 ,   /),&
                &   a21,a,1x,a,  / )' )                    &
            'pe=',dace% pe,       '  spti ',               &
            'spti% corme        = ', spti% corme        ,  &
            'spti% col% c% dlat = ', spti% col% c% dlat ,  &
            'spti% col% c% dlon = ', spti% col% c% dlon ,  &
            'spti% actual_time  = ', cdate(1:8), cdate(9:)
    endif

!  pressure (vert.location)
    call set_datum (s% p  ,mpn(ilev,is)  ,spti% corme)

!  height or altitude (vert.location) [m]
    call set_datum (s% z  ,nflev(ilev,is),spti% corme)

!  wind direction/ speed
    call set_datum (s% dd ,ndndn(ilev,is) ,spti% corme)
    call set_datum (s% ff ,nfnfn(ilev,is) ,spti% corme)

!  temperature/dry bulb temperature
    call set_datum (s% t  ,mtdbt(ilev,is) ,spti% corme)

!  relative humidity / mixing ratio
    if (l_muuu) then
      if (muuu(ilev,is) /= rvind) then
          muuu(ilev,is) = muuu(ilev,is) / 100.0_sp      ! % -> (0..1)
      endif
      call set_datum (s% rh ,muuu(ilev,is) ,spti% corme)
      if (reject_q_0 .and. s% rh% o == 0.) s% rh% qc = QC_NOUSE
    end if

    if (l_mmixr) then
      call set_datum (s% mr ,mmixr(ilev,is) ,spti% corme)
      if (reject_q_0 .and. s% mr% o == 0.) s% mr% qc = QC_NOUSE
    end if

    if (l_mmrq) then
      select case (mmrq (is))
      case (-1,0,1,2)         ! -1 for missing
      case default
        s% rh% qc = QC_NOUSE
        s% mr% qc = QC_NOUSE
      end select
    end if

    !---------------------
    ! TAMDAR rel. humidity
    !---------------------
    if (l_muuuo) then
      if (muuuo(ilev,is) /= rvind) then
          muuuo(ilev,is) = muuuo(ilev,is) / 100.0_sp    ! % -> (0..1)
      endif
      call set_datum (s% rh ,muuuo(ilev,is) ,spti% corme)
      if (reject_q_0 .and. s% rh% o == 0.) s% rh% qc = QC_NOUSE

      select case (muuuoq (is))
      case (100)              ! 100 good data quality
      case default
        s% rh% qc = QC_NOUSE
      end select
    end if

!  phase of aircraft flight
    spti% phase    = mphai(     is)

!  aircraft roll angle quality
    s% ff% lev_sig = mqara(ilev,is)
    !------------------------------------------
    ! Check flight phase, quality of roll-angle
    !------------------------------------------
    lbadqi = .false.
    select case (chk_phase)
    case (0)
       ! Accept any
    case (1)
       ! Reject 'unsteady'
       if (spti% phase    == 2) lbadqi = .true.
    end select
    select case (chk_rollangle)
    case (0)
       ! Accept any
    case (1)
       ! Reject 'bad'
       if (s% ff% lev_sig == 1) lbadqi = .true.
    case (2)
       ! Require 'good'
       if (s% ff% lev_sig /= 0) lbadqi = .true.
    end select
    if (lbadqi) call decr_rpt_use (spti, CHK_QI)

    if ( lpr_airep .and. (is + ilev) < npr_extd ) then
      write (6,'(   a15, i3  ,a, /,                               &
         &        5(a15,f10.3,   /),                              &
         &        2(a15,e10.3,   /),                              &
         &          a15, i8  ,   / ,                              &
         &          a15, i8  ,a20, / )' )                         &
                  'pe='    ,    dace% pe,'  s    ',               &
                  's% p='  ,    s% p  %o,                         &
                  's% z='  ,    s% z  %o,                         &
                  's% dd=' ,    s% dd %o,                         &
                  's% ff=' ,    s% ff %o,                         &
                  's% t='  ,    s% t  %o,                         &
                  's% rh=' ,    s% rh %o,                         &
                  's% mr=' ,    s% mr %o,                         &
                  'spti% phase=', spti% phase,                    &
                  's% ff% lev_sig=',s% ff% lev_sig, ' aircraft roll angle quality'
    endif

    !-------------------------------------------------------------------------------------
    ! Meaning of some BUFR code table entries:
    !
    ! MPHAI (8004) : PHASE OF AIRCRAFT FLIGHT
    ! Code figure
    !   0 -   1   Reserved
    !      2      Unsteady (UNS)
    !      3      Level flight, routine observation (LVR)
    !      4      Level flight, highest wind encountered (LVW)
    !      5      Ascending (ASC)
    !      6      Descending (DES)
    !      7      Missing value
    !
    ! MQARA (2064) : AIRCRAFT ROLL ANGLE QUALITY
    ! Code figure
    !      0      Good
    !      1      Bad
    !      2      Reserved
    !      3      Missing value
    ! Note:       Bad is currently defined as a roll angle > 5
    !                 degrees from vertical
    !-------------------------------------------------------------------------------------
    ! NDEPF (8009) : DETAILED PHASE OF AIRCRAFT FLIGHT
    !
    ! Code figure
    !      0 -> 2 Level flight, routine observation, unsteady
    !      1 -> 2 Level flight, highest wind encountered, unsteady
    !      2 -> 2 Unsteady (UNS)
    !      3 -> 3 Level flight, routine observation (LVR)
    !      4 -> 4 Level flight, highest wind encountered (LVW)
    !      5 -> 5 Ascending (ASC)
    !      6 -> 6 Descending (DES)
    !      7 -> 5 Ascending, observation intervals selected by time increments
    !      8 -> 2 Ascending, observation intervals selected by time increments, unsteady
    !      9 -> 5 Ascending, observation intervals selected by pressure increments
    !     10 -> 2 Ascending, observation intervals selected by pressure increments, unsteady
    !     11 -> 6 Descending, observation intervals selected by time increments
    !     12 -> 2 Descending, observation intervals selected by time increments, unsteady
    !     13 -> 6 Descending, observation intervals selected by pressure increments
    !     14 -> 2 Descending, observation intervals selected by pressure increments, unsteady
    !     15 -> 7 Missing value
    !-------------------------------------------------------------------------------------
    ! MMRQI (033026) : Moisture quality
    !
    ! Code figure
    !      0     Normal operations      measurement mode
    !      1     Normal operations non-measurement mode
    !      2     small rh
    !      3     Humidity element is wet
    !      4     Humidity element contaminated
    !      5     Heater fail
    !      6     Heater fail and wet/contaminated humidity element
    !      7     At least one of the input parameters used in the calculation of mixing ratio is invalid
    !      8     Numeric error
    !      9     Sensor not installed
    !     10-62  Reserved
    !     63     Missing value
    !-------------------------------------------------------------------------------------

#ifndef CHECKCODES
!       case default
          !-----------------------------
          ! ignore the following entries
          !-----------------------------
#else
!       case ('MGGQ')    ! Q-BITS FOR FOLLOWING VALUE
        !---------------------------------------------
        ! ignore the following entries in AMDAR, AIREP
        !---------------------------------------------
!       case ('MOBITQ')  ! OVERALL QUALITY BITS
!       case ('MADDF')   ! ASSOCIATED FIELD SIGNIFICANCE
!       case ('YCCCC')   ! ICAO LOCATION INDICATOR
!       case ('YDDDD')   ! SHIP OR MOBILE LAND STATION IDENTIFIER
!       case ('MADDF0')  ! ASSOCIATED FIELD SIGNIFICANCE
!       case ('Loop000') ! Start of Loop - 101000
!       case ('MDREP')   ! DELAYED DESCRIPTOR REPLICATION FACTOR
!       case ('Lcnt000') ! Loop Counter
!       case ('YSUPL')   ! 008 CHARACTERS
        !------------------------------------------------------
        ! in addition ignore the following entries in ACARS_USA
        !------------------------------------------------------
!       case ('NIX')     ! TYPE OF STATION
!       case ('NIW')     ! TYPE OF INSTRUMENT.FOR WIND MEASUREMENT
!       case ('MPOTO')   ! PRECISION OF TEMPERATURE OBSERVATION
!       case ('NADRS')   ! TYPE OF AIRCRAFT DATA RELAY SYSTEM
!       case ('NOSLL')   ! ORIGINAL SPECIF. OF LATITUDE/LONGITUDE
!       case ('YAGRS')   ! ACARS GROUND RECEIVING STATION
!       case ('MMRQ')    ! MIXING RATIO QUALITY
!       case ('NGGTI')   ! TIME INCREMENT
!       case ('NGGTM')   ! DURAT.OF TIME RELAT.TO FOLLOWING VALUE
!       case ('011235')  ! Unknown Descriptor
!       case ('MAIV')    ! ACARS INTERPOLATED VALUES
!       case ('NDNDNQ')  ! Q-BITS FOR FOLLOWING VALUE
!       case ('MYYQ')    ! Q-BITS FOR FOLLOWING VALUE
        !-----------------------------------------------------
        ! in addition ignore the following entries in ACARS_EU
        !-----------------------------------------------------
!       case ('NS1')     ! AIRCRAFT NAVIGATIONAL SYSTEM
        !-----------------------------------------------------
        ! in addition ignore the following entries in ACARS_LH
        !-----------------------------------------------------
!       case ('MSREP')   ! SHORT DELAYED DESCRIPTOR REPLIC. FACTOR
!       case ('MSREP0':'MSREP4')
!       case ('Loop001':'Loop005')  ! Start of Loop
!       case ('Lcnt001':'Lcnt005')  ! Loop Counter
!       !----------------------------------------------
        ! since 2007062318 ignore the following entries
        !----------------------------------------------
!       case ('NFNFNQ')  ! Q-BITS FOR FOLLOWING VALUE (wind speed)
        !------------------------
        ! check for unknown codes
        !------------------------
!       case default
!         call bufr_get_entry_texts (bufr)
!         call bufr_get_entry_units (bufr)
!         write(0,'(a,a,a)') ymnem, bufr% ytext(id),      &
!                       '['//trim(bufr% yunit(id))//']'
!         unused = unused + 1
!         if (unused > 50) exit readdata
#endif
!       end select
!     end do readdata
      !--------------------------------
      ! abort if unkown codes are found
      !--------------------------------
!     if (unused /= 0) call finish('read_airep_netcdf','code(s) not implemented')

    !----------------
    ! standard checks
    !----------------
    lkeep    = .true.
    spti% ps = s% p% o
    call check_report_1 (spti)
    if (spti% use% state <= STAT_DISMISS) lkeep        = .false.
    if (spti% statid     == ' '         ) spti% statid = 'AIREP'

    if (lkeep) then
      !-------------------------------------------------------
      ! check for double entries in case of correction message
      !-------------------------------------------------------
!      if (spti% corme/=0) then
!        write(6,*) '---------------------------------------------------------'
!        write(6,*) 'AIREP correction message:'
!        write(6,*) '  statid,corme,file,rec,dbkz,spot,time,lat,lon,p=',   &
!          spti% statid, spti% corme, spti% hd% source, spti% hd% record,  &
!          spti% hd% dbkz,obs% n_spot+1,cyyyymmddhhmm(spti%  actual_time), &
!          real(spti% col% c% dlon), real(spti% col% c% dlat), real(s% p% o)
!        lrepl = .false.
!        write(6,*) '  replace messages:'
!        do i=obs% n_spot,1,-1
!          if (      obs% spot(i)% statid      == spti% statid      &
!              .and. obs% spot(i)% actual_time == spti% actual_time &
!              .and. obs% spot(i)% hd% dbkz    == spti% hd% dbkz    &
!              .and. obs% spot(i)% corme       == spti% corme - 1   &
!                                                                   ) then
!            lrepl = .true.
!            call decr_rpt_use (obs% spot(i), CHK_CORR)
!            write(6,*) '  statid,corme,file,rec,dbkz,spot,time,lat,lon,p=',&
!              obs% spot(i)% statid, obs% spot(i)% corme,                   &
!              obs% spot(i)% hd% source, obs% spot(i)% hd% record,          &
!              obs% spot(i)% hd% dbkz,obs% n_spot+1,                        &
!              cyyyymmddhhmm(obs% spot(i)%  actual_time),                   &
!              real(obs% spot(i)% col% c% dlon),                            &
!              real(obs% spot(i)% col% c% dlat),                            &
!              real(obs%  olev (obs% spot(i)%o%i+1))
!          endif
!        end do
!        if (.not.lrepl) then
!          lkeep = .false.
!          call decr_rpt_use (spti ,CHK_CORRERR ,STAT_DISMISS)
!          write(6,*) '  dismiss this and all messages:'
!          do i=obs% n_spot,1,-1
!            if (      obs% spot(i)% statid      == spti% statid      &
!                .and. obs% spot(i)% hd% dbkz    == spti% hd% dbkz    &
!                                                                   ) then
!              call decr_rpt_use (obs% spot(i) ,CHK_CORRERR ,STAT_DISMISS)
!              write(6,*) '  statid,corme,file,rec,dbkz,spot,time,lat,lon,p=',&
!                obs% spot(i)% statid, obs% spot(i)% corme,                   &
!                obs% spot(i)% hd% source, obs% spot(i)% hd% record,          &
!                obs% spot(i)% hd% dbkz,obs% n_spot+1,                        &
!                cyyyymmddhhmm(obs% spot(i)%  actual_time),                   &
!                real(obs% spot(i)% col% c% dlon),                            &
!                real(obs% spot(i)% col% c% dlat),                            &
!                real(obs%  olev (obs% spot(i)%o%i+1))
!            endif
!          end do
!        endif
!      endif
       !--------------------------------------------------------
       ! insert in list if most important information is present
       !--------------------------------------------------------
       if (lkeep) call check_store_airep (s, spti, obs, lkeep)
       if (lkeep) then
         nkeep = nkeep + 1
       endif
    endif
   !------------------------
   ! end loop over levels
   !------------------------
   end do
  !------------------------
  ! end loop over reports
  !------------------------
  end do
  lkeep = nkeep > 0

  !------------
  ! deallocate
  !------------

  if (allocated(mpn     )) deallocate (mpn     )
  if (allocated(nflev   )) deallocate (nflev   )
  if (allocated(ndndn   )) deallocate (ndndn   )
  if (allocated(nfnfn   )) deallocate (nfnfn   )
  if (allocated(mtdbt   )) deallocate (mtdbt   )
  if (allocated(mmixr   )) deallocate (mmixr   )
  if (allocated(muuu    )) deallocate (muuu    )
  if (allocated(muuuo   )) deallocate (muuuo   )
  if (allocated(muuuoq  )) deallocate (muuuoq  )

  if (allocated(ndepf   )) deallocate (ndepf   )
  if (allocated(mphai   )) deallocate (mphai   )
  if (allocated(mqara   )) deallocate (mqara   )
  if (allocated(mmrq    )) deallocate (mmrq    )
  if (allocated(mpoto   )) deallocate (mpoto   )
  if (allocated(mdrep   )) deallocate (mdrep   )
  if (allocated(mlah0   )) deallocate (mlah0   )
  if (allocated(mloh0   )) deallocate (mloh0   )

  if (allocated(ifield  )) deallocate (ifield  )
  if (allocated(rfield  )) deallocate (rfield  )
  if (allocated(ifield1 )) deallocate (ifield1 )
  if (allocated(rfield1 )) deallocate (rfield1 )

  if ( lpr_airep ) then
     write (6,'()')
     write (6,'(a,i3,a)' )'pe=',dace% pe,' ENDE    MO_AIREP.F90 xxxxxxxxxxxxxxxxxxxxxxx'
     write (6,'()')
  endif
  contains
  !----------------------------------------------------------------------------
    subroutine andbits (target, source, mask, shift)
    integer(i1) ,intent(inout) :: target
    integer     ,intent(in)    :: source
    integer     ,intent(in)    :: mask
    integer     ,intent(in)    :: shift
!#if defined(__SX__)
!     integer(i1) :: k              ! Temporary for workaround of bug on SX-6
!     k = int (source+mask, i1)
!     target = iand (target, ishftc (k,shift))
!#else
      target = iand (target, ishftc (int(source+mask,i1),shift))
!#endif
    end subroutine andbits
  end subroutine read_airep_netcdf
!==============================================================================
end module mo_airep
